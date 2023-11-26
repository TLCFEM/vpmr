/*******************************************************************************
 * Copyright (C) 2023 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include <iomanip>
#include <mpfr/exprtk_mpfr_adaptor.hpp>
#include <mutex>
#include <string>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_for.h>
#include "BigInt.hpp"
#include "Cache.hpp"
#include "Eigen/Core"
#include "Eigen/MatrixFunctions"
#include "Gauss.hpp"
#ifdef HAVE_WINDOWS_H
#include <Windows.h>
#endif

using namespace mpfr;

bool OUTPUT_W = false;              // output W
bool OUTPUT_EIG = false;            // output eigenvalues
int NC = 4;                         // maximum exponent
int N = 10;                         // number of terms
int DIGIT = 512;                    // number of digits
int QUAD_ORDER = 500;               // number of quadrature points
int SCALE = 6;                      // scaling
mpreal TOL = mpreal("1E-8", DIGIT); // tolerance
std::string KERNEL = "exp(-t*t/4)"; // default kernel

BigInt comb_impl(unsigned, unsigned);

// cached function for combinatorial number calculation
Cache<BigInt, unsigned, unsigned> comb(comb_impl);

// recursive implementation of combinatorial number calculation
BigInt comb_impl(unsigned n, unsigned k) {
    if(k == 0 || k == n) return {1};
    if(k == 1 || k == n - 1) return {n};
    if(k > n) return {0};
    if(2 * k > n) k = n - k;

    return comb(n - 1, k - 1) + comb(n - 1, k);
}

// kernel function
class Expression {
    mpreal time;

    exprtk::symbol_table<mpreal> symbol_table;
    exprtk::expression<mpreal> expression;
    exprtk::parser<mpreal> parser;

    std::mutex eval_mutex;

public:
    bool compile() {
        symbol_table.add_variable("t", time);
        symbol_table.add_constants();
        expression.register_symbol_table(symbol_table);
        return parser.compile(KERNEL, expression);
    }

    mpreal value(const mpreal& t) {
        std::scoped_lock lock(eval_mutex);
        time = t;
        return expression.value();
    }
};

// integrand in eq. 2.2 --- K(t)cos(jt)
class Integrand {
    static Expression* expression;
    static tbb::concurrent_unordered_map<int, mpreal> value;

    const int j;

public:
    explicit Integrand(const int J)
        : j(J) {}

    static void set_expression(Expression* E) { expression = E; }

    mpreal operator()(const int i, const mpreal& r) const {
        if(value.find(i) == value.end())
            value.insert(std::make_pair(i, expression->value(-NC * log((mpreal(1, DIGIT) + cos(r)) >> 1))));
        return value[i] * cos(j * r);
    }
};

Expression* Integrand::expression = nullptr;
tbb::concurrent_unordered_map<int, mpreal> Integrand::value(QUAD_ORDER);

using mat = Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic>;
using vec = Eigen::Matrix<mpreal, Eigen::Dynamic, 1>;
using cx_vec = Eigen::Matrix<std::complex<mpreal>, Eigen::Dynamic, 1>;

// step 2
mat lyap(const vec& A, mat&& C) {
    tbb::parallel_for(0, static_cast<int>(A.size()), [&](const int i) {
        for(auto j = 0; j < A.size(); ++j) C(i, j) /= (A(i) + A(j));
    });

    return std::move(C);
}

long pos(const vec& A) {
    auto sum = mpreal(0, DIGIT);
    long P = 0;
    for(auto I = A.size() - 1; I >= 0; --I) {
        sum += A(I);
        if(sum > TOL) {
            P = I;
            break;
        }
    }
    return P + 1;
}

mat lyap_rhs(const vec& M) {
    mat A(M.size(), M.size());
    tbb::parallel_for(0, static_cast<int>(M.size()), [&](const int I) {
        A(I, I) = -M(I) * M(I);
        for(auto J = I + 1; J < M.size(); ++J) A(I, J) = A(J, I) = -M(I) * M(J);
    });
    return A;
}

std::tuple<cx_vec, cx_vec> model_reduction(const vec& A, const vec& B, const vec& C) {
    // account for 4^j here separately
    vec W = A;
    auto ONE = mpreal(1, DIGIT);
    for(auto I = 0; I < W.size(); ++I) W(I) = (ONE <<= 2);

    // step 3
    std::cout << "[2/6] Solving Lyapunov equation...\n";
    const mat S = lyap(A, lyap_rhs(B)).llt().matrixL();
    const mat L = lyap(A, lyap_rhs(C)).llt().matrixL();

    // step 4
    std::cout << "[3/6] Solving SVD...\n";
    const auto SVD = Eigen::BDCSVD<mat, Eigen::ComputeFullU>(S.transpose() * W.asDiagonal() * L);
    const auto& U = SVD.matrixU();
    const auto& SIG = SVD.singularValues();
    if(OUTPUT_EIG) {
        std::cout << "SIGMA = \n";
        for(auto I = 0; I < SIG.size(); ++I) std::cout << SIG(I).toString() << '\n';
    }

    // step 7
    const auto P = pos(SIG);
    std::cout << "[4/6] Transforming (P=" << P << ")...\n";
    if(P == SIG.size())
        std::cout << "WARNING: No singular value is smaller than the given tolerance.\n";

    // step 5
    auto SS = SIG;
    for(auto I = 0; I < SS.size(); ++I)
        if(SS(I) > std::numeric_limits<mpreal>::epsilon())
            SS(I) = pow(SS(I), -.5);
        else {
            std::cout << "WARNING: Need to increase digits.\n";
            SS(I) = pow(std::numeric_limits<mpreal>::epsilon(), -.5);
        }

    // step 5
    const mat T = S * U * SS.asDiagonal();
    const auto LUT = T.fullPivLu();

    // step 6
    std::cout << "[5/6] Solving eigen decomposition...\n";
    const auto EIGEN = Eigen::EigenSolver<mat>(LUT.solve(A.asDiagonal() * T).block(0, 0, P, P));
    // step 8
    const auto X = EIGEN.eigenvectors();

    // step 9
    const cx_vec BB = LUT.solve(B).head(P);
    const cx_vec CC = ((C.cwiseProduct(W).transpose() * T).head(P) * X).transpose();

    return std::make_tuple(X.fullPivLu().solve(BB).cwiseProduct(CC), EIGEN.eigenvalues());
}

mpreal weight(const Quadrature& quad, const int j) {
    auto k_hat = [&quad](const int l) { return quad.integrate(Integrand(l)); };
    auto sgn_bit = [](const int l) { return l % 2 == 0 ? 1 : -1; };

    if(0 == j) {
        mpreal result = k_hat(0) / 2;
        for(auto l = 1; l <= N; ++l) result += sgn_bit(l) * k_hat(l);
        for(auto l = 1; l < N; ++l) result += sgn_bit(N + l) * mpreal(N - l, DIGIT) / mpreal(N, DIGIT) * k_hat(N + l);
        return result;
    }

    mpreal result{0, DIGIT};

    if(j > N) {
        for(auto l = j - N; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal((N - l) * (N + l), DIGIT) / mpreal(N * (N + l + j), DIGIT) *
                mpreal(comb(N + l + j, N + l - j).to_string(), DIGIT) * k_hat(N + l);
    }
    else {
        for(auto l = j; l <= N; ++l)
            result += sgn_bit(l - j) * mpreal(l, DIGIT) / mpreal(l + j, DIGIT) *
                mpreal(comb(l + j, l - j).to_string(), DIGIT) * k_hat(l);
        for(auto l = 1; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal((N - l) * (N + l), DIGIT) / mpreal(N * (N + l + j), DIGIT) *
                mpreal(comb(N + l + j, N + l - j).to_string(), DIGIT) * k_hat(N + l);
    }

    // unlike the original paper, we do not multiply 4^j here
    // to avoid multiplication of large complex numbers and large powers of 4
    return result;
}

std::vector<long> sort_index(const cx_vec& v) {
    std::vector<long> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(), [&v](const long i, const long j) {
        const auto value_i = norm(v(i)).toDouble(), value_j = norm(v(j)).toDouble();
        if(value_i > value_j) return true;
        if(value_i < value_j) return false;
        return v(i).imag() > 0;
    });

    return idx;
}

std::tuple<cx_vec, cx_vec> vpmr() {
    const auto quad = Quadrature(QUAD_ORDER);

    vec W = vec::Zero(2 * N);
    for(auto I = 0; I < W.size();) {
        W(I) = weight(quad, I);
        std::cout << "\r[1/6] Computing weights... [" << ++I << '/' << W.size() << ']' << std::flush;
    }
    std::cout << std::showpos << std::setprecision(16) << '\n';

    if(OUTPUT_W) {
        std::cout << "W = \n";
        for(auto I = 0; I < W.size(); ++I) std::cout << W(I).toString() << '\n';
    }

    // step 1
    vec A = vec::Zero(W.size() - 1), B = vec::Zero(A.size()), C = vec::Zero(A.size());
    for(decltype(W.size()) I = 0, J = 1; I < A.size(); ++I, ++J) {
        A(I) = -mpreal(J, DIGIT) / NC;
        B(I) = sqrt(abs(W(J)));
        C(I) = sgn(W(J)) * B(I);
    }

    auto [M, S] = model_reduction(A, B, C);

    auto ID = sort_index(M);
    for(int I = int(ID.size()) - 1; I >= 0; --I)
        if(norm(M(ID[I])) < TOL)
            ID.erase(ID.begin() + I);
        else
            break;

    std::cout << "[6/6] Done.\n\n";

    if(abs(W(0)) < TOL) return std::make_tuple(M(ID), -S(ID));

    cx_vec MM = cx_vec::Zero(long(ID.size()) + 1), SS = cx_vec::Zero(long(ID.size()) + 1);
    MM(0) = W(0);
    SS(0) = mpreal(0, DIGIT);
    MM.tail(ID.size()) = M(ID);
    SS.tail(ID.size()) = -S(ID);

    return std::make_tuple(MM, SS);
}

int print_helper() {
    std::cout << "--> \xF0\x9F\xA5\xB7 VPMR C++ Implementation <--\n\n";
    std::cout << "Usage: vpmr [options]\n\n";
    std::cout << "Options:\n\n";
    std::cout << "   -n <int>     number of terms (default: 10)\n";
    std::cout << "   -d <int>     number of precision bits (default: 512)\n";
    std::cout << "   -q <int>     quadrature order (default: 500)\n";
    std::cout << "   -m <int>     precision multiplier (default: 6)\n";
    std::cout << "   -nc <int>    controls the maximum exponent (default: 4)\n";
    std::cout << "   -e <float>   tolerance (default: 1E-8)\n";
    std::cout << "   -k <string>  file name of kernel function (default: exp(-t^2/4))\n";
    std::cout << "   -s           print singular values\n";
    std::cout << "   -w           print weights\n";
    std::cout << "   -h           print this help message\n";
    return 0;
}

#ifndef PYVPMR
int main(const int argc, const char** argv) {
#ifdef HAVE_WINDOWS_H
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif

    const auto start = std::chrono::high_resolution_clock::now();

    bool has_digit = false;
    for(auto I = 1; I < argc; ++I) {
        const auto token = std::string(argv[I]);

        if(token == "-nc")
            NC = std::stoi(argv[++I]);
        else if(token == "-n")
            N = std::stoi(argv[++I]);
        else if(token == "-d") {
            DIGIT = std::stoi(argv[++I]);
            has_digit = true;
        }
        else if(token == "-q")
            QUAD_ORDER = std::stoi(argv[++I]);
        else if(token == "-m")
            SCALE = std::stoi(argv[++I]);
        else if(token == "-e")
            TOL = mpreal(argv[++I]);
        else if(token == "-w")
            OUTPUT_W = true;
        else if(token == "-s")
            OUTPUT_EIG = true;
        else if(token == "-h")
            return print_helper();
        else if(token == "-k") {
            const std::string file_name(argv[++I]);
            std::ifstream file(file_name);
            if(!file.is_open()) {
                std::cerr << "Cannot open file: " << file_name << std::endl;
                return 1;
            }
            KERNEL = "";
            std::string line;
            while(std::getline(file, line)) KERNEL += line;
            file.close();
        }
        else {
            std::cerr << "Unknown option: " << token << std::endl;
            return 1;
        }
    }

    // check size
    BigInt comb_max = comb(2 * N, N);
    int comb_digit = 0;
    while((comb_max /= 2) > 0) ++comb_digit;
    comb_digit = std::max(5 * N, SCALE * comb_digit);

    if(!has_digit)
        DIGIT = comb_digit;
    else if(comb_digit >= DIGIT) {
        std::cout << "WARNING: Too few digits to hold combinatorial number, resetting digits to " << comb_digit << ".\n";
        DIGIT = comb_digit;
    }

    // set precision
    mpreal::set_default_prec(DIGIT);

    TOL.setPrecision(DIGIT);

    MP_PI = const_pi(2 * DIGIT);
    MP_PI_HALF = MP_PI / 2;

    TOL /= 2;

    Expression kernel;
    if(!kernel.compile()) {
        std::cerr << "Cannot compile kernel function: " << KERNEL << ".\n";
        return 1;
    }
    Integrand::set_expression(&kernel);

    std::cout << std::scientific << std::setprecision(4);

    std::cout << "Using the following parameters:\n";
    std::cout << "        nc = " << NC << ".\n";
    std::cout << "         n = " << N << ".\n";
    std::cout << "     order = " << QUAD_ORDER << ".\n";
    std::cout << " precision = " << DIGIT << ".\n";
    std::cout << " tolerance = " << (2 * TOL).toDouble() << ".\n";
    std::cout << "    kernel = " << KERNEL << ".\n\n";

    try {
        // run VPMR algorithm
        const auto [M, S] = vpmr();

        std::cout << "M = \n";
        for(auto I = 0; I < M.size(); ++I) std::cout << M(I).real().toDouble() << M(I).imag().toDouble() << "j\n";

        std::cout << "S = \n";
        for(auto I = 0; I < S.size(); ++I) std::cout << S(I).real().toDouble() << S(I).imag().toDouble() << "j\n";
    }
    catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << std::noshowpos << "\nRunning time: " << duration << " ms.\n";

    return 0;
}
#else
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>> vpmr_wrapper(
    const int n, const int d, const int q, const int m, const int nc, const double e, const std::string& k) {
    N = n;
    DIGIT = d;
    QUAD_ORDER = q;
    SCALE = m;
    NC = nc;
    TOL = mpreal(e);
    if(!k.empty()) KERNEL = k;

    // check size
    BigInt comb_max = comb(2 * N, N);
    int comb_digit = 0;
    while((comb_max /= 2) > 0) ++comb_digit;
    comb_digit = std::max(5 * N, SCALE * comb_digit);

    if(0 == d)
        DIGIT = comb_digit;
    else if(comb_digit >= DIGIT) {
        std::cout << "WARNING: Too few digits to hold combinatorial number, resetting digits to " << comb_digit << ".\n";
        DIGIT = comb_digit;
    }

    // set precision
    mpreal::set_default_prec(DIGIT);

    TOL.setPrecision(DIGIT);

    MP_PI = const_pi(2 * DIGIT);
    MP_PI_HALF = MP_PI / 2;

    TOL /= 2;

    std::vector<std::complex<double>> mm, ss;

    Expression kernel;
    if(!kernel.compile()) {
        std::cerr << "Cannot compile kernel function: " << KERNEL << ".\n";
        return {mm, ss};
    }
    Integrand::set_expression(&kernel);

    const auto [M, S] = vpmr();

    for(const auto& I : M) mm.emplace_back(I.real().toDouble(), I.imag().toDouble());
    for(const auto& I : S) ss.emplace_back(I.real().toDouble(), I.imag().toDouble());

    return {mm, ss};
}

PYBIND11_MODULE(_pyvpmr, m) {
    m.doc() = "The VPMR Algorithm";

    m.def(
        "vpmr", &vpmr_wrapper, pybind11::call_guard<pybind11::gil_scoped_release>(),
        pybind11::kw_only(), pybind11::arg("n") = 10, pybind11::arg("d") = 0, pybind11::arg("q") = 500, pybind11::arg("m") = 6, pybind11::arg("nc") = 4, pybind11::arg("e") = 1E-8, pybind11::arg("k") = "",
        "The VPMR Algorithm.\n\n"
        ":param n: number of terms (default: 10)\n"
        ":param d: number of precision bits (default: 512)\n"
        ":param q: quadrature order (default: 500)\n"
        ":param m: precision multiplier (default: 6)\n"
        ":param nc: maximum exponent (default: 4)\n"
        ":param e: tolerance (default: 1E-8)\n"
        ":param k: kernel function (default: exp(-t^2/4))\n"
        ":return: M, S\n");
}
#endif

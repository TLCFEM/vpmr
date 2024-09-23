/*******************************************************************************
 * Copyright (C) 2024 Theodore Chang
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

#include "VPMR.h"
#include <mpfr/exprtk_mpfr_adaptor.hpp>
#include <mutex>
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

Config config;

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
        return parser.compile(config.kernel, expression);
    }

    mpreal value(const mpreal& t) {
        std::scoped_lock lock(eval_mutex);
        time = t;
        return expression.value();
    }
};

Expression kernel{};

// integrand in eq. 2.2 --- K(t)cos(jt)
class Integrand {
    static tbb::concurrent_unordered_map<int, mpreal> value;

    const int j;

public:
    explicit Integrand(const int J)
        : j(J) {}

    mpreal operator()(const int i, const mpreal& r) const {
        if(value.find(i) == value.end()) value.insert(std::make_pair(i, kernel.value(-config.max_exponent * log((mpreal(1, config.precision_bits) + cos(r)) >> 1))));
        return value[i] * cos(j * r);
    }
};

tbb::concurrent_unordered_map<int, mpreal> Integrand::value(100);

using mat = Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic>;
using vec = Eigen::Matrix<mpreal, Eigen::Dynamic, 1>;
using cx_vec = Eigen::Matrix<std::complex<mpreal>, Eigen::Dynamic, 1>;

// step 2
mat lyap(const vec& A, mat&& C) {
    tbb::parallel_for(0, static_cast<int>(A.size()), [&](const int i) { for(auto j = 0; j < A.size(); ++j) C(i, j) /= (A(i) + A(j)); });

    return std::move(C);
}

long pos(const vec& A) {
    auto sum = mpreal(0, config.precision_bits);
    long P = 0;
    for(auto I = A.size() - 1; I >= 0; --I) {
        sum += A(I);
        if(sum > config.tolerance) {
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
    auto ONE = mpreal(1, config.precision_bits);
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
    if(config.print_singular_value) {
        std::cout << "SIGMA = \n";
        for(auto I = 0; I < SIG.size(); ++I) std::cout << SIG(I).toString() << '\n';
    }

    // step 7
    const auto P = pos(SIG);
    std::cout << "[4/6] Transforming (P=" << P << ")...\n";
    if(P == SIG.size()) std::cout << "WARNING: No singular value is smaller than the given tolerance.\n";

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

    const auto N = config.max_terms;
    const auto D = config.precision_bits;

    if(0 == j) {
        mpreal result = k_hat(0) / 2;
        for(auto l = 1; l <= N; ++l) result += sgn_bit(l) * k_hat(l);
        for(auto l = 1; l < N; ++l) result += sgn_bit(N + l) * mpreal(N - l, D) / mpreal(N, D) * k_hat(N + l);
        return result;
    }

    mpreal result{0, D};

    if(j > N) {
        for(auto l = j - N; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal((N - l) * (N + l), D) / mpreal(N * (N + l + j), D) *
                mpreal(comb(N + l + j, N + l - j).to_string(), D) * k_hat(N + l);
    }
    else {
        for(auto l = j; l <= N; ++l)
            result += sgn_bit(l - j) * mpreal(l, D) / mpreal(l + j, D) *
                mpreal(comb(l + j, l - j).to_string(), D) * k_hat(l);
        for(auto l = 1; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal((N - l) * (N + l), D) / mpreal(N * (N + l + j), D) *
                mpreal(comb(N + l + j, N + l - j).to_string(), D) * k_hat(N + l);
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
    const auto quad = Quadrature(config.quadrature_order);

    vec W = vec::Zero(2 * config.max_terms);
    for(auto I = 0; I < W.size();) {
        W(I) = weight(quad, I);
        std::cout << "\r[1/6] Computing weights... [" << ++I << '/' << W.size() << ']' << std::flush;
    }
    std::cout << std::showpos << std::setprecision(16) << '\n';

    if(config.print_weight) {
        std::cout << "W = \n";
        for(auto I = 0; I < W.size(); ++I) std::cout << W(I).toString() << '\n';
    }

    // step 1
    vec A = vec::Zero(W.size() - 1), B = vec::Zero(A.size()), C = vec::Zero(A.size());
    for(decltype(W.size()) I = 0, J = 1; I < A.size(); ++I, ++J) {
        A(I) = -mpreal(J, config.precision_bits) / config.max_exponent;
        B(I) = sqrt(abs(W(J)));
        C(I) = sgn(W(J)) * B(I);
    }

    auto [M, S] = model_reduction(A, B, C);

    auto ID = sort_index(M);
    for(int I = static_cast<int>(ID.size()) - 1; I >= 0; --I)
        if(norm(M(ID[I])) < config.tolerance)
            ID.erase(ID.begin() + I);
        else
            break;

    std::cout << "[6/6] Done.\n\n";

    if(abs(W(0)) < config.tolerance) return std::make_tuple(M(ID), -S(ID));

    cx_vec MM = cx_vec::Zero(static_cast<long>(ID.size()) + 1), SS = cx_vec::Zero(static_cast<long>(ID.size()) + 1);
    MM(0) = W(0);
    SS(0) = mpreal(0, config.precision_bits);
    MM.tail(ID.size()) = M(ID);
    SS.tail(ID.size()) = -S(ID);

    return std::make_tuple(MM, SS);
}

int print_helper() {
    std::cout << "--> \xF0\x9F\xA5\xB7 VPMR C++ Implementation <--\n\n";
    std::cout << "Compiled on: " << __DATE__ << "\n\n";
    std::cout << "Usage: vpmr [options]\n\n";
    std::cout << "Options:\n\n";
    std::cout << "    -n, --max-terms             <int>     number of terms (default: 10)\n";
    std::cout << "    -c, --max-exponent          <int>     maximum exponent (default: 4)\n";
    std::cout << "    -d, --precision-bits        <int>     number of precision bits (default: 512)\n";
    std::cout << "    -q, --quadrature-order      <int>     quadrature order (default: 500)\n";
    std::cout << "    -m, --precision-multiplier  <float>   precision multiplier (default: 1.5)\n";
    std::cout << "    -e, --tolerance             <float>   tolerance (default: 1E-8)\n";
    std::cout << "    -k, --kernel                <string>  file name of kernel function (default uses: exp(-t^2/4))\n";
    std::cout << "    -s, --singular-values                 print singular values\n";
    std::cout << "    -w, --weights                         print weights\n";
    std::cout << "    -h, --help                            print this help message\n\n";
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
        if(const auto token = std::string(argv[I]); token == "-n" || token == "--max-terms")
            config.max_terms = std::max(1, std::stoi(argv[++I]));
        else if(token == "-c" || token == "--max-exponent")
            config.max_exponent = std::max(1, std::stoi(argv[++I]));
        else if(token == "-d" || token == "--precision-bits") {
            config.precision_bits = std::max(1, std::stoi(argv[++I]));
            has_digit = true;
        }
        else if(token == "-q" || token == "--quadrature-order")
            config.quadrature_order = std::max(1, std::stoi(argv[++I]));
        else if(token == "-m" || token == "--precision-multiplier")
            config.precision_multiplier = std::max(1.5, std::stod(argv[++I]));
        else if(token == "-e" || token == "--tolerance")
            config.tolerance = mpreal(argv[++I]);
        else if(token == "-k" || token == "--kernel") {
            const std::string file_name(argv[++I]);
            std::ifstream file(file_name);
            if(!file.is_open()) {
                std::cerr << "Cannot open file: " << file_name << std::endl;
                return 1;
            }
            config.kernel = "";
            std::string line;
            while(std::getline(file, line)) config.kernel += line;
            file.close();
        }
        else if(token == "-w" || token == "--weights")
            config.print_weight = true;
        else if(token == "-s" || token == "--singular-values")
            config.print_singular_value = true;
        else if(token == "-h" || token == "--help")
            return print_helper();
        else {
            std::cerr << "Unknown option: " << token << std::endl;
            return 1;
        }
    }

    // check size
    BigInt comb_max = comb(4 * config.max_terms, 2 * config.max_terms);
    int comb_digit = 1;
    while((comb_max /= 2) > 0) ++comb_digit;
    comb_digit = std::max(64, static_cast<int>(std::max(1.5, config.precision_multiplier) * static_cast<double>(comb_digit + 4 * config.max_terms)));

    if(!has_digit)
        config.precision_bits = comb_digit;
    else if(comb_digit >= config.precision_bits) {
        std::cout << "WARNING: Too few digits to hold combinatorial number, resetting digits to " << comb_digit << ".\n";
        config.precision_bits = comb_digit;
    }

    // set precision
    mpreal::set_default_prec(config.precision_bits);

    config.tolerance.setPrecision(config.precision_bits);

    MP_PI = const_pi(config.precision_bits);
    MP_PI_HALF = MP_PI / 2;

    config.tolerance /= 2;

    if(!kernel.compile()) {
        std::cerr << "Cannot compile kernel function: " << config.kernel << ".\n";
        return 1;
    }

    config.print();

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
    const int n, const int c, const int d, const int q, const double m, const double e, const std::string& k) {
    config.max_terms = std::max(1, n);
    config.max_exponent = std::max(1, c);
    config.precision_bits = std::max(1, d);
    config.quadrature_order = std::max(1, q);
    config.precision_multiplier = std::max(1.5, m);
    config.tolerance = mpreal(e);
    if(!k.empty()) config.kernel = k;

    // check size
    BigInt comb_max = comb(4 * config.max_terms, 2 * config.max_terms);
    int comb_digit = 1;
    while((comb_max /= 2) > 0) ++comb_digit;
    comb_digit = std::max(64, static_cast<int>(std::max(1.5, config.precision_multiplier) * static_cast<double>(comb_digit + 4 * config.max_terms)));

    if(0 == d)
        config.precision_bits = comb_digit;
    else if(comb_digit >= config.precision_bits) {
        std::cout << "WARNING: Too few digits to hold combinatorial number, resetting digits to " << comb_digit << ".\n";
        config.precision_bits = comb_digit;
    }

    // set precision
    mpreal::set_default_prec(config.precision_bits);

    config.tolerance.setPrecision(config.precision_bits);

    MP_PI = const_pi(config.precision_bits);
    MP_PI_HALF = MP_PI / 2;

    config.tolerance /= 2;

    std::vector<std::complex<double>> mm, ss;

    if(!kernel.compile()) {
        std::cerr << "Cannot compile kernel function: " << config.kernel << ".\n";
        return {mm, ss};
    }

    const auto [M, S] = vpmr();

    for(const auto& I : M) mm.emplace_back(I.real().toDouble(), I.imag().toDouble());
    for(const auto& I : S) ss.emplace_back(I.real().toDouble(), I.imag().toDouble());

    return {mm, ss};
}

PYBIND11_MODULE(_pyvpmr, m) {
    m.doc() = "The VPMR Algorithm";

    m.def(
        "vpmr", &vpmr_wrapper, pybind11::call_guard<pybind11::gil_scoped_release>(),
        pybind11::kw_only(), pybind11::arg("n") = 10, pybind11::arg("c") = 4, pybind11::arg("d") = 0, pybind11::arg("q") = 500, pybind11::arg("m") = 1.5, pybind11::arg("e") = 1E-8, pybind11::arg("k") = "",
        "The VPMR Algorithm.\n\n"
        ":param n: number of terms (default: 10)\n"
        ":param c: maximum exponent (default: 4)\n"
        ":param d: number of precision bits (default: 512)\n"
        ":param q: quadrature order (default: 500)\n"
        ":param m: precision multiplier (default: 1.5)\n"
        ":param e: tolerance (default: 1E-8)\n"
        ":param k: kernel function (default: exp(-t^2/4))\n"
        ":return: M, S\n");
}
#endif

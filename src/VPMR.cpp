/*******************************************************************************
 * Copyright (C) 2024-2026 Theodore Chang
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

#include <algorithm>
#include <gmpxx.h>
#include <iomanip>
#include <memory>
#include <mpfr/exprtk_mpfr_adaptor.hpp>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include "Eigen/Core"
#include "Eigen/MatrixFunctions"
#ifdef HAVE_WINDOWS_H
#include <Windows.h>
#endif

using namespace mpfr;

struct Config {
    bool print_weight = false;
    bool print_singular_value = false;
    bool omit_trivial_terms = true;
    bool verbose = true;
    int max_terms = 20;
    int max_exponent = 4;
    int precision_bits = 512;
    int quadrature_order = 500;
    double precision_multiplier = 1.05;
    mpreal tolerance{"1E-8", precision_bits};
    std::string kernel = "exp(-t*t/4)";

    void print() const {
        if(!verbose) return;

        std::cout << std::scientific << std::setprecision(4);

        std::cout << "Using the following parameters:\n";
        std::cout << "       terms = " << max_terms << ".\n";
        std::cout << "    exponent = " << max_exponent << ".\n";
        std::cout << "   precision = " << precision_bits << ".\n";
        std::cout << " quad. order = " << quadrature_order << ".\n";
        std::cout << "  multiplier = " << precision_multiplier << ".\n";
        std::cout << "   tolerance = " << (2 * tolerance).toDouble() << ".\n";
        std::cout << "      kernel = " << kernel << ".\n\n";
    }
};

Config config;

mpreal MP_PI = const_pi();
mpreal MP_PI_HALF{MP_PI / 2};

class Quadrature {
    class LegendrePolynomial {
        class Evaluation {
            const size_t degree;

            mpreal _x, _v, _d;
            // Reused work buffers avoid allocating several MPFR temporaries
            // for every term of the Legendre recurrence.
            mpreal _p0, _p1, _p2, _tmp;

        public:
            explicit Evaluation(const mpreal& X, const size_t D)
                : degree(D)
                , _x(X)
                , _v(1, config.precision_bits)
                , _d(0, config.precision_bits)
                , _p0(1, config.precision_bits)
                , _p1(X)
                , _p2(0, config.precision_bits)
                , _tmp(0, config.precision_bits) { this->evaluate(X); }

            void evaluate(const mpreal& x) {
                this->_x = x;

                // P_0 = 1, P_1 = x, (i)P_i = (2i - 1)xP_{i-1} - (i - 1)P_{i-2}
                this->_p0 = 1;
                this->_p1 = x;

                for(size_t i = 2; i <= degree; ++i) {
                    this->_p2 = this->_p1;
                    this->_p2 *= x;
                    this->_p2 *= 2 * i - 1;

                    this->_tmp = this->_p0;
                    this->_tmp *= i - 1;

                    this->_p2 -= this->_tmp;
                    this->_p2 /= i;

                    this->_p0 = this->_p1;
                    this->_p1 = this->_p2;
                }

                this->_v = this->_p1;

                this->_d = this->_v;
                this->_d *= x;
                this->_d -= this->_p0;
                this->_d *= degree;
                this->_tmp = x;
                this->_tmp *= x;
                this->_tmp -= 1;
                this->_d /= this->_tmp;
            }

            [[nodiscard]] const mpreal& v() const { return this->_v; }
            [[nodiscard]] const mpreal& d() const { return this->_d; }
            [[nodiscard]] const mpreal& x() const { return this->_x; }
        };

        static tbb::affinity_partitioner ap;

        const mpreal ONE = mpreal(1, config.precision_bits);
        const mpreal TWO = mpreal(2, config.precision_bits);

    public:
        const size_t degree;

    private:
        std::unique_ptr<mpreal[]> _r, _w;

    public:
        explicit LegendrePolynomial(const size_t D)
            : degree(D > 2 ? D : 2)
            , _r(std::make_unique<mpreal[]>(degree))
            , _w(std::make_unique<mpreal[]>(degree)) {
            tbb::parallel_for(static_cast<size_t>(0), degree / 2 + 1, [&](const size_t i) {
            mpreal dr{1, config.precision_bits};

            Evaluation eval(cos(MP_PI * mpreal(4 * i + 3, config.precision_bits) / mpreal(4 * degree + 2, config.precision_bits)), degree);
            do {
                dr = eval.v() / eval.d();
                eval.evaluate(eval.x() - dr);
            }
            while(abs(dr) > std::numeric_limits<mpreal>::epsilon());

            this->_r[i] = eval.x();
            this->_w[i] = TWO / ((ONE - eval.x() * eval.x()) * eval.d() * eval.d()); }, ap);

            tbb::parallel_for(degree / 2, degree, [&](const size_t i) {
            this->_r[i] = -this->_r[degree - i - 1];
            this->_w[i] = this->_w[degree - i - 1]; }, ap);

            // Map the roots from [-1, 1] to the quadrature interval [0, pi]
            // once instead of repeating the affine transform for every moment.
            for(size_t i = 0; i < degree; ++i) this->_r[i] = MP_PI_HALF * this->_r[i] + MP_PI_HALF;
        }

        [[nodiscard]] const mpreal& node(const int i) const { return this->_r[i]; }
        [[nodiscard]] const mpreal& weight(const int i) const { return this->_w[i]; }
    };

    LegendrePolynomial poly;

public:
    explicit Quadrature(const size_t D)
        : poly(D) {}

    [[nodiscard]] size_t size() const { return poly.degree; }
    [[nodiscard]] const mpreal& node(const int i) const { return poly.node(i); }

    template<typename Function> mpreal integrate(Function&& f) const {
        return tbb::parallel_deterministic_reduce(
            tbb::blocked_range(0, static_cast<int>(poly.degree)), mpreal(0, config.precision_bits),
            [&](const tbb::blocked_range<int>& r, mpreal running_total) {
                for(auto i = r.begin(); i < r.end(); ++i) running_total += poly.weight(i) * f(i, poly.node(i));
                return running_total;
            },
            std::plus<>());
    }
};

tbb::affinity_partitioner Quadrature::LegendrePolynomial::ap{};

template<typename R, typename... A> class Cache {
public:
    Cache() = default;

    explicit Cache(std::function<R(A...)> f)
        : f_(std::move(f)) {}

    R operator()(A... a) {
        std::tuple<A...> key(a...);
        if(auto search = map_.find(key); search != map_.end()) return search->second;
        auto result = f_(a...);
        map_[key] = result;
        return result;
    }

private:
    std::function<R(A...)> f_;
    std::map<std::tuple<A...>, R> map_;
};

using big_int = mpz_class;

big_int comb_impl(unsigned, unsigned);

// cached function for combinatorial number calculation
Cache<big_int, unsigned, unsigned> comb(comb_impl);

// recursive implementation of combinatorial number calculation
big_int comb_impl(unsigned n, unsigned k) {
    if(k > n) return {0};
    if(k == 0 || k == n) return {1};
    if(k == 1 || k == n - 1) return {n};
    if(2 * k > n) k = n - k;

    return comb(n - 1, k - 1) + comb(n - 1, k);
}

// Estimate the number of bits needed to represent the largest binomial
// coefficient exactly. Use GMP's bit-length query instead of repeatedly
// dividing a big integer by two.
int estimate_precision_bits(const int max_terms, const double multiplier) {
    const auto comb_max = comb(4 * max_terms, 2 * max_terms);
    const auto comb_bits = static_cast<int>(mpz_sizeinbase(comb_max.get_mpz_t(), 2));
    return std::max(64, static_cast<int>(std::max(1.5, multiplier) * static_cast<double>(comb_bits + 4 * max_terms)));
}

// kernel function
class Expression {
    mpreal time;

    // The expression owns references to the symbol table. Rebuilding both on
    // every compile keeps repeated calls with different kernels independent.
    std::unique_ptr<exprtk::symbol_table<mpreal>> symbol_table;
    std::unique_ptr<exprtk::expression<mpreal>> expression;
    std::unique_ptr<exprtk::parser<mpreal>> parser;

    std::mutex eval_mutex;

public:
    bool compile() {
        symbol_table = std::make_unique<exprtk::symbol_table<mpreal>>();
        expression = std::make_unique<exprtk::expression<mpreal>>();
        parser = std::make_unique<exprtk::parser<mpreal>>();

        symbol_table->add_variable("t", time);
        symbol_table->add_constants();
        expression->register_symbol_table(*symbol_table);
        return parser->compile(config.kernel, *expression);
    }

    mpreal value(const mpreal& t) {
        std::scoped_lock lock(eval_mutex);
        time = t;
        return expression->value();
    }
};

Expression kernel{};

using mat = Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic>;
using vec = Eigen::Matrix<mpreal, Eigen::Dynamic, 1>;
using cx_vec = Eigen::Matrix<std::complex<mpreal>, Eigen::Dynamic, 1>;

// step 2
mat lyap(const vec& A, mat&& C) {
    tbb::parallel_for(0, static_cast<int>(A.size()), [&](const int i) { for(auto j = 0; j < A.size(); ++j) C(i, j) /= (A(i) + A(j)); });

    return std::move(C);
}

long pos(const vec& A) {
    const auto target_tolerance = .1 * config.tolerance;
    auto sum = mpreal(0, config.precision_bits);
    long P = 0;
    for(auto I = A.size() - 1; I >= 0; --I) {
        sum += A(I);
        if(sum > target_tolerance) {
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
    if(config.verbose) std::cout << "[2/6] Solving Lyapunov equation...\n";
    const mat S = lyap(A, lyap_rhs(B)).llt().matrixL();
    const mat L = lyap(A, lyap_rhs(C)).llt().matrixL();

    // step 4
    if(config.verbose) std::cout << "[3/6] Solving SVD...\n";
    const auto SVD = Eigen::BDCSVD<mat, Eigen::ComputeFullU>(S.transpose() * W.asDiagonal() * L);
    const auto& U = SVD.matrixU();
    const auto& SIG = SVD.singularValues();
    if(config.print_singular_value) {
        std::cout << "SIGMA = \n";
        for(auto I = 0; I < SIG.size(); ++I) std::cout << SIG(I).toString() << '\n';
    }

    // step 7
    const auto P = pos(SIG);
    if(config.verbose) std::cout << "[4/6] Transforming (P=" << P << ")...\n";
    if(P == SIG.size()) std::cerr << "WARNING: No singular value is smaller than the given tolerance.\n";

    // step 5
    auto SS = SIG;
    for(auto I = 0; I < SS.size(); ++I)
        if(SS(I) > std::numeric_limits<mpreal>::epsilon())
            SS(I) = pow(SS(I), -.5);
        else {
            std::cerr << "WARNING: Need to increase digits.\n";
            SS(I) = pow(std::numeric_limits<mpreal>::epsilon(), -.5);
        }

    // step 5
    const mat T = S * U * SS.asDiagonal();
    const auto LUT = T.fullPivLu();

    // step 6
    if(config.verbose) std::cout << "[5/6] Solving eigen decomposition...\n";
    const auto EIGEN = Eigen::EigenSolver<mat>(LUT.solve(A.asDiagonal() * T).block(0, 0, P, P));
    // step 8
    const auto X = EIGEN.eigenvectors();

    // step 9
    const cx_vec BB = LUT.solve(B).head(P);
    const cx_vec CC = ((C.cwiseProduct(W).transpose() * T).head(P) * X).transpose();

    return std::make_tuple(X.fullPivLu().solve(BB).cwiseProduct(CC), EIGEN.eigenvalues());
}

mpreal weight(const std::vector<mpreal>& k_hat, const int j) {
    auto sgn_bit = [](const int l) { return l % 2 == 0 ? 1 : -1; };

    const auto N = config.max_terms;
    const auto D = config.precision_bits;

    // Convert exact binomial coefficients directly from GMP integers. Going
    // through decimal strings loses the benefit of both libraries.
    auto binom = [D](const unsigned n, const unsigned k) {
        const auto value = comb(n, k);
        return mpreal(value.get_mpz_t(), D);
    };

    if(0 == j) {
        mpreal result = k_hat[0] / 2;
        for(auto l = 1; l <= N; ++l) result += sgn_bit(l) * k_hat[l];
        for(auto l = 1; l < N; ++l) result += sgn_bit(N + l) * mpreal(N - l, D) / mpreal(N, D) * k_hat[N + l];
        return result;
    }

    mpreal result{0, D};

    if(j > N) {
        for(auto l = j - N; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal((N - l) * (N + l), D) / mpreal(N * (N + l + j), D) *
                binom(N + l + j, N + l - j) * k_hat[N + l];
    }
    else {
        for(auto l = j; l <= N; ++l)
            result += sgn_bit(l - j) * mpreal(l, D) / mpreal(l + j, D) * binom(l + j, l - j) * k_hat[l];
        for(auto l = 1; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal((N - l) * (N + l), D) / mpreal(N * (N + l + j), D) *
                binom(N + l + j, N + l - j) * k_hat[N + l];
    }

    // unlike the original paper, we do not multiply 4^j here
    // to avoid multiplication of large complex numbers and large powers of 4
    return result;
}

std::vector<long> sort_index(const cx_vec& v) {
    std::vector<long> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // Keep the comparison in MPFR. Converting to double here discards the
    // working precision and can turn large coefficients into infinities.
    std::vector<mpreal> magnitude(v.size());
    for(auto I = 0; I < v.size(); ++I) magnitude[I] = norm(v(I));

    std::stable_sort(idx.begin(), idx.end(), [&v, &magnitude](const long i, const long j) {
        if(magnitude[i] > magnitude[j]) return true;
        if(magnitude[i] < magnitude[j]) return false;
        return v(i).imag() > v(j).imag();
    });

    return idx;
}

std::tuple<cx_vec, cx_vec> vpmr() {
    const auto quad = Quadrature(config.quadrature_order);

    // The kernel only depends on the quadrature node, so sample it once
    // instead of maintaining a global cache across calls and moments.
    const mpreal ONE = mpreal(1, config.precision_bits);
    std::vector<mpreal> sample(quad.size());
    for(size_t I = 0; I < sample.size(); ++I) sample[I] = kernel.value(-config.max_exponent * log((ONE + cos(quad.node(I))) >> 1));

    // The moments only depend on the kernel and quadrature rule. Computing
    // them once avoids repeating the same O(q) quadrature for every weight.
    const auto moment_count = 2 * config.max_terms;
    std::vector<mpreal> k_hat(moment_count);
    for(auto I = 0; I < moment_count; ++I) {
        k_hat[I] = quad.integrate([&](const int i, const mpreal& r) { return sample[i] * cos(mpreal(I, config.precision_bits) * r); });
        if(config.verbose) std::cout << "\r[1/6] Computing weights... [" << I + 1 << '/' << moment_count << ']' << std::flush;
    }
    if(config.verbose) std::cout << std::showpos << std::setprecision(16) << '\n';

    vec W = vec::Zero(moment_count);
    for(auto I = 0; I < W.size(); ++I) W(I) = weight(k_hat, I);

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
    if(config.omit_trivial_terms) {
        // The coefficients are sorted by descending magnitude. The old code
        // compared the squared magnitude with a linear tolerance, which removed
        // terms that were actually larger than the tolerance.
        const auto first = std::find_if(ID.begin(), ID.end(), [&M](const long I) { return abs(M(I)) <= config.tolerance; });
        ID.erase(first, ID.end());
    }

    if(config.verbose) std::cout << "[6/6] Done with " << ID.size() << " final terms.\n\n";

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
    std::cout << "    -n, --max-terms             <int>     number of terms (default: 20)\n";
    std::cout << "    -c, --max-exponent          <int>     maximum exponent (default: 4)\n";
    std::cout << "    -d, --precision-bits        <int>     number of precision bits (default: automatic)\n";
    std::cout << "    -q, --quadrature-order      <int>     quadrature order (default: 500)\n";
    std::cout << "    -m, --precision-multiplier  <float>   precision multiplier (default: 1.05)\n";
    std::cout << "    -e, --tolerance             <float>   tolerance (default: 1E-8)\n";
    std::cout << "    -k, --kernel                <string>  file name of kernel function (default uses: exp(-t*t/4))\n";
    std::cout << "        --no-omit                         do not omit trivial terms\n";
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
    std::string tolerance_string{"1E-8"};
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
            tolerance_string = argv[++I];
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
        else if(token == "--no-omit")
            config.omit_trivial_terms = false;
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
    const auto comb_digit = estimate_precision_bits(config.max_terms, config.precision_multiplier);

    if(!has_digit)
        config.precision_bits = comb_digit;
    else if(comb_digit >= config.precision_bits) {
        std::cout << "WARNING: Too few digits to hold combinatorial number, resetting digits to " << comb_digit << ".\n";
        config.precision_bits = comb_digit;
    }

    // set precision
    mpreal::set_default_prec(config.precision_bits);

    // Parse the tolerance only once the working precision is known.
    config.tolerance = mpreal(tolerance_string, config.precision_bits);

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
    const int n, const int c, const int d, const int q, const double m, const double e, const std::string& k, const bool omit, const bool verbose) {
    // The algorithm keeps its configuration and caches in globals. Serialize
    // calls and reset the configuration so successive Python calls cannot
    // influence one another.
    static std::mutex run_mutex;
    const std::scoped_lock lock(run_mutex);

    config = Config{};
    config.max_terms = std::max(1, n);
    config.max_exponent = std::max(1, c);
    config.precision_bits = std::max(1, d);
    config.quadrature_order = std::max(1, q);
    config.precision_multiplier = std::max(1.05, m);
    config.omit_trivial_terms = omit;
    config.verbose = verbose;
    if(!k.empty()) config.kernel = k;

    // check size
    const auto comb_digit = estimate_precision_bits(config.max_terms, config.precision_multiplier);

    if(0 == d)
        config.precision_bits = comb_digit;
    else if(comb_digit >= config.precision_bits) {
        if(config.verbose)
            std::cerr << "WARNING: Too few digits to hold combinatorial number, resetting digits to " << comb_digit << ".\n";
        config.precision_bits = comb_digit;
    }

    // set precision
    mpreal::set_default_prec(config.precision_bits);

    // Parse the tolerance only once the working precision is known.
    config.tolerance = mpreal(e, config.precision_bits);

    MP_PI = const_pi(config.precision_bits);
    MP_PI_HALF = MP_PI / 2;

    config.tolerance /= 2;

    if(!kernel.compile()) throw std::runtime_error("Cannot compile kernel function: " + config.kernel);

    config.print();

    const auto [M, S] = vpmr();

    std::vector<std::complex<double>> mm, ss;
    for(const auto& I : M) mm.emplace_back(I.real().toDouble(), I.imag().toDouble());
    for(const auto& I : S) ss.emplace_back(I.real().toDouble(), I.imag().toDouble());

    return {mm, ss};
}

PYBIND11_MODULE(_pyvpmr, m) {
    m.doc() = "The VPMR Algorithm";

    m.def(
        "vpmr", &vpmr_wrapper,
        "The VPMR Algorithm.\n\n"
        ":param n: number of terms (default: 20)\n"
        ":param c: maximum exponent (default: 4)\n"
        ":param d: number of precision bits (default: 512)\n"
        ":param q: quadrature order (default: 500)\n"
        ":param m: precision multiplier (default: 1.05)\n"
        ":param e: tolerance (default: 1E-8)\n"
        ":param k: kernel expression (default: exp(-t*t/4))\n"
        ":param omit: omit trivial terms (default: True)\n"
        ":param verbose: print progress information (default: True)\n"
        ":return: (M, S)\n",
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        pybind11::kw_only(), pybind11::arg("n") = 20, pybind11::arg("c") = 4, pybind11::arg("d") = 0, pybind11::arg("q") = 500, pybind11::arg("m") = 1.05, pybind11::arg("e") = 1E-8, pybind11::arg("k") = "", pybind11::arg("omit") = true, pybind11::arg("verbose") = true);
}
#endif

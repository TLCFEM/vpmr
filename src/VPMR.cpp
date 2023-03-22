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

#include <iostream>
#include "Cache.hpp"
#include "BigInt.hpp"
#include "Gauss.hpp"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "unsupported/Eigen/MPRealSupport"
#include "unsupported/Eigen/MatrixFunctions"
#include <chrono>

using namespace mpfr;
using namespace Eigen;

int NC = 4; // maximum exponent
int N = 10; // number of terms
int DIGIT = 512; // number of digits
int QUAD_ORDER = 500; // number of quadrature points
mpreal TOL = mpreal("1E-8"); // tolerance

BigInt comb_impl(unsigned, unsigned);

// cached function for combinatorial number calculation
Cache<BigInt, unsigned, unsigned> comb(comb_impl);

// recursive implementation of combinatorial number calculation
BigInt comb_impl(unsigned n, unsigned k) {
    if (k == 0 || k == n) return {1};
    if (k == 1 || k == n - 1) return {n};
    if (k > n) return {0};
    if (2 * k > n) k = n - k;

    return comb(n - 1, k - 1) + comb(n - 1, k);
}

// integrand in eq. 2.2 --- K(t)cos(jt)
class Kernel {
    const int j;
public:
    explicit Kernel(const int J) : j(J) {}

    mpreal operator()(const mpreal &r) const {
        return smooth_function(-NC * log((mpreal(1) + cos(r)) >> 1)) * cos(j * r);
    }

    // change this function to the target kernel function in convolution
    static mpreal smooth_function(const mpreal &x) {
        return exp(-x * x / 4);
    }
};

using mat = Eigen::Matrix<mpreal, Dynamic, Dynamic>;
using vec = Eigen::Matrix<mpreal, Dynamic, 1>;
using cx_vec = Eigen::Matrix<std::complex<mpreal>, Dynamic, 1>;

// step 2
mat lyap(const mat &A, const mat &C) {
    return Eigen::internal::matrix_function_solve_triangular_sylvester(A, A, C);
}

Index pos(const vec &A) {
    const auto target = TOL * .5;
    auto sum = mpreal(0);
    Index P = 0;
    for (auto I = A.size() - 1; I >= 0; --I) {
        sum += A(I);
        if (sum > target) {
            P = I + 1;
            break;
        }
    }
    return P;
}

std::tuple<cx_vec, cx_vec> model_reduction(const vec &A, const vec &B, const vec &C) {
    // step 3
    const mat S = lyap(A.asDiagonal(), -B * B.transpose()).llt().matrixL();
    const mat L = lyap(A.asDiagonal(), -C * C.transpose()).llt().matrixL();

    // step 4
    const auto SVD = BDCSVD<mat, ComputeFullV | ComputeFullU>(S.transpose() * L);
    const auto &U = SVD.matrixU();
    const auto &SIG = SVD.singularValues();

    // step 7
    const auto P = pos(SIG);
    if (P == SIG.size())
        throw std::invalid_argument("No solution can be found, try to increase tolerance (e) or number of terms (n).");

    auto SS = SIG;
    for (auto I = 0; I < SS.size(); ++I) SS(I) = pow(SS(I), -.5);

    // step 5
    const mat T = S * U * SS.asDiagonal();
    const auto LUT = T.lu();
    // step 6
    const auto EIGEN = EigenSolver<mat>(LUT.solve(A.asDiagonal() * T).block(0, 0, P, P));
    // step 8
    const auto X = EIGEN.eigenvectors();

    // step 9
    const cx_vec BB = LUT.solve(B).head(P);
    const cx_vec CC = ((C.transpose() * T).head(P) * X).transpose();

    return std::make_tuple(X.fullPivLu().solve(BB).cwiseProduct(CC), EIGEN.eigenvalues());
}

mpreal weight(const Quadrature &quad, const int j) {
    auto k_hat = [&quad](const int l) { return 2 / MP_PI * quad.integrate(Kernel(l)); };
    auto sgn_bit = [](const int l) { return l % 2 == 0 ? 1 : -1; };

    if (0 == j) {
        mpreal result = k_hat(0) / 2;
        for (auto l = 1; l <= N; ++l) result += sgn_bit(l) * k_hat(l);
        for (auto l = 1; l < N; ++l) result += sgn_bit(N + l) * mpreal(N - l) / mpreal(N) * k_hat(N + l);
        return result;
    }

    mpreal result{0};

    if (j > N) {
        for (auto l = j - N; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal(N - l) / mpreal(N) * mpreal(N + l) / mpreal(N + l + j) *
                      mpreal(comb(N + l + j, N + l - j).to_string()) * k_hat(N + l);
    } else {
        for (auto l = j; l <= N; ++l)
            result += sgn_bit(l - j) * mpreal(l) / mpreal(l + j) *
                      mpreal(comb(l + j, l - j).to_string()) * k_hat(l);
        for (auto l = 1; l < N; ++l)
            result += sgn_bit(N + l - j) * mpreal(N - l) / mpreal(N) * mpreal(N + l) / mpreal(N + l + j) *
                      mpreal(comb(N + l + j, N + l - j).to_string()) * k_hat(N + l);
    }

    return result * pow(mpreal(4), j);
}

std::tuple<cx_vec, cx_vec> vpmr() {
    const auto quad = Quadrature(QUAD_ORDER);

    vec W = vec::Zero(2 * N);
    for (auto I = 0; I < W.size(); ++I) W(I) = weight(quad, I);

    // step 1
    vec A = vec::Zero(W.size() - 1), B = vec::Zero(A.size()), C = vec::Zero(A.size());
    for (decltype(W.size()) I = 0, J = 1; I < A.size(); ++I, ++J) {
        A(I) = -mpreal(J) / NC;
        B(I) = sqrt(abs(W(J)));
        C(I) = sgn(W(J)) * B(I);
    }

    auto [M, S] = model_reduction(A, B, C);

    cx_vec MM = cx_vec::Zero(M.size() + 1);
    MM(0) = W(0);
    MM.tail(M.size()) = M;

    cx_vec SS = cx_vec::Zero(S.size() + 1);
    SS(0) = mpreal(0);
    SS.tail(S.size()) = -S;

    return std::make_tuple(MM, SS);
}

int print_helper() {
    std::cout << "Usage: vpmr [options]\n" << std::endl;
    std::cout << "Options:\n" << std::endl;
    std::cout << "   -nc <int>    number of exponent (default: 4)\n";
    std::cout << "   -n <int>     number of terms (default: 10)\n";
    std::cout << "   -d <int>     number of digits (default: 512)\n";
    std::cout << "   -q <int>     quadrature order (default: 500)\n";
    std::cout << "   -e <float>   tolerance (default: 1E-8)\n";
    std::cout << "   -h           print this help message\n";
    return 0;
}

int main(const int argc, const char **argv) {
    const auto start = std::chrono::high_resolution_clock::now();

    // set precision
    mpreal::set_default_prec(DIGIT);

    for (auto I = 1; I < argc; ++I) {
        const auto token = std::string(argv[I]);

        if (token == "-nc") NC = std::stoi(argv[++I]);
        else if (token == "-n") N = std::stoi(argv[++I]);
        else if (token == "-d") DIGIT = std::stoi(argv[++I]);
        else if (token == "-q") QUAD_ORDER = std::stoi(argv[++I]);
        else if (token == "-e") TOL = mpreal(argv[++I]);
        else if (token == "-h") return print_helper();
        else {
            std::cerr << "Unknown option: " << token << std::endl;
            return 1;
        }
    }

    std::cout << std::scientific << std::setprecision(4);

    std::cout << "Using the following parameters:\n";
    std::cout << "        nc = " << NC << ".\n";
    std::cout << "         n = " << N << ".\n";
    std::cout << "     order = " << QUAD_ORDER << ".\n";
    std::cout << " precision = " << DIGIT << ".\n";
    std::cout << " tolerance = " << TOL.toDouble() << ".\n\n";

    std::cout << std::scientific << std::showpos << std::setprecision(16);

    try {
        // run VPMR algorithm
        const auto [M, S] = vpmr();

        std::cout << "M = \n";
        for (auto I = 0; I < M.size(); ++I) std::cout << M(I).real().toDouble() << M(I).imag().toDouble() << "j\n";

        std::cout << "S = \n";
        for (auto I = 0; I < S.size(); ++I) std::cout << S(I).real().toDouble() << S(I).imag().toDouble() << "j\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << std::noshowpos << "\nRunning time: " << duration << " s.\n";

    return 0;
}
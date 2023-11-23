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

#include "VPMR.h"
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using cx_vec = Eigen::Matrix<std::complex<mpreal>, Eigen::Dynamic, 1>;

std::tuple<cx_vec, cx_vec> vpmr();

std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>> vpmr_wrapper(
    const int n, const int d, const int q, const int m, const int nc, const double e, const std::string& k) {
    N = n;
    DIGIT = d;
    QUAD_ORDER = q;
    SCALE = m;
    NC = nc;
    TOL = mpreal(e);
    if(!k.empty()) KERNEL = k;

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

PYBIND11_MODULE(pyvpmr, m) {
    m.def("vpmr", &vpmr_wrapper, pybind11::call_guard<pybind11::gil_scoped_release>(), pybind11::kw_only(), pybind11::arg("n") = 10, pybind11::arg("d") = 512, pybind11::arg("q") = 500, pybind11::arg("m") = 6, pybind11::arg("nc") = 4, pybind11::arg("e") = 1E-8, pybind11::arg("k") = "");
}
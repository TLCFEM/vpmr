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

std::tuple<std::vector<std::complex<double>>, std::vector<std::complex<double>>> vpmr_wrapper() {
    Expression kernel;
    if(!kernel.compile()) {
        std::cerr << "Cannot compile kernel function: " << KERNEL << ".\n";
        throw std::runtime_error("Cannot compile kernel function.");
    }
    Integrand::set_expression(&kernel);

    const auto [M, S] = vpmr();

    std::vector<std::complex<double>> m(M.size()), s(S.size());
    for(const auto& I : M) m.emplace_back(I.real().toDouble(), I.imag().toDouble());
    for(const auto& I : S) s.emplace_back(I.real().toDouble(), I.imag().toDouble());

    return std::make_tuple(m, s);
}

PYBIND11_MODULE(pyvpmr, m) {
    m.def("vpmr", &vpmr_wrapper, pybind11::call_guard<pybind11::gil_scoped_release>());
}
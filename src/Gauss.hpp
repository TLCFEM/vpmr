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

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include "mpreal.h"
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

using namespace mpfr;

mpreal MP_PI = const_pi(1000);
mpreal MP_PI_HALF{MP_PI / 2};

class Evaluation {
    const size_t degree;

    mpreal _x, _v, _d;
public:
    explicit Evaluation(const mpreal &X, const size_t D) : degree(D), _x(X), _v(1), _d(0) { this->evaluate(X); }

    void evaluate(const mpreal &x) {
        this->_x = x;

        mpreal left{x}, right{1};

        for (int i = 2; i <= degree; ++i) {
            this->_v = ((mpreal(2) * i - mpreal(1)) * x * left - (i - mpreal(1)) * right) / i;

            right = left;
            left = this->_v;
        }

        this->_d = degree / (x * x - mpreal(1)) * (x * this->_v - right);
    }

    [[nodiscard]] mpreal v() const { return this->_v; }

    [[nodiscard]] mpreal d() const { return this->_d; }

    [[nodiscard]] mpreal x() const { return this->_x; }
};

class LegendrePolynomial {
public:
    const size_t degree;
private:
    std::unique_ptr<mpreal[]> _r, _w;
public:
    explicit LegendrePolynomial(const size_t D)
            : degree(D > 2 ? D : 2),
              _r(std::make_unique<mpreal[]>(degree)),
              _w(std::make_unique<mpreal[]>(degree)) {
        for (size_t i = 0; i < degree / 2 + 1; ++i) {
            mpreal dr{1};

            Evaluation eval(cos(MP_PI * mpreal(4 * i + 3) / mpreal(4 * degree + 2)), degree);
            do {
                dr = eval.v() / eval.d();
                eval.evaluate(eval.x() - dr);
            } while (abs(dr) > std::numeric_limits<mpreal>::epsilon());

            this->_r[i] = eval.x();
            this->_w[i] = mpreal(2) / ((mpreal(1) - eval.x() * eval.x()) * eval.d() * eval.d());
        }

        for (size_t i = degree - 1; i >= degree / 2; --i) {
            this->_r[i] = -this->_r[degree - i - 1];
            this->_w[i] = this->_w[degree - i - 1];
        }
    }

    [[nodiscard]] mpreal root(int i) const { return this->_r[i]; }

    [[nodiscard]] mpreal weight(int i) const { return this->_w[i]; }
};

class Quadrature {
    LegendrePolynomial poly;
public:
    explicit Quadrature(const size_t D) : poly(D) {}

    template<typename Function>
    mpreal integrate(Function &&f) const {
//        mpreal total(0);
        auto total = tbb::parallel_reduce(
                tbb::blocked_range<int>(0, int(poly.degree)), mpreal(0),
                [&](const tbb::blocked_range<int> &r, mpreal running_total) {
                    for (auto i = r.begin(); i < r.end(); ++i)
                        running_total += poly.weight(i) * f(i, MP_PI_HALF * poly.root(i) + MP_PI_HALF);
                    return running_total;
                }, std::plus<>());
//        for (int i = 0; i < poly.degree; ++i) total += poly.weight(i) * f(i, MP_PI_HALF * poly.root(i) + MP_PI_HALF);
        return total;
    }
};

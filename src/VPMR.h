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

#include <mpfr/exprtk_mpfr_adaptor.hpp>
#include <mutex>
#include <string>
#include <tbb/concurrent_unordered_map.h>
#include "Eigen/Core"

using namespace mpfr;

extern bool OUTPUT_W;      // output W
extern bool OUTPUT_EIG;    // output eigenvalues
extern int NC;             // maximum exponent
extern int N;              // number of terms
extern int DIGIT;          // number of digits
extern int QUAD_ORDER;     // number of quadrature points
extern int SCALE;          // scaling
extern mpreal TOL;         // tolerance
extern std::string KERNEL; // default kernel

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
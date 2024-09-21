#pragma once

#include <iomanip>
#include <string>
#include "mpreal.h"

using namespace mpfr;

struct Config {
    bool print_weight = false;
    bool print_singular_value = false;
    int max_terms = 10;
    int max_exponent = 4;
    int precision_bits = 512;
    int quadrature_order = 500;
    double precision_multiplier = 1.5;
    mpreal tolerance = mpreal("1E-8", precision_bits);
    std::string kernel = "exp(-t*t/4)";

    void print() const {
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
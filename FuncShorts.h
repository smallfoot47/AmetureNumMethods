#pragma once
#include <iostream>
#include <vector>

#include "NumMethods.h"

#define MAKE_FUNC(def) [&](double x) { def } //eg:- MAKE_FUNCTION( return pow(x, 2) - 1; )
#define PRINT_AT_VALUE(function, value) std::cout << #function << "(" << #value << ") = " << function(value) << std::endl //eg:- PRINT_AT_VALUE(abs, -2);
#define PRINT_FUNC_COORD(answer) std::cout << #answer << "(" << answer.first << ") = " << answer.second << std::endl //eg:- PRINT_FUNC_COORD({ 2.5, 1 });

typedef std::function<double(double)> func_t;//double->double function storage object

func_t MAKE_POLY_FUNC(std::vector<double> coeffs); //eg:- MAKE_POLY_FUNC({2, 3, 1}); => 2x^2 + 3x + 1
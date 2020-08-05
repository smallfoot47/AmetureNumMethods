#pragma once
#include <functional>
#include <limits>

std::function<double(double)> derive(std::function<double(double)>& f); //Numerical approxiamation of derivative
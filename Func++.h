#pragma once
#include <functional>
#include <limits>

template <class T> struct Interval;

template <class N_t = double>
N_t gamma(N_t x) {//currently restricted to real numbers
	if (x > (N_t)1)
		return x * gamma(x - (N_t)1);
	return (N_t)1;
}

template <class N_t = double>
inline std::function<N_t(N_t)> derive(std::function<N_t(N_t)>& f, const double lifting_param = 10000.0) //Numerical approximation of derivative
{
	N_t dx = lifting_param * std::numeric_limits<N_t>::epsilon();
	N_t denom = (N_t)2 * dx;

	return [=](N_t x)
	{
		return (f(x + dx) - f(x - dx)) / denom;
	};
}

#pragma once
#include <functional>
#include <limits>

std::function<double(double)> derive(std::function<double(double)>& f); //Numerical approxiamation of derivative
/*
template <class N_t>
inline std::function<N_t(N_t)> derive(std::function<N_t(N_t)>& f)
{
	N_t dx = std::numeric_limits<N_t>::epsilon();
	N_t denom = 2 * dx;

	return [&](N_t x)
	{
		return (f(x + dx) - f(x - dx)) / denom;
	};
}
*/
template <class N_t>
inline N_t sigmoid(N_t x)//Mathematical sigmoid function
{
	return 1 / (exp(-x) + 1);
}

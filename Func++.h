#pragma once
#include <functional>
#include <limits>

template <class N_t = double>
inline std::function<N_t(N_t)> derive(std::function<N_t(N_t)>& f) //Numerical approximation of derivative
{
	N_t dx = std::numeric_limits<N_t>::epsilon();
	N_t denom = 2 * dx;

	return [&](N_t x)
	{
		return (f(x + dx) - f(x - dx)) / denom;
	};
}

template <class N_t = double>
inline N_t sigmoid(N_t x)//Mathematical sigmoid function
{
	return 1.0 / (exp(-x) + 1.0);
}

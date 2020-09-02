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

template <class N_t = double, class C_t = size_t>
N_t permutations(N_t space_size, C_t groups) {
	if (groups > (C_t)0)
		return space_size * permutations(space_size - (N_t)1, groups - (C_t)1);
	return (N_t)1;
}

template <class N_t = double, class C_t = size_t>
N_t combinations(N_t space_size, C_t groups)
{
	return permutations(space_size, groups) / gamma(groups);
}

template<class N_t = double>
N_t clamp(N_t x, Interval<N_t> bounds = { (N_t)0, (N_t)1 }) {
	if (x < bounds)
		return bounds.lower_bound;
	if (x > bounds)
		return bounds.upper_bound;
	return x;
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

template <class N_t = double>
inline N_t sigmoid(N_t x)//Mathematical sigmoid function
{
	return (N_t)1 / (exp<N_t>(-x) + (N_t)1);
}

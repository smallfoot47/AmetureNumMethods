#include "Func++.h"

std::function<double(double)> derive(std::function<double(double)>& f)
{
	double dx = std::numeric_limits<double>::epsilon();
	double denom = 2 * dx;
	return [&](double x)
	{
		return (f(x + dx) - f(x - dx)) / denom;
	};
}
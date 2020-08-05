#include "NumMethods.h"
#include "Func++.h"

Coord2D IntermediateValue(std::function<double(double)> f, std::pair<double, double> interval, uint8_t n)
{
	double x = (interval.first + interval.second) / 2;
	double y = f(x);

	if (n)
	{
		if (y * f(interval.first) < 0)
			return IntermediateValue(f, { interval.first, x }, n - 1);
		else if (y * f(interval.second) < 0)
			return IntermediateValue(f, { x, interval.second }, n - 1);
	}

	return { x, y };
}

inline Coord2D NewtonsMethod(std::function<double(double)> f, double x0, uint8_t n)
{
	return NewtonsMethod(f, derive(f), x0, n);
}

Coord2D NewtonsMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, uint8_t n)
{
	if (n && df(x0))
		return NewtonsMethod(f, df, x0 - f(x0) / df(x0), n - 1);
			
	return { x0, f(x0) };
}

inline Coord2DMap Function_to_Map(std::function<double(double)> f, std::pair<double, double> interval)
{
	return Function_to_Map(f, interval, std::numeric_limits<double>::epsilon());
}

Coord2DMap Function_to_Map(std::function<double(double)> f, std::pair<double, double> interval, uint32_t n)
{
	Coord2DMap F(n);

	double epsilon = (interval.second - interval.first) / n;

	for (uint32_t i = 0; i < n; ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

	return F;
}

Coord2DMap Function_to_Map(std::function<double(double)> f, std::pair<double, double> interval, double epsilon)
{
	Coord2DMap F((size_t)((interval.second - interval.first) / epsilon));

	for (uint32_t i = 0; i < F.size(); ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

	return F;
}

std::function<double(double)> Interpolate(Coord2DMap map, std::pair<double, double> interval)
{
	const double domain_width = (interval.second - interval.first);
	const double internodal_width = domain_width / map.size();
	const double bias = -interval.first;
	
	return [&](double x)
	{
		x = (x + bias) / domain_width;
		return map[floor(x) * domain_width + remainder(x, 1) * internodal_width];
	};
}

std::function<double(double)> Interpolate(Coord2DMap map, std::pair<double, double> interval, uint8_t poly_factor, uint32_t disloc_factor)
{
	double domain_width = (interval.second - interval.first);
	double internodal_width = domain_width / map.size();
	double bias = disloc_factor * internodal_width - interval.first;
	internodal_width *= poly_factor;

	return [&](double x)
	{
		x = (x + bias) / domain_width;
		return map[floor(x) * domain_width + remainder(x, 1) * internodal_width];
	};
}

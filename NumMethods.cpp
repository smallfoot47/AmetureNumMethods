#include "NumMethods.h"

Coord2D BisectionMethod(std::function<double(double)> f, Interval<double> interval, uint8_t n)
{
	double x = interval.mean();
	double y = f(x);

	if (n > 1u)
	{
		if (y * f(interval.lower) < 0.0)
			return BisectionMethod(f, { interval.lower, x }, n - 1u);
		else if (y * f(interval.upper) < 0.0)
			return BisectionMethod(f, { x, interval.upper }, n - 1u);
	}

	return { x, y };
}

inline Coord2D NewtonsMethod(std::function<double(double)> f, double x0, uint8_t n)
{
	return NewtonsMethod(f, derive(f), x0, n);
}

Coord2D NewtonsMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, uint8_t n)
{
	if (df(x0))
	{
		x0 -= f(x0) / df(x0);

		if (n > 1u)
			return NewtonsMethod(f, df, x0, n - 1u);
	}
	return { x0, f(x0) };
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval)
{
	return Function_to_Map(f, interval, 1E-5);
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, uint32_t n)
{
	Coord2DMap F(n--);

	double epsilon = interval / n;

	for (uint32_t i = 0u; i <= n; ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

#include "NumMethods.h"

Coord2D BisectionMethod(std::function<double(double)> f, Interval<double> interval, uint8_t n)
{
	double x = interval.mean();
	double y = f(x);

	if (n > 1u)
	{
		if (y * f(interval.lower) < 0.0)
			return BisectionMethod(f, { interval.lower, x }, n - 1u);
		else if (y * f(interval.upper) < 0.0)
			return BisectionMethod(f, { x, interval.upper }, n - 1u);
	}

	return { x, y };
}

inline Coord2D NewtonsMethod(std::function<double(double)> f, only<double> x0, uint8_t n)
{
	return NewtonsMethod(f, derive(f), x0, n);
}

Coord2D NewtonsMethod(std::function<double(double)> f, std::function<double(double)> df, only<double> x0, uint8_t n)
{
	if (df(x0))
	{
		x0 -= f(x0) / df(x0);

		if (n > 1u)
			return NewtonsMethod(f, df, x0, n - 1u);
	}
	return { x0, f(x0) };
}

Coord2D SecantMethod(std::function<double(double)> f, only<double> x0, only<double> x1, uint8_t n)
{
	double F0 = f(x0), F1 = f(x1);

	if(F1 != F0)
	{
		double x2 = x1 - F1 * (x1 - x0) / (F1 - F0);
		
		if (n > 2)
			return SecantMethod(f, x1, x2, n - 1);
		
		return { x2, f(x2) };
	}

	return { x1, F1 };
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval)
{
	return Function_to_Map(f, interval, 1E-5);
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, uint32_t n)
{
	Coord2DMap F(n--);

	double epsilon = interval / n;

	for (uint32_t i = 0u; i <= n; ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

	return F;
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, double epsilon)
{
	uint32_t size = (uint32_t)(interval / epsilon) + 1u;
	Coord2DMap F; F.reserve(size);
	
	for (uint32_t i = 0u; i < size; ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

	return F;
}

std::function<double(double)> Interpolate(Coord2DMap& map, const Interval<double>& interval)
{
	return Interpolate(map, interval, map.size());
}

std::function<double(double)> Interpolate(Coord2DMap& map, const Interval<double>& interval, uint32_t n)
{
	const double internodal_width = interval / (n - 1u);
	const double bias = -interval.first;

	return [=](double x) mutable
	{
		if (interval[x])
		{
			x = (x + bias) / internodal_width;
			double foot = floor(x) * internodal_width - bias;
			return LinearInterpolate(map[foot], map[foot + internodal_width], abs(remainder(x, 1.0)));
		}

		return 0.0;
	};
}

std::function<double(double)> Interpolate(Coord2DMap& map, const Interval<double>& interval, uint16_t poly_factor, uint32_t disloc_factor)
{
	double internodal_width = interval / (map.size() - 1u);
	double bias = disloc_factor * internodal_width - interval.first;
	internodal_width *= poly_factor;

	return [=](double x) mutable
	{
		if (interval[x])
		{
			x = (x + bias) / internodal_width;
			double foot = floor(x) * internodal_width + interval.first;
			return LinearInterpolate(map[foot], map[foot + internodal_width], abs(remainder(x, 1.0)));
		}

		return 0.0;
	};
}

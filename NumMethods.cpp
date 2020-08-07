#include "NumMethods.h"

Coord2D IntermediateValue(std::function<double(double)> f, Interval<double> interval, uint8_t n)
{
	double x = interval / 2;
	double y = f(x);

	if (n > 1)
	{
		if (y * f(interval.lower) < 0)
			return IntermediateValue(f, { interval.lower, x }, n - 1);
		else if (y * f(interval.upper) < 0)
			return IntermediateValue(f, { x, interval.upper }, n - 1);
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

		if (n > 1)
			return NewtonsMethod(f, df, x0, n - 1);
	}
	return { x0, f(x0) };
}

inline Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval)
{
	return Function_to_Map(f, interval, std::numeric_limits<double>::epsilon());
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, uint32_t n)
{
	Coord2DMap F(n--);

	double epsilon = interval / n;

	for (uint32_t i = 0; i <= n; ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

	return F;
}

Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, double epsilon)
{
	Coord2DMap F((size_t)(interval / epsilon + 1));

	for (uint32_t i = 0; i < F.size(); ++i)
	{
		double x = interval.first + i * epsilon;
		F[x] = f(x);
	}

	return F;
}

std::function<double(double)> Interpolate(Coord2DMap map, Interval<double> interval)
{
	const double internodal_width = interval / map.size();
	const double bias = -interval.first;

	return [&](double x)->double
	{
		if (interval[x])
		{
			x = (x + bias) / internodal_width;
			double foot = floor(x) * internodal_width - bias;
			return LinearInterpolate(map[foot], map[foot + internodal_width], remainder(x, 1));
		}

		return 0;
	};
}

std::function<double(double)> Interpolate(Coord2DMap map, Interval<double> interval, uint8_t poly_factor, uint32_t disloc_factor)
{
	double internodal_width = interval / map.size();
	double bias = disloc_factor * internodal_width - interval.first;
	internodal_width *= poly_factor;

	return [&](double x)->double
	{
		if (interval[x])
		{
			x = (x + bias) / internodal_width;
			double foot = floor(x) * internodal_width + interval.first;
			return LinearInterpolate(map[foot], map[foot + internodal_width], remainder(x, 1));
		}

		return 0;
	};
}

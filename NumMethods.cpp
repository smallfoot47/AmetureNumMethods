#include "NumMethods.h"

namespace IterationBasedApproximation
{
	Coord2D BisectionMethod(func_t f, Interval<double> interval, uint8_t n)
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

	inline Coord2D NewtonsMethod(func_t f, only<double> x0, uint8_t n)
	{
		return NewtonsMethod(f, derive(f), x0, n);
	}

	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, uint8_t n)
	{
		/*
		  df is assumed to be the explicitly stated derivative of f w.r.t it's input variable
		*/
		if (df(x0))
		{
			x0 -= f(x0) / df(x0);

			if (n > 1u)
				return NewtonsMethod(f, df, x0, n - 1u);
		}
		return { x0, f(x0) };
	}

	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, uint8_t n)
	{
		double F0 = f(x0), F1 = f(x1);

		if (F1 != F0)
		{
			double x2 = x1 - F1 * (x1 - x0) / (F1 - F0);

			if (n > 2)
				return SecantMethod(f, x1, x2, n - 1u);

			return { x2, f(x2) };
		}

		return { x1, F1 };
	}

	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, uint8_t n)
	{
		/*
		  Accepting in interval as { ai, bi }
		*/
		auto SMApprox = SecantMethod(f, interval.left, interval.right, 2u);

		if (SMApprox.dependent * f(interval.left) < 0.0)       //f(x[i+1]) * f(ai) < 0: root probably lies in { ai, x[i+1] }
			interval.upper_bound = SMApprox.independent;     //interval => { ai, x[i+1] }
		else if (SMApprox.dependent * f(interval.right) < 0.0) //f(x[i+1]) * f(bi) < 0: root probably lies in { x[i+1], bi }
			interval.lower_bound = SMApprox.independent;     //interval => { x[i+1], bi }
		else
			return SMApprox;

		if (n > 1u)
			return RegulaFalsiMethod(f, interval, n - 1u);

		return SMApprox;
	}

	Coord2D FixedPointMethod(func_t f, only<double> x0, uint8_t n)
	{
		double F0 = f(x0);

		if (n > 1)
			return FixedPointMethod(f, F0, n - 1);

		return { F0, f(F0) };
	}
}

namespace ErrorBasedApproximation
{
	Coord2D BisectionMethod(func_t f, Interval<double> interval, only<double> accuracy)
	{

		Coord2D Approx = IterationBasedApproximation::BisectionMethod(f, interval, 1u);
		double last_x = Approx.independent;
		double magnitude = 0.0;

		do {
			if (Approx.dependent * f(interval.left) < 0.0)
				interval.upper_bound = Approx.independent;
			else
				if (Approx.dependent * f(interval.right) < 0.0)
					interval.lower_bound = Approx.independent;
				else
					return Approx;

			magnitude = -log10(abs(Approx.independent - last_x));

			last_x = Approx.independent;
			Approx = IterationBasedApproximation::BisectionMethod(f, interval, 1u);

		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D NewtonsMethod(func_t f, only<double> x0, only<double> accuracy)
	{

		Coord2D Approx = IterationBasedApproximation::NewtonsMethod(f, x0, 1u);
		double last_x = Approx.independent;
		double magnitude = 0.0;

		do {
			magnitude = -log10(abs(Approx.independent - last_x));

			last_x = Approx.independent;
			Approx = IterationBasedApproximation::NewtonsMethod(f, last_x, 1u);

		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, only<double> accuracy)
	{
		Coord2D Approx = IterationBasedApproximation::NewtonsMethod(f, df, x0, 1u);
		double last_x = Approx.independent;
		double magnitude = 0.0;

		do {
			magnitude = -log10(abs(Approx.independent - last_x));

			last_x = Approx.independent;
			Approx = IterationBasedApproximation::NewtonsMethod(f, df, last_x, 1u);

		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, only<double> accuracy)
	{
		Coord2D Approx = IterationBasedApproximation::SecantMethod(f, x0, x1, 2u);
		double last_x = Approx.independent;
		double magnitude = 0.0;

		do {
			magnitude = -log10(abs(Approx.independent - last_x));

			last_x = Approx.independent;
			Approx = IterationBasedApproximation::SecantMethod(f, last_x, Approx.x, 2u);

		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, only<double> accuracy)
	{
		Coord2D Approx = IterationBasedApproximation::SecantMethod(f, interval.left, interval.right, 2u);
		double magnitude = 0.0;
		double last_x = Approx.independent;

		do {
			if (Approx.dependent * f(interval.left) < 0.0)       //f(x[i+1]) * f(ai) < 0: root probably lies in { ai, x[i+1] }
				interval.upper_bound = Approx.independent;       //interval => { ai, x[i+1] }
			else
				if (Approx.dependent * f(interval.right) < 0.0)      //f(x[i+1]) * f(bi) < 0: root probably lies in { x[i+1], bi }
					interval.lower_bound = Approx.independent;       //interval => { x[i+1], bi }
				else
					return Approx;

			magnitude = -log10(abs(Approx.independent - last_x));

			last_x = Approx.independent;
			Approx = IterationBasedApproximation::SecantMethod(f, last_x, Approx.x, 2u);

		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D FixedPointMethod(func_t f, only<double> x0, only<double> accuracy)
	{
		Coord2D Approx; Approx.dependent = x0;
		double magnitude = 0.0;

		do {
			Approx.independent = Approx.dependent;
			Approx.dependent = f(Approx.independent);

			magnitude = -log10(abs(Approx.dependent - Approx.independent));
		
		} while (magnitude < accuracy);
	}
}

Coord2DMap Function_to_Map(func_t f, Interval<double> interval)
{
	return Function_to_Map(f, interval, 1E-5);
}

Coord2DMap Function_to_Map(func_t f, Interval<double> interval, uint32_t n)
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

Coord2DMap Function_to_Map(func_t f, Interval<double> interval, double epsilon)
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

func_t Interpolate(Coord2DMap& map, const Interval<double>& interval)
{
	/*
	  Linearly interpolates between each consecutive node in map
	*/
	return Interpolate(map, interval, map.size());
}

func_t Interpolate(Coord2DMap& map, const Interval<double>& interval, uint32_t n)
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

func_t Interpolate(Coord2DMap& map, const Interval<double>& interval, uint16_t poly_factor, uint32_t disloc_factor)
{
	/*
	  poly_factor determines the number of nodes to skip over to reach the next relavent node for interpolation
	  disloc_factor determines the number of nodes the initial node is displaced by(and this dislocation carries over to all subsequent nodes)
	*/
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

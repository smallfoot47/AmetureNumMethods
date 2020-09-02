#include "NumMethods.h"
#include "StepTracking.h"

namespace IterationBasedApproximation
{
	Coord2D BisectionMethod(func_t f, Interval<double> interval, uint16_t n)
	{
		double x = interval.mean();
		double y = f(x);

		{
			PUSH_STEP_LOG(
				std::cout << std::endl << "f({" << interval.lower_bound << " + " << interval.upper_bound << "}/2) = f(" << x << ") = " << y << std::endl;
			);
		}
		if (n > 1u)
		{
			if (y * f(interval.lower) < 0.0) {
				{
					PUSH_STEP_LOG(
						std::cout << "f(" << x << ").f(" << interval.lower << ") < " << 0 << std::endl;
						y * f(interval.upper_bound) > 0.0 ?
						std::cout << "f(" << x << ").f(" << interval.upper << ") > " << 0 << std::endl :
						std::cout << "f(" << x << ").f(" << interval.upper << ") < " << 0 << std::endl;
					)
				}
				return BisectionMethod(f, { interval.lower, x }, n - 1u);
			}
			if (y * f(interval.upper) < 0.0) {
				{
					PUSH_STEP_LOG(
						std::cout << "f(" << x << ").f(" << interval.upper << ") < " << 0 << std::endl;
						y * f(interval.lower_bound) > 0.0 ?
						std::cout << "f(" << x << ").f(" << interval.lower << ") > " << 0 << std::endl :
						std::cout << "f(" << x << ").f(" << interval.lower << ") < " << 0 << std::endl;
					)
				}
				return BisectionMethod(f, { x, interval.upper }, n - 1u);
			}
		}

		return { x, y };
	}

	Coord2D NewtonsMethod(func_t f, only<double> x0, uint16_t n)
	{
		return NewtonsMethod(f, derive(f), x0, n);
	}

	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, uint16_t n)
	{
		/*
		  df is assumed to be the explicitly stated derivative of f w.r.t it's input variable
		*/
		if (df(x0))
		{
			{
				PUSH_STEP_LOG(
					std::cout << std::endl << "x{next} = " << x0 << " - f(" << x0 << ")/f'(" << x0 << ")" << std::endl <<
					"x{next} = " << x0 << " - " << f(x0) << "/" << df(x0) << std::endl <<
					"x{next} = " << x0 << " - " << f(x0) / df(x0) << " = " <<
					x0 - f(x0) / df(x0) << std::endl;
				)
			}
			x0 -= f(x0) / df(x0);

			if (n > 1u)
				return NewtonsMethod(f, df, x0, n - 1u);
		}
		else {
			PUSH_STEP_LOG(
				std::cout << std::endl << "df(" << x0 << ") = 0" << std::endl;
			)
		}
		
		return { x0, f(x0) };
	}

	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, uint16_t n)
	{
		double F0 = f(x0), F1 = f(x1);
		if (F1 != F0)
		{
			{
				PUSH_STEP_LOG(
					std::cout << "x{next} = " << x1 << " - " << F1 << "(" << x1 << " - " << x0 << ")/(" << F1 << " - " << F0 << ")" << std::endl <<
					"x{next} = " << x1 << " - " << F1 << " * " << x1 - x0 << "/" << F1 - F0 << std::endl <<
					"x{next} = " << x1 << " - " << F1 * (x1 - x0) / (F1 - F0) << " = " << x1 - F1 * (x1 - x0) / (F1 - F0) << std::endl;
				)
			}
			double x2 = x1 - F1 * (x1 - x0) / (F1 - F0);

			if (n > 2u)
				return SecantMethod(f, x1, x2, n - 1u);

			return { x2, f(x2) };
		}
		else {
			PUSH_STEP_LOG(
				std::cout << std::endl << "f(" << x0 << ") = f(" << x1 << ") = " << F0 << std::endl;
			)
		}

		return { x1, F1 };
	}

	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, uint16_t n)
	{
		/*
		  Accepting in interval as { ai, bi }
		*/
		auto SMApprox = SecantMethod(f, interval.left, interval.right, 2u);

		if (SMApprox.dependent * f(interval.left) < 0.0) {       //f(x[i+1]) * f(ai) < 0: root probably lies in { ai, x[i+1] }
			{
				PUSH_STEP_LOG(
					std::cout << std::endl << "f(" << SMApprox.dependent << ").f(" << interval.left << ") < 0" << std::endl;
				if (SMApprox.dependent * f(interval.right) < 0.0)
					std::cout << "f(" << SMApprox.dependent << ").f(" << interval.right << ") < 0" << std::endl;
				else if (SMApprox.dependent * f(interval.right) > 0.0)
					std::cout << "f(" << SMApprox.dependent << ").f(" << interval.right << ") > 0" << std::endl;
				)
			}
			interval.upper_bound = SMApprox.independent;     //interval => { ai, x[i+1] }
		}
		else if (SMApprox.dependent * f(interval.right) < 0.0) { //f(x[i+1]) * f(bi) < 0: root probably lies in { x[i+1], bi }
			{
				PUSH_STEP_LOG(
					std::cout << std::endl << "f(" << SMApprox.dependent << ").f(" << interval.right << ") < 0" << std::endl;
				if (SMApprox.dependent * f(interval.left) < 0.0)
					std::cout << "f(" << SMApprox.dependent << ").f(" << interval.left << ") < 0" << std::endl;
				else if (SMApprox.dependent * f(interval.left) > 0.0)
					std::cout << "f(" << SMApprox.dependent << ").f(" << interval.left << ") > 0" << std::endl;
				)
			}
			interval.lower_bound = SMApprox.independent;     //interval => { x[i+1], bi }
		}
		else {
			{
				PUSH_STEP_LOG(
					std::cout << "f(" << SMApprox.independent << ") = " << SMApprox.dependent << " in (" << interval.lower << ", " << interval.upper << ")" << std::endl;
				)
			}
			return SMApprox;
		}
		{
			PUSH_STEP_LOG(
				std::cout << "f(" << SMApprox.independent << ") = " << SMApprox.dependent << " in (" << interval.lower << ", " << interval.upper << ")" << std::endl;
			)
		}

		if (n > 1u)
			return RegulaFalsiMethod(f, interval, n - 1u);

		return SMApprox;
	}

	Coord2D FixedPointMethod(func_t f, only<double> x0, uint16_t n)
	{
		double F0 = f(x0);

		{
			PUSH_STEP_LOG(
				std::cout << std::endl << "f(" << x0 << ") = " << F0 << std::endl;
			)
		}
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
		double last_x;
		double magnitude = 0.0;

		do {
			if (Approx.dependent * f(interval.left) < 0.0) {
				{
					PUSH_STEP_LOG(
						std::cout << "f(" << Approx.x << ").f(" << interval.lower << ") < " << 0 << std::endl;
     					Approx.y * f(interval.upper_bound) > 0.0 ?
						std::cout << "f(" << Approx.x << ").f(" << interval.upper << ") > " << 0 << std::endl :
						std::cout << "f(" << Approx.x << ").f(" << interval.upper << ") < " << 0 << std::endl;
					)
				}
				interval.upper_bound = Approx.independent;
			}
			else if (Approx.dependent * f(interval.right) < 0.0) {
				{
					PUSH_STEP_LOG(
						std::cout << "f(" << Approx.x << ").f(" << interval.upper << ") < " << 0 << std::endl;
					    Approx.y * f(interval.lower_bound) > 0.0 ?
						std::cout << "f(" << Approx.x << ").f(" << interval.lower << ") > " << 0 << std::endl :
						std::cout << "f(" << Approx.x << ").f(" << interval.lower << ") < " << 0 << std::endl;
					)
				}
				interval.lower_bound = Approx.independent;
				}
				else {
					{
						PUSH_STEP_LOG(
							std::cout << std::endl << "No apparent Error" << std::endl;
						)
					}
					return Approx;
				}

			last_x = Approx.independent;
			Approx = IterationBasedApproximation::BisectionMethod(f, interval, 1u);

			magnitude = -log10(abs(Approx.independent - last_x));
			{
				PUSH_STEP_LOG(
					Approx.independent == last_x ? std::cout << std::endl << "No apparent Error" << std::endl
					: std::cout << std::endl << "Error lies in the order of 10^(-" << ceil(magnitude) << ")" << std::endl;
				)
			}
		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D NewtonsMethod(func_t f, only<double> x0, only<double> accuracy)
	{
		return NewtonsMethod(f, derive(f), x0, accuracy);
	}

	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, only<double> accuracy)
	{
		Coord2D Approx = { x0 };
		double last_x;
		double magnitude = 0.0;

		do {
			last_x = Approx.independent;
			Approx = IterationBasedApproximation::NewtonsMethod(f, df, last_x, 1u);

			magnitude = -log10(abs(Approx.independent - last_x));
			{
				PUSH_STEP_LOG(
					Approx.independent == last_x ? std::cout << std::endl << "no apparent Error" << std::endl
					: std::cout << std::endl << "Error lies in the order of 10^(-" << ceil(magnitude) << ")" << std::endl;
				)
			}
		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, only<double> accuracy)
	{
		Coord2D Approx = { x1 };
		Interval<double> last_xs = { x0 };
		double magnitude = 0.0;

		do {
			last_xs.upper = Approx.x;
			Approx = IterationBasedApproximation::SecantMethod(f, Approx.x, last_xs.lower, 2u);
			last_xs.lower = last_xs.upper;

			magnitude = -log10(abs(Approx.independent - last_xs.upper));
			{
				PUSH_STEP_LOG(
					Approx.independent == last_xs.upper ? std::cout << std::endl << "no apparent Error" << std::endl
					: std::cout << std::endl << "Error lies in the order of 10^(-" << ceil(magnitude) << ")\n" << std::endl;
				)
			}
		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, only<double> accuracy)
	{
		Coord2D Approx = IterationBasedApproximation::SecantMethod(f, interval.left, interval.right, 2u);
		double last_x;
		double magnitude = 0.0;

		do {
			if (Approx.dependent * f(interval.left) < 0.0) {      //f(x[i+1]) * f(ai) < 0: root probably lies in { ai, x[i+1] }
				{
					PUSH_STEP_LOG(
						std::cout << "f(" << Approx.x << ").f(" << interval.left << ") < 0";
		    			if (Approx.dependent * f(interval.right) < 0.0)
	    					std::cout << std::endl << "f(" << Approx.x << ").f(" << interval.right << ") < 0\n" << std::endl;
			    	 	else if (Approx.dependent * f(interval.right) > 0.0)
			    			std::cout << "f(" << Approx.x << ").f(" << interval.right << ") > 0\n" << std::endl;
						else
							std::cout << std::endl;
					);
				}
				interval.upper_bound = Approx.independent;       //interval => { ai, x[i+1] }
			}
			else if (Approx.dependent * f(interval.right) < 0.0) {      //f(x[i+1]) * f(bi) < 0: root probably lies in { x[i+1], bi }
				{
					PUSH_STEP_LOG(
						std::cout << std::endl << "f(" << Approx.x << ").f(" << interval.right << ") < 0" << std::endl;
			     		if(Approx.dependent * f(interval.left) < 0.0)
			    			std::cout << "f(" << Approx.x << ").f(" << interval.left << ") < 0\n" << std::endl;
						else if(Approx.dependent * f(interval.left) > 0.0)
            			    std::cout << "f(" << Approx.x << ").f(" << interval.left << ") > 0\n" << std::endl;
						else 
							std::cout << std::endl;
					);
				}
				interval.lower_bound = Approx.independent;       //interval => { x[i+1], bi }
			}
			else {
				{
					PUSH_STEP_LOG(
						std::cout << std::endl << "no apparent Error" << std::endl;
					)
				}
				return Approx;
			}
			last_x = Approx.independent;
			Approx = IterationBasedApproximation::SecantMethod(f, interval.left, interval.right, 2u);

			magnitude = -log10(abs(Approx.independent - last_x));
			{
				PUSH_STEP_LOG(
					Approx.independent == last_x ? std::cout << std::endl << "no apparent Error" << std::endl
					: std::cout << std::endl << "Error lies in the order of 10^(-" << ceil(magnitude) << ")" << std::endl;
				)
			}
		} while (magnitude < accuracy);

		return Approx;
	}

	Coord2D FixedPointMethod(func_t f, only<double> x0, only<double> accuracy)
	{
		Coord2D Approx; Approx.dependent = x0;
		double magnitude = 0.0;

		do {
			Approx.independent = Approx.dependent;
			{
				PUSH_STEP_LOG(
					std::cout << std::endl << "f(" << Approx.independent << ") = " << f(Approx.independent) << std::endl;
				)
			}
			Approx.dependent = f(Approx.independent);

			magnitude = -log10(abs(Approx.dependent - Approx.independent));
			{
				PUSH_STEP_LOG(
					Approx.dependent == Approx.independent ? std::cout << std::endl << "no apparent Error" << std::endl
					:std::cout << std::endl << "Error lies in the order of 10^(-" << ceil(magnitude) << ")" << std::endl;
				)
			}
		} while (magnitude < accuracy);

		return Approx;
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

std::vector<double> DifferenceLadder(std::vector<double> Ypoints)
{
	{
		PUSH_STEP_LOG(
			std::cout << std::endl << "delta 0:{";

			for (int i = Ypoints.size() - 1; i >= 0; --i) {
				std::cout << "{Y" << i << ":" << Ypoints[i] << "}";
				if (i > 0) std::cout << ", ";
			}

			std::cout << "}";
		)
	}
	/*
		generating deltas of Y
	*/
	for (size_t i = 0; i < Ypoints.size() - 1; ++i) {
		{
			PUSH_STEP_LOG(
				std::cout << std::endl << "delta " << i + 1 << ":{";
			);
		}
		for (size_t j = Ypoints.size() - 1; j > i; --j) {
			Ypoints[j] -= Ypoints[j - 1];
			{
				PUSH_STEP_LOG(
					std::cout << "{Y" << j - i - 1 << ":" << Ypoints[j] << "}";
				    if (j > i + 1) std::cout << ", ";
				)
			}
		}
		{
			PUSH_STEP_LOG(
				std::cout << "}";
			)
		}
	}

	return Ypoints;
}

std::vector<double> DevidedDifferenceLadder(std::vector<Coord2D> points)
{
	std::vector<Interval<double>> neighbourhood(points.size());
	std::vector<double> ddiff(points.size());
	
	/*
	   collecting Ys and their bounded regions
	*/
	for (size_t i = 0; i < points.size(); ++i) {
		ddiff[i] = points[i].y;
		neighbourhood[i] = { points[i].x, points[i].x };
	}
	{
		PUSH_STEP_LOG(
			std::cout << std::endl << "delta 0:{";
		    
		    for (int i = points.size() - 1; i >= 0; --i) {
			    std::cout << "{Y" << i << ":" << points[i].y << "}";
			    if (i > 0) std::cout << ", ";
		    }
			
			std::cout << "}";
		)
	}
	/*
		generating devided differences of Y in their neighbourhood
	*/
	for (size_t i = 0; i < points.size() - 1; ++i)
	{
		{
			PUSH_STEP_LOG(
				std::cout << std::endl << "delta " << i + 1 << ":{";
			);
		}
		for (size_t j = points.size() - 1; j > i; --j) {
			ddiff[j] = (ddiff[j] - ddiff[j - 1]) / (neighbourhood[j].upper - neighbourhood[j - 1].lower);
			neighbourhood[j] = { neighbourhood[j - 1].lower, neighbourhood[j].upper };
			{
				PUSH_STEP_LOG(
					std::cout << "{Y" << j - i - 1 << ":" << ddiff[j] << "}";
    				if (j > i + 1) std::cout << ", ";
				)
			}
		}
		{
			PUSH_STEP_LOG(
				std::cout << "}";
			)
		}
	}

	return ddiff;
}

func_t LinearInterpolation(std::vector<Coord2D> points)
{
	std::function<Interval<double>(double)> get_interval_of =
		[=](double x) {
		
		Interval<double> I = { 0, 0 };

		if (x < points.back().x)
			return I;

		for (size_t i = 0; i < points.size() - 1; ++i)
		{
			I = { points[i].x, points[i + 1].x };

			if (I[x])
				return I;
		}

		return I;
	};

	return [=](double x) { 
		auto subsection = get_interval_of(x);
		return LinearInterpolate(subsection.lower_bound, subsection.upper_bound, remainder((x - subsection.lower_bound) / subsection, 1.0));
	};
}

func_t LinearInterpolation(Coord2DMap& map, const Interval<double>& interval)
{
	/*
	  Linearly interpolates between each consecutive node in map
	*/
	return LinearInterpolation(map, interval, map.size());
}

func_t LinearInterpolation(Coord2DMap& map, const Interval<double>& interval, uint32_t n)
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

func_t LinearInterpolation(Coord2DMap& map, const Interval<double>& interval, uint16_t poly_factor, uint32_t disloc_factor)
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

polynomial::poly_func<double> NewtonsInterpolation(std::vector<double> points, Interval<double> input_domain) {
	auto Y = DifferenceLadder(points);

	double scale = input_domain / (Y.size() - 1.0), 
		   bias = input_domain.lower_bound;
	polynomial::poly_func<double> poly = { Y.front() };
	{
		PUSH_STEP_LOG(
			std::cout << std::endl << "\nterm0 = " << poly.String("x") << std::endl;
		)
	}
	polynomial::poly_func<double> term = { 1.0 };
	for (size_t i = 1u; i < Y.size(); ++i)
	{
		{
			PUSH_STEP_LOG(
				std::cout << std::endl << "term" << i << " = (" << Y[i] << "/" << gamma<double>(i) << ")(" << term.String("x") << ")(" << polynomial::FromRoots({ (i - 1.0) * scale + bias }).String("x") << ")" <<
				std::endl << "term" << i << " = " << Y[i] / gamma<double>(i) << "(" << term.AndRoots({ (i - 1.0) * scale + bias }).String("x") << ")" << std::endl;
			)
		}
		term.addRoots({ (i - 1.0) * scale + bias });
		{
			PUSH_STEP_LOG(
				std::cout << std::endl << "p(x) = (" << poly.String("x") << ") + (" << (term * (Y[i]/gamma<double>(i))).String("x") << ")" << std::endl <<
				"p(x) = " << (poly + term * (Y[i]/gamma<double>(i))).String("x") << std::endl;
			)
		}
		poly = poly + term * (Y[i] / gamma<double>(i));
	}

	return poly;
}

polynomial::poly_func<double> NewtonsInterpolation(std::vector<Coord2D> points) {
	auto Y = DevidedDifferenceLadder(points);

	polynomial::poly_func<double> poly = { Y.front() };
	{
		PUSH_STEP_LOG(
			std::cout << std::endl << "\np(x) = " << poly.String("x") << std::endl;
		);
	}

	auto term  = polynomial::one<double>();

	for (size_t i = 1u; i < Y.size(); ++i)
	{
		{
			PUSH_STEP_LOG(
				std::cout << "term" << i << " = " << Y[i] << "(" << term.String<double>("x") << ") * (" << polynomial::FromRoots({ points[i - 1u].x }).String("x") << ")" << std::endl <<
				"term" << i << " = " << Y[i] << "(" << term.AndRoots({ points[i - 1u].x }).String("x") << ")" << std::endl;
			)
		}
		term.addRoots({ points[i - 1u].x });
		{
			PUSH_STEP_LOG(
				std::cout << "p(x) = (" << poly.String("x") << ") + (" << (term * (Y[i]/gamma<double>(i))).String("x") << ")" << std::endl <<
				"p(x) = " << (poly + term * (Y[i]/gamma<double>(i))).String("x") << std::endl;
			)
		}
		poly = poly + term * Y[i];
	}

	return poly;
}

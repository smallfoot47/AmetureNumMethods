#pragma once
#define ENABLE_STEP_TRACKING
#include "Skintainers.h"
#include "Containers.h"
#include "FuncShorts.h"
#include "Func++.h"

/*
	  Limited Iterative Models: returns { x, f(x) }. Process is repeated a given number of times
*/
namespace IterationBasedApproximation
{
	/*
	  Closest root within iterative limit
	  -Accepts {function}, {interval to be bisected}, {number of iterations left}
	  [minimum n = 1]
	*/
	Coord2D BisectionMethod(func_t f, Interval<double> interval, uint16_t n);
	/*
	  Approximate Newton's Method for approaching roots within iterative limit
	  (numerically estimates the derivative at a point)
	  -Accepts {function}, {starting point}, {number of iterations left}
	  [minimum n = 1]
	*/
	Coord2D NewtonsMethod(func_t f, only<double> x0, uint16_t n);
	/*
	  Newton's Method for approaching roots within iterative limit
	  -Accepts {function}, {function derivative}, {starting point}, {number of iterations left}
	  [minimum n = 1]
	*/
	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, uint16_t n);
	/*
	  Closest root within iterative limit
	  -Accepts {function}, {starting left value}, {starting right value}, {number of iterations left}
	  [minimum n = 2]
	*/
	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, uint16_t n);
	/*
	  Closest root within iterative limit
	  -Accepts {function}, {initial interval}, {number of iterations left}
	  [minimum n = 1]
	*/
	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, uint16_t n);
	/*
	  Closest root within iterative limit
	  -Accepts {function}, {starting point}, {number of iterations left}
	  [minimum n = 1]
	*/
	Coord2D FixedPointMethod(func_t f, only<double> x0, uint16_t n);
}

/*
	  Error-bound Iterative Models: returns a function of the magnitude of maximum error to the base 10 which returns { x, f(x) }
	 [
	  Here error is vaugely used. A more appropriate wording would be variation as these functions measure the variation in outputs per iteration.
	 ]
*/
namespace ErrorBasedApproximation
{
	/*
	  Closest root within given error
	  -Accepts {function}, {interval to be bisected}, {order in base 1/10 of maximum error}
	  [minimum n = 1]
	*/
	Coord2D BisectionMethod(func_t f, Interval<double> interval, only<double> accuracy);
	/*
	  Approximate Newton's Method for approaching roots within given error
	  (numerically estimates the derivative at a point)
	  -Accepts {function}, {starting point}, {order in base 1/10 of maximum error}
	  [minimum n = 1]
	*/
	Coord2D NewtonsMethod(func_t f, only<double> x0, only<double> accuracy);
	/*
	  Newton's Method for approaching roots within given error
	  -Accepts {function}, {function derivative}, {starting point}, {order in base 1/10 of maximum error}
	  [minimum n = 1]
	*/
	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, only<double> accuracy);
	/*
	  Closest root within given error
	  -Accepts {function}, {starting left value}, {starting right value}, {order in base 1/10 of maximum error}
	  [minimum n = 2]
	*/
	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, only<double> accuracy);
	/*
	  Closest root within given error
	  -Accepts {function}, {initial interval}, {order in base 1/10 of maximum error}
	  [minimum n = 1]
	*/
	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, only<double> accuracy);
	/*
	  Closest root within given error
	  -Accepts {function}, {starting point}, {order in base 1/10 of maximum error}
	  [minimum n = 1]
	*/
	Coord2D FixedPointMethod(func_t f, only<double> x0, only<double> accuracy);
}

//converts a bounded function to a Coord2DMap, disregarding memory optimization
Coord2DMap Function_to_Map(func_t f, Interval<double> interval);
//converts a bounded function to a Coord2DMap, taking n discrete points from f
Coord2DMap Function_to_Map(func_t f, Interval<double> interval, uint32_t n);
//converts a bounded function to a Coord2DMap, seperated by epsilon amoung adjacent points along X
Coord2DMap Function_to_Map(func_t f, Interval<double> interval, double epsilon);

/*
  generates a row of differences(in delta form) for first point
  -Accepts {collection of Y values}
*/
std::vector<double> DifferenceLadder(std::vector<double> points);
/*
  generates a row of devided differences(in delta form) for first point
  -Accepts {collection of Coord2Ds(x, y)}
*/
std::vector<double> DevidedDifferenceLadder(std::vector<Coord2D> points);

//makes mapping linearly continous
func_t LinearInterpolation(std::vector<Coord2D> points);
//makes mapping linearly continous
func_t LinearInterpolation(Coord2DMap& map, const Interval<double>& interval);
//makes mapping linearly continous with n nodes
func_t LinearInterpolation(Coord2DMap& map, const Interval<double>& interval, uint32_t n);
//makes mapping linearly continous with vertice-count controlled by poly_factor and input bias controlled by disloc_factor
func_t LinearInterpolation(Coord2DMap& map, const Interval<double>& interval, uint16_t poly_factor, uint32_t disloc_factor = 0);

/*
  estimates a polynomial that passes through the given points
  provided the X values are equally spaced
  -Accepts {Y values per each point}, {open interval where X values lie(Ascending order)}
*/
polynomial::poly_func<double> NewtonsInterpolation(std::vector<double> points, Interval<double> input_domain);
/*
  estimates a polynomial that passes through the given points
  -Accepts {collection of Coord2Ds(x, y) (Ascending order of x)}
*/
polynomial::poly_func<double> NewtonsInterpolation(std::vector<Coord2D> points);
/*
  estimates a polynomial that passes through the given points
  -Accepts {collection of Coord2Ds(x, y) (Ascending order of x)}
*/
polynomial::poly_func<double> LegrangesInterpolation(std::vector<Coord2D> points);

/*
converts from [0, 1] space to [foot, head] space linearly
*/
template <class V_t, class P_t>
inline V_t LinearInterpolate(V_t foot, V_t head, P_t percent)
{
	return foot + (head - foot) * clamp(percent);
}

/*
  converts from (-inf, inf) space to [foot, head] space wrt the sigmoid function
*/
template <class V_t, class P_t>
inline V_t SigmoidInterpolate(V_t foot, V_t head, P_t affinity){
	return foot + (head - foot) * sigmoid(affinity);
}

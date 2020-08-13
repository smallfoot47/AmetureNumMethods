#pragma once
#include <vector>

#include "Skintainers.h"
#include "Containers.h"
#include "FuncShorts.h"
#include "Func++.h"

namespace IterationBasedApproximation
{
	/*
      Limited Iterative Models: returns { x, f(x) }. Process is repeated a given number of times
    */
	Coord2D BisectionMethod(func_t f, Interval<double> interval, uint8_t n);//Closest root within iterative limit: starts at n = 1
	Coord2D NewtonsMethod(func_t f, only<double> x0, uint8_t n);//Approximate Newton's Method for approaching roots within iterative limit: starts at n = 1
	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, uint8_t n);//Newton's Method for approaching roots within iterative limit: starts at n = 1
	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, uint8_t n);//Closest root within iterative limit: starts at n = 2
	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, uint8_t n);//Closest root within iterative limit: starts at n = 1
}

namespace ErrorBasedApproximation
{
	/*
      Error-bound Iterative Models: returns a function of the magnitude of maximum error to the base 10 which returns { x, f(x) }
     [
      Here error is vaugely used. A more appropriate wording would be variation as these functions measure the variation in outputs per iteration.
     ]
    */
	Coord2D BisectionMethod(func_t f, Interval<double> interval, only<double> accuracy);//Closest root within given error
	Coord2D NewtonsMethod(func_t f, only<double> x0, only<double> accuracy);//Approximate Newton's Method for approaching roots within given error
	Coord2D NewtonsMethod(func_t f, func_t df, only<double> x0, only<double> accuracy);//Newton's Method for approaching roots within given error
	Coord2D SecantMethod(func_t f, only<double> x0, only<double> x1, only<double> accuracy);//Closest root within given error
	Coord2D RegulaFalsiMethod(func_t f, Interval<double> interval, only<double> accuracy);//Closest root within given error
}

Coord2DMap Function_to_Map(func_t f, Interval<double> interval);//converts a bounded function to a Coord2DMap, disregarding memory optimization
Coord2DMap Function_to_Map(func_t f, Interval<double> interval, uint32_t n);//converts a bounded function to a Coord2DMap, taking n discrete points from f
Coord2DMap Function_to_Map(func_t f, Interval<double> interval, double epsilon);//converts a bounded function to a Coord2DMap, seperated by epsilon amoung adjacent points along X

func_t Interpolate(Coord2DMap& map, const Interval<double>& interval);//makes mapping linearly continous
func_t Interpolate(Coord2DMap& map, const Interval<double>& interval, uint32_t n);//makes mapping linearly continous with n nodes
func_t Interpolate(Coord2DMap& map, const Interval<double>& interval, uint16_t poly_factor, uint32_t disloc_factor = 0);//makes mapping linearly continous with vertice-count controlled by poly_factor and input bias controlled by disloc_factor

template <class V_t, class P_t>
inline V_t LinearInterpolate(V_t foot, V_t head, P_t percent)//converts from [0, 1] space to [foot, head] space linearly
{
	return foot + (head - foot) * percent;
}

template <class V_t, class P_t>
inline V_t SigmoidInterpolate(V_t foot, V_t head, P_t affinity)//converts from (-inf, inf) space to [foot, head] space wrt the sigmoid function
{
	return foot + (head - foot) * sigmoid(affinity);
}

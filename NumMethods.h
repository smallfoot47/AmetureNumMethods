#pragma once
#include <vector>

#include "Containers.h"
#include "Func++.h"

Coord2D BisectionMethod(std::function<double(double)> f, Interval<double> interval, uint8_t n);//Closest zero within iterative limits
inline Coord2D NewtonsMethod(std::function<double(double)> f, double x0, uint8_t n);//Approximate Newton's Method for approaching roots
Coord2D NewtonsMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, uint8_t n);//Explicitly defined Newton's Method for approaching roots

inline Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval);//converts a bounded function to a Coord2DMap, disregarding memory optimization
Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, uint32_t n);//converts a bounded function to a Coord2DMap, taking n discrete points from f
Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, double epsilon);//converts a bounded function to a Coord2DMap, seperated by epsilon amoung adjacent points along X

std::function<double(double)> Interpolate(Coord2DMap map, Interval<double> interval);//makes mapping linearly continous
std::function<double(double)> Interpolate(Coord2DMap map, Interval<double> interval, uint8_t poly_factor, uint32_t disloc_factor = 0);//makes mapping linearly continous with vertice-count controlled by poly_factor and input bias controlled by disloc_factor

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

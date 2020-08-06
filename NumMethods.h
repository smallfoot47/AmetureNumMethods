#pragma once
#include <functional>
#include <vector>
#include <unordered_map>

#include "Containers.h"

Coord2D IntermediateValue(std::function<double(double)> f, Interval<double> interval, uint8_t n);//Closest zero within iterative limits
inline Coord2D NewtonsMethod(std::function<double(double)> f, double x0, uint8_t n);//Approximate Newton's Method for approaching roots
Coord2D NewtonsMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, uint8_t n);//Explicitly defined Newton's Method for approaching roots

inline Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval);//converts a bounded function to a Coord2DMap, disregarding memory optimization
Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, uint32_t n);//converts a bounded function to a Coord2DMap, taking n discrete points from f
Coord2DMap Function_to_Map(std::function<double(double)> f, Interval<double> interval, double epsilon);//converts a bounded function to a Coord2DMap, seperated by epsilon amoung adjacent points along X

template <class V_t, class P_t>
V_t Interpolate(V_t foot, V_t head, P_t percent);

std::function<double(double)> Interpolate(Coord2DMap map, Interval<double> interval);//makes mapping linearly continous
std::function<double(double)> Interpolate(Coord2DMap map, Interval<double> interval, uint8_t poly_factor, uint32_t disloc_factor = 0);//makes mapping linearly continous with vertice-count controlled by poly_factor and input bias controlled by disloc_factor

template<class V_t, class P_t>
inline V_t Interpolate(V_t foot, V_t head, P_t percent)
{
	return foot + (head - foot) * percent;
}

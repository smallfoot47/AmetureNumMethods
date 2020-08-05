#pragma once
#include <functional>
#include <vector>
#include <unordered_map>

typedef std::unordered_map<double, double> Coord2DMap;

struct Coord2D //dual coordinate containment class
{
	union { double x, h, horizontal, right, abscissa, first; };
	union { double y, k, vertical, up, ordinate, second; };

	operator bool()
	{
		return x && y;
	}

	operator double()
	{
		return sqrt(x * x + y * y);
	}
	
	operator float()
	{
		float X = (float)x;
		float Y = (float)y;

		return sqrtf(X * X + Y * Y);
	}
	
	operator long double()
	{
		long double X = (long double)x;
		long double Y = (long double)y;

		return sqrtl(X * X + Y * Y);
	}
};

Coord2D IntermediateValue(std::function<double(double)> f, std::pair<double, double> interval, uint8_t n);//Closest zero within iterative limits
inline Coord2D NewtonsMethod(std::function<double(double)> f, double x0, uint8_t n);//Approximate Newton's Method for approaching roots
Coord2D NewtonsMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, uint8_t n);//Explicitly defined Newton's Method for approaching roots

Coord2DMap Function_to_Map(std::function<double(double)> f, std::pair<double, double> interval);//converts a bounded function to a Coord2DMap, disregarding memory optimization
Coord2DMap Function_to_Map(std::function<double(double)> f, std::pair<double, double> interval, uint32_t n);//converts a bounded function to a Coord2DMap, taking n discrete points from f
Coord2DMap Function_to_Map(std::function<double(double)> f, std::pair<double, double> interval, double epsilon);//converts a bounded function to a Coord2DMap, seperated by epsilon amoung adjacent points along X

std::function<double(double)> Interpolate(Coord2DMap map, std::pair<double, double> interval);//makes mapping linearly continous
std::function<double(double)> Interpolate(Coord2DMap map, std::pair<double, double> interval, uint8_t poly_factor, uint32_t disloc_factor = 0);//makes mapping linearly continous with vertice-count controlled by poly_factor and input bias controlled by disloc_factor
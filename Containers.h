#pragma once
#include <unordered_map>
#include <math.h>

typedef std::unordered_map<double, double> Coord2DMap;

struct Coord2D //dual coordinate containment class
{
	union { double x, h, horizontal, independent, abscissa, right,  first; };
	union { double y, k, vertical, dependent, ordinate, up, second; };

	operator bool()
	{
		return x && y;
	}

	operator double()
	{
		return sqrt(x * x + y * y);
	}
	
	Coord2D const operator +(Coord2D C2D)
	{
		return { x + C2D.x, y + C2D.y };
	}
	
	Coord2D const operator -(Coord2D C2D)
	{
		return { x - C2D.x, y - C2D.y };
	}
	
	void operator +=(Coord2D C2D)
	{
		x += C2D.x;
		y += C2D.y;
	}
	
	void operator -=(Coord2D C2D)
	{
		x -= C2D.x;
		y -= C2D.y;
	}

	Coord2D const operator *(double scale)
	{
		return { x * scale, y * scale };
	}

	Coord2D const operator /(double scale)
	{
		return { x / scale, y / scale };
	}
	
	void operator *=(double scale)
	{
		x *= scale;
		y *= scale;
	}

	void operator /=(double scale)
	{
		x *= scale;
		y *= scale;
	}
};

template <class D_t>
struct Interval
{
	union { D_t left, lower, lower_bound, low, min, first; };
	union { D_t right, upper, upper_bound, up, high, max, second; };

	template <class T_t>
	operator Interval<T_t>()
	{
		return { (T_t)left, (T_t)right };
	}

	bool operator ()(D_t x)//closed interval
	{
		return left < x && x < right;
	}

	bool operator [](D_t x)//open interval
	{
		return left <= x && x <= right;
	}

	bool operator <<(D_t x)//left closed interval
	{
		return left < x && x <= right;
	}

	bool operator >>(D_t x)//right closed interval
	{
		return left <= x && x < right;
	}

	bool operator <(D_t x)//beyond upper bound
	{
		return x > right;
	}

	bool operator <=(D_t x)//beyond or on upper bound
	{
		return x >= right;
	}

	bool operator >(D_t x)//below lower bound
	{
		return x > left;
	}

	bool operator >=(D_t x)//below or on lower bound
	{
		return x >= left;
	}

	operator D_t()
	{
		return right - left;
	}
};

template <class D_t> bool operator >(D_t x, Interval<D_t> I)
{
	return I < x;
}

template <class D_t> bool operator >=(D_t x, Interval<D_t> I)
{
	return I <= x;
}

template <class D_t> bool operator <(D_t x, Interval<D_t> I)
{
	return I > x;
}

template <class D_t> bool operator <=(D_t x, Interval<D_t> I)
{
	return I >= x;
}

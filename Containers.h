#pragma once
#include <unordered_map>
#include <math.h>

typedef std::unordered_map<double, double> Coord2DMap;

struct Coord2D //dual coordinate containment class
{
	union { double x, h, horizontal, independent, abscissa, right,  first; };
	union { double y, k, vertical, dependent, ordinate, up, second; };

	operator bool() const
	{
		return x && y;
	}

	operator double () const
	{
		return sqrt(x * x + y * y);
	}
	
	Coord2D operator +(const Coord2D C2D) const
	{
		return { x + C2D.x, y + C2D.y };
	}
	
	Coord2D operator -(const Coord2D C2D) const
	{
		return { x - C2D.x, y - C2D.y };
	}
	
	void operator +=(const Coord2D C2D)
	{
		x += C2D.x;
		y += C2D.y;
	}
	
	void operator -=(const Coord2D C2D)
	{
		x -= C2D.x;
		y -= C2D.y;
	}

	Coord2D operator *(const double scale) const
	{
		return { x * scale, y * scale };
	}

	Coord2D operator /(const double scale) const
	{
		return { x / scale, y / scale };
	}
	
	void operator *=(const double scale)
	{
		x *= scale;
		y *= scale;
	}

	void operator /=(const double scale)
	{
		x *= scale;
		y *= scale;
	}
};

template <class D_t>
struct Interval
{
	union { D_t left, lower, lower_bound, low, first; };
	union { D_t right, upper, upper_bound, up, high, second; };

	D_t mean() const
	{
		return (low + high) / 2.0;
	}

	template <class T_t>
	operator Interval<T_t>() const
	{
		return { (T_t)left, (T_t)right };
	}
	
	bool operator ()(const D_t& x) const//closed interval
	{
		return left < x && x < right;
	}

	bool operator [](const D_t& x) const//open interval
	{
		return left <= x && x <= right;
	}

	bool operator << (const D_t& x) const//left closed interval
	{
		return left < x && x <= right;
	}

	bool operator >> (const D_t& x) const//right closed interval
	{
		return left <= x && x < right;
	}

	bool operator < (const D_t& x) const//beyond upper bound
	{
		return x > right;
	}

	bool operator <= (const D_t& x) const//beyond or on upper bound
	{
		return x >= right;
	}

	bool operator > (const D_t& x) const//below lower bound
	{
		return x > left;
	}

	bool operator >= (const D_t& x) const//below or on lower bound
	{
		return x >= left;
	}

	operator D_t() const
	{
		return right - left;
	}
};

template <class D_t> bool operator >(const D_t& x, const Interval<D_t>& I)
{
	return I < x;
}

template <class D_t> bool operator >=(const D_t& x, const Interval<D_t>& I)
{
	return I <= x;
}

template <class D_t> bool operator <(const D_t& x, const Interval<D_t>& I)
{
	return I > x;
}

template <class D_t> bool operator <=(const D_t& x, const Interval<D_t>& I)
{
	return I >= x;
}

template <class D_t> Interval<D_t> Ascending(const Interval<D_t> I)
{
	D_t min = I.left < I.right ? I.left : I.right;
	return { min, I.left + I.right - min };
}

template <class D_t> Interval<D_t> Descending(const Interval<D_t> I)
{
	D_t min = I.left < I.right ? I.left : I.right;
	return { I.left + I.right - min, min };
}

#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include <string>

/*
   Makes a lambda with double x as a parameter and gives it the defenition used
*/
#define MAKE_FUNC(def) [&](double x) { def } 
/*
   Makes an initializer_list with value and function(value) as elements
*/
#define MAKE_FUNC_COORD(function, value) { value, function(value) }
/*
   Prints "<function name>(<value name>) = function(value)" and then moves to the cursor to the next line
*/
#define PRINT_AT_VAR(function, value) std::cout << #function << "(" << #value << ") = " << function(value) << std::endl
/*
   Prints "<function name>(value) = function(value)" and then moves the cursor to the next line
*/
#define PRINT_AT_VALUE(function, value) std::cout << #function << "(" << value  << ") = " << function(value) << std::endl
/*
   Prints "<container name>(<container>.first) = <container>.second" and then moves the cursor to the next line
*/
#define PRINT_FUNC_COORD(answer) std::cout << #answer << "(" << answer.first << ") = " << answer.second << std::endl

typedef std::function<double(double)> func_t;//double->double function storage object

namespace polynomial
{
	template <class D_t = double>
	struct poly_func
	{
	private:
		//coefficents in decreasing order of significance
		std::vector<D_t> Coeffs;

	public:
		//constructors----------------------------------------
		/*
		  default to p(x) = 0
		*/
		poly_func() {};
		/*
		   constructs a polynomail container wrt the coefficients in order of significance
		*/
		poly_func(const std::vector<D_t> coefficients) : Coeffs(coefficients) {}
		/*
		   constructs a polynomail container wrt the coefficients in order of significance
		*/
		poly_func(const std::initializer_list<D_t> coefficents) : Coeffs(coefficents) {}
		//getters---------------------------------------------
		/*
		  returns the degree of the polynomial stored in memory
		*/
		size_t Degree() const 
		{ 
			return Coeffs.empty() ? 0u : Coeffs.size() - 1u; 
		}
		/*
		   returns a reference to the coefficient of power: power % (max degree + 1)
		   NOTE: -1 corresponds to the most significant coefficient
		*/
		D_t& Coefficient(const int power) 
		{
			int cap = Coeffs.size();
			return Coeffs[(cap - 1 - power) % cap];
		}
		/*
		   returns a constant reference to the coefficient of power: power % (max degree + 1)
		   NOTE: -1 corresponds to the most significant coefficient
		*/
		const D_t& Coefficient(const int power) const
		{
			int cap = Coeffs.size();
			return Coeffs[(cap - 1 - power) % cap];
		}
		/*
		   returns function container that estimates the value of the polynomail at a point
		   return type can be varied as such: <polynomial>.Function<<new type>>() [ default type : polynomial type ]
		*/
		template<class R_t = D_t>
		std::function<R_t(D_t)> Function() const;
		/*
		  Normalized: most significant coefficient is 1 OR p(x) = 0
		*/
		bool isNormalized() const { 
			if (*this)
				return Coeffs.front() == 1;
			return true;
		}
		/*
		  Simplified: least significant term is non-zero OR p(x) = 0
		*/
		bool isSimplified() const { 
			if (*this)
				return Coeffs.back() != 0;
			return true;
		}
		/*
		  returns the normalized form of this polynomial
		*/
		poly_func Normalized() const {
			poly_func<D_t> p(*this);
			p.normalize();
			return p;
		}
		/*
		  returns the simplified form of this polynomial
        */
		poly_func Simplified() const {
			poly_func<D_t> p(*this);
			p.simplify();
			return p;
		}
		/*
		  returns the polynomial passing through the roots of this polynomial and the roots specified
		*/
		poly_func AndRoots(const std::vector<D_t> roots) const;
		/*
		   returns the polynomial passing through the roots of this polynomial and the roots specified
		*/
		poly_func AndRoots(const std::initializer_list<D_t> roots) const;
		//mutators--------------------------------------------
		/*
		   scales the polynomail such that the most significant coefficient is 1 if p(x) is non zero
		*/
		void normalize();
		/*
		   reduces the polynomial such that the least significant coefficient in non zero if p(x) is non zero
		*/
		void simplify();
		/*
		  multiplies this polynomial by (1 - root) for each root
		*/
		void addRoots(const std::vector<D_t> roots);
		/*
		  multiplies this polynomial by (1 - root) for each root
		*/
		void addRoots(const std::initializer_list<D_t> roots);
		//type-conversions------------------------------------
		operator std::function<D_t(D_t)>() const;
		operator std::string() const;
		/*
		  Converts to the form a0 x^n + a1 x^n-1 ...
		  where x will be replaced by variableName
		*/
		template<class R_t = D_t>
		std::string String(const std::string variableName) const;
		operator bool() const {
			if (Coeffs.empty())
				return false;

			for (auto itr = Coeffs.begin(); itr != Coeffs.end(); ++itr)
				if (*itr) return true;
			return false;
		}
		//operators-------------------------------------------
		/*
		  polynomial at value
		*/
		D_t operator()(const D_t x) const;
		poly_func operator +(const poly_func p) const {
			if (!*this)
				return p;
			if (!p)
				return *this;

			std::vector<D_t> coeffs(std::max<size_t>(p.Coeffs.size(), Coeffs.size()));
			int delta = this->Degree() - p.Degree();
			
			auto itrSum = coeffs.begin(); //coeffs iterator
			auto itrThis = Coeffs.begin(); //Coeffs iterator
			auto itrP = p.Coeffs.begin(); //p.Coeffs iterator

			if (delta >= 0)
			{
				auto blockThis = itrThis + delta;
				while (itrThis != blockThis){
					*itrSum++ = *itrThis++;
				}

				while (itrP != p.Coeffs.end()) {
					*itrSum++ = *itrThis++ + *itrP++;
				}
			}
			else
			{
				auto blockP = itrP - delta;
				while (itrP != blockP) {
					*itrSum++ = *itrP++;
				}

				while (itrThis != Coeffs.end()) {
					*itrSum++ = *itrThis++ + *itrP++;
				}
			}

			return poly_func(coeffs);
		}
		poly_func operator -(const poly_func p) const {
			if (!*this)
				return p;
			if (!p)
				return *this;

			std::vector<D_t> coeffs(std::max<size_t>(p.Coeffs.size(), Coeffs.size()));
			int delta = this->Degree() - p.Degree();

			auto itrSum = coeffs.begin(); //coeffs iterator
			auto itrThis = Coeffs.begin(); //Coeffs iterator
			auto itrP = p.Coeffs.begin(); //p.Coeffs iterator

			if (delta >= 0)//*this has more terms
			{
				auto blockThis = itrThis + delta;
				while (itrThis != blockThis) {
					*itrSum++ = *itrThis++;
				}

				while (itrP != p.Coeffs.end()) {
					*itrSum++ = *itrThis++ - *itrP++;
				}
			}
			else//p has more terms
			{
				auto blockP = itrP - delta;
				while (itrP != blockP) {
					*itrSum++ = -*itrP++;
				}

				while (itrThis != Coeffs.end()) {
					*itrSum++ = *itrThis++ - *itrP++;
				}
			}

			return poly_func(coeffs);
		}
		poly_func operator *(const poly_func p) const {
			std::vector<D_t> coeffs(Degree() + p.Degree() + 1u);//multiplied product will have a degree equal to the sum of degrees of the operands

			auto Coefficient_itr = coeffs.begin();//for writing into coeffs faster

			auto anchor_This = Coeffs.begin();//most significant coefficient under the constraint of constant degree
		    auto anchor_P = p.Coeffs.begin(); //least significant coefficient under the constraint of constant degree

			const auto stop_This = Coeffs.end() - 1;
			const auto stop_P = p.Coeffs.begin();

			do {
				auto slide_This = anchor_This;//decreases from anchor_This as much as possible under the contraint of constant degree
				auto slide_P = anchor_P;//increases from anchor_P as much as possible under the contraint of constant degree

				//linearly picks out from Coeffs and coeffs such that the degree is contant
				while (slide_P != stop_P && slide_This != stop_This) {
					*Coefficient_itr += (*slide_This++) * (*slide_P--);
				};
				//sums in the excluded boundary term
				*Coefficient_itr += (*slide_This) * (*slide_P);

				if (anchor_P + 1 != p.Coeffs.end())
					++anchor_P;//reduce significance of least significant by 1
				else
					++anchor_This;//reduce significance of mosst significant by 1
			} while (++Coefficient_itr != coeffs.end());

			return poly_func(coeffs);
		}
		poly_func operator /(const poly_func p) const {
			size_t Deg[] = { this->Degree(), p.Degree() };
			if (Deg[0] < Deg[1]) return {};

			std::vector<D_t> coeffs(Deg[0] - Deg[1] + 1u);

			auto itrDiv = coeffs.begin();
			for (auto itr = Coeffs.begin(); itrDiv != coeffs.end(); ++itr) {
				*itrDiv = *itr;
				/*
				  evaluates for x ^ { delta - (itrDiv - coeffs.begin()) } coefficient
				*/
				{
					auto  PrevItrDiv = coeffs.begin(); //to access all previous divident terms upto but not including current term
					auto itrP = p.Coeffs.begin() + (itrDiv - PrevItrDiv); //positon in 'p' polynomail

					while (PrevItrDiv != itrDiv) {
						*itrDiv -= *(PrevItrDiv++) * *(itrP--);// coeff[i] * p.Coeff[current_term_index - i]
					}
				}
				*(itrDiv++) /= p.Coeffs.front();
			}

			return poly_func(coeffs);
		}
		poly_func operator %(const poly_func p) const {
			return *this - ((*this / p) * p);
		}
		poly_func operator +(const D_t constant) const {
			poly_func p(*this);
			if (p.Coeffs.empty()) p.Coeffs.push_back(0);
			p.Coeffs.back() += constant;

			return p;
		}
		poly_func operator -(const D_t constant) const {
			poly_func p(*this);
			p.Coeffs.back() -= constant;

			return p;
		}
		poly_func operator *(const D_t scale) const {
			if (!*this)
				return {};

			std::vector<D_t> coeffs(Coeffs.size());
			
			auto itrScaled = coeffs.begin();
			auto itrThis = Coeffs.begin();
			while (itrThis != Coeffs.end()) {
				*itrScaled++ = *itrThis++ * scale;
			}

			return poly_func(coeffs);
		}
		poly_func operator /(const D_t scale) const {
			std::vector<D_t> coeffs(Coeffs.size());

			auto itrScaled = coeffs.begin();

			for (auto itrThis = Coeffs.begin(); itrThis != Coeffs.end(); ++itrThis)
				*(itrScaled++) = *itrThis / scale;

			return poly_func(coeffs);
		}
		poly_func operator %(const D_t divident) const {
			std::vector<D_t> coeffs(Coeffs.size());

			auto itrScaled = coeffs.begin();

			for (auto itrThis = Coeffs.begin(); itrThis != Coeffs.end(); ++itrThis)
				*(itrScaled++) = remainder(*itrThis, divident);

			return poly_func(coeffs);
		}
		poly_func derivative() const {
			std::vector<D_t> coeffs(Degree());

			auto itrP = Coeffs.rbegin() + 1;

			size_t n = 1;
			for (auto itrD = coeffs.rbegin(); itrD != coeffs.rend(); ++itrD)
				*itrD = n++ * *itrP++;

			return poly_func(coeffs);
		}
		poly_func antiderivative(const D_t constant = 0) const {
			std::vector<D_t> coeffs(Coeffs.size() + 1);
			coeffs.back() = constant;

			auto itrP = Coeffs.rbegin();

			size_t n = 1;
			for (auto itrD = coeffs.rbegin() + 1; itrD != coeffs.rend(); ++itrD)
				*itrD = n++ * *itrP++;

			return poly_func(coeffs);
		}
		void operator +=(const poly_func p) {	
			if (!p)
				return;
			else if (!(*this)) {
				this->Coeffs = p.Coeffs;
				return;
			}

			auto itrThis = Coeffs.begin(); //Coeffs iterator
			auto itrP = p.Coeffs.begin(); //p iterator

			int delta = this->Degree() - p.Degree();
			if (delta >= 0)
				itrThis += delta;
			else {
				auto Pfront = p.Coeffs.begin();
				itrP = itrP - delta + 1;
				itrThis = Coeffs.insert(itrThis, Pfront, itrP) - delta;
			}

			if(itrP != p.Coeffs.end())
			while (itrThis != Coeffs.end()) {
				*(itrThis++) += *(itrP++);
			}
		}
		void operator -=(const poly_func p) {
			if (!p)
				return;
			else if (!(*this)) {
				this->Coeffs = p.Coeffs * (D_t)-1;
				return;
			}

			
			if (this->Degree() >= p.Degree())
			{
				Coeffs.resize(p.Coeffs.size());
				Coeffs.reserve(Coeffs.size());
			}
			else
				Coeffs.reserve(p.Coeffs.size());

			auto itrThis = Coeffs.begin(); //Coeffs iterator

			for (auto itrP = p.Coeffs.begin(); itrP != p.Coeffs.end();)//p iterator
				*(itrThis++) -= *(itrP++);
		}
		void operator *=(const poly_func p) {
			if (!p) {
				this->Coeffs = {};
				return;
			}
			else if (!*this) 
				return;

			*this = *this * p;
		}
		void operator /=(const poly_func p) {
			if (!*this)
				return;
			*this = *this / p;
		}
		void operator %=(const poly_func p) {
			if (!*this)
				return;

			*this -= (*this / p) * p;
		}
		void operator +=(const D_t constant) {
			if (Coeffs.empty())
				Coeffs.push_back(constant);
			else
				Coeffs.back() += constant;
		}
		void operator -=(const D_t constant) {
			if (Coeffs.empty())
				Coeffs.push_back(-constant);
			else
				Coeffs.back() -= constant;
		}
		void operator *=(const D_t scale) {
			for (auto& coefficient : Coeffs)
				coefficient *= scale;
		}
		void operator /=(const D_t scale) {
			for (auto& coefficient : Coeffs)
				coefficient /= scale;
		}
		void operator %=(const D_t divident) {
			for (auto& coefficient : Coeffs)
				coefficient = remainder<D_t, D_t>(coefficient, divident);
		}
		void derive() {
			if (Coeffs.empty()) 
				return;
		
			auto itrThis = Coeffs.rbegin();
			auto itrStop = Coeffs.rend() - 1;

			D_t scale = 1.0;
			while (itrThis != itrStop) {
				auto prevThis = itrThis++;
				prevThis = scale * *itrThis;
			}

			Coeffs.erase(itrStop);
		}
		void integrate(const D_t constant = 0) {
			if (Coeffs.empty()) {
				Coeffs.assign({ (D_t)1, constant });
				return;
			}

			Coeffs.push_back(constant);

			if (Coeffs.size() > 2u)
			{
				auto itrThis = Coeffs.rbegin() + 2;
				auto itrStop = Coeffs.rend();

				double scale = (D_t)2;
				while (itrThis != itrStop) {
					auto prevThis = itrThis++;
					prevThis = itrThis / (scale++);
				}
			}
		}
	};

	template<class D_t>
	poly_func<D_t> zero() { 
		return {};
	}

	template<class D_t>
	poly_func<D_t> one() {
		return { (D_t)1 };
	}

	template<class D_t>
	poly_func<D_t> signif(const size_t degree) {
		std::vector<D_t> dominant(degree + 1u);
		dominant.front() = (D_t)1;

		return poly_func<D_t>(dominant);
	}

	template<class D_t>
	template<class R_t>
	inline std::function<R_t(D_t)> poly_func<D_t>::Function() const
	{
		return [=](const D_t x)->R_t {
			return (*this)(x);
		};
	}

	template<class D_t>
	inline poly_func<D_t> poly_func<D_t>::AndRoots(const std::vector<D_t> roots) const
	{
		poly_func<D_t> p = *this;
		p.addRoots(roots);

		return p;
	}

	template<class D_t>
	inline poly_func<D_t> poly_func<D_t>::AndRoots(const std::initializer_list<D_t> roots) const
	{
		poly_func<D_t> p = *this;
		p.addRoots(roots);

		return p;
	}

	template<class D_t>
	inline void poly_func<D_t>::normalize()
	{
		if (!*this) return;

		while (Coeffs.front() == 0) { Coeffs.erase(Coeffs.begin()); }

		auto denominator = Coeffs.front();
		for (auto& coefficient : Coeffs)
			coefficient /= denominator;
	}

	template<class D_t>
	inline void poly_func<D_t>::simplify()
	{
		if (!*this) return;

		while (Coeffs.back() == 0) { Coeffs.pop_back(); }
		while (Coeffs.front() == 0) { Coeffs.erase(Coeffs.begin()); }
	}

	template<class D_t>
	inline void poly_func<D_t>::addRoots(const std::vector<D_t> roots)
	{
		for (auto& root : roots)
			*this *= poly_func({ (D_t)1, -root });
	}

	template<class D_t>
	inline void poly_func<D_t>::addRoots(const std::initializer_list<D_t> roots)
	{
		for (auto& root : roots)
			*this *= poly_func({ (D_t)1, -root });
	}

	template<class D_t>
	inline D_t poly_func<D_t>::operator()(const D_t x) const
	{
		D_t y = (D_t)0;
		size_t power = this->Degree();
		auto itrCoeff = Coeffs.begin();

		while (itrCoeff != Coeffs.end()) {
			if (*itrCoeff)
				y += *(itrCoeff++) * pow<D_t, size_t>(x, power--);
			else {
				++itrCoeff;
				--power;
			}
		}

		return y;
	}

	template<class D_t>
	inline poly_func<D_t>::operator std::function<D_t(D_t)>() const{
		return this->Function<D_t>();
	}

	template<class D_t>
	inline poly_func<D_t>::operator std::string() const
	{
		return String("x");
	}

	template<class D_t>
	template<class R_t>
	inline std::string poly_func<D_t>::String(const std::string variableName) const
	{
		if (*this) {
			size_t power = this->Degree();//power of current term
			std::string expression;//higher level form

			if (power != 0) {
				/*
				   generates [ variableName^pow] or [ variableName](for pow = 1) or [](for pow = 0) depending on pow
				*/
				auto var_power_string = [&](size_t pow, bool gap = true)->std::string {
					if (pow > 1)
						return (gap ? " " : "") + variableName + "^" + std::to_string(pow);
					if (pow) return (gap ? " " : "") + variableName;

					return "";
				};
				/*
				  generates [|coeff|] or [](for coeff = 0 or 1) depending on coeff
				*/
				auto coeff_string = [&](R_t coeff)->std::string{
					if (abs(coeff) != (R_t)1)
						return std::to_string(abs(coeff));
					return "";
				};
				/*
				  appends first term depending on it's coefficient
				*/
				{
					R_t first_coeff = (R_t)Coeffs.front();
					if (first_coeff)
						expression = (first_coeff < (R_t)0 ? "-" : "") + coeff_string(first_coeff) + var_power_string(power--, abs(first_coeff) != (R_t)1);//([-]|[])[|coefficient|]([ ]|[])[variableName^degree]
					else
						--power;
				}
				auto permCoeff = Coeffs.end() - 1;
				for (auto itrCoeff = Coeffs.begin() + 1; itrCoeff != permCoeff; ++itrCoeff)
				{
					R_t current_coeff = (R_t)*itrCoeff;
					if (current_coeff)
						expression += (current_coeff >= 0 ? " + " : " - ") + coeff_string(current_coeff) + var_power_string(power--, abs(current_coeff) != (R_t)1);//[ sign ][|coefficient|]([ ]|[])[variableName^degree]
					else
						--power;
				}
				if ((R_t)Coeffs.back() != 0)
					if (expression == "")
						expression = ((R_t)Coeffs.back() < 0 ? "-" : "") + std::to_string(abs((R_t)Coeffs.back()));
					else
						expression += ((R_t)Coeffs.back() >= 0 ? " + " : " - ") + std::to_string(abs((R_t)Coeffs.back()));


				return expression;
			}
			else
				return std::to_string((R_t)Coeffs.back());
		}
		return "0";
	}

	template<class D_t>
	inline poly_func<D_t> FromRoots(const std::vector<D_t> roots)
	{
		return one<D_t>().AndRoots(roots);
	}
	
	template<class D_t>
	inline poly_func<D_t> FromRoots(const std::initializer_list<D_t> roots)
	{
		return one<D_t>().AndRoots(roots);
	}
}

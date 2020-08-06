#include "FuncShorts.h"

func_t PolynomialFunction(std::vector<double> coeffs)
{
	return
		[&coeffs](double x) {
		double y = 0.0;

		uint8_t i = 0;
		for (auto& each : coeffs)
			y += pow(x, coeffs.size() - 1 - i++) * each;

		return y;
	};
}

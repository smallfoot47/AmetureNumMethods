#include "NumMethods.h"
#include "FuncShorts.h"

#define ENABLE_LOG//this segment is temporary. Will try to automate this via a function call; comment off to disable log
#include "log.h"

void NewtonsMethodApp()//Example on applying NewtonsMethod
{
	func_t Func = MAKE_FUNC(
		return cos(x) - x * exp(x);  //define function here: in terms of x
	);

	func_t Delta = MAKE_FUNC(
		return -sin(x) - (x + 1) * exp(x);  //define derivative here: in terms of x
	);

	double last_x = 1;  //initialize to x0
	std::cout << "x0>> "; PRINT_AT_VALUE(Func, 0.5);

	for (int i = 0; i < 20; ++i)
	{
		/*
		  auto BaseFunc = NewtonsMethod(Func, 0.5, i)//function, x0, iteration count: returns (x, y)
		  #Alternative for functions that are harder to differentiate, at the cost of the rate of convergence
		*/
		auto NewtonsFunc = NewtonsMethod(Func, Delta, last_x, i + 1);//function, derivative, x0, iteration count: returns (x, y)
		std::cout << "x" << i + 1 << ">> "; PRINT_FUNC_COORD(NewtonsFunc);
		LOG("Variation: " << abs(1 - NewtonsFunc.independent / last_x) * 100 << "%");

		last_x = NewtonsFunc.x;
	}
}

void BisectionMethodApp()//Example on applying BisectionMethod AKA IntermediateValueMethod
{
	func_t Func = MAKE_FUNC(
		return cos(x) - x * exp(x);
	);

	double last_x = 0.5;  //initialize to x0
	std::cout << "x0>> "; PRINT_AT_VALUE(Func, 0.5);

	for (int i = 0; i < 20; ++i)
	{
		auto BisectionFunc = BisectionMethod(Func, { 0, 1 }, i + 1);//function, interval, iteration count: returns (x, y); uses Intermediate Value Theorem
		std::cout << "x" << i + 1 << ">> "; PRINT_FUNC_COORD(BisectionFunc);
		LOG("Variation: " << abs(1 - BisectionFunc.independent / last_x) * 100 << "%");

		last_x = BisectionFunc.independent;
	}
}

void FunctionAbstractionApp()
{
	func_t Func = MAKE_FUNC(
		return cos(x) - x * exp(x);
	);
	
	Interval<double> domain = { 0, 1 };
	auto Mapping = Function_to_Map(Func, domain, 0.1);
	auto Approx = Interpolate(Mapping, domain);

	for (double x = 0; x <= 1; x += 0.25)
		std::cout << " at x = " << x << ", Approx(x) = " << Approx(x) << std::endl;
}

int main()
{
	LOG("LOG initailized");

	NewtonsMethodApp();
	//BisectionMethodApp();
	//FunctionAbstractionApp();
}

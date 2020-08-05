#include "NumMethods.h"
#include "FuncShorts.h"

#define ENABLE_LOG//this segment is temporary. Will try to automate this via a function call; comment off to disable log
#include "log.h"

int main()
{
	LOG("LOG is working");

	func_t Func = MAKE_FUNC(
		return cos(x) - x * exp(x);  //define function here: in terms of x
	);

	auto BaseFunc = IntermediateValue(Func, { 0, 1 }, 100);//function, interval, iteration count: returns (x, y)
	PRINT_FUNC_COORD(BaseFunc);
}
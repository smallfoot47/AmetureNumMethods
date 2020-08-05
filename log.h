#pragma once
#include <iostream>
//define ENABLE_LOG before inclusion to enable
#ifdef ENABLE_LOG
#define LOG(x) std::cout << "LOG::" << x << std::endl //logs statements concatinated by cascades
#define LOG_COMMAND(x) x //logs statement executions
#else
#define LOG(x) //#define ENABLE_LOG before #include "...log.h" to enable
#define LOG_COMMAND(x) //#define ENABLE_LOG before #include "...log.h" to enable
#endif
#pragma once
#include <iostream>
//define ENABLE_LOG before inclusion to enable
#ifdef ENABLE_LOG
 #define LOG(x) std::cout << "LINE:" << __LINE__ << ":: " << x << std::endl //logs statements concatinated by cascades
 #define FIXED_LOG(x) std::cout <<"LINE:" << __LINE__ << ":: " << x << " "; //logs statements concatinated by cascades and fixes the cursor at the end 
 #define LOG_COMMAND(x) x //logs statement executions
#else
 #define LOG(x) //#define ENABLE_LOG before #include "...log.h" to enable
 #define FIXED_LOG(x) //#define ENABLE_LOG before #include "...log.h" to enable
 #define LOG_COMMAND(x) //#define ENABLE_LOG before #include "...log.h" to enable
#endif

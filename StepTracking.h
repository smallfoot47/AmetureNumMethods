#pragma once

#ifdef ENABLE_STEP_TRACKING
  #define PUSH_STEP_LOG(output_command) { output_command }
#else
 #define PUSH_STEP_LOG(ignored_command)
#endif // ENABLE_STEP_TRACKING
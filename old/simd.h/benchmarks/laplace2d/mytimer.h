#include <sys/time.h>

#ifndef __TIMER__H__
#define __TIMER__H__

inline double get_wtime(){
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#endif // __TIMER__H__

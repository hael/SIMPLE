/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 1st May 2011
 *                           Modified: 9th of October 2014
 *      October 2014
 *
 *      c code for the MACOSX code this wrapps around the timming_c.cpp 
 *      methods for the timming methods. It also includes the precise nano
 *      second testing for both LINUX and MACOSX where the
 *      clock_gettime(CLOCK_REALTIME, &time) is not defined in the time.h 
 *      Header file.
 *
 * @precisions normal z -> s d c
 */

/* The Simple header */
#include "simple.h"

double get_cpu_time_c_(tt)
float tt[2];
{
  double cpu_time = 0.0;
  double *elapsed_cpu_time = NULL;
  double *tt_c = (double*)tt;

  elapsed_cpu_time = (double*)malloc(sizeof(double));

  GETCPUTIME_C(tt_c,elapsed_cpu_time);
  cpu_time = *elapsed_cpu_time;

  free(elapsed_cpu_time);

  return(cpu_time);

}
void GET_CPU_TIME_C_() __attribute__((weak,alias("get_cpu_time_c_")));
void get_cpu_time_c__() __attribute__((weak,alias("get_cpu_time_c_")));
void GET_CPU_TIME_C__() __attribute__((weak,alias("get_cpu_time_c_")));

#if defined (BENCH)
unsigned long long int get_cpu_nano_time_c_() {
  long long int t = 0;

  t = GETCPUNANOTIME_C(t);

  return(t);
}

void GET_CPU_NANO_TIME_C_() __attribute__((weak,alias("get_cpu_nano_time_c_")));
void get_cpu_nano_time_c__() __attribute__((weak,alias("get_cpu_nano_time_c_")));
void GET_CPU_NANO_TIME_C__() __attribute__((weak,alias("get_cpu_nano_time_c_")));
#endif

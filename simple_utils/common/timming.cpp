/*
 *   -- MAGMA addon
 *      Author: Frederic Bonnet, Date: 5th of September 2011
 *                           Modified: 10th of October 2014
 *                           Modified: 19th of May 2015
 *      May 2015
 *
 *   code to perform timing tests.
 *
 * @precisions normal z -> s d c
 */

#include "timming.h"
/*Header for the MACOSX clock environment */
#if defined (MACOSX)
#include <mach/clock.h>
#include <mach/mach.h>
#endif
#define M_LN10 2.30258509299404568402  /* log_e 10 */
/* 
 *
 * getcpunanotime_c will get the cpu time in nano seconds.
 *
 */

#ifdef __cplusplus
extern "C" {
#endif
/******************************************************************************/
/**
 *  FORTRAN API - timming functions (simple interface)
 **/

/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions / Data on CPU to get the time of instance
    void gettimeofday_c(struct timeval *restrict tp, void *restrict tzp);
*/
  void GETTIMEOFDAY_C( double Time[2] ) {
    struct timeval tp;
    struct timezone tz;

    gettimeofday(&tp, &tz);

    Time[0] = tp.tv_sec;
    Time[1] = tp.tv_usec;
  }
/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions / Data on CPU calculates the time diffrence between
    time instances
    void elapsed_time_c( double start_Time[2], double end_Time[2],
                         double *elapsed_time)
*/ 
  void ELAPSED_TIME_C( double start_Time[2], double end_Time[2], double *elapsed_time) {
    *elapsed_time = (end_Time[0] - start_Time[0]) + (end_Time[1] - start_Time[1])/1000000.0;
  }
/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions / Data on CPU gets the cpu time, true time
    void GETCPUTIME_C( double cpu_time[2], double *elapsed_cpu_time)
*/
  void GETCPUTIME_C( double cpu_time[2], double *elapsed_cpu_time) {
    int process;
    struct rusage cpu_usage;
    process = 0;
    getrusage(process, &cpu_usage);
    cpu_time[0] = cpu_usage.ru_utime.tv_sec+((cpu_usage.ru_utime.tv_usec)/1000000.);
    cpu_time[1] = cpu_usage.ru_stime.tv_sec+((cpu_usage.ru_stime.tv_usec)/1000000.);
    *elapsed_cpu_time = cpu_time[0]+cpu_time[1];
  }
/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions / Data on CPU gets the cpu time, true time nano secs
    unsigned long long int GETCPUNANOTIME_C(long long int t)
*/
  unsigned long long int GETCPUNANOTIME_C(long long int t) {
    struct timespec time;

#if defined (LINUX)
    clock_gettime(CLOCK_REALTIME, &time);
    t  = time.tv_sec;
    t *= 1000000000;
    t += time.tv_nsec;
#elif defined (MACOSX)
    clock_serv_t cclock;
    mach_timespec_t mtime;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mtime);
    mach_port_deallocate(mach_host_self(), cclock);
    t  = mtime.tv_sec;
    t *= 1000000000;
    t += mtime.tv_nsec;
#endif
    return(t);
  }
/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions that gets the length in terms of digits a positive
 integer, float and double
*/
  int get_length_of_string(int *int_in){
    int length = 1;
    int num = *int_in;
    if (num > 0) length = floor(log10(abs(num))) + 1;
    return length;
  }
  int get_length_of_stringf(float *float_in){
    int length = 1;
    float num = *float_in;
    if (num > 0) length = floor(log10(abs(num))) + 1;
    return length;
  }
  int get_length_of_stringd(double *double_in){
    int length = 1;
    double num = *double_in;
    if (num > 0) length = floor(log10(abs(num))) + 1;
    return length;
  }
/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions that copnverts a integer to a string when input
 int is greater or equal to 0 with the input already setup with input integer
*/
  int convert_int2char_pos(char *char_out, int *int_in) {
    int rc = RC_SUCCESS; //return code
    int length = get_length_of_string(int_in);
    int size_int_in = (length) * sizeof(char);
    char *buff = (char*)malloc(size_int_in);
    sprintf(buff,"%i",*int_in);
    //copying over the buffer into the char_out
    for (int il=0 ; il < length ; il++) {
      char_out[il] = buff[il];
    }
    free(buff);
    return rc;
  }
/* ////////////////////////////////////////////////////////////////////////////
 -- C function definitions that converts a real to a string when input
 real is greater or equal to 0 with the input already setup with input integer
*/

  int convert_float2char_pos(char *char_out, float *float_in) {
    int rc = RC_SUCCESS; //return code
    int length = get_length_of_stringf(float_in);
    int size_float_in = (length) * sizeof(char);
    char *buff = (char*)malloc(size_float_in);
    sprintf(buff,"%f",*float_in);
    //copying over the buffer into the char_out
    for (int il=0 ; il < length ; il++) {
      char_out[il] = buff[il];
    }
    free(buff);
    return rc;
  }

/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions that converts a real to a string when input
 real is greater or equal to 0 with the input already setup with input integer
*/
  int convert_double2char_pos(char *char_out, double *double_in) {
    int rc = RC_SUCCESS; //return code
    int length = get_length_of_stringd(double_in);
    int size_double_in = (length) * sizeof(char);
    char *buff = (char*)malloc(size_double_in);
    sprintf(buff,"%d",*double_in);
    //copying over the buffer into the char_out
    for (int il=0 ; il < length ; il++) {
      char_out[il] = buff[il];
    }
    free(buff);
    return rc;
  }

/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions that copnverts a integer to a string when input
 int is greater or equal to 0 with the input not setup with input integer
*/
  int convert_int2char_indexed(char *char_out, int *iter, int *iint,
                               int *max_iter) {
    int rc = RC_SUCCESS; //return code
    int length = get_length_of_string(max_iter);
    int indexed_int = get_indexed_int(length);
    int num = indexed_int + *iint + *iter;
    int length_num = get_length_of_string(&num);
    int size_iter = (indexed_int) * sizeof(char);
    char *buff = (char*)malloc(size_iter);
    sprintf(buff,"%i",num);
    //copying over the buffer into the char_out
    for (int il=0 ; il < length_num ; il++) {
      char_out[il] = buff[il];
    }
    free(buff);
    return rc;
  }
/* /////////////////////////////////////////////////////////////////////////////
 -- C function definitions that increases the indexing by a power based 10 the
 index of the number
*/
  int get_indexed_int(int length) {
#if defined (LINUX)
    int indexed_int = exp10(length+1);
#else if defined (MACOSX)
    int indexed_int = exp(M_LN10 * length +1);
#endif
    return indexed_int;
  }
  
  /* the aliases for external access */
  extern "C" int convert_int2char_pos_() __attribute__((weak,alias("convert_int2char_pos")));
  extern "C" int convert_float2char_pos_() __attribute__((weak,alias("convert_float2char_pos")));
  extern "C" int convert_double2char_pos_() __attribute__((weak,alias("convert_double2char_pos")));
  extern "C" int convert_int2char_indexed_() __attribute__((weak,alias("convert_int2char_indexed")));
  extern "C" int get_length_of_string_() __attribute__((weak,alias("get_length_of_string")));
  extern "C" int get_length_of_stringf_() __attribute__((weak,alias("get_length_of_stringf")));
  extern "C" int get_length_of_stringd_() __attribute__((weak,alias("get_length_of_stringd")));
  extern "C" int get_indexed_int_() __attribute__((weak,alias("get_indexed_int")));
  extern "C" void gettimeofday_c_() __attribute__((weak,alias("GETTIMEOFDAY_C")));
  extern "C" void elapsed_time_c_() __attribute__((weak,alias("ELAPSED_TIME_C")));
  extern "C" void getcputime_c_() __attribute__((weak,alias("GETCPUTIME_C")));
  extern "C" void getcpunanotime_c_() __attribute__((weak,alias("GETCPUNANOTIME_C")));

#ifdef __cplusplus
}
#endif

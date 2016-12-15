/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 5th of September 2011
 *                           Modified: 10th of October 2014
 *                           Modified: 19th of May 2015
 *      March 2016
 *
 *   code to perform timing tests.
 *
 * @precisions normal z -> s d c
 */

#include "timming.h"

/* 
 *
 * getcpunanotime_c will get the cpu time in nano seconds.
 *
 */

/***************************************************************************/
/**
 *  FORTRAN API - timming functions (simple interface) for MacOSX
 **/


/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions that gets the length in terms of digits a positive
 integer, float and double
*/
/*int*/
int get_length_of_string_c_(int *int_in){
  int length = 1;
  length = get_length_of_string(int_in);
  return length;
}
void GET_LENGTH_OF_STRING_C_() __attribute__((weak,alias("get_length_of_string_c_")));
void get_length_of_string_c__() __attribute__((weak,alias("get_length_of_string_c_")));
void GET_LENGTH_OF_STRING_C__() __attribute__((weak,alias("get_length_of_string_c_")));
/*float*/
int get_length_of_stringf_c_(float *float_in){
  int length = 1;
  length = get_length_of_stringf(float_in);
  return length;
}
void GET_LENGTH_OF_STRINGF_C_() __attribute__((weak,alias("get_length_of_stringf_c_")));
void get_length_of_stringf_c__() __attribute__((weak,alias("get_length_of_stringf_c_")));
void GET_LENGTH_OF_STRINGF_C__() __attribute__((weak,alias("get_length_of_stringf_c_")));
/*double*/
int get_length_of_stringd_c_(double *double_in){
  int length = 1;
  length = get_length_of_stringd(double_in);
  return length;
}
void GET_LENGTH_OF_STRINGD_C_() __attribute__((weak,alias("get_length_of_stringd_c_")));
void get_length_of_stringd_c__() __attribute__((weak,alias("get_length_of_stringd_c_")));
void GET_LENGTH_OF_STRINGD_C__() __attribute__((weak,alias("get_length_of_stringd_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions that gets the length in terms of digits a positive
 integer
*/
int get_indexed_int_c_(int *length){
  int length_in = *length;
  int indexed_int;
  indexed_int = get_indexed_int(length_in);
  return indexed_int;
}
void GET_INDEXED_INT_C_() __attribute__((weak,alias("get_indexed_int_c_")));
void get_indexed_int_c__() __attribute__((weak,alias("get_indexed_int_c_")));
void GET_INDEXED_INT_C__() __attribute__((weak,alias("get_indexed_int_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions that copnverts a integer to a string when input
 int,float,double is greater or equal to 0
*/
/*int*/
int convert_int2char_pos_c_(char *char_out, int *int_in) {
  int rc = RC_SUCCESS; //return code
  rc = convert_int2char_pos(char_out, int_in);
  return rc;
}
void CONVERT_INT2CHAR_POS_C_() __attribute__((weak,alias("convert_int2char_pos_c_")));
void convert_int2char_pos_c__() __attribute__((weak,alias("convert_int2char_pos_c_")));
void CONVERT_INT2CHAR_POS_C__() __attribute__((weak,alias("convert_int2char_pos_c_")));
/*float*/
int convert_float2char_pos_c_(char *char_out, float *float_in) {
  int rc = RC_SUCCESS; //return code
  rc = convert_float2char_pos(char_out, float_in);
  return rc;
}
void CONVERT_FLOAT2CHAR_POS_C_() __attribute__((weak,alias("convert_float2char_pos_c_")));
void convert_float2char_pos_c__() __attribute__((weak,alias("convert_float2char_pos_c_")));
void CONVERT_FLOAT2CHAR_POS_C__() __attribute__((weak,alias("convert_float2char_pos_c_")));
/*double*/
int convert_double2char_pos_c_(char *char_out, double *double_in) {
  int rc = RC_SUCCESS; //return code
  rc = convert_double2char_pos(char_out, double_in);
  return rc;
}
void CONVERT_DOUBLE2CHAR_POS_C_() __attribute__((weak,alias("convert_double2char_pos_c_")));
void convert_double2char_pos_c__() __attribute__((weak,alias("convert_double2char_pos_c_")));
void CONVERT_DOUBLE2CHAR_POS_C__() __attribute__((weak,alias("convert_double2char_pos_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions that copnverts a integer to a string when input
 int is greater or equal to 0
*/
int convert_int2char_indexed_c_(char *char_out, int *iter, int *iint,
                                int *max_iter) {
  int rc = RC_SUCCESS; //return code
  rc = convert_int2char_indexed(char_out, iter, iint, max_iter);
  return rc;
}
void CONVERT_INT2CHAR_INDEXED_C_() __attribute__((weak,alias("convert_int2char_indexed_c_")));
void convert_int2char_indexed_c__() __attribute__((weak,alias("convert_int2char_indexed_c_")));
void CONVERT_INT2CHAR_INDEXED_C__() __attribute__((weak,alias("convert_int2char_indexed_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions / Data on CPU to get the time of instance
    void gettimeofday_c_(struct timeval *restrict tp, void *restrict tzp);
*/
void gettimeofday_c_( double Time[2] ) {
  GETTIMEOFDAY_C(Time);
}
void GETTIMEOFDAY_C_() __attribute__((weak,alias("gettimeofday_c_")));
void gettimeofday_c__() __attribute__((weak,alias("gettimeofday_c_")));
void GETTIMEOFDAY_C__() __attribute__((weak,alias("gettimeofday_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions / Data on CPU calculates the time diffrence between
    time instances
    void elapsed_time_c_( double start_Time[2], double end_Time[2],
                         double *elapsed_time)
*/ 
 void elapsed_time_c_( double start_Time[2], double end_Time[2], double *elapsed_time) {
   ELAPSED_TIME_C( start_Time, end_Time, elapsed_time);
 }
void ELAPSED_TIME_C_() __attribute__((weak,alias("elapsed_time_c_")));
void elapsed_time_c__() __attribute__((weak,alias("elapsed_time_c_")));
void ELAPSED_TIME_C__() __attribute__((weak,alias("elapsed_time_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions / Data on CPU gets the cpu time, true time
    void GETCPUTIME_C( double cpu_time[2], double *elapsed_cpu_time)
*/
void getcputime_c_( double cpu_time[2], double *elapsed_cpu_time) {
  GETCPUTIME_C(cpu_time, elapsed_cpu_time);
}
void GETCPUTIME_C_() __attribute__((weak,alias("getcputime_c_")));
void getcputime_c__() __attribute__((weak,alias("getcputime_c_")));
void GETCPUTIME_C__() __attribute__((weak,alias("getcputime_c_")));
/* //////////////////////////////////////////////////////////////////////////// 
 -- C function definitions / Data on CPU gets the cpu time, true time nano secs
    unsigned long long int GETCPUNANOTIME_C(long long int t)
*/
unsigned long long int getcpunanotime_c_(long long int t) {
  GETCPUNANOTIME_C(t);
  return(t);
}
void GETCPUNANOTIME_C_() __attribute__((weak,alias("getcpunanotime_c_")));
void getcpunanotime_c__() __attribute__((weak,alias("getcpunanotime_c_")));
void GETCPUNANOTIME_C__() __attribute__((weak,alias("getcpunanotime_c_")));


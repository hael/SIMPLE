\
/*
  C/FORTRAN  preprocessor macros for timing module blocks
  Note: __PRETTY_FUNCTION__ is GNU specific and not standard
        __FILE__ is standard but is different some compilers
        __LINE__ is defined in the standard


     2017-04-08        Michael Eager (michael.eager@monash.edu)
*/

#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H

/*  getting into ## preprocessor magic 
#define CAT(prefix, suffix)            prefix ## suffix
#define _UNIQUE_LABEL(prefix, suffix)  CAT(prefix, suffix)
#define UNIQUE_LABEL(prefix)           _UNIQUE_LABEL(prefix, __LINE__)
*/
/* calculate the number of arguments in macro - max 5 */
/*
#define VA_ARGS_NUM_PRIV(P1, P2, P3, P4, P5, P6, Pn, ...)  Pn
#define VA_ARGS_NUM(...) VA_ARGS_NUM_PRIV(-1, ##__VA_ARGS__, 5, 4, 3, 2, 1, 0)
*/



/*
#define VA_NARGS_IMPL(_1, _2, _3, _4, _5, N, ...) N
#define VA_NARGS(...) VA_NARGS_IMPL(__VA_ARGS__, 5, 4, 3, 2, 1)
#define VA_NUM_ARGS(...)                            \
  (sizeof(#__VA_ARGS__) == sizeof("")               \
   ? 0 : VA_NUM_ARGS_IMPL(__VA_ARGS__, 5,4,3,2,1))
*/
/*

#ifdef _MSC_VER // Microsoft compilers
#define GET_ARG_COUNT(...)  INTERNAL_EXPAND_ARGS_PRIVATE(INTERNAL_ARGS_AUGMENTER(__VA_ARGS__))
#define INTERNAL_ARGS_AUGMENTER(...) unused, __VA_ARGS__
#define INTERNAL_EXPAND(x) x
#define INTERNAL_EXPAND_ARGS_PRIVATE(...) INTERNAL_EXPAND(INTERNAL_GET_ARG_COUNT_PRIVATE(__VA_ARGS__, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0))
#define INTERNAL_GET_ARG_COUNT_PRIVATE(_1_, _2_, _3_, _4_, _5_, _6_, _7_, _8_, _9_, _10_, _11_, _12_, _13_, _14_, _15_, _16_, _17_, _18_, _19_, _20_, _21_, _22_, _23_, _24_, _25_, _26_, _27_, _28_, _29_, _30_, _31_, _32_, _33_, _34_, _35_, _36, _37, _38, _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62, _63, _64, _65, _66, _67, _68, _69, _70, count, ...) count

#else // Non-Microsoft compilers
#if defined(__GFORTRAN__)
#define GET_ARG_COUNT(...) INTERNAL_GET_ARG_COUNT_PRIVATE(0, ##__VA_ARGS__, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define INTERNAL_GET_ARG_COUNT_PRIVATE(_0, _1_, _2_, _3_, _4_, _5_, _6_, _7_, _8_, _9_, _10_, _11_, _12_, _13_, _14_, _15_, _16_, _17_, _18_, _19_, _20_, _21_, _22_, _23_, _24_, _25_, _26_, _27_, _28_, _29_, _30_, _31_, _32_, _33_, _34_, _35_, _36, _37, _38, _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62, _63, _64, _65, _66, _67, _68, _69, _70, count, ...) count





#elif defined(__INTEL_COMPILER)

#else
    #error "Name-mangling macro not set for your compiler."
#endif


#endif
*/

#define TBLOCK()                                \
  print *,"TBLOCK:  Start timer: ", tic()


#define TSTOP()                                                         \
  write(*,'(A,A,1i4,A,F20.10)') __FILE__,":",__LINE__,":Elapsed time (sec) ", toc()

#define TBREAK(TSTRING)                         \
  write(*,'(A,A,1i4,A,A,F20.10)') __FILE__,":",__LINE__,TSTRING," time (sec)" toc()


#define TPROFILE5(NLOOP,TOK1,TOK2,TOK3,TOK4,TOK5) \
  call timer_profile_setup(NUMLOOPS,TOK1,TOK2,TOK3,TOK4,TOK5)
#define TPROFILE4(NLOOP,TOK1,TOK2,TOK3,TOK4)      \
  call timer_profile_setup(NUMLOOPS,TOK1,TOK2,TOK3,TOK4)
#define TPROFILE3(NLOOP,TOK1,TOK2,TOK3)           \
  call timer_profile_setup(NUMLOOPS,TOK1,TOK2,TOK3)
#define TPROFILE2(NLOOP,TOK1,TOK2)                \
  call timer_profile_setup(NUMLOOPS,TOK1,TOK2)
#define TPROFILE1(NLOOP,TOK1)                     \
  call timer_profile_setup(NUMLOOPS,TOK1)

/*#define TPROFILE(NUMLOOPS, ...)  TPROFILE_AUX##VA_NUM_ARGS(__VA_ARGS__)(NUMLOOPS,__VA_ARGS__)
 */
#define TBEG(TOKEN) call timer_profile_start(TOKEN)

#define TEND(TOKEN) call timer_profile_break(TOKEN)

#define TREPORT(COMMENT) call timer_profile_report(COMMENT)


#define START_TIMER_LOOP(N)                     \
  call timer_loop_start(N);                     \
  do

#define STOP_TIMER_LOOP()                       \
  if (.not.in_timer_loop()) exit;               \
  end do;                                       \
  call timer_loop_end();



#if defined(OPENMP)
#define TBLOCKOMP()                             \
  print *,"Start timer: ", ticOMP()

#define TSTOPOMP()                                                      \
  write(*,'(A,A,1i4,A,1d20.10)') __FILE__,":",__LINE__,": Elapsed time (s) ", tocOMP()

#endif /*OPENMP*/


/*

write(*,'(A,1d20.10)') "__PRETTY_FUNCTION__:__FILE__:__LINE__: Elapsed time (sec): "  REAL(etime,dp)

#define TIME_LOOP(...) TIMER_LOOP_START(integer(dp)) (__VA_ARGS__)

#define TIMED_LOOP(XX) TIMER_LOOP_START(integer(dp)) \
  { XX } \
  TIMER_LOOP_END()

#define TIMER_LOOP_START(typ)                                    \
  typ :: UNIQUE_LABEL(timestamp);UNIQUE_LABEL(timestamp) = tic(); \
  real(dp),dimension(3)::UNIQUE_LABEL(elapsed); \
          integer::loopindex=1; do loopindex=1,3


#define TIMER_LOOP_END()                                           \
  UNIQUE_LABEL(elapsed)(loopindex) = toc(UNIQUE_LABEL(timestamp)); \
  end do;\
write(*,'(A,A,1i4,A,1d20.10,A,1d20.10,A)') __FILE__,":",__LINE__,":Elapsed time(s): average ",  SUM(UNIQUE_LABEL(elapsed))/REAL(3.,dp) , " (sec)   Min ", MINVAL(UNIQUE_LABEL(elapsed)), " (sec)"

  */




/*
#define TBLOCK( EVALBLOCK )                \
use simple_timer;                            \
integer, precision :: __t1;                  \
real,precision :: elapsed;                   \
starttime = tic();                           \
EVALBLOCK  ;                                 \
write(*,"('__PRETTY_FUNCTION__:__FILE__:__LINE__: \
Elapsed time (sec): ' , 1d20.10,/)" ) toc(__t1)



#define TBLOCK2( EVALBLOCK , SUBSTR )        \
use simple_timer;                            \
integer, precision :: __t1;                  \
real,precision :: elapsed;                   \
starttime = tic();                           \
EVALBLOCK  ;                                 \
print *, "SUBSTR:__LINE__: Elapsed time (sec): ", toc(__t1)

*/

/* now your watch has ended */


#ifdef CUDA
#define GPU_STOPWATCH( EVALBLOCK )              \
use simple_timer_cuda;                          \
type(cudaEvent) :: __t1;                        \
real(fp_kind) :: elapsed;                       \
starttime = tic();                              \
EVALBLOCK  ;                                    \
  print *, "__FILE__:__LINE__: CUDA Event timer (sec): ", toc(__t1)

#endif

#define CPU_STOPWATCH( EVALBLOCK )              \
double, precision :: __t1,__t2;                 \
call cpu_time(__t);                             \
EVALBLOCK  ;                                    \
call cpu_time(__t2);                            \
print *, "__FILE__:__LINE__: CPU wall time (sec): ", __t2-__t1



#endif/*SIMPLE_TIMER_H*/

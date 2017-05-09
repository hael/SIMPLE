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

/*  getting into ## preprocessor magic */
#define CAT(prefix, suffix)            prefix ## suffix
#define _UNIQUE_LABEL(prefix, suffix)  CAT(prefix, suffix)
#define UNIQUE_LABEL(prefix)           _UNIQUE_LABEL(prefix, __LINE__)
/* calculate the number of arguments in macro - max 5 */
#define VA_ARGS_NUM_PRIV(P1, P2, P3, P4, P5, P6, Pn, ...) Pn
#define VA_ARGS_NUM(...) VA_ARGS_NUM_PRIV(-1, ##__VA_ARGS__, 5, 4, 3, 2, 1, 0)


#define TBLOCK()                                \
  print *,"TBLOCK:  Start timer: ", tic()


#define TSTOP()                                                         \
  write(*,'(A,A,1i4,A,F20.10)') __FILE__,":",__LINE__,":Elapsed time (sec) ", toc()

#define TBREAK(TSTRING)                         \
  write(*,'(A,A,1i4,A,A,F20.10)') __FILE__,":",__LINE__,TSTRING," time (sec)" toc()


#define TPROFILE(NUMLOOPS, ... )  call timer_profile_setup(NUMLOOPS, VA_ARGS_NUM(__VA_ARGS__), __VA_ARGS__)
#define TBEG(LABEL) call timer_profile_start(LABEL)
#define TEND(LABEL) call timer_profile_break(LABEL)
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


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

/* calculate the number of arguments in macro - max 10 */

!    ifdef PROFILER
/* C99 __VA_ARGS__ versions */ /* If only ## worked.*/
#define c99_count(...)    _c99_count1 ( , ##__VA_ARGS__)/* */
#define _c99_count1(...)  _c99_count2 (__VA_ARGS__,10,9,8,7,6,5,4,3,2,1,0)
#define _c99_count2(_,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,n,...) n



#define __VA_NARG__(...) \
(__VA_NARG_(_0, ## __VA_ARGS__, __RSEQ_N()))
#define __VA_NARG_(...) \
  __VA_ARG_N(__VA_ARGS__)
#define __VA_ARG_N( \
                   _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
                   _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
                   _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
                   _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
                   _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
                   _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
                   _61,_62,_63,N,...) N
#define __RSEQ_N() \
  62, 61, 60,                             \
    59, 58, 57, 56, 55, 54, 53, 52, 51, 50, \
    49, 48, 47, 46, 45, 44, 43, 42, 41, 40, \
    39, 38, 37, 36, 35, 34, 33, 32, 31, 30, \
    29, 28, 27, 26, 25, 24, 23, 22, 21, 20, \
    19, 18, 17, 16, 15, 14, 13, 12, 11, 10, \
    9,  8,  7,  6,  5,  4,  3,  2,  1,  0 




#define TPROFILER(NLOOPS,IDX,...)  block;\
  use simple_timer; \
  character(len=20)::p_tmp;        \
  character(len=80)::p_comment;    \
  integer(dp) :: IDX,tn;        \
  integer,parameter :: nv=c99_count (__VA_ARGS__);  \
  character(len=20),dimension(nv)::p_tokens=(/ __VA_ARGS__ /);  \
  tn=tic();\
  call timer_profile_setup(NLOOPS,nv,p_tokens);


#define TBEG(TOKEN) p_tmp = TOKEN; \
 call timer_profile_start(trim(p_tmp))

#define TEND(TOKEN) p_tmp = TOKEN; \
 call timer_profile_break(trim(p_tmp))

#define TREPORT(COMMENT) p_comment = COMMENT;         \
  call timer_profile_report(trim(adjustl(p_comment)),toc(tn));  \
 end block

! endif   /* PROFILER */

#define TBLOCK() \
  print *,"TBLOCK:  Start timer: ", tic()

#define TIMER_BLOCK(YBLOCK,COMMENT) \
block; \
 use simple_timer;\
character(len=80) :: comment_,srcname; \
integer(dp) :: t1,srcline; \
t1=tic();srcline=__LINE__; \
 comment_=trim(adjustl(COMMENT)); \
srcname=trim(__FILE__); \
YBLOCK; \
 write(*,'(A,A,A,1i4,A,A,A,1ES20.6)') "TIMER_BLOCK:",trim(srcname),":",srcline,":",trim(adjustl(comment_)),': Elapsed time (sec): ', toc(t1); \
end block


#define TSTOP()                                 \
  write(*,'(A,A,1i4,A,F20.10)') __FILE__,":",__LINE__,": Elapsed time (sec) ", toc()

#define TBREAK(TSTRING)                         \
  write(*,'(A,A,1i4,A,A,F20.10)') __FILE__,":",__LINE__,TSTRING," time (sec)" toc()


#define START_TIMER_LOOP(NLOOP)  \
  block;character(len=80) :: cblock;\
  call timer_loop_start(NLOOP);     \
do


#define STOP_TIMER_LOOP_(COMMENT) \
  if (.not.in_timer_loop()) exit; \
  end do; \
  cblock=trim(COMMENT); \
  call timer_loop_end(trim(cblock));end block


/*   write(*,'(A,A,1i4)') __FILE__,":",__LINE__;   */

#define TIMER_BLOCK_LOOP(NLOOP,COMMENT,YBLOCK) block;      \
  character(len=*) :: cblock=trim(COMMENT);                \
  call timer_loop_start(NLOOP);do i=1,NLOOP; YBLOCK;\
  if (.not.in_timer_loop()) exit; \
end do;\
call timer_loop_end(cblock);end block

/*  write(*,'(A,A,1i4)') __FILE__,":",__LINE__; */



#ifdef OPENMP

#define TBLOCKOMP()                             \
  print *,"Start timer: ", ticOMP()

#define TSTOPOMP()                              \
  write(*,'(A,A,1i4,A,1d20.10)') __FILE__,":",__LINE__,": Elapsed time (s) ", tocOMP()

#endif /*OPENMP*/



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



#endif  /* SIMPLE_TIMER_H */

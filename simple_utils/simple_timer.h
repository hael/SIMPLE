
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


#define TPROFILER(NUMLOOPS,IDX,...)  block;\
  character(len=20)::p_tmp;        \
  character(len=80)::p_comment;    \
  integer(dp) :: IDX,tn,np;        \
  integer,parameter :: nv=c99_count (__VA_ARGS__);  \
  character(len=20)::p_tokens(nv);                  \
  np=NUMLOOPS;tn=tic();p_tokens=(/ __VA_ARGS__ /);  \
  call timer_profile_setup(np,nv,p_tokens);


#define TBEG(TOKEN) p_tmp = TOKEN; \
 call timer_profile_start(trim(p_tmp))

#define TEND(TOKEN) p_tmp = TOKEN; \
 call timer_profile_break(trim(p_tmp))

#define TREPORT(COMMENT) p_comment = COMMENT;         \
 call timer_profile_report(trim(p_comment),toc(tn));  \
 end block

! endif   /* PROFILER */

#define TBLOCK() \
  print *,"TBLOCK:  Start timer: ", tic()

#define TIMER_BLOCK(YBLOCK,COMMENT) \
block; \
character(len=80) :: cblock,srcname; \
integer(dp) :: t1,srcline; \
t1=tic();srcline=__LINE__; \
cblock=trim(COMMENT); \
srcname=trim(__FILE__); \
YBLOCK; \
 write(*,'(A,A,A,1i4,A ,A)') "TIMER_BLOCK:",trim(srcname),":",srcline,":",trim(cblock); \
 write(*,'(A,1F20.10)') 'Elapsed time (sec): ', toc(t1); \
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

/*-----------------------------------------------------------------------------*/
/* SIMPLE v3.0          Elmlund & Elmlund Lab         simplecryoem.com         */
/*-----------------------------------------------------------------------------*/
/** Simple timer macros
* \brief C/FORTRAN preprocessor macros for timing module blocks
* simple timer macros with variadic arg support
* Note: __FILENAME__ is used instead of __FILE__.
* __FILE__ is standard but is different some compilers
* __LINE__ is defined in the standard
* Intel preprocessor fpp does not support variadic macros
* PGI does not like block command statements
* See simple_test_timer and the basic and profile test modules
*
* \author
*       Michael Eager
*
* Revisions:
* Original TBLOCK TSTOP  TIMER_LOOP
* Profiling required variadic macros and C preprocessor
**/

#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H

/*  getting into ## preprocessor magic */
#define CAT(prefix, suffix)            prefix ## suffix
#define _UNIQUE_LABEL(prefix, suffix)  CAT(prefix, suffix)
#define UNIQUE_LABEL(prefix)           _UNIQUE_LABEL(prefix, __LINE__)


#if 0
/*calculate the number of arguments in macro - max 10*/
/*C99 __VA_ARGS__ versions */
#define c99_count(...)    _c99_count1 ( , ##__VA_ARGS__)
#define _c99_count1(...)  _c99_count2 (__VA_ARGS__,10,9,8,7,6,5,4,3,2,1,0)
#define _c99_count2(_,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,n,...) n


#define TPROFILER(NLOOPS,IDX,...)  block;\
  use simple_timer;\
  character(len=20)::p_tmp;\
  character(len=255)::p_comment;\
  integer :: IDX,np;\
  integer(timer_int_kind):: tn;\
  integer,parameter :: nv=c99_count (__VA_ARGS__);\
  character(255)::p_tokens= #__VA_ARGS__ ; \
  tn=tic();\
  np=NLOOPS;\
  call timer_profile_setup(np,nv,p_tokens);

#else

#define TPROFILER(NLOOPS,IDX,TOKENS) block;\
 use simple_timer;\
 character(len=20)::p_tmp;\
 character(len=255)::p_comment;\
 integer :: IDX,np;\
 integer(timer_int_kind):: tn;\
 character(len=255)::p_tokens= #TOKENS; \
 tn=tic();\
 np=NLOOPS;\
 call timer_profile_setup(np,-1,p_tokens);

#endif

#define TBEG(TOKEN) p_tmp = #TOKEN; \
 call timer_profile_start(trim(p_tmp))

#define TEND(TOKEN) p_tmp = #TOKEN; \
 call timer_profile_break(trim(p_tmp))

#define TREPORT(COMMENT) p_comment = #COMMENT; \
  call timer_profile_report(trim(adjustl(p_comment)),toc(tn));  \
 end block


#define TBLOCK() \
  print *,"TBLOCK:  Start timer: ", tic()

#define TIMER_BLOCK(YBLOCK,COMMENT) \
block; \
 use simple_timer;\
character(len=80) :: comment_,srcname; \
 integer(timer_int_kind) :: t1;integer :: srcline;  \
t1=tic();srcline=__LINE__; \
 comment_=trim(adjustl(COMMENT)); \
srcname=trim(__FILENAME__); \
YBLOCK; \
 write(*,'(A,A,A,1i4,A,A,A,1ES20.6)') "TIMER_BLOCK:",trim(srcname),":",srcline,":",trim(adjustl(comment_)),': Elapsed time (sec): ', toc(t1); \
end block


#define TSTOP()                                 \
  write(*,'(A,A,1i4,A,F20.10)') __FILENAME__,":",__LINE__,": Elapsed time (sec) ", toc()

#define TBREAK(TSTRING)                         \
  write(*,'(A,A,1i4,A,A,F20.10)') __FILENAME__,":",__LINE__,TSTRING," time (sec)" toc()

#ifndef PGI
#define START_TIMER_LOOP(NLOOP)  \
  block;\
  character(len=80) :: cblock;      \
  call timer_loop_start(NLOOP);     \
do

#define STOP_TIMER_LOOP_(COMMENT) \
  if (.not.in_timer_loop()) exit; \
  end do; \
  cblock=trim(COMMENT); \
  call timer_loop_end(trim(cblock));end block
#else
#define START_TIMER_LOOP(NLOOP)                 \
  call timer_loop_start(NLOOP);                 \
  do

#define STOP_TIMER_LOOP_(COMMENT)               \
  if (.not.in_timer_loop()) exit;               \
  end do;                                       \
  call timer_loop_end(#COMMENT)
#endif



#define TIMER_BLOCK_LOOP(NLOOP,COMMENT,YBLOCK) block;      \
  character(len=*) :: cblock=trim(COMMENT);                \
  call timer_loop_start(NLOOP);do i=1,NLOOP; YBLOCK;\
  if (.not.in_timer_loop()) exit; \
end do;\
call timer_loop_end(cblock);end block



#if defined(OPENMP) || defined(_OPENMP)

#define TBLOCK_omp()                             \
  print *,"Start timer: ", tic_omp()

#define TSTOP_omp()                              \
  write(*,'(A,A,1i4,A,1d20.10)') __FILENAME__,":",__LINE__,": Elapsed time (s) ", toc_omp()

#endif
  /*OPENMP*/


#if defined(CUDA)
#define GPU_STOPWATCH( EVALBLOCK )              \
use simple_timer_cuda;                          \
type(cudaEvent) :: __t1;                        \
real(timer_int_kind) :: elapsed;  \
starttime = tic();                              \
EVALBLOCK  ;                                    \
  print *, "__FILENAME__:__LINE__: CUDA Event timer (sec): ", toc(__t1)

#endif

#define CPU_STOPWATCH( EVALBLOCK )              \
double, precision :: __t1,__t2;                 \
call cpu_time(__t);                             \
EVALBLOCK  ;                                    \
call cpu_time(__t2);                            \
print *, "__FILENAME__:__LINE__: CPU wall time (sec): ", __t2-__t1


#endif

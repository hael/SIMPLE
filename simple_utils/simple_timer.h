
/*
  C/FORTRAN  preprocessor macros for timing module blocks
  Note: __PRETTY_FUNCTION__ is GNU specific and not standard
        __FILE__ is standard but is different some compilers
        __LINE__ is defined in the standard


     2017-04-08        Michael Eager (michael.eager@monash.edu)
*/

#if !defined(SIMPLE_TIMER_H)
#define SIMPLE_TIMER_H


#define TBLOCK(XX)                              \
  XX = tic();
!  write(*,'(A)') "Begin TBLOCK"

#define TSTOP(XX)                                                       \
  write(*,'(A,A,1i4,A,1d20.10)') __FILE__,":",__LINE__,":Elapsed time (s) ", toc()


/*write(*,'(A,1d20.10)') "__PRETTY_FUNCTION__:__FILE__:__LINE__: Elapsed time (sec): "  REAL(etime,dp)*/



/*
#define TBLOCK( EVALBLOCK )                \
use simple_timer;                            \
integer, precision :: __t1;                  \
real,precision :: elapsed;                   \
starttime = tic();                           \
EVALBLOCK  ;                                 \
write(*,"('__PRETTY_FUNCTION__:__FILE__:__LINE__: \
Elapsed time (sec): ' , 1d20.10,/)" ) toc(__t1)
 now your watch has ended */


#define TBLOCK2( EVALBLOCK , SUBSTR )        \
use simple_timer;                            \
integer, precision :: __t1;                  \
real,precision :: elapsed;                   \
starttime = tic();                           \
EVALBLOCK  ;                                 \
print *, "SUBSTR:__LINE__: Elapsed time (sec): ", toc(__t1)
/* now your watch has ended */




#if defined CUDA
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



#if defined OPENMP
#define OMP_STOPWATCH( EVALBLOCK )                 \
use omp_lib;                                        \
double, precision :: __t1,__t2;                     \
__t1 = omp_get_wtime();                             \
EVALBLOCK  ;                                        \
__t2 = omp_get_wtime();                             \
print *, "__FILE__:__LINE__: OpenMP wall time (sec): ", __t2-__t1
#endif

#endif

program simple_test_new_corrcalc
use simple_defs
use simple_cmdline,            only: cmdline
use simple_stat,               only: pearsn
use simple_math,               only: euclid,sum_diff,rel_sum_diff
use simple_hadamard3D_matcher, only: preppftcc4align, pftcc
use simple_polarft_corrcalc,   only: polarft_corrcalc
use simple_params,             only: params
use simple_build,              only: build
use simple_timing
use invert_cpu_utils
use matmul_cpu_utils
!gpu modules and handlers
use simple_cuda_defs
use simple_cuda
use matmul_gpu_utils
use invert_gpu_utils
use simple_math_gpu
implicit none
!----------------------------------------------------------
type(params)                  :: p
type(build)                   :: b
type(cmdline)                 :: cline
!type(cuda)                   :: cudaQ !potential class replacement
!CUDA err variable for the return function calls
integer                       :: err
!correlation variables
real(sp)                      :: corr, sdiff, rel_sdiff, dist
integer                       :: iptcl, cnt, iref
real, allocatable             :: corrmat3d_new(:,:,:), corrmat3d_old(:,:,:)
real, allocatable             :: corrmat2d_new(:,:), corrmat2d_old(:,:), corrmat2d_tmp(:,:)
integer, allocatable          :: inplmat2d_new(:,:), inplmat2d_old(:,:), inplmat2d_tmp(:,:)
integer, allocatable          :: indices(:)
!timer variables
real(sp)                      :: speedup
double precision              :: elps_T
double precision,dimension(2) :: st_T_r, et_T_r
double precision              :: elps_D
double precision,dimension(2) :: st_D_r, et_D_r
double precision              :: elps_A
double precision,dimension(2) :: st_A_r, et_A_r
double precision              :: elps_corr_cpu
double precision,dimension(2) :: st_corr_cpu, et_corr_cpu
double precision              :: elps_L
double precision,dimension(2) :: st_L_r, et_L_r
!indexers
integer                       :: ipart,jpart
!----------------------------------------------------------
call timestamp()
call start_Alltimers_cpu()
call simple_cuda_init(err)
if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'

write(*,*) "has_gpu: ",has_gpu

if( command_argument_count() < 3 )then
   write(*,'(a)',advance='no') 'simple_test_new_corrcalc stk=<stack.ext> vol1=<invol.ext>'
    write(*,'(a)',advance='no') ' nspace=<nr of ptcls> lp=<low-pass limit{20}> smpd=<sampling'
    write(*,'(a)') ' distance(in A)> msk=<mask radius(in pixels)> use_gpu=<yes|no{no}>'
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('vol1',   2)
call cline%checkvar('nspace', 3)
call cline%checkvar('lp',     4)
call cline%checkvar('smpd',   5)
call cline%checkvar('msk',    6)
call cline%checkvar('use_gpu',7)
call cline%check
call cline%set('prg', 'test_new_corrcalc')
p = params(cline)                   ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
p%eo = 'no'                  ! default
!write(*,*) "Allocating memory "
allocate( corrmat3d_new(p%nptcls,p%nptcls,p%nrots),&
          corrmat3d_old(p%nptcls,p%nptcls,p%nrots),&
          corrmat2d_tmp(p%nptcls,p%nptcls),&
          corrmat2d_new(p%nptcls,p%nptcls),&
          corrmat2d_old(p%nptcls,p%nptcls),&
          inplmat2d_tmp(p%nptcls,p%nptcls),&
          inplmat2d_new(p%nptcls,p%nptcls),&
          inplmat2d_old(p%nptcls,p%nptcls),&
          indices(p%nptcls))
          
!write(*,*)"Before build Hadamard prime tbox" 
call b%build_hadamard_prime3D_tbox(p)
!write(*,*)"After build Hadamard prime tbox" 
!write(*,*)"Before preppftcc4align" 

call start_timer_cpu("preppftcc4align")
call preppftcc4align(b, p, cline)
call stop_timer_cpu("preppftcc4align")
write(*,*)'*****************CPU corr****************************************'
write(*,*)"Before gencorrs_all_tester CPU"
write(*,'(x,a,x,i3)')'Number of OpenMP threads: ',p%nthr
write(*,*)"-----------------------------------------------------"
#if defined (BENCH)
    call gettimeofday_c(st_T_r)
#endif
call pftcc%gencorrs_all_tester(corrmat3d_old)
#if defined (BENCH)
call gettimeofday_c(et_T_r)
call elapsed_time_c(st_T_r,et_T_r,elps_T)
#endif

!write(*,*)"Before expand_dim"

#if defined (BENCH)
call gettimeofday_c(st_D_r)
#endif
call pftcc%expand_dim
#if defined (BENCH)
call gettimeofday_c(et_D_r)
call elapsed_time_c(st_D_r,et_D_r,elps_D)
#endif

call pftcc%print()

write(*,*)'*****************GPU corr****************************************'
write(*,*)"Before gencorrs_all GPU if has_gpu: ",has_gpu,", use_gpu: ",use_gpu
write(*,*)"-----------------------------------------------------"
#if defined (BENCH)
call start_timer_cpu("genAll")
call gettimeofday_c(st_A_r)
#endif
if( p%use_gpu .eq. 'yes' )then
  call pftcc%gencorrs_all_gpu(p,corrmat3d_new)
else
  call pftcc%gencorrs_all_cpu(corrmat3d_new)
endif
#if defined (BENCH)
    call gettimeofday_c(et_A_r)
    call elapsed_time_c(st_A_r,et_A_r,elps_A)
    call stop_timer_cpu("genAll")
#endif
write(*,*)'*****************************************************************'

! reverting to the old mode of comparison

!$omp parallel default(shared) private(iptcl)
!$omp workshare 
corrmat2d_tmp = maxval(corrmat3d_new, dim=3)
inplmat2d_tmp = maxloc(corrmat3d_new, dim=3)
corrmat2d_old = maxval(corrmat3d_old, dim=3)
inplmat2d_old = maxloc(corrmat3d_old, dim=3)
!$omp end workshare
!$omp end parallel

! initialise validation variables
! corr      = 0.
! sdiff     = 0.
! rel_sdiff = 0.
! dist      = 0.

call get_corrmat_inplmat_2D(p,corrmat3d_new, &
                            corrmat2d_tmp,inplmat2d_tmp, &
                            corrmat2d_new,inplmat2d_new)

do ipart=1,p%nptcls
   do jpart=1,p%nptcls
      write(2000,*)corrmat2d_new(ipart,jpart),corrmat2d_old(ipart,jpart)
   end do
end do

#if defined (CUDA) && defined (MAGMA)
write(*,*)"switching to difference and relative difference comparisons"
write(*,*) 'corr,         corr=', &
     pearsn(reshape(corrmat2d_new, shape=[p%nptcls**2]),             &
     reshape(corrmat2d_old, shape=[p%nptcls**2]))
write(*,*) 'corr,     sum_diff=', &
     sum_diff(reshape(corrmat2d_new, shape=[p%nptcls**2]),           &
     reshape(corrmat2d_old, shape=[p%nptcls**2]))
write(*,*) 'corr, rel_sum_diff=', &
     rel_sum_diff(reshape(corrmat2d_new, shape=[p%nptcls**2]),       &
     reshape(corrmat2d_old, shape=[p%nptcls**2]))
write(*,*) 'inpl,         corr=', &
     pearsn(real(reshape(inplmat2d_new, shape=[p%nptcls**2])),       &
     real(reshape(inplmat2d_old, shape=[p%nptcls**2])))
write(*,*) 'inpl,     sum_diff=', &
     sum_diff(real(reshape(inplmat2d_new, shape=[p%nptcls**2])),     &
     real(reshape(inplmat2d_old, shape=[p%nptcls**2])))
write(*,*) 'inpl, rel_sum_diff=', &
     rel_sum_diff(real(reshape(inplmat2d_new, shape=[p%nptcls**2])), &
     real(reshape(inplmat2d_old, shape=[p%nptcls**2])))
#else
write(*,*)"switching to Euclidean distance comparisons"
write(*,*) 'corr, corr=', pearsn(reshape(corrmat2d_new, shape=[p%nptcls**2]),&
reshape(corrmat2d_old, shape=[p%nptcls**2]))
write(*,*) 'corr, euclid=', euclid(reshape(corrmat2d_new, shape=[p%nptcls**2]),&
reshape(corrmat2d_old, shape=[p%nptcls**2]))
write(*,*) 'inpl, corr=', pearsn(real(reshape(inplmat2d_new, shape=[p%nptcls**2])),&
real(reshape(inplmat2d_old, shape=[p%nptcls**2])))
write(*,*) 'inpl, euclid=', euclid(real(reshape(inplmat2d_new, shape=[p%nptcls**2])),&
real(reshape(inplmat2d_old, shape=[p%nptcls**2])))
#endif

#if defined (BENCH) && defined (LINUX)
speedup = real(elps_T) / real(elps_A)
write(*,*)"gencorrs_all_tester     expand_dim     gencorrs_all"
write(*,*)"---------------------------------------------------"
write(*,'(f15.8,5x,f15.8,x,f15.8)') elps_T,elps_D,elps_A
write(*,*)'*****************************************************************'
write(*,*) " Speed up from CPU to GPU     : ",speedup
if (speedup<1.0) write(*,*)"speedup is < 1, try increasing number of particles"
write(*,*)'*****************************************************************'
#endif

! corr = corr/real(cnt)
! #if defined (CUDA) && defined (MAGMA)
!     write(*,*)"switching to difference and relative difference comparisons"
!     write(*,*) 'corr,         corr=', corr
!     write(*,*) 'corr,     sum_diff=', sdiff
!     write(*,*) 'corr, rel_sum_diff=', rel_sdiff
! #else
!     write(*,*)"switching to Euclidean distance comparisons"
!     write(*,*) 'corr, corr=',   corr
!     write(*,*) 'corr, euclid=', dist
! #endif
! #if defined (BENCH) && defined (LINUX)
!     speedup = real(elps_T)/real(elps_A)
!     write(*,*)"gencorrs_all_tester     expand_dim     gencorrs_all"
!     write(*,*)"---------------------------------------------------"
!     write(*,'(f15.8,5x,f15.8,x,f15.8)') elps_T,elps_D,elps_A
!     write(*,*)'***************************************************************'
!     write(*,*) " Speed up from CPU to GPU     : ",speedup
!     write(*,*)'***************************************************************'
! #endif

!freeing ressources on host
deallocate(corrmat3d_new)
deallocate(corrmat3d_old)
deallocate(corrmat2d_tmp)
deallocate(corrmat2d_new)
deallocate(corrmat2d_old)
deallocate(inplmat2d_tmp)
deallocate(inplmat2d_new)
deallocate(inplmat2d_old)
deallocate(indices)

!shutting down the environment
call simple_cuda_shutdown()
!shutting down the timers
call stop_Alltimers_cpu()

! contains
!
!     subroutine metrics_pft_corrs( a_new, a_old, corr, sdiff, rel_sdiff, dist )
!         real, intent(in)    :: a_new(:), a_old(:)
!         real, intent(inout) :: corr, sdiff, rel_sdiff, dist
!         if( size(a_new) /= size(a_old) )then
!             stop 'arrays 4 comparison not equal dims'
!         endif
!         corr      = corr      + pearsn(       a_new, a_old )
!         sdiff     = sdiff     + sum_diff(     a_new, a_old )
!         rel_sdiff = rel_sdiff + rel_sum_diff( a_new, a_old )
!         dist      = dist      + euclid(       a_new, a_old )
!     end subroutine

end program simple_test_new_corrcalc
!*******************************************************************************
!    Subroutine to construct the corrmat2d inherited from the GPU gencorrAll
!
!*******************************************************************************
!
subroutine get_corrmat_inplmat_2D(p,corrmat3d_in,corrmat2d_in,inplmat2d_in, &
                                                 corrmat2d_out,inplmat2d_out)
  use simple_defs
  use simple_params
  implicit none
  type(params) :: p
  integer      :: inplmat2d_out(p%nptcls,*)
  real(sp)     :: corrmat2d_out(p%nptcls,*)
  integer      :: inplmat2d_in(p%nptcls,*)
  real(sp)     :: corrmat2d_in(p%nptcls,*)
  real(sp)     :: corrmat3d_in(p%nptcls,p%nptcls,*)
  integer      :: indices(p%nptcls)
  !timer variables
  double precision              :: elps_L
  double precision,dimension(2) :: st_L_r, et_L_r
  !indexers
  integer                       :: iptcl, cnt, iref
  !start of the execution commands

  !$omp parallel default(shared) private(iptcl)
  !$omp do schedule(auto)
  do iptcl=1,p%nptcls
     indices(iptcl) = iptcl ! generate indices for cshifting
  end do
  !$omp end do nowait

  ! initialise validation variables
  ! corr      = 0.
  ! sdiff     = 0.
  ! rel_sdiff = 0.
  ! dist      = 0.

  !$omp end parallel

  cnt       = 0
#if defined (BENCH)
  call gettimeofday_c(st_L_r)
#endif
  do iref=1,p%nptcls
     if( iref /= 1 )then
        indices = cshift(indices, shift=1)
     endif
     !$omp parallel do schedule(auto) private(iptcl)
     do iptcl=1,p%nptcls
        corrmat2d_out(indices(iptcl),iref) = corrmat2d_in(iref,iptcl)
        inplmat2d_out(indices(iptcl),iref) = inplmat2d_in(iref,iptcl)
        cnt = cnt+1
     end do
     !$omp end parallel do
  end do
#if defined (BENCH)
  call gettimeofday_c(et_L_r)
  call elapsed_time_c(st_L_r,et_L_r,elps_L)
  write(*,*)"---------------------------------------------------"
  write(*,'(a,x,f15.8)')"Matrix mapping (secs): ", elps_L
  write(*,*)"---------------------------------------------------"
#endif

  return
end subroutine get_corrmat_inplmat_2D

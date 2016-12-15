program simple_test_cartesian_corr
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use simple_timing
  use simple_cmdline,     only: cmdline
  use simple_params,      only: params
  use simple_build,       only: build
  use simple_image,       only: image
  use simple_ft_expanded
  use simple_syscalls
  implicit none

  integer, parameter   :: NITS=1
  integer, parameter   :: vx = 4096, vy = 4096, vz = 1

  type(params)         :: p
  type(build)          :: b
  type(cmdline)        :: cline
  type(image)          :: img1, img2
  type(ft_expanded)    :: ftexp1, ftexp2
  integer, allocatable :: phys(:,:,:,:)
  logical, allocatable :: lmsk(:,:,:)
  integer              :: lims(3,2), i
  real, parameter      :: hp=100.0, lp=8.0, shvec(3)=[0.,0.,0.]
  real                 :: corr
  real(sp)             :: speedup
  !CUDA err variable for the return function calls
  integer                       :: err
  !timer variables
  double precision              :: elps_T
  double precision,dimension(2) :: st_T_r, et_T_r
  double precision              :: elps_D
  double precision,dimension(2) :: st_D_r, et_D_r
  double precision              :: elps_S
  double precision,dimension(2) :: st_S_r, et_S_r
  
  !start of the execution commands
  call timestamp()
  call start_Alltimers_cpu()
  call simple_cuda_init(err)
  if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'

  write(*,*) "has_gpu: ",has_gpu

  if( command_argument_count() < 2 )then
     write(*,'(a)',advance='no') 'simple_test_cartesian_corr '
     write(*,'(a)',advance='no') 'use_gpu=<yes|no{no}> '
     write(*,'(a)')              'nthr=<nr of OpenMP threads{1}>'
     stop
  endif
  call cline%parse
  call cline%checkvar('use_gpu',1)
  call cline%checkvar('nthr',   2)
  call cline%check
  p = params(cline)                   ! parameters generated
  call b%build_general_tbox(p, cline) ! general objects built
  p%eo = 'no'                         ! default

  ! make random images
  call img1%new([vx,vy,vz], 1.77)
  call img1%ran
  call img2%new([vx,vy,vz], 1.77)
  call img2%ran
  ! prepare objects for corrcalc
  call ftexp1%new(img1, hp, lp)
  call ftexp2%new(img2, hp, lp)
  lims = img1%loop_lims(1,lp)
  ! time the different routines
  call start_timer_cpu("old")
  call gettimeofday_c(st_T_r)
  do i=1,NITS
     corr  = img1%corr(img2, lp) 
  end do
  call gettimeofday_c(et_T_r)
  call elapsed_time_c(st_T_r,et_T_r,elps_T)
  call stop_timer_cpu("old")

  call start_timer_cpu("recast")
  call gettimeofday_c(st_D_r)
  do i=1,NITS
     corr  = ftexp1%corr(ftexp2)
  end do
  call gettimeofday_c(et_D_r)
  call elapsed_time_c(st_D_r,et_D_r,elps_D)
  call stop_timer_cpu("recast")

  call start_timer_cpu("shifted")
  call gettimeofday_c(st_S_r)
  do i=1,NITS
     corr = ftexp1%corr_shifted(ftexp2, shvec)
  end do
  call gettimeofday_c(et_S_r)
  call elapsed_time_c(st_S_r,et_S_r,elps_S)
  call stop_timer_cpu("shifted")
  
  speedup = elps_T/elps_D
  write(*,*)'*****************************************************************'
  write(*,*)"Correlations            old            re-cast     shifted"
  write(*,*)"-----------------------------------------------------------------"
  write(*,'(18x,f15.8,5x,f15.8)') elps_T,elps_D,elps_S
  write(*,*)'*****************************************************************'
  write(*,*) " Speed up from old to re-cast(GPU)   : ",speedup
  if (speedup<1.0) write(*,*)"speedup is < 1, try increasing nmbr of particles"
  write(*,*)'*****************************************************************'

  !shutting down the environment
  call simple_cuda_shutdown()
  !shutting down the timers
  call stop_Alltimers_cpu()
  
end program simple_test_cartesian_corr

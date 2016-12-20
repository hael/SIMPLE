! ============================================================================
! Name        : testing_congruence.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the corr calc calculation with CUDA
!             :
! ============================================================================
!
program testing_congruence
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu

  implicit none

#define devptr_t integer*8
#define verbose .true.
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD

  integer, parameter            :: nptcls = 63119
  real(sp), parameter           :: tol = 0.1 !%tage of tolerance in the level
  
  real(sp), parameter           :: range_chunksz_init = 45.e6
  real(sp), parameter           :: range_chunksz_fnit = 200.e6

  integer, parameter            :: chunksz=3000
  integer, parameter            :: start_chunksz=2000, istep_chunksz=16

  integer, parameter            :: nrot =478
  integer, parameter            :: start_nrot=439, istep_nrot=1

  integer, parameter            :: nk = 100
  integer, parameter            :: start_nk=100, istep_nk=1
  !GPU optimastion problem
  integer, parameter            :: nthreads = 256 !number of threads
  !ressources available
  integer, parameter            :: nnodes = 14
  integer, parameter            :: ncores = 16

  integer                       :: err
  integer                       :: lda
  integer                       :: ldb,ldc
  !local variables
  integer                       :: nCPU_cores
  integer                       :: ndev
  integer                       :: N !total number of elemts chunksz*nrot*nk 
  integer                       :: tcores = ncores*nnodes

  !variables for the function 

  integer                       :: klp         !< low-pass frequency limit
  integer                       :: khp         !< high-pass frequency limit 

  real(sp)                      :: alpha,beta  
  real(sp)                      :: r, sumasq, sumbsq
  integer, allocatable          :: mask_pft(:,:,:)

  integer, allocatable          :: mask_chunksz(:)
  integer, allocatable          :: mask_nrot(:)
  integer, allocatable          :: mask_nk(:)
  integer, allocatable          :: mask_size(:)
  integer, allocatable          :: mask_required_cores(:)
  real(sp), allocatable         :: mask_rejct(:)

  complex(sp), allocatable      :: pft1(:,:,:)
  complex(sp), allocatable      :: pft2(:,:,:)

  !gpu variables for the calculation of the r, sumasq and sumbsq

  real(sp)                      :: r_gpu, sumasq_gpu, sumbsq_gpu
  real(sp)                      :: dble_alpha,dble_beta

  !timer variables

  double precision              :: elps_corr_FN_gpu
  double precision,dimension(2) :: st_FN_r, et_FN_r

  double precision              :: elps_corr_NN_gpu
  double precision,dimension(2) :: st_NN_r, et_NN_r

  double precision              :: elps_corr_cpu
  double precision,dimension(2) :: st_corr_cpu, et_corr_cpu

  !indexers
  integer                       :: ifound_range
  !counters
  integer                       :: nfound, nfound_range,count_min

  !interfaces
  interface
     subroutine get_nfound_range(nfound, nfound_range,                 &
                                 start_chunksz, chunksz, istep_chunksz,      &
                                 start_nrot, nrot, istep_nrot,         &
                                 start_nk, nk, istep_nk,               &
                                 nthreads, nptcls,                     &
                                 range_chunksz_init,range_chunksz_fnit,    &
                                 mask_chunksz,mask_nrot,mask_nk,         &
                                 mask_size, mask_required_cores,       &
                                 mask_rejct)
       use simple_defs
       implicit none
       integer,intent(out)            :: nfound, nfound_range
       integer,intent(in)             :: chunksz, start_chunksz, istep_chunksz
       integer,intent(in)             :: nrot, start_nrot, istep_nrot
       integer,intent(in)             :: nk, start_nk, istep_nk
       integer,intent(in)             :: nthreads
       integer,intent(in)             :: nptcls
       real(sp),intent(in)            :: range_chunksz_init
       real(sp),intent(in)            :: range_chunksz_fnit
       integer,intent(inout),optional :: mask_chunksz(*)
       integer,intent(inout),optional :: mask_nrot(*)
       integer,intent(inout),optional :: mask_nk(*)
       integer,intent(inout),optional :: mask_size(*)
       integer,intent(inout),optional :: mask_required_cores(*)
       real(sp),intent(inout),optional :: mask_rejct(*)
     end subroutine get_nfound_range
  end interface
  
  !start of the execution commands
  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  call sysQ%new_systemQuery_cpu(hstD)

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(sysQ) and the data structure(hstD) '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  write(*,*)'Number of cores on the system                    : ',nCPU_cores
  !writting the number of cores avaible on the system
  write(*,*)'Number of cores returned from data structure hstD: ',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()
  write(*,*)'******************************************************************'

!*******************************************************************************
!     now testing for congruence
!
!*******************************************************************************

  write(*,*)'                                                               '
  write(*,*)'***************************************************************'
  write(*,*)'     now scanning for size of optimal size and factor of 256   '
  write(*,*)'***************************************************************'
  write(*,*)
  write(*,*)"Sumstack n particlesSize : ",nptcls
  write(*,*)"Number of nodes          : ",nnodes
  write(*,*)"Number of cores          : ",ncores
  write(*,*)"Total number of cores    : ",tcores
  write(*,*)"Size of blocks in threads: ",nthreads
  write(*,*)"Range for chunksz          : [1,",chunksz,"]"
  write(*,*)"Range for nrot           : [1,",nrot,"]"
  write(*,*)"Range for nk             : [1,",nk,"]"
  write(*,*)"Floating point window    : [",range_chunksz_init,",",range_chunksz_fnit,"]"
  write(*,*)"In steps of istep_chunksz  :  ",istep_chunksz
  write(*,*)"In steps of istep_nrot   :  ",istep_nrot
  write(*,*)"In steps of istep_nk     :  ",istep_nk
  write(*,*)
  write(*,*) "scanning for value..."
  write(*,'(9x,a,4x,a,4x,a,5x,a,4x,a,5x,a,5x,a,2x,a)') &
       "ifound","ipart","irot","nk","# elms in pft(ipart,irot,nk)", &
       "(Millions)","cores needed","Thrown particles"

  !scanning for values

  call start_timer_cpu("scan_vals")
  call gettimeofday_c(st_FN_r)

  call get_nfound_range(nfound, nfound_range,            &
                        start_chunksz, chunksz, istep_chunksz, &
                        start_nrot, nrot, istep_nrot,    &
                        start_nk, nk, istep_nk,          &
                        nthreads, nptcls,                &
                        range_chunksz_init,range_chunksz_fnit)

  call gettimeofday_c(et_FN_r)
  call elapsed_time_c(st_FN_r,et_FN_r,elps_corr_FN_gpu)
  call stop_timer_cpu("scan_vals")

  write(*,*)"Total found in the range: ",nfound_range
  write(*,*)"Total found             : ",nfound
  write(*,*)"Elapsed time for scan   : ",real(elps_corr_FN_gpu),"(seconds)"
  write(*,*)
  
  !alocating the mask value
  allocate(mask_chunksz(nfound_range))
  allocate(mask_nrot(nfound_range))
  allocate(mask_nk(nfound_range))
  allocate(mask_size(nfound_range))
  allocate(mask_required_cores(nfound_range))
  allocate(mask_rejct(nfound_range))
  
  !Initialising the mask
  mask_chunksz = 0
  mask_nrot = 0
  mask_nk = 0
  mask_size = 0
  
  call start_timer_cpu("scan_range")
  call get_nfound_range(nfound, nfound_range,                 &
                        start_chunksz, chunksz, istep_chunksz,      &
                        start_nrot, nrot, istep_nrot,         &
                        start_nk, nk, istep_nk,               &
                        nthreads, nptcls,                     &
                        range_chunksz_init,range_chunksz_fnit,    &
                        mask_chunksz,mask_nrot,mask_nk,         &
                        mask_size, mask_required_cores,       &
                        mask_rejct)
  call stop_timer_cpu("scan_range")

!  write(*,'(3x,i5,5x,i5,4x,i5,2x,i5,10x,i10,13x,f15.8)') &
!       nfound, ipart, irot, nk, ipart*nrot*nk, ipart*nrot*nk/1.e6
  write(*,*) "Found min reject..."
  write(*,'(9x,a,4x,a,4x,a,5x,a,4x,a,5x,a,5x,a,2x,a)') &
       "ifound","ipart","irot","nk","# elms in pft(ipart,irot,nk)", &
       "(Millions)","cores needed","Thrown particles"

  count_min = 0
  do ifound_range=1,nfound_range
     if (mask_rejct(ifound_range) <= minval(mask_rejct)+(nptcls*tol/100.0) ) then
        count_min = count_min + 1
     end if
  end do

  count_min = 0
  do ifound_range=1,nfound_range

     if (mask_rejct(ifound_range) <= minval(mask_rejct)+(nptcls*tol/100.0) ) then
        count_min = count_min + 1
        write(*,'(3x,i10,5x,i8,4x,i5,2x,i5,10x,i10,13x,f15.8,5x,i6,7x,f15.4)') &
             count_min,                          &
             mask_chunksz(ifound_range),           &
             mask_nrot(ifound_range),            &
             mask_nk(ifound_range),              &
             mask_size(ifound_range),            &
             mask_size(ifound_range)/1.e6,       &
             mask_required_cores(ifound_range),  &
             mask_rejct(ifound_range)
     end if

  end do
  
  !Freeing the masks on CPU memory 
  deallocate(mask_chunksz)
  deallocate(mask_nrot)
  deallocate(mask_nk)
  deallocate(mask_size)
  deallocate(mask_required_cores)
  deallocate(mask_rejct)
  
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_congruence
!*******************************************************************************
! DESCRIPTION
! Subroutine to get the nfound and nfound_range values by scaning possibilities
!
!*******************************************************************************
! SOURCE
subroutine get_nfound_range(nfound, nfound_range,                 &
                            start_chunksz, chunksz, istep_chunksz,      &
                            start_nrot, nrot, istep_nrot,         &
                            start_nk, nk, istep_nk,               &
                            nthreads,nptcls,                      &
                            range_chunksz_init,range_chunksz_fnit,    &
                            mask_chunksz,mask_nrot,mask_nk,         &
                            mask_size,mask_required_cores,        &
                            mask_rejct)
  use simple_defs
  implicit none
  integer,intent(out)            :: nfound, nfound_range
  integer,intent(in)             :: chunksz, start_chunksz, istep_chunksz
  integer,intent(in)             :: nrot, start_nrot, istep_nrot
  integer,intent(in)             :: nk, start_nk, istep_nk
  integer,intent(in)             :: nthreads
  integer,intent(in)             :: nptcls
  real(sp),intent(in)            :: range_chunksz_init
  real(sp),intent(in)            :: range_chunksz_fnit
  integer,intent(inout),optional :: mask_chunksz(*)
  integer,intent(inout),optional :: mask_nrot(*)
  integer,intent(inout),optional :: mask_nk(*)
  integer,intent(inout),optional :: mask_size(*)
  integer,intent(inout),optional :: mask_required_cores(*)
  real(sp),intent(inout),optional :: mask_rejct(*)
  !local variable
  integer                       :: N
  !indexers
  integer                       :: ipart,irot, ik, ithreads

  !scanning for values
  nfound = 0
  nfound_range = 0
  do irot = start_nrot, nrot, istep_nrot
     do ipart = start_chunksz, chunksz, istep_chunksz
        do ik = start_nk, nk, istep_nk

           N = ipart * irot * nk
           if ( modulo(N,nthreads) == 0 ) then
              nfound = nfound + 1

              if ( range_chunksz_init < N .and. N < range_chunksz_fnit ) then
                 nfound_range = nfound_range + 1
                 if ( verbose ) then
                    if ( .not. ( present(mask_chunksz) .and. &
                         present(mask_nrot) .and.  &
                         present(mask_nk) ) ) then
                       write(*,'(3x,i10,5x,i8,4x,i5,2x,i5,10x,i10,13x,f15.8,5x,i6,7x,f15.4)') &
                            nfound_range, ipart, irot, ik, N, N/1.e6, nptcls/ipart,(real(nptcls)/real(ipart)-nptcls/ipart)*ipart

                    end if
                 end if

                 if ( present(mask_chunksz) .and. &
                      present(mask_nrot)  .and. &
                      present(mask_nk)          ) then

                    mask_chunksz(nfound_range) = ipart
                    mask_nrot(nfound_range)  = irot
                    mask_nk(nfound_range)    = ik
                    mask_size(nfound_range)  = N
                    mask_required_cores(nfound_range) = nptcls/ipart
                    mask_rejct(nfound_range) = (real(nptcls)/real(ipart)-nptcls/ipart)*ipart

                 end if

              end if

           end if
        end do
     end do
  end do

  return
end subroutine get_nfound_range

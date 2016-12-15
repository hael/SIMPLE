! ============================================================================
! Name        : testing_discrepancy.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 07th of October 2015
! Description : tests the discrepancy between a CPU and GPU calculation
!             :
! ============================================================================
!
program testing_discrepancy
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use greeting_version
  use simple_timing

  implicit none

  real(sp), parameter           :: epsilon=FTOL/10 !0.00001

  character(len=30)             :: filename1 = 'AAstar_cpu_GNUf.log'
  character(len=30)             :: filename2 = 'AAstar_gpu_CUDA.log'

  integer, parameter            :: npart=2141
  integer, parameter            :: start_npart=2141, istep=125
  integer, parameter            :: nrot = 360
  integer, parameter            :: nk = 59

  integer                       :: ix,iy,iz

  real(sp), allocatable         :: a_cpu(:,:,:)
  real(sp), allocatable         :: a_gpu(:,:,:)
  real(sp), allocatable         :: rel_diff(:,:,:)

  !indexers
  integer                       :: ipart,irot,ik
  !counters
  integer                       :: count, count_diff
  !start of the execution commands
  !start of the greeting message
  call hello_cpu_discrepancy_checker()
  call timestamp()
  call start_Alltimers_cpu()

  write(*,*)'                                                               '
  write(*,*)'***************CPU vs GPU discrepancy analyser******************'
  write(*,*)'                                                               '

  !allocating the complex matrices
  allocate(a_cpu(npart,nrot,nk))
  allocate(a_gpu(npart,nrot,nk))
  allocate(rel_diff(npart,nrot,nk))

  !openeing the files in question
  open(1,file=filename1,action='read')
  open(2,file=filename2,action='read')

  !reading in the files in questions

  write(*,'(4x,a,4x,a,5x,a)') "inpart","nrot","nk"
  write(*,'(7(3x,i5))') npart,nrot,nk
  write(*,'(x,a,x,f15.8)') "Tolerance in the relative difference  : ",epsilon
  write(*,*) "Number of elements in a(npart,nrot,nk): ",npart*nrot*nk
  write(*,*) "Number of elements in a(npart,nrot,nk): ",npart*nrot*nk/1.e6,"(Millions)"
  write(*,*) 
  write(*,*)'               Checking difference....                          '
  write(*,*) 
  write(*,'(6x,a,11x,a,10x,a,10x,a,6x,a,12x,a,12x,a)') &
       "entry","ix", "iy", "iz","a_cpu", "a_gpu", "rel_diff"

  count = 0
  count_diff = 0
  do ipart=1,npart!-(npart-2)
     do irot=1,nrot!-(nrot-2)
        do ik=1,nk!-(nk-2)
           read(1,*)ix, iy, iz,a_cpu(ipart,irot,ik)
           read(2,*)ix, iy, iz,a_gpu(ipart,irot,ik)
           
           rel_diff(ipart,irot,ik) = a_cpu(ipart,irot,ik) - a_gpu(ipart,irot,ik)

           if ( rel_diff(ipart,irot,ik) > epsilon ) then
              write(*,*) count,ix, iy, iz,      &
                   a_cpu(ipart,irot,ik),  &
                   a_gpu(ipart,irot,ik),  &
                   rel_diff(ipart,irot,ik)
              count_diff = count_diff + 1
           end if
           count = count + 1
        end do
     end do
  end do

  if ( count_diff == 0 ) then
     write(*,*)"The relative diff is within the tolerance of: ",epsilon
     write(*,*)"Number of elements found outside tolerance  : ",count_diff
  end if
  write(*,*)
  write(*,*)'           ....Checked difference.                       '
  write(*,*)
  write(*,*)"Full sum of a_cpu is    : ",sum(a_cpu)
  write(*,*)"Full sum of a_gpu is    : ",sum(a_gpu)

  rel_diff(:,:,:) = a_cpu(:,:,:) - a_gpu(:,:,:)

  close(1)
  close(2)

  !freeing the ressources
  deallocate(a_cpu)
  deallocate(a_gpu)
  deallocate(rel_diff)

  !shutting down the timers
  call stop_Alltimers_cpu()
  !end of greeting message
  call bye_cpu_discrepancy_checker()

end program testing_discrepancy



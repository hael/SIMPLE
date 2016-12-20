! ============================================================================
! Name        : testing_yaml.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 11th of February 2016
! Description : tests the yaml output and so on
!             :
! ============================================================================
!
program testing_DataParallelIO
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu
  use simple_yaml_output
  use simple_yaml_strings
  !file and resources data structure defs
  use simple_file_defs
  use simple_resources_defs
  !parralel IO module 
  use simple_DataParallel_IO
  
  implicit none

#define devptr_t integer*8
#define verbose .true.
  type(systemQuery_cpu)         :: sysQ            !system details query object
  type(systemDetails)           :: hstD            !systemDetails data structure
  type(polar_corr_calc)         :: s_polar         !polar ft data structure
  type(fileDetails_Parallel_IO) :: fpioD           !fle details for parrellel IO
  type(resources_avail)         :: resD            !ressources for the system
  integer, parameter            :: nptcls = 100000

  integer, parameter            :: npart=32000
  integer, parameter            :: start_npart=32000, istep_npart=64

  integer, parameter            :: nrot =3140
  integer, parameter            :: start_nrot=3140, istep_nrot=1

  integer                       :: ikernel = 6           !Kernel option
  integer                       :: threadsPerBlock = 512 !threadsPerBlock
  integer                       :: nx_3D=16, ny_3D=16, nz_3D=4 !3D mesh
  
  integer                       :: omp_nthr=12           !openMP variables

  !ressources available
  integer, parameter            :: nnodes = 10
  integer, parameter            :: size_socket_logi = 12
  integer, parameter            :: size_socket_phys = 6

  integer                       :: lda
  integer                       :: ldb,ldc
  integer                       :: size_pft1, size_pft2  
  complex(sp), allocatable      :: pft1(:,:)
  complex(sp), allocatable      :: pft2(:,:)

  !Yaml imput parameters for the fucntion code
  integer                       :: unit = 1
!  character(len=*)              :: mapname = 'mapname'
!  character(len=*)              :: mapvalue
!  character(len=*)              :: label
!  character(len=*)              :: tag
!  character(len=*)              :: advance
  !index variable major
  integer                       :: ipart
  !indexers
  integer                       :: i,j,jpart,ik,irot,iptcl
  !file details
  !file_open(file,unit,status,position,action,binary)
  character(len=80) :: file1='pft1'
  character(len=80) :: filename

!*******************************************************************************
!     start of the execution commands
!*******************************************************************************

  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  call Sanity_check_cpu(sysQ, hstD)
  write(*,*)'******************************************************************'
  
!*******************************************************************************
!     now testing for the YAML output
!
!*******************************************************************************

  write(*,*)'                                                                  '
  write(*,*)'********************* YAML output ********************************'
  write(*,*)'                                                                  '
  
  call yaml_map('mapname',hstD%nCPUcores,fmt='(es24.17)')
  
  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'

!*******************************************************************************
!     now let's generates some data for using as input arrays
!
!*******************************************************************************
  do ipart = start_npart, npart, istep_npart

     lda = ipart
     ldb = nrot
     ldc = ipart

!*******************************************************************************
!    Data allocation and initialisation of the matrices for testing purposes.
!
!*******************************************************************************

     !cpu allocation
     allocate( pft1(ipart,nrot) ) !pfts_refs
     allocate( pft2(ipart,nrot) ) !pfts_ptcls

     pft1 = 0.0
     pft2 = 0.0

!*******************************************************************************
!    Filling and printing PFTs array with random numbers using test functions
!
!*******************************************************************************
     call getCRandomGaussianDistr_2D( ipart, nrot, pft1)
     call getCRandomGaussianDistr_2D( ipart, nrot, pft2)

     write(*,'(10x,a,5x,a,10x,a,16x,a,26x,a)')"ipart","irot", &
                                              "(pft1)","(pft2)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
              write(*,*)i, irot,ik,pft1(i,irot), pft2(i,irot)
        end do
     end do
     write(*,'(50x,a,12x,a)')"conjg(pft2)","real(pft1(i1,irot,ik)*conjg(pft2(i2,irot,ik)))"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
              write(*,*)i, irot, ik, conjg(pft2(i,irot)) , &
                   real( pft1(i,irot)*conjg(pft2(i,irot)) ), &
                   real(pft1(i,irot)) * real(pft2(i,irot)) + &
                   imag(pft1(i,irot)) * imag(pft2(i,irot)) 
        end do
     end do
!*******************************************************************************
!     now testign the corr function on both CPU(OpenMP) vs GPU single thread
!
!*******************************************************************************
     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     Now testing the Parrallel IO                              '
     write(*,*)'***************************************************************'
     size_pft1 = ipart*nrot
     size_pft2 = ipart*nrot
     write(*,*)
     write(*,'(3x,a,3x,a,4x,a,6x,a,2x,a)') &
          "npart","ipart","nrot","OpenMP thrd"
     write(*,'(8(3x,i5))') npart,  ipart,  nrot, omp_nthr
     write(*,*) "N (Ref)   pft1(npart,nrot/2,nk): ",size_pft1
     write(*,*) "N (Ref)   pft1(npart,nrot/2,nk): ",size_pft1/1.e6,"(Millions)"
     write(*,*) "N (ptc) pft2(npart*2,nrot*2,nk): ",size_pft2
     write(*,*) "N (ptc) pft2(npart*2,nrot*2,nk): ",size_pft2/1.e6,"(Millions)"
     write(*,*) " Ratio of size  (PFT2 and PFT1): ",size_pft2 / real(size_pft1)
     write(*,*)'                                                               '
     write(*,*)'************** CPU parralel IO ********************************'
     write(*,*)'                                                               '

     call hello_DataParallel_IO()

     unit = 2 !for data file
     filename = trim(adjustl(file1))//'.dat'

     fpioD%file = filename
     fpioD%unit = unit
     fpioD%status = 'unknown'
     fpioD%position = 'asis'
     fpioD%action  = 'write'
     fpioD%binary = .true.

     call writeData(ipart,nrot,pft1,nnodes,hstD,fpioD)!filename,unit,'unknown','asis','write',.true.)

!*******************************************************************************
!    Data deallocation
!
!*******************************************************************************

     deallocate(pft1)
     deallocate(pft2)

  end do!end of do ipart = start_npart, npart, istep
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
     
  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_DataParallelIO
!
!*******************************************************************************
!    Unused code for checks of the input parameters later used for the Yaml-out
!
!*******************************************************************************
!
!GPU optimastion problem
!integer, parameter            :: nthreads = 256 !number of threads
!ressources available
!integer, parameter            :: nnodes = 14
!integer, parameter            :: ncores = 16
!local variables
!integer                       :: nCPU_cores
!integer                       :: ndev
!integer                       :: N !total number of elemts npart*nrot*nk 
!integer                       :: tcores = ncores*nnodes
!write(*,*)'                                                                  !'
!write(*,*)'******************************************************************'
!write(*,*)'     now scanning for size of optimal size and factor of 256      '
!write(*,*)'******************************************************************'
!write(*,*)
!write(*,*)"Sumstack n particlesSize : ",nptcls
!write(*,*)"Number of nodes          : ",nnodes
!write(*,*)"Number of cores          : ",ncores
!write(*,*)"Total number of cores    : ",tcores
!write(*,*)"Size of blocks in threads: ",nthreads
!write(*,*)"Range for npart          : [1,",npart,"]"
!write(*,*)"Range for nrot           : [1,",nrot,"]"
!write(*,*)"Range for nk             : [1,",nk,"]"
!write(*,*)"In steps of istep_npart  :  ",istep_npart
!write(*,*)"In steps of istep_nrot   :  ",istep_nrot
!write(*,*)"In steps of istep_nk     :  ",istep_nk
!
!
!*******************************************************************************
!    Subroutine to run sanity checks on the data structure passed CPU
!
!*******************************************************************************
!
subroutine Sanity_check_cpu(sysQ, hstD)
  use simple_defs
  use greeting_version
  use simple_systemQuery_cpu
  implicit none
  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  !local variables
  !cpu gear
  integer                         :: nCPU_cores
  integer*8                       :: h_TotalMemSize
  integer*8                       :: h_AvailMemSize

  !start of the execution commands

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'  Sanity checks on the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  h_TotalMemSize = sysQ%get_SysTotalMem_size()
  h_AvailMemSize = sysQ%get_SysAvailMem_size()
  write(*,*)'Number of cores on system     (sysQ):',nCPU_cores
  write(*,*)'Number of cores               (hstD):',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()
  write(*,*)'Total Mem on system           (sysQ):',h_TotalMemSize
  write(*,*)'Total Mem on system           (hstD):',hstD%mem_Size
  if(h_TotalMemSize /= hstD%mem_Size)call sysQ%get_warning_dataStruct()
  write(*,*)'Total Available Mem on system (sysQ):',h_AvailMemSize
#if defined (MACOSX)
  write(*,*)'Total Mem on system           (hstD):',hstD%mem_User
  if(h_AvailMemSize /= hstD%mem_User)call sysQ%get_warning_dataStruct()
#elif defined (LINUX)
  write(*,*)'Total Mem on system           (hstD):',hstD%avail_Mem
  if(h_AvailMemSize /= hstD%avail_Mem)call sysQ%get_warning_dataStruct()
#endif
  
  return
end subroutine Sanity_check_cpu

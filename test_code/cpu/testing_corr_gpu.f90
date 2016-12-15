program testing_gen_polar_coords

  use simple_defs
  use simple_timing
  use greeting_version
  use matrixGetter
  use simple_math, only: calc_corr, csq

  implicit none
#define devptr_t integer*8

  integer                      :: nradial1           !< #radial vectors (angle index) (= # of components in each shell)
  integer                      :: nradial2           !< #radial vectors (angle index) (= # of components in each shell)
  integer                      :: klp               !< low-pass frequency limit
  integer                      :: khp               !< high-pass frequency limit 

  integer                      :: rot1, rot2
  integer                      :: ring2
  integer                      :: lda
  integer, parameter           :: n=3

  real(sp)                     :: r, sumasq, sumbsq
  complex(sp), allocatable     :: pft1(:,:)
  complex(sp), allocatable     :: pft2(:,:)

  double precision,allocatable :: re_matZA(:,:)
  double precision,allocatable :: im_matZA(:,:)

  !counters

  integer                      :: i,j,k,in
  integer                      :: i1,i2

  !start of the execution commands
  !start of the greeting message
  call hello_gpu_math()
  call timestamp()
  call start_Alltimers_cpu()

  ring2 = 2
 
  rot1 = 1
  rot2 = 1
  nradial1 = n
  nradial2 = n
  klp = n
  khp = 1

  write(*,*) 
  write(*,'(7(5x,a))') "ring2","rot1","rot2","nradial1","nradial2","klp","khp"
  write(*,'(7(7x,i3))') ring2, rot1, rot2, nradial1, nradial2, klp, khp
  write(*,*) 

  lda = n

  !allocating the complex matrices
  allocate(pft1(n,n))
  allocate(pft2(n,n))

  call get1to9RowMajZspMat_cpu(n,lda,pft1)
  call get1to9ColMajZspMat_cpu(n,lda,pft2)

  write(*,'(34x,a,32x,a)')"(pft1) row major","(pft2) column major"
  do i=1,n
     do j=1,n
        write(*,*)i, j,pft1(i,j), pft2(i,j)
     end do
  end do
  write(*,'(34x,a)')"conjg(pft2)"
  do i=1,n
     do j=1,n
        write(*,*)i, j,conjg(pft2(i,j))
     end do
  end do

  i1     = rot1
  i2     = rot2
  sumasq = 0.
  sumbsq = 0.
  r      = 0.

  write(*,*)"the correlator"
!  do i=1,nradial1/2
  do i=1,nradial1
     do j=khp,klp
        r = r+real(pft1(i1,j)*conjg(pft2(i2,j)))
        sumasq = sumasq+csq(pft1(i1,j))
        sumbsq = sumbsq+csq(pft2(i2,j))
        write(*,*)i, j, r
     end do
     i1 = i1+1
     if( i1 > nradial1 ) i1 = 1
     i2 = i2+1
     if( i2 > nradial2 ) i2 = 1
     write(*,*)i1,i2,i,j,r
  end do
  r = calc_corr(r,sumasq*sumbsq)

  write(*,*)"after loops r= ",r

  deallocate(pft1)  
  deallocate(pft2)

  !shutting down the environment

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_math()



end program testing_gen_polar_coords


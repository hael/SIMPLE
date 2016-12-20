!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 18th of June 2015.
!
! Name:
! testfunction - Various utilities and testfunctions
!
! Description:
! testfunctions provides test fcuntions for test codes and other modules:
!*******************************************************************************
!
module simple_testfunction
  use simple_defs

implicit none

contains
!*******************************************************************************
!    3D atom_like
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  !Generating 3D double complex functional for testign purposes.
  subroutine atom_like_3D_complex(ui,uf,vi,vf,wi,wf,nx,ny,nz,a,b,c,u,v,w,fuvw,n)
    use simple_defs
    implicit none
    integer  :: nx,ny,nz
    real(dp) :: a,b,c,n
    complex(dp) :: ui,uf,vi,vf,wi,wf
    complex(dp) :: u(1:nx),v(1:ny),w(1:nz)
    complex(dp) :: fuvw(nx,ny,*)
    !local variables
    real(dp) :: delta_ur,delta_ui
    real(dp) :: delta_vr,delta_vi
    real(dp) :: delta_wr,delta_wi
    real(dp),dimension(1:nx) :: re_u,im_u
    real(dp),dimension(1:ny) :: re_v,im_v
    real(dp),dimension(1:nz) :: re_w,im_w
    real(dp),dimension(1:nx,1:ny,1:nz) :: re_fuvw,im_fuvw
    !counters
    integer :: i,j,k

    delta_ur = ( real(uf) - real(ui) ) / (nx-1)
    delta_ui = ( imag(uf) - imag(ui) ) / (nx-1)
    do i=1,nx
       re_u(i) = real(ui) + (i-1) * delta_ur
       im_u(i) = imag(ui) + (i-1) * delta_ui
    end do

    delta_vr = ( real(vf) - real(vi) ) / (ny-1)
    delta_vi = ( imag(vf) - imag(vi) ) / (ny-1)
    do j=1,ny
       re_v(j) = real(vi) + (j-1) * delta_vr
       im_v(j) = imag(vi) + (j-1) * delta_vi
    end do

    delta_wr = ( real(wf) - real(wi) ) / (nz-1)
    delta_wi = ( imag(wf) - imag(wi) ) / (nz-1)
    do k=1,nz
       re_w(k) = real(wi) + (k-1) * delta_wr
       im_w(k) = imag(wi) + (k-1) * delta_wi
    end do

    u(1:nx) = cmplx(re_u(1:nx),im_u(1:nx),dp)
    v(1:ny) = cmplx(re_v(1:ny),im_v(1:ny),dp)
    w(1:nz) = cmplx(re_w(1:nz),im_w(1:nz),dp)

    fuvw(1:nx,1:ny,1:nz) = 0.0d0
    
    do i=1,nx
       do j=1,ny
          do k=1,nz
             re_fuvw(i,j,k) = cos(( real(u(i)) + real(v(j)) + real(w(k)) )**b ) / &
                  ( (real(u(i))**a + real(v(j))**a + real(w(k))**a + c )**n)
             im_fuvw(i,j,k) = sin(( imag(u(i)) + imag(v(j)) + imag(w(k)) )**b ) / &
                  ( (imag(u(i))**a + imag(v(j))**a + imag(v(j))**a + c )**n)
          end do
       end do
    end do

    fuvw(1:nx,1:ny,1:nz) = cmplx(re_fuvw(1:nx,1:ny,1:nz),im_fuvw(1:nx,1:ny,1:nz),dp)

    return
  end subroutine atom_like_3D_complex
!*******************************************************************************
!    3D Sequencer
!    Filling and printing array with sequencing/1000.0
!
!*******************************************************************************
!
  !Generating 3D double complex precision gaussian distributed matrix.
  subroutine getCSeq_3D(nx,ny,nz,fuvw)
    use simple_defs
    implicit none
    integer  :: nx,ny,nz
    complex(sp) :: fuvw(nx,ny,*)
    !local variables
    integer  :: ix,iy,iz

    fuvw(:nx,:ny,:nz) = 0.0
    do ix=1,nx
       do iy=1,ny
          do iz=1,nz
             fuvw(ix,iy,iz) = cmplx(ix+iy+iz,sp) / 1000.0
          end do
       end do
    end do
    
    return
  end subroutine getCSeq_3D
!*******************************************************************************
!    3D RandomGaussianDistr
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  !Generating 3D double complex precision gaussian distributed matrix.
  subroutine getZRandomGaussianDistr_3D(nx,ny,nz,fuvw)
    use simple_defs
    implicit none
    integer  :: nx,ny,nz
    complex(dp) :: fuvw(nx,ny,*)
    !local variables
    real(dp),dimension(nx*ny*nz) :: harvestr
    real(dp),dimension(nx,ny,nz) :: r
    real(dp),dimension(nx,ny,nz) :: theta
    integer,dimension(3)         :: shaper

    shaper(1) = nx
    shaper(2) = ny
    shaper(3) = nz

    call random_number (harvestr)
    where (harvestr == 0.0d0) harvestr = 1.0d0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0d0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0d0) harvestr = 1.0d0
    theta = reshape(harvestr, shaper)
    theta = 2.0d0 * 4.0 * atan(1.0d0) * theta
    
    fuvw(1:nx,1:ny,1:nz) = cmplx(r(:,:,:) * cos(theta(:,:,:)), &
                                 r(:,:,:) * sin(theta(:,:,:)), dp)

    return
  end subroutine getZRandomGaussianDistr_3D

  !Generating 3D double complex precision gaussian distributed matrix.
  subroutine getCRandomGaussianDistr_3D(nx,ny,nz,fuvw)
    use simple_defs
    implicit none
    integer  :: nx,ny,nz
    complex(sp) :: fuvw(nx,ny,*)
    !local variables
    real(sp),dimension(nx*ny*nz) :: harvestr
    real(sp),dimension(nx,ny,nz) :: r
    real(sp),dimension(nx,ny,nz) :: theta
    integer,dimension(3)         :: shaper

    shaper(1) = nx
    shaper(2) = ny
    shaper(3) = nz

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fuvw(1:nx,1:ny,1:nz) = cmplx(r(:,:,:) * cos(theta(:,:,:)), &
                                 r(:,:,:) * sin(theta(:,:,:)), sp)

    return
  end subroutine getCRandomGaussianDistr_3D
  
  !Generating 3D double precision gaussian distributed matrix.
  subroutine getDRandomGaussianDistr_3D(nx,ny,nz,fxyz)
    use simple_defs
    implicit none
    integer  :: nx,ny,nz
    real(dp) :: fxyz(nx,ny,*)
    !local variables
    real(dp),dimension(nx*ny*nz) :: harvestr
    real(dp),dimension(nx,ny,nz) :: r
    real(dp),dimension(nx,ny,nz) :: theta
    integer,dimension(3)         :: shaper

    shaper(1) = nx
    shaper(2) = ny
    shaper(3) = nz

    call random_number (harvestr)
    where (harvestr == 0.0d0) harvestr = 1.0d0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0d0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0d0) harvestr = 1.0d0
    theta = reshape(harvestr, shaper)
    theta = 2.0d0 * 4.0 * atan(1.0d0) * theta

    fxyz(1:nx,1:ny,1:nz) = r(:,:,:) * ( cos(theta(:,:,:)) + sin(theta(:,:,:)) )

    return
  end subroutine getDRandomGaussianDistr_3D

  !Generating 3D single precision gaussian distributed matrix.
  subroutine getSRandomGaussianDistr_3D(nx,ny,nz,fxyz)
    use simple_defs
    implicit none
    integer  :: nx,ny,nz
    real(sp) :: fxyz(nx,ny,*)
    !local variables
    real(sp),dimension(nx*ny*nz) :: harvestr
    real(sp),dimension(nx,ny,nz) :: r
    real(sp),dimension(nx,ny,nz) :: theta
    integer,dimension(3)         :: shaper

    shaper(1) = nx
    shaper(2) = ny
    shaper(3) = nz

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fxyz(1:nx,1:ny,1:nz) = r(:,:,:) * ( cos(theta(:,:,:)) + sin(theta(:,:,:)) )

    return
  end subroutine getSRandomGaussianDistr_3D
!*******************************************************************************
!    2D RandomGaussianDistr
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  !Generating 2D single precision gaussian distributed matrix.
  subroutine getSRandomGaussianDistr_2D(nx,ny,fxy)
    use simple_defs
    implicit none
    integer  :: nx,ny
    real(sp) :: fxy(nx,*)
    !local variables
    real(sp),dimension(nx*ny) :: harvestr
    real(sp),dimension(nx,ny) :: r
    real(sp),dimension(nx,ny) :: theta
    integer,dimension(2)      :: shaper

    shaper(1) = nx
    shaper(2) = ny

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fxy(1:nx,1:ny) = r(:,:) * ( cos(theta(:,:)) + sin(theta(:,:)) )

    return
  end subroutine getSRandomGaussianDistr_2D
  !Generating 2D double precision gaussian distributed matrix.
  subroutine getDRandomGaussianDistr_2D(nx,ny,fxy)
    use simple_defs
    implicit none
    integer  :: nx,ny
    real(dp) :: fxy(nx,*)
    !local variables
    real(dp),dimension(nx*ny) :: harvestr
    real(dp),dimension(nx,ny) :: r
    real(dp),dimension(nx,ny) :: theta
    integer,dimension(2)      :: shaper

    shaper(1) = nx
    shaper(2) = ny

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fxy(1:nx,1:ny) = r(:,:) * ( cos(theta(:,:)) + sin(theta(:,:)) )

    return
  end subroutine getDRandomGaussianDistr_2D
  !Generating 2D complex precision gaussian distributed matrix.
  subroutine getCRandomGaussianDistr_2D(nx,ny,fuvw)
    use simple_defs
    implicit none
    integer  :: nx,ny
    complex(sp) :: fuvw(nx,*)
    !local variables
    real(sp),dimension(nx*ny) :: harvestr
    real(sp),dimension(nx,ny) :: r
    real(sp),dimension(nx,ny) :: theta
    integer,dimension(2)      :: shaper

    shaper(1) = nx
    shaper(2) = ny

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fuvw(1:nx,1:ny) = cmplx(r(:,:) * cos(theta(:,:)), &
                            r(:,:) * sin(theta(:,:)), sp)

    return
  end subroutine getCRandomGaussianDistr_2D
  !Generating 2D complex precision gaussian distributed matrix.
  subroutine getZRandomGaussianDistr_2D(nx,ny,fuvw)
    use simple_defs
    implicit none
    integer  :: nx,ny
    complex(dp) :: fuvw(nx,*)
    !local variables
    real(dp),dimension(nx*ny) :: harvestr
    real(dp),dimension(nx,ny) :: r
    real(dp),dimension(nx,ny) :: theta
    integer,dimension(2)      :: shaper

    shaper(1) = nx
    shaper(2) = ny

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fuvw(1:nx,1:ny) = cmplx(r(:,:) * cos(theta(:,:)), &
                            r(:,:) * sin(theta(:,:)), dp)

    return
  end subroutine getZRandomGaussianDistr_2D

!*******************************************************************************
!    2D atom_like
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  !Generating 2D double complex functional for testign purposes.
  subroutine atom_like_2D_complex(ui,uf,vi,vf,nx,ny,a,b,c,u,v,fuv,n)! result(fuv)
    use simple_defs
    implicit none
    integer  :: nx,ny
    real(dp) :: a,b,c,n
    complex(dp) :: ui,uf,vi,vf
    complex(dp) :: u(*),v(*)
    complex(dp) :: fuv(nx,*)
    !local variables
    real(dp) :: delta_ur,delta_ui
    real(dp) :: delta_vr,delta_vi
    real(dp),dimension(nx) :: re_u,im_u
    real(dp),dimension(ny) :: re_v,im_v
    real(dp),dimension(nx,ny) :: re_fuv,im_fuv
    !counters
    integer :: i,j
    !here insert the function code
    !TODO: check the nx and ny indexing inb the construction

    delta_ur = ( real(uf) - real(ui) ) / (nx-1)
    delta_ui = ( imag(uf) - imag(ui) ) / (nx-1)
    do i=1,nx
       re_u(i) = real(ui) + (i-1) * delta_ur
       im_u(i) = imag(ui) + (i-1) * delta_ui
    end do

    delta_vr = ( real(vf) - real(vi) ) / (ny-1)
    delta_vi = ( imag(vf) - imag(vi) ) / (ny-1)
    do j=1,ny
       re_v(j) = real(vi) + (j-1) * delta_vr
       im_v(j) = imag(vi) + (j-1) * delta_vi
    end do

    u(1:nx) = cmplx(re_u(1:nx),im_u(1:nx),dp)
    v(1:ny) = cmplx(re_v(1:ny),im_v(1:ny),dp)

    fuv(1:nx,1:ny) = 0.0d0

    do i=1,nx
       do j=1,ny
          re_fuv(i,j) = cos(real(u(i))**b) / ( (real(u(i))**a + real(v(j))**a + c )**n)
          im_fuv(i,j) = sin(imag(v(j))**b) / ( (imag(u(i))**a + imag(v(j))**a + c )**n)
       end do
    end do

    fuv(1:nx,1:ny) = cmplx(re_fuv(1:nx,1:ny),im_fuv(1:nx,1:ny),dp)

    return
  end subroutine atom_like_2D_complex
  !Generating 2D single complex functional for testing purposes.
  subroutine atom_like_2D_single_complex(ui,uf,vi,vf,nx,ny,a,b,c,u,v,fuv,n)! result(fuv)
    use simple_defs
    implicit none
    integer  :: nx,ny
    real(sp) :: a,b,c,n
    complex(sp) :: ui,uf,vi,vf
    complex(sp) :: u(*),v(*)
    complex(sp) :: fuv(nx,*)
    !local variables
    real(sp) :: delta_ur,delta_ui
    real(sp) :: delta_vr,delta_vi
    real(sp),dimension(nx) :: re_u,im_u
    real(sp),dimension(ny) :: re_v,im_v
    real(sp),dimension(nx,ny) :: re_fuv,im_fuv
    !counters
    integer :: i,j
    !here insert the function code
    !TODO: check the nx and ny indexing inb the construction

    delta_ur = ( real(uf) - real(ui) ) / (nx-1)
    delta_ui = ( imag(uf) - imag(ui) ) / (nx-1)
    do i=1,nx
       re_u(i) = real(ui) + (i-1) * delta_ur
       im_u(i) = imag(ui) + (i-1) * delta_ui
    end do

    delta_vr = ( real(vf) - real(vi) ) / (ny-1)
    delta_vi = ( imag(vf) - imag(vi) ) / (ny-1)
    do j=1,ny
       re_v(j) = real(vi) + (j-1) * delta_vr
       im_v(j) = imag(vi) + (j-1) * delta_vi
    end do

    u(1:nx) = cmplx(re_u(1:nx),im_u(1:nx),dp)
    v(1:ny) = cmplx(re_v(1:ny),im_v(1:ny),dp)

    fuv(1:nx,1:ny) = 0.0

    do i=1,nx
       do j=1,ny
          re_fuv(i,j) = cos(real(u(i))**b) / ( (real(u(i))**a + real(v(j))**a + c )**n)
          im_fuv(i,j) = sin(imag(v(j))**b) / ( (imag(u(i))**a + imag(v(j))**a + c )**n)
       end do
    end do

    fuv(1:nx,1:ny) = cmplx(re_fuv(1:nx,1:ny),im_fuv(1:nx,1:ny),dp)

    return
  end subroutine atom_like_2D_single_complex
  !Generating 2D double functional for testign purposes.
  subroutine atom_like_2D_double(ui,uf,nx,ny,a,b,c,x_in,y_in,fxy,n)! result(fxy)
    use simple_defs
    implicit none
    integer  :: nx,ny
    real(dp) :: a,b,c,n
    complex(dp) :: ui,uf
    real(dp) :: deltax,deltay
    real(dp) :: fxy(nx,*)
    real(dp) :: x_in(*)
    real(dp) :: y_in(*)
    !local variables
    real(dp),dimension(nx) :: x
    real(dp),dimension(ny) :: y
    real(dp),dimension(nx,ny) :: re_fxy
    !counters 
    integer :: i,j

    deltax = ( real(uf) - real(ui) ) / (nx-1)
    do i=1,nx
       x(i) = real(ui) + (i-1) * deltax
    end do

    deltay = ( imag(uf) - imag(ui) ) / (ny-1)
    do j=1,ny
       y(j) = imag(ui) + (j-1) * deltay
    end do

    fxy(1:nx,1:ny) = 0.0d0

    do i = 1,nx
       do j = 1,ny
          re_fxy(i,j) = cos(x(i)**b)     / (  x(i)**a + y(j)**a + c )**n + &
               cos((x(i)-5)**b) / ( (x(i)-5)**a + (y(j)-5)**a + c )**n + &
               cos((x(i)+5)**b) / ( (x(i)+5)**a + (y(j)+5)**a + c )**n;  
       end do
    end do

    x_in(1:nx) = x(1:nx) 
    y_in(1:nx) = y(1:ny) 
    fxy(1:nx,1:ny) = re_fxy(1:nx,1:ny)

    return
  end subroutine atom_like_2D_double
  !Generating 2D double functional for testign purposes.
  subroutine atom_like_2D_single(ui,uf,nx,ny,a,b,c,x_in,y_in,fxy,n)! result(fxy)
    use simple_defs
    implicit none
    integer  :: nx,ny
    real(sp) :: a,b,c,n
    complex(sp) :: ui,uf
    real(sp) :: deltax,deltay
    real(sp) :: fxy(nx,*)
    real(sp) :: x_in(*)
    real(sp) :: y_in(*)
    !local variables
    real(sp),dimension(nx) :: x
    real(sp),dimension(ny) :: y
    real(sp),dimension(nx,ny) :: re_fxy
    !counters 
    integer :: i,j

    deltax = ( real(uf) - real(ui) ) / (nx-1)
    do i=1,nx
       x(i) = real(ui) + (i-1) * deltax
    end do

    deltay = ( imag(uf) - imag(ui) ) / (ny-1)
    do j=1,ny
       y(j) = imag(ui) + (j-1) * deltay
    end do

    fxy(1:nx,1:ny) = 0.0d0

    do i = 1,nx
       do j = 1,ny
          re_fxy(i,j) = cos(x(i)**b) / (  x(i)**a + y(j)**a + c )**n
          !+ &
          !cos((x(i)-5)**b) / ( (x(i)-5)**a + (y(j)-5)**a + c )**n + &
          !cos((x(i)+5)**b) / ( (x(i)+5)**a + (y(j)+5)**a + c )**n;  
       end do
    end do

    x_in(1:nx) = x(1:nx) 
    y_in(1:nx) = y(1:ny) 
    fxy(1:nx,1:ny) = re_fxy(1:nx,1:ny)

    return
  end subroutine atom_like_2D_single

  !Generating 2D double functional ptwise for testign purposes.
  function atom_like_2D_single_ptwise(a,b,c,x_in,y_in,n) result(fxy_ptwise)
    use simple_defs
    implicit none
    real(sp) :: a,b,c,n
    real(sp) :: deltax,deltay
    real(sp) :: fxy_ptwise
    real(sp) :: x_in
    real(sp) :: y_in

    fxy_ptwise = 0.0d0

    fxy_ptwise = cos(x_in**b) / (  x_in**a + y_in**a + c )**n
    ! + &
    !cos((x_in-5)**b) / ( (x_in-5)**a + (y_in-5)**a + c )**n + &
    !cos((x_in+5)**b) / ( (x_in+5)**a + (y_in+5)**a + c )**n

    return
  end function atom_like_2D_single_ptwise

  !Generating 2D double functional ptwise for testign purposes.
  function atom_like_2D_single_ptwise_symbolicDiff(a,b,c,n,x_D,D) result(grad_fxy_ptwise)
    real(sp) :: a,b,c,n
    integer, intent(in) :: D
    real, intent(inout) :: x_D(D)
    real(sp) :: grad_fxy_ptwise(D)

    grad_fxy_ptwise(1) = (-sin(x_D(1)**b)*b*x_D(1)**(b-1.0) - &
         cos(x_D(1)**b)*n*a*x_D(1)**(a-1.0) * (x_D(1)**a + x_D(2)**a + c)**(-1.0) ) / &
         (  x_D(1)**a + x_D(2)**a + c )**n

    grad_fxy_ptwise(2) = atom_like_2D_single_ptwise(a,b,c,x_D(1),x_D(2),n) * &
         ( (-n / ( x_D(1)**a + x_D(2)**a + c )**n ) * a * x_D(2) * (a-1.0) )

  end function atom_like_2D_single_ptwise_symbolicDiff
!*******************************************************************************
!    2D RandomGaussianDistr
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  !Generating 1D single precision gaussian distributed matrix.
  subroutine getSRandomGaussianDistr_1D(nx,fx)
    use simple_defs
    implicit none
    integer  :: nx
    real(sp) :: fx(*)
    !local variables
    real(sp),dimension(nx) :: harvestr
    real(sp),dimension(nx) :: r
    real(sp),dimension(nx) :: theta
    integer,dimension(1)   :: shaper

    shaper(1) = nx

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fx(1:nx) = r(:) * ( cos(theta(:)) + sin(theta(:)) )

    return
  end subroutine getSRandomGaussianDistr_1D
  !Generating 1D double precision gaussian distributed matrix.
  subroutine getDRandomGaussianDistr_1D(nx,fx)
    use simple_defs
    implicit none
    integer  :: nx
    real(dp) :: fx(*)
    !local variables
    real(dp),dimension(nx) :: harvestr
    real(dp),dimension(nx) :: r
    real(dp),dimension(nx) :: theta
    integer,dimension(1)   :: shaper

    shaper(1) = nx

    call random_number (harvestr)
    where (harvestr == 0.0d0) harvestr = 1.0d0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0d0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0d0) harvestr = 1.0d0
    theta = reshape(harvestr, shaper)
    theta = 2.0d0 * 4.0d0 * atan(1.0d0) * theta

    fx(1:nx) = r(:) * ( cos(theta(:)) + sin(theta(:)) )

    return
  end subroutine getDRandomGaussianDistr_1D
  !Generating 1D double complex precision gaussian distributed matrix.
  subroutine getCRandomGaussianDistr_1D(nx,fu)
    use simple_defs
    implicit none
    integer  :: nx
    complex(sp) :: fu(*)
    !local variables
    real(sp),dimension(nx) :: harvestr
    real(sp),dimension(nx) :: r
    real(sp),dimension(nx) :: theta
    integer,dimension(1)   :: shaper

    shaper(1) = nx

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fu(1:nx) = cmplx(r(:) * cos(theta(:)), &
                     r(:) * sin(theta(:)), sp)

    return
  end subroutine getCRandomGaussianDistr_1D
  !Generating 1D double complex precision gaussian distributed matrix.
  subroutine getZRandomGaussianDistr_1D(nx,fu)
    use simple_defs
    implicit none
    integer  :: nx
    complex(dp) :: fu(*)
    !local variables
    real(dp),dimension(nx) :: harvestr
    real(dp),dimension(nx) :: r
    real(dp),dimension(nx) :: theta
    integer,dimension(1)   :: shaper

    shaper(1) = nx

    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    r = reshape(harvestr,shaper)
    r = sqrt(-2.0*log(r) )
    call random_number (harvestr)
    where (harvestr == 0.0) harvestr = 1.0
    theta = reshape(harvestr, shaper)
    theta = 2.0 * 4.0 * atan(1.0) * theta

    fu(1:nx) = cmplx(r(:) * cos(theta(:)), &
                     r(:) * sin(theta(:)), dp)

    return
  end subroutine getZRandomGaussianDistr_1D
!*******************************************************************************
!    1D atom_like
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  !Generating 1D functional for testign purposes.
  function atom_like_complex(ui,uf,nstep,a,b,c,u,n) result(fu)
    use simple_defs
    integer  :: nstep
    real(dp) :: a,b,c,n
    complex(dp) :: ui,uf
    complex(dp),dimension(nstep) :: u,fu
    !local variables
    real(dp) :: re_delta,im_delta
    real(dp) :: re_u,im_u
    !counters
    integer :: i

    re_delta = ( real(uf) - real(ui) ) / (nstep-1)
    im_delta = ( imag(uf) - imag(ui) ) / (nstep-1)
    do i=1,nstep
       re_u = real(ui) + (i-1) * re_delta
       im_u = imag(ui) + (i-1) * im_delta
       u(i) = cmplx(re_u , im_u, dp)
    end do

    fu(:) = cmplx( &
         cos(real(u(:))**b) / ( (real(u(:))**a) + c )**n, &
         cos(imag(u(:))**b) / ( (imag(u(:))**a) + c )**n, dp)

  end function atom_like_complex

  !Generating 1D functional for testign purposes.
  function atom_like_real(xi,xf,nx,a,b,c,x,n) result(fx)
    use simple_defs
    integer  :: nx
    real(sp) :: xi,xf
    real(sp) :: a,b,c,n
    real(sp) :: delta
    real(sp),dimension(nx) :: x,fx
    real(sp) :: x_i
    !counters
    integer :: i

    delta = ( xf - xi ) / (nx-1)
    x_i = xi
    do i=1,nx
       x(i) = xi + (i-1) * delta
    end do

    fx(:) = cos(x(:)**b) / ( (x(:)**a) + c )**n

  end function atom_like_real
!*******************************************************************************
!    1D atom_like element wise
!    Filling and printing array with random numbers using test functions
!
!*******************************************************************************
!
  function atom_like_real_ptwise(a,b,c,x,n) result(fx_ptwise)
    use simple_defs
    real(sp) :: a,b,c,n
    real(sp) :: x,fx_ptwise

    fx_ptwise = cos(x**b) / ( (x**a) + c )**n

  end function atom_like_real_ptwise

  function atom_like_real_ptwise_symbolicDiff(x) result(fx_ptwise)
    real(sp) :: x,fx_ptwise

    fx_ptwise = ( -2.0*x/(x**2.0+1.0)**2.0 ) * ( cos(x**2.0)+sin(x**2.0)*(x**2.0 + 1) )

  end function atom_like_real_ptwise_symbolicDiff

end module simple_testfunction

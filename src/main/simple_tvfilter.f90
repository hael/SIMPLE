! Filter based on total variation in real space
module simple_tvfilter
include 'simple_lib.f08'
use simple_image, only: image
implicit none
private
public :: tvfilter, test_tvfilter, test_tvfilter_3d
#include "simple_local_flags.inc"

type :: tvfilter
    logical               :: existence
    type(image)           :: interpolate_coeffs
    real, pointer         :: interpolate_coeffs_rmat(:,:,:) => null()
    integer               :: img_dims(2), img_dims_3d(3)
    type(image)           :: r_img, b_img
    integer               :: ldim(2)
contains
    procedure :: new               => new_tvfilter
    procedure :: apply_filter     ! 2D by default
    procedure :: apply_filter_3d
    procedure :: prepare_interpolation
    procedure :: interp_val
    procedure :: kill              => kill_tvfilter
    procedure, private :: fill_r
    procedure, private :: fill_b
    procedure, private :: fill_r_3d
    procedure, private :: fill_b_3d
end type tvfilter

contains

    subroutine new_tvfilter( self )
        class(tvfilter), intent(inout) :: self
        self%existence   = .true.
    end subroutine new_tvfilter

    subroutine apply_filter( self, img, lambda, idx )
        class(tvfilter),   intent(inout) :: self
        class(image),      intent(inout) :: img
        real,              intent(in)    :: lambda ! >0.; 0.1 is a starting point
        integer, optional, intent(in)    :: idx
        integer :: idx_here
        integer :: img_ldim(3), rb_ldim(3)
        real    :: img_smpd
        logical :: img_ft_prev
        complex(kind=c_float_complex), pointer :: cmat_b(:,:,:), cmat_r(:,:,:), cmat_img(:,:,:)
        logical :: do_alloc
        integer :: dims1
        if (.not. present(idx)) then
            idx_here = 1
        else
            idx_here = idx
        end if
        img_ldim = img%get_ldim()
        self%img_dims(1:2) = img_ldim(1:2)
        if (idx_here > img_ldim(3)) then
            THROW_HARD('tvfilter::apply_filter : idx greater than stack size')
        end if
        img_smpd = img%get_smpd()
        do_alloc = .true.
        if (self%r_img%exists()) then
            rb_ldim  = self%r_img%get_ldim()
            if (all(img_ldim(1:2) == rb_ldim(1:2))) then
                do_alloc = .false.
            end if
        end if
        if (do_alloc) then
            rb_ldim(1:2) = img_ldim(1:2)
            rb_ldim(3)   = 1
            call self%r_img%new(rb_ldim, img_smpd)
            call self%b_img%new(rb_ldim, img_smpd)
            call self%fill_b()
            call self%fill_r()
            call self%r_img%fft_noshift()
            call self%b_img%fft_noshift()
        end if
        img_ft_prev = img%is_ft()
        if (.not. img_ft_prev) call img%fft()
        call self%b_img%get_cmat_ptr(cmat_b)
        call self%r_img%get_cmat_ptr(cmat_r)
        call img%get_cmat_ptr(cmat_img)
        dims1 = int(img_ldim(1)/2)+1
        cmat_img(1:dims1,:,idx_here) = cmat_img(1:dims1,:,idx_here) * (real(cmat_b(1:dims1,:,1))**2 + aimag(cmat_b(1:dims1,:,1))**2) / &
            (real(cmat_b(1:dims1,:,1))**2 + aimag(cmat_b(1:dims1,:,1))**2 + lambda * cmat_r(1:dims1,:,1))
        if (.not. img_ft_prev) call img%ifft()
    end subroutine apply_filter

    subroutine apply_filter_3d( self, img, lambda )
        class(tvfilter),   intent(inout) :: self
        class(image),      intent(inout) :: img
        real,              intent(in)    :: lambda ! >0.; 0.1 is a starting point
        integer :: img_ldim(3), rb_ldim(3)
        real    :: img_smpd
        logical :: img_ft_prev
        complex(kind=c_float_complex), pointer :: cmat_b(:,:,:), cmat_r(:,:,:), cmat_img(:,:,:)
        logical :: do_alloc
        integer :: dims1
        img_ldim = img%get_ldim()
        self%img_dims_3d(1:3) = img_ldim(1:3)
        img_smpd = img%get_smpd()
        do_alloc = .true.
        if (self%r_img%exists()) then
            rb_ldim  = self%r_img%get_ldim()
            if (all(img_ldim(1:3) == rb_ldim(1:3))) then
                do_alloc = .false.
            end if
        end if
        if (do_alloc) then
            rb_ldim(1:3) = img_ldim(1:3)
            call self%r_img%new(rb_ldim, img_smpd)
            call self%b_img%new(rb_ldim, img_smpd)
            call self%fill_b_3d()
            call self%fill_r_3d()
            call self%r_img%fft_noshift()
            call self%b_img%fft_noshift()
        end if
        img_ft_prev = img%is_ft()
        if (.not. img_ft_prev) call img%fft()
        call self%b_img%get_cmat_ptr(cmat_b)
        call self%r_img%get_cmat_ptr(cmat_r)
        call img%get_cmat_ptr(cmat_img)
        dims1 = int(img_ldim(1)/2)+1
        cmat_img(1:dims1,:,:) = cmat_img(1:dims1,:,:) * (real(cmat_b(1:dims1,:,:))**2 + aimag(cmat_b(1:dims1,:,:))**2) / &
            (real(cmat_b(1:dims1,:,:))**2 + aimag(cmat_b(1:dims1,:,:))**2 + lambda * cmat_r(1:dims1,:,:))
        if (.not. img_ft_prev) call img%ifft()
    end subroutine apply_filter_3d

    subroutine prepare_interpolation( self, img, lambda, idx )
        class(tvfilter),   intent(inout) :: self
        class(image),      intent(inout) :: img
        real,              intent(in)    :: lambda ! >0.; 0.1 is a starting point
        integer, optional, intent(in)    :: idx
        integer :: idx_here
        integer :: img_ldim(3), rb_ldim(3)
        real    :: img_smpd
        logical :: img_ft_prev
        complex(kind=c_float_complex), pointer :: cmat_b(:,:,:), cmat_r(:,:,:), cmat_img(:,:,:)
        logical :: do_alloc
        integer :: dims1
        call self%interpolate_coeffs%kill
        self%interpolate_coeffs = img
        if (.not. present(idx)) then
            idx_here = 1
        else
            idx_here = idx
        end if
        img_ldim = img%get_ldim()
        self%img_dims(1:2) = img_ldim(1:2)
        if (idx_here > img_ldim(3)) then
            THROW_HARD('tvfilter::apply_filter : idx greater than stack size')
        end if
        self%ldim(1:2) = img_ldim(1:2)
        img_smpd = img%get_smpd()
        do_alloc = .true.
        if (self%r_img%exists()) then
            rb_ldim  = self%r_img%get_ldim()
            if (all(img_ldim(1:2) == rb_ldim(1:2))) then
                do_alloc = .false.
            end if
        end if
        if (do_alloc) then
            rb_ldim(1:2) = img_ldim(1:2)
            rb_ldim(3)   = 1
            call self%r_img%new(rb_ldim, img_smpd)
            call self%b_img%new(rb_ldim, img_smpd)
            call self%fill_b()
            call self%fill_r()
            call self%r_img%fft_noshift()
            call self%b_img%fft_noshift()
        end if
        img_ft_prev = self%interpolate_coeffs%is_ft()
        if (.not. img_ft_prev) call self%interpolate_coeffs%fft()
        call self%b_img%get_cmat_ptr(cmat_b)
        call self%r_img%get_cmat_ptr(cmat_r)
        call self%interpolate_coeffs%get_cmat_ptr(cmat_img)
        dims1 = int(img_ldim(1)/2)+1
        cmat_img(1:dims1,:,idx_here) = cmat_img(1:dims1,:,idx_here) * conjg(cmat_b(1:dims1,:,1)) / &
            (real(cmat_b(1:dims1,:,1))**2 + aimag(cmat_b(1:dims1,:,1))**2 + lambda * cmat_r(1:dims1,:,1))
        call self%interpolate_coeffs%ifft()
        call self%interpolate_coeffs%get_rmat_ptr(self%interpolate_coeffs_rmat)
    end subroutine prepare_interpolation

    function interp_val(self, x, y) result( val )
        class(tvfilter), intent(in) :: self
        real,            intent(in) :: x, y
        real :: val
        integer :: xi(5), yi(5) ! 3=floor'ed, 1=-2, 2=-1, 4=+1, 5=+2
        integer :: xi_h(5), yi_h(5)  ! TODO: reduce to 4x4 points
        integer :: i, j
        real    :: dx, dy, b3x, b3y
        xi(3) = floor(x)
        yi(3) = floor(y)
        xi(1) = xi(3) - 2
        yi(1) = yi(3) - 2
        xi(2) = xi(3) - 1
        yi(2) = yi(3) - 1
        xi(4) = xi(3) + 1
        yi(4) = yi(3) + 1
        xi(5) = xi(3) + 2
        yi(5) = yi(3) + 2
        xi_h = xi
        yi_h = yi
        do i = 1, 5
            xi_h(i) = max(xi(i), 1)
            yi_h(i) = max(yi(i), 1)
            xi_h(i) = min(xi_h(i), self%ldim(1))
            yi_h(i) = min(yi_h(i), self%ldim(2))
        end do
        val = 0.
        do i = 1, 5
            dx  = x - real(xi_h(i))
            b3x = b3(dx)
            if (b3x < EPSILON(b3x)) then
                cycle
            end if
            do j = 1, 5
                dy = y - real(yi_h(j))
                b3y = b3(dy)
                val = val + b3x * b3y * self%interpolate_coeffs_rmat(xi_h(i), yi_h(i), 1)
            end do
        end do
    end function interp_val

    subroutine fill_b(self)
        class(tvfilter), intent(inout) :: self
        integer :: nonzero_x(3), nonzero_y(3)         ! indices of non-zero entries in x- and y-component
        real    :: nonzero_pt_x(3),  nonzero_pt_y(3)  ! points of non-zero entries
        real    :: nonzero_val_x(3), nonzero_val_y(3) ! values of non-zero entries
        integer :: i, j, x, y
        real, pointer :: b(:,:,:)
        call self%b_img%get_rmat_ptr(b)
        nonzero_x   (1) = 1
        nonzero_pt_x(1) = 0.
        nonzero_x   (2) = 2
        nonzero_pt_x(2) = 1.
        nonzero_x   (3) = self%img_dims(1)
        nonzero_pt_x(3) = -1.
        nonzero_y   (1) = 1
        nonzero_pt_y(1) = 0.
        nonzero_y   (2) = 2
        nonzero_pt_y(2) = 1.
        nonzero_y   (3) = self%img_dims(2)
        nonzero_pt_y(3) = -1.
        do i = 1,size(nonzero_x)
            nonzero_val_x(i) = b3(nonzero_pt_x(i))
        end do
        do i = 1,size(nonzero_y)
            nonzero_val_y(i) = b3(nonzero_pt_y(i))
        end do
        b = 0.
        do j = 1,size(nonzero_y)
            y = nonzero_y(j)
            do i = 1,size(nonzero_x)
                x = nonzero_x(i)
                b(x,y,1) = nonzero_val_x(i) * nonzero_val_y(j)
            end do
        end do
    end subroutine fill_b

    subroutine fill_b_3d(self)
        class(tvfilter), intent(inout) :: self
        integer :: nonzero_x(3),     nonzero_y(3),     nonzero_z(3)     ! indices of non-zero entries in x-, y- and z-component
        real    :: nonzero_pt_x(3),  nonzero_pt_y(3),  nonzero_pt_z(3)  ! points of non-zero entries
        real    :: nonzero_val_x(3), nonzero_val_y(3), nonzero_val_z(3) ! values of non-zero entries
        integer :: i, j, k, x, y, z
        real, pointer :: b(:,:,:)
        call self%b_img%get_rmat_ptr(b)
        nonzero_x   (1) = 1
        nonzero_pt_x(1) = 0.
        nonzero_x   (2) = 2
        nonzero_pt_x(2) = 1.
        nonzero_x   (3) = self%img_dims_3d(1)
        nonzero_pt_x(3) = -1.
        nonzero_y   (1) = 1
        nonzero_pt_y(1) = 0.
        nonzero_y   (2) = 2
        nonzero_pt_y(2) = 1.
        nonzero_y   (3) = self%img_dims_3d(2)
        nonzero_pt_y(3) = -1.
        nonzero_z   (1) = 1
        nonzero_pt_z(1) = 0.
        nonzero_z   (2) = 2
        nonzero_pt_z(2) = 1.
        nonzero_z   (3) = self%img_dims_3d(3)
        nonzero_pt_z(3) = -1.
        do i = 1,size(nonzero_x)
            nonzero_val_x(i) = b3(nonzero_pt_x(i))
        end do
        do i = 1,size(nonzero_y)
            nonzero_val_y(i) = b3(nonzero_pt_y(i))
        end do
        do i = 1,size(nonzero_z)
            nonzero_val_z(i) = b3(nonzero_pt_z(i))
        end do
        b = 0.
        do k = 1,size(nonzero_z)
            z = nonzero_z(k)
            do j = 1,size(nonzero_y)
                y = nonzero_y(j)
                do i = 1,size(nonzero_x)
                    x = nonzero_x(i)
                    b(x,y,z) = nonzero_val_x(i) * nonzero_val_y(j) * nonzero_val_z(k)
                end do
            end do
        end do
    end subroutine fill_b_3d

    pure elemental function b3(x) result(y)
        real, intent(in) :: x
        real :: y
        y = b3_help(x + 1.5)
    contains
        pure elemental function b3_help(xx) result(yy)
            real, intent(in) :: xx
            real :: yy
            if (xx < 0.) then
                yy = 0.
                return
            end if
            if (xx < 1.) then
                yy = 0.5*xx**2
                return
            end if
            if (xx < 2.) then
                yy = -xx**2+3.*xx-1.5
                return
            end if
            if (xx < 3.) then
                yy = 0.5*xx**2-3.*xx+4.5
                return
            end if
            yy = 0.
        end function b3_help
    end function b3

    subroutine fill_r(self)
        class(tvfilter), intent(inout) :: self
        integer :: nonzero_x(5), nonzero_y(5)   ! indices of non-zero entries in x- and y-component
        real    :: nonzero_val_x_a0(5), nonzero_val_x_a2(5) ! values of non-zero entries
        real    :: nonzero_val_y_a0(5), nonzero_val_y_a2(5) ! values of non-zero entries
        integer :: i, j, x, y
        real, pointer :: r(:,:,:)
        call self%r_img%get_rmat_ptr(r)
        nonzero_x(1) = 1
        nonzero_x(2) = 2
        nonzero_x(3) = 3
        nonzero_x(4) = self%img_dims(1)-1
        nonzero_x(5) = self%img_dims(1)
        nonzero_y(1) = 1
        nonzero_y(2) = 2
        nonzero_y(3) = 3
        nonzero_y(4) = self%img_dims(2)-1
        nonzero_y(5) = self%img_dims(2)
        do i = 1, size(nonzero_x)
            x = nonzero_x(i)
            nonzero_val_x_a0(i) = a0xy(x, self%img_dims(1))
            nonzero_val_x_a2(i) = a2xy(x, self%img_dims(1))
        end do
        do i = 1, size(nonzero_y)
            y = nonzero_y(i)
            nonzero_val_y_a0(i) = a0xy(y, self%img_dims(2))
            nonzero_val_y_a2(i) = a2xy(y, self%img_dims(2))
        end do
        r = 0.
        do j = 1, size(nonzero_x)
            y = nonzero_y(j)
            do i = 1, size(nonzero_x)
                x = nonzero_x(i)
                r(x,y,1) = nonzero_val_x_a0(i) * nonzero_val_y_a2(j) + nonzero_val_x_a2(i) * nonzero_val_y_a0(j)
            end do
        end do
    end subroutine fill_r

    pure function a0xy(i, this_ldim) result(y)
        integer, intent(in) :: i, this_ldim
        real :: y
        if (i == 1) then
            y = 11. / 20.
            return
        end if
        if (i == 2) then
            y = 13. / 60.
            return
        end if
        if (i == 3) then
            y = 1. / 120.
            return
        end if
        if (i == this_ldim-1) then
            y = 1. / 120.
            return
        end if
        if (i == this_ldim) then
            y = 13. / 60.
            return
        end if
        y = 0.
    end function a0xy

    pure function a2xy(i, this_ldim) result(y)
        integer, intent(in) :: i, this_ldim
        real :: y
        if (i == 1) then
            y = 1.
            return
        end if
        if (i == 2) then
            y = -1. / 3.
            return
        end if
        if (i == 3) then
            y = -1. / 6.
            return
        end if
        if (i == this_ldim-1) then
            y = -1. / 6.
            return
        end if
        if (i == this_ldim) then
            y = -1. / 3.
            return
        end if
        y = 0.
    end function a2xy

    subroutine fill_r_3d(self)
        class(tvfilter), intent(inout) :: self
        integer :: nonzero_x(5), nonzero_y(5), nonzero_z(5)   ! indices of non-zero entries in x-, y- and z-component
        real    :: nonzero_val_x_a0(5), nonzero_val_x_a2(5) ! values of non-zero entries
        real    :: nonzero_val_y_a0(5), nonzero_val_y_a2(5) ! values of non-zero entries
        real    :: nonzero_val_z_a0(5), nonzero_val_z_a2(5) ! values of non-zero entries
        integer :: i, j, k, x, y, z
        real, pointer :: r(:,:,:)
        call self%r_img%get_rmat_ptr(r)
        nonzero_x(1) = 1
        nonzero_x(2) = 2
        nonzero_x(3) = 3
        nonzero_x(4) = self%img_dims_3d(1)-1
        nonzero_x(5) = self%img_dims_3d(1)
        nonzero_y(1) = 1
        nonzero_y(2) = 2
        nonzero_y(3) = 3
        nonzero_y(4) = self%img_dims_3d(2)-1
        nonzero_y(5) = self%img_dims_3d(2)
        nonzero_z(1) = 1
        nonzero_z(2) = 2
        nonzero_z(3) = 3
        nonzero_z(4) = self%img_dims_3d(3)-1
        nonzero_z(5) = self%img_dims_3d(3)
        do i = 1, size(nonzero_x)
            x = nonzero_x(i)
            nonzero_val_x_a0(i) = a0xy(x, self%img_dims_3d(1))
            nonzero_val_x_a2(i) = a2xy(x, self%img_dims_3d(1))
        end do
        do i = 1, size(nonzero_y)
            y = nonzero_y(i)
            nonzero_val_y_a0(i) = a0xy(y, self%img_dims_3d(2))
            nonzero_val_y_a2(i) = a2xy(y, self%img_dims_3d(2))
        end do
        do i = 1, size(nonzero_y)
            z = nonzero_z(i)
            nonzero_val_z_a0(i) = a0xy(z, self%img_dims_3d(3))
            nonzero_val_z_a2(i) = a2xy(z, self%img_dims_3d(3))
        end do
        r = 0.
        do k = 1, size(nonzero_z)
            z = nonzero_z(k)
            do j = 1, size(nonzero_y)
                y = nonzero_y(j)
                do i = 1, size(nonzero_x)
                    x = nonzero_x(i)
                    r(x,y,z) = nonzero_val_x_a0(i) * nonzero_val_y_a2(j) * nonzero_val_z_a2(k) + &
                        nonzero_val_x_a2(i) * nonzero_val_y_a0(j) * nonzero_val_z_a0(k)
                end do
            end do
        end do
    end subroutine fill_r_3d


    !>  TEST ROUTINE
    subroutine test_tvfilter( ldim, smpd, lambda )
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, lambda
        type(tvfilter) :: tv
        type(image)    :: img_ref, img_filt, img_tmp
        real           :: cc
        integer        :: i, j, in, jn
        write(logfhandle,*)'>>> REPRODUCIBILITY'
        call tv%new()
        call img_ref%ring(ldim, smpd, real(minval(ldim(1:2)))/4., real(minval(ldim(1:2)))/6.)
        call img_ref%add_gauran(1.)
        call img_ref%shift([32.,48.,0.])
        call img_tmp%ring(ldim, smpd, real(minval(ldim(1:2)))/7., real(minval(ldim(1:2)))/8.)
        call img_tmp%mul(4.)
        call img_tmp%add_gauran(0.5)
        call img_ref%add(img_tmp)
        call img_ref%add_gauran(0.5)
        call img_tmp%kill
        call img_ref%write('test_tvfilter.mrc')
        img_filt = img_ref
        call tv%apply_filter(img_filt, lambda)
        call img_filt%write('test_tvfilter_filt.mrc')
        in = 10
        jn = 10
        do j=1,jn
            call tv%new()
            call progress(j,jn)
            do i = 1,in
                img_tmp = img_ref
                call tv%apply_filter(img_tmp, lambda)
                cc = img_filt%real_corr(img_tmp)
                if( cc < 0.99 )then
                    write(logfhandle,*)'*** Fail at trial: ',i,' - ',j
                    call img_tmp%write('test_tvfilter_filt_'//int2str(i)//'_'//int2str(j)//'.mrc')
                endif
            enddo
            call tv%kill
            call img_tmp%kill
        enddo
    end subroutine test_tvfilter

    subroutine test_tvfilter_3d( ldim, smpd, lambda )
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, lambda
        type(tvfilter) :: tv
        type(image)    :: img_ref, img_filt, img_tmp
        real           :: cc
        integer        :: i, j, in, jn
        write(logfhandle,*)'>>> REPRODUCIBILITY'
        call tv%new()
        call img_ref%ring(ldim, smpd, real(minval(ldim(1:2)))/4., real(minval(ldim(1:2)))/6.)
        call img_ref%add_gauran(1.)
        call img_ref%shift([32.,48.,10.])
        call img_tmp%ring(ldim, smpd, real(minval(ldim(1:2)))/7., real(minval(ldim(1:2)))/8.)
        call img_tmp%mul(4.)
        call img_tmp%add_gauran(0.5)
        call img_ref%add(img_tmp)
        call img_ref%add_gauran(0.5)
        call img_tmp%kill
        call img_ref%write('test_tvfilter.mrc')
        img_filt = img_ref
        call tv%apply_filter_3d(img_filt, lambda)
        call img_filt%write('test_tvfilter_filt_3d.mrc')
        in = 10
        jn = 10
        do j=1,jn
            call tv%new()
            call progress(j,jn)
            do i = 1,in
                img_tmp = img_ref
                call tv%apply_filter_3d(img_tmp, lambda)
                cc = img_filt%real_corr(img_tmp)
                if( cc < 0.99 )then
                    write(logfhandle,*)'*** Fail at trial: ',i,' - ',j
                    call img_tmp%write('test_tvfilter_filt_3d_'//int2str(i)//'_'//int2str(j)//'.mrc')
                endif
            enddo
            call tv%kill
            call img_tmp%kill
        enddo
    end subroutine test_tvfilter_3d

    !>  DESTRUCTOR
    subroutine kill_tvfilter( self )
        class(tvfilter), intent(inout) :: self
        call self%r_img%kill()
        call self%b_img%kill()
        self%interpolate_coeffs_rmat => null()
        call self%interpolate_coeffs%kill
        self%existence = .false.
    end subroutine kill_tvfilter


end module simple_tvfilter

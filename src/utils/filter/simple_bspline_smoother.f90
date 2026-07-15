!@descr: quadratic B-spline Laplacian smoother (Tikhonov regularization, Fourier-domain solve).
! "Generalized smoothing splines and the optimal discretization of the Wiener filter."
! Unser & Blu, IEEE Trans. Signal Processing 53(6), 2146–2159, 2005.
module simple_bspline_smoother
use simple_core_module_api
use simple_image, only: image
implicit none
private
public :: bspline_smoother, test_bspline_smoother, test_bspline_smoother_3d
#include "simple_local_flags.inc"

type :: bspline_smoother
    type(image)   :: r_img, b_img
    type(image)   :: interpolate_coeffs
    real, pointer :: interpolate_coeffs_rmat(:,:,:) => null()
    integer       :: img_dims(2), img_dims_3d(3), ldim(2)
    logical       :: existence
contains
    procedure, private :: new_bspline_smoother
    procedure, private :: new_bspline_smoother_img
    generic            :: new => new_bspline_smoother, new_bspline_smoother_img
    procedure          :: smooth
    procedure          :: smooth_3d
    procedure          :: kill => kill_bspline_smoother
    procedure, private :: fill_r
    procedure, private :: fill_b
    procedure, private :: fill_r_3d
    procedure, private :: fill_b_3d
end type bspline_smoother

contains

    subroutine new_bspline_smoother( self )
        class(bspline_smoother), intent(inout) :: self
        self%existence   = .true.
    end subroutine new_bspline_smoother
    
    subroutine new_bspline_smoother_img( self, img )
        class(bspline_smoother), intent(inout) :: self
        class(image),    intent(inout) :: img
        integer :: img_ldim(3), rb_ldim(3)
        real    :: img_smpd
        self%existence   = .true.
        img_ldim = img%get_ldim()
        if ( img_ldim(3) /= 1 )then
            write(logfhandle,*) 'ldim in bspline_smoother smooth: ', img_ldim(1), img_ldim(2), img_ldim(3)
            THROW_HARD('only for 2D images; bspline_smoother::smooth')
        endif
        self%img_dims(1:2) = img_ldim(1:2)
        img_smpd = img%get_smpd()
        rb_ldim(1:2) = img_ldim(1:2)
        rb_ldim(3)   = 1
        call self%r_img%new(rb_ldim, img_smpd)
        call self%b_img%new(rb_ldim, img_smpd)
        call self%fill_b()
        call self%fill_r()
        call self%r_img%fft_noshift()
        call self%b_img%fft_noshift()
    end subroutine new_bspline_smoother_img

    subroutine smooth( self, img, lambda )
        class(bspline_smoother),   intent(inout) :: self
        class(image),      intent(inout) :: img
        real,              intent(in)    :: lambda ! >0.; 0.1 is a starting point
        integer :: img_ldim(3), rb_ldim(3)
        real    :: img_smpd, lambda_here
        logical :: img_ft_prev
        logical :: do_alloc
        img_ldim = img%get_ldim()
        if ( img_ldim(3) /= 1 )then
            write(logfhandle,*) 'ldim in bspline_smoother smooth: ', img_ldim(1), img_ldim(2), img_ldim(3)
            THROW_HARD('only for 2D images; bspline_smoother::smooth')
        endif
        self%img_dims(1:2) = img_ldim(1:2)
        lambda_here = lambda / real(product(self%img_dims(1:2)))
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
        call img%bs_smooth(self%b_img, self%r_img, lambda_here)
        if (.not. img_ft_prev) call img%ifft()
    end subroutine smooth

    subroutine smooth_3d( self, img, lambda )
        class(bspline_smoother),   intent(inout) :: self
        class(image),      intent(inout) :: img
        real,              intent(in)    :: lambda ! >0.; 0.1 is a starting point
        integer :: img_ldim(3), rb_ldim(3)
        real    :: img_smpd, lambda_here
        logical :: img_ft_prev
        logical :: do_alloc
        img_ldim    = img%get_ldim()
        self%img_dims_3d(1:3) = img_ldim(1:3)
        lambda_here = lambda / real(product(self%img_dims_3d(1:3)))
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
        call img%bs_smooth(self%b_img, self%r_img, lambda_here)
        if (.not. img_ft_prev) call img%ifft()
    end subroutine smooth_3d

    subroutine fill_b(self)
        class(bspline_smoother), intent(inout) :: self
        integer :: nonzero_x(3), nonzero_y(3)         ! indices of non-zero entries in x- and y-component
        real    :: nonzero_pt_x(3),  nonzero_pt_y(3)  ! points of non-zero entries
        real    :: nonzero_val_x(3), nonzero_val_y(3) ! values of non-zero entries
        integer :: i, j, x, y
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
        call self%b_img%zero_and_unflag_ft
        do j = 1,size(nonzero_y)
            y = nonzero_y(j)
            do i = 1,size(nonzero_x)
                x = nonzero_x(i)
                call self%b_img%set([x,y,1], nonzero_val_x(i) * nonzero_val_y(j))
            end do
        end do
    end subroutine fill_b

    subroutine fill_b_3d(self)
        class(bspline_smoother), intent(inout) :: self
        integer :: nonzero_x(3),     nonzero_y(3),     nonzero_z(3)     ! indices of non-zero entries in x-, y- and z-component
        real    :: nonzero_pt_x(3),  nonzero_pt_y(3),  nonzero_pt_z(3)  ! points of non-zero entries
        real    :: nonzero_val_x(3), nonzero_val_y(3), nonzero_val_z(3) ! values of non-zero entries
        integer :: i, j, k, x, y, z
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
        call self%b_img%zero_and_unflag_ft
        do k = 1,size(nonzero_z)
            z = nonzero_z(k)
            do j = 1,size(nonzero_y)
                y = nonzero_y(j)
                do i = 1,size(nonzero_x)
                    x = nonzero_x(i)
                    call self%b_img%set([x,y,z], nonzero_val_x(i) * nonzero_val_y(j) * nonzero_val_z(k))
                end do
            end do
        end do
    end subroutine fill_b_3d

    elemental function b3(x) result(y)
        real, intent(in) :: x
        real :: y
        y = b3_help(x + 1.5)
    contains
        elemental function b3_help(xx) result(yy)
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
        class(bspline_smoother), intent(inout) :: self
        integer :: nonzero_x(5), nonzero_y(5)   ! indices of non-zero entries in x- and y-component
        real    :: nonzero_val_x_a0(5), nonzero_val_x_a2(5) ! values of non-zero entries
        real    :: nonzero_val_y_a0(5), nonzero_val_y_a2(5) ! values of non-zero entries
        integer :: i, j, x, y
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
        call self%r_img%zero_and_unflag_ft
        do j = 1, size(nonzero_x)
            y = nonzero_y(j)
            do i = 1, size(nonzero_x)
                x = nonzero_x(i)
                call self%r_img%set([x,y,1], nonzero_val_x_a0(i) * nonzero_val_y_a2(j) + nonzero_val_x_a2(i) * nonzero_val_y_a0(j))
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
        class(bspline_smoother), intent(inout) :: self
        integer :: nonzero_x(5), nonzero_y(5), nonzero_z(5)   ! indices of non-zero entries in x-, y- and z-component
        real    :: nonzero_val_x_a0(5), nonzero_val_x_a2(5) ! values of non-zero entries
        real    :: nonzero_val_y_a0(5), nonzero_val_y_a2(5) ! values of non-zero entries
        real    :: nonzero_val_z_a0(5), nonzero_val_z_a2(5) ! values of non-zero entries
        real    :: r
        integer :: i, j, k, x, y, z
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
        do i = 1, size(nonzero_z)
            z = nonzero_z(i)
            nonzero_val_z_a0(i) = a0xy(z, self%img_dims_3d(3))
            nonzero_val_z_a2(i) = a2xy(z, self%img_dims_3d(3))
        end do
        call self%r_img%zero_and_unflag_ft
        do k = 1, size(nonzero_z)
            z = nonzero_z(k)
            do j = 1, size(nonzero_y)
                y = nonzero_y(j)
                do i = 1, size(nonzero_x)
                    x = nonzero_x(i)
                    r = nonzero_val_x_a2(i) * nonzero_val_y_a0(j) * nonzero_val_z_a0(k) + &
                        nonzero_val_x_a0(i) * nonzero_val_y_a2(j) * nonzero_val_z_a0(k) + &
                        nonzero_val_x_a0(i) * nonzero_val_y_a0(j) * nonzero_val_z_a2(k)
                    call self%r_img%set([x,y,z], r)
                end do
            end do
        end do
    end subroutine fill_r_3d

    !>  TEST ROUTINE
    subroutine test_bspline_smoother( ldim, smpd, lambda )
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, lambda
        type(bspline_smoother) :: bs
        type(image)    :: img_a, img_b, img_ref_filt, img_lhs, img_rhs, img_tmp
        real           :: cc, cc_mild, cc_strong
        integer        :: i
        write(*,*) 'jhsfkgfjhwgfjhwegfjhkwgefjkhgwrgf'
        ! setup: two independent Gaussian noise images
        call img_a%new(ldim, smpd)
        call img_a%gauran(0.,1.)
        call img_b%new(ldim, smpd)
        call img_b%gauran(0.,1.)
        ! 1. reproducibility: repeated calls with same input give identical output
        write(logfhandle,'(a)') '**info(test_bspsmooth, 2D part 1): reproducibility'
        call bs%new()
        img_ref_filt = img_a
        call bs%smooth(img_ref_filt, lambda)
        call img_ref_filt%write(string('test_bspsmooth_filt.mrc'))
        do i = 1, 5
            img_tmp = img_a
            call bs%smooth(img_tmp, lambda)
            cc = img_ref_filt%real_corr(img_tmp)
            if( cc < 0.9999 ) THROW_HARD('2D reproducibility test failed; test_bspsmooth')
        end do
        call bs%kill
        ! 2. linearity: filter(2*a + 3*b) == 2*filter(a) + 3*filter(b)
        write(logfhandle,'(a)') '**info(test_bspsmooth, 2D part 2): linearity'
        call bs%new()
        img_lhs = img_a
        call img_lhs%mul(2.)
        img_tmp = img_b
        call img_tmp%mul(3.)
        call img_lhs%add(img_tmp)           ! 2*a + 3*b
        call bs%smooth(img_lhs, lambda) ! filter(2*a + 3*b)
        img_rhs = img_a
        call bs%smooth(img_rhs, lambda)
        call img_rhs%mul(2.)                ! 2*filter(a)
        img_tmp = img_b
        call bs%smooth(img_tmp, lambda)
        call img_tmp%mul(3.)                ! 3*filter(b)
        call img_rhs%add(img_tmp)           ! 2*filter(a) + 3*filter(b)
        cc = img_lhs%real_corr(img_rhs)
        if( cc < 0.9999 ) THROW_HARD('2D linearity test failed; test_bspsmooth')
        call bs%kill
        ! 3. monotone smoothing: larger lambda yields output less correlated with input
        write(logfhandle,'(a)') '**info(test_bspsmooth, 2D part 3): monotone smoothing'
        call bs%new()
        img_lhs = img_a
        call bs%smooth(img_lhs, lambda)
        cc_mild = img_a%real_corr(img_lhs)
        img_rhs = img_a
        call bs%smooth(img_rhs, max(lambda, 1.) * 100.)
        cc_strong = img_a%real_corr(img_rhs)
        if( cc_mild <= cc_strong ) THROW_HARD('2D monotone smoothing test failed; test_bspsmooth')
        call bs%kill
        ! 4. pre-allocation consistency: tv%new() and tv%new(img) give identical output
        write(logfhandle,'(a)') '**info(test_bspsmooth, 2D part 4): pre-allocation consistency'
        img_lhs = img_a
        call bs%new()
        call bs%smooth(img_lhs, lambda)
        call bs%kill
        img_rhs = img_a
        call bs%new(img_rhs)
        call bs%smooth(img_rhs, lambda)
        call bs%kill
        cc = img_lhs%real_corr(img_rhs)
        if( cc < 0.9999 ) THROW_HARD('2D pre-allocation consistency test failed; test_bspsmooth')
        ! 5. FT-state agnosticism: real-space and FT-space input give identical output
        write(logfhandle,'(a)') '**info(test_bspsmooth, 2D part 5): FT-state agnosticism'
        call bs%new()
        img_lhs = img_a
        call bs%smooth(img_lhs, lambda)   ! input: real space -> result: real space
        img_rhs = img_a
        call img_rhs%fft()
        call bs%smooth(img_rhs, lambda)   ! input: FT space -> result: FT space
        call img_rhs%ifft()
        cc = img_lhs%real_corr(img_rhs)
        if( cc < 0.9999 ) THROW_HARD('2D FT-state agnosticism test failed; test_bspsmooth')
        call bs%kill
        write(logfhandle,'(a)') 'TEST_BSPSMOOTH 2D COMPLETED SUCCESSFULLY'
        call img_a%kill
        call img_b%kill
        call img_ref_filt%kill
        call img_lhs%kill
        call img_rhs%kill
        call img_tmp%kill
    end subroutine test_bspline_smoother

    subroutine test_bspline_smoother_3d( ldim, smpd, lambda )
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, lambda
        type(bspline_smoother) :: bs
        type(image)    :: img_a, img_b, img_ref_filt, img_lhs, img_rhs, img_tmp
        real           :: cc, cc_mild, cc_strong
        integer        :: i
        ! setup: two independent Gaussian noise volumes
        call img_a%new(ldim, smpd)
        call img_a%zero_and_unflag_ft
        call img_a%gauran(0.,1.)
        call img_b%new(ldim, smpd)
        call img_b%zero_and_unflag_ft
        call img_b%gauran(0.,1.)
        ! 1. reproducibility
        write(logfhandle,'(a)') '**info(test_bspsmooth, 3D part 1): reproducibility'
        call bs%new()
        img_ref_filt = img_a
        call bs%smooth_3d(img_ref_filt, lambda)
        call img_ref_filt%write(string('test_bspsmooth_filt_3d.mrc'))
        do i = 1, 5
            img_tmp = img_a
            call bs%smooth_3d(img_tmp, lambda)
            cc = img_ref_filt%real_corr(img_tmp)
            if( cc < 0.9999 ) THROW_HARD('3D reproducibility test failed; test_bspsmooth_3d')
        end do
        call bs%kill
        ! 2. linearity: filter(2*a + 3*b) == 2*filter(a) + 3*filter(b)
        write(logfhandle,'(a)') '**info(test_bspsmooth, 3D part 2): linearity'
        call bs%new()
        img_lhs = img_a
        call img_lhs%mul(2.)
        img_tmp = img_b
        call img_tmp%mul(3.)
        call img_lhs%add(img_tmp)             ! 2*a + 3*b
        call bs%smooth_3d(img_lhs, lambda) ! filter(2*a + 3*b)
        img_rhs = img_a
        call bs%smooth_3d(img_rhs, lambda)
        call img_rhs%mul(2.)                  ! 2*filter(a)
        img_tmp = img_b
        call bs%smooth_3d(img_tmp, lambda)
        call img_tmp%mul(3.)                  ! 3*filter(b)
        call img_rhs%add(img_tmp)             ! 2*filter(a) + 3*filter(b)
        cc = img_lhs%real_corr(img_rhs)
        if( cc < 0.9999 ) THROW_HARD('3D linearity test failed; test_bspsmooth_3d')
        call bs%kill
        ! 3. monotone smoothing: larger lambda yields output less correlated with input
        write(logfhandle,'(a)') '**info(test_bspsmooth, 3D part 3): monotone smoothing'
        call bs%new()
        img_lhs = img_a
        call bs%smooth_3d(img_lhs, lambda)
        cc_mild = img_a%real_corr(img_lhs)
        img_rhs = img_a
        call bs%smooth_3d(img_rhs, max(lambda, 1.) * 100.)
        cc_strong = img_a%real_corr(img_rhs)
        if( cc_mild <= cc_strong ) THROW_HARD('3D monotone smoothing test failed; test_bspsmooth_3d')
        call bs%kill
        ! 4. FT-state agnosticism: real-space and FT-space input give identical output
        write(logfhandle,'(a)') '**info(test_bspsmooth, 3D part 4): FT-state agnosticism'
        call bs%new()
        img_lhs = img_a
        call bs%smooth_3d(img_lhs, lambda)   ! input: real space -> result: real space
        img_rhs = img_a
        call img_rhs%fft()
        call bs%smooth_3d(img_rhs, lambda)   ! input: FT space -> result: FT space
        call img_rhs%ifft()
        cc = img_lhs%real_corr(img_rhs)
        if( cc < 0.9999 ) THROW_HARD('3D FT-state agnosticism test failed; test_bspsmooth_3d')
        call bs%kill
        write(logfhandle,'(a)') 'TEST_BSPSMOOTH 3D COMPLETED SUCCESSFULLY'
        call img_a%kill
        call img_b%kill
        call img_ref_filt%kill
        call img_lhs%kill
        call img_rhs%kill
        call img_tmp%kill
    end subroutine test_bspline_smoother_3d

    !>  DESTRUCTOR
    subroutine kill_bspline_smoother( self )
        class(bspline_smoother), intent(inout) :: self
        call self%r_img%kill()
        call self%b_img%kill()
        self%interpolate_coeffs_rmat => null()
        call self%interpolate_coeffs%kill
        self%existence = .false.
    end subroutine kill_bspline_smoother

end module simple_bspline_smoother

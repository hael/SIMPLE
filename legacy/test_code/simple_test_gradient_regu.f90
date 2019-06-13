! Test file to numerically assess the validity of the gradient of the regularization term in the anisotropic motion correction
! Use simple_test_gradient_regu prg=dummy smpd=1.
module mtest_gradient_regu
include 'simple_lib.f08'
use simple_cmdline,             only: cmdline
use simple_image,               only: image
use simple_parameters,          only: parameters
use simple_strings,             only: int2str
use simple_motion_anisocor_dbl, only: motion_anisocor_dbl
use simple_defs
implicit none
#include "simple_local_flags.inc"


type :: test_gradient_regu
    type(cmdline)    :: cline
    type(parameters) :: p
    type(image)      :: aref, aframe
    type(motion_anisocor_dbl) :: motion_aniso
    integer          :: ldim(3)
contains
    procedure :: run
    procedure :: run_tests
    procedure :: alloc_imgs
    procedure :: one_imgs
    procedure :: alloc_motion_aniso
    procedure :: kill
end type test_gradient_regu

contains

    subroutine run(self)
        class(test_gradient_regu), intent(inout) :: self
        call self%cline%parse_oldschool()
        call self%cline%checkvar('smpd', 1)
        call self%cline%check()
        call self%p%new(self%cline)
        call self%alloc_imgs()
        call self%one_imgs()
        call self%alloc_motion_aniso()
        call self%run_tests()
        call self%kill()
    end subroutine run

    subroutine one_imgs(self)
        class(test_gradient_regu), intent(inout) :: self
        real, pointer :: ref_rmat(:,:,:), frame_rmat(:,:,:)
        call self%aref%get_rmat_ptr(ref_rmat)
        call self%aframe%get_rmat_ptr(frame_rmat)
        ref_rmat = 1.
        frame_rmat = 1.
    end subroutine one_imgs

    subroutine run_tests(self)
        class(test_gradient_regu), intent(inout) :: self
        real(dp) :: deltas(5) = (/0.001,0.0005,0.0001,0.00001, 0.000001/)
        integer :: idelta
        real(dp) :: delta
        real(dp) :: a(12), atmp(12)
        real(dp) :: grad(12), grad_tmp(12), f, ftmp
        real(dp) :: val
        integer  :: grad_idx
        a = 0._dp
        a(1) = 0.01
        a(2) = 0.02
        a(3) = 0.01
        a(4) = 0.005
        a(5) = -0.01
        a(6) = -0.005
        a(7) = 0.01
        a(8) = 0.002
        a(9) = 0.005
        a(10) = -0.01
        a(11) = 0.005
        a(12) = -0.001
        call self%motion_aniso%eval_fdf_foo(self%aref, self%aframe, a, f, grad)
        do idelta = 1, size(deltas)
            delta = deltas(idelta)
            write (*,*) 'delta=', delta
            do grad_idx = 1,12
                atmp = a
                atmp(grad_idx) = atmp(grad_idx) + delta
                call self%motion_aniso%eval_fdf_foo(self%aref, self%aframe, atmp, ftmp, grad_tmp)
                val = (ftmp - f) / delta
                write (*,*) 'grad(', grad_idx, ') = ', val, '    theor: ', grad(grad_idx)
            end do
        end do
    end subroutine run_tests

    subroutine alloc_imgs(self)
        class(test_gradient_regu), intent(inout) :: self

        self%ldim(1) = 100
        self%ldim(2) = 100
        self%ldim(3) = 1
        call self%aref%new(self%ldim, 1.)
        call self%aframe%new(self%ldim, 1.)
    end subroutine alloc_imgs

    subroutine alloc_motion_aniso(self)
        class(test_gradient_regu), intent(inout) :: self
        call self%motion_aniso%new()
    end subroutine alloc_motion_aniso

    subroutine kill(self)
        class(test_gradient_regu), intent(inout) :: self
        call self%aframe%kill()
        call self%aref%kill()
    end subroutine kill

end module mtest_gradient_regu

program simple_test_gradient_regu
    use mtest_gradient_regu
    implicit none
    type(test_gradient_regu) :: atest_gradient_regu
    call atest_gradient_regu%run()
end program simple_test_gradient_regu

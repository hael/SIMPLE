!@descr: unit tests for class-average registration on the polar Fourier transform (match_imgs, match_imgs2ref)
! Moved out of simple_strategy2D_utils (its self-test test_cavg_registration, which the
! cavg_registration test program ran) by the utils review (plan, section 9.7); the two matchers
! are exported for it. Five copies of an asymmetric image of three Gaussians, rotated in 30 degree
! steps, are registered all against all and against the first; then the copies are also shifted
! by 0.25 (i-1) pixels in x and y and registered against the first, and the rotation and shift
! found are compared with the ones applied. The sign convention of e3 is not pinned (the found
! angle may be the applied one or its negative), and the shift is compared by length, which
! does not depend on whether the rotation acts before or after the shift.
module simple_cavg_registration_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image,            only: image
use simple_imgarr_utils,     only: dealloc_imgarr
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_strategy2D_utils, only: match_imgs, match_imgs2ref
use simple_test_utils
implicit none
private
public :: run_all_cavg_registration_tests

contains

    subroutine run_all_cavg_registration_tests()
        write(*,'(A)') '**** running all class-average registration tests ****'
        call test_rotated_copies()
    end subroutine run_all_cavg_registration_tests

    subroutine test_rotated_copies()
        integer, parameter :: NIMGS     = 5, BOX = 64
        real,    parameter :: SMPD      = 1.5
        real,    parameter :: CORR_TOL  = 0.95
        real,    parameter :: ANG_TOL   = 3.6   ! one polar step: mskdiam 48 A at 1.5 A is a radius of 16 px, 100 rotations
        real,    parameter :: SHIFT_TOL = 0.5   ! pixels
        type(inpl_struct), allocatable :: alg_info1(:), alg_info2(:,:)
        type(image),       allocatable :: imgs_ref(:), imgs_targ(:)
        type(parameters)   :: params
        type(cmdline)      :: cline
        real, allocatable  :: rmat(:,:,:)
        real               :: rx, ry, ang, dang, shift_applied, shift_found
        integer            :: i, x, y
        character(len=1)   :: ci
        write(*,'(A)') 'test_rotated_copies'
        call cline%set('smpd',    SMPD)
        call cline%set('lp',      6.)
        call cline%set('hp',      20.)
        call cline%set('nthr',    1)
        call cline%set('trs',     8.)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        call cline%set('box',     BOX)
        call cline%set('mskdiam', 48.)
        call params%new(cline)
        allocate(imgs_ref(NIMGS), imgs_targ(NIMGS), rmat(BOX,BOX,1))
        do y = 1, BOX
            ry = real(y - (BOX / 2 + 1))
            do x = 1, BOX
                rx = real(x - (BOX / 2 + 1))
                rmat(x,y,1) = exp(-((rx - 8.)**2 + (ry + 6.)**2) / 18.) + &
                    &0.7 * exp(-((rx + 7.)**2 + (ry - 4.)**2) / 32.) + &
                    &0.4 * exp(-((rx - 2.)**2 + (ry - 10.)**2) / 8.)
            enddo
        enddo
        call imgs_ref(1)%new([BOX,BOX,1], SMPD)
        call imgs_ref(1)%set_rmat(rmat, .false.)
        call imgs_ref(1)%norm
        call imgs_targ(1)%copy(imgs_ref(1))
        do i = 2, NIMGS
            call imgs_ref(i)%copy(imgs_ref(1))
            call imgs_ref(i)%rtsq(real(i - 1) * 30., 0., 0.)
            call imgs_targ(i)%copy(imgs_ref(i))
        enddo
        ! all against all: every pair of copies registers to a high correlation
        alg_info2 = match_imgs(params, params%hp, params%lp, params%trs, imgs_ref, imgs_targ)
        call assert_true(all(ieee_is_finite(alg_info2%corr)), 'match_imgs: every correlation is finite')
        call assert_true(all(alg_info2%corr >= CORR_TOL), 'match_imgs: every pair of rotated copies registers at 0.95 or more')
        ! rotated and shifted copies against the first: correlation, rotation and shift
        do i = 2, NIMGS
            call imgs_targ(i)%copy(imgs_ref(1))
            call imgs_targ(i)%rtsq(real(i - 1) * 30., 0.25 * real(i - 1), 0.25 * real(i - 1))
        enddo
        alg_info1 = match_imgs2ref(params, params%hp, params%lp, params%trs, imgs_ref(1), imgs_targ)
        do i = 1, NIMGS
            write(ci,'(I1)') i
            call assert_true(ieee_is_finite(alg_info1(i)%corr) .and. alg_info1(i)%corr >= CORR_TOL, &
                &'match_imgs2ref: copy '//ci//' registers at 0.95 or more')
            ang  = real(i - 1) * 30.
            dang = min(angdist(alg_info1(i)%e3, ang), angdist(alg_info1(i)%e3, -ang))
            call assert_real(0., dang, ANG_TOL, 'match_imgs2ref: copy '//ci//' rotation within one polar step of the applied one')
            shift_applied = 0.25 * real(i - 1) * sqrt(2.)
            shift_found   = sqrt(alg_info1(i)%x**2 + alg_info1(i)%y**2)
            call assert_real(shift_applied, shift_found, SHIFT_TOL, 'match_imgs2ref: copy '//ci//' shift length within 0.5 px')
            call assert_false(alg_info1(i)%l_mirr, 'match_imgs2ref: copy '//ci//' is not taken for a mirror')
        enddo
        call dealloc_imgarr(imgs_ref)
        call dealloc_imgarr(imgs_targ)
        deallocate(alg_info1, alg_info2, rmat)
        call cline%kill

      contains

        !> angular distance in degrees, in [0, 180]
        real function angdist( a, b )
            real, intent(in) :: a, b
            angdist = abs(modulo(a - b + 180., 360.) - 180.)
        end function angdist

    end subroutine test_rotated_copies

end module simple_cavg_registration_tester

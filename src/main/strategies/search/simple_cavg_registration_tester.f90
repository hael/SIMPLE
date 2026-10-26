!@descr: unit tests for class-average registration on the polar Fourier transform (match_imgs, match_imgs2ref)
! Five copies of an asymmetric three-Gaussian image, rotated in 30 degree steps, are registered all
! against all; shifted also by 0.25 (i-1) px in x and y, against the first, where the found rotation
! and shift are compared with the applied ones. The sign of e3 is not pinned, and the shift is compared
! by length, which does not depend on whether the rotation acts before or after the shift.
module simple_cavg_registration_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image,            only: image
use simple_imgarr_utils,     only: alloc_imgarr, dealloc_imgarr
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_strategy2D_utils, only: flag_non_junk_cavgs, match_imgs, match_imgs2ref, rtsq_imgs
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
        !$ use omp_lib, only: omp_set_num_threads
        integer, parameter :: NIMGS     = 5, BOX = 64
        real,    parameter :: SMPD      = 1.5
        real,    parameter :: CORR_TOL  = 0.95
        real,    parameter :: SPATIAL_CORR_TOL = 0.90
        real,    parameter :: ANG_TOL   = 3.6   ! one polar step: mskdiam 48 A at 1.5 A is a radius of 16 px, 100 rotations
        real,    parameter :: SHIFT_TOL = 0.5   ! pixels
        type(inpl_struct), allocatable :: alg_info1(:), alg_info1_serial(:), alg_info2(:,:), alg_info2_serial(:,:)
        type(inpl_struct) :: one_info(1)
        type(image),       allocatable :: imgs_ref(:), imgs_targ(:), imgs_ref_serial(:), imgs_targ_serial(:), aligned(:)
        type(image) :: one_aligned(1)
        type(parameters)   :: params
        type(cmdline)      :: cline
        real, allocatable  :: rmat(:,:,:)
        real               :: rx, ry, ang, dang, shift_applied, shift_found, spatial_corr
        integer            :: i, j, x, y, nthr_saved
        logical, allocatable :: non_junk(:), non_junk_serial(:)
        character(len=1)   :: ci
        character(len=128) :: message
        write(*,'(A)') 'test_rotated_copies'
        nthr_saved = nthr_glob
        call cline%set('smpd',    SMPD)
        call cline%set('lp',      6.)
        call cline%set('hp',      20.)
        call cline%set('nthr',    3)
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
        call assert_true(imgs_ref(1)%is_wthreads(), 'registration fixtures use threaded source images')
        call alloc_imgarr(NIMGS, [BOX,BOX,1], SMPD, imgs_ref_serial, wthreads=.false.)
        call alloc_imgarr(NIMGS, [BOX,BOX,1], SMPD, imgs_targ_serial, wthreads=.false.)
        do i = 1,NIMGS
            call imgs_ref_serial(i)%copy_fast(imgs_ref(i))
            call imgs_targ_serial(i)%copy_fast(imgs_targ(i))
        enddo
        call flag_non_junk_cavgs(imgs_ref, 20., 16., non_junk)
        call flag_non_junk_cavgs(imgs_ref_serial, 20., 16., non_junk_serial)
        call assert_true(all(non_junk .eqv. non_junk_serial), &
            &'flag_non_junk_cavgs gives the same decisions for threaded and serial sources')
        ! all against all: every pair of copies registers to a high correlation
        alg_info2 = match_imgs(params, params%hp, params%lp, params%trs, imgs_ref, imgs_targ)
        alg_info2_serial = match_imgs(params, params%hp, params%lp, params%trs, imgs_ref_serial, imgs_targ_serial)
        call assert_true(all(abs(alg_info2%corr - alg_info2_serial%corr) < 1.e-6), &
            &'match_imgs gives the same correlations for threaded and serial sources')
        call assert_true(all(abs(alg_info2%e3 - alg_info2_serial%e3) < 1.e-6), &
            &'match_imgs gives the same rotations for threaded and serial sources')
        call assert_true(all(abs(alg_info2%x - alg_info2_serial%x) < 1.e-6) .and. &
            &all(abs(alg_info2%y - alg_info2_serial%y) < 1.e-6), &
            &'match_imgs gives the same shifts for threaded and serial sources')
        call assert_true(all(alg_info2%l_mirr .eqv. alg_info2_serial%l_mirr), &
            &'match_imgs gives the same mirror decisions for threaded and serial sources')
        call assert_true(all(ieee_is_finite(alg_info2%corr)), 'match_imgs: every correlation is finite')
        call assert_true(all(alg_info2%corr >= CORR_TOL), 'match_imgs: every pair of rotated copies registers at 0.95 or more')
        do i = 1, NIMGS
            do j = 1, NIMGS
                write(message,'(A,I0,A,I0,A,F7.4,A,F5.2)') 'match_imgs reference ', i, ', target ', j, &
                    &' search correlation=', alg_info2(i,j)%corr, ', minimum=', CORR_TOL
                call assert_true(alg_info2(i,j)%corr >= CORR_TOL, trim(message))
                call one_aligned(1)%copy(imgs_targ(j))
                one_info(1) = alg_info2(i,j)
                call rtsq_imgs(1, one_info, one_aligned)
                spatial_corr = imgs_ref(i)%real_corr(one_aligned(1))
                write(message,'(A,I0,A,I0,A,F7.4,A,F5.2)') 'match_imgs reference ', i, ', target ', j, &
                    &' pixel correlation=', spatial_corr, ', minimum=', SPATIAL_CORR_TOL
                call assert_true(spatial_corr >= SPATIAL_CORR_TOL, trim(message))
                call one_aligned(1)%kill
            enddo
        enddo
        ! rotated and shifted copies against the first: correlation, rotation and shift
        do i = 2, NIMGS
            call imgs_targ(i)%copy(imgs_ref(1))
            call imgs_targ(i)%rtsq(real(i - 1) * 30., 0.25 * real(i - 1), 0.25 * real(i - 1))
            call imgs_targ_serial(i)%copy_fast(imgs_targ(i))
        enddo
        alg_info1 = match_imgs2ref(params, params%hp, params%lp, params%trs, imgs_ref(1), imgs_targ)
        alg_info1_serial = match_imgs2ref(params, params%hp, params%lp, params%trs, imgs_ref_serial(1), imgs_targ_serial)
        call assert_true(all(abs(alg_info1%corr - alg_info1_serial%corr) < 1.e-6), &
            &'match_imgs2ref gives the same correlations for threaded and serial sources')
        call assert_true(all(abs(alg_info1%e3 - alg_info1_serial%e3) < 1.e-6), &
            &'match_imgs2ref gives the same rotations for threaded and serial sources')
        call assert_true(all(abs(alg_info1%x - alg_info1_serial%x) < 1.e-6) .and. &
            &all(abs(alg_info1%y - alg_info1_serial%y) < 1.e-6), &
            &'match_imgs2ref gives the same shifts for threaded and serial sources')
        call assert_true(all(alg_info1%l_mirr .eqv. alg_info1_serial%l_mirr), &
            &'match_imgs2ref gives the same mirror decisions for threaded and serial sources')
        allocate(aligned(NIMGS))
        do i = 1, NIMGS
            call aligned(i)%copy(imgs_targ(i))
        enddo
        call rtsq_imgs(NIMGS, alg_info1, aligned)
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
            spatial_corr = imgs_ref(1)%real_corr(aligned(i))
            write(message,'(A,I0,A,F7.4,A,F5.2)') 'match_imgs2ref target ', i, &
                &' pixel correlation=', spatial_corr, ', minimum=', SPATIAL_CORR_TOL
            call assert_true(spatial_corr >= SPATIAL_CORR_TOL, trim(message))
        enddo
        call dealloc_imgarr(imgs_ref)
        call dealloc_imgarr(imgs_targ)
        call dealloc_imgarr(imgs_ref_serial)
        call dealloc_imgarr(imgs_targ_serial)
        call dealloc_imgarr(aligned)
        deallocate(alg_info1, alg_info1_serial, alg_info2, alg_info2_serial, non_junk, non_junk_serial, rmat)
        call cline%kill
        nthr_glob = nthr_saved
        !$ call omp_set_num_threads(nthr_saved)

      contains

        !> angular distance in degrees, in [0, 180]
        real function angdist( a, b )
            real, intent(in) :: a, b
            angdist = abs(modulo(a - b + 180., 360.) - 180.)
        end function angdist

    end subroutine test_rotated_copies

end module simple_cavg_registration_tester

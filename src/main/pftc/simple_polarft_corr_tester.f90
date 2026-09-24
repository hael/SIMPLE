!@descr: unit tests of the polar-Fourier correlations (gen_objfun_vals) on generated images, through the production polarisation path
! Moved from the gencorrs_fft test program by the review of the remaining tests (plan, section 9.7):
! three smooth zero-mean random images are polarised into a polarft_calc as references and, in
! rotated or unrelated form, as particles; gen_objfun_vals (objfun=cc) must then peak at rotation 1
! with correlation ~1 for an image against itself, at the applied rotation (a sixth of a turn) for
! a rotated copy, and stay low for an unrelated image. The noise is zero-mean so that the shared
! mask envelope carries no correlation.
module simple_polarft_corr_tester
use simple_core_module_api
use simple_cmdline,      only: cmdline
use simple_parameters,   only: parameters
use simple_polarft_calc, only: polarft_calc
use simple_image,        only: image
use simple_test_utils
implicit none
private
public :: run_all_polarft_corr_tests

contains

    subroutine run_all_polarft_corr_tests()
        write(*,'(A)') '**** running all polar correlation tests ****'
        call test_generated_images()
    end subroutine run_all_polarft_corr_tests

    subroutine test_generated_images()
        integer, parameter :: BOX = 64, NIMGS = 3
        real,    parameter :: SMPD = 1.0, MSKDIAM = 48.0, LP = 6.0, MSKRAD = 22.0
        type(parameters), target :: p
        type(cmdline)       :: cline
        type(polarft_calc)  :: pftc
        type(image)         :: imgs(NIMGS), ptcls(NIMGS)
        real, allocatable   :: cc(:)
        real    :: ang, dang
        integer :: pdim(3), kfromto(2), nrots, rotstep, i, loc, loc_fwd, loc_bwd
        logical :: peak_ok
        write(*,'(A)') 'test_generated_images'
        ! a polarft_calc from parameters alone (no stack on disk)
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nptcls',  real(NIMGS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        call p%new(cline)
        ! fixed seed after parameters%new, which reseeds (seed_rnd): the fixture is generated, so the run is reproducible
        call set_fixed_seed(20260922)
        kfromto = [2, calc_fourier_index(LP, BOX, SMPD)]
        call pftc%new(p, NIMGS, [1,NIMGS], kfromto)
        pdim  = pftc%get_pdim_srch()
        nrots = pftc%get_nrots()
        dang    = 360.0 / real(nrots)
        rotstep = nrots / 6
        call assert_true(rotstep >= 3, 'enough in-plane rotations for the probe step')
        ! smooth, asymmetric, masked test images: Gaussian low-passed zero-mean Gaussian noise
        do i = 1,NIMGS
            call imgs(i)%new([BOX,BOX,1], SMPD)
            call imgs(i)%gauran(0.0, 1.0)
            call imgs(i)%fft
            call imgs(i)%bpgau2D(0.0, LP)
            call imgs(i)%ifft
            call imgs(i)%memoize_mask_coords
            call imgs(i)%mask2D_soft(MSKRAD, backgr=0.0)
        end do
        ! particles: 1 = image 1 itself, 2 = image 2 rotated by rotstep polar steps, 3 = image 3 (unrelated to reference 1)
        ang = real(rotstep) * dang
        call ptcls(1)%copy(imgs(1))
        call imgs(2)%rtsq(ang, 0.0, 0.0, ptcls(2))
        call ptcls(3)%copy(imgs(3))
        ! image -> polar transform through the production path
        call imgs(1)%memoize4polarize(pdim)
        do i = 1,NIMGS
            call imgs(i)%fft
            call ptcls(i)%fft
            call pftc%polarize_ref_pft(imgs(i), i, iseven=.true., pdim=pdim, oversamp=.false.)
            call pftc%polarize_ptcl_pft(ptcls(i), i, pdim=pdim, oversamp=.false.)
            call pftc%set_eo(i, .true.)
        end do
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        allocate(cc(nrots))
        ! an image against itself: correlation 1 at rotation 1, nowhere higher
        call pftc%gen_objfun_vals(1, 1, [0.0,0.0], cc)
        loc = maxloc(cc, dim=1)
        call assert_int(1, loc,                    'self-correlation peaks at rotation 1')
        call assert_real(1.0, cc(1), 1.0e-3,       'self-correlation at rotation 1 is 1')
        call assert_true(all(cc <= cc(1) + 1.0e-5), 'no rotation correlates better than the identity')
        ! a rotated copy: the peak sits rotstep polar steps from the identity, in either
        ! angular convention (the sign of the real-space rotation is not what is tested here)
        call pftc%gen_objfun_vals(2, 2, [0.0,0.0], cc)
        loc     = maxloc(cc, dim=1)
        loc_fwd = rotstep + 1
        loc_bwd = nrots - rotstep + 1
        peak_ok = abs(loc - loc_fwd) <= 1 .or. abs(loc - loc_bwd) <= 1
        write(logfhandle,'(A,I4,A,F7.2,A,I4,A,I4,A,F6.3)') 'rotated copy: peak index', loc, ' (', pftc%get_rot(loc),&
            &' deg), expected', loc_fwd, ' or', loc_bwd, ', cc=', cc(loc)
        call assert_true(peak_ok,                  'rotated copy peaks at the applied rotation (within one step)')
        call assert_true(cc(loc) > 0.9,            'rotated copy correlates above 0.9 at its peak')
        call assert_true(cc(1) < cc(loc) - 0.3,    'rotated copy does not peak at the identity')
        ! an unrelated image: no rotation correlates strongly
        call pftc%gen_objfun_vals(1, 3, [0.0,0.0], cc)
        call assert_true(maxval(cc) < 0.5,         'unrelated image stays below 0.5 at every rotation')
        ! cleanup
        do i = 1,NIMGS
            call imgs(i)%kill
            call ptcls(i)%kill
        end do
        call pftc%kill
        deallocate(cc)
    end subroutine test_generated_images

end module simple_polarft_corr_tester

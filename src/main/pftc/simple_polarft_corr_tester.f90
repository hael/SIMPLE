!@descr: unit tests of the polar-Fourier correlations (gen_objfun_vals, calc_frc) on generated images, through the production polarisation path
! Moved from the gencorrs_fft test program by the review of the remaining tests (plan, section 9.7):
! three smooth zero-mean random images are polarised into a polarft_calc as references and, in
! rotated or unrelated form, as particles; gen_objfun_vals (objfun=cc) must then peak at rotation 1
! with correlation ~1 for an image against itself, at the applied rotation (a sixth of a turn) for
! a rotated copy, and stay low for an unrelated image. The noise is zero-mean so that the shared
! mask envelope carries no correlation. calc_frc rotates the reference with rotate_ref_8, not the
! FFT path of gen_objfun_vals, so its shell correlations must peak at the same rotation; the probe
! rotations (+60, -60 and 180 degrees) put the peaks in each branch of rotate_ref_8.
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
        integer, parameter :: BOX = 64, NREFS = 3, NPTCLS = 5
        real,    parameter :: SMPD = 1.0, MSKDIAM = 48.0, LP = 6.0, MSKRAD = 22.0
        type(parameters), target :: p
        type(cmdline)       :: cline
        type(polarft_calc)  :: pftc
        type(image)         :: imgs(NREFS), ptcls(NPTCLS)
        real, allocatable   :: cc(:)
        real    :: ang, dang
        integer :: pdim(3), kfromto(2), nrots, pftsz, rotstep, i, loc, loc_fwd, loc_bwd, loc_plus, loc_minus
        logical :: peak_ok
        write(*,'(A)') 'test_generated_images'
        ! a polarft_calc from parameters alone (no stack on disk)
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nptcls',  real(NPTCLS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        call p%new(cline)
        ! fixed seed after parameters%new, which reseeds (seed_rnd): the fixture is generated, so the run is reproducible
        call set_fixed_seed(20260922)
        kfromto = [2, calc_fourier_index(LP, BOX, SMPD)]
        call pftc%new(p, NREFS, [1,NPTCLS], kfromto)
        pdim  = pftc%get_pdim_srch()
        nrots = pftc%get_nrots()
        pftsz = pftc%get_pftsz()
        dang    = 360.0 / real(nrots)
        rotstep = nrots / 6
        call assert_true(rotstep >= 3, 'enough in-plane rotations for the probe step')
        ! smooth, asymmetric, masked test images: Gaussian low-passed zero-mean Gaussian noise
        do i = 1,NREFS
            call imgs(i)%new([BOX,BOX,1], SMPD)
            call imgs(i)%gauran(0.0, 1.0)
            call imgs(i)%fft
            call imgs(i)%bpgau2D(0.0, LP)
            call imgs(i)%ifft
            call imgs(i)%memoize_mask_coords
            call imgs(i)%mask2D_soft(MSKRAD, backgr=0.0)
        end do
        ! particles: 1 = image 1 itself, 2 = image 2 rotated by +rotstep polar steps, 3 = image 3
        ! (unrelated to reference 1), 4 = image 2 rotated by -rotstep polar steps, 5 = image 1 turned
        ! by 180 degrees
        ang = real(rotstep) * dang
        call ptcls(1)%copy(imgs(1))
        call imgs(2)%rtsq( ang,   0.0, 0.0, ptcls(2))
        call ptcls(3)%copy(imgs(3))
        call imgs(2)%rtsq(-ang,   0.0, 0.0, ptcls(4))
        call imgs(1)%rtsq(180.0,  0.0, 0.0, ptcls(5))
        ! image -> polar transform through the production path
        call imgs(1)%memoize4polarize(pdim)
        do i = 1,NREFS
            call imgs(i)%fft
            call pftc%polarize_ref_pft(imgs(i), i, iseven=.true., pdim=pdim, oversamp=.false.)
        end do
        do i = 1,NPTCLS
            call ptcls(i)%fft
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
        ! calc_frc (rotate_ref_8) against gen_objfun_vals (FFT): the identity branch, both halves of
        ! the in-plane range and the half-turn branch
        call check_frc_identity()
        call check_frc_peak(2, 2, 'rotated by +60 degrees', loc_plus)
        call check_frc_peak(2, 4, 'rotated by -60 degrees', loc_minus)
        call assert_true((loc_plus  >= 2 .and. loc_plus  <= pftsz .and. loc_minus > pftsz + 1) .or.&
                        &(loc_minus >= 2 .and. loc_minus <= pftsz .and. loc_plus  > pftsz + 1),&
                        &'the +60 and -60 degree peaks fall in opposite halves of the in-plane range')
        call check_frc_peak(1, 5, 'turned by 180 degrees', loc)
        call assert_true(cyclic_dist(loc, pftsz + 1) <= 1, 'the half turn peaks at rotation pftsz+1 (within one step)')
        ! cleanup
        do i = 1,NREFS
            call imgs(i)%kill
        end do
        do i = 1,NPTCLS
            call ptcls(i)%kill
        end do
        call pftc%kill
        deallocate(cc)

      contains

        ! the reference against its own image, unrotated: every shell correlates perfectly
        subroutine check_frc_identity()
            real :: frc(kfromto(1):kfromto(2))
            call pftc%calc_frc(1, 1, 1, [0.0,0.0], frc)
            call assert_true(all(abs(frc - 1.0) < 1.0e-4), 'calc_frc: an image against itself at rotation 1 is 1 in every shell')
        end subroutine check_frc_identity

        ! the rotation that maximises the mean shell correlation of calc_frc is the peak of
        ! gen_objfun_vals (within one step), and the rotated copy correlates there in every shell
        subroutine check_frc_peak( iref, iptcl, label, loc_frc )
            integer,          intent(in)  :: iref, iptcl
            character(len=*), intent(in)  :: label
            integer,          intent(out) :: loc_frc
            real    :: frc(kfromto(1):kfromto(2)), frc_mean(nrots), frc_at_peak(kfromto(1):kfromto(2))
            integer :: irot, loc_cc
            call pftc%gen_objfun_vals(iref, iptcl, [0.0,0.0], cc)
            loc_cc = maxloc(cc, dim=1)
            do irot = 1,nrots
                call pftc%calc_frc(iref, iptcl, irot, [0.0,0.0], frc)
                frc_mean(irot) = sum(frc) / real(size(frc))
            end do
            loc_frc = maxloc(frc_mean, dim=1)
            call pftc%calc_frc(iref, iptcl, loc_frc, [0.0,0.0], frc_at_peak)
            write(logfhandle,'(A,I4,A,I4,A,F6.3)') 'calc_frc '//label//': peak index', loc_frc, ', gen_objfun_vals', loc_cc,&
                &', mean FRC', frc_mean(loc_frc)
            call assert_true(cyclic_dist(loc_frc, loc_cc) <= 1, 'calc_frc '//label//': peaks where gen_objfun_vals does (within one step)')
            call assert_true(frc_mean(loc_frc) > 0.9,           'calc_frc '//label//': mean shell correlation above 0.9 at the peak')
            call assert_true(minval(frc_at_peak) > 0.5,         'calc_frc '//label//': every shell correlates at the peak')
        end subroutine check_frc_peak

        integer function cyclic_dist( a, b )
            integer, intent(in) :: a, b
            cyclic_dist = min(abs(a - b), nrots - abs(a - b))
        end function cyclic_dist

    end subroutine test_generated_images

end module simple_polarft_corr_tester

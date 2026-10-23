!@descr: unit tests of the polar-Fourier correlations (gen_objfun_vals, calc_frc, the cross-half, at-pose and seed-bank Euclidean losses, the in-plane posterior draw, the slot detector, the reference model loaded at an offset) on generated images
! Three smooth zero-mean images (so the shared mask envelope carries no correlation) are references
! and, rotated or unrelated, particles. gen_objfun_vals (objfun=cc) must peak at rotation 1 with cc ~1
! for an image against itself, at the applied sixth of a turn for a rotated copy, and stay low for an
! unrelated image. calc_frc (rotate_ref_8, not the FFT path) must peak at the same rotations; the probes
! (+60, -60 and 180 degrees) cover each branch of rotate_ref_8.
module simple_polarft_corr_tester
use simple_core_module_api
use simple_cmdline,      only: cmdline
use simple_parameters,   only: parameters
use simple_polarft_calc, only: polarft_calc
use simple_pftc_shsrch_grad, only: shift_ptcl2ref
use simple_image,        only: image
use simple_test_utils
implicit none
private
public :: run_all_polarft_corr_tests

contains

    subroutine run_all_polarft_corr_tests()
        write(*,'(A)') '**** running all polar correlation tests ****'
        call test_generated_images()
        call test_cross_half_euclid()
        call test_ref_pfts_offset_io()
        call test_bank_and_at_pose_kernels()
        call test_posterior_inpl_draw()
    end subroutine run_all_polarft_corr_tests

    !> The in-plane posterior draw (gen_posterior_val and its bank form, prob_assign=posterior): at a huge
    !! temperature the argmin of the same kernel (gen_best_objfun_val, gen_best_euclid_val_bank), with its loss;
    !! at temperature 0 uniform over the grid (many distinct rotations, the argmin a minority), the returned loss
    !! the at-pose loss of the drawn rotation; at the particle's own temperature (posterior_temperature) the
    !! draw's distribution is the posterior computed from the row (gen_objfun_vals): the mode's frequency and the
    !! mass within two grid steps of it over 200 draws agree with the computed masses to 0.1, in two regimes. At
    !! unit sigma2 the fixture's temperature is near zero (sigma2 = 1 is enormous against the image's power in
    !! the calculator's units), the posterior is uniform and so are the draws; with sigma2 at the image's own
    !! power scale, chosen so that the temperature is 100 (a constant sigma2 leaves the normalized loss
    !! unchanged and scales the temperature as its inverse), the posterior and the draws concentrate about the
    !! mode.
    subroutine test_posterior_inpl_draw()
        use simple_eul_prob_tab_utils, only: posterior_temperature
        integer, parameter :: BOX = 64, NREFS = 2, NPTCLS = 1, NDRAWS = 200
        real,    parameter :: SMPD = 1.0, MSKDIAM = 48.0, LP = 6.0, MSKRAD = 22.0
        type(parameters), target :: p
        type(cmdline)       :: cline
        type(polarft_calc)  :: pftc
        type(image)         :: img_self, img_other
        real,    allocatable, target :: sigma2(:,:)
        real,    allocatable :: vals(:)
        real(dp), allocatable :: w(:)
        complex, allocatable :: bank(:,:)
        logical, allocatable :: seen(:)
        real     :: s(2), u(2), dist, dist_b, dist_at, temp, max_err, p_near, f_near
        integer  :: pdim(3), kfromto(2), irot, irot_b, nrots, idraw, nmode, ndistinct, noff
        write(*,'(A)') 'test_posterior_inpl_draw'
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nptcls',  real(NPTCLS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline)
        call set_fixed_seed(20261025)
        kfromto = [2, calc_fourier_index(LP, BOX, SMPD)]
        call pftc%new(p, NREFS, [1,NPTCLS], kfromto)
        pdim  = pftc%get_pdim_srch()
        nrots = pftc%get_nrots()
        call img_self%new([BOX,BOX,1], SMPD)
        call img_self%gauran(0.0, 1.0)
        call img_self%fft
        call img_self%bpgau2D(0.0, LP)
        call img_self%ifft
        call img_self%memoize_mask_coords
        call img_self%mask2D_soft(MSKRAD, backgr=0.0)
        call img_other%new([BOX,BOX,1], SMPD)
        call img_other%gauran(0.0, 1.0)
        call img_other%fft
        call img_other%bpgau2D(0.0, LP)
        call img_other%ifft
        call img_other%memoize_mask_coords
        call img_other%mask2D_soft(MSKRAD, backgr=0.0)
        call img_self%memoize4polarize(pdim)
        call img_self%fft
        call img_other%fft
        ! reference 1 is the particle itself in both slots, reference 2 an unrelated image
        call pftc%polarize_ref_pft(img_self,  1, iseven=.true.,  pdim=pdim, oversamp=.false.)
        call pftc%polarize_ref_pft(img_self,  1, iseven=.false., pdim=pdim, oversamp=.false.)
        call pftc%polarize_ref_pft(img_other, 2, iseven=.true.,  pdim=pdim, oversamp=.false.)
        call pftc%polarize_ref_pft(img_other, 2, iseven=.false., pdim=pdim, oversamp=.false.)
        call pftc%polarize_ptcl_pft(img_self, 1, pdim=pdim, oversamp=.false.)
        call pftc%set_eo(1, .true.)
        allocate(sigma2(kfromto(1):kfromto(2),2), source=1.0)
        call pftc%assign_sigma2_noise(sigma2)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        allocate(bank(pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1), seen(nrots), vals(nrots), w(nrots))
        ! a huge temperature: the argmin of the same kernel, at a nonzero reference shift
        s = [1.3, -0.6]
        call pftc%gen_best_objfun_val(2, 1, s, dist_b, irot_b)
        call pftc%gen_posterior_val(2, 1, s, huge(temp), dist, irot)
        call assert_int(irot_b, irot, 'a huge temperature draws the argmin rotation')
        call assert_real(dist_b, dist, 1.e-6, 'a huge temperature returns the argmin loss')
        ! the bank form against gen_best_euclid_val_bank
        u = [0.9, 1.4]
        call pftc%memoize_shifted_ptcl(1, u, bank)
        call pftc%gen_best_euclid_val_bank(2, 1, bank, dist_b, irot_b)
        call pftc%gen_posterior_val_bank(2, 1, bank, huge(temp), dist, irot)
        call assert_int(irot_b, irot, 'bank form: a huge temperature draws the argmin rotation')
        call assert_real(dist_b, dist, 1.e-6, 'bank form: a huge temperature returns the argmin loss')
        ! temperature 0: uniform over the grid; the returned loss is the drawn rotation's at-pose loss
        call pftc%gen_best_objfun_val(1, 1, s, dist_b, irot_b)
        seen    = .false.
        nmode   = 0
        noff    = 0
        max_err = 0.
        do idraw = 1,NDRAWS
            call pftc%gen_posterior_val(1, 1, s, 0., dist, irot)
            if( irot < 1 .or. irot > nrots )then
                noff = noff + 1
                cycle
            endif
            seen(irot) = .true.
            if( irot == irot_b ) nmode = nmode + 1
            if( idraw <= 10 )then
                call pftc%gen_euclid_val_at(1, 1, s, irot, dist_at)
                max_err = max(max_err, abs(dist_at - dist))
            endif
        enddo
        call assert_int(0, noff, 'temperature 0: every drawn rotation is on the grid')
        call assert_true(max_err <= 1.e-5, 'the returned loss is the at-pose loss of the drawn rotation')
        ndistinct = count(seen)
        call assert_true(ndistinct >= min(nrots, NDRAWS) / 4, 'temperature 0 spreads the draws over the grid')
        call assert_true(nmode < NDRAWS / 2, 'temperature 0 draws the argmin a minority of the time')
        ! the particle's own temperature: the draw's distribution is the posterior computed from the row. At unit
        ! sigma2 the temperature is near zero and the posterior uniform
        temp = posterior_temperature(pftc%wsqsum_ptcl_now(1), pftc%get_pftsz())
        call assert_true(temp > 0. .and. temp < huge(temp), 'the posterior temperature is positive and finite')
        call pftc%gen_best_objfun_val(1, 1, [0.,0.], dist_b, irot_b)
        call assert_int(1, irot_b, 'the image against itself peaks at rotation 1')
        call draws_against_posterior(temp, p_near, f_near)
        call assert_true(p_near < 0.5, 'at unit sigma2 the fixture''s posterior over rotations is wide')
        ! sigma2 at the image's own power scale, so that the temperature is 100: the normalized loss is unchanged
        ! (a constant sigma2 cancels), the temperature scales as 1/sigma2, and the draws concentrate about the mode
        sigma2 = temp / 100.
        call pftc%memoize_sigma2_weights
        temp = posterior_temperature(pftc%wsqsum_ptcl_now(1), pftc%get_pftsz())
        call assert_real(100., temp, 0.1, 'the temperature scales as the inverse of a constant sigma2')
        call draws_against_posterior(temp, p_near, f_near)
        call assert_true(p_near > 0.5 .and. f_near > 0.5, 'at the image''s power scale the draws concentrate about the mode')
        call img_self%kill
        call img_other%kill
        call pftc%kill
        deallocate(sigma2, bank, seen, vals, w)

    contains

        !> the posterior over the rotations of reference 1 against the particle at temperature temp_loc, from
        !! every rotation's loss (gen_objfun_vals), against NDRAWS draws: the mode's mass and frequency, and the
        !! mass and frequency within two grid steps of the mode, must agree to 0.1
        subroutine draws_against_posterior( temp_loc, p_near_loc, f_near_loc )
            real, intent(in)  :: temp_loc
            real, intent(out) :: p_near_loc, f_near_loc
            real(dp) :: wsum
            real     :: lmin, p_mode, f_mode, dist_loc
            integer  :: irot_loc, idraw_loc, dstep, nmode_loc, nnear_loc
            call pftc%gen_objfun_vals(1, 1, [0.,0.], vals)  ! exp(-L_r)
            lmin = -log(max(maxval(vals), 1.e-30))
            wsum = 0.d0
            do irot_loc = 1,nrots
                w(irot_loc) = exp(-real(temp_loc,dp) * (real(-log(max(vals(irot_loc), 1.e-30)),dp) - real(lmin,dp)))
                wsum        = wsum + w(irot_loc)
            enddo
            p_mode     = real(w(irot_b) / wsum)
            p_near_loc = 0.
            do irot_loc = 1,nrots
                dstep = min(abs(irot_loc - irot_b), nrots - abs(irot_loc - irot_b))
                if( dstep <= 2 ) p_near_loc = p_near_loc + real(w(irot_loc) / wsum)
            enddo
            nmode_loc = 0
            nnear_loc = 0
            do idraw_loc = 1,NDRAWS
                call pftc%gen_posterior_val(1, 1, [0.,0.], temp_loc, dist_loc, irot_loc)
                dstep = min(abs(irot_loc - irot_b), nrots - abs(irot_loc - irot_b))
                if( irot_loc == irot_b ) nmode_loc = nmode_loc + 1
                if( dstep <= 2 )         nnear_loc = nnear_loc + 1
            enddo
            f_mode     = real(nmode_loc) / real(NDRAWS)
            f_near_loc = real(nnear_loc) / real(NDRAWS)
            write(*,'(A,F10.3,A,F6.3,A,F6.3,A,F6.3,A,F6.3)') '  posterior temperature ', temp_loc, ': mode mass ', p_mode,&
                &' drawn ', f_mode, '; mass within 2 steps ', p_near_loc, ' drawn ', f_near_loc
            call assert_true(abs(f_mode - p_mode) <= 0.1, 'the mode is drawn at the posterior''s mode mass')
            call assert_true(abs(f_near_loc - p_near_loc) <= 0.1, 'two grid steps about the mode are drawn at the posterior''s mass')
        end subroutine draws_against_posterior

    end subroutine test_posterior_inpl_draw

    !> The at-pose loss at a nonzero reference shift is the shift search's objective (gen_objfun_vals, exp(-L))
    !! at the same shift and rotation: the identical-pose null test of prob_post_diag re-scores a refined mode
    !! through it. A seed-bank entry (the particle shifted by a particle-frame shift u) scores at rotation irot
    !! as the reference shifted by shift_ptcl2ref(u, rot(irot)), up to the angular Nyquist bin the
    !! shifted-reference kernel drops; the bank's best-rotation and at-pose losses agree with each other
    !! exactly, and cross swaps the half as for the unshifted kernels.
    subroutine test_bank_and_at_pose_kernels()
        integer, parameter :: BOX = 64, NREFS = 1, NPTCLS = 1
        real,    parameter :: SMPD = 1.0, MSKDIAM = 48.0, LP = 6.0, MSKRAD = 22.0
        real,    parameter :: NYQUIST_ATOL = 5.e-2 ! the dropped bin separates the kernels; a frame error is O(1)
        type(parameters), target :: p
        type(cmdline)       :: cline
        type(polarft_calc)  :: pftc
        type(image)         :: img_even, img_odd
        real,    allocatable, target :: sigma2(:,:)
        real,    allocatable :: vals(:)
        complex, allocatable :: bank(:,:)
        real    :: s(2), u(2), dist_at, dist_b, dist_bx, dist_r, dist_rx
        integer :: pdim(3), kfromto(2), irot, irot_b, irot_bx, irot_r, irot_rx, nrots
        write(*,'(A)') 'test_bank_and_at_pose_kernels'
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nptcls',  real(NPTCLS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline)
        call set_fixed_seed(20261024)
        kfromto = [2, calc_fourier_index(LP, BOX, SMPD)]
        call pftc%new(p, NREFS, [1,NPTCLS], kfromto)
        pdim  = pftc%get_pdim_srch()
        nrots = pftc%get_nrots()
        call img_even%new([BOX,BOX,1], SMPD)
        call img_even%gauran(0.0, 1.0)
        call img_even%fft
        call img_even%bpgau2D(0.0, LP)
        call img_even%ifft
        call img_even%memoize_mask_coords
        call img_even%mask2D_soft(MSKRAD, backgr=0.0)
        call img_odd%new([BOX,BOX,1], SMPD)
        call img_odd%gauran(0.0, 1.0)
        call img_odd%fft
        call img_odd%bpgau2D(0.0, LP)
        call img_odd%ifft
        call img_odd%memoize_mask_coords
        call img_odd%mask2D_soft(MSKRAD, backgr=0.0)
        call img_even%memoize4polarize(pdim)
        call img_even%fft
        call img_odd%fft
        call pftc%polarize_ref_pft(img_even, 1, iseven=.true.,  pdim=pdim, oversamp=.false.)
        call pftc%polarize_ref_pft(img_odd,  1, iseven=.false., pdim=pdim, oversamp=.false.)
        call pftc%polarize_ptcl_pft(img_even, 1, pdim=pdim, oversamp=.false.)
        call pftc%set_eo(1, .true.)
        allocate(sigma2(kfromto(1):kfromto(2),2), source=1.0)
        call pftc%assign_sigma2_noise(sigma2)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        allocate(vals(nrots), bank(pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1))
        ! the at-pose loss at a nonzero reference shift is the shift search's objective, every rotation
        s = [1.3, -0.6]
        call pftc%gen_objfun_vals(1, 1, s, vals)
        do irot = 1,nrots,max(1,nrots/7)
            call pftc%gen_euclid_val_at(1, 1, s, irot, dist_at)
            call assert_real(-log(max(vals(irot), 1.e-30)), dist_at, 1.e-4, &
                &'the at-pose loss at a reference shift is the shift search''s objective at that pose')
        end do
        ! a bank entry: the particle shifted by u scores as the reference shifted by shift_ptcl2ref(u, rot)
        u = [0.9, 1.4]
        call pftc%memoize_shifted_ptcl(1, u, bank)
        call pftc%gen_best_euclid_val_bank(1, 1, bank, dist_b, irot_b)
        call pftc%gen_euclid_val_at_bank(1, 1, bank, irot_b, dist_at)
        call assert_real(dist_b, dist_at, 1.e-6, 'the bank at-pose loss reproduces the bank best value at its rotation')
        call pftc%gen_euclid_val_at_bank(1, 1, bank, 0, dist_at)
        call assert_true(dist_at > 1.e30, 'the bank at-pose loss is huge off the rotation grid')
        call pftc%gen_euclid_val_at(1, 1, shift_ptcl2ref(u, pftc%get_rot(irot_b)), irot_b, dist_r)
        call assert_true(abs(dist_r - dist_b) <= NYQUIST_ATOL, &
            &'the bank loss at its best rotation is the shifted-reference loss at the frame-converted seed, up to the Nyquist bin')
        call pftc%gen_euclid_val_at(1, 1, u, irot_b, dist_rx)
        ! the unconverted seed is the wrong frame unless the rotation is near zero or the shift along the axis
        if( sum(abs(shift_ptcl2ref(u, pftc%get_rot(irot_b)) - u)) > 0.5 ) call assert_true(abs(dist_rx - dist_b) > abs(dist_r - dist_b), &
            &'the frame-converted seed is closer to the bank loss than the unconverted one')
        ! cross swaps the half for the bank kernels as for the unshifted ones
        call pftc%gen_best_euclid_val_bank(1, 1, bank, dist_bx, irot_bx, cross=.true.)
        call pftc%gen_best_euclid_val_cross(1, 1, shift_ptcl2ref(u, pftc%get_rot(irot_bx)), dist_rx, irot_rx)
        call assert_true(abs(dist_bx - dist_rx) <= NYQUIST_ATOL, 'the cross bank loss is the cross shifted-reference loss at the converted seed')
        call assert_true(dist_bx > dist_b + 0.3, 'the other half (an unrelated image) loses to the own half for the shifted particle too')
        call pftc%set_eo(1, .false.)
        call pftc%memoize_sigma2_weights
        call pftc%gen_best_euclid_val_bank(1, 1, bank, dist_r, irot_r, cross=.true.)
        call assert_real(dist_b, dist_r, 1.e-5, 'declared odd: the cross bank loss is the former own one')
        call assert_int(irot_b, irot_r, 'declared odd: at the same rotation')
        call img_even%kill
        call img_odd%kill
        call pftc%kill
        deallocate(sigma2, vals, bank)
    end subroutine test_bank_and_at_pose_kernels

    !> A reprojection model written by a two-reference object (write_ref_pfts) and read back by a
    !! four-reference one at iref_from 1 and 3 (read_ref_pfts): the sections of both halves land at
    !! the offset unchanged, as the pose-posterior diagnostic loads the real halves after the matcher's
    !! references; the plain read of a matching object is the iref_from 1 case.
    subroutine test_ref_pfts_offset_io()
        integer, parameter :: BOX = 64, NREFS = 2, NPTCLS = 1
        real,    parameter :: SMPD = 1.0, MSKDIAM = 48.0, LP = 6.0, MSKRAD = 22.0
        type(parameters), target :: p
        type(cmdline)       :: cline
        type(polarft_calc)  :: pftc_src, pftc_dst, pftc_same
        type(image)         :: imgs(2*NREFS)
        type(string)        :: fname_even, fname_odd
        complex, allocatable :: pft_src(:,:), pft_dst(:,:)
        integer :: pdim(3), kfromto(2), i, iref, ioff
        logical :: iseven
        write(*,'(A)') 'test_ref_pfts_offset_io'
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nptcls',  real(NPTCLS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline)
        call set_fixed_seed(20261023)
        kfromto = [2, calc_fourier_index(LP, BOX, SMPD)]
        call pftc_src%new(p, NREFS, [1,NPTCLS], kfromto)
        pdim = pftc_src%get_pdim_srch()
        do i = 1,2*NREFS
            call imgs(i)%new([BOX,BOX,1], SMPD)
            call imgs(i)%gauran(0.0, 1.0)
            call imgs(i)%fft
            call imgs(i)%bpgau2D(0.0, LP)
            call imgs(i)%ifft
            call imgs(i)%memoize_mask_coords
            call imgs(i)%mask2D_soft(MSKRAD, backgr=0.0)
        end do
        call imgs(1)%memoize4polarize(pdim)
        ! four distinct sections: reference 1 even/odd from images 1/2, reference 2 from images 3/4
        do i = 1,2*NREFS
            call imgs(i)%fft
            iref   = (i + 1) / 2
            iseven = mod(i, 2) == 1
            call pftc_src%polarize_ref_pft(imgs(i), iref, iseven=iseven, pdim=pdim, oversamp=.false.)
        end do
        fname_even = string('test_reproj_model_offset_even.bin')
        fname_odd  = string('test_reproj_model_offset_odd.bin')
        call pftc_src%write_ref_pfts(fname_even, .true.)
        call pftc_src%write_ref_pfts(fname_odd,  .false.)
        call assert_true(file_exists(fname_even) .and. file_exists(fname_odd), 'both halves of the model are written')
        ! the same model twice in a four-reference object: at 1 and at NREFS+1
        call pftc_dst%new(p, 2*NREFS, [1,NPTCLS], kfromto)
        call pftc_dst%read_ref_pfts(fname_even, .true.,  iref_from=1)
        call pftc_dst%read_ref_pfts(fname_odd,  .false., iref_from=1)
        call pftc_dst%read_ref_pfts(fname_even, .true.,  iref_from=NREFS+1)
        call pftc_dst%read_ref_pfts(fname_odd,  .false., iref_from=NREFS+1)
        allocate(pft_src(pftc_src%get_pftsz(),kfromto(1):kfromto(2)), pft_dst(pftc_src%get_pftsz(),kfromto(1):kfromto(2)))
        do ioff = 0,NREFS,NREFS
            do iref = 1,NREFS
                call pftc_src%get_ref_pft(iref, .true., pft_src)
                call pftc_dst%get_ref_pft(iref+ioff, .true., pft_dst)
                call assert_true(all(pft_src == pft_dst), 'even section '//int2str(iref)//' lands at offset '//int2str(ioff)//' unchanged')
                call pftc_src%get_ref_pft(iref, .false., pft_src)
                call pftc_dst%get_ref_pft(iref+ioff, .false., pft_dst)
                call assert_true(all(pft_src == pft_dst), 'odd section '//int2str(iref)//' lands at offset '//int2str(ioff)//' unchanged')
            end do
        end do
        call assert_false(pftc_dst%ref_slots_identical(NREFS+1), 'the loaded halves keep their distinct slots')
        ! the plain read of an object of the model's size
        call pftc_same%new(p, NREFS, [1,NPTCLS], kfromto)
        call pftc_same%read_ref_pfts(fname_even, .true.)
        call pftc_same%read_ref_pfts(fname_odd,  .false.)
        call pftc_src%get_ref_pft(NREFS, .false., pft_src)
        call pftc_same%get_ref_pft(NREFS, .false., pft_dst)
        call assert_true(all(pft_src == pft_dst), 'the plain read of a matching object reproduces the sections')
        call del_file(fname_even)
        call del_file(fname_odd)
        do i = 1,2*NREFS
            call imgs(i)%kill
        end do
        call pftc_src%kill
        call pftc_dst%kill
        call pftc_same%kill
        deallocate(pft_src, pft_dst)
    end subroutine test_ref_pfts_offset_io

    !> The Euclidean loss against the other half's reference (gen_best_euclid_val_cross): with image 1 as
    !! the even reference and an unrelated image 3 as the odd one, an even particle that is image 1 has
    !! loss 0 at rotation 1 against its own half (no CTF, unit sigma2: L = |X - R|^2 / |X|^2) and a loss
    !! near 2 (two uncorrelated images of equal power) against the other; flipping the particle's half
    !! swaps the two, so the switch flips the reference half alone. gen_euclid_val_at reproduces the best
    !! value at its rotation and is huge off the grid; ref_slots_identical tells the fixture's distinct
    !! slots from a copy of one section into both.
    subroutine test_cross_half_euclid()
        use simple_eul_prob_tab_utils, only: posterior_temperature
        integer, parameter :: BOX = 64, NREFS = 1, NPTCLS = 1
        real,    parameter :: SMPD = 1.0, MSKDIAM = 48.0, LP = 6.0, MSKRAD = 22.0
        type(parameters), target :: p
        type(cmdline)       :: cline
        type(polarft_calc)  :: pftc
        type(image)         :: img_even, img_odd
        real, allocatable, target :: sigma2(:,:)
        real    :: dist_own, dist_cross, dist_at, temp
        integer :: pdim(3), kfromto(2), irot_own, irot_cross
        write(*,'(A)') 'test_cross_half_euclid'
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('nptcls',  real(NPTCLS))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline)
        call set_fixed_seed(20261022)
        kfromto = [2, calc_fourier_index(LP, BOX, SMPD)]
        call pftc%new(p, NREFS, [1,NPTCLS], kfromto)
        pdim = pftc%get_pdim_srch()
        call img_even%new([BOX,BOX,1], SMPD)
        call img_even%gauran(0.0, 1.0)
        call img_even%fft
        call img_even%bpgau2D(0.0, LP)
        call img_even%ifft
        call img_even%memoize_mask_coords
        call img_even%mask2D_soft(MSKRAD, backgr=0.0)
        call img_odd%new([BOX,BOX,1], SMPD)
        call img_odd%gauran(0.0, 1.0)
        call img_odd%fft
        call img_odd%bpgau2D(0.0, LP)
        call img_odd%ifft
        call img_odd%memoize_mask_coords
        call img_odd%mask2D_soft(MSKRAD, backgr=0.0)
        call img_even%memoize4polarize(pdim)
        call img_even%fft
        call img_odd%fft
        call pftc%polarize_ref_pft(img_even, 1, iseven=.true.,  pdim=pdim, oversamp=.false.)
        call pftc%polarize_ref_pft(img_odd,  1, iseven=.false., pdim=pdim, oversamp=.false.)
        call pftc%polarize_ptcl_pft(img_even, 1, pdim=pdim, oversamp=.false.)
        call pftc%set_eo(1, .true.)
        allocate(sigma2(kfromto(1):kfromto(2),2), source=1.0)
        call pftc%assign_sigma2_noise(sigma2)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        ! an even particle: its own half is the even reference (itself), the other half the unrelated image
        call pftc%gen_best_objfun_val(1, 1, [0.0,0.0], dist_own, irot_own)
        call pftc%gen_best_euclid_val_cross(1, 1, [0.0,0.0], dist_cross, irot_cross)
        call assert_int(1, irot_own, 'own half: the image against itself peaks at rotation 1')
        call assert_true(dist_own < 1.e-3, 'own half: the image against itself has loss 0')
        call assert_true(dist_cross > 0.5, 'other half: an unrelated reference has a loss of order 1')
        call assert_true(dist_cross > dist_own + 0.5, 'the cross-half loss is the other reference, not the own one')
        ! the loss at a given pose: the best value at its rotation, the other slot with cross, huge off the grid
        call pftc%gen_euclid_val_at(1, 1, [0.0,0.0], irot_own, dist_at)
        call assert_real(dist_own, dist_at, 1.e-5, 'at-pose loss reproduces the best value at its rotation')
        call pftc%gen_euclid_val_at(1, 1, [0.0,0.0], irot_cross, dist_at, cross=.true.)
        call assert_real(dist_cross, dist_at, 1.e-5, 'at-pose loss on the other slot reproduces the cross best value')
        call pftc%gen_euclid_val_at(1, 1, [0.0,0.0], 0, dist_at)
        call assert_true(dist_at > 1.e30, 'at-pose loss is huge off the rotation grid')
        call pftc%gen_euclid_val_at(1, 1, [0.0,0.0], pftc%get_pftsz()+1, dist_at)
        call assert_true(dist_at > 0.5, 'at-pose loss a half turn away is large')
        call assert_false(pftc%ref_slots_identical(1), 'distinct even and odd sections: the slots differ')
        call assert_false(pftc%ref_slots_identical(0), 'slot detector: an invalid reference is not identical')
        ! the same particle declared odd: own and cross swap
        call pftc%set_eo(1, .false.)
        call pftc%memoize_sigma2_weights
        call pftc%gen_best_objfun_val(1, 1, [0.0,0.0], dist_own, irot_own)
        call pftc%gen_best_euclid_val_cross(1, 1, [0.0,0.0], dist_cross, irot_cross)
        call assert_true(dist_own > 0.5, 'declared odd: its own half is now the unrelated reference')
        call assert_int(1, irot_cross, 'declared odd: the cross half is the image itself, rotation 1')
        call assert_true(dist_cross < 1.e-3, 'declared odd: the cross-half loss is 0')
        ! the posterior temperature of this particle, pi wsqsum / (2 pftsz): wsqsum = sum_k k sum_p |X_kp|^2 at
        ! unit sigma2 is the self-term the calculator memoized, so the temperature is positive and finite
        temp = posterior_temperature(pftc%wsqsum_ptcl_now(1), pftc%get_pftsz())
        call assert_true(temp > 0. .and. temp < huge(temp), 'posterior temperature from the memoized self-term is positive and finite')
        ! one section in both slots, as solve3D's LP-set stages leave them: identical, and cross equals own
        call pftc%polarize_ref_pft(img_even, 1, iseven=.false., pdim=pdim, oversamp=.false.)
        call pftc%memoize_refs
        call assert_true(pftc%ref_slots_identical(1), 'one section in both slots: the slots are identical')
        call pftc%set_eo(1, .true.)
        call pftc%memoize_sigma2_weights
        call pftc%gen_best_objfun_val(1, 1, [0.0,0.0], dist_own, irot_own)
        call pftc%gen_best_euclid_val_cross(1, 1, [0.0,0.0], dist_cross, irot_cross)
        call assert_real(dist_own, dist_cross, 1.e-6, 'identical slots: the cross-half loss is the own loss')
        call assert_int(irot_own, irot_cross, 'identical slots: at the same rotation')
        call img_even%kill
        call img_odd%kill
        call pftc%kill
        deallocate(sigma2)
    end subroutine test_cross_half_euclid

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

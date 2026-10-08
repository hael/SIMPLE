!@descr: PCG fractional-update gate of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_frac
implicit none
#include "simple_local_flags.inc"

contains

module subroutine exec_test_pcg_frac_update( self, cline )
    use simple_builder,    only: builder
    use simple_parameters, only: parameters
    class(commander_test_pcg_frac_update), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    type(parameters) :: params
    type(builder)    :: build

    if( .not. cline%defined('trs') ) call cline%set('trs', 5.)
    call cline%set('oritype', 'ptcl3D')
    call cline%set('mkdir', 'no')
    call build%init_params_and_build_general_tbox(cline, params)
    call build%build_strategy3D_tbox(params)
    if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
    call build%spproj%write_segment_inside(params%oritype)
    call validate_rec3D_pcg_fractional_updates(params, build, cline)
    call build%esig%kill
    call build%kill_strategy3D_tbox
    call build%kill_general_tbox
    call simple_end('**** SIMPLE_TEST_PCG_FRAC_UPDATE NORMAL STOP ****', print_simple=.false.)
end subroutine exec_test_pcg_frac_update

!> Project-backed pre-integration gate for accumulator-domain fractional
!! updates. The caller supplies exactly the reconstruction project and
!! geometry; complementary subsets, realized fractions and continuation
!! artifacts are generated deterministically here.
subroutine validate_rec3D_pcg_fractional_updates( params, build, cline )
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters
    use simple_image,              only: image
    use simple_reconstructor_pcg,  only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, PCG_STOP_INDEFINITE, PCG_LAMBDA
    use simple_rec3D_pcg_strategy, only: validate_pcg_common
    use simple_matcher_ptcl_io,    only: prepimgbatch, discrete_read_imgbatch, killimgbatch, prep_rec_observation
    use simple_sigma2_files,       only: load_sigma2_groups
    use simple_math_ft,            only: resample_sigma2
    type(parameters), intent(inout) :: params
    type(builder),    intent(inout) :: build
    class(cmdline),   intent(inout) :: cline
    real, parameter :: RAW_TOL = 2.0e-5, REPLAY_TOL = 2.0e-6
    integer, allocatable :: selected_pinds(:), state_pinds(:), state_subset1(:), state_subset2(:)
    integer, allocatable :: half_pinds(:), subset1(:), subset2(:)
    logical              :: l_sigma_loaded
    integer              :: nselected, state, eo, n_half, nhalves_tested
    real                 :: state_f1, state_f2
    character(len=4)     :: half
    integer              :: funit
    type(string)         :: diag_fname

    call validate_pcg_common(params)
    if( cline%defined('part') ) THROW_HARD('pcg_frac_update cannot run as a worker part')
    nselected = 0
    call build%spproj_field%sample4rec([params%fromp,params%top], nselected, selected_pinds)
    if( nselected < 1 ) THROW_HARD('no active particles selected for PCG fractional-update validation')
    if( params%cc_objfun == OBJFUN_EUCLID )then
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, &
            &build%spproj_field, l_sigma_loaded)
        if( .not. l_sigma_loaded ) THROW_HARD('PCG fractional-update validation requires sigma2 for objfun=euclid')
    endif
    call prepimgbatch(params, build, MAXIMGBATCHSZ)
    diag_fname = 'pcg_fractional_update_validation.txt'
    call fopen(funit, file=diag_fname, status='replace', action='write')
    write(funit,'(A)') 'test=pcg_frac_update'
    write(funit,'(A,I0)') 'selected_particles=', nselected
    write(funit,'(A,I0)') 'box_crop=', params%box_crop
    write(funit,'(A,F12.6)') 'smpd_crop=', params%smpd_crop
    write(funit,'(A,I0)') 'maxits_pcg=', params%maxits_pcg
    write(funit,'(A,ES14.6)') 'rtol=', params%rtol
    nhalves_tested = 0

    do state = 1, params%nstates
        call collect_state(state, selected_pinds, state_pinds)
        if( size(state_pinds) == 0 )then
            deallocate(state_pinds)
            cycle
        endif
        call split_complementary(state_pinds, state_subset1, state_subset2)
        state_f1 = real(size(state_subset1)) / real(size(state_pinds))
        state_f2 = real(size(state_subset2)) / real(size(state_pinds))
        do eo = 0, 1
            call collect_state_half(state, eo, state_pinds, half_pinds)
            call collect_state_half(state, eo, state_subset1, subset1)
            call collect_state_half(state, eo, state_subset2, subset2)
            n_half = size(half_pinds)
            if( n_half == 0 )then
                deallocate(half_pinds, subset1, subset2)
                cycle
            endif
            if( n_half < 2 ) THROW_HARD('PCG fractional test needs two particles per state/half')
            if( size(subset1) < 1 .or. size(subset2) < 1 )then
                THROW_HARD('PCG fractional test needs both state subsets represented in each half')
            endif
            half = merge('odd ', 'even', eo == 1)
            call validate_half(state, eo, trim(half), half_pinds, subset1, subset2, &
                &state_f1, state_f2, funit)
            nhalves_tested = nhalves_tested + 1
            deallocate(half_pinds, subset1, subset2)
        enddo
        deallocate(state_pinds, state_subset1, state_subset2)
    enddo
    if( nhalves_tested < 1 ) THROW_HARD('PCG fractional-update validation found no populated state/half')
    call fclose(funit)
    call killimgbatch(build)
    deallocate(selected_pinds)
    write(logfhandle,'(A,1X,A)') '>>> PCG FRACTIONAL-UPDATE VALIDATION: PASS; DIAGNOSTICS:', diag_fname%to_char()
    call diag_fname%kill

contains

    subroutine collect_state( state_here, pinds, selected )
        integer,              intent(in)  :: state_here, pinds(:)
        integer, allocatable, intent(out) :: selected(:)
        integer :: i, n
        n = 0
        do i = 1, size(pinds)
            if( build%spproj_field%get_state(pinds(i)) == state_here ) n = n + 1
        enddo
        allocate(selected(n))
        n = 0
        do i = 1, size(pinds)
            if( build%spproj_field%get_state(pinds(i)) /= state_here ) cycle
            n = n + 1
            selected(n) = pinds(i)
        enddo
    end subroutine collect_state

    subroutine collect_state_half( state_here, eo_here, pinds, selected )
        integer,              intent(in)  :: state_here, eo_here, pinds(:)
        integer, allocatable, intent(out) :: selected(:)
        integer :: i, n, p
        n = 0
        do i = 1, size(pinds)
            p = pinds(i)
            if( build%spproj_field%get_state(p) == state_here .and. build%spproj_field%get_eo(p) == eo_here ) n = n + 1
        enddo
        allocate(selected(n))
        n = 0
        do i = 1, size(pinds)
            p = pinds(i)
            if( build%spproj_field%get_state(p) /= state_here .or. build%spproj_field%get_eo(p) /= eo_here ) cycle
            n = n + 1
            selected(n) = p
        enddo
    end subroutine collect_state_half

    subroutine split_complementary( pinds, first, second )
        integer,              intent(in)  :: pinds(:)
        integer, allocatable, intent(out) :: first(:), second(:)
        integer :: i, i1, i2
        ! Deliberately use an unequal one-third/two-thirds split. A 50/50
        ! split with the default u=0.5 would make u/f=1 and fail to exercise
        ! the mass-renormalization step at all.
        allocate(first((size(pinds)+2)/3), second(size(pinds)-(size(pinds)+2)/3))
        i1 = 0
        i2 = 0
        do i = 1, size(pinds)
            if( mod(i-1,3) == 0 )then
                i1 = i1 + 1
                first(i1) = pinds(i)
            else
                i2 = i2 + 1
                second(i2) = pinds(i)
            endif
        enddo
    end subroutine split_complementary

    subroutine validate_half( state_here, eo_here, half_here, full_pinds, pinds1, pinds2, f1, f2, unit )
        integer,          intent(in) :: state_here, eo_here, full_pinds(:), pinds1(:), pinds2(:), unit
        character(len=*), intent(in) :: half_here
        real,             intent(in) :: f1, f2
        type(reconstructor_pcg) :: op_full, op_subset, op_sum, op_ref
        type(reconstructor_pcg) :: op_blend, op_oracle, op_replay, op_ensemble
        type(string)            :: f_full, f_sub1, f_sub2, f_chain1, f_chain2
        real, allocatable       :: x_direct(:,:,:), x_replay(:,:,:), x_full(:,:,:), hist(:)
        character(len=256)      :: provenance
        real                    :: u, berr, derr, solve_err, sample_map_err
        integer                 :: nraw, nraw_total, niters

        provenance = 'pcg-frac-update-v1'
        u  = 0.5
        if( params%l_ufrac_trec_defined ) u = params%ufrac_trec
        if( u <= 0.0 .or. u > 1.0 ) THROW_HARD('pcg_frac_update requires 0 < ufrac_trec <= 1')
        f_full   = frac_fname(state_here, half_here, 'full')
        f_sub1   = frac_fname(state_here, half_here, 'subset1')
        f_sub2   = frac_fname(state_here, half_here, 'subset2')
        f_chain1 = frac_fname(state_here, half_here, 'chain1')
        f_chain2 = frac_fname(state_here, half_here, 'chain2')
        write(unit,'(A,I0,A,A)') 'state=', state_here, ' half=', trim(half_here)
        write(unit,'(A,I0,A,I0,A,I0)') 'n_full=', size(full_pinds), &
            &' n_subset1=', size(pinds1), ' n_subset2=', size(pinds2)
        write(unit,'(A,F10.6,A,F10.6,A,F10.6)') 'state_f1=', f1, ' state_f2=', f2, ' u=', u

        call accumulate_raw(full_pinds, op_full)
        call op_full%write_raw_accum(f_full, state_here, eo_here, 1, 1, size(full_pinds), provenance)
        call op_full%kill
        call accumulate_raw(pinds1, op_subset)
        call op_subset%write_raw_accum(f_sub1, state_here, eo_here, 1, 2, size(pinds1), provenance)
        call op_subset%kill
        call accumulate_raw(pinds2, op_subset)
        call op_subset%write_raw_accum(f_sub2, state_here, eo_here, 2, 2, size(pinds2), provenance)
        call op_subset%kill

        call load_weighted(op_blend, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
        call op_blend%scale_raw_accum(1.0/f1)
        call op_blend%write_raw_accum(f_chain1, state_here, eo_here, 1, 1, size(full_pinds), provenance)
        call op_blend%scale_raw_accum(f1)
        call load_weighted(op_ref, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
        call op_blend%compare_raw_accum(op_ref, berr, derr)
        call require_raw('bootstrap working-mass restore', berr, derr, REPLAY_TOL)
        call op_blend%kill
        call op_ref%kill
        call load_weighted(op_replay, f_chain1, state_here, eo_here, 1, 1, provenance, 1.0)
        call load_weighted(op_oracle, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
        call op_oracle%scale_raw_accum(1.0/f1)
        call op_replay%compare_raw_accum(op_oracle, berr, derr)
        call require_raw('bootstrap full-mass chain seed', berr, derr, REPLAY_TOL)
        call op_replay%kill
        call op_oracle%kill

        call new_reduction(op_sum)
        call op_sum%add_raw_accum(f_sub1, state_here, eo_here, 1, 2, provenance, nraw)
        nraw_total = nraw
        call op_sum%add_raw_accum(f_sub2, state_here, eo_here, 2, 2, provenance, nraw)
        nraw_total = nraw_total + nraw
        if( nraw_total /= size(full_pinds) ) THROW_HARD('complementary raw PCG particle counts do not close')
        call load_weighted(op_ref, f_full, state_here, eo_here, 1, 1, provenance, 1.0)
        call op_sum%compare_raw_accum(op_ref, berr, derr)
        call require_raw('complementary additivity', berr, derr, RAW_TOL)
        call op_sum%kill
        call op_ref%kill

        call build_blend(f_sub1, 1, f1, f_full, state_here, eo_here, provenance, u, op_blend)
        call load_weighted(op_oracle, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
        call op_oracle%scale_raw_accum(u/f1)
        call add_weighted(op_oracle, f_full, state_here, eo_here, 1, 1, provenance, 1.0-u)
        call op_blend%compare_raw_accum(op_oracle, berr, derr)
        call require_raw('subset1 u/f blend', berr, derr, REPLAY_TOL)
        call op_blend%write_raw_accum(f_chain1, state_here, eo_here, 1, 1, size(full_pinds), provenance)
        call op_oracle%kill

        call load_weighted(op_replay, f_chain1, state_here, eo_here, 1, 1, provenance, 1.0)
        call op_blend%compare_raw_accum(op_replay, berr, derr)
        call require_raw('continuation write/read', berr, derr, REPLAY_TOL)
        call finalize_and_solve(op_blend, x_direct, hist, niters)
        deallocate(hist)
        call finalize_and_solve(op_replay, x_replay, hist, niters)
        deallocate(hist)
        solve_err = volume_rel_error(x_direct, x_replay)
        if( solve_err > REPLAY_TOL ) THROW_HARD('PCG continuation replay changed the reconstructed solution')
        call op_blend%kill
        call op_replay%kill
        deallocate(x_replay)

        call build_blend(f_sub2, 2, f2, f_full, state_here, eo_here, provenance, u, op_blend)
        call load_weighted(op_oracle, f_sub2, state_here, eo_here, 2, 2, provenance, 1.0)
        call op_oracle%scale_raw_accum(u/f2)
        call add_weighted(op_oracle, f_full, state_here, eo_here, 1, 1, provenance, 1.0-u)
        call op_blend%compare_raw_accum(op_oracle, berr, derr)
        call require_raw('subset2 u/f blend', berr, derr, REPLAY_TOL)
        call op_blend%write_raw_accum(f_chain2, state_here, eo_here, 1, 1, size(full_pinds), provenance)
        call op_blend%kill
        call op_oracle%kill

        call load_weighted(op_ensemble, f_chain1, state_here, eo_here, 1, 1, provenance, f1)
        call add_weighted(op_ensemble, f_chain2, state_here, eo_here, 1, 1, provenance, f2)
        call load_weighted(op_ref, f_full, state_here, eo_here, 1, 1, provenance, 1.0)
        call op_ensemble%compare_raw_accum(op_ref, berr, derr)
        call require_raw('full-mass weighted ensemble', berr, derr, RAW_TOL)
        call op_ensemble%kill
        call finalize_and_solve(op_ref, x_full, hist, niters)
        deallocate(hist)
        sample_map_err = volume_rel_error(x_direct, x_full)
        call op_ref%kill

        write(unit,'(A,ES14.6)') 'continuation_replay_solution_relerr=', solve_err
        write(unit,'(A,ES14.6)') 'single_subset_vs_full_map_relerr_diagnostic_only=', sample_map_err
        write(unit,'(A)') 'status=PASS'
        write(logfhandle,'(A,I0,A,A,A,F7.4,A,ES10.3)') '>>> PCG FRAC | STATE=', state_here, &
            &' | HALF=', trim(half_here), ' | U=', u, ' | REPLAY=', solve_err

        deallocate(x_direct, x_full)
        call del_file(f_full)
        call del_file(f_sub1)
        call del_file(f_sub2)
        call del_file(f_chain1)
        call del_file(f_chain2)
        call f_full%kill
        call f_sub1%kill
        call f_sub2%kill
        call f_chain1%kill
        call f_chain2%kill
    end subroutine validate_half

    subroutine accumulate_raw( pinds, op )
        integer,                 intent(in)    :: pinds(:)
        type(reconstructor_pcg), intent(inout) :: op
        type(oris)           :: selection
        type(ori)            :: orientation
        type(ctfparams)      :: ctfparms
        type(image)          :: obs
        complex, allocatable :: y_batch(:,:,:)
        real, allocatable    :: sig2(:,:)
        integer              :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
        integer              :: i, ii, iptcl, ibatch
        real                 :: shift(2), crop_factor
        call op%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
        call op%set_sym(build%pgrpsyms)
        call op%set_window_sphere(params%msk_crop)
        lims2 = op%get_lims2()
        R = lims2(1,2)
        allocate(sig2(0:R,size(pinds)), source=1.0)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            kfromto = build%esig%get_kfromto()
            do i = 1, size(pinds)
                call resample_sigma2(kfromto(1), kfromto(2), &
                    &build%esig%sigma2_noise(kfromto(1):kfromto(2),pinds(i)), R, 1.0, sig2(0:R,i))
            enddo
        endif
        call selection%new(size(pinds), .true.)
        call orientation%new(.false.)
        crop_factor = real(params%box_crop) / real(params%box)
        do i = 1, size(pinds)
            iptcl = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
            ctfparms%smpd = params%smpd_crop
            shift = build%spproj_field%get_2Dshift(iptcl) * crop_factor
            call orientation%set_ctfvars(ctfparms)
            call orientation%set_shift(shift)
            call selection%set_ori(i, orientation)
        enddo
        call op%prep_particles(selection, use_ctf=.true., sig2=sig2)
        allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
        call op%begin_accum
        call obs%new([params%box_crop,params%box_crop,1], params%smpd_crop)
        do ibatch = 1, size(pinds), MAXIMGBATCHSZ
            batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
            batchsz = batchlims(2)-batchlims(1)+1
            call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
            do ii = 1, batchsz
                ! the backend-neutral observation (normalize, crop, taper), see prep_rec_observation
                call prep_rec_observation(build%imgbatch(ii), build%lmsk, obs, .true.)
                call obs%fft
                y_batch(:,:,ii) = op%extract_native_plane(obs)
            enddo
            call op%accumulate_batch(y_batch, batchsz, batchlims(1))
        enddo
        call obs%kill
        call selection%kill
        call orientation%kill
        deallocate(y_batch, sig2)
    end subroutine accumulate_raw

    subroutine new_reduction( op )
        type(reconstructor_pcg), intent(inout) :: op
        call op%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
        call op%set_window_sphere(params%msk_crop)
        call op%begin_reduction
    end subroutine new_reduction

    subroutine load_weighted( op, fname, state_here, eo_here, part, nparts, provenance, weight )
        type(reconstructor_pcg), intent(inout) :: op
        class(string),           intent(in)    :: fname
        integer,                 intent(in)    :: state_here, eo_here, part, nparts
        character(len=*),        intent(in)    :: provenance
        real,                    intent(in)    :: weight
        call new_reduction(op)
        call add_weighted(op, fname, state_here, eo_here, part, nparts, provenance, weight)
    end subroutine load_weighted

    subroutine add_weighted( op, fname, state_here, eo_here, part, nparts, provenance, weight )
        type(reconstructor_pcg), intent(inout) :: op
        class(string),           intent(in)    :: fname
        integer,                 intent(in)    :: state_here, eo_here, part, nparts
        character(len=*),        intent(in)    :: provenance
        real,                    intent(in)    :: weight
        integer :: nraw
        call op%add_raw_accum_weighted(fname, state_here, eo_here, part, nparts, provenance, weight, nraw)
    end subroutine add_weighted

    subroutine build_blend( f_subset, subset_part, frac, f_full, state_here, eo_here, provenance, u, op )
        class(string),           intent(in)    :: f_subset, f_full
        integer,                 intent(in)    :: subset_part, state_here, eo_here
        real,                    intent(in)    :: frac, u
        character(len=*),        intent(in)    :: provenance
        type(reconstructor_pcg), intent(inout) :: op
        call load_weighted(op, f_subset, state_here, eo_here, subset_part, 2, provenance, u/frac)
        call add_weighted(op, f_full, state_here, eo_here, 1, 1, provenance, 1.0-u)
    end subroutine build_blend

    subroutine finalize_and_solve( op, x, history, niters )
        type(reconstructor_pcg), intent(inout) :: op
        real, allocatable,       intent(out)   :: x(:,:,:), history(:)
        integer,                 intent(out)   :: niters
        type(pcg_solver_outcome) :: outcome
        allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
        call op%end_accum(.true.)
        call op%set_op_mode(PCG_OP_KERNEL)
        call op%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, rel_res_hist=history, &
            &niters=niters, outcome=outcome)
        if( trim(outcome%stop_reason) == PCG_STOP_INDEFINITE ) &
            &THROW_HARD('PCG lost positive-definiteness in the harness solve')
    end subroutine finalize_and_solve

    subroutine require_raw( label, b_err, d_err, tolerance )
        character(len=*), intent(in) :: label
        real,             intent(in) :: b_err, d_err, tolerance
        write(funit,'(A,A)') 'gate=', trim(label)
        write(funit,'(A,ES14.6)') 'B_relerr=', b_err
        write(funit,'(A,ES14.6)') 'D_relerr=', d_err
        if( b_err > tolerance .or. d_err > tolerance )then
            write(logfhandle,'(A,A,A,ES12.4,A,ES12.4)') '>>> PCG FRAC FAIL | ', trim(label), &
                &' | B=', b_err, ' | D=', d_err
            THROW_HARD('PCG fractional-update raw accumulator validation failed')
        endif
    end subroutine require_raw

    pure real function volume_rel_error( lhs, rhs ) result(err)
        real, intent(in) :: lhs(:,:,:), rhs(:,:,:)
        real(dp) :: numerator, denominator
        integer  :: i, j, k
        numerator = 0.0_dp
        denominator = 0.0_dp
        do k = 1, size(lhs,3)
            do j = 1, size(lhs,2)
                do i = 1, size(lhs,1)
                    numerator = numerator + real(lhs(i,j,k)-rhs(i,j,k),dp)**2
                    denominator = denominator + real(rhs(i,j,k),dp)**2
                enddo
            enddo
        enddo
        err = real(sqrt(numerator) / max(1.0_dp, sqrt(denominator)))
    end function volume_rel_error

    function frac_fname( state_here, half_here, tag ) result(fname)
        integer,          intent(in) :: state_here
        character(len=*), intent(in) :: half_here, tag
        type(string) :: fname
        fname = 'pcg_frac_state'//int2str_pad(state_here,2)//'_'//trim(half_here)//'_'//trim(tag)//'.raw'
    end function frac_fname

end subroutine validate_rec3D_pcg_fractional_updates

end submodule simple_commanders_test_highlevel_frac

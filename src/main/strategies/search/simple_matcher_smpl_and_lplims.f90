!@descr: search-space and particle-selection policy routines for matcher workflows
module simple_matcher_smpl_and_lplims
use simple_pftc_srch_api
use simple_builder, only: builder
use simple_solve2D_controller, only: NSAMPLE_DEFAULT_2D
use simple_refine3D_fnames, only: refine3D_fsc_fname
use simple_rnd, only: rnd_key, rnd_scope_key, rnd_draw_key, rnd_run_seed, RND_PROG_REFINE3D, RND_PROG_REFINE2D, &
    &RND_PASS_ALIGN, RND_FAM_SAMPLE, ensure_run_seed, set_run_seed, RND_KEYED_VERSION
use, intrinsic :: iso_fortran_env, only: int64
implicit none

public :: set_bp_range3D, set_bp_range2D
public :: sample_ptcls4update3D, sample_ptcls4fillin, sample_ptcls4missing3D, sample_ptcls4update2D
public :: sample_ptcls4iter3D, master_draws_sample3D, master_draws_sample2D, adopt_run_seed
public :: sample_ptcls4prob3D, sample_ptcls4prob2D
public :: refine2D_requires_full_assignment, all_active_ptcls_2D_assigned, refine2D_blends_carryover
private
#include "simple_local_flags.inc"

contains

    subroutine set_bp_range3D( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        real, allocatable :: resarr(:), fsc_arr(:)
        real              :: fsc0143, fsc05, nu_align_lp
        real              :: mapres(params%nstates)
        integer           :: s, loc, lp_ind, arr_sz, fsc_sz, nyqcrop_ind
        type(string)      :: fsc_fname
        logical           :: fsc_bin_exists(params%nstates), all_fsc_bin_exist, l_nu_align_lp
        l_nu_align_lp  = .false.
        mapres(:)      = huge(0.0)
        build%fsc(:,:) = 0.0
        if( params%l_lpset )then
            ! set Fourier index range
            params%kfromto(2) = calc_fourier_index(params%lp, params%box, params%smpd)
            if( cline%defined('lpstop') )then
                params%kfromto(2) = min(params%kfromto(2),&
                    &calc_fourier_index(params%lpstop, params%box, params%smpd))
            endif
            ! FSC values are read anyway
            do s=1,params%nstates
                fsc_fname = refine3D_fsc_fname(s)
                if( file_exists(fsc_fname) )then
                    fsc_arr = file2rarr(fsc_fname)
                    fsc_sz  = size(build%fsc(s,:))
                    arr_sz  = size(fsc_arr)
                    if( fsc_sz == arr_sz )then
                        build%fsc(s,:) = fsc_arr(:)
                    else if( fsc_sz > arr_sz )then
                        ! padding
                        build%fsc(s,:arr_sz)   = fsc_arr(:)
                        build%fsc(s,arr_sz+1:) = 0.
                    else
                        ! clipping
                        build%fsc(s,:fsc_sz)   = fsc_arr(:fsc_sz)
                    endif
                    deallocate(fsc_arr)
                endif
            enddo
        else
            ! check all fsc_state*.bin exist
            all_fsc_bin_exist = .true.
            fsc_bin_exists    = .false.
            do s=1,params%nstates
                fsc_fname = refine3D_fsc_fname(s)
                fsc_bin_exists( s ) = file_exists(fsc_fname)
                if( build%spproj_field%get_pop(s, 'state') > 0 .and. .not.fsc_bin_exists(s))&
                    & all_fsc_bin_exist = .false.
            enddo
            if(build%spproj%is_virgin_field(params%oritype)) &
                all_fsc_bin_exist = (count(fsc_bin_exists)==params%nstates)
            l_nu_align_lp = try_nu_project_working_lp(params, build, nu_align_lp)
            ! set low-pass Fourier index limit
            if( all_fsc_bin_exist )then
                resarr = build%img%get_res()
                do s=1,params%nstates
                    if( fsc_bin_exists(s) )then
                        fsc_fname = refine3D_fsc_fname(s)
                        fsc_arr = file2rarr(fsc_fname)
                        fsc_sz  = size(build%fsc(s,:))
                        arr_sz  = size(fsc_arr)
                        if( fsc_sz == arr_sz )then
                            build%fsc(s,:) = fsc_arr(:)
                        else if( fsc_sz > arr_sz )then
                            ! padding
                            build%fsc(s,:arr_sz)   = fsc_arr(:)
                            build%fsc(s,arr_sz+1:) = 0.
                        else
                            ! clipping
                            build%fsc(s,:fsc_sz)   = fsc_arr(:fsc_sz)
                        endif
                        deallocate(fsc_arr)
                        call get_resolution(build%fsc(s,:), resarr, fsc05, fsc0143)
                        mapres(s) = fsc0143
                    endif
                end do
            endif
            if( l_nu_align_lp )then
                params%kfromto(2) = calc_fourier_index(nu_align_lp, params%box, params%smpd)
            else if( all_fsc_bin_exist )then
                loc = minloc(mapres, dim=1) ! best resolved
                if( params%nstates == 1 )then
                    lp_ind = get_find_at_crit(build%fsc(1,:), params%lplim_crit, incrreslim=params%l_incrreslim)
                else
                    lp_ind = get_find_at_crit(build%fsc(loc,:), params%lplim_crit)
                endif
                ! interpolation limit is NOT Nyqvist in correlation search
                params%kfromto(2) = calc_fourier_index(resarr(lp_ind), params%box, params%smpd)
            else if( build%spproj_field%isthere(params%fromp,'lp') )then
                params%kfromto(2) = calc_fourier_index(&
                    build%spproj_field%get(params%fromp,'lp'), params%box, params%smpd)
            else
                THROW_HARD('no method for low-pass limit: need fsc file or lp field; set_bp_range3D')
            endif
            ! lpstop overrides any other method for setting the low-pass limit
            if( cline%defined('lpstop') )then
                params%kfromto(2) = min(params%kfromto(2), &
                    calc_fourier_index(params%lpstop, params%box, params%smpd))
            endif
            ! re-set the low-pass limit
            params%lp = calc_lowpass_lim(params%kfromto(2), params%box, params%smpd)
        endif
        ! making sure the resolution limit does not exceed limits of box_crop
        nyqcrop_ind       = calc_fourier_index(2.*params%smpd_crop, params%box_crop, params%smpd_crop)
        params%kfromto(2) = min(params%kfromto(2), nyqcrop_ind)
        params%lp         = calc_lowpass_lim(params%kfromto(2), params%box, params%smpd)
        if( l_nu_align_lp )then
            write(logfhandle,'(A,F8.3,A)') '>>> NU filter matching low-pass limit: ', params%lp, ' A'
            if( params%l_nonuniform_lpset )then
                params%l_lpset = .true.
                call cline%set('lp', params%lp)
                write(logfhandle,'(A,F8.3,A)') &
                    &'>>> NU filter promoted matching low-pass to command line: ', params%lp, ' A'
            endif
        endif
        ! update low-pass limit in project
        call build%spproj_field%set_all2single('lp',params%lp)
    end subroutine set_bp_range3D

    logical function try_nu_project_working_lp( params, build, align_lp ) result( l_use )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        real,              intent(out)   :: align_lp
        logical :: l_fresh_start
        align_lp = 0.
        l_use    = .false.
        if( .not. params%l_nonuniform ) return
        l_fresh_start = params%which_iter <= params%startit .and. trim(params%continue).ne.'yes'
        if( l_fresh_start )then
            ! A fresh start normally cannot trust an lp left in the project by
            ! some earlier run. It can when the NU-filtered references that lp
            ! was derived from are on disk: the same assembly wrote both, so
            ! the handoff describes the references this iteration will match
            ! (refine3D_auto's startup bootstrap). Without this the first
            ! iteration falls back to the FSC band while every later one uses
            ! the evidence-derived band -- the discontinuity the bootstrap
            ! exists to remove.
            if( .not. file_exists(add2fbody(params%vols_even(1), MRC_EXT, NUFILT_SUFFIX)) ) return
        endif
        if( .not. build%spproj_field%isthere(params%fromp, 'lp') ) return
        align_lp = build%spproj_field%get(params%fromp, 'lp')
        l_use = align_lp > TINY
        if( align_lp <= TINY ) l_use = .false.
    end function try_nu_project_working_lp

    subroutine set_bp_range2D( params, build, cline, which_iter, frac_srch_space )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: which_iter
        real,              intent(in)    :: frac_srch_space
        real    :: lplim
        integer :: lpstart_find
        params%kfromto(1) = max(2,calc_fourier_index(params%hp, params%box, params%smpd))
        if( params%l_lpset )then
            lplim = params%lp
            params%kfromto(2) = calc_fourier_index(lplim, params%box_crop, params%smpd_crop)
        else
            if( trim(params%stream2d).eq.'yes' )then
                if( file_exists(params%frcs) )then
                    lplim = build%clsfrcs%estimate_lp_for_align()
                else
                    lplim = params%lplims2D(3)
                endif
                if( cline%defined('lpstop') ) lplim = max(lplim, params%lpstop)
            else
                if( file_exists(params%frcs) .and. which_iter >= LPLIM1ITERBOUND )then
                    lplim = build%clsfrcs%estimate_lp_for_align()
                else
                    if( which_iter < LPLIM1ITERBOUND )then
                        lplim = params%lplims2D(1)
                    else if( frac_srch_space >= FRAC_SH_LIM .and. which_iter > LPLIM3ITERBOUND )then
                        lplim = params%lplims2D(3)
                    else
                        lplim = params%lplims2D(2)
                    endif
                endif
            endif
            params%kfromto(2) = calc_fourier_index(lplim, params%box_crop, params%smpd_crop)
            ! to avoid pathological cases, fall-back on lpstart
            lpstart_find = calc_fourier_index(params%lpstart, params%box_crop, params%smpd_crop)
            if( lpstart_find > params%kfromto(2) ) params%kfromto(2) = lpstart_find
            lplim     = calc_lowpass_lim(params%kfromto(2), params%box_crop, params%smpd_crop)
            params%lp = lplim
        endif
        ! update low-pas limit in project
        call build%spproj_field%set_all2single('lp',lplim)
    end subroutine set_bp_range2D

    !> Fractional-update selection by balance: none takes the global lowest-updatecnt tiers, class and cavg
    !! the nested quota over the sampling units of CLASS_SAMPLING_FILE. Under cohort_sampling a stage
    !! draws at its first iteration and rescores that cohort at the others.
    subroutine sample_ptcls4update3D( params, build, pfromto, l_incr_sampl, nptcls2update, pinds, allow_empty )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        logical, optional,    intent(in)    :: allow_empty  !< a distributed partition may sample nothing
        type(class_sample),   allocatable   :: clssmp(:)
        type(rnd_key),        allocatable   :: key
        type(string) :: fname
        if( .not. params%l_update_frac )then
            ! we sample all state > 0
            call build%spproj_field%sample4update_all(pfromto, nptcls2update, pinds, l_incr_sampl)
            return
        endif
        ! the probabilistic routes draw the sample on the keyed generator (an absent key: the intrinsic stream)
        if( params%l_prob_align_mode ) key = sample_rnd_key(params, RND_PROG_REFINE3D)
        if( params%l_cohort_sampling .and. params%which_iter /= params%startit )then
            call build%spproj_field%sample4update_rescore(pfromto, nptcls2update, pinds, allow_empty=allow_empty)
            return
        endif
        select case(trim(params%balance))
            case('class', 'cavg')
                fname = CLASS_SAMPLING_FILE
                if( .not. file_exists(fname) )then
                    THROW_HARD('File for class-balanced sampling in fractional update: '//CLASS_SAMPLING_FILE//' does not exist!')
                endif
                call read_class_samples(clssmp, fname)
                if( params%l_frac_best )then
                    call build%spproj_field%sample4update_class(clssmp, pfromto, params%update_frac, &
                        &nptcls2update, pinds, l_incr_sampl, params%l_greedy_smpl, frac_best=params%frac_best, &
                        &allow_empty=allow_empty, key=key)
                else
                    call build%spproj_field%sample4update_class(clssmp, pfromto, params%update_frac, &
                        &nptcls2update, pinds, l_incr_sampl, params%l_greedy_smpl, &
                        &allow_empty=allow_empty, key=key)
                endif
                call deallocate_class_samples(clssmp)
                call fname%kill
            case DEFAULT
                call build%spproj_field%sample4update_cnt(pfromto, params%update_frac, &
                    &nptcls2update, pinds, l_incr_sampl, allow_empty=allow_empty, key=key)
        end select
    end subroutine sample_ptcls4update3D

    !> The project owner's seed at the start of a refine3D or refine2D run, before any part is scheduled: a run
    !! that continues (continue=yes) reuses the seed the project records (a SIMPLE_SEED set to another value is
    !! inconsistent restart metadata and stops the run); any other run keeps a set SIMPLE_SEED or draws one and
    !! exports it. The seed and the generator version are recorded in the project's info segment and written.
    subroutine adopt_run_seed( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer(int64) :: seed
        integer        :: version
        logical        :: found
        character(len=96) :: msg
        found = .false.
        if( trim(params%continue) == 'yes' ) call build%spproj%get_rnd_seed(seed, version, found)
        if( found )then
            if( version /= RND_KEYED_VERSION )then
                write(msg,'(a,i0,a,i0)') 'the project records generator version ', version, ', this build is ', RND_KEYED_VERSION
                THROW_HARD(msg)
            endif
            call set_run_seed(seed)
        else
            seed = ensure_run_seed()
        endif
        call build%spproj%set_rnd_seed(seed, RND_KEYED_VERSION)
        call build%spproj%write_segment_inside('projinfo', params%projfile)
    end subroutine adopt_run_seed

    !> the key of the probabilistic routes' sample draw (family SAMPLE, the project owner's particle 0): the
    !! run's seed, the program, the discrete alignment pass and the iteration, so that the sample depends on
    !! the run and the field alone
    function sample_rnd_key( params, program ) result( key )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: program
        type(rnd_key) :: key
        key = rnd_draw_key(rnd_scope_key(rnd_run_seed(), program, RND_PASS_ALIGN, params%which_iter), 0, RND_FAM_SAMPLE)
    end function sample_rnd_key

    !> The particle sample of one refine3D search iteration over pfromto, advancing sampled and updatecnt:
    !! the missing-update pass, the periodic fill-in, or the fractional-update policy of sample_ptcls4update3D
    subroutine sample_ptcls4iter3D( params, build, which_iter, pfromto, nptcls2update, pinds, allow_empty )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: which_iter
        integer,              intent(in)    :: pfromto(2)
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        logical, optional,    intent(in)    :: allow_empty  !< a distributed partition may sample nothing
        if( params%l_update_missing )then
            call sample_ptcls4missing3D(build, pfromto, .true., nptcls2update, pinds)
        else if( params%l_fillin .and. mod(which_iter,5) == 0 )then
            call sample_ptcls4fillin(params, build, pfromto, .true., nptcls2update, pinds, allow_empty=allow_empty)
        else
            call sample_ptcls4update3D(params, build, pfromto, .true., nptcls2update, pinds, allow_empty=allow_empty)
        endif
    end subroutine sample_ptcls4iter3D

    !> Whether the distributed refine3D master draws the iteration's sample once, over the whole project,
    !! for the workers to reproduce. Each worker would otherwise draw its own sample over the whole project
    !! with its own random seed and keep the part in its range, so the union is not one sample. The
    !! probabilistic modes always do (their parts run the pipeline on the master's sample, whatever the
    !! update fraction, so the field's sample stamps must be the master's); otherwise under a fractional
    !! update, except that the polish pass reproduces the discrete pass it follows, the missing-update
    !! pass and full sampling are deterministic, and the eval and sigma passes keep their own selection.
    !! Master and workers evaluate this on the same command line, so they agree.
    logical function master_draws_sample3D( params ) result( l_master )
        class(parameters), intent(in) :: params
        l_master = .false.
        if( params%l_cont_polish .or. params%l_update_missing ) return
        if( params%l_prob_align_mode )then
            l_master = .true.
            return
        endif
        if( .not. params%l_update_frac ) return
        select case(trim(params%refine))
            case('eval','sigma')
                return
        end select
        l_master = .true.
    end function master_draws_sample3D

    !> Whether the distributed refine2D master draws the iteration's sample once, over the whole project, for
    !! the workers to reproduce, as master_draws_sample3D does: otherwise each worker draws within its own
    !! particles, and the sample depends on the number of parts. The probabilistic modes always do (their parts
    !! run the pipeline on the master's sample); the stream pool keeps its own sampling. Master and workers
    !! evaluate this on the same command line, so they agree.
    logical function master_draws_sample2D( params ) result( l_master )
        class(parameters), intent(in) :: params
        l_master = .false.
        if( trim(params%stream2d) == 'yes' ) return
        if( params%l_prob_align_mode )then
            l_master = .true.
            return
        endif
        if( .not. params%l_update_frac ) return
        l_master = .true.
    end function master_draws_sample2D

    !> The probabilistic routes' sample draw of one refine3D iteration over the whole project, by the project
    !! owner (the shared-memory matcher, or the distributed master whose parts reproduce it): the sampling
    !! history is reset at the run's first iteration only, the fill-in sample every fifth iteration under
    !! fillin, else the fractional-update draw; then the disk-set readout, whose inputs are the owner's
    subroutine sample_ptcls4prob3D( params, build, nptcls2update, pinds )
        use simple_eul_prob_tab_utils, only: print_disk_set
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(out)   :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        if( params%l_fillin .and. mod(params%which_iter,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls2update, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls2update, pinds)
        endif
        call print_disk_set(build%spproj_field, params%l_centring, params%trs, pinds, nptcls2update)
    end subroutine sample_ptcls4prob3D

    !> The same for refine2D: the history reset at the run's first iteration, the draw of sample_ptcls4update2D
    !! (l_updatefrac the strategy's rule: false on a full update), the disk-set readout
    subroutine sample_ptcls4prob2D( params, build, l_updatefrac, nptcls2update, pinds )
        use simple_eul_prob_tab_utils, only: print_disk_set
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        logical,              intent(in)    :: l_updatefrac
        integer,              intent(out)   :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        call sample_ptcls4update2D(params, build, [1,params%nptcls], l_updatefrac, nptcls2update, pinds)
        call print_disk_set(build%spproj_field, params%l_centring, params%trs, pinds, nptcls2update)
    end subroutine sample_ptcls4prob2D

    subroutine sample_ptcls4fillin( params, build, pfromto, l_incr_sampl, nptcls2update, pinds, allow_empty )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        logical, optional,    intent(in)    :: allow_empty  !< a distributed partition may sample nothing
        type(rnd_key), allocatable :: key
        if( params%l_prob_align_mode ) key = sample_rnd_key(params, RND_PROG_REFINE3D)
        call build%spproj_field%sample4update_fillin(pfromto, params%update_frac, nptcls2update, pinds, l_incr_sampl, &
            &allow_empty=allow_empty, key=key)
    end subroutine sample_ptcls4fillin

    subroutine sample_ptcls4missing3D( build, pfromto, l_incr_sampl, nptcls2update, pinds )
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        call build%spproj_field%sample4update_missing(pfromto, nptcls2update, pinds, l_incr_sampl)
    end subroutine sample_ptcls4missing3D

    subroutine sample_ptcls4update2D( params, build, pfromto, l_updatefrac, nptcls, pinds )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        logical,              intent(in)    :: l_updatefrac
        integer,              intent(in)    :: pfromto(2)
        integer,              intent(inout) :: nptcls
        integer, allocatable, intent(inout) :: pinds(:)
        type(rnd_key), allocatable :: key
        logical :: l_has_been_sampled
        l_has_been_sampled = build%spproj_field%has_been_sampled()
        if( params%l_prob_align_mode ) key = sample_rnd_key(params, RND_PROG_REFINE2D)
        if( l_updatefrac )then
            ! solve2D controller policy:
            ! startit==1  -> sticky sampled subset stage
            ! startit>1   -> stochastic resampling biased toward low updatecnt
            ! assign_all waits for every particle to be assigned, which a sticky subset never reaches,
            ! so a fresh start with assign_all resamples by update count as well
            if( params%startit == 1 .and. .not. params%l_assign_all )then
                if( l_has_been_sampled )then
                    call build%spproj_field%sample4update_reprod(pfromto, nptcls, pinds)
                else
                    call build%spproj_field%sample4update_rnd(pfromto, params%update_frac, nptcls, pinds, .true., key=key)
                    call build%spproj_field%set_updatecnt(1, pinds)
                endif
            else
                call build%spproj_field%sample4update_cnt(pfromto, params%update_frac, nptcls, pinds, .true., key=key)
            endif
        else
            call build%spproj_field%sample4update_all(pfromto, nptcls, pinds, .true.)
        endif
    end subroutine sample_ptcls4update2D

    !> Whether a refine2D iteration blends carried class sums (fractional class-average
    !! restoration): a fractional update past a fresh start, or, in streaming, any iteration
    !! after the first with an update fraction below 0.99. Shared by the matcher and by the
    !! masters that check the carried set before the iteration.
    logical function refine2D_blends_carryover( params, which_iter ) result( l_blend )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: which_iter
        if( trim(params%stream2d) == 'yes' )then
            l_blend = (which_iter > 1) .and. (params%update_frac < 0.99)
        else
            l_blend = params%l_update_frac .and. (params%startit /= 1)
        endif
    end function refine2D_blends_carryover

    !> assign_all=yes (under fractional updates): refine2D converges only once every active particle has a class.
    !! A coverage guard only; fillin keeps the sampling it has in 3D.
    logical function refine2D_requires_full_assignment( params ) result( l_required )
        class(parameters), intent(in) :: params
        l_required = params%l_assign_all
    end function refine2D_requires_full_assignment

    logical function all_active_ptcls_2D_assigned( os, pfromto, n_missing ) result( l_assigned )
        class(oris),       intent(in)  :: os
        integer,           intent(in)  :: pfromto(2)
        integer, optional, intent(out) :: n_missing
        integer :: iptcl, missing
        missing = 0
        do iptcl = pfromto(1), pfromto(2)
            if( os%get_state(iptcl) == 0 ) cycle
            if( os%get_updatecnt(iptcl) <= 0 .or. os%get_class(iptcl) <= 0 .or. nint(os%get(iptcl,'inpl')) <= 0 )then
                missing = missing + 1
            endif
        enddo
        if( present(n_missing) ) n_missing = missing
        l_assigned = missing == 0
    end function all_active_ptcls_2D_assigned

end module simple_matcher_smpl_and_lplims

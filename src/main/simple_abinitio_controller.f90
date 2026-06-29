submodule(simple_abinitio_utils) simple_abinitio_controller
implicit none
#include "simple_local_flags.inc"

! Output naming
character(len=*), parameter :: REC_FBODY               = 'rec_final_state'

! Stage schedule and search-space sizes
integer,          parameter :: NSTAGES                 = 8
integer,          parameter :: NSTAGES_INI3D           = 4    ! # of ini3D stages used for initialization
integer,          parameter :: NSTAGES_INI3D_MAX       = 7
integer,          parameter :: MAXITS(8)               = [20,20,17,15,12,12,12,25]
integer,          parameter :: NSPACE(8)               = [1000,1000,1000,1000,2500,2500,5000,5000]
integer,          parameter :: NSPACE_SUB              = 126
integer,          parameter :: NSPACE_SUB_BASE         = 2500

! Stage transition policy
integer,          parameter :: TURNED_OFF              = NSTAGES + 1 ! value for stage-based policies to indicate "turned off" 
integer,          parameter :: GAUREF_LAST_STAGE       = 2           ! stop gaussian filtering after early stages
integer,          parameter :: SYMSRCH_STAGE           = 3           ! search symmetry axis
integer,          parameter :: PROB_REFINE_STAGE       = 3           ! prob refinement stages 3-5
integer,          parameter :: TRAILREC_STAGE_SINGLE   = 5           ! first stage where trail_rec behavior changes
integer,          parameter :: STOCH_SAMPL_STAGE       = 5           ! switch from greedy to stochastic sampling
integer,          parameter :: STOCH_SAMPL_STAGE_INDEP = 4           ! independent multi-state needs earlier stochastic coverage
integer,          parameter :: NU_FILTER_STAGE         = 6           ! switch on staged NU filtering
integer,          parameter :: PROB_NEIGH_REFINE_STAGE = 6           ! prob_neigh refinement stages 6-8
integer,          parameter :: NSTAGES_INDEPENDENT     = PROB_NEIGH_REFINE_STAGE - 1
integer,          parameter :: GOLD_STD_STAGE          = TURNED_OFF  ! gold-standard doesn't work for abinitio 3D 
integer,          parameter :: AUTOMSK_STAGE           = NSTAGES     ! switch on automasking
integer,          parameter :: TRAILREC_STAGE_MULTI    = NSTAGES
integer,          parameter :: HET_DOCKED_STAGE        = 6           ! split after stage 5; stage 6 stabilizes split states
character(len=*), parameter :: PROB_NEIGH_MODE_EARLY   = 'shc'
character(len=*), parameter :: PROB_NEIGH_MODE_LATE    = 'state'
character(len=*), parameter :: PROB_NEIGH_MODE_MULTI   = 'sum'
character(len=*), parameter :: PROB_NEIGH_MODE_DOCKED  = 'geom'

! Filtering and low-pass defaults
real,             parameter :: LPSTOP_BOUNDS(2)        = [4.5,6.0]
real,             parameter :: LPSTART_BOUNDS(2)       = [10.,20.]
real,             parameter :: CENLP_DEFAULT           = 30.
real,             parameter :: LPSYMSRCH_LB            = 12.
real,             parameter :: LPSTART_INI3D           = 20.  ! default lpstart for abinitio3D_cavgs/cavgs_ini
real,             parameter :: LPSTOP_INI3D            = 8.   ! default lpstop for abinitio3D_cavgs/cavgs_ini
real,             parameter :: LPSTOP_INDEPENDENT      = 6.   ! conservative default for independent multi-state abinitio3D

! Sampling and update defaults
real,             parameter :: UPDATE_FRAC_MAX            = 0.9  ! ensures fractional update remains on
real,             parameter :: FULL_SAMPLE_SWITCH_FRAC    = 0.9  ! force all-active sampling once nsample/active reaches this fraction
integer,          parameter :: NSAMPLE_ABINITIO3D_DEFAULT = 10000

type :: refine3D_stage_cfg
    type(string) :: ml_reg, fillin, conical_fsc
    type(string) :: refine, trail_rec, pgrp, balance, filt_mode, automsk, nu_refine, greedy_sampling, prob_neigh_mode
    integer :: iter, inspace, inspace_sub, imaxits
    real    :: trs, frac_best, overlap, fracsrch
    real    :: snr_noise_reg, gaufreq, update_frac_dyn
end type refine3D_stage_cfg

contains

    module function abinitio_rec_fbody() result(fbody)
        character(len=15) :: fbody
        fbody = REC_FBODY
    end function abinitio_rec_fbody

    module function abinitio_lpstop_bounds() result(bounds)
        real :: bounds(2)
        bounds = LPSTOP_BOUNDS
    end function abinitio_lpstop_bounds

    module function abinitio_lpstart_bounds() result(bounds)
        real :: bounds(2)
        bounds = LPSTART_BOUNDS
    end function abinitio_lpstart_bounds

    module function abinitio_cenlp_default() result(cenlp)
        real :: cenlp
        cenlp = CENLP_DEFAULT
    end function abinitio_cenlp_default

    module function abinitio_lpsymsrch_lb() result(lp)
        real :: lp
        lp = LPSYMSRCH_LB
    end function abinitio_lpsymsrch_lb

    module function abinitio_update_frac_max() result(frac)
        real :: frac
        frac = UPDATE_FRAC_MAX
    end function abinitio_update_frac_max

    module function abinitio_full_sample_switch_frac() result(frac)
        real :: frac
        frac = FULL_SAMPLE_SWITCH_FRAC
    end function abinitio_full_sample_switch_frac

    module function abinitio_lpstart_ini3D() result(lp)
        real :: lp
        lp = LPSTART_INI3D
    end function abinitio_lpstart_ini3D

    module function abinitio_lpstop_ini3D() result(lp)
        real :: lp
        lp = LPSTOP_INI3D
    end function abinitio_lpstop_ini3D

    module function abinitio_independent_lpstop_default() result(lp)
        real :: lp
        lp = LPSTOP_INDEPENDENT
    end function abinitio_independent_lpstop_default

    module function abinitio_nstages() result(nstages_out)
        integer :: nstages_out
        nstages_out = NSTAGES
    end function abinitio_nstages

    module function abinitio_independent_nstages_default() result(nstages_out)
        integer :: nstages_out
        nstages_out = NSTAGES_INDEPENDENT
    end function abinitio_independent_nstages_default

    module function abinitio_nstages_ini3D() result(nstages_out)
        integer :: nstages_out
        nstages_out = NSTAGES_INI3D
    end function abinitio_nstages_ini3D

    module function abinitio_nstages_ini3D_max() result(nstages_out)
        integer :: nstages_out
        nstages_out = NSTAGES_INI3D_MAX
    end function abinitio_nstages_ini3D_max

    module function abinitio_symsrch_stage() result(istage)
        integer :: istage
        istage = SYMSRCH_STAGE
    end function abinitio_symsrch_stage

    module function abinitio_het_docked_stage() result(istage)
        integer :: istage
        istage = HET_DOCKED_STAGE
    end function abinitio_het_docked_stage

    module function abinitio_stoch_sampl_stage(params) result(istage)
        class(parameters), intent(in) :: params
        integer :: istage
        istage = STOCH_SAMPL_STAGE
        if( trim(params%multivol_mode).eq.'independent' ) istage = STOCH_SAMPL_STAGE_INDEP
    end function abinitio_stoch_sampl_stage

    module function abinitio_nsample_default() result(nsample)
        integer :: nsample
        nsample = NSAMPLE_ABINITIO3D_DEFAULT
    end function abinitio_nsample_default

    integer function active_refine3D_nstages() result(nstages_active)
        nstages_active = nstages_refine3D
        if( nstages_active <= 0 ) nstages_active = NSTAGES
        nstages_active = min(nstages_active, size(NSPACE), size(MAXITS))
        if( allocated(lpinfo) ) nstages_active = min(nstages_active, size(lpinfo))
    end function active_refine3D_nstages

    module procedure set_cline_refine3D
        type(refine3D_stage_cfg) :: cfg
        call build_refine3D_stage_cfg( cfg, params, istage, l_cavgs )
        call emit_refine3D_stage_cfg( cfg, params, istage, l_cavgs, l_refine3D_lp_override )
    end procedure set_cline_refine3D

    subroutine build_refine3D_stage_cfg( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        call init_refine3D_iteration( cfg )
        call set_refine3D_update_policy( cfg, params, istage )
        call set_refine3D_symmetry_policy( cfg, params, istage )
        call set_refine3D_mode_policy( cfg, params, istage )
        call set_refine3D_balance_policy( cfg )
        call set_refine3D_gauref_policy( cfg, params, istage, l_cavgs )
        call set_refine3D_trailrec_policy( cfg, params, istage )
        call set_refine3D_filtering_policy( cfg, params, istage, l_cavgs )
        call set_refine3D_automsk_policy( cfg, params, istage, l_cavgs )
        call set_refine3D_stage_controls( cfg, params, istage )
        call apply_refine3D_search_overrides( cfg )
    end subroutine build_refine3D_stage_cfg

    subroutine init_refine3D_iteration( cfg )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        cfg%iter = 0
        if( cline_refine3D%defined('endit') )then
            cfg%iter = cline_refine3D%get_iarg('endit')
        endif
        cfg%iter = cfg%iter + 1
        cfg%inspace_sub = 0
    end subroutine init_refine3D_iteration

    subroutine set_refine3D_update_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        real :: update_frac_stage
        if( docked_split_stage(params, istage) .or. force_full_sampling_mode(params) )then
            cfg%fillin          = 'no'
            cfg%update_frac_dyn = 1.0
            return
        endif
        update_frac_stage = update_frac
        if( istage == active_refine3D_nstages() )then
            cfg%fillin = 'yes'
            if( params%nstates > 1 ) cfg%fillin = 'no'
            cfg%update_frac_dyn = update_frac_stage
        else
            cfg%fillin = 'no'
            cfg%update_frac_dyn = update_frac_stage
        endif
        cfg%update_frac_dyn = min(UPDATE_FRAC_MAX, cfg%update_frac_dyn)
    end subroutine set_refine3D_update_policy

    subroutine set_refine3D_symmetry_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%pgrp = trim(params%pgrp)
        if( l_srch4symaxis )then
            if( istage <= SYMSRCH_STAGE )then
                cfg%pgrp = trim(params%pgrp_start)
            endif
        endif
    end subroutine set_refine3D_symmetry_policy

    subroutine set_refine3D_mode_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%prob_neigh_mode = ''
        if( l_refine3D_mode_override )then
            cfg%refine = refine3D_mode_override
            if( cfg%refine.eq.'prob_neigh' )then
                cfg%prob_neigh_mode = trim(params%prob_neigh_mode)
                if( params%nstates == 1 .and. cfg%prob_neigh_mode.eq.'sum' ) cfg%prob_neigh_mode = PROB_NEIGH_MODE_LATE
            endif
        else if( istage <  PROB_REFINE_STAGE )then
            cfg%refine           = 'prob_neigh'
            cfg%prob_neigh_mode  = PROB_NEIGH_MODE_EARLY
        else if( istage < PROB_NEIGH_REFINE_STAGE )then
            cfg%refine = 'prob'
        else
            cfg%refine           = 'prob_neigh'
            if( params%nstates > 1 )then
                if( trim(params%multivol_mode).eq.'docked' )then
                    cfg%prob_neigh_mode  = PROB_NEIGH_MODE_DOCKED
                else
                    cfg%prob_neigh_mode  = PROB_NEIGH_MODE_MULTI
                endif
            else
                cfg%prob_neigh_mode  = PROB_NEIGH_MODE_LATE
            endif
        endif
        if( docked_split_stage(params, istage) )then
            cfg%refine           = 'prob_state'
            cfg%prob_neigh_mode  = ''
        endif
    end subroutine set_refine3D_mode_policy

    subroutine set_refine3D_balance_policy( cfg )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        cfg%balance = 'yes'
    end subroutine set_refine3D_balance_policy

    subroutine set_refine3D_gauref_policy( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        cfg%gaufreq = -1.
        if( l_cavgs .or. istage <= GAUREF_LAST_STAGE )then
            cfg%gaufreq = lpinfo(istage)%lp
        endif
    end subroutine set_refine3D_gauref_policy

    subroutine set_refine3D_trailrec_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%trail_rec = 'no'
        if( force_full_sampling_mode(params) ) return
        select case(trim(params%multivol_mode))
            case('single')
                if( istage >= TRAILREC_STAGE_SINGLE ) cfg%trail_rec = 'yes'
            case('independent')
                if( istage >= TRAILREC_STAGE_MULTI  ) cfg%trail_rec = 'yes'
            case('docked')
                if( istage >= TRAILREC_STAGE_SINGLE .and. istage /= params%split_stage )then
                    cfg%trail_rec = 'yes'
                endif
            case default
                cfg%trail_rec = 'no'
        end select
    end subroutine set_refine3D_trailrec_policy

    subroutine set_refine3D_filtering_policy( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        cfg%filt_mode  = 'none'
        cfg%nu_refine  = 'no'
        if( l_cavgs ) return
        if( l_nonuniform .and. &
            &(istage >= NU_FILTER_STAGE .or. (l_state_continue_mode .and. istage >= TRAILREC_STAGE_SINGLE)) )then
            cfg%filt_mode = trim(params%filt_mode)
            if( cfg%filt_mode.eq.'nonuniform' .and. &
                &(istage < GOLD_STD_STAGE .or. params%nstates > 1) ) cfg%filt_mode = 'nonuniform_lpset'
            if( cfg%filt_mode.eq.'nonuniform_lpset' .and. &
                &params%nstates == 1 .and. istage >= GOLD_STD_STAGE ) cfg%filt_mode = 'nonuniform'
            cfg%nu_refine = 'no'
        endif
    end subroutine set_refine3D_filtering_policy

    subroutine set_refine3D_automsk_policy( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        cfg%automsk = 'no'
        if( l_cavgs ) return
        if( istage >= AUTOMSK_STAGE .and. l_automsk ) cfg%automsk = trim(params%automsk)
    end subroutine set_refine3D_automsk_policy

    subroutine set_refine3D_stage_controls( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        integer :: stoch_stage
        stoch_stage = abinitio_stoch_sampl_stage(params)
        cfg%inspace = NSPACE(istage)
        cfg%greedy_sampling = 'yes'
        select case(istage)
            case(1,2)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = 0.
                cfg%ml_reg        = 'no'
                cfg%conical_fsc   = 'no'
                cfg%frac_best     = 1.0
                cfg%overlap       = 0.99
                cfg%fracsrch      = 99.
                cfg%snr_noise_reg = 2.0
            case(3,4,5,6)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = lpinfo(istage)%trslim
                cfg%ml_reg        = 'yes'
                cfg%conical_fsc   = params%conical_fsc
                cfg%frac_best     = 1.0
                if( trim(params%multivol_mode).eq.'independent' .and. istage >= stoch_stage )then
                    cfg%greedy_sampling = 'no'
                else if( istage >= stoch_stage )then
                    cfg%frac_best = 0.5
                endif
                if( istage > SYMSRCH_STAGE )then
                    cfg%overlap  = 0.9 ! early stopping
                    cfg%fracsrch = 90. ! early stopping
                else
                    cfg%overlap  = 0.99
                    cfg%fracsrch = 99.
                endif
                cfg%snr_noise_reg = 4.0
            case(7,8)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = lpinfo(istage)%trslim
                cfg%ml_reg        = 'yes'
                cfg%conical_fsc   = params%conical_fsc
                if( trim(params%multivol_mode).eq.'independent' )then
                    cfg%frac_best       = 1.0
                    cfg%greedy_sampling = 'no'
                else if( params%nstates > 1 )then
                    cfg%frac_best = 0.98
                else
                    cfg%frac_best = 0.85
                endif
                cfg%overlap       = 0.95 ! early stopping
                cfg%fracsrch      = 90.  ! doesn't affect later stages
                cfg%snr_noise_reg = 6.0
            case default
                THROW_HARD('Invalid istage index in set_refine3D_stage_controls')
        end select
    end subroutine set_refine3D_stage_controls

    subroutine apply_refine3D_search_overrides( cfg )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        select case(cfg%refine%to_char())
            case('prob_neigh')
                select case(cfg%prob_neigh_mode%to_char())
                    case('shc','snhc')
                        cfg%inspace_sub = 0
                    case DEFAULT
                        ! Keep neighborhood granularity roughly constant when
                        ! late stages increase nspace (e.g. 2500 -> 5000).
                        cfg%inspace_sub = max(1, nint(real(NSPACE_SUB) * real(cfg%inspace) / real(NSPACE_SUB_BASE)))
                end select
        end select
    end subroutine apply_refine3D_search_overrides

    subroutine emit_refine3D_stage_cfg( cfg, params, istage, l_cavgs, l_cmdline_lp_override )
        type(refine3D_stage_cfg), intent(in) :: cfg
        class(parameters),        intent(in) :: params
        integer,                  intent(in) :: istage
        logical,                  intent(in) :: l_cavgs
        logical,                  intent(in) :: l_cmdline_lp_override
        character(len=STDLEN) :: ptcl_src_eff
        real :: lp_eff
        logical :: l_den_src
        logical :: l_full_update_stage
        l_full_update_stage = docked_split_stage(params, istage) .or. force_full_sampling_mode(params)
        ptcl_src_eff        = stage_ptcl_src(cfg, params)
        l_den_src           = trim(ptcl_src_eff) == 'den'
        lp_eff              = stage_matching_lp(cfg, params, istage, l_cmdline_lp_override)
        call cline_refine3D%set('prg',                     'refine3D')
        if( l_cavgs )then
            call cline_refine3D%set('envfsc',              'no')
        else if( params%nstates == 1 .and. istage >= GOLD_STD_STAGE )then
            call cline_refine3D%set('envfsc',              'yes')
        else
            call cline_refine3D%set('envfsc',              'no')
        endif
        if( l_full_update_stage )then
            call cline_refine3D%delete('update_frac')
            call cline_refine3D%delete('fillin')
        else if( istage == active_refine3D_nstages() )then
            call cline_refine3D%set('update_frac',        cfg%update_frac_dyn)
            call cline_refine3D%set('fillin',             cfg%fillin)
        else
            call cline_refine3D%set('update_frac',        cfg%update_frac_dyn)
            call cline_refine3D%delete('fillin')
        endif
        if( .not. l_cavgs )then
            if( l_full_update_stage )then
                call cline_refine3D%delete('nsample')
            else
                call cline_refine3D%set('nsample', params%nsample)
            endif
        endif
        call cline_refine3D%set('box_crop',               lpinfo(istage)%box_crop)
        call cline_refine3D%set('startit',                cfg%iter)
        call cline_refine3D%set('which_iter',             cfg%iter)
        call cline_refine3D%set('pgrp',                   cfg%pgrp)
        call cline_refine3D%set('refine',                 cfg%refine)
        if( cfg%refine.eq.'prob_neigh' )then
            call cline_refine3D%set('prob_neigh_mode',    cfg%prob_neigh_mode)
        else
            call cline_refine3D%delete('prob_neigh_mode')
        endif
        call cline_refine3D%set('balance',                cfg%balance)
        call cline_refine3D%set('trail_rec',              cfg%trail_rec)
        call cline_refine3D%set('filt_mode',              cfg%filt_mode)
        call cline_refine3D%set('ptcl_src',               ptcl_src_eff)
        call cline_refine3D%set('nu_refine',              cfg%nu_refine)
        call cline_refine3D%delete('lpstart')
        call cline_refine3D%delete('lpstop')
        call cline_refine3D%set('automsk',                cfg%automsk)
        if( params%nstates == 1 .and. istage >= GOLD_STD_STAGE )then
            ! Past this point, NU filtering promotes the selected matching
            ! bandwidth; the controller no longer injects schedule LP.
            call cline_refine3D%delete('lp')
        else
            call cline_refine3D%set('lp',                 lp_eff)
        endif
        call cline_refine3D%set('nspace',                 cfg%inspace)
        if( cfg%inspace_sub > 0 )then
            call cline_refine3D%set('nspace_sub',         cfg%inspace_sub)
        else
            call cline_refine3D%delete('nspace_sub')
        endif
        call cline_refine3D%set('maxits',                 cfg%imaxits)
        call cline_refine3D%set('trs',                    cfg%trs)
        ! if( l_den_src )then
        !     call cline_refine3D%set('ml_reg',             'no')
        ! else
            call cline_refine3D%set('ml_reg',             cfg%ml_reg)
        ! endif
        call cline_refine3D%set('conical_fsc',            cfg%conical_fsc)
        call cline_refine3D%set('greedy_sampling',        cfg%greedy_sampling)
        call cline_refine3D%set('frac_best',              cfg%frac_best)
        call cline_refine3D%set('overlap',                cfg%overlap)
        call cline_refine3D%set('fracsrch',               cfg%fracsrch)
        if( l_cavgs )then
            call cline_refine3D%set('snr_noise_reg',      cfg%snr_noise_reg)
            call cline_refine3D%delete('update_frac')
        else
            call cline_refine3D%delete('snr_noise_reg')
        endif
        ! if( l_den_src )then
        !     call cline_refine3D%set('gauref',             'yes')
        !     call cline_refine3D%set('gaufreq',            lp_eff)
        ! else 
        if( cfg%gaufreq > 0. )then
            call cline_refine3D%set('gauref',             'yes')
            call cline_refine3D%set('gaufreq',            cfg%gaufreq)
        else
            call cline_refine3D%delete('gauref')
            call cline_refine3D%delete('gaufreq')
        endif
    end subroutine emit_refine3D_stage_cfg

    logical function docked_split_stage( params, istage )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        docked_split_stage = trim(params%multivol_mode).eq.'docked' .and. istage == params%split_stage
    end function docked_split_stage

    logical function force_full_sampling_mode( params ) result( l_force_full )
        class(parameters), intent(in) :: params
        real :: sample_frac
        l_force_full = .false.
        if( nptcls_eff <= 0 ) return
        if( params%nsample <= 0 ) return
        sample_frac  = real(params%nsample) / real(nptcls_eff)
        l_force_full = sample_frac > abinitio_full_sample_switch_frac()
    end function force_full_sampling_mode

    real function stage_matching_lp( cfg, params, istage, l_cmdline_lp_override ) result( lp )
        type(refine3D_stage_cfg), intent(in) :: cfg
        class(parameters),        intent(in) :: params
        integer,                  intent(in) :: istage
        logical,                  intent(in) :: l_cmdline_lp_override
        lp = lpinfo(istage)%lp
        if( l_cmdline_lp_override .and. cfg%ml_reg.eq.'yes' ) lp = params%lp
    end function stage_matching_lp

    character(len=STDLEN) function stage_ptcl_src( cfg, params ) result( ptcl_src )
        type(refine3D_stage_cfg), intent(in) :: cfg
        class(parameters),        intent(in) :: params
        ptcl_src = trim(params%ptcl_src)
        ! if( ptcl_src == 'den' )then
        !     select case(cfg%filt_mode%to_char())
        !         case('nonuniform','nonuniform_lpset')
        !             ptcl_src = 'raw'
        !     end select
        ! endif
    end function stage_ptcl_src

end submodule simple_abinitio_controller

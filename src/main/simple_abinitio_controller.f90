submodule(simple_abinitio_utils) simple_abinitio_controller
implicit none
#include "simple_local_flags.inc"

real,    parameter :: UPDATE_FRAC_MIN                 = 0.1                  ! 10% of the particles updated each iteration
integer, parameter :: NSPACE(8)                       = [500,500,1000,1000,2500,2500,5000,5000]
integer, parameter :: SHC_REFINE_STAGE                = 1                    ! shc-style refinement stages 1-4
integer, parameter :: PROB_REFINE_STAGE               = 5                    ! prob refinement stages 5-6
integer, parameter :: PROB_NEIGH_REFINE_STAGE         = 7                    ! prob_neigh refinement stages 7-8
integer, parameter :: STOCH_SAMPL_STAGE               = PROB_REFINE_STAGE    ! we switch from greedy to stochastic balanced class sampling when prob is switched on
integer, parameter :: LPAUTO_STAGE                    = PROB_REFINE_STAGE + 1 ! we switch on automatic low-pass limit estimation after one prob stage, but only if lp_auto is switched on by user
integer, parameter :: REFINE3D_ROUTE_STD              = 1
integer, parameter :: REFINE3D_ROUTE_CAVGS            = 2
integer, parameter :: REFINE3D_ROUTE_POLAR            = 3
integer, parameter :: REFINE3D_ROUTE_POLAR_CAVGS      = 4
integer, parameter :: NEIGH_NSPACES(2)                = [126,5000]

type :: refine3D_stage_cfg
    type(string) :: ml_reg, fillin
    type(string) :: refine, trail_rec, pgrp, balance, lp_auto, automsk
    integer :: iter, inspace, inspace_sub, imaxits, nsample_dyn, ipftsz
    real    :: trs, frac_best, overlap, fracsrch, lpstart, lpstop
    real    :: snr_noise_reg, gaufreq, update_frac_dyn
end type refine3D_stage_cfg

contains

    integer function classify_refine3D_route( l_cavgs ) result(route)
        logical, intent(in) :: l_cavgs
        if( l_polar )then
            if( l_cavgs )then
                route = REFINE3D_ROUTE_POLAR_CAVGS
            else
                route = REFINE3D_ROUTE_POLAR
            endif
        else
            if( l_cavgs )then
                route = REFINE3D_ROUTE_CAVGS
            else
                route = REFINE3D_ROUTE_STD
            endif
        endif
    end function classify_refine3D_route

    module procedure set_cline_refine3D
        type(refine3D_stage_cfg) :: cfg
        integer :: route
        route = classify_refine3D_route(l_cavgs)
        call build_refine3D_stage_cfg( cfg, params, istage, l_cavgs, route )
        call emit_refine3D_stage_cfg( cfg, params, istage, l_cavgs, route )
    end procedure set_cline_refine3D

    subroutine build_refine3D_stage_cfg( cfg, params, istage, l_cavgs, route )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        integer,                  intent(in)    :: route
        call init_refine3D_iteration( cfg )
        call set_refine3D_update_policy( cfg, params, istage )
        call set_refine3D_symmetry_policy( cfg, params, istage )
        call set_refine3D_mode_policy( cfg, params, istage, route )
        call set_refine3D_balance_policy( cfg )
        call set_refine3D_gauref_policy( cfg, params, istage )
        call set_refine3D_trailrec_policy( cfg, params, istage )
        call set_refine3D_lp_auto_policy( cfg, params, istage )
        call set_refine3D_automsk_policy( cfg, istage, route )
        call set_refine3D_stage_controls( cfg, params, istage )
        call apply_refine3D_route_overrides( cfg, params, istage, route )
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
        if( istage == NSTAGES )then
            cfg%fillin = 'yes'
            if( params%nstates > 1 ) cfg%fillin = 'no'
            if( l_nsample_stop_given )then
                cfg%update_frac_dyn = real(nsample_minmax(2)) / real(nptcls_eff)
            else if( l_nsample_given )then
                cfg%update_frac_dyn = update_frac
            else
                cfg%nsample_dyn     = nint(UPDATE_FRAC_MIN * real(nptcls_eff) / real(params%nstates))
                cfg%nsample_dyn     = max(NSAMPLE_MINMAX_DEFAULT(1), min(NSAMPLE_MINMAX_DEFAULT(2), cfg%nsample_dyn))
                cfg%update_frac_dyn = real(cfg%nsample_dyn * params%nstates) / real(nptcls_eff)
            endif
        else
            cfg%fillin = 'no'
            cfg%update_frac_dyn = calc_update_frac_dyn(nptcls_eff, params%nstates, nsample_minmax, cfg%iter, maxits_dyn)
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

    subroutine set_refine3D_mode_policy( cfg, params, istage, route )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        integer,                  intent(in)    :: route
        if( istage <  PROB_REFINE_STAGE )then
            cfg%refine = 'shc_smpl'
        else if( istage < PROB_NEIGH_REFINE_STAGE )then
            cfg%refine = 'prob'
        else
            cfg%refine = 'prob_neigh'
        endif
        if( trim(params%multivol_mode).eq.'input_oris_fixed' )then
            cfg%refine = 'prob_state'
        endif  
    end subroutine set_refine3D_mode_policy

    subroutine set_refine3D_balance_policy( cfg )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        cfg%balance = 'yes'
    end subroutine set_refine3D_balance_policy

    subroutine set_refine3D_gauref_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%gaufreq = -1.
        if( istage <= GAUREF_LAST_STAGE )then
            cfg%gaufreq = lpinfo(istage)%lp
        endif
    end subroutine set_refine3D_gauref_policy

    subroutine set_refine3D_trailrec_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%trail_rec = 'no'
        select case(trim(params%multivol_mode))
            case('single')
                if( istage >= TRAILREC_STAGE_SINGLE ) cfg%trail_rec = 'yes'
            case('independent')
                if( istage >= TRAILREC_STAGE_MULTI  ) cfg%trail_rec = 'yes'
            case('docked')
                if( istage == NSTAGES )then
                    cfg%trail_rec = 'no'
                else if( istage >= TRAILREC_STAGE_SINGLE )then
                    cfg%trail_rec = 'yes'
                endif
            case('input_oris_fixed')
                cfg%trail_rec = 'no'
            case('input_oris_start')
                cfg%trail_rec = 'no'
            case default
                cfg%trail_rec = 'no'
        end select
    end subroutine set_refine3D_trailrec_policy

    subroutine set_refine3D_lp_auto_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%lp_auto = 'no'
        cfg%lpstart = 0.
        cfg%lpstop  = 0.
        if( istage >= LPAUTO_STAGE .and. l_lpauto )then
            cfg%lp_auto = trim(params%lp_auto)
            cfg%lpstart = lpinfo(istage - 1)%lp
            if( istage == NSTAGES )then
                cfg%lpstop = lpinfo(istage)%smpd_crop * 2.
            else
                cfg%lpstop = lpinfo(istage + 1)%lp
            endif
        endif
    end subroutine set_refine3D_lp_auto_policy

    subroutine set_refine3D_automsk_policy( cfg, istage, route )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        integer,                  intent(in)    :: istage
        integer,                  intent(in)    :: route
        cfg%automsk = 'no'
        if( route == REFINE3D_ROUTE_STD )then
            if( istage >= AUTOMSK_STAGE .and. l_automsk ) cfg%automsk = 'yes'
        endif
    end subroutine set_refine3D_automsk_policy

    subroutine set_refine3D_stage_controls( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%inspace = NSPACE(istage)
        select case(istage)
            case(1,2)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = 0.
                cfg%ml_reg        = 'no'
                cfg%frac_best     = 1.0
                cfg%overlap       = 0.99
                cfg%fracsrch      = 99.
                cfg%snr_noise_reg = 2.0
            case(3,4,5,6)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = lpinfo(istage)%trslim
                cfg%ml_reg        = 'yes'
                if( istage >= STOCH_SAMPL_STAGE )then
                    cfg%frac_best = 0.5
                else
                    cfg%frac_best = 1.0
                endif
                if( istage > SYMSRCH_STAGE )then
                    cfg%overlap  = 0.9
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
                if( params%nstates > 1 )then
                    cfg%frac_best = 0.98
                else
                    cfg%frac_best = 0.85
                endif
                cfg%overlap       = 0.95 ! early stopping
                cfg%fracsrch      = 90.
                cfg%snr_noise_reg = 6.0
            case default
                THROW_HARD('Invalid istage index in set_refine3D_stage_controls')
        end select
    end subroutine set_refine3D_stage_controls

    subroutine apply_refine3D_route_overrides( cfg, params, istage, route )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        integer,                  intent(in)    :: route
        cfg%ipftsz = 0
        select case(route)
            case(REFINE3D_ROUTE_STD, REFINE3D_ROUTE_CAVGS)
                ! nothing to override for standard and cavgs routes
            case(REFINE3D_ROUTE_POLAR, REFINE3D_ROUTE_POLAR_CAVGS)
                if( trim(params%gauref).eq.'no' ) cfg%gaufreq = -1.
                if( cfg%ml_reg.eq.'yes' )         cfg%gaufreq = -1.
                if( cfg%trail_rec=='yes' )then
                    if( trim(params%multivol_mode)=='docked' )then
                        cfg%ipftsz = magic_pftsz(params%msk, params%box, lpinfo(NSTAGES-1)%box_crop)
                    else
                        cfg%ipftsz = magic_pftsz(params%msk, params%box, lpinfo(NSTAGES)%box_crop)
                    endif
                    cfg%inspace = NSPACE(NSTAGES)
                endif
        end select
        select case(cfg%refine%to_char())
            case('prob_neigh')
                cfg%inspace_sub = NEIGH_NSPACES(1)
                cfg%inspace     = NEIGH_NSPACES(2)
        end select
        ! This is a temporary fix because nspace can only be increased after
        ! a cartesian reconstruction that cannot happen with trail_rec=yes for now
        if( params%l_polar )then
            select case(cfg%refine%to_char())
                case('prob_neigh')
                    cfg%inspace = NSPACE(NSTAGES)
            end select
        endif
    end subroutine apply_refine3D_route_overrides

    subroutine emit_refine3D_stage_cfg( cfg, params, istage, l_cavgs, route )
        type(refine3D_stage_cfg), intent(in) :: cfg
        class(parameters),        intent(in) :: params
        integer,                  intent(in) :: istage
        logical,                  intent(in) :: l_cavgs
        integer,                  intent(in) :: route
        call cline_refine3D%set('prg',                     'refine3D')
        if( l_update_frac_dyn .or. istage == NSTAGES )then
            call cline_refine3D%set('update_frac',        cfg%update_frac_dyn)
            call cline_refine3D%set('fillin',             cfg%fillin)
        else
            call cline_refine3D%set('update_frac',        update_frac)
            call cline_refine3D%delete('fillin')
        endif
        call cline_refine3D%set('lp',                     lpinfo(istage)%lp)
        call cline_refine3D%set('smpd_crop',              lpinfo(istage)%smpd_crop)
        call cline_refine3D%set('box_crop',               lpinfo(istage)%box_crop)
        call cline_refine3D%set('startit',                cfg%iter)
        call cline_refine3D%set('which_iter',             cfg%iter)
        call cline_refine3D%set('pgrp',                   cfg%pgrp)
        call cline_refine3D%set('refine',                 cfg%refine)
        call cline_refine3D%set('balance',                cfg%balance)
        call cline_refine3D%set('trail_rec',              cfg%trail_rec)
        call cline_refine3D%set('lp_auto',                cfg%lp_auto)
        if( cfg%lp_auto.eq.'yes' )then
            call cline_refine3D%set('lpstart',            cfg%lpstart)
            call cline_refine3D%set('lpstop',             cfg%lpstop)
        else
            call cline_refine3D%delete('lpstart')
            call cline_refine3D%delete('lpstop')
        endif
        call cline_refine3D%set('automsk',                cfg%automsk)
        call cline_refine3D%set('nspace',                 cfg%inspace)
        if( cfg%inspace_sub > 0 )then
            call cline_refine3D%set('nspace_sub',         cfg%inspace_sub)
        else
            call cline_refine3D%delete('nspace_sub')
        endif
        call cline_refine3D%set('maxits',                 cfg%imaxits)
        call cline_refine3D%set('trs',                    cfg%trs)
        call cline_refine3D%set('ml_reg',                 cfg%ml_reg)
        call cline_refine3D%set('frac_best',              cfg%frac_best)
        call cline_refine3D%set('overlap',                cfg%overlap)
        call cline_refine3D%set('fracsrch',               cfg%fracsrch)
        if( l_cavgs )then
            call cline_refine3D%set('snr_noise_reg',      cfg%snr_noise_reg)
            call cline_refine3D%delete('update_frac')
        else
            call cline_refine3D%delete('snr_noise_reg')
        endif
        if( cfg%gaufreq > 0. )then
            call cline_refine3D%set('gauref',             'yes')
            call cline_refine3D%set('gaufreq',            cfg%gaufreq)
        else
            call cline_refine3D%delete('gauref')
            call cline_refine3D%delete('gaufreq')
        endif
        call cline_refine3D%delete('ipftsz')
        select case(route)
            case(REFINE3D_ROUTE_POLAR, REFINE3D_ROUTE_POLAR_CAVGS)
                call cline_refine3D%set('center_type',    'params')
                if( cfg%ipftsz > 0 )then
                    call cline_refine3D%set('pftsz',      cfg%ipftsz)
                endif
        end select
    end subroutine emit_refine3D_stage_cfg

end submodule simple_abinitio_controller

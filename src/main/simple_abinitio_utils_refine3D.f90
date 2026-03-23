submodule(simple_abinitio_utils) simple_abinitio_utils_refine3D
implicit none
#include "simple_local_flags.inc"

integer, parameter :: REFINE3D_ROUTE_STD              = 1
integer, parameter :: REFINE3D_ROUTE_CAVGS            = 2
integer, parameter :: REFINE3D_ROUTE_POLAR            = 3
integer, parameter :: REFINE3D_ROUTE_POLAR_CAVGS      = 4
integer, parameter :: REFINE3D_ROUTE_TREE_STD         = 5
integer, parameter :: REFINE3D_ROUTE_TREE_CAVGS       = 6
integer, parameter :: REFINE3D_ROUTE_TREE_POLAR       = 7
integer, parameter :: REFINE3D_ROUTE_TREE_POLAR_CAVGS = 8
integer, parameter :: SHC_PTREE_NSPACE                = 2000
integer, parameter :: SHC_PTREE_NSPACE_SUB            = 50
integer, parameter :: PROB_NEIGH_NSPACE               = 4000
integer, parameter :: PROB_NEIGH_NSPACE_SUB           = 100

type :: refine3D_stage_cfg
    type(string) :: ml_reg, fillin, cavgw, neigh_type
    type(string) :: refine, icm, trail_rec, pgrp, balance, lp_auto, automsk
    integer :: iphase, iter, inspace, inspace_sub, imaxits, nsample_dyn, ipftsz
    real    :: trs, frac_best, overlap, fracsrch, lpstart, lpstop
    real    :: snr_noise_reg, gaufreq, update_frac_dyn
end type refine3D_stage_cfg

contains

    integer function classify_refine3D_route( l_cavgs ) result(route)
        logical, intent(in) :: l_cavgs
        if( l_polar )then
            if( l_cavgs )then
                if( l_tree )then
                    route = REFINE3D_ROUTE_TREE_POLAR_CAVGS
                else
                    route = REFINE3D_ROUTE_POLAR_CAVGS
                endif
            else
                if( l_tree )then
                    route = REFINE3D_ROUTE_TREE_POLAR
                else
                    route = REFINE3D_ROUTE_POLAR
                endif
            endif
        else
            if( l_cavgs )then
                if( l_tree )then
                    route = REFINE3D_ROUTE_TREE_CAVGS
                else
                    route = REFINE3D_ROUTE_CAVGS
                endif
            else
                if( l_tree )then
                    route = REFINE3D_ROUTE_TREE_STD
                else
                    route = REFINE3D_ROUTE_STD
                endif
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
        call set_refine3D_icm_policy( cfg, params, istage )
        call set_refine3D_balance_policy( cfg )
        call set_refine3D_gauref_policy( cfg, params, istage )
        call set_refine3D_trailrec_policy( cfg, params, istage )
        call set_refine3D_lp_auto_policy( cfg, params, istage )
        call set_refine3D_automsk_policy( cfg, istage, route )
        call set_refine3D_cavgw_policy( cfg, params, istage, l_cavgs )
        call set_refine3D_phase_index( cfg, istage )
        call set_refine3D_phase_controls( cfg, params, istage, route )
        call apply_refine3D_route_overrides( cfg, params, istage, route )
        call normalize_refine3D_stage_cfg( cfg )
        call cap_refine3D_nspace( cfg, params )
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
        cfg%refine  = 'shc_smpl'
        cfg%neigh_type = ''
        if( is_tree_route(route) ) cfg%refine = 'shc_ptree'
        if( istage >= PROBREFINE_STAGE )then
            if( is_tree_route(route) )then
                ! cfg%refine     = 'prob_neigh'
                ! cfg%neigh_type = 'subspace_srch'
                cfg%refine     = 'ptree'
                cfg%neigh_type = 'subspace_srch'
            else
                cfg%refine  = 'prob'
            endif
        endif
        if( trim(params%multivol_mode).eq.'input_oris_fixed' )then
            cfg%refine = 'prob_state'
        endif
    end subroutine set_refine3D_mode_policy

    subroutine set_refine3D_icm_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%icm = 'no'
        if( istage >= ICM_STAGE ) cfg%icm = 'yes'
        if( cline_refine3D%defined('icm') )then
            if( trim(params%icm).eq.'no' ) cfg%icm = 'no'
        endif
    end subroutine set_refine3D_icm_policy

    subroutine set_refine3D_balance_policy( cfg )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        cfg%balance = 'yes'
    end subroutine set_refine3D_balance_policy

    subroutine set_refine3D_gauref_policy( cfg, params, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        cfg%gaufreq = -1.
        if( istage <= params%gauref_last_stage )then
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
        if( route == REFINE3D_ROUTE_STD .or. route == REFINE3D_ROUTE_TREE_STD )then
            if( istage >= AUTOMSK_STAGE .and. l_automsk ) cfg%automsk = 'yes'
        endif
    end subroutine set_refine3D_automsk_policy

    subroutine set_refine3D_cavgw_policy( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        cfg%cavgw = 'no'
        if( l_cavgs )then
            if( (trim(params%cavgw).eq.'yes') .and. (istage>=CAVGWEIGHTS_STAGE) )then
                cfg%cavgw = 'yes'
            endif
        endif
    end subroutine set_refine3D_cavgw_policy

    subroutine set_refine3D_phase_index( cfg, istage )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        integer,                  intent(in)    :: istage
        cfg%iphase = 0
        if(      istage <= PHASES(1) )then
            cfg%iphase = 1
        else if( istage <= PHASES(2) )then
            cfg%iphase = 2
        else if( istage <= PHASES(3) )then
            cfg%iphase = 3
        else
            THROW_HARD('Invalid istage index')
        endif
    end subroutine set_refine3D_phase_index

    subroutine set_refine3D_phase_controls( cfg, params, istage, route )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        integer,                  intent(in)    :: route
        select case(cfg%iphase)
            case(1)
                cfg%inspace       = NSPACE(1)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = 0.
                cfg%ml_reg        = 'no'
                cfg%frac_best     = 1.0
                cfg%overlap       = 0.99
                cfg%fracsrch      = 99.
                cfg%snr_noise_reg = 2.0
            case(2)
                cfg%inspace       = NSPACE(2)
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
                    cfg%fracsrch = 90.
                else
                    cfg%overlap  = 0.99
                    cfg%fracsrch = 99.
                endif
                cfg%snr_noise_reg = 4.0
            case(3)
                cfg%inspace       = NSPACE(3)
                cfg%imaxits       = MAXITS(istage)
                cfg%trs           = lpinfo(istage)%trslim
                cfg%ml_reg        = 'yes'
                if( params%nstates > 1 )then
                    cfg%frac_best = 0.98
                else
                    cfg%frac_best = 0.85
                endif
                cfg%overlap       = 0.9
                cfg%fracsrch      = 90.
                cfg%snr_noise_reg = 6.0
        end select
        if( is_tree_route(route) )then
            if( cfg%refine == 'shc_ptree' )then
                cfg%inspace     = SHC_PTREE_NSPACE
                cfg%inspace_sub = SHC_PTREE_NSPACE_SUB
            else if( cfg%refine == 'prob_neigh' )then
                cfg%inspace     = PROB_NEIGH_NSPACE
                cfg%inspace_sub = PROB_NEIGH_NSPACE_SUB
            endif
        endif
    end subroutine set_refine3D_phase_controls

    subroutine apply_refine3D_route_overrides( cfg, params, istage, route )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        integer,                  intent(in)    :: route
        cfg%ipftsz = 0
        select case(route)
            case(REFINE3D_ROUTE_STD, REFINE3D_ROUTE_CAVGS)
                ! nothing to override for standard and cavgs routes
            case(REFINE3D_ROUTE_POLAR, REFINE3D_ROUTE_POLAR_CAVGS, REFINE3D_ROUTE_TREE_POLAR, REFINE3D_ROUTE_TREE_POLAR_CAVGS)
                cfg%icm = 'no'
                if( trim(params%gauref).eq.'no' ) cfg%gaufreq = -1.
                if( cfg%ml_reg.eq.'yes' )         cfg%gaufreq = -1.
                if( cfg%trail_rec=='yes' )then
                    if( trim(params%multivol_mode)=='docked' )then
                        cfg%ipftsz = magic_pftsz(params%msk, params%box, lpinfo(NSTAGES-1)%box_crop)
                    else
                        cfg%ipftsz = magic_pftsz(params%msk, params%box, lpinfo(NSTAGES)%box_crop)
                    endif
                    cfg%inspace = NSPACE(3)
                endif
        end select
    end subroutine apply_refine3D_route_overrides

    subroutine normalize_refine3D_stage_cfg( cfg )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        if( cfg%icm.eq.'yes' ) cfg%ml_reg = 'no'
    end subroutine normalize_refine3D_stage_cfg

    subroutine cap_refine3D_nspace( cfg, params )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        if( cfg%refine == 'shc_ptree' ) return
        if( cline_refine3D%defined('nspace_max') )then
            cfg%inspace = min(cfg%inspace, params%nspace_max)
        endif
    end subroutine cap_refine3D_nspace

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
        if( cfg%neigh_type%strlen_trim() > 0 )then
            call cline_refine3D%set('neigh_type',         cfg%neigh_type)
        else
            call cline_refine3D%delete('neigh_type')
        endif
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
        if( l_cavgs )then
            call cline_refine3D%set('cavgw',              cfg%cavgw)
        endif
        call cline_refine3D%set('nspace',                 cfg%inspace)
        if( cfg%inspace_sub > 0 )then
            call cline_refine3D%set('nspace_sub',         cfg%inspace_sub)
        else
            call cline_refine3D%delete('nspace_sub')
        endif
        call cline_refine3D%set('maxits',                 cfg%imaxits)
        call cline_refine3D%set('trs',                    cfg%trs)
        call cline_refine3D%set('ml_reg',                 cfg%ml_reg)
        call cline_refine3D%set('icm',                    cfg%icm)
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
            case(REFINE3D_ROUTE_POLAR, REFINE3D_ROUTE_POLAR_CAVGS, REFINE3D_ROUTE_TREE_POLAR, REFINE3D_ROUTE_TREE_POLAR_CAVGS)
                call cline_refine3D%set('center_type',    'params')
                if( cfg%ipftsz > 0 )then
                    call cline_refine3D%set('pftsz',      cfg%ipftsz)
                endif
        end select
    end subroutine emit_refine3D_stage_cfg

    logical function is_tree_route( route )
        integer, intent(in) :: route
        is_tree_route = route == REFINE3D_ROUTE_TREE_STD .or. route == REFINE3D_ROUTE_TREE_CAVGS .or. &
                        route == REFINE3D_ROUTE_TREE_POLAR .or. route == REFINE3D_ROUTE_TREE_POLAR_CAVGS
    end function is_tree_route

end submodule simple_abinitio_utils_refine3D

submodule(simple_abinitio_utils) simple_abinitio_controller
implicit none
#include "simple_local_flags.inc"

integer, parameter :: NSPACE(8)                  = [500,1000,1000,1000,2500,2500,5000,5000]
integer, parameter :: SHC_REFINE_STAGE           = 1     ! shc-style refinement  stages 1-2
integer, parameter :: PROB_REFINE_STAGE          = 3     ! prob refinement       stages 3-5
integer, parameter :: PROB_NEIGH_REFINE_STAGE    = 6     ! prob_neigh refinement stages 6-8
integer, parameter :: STOCH_SAMPL_STAGE          = 5     ! we switch from greedy to stochastic balanced class sampling
integer, parameter :: LPAUTO_STAGE               = 6     ! we switch on automatic low-pass limit
integer, parameter :: NSPACE_SUB                 = 126

type :: refine3D_stage_cfg
    type(string) :: ml_reg, fillin
    type(string) :: refine, trail_rec, pgrp, balance, filt_mode, automsk, nu_refine
    integer :: iter, inspace, inspace_sub, imaxits
    real    :: trs, frac_best, overlap, fracsrch, lpstart, lpstop
    real    :: snr_noise_reg, gaufreq, update_frac_dyn
end type refine3D_stage_cfg

contains

    integer function active_refine3D_nstages() result(nstages)
        nstages = min(nstages_refine3D, size(NSPACE), size(MAXITS))
        if( allocated(lpinfo) ) nstages = min(nstages, size(lpinfo))
    end function active_refine3D_nstages

    module procedure set_cline_refine3D
        type(refine3D_stage_cfg) :: cfg
        call build_refine3D_stage_cfg( cfg, params, istage, l_cavgs )
        call emit_refine3D_stage_cfg( cfg, params, istage, l_cavgs )
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
        if( istage == active_refine3D_nstages() )then
            cfg%fillin = 'yes'
            if( params%nstates > 1 ) cfg%fillin = 'no'
            cfg%update_frac_dyn = update_frac
        else
            cfg%fillin = 'no'
            cfg%update_frac_dyn = update_frac
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
        select case(trim(params%multivol_mode))
            case('single')
                if( istage >= TRAILREC_STAGE_SINGLE ) cfg%trail_rec = 'yes'
            case('independent')
                if( istage >= TRAILREC_STAGE_MULTI  ) cfg%trail_rec = 'yes'
            case('docked')
                if( istage == active_refine3D_nstages() )then
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

    subroutine set_refine3D_filtering_policy( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(inout) :: cfg
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: istage
        logical,                  intent(in)    :: l_cavgs
        cfg%filt_mode  = 'none'
        cfg%nu_refine  = 'no'
        cfg%lpstart    = 0.
        cfg%lpstop     = 0.
        if( l_cavgs ) return
        if( istage >= LPAUTO_STAGE .and. (l_lpauto .or. l_nonuniform) )then
            cfg%filt_mode = trim(params%filt_mode)
            if( cfg%filt_mode.eq.'uniform' )then
                cfg%lpstart = lpinfo(istage - 1)%lp
                if( istage == active_refine3D_nstages() )then
                    cfg%lpstop = lpinfo(istage)%smpd_crop * 2.
                else
                    cfg%lpstop = lpinfo(istage + 1)%lp
                endif
            endif
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
                if( params%nstates > 1 )then
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
                cfg%inspace_sub = NSPACE_SUB
        end select
    end subroutine apply_refine3D_search_overrides

    subroutine emit_refine3D_stage_cfg( cfg, params, istage, l_cavgs )
        type(refine3D_stage_cfg), intent(in) :: cfg
        class(parameters),        intent(in) :: params
        integer,                  intent(in) :: istage
        logical,                  intent(in) :: l_cavgs
        call cline_refine3D%set('prg',                     'refine3D')
        if( istage == active_refine3D_nstages() )then
            call cline_refine3D%set('update_frac',        cfg%update_frac_dyn)
            call cline_refine3D%set('fillin',             cfg%fillin)
        else
            call cline_refine3D%set('update_frac',        update_frac)
            call cline_refine3D%delete('fillin')
        endif
        call cline_refine3D%set('box_crop',               lpinfo(istage)%box_crop)
        call cline_refine3D%set('startit',                cfg%iter)
        call cline_refine3D%set('which_iter',             cfg%iter)
        call cline_refine3D%set('pgrp',                   cfg%pgrp)
        call cline_refine3D%set('refine',                 cfg%refine)
        call cline_refine3D%set('balance',                cfg%balance)
        call cline_refine3D%set('trail_rec',              cfg%trail_rec)
        call cline_refine3D%set('filt_mode',              cfg%filt_mode)
        call cline_refine3D%set('nu_refine',              cfg%nu_refine)
        if( cfg%filt_mode.eq.'uniform' )then
            call cline_refine3D%set('lpstart',            cfg%lpstart)
            call cline_refine3D%set('lpstop',             cfg%lpstop)
        else
            call cline_refine3D%delete('lpstart')
            call cline_refine3D%delete('lpstop')
        endif
        call cline_refine3D%set('automsk',                cfg%automsk)
        if( cfg%automsk.eq.'no' )then
            ! frequency-limited refinement
            call cline_refine3D%set('lp',                 lpinfo(istage)%lp)
        else
            ! gold-standard refinement
            call cline_refine3D%delete('lp')
            write(logfhandle,'(A,I0,A)') 'emit_refine3D_stage_cfg: stage=', istage, ' automsk active, deleting lp for gold-standard refinement'
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
    end subroutine emit_refine3D_stage_cfg

end submodule simple_abinitio_controller

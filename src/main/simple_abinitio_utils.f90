! concrete commander: high-level workflows
module simple_abinitio_utils
use simple_commander_module_api
use simple_abinitio_config
use simple_commanders_volops,   only: commander_postprocess, commander_symmetrize_map
use simple_commanders_mask,     only: commander_automask
use simple_cluster_seed,        only: gen_labelling
use simple_class_frcs,          only: class_frcs
use simple_decay_funs,          only: calc_update_frac_dyn
implicit none

public :: prep_class_command_lines
public :: set_symmetry_class_vars
public :: set_lplims_from_frcs
public :: set_cline_refine3D
public :: exec_refine3D
public :: symmetrize
public :: calc_start_rec
public :: calc_rec4polar
public :: randomize_states
public :: gen_ortho_reprojs4viz
public :: calc_final_rec
public :: postprocess_final_rec
public :: generate_random_volumes
private
#include "simple_local_flags.inc"

contains

    subroutine prep_class_command_lines( cline, projfile )
        class(cmdline), intent(in) :: cline
        class(string),  intent(in) :: projfile
        cline_refine3D      = cline
        cline_symmap        = cline
        cline_reconstruct3D = cline
        cline_postprocess   = cline
        cline_reproject     = cline
        ! refine3D
        call cline_refine3D%set('prg',                         'refine3D')
        call cline_refine3D%set('pgrp',                  params_glob%pgrp)
        call cline_refine3D%set('projfile',                      projfile)
        ! symmetrization
        call cline_symmap%set('prg',                     'symmetrize_map')
        call cline_symmap%set('pgrp',                    params_glob%pgrp)
        call cline_symmap%set('projfile',                        projfile)
        call cline_symmap%set('center',                             'yes')
        if( .not. cline_symmap%defined('cenlp') )then
        call cline_symmap%set('cenlp',                      CENLP_DEFAULT)
        endif
        call cline_symmap%set('hp',                        params_glob%hp)
        ! re-reconstruct volume
        call cline_reconstruct3D%set('prg',               'reconstruct3D')
        call cline_reconstruct3D%set('box',               params_glob%box)
        call cline_reconstruct3D%set('smpd',             params_glob%smpd)
        call cline_reconstruct3D%set('projfile',                 projfile)
        call cline_reconstruct3D%set('pgrp',             params_glob%pgrp)
        call cline_reconstruct3D%set('ml_reg',                       'no')
        call cline_reconstruct3D%set('needs_sigma',                  'no')
        call cline_reconstruct3D%set('objfun',                       'cc')
        call cline_reconstruct3D%delete('polar')
        ! no fractional update
        call cline_reconstruct3D%delete('update_frac')
        ! individual particles reconstruction
        call cline_reconstruct3D%set('projrec', 'no')
        ! postprocess volume
        call cline_postprocess%set('prg',                   'postprocess')
        call cline_postprocess%set('projfile',                   projfile)
        call cline_postprocess%set('mkdir',                          'no')
        call cline_postprocess%set('imgkind',                       'vol')
        call cline_postprocess%delete('lp')   ! to obtain optimal filtration
        ! re-project volume, only with cavgs
        call cline_reproject%set('prg',                       'reproject')
        call cline_reproject%set('pgrp',                 params_glob%pgrp)
        call cline_reproject%set('outstk', 'reprojs'//params_glob%ext%to_char())
        call cline_reproject%set('smpd',                 params_glob%smpd)
        call cline_reproject%set('box',                   params_glob%box)
        call cline_reproject%set('oritab',               'final_oris.txt')
        call cline_reproject%set('nstates',           params_glob%nstates)
        call cline_reproject%delete('projfile')
    end subroutine prep_class_command_lines

    subroutine set_symmetry_class_vars
        type(string) :: pgrp, pgrp_start
        pgrp           = lowercase(trim(params_glob%pgrp))
        pgrp_start     = lowercase(trim(params_glob%pgrp_start))
        l_srch4symaxis = pgrp .ne. pgrp_start
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( pgrp_start.ne.'c1' .or. pgrp.ne.'c1' )then
            se1 = sym(pgrp_start%to_char())
            se2 = sym(pgrp%to_char())
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(pgrp%to_char()) )THROW_HARD('Incompatible symmetry groups; exec_abinitio3D')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(pgrp_start%to_char()) )THROW_HARD('Incompatible symmetry groups; exec_abinitio3D')
            endif
        endif
    end subroutine set_symmetry_class_vars

    subroutine set_lplims_from_frcs( spproj, l_cavgs, lpstart, lpstop )
        class(sp_project), intent(inout) :: spproj
        logical,           intent(in)    :: l_cavgs
        real, optional,    intent(in)    :: lpstart, lpstop
        real,         allocatable :: frcs_avg(:)
        integer,      allocatable :: states(:)
        type(string)     :: frcs_fname
        type(class_frcs) :: clsfrcs
        real             :: lpfinal
        integer          :: filtsz
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        ! work out low-pass limits and downscaling parameters
        params_glob%frcs = frcs_fname
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        states = nint(spproj%os_cls2D%get_all('state'))
        call clsfrcs%avg_frc_getter(frcs_avg, states)
        if( allocated(lpinfo) ) deallocate(lpinfo)
        allocate(lpinfo(NSTAGES))
        lpfinal = max(LPSTOP_BOUNDS(1),calc_lplim_final_stage(3))
        lpfinal = min(LPSTOP_BOUNDS(2),lpfinal)
        if( present(lpstop) ) lpfinal = max(lpstop,lpfinal)
        if( present(lpstart) )then
            call lpstages(params_glob%box, NSTAGES, frcs_avg, params_glob%smpd,&
            &lpstart, lpstart, lpfinal, lpinfo, l_cavgs, verbose=.not.l_polar)
        else
            call lpstages(params_glob%box, NSTAGES, frcs_avg, params_glob%smpd,&
            &LPSTART_BOUNDS(1), LPSTART_BOUNDS(2), lpfinal, lpinfo, l_cavgs, verbose=.not.l_polar)
        endif
        ! Stepped increase of dimensions with polar representation
        call edit_lpstages4polar( NSPACE_PHASE_POLAR, NSTAGES, lpinfo )
        ! cleanup
        call clsfrcs%kill
        contains

            function calc_lplim_final_stage( nbest ) result( lplim )
                integer, intent(in)  :: nbest
                real,    allocatable :: res(:), tmp_rarr(:)
                integer, allocatable :: tmp_iarr(:)
                real :: lplim
                tmp_rarr  = spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_final_stage

    end subroutine set_lplims_from_frcs

    subroutine set_cline_refine3D( istage, l_cavgs )
        integer,          intent(in)  :: istage
        logical,          intent(in)  :: l_cavgs
        type(string) :: sh_first, prob_sh, ml_reg, fillin, cavgw, ref_type
        type(string) :: refine, icm, trail_rec, pgrp, balance, lp_auto, automsk
        integer :: iphase, iter, inspace, imaxits, nsample_dyn, nspace_phase
        real    :: trs, frac_best, overlap, fracsrch, lpstart, lpstop, snr_noise_reg, gaufreq
        ! iteration number bookkeeping
        iter = 0
        if( cline_refine3D%defined('endit') )then
            iter = cline_refine3D%get_iarg('endit')
        endif
        iter = iter + 1
        ! dynamic update frac
        if( istage == NSTAGES )then
            fillin = 'yes'
            if( params_glob%nstates > 1 ) fillin = 'no' ! fill-in doesn't work with multi-state
            if( l_nsample_stop_given )then
                update_frac_dyn = real(nsample_minmax(2)) / real(nptcls_eff)
            else if( l_nsample_given )then
                update_frac_dyn = update_frac
            else
                ! we change the sampling method for the last stage (accelerated refinement)
                nsample_dyn     = nint(UPDATE_FRAC_MIN * real(nptcls_eff) / real(params_glob%nstates))
                nsample_dyn     = max(NSAMPLE_MINMAX_DEFAULT(1), min(NSAMPLE_MINMAX_DEFAULT(2), nsample_dyn))
                update_frac_dyn = real(nsample_dyn * params_glob%nstates) / real(nptcls_eff)
            endif
        else
            fillin = 'no'
            update_frac_dyn = calc_update_frac_dyn(nptcls_eff, params_glob%nstates, nsample_minmax, iter, maxits_dyn)
        endif
        update_frac_dyn = min(UPDATE_FRAC_MAX, update_frac_dyn) ! to ensure fractional update is always on
        ! symmetry
        pgrp = trim(params_glob%pgrp)
        if( l_srch4symaxis )then
            if( istage <= SYMSRCH_STAGE )then
                ! need to replace original point-group flag with c1/pgrp_start
                pgrp = trim(params_glob%pgrp_start)
            endif
        endif
        ! refinement mode
        refine  = 'shc_smpl'
        prob_sh = 'no'
        if( istage >= PROBREFINE_STAGE )then
            refine  = 'prob'
            prob_sh = 'yes'
        endif
        if( trim(params_glob%multivol_mode).eq.'input_oris_fixed' )then ! only state sorting, no 3D ori refinement
            refine = 'prob_state'
        endif
        ! ICM regularization
        icm = 'no'
        if( istage >= ICM_STAGE ) icm = 'yes'
        ! forcing icm to no if user sets icm to no
        if( cline_refine3D%defined('icm') )then
            if( trim(params_glob%icm).eq.'no')then
                icm = 'no'
            endif
        endif
        ! balance
        balance = 'yes'
        ! Gaussian filtering of reference volume
        gaufreq = -1.
        if( istage <= params_glob%gauref_last_stage )then
            gaufreq = lpinfo(istage)%lp
        endif
        ! trailing reconstruction
        trail_rec = 'no'
        select case(trim(params_glob%multivol_mode))
            case('single')
                if( istage >= TRAILREC_STAGE_SINGLE ) trail_rec = 'yes'
            case('independent')
                if( istage >= TRAILREC_STAGE_MULTI  ) trail_rec = 'yes'
            case('docked')
                if( istage == NSTAGES )then
                    trail_rec = 'no'
                else if( istage >= TRAILREC_STAGE_SINGLE )then
                    trail_rec = 'yes'
                endif
            case('input_oris_fixed')
                trail_rec = 'no'
            case('input_oris_start')
                trail_rec = 'no'
            case DEFAULT
                trail_rec = 'no' ! defaults to 'no' for safety
        end select
        ! automatic low-pass limit estimation
        lp_auto = 'no'
        if( istage >= LPAUTO_STAGE .and. l_lpauto )then
            lp_auto = trim(params_glob%lp_auto)
            lpstart = lpinfo(istage - 1)%lp
            if( istage == NSTAGES )then
                lpstop = lpinfo(istage)%smpd_crop * 2. ! Nyqvist limit
            else
                lpstop = lpinfo(istage + 1)%lp
            endif
        endif
        ! automasking
        automsk = 'no'
        if( (.not.l_cavgs) .and. (.not.l_polar) )then
            if( istage >= AUTOMSK_STAGE .and. l_automsk )then
                automsk = 'yes'
            endif
        endif
        ! cavgs weights, not supported for particles
        cavgw = 'no'
        if( l_cavgs )then
            if( (trim(params_glob%cavgw).eq.'yes') .and. (istage>=CAVGWEIGHTS_STAGE))then
                cavgw = 'yes'
            endif
        endif
        ! phase logics
        iphase = 0
        if(      istage <= PHASES(1) )then
            iphase = 1
        else if( istage <= PHASES(2) )then
            iphase = 2
        else if( istage <= PHASES(3) )then
            iphase = 3
        else 
            THROW_HARD('Invalid istage index')
        endif
        ! phase control parameters
        select case(iphase)
            case(1)
                inspace       = NSPACE(1)
                imaxits       = MAXITS(istage)
                trs           = 0.
                sh_first      = 'no'
                ml_reg        = 'no'
                frac_best     = 1.0 ! means it does not control sampling, greedy selection
                overlap       = 0.90
                fracsrch      = 90.
                snr_noise_reg = 2.0
            case(2)
                inspace       = NSPACE(2)
                imaxits       = MAXITS(istage)
                trs           = lpinfo(istage)%trslim
                sh_first      = 'yes'
                ml_reg        = 'yes'
                if( istage >= STOCH_SAMPL_STAGE )then
                frac_best     = 0.5 ! means sampling is done from top-ranking 50% particles in class
                else
                frac_best     = 1.0 ! means it does not control sampling, greedy selection
                endif
                overlap       = 0.90
                fracsrch      = 90.
                snr_noise_reg = 4.0
            case(3)
                inspace       = NSPACE(3)
                imaxits       = MAXITS(istage)
                trs           = lpinfo(istage)%trslim
                sh_first      = 'yes'
                ml_reg        = 'yes'
                if( params_glob%nstates > 1 )then
                ! turn off balancing
                frac_best     = 0.98 ! max out balanced sampling
                else
                frac_best     = 0.85 ! means sampling is done from top-ranking 85% particles in class
                endif
                if( istage == NSTAGES )then
                overlap       = 0.95
                fracsrch      = 95.
                else
                overlap       = 0.90
                fracsrch      = 90.
                endif
                snr_noise_reg = 6.0
        end select
        ! Specific options for polar representation
        if( l_polar )then
            ! # of projection directions
            nspace_phase = 3
            if( istage <= NSPACE_PHASE_POLAR(1) )then
                nspace_phase = 1
            else if( istage <= NSPACE_PHASE_POLAR(2) )then
                nspace_phase = 2
            endif
            inspace  = NSPACE_POLAR(nspace_phase)
            ! volume filtering
            icm = 'no'
            if( trim(params_glob%gauref).eq.'no' ) gaufreq = -1.
            if( ml_reg.eq.'yes' )                  gaufreq = -1.
            ! CL-based approach
            ref_type = trim(params_glob%ref_type)
        endif
        ! turn off ML-regularization when icm is on
        if( icm.eq.'yes' ) ml_reg = 'no'
        ! projection directions
        if( cline_refine3D%defined('nspace_max') )then
            inspace = min(inspace, params_glob%nspace_max)
        endif
        ! Development
        if( trim(params_glob%algorithm).eq.'mimic2D' )then
            imaxits  = params_glob%maxits_between   ! # of iterations
            inspace  = min(inspace, 2000)           ! nspace
            refine   = 'snhc_smpl'                  ! refinement
            sh_first = 'no'                         ! not implemented
            prob_sh  = 'no'
            lp_auto  = 'no'
            icm      = 'no'                         ! icm filter off
            if( istage <= params_glob%gauref_last_stage )then
             gaufreq = lpinfo(istage)%lp            ! gaussian on, ML off
             ml_reg  = 'no'
            else
             gaufreq = -1.0                         ! gaussian off, ML on
             ml_reg  = 'yes'
            endif
            call cline_refine3D%set('extr_lim', (NSTAGES-1)*imaxits)
            if( istage <= 1 )then
                call cline_refine3D%set('extr_iter', 1)
            elseif( istage < NSTAGES )then
                call cline_refine3D%set('extr_iter', (istage-1)*imaxits+1)
            else
                call cline_refine3D%delete('extr_iter')
                call cline_refine3D%delete('extr_lim')
            endif
        endif
        ! command line update
        call cline_refine3D%set('prg',                     'refine3D')
        ! class global control parameters
        if( l_update_frac_dyn .or. istage == NSTAGES )then
        call cline_refine3D%set('update_frac',        update_frac_dyn)
        call cline_refine3D%set('fillin',                      fillin)
        else
        call cline_refine3D%set('update_frac',            update_frac)
        call cline_refine3D%delete('fillin')
        endif
        call cline_refine3D%set('lp',             lpinfo(istage  )%lp)
        call cline_refine3D%set('smpd_crop', lpinfo(istage)%smpd_crop)
        call cline_refine3D%set('box_crop',   lpinfo(istage)%box_crop)
        ! iteration number
        call cline_refine3D%set('startit',                       iter)
        call cline_refine3D%set('which_iter',                    iter)
        ! dynamic control parameters
        call cline_refine3D%set('pgrp',                          pgrp)
        call cline_refine3D%set('refine',                      refine)
        call cline_refine3D%set('balance',                    balance)
        call cline_refine3D%set('trail_rec',                trail_rec)
        call cline_refine3D%set('lp_auto',                    lp_auto)
        if( lp_auto.eq.'yes' )then
        call cline_refine3D%set('lpstart',                    lpstart)
        call cline_refine3D%set('lpstop',                      lpstop)
        else
        call cline_refine3D%delete('lpstart')
        call cline_refine3D%delete('lpstop')
        endif
        call cline_refine3D%set('automsk',                    automsk)
        if( l_cavgs )then
        call cline_refine3D%set('cavgw',                        cavgw)
        endif
        ! phase control parameters
        call cline_refine3D%set('nspace',                     inspace)
        call cline_refine3D%set('maxits',                     imaxits)
        call cline_refine3D%set('trs',                            trs)
        call cline_refine3D%set('sh_first',                  sh_first)
        call cline_refine3D%set('prob_sh',                    prob_sh)
        call cline_refine3D%set('ml_reg',                      ml_reg)
        call cline_refine3D%set('icm',                            icm)
        call cline_refine3D%set('frac_best',                frac_best)
        call cline_refine3D%set('overlap',                    overlap)
        call cline_refine3D%set('fracsrch',                  fracsrch)
        if( l_cavgs )then
        call cline_refine3D%set('snr_noise_reg',        snr_noise_reg)
        call cline_refine3D%delete('update_frac') ! never on cavgs
        else
        call cline_refine3D%delete('snr_noise_reg')
        endif
        if( gaufreq > 0.)then
        call cline_refine3D%set('gauref',                      'yes')
        call cline_refine3D%set('gaufreq',                   gaufreq)
        else
        call cline_refine3D%delete('gauref')
        call cline_refine3D%delete('gaufreq')
        endif
        if( l_polar )then
            call cline_refine3D%set('center_type',           'params')
            call cline_refine3D%set('ref_type',              ref_type)
        endif
    end subroutine set_cline_refine3D

    subroutine exec_refine3D( istage, xrefine3D )
        integer,               intent(in)    :: istage
        class(commander_base), intent(inout) :: xrefine3D
        type(string) :: stage, str_state, vol_name, vol_pproc
        integer      :: state
        call cline_refine3D%delete('endit')
        call xrefine3D%execute_safe(cline_refine3D)
        call del_files(DIST_FBODY,      params_glob%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params_glob%nparts,ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        stage = '_stage_'//int2str(istage)
        if( .not. l_polar )then
            do state = 1, params_glob%nstates
                str_state = int2str_pad(state,2)
                vol_name  = string(VOL_FBODY)//str_state//params_glob%ext
                vol_pproc = add2fbody(vol_name, params_glob%ext, PPROC_SUFFIX)
                if( file_exists(vol_name) ) call simple_copy_file(vol_name,  add2fbody(vol_name, params_glob%ext,stage))
                if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params_glob%ext,stage))
            enddo
        endif
    end subroutine exec_refine3D

    subroutine symmetrize( istage, spproj, projfile, xreconstruct3D )
        integer,                         intent(in)    :: istage
        class(sp_project),               intent(inout) :: spproj
        class(string),                   intent(in)    :: projfile
        class(commander_base), optional, intent(inout) :: xreconstruct3D
        type(commander_symmetrize_map) :: xsymmap
        type(cmdline)                  :: cline_asymrec, cline_symrec
        type(string) :: vol_iter, vol_sym
        real :: lpsym
        if( l_symran )then
            call se1%symrandomize(spproj%os_ptcl3D)
            call spproj%write_segment_inside('ptcl3D', projfile)
        endif
        if( l_srch4symaxis )then
            ! asymmetric/low symmetry reconstruction
            if( l_polar )then
                if( .not.present(xreconstruct3D) )then
                    THROW_HARD('Reconstructor required with polar=yes')
                endif
                cline_asymrec = cline_refine3D
                call cline_asymrec%set('prg',        'reconstruct3D')
                call cline_asymrec%set('mkdir',      'no')
                call cline_asymrec%set('projfile',   projfile)
                call cline_asymrec%set('pgrp',       params_glob%pgrp_start)
                call cline_asymrec%set('ml_reg',     'no') ! no ml reg for now
                call cline_asymrec%set('objfun',     'cc')
                call cline_asymrec%set('needs_sigma','no')
                call cline_asymrec%delete('which_iter')
                call cline_asymrec%delete('endit')
                call xreconstruct3D%execute_safe(cline_asymrec)
                vol_iter = 'asymmetric_map'//params_glob%ext%to_char()
                call simple_copy_file(string(VOL_FBODY)//int2str_pad(1,2)//params_glob%ext%to_char(), vol_iter)
                call cline_asymrec%kill
            else
                ! Volume from previous stage
                vol_iter = string(VOL_FBODY)//STR_STATE_GLOB//params_glob%ext%to_char()
            endif
            ! symmetry determination & map symmetrization
            if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to map symmetrization does not exist')
            call cline_symmap%set('vol1', vol_iter)
            call cline_symmap%set('smpd', lpinfo(istage)%smpd_crop)
            call cline_symmap%set('box',  lpinfo(istage)%box_crop)
            vol_sym = 'symmetrized_map'//params_glob%ext%to_char()
            call cline_symmap%set('outvol', vol_sym)
            lpsym = max(LPSYMSRCH_LB,lpinfo(SYMSRCH_STAGE)%lp)
            call cline_symmap%set('lp', lpsym)
            write(logfhandle,'(A,F5.1)') '>>> DID SET MAP SYMMETRIZATION LOW-PASS LIMIT (IN A) TO: ', lpsym
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> MAP SYMMETRIZATION'
            write(logfhandle,'(A)') '>>>'
            call xsymmap%execute_safe(cline_symmap)
            call del_file('SYMAXIS_SEARCH_FINISHED')
            if( present(xreconstruct3D) )then
                ! symmetric reconstruction
                cline_symrec = cline_refine3D
                call cline_symrec%set('prg',        'reconstruct3D')
                call cline_symrec%set('mkdir',      'no')
                call cline_symrec%set('projfile',   projfile)
                call cline_symrec%set('pgrp',       params_glob%pgrp)
                call cline_symrec%set('which_iter', cline_refine3D%get_iarg('endit'))
                call cline_symrec%delete('endit')
                call xreconstruct3D%execute_safe(cline_symrec)
                vol_sym = VOL_FBODY//int2str_pad(1,2)//params_glob%ext%to_char()
                call simple_copy_file(vol_sym, string('symmetric_map')//params_glob%ext)
                call cline_symrec%kill
            endif
            call cline_refine3D%set('vol1', vol_sym)
        endif
    end subroutine symmetrize

    subroutine calc_start_rec( projfile, xreconstruct3D, istage )
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        integer,               intent(in)    :: istage
        type(commander_automask)       :: xautomask
        type(string)  :: str_state, vol, str, vol_even, vol_odd
        type(cmdline) :: cline_startrec, cline_automask
        integer       :: state
        cline_startrec = cline_refine3D
        call cline_startrec%set('prg',         'reconstruct3D')
        call cline_startrec%set('mkdir',       'no')
        call cline_startrec%set('projfile',    projfile)
        call cline_startrec%set('pgrp',        params_glob%pgrp)
        call cline_startrec%set('objfun',      'cc') ! ugly, but this is how it works in parameters 
        call cline_startrec%set('box_crop',    lpinfo(istage)%box_crop)
        call cline_startrec%set('projrec',     'no')
        call cline_startrec%delete('update_frac')    ! use all particles that have been updated
        call cline_startrec%delete('which_iter')
        call cline_startrec%delete('endit')
        call cline_startrec%delete('needs_sigma')
        call cline_startrec%delete('sigma_est')
        call cline_startrec%delete('automsk') ! no automask generated
        call cline_startrec%delete('mskfile') ! no masked FSC
        ! endif
        call xreconstruct3D%execute_safe(cline_startrec)
        do state = 1,params_glob%nstates
            ! rename volumes and update cline
            str_state = int2str_pad(state,2)
            vol       = string(VOL_FBODY)//str_state//params_glob%ext
            str       = string(STARTVOL_FBODY)//str_state//params_glob%ext
            call      simple_rename(vol, str)
            params_glob%vols(state) = str
            vol       = 'vol'//int2str(state)
            call      cline_refine3D%set(vol, str)
            vol_even  = string(VOL_FBODY)//str_state//'_even'//params_glob%ext%to_char()
            str       = string(STARTVOL_FBODY)//str_state//'_even_unfil'//params_glob%ext%to_char()
            call      simple_copy_file(vol_even, str)
            str       = string(STARTVOL_FBODY)//str_state//'_even'//params_glob%ext%to_char()
            call      simple_rename(vol_even, str)
            vol_odd   = string(VOL_FBODY)//str_state//'_odd' //params_glob%ext%to_char()
            str       = string(STARTVOL_FBODY)//str_state//'_odd_unfil'//params_glob%ext%to_char()
            call      simple_copy_file(vol_odd, str)
            str       = string(STARTVOL_FBODY)//str_state//'_odd'//params_glob%ext%to_char()
            call      simple_rename(vol_odd, str)
        enddo
        if( istage >= AUTOMSK_STAGE .and. l_automsk )then
            str_state = int2str_pad(1,2)
            vol_even  = string(STARTVOL_FBODY)//str_state//'_even'//params_glob%ext%to_char()
            vol_odd   = string(STARTVOL_FBODY)//str_state//'_odd'//params_glob%ext%to_char()
            call cline_automask%set('vol1', vol_odd)
            call cline_automask%set('vol2', vol_even)
            call cline_automask%set('smpd', lpinfo(istage)%smpd_crop)
            call cline_automask%set('amsklp', params_glob%amsklp)
            call cline_automask%set('automsk', 'yes')
            call cline_automask%set('mkdir',    'no')
            call cline_automask%set('nthr', params_glob%nthr)
            call xautomask%execute_safe(cline_automask)
            params_glob%mskfile = MSKVOL_FILE
            call cline_refine3D%set('mskfile', MSKVOL_FILE)
        endif
        call cline_startrec%kill
    end subroutine calc_start_rec

    ! Performs reconstruction at some set stages when polar=yes
    subroutine calc_rec4polar( xreconstruct3D, istage )
        use simple_class_frcs, only: class_frcs
        class(commander_base), intent(inout) :: xreconstruct3D
        integer,               intent(in)    :: istage
        real, allocatable :: fsc(:)
        type(cmdline)     :: cline_rec
        type(class_frcs)  :: frcs
        type(string)      :: src, dest, sstate, sstage, ext, pgrp
        integer :: i, inspace
        if( trim(params_glob%polar) /= 'yes' ) return
        if( (istage /= NSPACE_PHASE_POLAR(1)+1) .and. (istage /= NSPACE_PHASE_POLAR(2)+1) ) return
        ! Reconstruction
        pgrp = trim(params_glob%pgrp)
        if( istage <= SYMSRCH_STAGE ) pgrp = trim(params_glob%pgrp_start)
        cline_rec = cline_refine3D
        call cline_rec%set('prg',       'reconstruct3D')
        call cline_rec%set('mkdir',     'no')
        call cline_rec%set('projfile',  params_glob%projfile)
        call cline_rec%set('pgrp',      pgrp)
        call cline_rec%set('box_crop',  lpinfo(istage)%box_crop)
        call cline_rec%set('projrec',   'yes')
        call cline_rec%set('trail_rec', 'no')
        if( cline_rec%get_carg('ml_reg').ne.'yes' )then
            call cline_rec%set('objfun',    'cc')
        endif
        call cline_rec%delete('update_frac')
        call cline_rec%delete('endit')
        call cline_rec%delete('needs_sigma')
        call cline_rec%delete('automsk')
        call cline_rec%delete('mskfile')
        call xreconstruct3D%execute_safe(cline_rec)
        ! preserve volume (e/o will be overwritten in next iteration)
        sstate = int2str_pad(1,2)
        sstage = int2str_pad(istage-1,2)
        ext    = params_glob%ext
        src    = string(VOL_FBODY)//sstate//ext
        dest   = string(VOL_FBODY)//sstate//'_stage_'//sstage//ext
        call simple_rename(src, dest)
        ! Update refine3D command line
        call cline_refine3D%set('vol1', dest)
        ! Update FRCS
        fsc     = file2rarr(string(FSC_FBODY)//int2str_pad(1,2)//BIN_EXT)
        inspace = cline_refine3D%get_iarg('nspace')
        call frcs%new(inspace, lpinfo(istage)%box_crop, lpinfo(istage)%smpd_crop, 1)
        do i = 1,inspace
            call frcs%set_frc(i,fsc)
        enddo
        call frcs%write(string(FRCS_FILE))
        ! cleanup
        deallocate(fsc)
        call frcs%kill
        call cline_rec%kill
    end subroutine calc_rec4polar

    subroutine randomize_states( spproj, projfile, xreconstruct3D, istage )
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        integer,               intent(in)    :: istage
        call spproj%read_segment('ptcl3D', projfile)
        call gen_labelling(spproj%os_ptcl3D, params_glob%nstates, 'squared_uniform')
        call spproj%write_segment_inside(params_glob%oritype, projfile)
        call cline_refine3D%set(     'nstates', params_glob%nstates)
        call cline_reconstruct3D%set('nstates', params_glob%nstates)
        call cline_postprocess%set(  'nstates', params_glob%nstates)
        call cline_reproject%set(    'nstates', params_glob%nstates)
        call calc_start_rec(projfile, xreconstruct3D, istage=istage)
    end subroutine randomize_states

    subroutine gen_ortho_reprojs4viz( spproj )
        type(sp_project), intent(inout) :: spproj
        type(string) :: str_state, fname
        type(image)  :: final_vol, reprojs
        integer      :: state, ifoo, ldim(3)
        real         :: smpd
        call spproj%read_segment('out', params_glob%projfile)
        do state = 1, params_glob%nstates
            if( .not.spproj%isthere_in_osout('vol', state) )cycle   ! empty-state case
            str_state = int2str_pad(state,2)
            fname = string(VOL_FBODY)//str_state//params_glob%ext
            if( .not. file_exists(fname) )cycle
            exit
        enddo
        call find_ldim_nptcls(fname, ldim, ifoo, smpd)
        call final_vol%new(ldim, smpd)
        do state = 1, params_glob%nstates
            str_state = int2str_pad(state,2)
            if( spproj%isthere_in_osout('vol', state) )then
                str_state = int2str_pad(state,2)
                fname     = VOL_FBODY//str_state%to_char()//params_glob%ext%to_char()
                if( .not. file_exists(fname) )cycle
                call final_vol%read(fname)
                call final_vol%generate_orthogonal_reprojs(reprojs)
                call reprojs%write_jpg(string('orthogonal_reprojs_state')//str_state//'.jpg')
                call reprojs%kill
            endif
        enddo
        call final_vol%kill
    end subroutine gen_ortho_reprojs4viz

    subroutine calc_final_rec( spproj, projfile, xreconstruct3D )
        class(sp_project),     intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(commander_base), intent(inout) :: xreconstruct3D
        type(string) :: str_state, vol_name
        integer      :: state, pop
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call xreconstruct3D%execute_safe(cline_reconstruct3D)
        call spproj%read_segment('out', projfile)
        call spproj%read_segment('ptcl3D', projfile)
        do state = 1, params_glob%nstates
            pop = spproj%os_ptcl3D%get_pop(state, 'state')
            if( pop == 0 )cycle     ! empty-state case
            str_state = int2str_pad(state,2)
            vol_name  = string(VOL_FBODY)//str_state//params_glob%ext
            if( .not. file_exists(vol_name) )cycle
            call spproj%add_vol2os_out(vol_name, params_glob%smpd, state, 'vol', pop=pop)
            call spproj%add_fsc2os_out(string(FSC_FBODY)//str_state//BIN_EXT, state, params_glob%box)
        enddo
        call spproj%write_segment_inside('out', projfile)
    end subroutine calc_final_rec

    subroutine postprocess_final_rec( spproj )
        class(sp_project), intent(in) :: spproj
        type(commander_postprocess)   :: xpostprocess
        type(string) :: str_state, vol_name, vol_pproc, vol_pproc_mirr, vol_final
        integer :: state
        do state = 1, params_glob%nstates
            if( .not.spproj%isthere_in_osout('vol', state) )cycle ! empty-state case
            str_state      = int2str_pad(state,2)
            vol_name       = string(VOL_FBODY)//str_state//params_glob%ext  ! reconstruction from particles stored in project
            if( .not. file_exists(vol_name) )cycle
            call cline_postprocess%set('state', state)
            call xpostprocess%execute_safe(cline_postprocess)
            vol_pproc      = add2fbody(vol_name,params_glob%ext, PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol_name,params_glob%ext, PPROC_SUFFIX//MIRR_SUFFIX)
            vol_final      = string(REC_FBODY)//str_state//params_glob%ext
            if( file_exists(vol_name)       ) call simple_copy_file(vol_name,    vol_final)
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      add2fbody(vol_final,params_glob%ext, PPROC_SUFFIX))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, add2fbody(vol_final,params_glob%ext, PPROC_SUFFIX//MIRR_SUFFIX))
        enddo
    end subroutine postprocess_final_rec

    ! create starting noise volume(s)
    subroutine generate_random_volumes( box, smpd, cline )
        integer,        intent(in)    :: box
        real,           intent(in)    :: smpd
        type(cmdline),  intent(inout) :: cline
        type(string) :: vol_name
        type(image)  :: noisevol, signal
        real         :: b
        integer      :: s
        call noisevol%new([box,box,box], smpd)
        select case(trim(params_glob%inivol))
            case('rand')
                ! random uniform distribution [0,1]
                ! resulting std dev in projections is ~box/10
                do s = 1, params_glob%nstates
                    call noisevol%ran
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline%set('vol'//int2str(s), vol_name)
                    params_glob%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%ran
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%ran
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
            case('rand_scaled')
                ! random uniform distribution [0,5/box]
                ! resulting std dev in projections is ~0.2
                b = 5.0/real(box)
                do s = 1, params_glob%nstates
                    call noisevol%ran(b=b)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline%set('vol'//int2str(s), vol_name)
                    params_glob%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%ran(b=b)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%ran(b=b)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
            case('sphere')
                ! random normal sphere N(0,5/box) + normal background N(0,5/box)
                b = 5.0/real(box)
                call signal%new([box,box,box], smpd)
                call signal%gauran(0.0, b)
                call signal%mask(0.25*real(box), 'soft', backgr=0.)
                do s = 1, params_glob%nstates
                    call noisevol%gauran(0., b)
                    call noisevol%add(signal)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                    call cline%set('vol'//int2str(s), vol_name)
                    params_glob%vols(s) = vol_name
                    call noisevol%write(vol_name)
                    call noisevol%gauran(0., b)
                    call noisevol%add(signal)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                    call noisevol%write(vol_name)
                    call noisevol%gauran(0., b)
                    call noisevol%add(signal)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                    call noisevol%write(vol_name)
                    vol_name = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                    call noisevol%write(vol_name)
                end do
                call signal%kill
        end select
        call noisevol%kill
    end subroutine generate_random_volumes

end module simple_abinitio_utils

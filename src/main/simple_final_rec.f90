!@descr: the shared final all-particle reconstruction at native sampling
module simple_final_rec
use, intrinsic :: iso_fortran_env, only: int64
use simple_commanders_api
use simple_parameters,        only: parameters
use simple_solve3D_utils,     only: configure_final_pcg_solve_budget, write_final_rec_outputs, &
    &gen_ortho_reprojs4viz
use simple_refine3D_fnames,   only: refine3D_fsc_fname, refine3D_state_vol_fname
use simple_vol_pproc_policy,  only: state_mask_is_compatible
use simple_sigma2_files,      only: canonical_sigma2_consumable
use simple_state_weight_set,  only: state_weight_set
use simple_reconstructor_pcg, only: PCG_WEIGHT_THRESHOLD
implicit none
#include "simple_local_flags.inc"

public :: calc_final_rec, calc_final_rec_addon
private

abstract interface
    !> solve3D_addon membership switch of the working project: l_union restores the union (frozen
    !! rows back), and with l_means republishes the sigma2 state with both global means over it;
    !! otherwise the frozen and retired rows are masked again, so a reconstruction reads the cohort only
    subroutine addon_membership_switch( l_union, l_means )
        logical, intent(in) :: l_union, l_means
    end subroutine addon_membership_switch
end interface

contains

    !> The final all-particle reconstruction at native sampling shared by solve3D, refine3D_auto,
    !! refine3D_states and classify3D_refs, driven by the refinement command line it continues from. Sigmas
    !! not consumable at the native box are rebuilt by bootstrap_rec3D; without ML regularization it ships a cc map.
    !! With m_estimator=flex both routes reconstruct from the project's state weight set, whose identity must
    !! not change while they run, and every weighted map is registered with its hard population, applied
    !! mass and effective sample size.
    subroutine calc_final_rec( params, spproj, projfile, cline_refine, xrec3D, xbootstrap_rec3D, &
        &l_postprocess, lp_snapshot )
        class(parameters),     intent(in)    :: params
        type(sp_project),      intent(inout) :: spproj
        class(string),         intent(in)    :: projfile
        class(cmdline),        intent(in)    :: cline_refine  !< refinement command line the final map continues from
        class(commander_base), intent(inout) :: xrec3D
        class(commander_base), intent(inout) :: xbootstrap_rec3D
        logical,               intent(in)    :: l_postprocess
        real,                  intent(in)    :: lp_snapshot   !< planned low-pass fallback for the diagnostic snapshot
        type(cmdline)          :: cline_final
        type(state_weight_set) :: wset
        type(string)           :: vol_name, stkname, vol_pproc, vol_mirr
        integer, allocatable   :: rec_rows(:)
        integer(int64)         :: wset_generation, wset_digest
        real                   :: smpd, wthres, mass, ess
        integer                :: ldim(3), state, pop, stkind, ind_in_stk, nptcls, bootstrap_sigma_iter, iptcl
        logical                :: l_bootstrap_sigmas
        wset_generation = 0_int64
        wset_digest     = 0_int64
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call spproj%read(projfile) ! ensure we have the latest project info
        call spproj%map_ptcl_ind2stk_ind('ptcl3D', 1, stkind, ind_in_stk)
        stkname = spproj%os_stk%get_str(stkind, 'stk')
        call find_ldim_nptcls(stkname, ldim, nptcls)
        smpd = spproj%os_stk%get(stkind, 'smpd')
        write(logfhandle,'(A,I0,A,F8.4)') '>>> FINAL RECONSTRUCTION SAMPLING: box=', ldim(1), ' smpd=', smpd
        if( params%l_m_estimator_flex )then
            ! the frozen weights the final map is reconstructed from
            call wset%new(spproj, spproj%os_ptcl3D)
            wset_generation = wset%get_generation()
            wset_digest     = wset%get_layout_digest()
            write(logfhandle,'(A,I0,A,I0,A,A,A,I0)') '>>> FINAL RECONSTRUCTION WEIGHTS: state weight set generation ', &
                &wset_generation, ', layout digest ', wset_digest, ', producer ', wset%get_producer(), &
                &', parent state ', wset%get_parent_state()
            call wset%kill
        endif
        call prep_final_rec_cline(params, projfile, cline_refine, l_postprocess, cline_final, 'reconstruct3D')
        l_bootstrap_sigmas = .false.
        if( final_stage_uses_ml_reg(cline_refine) )then
            ! a valid committed state is reused directly; a missing, stale,
            ! wrong-grid or wrong-layout state, or one whose header the current
            ! format refuses, is rebuilt from particle power and
            ! residual-upgraded by bootstrap_rec3D.
            ! Sigmas estimated at a cropped registration box are refreshed at
            ! native sampling before the shipped map.
            l_bootstrap_sigmas = canonical_final_rec_needs_bootstrap()
            if( .not. l_bootstrap_sigmas ) l_bootstrap_sigmas = final_rec_box_changed()
            if( .not. l_bootstrap_sigmas ) write(logfhandle,'(A)') &
                &'>>> FINAL RECONSTRUCTION: reusing committed canonical sigmas'
        endif
        if( l_bootstrap_sigmas )then
            bootstrap_sigma_iter = final_rec_bootstrap_sigma_iter(cline_refine)
            call prep_final_rec_cline(params, projfile, cline_refine, l_postprocess, cline_final, 'bootstrap_rec3D')
            call cline_final%set('which_iter', bootstrap_sigma_iter)
            write(logfhandle,'(A,I0)') '>>> FINAL RECONSTRUCTION BOOTSTRAP SIGMA ITERATION: ', bootstrap_sigma_iter
            if( trim(params%rec_backend) == 'pcg' ) write(logfhandle,'(A,I0)') &
                &'>>> FINAL PCG COLD-SOLVE ITERATION BUDGET: ', cline_final%get_iarg('maxits_pcg')
            ! bootstrap_rec3D owns the whole sequence (image-power seed, euclid ML bootstrap map, one residual
            ! sigma2 pass at the final sampling, shipped euclid ML map); prg=bootstrap_rec3D runs it standalone.
            call xbootstrap_rec3D%execute(cline_final)
        else
            if( trim(params%rec_backend) == 'pcg' ) write(logfhandle,'(A,I0)') &
                &'>>> FINAL PCG COLD-SOLVE ITERATION BUDGET: ', cline_final%get_iarg('maxits_pcg')
            call xrec3D%execute(cline_final)
        endif
        if( .not. l_postprocess )then
            do state = 1, params%nstates
                vol_name  = refine3D_state_vol_fname(state)
                vol_pproc = add2fbody(vol_name, MRC_EXT, PPROC_SUFFIX)
                if( file_exists(vol_pproc) ) call del_file(vol_pproc)
                vol_mirr = add2fbody(vol_pproc, MRC_EXT, MIRR_SUFFIX)
                if( file_exists(vol_mirr) ) call del_file(vol_mirr)
                call vol_name%kill
                call vol_pproc%kill
                call vol_mirr%kill
            enddo
        endif
        call spproj%read_segment('out', projfile)
        call spproj%read_segment('ptcl3D', projfile)
        if( params%l_m_estimator_flex )then
            call wset%new(spproj, spproj%os_ptcl3D)
            if( wset%get_generation() /= wset_generation .or. wset%get_layout_digest() /= wset_digest ) &
                &THROW_HARD('the state weight set changed during the final reconstruction; calc_final_rec')
            if( wset%get_nstates() /= params%nstates ) THROW_HARD('state weight set and nstates disagree; calc_final_rec')
            ! the weights as the backend applied them, over the reconstructed rows
            wthres = 0.
            if( trim(params%rec_backend) == 'pcg' ) wthres = PCG_WEIGHT_THRESHOLD
            rec_rows = pack([(iptcl, iptcl = 1, spproj%os_ptcl3D%get_noris())], spproj%os_ptcl3D%get_all_asint('state') > 0)
        endif
        do state = 1, params%nstates
            if( params%l_m_estimator_flex )then
                pop  = wset%get_pop(state)
                mass = real(wset%applied_mass(state, wthres, rec_rows))
                ess  = real(wset%applied_ess(state, wthres, rec_rows))
                if( mass <= 0. )then
                    call spproj%remove_state_artifacts_from_osout(state)
                    cycle
                endif
            else
                pop  = spproj%os_ptcl3D%get_pop(state, 'state')
                mass = 0.
                ess  = 0.
                if( pop == 0 )then
                    call spproj%remove_state_artifacts_from_osout(state)
                    cycle
                endif
            endif
            call register_state_map(params, spproj, state, smpd, ldim(1), pop, params%l_m_estimator_flex, mass, ess)
        enddo
        call spproj%write_segment_inside('out', projfile)
        call wset%kill
        if( l_postprocess )then
            ! final raw and low-pass diagnostic outputs, then the orthogonal
            ! reprojections for visualization
            call write_final_rec_outputs(params, spproj, lp_snapshot)
            call gen_ortho_reprojs4viz(params, spproj)
        endif
        call cline_final%kill
        call stkname%kill

        contains

            !> Sigmas estimated at a cropped registration box are refreshed at
            !! the final (native) box before the shipped map.
            logical function final_rec_box_changed() result( l_changed )
                integer :: reg_box
                l_changed = .false.
                reg_box   = params%box_crop
                if( cline_refine%defined('box_crop') ) reg_box = cline_refine%get_iarg('box_crop')
                if( reg_box > 0 .and. reg_box /= ldim(1) )then
                    l_changed = .true.
                    write(logfhandle,'(A,I0,A,I0)') &
                        &'>>> FINAL RECONSTRUCTION: registration/final boxes differ; bootstrapping sigmas: ', &
                        &reg_box, ' -> ', ldim(1)
                endif
            end function final_rec_box_changed

            logical function canonical_final_rec_needs_bootstrap() result( l_bootstrap )
                character(len=STDLEN) :: message
                ! one validation boundary for every canonical consumer, at the
                ! final (original) sampling
                l_bootstrap = .not. canonical_sigma2_consumable(spproj, spproj%os_ptcl3D, ldim(1), smpd, message)
                if( l_bootstrap ) write(logfhandle,'(A)') &
                    &'>>> FINAL RECONSTRUCTION: rebuilding canonical sigmas: '//trim(message)
            end function canonical_final_rec_needs_bootstrap




    end subroutine calc_final_rec

    !> solve3D_addon's incremental final reconstruction at native sampling over the run's own particles
    !! (the cohort; entered and left with the working project masked to it by switch_membership) and
    !! the carried frozen sums (frozen_rec, the run context, at the native and the stage box), with no
    !! image-power seed: the bootstrap map on the stage's sigmas, one residual pass that commits the
    !! cohort's sigmas, then the shipped map with the new carried sums (frozen_carry), at the stage box
    !! first when it differs and at the native box last. Between them the union's global means are
    !! republished (copy-on-write sigma2 transaction). No frozen particle is read. The maps are
    !! registered over the union. Contract: doc/policies/3D/solve3D_addon_policy.md sec. 8.
    subroutine calc_final_rec_addon( params, projfile, cline_refine, xrec3D, xrefine3D, &
        &frozen_rec, frozen_carry, stage_box, switch_membership, lp_snapshot )
        use simple_sigma2_bootstrap, only: prepare_residual_sigma2_pass_cline
        class(parameters),     intent(in)    :: params
        class(string),         intent(in)    :: projfile
        class(cmdline),        intent(in)    :: cline_refine
        class(commander_base), intent(inout) :: xrec3D, xrefine3D
        class(string),         intent(in)    :: frozen_rec, frozen_carry
        integer,               intent(in)    :: stage_box
        procedure(addon_membership_switch)   :: switch_membership
        real,                  intent(in)    :: lp_snapshot
        type(sp_project), allocatable :: spproj    !< read from projfile after every membership switch
        type(cmdline) :: cline_rec, cline_sigma
        type(string), allocatable :: seed_vols(:)
        type(string)  :: stkname
        real          :: smpd
        integer       :: ldim(3), nptcls, stkind, ind_in_stk, iter, state, pop
        if( trim(params%rec_backend) /= 'gridding' ) THROW_HARD('solve3D_addon reconstructs by gridding only')
        if( params%l_m_estimator_flex ) THROW_HARD('solve3D_addon carries hard-label sums only; no state weight set')
        if( .not. final_stage_uses_ml_reg(cline_refine) ) &
            &THROW_HARD('solve3D_addon carries euclid ML-regularized sums; the stage is not')
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> SOLVE3D_ADDON: INCREMENTAL RECONSTRUCTION AT ORIGINAL SAMPLING (COHORT + CARRIED SUMS)'
        write(logfhandle,'(A)') '>>>'
        allocate(spproj)
        call spproj%read(projfile)
        call spproj%map_ptcl_ind2stk_ind('ptcl3D', 1, stkind, ind_in_stk)
        stkname = spproj%os_stk%get_str(stkind, 'stk')
        call find_ldim_nptcls(stkname, ldim, nptcls)
        smpd = spproj%os_stk%get(stkind, 'smpd')
        iter = final_rec_bootstrap_sigma_iter(cline_refine)
        ! the stage's sigmas with both means over the union
        call switch_membership(.true., .true.)
        call switch_membership(.false., .false.)
        ! 1. bootstrap map: the carried native sums plus the cohort at its final poses, on the stage's sigmas
        call prep_final_rec_cline(params, projfile, cline_refine, .false., cline_rec, 'bootstrap_rec3D')
        call set_addon_rec_cline(cline_rec, iter)
        write(logfhandle,'(A,I0)') '>>> SOLVE3D_ADDON: BOOTSTRAP MAP (COHORT + CARRIED SUMS), SIGMA ITERATION ', iter
        call xrec3D%execute(cline_rec)
        ! 2. residual pass: the cohort against the bootstrap map at native sampling; commits its sigmas
        allocate(seed_vols(params%nstates))
        do state = 1, params%nstates
            seed_vols(state) = refine3D_state_vol_fname(state)
        enddo
        call prepare_residual_sigma2_pass_cline(cline_rec, iter, params%nstates, seed_vols, cline_sigma)
        write(logfhandle,'(A,I0)') '>>> SOLVE3D_ADDON: RESIDUAL SIGMA2 PASS OF THE COHORT, ITERATION ', iter
        call xrefine3D%execute(cline_sigma)
        call cline_sigma%kill
        call cline_rec%kill
        call switch_membership(.true., .true.)
        call switch_membership(.false., .false.)
        ! 3. the new carried sums: the cohort again on its new sigmas plus the carried sums, at the stage
        ! box when it differs, then the shipped map at the native box last (its outputs are the run's)
        if( stage_box /= ldim(1) )then
            call prep_final_rec_cline(params, projfile, cline_refine, .false., cline_rec, 'reconstruct3D')
            call set_addon_rec_cline(cline_rec, iter + 1, l_carry=.true.)
            call cline_rec%set('box_crop', stage_box)
            write(logfhandle,'(A,I0)') '>>> SOLVE3D_ADDON: CARRIED SUMS AT THE STAGE BOX ', stage_box
            call xrec3D%execute(cline_rec)
            call cline_rec%kill
        endif
        call prep_final_rec_cline(params, projfile, cline_refine, .true., cline_rec, 'reconstruct3D')
        call set_addon_rec_cline(cline_rec, iter + 1, l_carry=.true.)
        write(logfhandle,'(A,I0)') '>>> SOLVE3D_ADDON: SHIPPED MAP AND CARRIED SUMS AT THE NATIVE BOX ', ldim(1)
        call xrec3D%execute(cline_rec)
        call cline_rec%kill
        ! the union, its maps registered
        call switch_membership(.true., .false.)
        call spproj%read(projfile)
        do state = 1, params%nstates
            pop = spproj%os_ptcl3D%get_pop(state, 'state')
            if( pop == 0 )then
                call spproj%remove_state_artifacts_from_osout(state)
                cycle
            endif
            call register_state_map(params, spproj, state, smpd, ldim(1), pop, .false., 0., 0.)
        enddo
        call spproj%write_segment_inside('out', projfile)
        call write_final_rec_outputs(params, spproj, lp_snapshot)
        call gen_ortho_reprojs4viz(params, spproj)
        do state = 1, params%nstates
            call seed_vols(state)%kill
        enddo
        deallocate(seed_vols)
        call stkname%kill
        call spproj%kill
        deallocate(spproj)

    contains

        !> a euclid ML-regularized gridding reconstruct3D with the frozen term, and the carried sums
        subroutine set_addon_rec_cline( cline_here, iter_here, l_carry )
            class(cmdline),    intent(inout) :: cline_here
            integer,           intent(in)    :: iter_here
            logical, optional, intent(in)    :: l_carry
            call cline_here%set('prg',         'reconstruct3D')
            call cline_here%set('oritype',     'ptcl3D')
            call cline_here%set('nstates',     params%nstates)
            call cline_here%set('which_iter',  iter_here)
            call cline_here%set('trail_rec',   'no')
            call cline_here%set('combine_eo',  'no')
            call cline_here%set('objfun',      'euclid')
            call cline_here%set('ml_reg',      'yes')
            call cline_here%set('rec_backend', 'gridding')
            call cline_here%set('frozen_rec',  frozen_rec)
            if( present(l_carry) )then
                if( l_carry ) call cline_here%set('frozen_carry', frozen_carry)
            endif
        end subroutine set_addon_rec_cline

    end subroutine calc_final_rec_addon

    !> one state's final map, FSC and (envfsc) mask into the out segment; under a state weight set
    !! with its hard population, applied mass and effective sample size
    subroutine register_state_map( params, spproj, state, smpd, box, pop, l_weighted, mass, ess )
        class(parameters), intent(in)    :: params
        type(sp_project),  intent(inout) :: spproj
        integer,           intent(in)    :: state, box, pop
        real,              intent(in)    :: smpd, mass, ess
        logical,           intent(in)    :: l_weighted
        type(string) :: str_state, vol_name, vol_envmsk
        integer :: ldim(3)
        logical :: l_mask_exists, l_mask_compatible
        ldim = box
        str_state = int2str_pad(state,2)
        vol_name  = refine3D_state_vol_fname(state)
        if( .not. file_exists(vol_name) )return
        if( l_weighted )then
            call spproj%add_vol2os_out(vol_name, smpd, state, 'vol', pop=pop, mass=mass, ess=ess)
            write(logfhandle,'(A,I0,A,I0,A,F12.3,A,F12.3)') '>>> FINAL RECONSTRUCTION REGISTERED STATE ', state, &
                &': hard population ', pop, ', applied mass ', mass, ', effective sample size ', ess
        else
            call spproj%add_vol2os_out(vol_name, smpd, state, 'vol', pop=pop)
        endif
        call spproj%add_fsc2os_out(refine3D_fsc_fname(state), state, ldim(1))
        if( params%l_envfsc )then
            if( trim(params%automsk) == 'nu' )then
                vol_envmsk = NU_ENVMASK_FBODY//trim(str_state%to_char())//MRC_EXT
                call state_mask_is_compatible(vol_envmsk, ldim(1), smpd, l_mask_exists, l_mask_compatible)
                if( .not. l_mask_compatible )then
                    vol_envmsk = AUTOMASK_FBODY//trim(str_state%to_char())//MRC_EXT
                endif
            else
                vol_envmsk = AUTOMASK_FBODY//trim(str_state%to_char())//MRC_EXT
            endif
            call state_mask_is_compatible(vol_envmsk, ldim(1), smpd, l_mask_exists, l_mask_compatible)
            if( l_mask_compatible )then
                call spproj%add_vol2os_out(vol_envmsk, smpd, state, 'vol_msk', ldim(1))
            else if( l_mask_exists )then
                THROW_WARN('>>> FINAL RECONSTRUCTION: envfsc mask has incompatible dimensions or sampling for state '//str_state%to_char())
            else
                THROW_WARN('>>> FINAL RECONSTRUCTION: expected envfsc mask file does not exist for state '//str_state%to_char())
            endif
        endif
    end subroutine register_state_map

    !> The child command line is built from scratch: only the controls
    !! that define the final reconstruction are copied from params and
    !! from the refinement command line.
    subroutine prep_final_rec_cline( params, projfile, cline_refine, l_postprocess, child_cline, prg )
        class(parameters), intent(in)    :: params
        class(string),     intent(in)    :: projfile
        class(cmdline),    intent(in)    :: cline_refine
        logical,           intent(in)    :: l_postprocess
        class(cmdline),    intent(inout) :: child_cline
        character(len=*),  intent(in)    :: prg
        integer :: istate
        call child_cline%kill
        call child_cline%set('prg',      prg)
        call child_cline%set('mkdir',    'no')
        call child_cline%set('projfile', projfile)
        ! volassemble appends _STATENN and writes the extension-less
        ! resolution document next to rec_final_stateNN.mrc.
        call child_cline%set('outfile', 'RESOLUTION_FINAL.txt')
        call child_cline%set('pgrp',    params%pgrp)
        call child_cline%set('envfsc',   params%envfsc)
        call child_cline%set('envmsklp', params%envmsklp)
        call child_cline%set('binwidth', params%binwidth)
        if( params%nthr    > 1  ) call child_cline%set('nthr',    params%nthr)
        if( params%mskdiam > 0. ) call child_cline%set('mskdiam', params%mskdiam)
        if( params%nparts  > 1  ) call child_cline%set('nparts',  params%nparts)
        if( params%nstates > 1  ) call child_cline%set('nstates', params%nstates)
        ! the reconstruction's particle weighting and state restriction follow the refinement
        if( cline_refine%defined('m_estimator') ) &
            &call child_cline%set('m_estimator', cline_refine%get_carg('m_estimator'))
        if( cline_refine%defined('state') ) call child_cline%set('state', cline_refine%get_iarg('state'))
        ! automsk is inherited on BOTH routes: the shipped PCG map is estimated on the same
        ! density/NU support, with the same estimator-constrained FSC, as every refinement iteration.
        if( cline_refine%defined('automsk') ) &
            &call child_cline%set('automsk',   cline_refine%get_carg('automsk'))
        if( trim(params%automsk) == 'nu' )then
            ! automsk=nu is valid only in an NU filter mode, including
            ! internal final reconstruction and assembly commands.
            call child_cline%set('filt_mode', params%filt_mode)
            if( cline_refine%defined('nu_msk_sig') ) &
                &call child_cline%set('nu_msk_sig', cline_refine%get_rarg('nu_msk_sig'))
            if( cline_refine%defined('amsklp') ) &
                &call child_cline%set('amsklp', cline_refine%get_rarg('amsklp'))
        endif
        if( prg.eq.'bootstrap_rec3D' )then
            ! the residual sigmas depend on the regularization of the
            ! reference they are scored against: the bootstrap map is
            ! regularized exactly as the refinement's matching
            ! references were; the shipped map is made classical by
            ! bootstrap_rec3D itself
            if( cline_refine%defined('filt_mode') ) &
                &call child_cline%set('filt_mode', cline_refine%get_carg('filt_mode'))
            if( trim(params%automsk) /= 'nu' .and. cline_refine%defined('nu_msk_sig') ) &
                &call child_cline%set('nu_msk_sig', cline_refine%get_rarg('nu_msk_sig'))
        else
            ! the direct route ships a classical map (no nonuniform
            ! filtering, a matching-reference feature); on PCG the
            ! last refinement volume at this box is the lag-one
            ! reference the solve support derives from
            if( trim(params%automsk) /= 'nu' ) call child_cline%set('filt_mode', 'none')
            if( trim(params%rec_backend) == 'pcg' )then
                do istate = 1, params%nstates
                    if( file_exists(refine3D_state_vol_fname(istate)) ) &
                        &call child_cline%set('vol'//int2str(istate), refine3D_state_vol_fname(istate))
                enddo
            endif
        endif
        if( .not. l_postprocess )then
            call child_cline%set('postprocess', 'no')
        endif
        if( prg.eq.'reconstruct3D' .and. .not. final_stage_uses_ml_reg(cline_refine) )then
            call child_cline%set('objfun', 'cc')
            call child_cline%set('ml_reg', 'no')
        endif
        ! the final reconstruction runs on the refinement's backend;
        ! final-map postprocessing is classical on both backends
        if( trim(params%rec_backend) == 'pcg' )then
            call child_cline%set('rec_backend', 'pcg')
            call configure_final_pcg_solve_budget(cline_refine, child_cline)
            if( cline_refine%defined('rtol') )&
                &call child_cline%set('rtol', cline_refine%get_rarg('rtol'))
            ! the soft solvent prior of the refinement's PCG stage
            ! continues into the shipped map
            if( cline_refine%defined('pcg_solvent') )then
                call child_cline%set('pcg_solvent', cline_refine%get_carg('pcg_solvent'))
                if( cline_refine%defined('pcg_solvent_lambda') ) &
                    &call child_cline%set('pcg_solvent_lambda', cline_refine%get_rarg('pcg_solvent_lambda'))
                if( cline_refine%defined('pcg_solvent_check') ) &
                    &call child_cline%set('pcg_solvent_check', cline_refine%get_carg('pcg_solvent_check'))
            endif
        endif
    end subroutine prep_final_rec_cline

    logical function final_stage_uses_ml_reg( cline_refine ) result( l_ml_reg )
        class(cmdline), intent(in) :: cline_refine
        l_ml_reg = .false.
        if( .not. cline_refine%defined('ml_reg') ) return
        if( cline_refine%get_carg('ml_reg').ne.'yes' ) return
        if( cline_refine%defined('objfun') )then
            l_ml_reg = cline_refine%get_carg('objfun').eq.'euclid'
        else
            l_ml_reg = .true.
        endif
    end function final_stage_uses_ml_reg

    !> bootstrap_rec3D's residual sigma pass is a refine3D iteration;
    !! number it beyond the refinement's own iterations so its
    !! iteration files never collide with theirs. The canonical sigma
    !! state is one committed file and carries no iteration number.
    integer function final_rec_bootstrap_sigma_iter( cline_refine ) result( iter )
        class(cmdline), intent(in) :: cline_refine
        iter = 1
        if( cline_refine%defined('endit') )then
            iter = cline_refine%get_iarg('endit') + 2
        else if( cline_refine%defined('which_iter') )then
            iter = cline_refine%get_iarg('which_iter') + 2
        endif
        iter = max(1, iter)
    end function final_rec_bootstrap_sigma_iter

end module simple_final_rec

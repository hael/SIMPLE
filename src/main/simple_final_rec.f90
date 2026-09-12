!@descr: the shared final all-particle reconstruction at native sampling
module simple_final_rec
use simple_commanders_api
use simple_parameters,       only: parameters
use simple_abinitio_utils,   only: configure_final_pcg_solve_budget, write_final_rec_outputs, &
    &gen_ortho_reprojs4viz
use simple_refine3D_fnames,  only: refine3D_fsc_fname, refine3D_state_vol_fname
use simple_vol_pproc_policy, only: state_mask_is_compatible
use simple_sigma2_files,     only: canonical_sigma2_consumable
implicit none
#include "simple_local_flags.inc"

public :: calc_final_rec
private

contains

    !> The one ending shared by abinitio3D, refine3D_auto, refine3D_states and
    !! classify3D_refs: the final all-particle reconstruction at the native
    !! particle sampling, driven by the refinement command line the map
    !! continues from (filtering, automask, backend, crop and iteration
    !! provenance). Valid committed canonical sigmas are reused; a missing,
    !! stale, wrong-grid or crop-box-registered state is rebuilt by
    !! bootstrap_rec3D (image-power seed, euclid ML bootstrap map, one residual
    !! sigma pass, shipped euclid ML map). A refinement without ML
    !! regularization ships a classical correlation map. State maps, FSCs and
    !! envfsc masks are registered in the project; with l_postprocess the final
    !! products (rec_final maps, low-pass snapshots, postprocessed and mirrored
    !! maps, orthogonal reprojections) are written as well.
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
        type(cmdline) :: cline_final
        type(string)  :: str_state, vol_name, stkname, vol_pproc, vol_mirr, vol_envmsk
        integer       :: ldim(3), state, pop, stkind, ind_in_stk, nptcls, bootstrap_sigma_iter
        real          :: smpd
        logical       :: l_bootstrap_sigmas, l_mask_exists, l_mask_compatible
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call spproj%read(projfile) ! ensure we have the latest project info
        call spproj%map_ptcl_ind2stk_ind('ptcl3D', 1, stkind, ind_in_stk)
        stkname = spproj%os_stk%get_str(stkind, 'stk')
        call find_ldim_nptcls(stkname, ldim, nptcls)
        smpd = spproj%os_stk%get(stkind, 'smpd')
        write(logfhandle,'(A,I0,A,F8.4)') '>>> FINAL RECONSTRUCTION SAMPLING: box=', ldim(1), ' smpd=', smpd
        call prep_final_rec_cline(cline_final, 'reconstruct3D')
        l_bootstrap_sigmas = .false.
        if( final_stage_uses_ml_reg() )then
            ! a valid committed state is reused directly; a missing, stale,
            ! wrong-grid, wrong-layout or wrong-grouping state is rebuilt from
            ! particle power and residual-upgraded by bootstrap_rec3D.
            ! Sigmas estimated at a cropped registration box are refreshed at
            ! native sampling before the shipped map (2026-09-07).
            l_bootstrap_sigmas = canonical_final_rec_needs_bootstrap()
            if( .not. l_bootstrap_sigmas ) l_bootstrap_sigmas = final_rec_box_changed()
            if( .not. l_bootstrap_sigmas ) write(logfhandle,'(A)') &
                &'>>> FINAL RECONSTRUCTION: reusing committed canonical sigmas'
        endif
        if( l_bootstrap_sigmas )then
            bootstrap_sigma_iter = final_rec_bootstrap_sigma_iter()
            call prep_final_rec_cline(cline_final, 'bootstrap_rec3D')
            call cline_final%set('which_iter', bootstrap_sigma_iter)
            write(logfhandle,'(A,I0)') '>>> FINAL RECONSTRUCTION BOOTSTRAP SIGMA ITERATION: ', bootstrap_sigma_iter
            if( trim(params%rec_backend) == 'pcg' ) write(logfhandle,'(A,I0)') &
                &'>>> FINAL PCG COLD-SOLVE ITERATION BUDGET: ', cline_final%get_iarg('maxits_pcg')
            ! bootstrap_rec3D owns the complete sequence: image-power seed,
            ! euclid ML bootstrap map, one residual sigma2 pass (refine=sigma)
            ! against it at the final sampling, group consolidation as the
            ! next iteration and the shipped euclid ML reconstruction on the
            ! residual sigmas. The same program is the standalone test entry
            ! point for this stage on any project with 3D orientations
            ! (simple_exec prg=bootstrap_rec3D), 2026-09-07.
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
        do state = 1, params%nstates
            pop = spproj%os_ptcl3D%get_pop(state, 'state')
            if( pop == 0 )cycle     ! empty-state case
            str_state = int2str_pad(state,2)
            vol_name  = refine3D_state_vol_fname(state)
            if( .not. file_exists(vol_name) )cycle
            call spproj%add_vol2os_out(vol_name, smpd, state, 'vol', pop=pop)
            call spproj%add_fsc2os_out(refine3D_fsc_fname(state), state, ldim(1))
            if( params%l_envfsc )then
                vol_envmsk = AUTOMASK_FBODY//trim(str_state%to_char())//MRC_EXT
                call state_mask_is_compatible(vol_envmsk, ldim(1), smpd, l_mask_exists, l_mask_compatible)
                if( l_mask_compatible )then
                    call spproj%add_vol2os_out(vol_envmsk, smpd, state, 'vol_msk', ldim(1))
                else if( l_mask_exists )then
                    THROW_WARN('>>> FINAL RECONSTRUCTION: envfsc mask has incompatible dimensions or sampling for state '//str_state%to_char())
                else
                    THROW_WARN('>>> FINAL RECONSTRUCTION: expected envfsc mask file does not exist for state '//str_state%to_char())
                endif
            endif
        enddo
        call spproj%write_segment_inside('out', projfile)
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
                l_bootstrap = .not. canonical_sigma2_consumable(spproj, spproj%os_ptcl3D, ldim(1), smpd, &
                    &params%l_sigma_glob, message)
                if( l_bootstrap ) write(logfhandle,'(A)') &
                    &'>>> FINAL RECONSTRUCTION: rebuilding canonical sigmas: '//trim(message)
            end function canonical_final_rec_needs_bootstrap

            !> bootstrap_rec3D's residual sigma pass is a refine3D iteration;
            !! number it beyond the refinement's own iterations so its
            !! iteration files never collide with theirs. The canonical sigma
            !! state is one committed file and carries no iteration number.
            integer function final_rec_bootstrap_sigma_iter() result( iter )
                iter = 1
                if( cline_refine%defined('endit') )then
                    iter = cline_refine%get_iarg('endit') + 2
                else if( cline_refine%defined('which_iter') )then
                    iter = cline_refine%get_iarg('which_iter') + 2
                endif
                iter = max(1, iter)
            end function final_rec_bootstrap_sigma_iter

            !> The child command line is built from scratch: only the controls
            !! that define the final reconstruction are copied from params and
            !! from the refinement command line.
            subroutine prep_final_rec_cline( child_cline, prg )
                class(cmdline), intent(inout) :: child_cline
                character(len=*), intent(in)  :: prg
                call child_cline%kill
                call child_cline%set('prg',      prg)
                call child_cline%set('mkdir',    'no')
                call child_cline%set('projfile', projfile)
                call child_cline%set('sigma_est',   params%sigma_est)
                ! volassemble appends _STATENN and writes the extension-less
                ! resolution document next to rec_final_stateNN.mrc.
                call child_cline%set('outfile', 'RESOLUTION_FINAL.txt')
                call child_cline%set('pgrp',    params%pgrp)
                call child_cline%set('ptcl_src', params%ptcl_src)
                call child_cline%set('envfsc',   params%envfsc)
                call child_cline%set('envmsklp', params%envmsklp)
                call child_cline%set('binwidth', params%binwidth)
                if( params%nthr    > 1  ) call child_cline%set('nthr',    params%nthr)
                if( params%mskdiam > 0. ) call child_cline%set('mskdiam', params%mskdiam)
                if( params%nparts  > 1  ) call child_cline%set('nparts',  params%nparts)
                if( params%nstates > 1  ) call child_cline%set('nstates', params%nstates)
                if( final_stage_uses_ml_reg() ) call child_cline%set('conical_fsc', params%conical_fsc)
                if( prg.eq.'bootstrap_rec3D' )then
                    ! the residual sigmas depend on the regularization of the
                    ! reference they are scored against: the bootstrap map is
                    ! regularized exactly as the refinement's matching
                    ! references were; the shipped map is made classical by
                    ! bootstrap_rec3D itself (2026-09-07)
                    if( cline_refine%defined('filt_mode') ) &
                        &call child_cline%set('filt_mode', cline_refine%get_carg('filt_mode'))
                    if( cline_refine%defined('nu_refine') ) &
                        &call child_cline%set('nu_refine', cline_refine%get_carg('nu_refine'))
                    if( cline_refine%defined('automsk') ) &
                        &call child_cline%set('automsk',   cline_refine%get_carg('automsk'))
                    if( cline_refine%defined('nu_msk_sig') ) &
                        &call child_cline%set('nu_msk_sig', cline_refine%get_rarg('nu_msk_sig'))
                endif
                if( .not. l_postprocess )then
                    call child_cline%set('postprocess', 'no')
                endif
                if( prg.eq.'reconstruct3D' .and. .not. final_stage_uses_ml_reg() )then
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
                endif
            end subroutine prep_final_rec_cline

            logical function final_stage_uses_ml_reg() result( l_ml_reg )
                l_ml_reg = .false.
                if( .not. cline_refine%defined('ml_reg') ) return
                if( cline_refine%get_carg('ml_reg').ne.'yes' ) return
                if( cline_refine%defined('objfun') )then
                    l_ml_reg = cline_refine%get_carg('objfun').eq.'euclid'
                else
                    l_ml_reg = .true.
                endif
            end function final_stage_uses_ml_reg

    end subroutine calc_final_rec

end module simple_final_rec

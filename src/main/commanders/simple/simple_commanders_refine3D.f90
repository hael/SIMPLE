!@descr: supporting 3D orientation search
module simple_commanders_refine3D
use simple_commanders_api
use simple_pftc_srch_api
use simple_refine3D_fnames,   only: refine3D_state_vol_fname, refine3D_fsc_fname
use simple_refine3D_stage_plan, only: refine3D_stage_plan_entry, plan_refine3D_frequency_stages
use simple_external_reference_pose_initialization, only: initialize_poses_against_external_references
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander

type, extends(commander_base) :: commander_refine3D_auto
  contains
    procedure :: execute      => exec_refine3D_auto
end type commander_refine3D_auto

type, extends(commander_base) :: commander_refine3D_states
  contains
    procedure :: execute      => exec_refine3D_states
end type commander_refine3D_states

type, extends(commander_base) :: commander_classify3D_refs
  contains
    procedure :: execute      => exec_classify3D_refs
end type commander_classify3D_refs

type, extends(commander_base) :: commander_refine3D
  contains
    procedure :: execute      => exec_refine3D
end type commander_refine3D

type, extends(commander_base) :: commander_refine3D_distr_worker
  contains
    procedure :: execute      => exec_refine3D_distr_worker
end type commander_refine3D_distr_worker

contains

    subroutine exec_nspace(self,cline)
        class(nspace_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(oris)       :: o
        real             :: ares
        integer          :: i
        call params%new(cline)
        do i=500,5000,500
            o = oris(i, is_ptcl=.false.)
            call o%spiral
            ares = o%find_angres()
            write(logfhandle,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, params%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_refine3D_auto( self, cline )
        use simple_abinitio_utils, only: gen_ortho_reprojs4viz, write_final_rec_outputs
        use simple_commanders_rec, only: commander_rec3D
        use simple_estimate_ssnr,     only: lpstages_setlims
        use simple_commanders_rec,    only: commander_bootstrap_rec3D
        use simple_commanders_euclid, only: commander_calc_pspec
        use simple_refine3D_strategy, only: strip_refine3D_search_only_args
        class(commander_refine3D_auto), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(cmdline)               :: cline_rec3D, cline_boot
        type(parameters)            :: params, params_final_rec
        type(sp_project)            :: spproj
        type(string)                :: init_vol
        type(string)                :: pose_init_refs(1), pose_init_checkpoint(1)
        integer, parameter :: NSAMPLE_REFINE3D_AUTO = 25000
        real,    parameter :: SMPD_TARGET_DEFAULT = 1.3
        real,    parameter :: TARGET_UPDATES_PER_PARTICLE_REFINE3D_AUTO = 4.0
        character(len=*), parameter :: WORKFLOW_LABEL = 'REFINE3D_AUTO'
        logical, parameter :: DEBUG  = .true.
        integer, parameter :: MINBOX = 256
        integer, parameter :: MINITS_REFINE3D_AUTO = 10
        integer, parameter :: MAXITS_REFINE3D_AUTO_CAP = 50
        real    :: smpd_target, smpd_crop, scale, trslim, init_smpd, update_frac_auto
        integer :: box_crop, init_box, nptcls_eff, nsample_target, maxits_user
        logical :: l_autoscale, l_have_init_vol, l_maxits_defined
        logical :: l_external_input, l_ref_pose_init_requested
        ! commanders
        type(commander_rec3D)           :: xrec3D
        type(commander_bootstrap_rec3D) :: xbootstrap_rec3D
        type(commander_calc_pspec)      :: xcalc_pspec
        type(commander_refine3D)        :: xrefine3D
        maxits_user      = 0
        l_external_input = cline%defined('vol1')
        ! hard defaults
        call cline%set('balance',         'no') ! no balancing based on 2D clustering
        call cline%set('greedy_sampling', 'no') ! only active when balance is 'yes'`
        call cline%set('trail_rec',      'yes') ! trailing average 3D reconstruction
        call cline%set('refine',  'prob_neigh') ! probabilistic neighborhood 3D refinement
        call cline%set('ml_reg',         'yes') ! ML regularization is on
        call cline%set('overlap',         0.99) ! convergence if overlap > 99%
        call cline%set('nstates',            1) ! only single-state refinement is supported
        call cline%set('objfun',      'euclid') ! the objective function is noise-normalized Euclidean distance
        call cline%set('lplim_crit',     0.143) ! we use the 0.143 criterion for low-pass limitation
        call cline%set('incrreslim',      'no') ! if anything 'yes' makes it slightly worse, but no real difference right now
        ! overridable defaults
        if( .not. cline%defined('envfsc') )then
            ! Density-envelope/phase-randomization-corrected FSC steers the
            ! matching LP, NU gating, and convergence.
            call cline%set('envfsc', 'yes')
        endif
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',            'no') ! 4 now, probably fine
        if( cline%defined('ref_pose_init') .and. cline%get_carg('ref_pose_init').eq.'cc' )then
            call cline%set('sigma_est', 'global')
        else if( .not. cline%defined('sigma_est') )then
            call cline%set('sigma_est', 'global') ! 4 now, probably fine
        endif
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',        'no') ! 4 now, to allow more rapid testing
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',        'yes') ! no difference at this stage, so prefer 'yes'
        if( .not. cline%defined('nsample')     ) call cline%set('nsample', NSAMPLE_REFINE3D_AUTO)
        if( .not. cline%defined('autoscale')   ) call cline%set('autoscale',        'yes')
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode', 'nonuniform') ! obvioulsy
        ! nu_refine=yes means conservative resolution-bank expansion on BOTH
        ! backends (Stage 6.6): the gridding filter challenger, or its mirror
        ! on the pcg path -- Q_NU evidence-bank shell extension with the same
        ! proven win-fraction acceptance. abinitio3D keeps the discrete
        ! static ladder via its stage policy.
        if( .not. cline%defined('nu_refine')   ) call cline%set('nu_refine',        'yes') ! allow conservative NU resolution-bank expansion
        if( .not. cline%defined('automsk') )then
            call cline%set('automsk', 'yes') ! evidence-constrained background filtering
        endif
        l_maxits_defined = cline%defined('maxits')
        if( l_maxits_defined )then
            maxits_user = cline%get_iarg('maxits')
            if( maxits_user < 1 ) THROW_HARD('maxits must be >= 1 for '//WORKFLOW_LABEL)
            call cline%set('minits', maxits_user)
        else if( cline%defined('minits') )then
            call cline%set('minits', max(MINITS_REFINE3D_AUTO, cline%get_iarg('minits')))
        else
            call cline%set('minits', MINITS_REFINE3D_AUTO)
        endif
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol', 'no') ! we do not keep volumes for each iteration by deafult
        call params%new(cline)
        l_ref_pose_init_requested = trim(params%ref_pose_init).eq.'cc'
        if( l_ref_pose_init_requested .and. .not. l_external_input )then
            THROW_HARD(WORKFLOW_LABEL//' ref_pose_init=cc requires an external vol1 input')
        endif
        call cline%set('mkdir', 'no') ! to avoid nested directory structure
        call set_refine3D_auto_sampling()
        l_have_init_vol = .false.
        init_box        = 0
        init_smpd       = 0.
        if( cline%defined('vol1') )then
            init_vol = cline%get_carg('vol1')
            if( .not. file_exists(init_vol) )then
                THROW_HARD('File: '//init_vol%to_char()//' does not exist! '//WORKFLOW_LABEL)
            endif
            call prepare_external_init_vol(init_vol)
            l_have_init_vol = .true.
            call cline%delete('vol1') ! because now part of the project
            write(logfhandle,'(A,1X,A)') '>>> '//WORKFLOW_LABEL//' USING INPUT VOLUME:', init_vol%to_char()
            if( .not. l_ref_pose_init_requested )then
                write(logfhandle,'(A)') &
                    &'>>> WARNING: EXTERNAL VOL1 IS TRUSTED WITHOUT CC POSE INITIALIZATION (REF_POSE_INIT=NONE)'
            endif
        else
            call spproj%read_segment('out', params%projfile)
            if( spproj%isthere_in_osout('vol', 1) )then
                call spproj%get_vol('vol', 1, init_vol, init_smpd, init_box)
                if( file_exists(init_vol) )then
                    if( project_init_vol_compatible() )then
                        l_have_init_vol = .true.
                        write(logfhandle,'(A,1X,A)') '>>> '//WORKFLOW_LABEL//' USING PROJECT VOLUME:', init_vol%to_char()
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> PROJECT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
                    else
                        write(logfhandle,'(A,1X,A)') &
                            &'>>> '//WORKFLOW_LABEL//' PROJECT VOLUME SAMPLING MISMATCH, RECONSTRUCTING:', &
                            &init_vol%to_char()
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> PROJECT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> CURRENT RUN BOX/SMPD:    ', params%box, '/', params%smpd
                    endif
                else
                    write(logfhandle,'(A,1X,A)') &
                        &'>>> '//WORKFLOW_LABEL//' PROJECT VOLUME MISSING, RECONSTRUCTING:', init_vol%to_char()
                endif
            endif
            call spproj%kill
        endif
        smpd_target = SMPD_TARGET_DEFAULT
        if( .not. params%l_autoscale .or. params%box <= MINBOX )then
            smpd_target = params%smpd
            smpd_crop   = params%smpd
            box_crop    = params%box
            scale       = 1.0
            l_autoscale = .false.
        else
            call autoscale(params%box, params%smpd, smpd_target, box_crop, smpd_crop, scale, minbox=MINBOX)
            l_autoscale = box_crop < params%box
        endif
        trslim = min(8.,max(2.0, AHELIX_WIDTH / smpd_crop))
        if( DEBUG )then
            print *, 'smpd_target: ', smpd_target
            print *, 'box:         ', params%box
            print *, 'box_crop:    ', box_crop
            print *, 'smpd:        ', params%smpd
            print *, 'smpd_crop:   ', smpd_crop
            print *, 'scale:       ', scale
            print *, 'trslim:      ', trslim
            print *, 'l_autoscale: ', l_autoscale
        endif
        call cline%set('trs', trslim)
        if( l_autoscale )then
            call cline%set('box_crop',  box_crop)
            call cline%set('smpd_crop', smpd_crop)
        else
            call cline%delete('box_crop')
            call cline%delete('smpd_crop')
        endif
        if( l_ref_pose_init_requested ) call initialize_external_reference_poses
        ! generate an initial 3D reconstruction
        cline_rec3D = cline
        call cline_rec3D%set('prg', 'reconstruct3D') ! required for distributed call
        call cline_rec3D%delete('trail_rec')
        call cline_rec3D%delete('objfun')
        call cline_rec3D%delete('objfun_den')
        call cline_rec3D%delete('objfun_den_w')
        call cline_rec3D%delete('sigma_est')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%delete('box_crop')           ! original image dimensions
        call cline_rec3D%delete('smpd_crop')
        call cline_rec3D%set('objfun', 'cc') ! ugly, but this is how it works in parameters
        call cline_rec3D%set('postprocess', 'no')
        call cline_rec3D%set('nu_refine', 'no')
        ! the initial bootstrap reconstruction runs objfun=cc, so the Q_NU
        ! replay cannot engage there: strip the pinning keys or the PCG
        ! validator hard-errors on the explicit-activation contract
        ! (pcg_priors.md S6.2). The final reconstruction re-forwards them
        ! explicitly below, where its regularized bootstrap pass can honor
        ! them.
        call cline_rec3D%delete('pcg_nu_lambda_rel')
        call cline_rec3D%delete('pcg_nu_supp_target')
        ! STARTUP BOOTSTRAP: reconstruct -> build masks -> re-reconstruct with
        ! the masks and the NU prior, before any matching. Without it the
        ! first iteration matches raw, spherically masked references while
        ! every later iteration matches envelope-masked, NU-filtered ones --
        ! the reference convention changes underneath an already converged
        ! alignment, which is where refine3D_auto has been losing particles.
        ! Sigmas come FIRST, from the particle power spectra, so the startup
        ! reconstruction and every refinement iteration share one sigma
        ! basis. Deriving them instead from the startup half maps (what
        ! bootstrap_rec3D does, for the case where no box-compatible sigmas
        ! can exist) leaves the startup regularized against sigmas the
        ! refinement then discards, and its heavy rescaling conditions the
        ! euclid system markedly worse (bgal: base residual 0.23 vs 0.08).
        ! refine3D reuses these rather than re-deriving them, see
        ! sigma2_stage_needs_bootstrap.
        cline_boot = cline
        call strip_refine3D_search_only_args(cline_boot)
        call cline_boot%set('prg', 'calc_pspec')
        call cline_boot%set('which_iter', 1)
        call xcalc_pspec%execute(cline_boot)
        ! With the sigmas in hand the startup is a SINGLE regularized
        ! reconstruction carrying the refinement's own filtering settings; it
        ! leaves behind the automask, the NU evidence envelope, the _nu_filt
        ! matching references and the matching-lp handoff that iteration 1
        ! then consumes.
        cline_boot = cline
        call cline_boot%set('prg', 'reconstruct3D')
        call cline_boot%delete('trail_rec')
        call cline_boot%delete('objfun_den')
        call cline_boot%delete('objfun_den_w')
        call cline_boot%delete('update_frac')
        call cline_boot%delete('endit')
        ! alignment-loop controls have no meaning for a reconstruction
        call cline_boot%delete('refine')
        call cline_boot%delete('maxits')
        call cline_boot%delete('minits')
        call cline_boot%delete('continue')
        call cline_boot%set('postprocess', 'no')
        call cline_boot%set('which_iter',      1)
        call xrec3D%execute(cline_boot)
        ! the bootstrap reconstruction is now the starting reference, and its
        ! NU-filtered halves, masks and matching-lp handoff are on disk and in
        ! the project, so iteration 1 matches exactly what iteration N will
        call cline%set('vol1', refine3D_state_vol_fname(1))
        ! the NU replay controllers compare consecutive REFINEMENT iterations.
        ! The bootstrap is a reconstruction from the incoming orientations, so
        ! leaving its readout behind makes the first auto-target comparison
        ! bootstrap-vs-iteration-1, which is not like for like and produced a
        ! spurious setpoint step. Start the controllers clean.
        if( file_exists(PCG_NU_STATS_FILE) ) call del_file(PCG_NU_STATS_FILE)
        call cline_boot%kill
        call seed_refine3D_auto_nonuniform_lpset()
        ! 3D refinement iterations
        call cline%set('prg',                   'refine3D')
        call cline%set('maxits',             params%maxits)
        call xrefine3D%execute(cline)
        ! re-reconstruct from all particle images at original sampling: the
        ! refinement sigmas are crop-box incompatible, so bootstrap_rec3D
        ! derives compatible sigmas from an unregularized pass and ships a
        ! euclid ML-regularized final map -- with the Q_NU replay in-solve
        ! on the pcg backend
        call cline_rec3D%set('prg', 'bootstrap_rec3D')
        call cline_rec3D%set('outfile', 'RESOLUTION_FINAL.txt')
        call cline_rec3D%set('postprocess', 'yes')
        call cline_rec3D%delete('objfun') ! bootstrap_rec3D owns objfun/ml_reg per pass
        call cline_rec3D%delete('ml_reg')
        if( cline%defined('endit') )then
            ! write bootstrap sigmas beyond the refinement's own iterations
            ! so no crop-box sigma star is overwritten
            call cline_rec3D%set('which_iter', cline%get_iarg('endit') + 2)
        else
            call cline_rec3D%set('which_iter', MAXITS_REFINE3D_AUTO_CAP + 2)
        endif
        if( params%l_nonuniform )then
            if( trim(params%rec_backend) == 'pcg' )then
                ! keep the NU machinery: the Q_NU prior regularizes the final
                ! map in-solve (auto-lambda resumes from the persisted stats
                ! file; explicit keys pin)
                call cline_rec3D%set('filt_mode', 'nonuniform')
                if( cline%defined('pcg_nu_lambda_rel') )&
                    &call cline_rec3D%set('pcg_nu_lambda_rel', cline%get_rarg('pcg_nu_lambda_rel'))
                if( cline%defined('pcg_nu_supp_target') )&
                    &call cline_rec3D%set('pcg_nu_supp_target', cline%get_rarg('pcg_nu_supp_target'))
            else
                call cline_rec3D%set('filt_mode', 'none')
            endif
            call cline_rec3D%set('automsk', 'no')
        endif
        call cline_rec3D%set('nu_refine', 'no')
        call xbootstrap_rec3D%execute(cline_rec3D)
        call params_final_rec%new(cline_rec3D)
        params_final_rec%box  = params%box
        params_final_rec%smpd = params%smpd
        call spproj%read_segment('out', params_final_rec%projfile)
        call write_final_rec_outputs(params_final_rec, spproj, params_final_rec%res_target)
        call spproj%kill
        call init_vol%kill

    contains

        subroutine seed_refine3D_auto_nonuniform_lpset()
            type(sp_project) :: seed_proj
            type(string) :: fsc_fname
            real, allocatable :: fsc(:), res(:)
            real :: fsc05, fsc0143
            integer :: fsc_box
            if( .not. params%l_nonuniform_lpset ) return
            if( cline%defined('lp') ) return
            fsc_box = 0
            call seed_proj%read_segment('out', params%projfile)
            if( seed_proj%isthere_in_osout('fsc', 1) )then
                call seed_proj%get_fsc(1, fsc_fname, fsc_box)
            else
                fsc_fname = refine3D_fsc_fname(1)
                fsc_box = params%box
            endif
            call seed_proj%kill
            if( .not. file_exists(fsc_fname) )then
                THROW_HARD(WORKFLOW_LABEL//' filt_mode=nonuniform_lpset requires starting-volume FSC metadata or explicit lp')
            endif
            fsc = file2rarr(fsc_fname)
            if( size(fsc) < 2 ) THROW_HARD('starting-volume FSC has too few shells; '//WORKFLOW_LABEL//' filt_mode=nonuniform_lpset')
            if( fsc_box < 1 ) fsc_box = params%box
            res = get_resarr(fsc_box, params%smpd)
            if( size(res) < size(fsc) ) THROW_HARD('starting-volume FSC/box size mismatch; '//WORKFLOW_LABEL)
            call get_resolution(fsc, res, fsc05, fsc0143)
            if( fsc0143 <= TINY )then
                THROW_HARD(WORKFLOW_LABEL//' filt_mode=nonuniform_lpset could not derive a positive starting resolution; set lp explicitly')
            endif
            params%lp         = fsc0143
            params%kfromto(2) = calc_fourier_index(params%lp, params%box, params%smpd)
            params%l_lpset    = .true.
            call cline%set('lp', params%lp)
            write(logfhandle,'(A,F8.3,A)') &
                &'>>> '//WORKFLOW_LABEL//' nonuniform_lpset seeded matching low-pass from starting FSC 0.143: ', &
                &params%lp, ' A'
            if( allocated(fsc) ) deallocate(fsc)
            if( allocated(res) ) deallocate(res)
            call fsc_fname%kill
        end subroutine seed_refine3D_auto_nonuniform_lpset

        logical function project_init_vol_compatible() result( l_compatible )
            l_compatible = init_box == params%box .and. init_smpd > TINY .and. &
                &abs(init_smpd - params%smpd) <= 1.e-6
        end function project_init_vol_compatible

        subroutine set_refine3D_auto_sampling()
            type(sp_project) :: sampling_proj
            integer :: maxits_auto, nptcls_per_iter
            nsample_target = params%nsample
            if( nsample_target < 1 ) THROW_HARD('nsample must be >= 1 for '//WORKFLOW_LABEL)
            call sampling_proj%read(params%projfile)
            nptcls_eff = sampling_proj%count_state_gt_zero()
            call sampling_proj%kill
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//WORKFLOW_LABEL)
            nptcls_per_iter = min(nptcls_eff, nsample_target)
            if( nptcls_eff <= nsample_target )then
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)') &
                    &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    &nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
            else
                update_frac_auto = real(nsample_target) / real(nptcls_eff)
                if( update_frac_auto <= 0.99 )then
                    params%update_frac   = update_frac_auto
                    params%l_update_frac = .true.
                    params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                    call cline%set('update_frac', update_frac_auto)
                    write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                        &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                        &nptcls_eff, '/', nsample_target, '/', update_frac_auto
                else
                    params%update_frac   = 1.0
                    params%l_update_frac = .false.
                    params%l_trail_rec   = .false.
                    call cline%delete('update_frac')
                    write(logfhandle,'(A,I0,A,I0,A)') &
                        &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                        &nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
                endif
            endif
            if( .not. l_maxits_defined )then
                maxits_auto = ceiling((TARGET_UPDATES_PER_PARTICLE_REFINE3D_AUTO * real(nptcls_eff)) / real(nptcls_per_iter))
                maxits_auto = max(params%minits, min(MAXITS_REFINE3D_AUTO_CAP, max(2, maxits_auto)))
                params%maxits = maxits_auto
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0,A,F5.1,A,I0,A)') '>>> '//WORKFLOW_LABEL//' MAXITS: ', &
                    &params%maxits, ' FOR ~', TARGET_UPDATES_PER_PARTICLE_REFINE3D_AUTO, &
                    &' UPDATES/PARTICLE (MINIMUM: ', params%minits, ')'
            else
                write(logfhandle,'(A,I0)') &
                    &'>>> '//WORKFLOW_LABEL//' MAXITS COMMAND-LINE OVERRIDE: ', params%maxits
            endif
        end subroutine set_refine3D_auto_sampling

        ! Prepare external input e/o volumes & FSC for refinement and import into project
        subroutine prepare_external_init_vol( init_vol )
            use simple_refine3D_fnames, only: refine3D_startvol_fname, refine3D_fsc_fname
            type(string), intent(inout) :: init_vol
            type(string)      :: init_even, init_odd, new_vol, new_even, new_odd
            type(image)       :: vol, vol_even, vol_odd
            real, allocatable :: fsc(:)
            real    :: ave, stdev, maxv, minv, msk, v
            integer :: ldim(3), n
            init_even = add2fbody(init_vol, MRC_EXT, '_even')
            init_odd  = add2fbody(init_vol, MRC_EXT, '_odd')
            if( .not.(file_exists(init_even) .and. file_exists(init_odd)) )then
                THROW_HARD('Expected even/odd half-volumes for input volume '//init_vol%to_char()//' not found;&
                & looking for: '//init_even%to_char()//' and '//init_odd%to_char())
            endif
            call find_ldim_nptcls(init_vol, ldim, n)
            init_smpd = find_img_smpd(init_vol)
            init_box  = ldim(1)
            if( project_init_vol_compatible() )then
                write(logfhandle,'(A,1X,A)') '>>> '//WORKFLOW_LABEL//' USING E/O INPUT VOLUMES'
                write(logfhandle,'(A,I0,A,F8.4)') '>>> INPUT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
            else
                THROW_HARD('Input e/o volumes must have same dimensions as image')
            endif
            call vol%new(ldim, init_smpd)
            call vol%read(init_vol)
            call vol_even%new(ldim, init_smpd)
            call vol_even%read(init_even)
            if( (abs(vol_even%get_smpd() - init_smpd) > 1.e-6) .or. (vol_even%get_box() /= init_box) )then
                THROW_HARD('Even half-volumes for input volume '//init_vol%to_char()//' have different dimensions/sampling')
            endif
            call vol_odd%new(ldim, init_smpd)
            call vol_odd%read(init_odd)
            if( (abs(vol_odd%get_smpd() - init_smpd) > 1.e-6) .or. (vol_odd%get_box() /= init_box) )then
                THROW_HARD('Odd half-volumes for input volume '//init_vol%to_char()//' have different dimensions/sampling')
            endif
            ! calculate fsc and import into project
            new_vol  = refine3D_startvol_fname(1)
            new_even = add2fbody(new_vol, MRC_EXT, '_even')
            new_odd  = add2fbody(new_vol, MRC_EXT, '_odd')
            ! normalization of volumes
            msk = 0.5 * params%mskdiam / params%smpd
            call vol_even%stats('foreground', ave, stdev, maxv, minv, msk=msk)
            v = stdev * real(init_box) / real(params%box) ! foreground stdev to the data-quotient reference scale
            call vol_even%norm_ext(ave, v)
            call vol_odd%stats('foreground', ave, stdev, maxv, minv, msk=msk)
            v = stdev * real(init_box) / real(params%box) ! foreground stdev to the data-quotient reference scale
            call vol_odd%norm_ext(ave, v)
            call vol%stats('foreground', ave, stdev, maxv, minv, msk=msk)
            v = stdev * real(init_box) / real(params%box) ! foreground stdev to the data-quotient reference scale
            call vol%norm_ext(ave, v)
            call vol_even%write(new_even)
            call vol_odd%write(new_odd)
            call vol%write(new_vol)
            call vol_even%mask3D_soft(msk)
            call vol_odd%mask3D_soft(msk)
            call vol_even%fft
            call vol_odd%fft
            allocate(fsc(vol%get_lfny(1)), source=0.)
            call vol_even%fsc(vol_odd, fsc)
            call arr2file(fsc, refine3D_fsc_fname(1))
            if( all(fsc < 0.9) ) THROW_HARD('Calculated FSC is too low for refinement')
            call spproj%read_segment('out', params%projfile)
            init_vol = new_vol  ! renaming of global filename
            call spproj%add_vol2os_out(init_vol, init_smpd, 1, 'vol')
            call spproj%add_fsc2os_out(refine3D_fsc_fname(1), 1, init_box)
            call spproj%write_segment_inside('out', params%projfile)
            ! cleanup
            deallocate(fsc)
            call spproj%kill
            call vol%kill; call vol_even%kill; call vol_odd%kill
            call init_even%kill; call init_odd%kill
            call new_vol%kill; call new_even%kill; call new_odd%kill
        end subroutine prepare_external_init_vol

        subroutine initialize_external_reference_poses
            type(sp_project) :: pose_init_proj
            integer :: nactive
            call pose_init_proj%read_segment('ptcl3D', params%projfile)
            nactive = pose_init_proj%os_ptcl3D%count_state_gt_zero()
            call pose_init_proj%kill
            if( nactive < 1 ) THROW_HARD(WORKFLOW_LABEL//' pose initialization found no active particles')
            pose_init_refs(1) = init_vol
            call initialize_poses_against_external_references(params, cline, xrefine3D, xrec3D, &
                &nactive, pose_init_refs, pose_init_checkpoint, 1)
            init_vol = pose_init_checkpoint(1)
            call cline%set('startit',    2)
            call cline%set('which_iter', 2)
            call cline%set('extr_iter',  2)
            write(logfhandle,'(A)') &
                &'>>> REFINE3D_AUTO POSE INITIALIZATION CHECKPOINT RECONSTRUCTED; ENTERING EUCLIDEAN REFINEMENT'
        end subroutine initialize_external_reference_poses

    end subroutine exec_refine3D_auto

    subroutine exec_refine3D_states( self, cline )
        use simple_abinitio_utils, only: gen_ortho_reprojs4viz, write_final_rec_outputs
        use simple_commanders_rec, only: commander_rec3D
        use simple_estimate_ssnr,  only: lpstages_setlims
        class(commander_refine3D_states), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(cmdline)             :: cline_rec3D
        type(parameters)          :: params, params_final_rec
        type(sp_project)          :: spproj
        type(lp_crop_inf)         :: lpinfo_multi(2)
        type(string), allocatable :: init_vols(:)
        type(string)              :: multivol_mode, flex_arg, pose_policy_arg
        integer, parameter :: NSAMPLE_PER_STATE_REFINE3D_STATES = 10000
        integer, parameter :: NSAMPLE_REFINE3D_STATES_CAP       = 100000
        integer, parameter :: STAGE1_NSPACE                    = 2500
        integer, parameter :: STAGE2_NSPACE                    = 5000
        integer, parameter :: STAGE2_NSPACE_SUB                = 500
        integer, parameter :: FREQUENCY_BLOCK_NITS             = 3
        integer, parameter :: INIT_MAXITS_REFINE3D_STATES = 5
        integer, parameter :: STAGE2_MINITS              = 5
        integer, parameter :: MINITS_REFINE3D_STATES      = 10
        integer, parameter :: MAXITS_REFINE3D_STATES_CAP  = 50
        real,    parameter :: TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES = 4.0
        real,    parameter :: STATE_OVERLAP_EARLY_REFINE3D_STATES         = 0.95
        real,    parameter :: STATE_OVERLAP_NEIGH_REFINE3D_STATES         = 0.99
        real,    parameter :: LPSTART_REFINE3D_STATES = 10.0
        real,    parameter :: LPSTOP_REFINE3D_STATES  = 6.0
        character(len=*), parameter :: WORKFLOW_LABEL = 'REFINE3D_STATES'
        integer :: nstates_project, nptcls_eff, nsample_target, nptcls_per_iter, local_nspace_sub
        integer :: maxits_user, stage_cap, init_niters, stage2_niters, total_iter
        integer :: maxits_glob_multi, init_stage_cap, min_maxits_required
        real    :: update_frac_auto, state_overlap, local_ang_bound, local_inpl_bound, local_shift_bound
        logical :: l_maxits_defined, l_init_state_assignment, l_nstates_on_cline, l_flex_requested, l_nsample_auto
        logical :: l_has_project_multistates, l_run_init_stage, l_run_prob_neigh_stage
        ! commanders
        type(commander_rec3D)    :: xrec3D
        type(commander_refine3D) :: xrefine3D
        maxits_user    = 0
        init_niters    = 0
        stage2_niters  = 0
        init_stage_cap = INIT_MAXITS_REFINE3D_STATES
        l_init_state_assignment   = .false.
        l_has_project_multistates = .false.
        l_run_init_stage          = .false.
        l_run_prob_neigh_stage    = .false.
        l_flex_requested          = .false.
        min_maxits_required       = STAGE2_MINITS
        local_nspace_sub          = STAGE2_NSPACE_SUB
        local_ang_bound           = -1.
        local_inpl_bound          = -1.
        local_shift_bound         = -1.
        call cline%set('prg', 'refine3D_states')
        l_nstates_on_cline = cline%defined('nstates')
        if( cline%defined('multivol_mode') ) THROW_HARD(WORKFLOW_LABEL//' uses pose_policy, not multivol_mode')
        if( cline%defined('prob_neigh_mode') ) THROW_HARD(WORKFLOW_LABEL//' derives prob_neigh_mode from pose_policy')
        if( .not. cline%defined('pose_policy') ) call cline%set('pose_policy', 'global')
        pose_policy_arg = cline%get_carg('pose_policy')
        select case(trim(pose_policy_arg%to_char()))
            case('fixed')
                call cline%set('multivol_mode', 'input_oris_fixed')
                call cline%set('prob_neigh_mode', 'geom')
            case('local')
                call cline%set('multivol_mode', 'input_oris_refine')
                call cline%set('prob_neigh_mode', 'geom')
            case('global')
                call cline%set('multivol_mode', 'input_oris_refine')
                call cline%set('prob_neigh_mode', 'state')
            case default
                THROW_HARD(WORKFLOW_LABEL//' supports pose_policy=fixed|local|global')
        end select
        if( cline%defined('local_ang_bound') ) local_ang_bound = cline%get_rarg('local_ang_bound')
        if( cline%defined('local_inpl_bound') ) local_inpl_bound = cline%get_rarg('local_inpl_bound')
        if( cline%defined('local_shift_bound') ) local_shift_bound = cline%get_rarg('local_shift_bound')
        if( trim(pose_policy_arg%to_char()).ne.'local' .and. &
            &any([local_ang_bound, local_inpl_bound, local_shift_bound] >= 0.) )then
            THROW_HARD('local pose-bound overrides require pose_policy=local')
        endif
        if( local_ang_bound >= 0. )then
            if( local_ang_bound > 180. ) THROW_HARD('local_ang_bound must be between 0 and 180 degrees')
            if( local_ang_bound < 0.01 )then
                local_nspace_sub = STAGE2_NSPACE
            else
                local_nspace_sub = min(STAGE2_NSPACE, max(2, round2even(2. / &
                    &max(1.e-6, 1. - cos(deg2rad(local_ang_bound))))))
            endif
        endif
        if( local_inpl_bound >= 0. ) call cline%set('prob_athres', local_inpl_bound)
        if( local_shift_bound >= 0. ) call cline%set('trs', local_shift_bound)
        multivol_mode = cline%get_carg('multivol_mode')
        if( .not. cline%defined('filt_mode')     ) call cline%set('filt_mode', 'nonuniform_lpset')
        if( cline%defined('flex') )then
            ! Multi-state initialization with flex_pca
            flex_arg = cline%get_carg('flex')
            l_flex_requested = trim(flex_arg%to_char()).eq.'yes'
            call flex_arg%kill
        endif
        if( l_flex_requested )then
            if( .not. l_nstates_on_cline ) THROW_HARD(WORKFLOW_LABEL//' flex=yes requires nstates >= 3')
            nstates_project = cline%get_iarg('nstates')
            if( nstates_project < 3 ) THROW_HARD(WORKFLOW_LABEL//' flex=yes requires nstates >= 3')
        else
            call set_refine3D_states_nstates()
        endif
        ! hard defaults
        call cline%set('balance',        'yes')
        call cline%set('greedy_sampling', 'no')
        call cline%set('frac_best',       1.0)
        if( trim(multivol_mode%to_char()).eq.'input_oris_fixed' )then
            call cline%set('trail_rec', 'no')
        else
            call cline%set('trail_rec', 'yes')
        endif
        call cline%set('objfun',      'euclid')
        call cline%set('lplim_crit',       0.5)
        call cline%set('incrreslim',      'no')
        call cline%set('nu_refine',       'no')
        ! overridable defaults
        if( .not. cline%defined('combine_eo')      ) call cline%set('combine_eo',        'no')
        if( .not. cline%defined('envfsc')          ) call cline%set('envfsc',            'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('center')          ) call cline%set('center',            'no')
        if( .not. cline%defined('sigma_est')       ) call cline%set('sigma_est',     'global')
        if( .not. cline%defined('prob_inpl')       ) call cline%set('prob_inpl',        'yes')
        l_nsample_auto = .not. cline%defined('nsample')
        if( .not. l_nsample_auto ) l_nsample_auto = cline%get_iarg('nsample') <= 0
        if( l_nsample_auto )then
            call cline%set('nsample', default_refine3D_states_nsample())
        endif
        ! staged frequency marching always uses planner-driven downscaling
        if( cline%defined('autoscale') ) write(logfhandle,'(A)') &
            &'>>> '//WORKFLOW_LABEL//' IGNORES autoscale: frequency marching owns downscaling'
        call cline%set('autoscale', 'yes')
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',               'yes')
        if( .not. cline%defined('lpstart')     ) call cline%set('lpstart', LPSTART_REFINE3D_STATES)
        if( .not. cline%defined('lpstop')      ) call cline%set('lpstop',  LPSTOP_REFINE3D_STATES)
        if( .not. cline%defined('automsk')     ) call cline%set('automsk',               'no')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap', STATE_OVERLAP_NEIGH_REFINE3D_STATES)
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol',               'no')
        l_maxits_defined = cline%defined('maxits')
        if( l_maxits_defined )then
            if( trim(multivol_mode%to_char()).eq.'input_oris_fixed' ) min_maxits_required = 1
            maxits_user = cline%get_iarg('maxits')
            if( maxits_user < min_maxits_required )then
                THROW_HARD('maxits must be >= '//int2str(min_maxits_required)//' for '//WORKFLOW_LABEL)
            endif
        endif
        call params%new(cline)
        write(logfhandle,'(A,A)') '>>> '//WORKFLOW_LABEL//' POSE_POLICY: ', trim(params%pose_policy)
        call validate_refine3D_states_mode()
        call validate_refine3D_states_combine_eo()
        call validate_refine3D_states_filtering()
        call validate_refine3D_states_prob_neigh_mode()
        call cline%set('mkdir', 'no')
        if( l_flex_requested )then
            call run_flex_pca()
            call set_refine3D_states_nstates()
            params%nstates = nstates_project
            if( l_nsample_auto )then
                params%nsample = default_refine3D_states_nsample()
                call cline%set('nsample', params%nsample)
            endif
        endif
        call set_refine3D_states_sampling()
        call prepare_refine3D_states_class_sampling()
        call configure_refine3D_states_stages()
        call set_refine3D_states_downscaling()
        if( .not. l_flex_requested ) call initialize_state_volumes()
        call cline%set('prg', 'refine3D')
        maxits_glob_multi = 0
        if( l_run_init_stage       ) maxits_glob_multi = maxits_glob_multi + init_stage_cap
        if( l_run_prob_neigh_stage ) maxits_glob_multi = maxits_glob_multi + stage_cap
        total_iter = 0
        if( cline%defined('startit') ) total_iter = max(0, cline%get_iarg('startit') - 1)
        call cline%set('maxits_glob', max(1, total_iter + maxits_glob_multi))
        if( l_run_init_stage )then
            if( trim(params%pose_policy).eq.'fixed' )then
                call run_refine3D_states_frequency_march('prob_state', STAGE1_NSPACE, 0, 1, init_niters, &
                    &init_stage_cap, STATE_OVERLAP_EARLY_REFINE3D_STATES)
            else
                call cline%set('lp', params%lpstart)
                call cline%set('lpstop', params%lpstart)
                call run_refine3D_states_stage(0, 'prob_state', STAGE1_NSPACE, 0, 1, init_niters, &
                    &init_stage_cap, STATE_OVERLAP_EARLY_REFINE3D_STATES)
            endif
        endif
        if( l_run_prob_neigh_stage )then
            call run_refine3D_states_frequency_march('prob_neigh', STAGE2_NSPACE, local_nspace_sub, &
                &STAGE2_MINITS, stage2_niters, stage_cap, params%overlap)
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> '//WORKFLOW_LABEL//' STAGE ITERATIONS INIT/PROB_NEIGH/TOTAL: ', &
            &init_niters, '/', stage2_niters, '/', total_iter
        call ensure_all_active_particles_updated()
        ! re-reconstruct from all particle images
        cline_rec3D = cline
        call cline_rec3D%set('prg',            'reconstruct3D')
        call cline_rec3D%set('outfile', 'RESOLUTION_FINAL.txt')
        call cline_rec3D%set('postprocess',              'yes')
        call cline_rec3D%delete('trail_rec')
        call cline_rec3D%delete('refine')
        call cline_rec3D%delete('objfun_den')
        call cline_rec3D%delete('objfun_den_w')
        call cline_rec3D%delete('sigma_est')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%delete('ufrac_trec')
        call cline_rec3D%delete('endit')
        call cline_rec3D%delete('box_crop')
        call cline_rec3D%delete('smpd_crop')
        call cline_rec3D%set('objfun', 'cc')
        if( params%l_nonuniform )then
            call cline_rec3D%set('filt_mode', 'none')
            call cline_rec3D%set('automsk', 'no')
        endif
        call cline_rec3D%set('nu_refine', 'no')
        call xrec3D%execute(cline_rec3D)
        call params_final_rec%new(cline_rec3D)
        params_final_rec%box  = params_final_rec%box_crop
        params_final_rec%smpd = params_final_rec%smpd_crop
        call spproj%read_segment('out', params_final_rec%projfile)
        call write_final_rec_outputs(params_final_rec, spproj, params_final_rec%lpstop)
        call gen_ortho_reprojs4viz(params_final_rec, spproj)
        call spproj%kill
        call cleanup_init_vols()
        call pose_policy_arg%kill
        call multivol_mode%kill
        call simple_end('**** SIMPLE_REFINE3D_STATES NORMAL STOP ****')

    contains

        integer function default_refine3D_states_nsample() result(nsample_default)
            nsample_default = min(NSAMPLE_REFINE3D_STATES_CAP, &
                &NSAMPLE_PER_STATE_REFINE3D_STATES * nstates_project)
        end function default_refine3D_states_nsample

        subroutine validate_refine3D_states_mode()
            select case(trim(params%multivol_mode))
                case('input_oris_refine', 'input_oris_fixed')
                    ! supported
                case default
                    THROW_HARD('Unsupported multivol_mode for '//WORKFLOW_LABEL//': '//trim(params%multivol_mode))
            end select
            write(logfhandle,'(A,A)') '>>> '//WORKFLOW_LABEL//' MULTIVOL_MODE: ', trim(params%multivol_mode)
        end subroutine validate_refine3D_states_mode

        subroutine validate_refine3D_states_combine_eo()
            if( trim(params%combine_eo).ne.'no' )then
                THROW_HARD(WORKFLOW_LABEL//' does not support combine_eo; it is an LP-set multi-state workflow')
            endif
        end subroutine validate_refine3D_states_combine_eo

        subroutine validate_refine3D_states_filtering()
            select case(trim(params%filt_mode))
                case('fsc', 'nonuniform_lpset', 'none')
                    ! supported
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports filt_mode=fsc|nonuniform_lpset|none')
            end select
        end subroutine validate_refine3D_states_filtering

        subroutine validate_refine3D_states_prob_neigh_mode()
            select case(trim(params%prob_neigh_mode))
                case('geom', 'state')
                    ! supported
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports prob_neigh_mode=geom|state')
            end select
        end subroutine validate_refine3D_states_prob_neigh_mode

        subroutine configure_refine3D_states_stages()
            select case(trim(params%multivol_mode))
                case('input_oris_fixed')
                    l_run_init_stage       = .true.
                    l_run_prob_neigh_stage = .false.
                    init_stage_cap         = stage_cap
                case('input_oris_refine')
                    l_run_init_stage       = l_init_state_assignment
                    l_run_prob_neigh_stage = .true.
            end select
            write(logfhandle,'(A,L1,A,L1)') '>>> '//WORKFLOW_LABEL//' STAGES INIT/PROB_NEIGH: ', &
                &l_run_init_stage, '/', l_run_prob_neigh_stage
        end subroutine configure_refine3D_states_stages

        subroutine set_refine3D_states_nstates()
            type(sp_project) :: state_proj
            type(string)     :: projfile
            integer, allocatable :: pops(:)
            integer :: state, nactive_labels, nstates_cline, nstates_labels
            logical :: l_has_state_labels
            if( .not. cline%defined('projfile') )then
                THROW_HARD('projfile is required for '//WORKFLOW_LABEL)
            endif
            projfile = cline%get_carg('projfile')
            call state_proj%read_segment('ptcl3D', projfile)
            nactive_labels    = 0
            nstates_labels    = 1
            l_has_state_labels = state_proj%os_ptcl3D%isthere('state')
            if( l_has_state_labels )then
                nactive_labels = state_proj%os_ptcl3D%count_state_gt_zero()
                if( nactive_labels > 0 ) nstates_labels = state_proj%os_ptcl3D%get_n('state')
            endif
            l_has_project_multistates = nactive_labels > 0 .and. nstates_labels > 1
            if( l_has_project_multistates )then
                select case(trim(multivol_mode%to_char()))
                    case('input_oris_fixed')
                        if( .not. l_flex_requested ) THROW_HARD(WORKFLOW_LABEL//' input_oris_fixed expects state=0/1 input')
                end select
                nstates_project = nstates_labels
                if( l_nstates_on_cline )then
                    nstates_cline = cline%get_iarg('nstates')
                    if( nstates_cline /= nstates_project )then
                        THROW_HARD('command-line nstates does not match project state labels for '//WORKFLOW_LABEL)
                    endif
                endif
                call state_proj%os_ptcl3D%get_pops(pops, 'state', maxn=nstates_project)
                do state = 1,nstates_project
                    if( pops(state) < 1 )then
                        write(logfhandle,*) 'state, population: ', state, pops(state)
                        THROW_HARD(WORKFLOW_LABEL//' requires every state label to have at least one active particle')
                    endif
                enddo
                write(logfhandle,'(A,I0)') '>>> '//WORKFLOW_LABEL//' NSTATES FROM PROJECT: ', nstates_project
            else
                if( .not. l_nstates_on_cline )then
                    THROW_HARD(WORKFLOW_LABEL//' requires nstates > 1 for state=0/1 input')
                endif
                nstates_cline = cline%get_iarg('nstates')
                if( nstates_cline <= 1 )then
                    THROW_HARD('nstates must be > 1 for '//WORKFLOW_LABEL//' initial state assignment mode')
                endif
                nstates_project = nstates_cline
                l_init_state_assignment = .true.
                write(logfhandle,'(A,I0)') &
                    &'>>> '//WORKFLOW_LABEL//' NO PROJECT MULTI-STATE ASSIGNMENTS; INITIALIZING NSTATES: ', nstates_project
            endif
            call cline%set('nstates', nstates_project)
            if( allocated(pops) ) deallocate(pops)
            call state_proj%kill
            call projfile%kill
        end subroutine set_refine3D_states_nstates

        subroutine set_refine3D_states_sampling()
            type(sp_project) :: sampling_proj
            integer :: maxits_auto
            nsample_target = params%nsample
            if( nsample_target < 1 ) THROW_HARD('nsample must be >= 1 for '//WORKFLOW_LABEL)
            if( l_init_state_assignment )then
                call sampling_proj%read_segment('ptcl3D', params%projfile)
                nptcls_eff = sampling_proj%os_ptcl3D%get_noris(consider_state=.true.)
                if( nptcls_eff < 1 ) nptcls_eff = sampling_proj%os_ptcl3D%get_noris()
            else
                call sampling_proj%read(params%projfile)
                nptcls_eff = sampling_proj%count_state_gt_zero()
            endif
            call sampling_proj%kill
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//WORKFLOW_LABEL)
            nptcls_per_iter = min(nptcls_eff, nsample_target)
            if( trim(params%multivol_mode).eq.'input_oris_fixed' )then
                nptcls_per_iter       = nptcls_eff
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A)') &
                    &'>>> '//WORKFLOW_LABEL//' INPUT_ORIS_FIXED ACTIVE PARTICLES: ', &
                    &nptcls_eff, ' -> FULL UPDATE'
            else if( nptcls_eff <= nsample_target )then
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)') &
                    &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    &nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
            else
                update_frac_auto = real(nsample_target) / real(nptcls_eff)
                if( update_frac_auto <= 0.99 )then
                    params%update_frac   = update_frac_auto
                    params%l_update_frac = .true.
                    params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                    call cline%set('update_frac', update_frac_auto)
                    write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                        &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                        &nptcls_eff, '/', nsample_target, '/', update_frac_auto
                else
                    params%update_frac   = 1.0
                    params%l_update_frac = .false.
                    params%l_trail_rec   = .false.
                    call cline%delete('update_frac')
                    write(logfhandle,'(A,I0,A,I0,A)') &
                        &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                        &nptcls_eff, '/', nsample_target, ' -> FULL UPDATE'
                endif
            endif
            if( l_maxits_defined )then
                stage_cap = maxits_user
                params%maxits = stage_cap
                write(logfhandle,'(A,I0)') &
                    &'>>> '//WORKFLOW_LABEL//' STAGE MAXITS COMMAND-LINE OVERRIDE: ', stage_cap
            else
                maxits_auto = ceiling((TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES * real(nptcls_eff)) / real(nptcls_per_iter))
                stage_cap   = max(MINITS_REFINE3D_STATES, min(MAXITS_REFINE3D_STATES_CAP, max(STAGE2_MINITS, maxits_auto)))
                params%maxits = stage_cap
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0,A,F5.1,A)') '>>> '//WORKFLOW_LABEL//' STAGE MAXITS: ', &
                    &stage_cap, ' FOR ~', TARGET_UPDATES_PER_PARTICLE_REFINE3D_STATES, ' UPDATES/PARTICLE'
            endif
        end subroutine set_refine3D_states_sampling

        subroutine prepare_refine3D_states_class_sampling()
            type(builder)    :: sampling_build
            type(cmdline)    :: cline_sampling
            type(parameters) :: params_sampling
            type(class_sample), allocatable :: clssmp(:)
            integer :: nproj_bins
            if( trim(params%balance) /= 'yes' ) return
            if( .not. params%l_update_frac ) return
            cline_sampling = cline
            call cline_sampling%set('prg', 'refine3D')
            call cline_sampling%set('refine', 'prob')
            call cline_sampling%set('nspace', STAGE2_NSPACE)
            call sampling_build%init_params_and_build_general_tbox(cline_sampling, params_sampling, do3d=.true.)
            if( sampling_build%spproj%is_virgin_field(params_sampling%oritype) )then
                call sampling_build%kill_general_tbox
                call cline_sampling%kill
                THROW_HARD(WORKFLOW_LABEL//' projection-balanced fractional updates require prior 3D orientations')
            endif
            call make_projdir_class_samples(sampling_build, clssmp, nproj_bins)
            call write_class_samples(clssmp, string(CLASS_SAMPLING_FILE))
            write(logfhandle,'(A,I0,A,F5.2)') &
                &'>>> '//WORKFLOW_LABEL//' PROJDIR-BALANCED SAMPLING BINS/FRAC_BEST: ', &
                &nproj_bins, '/', params%frac_best
            call deallocate_class_samples(clssmp)
            call sampling_build%kill_general_tbox
            call cline_sampling%kill
        end subroutine prepare_refine3D_states_class_sampling

        subroutine make_projdir_class_samples( sampling_build, clssmp, nproj_bins )
            type(builder),                    intent(inout) :: sampling_build
            type(class_sample), allocatable,  intent(inout) :: clssmp(:)
            integer,                          intent(out)   :: nproj_bins
            integer, allocatable :: allinds(:), states(:), projs(:), pinds(:), proj_pops(:)
            real,    allocatable :: corrs(:)
            integer :: iproj, ibin, nparticles, nprojs
            if( allocated(clssmp) ) call deallocate_class_samples(clssmp)
            call sampling_build%spproj_field%set_projs(sampling_build%eulspace)
            nprojs     = sampling_build%eulspace%get_noris()
            states     = sampling_build%spproj_field%get_all_asint('state')
            projs      = sampling_build%spproj_field%get_all_asint('proj')
            corrs      = sampling_build%spproj_field%get_all('corr')
            nparticles = size(states)
            allinds    = (/(iproj, iproj=1,nparticles)/)
            allocate(proj_pops(nprojs), source=0)
            do iproj = 1,nprojs
                proj_pops(iproj) = count(states > 0 .and. projs == iproj)
            enddo
            nproj_bins = count(proj_pops > 0)
            if( nproj_bins < 1 )then
                if( allocated(allinds)   ) deallocate(allinds)
                if( allocated(states)    ) deallocate(states)
                if( allocated(projs)     ) deallocate(projs)
                if( allocated(corrs)     ) deallocate(corrs)
                if( allocated(proj_pops) ) deallocate(proj_pops)
                THROW_HARD(WORKFLOW_LABEL//' projection-balanced sampling found no active projection bins')
            endif
            allocate(clssmp(nproj_bins))
            ibin = 0
            do iproj = 1,nprojs
                if( proj_pops(iproj) < 1 ) cycle
                ibin = ibin + 1
                pinds = pack(allinds, mask=states > 0 .and. projs == iproj)
                clssmp(ibin)%clsind = iproj
                clssmp(ibin)%pop    = size(pinds)
                allocate(clssmp(ibin)%pinds(clssmp(ibin)%pop), source=pinds)
                allocate(clssmp(ibin)%ccs(clssmp(ibin)%pop), source=pack(corrs, mask=states > 0 .and. projs == iproj))
                call hpsort(clssmp(ibin)%ccs, clssmp(ibin)%pinds)
                call reverse(clssmp(ibin)%ccs)
                call reverse(clssmp(ibin)%pinds)
                if( allocated(pinds) ) deallocate(pinds)
            enddo
            if( allocated(allinds)   ) deallocate(allinds)
            if( allocated(states)    ) deallocate(states)
            if( allocated(projs)     ) deallocate(projs)
            if( allocated(corrs)     ) deallocate(corrs)
            if( allocated(proj_pops) ) deallocate(proj_pops)
        end subroutine make_projdir_class_samples

        subroutine set_refine3D_states_downscaling()
            lpinfo_multi(2)%trslim      = min(8., max(2.0, AHELIX_WIDTH / params%smpd))
            lpinfo_multi(2)%box_crop    = params%box
            lpinfo_multi(2)%smpd_crop   = params%smpd
            lpinfo_multi(2)%l_autoscale = .false.
            call lpstages_setlims(params%box, 2, params%smpd, params%lpstart, params%lpstop, lpinfo_multi)
            if( local_shift_bound >= 0. .and. trim(params%pose_policy).eq.'local' )then
                params%trs = local_shift_bound
                call cline%set('trs', params%trs)
            else
                params%trs = lpinfo_multi(2)%trslim
                call cline%set('trs', params%trs)
            endif
            if( lpinfo_multi(2)%l_autoscale )then
                params%box_crop  = lpinfo_multi(2)%box_crop
                params%smpd_crop = lpinfo_multi(2)%smpd_crop
                call cline%set('box_crop', params%box_crop)
                call cline%set('smpd_crop', params%smpd_crop)
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> '//WORKFLOW_LABEL//' DOWNSCALING BOX/SMPD_CROP: ', &
                    &params%box, '/', params%box_crop, '/', params%smpd_crop
            else
                call cline%delete('box_crop')
                call cline%delete('smpd_crop')
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' DOWNSCALING: native sampling retained'
            endif
        end subroutine set_refine3D_states_downscaling

        subroutine initialize_state_volumes()
            integer :: state
            allocate(init_vols(nstates_project))
            if( complete_input_volumes_defined() )then
                call validate_input_volumes()
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' USING INPUT STATE VOLUMES'
                return
            else if( any_input_volumes_defined() )then
                THROW_HARD(WORKFLOW_LABEL//' requires either all vol1..volN inputs or none')
            endif
            if( project_state_volumes_compatible() )then
                do state = 1,nstates_project
                    call cline%set('vol'//int2str(state), init_vols(state))
                enddo
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' USING PROJECT STATE VOLUMES'
                return
            endif
            if( l_init_state_assignment )then
                if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
                    write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' STARTUP STATE VOLUMES DELEGATED TO BASE REFINE3D'
                    return
                endif
                THROW_HARD(WORKFLOW_LABEL//' state=0/1 input without vol1..volN requires distributed startup; set nparts')
            endif
            call prepare_startup_reconstruct3D_cline()
            call xrec3D%execute(cline_rec3D)
            do state = 1,nstates_project
                call cline%set('vol'//int2str(state), refine3D_state_vol_fname(state))
            enddo
            write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' INITIALIZED STATE VOLUMES BY RECONSTRUCTION'
        end subroutine initialize_state_volumes

        subroutine run_flex_pca()
            use simple_commanders_flex_pca, only: commander_flex_pca
            type(commander_flex_pca) :: xflex
            type(cmdline)            :: cline_flex
            type(sp_project)         :: flex_proj
            type(string)             :: vol1
            integer :: state, flex_box, nstates_requested, nstates_flex
            integer :: nactive_labels, nstates_labels
            real    :: flex_smpd
            ! validate inputs
            nstates_requested = params%nstates
            if( nstates_requested < 3 ) THROW_HARD(WORKFLOW_LABEL//' flex=yes requires nstates >= 3')
            do state = 1,nstates_requested
                if( cline%defined('vol'//int2str(state)) )then
                    THROW_HARD(WORKFLOW_LABEL//' flex=yes does not support vol1..volN inputs')
                endif
            enddo
            ! validate project
            call flex_proj%read_segment('ptcl3D', params%projfile)
            nactive_labels = 0
            nstates_labels = 1
            if( flex_proj%os_ptcl3D%isthere('state') )then
                nactive_labels = flex_proj%os_ptcl3D%count_state_gt_zero()
                if( nactive_labels > 0 ) nstates_labels = flex_proj%os_ptcl3D%get_n('state')
            endif
            call flex_proj%kill
            if( nactive_labels > 0 .and. nstates_labels > 1 )then
                THROW_HARD(WORKFLOW_LABEL//' flex=yes requires an input project with a single state')
            endif
            call flex_proj%read_segment('out', params%projfile)
            if( .not. flex_proj%isthere_in_osout('vol', 1) )then
                THROW_HARD(WORKFLOW_LABEL//' flex=yes requires a consensus project volume')
            endif
            call flex_proj%get_vol('vol', 1, vol1, flex_smpd, flex_box)
            call flex_proj%kill
            if( .not. file_exists(vol1) )then
                THROW_HARD(WORKFLOW_LABEL//' flex=yes consensus project volume does not exist')
            endif
            if( flex_box /= params%box .or. flex_smpd <= 0. .or. abs(flex_smpd - params%smpd) > 1.e-6 )then
                THROW_HARD(WORKFLOW_LABEL//' flex=yes consensus volume must match the project particle sampling')
            endif
            vol1 = simple_abspath(vol1)
            ! execution prg=flex_pca
            cline_flex = cline
            call cline_flex%set('prg',        'flex_pca')
            call cline_flex%set('mkdir',      'no')
            call cline_flex%set('npreimages', nstates_requested)
            call cline_flex%set('vol1',       vol1)
            call cline_flex%delete('nstates')
            call xflex%execute(cline_flex)
            ! output parsing
            call flex_proj%read_segment('ptcl3D', params%projfile)
            if( .not. flex_proj%os_ptcl3D%isthere('state') )then
                THROW_HARD(WORKFLOW_LABEL//' flex_pca did not write particle state labels')
            endif
            nstates_flex = flex_proj%os_ptcl3D%get_n('state')
            call flex_proj%kill
            if( nstates_flex <= 1 ) THROW_HARD(WORKFLOW_LABEL//' flex_pca produced fewer than two states')
            call cline%set('nstates', nstates_flex)
            allocate(init_vols(nstates_flex))
            call flex_proj%read_segment('out', params%projfile)
            do state = 1,nstates_flex
                if( .not. flex_proj%isthere_in_osout('vol_flex', state) )then
                    THROW_HARD(WORKFLOW_LABEL//' flex_pca did not write every state volume')
                endif
                call flex_proj%get_vol('vol_flex', state, init_vols(state), flex_smpd, flex_box)
                if( .not. file_exists(init_vols(state)) )then
                    THROW_HARD(WORKFLOW_LABEL//' flex_pca state volume does not exist')
                endif
                if( flex_box /= params%box .or. flex_smpd <= 0. .or. abs(flex_smpd - params%smpd) > 1.e-6 )then
                    THROW_HARD(WORKFLOW_LABEL//' flex_pca state volumes must match the project particle sampling')
                endif
                ! update command line with flex volumes
                init_vols(state) = simple_abspath(init_vols(state))
                call cline%set('vol'//int2str(state), init_vols(state))
            enddo
            call cline%delete('neigs')
            call flex_proj%kill
            write(logfhandle,'(A,I0)') '>>> '//WORKFLOW_LABEL//' FLEX_PCA INITIALIZED NSTATES: ', nstates_flex
            call cline_flex%kill
            call vol1%kill
        end subroutine run_flex_pca

        subroutine prepare_startup_reconstruct3D_cline()
            cline_rec3D = cline
            call cline_rec3D%set('prg', 'reconstruct3D')
            call cline_rec3D%delete('trail_rec')
            call cline_rec3D%delete('refine')
            call cline_rec3D%delete('objfun_den')
            call cline_rec3D%delete('objfun_den_w')
            call cline_rec3D%delete('sigma_est')
            call cline_rec3D%delete('update_frac')
            call cline_rec3D%delete('ufrac_trec')
            call cline_rec3D%delete('endit')
            call cline_rec3D%delete('box_crop')
            call cline_rec3D%delete('smpd_crop')
            call cline_rec3D%set('objfun', 'cc')
            call cline_rec3D%set('postprocess', 'no')
            call cline_rec3D%set('nu_refine', 'no')
        end subroutine prepare_startup_reconstruct3D_cline

        logical function any_input_volumes_defined() result(l_any)
            integer :: state
            l_any = .false.
            do state = 1,nstates_project
                if( cline%defined('vol'//int2str(state)) )then
                    l_any = .true.
                    return
                endif
            enddo
        end function any_input_volumes_defined

        logical function complete_input_volumes_defined() result(l_complete)
            integer :: state
            l_complete = .true.
            do state = 1,nstates_project
                if( .not. cline%defined('vol'//int2str(state)) )then
                    l_complete = .false.
                    return
                endif
            enddo
        end function complete_input_volumes_defined

        subroutine validate_input_volumes()
            type(string) :: vol
            integer :: state, ldim(3), nptcls_dummy
            real    :: vol_smpd, extent_native, extent_vol
            do state = 1,nstates_project
                vol = cline%get_carg('vol'//int2str(state))
                if( .not. file_exists(vol) ) THROW_HARD('Input volume does not exist: '//vol%to_char())
                call find_ldim_nptcls(vol, ldim, nptcls_dummy)
                vol_smpd = find_img_smpd(vol)
                if( ldim(1) /= ldim(2) .or. ldim(1) /= ldim(3) )then
                    THROW_HARD('Input state volumes must be cubic')
                endif
                if( ldim(1) > params%box )then
                    THROW_HARD('Input state volumes cannot exceed the native project box')
                endif
                ! Any downscaled sampling of the native grid is acceptable — base
                ! refine3D rescales references to the stage crop, and the abinitio3D
                ! split checkpoint hands over volumes at the abinitio ladder crop,
                ! which need not match this workflow's own frequency-stage crop.
                extent_native = real(params%box) * params%smpd
                extent_vol    = real(ldim(1))   * vol_smpd
                if( abs(extent_vol - extent_native) > 0.01 * extent_native )then
                    THROW_HARD('Input state volumes must cover the native physical extent (box*smpd)')
                endif
            enddo
            call vol%kill
        end subroutine validate_input_volumes

        logical function project_state_volumes_compatible() result(l_compatible)
            real    :: init_smpd
            integer :: state, init_box
            l_compatible = .false.
            call spproj%read_segment('out', params%projfile)
            do state = 1,nstates_project
                if( .not. spproj%isthere_in_osout('vol', state) )then
                    call spproj%kill
                    return
                endif
                call spproj%get_vol('vol', state, init_vols(state), init_smpd, init_box)
                if( .not. file_exists(init_vols(state)) )then
                    call spproj%kill
                    return
                endif
                if( init_smpd <= 0. )then
                    call spproj%kill
                    return
                endif
                ! any downscaled sampling of the native grid is acceptable (see
                ! validate_input_volumes); base refine3D rescales to the stage crop
                if( init_box > params%box .or. &
                    &abs(real(init_box)*init_smpd - real(params%box)*params%smpd) > &
                    &0.01 * real(params%box)*params%smpd )then
                    call spproj%kill
                    return
                endif
            enddo
            l_compatible = .true.
            call spproj%kill
        end function project_state_volumes_compatible

        subroutine run_refine3D_states_frequency_march( refine_mode, nspace_stage, nspace_sub_stage, &
            &min_total_iters, niters, max_total_iters, overlap_target )
            character(len=*), intent(in)  :: refine_mode
            integer,          intent(in)  :: nspace_stage, nspace_sub_stage, min_total_iters, max_total_iters
            integer,          intent(out) :: niters
            real,             intent(in)  :: overlap_target
            type(refine3D_stage_plan_entry), allocatable :: stages(:)
            integer :: block, block_niters, block_minits
            niters = 0
            if( max_total_iters < min_total_iters )then
                THROW_HARD(WORKFLOW_LABEL//' frequency-stage cap is below its required minimum')
            endif
            call plan_refine3D_frequency_stages(params%box, params%smpd, params%lpstart, params%lpstop, &
                &max_total_iters, FREQUENCY_BLOCK_NITS, total_iter + 1, lpinfo_multi(2), stages)
            write(logfhandle,'(A,I0,A,F7.2,A,F7.2)') &
                &'>>> '//WORKFLOW_LABEL//' FREQUENCY MARCHING NSTAGES/LPSTART/LPSTOP: ', &
                &size(stages), '/', params%lpstart, '/', params%lpstop
            do block = 1,size(stages)
                call cline%set('trs',    stages(block)%lpinfo%trslim)
                if( local_shift_bound >= 0. .and. trim(params%pose_policy).eq.'local' )then
                    call cline%set('trs', local_shift_bound)
                endif
                call cline%set('lp',     stages(block)%lpinfo%lp)
                call cline%set('lpstop', stages(block)%lpinfo%lp)
                if( stages(block)%lpinfo%l_autoscale )then
                    call cline%set('box_crop',  stages(block)%lpinfo%box_crop)
                    call cline%set('smpd_crop', stages(block)%lpinfo%smpd_crop)
                else
                    call cline%delete('box_crop')
                    call cline%delete('smpd_crop')
                endif
                ! Each block may stop early on the state-overlap target once the
                ! march-wide iteration minimum is satisfied; the march still visits
                ! every frequency block so the low-pass schedule completes.
                block_minits = max(1, min(stages(block)%nits, min_total_iters - niters))
                call run_refine3D_states_stage(block, refine_mode, nspace_stage, nspace_sub_stage, &
                    &block_minits, block_niters, stages(block)%nits, overlap_target)
                niters = niters + block_niters
            enddo
            deallocate(stages)
        end subroutine run_refine3D_states_frequency_march

        subroutine run_refine3D_states_stage( stage, refine_mode, nspace_stage, nspace_sub_stage, min_stage_iters, &
            &niters, max_stage_iters, stage_overlap_target )
            integer,          intent(in)  :: stage, nspace_stage, nspace_sub_stage, min_stage_iters
            character(len=*), intent(in)  :: refine_mode
            integer,          intent(out) :: niters
            integer,          intent(in)  :: max_stage_iters
            real,             intent(in)  :: stage_overlap_target
            integer :: stage_start, stage_limit
            niters        = 0
            state_overlap = 0.
            stage_limit   = max_stage_iters
            stage_start   = total_iter + 1
            write(logfhandle,'(A,I0,A,A,A,I0,A,I0,A,F7.4)') '>>> '//WORKFLOW_LABEL//' ENTERING STAGE ', stage, &
                &' REFINE=', trim(refine_mode), ' NSPACE/NSPACE_SUB=', nspace_stage, '/', nspace_sub_stage, &
                &' OVERLAP_TARGET=', stage_overlap_target
            call cline%set('refine', refine_mode)
            call cline%set('nspace', nspace_stage)
            if( nspace_sub_stage > 0 )then
                call cline%set('nspace_sub', nspace_sub_stage)
            else
                call cline%delete('nspace_sub')
            endif
            call cline%set('maxits', stage_limit)
            call cline%set('minits', min_stage_iters)
            call cline%set('startit', stage_start)
            call cline%set('which_iter', stage_start)
            call cline%set('extr_iter', stage_start)
            call cline%set('overlap', stage_overlap_target)
            call cline%delete('endit')
            call xrefine3D%execute(cline)
            if( cline%defined('endit') )then
                total_iter = nint(cline%get_rarg('endit'))
            else
                total_iter = stage_start + stage_limit - 1
            endif
            niters = total_iter - stage_start + 1
            state_overlap = read_state_overlap()
            write(logfhandle,'(A,I0,A,I0,A,F7.4,A,F7.4)') '>>> '//WORKFLOW_LABEL//' STAGE ', stage, &
                &' NITERS ', niters, ' FINAL STATE_OVERLAP: ', state_overlap, ' TARGET: ', stage_overlap_target
            if( niters >= stage_limit .and. state_overlap < stage_overlap_target )then
                write(logfhandle,'(A,I0,A,F7.4)') '>>> '//WORKFLOW_LABEL//' STAGE ', stage, &
                    &' REACHED STAGE CAP BEFORE STATE_OVERLAP TARGET: ', state_overlap
            endif
        end subroutine run_refine3D_states_stage

        subroutine ensure_all_active_particles_updated()
            integer :: nactive, nupdated, nmissing
            call read_update_coverage(nactive, nupdated, nmissing)
            if( nactive < 1 )then
                THROW_HARD(WORKFLOW_LABEL//' has no active particles after staged refinement')
            endif
            if( nmissing > 0 )then
                call run_refine3D_states_missing_update(nmissing, nactive)
                call read_update_coverage(nactive, nupdated, nmissing)
                if( nmissing > 0 )then
                    THROW_HARD(WORKFLOW_LABEL//' final missing-update pass failed to update every active particle')
                endif
            endif
        end subroutine ensure_all_active_particles_updated

        subroutine read_update_coverage( nactive, nupdated, nmissing )
            type(sp_project) :: update_proj
            integer, allocatable :: states(:), updatecnts(:)
            integer, intent(out) :: nactive, nupdated, nmissing
            call update_proj%read_segment('ptcl3D', params%projfile)
            if( .not. update_proj%os_ptcl3D%isthere('updatecnt') )then
                call update_proj%kill
                THROW_HARD(WORKFLOW_LABEL//' cannot finish before active particles are updated')
            endif
            states     = update_proj%os_ptcl3D%get_all_asint('state')
            updatecnts = update_proj%os_ptcl3D%get_all_asint('updatecnt')
            nactive    = count(states > 0)
            nupdated   = count(states > 0 .and. updatecnts > 0)
            nmissing   = nactive - nupdated
            write(logfhandle,'(A,I0,A,I0,A,I0)') &
                &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLE UPDATE COVERAGE UPDATED/ACTIVE/MISSING: ', &
                &nupdated, '/', nactive, '/', nmissing
            if( allocated(states)     ) deallocate(states)
            if( allocated(updatecnts) ) deallocate(updatecnts)
            call update_proj%kill
        end subroutine read_update_coverage

        subroutine run_refine3D_states_missing_update( nmissing, nactive )
            integer, intent(in) :: nmissing, nactive
            type(cmdline) :: cline_missing
            integer       :: iter_missing
            iter_missing = next_refine3D_states_iteration()
            write(logfhandle,'(A,A,I0,A,I0,A,I0)') &
                &'>>> '//WORKFLOW_LABEL//' FINAL MISSING-UPDATE ASSIGNMENT', &
                &' MISSING/ACTIVE/ITER: ', nmissing, '/', nactive, '/', iter_missing
            call flush(logfhandle)
            cline_missing = cline
            call cline_missing%set('prg',           'refine3D')
            call cline_missing%set('mkdir',              'no')
            call cline_missing%set('balance',            'no')
            call cline_missing%set('frac_best',           1.0)
            call cline_missing%set('fillin',             'no')
            call cline_missing%set('update_frac',         1.0)
            call cline_missing%set('trail_rec',          'no')
            call cline_missing%set('volrec',             'no')
            call cline_missing%set('maxits',                1)
            call cline_missing%set('startit',    iter_missing)
            call cline_missing%set('which_iter', iter_missing)
            call cline_missing%set('extr_iter',  iter_missing)
            call cline_missing%delete('endit')
            call cline_missing%delete('partition')
            select case(trim(params%multivol_mode))
                case('input_oris_fixed')
                    call cline_missing%set('refine', 'prob_state')
                    call cline_missing%delete('update_missing')
                    call cline_missing%delete('greedy_sampling')
                case default
                    call cline_missing%set('refine',          'greedy')
                    call cline_missing%set('greedy_sampling',   'yes')
                    call cline_missing%set('update_missing',    'yes')
            end select
            call xrefine3D%execute(cline_missing)
            call del_files(DIST_FBODY,       params%nparts, ext='.dat')
            call del_files(ASSIGNMENT_FBODY, params%nparts, ext='.dat')
            call del_file(DIST_FBODY//'.dat')
            call del_file(ASSIGNMENT_FBODY//'.dat')
            call cline_missing%kill
        end subroutine run_refine3D_states_missing_update

        integer function next_refine3D_states_iteration() result(iter)
            iter = 1
            if( cline%defined('endit') )then
                iter = cline%get_iarg('endit') + 1
            else if( cline%defined('which_iter') )then
                iter = cline%get_iarg('which_iter') + 1
            endif
            iter = max(1, iter)
        end function next_refine3D_states_iteration

        real function read_state_overlap() result(overlap)
            type(sp_project) :: state_proj
            real, allocatable :: states(:), mi_state(:), sampled(:)
            logical, allocatable :: mask(:)
            real :: sampled_lb
            overlap = 0.
            call state_proj%read_segment('ptcl3D', params%projfile)
            if( .not. state_proj%os_ptcl3D%isthere('mi_state') )then
                call state_proj%kill
                return
            endif
            states   = state_proj%os_ptcl3D%get_all('state')
            mi_state = state_proj%os_ptcl3D%get_all('mi_state')
            if( state_proj%os_ptcl3D%isthere('sampled') )then
                sampled = state_proj%os_ptcl3D%get_all('sampled')
                sampled_lb = maxval(sampled) - 0.5
                allocate(mask(size(states)), source=sampled > sampled_lb .and. states > 0.5)
            else
                allocate(mask(size(states)), source=states > 0.5)
            endif
            if( count(mask) > 0 ) overlap = sum(mi_state, mask=mask) / real(count(mask))
            if( allocated(states)   ) deallocate(states)
            if( allocated(mi_state) ) deallocate(mi_state)
            if( allocated(sampled)  ) deallocate(sampled)
            if( allocated(mask)     ) deallocate(mask)
            call state_proj%kill
        end function read_state_overlap

        subroutine cleanup_init_vols()
            integer :: state
            if( allocated(init_vols) )then
                do state = 1,size(init_vols)
                    call init_vols(state)%kill
                enddo
                deallocate(init_vols)
            endif
        end subroutine cleanup_init_vols

    end subroutine exec_refine3D_states

    subroutine exec_classify3D_refs( self, cline )
        use simple_abinitio_utils, only: write_final_rec_outputs, gen_ortho_reprojs4viz
        use simple_commanders_rec, only: commander_rec3D
        use simple_estimate_ssnr,  only: lpstages_setlims
        class(commander_classify3D_refs), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer, parameter :: NSAMPLE_PER_STATE_CLASSIFY3D_REFS = 10000
        integer, parameter :: NSAMPLE_CLASSIFY3D_REFS_CAP = 100000
        integer, parameter :: NSPACE                      = 5000
        integer, parameter :: NSPACE_SUB                  = 500
        integer, parameter :: MINITS_CLASSIFY3D_REFS      = 10
        integer, parameter :: MAXITS_CLASSIFY3D_REFS      = 50
        real,    parameter :: TARGET_UPDATES_PER_PARTICLE = 4.0
        real,    parameter :: LPSTART_CLASSIFY3D_REFS     = 10.0
        real,    parameter :: LPSTOP_CLASSIFY3D_REFS      = 6.0
        character(len=*), parameter :: WORKFLOW_LABEL     = 'CLASSIFY3D_REFS'
        type(commander_rec3D)     :: xrec3D
        type(commander_refine3D)  :: xrefine3D
        type(cmdline)             :: cline_rec3D
        type(parameters)          :: params, params_final_rec
        type(sp_project)          :: spproj
        type(lp_crop_inf)         :: lpinfo_master(1)
        type(string), allocatable :: init_vols(:)
        integer,      allocatable :: state_pops(:)
        integer :: nptcls_eff, iter_glob, state
        call cline%set('prg', 'classify3D_refs')
        ! hard defaults
        call cline%set('balance',         'yes')
        call cline%set('greedy_sampling', 'no')
        call cline%set('frac_best',       1.0)
        call cline%set('trail_rec',       'yes')
        call cline%set('objfun',          'euclid')
        call cline%set('lplim_crit',      0.5)
        call cline%set('incrreslim',      'no')
        call cline%set('nu_refine',       'no')
        call cline%set('combine_eo',      'no')
        call cline%set('multivol_mode',   'independent')
        ! overridable defaults
        if( .not. cline%defined('filt_mode')       ) call cline%set('filt_mode',       'nonuniform_lpset')
        if( .not. cline%defined('envfsc')          ) call cline%set('envfsc',          'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',           'yes')
        if( .not. cline%defined('center')          ) call cline%set('center',          'no')
        call cline%set('sigma_est', 'global')
        if( .not. cline%defined('prob_inpl')       ) call cline%set('prob_inpl',       'yes')
        if( .not. cline%defined('refine')          ) call cline%set('refine',          'prob_neigh')
        if( .not. cline%defined('prob_neigh_mode') ) call cline%set('prob_neigh_mode', 'state')
        ! staged frequency marching always uses planner-driven downscaling
        if( cline%defined('autoscale') ) write(logfhandle,'(A)') &
            &'>>> '//WORKFLOW_LABEL//' IGNORES autoscale: frequency marching owns downscaling'
        call cline%set('autoscale', 'yes')
        if( .not. cline%defined('ml_reg')          ) call cline%set('ml_reg',          'yes')
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart', LPSTART_CLASSIFY3D_REFS)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',  LPSTOP_CLASSIFY3D_REFS)
        if( .not. cline%defined('automsk')         ) call cline%set('automsk',         'no')
        if( .not. cline%defined('nsample')         ) call cline%set('nsample', NSAMPLE_PER_STATE_CLASSIFY3D_REFS)
        if( .not. cline%defined('keepvol')         ) call cline%set('keepvol',         'no')
        if( .not. cline%defined('nspace')          ) call cline%set('nspace',          NSPACE)
        if( .not. cline%defined('nspace_sub')      ) call cline%set('nspace_sub',      NSPACE_SUB)
        call params%new(cline)
        if( params%nstates < 2 ) THROW_HARD('nstates must be >= 2 for '//WORKFLOW_LABEL)
        call cline%set('mkdir', 'no')
        call spproj%read( params%projfile )
        ! Search planning
        iter_glob = 0
        call set_classify3D_refs_nstates(nptcls_eff, state_pops)
        call set_classify3D_refs_nsample
        call validate_classify3D_refs_filtering
        call validate_classify3D_refs_search_mode
        call set_classify3D_refs_sampling
        call prepare_classify3D_refs_class_sampling
        call spproj%os_ptcl2D%kill
        call set_classify3D_refs_downscaling
        if( spproj%os_ptcl3D%has_been_sampled() )then
            call spproj%os_ptcl3D%clean_entry('sampled', 'updatecnt')
        endif
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        call initialize_state_volumes
        ! Search
        call cline%set('prg',         'refine3D')
        call cline%set('maxits_glob', max(1, params%maxits_glob))
        call cline%set('balance',     params%balance)
        call run_classify3D_refs
        ! Final particle mapping pass
        call ensure_all_active_particles_updated
        ! Final reconstruction
        call reconstruct_all_particles_volumes
        if( allocated(state_pops) ) deallocate(state_pops)
        if( allocated(init_vols) )then
            do state = 1,size(init_vols)
                call init_vols(state)%kill
            enddo
            deallocate(init_vols)
        endif
        call spproj%kill
        call simple_end('**** SIMPLE_CLASSIFY3D_REFS NORMAL STOP ****')
    contains

        subroutine set_classify3D_refs_nstates( nptcls_eff, state_pops )
            integer,              intent(inout) :: nptcls_eff
            integer, allocatable, intent(inout) :: state_pops(:)
            integer :: state, nstates_labels
            nptcls_eff = spproj%os_ptcl3D%count_state_gt_zero()
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//WORKFLOW_LABEL)
            nstates_labels = spproj%os_ptcl3D%get_n('state')
            if( nstates_labels == 1 )then
                allocate(state_pops(params%nstates), source=0)
                state_pops(1) = nptcls_eff
            elseif( nstates_labels == params%nstates )then
                call spproj%os_ptcl3D%get_pops(state_pops, 'state', maxn=params%nstates)
                write(logfhandle,'(A,I0)') '>>> '//WORKFLOW_LABEL//' NSTATES FROM PROJECT: ', nstates_labels
            else
                THROW_HARD('command-line nstates does not match project state labels for '//WORKFLOW_LABEL)
            endif
            do state = 1, size(state_pops)
                write(logfhandle,'(A,I3,A,I0)') '>>> '//WORKFLOW_LABEL//' STATE POPULATION: ', state, '/', state_pops(state)
            enddo
        end subroutine set_classify3D_refs_nstates

        subroutine set_classify3D_refs_nsample()
            params%nsample = min(NSAMPLE_CLASSIFY3D_REFS_CAP, params%nsample * params%nstates)
            params%nsample = min(nptcls_eff, params%nsample)
            if( params%nsample < 1 ) THROW_HARD('nsample must be >= 1 for '//WORKFLOW_LABEL)
            write(logfhandle,'(A,I0,A,I0)')'>>> '//WORKFLOW_LABEL//' NSAMPLE: ', params%nsample
        end subroutine set_classify3D_refs_nsample

        subroutine validate_classify3D_refs_filtering()
            select case(trim(params%filt_mode))
                case('nonuniform_lpset', 'none', 'nonuniform')
                    ! supported
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports filt_mode=nonuniform_lpset|none')
            end select
        end subroutine validate_classify3D_refs_filtering

        subroutine validate_classify3D_refs_search_mode()
            select case(trim(params%refine))
                case('prob_neigh', 'prob', 'shc')
                    call cline%set('refine', params%refine)
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports refine=prob_neigh|prob|shc')
            end select
            call cline%set('refine', params%refine)
            if( trim(params%refine).eq.'prob_neigh' )then
                select case(trim(params%prob_neigh_mode))
                case('geom', 'state')
                    call cline%set('prob_neigh_mode', params%prob_neigh_mode)
                case default
                    THROW_HARD(WORKFLOW_LABEL//' supports prob_neigh_mode=geom|state')
                end select
            endif
        end subroutine validate_classify3D_refs_search_mode

        subroutine set_classify3D_refs_sampling()
            real    :: update_frac_auto
            integer :: maxits_auto, nptcls_per_iter, stage_cap
            if( nptcls_eff < 1 ) THROW_HARD('no active particles available for '//WORKFLOW_LABEL)
            nptcls_per_iter = min(nptcls_eff, params%nsample)
            if( nptcls_eff <= params%nsample )then
                params%update_frac   = 1.0
                params%l_update_frac = .false.
                params%l_trail_rec   = .false.
                call cline%delete('update_frac')
                write(logfhandle,'(A,I0,A,I0,A)')'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                    &nptcls_eff, '/', params%nsample, ' -> FULL UPDATE'
            else
                update_frac_auto = real(params%nsample) / real(nptcls_eff)
                if( update_frac_auto <= 0.99 )then
                    params%update_frac   = update_frac_auto
                    params%l_update_frac = .true.
                    params%l_trail_rec   = trim(params%trail_rec).eq.'yes'
                    call cline%set('update_frac', update_frac_auto)
                    write(logfhandle,'(A,I0,A,I0,A,F8.4)')'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC: ', &
                        &nptcls_eff, '/', params%nsample, '/', update_frac_auto
                else
                    params%update_frac   = 1.0
                    params%l_update_frac = .false.
                    params%l_trail_rec   = .false.
                    call cline%delete('update_frac')
                    write(logfhandle,'(A,I0,A,I0,A)') &
                        &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLES/SAMPLE TARGET: ', &
                        &nptcls_eff, '/', params%nsample, ' -> FULL UPDATE'
                endif
            endif
            if( cline%defined('maxits') )then
                write(logfhandle,'(A,I0)') &
                    &'>>> '//WORKFLOW_LABEL//' STAGE MAXITS COMMAND-LINE OVERRIDE: ', params%maxits
            else
                maxits_auto   = ceiling((TARGET_UPDATES_PER_PARTICLE * real(nptcls_eff)) / real(nptcls_per_iter))
                stage_cap     = max(MINITS_CLASSIFY3D_REFS, min(MAXITS_CLASSIFY3D_REFS, maxits_auto))
                params%maxits = stage_cap
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0,A,F5.1,A)') '>>> '//WORKFLOW_LABEL//' STAGE MAXITS: ', &
                    &stage_cap, ' FOR ~', TARGET_UPDATES_PER_PARTICLE, ' UPDATES/PARTICLE'
            endif
        end subroutine set_classify3D_refs_sampling

        subroutine prepare_classify3D_refs_class_sampling()
            type(class_sample), allocatable :: clssmp(:)
            integer, allocatable :: tmpinds(:), clsinds(:), cls_states(:)
            integer              :: icls
            if( spproj%is_virgin_field('ptcl2D') )then
                params%balance = 'no'
            else
                if( params%update_frac > 0.99 )then
                    write(logfhandle,'(A)') '>>> FORCING FULL ACTIVE SAMPLING (NO FRACTIONAL OR TRAILING UPDATE)'
                    params%balance = 'no'
                else
                    ! generate a data structure for class sampling on disk
                    if( trim(params%balance).eq.'yes' )then
                        if( trim(params%partition).eq.'yes' )then
                            if( .not. spproj%os_cls2D%isthere('cluster') )then
                                THROW_HARD('Missing CLUSTER metadata in CLS2D field needed for PARTITION=YES')
                            endif
                            cls_states = nint(spproj%os_cls2D%get_all('state'))
                            tmpinds    = nint(spproj%os_cls2D%get_all('cluster'))
                            where( cls_states == 0 ) tmpinds = 0
                            clsinds = (/(icls,icls=1,maxval(tmpinds))/)
                            do icls = 1,size(clsinds)
                                if(count(tmpinds==icls) == 0) clsinds(icls) = 0
                            enddo
                            clsinds = pack(clsinds, mask=clsinds>0)
                            call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp, label='cluster')
                            deallocate(cls_states,tmpinds)
                        else
                            clsinds = spproj%get_selected_clsinds()
                            call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
                        endif
                        call write_class_samples(clssmp, string(CLASS_SAMPLING_FILE))
                        call deallocate_class_samples(clssmp)
                        deallocate(clsinds)
                        write(logfhandle,'(A)') '>>> SETUP 2D DERIVED CLASS SAMPLING'
                    endif
                endif
            endif
        end subroutine prepare_classify3D_refs_class_sampling

        subroutine set_classify3D_refs_downscaling()
            lpinfo_master(1)%trslim      = min(8.,max(2.0, AHELIX_WIDTH/params%smpd))
            lpinfo_master(1)%box_crop    = params%box
            lpinfo_master(1)%smpd_crop   = params%smpd
            lpinfo_master(1)%l_autoscale = .false.
            call lpstages_setlims(params%box, 1, params%smpd, params%lpstart, params%lpstop, lpinfo_master(1))
            params%trs = lpinfo_master(1)%trslim
            call cline%set('trs', params%trs)
            if( lpinfo_master(1)%l_autoscale )then
                params%box_crop  = lpinfo_master(1)%box_crop
                params%smpd_crop = lpinfo_master(1)%smpd_crop
                call cline%set('box_crop', params%box_crop)
                call cline%set('smpd_crop', params%smpd_crop)
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> '//WORKFLOW_LABEL//' DOWNSCALING BOX/SMPD_CROP: ', &
                    &params%box, '/', params%box_crop, '/', params%smpd_crop
            else
                call cline%delete('box_crop')
                call cline%delete('smpd_crop')
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' DOWNSCALING: native sampling retained'
            endif
        end subroutine set_classify3D_refs_downscaling

        subroutine initialize_state_volumes()
            type(string), allocatable :: checkpoint_vols(:)
            logical :: vols_defined(params%nstates)
            integer :: state
            do state = 1,params%nstates
                vols_defined(state) = cline%defined('vol'//int2str(state))
            enddo
            if( .not. all(vols_defined) )then
                THROW_HARD(WORKFLOW_LABEL//' requires a complete external vol1..volN reference set')
            endif
            allocate(init_vols(params%nstates), checkpoint_vols(params%nstates))
            call validate_input_volumes()
            init_vols(1:params%nstates) = params%vols(1:params%nstates)
            write(logfhandle,'(A)')'>>> '//WORKFLOW_LABEL//&
                &' EXTERNAL REFERENCES ARE UNTRUSTED UNTIL CC POSE INITIALIZATION COMPLETES'
            write(logfhandle,'(A,F7.2)') '>>> '//WORKFLOW_LABEL//' INITIALIZING POSES WITH LP:', params%lpstart
            call initialize_poses_against_external_references(params, cline, xrefine3D, xrec3D, &
                &nptcls_eff, init_vols, checkpoint_vols, 1, lp_init=params%lpstart)
            params%vols(1:params%nstates) = checkpoint_vols
            iter_glob = 1
            ! update command line
            do state = 1,params%nstates
                call cline%set('vol'//int2str(state), params%vols(state))
            enddo
            deallocate(checkpoint_vols)
        end subroutine initialize_state_volumes

        subroutine validate_input_volumes()
            type(string) :: vol
            integer :: state, ldim(3), nptcls_dummy
            real    :: vol_smpd
            do state = 1,params%nstates
                vol = params%vols(state)
                if( .not. file_exists(vol) ) THROW_HARD('Input volume does not exist: '//vol%to_char())
                call find_ldim_nptcls(vol, ldim, nptcls_dummy)
                vol_smpd = find_img_smpd(vol)
                if( any(ldim /= [params%box,params%box,params%box]) .or. abs(vol_smpd - params%smpd) > 1.e-6 )then
                    THROW_HARD('Input state volumes must have same dimensions/sampling as the project particles')
                endif
            enddo
            call vol%kill
        end subroutine validate_input_volumes

        subroutine run_classify3D_refs()
            integer, parameter  :: STEP_IT = 3
            type(refine3D_stage_plan_entry), allocatable :: stages(:)
            integer :: startit, lastit, stage
            call plan_refine3D_frequency_stages(params%box, params%smpd, params%lpstart, params%lpstop, &
                &params%maxits, STEP_IT, iter_glob + 1, lpinfo_master(1), stages)
            write(logfhandle,'(A,I0,A,F7.2,A,F7.2)') &
                &'>>> '//WORKFLOW_LABEL//' FREQUENCY MARCHING NSTAGES/LPSTART/LPSTOP: ',&
                &size(stages), '/', real(params%lpstart), '/', real(params%lpstop)
            do stage = 1,size(stages)
                write(logfhandle,'(A,I3,A,F7.2)')'>>> '//WORKFLOW_LABEL//' PLANNED STAGE/LP: ',&
                    &stage, ' / ', stages(stage)%lpinfo%lp
            enddo
            lastit = 0
            do stage = 1,size(stages)
                startit = iter_glob + 1
                call cline%set('minits',     1)
                call cline%set('maxits',     stages(stage)%nits)
                call cline%set('startit',    startit)
                call cline%set('which_iter', startit)
                call cline%set('extr_iter',  startit)
                call cline%set('trs',        stages(stage)%lpinfo%trslim)
                call cline%set('lp',         stages(stage)%lpinfo%lp)
                call cline%set('lpstop',     stages(stage)%lpinfo%lp)
                if( stages(stage)%lpinfo%l_autoscale )then
                    call cline%set('box_crop',  stages(stage)%lpinfo%box_crop)
                    call cline%set('smpd_crop', stages(stage)%lpinfo%smpd_crop)
                else
                    call cline%delete('box_crop')
                    call cline%delete('smpd_crop')
                endif
                write(logfhandle,'(A,I0,A,F7.2)')'>>> '//WORKFLOW_LABEL//' ENTERING STAGE ', &
                    &stage, ' LP: ', stages(stage)%lpinfo%lp
                call xrefine3D%execute(cline)
                lastit    = cline%get_iarg('endit')
                iter_glob = lastit
                call cline%delete('endit')
            enddo
            call del_files(DIST_FBODY,       params%nparts, ext='.dat')
            call del_files(ASSIGNMENT_FBODY, params%nparts, ext='.dat')
            call del_file(DIST_FBODY//'.dat')
            call del_file(ASSIGNMENT_FBODY//'.dat')
            write(logfhandle,'(A,I0)') '>>> '//WORKFLOW_LABEL//' EXITING FREQUENCY MARCHING AT ITERATION ', iter_glob
            deallocate(stages)
        end subroutine run_classify3D_refs

        subroutine ensure_all_active_particles_updated()
            integer :: nactive, nupdated, nmissing
            call read_update_coverage(nactive, nupdated, nmissing)
            if( nactive < 1 )then
                THROW_HARD(WORKFLOW_LABEL//' has no active particles after staged refinement')
            endif
            if( nmissing > 0 )then
                call run_classify3D_refs_missing_update(nmissing, nactive)
                call read_update_coverage(nactive, nupdated, nmissing)
                if( nmissing > 0 )then
                    THROW_HARD(WORKFLOW_LABEL//' final missing-update pass failed to update every active particle')
                endif
            endif
        end subroutine ensure_all_active_particles_updated

        subroutine read_update_coverage( nactive, nupdated, nmissing )
            type(sp_project) :: update_proj
            integer, allocatable :: states(:), updatecnts(:)
            integer, intent(out) :: nactive, nupdated, nmissing
            call update_proj%read_segment('ptcl3D', params%projfile)
            if( .not. update_proj%os_ptcl3D%isthere('updatecnt') )then
                call update_proj%kill
                THROW_HARD(WORKFLOW_LABEL//' cannot finish before active particles are updated')
            endif
            states     = update_proj%os_ptcl3D%get_all_asint('state')
            updatecnts = update_proj%os_ptcl3D%get_all_asint('updatecnt')
            nactive    = count(states > 0)
            nupdated   = count(states > 0 .and. updatecnts > 0)
            nmissing   = nactive - nupdated
            write(logfhandle,'(A,I0,A,I0,A,I0)') &
                &'>>> '//WORKFLOW_LABEL//' ACTIVE PARTICLE UPDATE COVERAGE UPDATED/ACTIVE/MISSING: ', &
                &nupdated, '/', nactive, '/', nmissing
            if( allocated(states)     ) deallocate(states)
            if( allocated(updatecnts) ) deallocate(updatecnts)
            call update_proj%kill
        end subroutine read_update_coverage

        subroutine run_classify3D_refs_missing_update( nmissing, nactive )
            integer, intent(in) :: nmissing, nactive
            type(cmdline) :: cline_missing
            iter_glob = iter_glob + 1
            write(logfhandle,'(A,A,I0,A,I0,A,I0)') &
                &'>>> '//WORKFLOW_LABEL//' FINAL MISSING-UPDATE ASSIGNMENT', &
                &' MISSING/ACTIVE/ITER: ', nmissing, '/', nactive, '/', iter_glob
            call flush(logfhandle)
            cline_missing = cline
            call cline_missing%set('prg',           'refine3D')
            call cline_missing%set('mkdir',              'no')
            call cline_missing%set('balance',            'no')
            call cline_missing%set('frac_best',           1.0)
            call cline_missing%set('fillin',             'no')
            call cline_missing%set('update_frac',         1.0)
            call cline_missing%set('trail_rec',          'no')
            call cline_missing%set('volrec',             'no')
            call cline_missing%set('maxits',                1)
            call cline_missing%set('startit',       iter_glob)
            call cline_missing%set('which_iter',    iter_glob)
            call cline_missing%set('extr_iter',     iter_glob)
            call cline_missing%set('refine',         'greedy')
            call cline_missing%set('greedy_sampling',   'yes')
            call cline_missing%set('update_missing',    'yes')
            call cline_missing%delete('endit')
            call cline_missing%delete('partition')
            call xrefine3D%execute(cline_missing)
            call cline_missing%kill
        end subroutine run_classify3D_refs_missing_update

        subroutine reconstruct_all_particles_volumes
            cline_rec3D = cline
            call cline_rec3D%set('prg',            'reconstruct3D')
            call cline_rec3D%set('outfile', 'RESOLUTION_FINAL.txt')
            call cline_rec3D%set('postprocess',              'yes')
            call cline_rec3D%delete('trail_rec')
            call cline_rec3D%delete('refine')
            call cline_rec3D%delete('objfun_den')
            call cline_rec3D%delete('objfun_den_w')
            call cline_rec3D%delete('sigma_est')
            call cline_rec3D%delete('update_frac')
            call cline_rec3D%delete('ufrac_trec')
            call cline_rec3D%delete('endit')
            call cline_rec3D%delete('box_crop')
            call cline_rec3D%delete('smpd_crop')
            call cline_rec3D%set('objfun', 'cc')
            if( params%l_nonuniform )then
                call cline_rec3D%set('filt_mode', 'none')
                call cline_rec3D%set('automsk', 'no')
            endif
            call cline_rec3D%set('nu_refine', 'no')
            call xrec3D%execute(cline_rec3D)
            call params_final_rec%new(cline_rec3D)
            params_final_rec%box  = params_final_rec%box_crop
            params_final_rec%smpd = params_final_rec%smpd_crop
            call spproj%read_segment('out', params_final_rec%projfile)
            call write_final_rec_outputs(params_final_rec, spproj, params_final_rec%lpstop)
            call gen_ortho_reprojs4viz(params_final_rec, spproj)
        end subroutine reconstruct_all_particles_volumes

    end subroutine exec_classify3D_refs

    !> Single entrypoint (shared-memory OR distributed master), driven by a strategy.
    subroutine exec_refine3D( self, cline )
        use simple_core_module_api
        use simple_refine3D_strategy
        class(commander_refine3D), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(refine3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        type(string)     :: filt_mode_arg
        logical          :: converged
        integer          :: niters
        ! sanity check: multiple input volumes require nstates > 1
        if( cline%defined('vol2') )then
            if( .not. cline%defined('nstates') )then
                THROW_HARD('Multiple volumes (vol1, vol2, ...) provided on command line but NSTATES is not set; set NSTATES to the number of states')
            else if( cline%get_iarg('nstates') <= 1 )then
                THROW_HARD('Multiple volumes (vol1, vol2, ...) provided on command line but NSTATES <= 1; set NSTATES to the number of states')
            endif
        endif
        if( cline%defined('filt_mode') .and. cline%defined('nstates') )then
            filt_mode_arg = cline%get_carg('filt_mode')
            if( cline%get_iarg('nstates') > 1 .and. trim(filt_mode_arg%to_char()).eq.'uniform' )then
                THROW_HARD('filt_mode=uniform is disabled for multi-state refine3D search')
            endif
            call filt_mode_arg%kill
        endif
        ! local defaults (kept consistent with previous distributed master)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('prg', 'refine3D')
        ! Select execution strategy (shared-memory vs distributed master)
        strategy = create_refine3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        if( params%nstates > 1 .and. trim(params%filt_mode).eq.'uniform' )then
            THROW_HARD('filt_mode=uniform is disabled for multi-state refine3D search')
        endif
        ! Main loop counter semantics:
        !   - params%maxits is the *number of iterations to run* in this invocation.
        !   - params%which_iter starts at params%startit.
        niters            = 0
        params%which_iter = params%startit - 1
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        do
            niters            = niters + 1
            params%which_iter = params%which_iter + 1
            params%extr_iter  = params%extr_iter  + 1
            call strategy%execute_iteration(params, build, cline, converged)
            call strategy%finalize_iteration(params, build)
            if( converged .or. niters >= params%maxits ) exit
        end do
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        ! Global teardown (strategies may have built different toolboxes)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call build%pftc%kill
        call simple_end('**** SIMPLE_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D

    !> Distributed worker (single-iteration execution). This should be the command
    !> invoked by the scheduler for each partition.
    subroutine exec_refine3D_distr_worker( self, cline )
        use simple_core_module_api
        use simple_strategy3D_matcher, only: refine3D_exec
        class(commander_refine3D_distr_worker), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        logical          :: converged
        logical          :: l_write_partial_recs
        ! Flags required for worker execution
        if( .not. cline%defined('part')    ) THROW_HARD('PART must be defined for distributed worker execution')
        if( .not. cline%defined('outfile') ) THROW_HARD('OUTFILE must be defined for distributed worker execution')
        ! Worker needs the alignment toolboxes
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        if( params%which_iter < 1 ) params%which_iter = max(1, params%startit)
        if( .not. cline%defined('extr_iter') ) params%extr_iter = params%which_iter
        call cline%set('which_iter', int2str(params%which_iter))
        l_write_partial_recs = trim(params%volrec) .eq. 'yes'
        call refine3D_exec(params, build, cline, params%which_iter, converged, l_write_partial_recs)
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
    end subroutine exec_refine3D_distr_worker

end module simple_commanders_refine3D

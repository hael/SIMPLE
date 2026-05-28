!@descr: supporting 3D orientation search
module simple_commanders_refine3D
use simple_commanders_api
use simple_pftc_srch_api
use simple_refine3D_fnames,   only: refine3D_state_vol_fname
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
        use simple_abinitio_utils, only: write_final_rec_outputs
        use simple_commanders_rec, only: commander_rec3D
        use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
            &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, &
            &extend_nu_filter_highres_shells, write_nu_local_resolution_map
        class(commander_refine3D_auto), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(cmdline)               :: cline_rec3D
        type(parameters)            :: params, params_final_rec
        type(sp_project)            :: spproj
        type(string)                :: init_vol
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
        ! commanders
        type(commander_rec3D)    :: xrec3D
        type(commander_refine3D) :: xrefine3D
        maxits_user      = 0
        ! hard defaults
        call cline%set('balance',         'no') ! no balancing based on 2D clustering
        call cline%set('greedy_sampling', 'no') ! only active when balance is 'yes'`
        call cline%set('trail_rec',      'yes') ! trailing average 3D reconstruction
        call cline%set('refine', 'prob_neigh') ! probabilistioc neighborhood 3D refinement
        call cline%set('ml_reg',         'yes') ! ML regularization is on
        call cline%set('overlap',         0.99) ! convergence if overlap > 99%
        call cline%set('nstates',            1) ! only single-state refinement is supported
        call cline%set('objfun',      'euclid') ! the objective function is noise-normalized Euclidean distance
        call cline%set('envfsc',         'yes') ! we use the envelope mask when calculating the FSC
        call cline%set('lplim_crit',     0.143) ! we use the 0.143 criterion for low-pass limitation
        call cline%set('incrreslim',      'no') ! if anything 'yes' makes it slightly worse, but no real difference right now
        ! overridable defaults
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',            'no') ! 4 now, probably fine
        if( .not. cline%defined('sigma_est')   ) call cline%set('sigma_est',     'global') ! 4 now, probably fine
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',        'no') ! 4 now, to allow more rapid testing
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',        'yes') ! no difference at this stage, so prefer 'yes'
        if( .not. cline%defined('nsample')     ) call cline%set('nsample', NSAMPLE_REFINE3D_AUTO)
        if( .not. cline%defined('autoscale')   ) call cline%set('autoscale',      'yes')
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',           'yes') ! better map with ml_reg='yes'
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode', 'nonuniform') ! obvioulsy
        if( .not. cline%defined('nu_refine')   ) call cline%set('nu_refine',        'yes') ! allow conservative NU resolution-bank expansion
        if( .not. cline%defined('automsk')     ) call cline%set('automsk',          'yes') ! envelope masking for background flattening
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
            l_have_init_vol = .true.
            write(logfhandle,'(A,1X,A)') '>>> '//WORKFLOW_LABEL//' USING INPUT VOLUME:', init_vol%to_char()
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
        ! generate an initial 3D reconstruction
        cline_rec3D = cline
        call cline_rec3D%set('prg', 'reconstruct3D') ! required for distributed call
        call cline_rec3D%delete('trail_rec')
        call cline_rec3D%delete('objfun')
        call cline_rec3D%delete('sigma_est')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%set('objfun', 'cc') ! ugly, but this is how it works in parameters
        call cline_rec3D%set('postprocess', 'no')
        call cline_rec3D%set('nu_refine', 'no')
        if( l_have_init_vol ) call prepare_nu_bootstrap_refs_from_raw_halves()
        if( l_have_init_vol )then
            call cline%set('vol1', init_vol)
        else
            call xrec3D%execute(cline_rec3D)
            call cline%set('vol1', refine3D_state_vol_fname(1))
        endif
        ! 3D refinement iterations
        call cline%set('prg',                   'refine3D')
        call cline%set('ufrac_trec',    params%update_frac)
        call cline%set('maxits',             params%maxits)
        call xrefine3D%execute(cline)
        ! re-reconstruct from all particle images
        call cline_rec3D%set('outfile', 'RESOLUTION_FINAL.txt')
        call cline_rec3D%set('postprocess', 'yes')
        if( params%l_nonuniform )then
            call cline_rec3D%set('filt_mode', 'none')
        endif
        call cline_rec3D%set('nu_refine', 'no')
        call xrec3D%execute(cline_rec3D)
        call params_final_rec%new(cline_rec3D)
        params_final_rec%box  = params_final_rec%box_crop
        params_final_rec%smpd = params_final_rec%smpd_crop
        call spproj%read_segment('out', params_final_rec%projfile)
        call write_final_rec_outputs(params_final_rec, spproj, params_final_rec%res_target)
        call spproj%kill
        call init_vol%kill

    contains

        logical function project_init_vol_compatible() result( l_compatible )
            l_compatible = init_box == params%box .and. init_smpd > TINY .and. &
                &abs(init_smpd - params%smpd) <= 1.e-6
        end function project_init_vol_compatible

        subroutine prepare_nu_bootstrap_refs_from_raw_halves()
            type(string)         :: init_even, init_odd, raw_even, raw_odd, candidate
            type(string)         :: out_even, out_odd, out_avg, out_locres
            type(image)          :: vol_even_raw, vol_odd_raw, vol_even_nu, vol_odd_nu, vol_msk
            type(image_msk)      :: envmsk
            logical, allocatable :: l_mask(:,:,:)
            integer, allocatable :: imat(:,:,:)
            integer              :: ldim_even(3), ldim_odd(3), ldim(3), nptcls_dummy, n_bootstrap_steps
            real                 :: mskrad_px
            logical              :: l_reconstruct_bootstrap
            if( .not. params%l_nonuniform ) return
            l_reconstruct_bootstrap = .false.
            init_even      = add2fbody(init_vol, MRC_EXT, '_even')
            init_odd       = add2fbody(init_vol, MRC_EXT, '_odd')
            candidate      = add2fbody(init_even, MRC_EXT, '_unfil')
            if( file_exists(candidate) )then
                raw_even = candidate
            else
                raw_even = init_even
            endif
            call candidate%kill
            candidate      = add2fbody(init_odd, MRC_EXT, '_unfil')
            if( file_exists(candidate) )then
                raw_odd = candidate
            else
                raw_odd = init_odd
            endif
            call candidate%kill
            if( .not. file_exists(raw_even) .or. .not. file_exists(raw_odd) )then
                write(logfhandle,'(A)') &
                    &'>>> '//WORKFLOW_LABEL//' BOOTSTRAP: raw native E/O pair missing; '//&
                    &'reconstructing startup references'
                l_reconstruct_bootstrap = .true.
            else
                call find_ldim_nptcls(raw_even, ldim_even, nptcls_dummy)
                call find_ldim_nptcls(raw_odd,  ldim_odd,  nptcls_dummy)
                if( any(ldim_even /= [params%box,params%box,params%box]) .or. any(ldim_odd /= ldim_even) )then
                    write(logfhandle,'(A)') &
                        &'>>> '//WORKFLOW_LABEL//' BOOTSTRAP: raw E/O dimensions incompatible; '//&
                        &'reconstructing startup references'
                    l_reconstruct_bootstrap = .true.
                else
                    ldim = ldim_even
                    call vol_even_raw%new(ldim, params%smpd)
                    call vol_odd_raw%new(ldim, params%smpd)
                    call vol_even_raw%read(raw_even)
                    call vol_odd_raw%read(raw_odd)
                    if( abs(vol_even_raw%get_smpd() - params%smpd) > 1.e-6 .or. &
                        &abs(vol_odd_raw%get_smpd()  - params%smpd) > 1.e-6 )then
                        write(logfhandle,'(A)') &
                            &'>>> '//WORKFLOW_LABEL//' BOOTSTRAP: raw E/O sampling incompatible; '//&
                            &'reconstructing startup references'
                        l_reconstruct_bootstrap = .true.
                    endif
                endif
            endif
            if( l_reconstruct_bootstrap )then
                l_have_init_vol = .false.
                call cline_rec3D%delete('vol1')
                call init_even%kill
                call init_odd%kill
                call raw_even%kill
                call raw_odd%kill
                call candidate%kill
                call vol_even_raw%kill
                call vol_odd_raw%kill
                call vol_msk%kill
                call envmsk%kill
                if( allocated(l_mask) ) deallocate(l_mask)
                if( allocated(imat)   ) deallocate(imat)
                return
            endif
            if( params%automsk .ne. 'no' )then
                call envmsk%automask3D(params, vol_even_raw, vol_odd_raw, l_tight=params%automsk.eq.'tight')
                call envmsk%set_imat
                call envmsk%get_imat(imat)
                allocate(l_mask(ldim(1),ldim(2),ldim(3)))
                l_mask = imat > 0
                deallocate(imat)
            else
                mskrad_px = 0.5 * params%mskdiam / params%smpd
                call vol_msk%disc(ldim, params%smpd, mskrad_px, l_mask)
            endif
            call setup_nu_dmats(vol_even_raw, vol_odd_raw, l_mask, [real ::])
            if( allocated(l_mask) ) deallocate(l_mask)
            call optimize_nu_cutoff_finds()
            if( params%l_nu_refine )then
                call extend_nu_filter_highres_shells(vol_even_raw, vol_odd_raw, nsteps=n_bootstrap_steps)
                write(logfhandle,'(A,I0)') '>>> NU bootstrap accepted high-resolution shell steps: ', n_bootstrap_steps
            endif
            call nu_filter_vols(vol_even_nu, vol_odd_nu)
            call print_nu_filtmap_lowpass_stats()
            call analyze_filtmap_neighbor_continuity()
            out_even = add2fbody(init_even, MRC_EXT, NUFILT_SUFFIX)
            out_odd  = add2fbody(init_odd,  MRC_EXT, NUFILT_SUFFIX)
            out_avg  = add2fbody(init_vol,  MRC_EXT, NUFILT_SUFFIX)
            out_locres = add2fbody(init_vol, MRC_EXT, NULOCRES_SUFFIX)
            call vol_even_nu%write(out_even, del_if_exists=.true.)
            call vol_odd_nu%write(out_odd, del_if_exists=.true.)
            call vol_even_nu%add(vol_odd_nu)
            call vol_even_nu%mul(0.5)
            call vol_even_nu%write(out_avg, del_if_exists=.true.)
            call write_nu_local_resolution_map(out_locres)
            call wait_for_closure(out_avg)
            call wait_for_closure(out_locres)
            write(logfhandle,'(A)') &
                &'>>> '//WORKFLOW_LABEL//' BOOTSTRAP: generated NU-filtered startup references from raw native E/O maps'
            call cleanup_nu_filter()
            call init_even%kill
            call init_odd%kill
            call raw_even%kill
            call raw_odd%kill
            call out_even%kill
            call out_odd%kill
            call out_avg%kill
            call out_locres%kill
            call candidate%kill
            call vol_even_raw%kill
            call vol_odd_raw%kill
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call vol_msk%kill
            call envmsk%kill
            if( allocated(l_mask) ) deallocate(l_mask)
            if( allocated(imat)   ) deallocate(imat)
        end subroutine prepare_nu_bootstrap_refs_from_raw_halves

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

    end subroutine exec_refine3D_auto

    !> Single entrypoint (shared-memory OR distributed master), driven by a strategy.
    subroutine exec_refine3D( self, cline )
        use simple_core_module_api
        use simple_refine3D_strategy
        class(commander_refine3D), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(refine3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
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
        ! local defaults (kept consistent with previous distributed master)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%set('prg', 'refine3D')
        ! Select execution strategy (shared-memory vs distributed master)
        strategy = create_refine3D_strategy(cline)
        call strategy%initialize(params, build, cline)
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

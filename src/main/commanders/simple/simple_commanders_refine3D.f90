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

type, extends(commander_base) :: commander_refine3D_multi
  contains
    procedure :: execute      => exec_refine3D_multi
end type commander_refine3D_multi

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
            call prepare_external_init_vol(init_vol)
            l_have_init_vol = .true.
            call cline%delete('vol1') ! because now part of the project
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
        call cline_rec3D%delete('box_crop')           ! original image dimensions
        call cline_rec3D%delete('smpd_crop')
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
        params_final_rec%box  = params%box
        params_final_rec%smpd = params%smpd
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
            call optimize_nu_cutoff_finds(histogram_potts=params%l_nu_hist_potts)
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
            v = stdev * real(init_box)
            call vol_even%norm_ext(ave, v)
            call vol_odd%stats('foreground', ave, stdev, maxv, minv, msk=msk)
            v = stdev * real(init_box)
            call vol_odd%norm_ext(ave, v)
            call vol%stats('foreground', ave, stdev, maxv, minv, msk=msk)
            v = stdev * real(init_box)
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
            call spproj%read_segment('out', params%projfile)
            ! cleanup
            deallocate(fsc)
            call spproj%kill
            call vol%kill; call vol_even%kill; call vol_odd%kill
            call init_even%kill; call init_odd%kill
            call new_vol%kill; call new_even%kill; call new_odd%kill
        end subroutine prepare_external_init_vol

    end subroutine exec_refine3D_auto

    subroutine exec_refine3D_multi( self, cline )
        use simple_abinitio_utils, only: write_final_rec_outputs
        use simple_commanders_rec, only: commander_rec3D
        use simple_estimate_ssnr, only: lpstages_setlims
        class(commander_refine3D_multi), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(cmdline)               :: cline_rec3D
        type(parameters)            :: params, params_final_rec
        type(sp_project)            :: spproj
        type(lp_crop_inf)           :: lpinfo_multi(2)
        type(string), allocatable   :: init_vols(:)
        type(string)                :: multivol_mode
        integer, parameter :: NSAMPLE_PER_STATE_REFINE3D_MULTI = 10000
        integer, parameter :: NSAMPLE_REFINE3D_MULTI_CAP       = 100000
        integer, parameter :: STAGE1_NSPACE                    = 2500
        integer, parameter :: STAGE2_NSPACE                    = 5000
        integer, parameter :: STAGE2_NSPACE_SUB                = 500
        integer, parameter :: INIT_MAXITS_REFINE3D_MULTI = 5
        integer, parameter :: STAGE2_MINITS = 5
        integer, parameter :: MINITS_REFINE3D_MULTI = 10
        integer, parameter :: MAXITS_REFINE3D_MULTI_CAP = 50
        real,    parameter :: TARGET_UPDATES_PER_PARTICLE_REFINE3D_MULTI = 4.0
        real,    parameter :: STATE_OVERLAP_EARLY_REFINE3D_MULTI = 0.95
        real,    parameter :: STATE_OVERLAP_NEIGH_REFINE3D_MULTI = 0.99
        real,    parameter :: LPSTART_REFINE3D_MULTI = 10.0
        real,    parameter :: LPSTOP_REFINE3D_MULTI = 6.0
        character(len=*), parameter :: WORKFLOW_LABEL = 'REFINE3D_MULTI'
        integer :: nstates_project, nptcls_eff, nsample_target, nptcls_per_iter
        integer :: maxits_user, stage_cap, init_niters, stage1_niters, stage2_niters, total_iter
        integer :: maxits_glob_multi, init_stage_cap, min_maxits_required
        real    :: update_frac_auto, state_overlap
        logical :: l_maxits_defined, l_init_state_assignment, l_nstates_on_cline
        logical :: l_has_project_multistates, l_run_init_stage, l_run_shc_stage, l_run_prob_neigh_stage
        ! commanders
        type(commander_rec3D)    :: xrec3D
        type(commander_refine3D) :: xrefine3D
        maxits_user = 0
        init_niters = 0
        stage1_niters = 0
        stage2_niters = 0
        init_stage_cap = INIT_MAXITS_REFINE3D_MULTI
        l_init_state_assignment = .false.
        l_has_project_multistates = .false.
        l_run_init_stage = .false.
        l_run_shc_stage = .false.
        l_run_prob_neigh_stage = .false.
        min_maxits_required = STAGE2_MINITS
        call cline%set('prg', 'refine3D_multi')
        l_nstates_on_cline = cline%defined('nstates')
        if( .not. cline%defined('multivol_mode') ) call cline%set('multivol_mode', 'input_oris_start')
        multivol_mode = cline%get_carg('multivol_mode')
        call validate_refine3D_multi_mode()
        call set_refine3D_multi_nstates()
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
        call cline%set('envfsc',          'no')
        call cline%set('lplim_crit',       0.5)
        call cline%set('incrreslim',      'no')
        call cline%set('nu_refine',       'no')
        ! overridable defaults
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',                  'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',                  'no')
        if( .not. cline%defined('sigma_est')   ) call cline%set('sigma_est',           'global')
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',              'no')
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',              'yes')
        if( .not. cline%defined('nsample')     ) call cline%set('nsample', default_refine3D_multi_nsample())
        if( .not. cline%defined('autoscale')   ) call cline%set('autoscale',              'yes')
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',                 'yes')
        if( .not. cline%defined('filt_mode')   ) call cline%set('filt_mode', 'nonuniform_lpset')
        if( .not. cline%defined('lpstop')      ) call cline%set('lpstop',  LPSTOP_REFINE3D_MULTI)
        if( .not. cline%defined('automsk')     ) call cline%set('automsk',                 'no')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap', STATE_OVERLAP_NEIGH_REFINE3D_MULTI)
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol',                 'no')
        l_maxits_defined = cline%defined('maxits')
        if( l_maxits_defined )then
            if( trim(multivol_mode%to_char()).eq.'input_oris_fixed' ) min_maxits_required = 1
            maxits_user = cline%get_iarg('maxits')
            if( maxits_user < min_maxits_required )then
                THROW_HARD('maxits must be >= '//int2str(min_maxits_required)//' for '//WORKFLOW_LABEL)
            endif
        endif
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call set_refine3D_multi_sampling()
        call prepare_refine3D_multi_class_sampling()
        call configure_refine3D_multi_stages()
        call set_refine3D_multi_downscaling()
        call initialize_state_volumes()
        call cline%set('prg', 'refine3D')
        call cline%set('ufrac_trec', params%update_frac)
        maxits_glob_multi = 0
        if( l_run_init_stage      ) maxits_glob_multi = maxits_glob_multi + init_stage_cap
        if( l_run_shc_stage       ) maxits_glob_multi = maxits_glob_multi + stage_cap
        if( l_run_prob_neigh_stage) maxits_glob_multi = maxits_glob_multi + stage_cap
        call cline%set('maxits_glob', max(1, maxits_glob_multi))
        total_iter = 0
        if( l_run_init_stage )then
            call run_refine3D_multi_stage(0, 'prob_state', STAGE1_NSPACE, 0, 1, init_niters, &
                &init_stage_cap, STATE_OVERLAP_EARLY_REFINE3D_MULTI)
        endif
        if( l_run_shc_stage )then
            call run_refine3D_multi_stage(1, 'shc_smpl', STAGE1_NSPACE, 0, 1, stage1_niters, &
                &stage_cap, STATE_OVERLAP_EARLY_REFINE3D_MULTI)
        endif
        if( l_run_prob_neigh_stage )then
            call run_refine3D_multi_stage(2, 'prob_neigh', STAGE2_NSPACE, STAGE2_NSPACE_SUB, STAGE2_MINITS, &
                &stage2_niters, stage_cap, params%overlap)
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> '//WORKFLOW_LABEL//' STAGE ITERATIONS INIT/STAGE1/STAGE2/TOTAL: ', &
            &init_niters, '/', stage1_niters, '/', stage2_niters, '/', total_iter
        call verify_all_active_particles_updated()
        ! re-reconstruct from all particle images
        cline_rec3D = cline
        call cline_rec3D%set('prg', 'reconstruct3D')
        call cline_rec3D%set('outfile', 'RESOLUTION_FINAL.txt')
        call cline_rec3D%set('postprocess', 'yes')
        call cline_rec3D%delete('trail_rec')
        call cline_rec3D%delete('refine')
        call cline_rec3D%delete('sigma_est')
        call cline_rec3D%delete('update_frac')
        call cline_rec3D%delete('ufrac_trec')
        call cline_rec3D%delete('box_crop')
        call cline_rec3D%delete('smpd_crop')
        call cline_rec3D%set('objfun', 'cc')
        if( params%l_nonuniform ) call cline_rec3D%set('filt_mode', 'none')
        call cline_rec3D%set('nu_refine', 'no')
        call xrec3D%execute(cline_rec3D)
        call params_final_rec%new(cline_rec3D)
        params_final_rec%box  = params_final_rec%box_crop
        params_final_rec%smpd = params_final_rec%smpd_crop
        call spproj%read_segment('out', params_final_rec%projfile)
        call write_final_rec_outputs(params_final_rec, spproj, params_final_rec%lpstop)
        call spproj%kill
        call cleanup_init_vols()
        call multivol_mode%kill
        call simple_end('**** SIMPLE_REFINE3D_MULTI NORMAL STOP ****')

    contains

        integer function default_refine3D_multi_nsample() result(nsample_default)
            nsample_default = min(NSAMPLE_REFINE3D_MULTI_CAP, &
                &NSAMPLE_PER_STATE_REFINE3D_MULTI * nstates_project)
        end function default_refine3D_multi_nsample

        subroutine validate_refine3D_multi_mode()
            select case(trim(multivol_mode%to_char()))
                case('input_oris_start', 'input_oris_refine', 'input_oris_fixed')
                    ! supported
                case default
                    THROW_HARD('Unsupported multivol_mode for '//WORKFLOW_LABEL//': '//trim(multivol_mode%to_char()))
            end select
            write(logfhandle,'(A,A)') '>>> '//WORKFLOW_LABEL//' MULTIVOL_MODE: ', trim(multivol_mode%to_char())
        end subroutine validate_refine3D_multi_mode

        subroutine configure_refine3D_multi_stages()
            select case(trim(multivol_mode%to_char()))
                case('input_oris_fixed')
                    l_run_init_stage       = .true.
                    l_run_shc_stage        = .false.
                    l_run_prob_neigh_stage = .false.
                    init_stage_cap         = stage_cap
                case('input_oris_refine')
                    l_run_init_stage       = l_init_state_assignment
                    l_run_shc_stage        = .false.
                    l_run_prob_neigh_stage = .true.
                case('input_oris_start')
                    l_run_init_stage       = l_init_state_assignment
                    l_run_shc_stage        = .true.
                    l_run_prob_neigh_stage = .true.
            end select
            write(logfhandle,'(A,L1,A,L1,A,L1)') '>>> '//WORKFLOW_LABEL//' STAGES INIT/SHC/PROB_NEIGH: ', &
                &l_run_init_stage, '/', l_run_shc_stage, '/', l_run_prob_neigh_stage
        end subroutine configure_refine3D_multi_stages

        subroutine set_refine3D_multi_nstates()
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
                        THROW_HARD(WORKFLOW_LABEL//' input_oris_fixed expects state=0/1 input')
                    case('input_oris_start')
                        if( l_nstates_on_cline )then
                            THROW_HARD(WORKFLOW_LABEL//' input_oris_start+nstates expects state=0/1 input')
                        endif
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
        end subroutine set_refine3D_multi_nstates

        subroutine set_refine3D_multi_sampling()
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
            if( trim(multivol_mode%to_char()).eq.'input_oris_fixed' )then
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
                maxits_auto = ceiling((TARGET_UPDATES_PER_PARTICLE_REFINE3D_MULTI * real(nptcls_eff)) / real(nptcls_per_iter))
                stage_cap   = max(MINITS_REFINE3D_MULTI, min(MAXITS_REFINE3D_MULTI_CAP, max(STAGE2_MINITS, maxits_auto)))
                params%maxits = stage_cap
                call cline%set('maxits', params%maxits)
                write(logfhandle,'(A,I0,A,F5.1,A)') '>>> '//WORKFLOW_LABEL//' STAGE MAXITS: ', &
                    &stage_cap, ' FOR ~', TARGET_UPDATES_PER_PARTICLE_REFINE3D_MULTI, ' UPDATES/PARTICLE'
            endif
        end subroutine set_refine3D_multi_sampling

        subroutine prepare_refine3D_multi_class_sampling()
            type(sp_project) :: cls_proj
            type(class_sample), allocatable :: clssmp(:)
            integer, allocatable :: clsinds(:)
            if( trim(params%balance) /= 'yes' ) return
            if( .not. params%l_update_frac ) return
            call cls_proj%read(params%projfile)
            clsinds = cls_proj%get_selected_clsinds()
            if( .not. allocated(clsinds) )then
                THROW_HARD(WORKFLOW_LABEL//' class-balanced fractional updates require selected cls2D entries')
            endif
            if( size(clsinds) < 1 )then
                THROW_HARD(WORKFLOW_LABEL//' class-balanced fractional updates require selected cls2D entries')
            endif
            call cls_proj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
            call write_class_samples(clssmp, string(CLASS_SAMPLING_FILE))
            write(logfhandle,'(A,I0,A,F5.2)') &
                &'>>> '//WORKFLOW_LABEL//' CLASS-BALANCED SAMPLING CLASSES/FRAC_BEST: ', &
                &size(clsinds), '/', params%frac_best
            call deallocate_class_samples(clssmp)
            if( allocated(clsinds) ) deallocate(clsinds)
            call cls_proj%kill
        end subroutine prepare_refine3D_multi_class_sampling

        subroutine set_refine3D_multi_downscaling()
            if( .not. params%l_autoscale )then
                call cline%delete('box_crop')
                call cline%delete('smpd_crop')
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' AUTOSCALE: off'
                return
            endif
            call lpstages_setlims(params%box, 2, params%smpd, LPSTART_REFINE3D_MULTI, params%lpstop, lpinfo_multi)
            params%trs = lpinfo_multi(2)%trslim
            call cline%set('trs', params%trs)
            if( lpinfo_multi(2)%l_autoscale )then
                params%box_crop  = lpinfo_multi(2)%box_crop
                params%smpd_crop = lpinfo_multi(2)%smpd_crop
                call cline%set('box_crop', params%box_crop)
                call cline%set('smpd_crop', params%smpd_crop)
                write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
                    &'>>> '//WORKFLOW_LABEL//' AUTOSCALE BOX/SMPD_CROP: ', &
                    &params%box, '/', params%box_crop, '/', params%smpd_crop
            else
                call cline%delete('box_crop')
                call cline%delete('smpd_crop')
                write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' AUTOSCALE: native sampling retained'
            endif
        end subroutine set_refine3D_multi_downscaling

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
                THROW_HARD(WORKFLOW_LABEL//' state=0/1 input requires complete vol1..volN or project state volumes')
            endif
            cline_rec3D = cline
            call cline_rec3D%set('prg', 'reconstruct3D')
            call cline_rec3D%delete('trail_rec')
            call cline_rec3D%delete('refine')
            call cline_rec3D%delete('sigma_est')
            call cline_rec3D%delete('update_frac')
            call cline_rec3D%delete('ufrac_trec')
            call cline_rec3D%delete('box_crop')
            call cline_rec3D%delete('smpd_crop')
            call cline_rec3D%set('objfun', 'cc')
            call cline_rec3D%set('postprocess', 'no')
            call cline_rec3D%set('nu_refine', 'no')
            call xrec3D%execute(cline_rec3D)
            do state = 1,nstates_project
                call cline%set('vol'//int2str(state), refine3D_state_vol_fname(state))
            enddo
            write(logfhandle,'(A)') '>>> '//WORKFLOW_LABEL//' INITIALIZED STATE VOLUMES BY RECONSTRUCTION'
        end subroutine initialize_state_volumes

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
            real    :: vol_smpd
            do state = 1,nstates_project
                vol = cline%get_carg('vol'//int2str(state))
                if( .not. file_exists(vol) ) THROW_HARD('Input volume does not exist: '//vol%to_char())
                call find_ldim_nptcls(vol, ldim, nptcls_dummy)
                vol_smpd = find_img_smpd(vol)
                if( any(ldim /= [params%box,params%box,params%box]) .or. abs(vol_smpd - params%smpd) > 1.e-6 )then
                    THROW_HARD('Input state volumes must have same dimensions/sampling as the project particles')
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
                if( init_box /= params%box .or. init_smpd <= 0. .or. abs(init_smpd - params%smpd) > 1.e-6 )then
                    call spproj%kill
                    return
                endif
            enddo
            l_compatible = .true.
            call spproj%kill
        end function project_state_volumes_compatible

        subroutine run_refine3D_multi_stage( stage, refine_mode, nspace_stage, nspace_sub_stage, min_stage_iters, &
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
        end subroutine run_refine3D_multi_stage

        subroutine verify_all_active_particles_updated()
            type(sp_project) :: update_proj
            integer, allocatable :: states(:), updatecnts(:)
            integer :: nactive, nupdated, nmissing
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
            if( nactive < 1 )then
                if( allocated(states)     ) deallocate(states)
                if( allocated(updatecnts) ) deallocate(updatecnts)
                call update_proj%kill
                THROW_HARD(WORKFLOW_LABEL//' has no active particles after staged refinement')
            endif
            if( nmissing > 0 )then
                if( allocated(states)     ) deallocate(states)
                if( allocated(updatecnts) ) deallocate(updatecnts)
                call update_proj%kill
                THROW_HARD(WORKFLOW_LABEL//' cannot finish before every active particle has been updated')
            endif
            if( allocated(states)     ) deallocate(states)
            if( allocated(updatecnts) ) deallocate(updatecnts)
            call update_proj%kill
        end subroutine verify_all_active_particles_updated

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

    end subroutine exec_refine3D_multi

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

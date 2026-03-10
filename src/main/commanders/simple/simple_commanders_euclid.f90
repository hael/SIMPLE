!@descr: for sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid
use simple_commanders_api
use simple_sigma2_binfile, only: sigma2_binfile
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec_distr
  contains
    procedure :: execute      => exec_calc_pspec_distr
end type commander_calc_pspec_distr

type, extends(commander_base) :: commander_calc_pspec
  contains
    procedure :: execute      => exec_calc_pspec
end type commander_calc_pspec

type, extends(commander_base) :: commander_calc_group_sigmas
  contains
    procedure :: execute      => exec_calc_group_sigmas
end type commander_calc_group_sigmas

type, extends(commander_base) :: estimate_first_sigmas_commander
  contains
    procedure :: execute      => exec_estimate_first_sigmas
end type estimate_first_sigmas_commander

contains

    subroutine exec_calc_pspec_distr( self, cline )
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble
        class(commander_calc_pspec_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! commanders
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        ! command lines
        type(cmdline)        :: cline_calc_pspec
        type(cmdline)        :: cline_calc_pspec_assemble
        ! other variables
        class(oris), pointer :: spproj_field => NULL()
        type(parameters)     :: params
        type(sp_project)     :: spproj
        type(qsys_env)       :: qenv
        type(chash)          :: job_descr
        logical              :: fall_over
        call cline%set('stream','no')
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('projfile') )then
            THROW_HARD('Missing project file entry; exec_calc_pspec_distr')
        endif
        ! init
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl2D','ptcl3D','cls3D')
                fall_over = spproj%get_nptcls() == 0
            case DEFAULT
                write(logfhandle,*)'Unsupported ORITYPE; simple_commanders_euclid :: exec_calc_pspec_distr'
        end select
        call spproj%ptr2oritype(params%oritype, spproj_field)
        if( fall_over )then
            THROW_HARD('no particles found! :exec_calc_pspec_distr')
        endif
        if( spproj_field%get_nevenodd() == 0 )then
            call spproj_field%partition_eo
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! prepare command lines from prototype master
        cline_calc_pspec          = cline
        cline_calc_pspec_assemble = cline
        ! initialise static command line parameters and static job description parameter
        call cline_calc_pspec%set('prg', 'calc_pspec' )                   ! required for distributed call
        call cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble' ) ! required for local call
        ! setup the environment for distributed execution
        call qenv%new(params, params%nparts)
        call cline_calc_pspec%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! assemble
        call xcalc_pspec_assemble%execute(cline_calc_pspec_assemble)
        ! end gracefully
        call spproj%kill
        call cline_calc_pspec%kill
        call cline_calc_pspec_assemble%kill
        call qenv%kill
        call job_descr%kill
        call qsys_cleanup(params)
        call simple_touch(CALCPSPEC_FINISHED)
        call simple_end('**** SIMPLE_DISTR_CALC_PSPEC NORMAL STOP ****')
    end subroutine exec_calc_pspec_distr

    subroutine exec_calc_pspec( self, cline )
        use simple_strategy2D3D_common, only: prepimgbatch, discrete_read_imgbatch, killimgbatch
        class(commander_calc_pspec), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)         :: params
        type(image)              :: sum_img
        type(builder)            :: build
        type(sigma2_binfile)     :: binfile
        type(string)             :: binfname
        complex(dp), allocatable :: cmat_thr_sum(:,:,:)
        complex,     allocatable :: cmat_sum(:,:,:)
        integer,     allocatable :: pinds(:)
        real,        allocatable :: pspec(:), sigma2(:,:)
        integer :: batchlims(2),kfromto(2)
        integer :: i,iptcl,imatch,nyq,nptcls_part_sel,batchsz_max,nbatch
        logical :: l_scale_update_frac
        call cline%set('mkdir', 'no')
        call cline%set('stream','no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! Sampling
        ! Because this is always run prior to reconstruction/search, sampling is not always informed
        ! or may change with workflows. Instead of setting a sampling for the following operations when
        ! l_update_frac, we sample uniformly AND do not write the corresponding field
        l_scale_update_frac = .false.
        if( params%l_update_frac )then
            call build%spproj_field%sample4update_rnd([params%fromp,params%top], params%update_frac, nptcls_part_sel, pinds, .false. )
            l_scale_update_frac = .true.
        else
            call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls_part_sel, pinds, .false.)
        endif
        ! init
        nyq = build%img%get_nyq()
        allocate(sigma2(nyq,params%fromp:params%top),pspec(nyq),source=0.)
        batchsz_max = min(nptcls_part_sel,50 * nthr_glob)
        call prepimgbatch(params, build, batchsz_max)
        call sum_img%new([params%box,params%box,1],params%smpd)
        call sum_img%zero_and_flag_ft
        cmat_sum = sum_img%allocate_cmat()
        allocate(cmat_thr_sum(size(cmat_sum,dim=1),size(cmat_sum,dim=2),1))
        ! mask memoization
        call build%imgbatch(1)%memoize_mask_coords
        do i = 1,nptcls_part_sel,batchsz_max
            batchlims = [i, min(i+batchsz_max-1,nptcls_part_sel)]
            nbatch    = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nbatch, pinds(batchlims(1):batchlims(2)), [1,nbatch])
            cmat_thr_sum = dcmplx(0.d0,0.d0)
            !$omp parallel do default(shared) private(iptcl,imatch,pspec)&
            !$omp schedule(static) proc_bind(close) reduction(+:cmat_thr_sum)
            do imatch = 1,nbatch
                iptcl = pinds(batchlims(1)+imatch-1)
                call build%imgbatch(imatch)%norm_noise_mask_fft_powspec(build%lmsk, params%msk, pspec)
                if( l_scale_update_frac )then
                    ! To account for spectra not included in sampling and yield the correct average
                    sigma2(:,iptcl) = pspec / (2.0 * params%update_frac)
                else
                    sigma2(:,iptcl) = pspec / 2.0
                endif
                ! thread average
                call build%imgbatch(imatch)%add_dble_cmat2mat(cmat_thr_sum(:,:,:))
            end do
            !$omp end parallel do
            ! global average
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmplx(cmat_thr_sum(:,:,:),kind=sp)
        end do
        call sum_img%set_cmat(cmat_sum)
        call sum_img%write(string('sum_img_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char())
        ! write to disk
        kfromto  = [1, nyq]
        binfname = 'init_pspec_part'//trim(int2str(params%part))//'.dat'
        call binfile%new(binfname,params%fromp,params%top,kfromto)
        call binfile%write(sigma2)
        ! destruct
        call build%kill_general_tbox
        call binfile%kill
        call killimgbatch(build)
        call sum_img%kill
        ! end gracefully
        call qsys_job_finished(params, string('simple_commanders_euclid :: exec_calc_pspec'))
        call simple_end('**** SIMPLE_CALC_PSPEC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_pspec

    subroutine exec_calc_group_sigmas( self, cline )
        class(commander_calc_group_sigmas), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        type(sigma2_binfile)  :: binfile
        type(sigma_array)     :: sigma2_array
        type(string)          :: starfile_fname
        real,     allocatable :: pspecs(:,:)
        real(dp), allocatable :: group_weights(:,:), group_pspecs(:,:,:)
        real(dp)              :: w
        integer               :: kfromto(2),iptcl,ipart,eo,ngroups,igroup,fromp,top
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! read sigmas from binfiles
        do ipart = 1,params%nparts
            sigma2_array%fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
            call binfile%new_from_file(sigma2_array%fname)
            call binfile%read(sigma2_array%sigma2)
            fromp = lbound(sigma2_array%sigma2,2)
            top   = ubound(sigma2_array%sigma2,2)
            if( (fromp<1).or.(top>params%nptcls) )then
                THROW_HARD('commander_euclid; exec_calc_group_sigmas; file ' // sigma2_array%fname%to_char() // ' has ptcl range ' // int2str(fromp) // '-' // int2str(top))
            end if
            if( ipart == 1 )then
                call binfile%get_resrange(kfromto)
                allocate(pspecs(kfromto(1):kfromto(2),params%nptcls))
            endif
            pspecs(:,fromp:top) = sigma2_array%sigma2(:,:)
            deallocate(sigma2_array%sigma2)
        end do
        call binfile%kill
        if( params%l_sigma_glob )then
            ngroups = 1
            allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.d0)
            !$omp parallel do default(shared) private(iptcl,eo,w)&
            !$omp schedule(static) proc_bind(close) reduction(+:group_pspecs,group_weights)
            do iptcl = 1,params%nptcls
                if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                eo = build%spproj_field%get_eo(iptcl) ! 0/1
                w  = real(build%spproj_field%get(iptcl,'w'),dp)
                if( w < TINY )cycle
                group_pspecs(eo+1,1,:) = group_pspecs (eo+1,1,:) + w * real(pspecs(:,iptcl),dp)
                group_weights(eo+1,1)  = group_weights(eo+1,1)   + w
            enddo
            !$omp end parallel do
        else
            ngroups = 0
            !$omp parallel do default(shared) private(iptcl,igroup)&
            !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
            do iptcl = 1,params%nptcls
                if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                igroup  = nint(build%spproj_field%get(iptcl,'stkind'))
                ngroups = max(igroup,ngroups)
            enddo
            !$omp end parallel do
            allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.d0)
            do iptcl = 1,params%nptcls
                if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                eo     = build%spproj_field%get_eo(iptcl) ! 0/1
                igroup = nint(build%spproj_field%get(iptcl,'stkind'))
                w      = real(build%spproj_field%get(iptcl,'w'),dp)
                if( w < TINY )cycle
                group_pspecs(eo+1,igroup,:) = group_pspecs (eo+1,igroup,:) + w * real(pspecs(:,iptcl),dp)
                group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)   + w
            enddo
        endif
        deallocate(pspecs)
        do eo = 1,2
            do igroup = 1,ngroups
                if( group_weights(eo,igroup) < TINY ) cycle
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / group_weights(eo,igroup)
            end do
        end do
        ! write group sigmas to starfile
        starfile_fname = SIGMA2_GROUP_FBODY//int2str(params%which_iter)//STAR_EXT
        call write_groups_starfile(starfile_fname, real(group_pspecs), ngroups)
        ! cleanup
        call build%kill_general_tbox
        call simple_touch('CALC_GROUP_SIGMAS_FINISHED')
        call simple_end('**** SIMPLE_CALC_GROUP_SIGMAS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_group_sigmas

    subroutine exec_estimate_first_sigmas( self, cline )
        use simple_strategy3D_matcher, only: refine3D_exec
        class(estimate_first_sigmas_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        ! command lines
        type(cmdline)    :: cline_first_sigmas, cline_calc_group_sigmas
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        logical          :: l_shmem, converged
        if( .not. cline%defined('pgrp')     ) THROW_HARD('point-group symmetry (pgrp) is needed for first sigma estimation')
        if( .not. cline%defined('mskdiam')  ) THROW_HARD('mask diameter (mskdiam) is needed for first sigma estimation')
        if( .not. cline%defined('nthr')     ) THROW_HARD('number of threads (nthr) is needed for first sigma estimation')
        if( .not. cline%defined('projfile') ) THROW_HARD('missing project file entry; exec_estimate_first_sigmas')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not.cline%defined('vol1') )then
            if( cline%defined('polar') )then
                if( cline%get_carg('polar').eq.'yes' )then
                    if( .not.file_exists(POLAR_REFS_FBODY//BIN_EXT) )then
                        THROW_HARD('starting polar references are needed for first sigma estimation')
                    endif
                else
                    THROW_HARD('starting volume is needed for first sigma estimation')
                endif
            else
                THROW_HARD('starting volume is needed for first sigma estimation')
            endif
        endif
        l_shmem = .not.cline%defined('nparts')
        cline_first_sigmas = cline
        call cline_first_sigmas%set('prg', 'refine3D')
        call cline_first_sigmas%set('center',    'no')
        call cline_first_sigmas%set('continue',  'no')
        call cline_first_sigmas%set('maxits',       1)
        call cline_first_sigmas%set('which_iter',   1)
        call cline_first_sigmas%set('objfun','euclid')
        call cline_first_sigmas%set('refine', 'sigma')
        call cline_first_sigmas%delete('update_frac') ! all particles neeed to contribute
        call cline_first_sigmas%delete('hp')
        call cline_first_sigmas%set('mkdir', 'no')    ! generate the sigma files in the root refine3D dir
        cline_calc_group_sigmas = cline_first_sigmas
        if( cline%defined('startit') )then
            call cline_calc_group_sigmas%set('which_iter', params%startit)
        else
            call cline_calc_group_sigmas%set('which_iter', params%which_iter)
        endif
        ! init
        if( l_shmem )then
            call build%init_params_and_build_strategy3D_tbox(cline_first_sigmas, params )
            call refine3D_exec(params, build, cline_first_sigmas, params%which_iter, converged)
            call build%kill_strategy3D_tbox
        else
            call build%init_params_and_build_spproj(cline_first_sigmas, params)
            ! setup the environment for distributed execution
            call qenv%new(params, params%nparts)
            ! prepare job description
            call cline_first_sigmas%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
            ! assemble
            call xcalc_group_sigmas%execute(cline_calc_group_sigmas)
            ! end gracefully
            call qsys_cleanup(params)
        endif
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call cline_first_sigmas%kill
        call cline_calc_group_sigmas%kill
        call simple_end('**** SIMPLE_ESTIMATE_FIRST_SIGMAS NORMAL STOP ****')
    end subroutine exec_estimate_first_sigmas

end module simple_commanders_euclid

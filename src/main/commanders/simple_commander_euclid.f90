module simple_commander_euclid
include 'simple_lib.f08'
use simple_builder,          only: builder, build_glob
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_parameters,       only: parameters, params_glob
use simple_sigma2_binfile,   only: sigma2_binfile
use simple_qsys_env,         only: qsys_env
use simple_euclid_sigma2,    only: write_groups_starfile
use simple_image,            only: image
use simple_qsys_funs
implicit none

public :: calc_pspec_commander_distr
public :: calc_pspec_commander
public :: calc_pspec_assemble_commander
public :: calc_group_sigmas_commander
public :: estimate_first_sigmas_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: calc_pspec_commander_distr
  contains
    procedure :: execute      => exec_calc_pspec_distr
end type calc_pspec_commander_distr

type, extends(commander_base) :: calc_pspec_commander
  contains
    procedure :: execute      => exec_calc_pspec
end type calc_pspec_commander

type, extends(commander_base) :: calc_pspec_assemble_commander
  contains
    procedure :: execute      => exec_calc_pspec_assemble
end type calc_pspec_assemble_commander

type, extends(commander_base) :: calc_group_sigmas_commander
  contains
    procedure :: execute      => exec_calc_group_sigmas
end type calc_group_sigmas_commander

type, extends(commander_base) :: estimate_first_sigmas_commander
  contains
    procedure :: execute      => exec_estimate_first_sigmas
end type estimate_first_sigmas_commander

type :: sigma_array
    character(len=:), allocatable :: fname
    real,             allocatable :: sigma2(:,:)
end type sigma_array

character(len=STDLEN), parameter :: PSPEC_FBODY = 'pspec_'

contains

    subroutine exec_calc_pspec_distr( self, cline )
        use simple_sp_project, only: sp_project
        class(calc_pspec_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
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
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('projfile') )then
            THROW_HARD('Missing project file entry; exec_calc_pspec_distr')
        endif
        ! init
        call params%new(cline)
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo')
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl2D','ptcl3D')
                fall_over = spproj%get_nptcls() == 0
            case DEFAULT
                write(logfhandle,*)'Unsupported ORITYPE; simple_commander_euclid :: exec_calc_pspec_distr'
        end select
        call spproj%ptr2oritype(params%oritype, spproj_field)
        if( fall_over )then
            THROW_HARD('no particles found! :exec_refine3D_distr')
        endif
        if( spproj_field%get_nevenodd() == 0 )then
            THROW_HARD('no even/odd flag found! :calc_pspec_distr')
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
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
        ! assemble
        call qenv%exec_simple_prg_in_queue(cline_calc_pspec_assemble, 'CALC_PSPEC_FINISHED')
        ! end gracefully
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_CALC_PSPEC NORMAL STOP ****')
    end subroutine exec_calc_pspec_distr

    subroutine exec_calc_pspec( self, cline )
        use simple_parameters,          only: params_glob
        use simple_strategy2D3D_common, only: prepimgbatch, read_imgbatch
        class(calc_pspec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: sum_img
        type(builder)        :: build
        complex, pointer     :: cmat(:,:,:), cmat_sum(:,:,:)
        real,    allocatable :: pspecs(:,:), pspec(:)
        logical, allocatable :: mask(:)
        real                 :: sdev_noise
        integer              :: batchlims(2),iptcl,iptcl_batch,imatch,nyq,nptcls_part,batchsz_max
        real, allocatable    :: sigma2(:,:)
        type(sigma2_binfile) :: binfile
        integer              :: kfromto(2)
        character(len=:), allocatable :: binfname
        call cline%set('mkdir', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! init
        nptcls_part = params%top-params%fromp+1
        nyq         = build%img%get_nyq()
        batchsz_max = 10 * nthr_glob
        allocate(mask(batchsz_max),source=.false.)
        allocate(pspecs(nyq,nptcls_part),source=0.)
        call prepimgbatch(batchsz_max)
        call sum_img%new([params%box,params%box,1],params%smpd)
        call sum_img%zero_and_flag_ft
        call sum_img%get_cmat_ptr(cmat_sum)
        do iptcl_batch=params_glob%fromp,params_glob%top,batchsz_max
            batchlims = [iptcl_batch, min(params_glob%top,iptcl_batch + batchsz_max - 1)]
            ! mask
            do iptcl = batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                mask(imatch) = .not. (build%spproj_field%get_state(iptcl) == 0)
            enddo
            ! read
            call read_imgbatch( batchlims )
            ! preprocess
            !$omp parallel do default(shared) private(iptcl,imatch,pspec)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                if( .not. mask(imatch) ) cycle
                ! normalize
                call build%imgbatch(imatch)%noise_norm(build%lmsk, sdev_noise)
                !  mask
                if( params%l_focusmsk )then
                    call build%imgbatch(imatch)%mask(params%focusmsk, 'softavg')
                else
                    call build%imgbatch(imatch)%mask(params%msk, 'softavg')
                endif
                call build%imgbatch(imatch)%fft
                ! power spectrum
                call build%imgbatch(imatch)%spectrum('power',pspec,norm=.true.)
                pspecs(:,iptcl-params_glob%fromp+1) = pspec / 2.0
            end do
            !$omp end parallel do
            ! global average
            do iptcl=batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                if( .not. mask(imatch) ) cycle
                call build%imgbatch(imatch)%get_cmat_ptr(cmat)
                !$omp workshare
                cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmat(:,:,:)
                !$omp end workshare
            enddo
        end do
        call sum_img%write('sum_img_part'//int2str_pad(params%part,params%numlen)//params%ext)
        ! write to disk
        kfromto(1) = 1
        kfromto(2) = nyq
        binfname = 'init_pspec_part'//trim(int2str(params%part))//'.dat'
        allocate(sigma2(nyq,params%fromp:params%top))
        do iptcl = params%fromp, params%top
            sigma2(:,iptcl) = pspecs(:,iptcl-params%fromp+1)
        end do
        ! taking account scaling of images
        if( cline%defined('scale') .and. abs(params%scale-1.0) > 0.01 )then
            !$omp workshare
            sigma2 = sigma2 * params%scale
            !$omp end workshare
        endif
        call binfile%new(binfname,params%fromp,params%top,kfromto)
        call binfile%write(sigma2)
        ! end gracefully
        call qsys_job_finished('simple_commander_euclid :: exec_calc_pspec')
        call simple_end('**** SIMPLE_CALC_PSPEC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_pspec

    subroutine exec_calc_pspec_assemble( self, cline )
        class(calc_pspec_assemble_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)                 :: params
        type(image)                      :: avg_img
        type(builder)                    :: build
        type(sigma2_binfile)             :: binfile
        type(sigma_array), allocatable   :: sigma2_arrays(:)
        character(len=:),  allocatable   :: part_fname,starfile_fname,outbin_fname
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup,nyq,pspec_l,pspec_u
        real,              allocatable   :: group_pspecs(:,:,:),pspec_ave(:),pspecs(:,:),sigma2_output(:,:)
        integer,           allocatable   :: group_weights(:,:)
        call cline%set('mkdir', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! set Fourier index range
        params%kfromto(1) = 1
        params%kfromto(2) = calc_fourier_index(2.*params%smpd, params%box, params%smpd)
        ! generate average power spectrum
        nptcls     = build%spproj_field%get_noris(consider_state=.false.)
        nptcls_sel = build%spproj_field%get_noris(consider_state=.true.)
        call avg_img%new([params%box,params%box,1], params%smpd)
        call avg_img%zero_and_flag_ft
        do ipart = 1,params%nparts
            call build%img%zero_and_flag_ft
            part_fname = 'sum_img_part'//int2str_pad(ipart,params%numlen)//params%ext
            call build%img%read(part_fname)
            call avg_img%add(build%img)
            call del_file(part_fname)
        enddo
        call avg_img%div(real(nptcls_sel))
        ! calculate power spectrum
        call avg_img%spectrum('power',pspec_ave,norm=.true.)
        pspec_ave = pspec_ave / 2.0
        nyq = avg_img%get_nyq()
        ! read power spectra of particles
        allocate(pspecs(nyq,params%nptcls),sigma2_arrays(params%nparts))
        do ipart = 1,params%nparts
            sigma2_arrays(ipart)%fname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call binfile%new_from_file(sigma2_arrays(ipart)%fname)
            call binfile%read(sigma2_arrays(ipart)%sigma2)
            pspec_l = lbound(sigma2_arrays(ipart)%sigma2,2)
            pspec_u = ubound(sigma2_arrays(ipart)%sigma2,2)
            if( (pspec_l<1).or.(pspec_u>params%nptcls) )then
                THROW_HARD('commander_euclid; exec_calc_pspec_assemble; file ' // sigma2_arrays(ipart)%fname // ' has ptcl range ' // int2str(pspec_l) // '-' // int2str(pspec_u))
            end if
            pspecs(:,pspec_l:pspec_u) = sigma2_arrays(ipart)%sigma2(:,:)
        end do
        ! generate group averages & write
        if( params_glob%l_sigma_glob )then
            ngroups = 1
        else
            ngroups = 0
            !$omp parallel do default(shared) private(iptcl,igroup)&
            !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
            do iptcl = 1,nptcls
                if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                igroup  = nint(build%spproj_field%get(iptcl,'stkind'))
                ngroups = max(igroup,ngroups)
            enddo
            !$omp end parallel do
        endif
        allocate(group_pspecs(2,ngroups,nyq),source=0.)
        allocate(group_weights(2,ngroups),source=0)
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build%spproj_field%get(iptcl,'eo')) ! 0/1
            if( params_glob%l_sigma_glob )then
                igroup = 1
            else
                igroup = nint(build%spproj_field%get(iptcl,'stkind'))
            endif
            group_pspecs(eo+1,igroup,:) = group_pspecs(eo+1,igroup,:) + pspecs(:, iptcl)
            group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)  + 1
        enddo
        do eo = 1,2
            do igroup = 1,ngroups
                if( group_weights(eo,igroup) < 1 ) cycle
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / real(group_weights(eo,igroup))
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) - pspec_ave(:)
                call remove_negative_sigmas(eo, igroup)
            end do
        end do
        ! write group sigmas to starfile
        if( cline%defined('which_iter') )then
            starfile_fname = trim(SIGMA2_GROUP_FBODY)//int2str(params%which_iter)//'.star'
        else
            starfile_fname = trim(SIGMA2_GROUP_FBODY)//'1.star'
        endif
        call write_groups_starfile(starfile_fname, group_pspecs, ngroups)
        ! update sigmas in binfiles to match averages
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build%spproj_field%get(iptcl,'eo')) ! 0/1
            if( params_glob%l_sigma_glob )then
                igroup = 1
            else
                igroup = nint(build%spproj_field%get(iptcl,'stkind'))
            endif
            pspecs(:,iptcl) = group_pspecs(eo+1,igroup,:)
        enddo
        ! write updated sigmas to disc
        do ipart = 1,params%nparts
            pspec_l = lbound(sigma2_arrays(ipart)%sigma2,2)
            pspec_u = ubound(sigma2_arrays(ipart)%sigma2,2)
            if( allocated(sigma2_output) ) deallocate(sigma2_output)
            allocate(sigma2_output(params%kfromto(1):params%kfromto(2),pspec_l:pspec_u))
            do iptcl = pspec_l, pspec_u
                sigma2_output(params%kfromto(1):params%kfromto(2),iptcl) = pspecs(params%kfromto(1):params%kfromto(2),iptcl)
            end do
            outbin_fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
            call binfile%new(outbin_fname, fromp=pspec_l, top=pspec_u, kfromto=(/params%kfromto(1), params%kfromto(2)/))
            call binfile%write(sigma2_output)
        end do
        ! end gracefully
        do ipart = 1,params%nparts
            deallocate(sigma2_arrays(ipart)%fname)
            deallocate(sigma2_arrays(ipart)%sigma2)
        end do
        call simple_touch('CALC_PSPEC_FINISHED',errmsg='In: commander_euclid::calc_pspec_assemble')
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine remove_negative_sigmas(eo, igroup)
            integer, intent(in) :: eo, igroup
            logical :: is_positive
            logical :: fixed_from_prev
            integer :: nn, idx
            ! remove any negative sigma2 noise values: replace by positive neighboring value
            do idx = 1, size(group_pspecs, 3)
                if( group_pspecs(eo,igroup,idx) < 0. )then
                    ! first try the previous value
                    fixed_from_prev = .false.
                    if( idx - 1 >= 1 )then
                        if( group_pspecs(eo,igroup,idx-1) > 0. )then
                            group_pspecs(eo,igroup,idx) = group_pspecs(eo,igroup,idx-1)
                            fixed_from_prev = .true.
                        end if
                    end if
                    if( .not. fixed_from_prev )then
                        is_positive = .false.
                        nn          = idx
                        do while (.not. is_positive)
                            nn = nn + 1
                            if( nn > size(group_pspecs,3) )then
                                THROW_HARD('BUG! Cannot find positive values in sigma2 noise spectrum; eo=' // trim(int2str(eo)) // ', igroup=' // trim(int2str(igroup)))
                            end if
                            if( group_pspecs(eo,igroup,nn) > 0. )then
                                is_positive = .true.
                                group_pspecs(eo,igroup,idx) = group_pspecs(eo,igroup,nn)
                            end if
                        end do
                    end if
                end if
            end do
        end subroutine remove_negative_sigmas

    end subroutine exec_calc_pspec_assemble

    subroutine exec_calc_group_sigmas( self, cline )
        class(calc_group_sigmas_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(sigma2_binfile)          :: binfile
        type(sigma_array)             :: sigma2_array
        character(len=:), allocatable :: starfile_fname
        real,             allocatable :: group_pspecs(:,:,:),glob_pspec(:,:,:),pspecs(:,:)
        real,             allocatable :: group_weights(:,:)
        real                          :: w
        integer                       :: kfromto(2),iptcl,ipart,eo,ngroups,igroup,fromp,top
        if( associated(build_glob) )then
            if( .not.associated(params_glob) )then
                THROW_HARD('Builder & parameters must be associated for shared memory execution!')
            endif
        else
            call cline%set('mkdir', 'no')
            if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
            call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        endif
        ! read sigmas from binfiles
        do ipart = 1,params_glob%nparts
            sigma2_array%fname = SIGMA2_FBODY//int2str_pad(ipart,params_glob%numlen)//'.dat'
            call binfile%new_from_file(sigma2_array%fname)
            call binfile%read(sigma2_array%sigma2)
            fromp = lbound(sigma2_array%sigma2,2)
            top   = ubound(sigma2_array%sigma2,2)
            if( (fromp<1).or.(top>params_glob%nptcls) )then
                THROW_HARD('commander_euclid; exec_calc_group_sigmas; file ' // sigma2_array%fname // ' has ptcl range ' // int2str(fromp) // '-' // int2str(top))
            end if
            if( ipart == 1 )then
                call binfile%get_resrange(kfromto)
                allocate(pspecs(kfromto(1):kfromto(2),params_glob%nptcls))
            endif
            pspecs(:,fromp:top) = sigma2_array%sigma2(:,:)
            deallocate(sigma2_array%sigma2)
        end do
        ngroups = 0
        !$omp parallel do default(shared) private(iptcl,igroup)&
        !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
        do iptcl = 1,params_glob%nptcls
            if( build_glob%spproj_field%get_state(iptcl) == 0 ) cycle
            igroup  = nint(build_glob%spproj_field%get(iptcl,'stkind'))
            ngroups = max(igroup,ngroups)
        enddo
        !$omp end parallel do
        allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)),glob_pspec(2,1,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.)
        do iptcl = 1,params_glob%nptcls
            if( build_glob%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build_glob%spproj_field%get(iptcl,'eo'    )) ! 0/1
            igroup = nint(build_glob%spproj_field%get(iptcl,'stkind'))
            w      = build_glob%spproj_field%get(iptcl,'w')
            if( w < TINY )cycle
            group_pspecs(eo+1,igroup,:) = group_pspecs (eo+1,igroup,:) + w * pspecs(:,iptcl)
            group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)   + w
        enddo
        if( params_glob%l_sigma_glob )then
            glob_pspec = 0.
            w          = 1./real(ngroups)
            do eo = 1,2
                do igroup = 1,ngroups
                    if( group_weights(eo,igroup) < TINY ) cycle
                    group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / group_weights(eo,igroup)
                    glob_pspec(eo,1,:)        = glob_pspec(eo,1,:) + w * group_pspecs(eo,igroup,:)
                end do
            end do
            ! write global sigma to starfile
            starfile_fname = trim(SIGMA2_GROUP_FBODY)//trim(int2str(params_glob%which_iter))//'.star'
            call write_groups_starfile(starfile_fname, glob_pspec, 1)
        else
            do eo = 1,2
                do igroup = 1,ngroups
                    if( group_weights(eo,igroup) < TINY ) cycle
                    group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / group_weights(eo,igroup)
                end do
            end do
            ! write group sigmas to starfile
            starfile_fname = trim(SIGMA2_GROUP_FBODY)//trim(int2str(params_glob%which_iter))//'.star'
            call write_groups_starfile(starfile_fname, group_pspecs, ngroups)
        endif
        call simple_touch('CALC_GROUP_SIGMAS_FINISHED',errmsg='In: commander_euclid::calc_group_sigmas')
        call simple_end('**** SIMPLE_CALC_GROUP_SIGMAS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_group_sigmas

    subroutine exec_estimate_first_sigmas( self, cline )
        class(estimate_first_sigmas_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        ! command lines
        type(cmdline)    :: cline_first_sigmas
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( .not. cline%defined('vol1')     ) THROW_HARD('starting volume is needed for first sigma estimation')
        if( .not. cline%defined('pgrp')     ) THROW_HARD('point-group symmetry (pgrp) is needed for first sigma estimation')
        if( .not. cline%defined('mskdiam')  ) THROW_HARD('mask diameter (mskdiam) is needed for first sigma estimation')
        if( .not. cline%defined('nthr')     ) THROW_HARD('number of threads (nthr) is needed for first sigma estimation')
        if( .not. cline%defined('nparts')   ) THROW_HARD('number of partitions (npart) is needed for first sigma estimation (distributed workflow)')
        if( .not. cline%defined('projfile') ) THROW_HARD('missing project file entry; exec_estimate_first_sigmas')
        cline_first_sigmas = cline
        call cline_first_sigmas%set('prg', 'refine3D')
        call cline_first_sigmas%set('center',    'no')
        call cline_first_sigmas%set('continue',  'no')
        call cline_first_sigmas%set('ptclw',     'no')
        call cline_first_sigmas%delete('lp_iters')
        call cline_first_sigmas%set('maxits',     1.0)
        call cline_first_sigmas%set('which_iter', 1.0)
        call cline_first_sigmas%set('objfun','euclid')
        call cline_first_sigmas%set('refine', 'sigma')
        call cline_first_sigmas%delete('update_frac') ! all particles neeed to contribute
        call cline_first_sigmas%delete('hp')
        call cline_first_sigmas%delete('lp')
        call cline_first_sigmas%delete('lpstop')
        call cline_first_sigmas%set('oritype', 'ptcl3D')
        call cline_first_sigmas%set('mkdir', 'no')    ! generate the sigma files in the root refine3D dir
        ! init
        call build%init_params_and_build_spproj(cline_first_sigmas, params)
        call build%spproj%update_projinfo(cline_first_sigmas)
        call build%spproj%write_segment_inside('projinfo')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline_first_sigmas%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! end gracefully
        call qsys_cleanup
        call build%spproj%kill
        call simple_end('**** SIMPLE_ESTIMATE_FIRST_SIGMAS NORMAL STOP ****')
    end subroutine exec_estimate_first_sigmas

end module simple_commander_euclid

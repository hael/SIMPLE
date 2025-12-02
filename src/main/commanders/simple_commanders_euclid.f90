module simple_commanders_euclid
include 'simple_lib.f08'
use simple_builder,        only: builder, build_glob
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_euclid_sigma2,  only: write_groups_starfile
use simple_image,          only: image
use simple_parameters,     only: parameters, params_glob
use simple_qsys_env,       only: qsys_env
use simple_sigma2_binfile, only: sigma2_binfile
use simple_qsys_funs
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

type, extends(commander_base) :: commander_calc_pspec_assemble
  contains
    procedure :: execute      => exec_calc_pspec_assemble
end type commander_calc_pspec_assemble

type, extends(commander_base) :: commander_calc_group_sigmas
  contains
    procedure :: execute      => exec_calc_group_sigmas
end type commander_calc_group_sigmas

type :: sigma_array
    type(string)      :: fname
    real, allocatable :: sigma2(:,:)
end type sigma_array

contains

    subroutine exec_calc_pspec_distr( self, cline )
        use simple_sp_project, only: sp_project
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
            call spproj%write_segment_inside(params_glob%oritype, params%projfile)
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
        call cline_calc_pspec%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! assemble
        call xcalc_pspec_assemble%execute_safe(cline_calc_pspec_assemble)
        ! end gracefully
        call spproj%kill
        call cline_calc_pspec%kill
        call cline_calc_pspec_assemble%kill
        call qenv%kill
        call job_descr%kill
        call qsys_cleanup
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
        complex,     pointer     :: cmat(:,:,:), cmat_sum(:,:,:)
        integer,     allocatable :: pinds(:)
        real,        allocatable :: pspec(:), sigma2(:,:)
        complex(dp), allocatable :: cmat_thr_sum(:,:,:)
        real    :: sdev_noise
        integer :: batchlims(2),kfromto(2)
        integer :: i,iptcl,imatch,nyq,nptcls_part_sel,batchsz_max,nbatch
        logical :: l_scale_update_frac
        call cline%set('mkdir', 'no')
        call cline%set('stream','no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! Sampling
        ! Because this is always run prior to reconstruction/search, sampling is not always informed
        ! or may change with workflows. Instead of setting a sampling for the following oprations when
        ! l_update_frac, we sample uniformly AND do not write the corresponding field
        l_scale_update_frac = .false.
        if( params%l_update_frac )then
            call build%spproj_field%sample4update_rnd([params%fromp,params%top], params_glob%update_frac, nptcls_part_sel, pinds, .false. )
            l_scale_update_frac = .true.
        else
            call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls_part_sel, pinds, .false.)
        endif
        ! init
        nyq = build%img%get_nyq()
        allocate(sigma2(nyq,params%fromp:params%top),pspec(nyq),source=0.)
        batchsz_max = min(nptcls_part_sel,50 * nthr_glob)
        call prepimgbatch(batchsz_max)
        call sum_img%new([params%box,params%box,1],params%smpd)
        call sum_img%zero_and_flag_ft
        call sum_img%get_cmat_ptr(cmat_sum)
        allocate(cmat_thr_sum(size(cmat_sum,dim=1),size(cmat_sum,dim=2),1))
        do i = 1,nptcls_part_sel,batchsz_max
            batchlims = [i, min(i+batchsz_max-1,nptcls_part_sel)]
            nbatch    = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(nbatch, pinds(batchlims(1):batchlims(2)), [1,nbatch])
            cmat_thr_sum = dcmplx(0.d0,0.d0)
            !$omp parallel do default(shared) private(iptcl,imatch,pspec,cmat)&
            !$omp schedule(static) proc_bind(close) reduction(+:cmat_thr_sum)
            do imatch = 1,nbatch
                iptcl = pinds(batchlims(1)+imatch-1)
                ! normalize
                if( params%l_noise_norm )then
                    call build%imgbatch(imatch)%norm_noise(build%lmsk, sdev_noise)
                else
                    call build%imgbatch(imatch)%norm_within(build%lmsk)
                endif
                !  mask
                if( params%l_focusmsk )then
                    call build%imgbatch(imatch)%mask(params%focusmsk, 'softavg')
                else
                    call build%imgbatch(imatch)%mask(params%msk, 'softavg')
                endif
                ! power spectrum
                call build%imgbatch(imatch)%fft
                call build%imgbatch(imatch)%power_spectrum(pspec)
                if( l_scale_update_frac )then
                    ! To account for spectra not included in sampling and yield the correct average
                    sigma2(:,iptcl) = pspec / (2.0 * params_glob%update_frac)
                else
                    sigma2(:,iptcl) = pspec / 2.0
                endif
                ! thread average
                call build%imgbatch(imatch)%get_cmat_ptr(cmat)
                cmat_thr_sum(:,:,:) = cmat_thr_sum(:,:,:) + cmplx(cmat(:,:,:),kind=dp)
            end do
            !$omp end parallel do
            ! global average
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmplx(cmat_thr_sum(:,:,:),kind=sp)
        end do
        call sum_img%write(string('sum_img_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char())
        ! write to disk
        kfromto  = [1, nyq]
        binfname = 'init_pspec_part'//trim(int2str(params%part))//'.dat'
        call binfile%new(binfname,params%fromp,params%top,kfromto)
        call binfile%write(sigma2)
        ! destruct
        call build%kill_general_tbox
        call binfile%kill
        call killimgbatch
        call sum_img%kill
        ! end gracefully
        call qsys_job_finished(string('simple_commanders_euclid :: exec_calc_pspec'))
        call simple_end('**** SIMPLE_CALC_PSPEC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_pspec

    subroutine exec_calc_pspec_assemble( self, cline )
        class(commander_calc_pspec_assemble), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)                 :: params
        type(image)                      :: avg_img
        type(builder)                    :: build
        type(sigma2_binfile)             :: binfile
        type(sigma_array), allocatable   :: sigma2_arrays(:)
        type(string)                     :: part_fname,starfile_fname,outbin_fname
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup,nyq,pspec_l,pspec_u
        real(dp),          allocatable   :: group_pspecs(:,:,:)
        real,              allocatable   :: pspec_ave(:),pspecs(:,:),sigma2_output(:,:)
        integer,           allocatable   :: group_weights(:,:)
        call cline%set('mkdir', 'no')
        call cline%set('stream','no')
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
            part_fname = 'sum_img_part'//int2str_pad(ipart,params%numlen)//params%ext%to_char()
            call build%img%read(part_fname)
            call avg_img%add(build%img)
            call del_file(part_fname)
        enddo
        call avg_img%div(real(nptcls_sel))
        ! calculate power spectrum
        call avg_img%spectrum('power',pspec_ave,norm=.true.)
        pspec_ave = pspec_ave / 2.0
        nyq = avg_img%get_nyq()
        call avg_img%kill
        ! read power spectra of particles
        allocate(pspecs(nyq,params%nptcls),sigma2_arrays(params%nparts))
        do ipart = 1,params%nparts
            sigma2_arrays(ipart)%fname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call binfile%new_from_file(sigma2_arrays(ipart)%fname)
            call binfile%read(sigma2_arrays(ipart)%sigma2)
            pspec_l = lbound(sigma2_arrays(ipart)%sigma2,2)
            pspec_u = ubound(sigma2_arrays(ipart)%sigma2,2)
            if( (pspec_l<1).or.(pspec_u>params%nptcls) )then
                THROW_HARD('commander_euclid; exec_calc_pspec_assemble; file ' // sigma2_arrays(ipart)%fname%to_char()// ' has ptcl range ' // int2str(pspec_l) // '-' // int2str(pspec_u))
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
                igroup  = build%spproj_field%get_int(iptcl,'stkind')
                ngroups = max(igroup,ngroups)
            enddo
            !$omp end parallel do
        endif
        allocate(group_pspecs(2,ngroups,nyq),source=0.d0)
        allocate(group_weights(2,ngroups),source=0)
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo = build%spproj_field%get_eo(iptcl) ! 0/1
            if( params_glob%l_sigma_glob )then
                igroup = 1
            else
                igroup = build%spproj_field%get_int(iptcl, 'stkind')
            endif
            group_pspecs(eo+1,igroup,:) = group_pspecs(eo+1,igroup,:) + real(pspecs(:, iptcl),dp)
            group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)  + 1
        enddo
        do eo = 1,2
            do igroup = 1,ngroups
                if( group_weights(eo,igroup) < 1 ) cycle
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / real(group_weights(eo,igroup),dp)
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) - real(pspec_ave(:),dp)
                call remove_negative_sigmas(eo, igroup)
            end do
        end do
        ! write group sigmas to starfile
        if( cline%defined('which_iter') )then
            starfile_fname = SIGMA2_GROUP_FBODY//int2str(params%which_iter)//STAR_EXT
        else
            starfile_fname = SIGMA2_GROUP_FBODY//'1'//STAR_EXT
        endif
        call write_groups_starfile(starfile_fname, real(group_pspecs), ngroups)
        ! update sigmas in binfiles to match averages
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build%spproj_field%get(iptcl,'eo')) ! 0/1
            if( params_glob%l_sigma_glob )then
                igroup = 1
            else
                igroup = nint(build%spproj_field%get(iptcl,'stkind'))
            endif
            pspecs(:,iptcl) = real(group_pspecs(eo+1,igroup,:))
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
            call binfile%new(outbin_fname, fromp=pspec_l, top=pspec_u, kfromto=[params%kfromto(1), params%kfromto(2)])
            call binfile%write(sigma2_output)
        end do
        ! end gracefully
        do ipart = 1,params%nparts
            call sigma2_arrays(ipart)%fname%kill
            deallocate(sigma2_arrays(ipart)%sigma2)
        end do
        deallocate(sigma2_arrays,group_pspecs,pspec_ave,pspecs,group_weights)
        call binfile%kill
        call build%kill_general_tbox
        call simple_touch('CALC_PSPEC_FINISHED')
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine remove_negative_sigmas(eo, igroup)
            integer, intent(in) :: eo, igroup
            logical :: is_positive
            logical :: fixed_from_prev
            integer :: nn, idx
            ! remove any negative sigma2 noise values: replace by positive neighboring value
            do idx = 1, size(group_pspecs, 3)
                if( group_pspecs(eo,igroup,idx) < 0.d0 )then
                    ! first try the previous value
                    fixed_from_prev = .false.
                    if( idx - 1 >= 1 )then
                        if( group_pspecs(eo,igroup,idx-1) > 0.d0 )then
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
                                THROW_HARD('BUG! Cannot find positive values in sigma2 noise spectrum; eo=' // int2str(eo) // ', igroup=' // int2str(igroup))
                            end if
                            if( group_pspecs(eo,igroup,nn) > 0.d0 )then
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
        if( associated(build_glob) )then
            if( .not.associated(params_glob) )then
                THROW_HARD('Builder & parameters must be associated for shared memory execution!')
            endif
        else
            call cline%set('mkdir',  'no')
            call cline%set('stream', 'no')
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
                THROW_HARD('commander_euclid; exec_calc_group_sigmas; file ' // sigma2_array%fname%to_char() // ' has ptcl range ' // int2str(fromp) // '-' // int2str(top))
            end if
            if( ipart == 1 )then
                call binfile%get_resrange(kfromto)
                allocate(pspecs(kfromto(1):kfromto(2),params_glob%nptcls))
            endif
            pspecs(:,fromp:top) = sigma2_array%sigma2(:,:)
            deallocate(sigma2_array%sigma2)
        end do
        call binfile%kill
        if( params_glob%l_sigma_glob )then
            ngroups = 1
            allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.d0)
            !$omp parallel do default(shared) private(iptcl,eo,w)&
            !$omp schedule(static) proc_bind(close) reduction(+:group_pspecs,group_weights)
            do iptcl = 1,params_glob%nptcls
                if( build_glob%spproj_field%get_state(iptcl) == 0 ) cycle
                eo = build_glob%spproj_field%get_eo(iptcl) ! 0/1
                w  = real(build_glob%spproj_field%get(iptcl,'w'),dp)
                if( w < TINY )cycle
                group_pspecs(eo+1,1,:) = group_pspecs (eo+1,1,:) + w * real(pspecs(:,iptcl),dp)
                group_weights(eo+1,1)  = group_weights(eo+1,1)   + w
            enddo
            !$omp end parallel do
        else
            ngroups = 0
            !$omp parallel do default(shared) private(iptcl,igroup)&
            !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
            do iptcl = 1,params_glob%nptcls
                if( build_glob%spproj_field%get_state(iptcl) == 0 ) cycle
                igroup  = nint(build_glob%spproj_field%get(iptcl,'stkind'))
                ngroups = max(igroup,ngroups)
            enddo
            !$omp end parallel do
            allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.d0)
            do iptcl = 1,params_glob%nptcls
                if( build_glob%spproj_field%get_state(iptcl) == 0 ) cycle
                eo     = build_glob%spproj_field%get_eo(iptcl) ! 0/1
                igroup = nint(build_glob%spproj_field%get(iptcl,'stkind'))
                w      = real(build_glob%spproj_field%get(iptcl,'w'),dp)
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
        starfile_fname = SIGMA2_GROUP_FBODY//int2str(params_glob%which_iter)//STAR_EXT
        call write_groups_starfile(starfile_fname, real(group_pspecs), ngroups)
        ! cleanup
        call build%kill_general_tbox
        call simple_touch('CALC_GROUP_SIGMAS_FINISHED')
        call simple_end('**** SIMPLE_CALC_GROUP_SIGMAS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_group_sigmas

end module simple_commanders_euclid

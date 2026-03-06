!@descr: for distributed sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid_distr
use simple_commanders_api
use simple_sigma2_binfile, only: sigma2_binfile
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec_distr_worker
  contains
    procedure :: execute      => exec_calc_pspec_distr_worker
end type commander_calc_pspec_distr_worker

type, extends(commander_base) :: commander_calc_pspec_assemble
  contains
    procedure :: execute      => exec_calc_pspec_assemble
end type commander_calc_pspec_assemble

contains

    !> Worker subroutine for distributed execution.
    !> This is the task executed by each worker node.
    subroutine exec_calc_pspec_distr_worker( self, cline )
        use simple_calc_pspec_common, only: calc_pspec_exec
        class(commander_calc_pspec_distr_worker), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        ! Initialize
        call cline%set('stream','no')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! The worker builds its own environment based on the command line arguments
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! Execute the core calculation
        call calc_pspec_exec(params, build)
        ! Cleanup
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_pspec :: exec_calc_pspec_distr_worker'))
        call simple_end('**** SIMPLE_CALC_PSPEC WORKER NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_pspec_distr_worker

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
        if( params%l_sigma_glob )then
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
            if( params%l_sigma_glob )then
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
            if( params%l_sigma_glob )then
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

end module simple_commanders_euclid_distr

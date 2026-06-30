!@descr: for distributed sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid_distr
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_commanders_api
use simple_sigma2_binfile, only: sigma2_binfile
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec_assemble
  contains
    procedure :: execute      => exec_calc_pspec_assemble
end type commander_calc_pspec_assemble

contains

    subroutine exec_calc_pspec_assemble( self, cline )
        class(commander_calc_pspec_assemble), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)                 :: params
        type(image)                      :: avg_img
        type(builder)                    :: build
        type(sigma2_binfile)             :: binfile
        type(string)                     :: part_fname,starfile_fname,outbin_fname
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup,nstks,nyq,pspec_l,pspec_u
        real(dp),          allocatable   :: group_pspecs(:,:,:)
        real,              allocatable   :: pspec_ave(:),sigma2_part(:,:),sigma2_output(:,:)
        integer,           allocatable   :: group_weights(:,:)
        logical,           allocatable   :: pspec_covered(:)
        call cline%set('mkdir', 'no')
        call cline%set('stream','no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! set Fourier index range
        params%kfromto(1) = 1
        params%kfromto(2) = calc_fourier_index(2.*params%smpd, params%box, params%smpd)
        ! generate average power spectrum for global sigma bootstrapping
        nptcls     = build%spproj_field%get_noris(consider_state=.false.)
        nptcls_sel = build%spproj_field%count_state_gt_zero()
        if( nptcls_sel < 1 ) THROW_HARD('No active particles found; exec_calc_pspec_assemble')
        nyq = build%img%get_nyq()
        if( params%l_sigma_glob )then
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
        else
            do ipart = 1,params%nparts
                part_fname = 'sum_img_part'//int2str_pad(ipart,params%numlen)//params%ext%to_char()
                call del_file(part_fname)
            enddo
        endif
        ! generate group averages & write
        if( params%l_sigma_glob )then
            ngroups = 1
        else
            nstks   = build%spproj%get_nstks()
            ngroups = 0
            do iptcl = 1,nptcls
                if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                igroup  = build%spproj_field%get_int(iptcl,'stkind')
                if( igroup < 1 .or. igroup > nstks )then
                    write(logfhandle,*) 'iptcl/stkind/nstks: ', iptcl, igroup, nstks
                    THROW_HARD('commander_euclid; exec_calc_pspec_assemble; particle stkind out of range')
                endif
                ngroups = max(igroup,ngroups)
            enddo
        endif
        allocate(group_pspecs(2,ngroups,nyq),source=0.d0)
        allocate(group_weights(2,ngroups),source=0)
        allocate(pspec_covered(params%nptcls), source=.false.)
        do ipart = 1,params%nparts
            part_fname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call binfile%new_from_file(part_fname)
            call binfile%read(sigma2_part)
            pspec_l = lbound(sigma2_part,2)
            pspec_u = ubound(sigma2_part,2)
            if( (pspec_l<1).or.(pspec_u>params%nptcls) )then
                write(logfhandle,*) 'file/ptcl range/nptcls: ', part_fname%to_char(), pspec_l, pspec_u, params%nptcls
                THROW_HARD('commander_euclid; exec_calc_pspec_assemble; invalid particle spectra range')
            end if
            if( any(pspec_covered(pspec_l:pspec_u)) )then
                THROW_HARD('commander_euclid; exec_calc_pspec_assemble; overlapping particle spectra ranges')
            endif
            pspec_covered(pspec_l:pspec_u) = .true.
            do iptcl = pspec_l,pspec_u
                if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                eo = build%spproj_field%get_eo(iptcl) ! 0/1
                if( params%l_sigma_glob )then
                    igroup = 1
                else
                    igroup = build%spproj_field%get_int(iptcl, 'stkind')
                endif
                if( (.not.all(ieee_is_finite(sigma2_part(:,iptcl)))) .or. &
                    ((.not.params%l_sigma_glob) .and. (.not.any(sigma2_part(:,iptcl) > real(DTINY)))) )then
                    write(logfhandle,*) 'iptcl/eo/igroup: ', iptcl, eo, igroup
                    write(logfhandle,*) 'finite/positive: ', &
                        all(ieee_is_finite(sigma2_part(:,iptcl))), any(sigma2_part(:,iptcl) > real(DTINY))
                    THROW_HARD('active particle sigma spectrum was not computed; exec_calc_pspec_assemble')
                endif
                group_pspecs(eo+1,igroup,:) = group_pspecs(eo+1,igroup,:) + real(sigma2_part(:, iptcl),dp)
                group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)  + 1
            end do
            deallocate(sigma2_part)
        end do
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
            if( .not.pspec_covered(iptcl) )then
                write(logfhandle,*) 'iptcl: ', iptcl
                THROW_HARD('active particle sigma spectrum was not covered; exec_calc_pspec_assemble')
            endif
        enddo
        do eo = 1,2
            do igroup = 1,ngroups
                if( group_weights(eo,igroup) < 1 ) cycle
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / real(group_weights(eo,igroup),dp)
                if( params%l_sigma_glob ) group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) - real(pspec_ave(:),dp)
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
        ! write updated sigmas to disc, one partition at a time
        do ipart = 1,params%nparts
            part_fname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call binfile%new_from_file(part_fname)
            call binfile%read(sigma2_part)
            pspec_l = lbound(sigma2_part,2)
            pspec_u = ubound(sigma2_part,2)
            if( allocated(sigma2_output) ) deallocate(sigma2_output)
            allocate(sigma2_output(params%kfromto(1):params%kfromto(2),pspec_l:pspec_u))
            sigma2_output = sigma2_part(params%kfromto(1):params%kfromto(2),pspec_l:pspec_u)
            do iptcl = pspec_l, pspec_u
                if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                eo = build%spproj_field%get_eo(iptcl) ! 0/1
                if( params%l_sigma_glob )then
                    igroup = 1
                else
                    igroup = build%spproj_field%get_int(iptcl, 'stkind')
                endif
                sigma2_output(params%kfromto(1):params%kfromto(2),iptcl) =&
                    &real(group_pspecs(eo+1,igroup,params%kfromto(1):params%kfromto(2)))
            end do
            outbin_fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
            call binfile%new(outbin_fname, fromp=pspec_l, top=pspec_u, kfromto=[params%kfromto(1), params%kfromto(2)])
            call binfile%write(sigma2_output)
            deallocate(sigma2_part)
        end do
        if( trim(params%match_src) == 'den' )then
            call assemble_pspec_channel('init_pspec_match_part', 'sum_img_match_part', &
                SIGMA2_MATCH_FBODY, SIGMA2_MATCH_GROUP_FBODY, 'match ')
        endif
        ! end gracefully
        deallocate(group_pspecs,group_weights,pspec_covered)
        if( allocated(sigma2_output) ) deallocate(sigma2_output)
        if( allocated(pspec_ave) ) deallocate(pspec_ave)
        call binfile%kill
        call build%kill_general_tbox
        call simple_touch('CALC_PSPEC_FINISHED')
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine assemble_pspec_channel(init_fbody, sum_fbody, sigma_fbody, group_fbody, label)
            character(len=*), intent(in) :: init_fbody, sum_fbody, sigma_fbody, group_fbody, label
            type(image)                    :: avg_img_l
            type(sigma2_binfile)           :: binfile_l
            type(string)                   :: part_fname_l, starfile_fname_l, outbin_fname_l
            integer                        :: iptcl_l, ipart_l, eo_l, ngroups_l, igroup_l, nstks_l
            integer                        :: pspec_l_l, pspec_u_l
            real(dp), allocatable          :: group_pspecs_l(:,:,:)
            real,     allocatable          :: pspec_ave_l(:), sigma2_part_l(:,:), sigma2_output_l(:,:)
            integer,  allocatable          :: group_weights_l(:,:)
            logical,  allocatable          :: pspec_covered_l(:)
            if( params%l_sigma_glob )then
                call avg_img_l%new([params%box,params%box,1], params%smpd)
                call avg_img_l%zero_and_flag_ft
                do ipart_l = 1,params%nparts
                    call build%img%zero_and_flag_ft
                    part_fname_l = trim(sum_fbody)//int2str_pad(ipart_l,params%numlen)//params%ext%to_char()
                    call build%img%read(part_fname_l)
                    call avg_img_l%add(build%img)
                    call del_file(part_fname_l)
                enddo
                call avg_img_l%div(real(nptcls_sel))
                call avg_img_l%spectrum('power',pspec_ave_l,norm=.true.)
                pspec_ave_l = pspec_ave_l / 2.0
                call avg_img_l%kill
            else
                do ipart_l = 1,params%nparts
                    part_fname_l = trim(sum_fbody)//int2str_pad(ipart_l,params%numlen)//params%ext%to_char()
                    call del_file(part_fname_l)
                enddo
            endif
            if( params%l_sigma_glob )then
                ngroups_l = 1
            else
                nstks_l   = build%spproj%get_nstks()
                ngroups_l = 0
                do iptcl_l = 1,nptcls
                    if( build%spproj_field%get_state(iptcl_l) <= 0 ) cycle
                    igroup_l  = build%spproj_field%get_int(iptcl_l,'stkind')
                    if( igroup_l < 1 .or. igroup_l > nstks_l )then
                        write(logfhandle,*) 'iptcl/stkind/nstks: ', iptcl_l, igroup_l, nstks_l
                        THROW_HARD('commander_euclid; exec_calc_pspec_assemble; match particle stkind out of range')
                    endif
                    ngroups_l = max(igroup_l,ngroups_l)
                enddo
            endif
            allocate(group_pspecs_l(2,ngroups_l,nyq),source=0.d0)
            allocate(group_weights_l(2,ngroups_l),source=0)
            allocate(pspec_covered_l(params%nptcls), source=.false.)
            do ipart_l = 1,params%nparts
                part_fname_l = trim(init_fbody)//trim(int2str(ipart_l))//'.dat'
                call binfile_l%new_from_file(part_fname_l)
                call binfile_l%read(sigma2_part_l)
                pspec_l_l = lbound(sigma2_part_l,2)
                pspec_u_l = ubound(sigma2_part_l,2)
                if( (pspec_l_l<1).or.(pspec_u_l>params%nptcls) )then
                    write(logfhandle,*) 'file/ptcl range/nptcls: ', part_fname_l%to_char(), pspec_l_l, pspec_u_l, params%nptcls
                    THROW_HARD('commander_euclid; exec_calc_pspec_assemble; invalid match particle spectra range')
                end if
                if( any(pspec_covered_l(pspec_l_l:pspec_u_l)) )then
                    THROW_HARD('commander_euclid; exec_calc_pspec_assemble; overlapping match particle spectra ranges')
                endif
                pspec_covered_l(pspec_l_l:pspec_u_l) = .true.
                do iptcl_l = pspec_l_l,pspec_u_l
                    if( build%spproj_field%get_state(iptcl_l) <= 0 ) cycle
                    eo_l = build%spproj_field%get_eo(iptcl_l)
                    if( params%l_sigma_glob )then
                        igroup_l = 1
                    else
                        igroup_l = build%spproj_field%get_int(iptcl_l, 'stkind')
                    endif
                    if( (.not.all(ieee_is_finite(sigma2_part_l(:,iptcl_l)))) .or. &
                        ((.not.params%l_sigma_glob) .and. (.not.any(sigma2_part_l(:,iptcl_l) > real(DTINY)))) )then
                        write(logfhandle,*) 'iptcl/eo/igroup: ', iptcl_l, eo_l, igroup_l
                        write(logfhandle,*) 'finite/positive: ', &
                            all(ieee_is_finite(sigma2_part_l(:,iptcl_l))), any(sigma2_part_l(:,iptcl_l) > real(DTINY))
                        THROW_HARD('active match particle sigma spectrum was not computed; exec_calc_pspec_assemble')
                    endif
                    group_pspecs_l(eo_l+1,igroup_l,:) = group_pspecs_l(eo_l+1,igroup_l,:) + real(sigma2_part_l(:, iptcl_l),dp)
                    group_weights_l(eo_l+1,igroup_l)  = group_weights_l(eo_l+1,igroup_l)  + 1
                end do
                deallocate(sigma2_part_l)
            end do
            do iptcl_l = 1,nptcls
                if( build%spproj_field%get_state(iptcl_l) <= 0 ) cycle
                if( .not.pspec_covered_l(iptcl_l) )then
                    write(logfhandle,*) 'iptcl: ', iptcl_l
                    THROW_HARD('active match particle sigma spectrum was not covered; exec_calc_pspec_assemble')
                endif
            enddo
            do eo_l = 1,2
                do igroup_l = 1,ngroups_l
                    if( group_weights_l(eo_l,igroup_l) < 1 ) cycle
                    group_pspecs_l(eo_l,igroup_l,:) = group_pspecs_l(eo_l,igroup_l,:) / real(group_weights_l(eo_l,igroup_l),dp)
                    if( params%l_sigma_glob )then
                        group_pspecs_l(eo_l,igroup_l,:) = group_pspecs_l(eo_l,igroup_l,:) - real(pspec_ave_l(:),dp)
                    endif
                    call remove_negative_sigmas_channel(group_pspecs_l, group_weights_l, eo_l, igroup_l, label)
                end do
            end do
            if( cline%defined('which_iter') )then
                starfile_fname_l = trim(group_fbody)//int2str(params%which_iter)//STAR_EXT
            else
                starfile_fname_l = trim(group_fbody)//'1'//STAR_EXT
            endif
            call write_groups_starfile(starfile_fname_l, real(group_pspecs_l), ngroups_l)
            do ipart_l = 1,params%nparts
                part_fname_l = trim(init_fbody)//trim(int2str(ipart_l))//'.dat'
                call binfile_l%new_from_file(part_fname_l)
                call binfile_l%read(sigma2_part_l)
                pspec_l_l = lbound(sigma2_part_l,2)
                pspec_u_l = ubound(sigma2_part_l,2)
                if( allocated(sigma2_output_l) ) deallocate(sigma2_output_l)
                allocate(sigma2_output_l(params%kfromto(1):params%kfromto(2),pspec_l_l:pspec_u_l))
                sigma2_output_l = sigma2_part_l(params%kfromto(1):params%kfromto(2),pspec_l_l:pspec_u_l)
                do iptcl_l = pspec_l_l, pspec_u_l
                    if( build%spproj_field%get_state(iptcl_l) <= 0 ) cycle
                    eo_l = build%spproj_field%get_eo(iptcl_l)
                    if( params%l_sigma_glob )then
                        igroup_l = 1
                    else
                        igroup_l = build%spproj_field%get_int(iptcl_l, 'stkind')
                    endif
                    sigma2_output_l(params%kfromto(1):params%kfromto(2),iptcl_l) =&
                        &real(group_pspecs_l(eo_l+1,igroup_l,params%kfromto(1):params%kfromto(2)))
                end do
                outbin_fname_l = trim(sigma_fbody)//int2str_pad(ipart_l,params%numlen)//'.dat'
                call binfile_l%new(outbin_fname_l, fromp=pspec_l_l, top=pspec_u_l, kfromto=[params%kfromto(1), params%kfromto(2)])
                call binfile_l%write(sigma2_output_l)
                deallocate(sigma2_part_l)
            end do
            deallocate(group_pspecs_l,group_weights_l,pspec_covered_l)
            if( allocated(sigma2_output_l) ) deallocate(sigma2_output_l)
            if( allocated(pspec_ave_l) ) deallocate(pspec_ave_l)
            call binfile_l%kill
        end subroutine assemble_pspec_channel

        subroutine remove_negative_sigmas_channel(group_pspecs_in, group_weights_in, eo, igroup, label)
            real(dp),         intent(inout) :: group_pspecs_in(:,:,:)
            integer,          intent(in)    :: group_weights_in(:,:)
            integer,          intent(in)    :: eo, igroup
            character(len=*), intent(in)    :: label
            logical :: is_positive
            logical :: fixed_from_prev
            integer :: nn, idx
            where( .not. ieee_is_finite(group_pspecs_in(eo,igroup,:)) ) group_pspecs_in(eo,igroup,:) = 1.d0
            if( .not. any(group_pspecs_in(eo,igroup,:) > DTINY) )then
                if( params%l_sigma_glob )then
                    write(logfhandle,*) '>>> WARNING: calc_pspec_assemble floored empty '//trim(label)//&
                        'global sigma2 spectrum to 1.0; eo/weight: ', eo, group_weights_in(eo,igroup)
                    group_pspecs_in(eo,igroup,:) = 1.d0
                    return
                else
                    write(logfhandle,*) 'eo/igroup/weight: ', eo, igroup, group_weights_in(eo,igroup)
                    THROW_HARD('BUG! Cannot find positive values in '//trim(label)//'group sigma2 noise spectrum')
                endif
            endif
            do idx = 1, size(group_pspecs_in, 3)
                if( group_pspecs_in(eo,igroup,idx) <= DTINY )then
                    fixed_from_prev = .false.
                    if( idx - 1 >= 1 )then
                        if( group_pspecs_in(eo,igroup,idx-1) > DTINY )then
                            group_pspecs_in(eo,igroup,idx) = group_pspecs_in(eo,igroup,idx-1)
                            fixed_from_prev = .true.
                        end if
                    end if
                    if( .not. fixed_from_prev )then
                        is_positive = .false.
                        nn          = idx
                        do while (.not. is_positive)
                            nn = nn + 1
                            if( nn > size(group_pspecs_in,3) )then
                                group_pspecs_in(eo,igroup,idx) = 1.d0
                                exit
                            end if
                            if( group_pspecs_in(eo,igroup,nn) > DTINY )then
                                is_positive = .true.
                                group_pspecs_in(eo,igroup,idx) = group_pspecs_in(eo,igroup,nn)
                            end if
                        end do
                    end if
                end if
            end do
        end subroutine remove_negative_sigmas_channel

        subroutine remove_negative_sigmas(eo, igroup)
            integer, intent(in) :: eo, igroup
            logical :: is_positive
            logical :: fixed_from_prev
            integer :: nn, idx
            where( .not. ieee_is_finite(group_pspecs(eo,igroup,:)) ) group_pspecs(eo,igroup,:) = 1.d0
            if( .not. any(group_pspecs(eo,igroup,:) > DTINY) )then
                if( params%l_sigma_glob )then
                    write(logfhandle,*) '>>> WARNING: calc_pspec_assemble floored empty global sigma2 spectrum to 1.0; eo/weight: ', &
                        eo, group_weights(eo,igroup)
                    group_pspecs(eo,igroup,:) = 1.d0
                    return
                else
                    write(logfhandle,*) 'eo/igroup/weight: ', eo, igroup, group_weights(eo,igroup)
                    THROW_HARD('BUG! Cannot find positive values in group sigma2 noise spectrum')
                endif
            endif
            ! remove any negative sigma2 noise values: replace by positive neighboring value
            do idx = 1, size(group_pspecs, 3)
                if( group_pspecs(eo,igroup,idx) <= DTINY )then
                    ! first try the previous value
                    fixed_from_prev = .false.
                    if( idx - 1 >= 1 )then
                        if( group_pspecs(eo,igroup,idx-1) > DTINY )then
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
                                group_pspecs(eo,igroup,idx) = 1.d0
                                exit
                            end if
                            if( group_pspecs(eo,igroup,nn) > DTINY )then
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

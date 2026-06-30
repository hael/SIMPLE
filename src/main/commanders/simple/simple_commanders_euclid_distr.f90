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
        ! end gracefully
        deallocate(group_pspecs,group_weights,pspec_covered)
        if( allocated(sigma2_output) ) deallocate(sigma2_output)
        if( allocated(pspec_ave) ) deallocate(pspec_ave)
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
            where( .not. ieee_is_finite(group_pspecs(eo,igroup,:)) ) group_pspecs(eo,igroup,:) = 1.d0
            if( .not. any(group_pspecs(eo,igroup,:) > DTINY) )then
                if( params%l_sigma_glob )then
                    write(logfhandle,*) '>>> WARNING: calc_pspec_assemble floored empty global sigma2 spectrum to 1.0; ', &
                        'eo/weight: ', eo, group_weights(eo,igroup)
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

!@descr: for sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid
use simple_commanders_api
use simple_sigma2_binfile, only: sigma2_binfile
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec
  contains
    procedure :: execute      => exec_calc_pspec
end type commander_calc_pspec

type, extends(commander_base) :: commander_calc_group_sigmas
  contains
    procedure :: execute      => exec_calc_group_sigmas
end type commander_calc_group_sigmas

contains

    subroutine exec_calc_pspec( self,cline )
        use simple_calc_pspec_strategy, only: calc_pspec_strategy, create_calc_pspec_strategy
        class(commander_calc_pspec), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        class(calc_pspec_strategy), allocatable :: strategy
        call cline%set('stream','no')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype','ptcl3D')
        ! Create and run strategy (shared-memory/worker vs distributed master)
        strategy = create_calc_pspec_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        ! Unified termination point
        call simple_end(strategy%end_msg, print_simple=strategy%end_print_simple)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_calc_pspec

    subroutine exec_calc_group_sigmas( self, cline )
        class(commander_calc_group_sigmas), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        call consolidate_channel(SIGMA2_FBODY, SIGMA2_GROUP_FBODY)
        if( trim(params%match_src) == 'den' )then
            call consolidate_channel(SIGMA2_MATCH_FBODY, SIGMA2_MATCH_GROUP_FBODY)
        endif
        ! cleanup
        call build%kill_general_tbox
        call simple_touch('CALC_GROUP_SIGMAS_FINISHED')
        call simple_end('**** SIMPLE_CALC_GROUP_SIGMAS NORMAL STOP ****', print_simple=.false.)
    contains

        subroutine consolidate_channel(part_fbody, group_fbody)
            character(len=*), intent(in) :: part_fbody, group_fbody
            type(sigma2_binfile)  :: binfile
            type(sigma_array)     :: sigma2_array
            type(string)          :: starfile_fname
            real,     allocatable :: pspecs(:,:)
            real(dp), allocatable :: group_weights(:,:), group_pspecs(:,:,:)
            integer               :: kfromto(2),iptcl,ipart,eo,ngroups,igroup,fromp,top
            do ipart = 1,params%nparts
                sigma2_array%fname = trim(part_fbody)//int2str_pad(ipart,params%numlen)//'.dat'
                call binfile%new_from_file(sigma2_array%fname)
                call binfile%read(sigma2_array%sigma2)
                fromp = lbound(sigma2_array%sigma2,2)
                top   = ubound(sigma2_array%sigma2,2)
                if( (fromp<1).or.(top>params%nptcls) )then
                    write(logfhandle,*) 'sigma2 file/ptcl range/nptcls: ', sigma2_array%fname%to_char(), fromp, top, params%nptcls
                    THROW_HARD('commander_euclid; exec_calc_group_sigmas; invalid sigma2 particle range')
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
                !$omp parallel do default(shared) private(iptcl,eo)&
                !$omp schedule(static) proc_bind(close) reduction(+:group_pspecs,group_weights)
                do iptcl = 1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                    eo = build%spproj_field%get_eo(iptcl) ! 0/1
                    group_pspecs(eo+1,1,:) = group_pspecs (eo+1,1,:) + real(pspecs(:,iptcl),dp)
                    group_weights(eo+1,1)  = group_weights(eo+1,1)   + 1.d0
                enddo
                !$omp end parallel do
            else
                ngroups = 0
                !$omp parallel do default(shared) private(iptcl,igroup)&
                !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
                do iptcl = 1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                    igroup  = nint(build%spproj_field%get(iptcl,'stkind'))
                    ngroups = max(igroup,ngroups)
                enddo
                !$omp end parallel do
                allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.d0)
                do iptcl = 1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                    eo     = build%spproj_field%get_eo(iptcl) ! 0/1
                    igroup = nint(build%spproj_field%get(iptcl,'stkind'))
                    group_pspecs(eo+1,igroup,:) = group_pspecs (eo+1,igroup,:) + real(pspecs(:,iptcl),dp)
                    group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)   + 1.d0
                enddo
            endif
            deallocate(pspecs)
            do eo = 1,2
                do igroup = 1,ngroups
                    if( group_weights(eo,igroup) < TINY ) cycle
                    group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / group_weights(eo,igroup)
                end do
            end do
            starfile_fname = trim(group_fbody)//int2str(params%which_iter)//STAR_EXT
            call write_groups_starfile(starfile_fname, real(group_pspecs), ngroups)
            deallocate(group_weights, group_pspecs)
        end subroutine consolidate_channel
    end subroutine exec_calc_group_sigmas

end module simple_commanders_euclid

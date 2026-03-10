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

type, extends(commander_base) :: estimate_first_sigmas_commander
  contains
    procedure :: execute      => exec_estimate_first_sigmas
end type estimate_first_sigmas_commander

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
        use simple_commanders_refine3D, only: commander_refine3D
        class(estimate_first_sigmas_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_refine3D)          :: xrefine3D
        ! command lines
        type(cmdline) :: cline_first_sigmas, cline_calc_group_sigmas
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
            call cline_calc_group_sigmas%set('which_iter', cline%get_iarg('startit'))
        else if( cline%defined('which_iter') )then
            call cline_calc_group_sigmas%set('which_iter', cline%get_iarg('which_iter'))
        else
            call cline_calc_group_sigmas%set('which_iter', 1)
        endif
        ! Execute refine3D using unified commander (handles both shared/distributed internally)
        call xrefine3D%execute(cline_first_sigmas)
        ! Execute group sigma calculation
        call xcalc_group_sigmas%execute(cline_calc_group_sigmas)
        ! Cleanup
        call cline_first_sigmas%kill
        call cline_calc_group_sigmas%kill
        call simple_end('**** SIMPLE_ESTIMATE_FIRST_SIGMAS NORMAL STOP ****')
    end subroutine exec_estimate_first_sigmas

end module simple_commanders_euclid

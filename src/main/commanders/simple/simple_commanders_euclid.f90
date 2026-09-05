!@descr: for sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_commanders_api
use simple_sigma2_binfile, only: sigma2_binfile
use simple_sigma2_state, only: sigma2_state_candidate_path, sigma2_state_commit, &
    &sigma2_state_merge_local_ranges, sigma2_state_project_layout_digest, &
    &sigma2_state_range_path, sigma2_state_reduce_groups, sigma2_state_validate_identity
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_create_candidate, &
    &sigma2_state_init_header, sigma2_state_read_groups, sigma2_state_read_header, &
    &sigma2_state_validate_file, sigma2_state_write_particles, SIGMA2_GROUP_GLOBAL, &
    &SIGMA2_GROUP_STACK, SIGMA2_PROV_LEGACY_PARTS, SIGMA2_PROV_STAR_SEED, &
    &SIGMA2_STATE_COMMITTED
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

type, extends(commander_base) :: commander_sigma2_convert
  contains
    procedure :: execute => exec_sigma2_convert
end type commander_sigma2_convert

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
        if( params%l_sigma_canonical )then
            call consolidate_canonical
        else
            call consolidate_channel(SIGMA2_FBODY, SIGMA2_GROUP_FBODY)
        endif
        ! cleanup
        call build%kill_general_tbox
        call simple_touch('CALC_GROUP_SIGMAS_FINISHED')
        call simple_end('**** SIMPLE_CALC_GROUP_SIGMAS NORMAL STOP ****', print_simple=.false.)
    contains

        subroutine consolidate_canonical
            type(string) :: state_path, candidate_path
            type(string), allocatable :: range_paths(:)
            logical, allocatable :: scheduled(:), active(:)
            integer, allocatable :: eo_ids(:), group_ids(:)
            integer :: ipart, iptcl, status
            logical :: found
            character(len=STDLEN) :: message
            call build%spproj%get_sigma2_state_path(state_path, found)
            if( .not. found ) THROW_HARD('particle project has no canonical sigma2 state path')
            candidate_path = sigma2_state_candidate_path(state_path%to_char())
            allocate(range_paths(params%nparts))
            do ipart = 1, params%nparts
                range_paths(ipart) = sigma2_state_range_path(state_path%to_char(), ipart, params%numlen)
            enddo
            allocate(scheduled(params%nptcls), source=.true.)
            call sigma2_state_merge_local_ranges(candidate_path%to_char(), range_paths, scheduled, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            allocate(active(params%nptcls), eo_ids(params%nptcls), group_ids(params%nptcls))
            do iptcl = 1, params%nptcls
                active(iptcl)    = build%spproj_field%get_state(iptcl) > 0
                eo_ids(iptcl)    = build%spproj_field%get_eo(iptcl)
                group_ids(iptcl) = build%spproj_field%get_int(iptcl, 'stkind')
            enddo
            call sigma2_state_reduce_groups(candidate_path%to_char(), active, eo_ids, group_ids, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call sigma2_state_commit(candidate_path%to_char(), state_path%to_char(), active, eo_ids, &
                &group_ids, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            do ipart = 1, params%nparts
                call del_file(range_paths(ipart))
            enddo
            call range_paths(:)%kill
            deallocate(range_paths, scheduled, active, eo_ids, group_ids)
            call state_path%kill
            call candidate_path%kill
        end subroutine consolidate_canonical

        subroutine consolidate_channel(part_fbody, group_fbody)
            character(len=*), intent(in) :: part_fbody, group_fbody
            type(sigma2_binfile)  :: binfile
            type(sigma_array)     :: sigma2_array
            type(string)          :: starfile_fname
            real,     allocatable :: pspecs(:,:)
            real(dp), allocatable :: group_weights(:,:), group_pspecs(:,:,:)
            integer               :: kfromto(2),iptcl,ipart,eo,ngroups,igroup,fromp,top
            logical               :: l_updated_only
            l_updated_only = trim(params%cc_emit_sigma) == 'yes'
            if( l_updated_only .and. (.not. build%spproj_field%isthere('updatecnt')) )then
                THROW_HARD('CC residual sigma consolidation requires update counts')
            endif
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
                    if( l_updated_only )then
                        if( build%spproj_field%get(iptcl,'updatecnt') <= 0. ) cycle
                    endif
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
                    if( l_updated_only )then
                        if( build%spproj_field%get(iptcl,'updatecnt') <= 0. ) cycle
                    endif
                    igroup  = nint(build%spproj_field%get(iptcl,'stkind'))
                    ngroups = max(igroup,ngroups)
                enddo
                !$omp end parallel do
                allocate(group_pspecs(2,ngroups,kfromto(1):kfromto(2)), group_weights(2,ngroups),source=0.d0)
                do iptcl = 1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                    if( l_updated_only )then
                        if( build%spproj_field%get(iptcl,'updatecnt') <= 0. ) cycle
                    endif
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

    subroutine exec_sigma2_convert( self, cline )
        use simple_euclid_sigma2, only: read_sigma2_groups_file, write_groups_starfile
        class(commander_sigma2_convert), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: project
        class(oris), pointer :: particles
        type(string) :: state_path, candidate_path, part_path
        type(sigma2_state_header) :: header
        type(sigma2_binfile) :: part_file
        real, allocatable :: star_groups(:,:,:), part_values(:,:)
        real(real32), allocatable :: spectra(:,:), stored_groups(:,:,:)
        logical, allocatable :: active(:), covered(:)
        integer, allocatable :: eo_ids(:), group_ids(:)
        integer(int64) :: layout_digest
        integer :: box, kfromto(2), ngroups, grouping, nptcls, project_ngroups
        integer :: ipart, iptcl, fromp, top, eo, igroup, status
        logical :: found
        character(len=STDLEN) :: message
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('sigma_action') ) THROW_HARD('sigma2_convert requires sigma_action')
        if( .not. cline%defined('outfile') ) THROW_HARD('sigma2_convert requires outfile')
        call params%new(cline)
        call project%read(params%projfile)
        call project%ptr2oritype(params%oritype, particles)
        nptcls = particles%get_noris(consider_state=.false.)
        if( nptcls < 1 ) THROW_HARD('sigma2_convert requires a populated particle field')
        box = project%get_box()
        if( box < 2 ) THROW_HARD('sigma2_convert requires valid project dimensions')
        layout_digest = sigma2_state_project_layout_digest(project, particles)
        if( layout_digest == 0_int64 ) THROW_HARD('cannot derive canonical sigma2 particle layout identity')
        project_ngroups = 0
        do iptcl = 1, nptcls
            if( particles%get_state(iptcl) <= 0 ) cycle
            project_ngroups = max(project_ngroups, particles%get_int(iptcl, 'stkind'))
        enddo
        select case(trim(params%sigma_action))
            case('star_export')
                call project%get_sigma2_state_path(state_path, found)
                if( .not. found ) THROW_HARD('project has no canonical sigma2 state path')
                call sigma2_state_validate_file(state_path%to_char(), status, message, deep=.true.)
                if( status /= 0 ) THROW_HARD(trim(message))
                call sigma2_state_read_header(state_path%to_char(), header, status, message)
                if( status /= 0 ) THROW_HARD(trim(message))
                if( header%state /= SIGMA2_STATE_COMMITTED ) THROW_HARD('sigma2 export requires committed state')
                call sigma2_state_validate_identity(state_path%to_char(), box, project%get_smpd(), 1, &
                    &fdim(box)-1, nptcls, layout_digest, status, message, &
                    &expected_state=SIGMA2_STATE_COMMITTED)
                if( status /= 0 ) THROW_HARD(trim(message))
                call sigma2_state_read_groups(state_path%to_char(), stored_groups, status, message)
                if( status /= 0 ) THROW_HARD(trim(message))
                allocate(star_groups(2,int(header%ngroups),int(header%kfrom):int(header%kto)))
                do igroup = 1, int(header%ngroups)
                    do eo = 1, 2
                        star_groups(eo,igroup,:) = real(stored_groups(:,eo,igroup))
                    enddo
                enddo
                call write_groups_starfile(params%outfile, star_groups, int(header%ngroups))
                deallocate(star_groups, stored_groups)
                write(logfhandle,'(A)') '>>> SIGMA2 CONVERTER EXPORTED GROUPS TO: '//params%outfile%to_char()
            case('star_import','parts_import')
                if( .not. cline%defined('infile') ) THROW_HARD('sigma2 import requires infile')
                state_path = params%outfile
                if( file_exists(state_path) ) THROW_HARD('sigma2 import refuses to replace an existing output state')
                kfromto = [1, fdim(box)-1]
                allocate(spectra(kfromto(1):kfromto(2),nptcls), source=0.0_real32)
                if( trim(params%sigma_action) == 'star_import' )then
                    call read_sigma2_groups_file(params%infile, star_groups, kfromto, ngroups)
                    if( kfromto(1) /= 1 .or. kfromto(2) /= fdim(box)-1 ) &
                        &THROW_HARD('STAR sigma2 seed does not cover the native project shell range')
                    if( params%l_sigma_glob .and. ngroups /= 1 ) &
                        &THROW_HARD('global sigma estimation requires a one-group STAR seed')
                    if( .not. params%l_sigma_glob )then
                        if( project_ngroups < 1 ) THROW_HARD('cannot derive sigma2 stack groups from project')
                        if( ngroups > 1 .and. ngroups /= project_ngroups ) &
                            &THROW_HARD('STAR sigma2 group count does not match the project')
                    endif
                    do iptcl = 1, nptcls
                        eo = particles%get_eo(iptcl) + 1
                        if( ngroups == 1 )then
                            igroup = 1
                        else if( particles%get_state(iptcl) <= 0 )then
                            ! Inactive rows are excluded from reduction and consumption.
                            ! Give them a valid seed without imposing active-group coverage.
                            igroup = 1
                        else
                            igroup = particles%get_int(iptcl, 'stkind')
                            if( igroup < 1 .or. igroup > ngroups ) &
                                &THROW_HARD('STAR sigma2 groups do not cover project stack membership')
                        endif
                        spectra(:,iptcl) = real(star_groups(eo,igroup,:), real32)
                    enddo
                    deallocate(star_groups)
                else
                    allocate(covered(nptcls), source=.false.)
                    do ipart = 1, params%nparts
                        part_path = params%infile//int2str_pad(ipart,params%numlen)//'.dat'
                        if( .not. file_exists(part_path) ) THROW_HARD('missing legacy sigma2 part: '//part_path%to_char())
                        call part_file%new_from_file(part_path)
                        call part_file%get_resrange(kfromto)
                        if( kfromto(1) /= 1 .or. kfromto(2) /= fdim(box)-1 ) &
                            &THROW_HARD('legacy sigma2 part does not cover the native project shell range')
                        call part_file%read(part_values)
                        fromp = lbound(part_values,2)
                        top   = ubound(part_values,2)
                        if( fromp < 1 .or. top > nptcls .or. any(covered(fromp:top)) ) &
                            &THROW_HARD('legacy sigma2 parts have invalid or overlapping particle ranges')
                        spectra(:,fromp:top) = real(part_values,real32)
                        covered(fromp:top) = .true.
                        deallocate(part_values)
                        call part_file%kill
                        call part_path%kill
                    enddo
                    if( .not. all(covered) ) THROW_HARD('legacy sigma2 parts do not exactly cover the project')
                    deallocate(covered)
                endif
                call publish_import(trim(params%sigma_action) == 'star_import')
            case default
                THROW_HARD('unsupported sigma2 conversion action')
        end select
        call project%kill
        call state_path%kill
        call candidate_path%kill
        nullify(particles)
        call simple_end('**** SIMPLE_SIGMA2_CONVERT NORMAL STOP ****')

      contains

        subroutine publish_import( lossy_seed )
            logical, intent(in) :: lossy_seed
            real :: smpd
            smpd = project%get_smpd()
            allocate(active(nptcls), eo_ids(nptcls), group_ids(nptcls))
            ngroups = 0
            do iptcl = 1, nptcls
                active(iptcl)    = particles%get_state(iptcl) > 0
                eo_ids(iptcl)    = particles%get_eo(iptcl)
                group_ids(iptcl) = particles%get_int(iptcl, 'stkind')
                if( active(iptcl) ) ngroups = max(ngroups, group_ids(iptcl))
            enddo
            if( params%l_sigma_glob )then
                grouping = SIGMA2_GROUP_GLOBAL
                ngroups = 1
                group_ids = 1
            else
                grouping = SIGMA2_GROUP_STACK
                if( ngroups < 1 ) THROW_HARD('cannot derive sigma2 stack groups from project')
            endif
            if( lossy_seed )then
                call sigma2_state_init_header(header, 1, fdim(box)-1, nptcls, box, smpd, ngroups, &
                    &grouping, 1_int64, layout_digest, SIGMA2_PROV_STAR_SEED)
            else
                call sigma2_state_init_header(header, 1, fdim(box)-1, nptcls, box, smpd, ngroups, &
                    &grouping, 1_int64, layout_digest, SIGMA2_PROV_LEGACY_PARTS)
            endif
            candidate_path = sigma2_state_candidate_path(state_path%to_char())
            call sigma2_state_create_candidate(candidate_path%to_char(), header, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call sigma2_state_write_particles(candidate_path%to_char(), 1, spectra, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call sigma2_state_reduce_groups(candidate_path%to_char(), active, eo_ids, group_ids, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call sigma2_state_commit(candidate_path%to_char(), state_path%to_char(), active, eo_ids, &
                &group_ids, status, message)
            if( status /= 0 ) THROW_HARD(trim(message))
            call project%set_sigma2_state_path(state_path)
            call project%write_segment_inside('projinfo', params%projfile)
            deallocate(spectra, active, eo_ids, group_ids)
            if( lossy_seed )then
                write(logfhandle,'(A)') '>>> SIGMA2 CONVERTER IMPORTED LOSSY GROUPED STAR SEED: '//state_path%to_char()
            else
                write(logfhandle,'(A)') '>>> SIGMA2 CONVERTER IMPORTED EXACT LEGACY PARTICLE PARTS: '//state_path%to_char()
            endif
        end subroutine publish_import

    end subroutine exec_sigma2_convert

end module simple_commanders_euclid

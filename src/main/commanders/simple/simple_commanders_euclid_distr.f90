!@descr: for distributed sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid_distr
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_commanders_api
use simple_sigma2_binfile, only: sigma2_binfile
use simple_sigma2_state, only: sigma2_state_candidate_path, sigma2_state_commit, &
    &sigma2_state_project_layout_digest, sigma2_state_reduce_groups
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_create_candidate, &
    &sigma2_state_init_header, sigma2_state_read_header, sigma2_state_read_particles, &
    &sigma2_state_write_particles, sigma2_state_validate_file, &
    &SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK, SIGMA2_PROV_PSPEC, SIGMA2_STATE_FNAME
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
        type(sigma2_state_header)        :: state_header, previous_header
        type(string)                     :: part_fname,starfile_fname,outbin_fname
        type(string)                     :: state_path,candidate_path,cwd
        integer(int64)                   :: generation,layout_digest,prefix_digest
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup,nstks,nyq,pspec_l,pspec_u
        integer                          :: prefix_first,prefix_last
        integer                          :: state_status
        real(dp),          allocatable   :: group_pspecs(:,:,:), bootstrap_pspecs(:,:)
        real,              allocatable   :: pspec_ave(:),sigma2_part(:,:),sigma2_output(:,:)
        real(real32),      allocatable   :: prefix_spectra(:,:)
        integer,           allocatable   :: bootstrap_weights(:),eo_ids(:),group_ids(:)
        logical,           allocatable   :: pspec_covered(:),active(:)
        logical                          :: state_path_found, preserve_prefix
        character(len=STDLEN)            :: state_message
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
        ! Bootstrap once from a global E/O sample. The resulting curves are
        ! replicated into the requested group layout; matcher residual updates
        ! establish group-specific spectra after the first search pass.
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
        allocate(bootstrap_pspecs(2,nyq), source=0.d0)
        allocate(bootstrap_weights(2), source=0)
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
                if( .not.all(ieee_is_finite(sigma2_part(:,iptcl)) ) )then
                    write(logfhandle,*) 'iptcl/eo: ', iptcl, eo
                    write(logfhandle,*) 'finite: ', all(ieee_is_finite(sigma2_part(:,iptcl)))
                    THROW_HARD('active particle sigma spectrum was not computed; exec_calc_pspec_assemble')
                endif
                bootstrap_pspecs(eo+1,:) = bootstrap_pspecs(eo+1,:) + real(sigma2_part(:, iptcl),dp)
                bootstrap_weights(eo+1)  = bootstrap_weights(eo+1) + 1
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
            if( bootstrap_weights(eo) > 0 )then
                bootstrap_pspecs(eo,:) = bootstrap_pspecs(eo,:) / real(bootstrap_weights(eo),dp)
                bootstrap_pspecs(eo,:) = bootstrap_pspecs(eo,:) - real(pspec_ave(:),dp)
                call remove_negative_sigmas(eo)
            else
                bootstrap_pspecs(eo,:) = 1.d0
            endif
            do igroup = 1,ngroups
                group_pspecs(eo,igroup,:) = bootstrap_pspecs(eo,:)
            enddo
        enddo
        if( params%l_sigma_canonical )then
            if( trim(params%oritype) /= 'ptcl2D' .and. trim(params%oritype) /= 'ptcl3D' ) &
                &THROW_HARD('canonical sigma2 bootstrap requires a particle orientation field')
            call build%spproj%get_sigma2_state_path(state_path, state_path_found)
            if( .not. state_path_found )then
                call simple_getcwd(cwd)
                state_path = filepath(cwd, string(SIGMA2_STATE_FNAME))
            endif
            candidate_path = sigma2_state_candidate_path(state_path%to_char())
            generation = 1_int64
            preserve_prefix = .false.
            if( file_exists(state_path) )then
                call sigma2_state_read_header(state_path%to_char(), previous_header, state_status, state_message)
                if( state_status == 0 )then
                    generation = previous_header%generation + 1_int64
                    if( previous_header%nptcls > 0 .and. previous_header%nptcls < nptcls )then
                        prefix_digest = sigma2_state_project_layout_digest(build%spproj, build%spproj_field, &
                            &int(previous_header%nptcls))
                        preserve_prefix = prefix_digest == previous_header%layout_digest .and. &
                            &previous_header%box == params%box .and. &
                            &previous_header%kfrom == params%kfromto(1) .and. &
                            &previous_header%kto == params%kfromto(2) .and. &
                            &abs(real(previous_header%smpd)-params%smpd) <= 1.e-5*max(1.,params%smpd)
                        if( preserve_prefix )then
                            call sigma2_state_validate_file(state_path%to_char(), state_status, state_message, deep=.true.)
                            preserve_prefix = state_status == 0
                        endif
                    endif
                endif
            endif
            layout_digest = sigma2_state_project_layout_digest(build%spproj, build%spproj_field)
            if( layout_digest == 0_int64 ) THROW_HARD('cannot derive canonical sigma2 particle layout identity')
            if( params%l_sigma_glob )then
                call sigma2_state_init_header(state_header, params%kfromto(1), params%kfromto(2), &
                    &nptcls, params%box, params%smpd, ngroups, SIGMA2_GROUP_GLOBAL, generation, &
                    &layout_digest, SIGMA2_PROV_PSPEC)
            else
                call sigma2_state_init_header(state_header, params%kfromto(1), params%kfromto(2), &
                    &nptcls, params%box, params%smpd, ngroups, SIGMA2_GROUP_STACK, generation, &
                    &layout_digest, SIGMA2_PROV_PSPEC)
            endif
            call sigma2_state_create_candidate(candidate_path%to_char(), state_header, state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            allocate(active(nptcls), eo_ids(nptcls), group_ids(nptcls))
            do iptcl = 1, nptcls
                active(iptcl)    = build%spproj_field%get_state(iptcl) > 0
                eo_ids(iptcl)    = build%spproj_field%get_eo(iptcl)
                group_ids(iptcl) = build%spproj_field%get_int(iptcl, 'stkind')
            enddo
        else
            ! write group sigmas to the legacy STAR history
            if( cline%defined('which_iter') )then
                starfile_fname = SIGMA2_GROUP_FBODY//int2str(params%which_iter)//STAR_EXT
            else
                starfile_fname = SIGMA2_GROUP_FBODY//'1'//STAR_EXT
            endif
            call write_groups_starfile(starfile_fname, real(group_pspecs), ngroups)
        endif
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
                sigma2_output(params%kfromto(1):params%kfromto(2),iptcl) =&
                    &real(bootstrap_pspecs(eo+1,params%kfromto(1):params%kfromto(2)))
            end do
            if( params%l_sigma_canonical )then
                call sigma2_state_write_particles(candidate_path%to_char(), pspec_l, &
                    &real(sigma2_output,real32), state_status, state_message)
                if( state_status /= 0 ) THROW_HARD(trim(state_message))
                call del_file(part_fname)
            else
                outbin_fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
                call binfile%new(outbin_fname, fromp=pspec_l, top=pspec_u, &
                    &kfromto=[params%kfromto(1), params%kfromto(2)])
                call binfile%write(sigma2_output)
            endif
            deallocate(sigma2_part)
        end do
        if( params%l_sigma_canonical )then
            if( preserve_prefix )then
                do prefix_first = 1, int(previous_header%nptcls), 4096
                    prefix_last = min(prefix_first+4095, int(previous_header%nptcls))
                    call sigma2_state_read_particles(state_path%to_char(), prefix_first, prefix_last, &
                        &prefix_spectra, state_status, state_message)
                    if( state_status /= 0 ) THROW_HARD(trim(state_message))
                    call sigma2_state_write_particles(candidate_path%to_char(), prefix_first, &
                        &prefix_spectra, state_status, state_message)
                    if( state_status /= 0 ) THROW_HARD(trim(state_message))
                    deallocate(prefix_spectra)
                enddo
                write(logfhandle,'(A,I0,A)') '>>> SIGMA2 APPEND: retained ', previous_header%nptcls, &
                    &' committed particle spectra and bootstrapped the appended suffix'
            endif
            call sigma2_state_reduce_groups(candidate_path%to_char(), active, eo_ids, group_ids, &
                &state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            call sigma2_state_commit(candidate_path%to_char(), state_path%to_char(), active, eo_ids, &
                &group_ids, state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            call build%spproj%set_sigma2_state_path(state_path)
            call build%spproj%write_segment_inside('projinfo', params%projfile)
            deallocate(active, eo_ids, group_ids)
        endif
        ! end gracefully
        deallocate(group_pspecs,bootstrap_pspecs,bootstrap_weights,pspec_covered)
        if( allocated(sigma2_output) ) deallocate(sigma2_output)
        if( allocated(pspec_ave) ) deallocate(pspec_ave)
        call binfile%kill
        call state_path%kill
        call candidate_path%kill
        call cwd%kill
        call build%kill_general_tbox
        call simple_touch('CALC_PSPEC_FINISHED')
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine remove_negative_sigmas(eo)
            integer, intent(in) :: eo
            logical :: is_positive
            logical :: fixed_from_prev
            integer :: nn, idx
            where( .not. ieee_is_finite(bootstrap_pspecs(eo,:)) ) bootstrap_pspecs(eo,:) = 1.d0
            if( .not. any(bootstrap_pspecs(eo,:) > DTINY) )then
                write(logfhandle,*) '>>> WARNING: calc_pspec_assemble floored empty global sigma2 spectrum to 1.0; ', &
                    'eo/weight: ', eo, bootstrap_weights(eo)
                bootstrap_pspecs(eo,:) = 1.d0
                return
            endif
            ! remove any negative sigma2 noise values: replace by positive neighboring value
            do idx = 1, size(bootstrap_pspecs, 2)
                if( bootstrap_pspecs(eo,idx) <= DTINY )then
                    ! first try the previous value
                    fixed_from_prev = .false.
                    if( idx - 1 >= 1 )then
                        if( bootstrap_pspecs(eo,idx-1) > DTINY )then
                            bootstrap_pspecs(eo,idx) = bootstrap_pspecs(eo,idx-1)
                            fixed_from_prev = .true.
                        end if
                    end if
                    if( .not. fixed_from_prev )then
                        is_positive = .false.
                        nn          = idx
                        do while (.not. is_positive)
                            nn = nn + 1
                            if( nn > size(bootstrap_pspecs,2) )then
                                bootstrap_pspecs(eo,idx) = 1.d0
                                exit
                            end if
                            if( bootstrap_pspecs(eo,nn) > DTINY )then
                                is_positive = .true.
                                bootstrap_pspecs(eo,idx) = bootstrap_pspecs(eo,nn)
                            end if
                        end do
                    end if
                end if
            end do
        end subroutine remove_negative_sigmas

    end subroutine exec_calc_pspec_assemble

end module simple_commanders_euclid_distr

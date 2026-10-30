!@descr: shared startup and continuation controller for refine3D_auto and refine3D_local
module simple_refine3D_workflow_controller
use simple_string,              only: string
use simple_refine3D_provenance, only: refine3D_provenance
implicit none
private
#include "simple_local_flags.inc"
public :: refine3D_workflow_controller
type :: refine3D_workflow_controller
    private
    type(string)              :: source_ref_vol, source_fsc_fname, source_sigma_state
    type(refine3D_provenance) :: source_reference
    character(len=14)         :: workflow_label = ''
    logical                   :: l_local = .false.
    logical                   :: l_external_input = .false.
    logical                   :: l_source_support_verified = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: preflight
    procedure :: load_project
    procedure :: publish_project
    procedure :: prepare_startup
    procedure, private :: preflight_cont_startup
    procedure, private :: stage_local_trailing_chain
    procedure, private :: prepare_external_init_vol
end type refine3D_workflow_controller

contains

    subroutine new( self, l_local, l_external_input )
        class(refine3D_workflow_controller), intent(inout) :: self
        logical, intent(in)                                :: l_local
        logical, intent(in)                                :: l_external_input

        call self%kill
        self%l_local = l_local
        self%l_external_input = l_external_input
        self%workflow_label = 'REFINE3D_AUTO'
        if( l_local ) self%workflow_label = 'REFINE3D_LOCAL'
    end subroutine new

    subroutine kill( self )
        class(refine3D_workflow_controller), intent(inout) :: self

        call self%source_reference%kill
        call self%source_ref_vol%kill
        call self%source_fsc_fname%kill
        call self%source_sigma_state%kill
        self%workflow_label = ''
        self%l_local = .false.
        self%l_external_input = .false.
        self%l_source_support_verified = .false.
    end subroutine kill

    subroutine preflight( self, cline )
        use simple_cmdline, only: cmdline
        class(refine3D_workflow_controller), intent(inout) :: self
        class(cmdline), intent(inout)                      :: cline

        self%l_source_support_verified = .false.
        if( self%l_local ) call self%preflight_cont_startup(cline)
    end subroutine preflight

    subroutine load_project( self, params, spproj )
        use simple_parameters, only: parameters
        use simple_sp_project, only: sp_project
        class(refine3D_workflow_controller), intent(inout) :: self
        class(parameters), intent(inout)                   :: params
        class(sp_project), intent(inout)                   :: spproj

        call spproj%read( params%projfile )
        call self%publish_project(params, spproj)
    end subroutine load_project

    subroutine publish_project( self, params, spproj )
        use simple_refine3D_strategy, only: stage_refine3D_local_sigma2
        use simple_parameters,        only: parameters
        use simple_sp_project,        only: sp_project
        class(refine3D_workflow_controller), intent(inout) :: self
        class(parameters), intent(inout)                   :: params
        class(sp_project), intent(inout)                   :: spproj

        call spproj%update_projinfo(params%projfile)
        if( self%l_local ) call stage_refine3D_local_sigma2(params, spproj, self%source_sigma_state)
    end subroutine publish_project

    subroutine prepare_startup( self, params, cline, spproj, init_vol, l_have_init_vol, box_crop, smpd_crop )
        use simple_refine3D_strategy, only: stage_refine3D_local_references
        use simple_refine3D_fnames,   only: refine3D_state_vol_fname
        use simple_syslib,            only: file_exists
        use simple_fileio,            only: add2fbody
        use simple_parameters,        only: parameters
        use simple_sp_project,        only: sp_project
        use simple_cmdline,           only: cmdline
        use simple_defs,              only: logfhandle
        use simple_defs_fname,        only: MRC_EXT
        use simple_error,             only: simple_exception
        class(refine3D_workflow_controller), intent(inout) :: self
        class(parameters), intent(inout)                   :: params
        class(cmdline), intent(inout)                      :: cline
        class(sp_project), intent(inout)                   :: spproj
        type(string), intent(inout)                        :: init_vol
        logical, intent(out)                               :: l_have_init_vol
        integer, intent(in)                                :: box_crop
        real, intent(in)                                   :: smpd_crop
        real    :: init_smpd
        integer :: init_box
        logical :: l_reused_derivatives

        l_have_init_vol = .false.
        init_box        = 0
        init_smpd       = 0.

        if( self%l_local )then
            if( .not. file_exists(self%source_ref_vol) )then
                THROW_HARD('REFINE3D_LOCAL starting volume does not exist: '//self%source_ref_vol%to_char())
            endif
            call stage_refine3D_local_references(params, spproj, self%source_reference, self%source_fsc_fname, &
                &self%l_source_support_verified, init_vol, l_reused_derivatives)
            call spproj%write_segment_inside('projinfo', params%projfile)
            if( l_reused_derivatives ) call spproj%write_segment_inside('ptcl3D', params%projfile)
            call spproj%write_segment_inside('out', params%projfile)
            params%vols(1)      = init_vol
            params%vols_even(1) = add2fbody(init_vol, MRC_EXT, '_even')
            params%vols_odd(1)  = add2fbody(init_vol, MRC_EXT, '_odd')
            l_have_init_vol = .true.
            call cline%set('vol1', init_vol)
            write(logfhandle,'(A,1X,A)') '>>> REFINE3D_LOCAL USING STAGED SAVED REFERENCES:', init_vol%to_char()
        else if( cline%defined('vol1') )then
            init_vol = params%vols(1)
            if( .not. file_exists(init_vol) )then
                THROW_HARD('File: '//init_vol%to_char()//' does not exist! '//self%workflow_label)
            endif
            call self%prepare_external_init_vol(params, spproj, init_vol)
            l_have_init_vol = .true.
            call cline%delete('vol1') ! because now part of the project
            write(logfhandle,'(A,1X,A)') '>>> '//self%workflow_label//' USING INPUT VOLUME:', init_vol%to_char()
            if( trim(params%ref_pose_init) /= 'cc' ) write(logfhandle,'(A)') &
                &'>>> WARNING: EXTERNAL VOL1 IS TRUSTED WITHOUT CC POSE INITIALIZATION (REF_POSE_INIT=NONE)'
        else
            if( spproj%isthere_in_osout('vol', 1) )then
                call spproj%get_vol('vol', 1, init_vol, init_smpd, init_box)
                if( file_exists(init_vol) )then
                    if( reference_grid_matches(params, init_box, init_smpd) )then
                        l_have_init_vol = .true.
                        write(logfhandle,'(A,1X,A)') '>>> '//self%workflow_label//' USING PROJECT VOLUME:', init_vol%to_char()
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> PROJECT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
                    else
                        write(logfhandle,'(A,1X,A)') &
                            &'>>> '//self%workflow_label//' PROJECT VOLUME SAMPLING MISMATCH, RECONSTRUCTING:', init_vol%to_char()
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> PROJECT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
                        write(logfhandle,'(A,I0,A,F8.4)') '>>> CURRENT RUN BOX/SMPD:    ', params%box, '/', params%smpd
                    endif
                else
                    write(logfhandle,'(A,1X,A)') &
                        &'>>> '//self%workflow_label//' PROJECT VOLUME MISSING, RECONSTRUCTING:', init_vol%to_char()
                endif
            endif
        endif
        if( self%l_local ) call self%stage_local_trailing_chain(params, cline, spproj, &
            &self%source_ref_vol, box_crop, smpd_crop)
    end subroutine prepare_startup

    ! Read source headers and project only; do not write maps or half assignments.
    subroutine preflight_cont_startup( self, cline )
        use simple_syslib,              only: file_exists, simple_abspath, resolve_project_files
        use simple_fileio,              only: get_fpath, add2fbody
        use simple_oris,                only: oris
        use simple_refine3D_strategy,   only: select_refine3D_state_rows
        use simple_strategy3D_cont,     only: check_cont_seeds, cont_reference_geometry_is_valid
        use simple_imghead,             only: get_mrcfile_info
        use simple_state_weight_set,    only: state_weight_set
        use simple_sigma2_files,        only: canonical_sigma2_consumable
        use simple_refine3D_provenance, only: REFSET_OK, REFSET_LEGACY
        use simple_sp_project,          only: sp_project
        use simple_cmdline,             only: cmdline
        use simple_defs,                only: STDLEN, logfhandle
        use simple_defs_fname,          only: MRC_EXT
        use simple_error,               only: simple_exception
        class(refine3D_workflow_controller), intent(inout) :: self
        class(cmdline), intent(inout)                      :: cline
        class(sp_project), allocatable       :: source
        class(state_weight_set), allocatable :: weights
        logical, allocatable                 :: selected(:)
        logical                              :: use_flex, source_sigma_found
        type(string)                         :: project_file, previous_dir, volume, names(3), project_dir, relative_volume
        class(oris), pointer                 :: source_os
        character(len=STDLEN)                :: validation_message
        type(string), allocatable            :: projects(:)
        integer                              :: next_dir, eo, state, i, dims(3,3), native_box, registered_box, halves(0:1), particle_state
        integer                              :: provenance_status
        real                                 :: sampling(3), native_smpd, registered_smpd

        if( cline%defined('projfile') )then
            project_file = cline%get_carg('projfile')
        else
            ! Match parameter-side project discovery, then pin that choice before directory creation.
            call resolve_project_files(string('.'), previous_dir, next_dir, projects)
            if( .not. allocated(projects) ) THROW_HARD('Continuous startup requires an existing project')
            if( size(projects) /= 1 ) THROW_HARD('Continuous startup requires an unambiguous projfile')
            project_file = projects(1)
            call cline%set('projfile', project_file)
        endif

        allocate(source)
        call source%read(project_file)
        call source%update_projinfo(project_file)
        self%source_fsc_fname = ''
        state = 1
        if( cline%defined('state') )then
            state = cline%get_iarg('state')
            if( state < 1 ) THROW_HARD('Continuous startup requires a positive state')
            use_flex = .false.
            if( cline%defined('m_estimator') ) use_flex = cline%get_carg('m_estimator') == 'flex'
            allocate(weights)
            call select_refine3D_state_rows(source, state, use_flex, selected, weights)
        else
            allocate(selected(source%os_ptcl3D%get_noris()), source=.true.)
            do i = 1, size(selected)
                if( source%os_ptcl3D%get_state(i) > 1 )then
                    THROW_HARD('refine3D_local needs state=X for a multistate project')
                endif
            enddo
        endif
        call check_cont_seeds(source%os_ptcl3D, selected)

        ! Seed validation above guarantees labels 0/1, which index these two counters.
        halves = 0
        do i = 1, size(selected)
            if( .not. selected(i) ) cycle
            particle_state = source%os_ptcl3D%get_state(i)
            if( particle_state <= 0 ) cycle
            eo = source%os_ptcl3D%get_eo(i)
            halves(eo) = halves(eo) + 1
        enddo
        if( any(halves == 0) ) THROW_HARD('Continuous startup requires particles in both half-sets')

        native_box = source%get_box()
        native_smpd = source%get_smpd()
        if( cline%defined('vol1') )then
            volume = cline%get_carg('vol1')
        else
            call source%get_vol('vol', state, volume, registered_smpd, registered_box)
            ! A stored relative map path may refer to the source project directory.
            if( .not. file_exists(volume) )then
                project_dir = get_fpath(project_file)
                relative_volume = project_dir%to_char()//volume%to_char()
                if( file_exists(relative_volume) ) volume = relative_volume
            endif
        endif
        self%source_ref_vol = simple_abspath(volume)
        ! Require the original pair before initialization can manufacture missing halves.
        names(1) = volume
        names(2) = add2fbody(volume, MRC_EXT, '_even')
        names(3) = add2fbody(volume, MRC_EXT, '_odd')
        do i = 1,3
            if( .not. file_exists(names(i)) )then
                THROW_HARD('Continuous startup requires an existing reference: '//names(i)%to_char())
            endif
            call get_mrcfile_info(names(i), dims(:,i), 'M', sampling(i), .false.)
        enddo
        if( .not. cont_reference_geometry_is_valid(dims, sampling, native_box, native_smpd) ) &
            THROW_HARD('Continuous startup reference dimensions, sampling, or physical extent are incompatible')

        call self%source_reference%new(self%source_ref_vol, 'reference', validation_message)
        provenance_status = self%source_reference%validate_reference_set(self%source_fsc_fname, &
            &self%l_source_support_verified, validation_message)
        select case(provenance_status)
            case(REFSET_OK)
                if( .not. self%source_reference%reference_matches_project(source, state, &
                    &self%l_external_input, validation_message) )then
                    THROW_HARD('Continuous startup reference provenance failed: '//trim(validation_message))
                endif
                if( len_trim(validation_message) > 0 ) write(logfhandle,'(A)') '>>> REFINE3D_LOCAL '//trim(validation_message)
                write(logfhandle,'(A)') '>>> REFINE3D_LOCAL VERIFIED THE SOURCE REFERENCE SET'
            case(REFSET_LEGACY)
                write(logfhandle,'(A)') '>>> REFINE3D_LOCAL LEGACY BASE PAIR: REBUILDING UNVERIFIED FSC AND DERIVED PRODUCTS'
            case default
                THROW_HARD('Continuous startup reference provenance failed: '//trim(validation_message))
        end select

        call source%ptr2oritype('ptcl3D', source_os)
        if( canonical_sigma2_consumable(source, source_os, native_box, native_smpd, validation_message) )then
            call source%get_sigma2_state_path(self%source_sigma_state, source_sigma_found)
            if( source_sigma_found ) self%source_sigma_state = simple_abspath(self%source_sigma_state)
        endif
        nullify(source_os)
        call source%kill
    end subroutine preflight_cont_startup

    ! Select chain reuse or a full first update using the caller's loaded project.
    subroutine stage_local_trailing_chain( self, params, cline, project, source_vol, box_crop, smpd_crop )
        use simple_refine3D_strategy, only: stage_refine3D_local_chain
        use simple_parameters,        only: parameters
        use simple_sp_project,        only: sp_project
        use simple_cmdline,           only: cmdline
        use simple_defs,              only: logfhandle
        class(refine3D_workflow_controller), intent(inout) :: self
        class(parameters), intent(in)                      :: params
        class(cmdline), intent(inout)                      :: cline
        class(sp_project), intent(inout)                   :: project
        integer, intent(in)                                :: box_crop
        real, intent(in)                                   :: smpd_crop
        type(string), intent(in)                           :: source_vol
        logical :: valid

        if( .not. params%l_update_frac .or. .not. params%l_trail_rec ) return
        valid = stage_refine3D_local_chain(params, project, source_vol, box_crop, smpd_crop)
        if( valid )then
            call cline%set('trail_seed', 'no')
            write(logfhandle,'(A)') '>>> REFINE3D_LOCAL REUSING A COMPATIBLE TRAILING ACCUMULATOR CHAIN'
        else
            call cline%set('trail_seed', 'yes')
            write(logfhandle,'(A)') '>>> REFINE3D_LOCAL WILL USE A FULL FIRST UPDATE TO SEED THE TRAILING CHAIN'
        endif
    end subroutine stage_local_trailing_chain

    ! Prepare external even/odd volumes and FSC for refinement and project registration.
    subroutine prepare_external_init_vol( self, params, spproj, init_vol )
        use simple_syslib,              only: file_exists
        use simple_fileio,              only: add2fbody, arr2file
        use simple_imghead,             only: find_ldim_nptcls, find_img_smpd
        use simple_image,               only: image
        use simple_refine3D_fnames,     only: refine3D_startvol_fname, refine3D_fsc_fname
        use simple_halfmap_diagnostics, only: remove_support_provenance
        use simple_parameters,          only: parameters
        use simple_sp_project,          only: sp_project
        use simple_defs,                only: logfhandle
        use simple_defs_fname,          only: MRC_EXT
        use simple_error,               only: simple_exception
        class(refine3D_workflow_controller), intent(inout) :: self
        class(parameters), intent(in)                      :: params
        class(sp_project), intent(inout)                   :: spproj
        type(string), intent(inout)                        :: init_vol
        type(string)      :: init_even, init_odd, new_vol, new_even, new_odd
        type(image)       :: vol, vol_even, vol_odd
        real, allocatable :: fsc(:)
        real              :: ave, stdev, maxv, minv, msk, v
        integer           :: ldim(3), n, init_box
        real              :: init_smpd
        init_even = add2fbody(init_vol, MRC_EXT, '_even')
        init_odd  = add2fbody(init_vol, MRC_EXT, '_odd')
        if( .not.(file_exists(init_even) .and. file_exists(init_odd)) )then
            THROW_HARD('Expected even/odd half-volumes for input volume '//init_vol%to_char()//' not found;&
            & looking for: '//init_even%to_char()//' and '//init_odd%to_char())
        endif
        call find_ldim_nptcls(init_vol, ldim, n)
        init_smpd = find_img_smpd(init_vol)
        init_box  = ldim(1)
        if( reference_grid_matches(params, init_box, init_smpd) )then
            write(logfhandle,'(A,1X,A)') '>>> '//self%workflow_label//' USING E/O INPUT VOLUMES'
            write(logfhandle,'(A,I0,A,F8.4)') '>>> INPUT VOLUME BOX/SMPD: ', init_box, '/', init_smpd
        else
            THROW_HARD('Input e/o volumes must have same dimensions as image')
        endif
        call vol%new(ldim, init_smpd)
        call vol%read(init_vol)
        call vol_even%new(ldim, init_smpd)
        call vol_even%read(init_even)
        if( (abs(vol_even%get_smpd() - init_smpd) > 1.e-6) .or. (vol_even%get_box() /= init_box) )then
            THROW_HARD('Even half-volumes for input volume '//init_vol%to_char()//' have different dimensions/sampling')
        endif
        call vol_odd%new(ldim, init_smpd)
        call vol_odd%read(init_odd)
        if( (abs(vol_odd%get_smpd() - init_smpd) > 1.e-6) .or. (vol_odd%get_box() /= init_box) )then
            THROW_HARD('Odd half-volumes for input volume '//init_vol%to_char()//' have different dimensions/sampling')
        endif
        ! calculate fsc and import into project
        new_vol  = refine3D_startvol_fname(1)
        new_even = add2fbody(new_vol, MRC_EXT, '_even')
        new_odd  = add2fbody(new_vol, MRC_EXT, '_odd')
        ! normalization of volumes
        msk = 0.5 * params%mskdiam / params%smpd
        call vol_even%stats('foreground', ave, stdev, maxv, minv, msk=msk)
        v = stdev * real(init_box) / real(params%box) ! foreground stdev to the data-quotient reference scale
        call vol_even%norm_ext(ave, v)
        call vol_odd%stats('foreground', ave, stdev, maxv, minv, msk=msk)
        v = stdev * real(init_box) / real(params%box) ! foreground stdev to the data-quotient reference scale
        call vol_odd%norm_ext(ave, v)
        call vol%stats('foreground', ave, stdev, maxv, minv, msk=msk)
        v = stdev * real(init_box) / real(params%box) ! foreground stdev to the data-quotient reference scale
        call vol%norm_ext(ave, v)
        call vol_even%write(new_even)
        call vol_odd%write(new_odd)
        call vol%write(new_vol)
        call remove_support_provenance(new_vol) ! an imported map carries no support record
        call vol_even%mask3D_soft(msk)
        call vol_odd%mask3D_soft(msk)
        call vol_even%fft
        call vol_odd%fft
        allocate(fsc(vol%get_lfny(1)), source=0.)
        call vol_even%fsc(vol_odd, fsc)
        call arr2file(fsc, refine3D_fsc_fname(1))
        if( all(fsc < 0.9) ) THROW_HARD('Calculated FSC is too low for refinement')
        init_vol = new_vol  ! renaming of global filename
        call spproj%add_vol2os_out(init_vol, init_smpd, 1, 'vol')
        call spproj%add_fsc2os_out(refine3D_fsc_fname(1), 1, init_box)
        call spproj%write_segment_inside('out', params%projfile)
        ! cleanup
        deallocate(fsc)
        call vol%kill; call vol_even%kill; call vol_odd%kill
        call init_even%kill; call init_odd%kill
        call new_vol%kill; call new_even%kill; call new_odd%kill
    end subroutine prepare_external_init_vol

    logical function reference_grid_matches( params, box, smpd ) result(compatible)
        use simple_parameters, only: parameters
        use simple_defs,       only: TINY
        class(parameters), intent(in) :: params
        integer, intent(in)           :: box
        real, intent(in)              :: smpd
        compatible = box == params%box .and. smpd > TINY .and. abs(smpd - params%smpd) <= 1.e-6
    end function reference_grid_matches

end module simple_refine3D_workflow_controller

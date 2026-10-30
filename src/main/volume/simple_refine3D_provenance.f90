!@descr: producer records and validation of local-refinement reference, sigma and accumulator handoffs
module simple_refine3D_provenance
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,                   only: STDLEN
use simple_defs_fname,             only: MRC_EXT, NUFILT_SUFFIX, NULOCRES_SUFFIX, NU_ENVMASK_FBODY, AUTOMASK_FBODY
use simple_error,                  only: simple_exception
use simple_string,                 only: string
use simple_string_utils,           only: int2str_pad
use simple_syslib,                 only: file_exists, simple_touch, del_file, simple_abspath
use simple_fileio,                 only: fopen, fclose, basename, get_fpath, add2fbody
use simple_parameters,             only: parameters
use simple_sp_project,             only: sp_project
use simple_oris,                   only: oris
use simple_ptcl_layout,            only: ptcl_layout_digest
use simple_halfmap_diagnostics,    only: support_provenance_fname
use simple_sigma2_state_file,      only: sigma2_state_digest_begin, sigma2_state_digest_integer, &
    &sigma2_state_digest_text, sigma2_state_digest_file
implicit none
private
public :: begin_reference_set, publish_reference_set, refine3D_provenance
public :: publish_reuse_provenance, reuse_provenance_matches, provenance_file, provenance_pending_file
public :: local_import_marker
public :: REFSET_OK, REFSET_LEGACY, REFSET_INVALID
#include "simple_local_flags.inc"

integer, parameter :: VERSION = 2, MAX_FILES = 16
integer, parameter :: REFSET_OK = 0, REFSET_LEGACY = 1, REFSET_INVALID = 2

type :: refine3D_provenance
    private
    type(string)              :: anchor
    type(string), allocatable :: files(:)
    integer(int64)            :: layout_digest = 0_int64, half_digest = 0_int64, population_digest = 0_int64
    integer(int64)            :: pose_digest = 0_int64, observation_digest = 0_int64, lineage_digest = 0_int64
    integer(int64)            :: settings_digest = 0_int64
    integer                   :: producer_state = 0
    real                      :: align_lp = 0.
    character(len=9)          :: kind = ''
    character(len=STDLEN)     :: validation_message = 'provenance record is not loaded'
    logical                   :: valid = .false., legacy = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: validate_reference_set
    procedure :: reference_matches_project
    procedure :: reference_derivatives_match
    procedure :: contains_file
    procedure :: get_anchor
    procedure, private :: read_record
    procedure, private :: reuse_matches
end type refine3D_provenance

contains

    subroutine new( self, anchor, kind, message )
        class(refine3D_provenance), intent(inout) :: self
        class(string),             intent(in)     :: anchor
        character(len=*),          intent(in)     :: kind
        character(len=*),          intent(out)    :: message
        call self%kill
        self%anchor = anchor
        self%kind = kind
        self%valid = self%read_record(message)
        self%validation_message = message
        self%legacy = .not. file_exists(provenance_file(anchor)) .and. &
            &.not. file_exists(provenance_pending_file(anchor))
    end subroutine new

    subroutine kill( self )
        class(refine3D_provenance), intent(inout) :: self
        call self%anchor%kill
        if( allocated(self%files) ) deallocate(self%files)
        self%layout_digest = 0_int64
        self%half_digest = 0_int64
        self%population_digest = 0_int64
        self%pose_digest = 0_int64
        self%observation_digest = 0_int64
        self%lineage_digest = 0_int64
        self%settings_digest = 0_int64
        self%producer_state = 0
        self%align_lp = 0.
        self%kind = ''
        self%validation_message = 'provenance record is not loaded'
        self%valid = .false.
        self%legacy = .false.
    end subroutine kill

    function get_anchor( self ) result(anchor)
        class(refine3D_provenance), intent(in) :: self
        type(string) :: anchor
        anchor = self%anchor
    end function get_anchor

    logical function contains_file( self, file ) result(found)
        class(refine3D_provenance), intent(in) :: self
        class(string),             intent(in)  :: file
        integer :: i
        found = .false.
        if( .not. self%valid ) return
        do i = 1,size(self%files)
            if( basename(self%files(i)) /= basename(file) ) cycle
            found = .true.
            return
        enddo
    end function contains_file

    function local_import_marker( volume ) result(fname)
        class(string), intent(in) :: volume
        type(string) :: fname
        fname = volume//'.local_import'
    end function local_import_marker

    function provenance_file( anchor ) result(fname)
        class(string), intent(in) :: anchor
        type(string) :: fname
        fname = anchor//'.provenance'
    end function provenance_file

    function provenance_pending_file( anchor ) result(fname)
        class(string), intent(in) :: anchor
        type(string) :: fname
        fname = provenance_file(anchor)//'.pending'
    end function provenance_pending_file

    ! The pending marker distinguishes interrupted publication from a legacy input.
    subroutine begin_reference_set( volume )
        class(string), intent(in) :: volume
        type(string) :: record, pending
        record = provenance_file(volume)
        pending = provenance_pending_file(volume)
        call simple_touch(pending)
        if( file_exists(record) ) call del_file(record)
        record = local_import_marker(volume)
        if( file_exists(record) ) call del_file(record)
    end subroutine begin_reference_set

    subroutine publish_reference_set( params, project, particles, state, volume, halves, fsc, align_lp )
        class(parameters), intent(in) :: params
        class(sp_project), intent(in) :: project
        class(oris),       intent(in) :: particles
        integer,           intent(in) :: state
        class(string),     intent(in) :: volume, halves(2), fsc
        real,              intent(in) :: align_lp
        type(string) :: files(MAX_FILES), name
        integer      :: n, i
        files(1) = volume
        files(2:3) = halves
        files(4) = fsc
        n = 4
        name = support_provenance_fname(volume)
        if( file_exists(name) ) call append(name)
        if( params%l_ml_reg )then
            do i = 1,2
                name = add2fbody(halves(i), MRC_EXT, '_unfil')
                if( file_exists(name) ) call append(name)
            enddo
        endif
        name = AUTOMASK_FBODY//int2str_pad(state,2)//MRC_EXT
        if( trim(params%automsk) /= 'no' .and. file_exists(name) ) call append(name)
        if( params%l_nonuniform )then
            name = add2fbody(volume, MRC_EXT, NUFILT_SUFFIX)
            call append(name)
            do i = 1,2
                name = add2fbody(halves(i), MRC_EXT, NUFILT_SUFFIX)
                call append(name)
            enddo
            name = add2fbody(volume, MRC_EXT, NULOCRES_SUFFIX)
            call append(name)
            name = NU_ENVMASK_FBODY//int2str_pad(state,2)//MRC_EXT
            if( file_exists(name) ) call append(name)
        endif
        call write_record(params, project, particles, 'reference', volume, files(:n), state, align_lp)
    contains
        subroutine append( file )
            class(string), intent(in) :: file
            n = n + 1
            files(n) = file
        end subroutine append
    end subroutine publish_reference_set

    integer function validate_reference_set( self, fsc, support_verified, message ) result(status)
        class(refine3D_provenance), intent(in)   :: self
        type(string),              intent(inout) :: fsc
        logical,                   intent(out)   :: support_verified
        character(len=*),          intent(out)   :: message
        type(string) :: expected
        integer      :: i
        call fsc%kill
        support_verified = .false.
        status = REFSET_INVALID
        message = self%validation_message
        if( self%kind /= 'reference' ) return
        if( self%legacy )then
            status = REFSET_LEGACY
            message = 'legacy base pair: derived-product provenance is unavailable'
            return
        endif
        if( .not. self%valid ) return
        if( size(self%files) < 4 ) return
        expected = add2fbody(self%anchor, MRC_EXT, '_even')
        message = 'reference constituent roles differ'
        if( basename(self%files(2)) /= basename(expected) ) return
        expected = add2fbody(self%anchor, MRC_EXT, '_odd')
        if( basename(self%files(3)) /= basename(expected) ) return
        fsc = self%files(4)
        expected = support_provenance_fname(self%anchor)
        do i = 5,size(self%files)
            if( basename(self%files(i)) == basename(expected) ) support_verified = .true.
        enddo
        message = ''
        status = REFSET_OK
    end function validate_reference_set

    ! A reference may lag poses, but its half partition must still match its source rows.
    logical function reference_matches_project( self, project, state, allow_external, message ) result(valid)
        class(refine3D_provenance), intent(in)  :: self
        class(sp_project),          intent(in)  :: project
        integer,                    intent(in)  :: state
        logical,                    intent(in)  :: allow_external
        character(len=*),           intent(out) :: message
        type(refine3D_provenance) :: current
        valid = .false.
        message = self%validation_message
        if( .not. self%valid ) return
        if( self%kind /= 'reference' ) return
        call particle_identity(project, project%os_ptcl3D, 'reference', current)
        if( self%layout_digest == 0_int64 .or. current%layout_digest == 0_int64 )then
            message = 'reference particle lineage is undefined'
            return
        endif
        if( self%layout_digest /= current%layout_digest )then
            message = 'reference ordered particle lineage differs from the source project'
            if( allow_external .and. self%lineage_digest /= current%lineage_digest )then
                message = 'explicit external reference: particle lineage differs; '//&
                    &'alignment and half relationship are caller supplied'
                valid = .true.
            endif
            return
        endif
        message = 'reference half assignments or active population differ from the source project'
        if( self%half_digest /= current%half_digest .or. self%population_digest /= current%population_digest ) return
        message = 'reference producer state differs from the selected state'
        if( self%producer_state /= state ) return
        message = ''
        valid = .true.
    end function reference_matches_project

    ! Validate source metadata; destination choices and file operations belong to startup staging.
    logical function reference_derivatives_match( self, params, producer_state, align_lp, message ) result(valid)
        class(refine3D_provenance), intent(in)  :: self
        class(parameters),          intent(in)  :: params
        integer,                    intent(out) :: producer_state
        real,                       intent(out) :: align_lp
        character(len=*),           intent(out) :: message
        valid = .false.
        producer_state = self%producer_state
        align_lp = self%align_lp
        message = self%validation_message
        if( .not. self%valid ) return
        if( self%kind /= 'reference' ) return
        message = 'reference preparation settings differ'
        if( self%settings_digest /= settings_identity(params, 'reference') ) return
        message = 'reference matching LP is unavailable'
        if( .not. ieee_is_finite(align_lp) ) return
        if( align_lp <= 0. ) return
        message = ''
        valid = .true.
    end function reference_derivatives_match

    subroutine publish_reuse_provenance( params, project, particles, kind, anchor, files )
        class(parameters), intent(in) :: params
        class(sp_project), intent(in) :: project
        class(oris),       intent(in) :: particles
        character(len=*),  intent(in) :: kind
        class(string),     intent(in) :: anchor, files(:)
        call write_record(params, project, particles, kind, anchor, files, 0, 0.)
    end subroutine publish_reuse_provenance

    logical function reuse_provenance_matches( params, project, particles, kind, anchor, expected, message ) result(valid)
        class(parameters), intent(in) :: params
        class(sp_project), intent(in) :: project
        class(oris),       intent(in) :: particles
        character(len=*),  intent(in) :: kind
        class(string),     intent(in) :: anchor, expected(:)
        character(len=*), intent(out) :: message
        type(refine3D_provenance) :: record
        call record%new(anchor, kind, message)
        valid = record%reuse_matches(params, project, particles, expected, message)
        call record%kill
    end function reuse_provenance_matches

    logical function reuse_matches( self, params, project, particles, expected, message ) result(valid)
        class(refine3D_provenance), intent(in)  :: self
        class(parameters),          intent(in)  :: params
        class(sp_project),          intent(in)  :: project
        class(oris),                intent(in)  :: particles
        class(string),              intent(in)  :: expected(:)
        character(len=*),           intent(out) :: message
        type(refine3D_provenance) :: current
        integer                   :: i
        valid = .false.
        message = self%validation_message
        if( .not. self%valid ) return
        message = 'handoff constituent list differs'
        if( size(self%files) /= size(expected) ) return
        do i = 1,size(self%files)
            if( basename(self%files(i)) /= basename(expected(i)) ) return
        enddo
        call particle_identity(project, particles, trim(self%kind), current)
        message = 'ordered rows, halves, population or observation identity differs'
        if( self%layout_digest /= current%layout_digest .or. self%half_digest /= current%half_digest .or. &
            &self%population_digest /= current%population_digest .or. &
            &self%observation_digest /= current%observation_digest ) return
        if( self%layout_digest == 0_int64 ) return
        if( self%kind == 'chain' )then
            message = 'chain producer poses differ from the saved project'
            if( self%pose_digest /= current%pose_digest ) return
        endif
        message = 'handoff producer settings differ'
        if( self%settings_digest /= settings_identity(params, trim(self%kind)) ) return
        message = ''
        valid = .true.
    end function reuse_matches

    ! Sigma estimates may lag the current poses; chain imports require the saved pose snapshot.
    subroutine particle_identity( project, particles, kind, identity )
        class(sp_project), intent(in)             :: project
        class(oris),       intent(in)             :: particles
        character(len=*),  intent(in)             :: kind
        class(refine3D_provenance), intent(inout) :: identity
        integer :: i
        identity%layout_digest = ptcl_layout_digest(project, particles)
        identity%half_digest = sigma2_state_digest_begin()
        identity%population_digest = sigma2_state_digest_begin()
        identity%pose_digest = sigma2_state_digest_begin()
        identity%observation_digest = sigma2_state_digest_begin()
        identity%lineage_digest = sigma2_state_digest_begin()
        ! Zero rows retain lineage identity without treating reordered rows as an external project.
        identity%lineage_digest = ptcl_layout_digest(project, particles, nrows=0)
        do i = 1,particles%get_noris()
            call sigma2_state_digest_integer(identity%half_digest, particles%get_eo(i))
            call sigma2_state_digest_integer(identity%population_digest, particles%get_state(i))
            if( kind == 'chain' )then
                call sigma2_state_digest_integer(identity%population_digest, &
                    &merge(1, 0, particles%get_updatecnt(i) > 0))
            endif
            call hash_reals(identity%pose_digest, particles%get_euler(i))
            call hash_reals(identity%pose_digest, particles%get_2Dshift(i))
        enddo
        call hash_observations(identity%observation_digest, particles)
        call hash_observations(identity%observation_digest, project%os_stk)
        call hash_observations(identity%observation_digest, project%os_optics)
    end subroutine particle_identity

    subroutine hash_observations( hash, os )
        integer(int64), intent(inout) :: hash
        class(oris), intent(in)       :: os
        character(len=8), parameter :: KEYS(10) = [character(len=8) :: &
            &'dfx','dfy','angast','phshift','kv','cs','fraca','smpd','w','ogid']
        type(string) :: flag
        real         :: value
        integer      :: i, j
        call sigma2_state_digest_integer(hash, os%get_noris())
        do i = 1,os%get_noris()
            do j = 1,size(KEYS)
                value = 0.
                if( os%isthere(i, trim(KEYS(j))) ) value = os%get(i, trim(KEYS(j)))
                call hash_reals(hash, [value])
            enddo
            flag = ''
            if( os%isthere(i, 'ctf') ) flag = os%get_str(i, 'ctf')
            call sigma2_state_digest_text(hash, flag%to_char())
        enddo
    end subroutine hash_observations

    integer(int64) function settings_identity( params, kind ) result(hash)
        class(parameters), intent(in) :: params
        character(len=*),  intent(in) :: kind
        hash = sigma2_state_digest_begin()
        call sigma2_state_digest_text(hash, params%noise_norm)
        call hash_reals(hash, [params%mskdiam])
        if( kind == 'sigma' ) return
        call sigma2_state_digest_text(hash, params%rec_backend)
        call sigma2_state_digest_text(hash, params%pgrp)
        call sigma2_state_digest_text(hash, params%m_estimator)
        if( kind == 'chain' ) return
        call sigma2_state_digest_text(hash, params%filt_mode)
        call sigma2_state_digest_text(hash, params%ml_reg)
        call sigma2_state_digest_text(hash, params%envfsc)
        call sigma2_state_digest_text(hash, params%automsk)
        call hash_reals(hash, [params%tau, params%amsklp, params%envmsklp, params%nu_msk_sig])
        call sigma2_state_digest_integer(hash, params%binwidth)
        call sigma2_state_digest_integer(hash, params%edge)
    end function settings_identity

    subroutine hash_reals( hash, values )
        integer(int64), intent(inout) :: hash
        real, intent(in)              :: values(:)
        integer :: i
        do i = 1,size(values)
            call sigma2_state_digest_integer(hash, transfer(values(i), 0))
        enddo
    end subroutine hash_reals

    subroutine write_record( params, project, particles, kind, anchor, files, state, align_lp )
        class(parameters), intent(in) :: params
        class(sp_project), intent(in) :: project
        class(oris),       intent(in) :: particles
        character(len=*),  intent(in) :: kind
        class(string),     intent(in) :: anchor, files(:)
        integer,           intent(in) :: state
        real,              intent(in) :: align_lp
        type(string)              :: record, pending, name
        integer(int64)            :: sizes(size(files)), digests(size(files)), settings
        type(refine3D_provenance) :: identity
        integer                   :: unit, status, i
        if( params%oritype /= 'ptcl3D' ) return
        record = provenance_file(anchor)
        pending = provenance_pending_file(anchor)
        call simple_touch(pending)
        call particle_identity(project, particles, kind, identity)
        settings = settings_identity(params, kind)
        do i = 1,size(files)
            sizes(i) = -1_int64
            if( file_exists(files(i)) ) inquire(file=files(i)%to_char(), size=sizes(i))
            digests(i) = sigma2_state_digest_file(files(i))
            if( sizes(i) <= 0 .or. digests(i) == 0_int64 )then
                THROW_HARD('Cannot publish an incomplete refinement handoff')
            endif
        enddo
        call fopen(unit, file=record, status='replace', action='write', iostat=status)
        if( status /= 0 ) THROW_HARD('Cannot open refinement handoff provenance')
        write(unit,*,iostat=status) VERSION, kind, size(files), state, params%which_iter, &
            &params%box_crop, params%smpd_crop, particles%get_noris(), &
            &identity%layout_digest, identity%half_digest, identity%population_digest, identity%pose_digest, &
            &identity%observation_digest, identity%lineage_digest, settings
        name = simple_abspath(params%projfile, check_exists=.false.)
        if( status == 0 ) write(unit,'(A)',iostat=status) name%to_char()
        if( status == 0 ) write(unit,*,iostat=status) trim(params%rec_backend), ' ', trim(params%pgrp), ' ', &
            &trim(params%filt_mode), ' ', params%ml_reg, ' ', params%envfsc, ' ', params%noise_norm, &
            &params%mskdiam, params%tau, params%amsklp, params%envmsklp, params%nu_msk_sig, &
            &params%binwidth, params%edge
        if( status == 0 ) write(unit,*,iostat=status) align_lp
        do i = 1,size(files)
            name = basename(files(i))
            if( status == 0 ) write(unit,'(A)',iostat=status) name%to_char()
            if( status == 0 ) write(unit,*,iostat=status) sizes(i), digests(i)
        enddo
        call fclose(unit)
        if( status /= 0 ) THROW_HARD('Cannot finish refinement handoff provenance')
        call del_file(pending)
    end subroutine write_record

    logical function read_record( self, message ) result(valid)
        class(refine3D_provenance), intent(inout) :: self
        character(len=*),          intent(out)    :: message
        type(string)          :: record, pending, dir
        character(len=STDLEN) :: name, producer, stored_kind, producer_settings
        integer(int64)        :: bytes, digest, actual_size
        integer               :: unit, status, version_read, n, i, state, iteration, box, nrows
        real                  :: smpd
        valid = .false.
        message = 'handoff provenance is absent, incomplete or unreadable'
        record = provenance_file(self%anchor)
        pending = provenance_pending_file(self%anchor)
        if( file_exists(pending) ) return
        if( .not. file_exists(record) ) return
        call fopen(unit, file=record, status='old', action='read', iostat=status)
        if( status /= 0 ) return
        read(unit,*,iostat=status) version_read, stored_kind, n, state, iteration, box, smpd, nrows, &
            &self%layout_digest, self%half_digest, self%population_digest, self%pose_digest, &
            &self%observation_digest, self%lineage_digest, self%settings_digest
        if( status == 0 )then
            if( version_read /= VERSION )then
                message = 'unsupported handoff provenance version; regenerate with the current build'
                status = 1
            endif
            if( trim(stored_kind) /= self%kind .or. n < 1 .or. n > MAX_FILES ) status = 1
            if( box <= 0 .or. nrows <= 0 .or. state < 0 .or. iteration < 0 ) status = 1
            if( .not. ieee_is_finite(smpd) ) status = 1
            if( smpd <= 0. ) status = 1
        endif
        if( status == 0 ) read(unit,'(A)',iostat=status) producer
        if( status == 0 ) read(unit,'(A)',iostat=status) producer_settings
        if( status == 0 ) read(unit,*,iostat=status) self%align_lp
        if( status == 0 ) self%producer_state = state
        if( status == 0 )then
            allocate(self%files(n))
            dir = get_fpath(self%anchor)
            if( .not. dir%is_allocated() ) dir = './'
            do i = 1,n
                read(unit,'(A)',iostat=status) name
                if( status /= 0 ) exit
                if( len_trim(name) == 0 .or. scan(trim(name), '/\') /= 0 )then
                    status = 1
                    exit
                endif
                self%files(i) = dir//trim(name)
                read(unit,*,iostat=status) bytes, digest
                if( status /= 0 ) exit
                actual_size = -1_int64
                if( file_exists(self%files(i)) ) inquire(file=self%files(i)%to_char(), size=actual_size)
                if( bytes <= 0 .or. actual_size /= bytes .or. digest == 0_int64 )then
                    status = 1
                    exit
                endif
                if( sigma2_state_digest_file(self%files(i)) /= digest )then
                    status = 1
                    exit
                endif
            enddo
        endif
        call fclose(unit)
        if( status /= 0 ) return
        if( basename(self%files(1)) /= basename(self%anchor) ) return
        valid = .true.
        message = ''
    end function read_record

end module simple_refine3D_provenance

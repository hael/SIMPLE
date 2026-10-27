!@descr: producer records and validation of local-refinement reference, sigma and accumulator handoffs
module simple_refine3D_provenance
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,              only: STDLEN
use simple_defs_fname,        only: MRC_EXT, NUFILT_SUFFIX, NULOCRES_SUFFIX, NU_ENVMASK_FBODY, AUTOMASK_FBODY
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_string_utils,      only: int2str_pad
use simple_syslib,            only: file_exists, simple_touch, del_file, simple_abspath
use simple_fileio,            only: fopen, fclose, basename, get_fpath, add2fbody, simple_copy_file
use simple_parameters,        only: parameters
use simple_sp_project,        only: sp_project
use simple_oris,              only: oris
use simple_ptcl_layout,       only: ptcl_layout_digest
use simple_halfmap_diagnostics, only: support_provenance_fname
use simple_sigma2_state_file, only: sigma2_state_digest_begin, sigma2_state_digest_integer, &
    &sigma2_state_digest_text, sigma2_state_digest_file
implicit none
private
public :: begin_reference_set, publish_reference_set, validate_reference_set
public :: publish_reuse_provenance, reuse_provenance_matches, provenance_file
public :: local_import_marker, reference_matches_project, stage_reference_derivatives
public :: REFSET_OK, REFSET_LEGACY, REFSET_INVALID
#include "simple_local_flags.inc"

integer, parameter :: VERSION = 2, MAX_FILES = 16
integer, parameter :: REFSET_OK = 0, REFSET_LEGACY = 1, REFSET_INVALID = 2

contains

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

    ! The pending marker distinguishes interrupted publication from a legacy input.
    subroutine begin_reference_set( volume )
        class(string), intent(in) :: volume
        type(string) :: record, pending
        record = provenance_file(volume)
        pending = record//'.pending'
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
        integer :: n, i
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

    integer function validate_reference_set( volume, fsc, support_verified, message ) result(status)
        class(string),     intent(in)    :: volume
        type(string),      intent(inout) :: fsc
        logical,           intent(out)   :: support_verified
        character(len=*),  intent(out)   :: message
        type(string), allocatable :: files(:)
        type(string) :: record, pending, expected
        integer(int64) :: ids(6), settings
        integer :: i, producer_state
        real :: align_lp
        logical :: valid
        call fsc%kill
        support_verified = .false.
        record = provenance_file(volume)
        pending = record//'.pending'
        message = 'reference set publication is incomplete'
        status = REFSET_INVALID
        if( file_exists(pending) ) return
        if( .not. file_exists(record) )then
            status = REFSET_LEGACY
            message = 'legacy base pair: derived-product provenance is unavailable'
            return
        endif
        valid = read_record(volume, 'reference', files, ids, settings, producer_state, align_lp, message)
        if( .not. valid ) return
        if( size(files) < 4 ) return
        expected = add2fbody(volume, MRC_EXT, '_even')
        message = 'reference constituent roles differ'
        if( basename(files(2)) /= basename(expected) ) return
        expected = add2fbody(volume, MRC_EXT, '_odd')
        if( basename(files(3)) /= basename(expected) ) return
        fsc = files(4)
        expected = support_provenance_fname(volume)
        do i = 5,size(files)
            if( basename(files(i)) == basename(expected) ) support_verified = .true.
        enddo
        message = ''
        status = REFSET_OK
    end function validate_reference_set

    ! A reference may lag poses, but its half partition must still match its source rows.
    logical function reference_matches_project( volume, project, state, allow_external, message ) result(valid)
        class(string),     intent(in) :: volume
        class(sp_project), intent(in) :: project
        integer,           intent(in) :: state
        logical,           intent(in) :: allow_external
        character(len=*), intent(out) :: message
        type(string), allocatable :: files(:)
        integer(int64) :: ids(6), current(6), settings
        integer :: producer_state
        real :: align_lp
        valid = read_record(volume, 'reference', files, ids, settings, producer_state, align_lp, message)
        if( .not. valid ) return
        call particle_identity(project, project%os_ptcl3D, 'reference', current)
        valid = .false.
        if( ids(1) == 0_int64 .or. current(1) == 0_int64 )then
            message = 'reference particle lineage is undefined'
            return
        endif
        if( ids(1) /= current(1) )then
            message = 'reference ordered particle lineage differs from the source project'
            if( allow_external .and. ids(6) /= current(6) )then
                message = 'explicit external reference: particle lineage differs; alignment and half relationship are caller supplied'
                valid = .true.
            endif
            return
        endif
        message = 'reference half assignments or active population differ from the source project'
        if( any(ids(2:3) /= current(2:3)) ) return
        message = 'reference producer state differs from the selected state'
        if( producer_state /= state ) return
        message = ''
        valid = .true.
    end function reference_matches_project

    ! Only recorded, compatible derivatives travel; embedded NU support needs no matcher remasking.
    logical function stage_reference_derivatives( params, volume, staged, align_lp, message ) result(reused)
        class(parameters), intent(in) :: params
        class(string),     intent(in) :: volume, staged
        real,             intent(out) :: align_lp
        character(len=*), intent(out) :: message
        type(string), allocatable :: files(:)
        type(string) :: sources(6), destinations(6), bases(3), targets(3), name, dir
        integer(int64) :: ids(6), settings
        integer :: producer_state, i, n
        logical :: valid
        reused = .false.
        align_lp = 0.
        message = 'NU reference bundle is not requested'
        if( .not. params%l_nonuniform ) return
        valid = read_record(volume, 'reference', files, ids, settings, producer_state, align_lp, message)
        if( .not. valid ) return
        message = 'reference preparation settings differ'
        if( settings /= settings_identity(params, 'reference') ) return
        message = 'reference matching LP is unavailable'
        if( .not. ieee_is_finite(align_lp) ) return
        if( align_lp <= 0. ) return
        bases(1) = volume
        bases(2) = add2fbody(volume, MRC_EXT, '_even')
        bases(3) = add2fbody(volume, MRC_EXT, '_odd')
        targets(1) = staged
        targets(2) = add2fbody(staged, MRC_EXT, '_even')
        targets(3) = add2fbody(staged, MRC_EXT, '_odd')
        n = 4
        do i = 1,3
            sources(i) = add2fbody(bases(i), MRC_EXT, NUFILT_SUFFIX)
            destinations(i) = add2fbody(targets(i), MRC_EXT, NUFILT_SUFFIX)
        enddo
        sources(4) = add2fbody(volume, MRC_EXT, NULOCRES_SUFFIX)
        destinations(4) = add2fbody(staged, MRC_EXT, NULOCRES_SUFFIX)
        message = 'NU reference constituents are not recorded'
        do i = 1,n
            if( .not. recorded(sources(i)) ) return
        enddo
        dir = get_fpath(volume)
        if( .not. dir%is_allocated() ) dir = './'
        name = dir//AUTOMASK_FBODY//int2str_pad(producer_state,2)//MRC_EXT
        if( recorded(name) )then
            n = n + 1
            sources(n) = name
            destinations(n) = AUTOMASK_FBODY//'01'//MRC_EXT
        else
            name = AUTOMASK_FBODY//'01'//MRC_EXT
            if( file_exists(name) ) call del_file(name)
        endif
        name = dir//NU_ENVMASK_FBODY//int2str_pad(producer_state,2)//MRC_EXT
        if( recorded(name) )then
            n = n + 1
            sources(n) = name
            destinations(n) = NU_ENVMASK_FBODY//'01'//MRC_EXT
        else
            name = NU_ENVMASK_FBODY//'01'//MRC_EXT
            if( file_exists(name) ) call del_file(name)
        endif
        do i = 1,n
            if( simple_abspath(sources(i)) /= simple_abspath(destinations(i), check_exists=.false.) ) &
                &call simple_copy_file(sources(i), destinations(i))
        enddo
        message = ''
        reused = .true.
    contains
        logical function recorded( file ) result(found)
            class(string), intent(in) :: file
            integer :: j
            found = .false.
            do j = 1,size(files)
                if( basename(files(j)) /= basename(file) ) cycle
                found = .true.
                return
            enddo
        end function recorded
    end function stage_reference_derivatives

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
        type(string), allocatable :: files(:)
        integer(int64) :: ids(6), current(6), settings
        integer :: i, producer_state
        real :: align_lp
        valid = read_record(anchor, kind, files, ids, settings, producer_state, align_lp, message)
        if( .not. valid ) return
        valid = .false.
        message = 'handoff constituent list differs'
        if( size(files) /= size(expected) ) return
        do i = 1,size(files)
            if( basename(files(i)) /= basename(expected(i)) ) return
        enddo
        call particle_identity(project, particles, kind, current)
        message = 'ordered rows, halves, population or observation identity differs'
        if( any(ids(1:3) /= current(1:3)) .or. ids(5) /= current(5) ) return
        if( ids(1) == 0_int64 ) return
        if( kind == 'chain' )then
            message = 'chain producer poses differ from the saved project'
            if( ids(4) /= current(4) ) return
        endif
        message = 'handoff producer settings differ'
        if( settings /= settings_identity(params, kind) ) return
        message = ''
        valid = .true.
    end function reuse_provenance_matches

    ! Sigma estimates may lag the current poses; chain imports require the saved pose snapshot.
    subroutine particle_identity( project, particles, kind, ids )
        class(sp_project), intent(in) :: project
        class(oris),       intent(in) :: particles
        character(len=*),  intent(in) :: kind
        integer(int64),   intent(out) :: ids(6)
        integer :: i
        ids(1) = ptcl_layout_digest(project, particles)
        ids(2:) = sigma2_state_digest_begin()
        ! Zero rows retain lineage identity without treating reordered rows as an external project.
        ids(6) = ptcl_layout_digest(project, particles, nrows=0)
        do i = 1,particles%get_noris()
            call sigma2_state_digest_integer(ids(2), particles%get_eo(i))
            call sigma2_state_digest_integer(ids(3), particles%get_state(i))
            if( kind == 'chain' ) call sigma2_state_digest_integer(ids(3), merge(1, 0, particles%get_updatecnt(i) > 0))
            call hash_reals(ids(4), particles%get_euler(i))
            call hash_reals(ids(4), particles%get_2Dshift(i))
        enddo
        call hash_observations(ids(5), particles)
        call hash_observations(ids(5), project%os_stk)
        call hash_observations(ids(5), project%os_optics)
    end subroutine particle_identity

    subroutine hash_observations( hash, os )
        integer(int64), intent(inout) :: hash
        class(oris), intent(in) :: os
        character(len=8), parameter :: KEYS(10) = [character(len=8) :: &
            &'dfx','dfy','angast','phshift','kv','cs','fraca','smpd','w','ogid']
        type(string) :: flag
        real :: value
        integer :: i, j
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
        real, intent(in) :: values(:)
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
        type(string) :: record, pending, name
        integer(int64) :: ids(6), sizes(size(files)), digests(size(files)), settings
        integer :: unit, status, i
        if( params%oritype /= 'ptcl3D' ) return
        record = provenance_file(anchor)
        pending = record//'.pending'
        call simple_touch(pending)
        call particle_identity(project, particles, kind, ids)
        settings = settings_identity(params, kind)
        do i = 1,size(files)
            sizes(i) = -1_int64
            if( file_exists(files(i)) ) inquire(file=files(i)%to_char(), size=sizes(i))
            digests(i) = sigma2_state_digest_file(files(i))
            if( sizes(i) <= 0 .or. digests(i) == 0_int64 ) THROW_HARD('Cannot publish an incomplete refinement handoff')
        enddo
        call fopen(unit, file=record, status='replace', action='write', iostat=status)
        if( status /= 0 ) THROW_HARD('Cannot open refinement handoff provenance')
        write(unit,*,iostat=status) VERSION, kind, size(files), state, params%which_iter, &
            &params%box_crop, params%smpd_crop, particles%get_noris(), ids, settings
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

    logical function read_record( anchor, kind, files, ids, settings, producer_state, align_lp, message ) result(valid)
        class(string),    intent(in)  :: anchor
        character(len=*), intent(in)  :: kind
        type(string), allocatable, intent(out) :: files(:)
        integer(int64), intent(out) :: ids(6), settings
        integer, intent(out) :: producer_state
        real, intent(out) :: align_lp
        character(len=*), intent(out) :: message
        type(string) :: record, pending, dir
        character(len=STDLEN) :: name, producer, stored_kind, producer_settings
        integer(int64) :: bytes, digest, actual_size
        integer :: unit, status, version_read, n, i, state, iteration, box, nrows
        real :: smpd
        valid = .false.
        producer_state = 0
        align_lp = 0.
        ids = 0_int64
        settings = 0_int64
        message = 'handoff provenance is absent, incomplete or unreadable'
        record = provenance_file(anchor)
        pending = record//'.pending'
        if( file_exists(pending) ) return
        if( .not. file_exists(record) ) return
        call fopen(unit, file=record, status='old', action='read', iostat=status)
        if( status /= 0 ) return
        read(unit,*,iostat=status) version_read, stored_kind, n, state, iteration, box, smpd, nrows, ids, settings
        if( status == 0 )then
            if( version_read /= VERSION )then
                message = 'unsupported handoff provenance version; regenerate with the current build'
                status = 1
            endif
            if( trim(stored_kind) /= kind .or. n < 1 .or. n > MAX_FILES ) status = 1
            if( box <= 0 .or. nrows <= 0 .or. state < 0 .or. iteration < 0 ) status = 1
            if( .not. ieee_is_finite(smpd) ) status = 1
            if( smpd <= 0. ) status = 1
        endif
        if( status == 0 ) read(unit,'(A)',iostat=status) producer
        if( status == 0 ) read(unit,'(A)',iostat=status) producer_settings
        if( status == 0 ) read(unit,*,iostat=status) align_lp
        if( status == 0 ) producer_state = state
        if( status == 0 )then
            allocate(files(n))
            dir = get_fpath(anchor)
            if( .not. dir%is_allocated() ) dir = './'
            do i = 1,n
                read(unit,'(A)',iostat=status) name
                if( status /= 0 ) exit
                if( len_trim(name) == 0 .or. scan(trim(name), '/\') /= 0 )then
                    status = 1
                    exit
                endif
                files(i) = dir//trim(name)
                read(unit,*,iostat=status) bytes, digest
                if( status /= 0 ) exit
                actual_size = -1_int64
                if( file_exists(files(i)) ) inquire(file=files(i)%to_char(), size=actual_size)
                if( bytes <= 0 .or. actual_size /= bytes .or. digest == 0_int64 )then
                    status = 1
                    exit
                endif
                if( sigma2_state_digest_file(files(i)) /= digest )then
                    status = 1
                    exit
                endif
            enddo
        endif
        call fclose(unit)
        if( status /= 0 ) return
        if( basename(files(1)) /= basename(anchor) ) return
        valid = .true.
        message = ''
    end function read_record

end module simple_refine3D_provenance

!@descr: the typed, versioned run manifest of solve3D: identity, solution, ladder, inputs and artifact digests of a completed run
! Written last, atomically and never fatally by exec_solve3D (write_run_manifest),
! then registered in projinfo by bare name. It is the only route into solve3D_addon:
! read_registered, validate_frozen against the frozen project, then replay. The replayed inputs are
! the add-on's defaults; the overridable ones (balance, nclust, mskdiam: how the cohort is sampled and
! masked, nothing the frozen term was built with) yield to the add-on's own command line.
! Plain key-value text with a schema line, a completion status and an FNV-1a checksum over
! every preceding line; read refuses unknown fields, a bad checksum, truncation and other versions.
module simple_solve3D_manifest
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_cmdline,           only: cmdline
use simple_sp_project,        only: sp_project
use simple_sigma2_state,      only: sigma2_state_project_layout_digest
use simple_sigma2_state_file, only: sigma2_state_digest_begin, sigma2_state_digest_text, &
    &sigma2_state_digest_integer, sigma2_state_digest_file
implicit none

public :: solve3D_manifest, solve3D_stage_record, MANIFEST_FNAME, manifest_records_input, manifest_overridable_input
private
#include "simple_local_flags.inc"

character(len=*), parameter :: MANIFEST_SCHEMA       = 'solve3D_manifest'
integer,          parameter :: MANIFEST_VERSION      = 1
character(len=*), parameter :: MANIFEST_FNAME        = 'solve3D_manifest.txt'
character(len=*), parameter :: MANIFEST_PROJINFO_KEY = 'solve3D_manifest'
character(len=*), parameter :: RUN_ID_PROJINFO_KEY   = 'solve3D_run_id'
integer,          parameter :: KLEN = 24

!> Keys of the base run's command line, as given at entry (before any default
!! is injected), that the manifest records: the solution, reconstruction and
!! search policy plus the entry routes (provenance)
character(len=KLEN), parameter :: MANIFEST_INPUT_KEYS(33) = [character(len=KLEN) :: &
    &'pgrp', 'mskdiam', 'nstates', 'nstages', 'rec_backend', 'maxits_pcg', &
    &'maxits_ml', 'pcg_solvent', 'pcg_solvent_lambda', 'filt_mode', 'automsk', 'envfsc', 'envmsklp', &
    &'objfun', 'sigma_est', 'hp', 'lp', 'lpstart', 'lpstop', 'force_lp_range', &
    &'prob_athres', 'bfac', 'gauref', 'balance', 'nclust', &
    &'lpstart_ini3D', 'lpstop_ini3D', 'center', 'cenlp', 'cavg_ini', 'cavg_ini_ext', 'pgrp_start', 'vol1']

!> The subset solve3D_addon replays as given: everything that describes the
!! model the cohort is aligned to, plus the overridable keys below as the
!! add-on's defaults. Entry routes and their controls, the state
!! layout (derived from the completed solution), the stage range
!! (from the ladder), centring (forced off) and the compute/convergence keys
!! the add-on accepts from its own command line are not replayed.
character(len=KLEN), parameter :: MANIFEST_REPLAY_KEYS(25) = [character(len=KLEN) :: &
    &'pgrp', 'mskdiam', 'rec_backend', 'maxits_pcg', 'maxits_ml', 'pcg_solvent', 'pcg_solvent_lambda', &
    &'filt_mode', 'automsk', 'envfsc', 'envmsklp', 'objfun', 'sigma_est', 'hp', &
    &'lp', 'lpstart', 'lpstop', 'force_lp_range', &
    &'prob_athres', 'bfac', 'gauref', 'balance', 'nclust', 'lpstart_ini3D', 'lpstop_ini3D']

!> Replayed keys the add-on may override from its own command line: they
!! govern how the cohort is sampled (balance, nclust) and masked (mskdiam),
!! which depends on the particles being added, not on the frozen solution.
!! The frozen accumulators are mask-free raw sums; the mask enters at
!! restoration, in the solve support and in matching, all on the union.
character(len=KLEN), parameter :: MANIFEST_OVERRIDABLE_KEYS(3) = [character(len=KLEN) :: &
    &'balance', 'nclust', 'mskdiam']

!> one stage of the ladder: the planned record and the limits actually
!! emitted (0 = not on the stage line, -1 = a stage the run never ran)
type :: solve3D_stage_record
    real    :: lp_planned = 0., lp_emitted = 0., lpstop_emitted = 0.
    real    :: smpd_crop = 0., scale = 1., trslim = 0., frc_crit = 0.
    integer :: box_crop = 0
    logical :: l_autoscale = .false., l_lpset = .false.
end type solve3D_stage_record

!> the manifest of one completed solve3D or solve3D_addon run
type :: solve3D_manifest
    private
    type(string)       :: fname              !< the file it was read from
    character(len=64)  :: run_id  = ''
    character(len=32)  :: program = ''
    logical            :: complete = .false.
    logical            :: eligible = .false. !< may serve as a frozen input of solve3D_addon
    ! project and particle-layout identity
    integer            :: nrows = 0
    integer(int64)     :: layout_digest = 0_int64, stack_digest = 0_int64, optics_digest = 0_int64
    ! the solution
    integer            :: nstates = 0, box = 0
    real               :: smpd = 0., mskdiam = 0.
    character(len=16)  :: pgrp = ''
    ! effective and provenance values of the base population
    integer            :: nsample = 0, nptcls_eff = 0
    real               :: update_frac = 1.
    logical            :: full_sampling = .false.
    ! stage command-line shape at planning time
    logical            :: lp_on_line = .false., lpstop_on_line = .false.
    real               :: lp_line = 0., lpstop_line = 0.
    ! the ladder
    integer            :: first_stage = 0, last_stage = 0
    type(solve3D_stage_record), allocatable :: stages(:)
    ! inputs as given
    character(len=KLEN), allocatable :: input_keys(:)
    type(string),        allocatable :: input_vals(:)
    ! artifacts: kind, state, digest, path
    character(len=16),   allocatable :: artifact_kind(:)
    integer,             allocatable :: artifact_state(:)
    integer(int64),      allocatable :: artifact_digest(:)
    type(string),        allocatable :: artifact_path(:)
    ! a record the text format cannot hold; write refuses to publish
    character(len=STDLEN) :: defect = ''
  contains
    ! construction by the completed run
    procedure          :: new
    procedure          :: set_solution
    procedure          :: set_sampling
    procedure          :: set_stage_line
    procedure          :: set_ladder
    procedure          :: record_inputs
    procedure          :: record_artifacts
    procedure, private :: add_input
    procedure, private :: add_artifact
    ! persistence and project registration
    procedure          :: write
    procedure          :: read
    procedure          :: register
    procedure          :: read_registered
    procedure, private :: build_lines
    ! use as a frozen input
    procedure          :: validate_frozen
    procedure          :: replay
    procedure          :: get_artifact
    procedure          :: matches_artifact
    procedure, private :: get_input
    procedure, private :: find_artifact
    ! getters
    procedure          :: get_fname
    procedure          :: get_run_id
    procedure          :: get_nstates
    procedure          :: get_box
    procedure          :: get_smpd
    procedure          :: get_last_stage
    procedure          :: get_nstages
    procedure          :: get_stage
    procedure          :: kill
end type solve3D_manifest

contains

    !> key is one of the run inputs a manifest records (MANIFEST_INPUT_KEYS):
    !! part of a completed run's settings, which solve3D_addon inherits
    logical function manifest_records_input( key ) result( l_recorded )
        character(len=*), intent(in) :: key
        l_recorded = any(MANIFEST_INPUT_KEYS == key)
    end function manifest_records_input

    !> key is a replayed input solve3D_addon may override from its own
    !! command line (MANIFEST_OVERRIDABLE_KEYS); the replayed value is its default
    logical function manifest_overridable_input( key ) result( l_overridable )
        character(len=*), intent(in) :: key
        l_overridable = any(MANIFEST_OVERRIDABLE_KEYS == key)
    end function manifest_overridable_input

    ! CONSTRUCTION

    !> A completed run of program_name over the particle layout of spproj:
    !! the identity of its rows, its stack table and its optics/CTF
    !! parameters. The records follow:
    !! set_solution (before record_artifacts), set_sampling, set_stage_line,
    !! set_ladder, record_inputs.
    subroutine new( self, run_id, program_name, eligible, spproj )
        class(solve3D_manifest), intent(inout) :: self
        character(len=*),           intent(in)    :: run_id, program_name
        logical,                    intent(in)    :: eligible
        class(sp_project),          intent(inout) :: spproj
        call self%kill
        self%run_id        = run_id
        self%program       = program_name
        self%complete      = .true.
        self%eligible      = eligible
        self%nrows         = spproj%os_ptcl3D%get_noris()
        self%layout_digest = sigma2_state_project_layout_digest(spproj, spproj%os_ptcl3D)
        self%stack_digest  = stack_table_digest(spproj)
        self%optics_digest = optics_ctf_digest(spproj)
    end subroutine new

    !> the completed solution: state layout, point group, native grid and mask
    subroutine set_solution( self, nstates, pgrp, box, smpd, mskdiam )
        class(solve3D_manifest), intent(inout) :: self
        integer,                    intent(in)    :: nstates, box
        character(len=*),           intent(in)    :: pgrp
        real,                       intent(in)    :: smpd, mskdiam
        self%nstates            = nstates
        self%pgrp               = trim(pgrp)
        self%box                = box
        self%smpd               = smpd
        self%mskdiam            = mskdiam
    end subroutine set_solution

    !> the effective sampling of the base population
    subroutine set_sampling( self, nsample, nptcls_eff, update_frac, full_sampling )
        class(solve3D_manifest), intent(inout) :: self
        integer,                    intent(in)    :: nsample, nptcls_eff
        real,                       intent(in)    :: update_frac
        logical,                    intent(in)    :: full_sampling
        self%nsample       = nsample
        self%nptcls_eff    = nptcls_eff
        self%update_frac   = update_frac
        self%full_sampling = full_sampling
    end subroutine set_sampling

    !> the stage command line's shape at planning time: which limits were on it
    subroutine set_stage_line( self, lp_on_line, lp, lpstop_on_line, lpstop )
        class(solve3D_manifest), intent(inout) :: self
        logical,                    intent(in)    :: lp_on_line, lpstop_on_line
        real,                       intent(in)    :: lp, lpstop
        self%lp_on_line     = lp_on_line
        self%lp_line        = lp
        self%lpstop_on_line = lpstop_on_line
        self%lpstop_line    = lpstop
    end subroutine set_stage_line

    !> the ladder as planned and emitted, and the stage range the run executed
    subroutine set_ladder( self, first_stage, last_stage, stages )
        class(solve3D_manifest),       intent(inout) :: self
        integer,                       intent(in)    :: first_stage, last_stage
        type(solve3D_stage_record), intent(in)       :: stages(:)
        self%first_stage = first_stage
        self%last_stage  = last_stage
        self%stages      = stages
    end subroutine set_ladder

    !> The MANIFEST_INPUT_KEYS given on the run's entry command line, as
    !! given. A value the text format cannot hold (empty or with blanks) makes
    !! the manifest defective: write then refuses to publish it, and the
    !! completed run is left as it is.
    subroutine record_inputs( self, cline )
        class(solve3D_manifest), intent(inout) :: self
        class(cmdline),             intent(in)    :: cline
        type(chash)  :: descr
        type(string) :: val
        character(len=:), allocatable :: key, raw
        integer :: i
        call cline%gen_job_descr(descr)
        do i = 1, size(MANIFEST_INPUT_KEYS)
            key = trim(MANIFEST_INPUT_KEYS(i))
            if( .not. descr%isthere(key) ) cycle
            val = descr%get(key)
            raw = trim(val%to_char())
            if( len(raw) == 0 .or. index(raw, ' ') > 0 )then
                if( len_trim(self%defect) == 0 ) self%defect = 'input '//key//' cannot be recorded (empty or with blanks)'
                cycle
            endif
            call self%add_input(key, raw)
        enddo
        call descr%kill
        call val%kill
    end subroutine record_inputs

    !> the run's registered products with their digests: every state's final
    !! map, its halves and FSC, and the committed residual sigma2 state
    subroutine record_artifacts( self, spproj )
        class(solve3D_manifest), intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        type(string) :: fname, halfvol
        integer :: state, box_vol
        real    :: smpd_vol
        logical :: found
        do state = 1, self%nstates
            if( .not. spproj%isthere_in_osout('vol', state) ) cycle
            call spproj%get_vol('vol', state, fname, smpd_vol, box_vol)
            call self%add_artifact('vol', state, sigma2_state_digest_file(fname), fname)
            halfvol = add2fbody(fname, MRC_EXT, '_even')
            if( file_exists(halfvol) ) &
                &call self%add_artifact('vol_even', state, sigma2_state_digest_file(halfvol), halfvol)
            halfvol = add2fbody(fname, MRC_EXT, '_odd')
            if( file_exists(halfvol) ) &
                &call self%add_artifact('vol_odd', state, sigma2_state_digest_file(halfvol), halfvol)
            if( spproj%isthere_in_osout('fsc', state) )then
                call spproj%get_fsc(state, fname, box_vol)
                call self%add_artifact('fsc', state, sigma2_state_digest_file(fname), fname)
            endif
        enddo
        call spproj%get_sigma2_state_path(fname, found)
        if( found .and. file_exists(fname) ) &
            &call self%add_artifact('sigma2_state', 0, sigma2_state_digest_file(fname), fname)
        call fname%kill
        call halfvol%kill
    end subroutine record_artifacts

    subroutine add_input( self, key, val )
        class(solve3D_manifest), intent(inout) :: self
        character(len=*),           intent(in)    :: key, val
        character(len=KLEN), allocatable :: keys(:)
        type(string),        allocatable :: vals(:)
        integer                          :: n
        if( .not. allocated(self%input_keys) )then
            allocate(self%input_keys(0), self%input_vals(0))
        endif
        n = size(self%input_keys)
        allocate(keys(n+1), vals(n+1))
        keys(1:n) = self%input_keys
        vals(1:n) = self%input_vals
        keys(n+1) = key
        vals(n+1) = val
        call move_alloc(keys, self%input_keys)
        call move_alloc(vals, self%input_vals)
    end subroutine add_input

    subroutine add_artifact( self, kind, state, digest, path )
        class(solve3D_manifest), intent(inout) :: self
        character(len=*),           intent(in)    :: kind
        integer,                    intent(in)    :: state
        integer(int64),             intent(in)    :: digest
        class(string),              intent(in)    :: path
        character(len=16), allocatable :: kinds(:)
        integer,           allocatable :: states(:)
        integer(int64),    allocatable :: digests(:)
        type(string),      allocatable :: paths(:)
        integer                        :: n
        if( .not. allocated(self%artifact_kind) )then
            allocate(self%artifact_kind(0), self%artifact_state(0), self%artifact_digest(0), self%artifact_path(0))
        endif
        n = size(self%artifact_kind)
        allocate(kinds(n+1), states(n+1), digests(n+1), paths(n+1))
        kinds(1:n) = self%artifact_kind;   kinds(n+1)   = kind
        states(1:n) = self%artifact_state; states(n+1)  = state
        digests(1:n) = self%artifact_digest; digests(n+1) = digest
        paths(1:n) = self%artifact_path;   paths(n+1)   = path
        call move_alloc(kinds,   self%artifact_kind)
        call move_alloc(states,  self%artifact_state)
        call move_alloc(digests, self%artifact_digest)
        call move_alloc(paths,   self%artifact_path)
    end subroutine add_artifact

    ! WRITER

    !> Publish the manifest atomically (.tmp then rename). Never fatal: returns
    !! status /= 0 with a message, and a failed write leaves no manifest.
    subroutine write( self, fname, status, msg )
        class(solve3D_manifest), intent(in)     :: self
        class(string),              intent(in)  :: fname
        integer,                    intent(out) :: status
        character(len=*),           intent(out) :: msg
        type(string), allocatable :: lines(:)
        type(string)              :: tmpname
        integer(int64)            :: checksum
        integer                   :: funit, io_stat, i, n
        status = 1
        msg    = ''
        if( len_trim(self%defect) > 0 )then
            msg = trim(self%defect)
            return
        endif
        call self%build_lines(lines)
        checksum = sigma2_state_digest_begin()
        do i = 1, size(lines)
            call sigma2_state_digest_text(checksum, lines(i)%to_char())
        enddo
        tmpname = fname//'.tmp'
        call del_file(tmpname)
        call fopen(funit, file=tmpname, status='REPLACE', action='WRITE', iostat=io_stat)
        if( io_stat /= 0 )then
            msg = 'cannot open '//tmpname%to_char()
            return
        endif
        n = size(lines)
        do i = 1, n
            write(funit,'(A)',iostat=io_stat) lines(i)%to_char()
            if( io_stat /= 0 ) exit
        enddo
        if( io_stat == 0 ) write(funit,'(A,1X,A)',iostat=io_stat) 'checksum', trim(int64_str(checksum))
        call fclose(funit)
        if( io_stat /= 0 )then
            call del_file(tmpname)
            msg = 'failed writing '//tmpname%to_char()
            return
        endif
        call simple_rename(tmpname, fname, overwrite=.true.)
        if( .not. file_exists(fname) )then
            msg = 'failed publishing '//fname%to_char()
            return
        endif
        status = 0
        call tmpname%kill
    end subroutine write

    !> the records of the manifest
    subroutine build_lines( self, lines )
        class(solve3D_manifest), intent(in)     :: self
        type(string), allocatable,  intent(out) :: lines(:)
        type(string), allocatable     :: buf(:)
        integer                       :: n, i
        type(solve3D_stage_record) :: st
        n = 0
        allocate(buf(512))
        call push(MANIFEST_SCHEMA//' '//int2str(MANIFEST_VERSION))
        call push('run_id '//trim(self%run_id))
        call push('program '//trim(self%program))
        call push('status '//trim(merge('complete  ', 'incomplete', self%complete)))
        call push('eligible '//trim(logical_str(self%eligible)))
        call push('nrows '//int2str(self%nrows))
        call push('layout_digest '//trim(int64_str(self%layout_digest)))
        call push('stack_digest '//trim(int64_str(self%stack_digest)))
        call push('optics_digest '//trim(int64_str(self%optics_digest)))
        call push('nstates '//int2str(self%nstates))
        call push('pgrp '//trim(self%pgrp))
        call push('box '//int2str(self%box))
        call push('smpd '//trim(fmt_real(self%smpd)))
        call push('mskdiam '//trim(fmt_real(self%mskdiam)))
        call push('nsample '//int2str(self%nsample))
        call push('nptcls_eff '//int2str(self%nptcls_eff))
        call push('update_frac '//trim(fmt_real(self%update_frac)))
        call push('full_sampling '//trim(logical_str(self%full_sampling)))
        call push('lp_line '//trim(logical_str(self%lp_on_line))//' '//trim(fmt_real(self%lp_line)))
        call push('lpstop_line '//trim(logical_str(self%lpstop_on_line))//' '//trim(fmt_real(self%lpstop_line)))
        call push('first_stage '//int2str(self%first_stage))
        call push('last_stage '//int2str(self%last_stage))
        if( allocated(self%stages) )then
            call push('nstages_ladder '//int2str(size(self%stages)))
            do i = 1, size(self%stages)
                st = self%stages(i)
                call push('stage '//int2str(i)//' '//trim(fmt_real(st%lp_planned))//' '//trim(fmt_real(st%lp_emitted))// &
                    &' '//trim(fmt_real(st%lpstop_emitted))//' '//int2str(st%box_crop)//' '//trim(fmt_real(st%smpd_crop))// &
                    &' '//trim(fmt_real(st%scale))//' '//trim(fmt_real(st%trslim))//' '//trim(fmt_real(st%frc_crit))// &
                    &' '//trim(logical_str(st%l_autoscale))//' '//trim(logical_str(st%l_lpset)))
            enddo
        else
            call push('nstages_ladder 0')
        endif
        if( allocated(self%input_keys) )then
            do i = 1, size(self%input_keys)
                call push('input '//trim(self%input_keys(i))//' '//self%input_vals(i)%to_char())
            enddo
        endif
        if( allocated(self%artifact_kind) )then
            do i = 1, size(self%artifact_kind)
                call push('artifact '//trim(self%artifact_kind(i))//' '//int2str(self%artifact_state(i))//' '// &
                    &trim(int64_str(self%artifact_digest(i)))//' '//self%artifact_path(i)%to_char())
            enddo
        endif
        call push('end')
        allocate(lines(n), source=buf(1:n))
        deallocate(buf)

    contains

        subroutine push( line )
            character(len=*), intent(in) :: line
            type(string), allocatable :: grown(:)
            if( n >= size(buf) )then
                allocate(grown(2*size(buf)))
                grown(1:n) = buf(1:n)
                call move_alloc(grown, buf)
            endif
            n = n + 1
            buf(n) = trim(line)
        end subroutine push

    end subroutine build_lines

    ! READER

    !> Parse and verify a manifest: schema and version, checksum over every
    !! preceding line, completion marker, no unknown field. Validation against
    !! a project is validate_frozen.
    subroutine read( self, fname, status, msg )
        class(solve3D_manifest), intent(inout) :: self
        class(string),              intent(in)    :: fname
        integer,                    intent(out)   :: status
        character(len=*),           intent(out)   :: msg
        type(string)                  :: rec
        character(len=:), allocatable :: line, rest, val
        character(len=64)             :: key, word, word2
        character(len=24)             :: flag
        integer(int64)                :: checksum, checksum_file
        integer                       :: funit, io_stat, io_stat_tail, version, istage, nstages, state, pos
        logical                       :: l_end, l_checksum
        call self%kill
        status     = 1
        msg        = ''
        l_end      = .false.
        l_checksum = .false.
        nstages    = -1
        if( .not. file_exists(fname) )then
            msg = 'solve3D manifest is missing: '//fname%to_char()
            return
        endif
        call fopen(funit, file=fname, status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 )then
            msg = 'solve3D manifest is unreadable'
            return
        endif
        checksum = sigma2_state_digest_begin()
        call rec%readline(funit, io_stat)
        if( io_stat == 0 )then
            line = rec%to_char()
            read(line,*,iostat=io_stat) word, version
        endif
        if( io_stat /= 0 .or. trim(word) /= MANIFEST_SCHEMA )then
            msg = 'not a solve3D manifest'
            call fclose(funit)
            return
        endif
        if( version /= MANIFEST_VERSION )then
            msg = 'unsupported solve3D manifest schema version'
            call fclose(funit)
            return
        endif
        call sigma2_state_digest_text(checksum, trim(line))
        do
            call rec%readline(funit, io_stat)
            if( io_stat /= 0 ) exit
            line = rec%to_char()
            if( len_trim(line) == 0 ) cycle
            read(line,*,iostat=io_stat) key
            if( io_stat /= 0 ) exit
            pos  = index(line, trim(key)) + len_trim(key)
            rest = adjustl(line(pos:))
            if( trim(key) == 'checksum' )then
                read(rest,*,iostat=io_stat) checksum_file
                l_checksum = io_stat == 0
                exit
            endif
            if( l_end )then
                msg = 'solve3D manifest has records after its end marker'
                io_stat = 0
                exit
            endif
            call sigma2_state_digest_text(checksum, trim(line))
            select case(trim(key))
                case('run_id');        read(rest,*,iostat=io_stat) self%run_id
                case('program');       read(rest,*,iostat=io_stat) self%program
                case('status')
                    read(rest,*,iostat=io_stat) word
                    self%complete = trim(word) == 'complete'
                case('eligible');      read(rest,*,iostat=io_stat) flag; self%eligible = trim(flag) == 'yes'
                case('nrows');         read(rest,*,iostat=io_stat) self%nrows
                case('layout_digest'); read(rest,*,iostat=io_stat) self%layout_digest
                case('stack_digest');  read(rest,*,iostat=io_stat) self%stack_digest
                case('optics_digest'); read(rest,*,iostat=io_stat) self%optics_digest
                case('nstates');       read(rest,*,iostat=io_stat) self%nstates
                case('pgrp');          read(rest,*,iostat=io_stat) self%pgrp
                case('box');           read(rest,*,iostat=io_stat) self%box
                case('smpd');          read(rest,*,iostat=io_stat) self%smpd
                case('mskdiam');       read(rest,*,iostat=io_stat) self%mskdiam
                case('nsample');       read(rest,*,iostat=io_stat) self%nsample
                case('nptcls_eff');    read(rest,*,iostat=io_stat) self%nptcls_eff
                case('update_frac');   read(rest,*,iostat=io_stat) self%update_frac
                case('full_sampling'); read(rest,*,iostat=io_stat) flag; self%full_sampling = trim(flag) == 'yes'
                case('lp_line')
                    read(rest,*,iostat=io_stat) flag, self%lp_line
                    self%lp_on_line = trim(flag) == 'yes'
                case('lpstop_line')
                    read(rest,*,iostat=io_stat) flag, self%lpstop_line
                    self%lpstop_on_line = trim(flag) == 'yes'
                case('first_stage');   read(rest,*,iostat=io_stat) self%first_stage
                case('last_stage');    read(rest,*,iostat=io_stat) self%last_stage
                case('nstages_ladder')
                    read(rest,*,iostat=io_stat) nstages
                    if( io_stat == 0 .and. nstages >= 0 .and. nstages <= 64 )then
                        allocate(self%stages(nstages))
                    else
                        io_stat = 1
                    endif
                case('stage')
                    read(rest,*,iostat=io_stat) istage
                    if( io_stat == 0 )then
                        if( .not. allocated(self%stages) )then
                            io_stat = 1
                        else if( istage < 1 .or. istage > size(self%stages) )then
                            io_stat = 1
                        else
                            associate(st => self%stages(istage))
                            read(rest,*,iostat=io_stat) istage, st%lp_planned, st%lp_emitted, st%lpstop_emitted, &
                                &st%box_crop, st%smpd_crop, st%scale, st%trslim, st%frc_crit, word, word2
                            st%l_autoscale = trim(word)  == 'yes'
                            st%l_lpset     = trim(word2) == 'yes'
                            end associate
                        endif
                    endif
                case('input')
                    ! the value is the rest of the record after its key: a
                    ! list-directed read would end a path at its first '/'
                    read(rest,*,iostat=io_stat) word
                    if( io_stat == 0 )then
                        if( .not. any(MANIFEST_INPUT_KEYS == word) )then
                            msg = 'unknown solve3D manifest input key: '//trim(word)
                            call fclose(funit)
                            call self%kill
                            return
                        endif
                        pos = index(rest, trim(word)) + len_trim(word)
                        val = adjustl(rest(pos:))
                        if( len_trim(val) == 0 .or. index(trim(val), ' ') > 0 )then
                            io_stat = 1
                        else
                            call self%add_input(trim(word), trim(val))
                        endif
                    endif
                case('artifact')
                    read(rest,*,iostat=io_stat) word, state, checksum_file
                    if( io_stat == 0 )then
                        ! the path is the rest of the record after the digest
                        pos = index(rest, trim(int64_str(checksum_file))) + len_trim(int64_str(checksum_file))
                        call self%add_artifact(trim(word), state, checksum_file, string(trim(adjustl(rest(pos:)))))
                    endif
                case('end')
                    l_end = .true.
                case DEFAULT
                    msg = 'unknown solve3D manifest field: '//trim(key)
                    call fclose(funit)
                    call self%kill
                    return
            end select
            if( io_stat /= 0 ) exit
        enddo
        if( l_checksum .and. len_trim(msg) == 0 )then
            ! the checksum line is the last record
            do
                call rec%readline(funit, io_stat_tail)
                if( io_stat_tail /= 0 ) exit
                if( rec%strlen_trim() > 0 )then
                    msg = 'solve3D manifest has records after its checksum'
                    exit
                endif
            enddo
        endif
        call fclose(funit)
        if( len_trim(msg) > 0 )then
            ! already set
        else if( io_stat /= 0 .and. .not. l_checksum )then
            if( l_end )then
                msg = 'solve3D manifest has no checksum'
            else
                msg = 'corrupt or truncated solve3D manifest'
            endif
        else if( .not. l_end )then
            msg = 'truncated solve3D manifest (no end marker)'
        else if( .not. l_checksum )then
            msg = 'solve3D manifest has no checksum'
        else if( checksum_file /= checksum )then
            msg = 'solve3D manifest checksum mismatch'
        else if( .not. self%complete )then
            msg = 'solve3D manifest does not describe a completed run'
        else if( len_trim(self%run_id) == 0 .or. self%nrows < 1 .or. self%nstates < 1 .or. self%box < 1 )then
            msg = 'solve3D manifest is incomplete'
        else if( .not. allocated(self%stages) )then
            msg = 'solve3D manifest has no ladder'
        else if( self%first_stage < 1 .or. self%last_stage < self%first_stage .or. self%last_stage > size(self%stages) )then
            msg = 'solve3D manifest stage range is inconsistent with its ladder'
        else
            status = 0
        endif
        if( status /= 0 )then
            call self%kill
        else
            self%fname = fname
        endif
    end subroutine read

    ! PROJECT REGISTRATION

    !> register the manifest in the project by name, with its run identifier
    subroutine register( self, spproj, name )
        class(solve3D_manifest), intent(in)       :: self
        class(sp_project),          intent(inout) :: spproj
        character(len=*),           intent(in)    :: name
        if( spproj%projinfo%get_noris() /= 1 ) call spproj%projinfo%new(1, is_ptcl=.false.)
        call spproj%projinfo%set(1, MANIFEST_PROJINFO_KEY, trim(name))
        call spproj%projinfo%set(1, RUN_ID_PROJINFO_KEY,   trim(self%run_id))
    end subroutine register

    !> Read the manifest a project registers. A bare name resolves against the
    !! directory of the project file itself (projinfo cwd is reset to the
    !! process directory by every builder read, so it is not a reliable
    !! project directory), and the manifest must carry the run identifier
    !! registered beside it.
    subroutine read_registered( self, spproj, projfile, status, msg )
        class(solve3D_manifest), intent(inout) :: self
        class(sp_project),          intent(in)    :: spproj
        class(string),              intent(in)    :: projfile
        integer,                    intent(out)   :: status
        character(len=*),           intent(out)   :: msg
        type(string) :: registered, projdir, path
        character(len=:), allocatable :: raw
        character(len=64) :: run_id_reg
        logical :: l_registered
        call self%kill
        status = 1
        msg    = ''
        l_registered = spproj%projinfo%get_noris() == 1
        if( l_registered ) l_registered = spproj%projinfo%isthere(1, MANIFEST_PROJINFO_KEY) .and. &
            &spproj%projinfo%isthere(1, RUN_ID_PROJINFO_KEY)
        if( .not. l_registered )then
            msg = 'the project registers no solve3D run manifest (only a completed solve3D run has one)'
            return
        endif
        call spproj%projinfo%getter(1, MANIFEST_PROJINFO_KEY, registered)
        raw = trim(registered%to_char())
        if( scan(raw, '/') > 0 )then
            path = raw
        else
            projdir = stemname(simple_abspath(projfile, check_exists=.false.))
            path    = filepath(projdir, registered)
        endif
        call spproj%projinfo%getter(1, RUN_ID_PROJINFO_KEY, registered)
        run_id_reg = trim(registered%to_char())
        call self%read(path, status, msg)
        if( status == 0 .and. trim(self%run_id) /= trim(run_id_reg) )then
            status = 1
            msg    = 'the manifest run identifier differs from the one registered in the project'
            call self%kill
        endif
        call registered%kill
        call projdir%kill
        call path%kill
    end subroutine read_registered

    ! USE AS A FROZEN INPUT

    !> A manifest may serve as a frozen input only for the project that
    !! registered it, unchanged: an eligible, completed solve3D or
    !! solve3D_addon run (an add-on output carries the union's sigma2
    !! state, so add-ons chain) with the same particle layout, stack table and
    !! optics/CTF identity, and every registered state map with the recorded
    !! digest. The registered run identifier is checked by read_registered.
    subroutine validate_frozen( self, spproj, status, msg )
        class(solve3D_manifest), intent(in)       :: self
        class(sp_project),          intent(inout) :: spproj
        integer,                    intent(out)   :: status
        character(len=*),           intent(out)   :: msg
        type(string) :: vol
        integer :: state, box
        real    :: smpd
        status = 1
        msg    = ''
        if( trim(self%program) /= 'solve3D' .and. trim(self%program) /= 'solve3D_addon' )then
            msg = 'the frozen project was not produced by solve3D or solve3D_addon'
            return
        endif
        if( .not. self%eligible )then
            msg = 'the frozen project is not eligible as a frozen input'
            return
        endif
        if( self%nrows /= spproj%os_ptcl3D%get_noris() )then
            msg = 'the frozen project particle count differs from its manifest'
            return
        endif
        if( self%layout_digest /= sigma2_state_project_layout_digest(spproj, spproj%os_ptcl3D) )then
            msg = 'the frozen project particle layout differs from its manifest'
            return
        endif
        if( self%stack_digest /= stack_table_digest(spproj) )then
            msg = 'the frozen project stack table differs from its manifest'
            return
        endif
        if( self%optics_digest /= optics_ctf_digest(spproj) )then
            msg = 'the frozen project optics/CTF parameters differ from its manifest'
            return
        endif
        do state = 1, self%nstates
            if( self%find_artifact('vol', state) == 0 )then
                msg = 'the manifest records no final map for state '//int2str(state)
                return
            endif
            if( .not. spproj%isthere_in_osout('vol', state) )then
                msg = 'the frozen project registers no map for state '//int2str(state)
                return
            endif
            call spproj%get_vol('vol', state, vol, smpd, box)
            if( .not. self%matches_artifact('vol', state, vol) )then
                msg = 'the registered map of state '//int2str(state)//' is not the map the base run produced'
                return
            endif
        enddo
        status = 0
        call vol%kill
    end subroutine validate_frozen

    !> The base run's settings onto a command line: the replayed inputs as
    !! given, the stage command line's shape at
    !! planning time (not the input text), and the completed solution (state
    !! layout, point group, mask, ladder end and effective nsample)
    subroutine replay( self, cline )
        class(solve3D_manifest), intent(in)       :: self
        class(cmdline),             intent(inout) :: cline
        type(string) :: val
        logical      :: found
        integer      :: i
        do i = 1, size(MANIFEST_REPLAY_KEYS)
            call self%get_input(trim(MANIFEST_REPLAY_KEYS(i)), val, found)
            if( found ) call cline%set_from_text(trim(MANIFEST_REPLAY_KEYS(i)), val%to_char())
        enddo
        call val%kill
        call cline%delete('lp')
        call cline%delete('lpstop')
        if( self%lp_on_line )     call cline%set('lp',     self%lp_line)
        if( self%lpstop_on_line ) call cline%set('lpstop', self%lpstop_line)
        call cline%set('pgrp',          trim(self%pgrp))
        call cline%set('pgrp_start',    trim(self%pgrp))
        call cline%set('mskdiam',       self%mskdiam)
        call cline%set('nstates',       self%nstates)
        call cline%set('nstages',       self%last_stage)
        call cline%set('nsample',       self%nsample)
    end subroutine replay

    !> the recorded path of an artifact; found=.false. when none is recorded
    subroutine get_artifact( self, kind, state, fname, found )
        class(solve3D_manifest), intent(in)       :: self
        character(len=*),           intent(in)    :: kind
        integer,                    intent(in)    :: state
        type(string),               intent(inout) :: fname
        logical,                    intent(out)   :: found
        integer :: ind
        call fname%kill
        ind   = self%find_artifact(kind, state)
        found = ind > 0
        if( found ) fname = self%artifact_path(ind)
    end subroutine get_artifact

    !> the file carries the digest recorded for the artifact (kind, state)
    logical function matches_artifact( self, kind, state, fname ) result( l_match )
        class(solve3D_manifest), intent(in) :: self
        character(len=*),           intent(in) :: kind
        integer,                    intent(in) :: state
        class(string),              intent(in) :: fname
        integer(int64) :: digest
        integer :: ind
        l_match = .false.
        ind = self%find_artifact(kind, state)
        if( ind == 0 ) return
        digest  = sigma2_state_digest_file(fname)
        l_match = digest /= 0_int64 .and. digest == self%artifact_digest(ind)
    end function matches_artifact

    !> the recorded value of an input key; found=.false. when the base run did not give it
    subroutine get_input( self, key, val, found )
        class(solve3D_manifest), intent(in)       :: self
        character(len=*),           intent(in)    :: key
        type(string),               intent(inout) :: val
        logical,                    intent(out)   :: found
        integer :: i
        call val%kill
        found = .false.
        if( .not. allocated(self%input_keys) ) return
        do i = 1, size(self%input_keys)
            if( trim(self%input_keys(i)) == trim(key) )then
                val   = self%input_vals(i)
                found = .true.
                return
            endif
        enddo
    end subroutine get_input

    integer function find_artifact( self, kind, state ) result( ind )
        class(solve3D_manifest), intent(in) :: self
        character(len=*),           intent(in) :: kind
        integer,                    intent(in) :: state
        integer :: i
        ind = 0
        if( .not. allocated(self%artifact_kind) ) return
        do i = 1, size(self%artifact_kind)
            if( trim(self%artifact_kind(i)) == trim(kind) .and. self%artifact_state(i) == state )then
                ind = i
                return
            endif
        enddo
    end function find_artifact

    ! GETTERS

    !> the file the manifest was read from
    function get_fname( self ) result( fname )
        class(solve3D_manifest), intent(in) :: self
        type(string) :: fname
        if( self%fname%is_allocated() ) fname = self%fname
    end function get_fname

    function get_run_id( self ) result( run_id )
        class(solve3D_manifest), intent(in) :: self
        character(len=64) :: run_id
        run_id = self%run_id
    end function get_run_id

    integer function get_nstates( self ) result( nstates )
        class(solve3D_manifest), intent(in) :: self
        nstates = self%nstates
    end function get_nstates

    !> native box of the solution
    integer function get_box( self ) result( box )
        class(solve3D_manifest), intent(in) :: self
        box = self%box
    end function get_box

    !> native sampling of the solution
    real function get_smpd( self ) result( smpd )
        class(solve3D_manifest), intent(in) :: self
        smpd = self%smpd
    end function get_smpd

    !> the last stage the run executed
    integer function get_last_stage( self ) result( last_stage )
        class(solve3D_manifest), intent(in) :: self
        last_stage = self%last_stage
    end function get_last_stage

    !> length of the ladder
    integer function get_nstages( self ) result( nstages )
        class(solve3D_manifest), intent(in) :: self
        nstages = 0
        if( allocated(self%stages) ) nstages = size(self%stages)
    end function get_nstages

    function get_stage( self, istage ) result( stage )
        class(solve3D_manifest), intent(in) :: self
        integer,                    intent(in) :: istage
        type(solve3D_stage_record) :: stage
        if( istage < 1 .or. istage > self%get_nstages() ) THROW_HARD('stage is outside the solve3D manifest ladder')
        stage = self%stages(istage)
    end function get_stage

    subroutine kill( self )
        class(solve3D_manifest), intent(inout) :: self
        call self%fname%kill
        self%run_id         = ''
        self%program        = ''
        self%complete       = .false.
        self%eligible       = .false.
        self%nrows          = 0
        self%layout_digest  = 0_int64
        self%stack_digest   = 0_int64
        self%optics_digest  = 0_int64
        self%nstates        = 0
        self%box            = 0
        self%smpd           = 0.
        self%mskdiam        = 0.
        self%pgrp           = ''
        self%nsample        = 0
        self%nptcls_eff     = 0
        self%update_frac    = 1.
        self%full_sampling  = .false.
        self%lp_on_line     = .false.
        self%lpstop_on_line = .false.
        self%lp_line        = 0.
        self%lpstop_line    = 0.
        self%first_stage    = 0
        self%last_stage     = 0
        if( allocated(self%stages)          ) deallocate(self%stages)
        if( allocated(self%input_keys)      ) deallocate(self%input_keys)
        if( allocated(self%input_vals)      ) deallocate(self%input_vals)
        if( allocated(self%artifact_kind)   ) deallocate(self%artifact_kind)
        if( allocated(self%artifact_state)  ) deallocate(self%artifact_state)
        if( allocated(self%artifact_digest) ) deallocate(self%artifact_digest)
        if( allocated(self%artifact_path)   ) deallocate(self%artifact_path)
        self%defect = ''
    end subroutine kill

    ! PRIVATE HELPERS

    !> identity of the stack table: every stack's name, particle range and physical size
    integer(int64) function stack_table_digest( spproj ) result( digest )
        class(sp_project), intent(in) :: spproj
        type(string) :: stk
        integer :: i
        digest = sigma2_state_digest_begin()
        do i = 1, spproj%os_stk%get_noris()
            stk = spproj%os_stk%get_str(i, 'stk')
            call sigma2_state_digest_text(digest, trim(stk%to_char()))
            call sigma2_state_digest_integer(digest, spproj%os_stk%get_int(i, 'fromp'))
            call sigma2_state_digest_integer(digest, spproj%os_stk%get_int(i, 'top'))
            call sigma2_state_digest_integer(digest, spproj%os_stk%get_int(i, 'box'))
            call sigma2_state_digest_text(digest, trim(real2str(spproj%os_stk%get(i, 'smpd'))))
        enddo
        call stk%kill
    end function stack_table_digest

    !> identity of the optics and CTF parameters of every particle row
    integer(int64) function optics_ctf_digest( spproj ) result( digest )
        class(sp_project), intent(inout) :: spproj
        type(ctfparams) :: ctfvars
        integer :: i
        digest = sigma2_state_digest_begin()
        do i = 1, spproj%os_ptcl3D%get_noris()
            ctfvars = spproj%get_ctfparams('ptcl3D', i)
            call sigma2_state_digest_integer(digest, int(ctfvars%ctfflag))
            call sigma2_state_digest_text(digest, trim(real2str(ctfvars%smpd))//' '//trim(real2str(ctfvars%kv))// &
                &' '//trim(real2str(ctfvars%cs))//' '//trim(real2str(ctfvars%fraca))//' '// &
                &trim(real2str(ctfvars%dfx))//' '//trim(real2str(ctfvars%dfy))//' '// &
                &trim(real2str(ctfvars%angast))//' '//trim(real2str(ctfvars%phshift)))
            if( spproj%os_ptcl3D%isthere(i, 'ogid') ) &
                &call sigma2_state_digest_integer(digest, spproj%os_ptcl3D%get_int(i, 'ogid'))
        enddo
    end function optics_ctf_digest

    function int64_str( i ) result( str )
        integer(int64), intent(in) :: i
        character(len=24) :: str
        write(str,'(I0)') i
    end function int64_str

    function logical_str( l ) result( str )
        logical, intent(in) :: l
        character(len=3) :: str
        str = merge('yes', 'no ', l)
    end function logical_str

    !> reals with full single-precision round trip
    function fmt_real( r ) result( str )
        real, intent(in) :: r
        character(len=24) :: str
        write(str,'(ES16.9)') r
        str = adjustl(str)
    end function fmt_real

end module simple_solve3D_manifest

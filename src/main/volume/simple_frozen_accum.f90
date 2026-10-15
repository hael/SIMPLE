!@descr: frozen accumulator store of one solve3D_addon run: run context, per-box sets, writers and validated adds into a reduction
! Raw frozen-particle accumulators (gridding even/odd S/rho), one set per state and box, added with
! coefficient one before any restoration or prior; nothing here scales or restores. Add-on chains
! reconstruct by gridding only, so the context and every set record the gridding backend.
! Sets are bound to the run context (run id, states, row counts, weighting, and the state
! weight set identity of the loading reconstruction, zero under hard labels) and
! validated against it and the consumer grid before any payload is read; defects are fatal.
! A carried store holds the union's raw sums an add-on publishes for the next add-on of its chain
! (carried_* names), with the rows they contain and a bundle manifest (carried_bundle), validated by
! headers, sizes and lineage, never by reading the payload.
! Contract: doc/policies/3D/solve3D_addon_policy.md sec. 5, 7 and 8.
module simple_frozen_accum
use simple_core_module_api
use simple_reconstructor,     only: reconstructor
use simple_refine3D_fnames,   only: refine3D_frozen_rec_fbody, refine3D_carried_rec_fbody, &
    &refine3D_carried_membership_fname
implicit none

public :: frozen_accum, carried_bundle, write_carried_membership, read_carried_membership
private
#include "simple_local_flags.inc"

character(len=*), parameter :: FROZEN_CONTEXT_SCHEMA = 'solve3D_addon_frozen_context'
character(len=*), parameter :: FROZEN_SET_SCHEMA     = 'solve3D_addon_frozen_set'
integer,          parameter :: FROZEN_SCHEMA_VERSION = 3   ! 3: sets record the state weight set identity
character(len=*), parameter :: FROZEN_BACKEND        = 'gridding'
real,             parameter :: SMPD_RELTOL           = 1.e-5
character(len=*), parameter :: BUNDLE_SCHEMA         = 'solve3D_addon_carried_bundle'
integer,          parameter :: BUNDLE_SCHEMA_VERSION = 1

!> one add-on run's frozen contribution: its identity (the run context) and
!! the per-box accumulator sets bound to it
type :: frozen_accum
    private
    character(len=64) :: run_id  = ''
    integer :: nstates      = 0 !< inherited state count
    integer :: nrows        = 0 !< rows of the add-on's working project (the consumers' index space)
    integer :: nrows_frozen = 0 !< rows of the frozen project (the producers' index space)
    integer :: nfrozen      = 0 !< frozen particles (state > 0 and updatecnt > 0 in the frozen project)
    integer, allocatable :: nfrozen_state(:)
    character(len=8)  :: weighting = '' !< the loading run's reconstruction weighting (euclid|cc)
    integer(kind=8)   :: wset_id(2) = 0_8 !< the loading run's state weight set identity (0 under hard labels)
    logical           :: l_carried = .false. !< a carried store: carried_* sets
    type(string)      :: dir                 !< directory of the sets, the working directory when empty
  contains
    ! the run context
    procedure          :: new
    procedure          :: write
    procedure          :: read
    procedure          :: validate
    procedure          :: load
    procedure          :: set_weighting
    procedure          :: set_dir
    procedure          :: get_nfrozen_state
    ! gridding sets
    procedure, private :: component
    procedure, private :: set_manifest
    procedure, private :: write_set_manifest
    procedure          :: write_gridding_set
    procedure          :: gridding_set_status
    procedure          :: add_gridding_set
    procedure          :: import_gridding_set
    procedure          :: kill
end type frozen_accum

!> The carried bundle of one add-on run: what its carried sets are bound to (run, rows, per-state
!! counts, weighting), where they are (states x boxes at one sampling), the run it continued (parent)
!! and the sigma2 generation it published with them; the membership file names the rows inside the sums
type :: carried_bundle
    character(len=64) :: run_id     = ''
    character(len=64) :: parent_run = ''
    character(len=8)  :: route      = ''   !< first: the frozen particles were read once; carried: none was read
    integer :: nstates = 0
    integer :: nrows   = 0
    integer, allocatable :: boxes(:)          !< the stage box and the native box (one when equal)
    real,    allocatable :: smpds(:)          !< the sampling at each box
    integer, allocatable :: ncarried_state(:)
    character(len=8) :: weighting = ''
    integer(kind=8)  :: wset_id(2) = 0_8
    integer(kind=8)  :: sigma2_generation = 0_8
  contains
    procedure :: write  => write_bundle
    procedure :: read   => read_bundle
    procedure :: store  => bundle_store
    procedure :: validate => validate_bundle
    procedure :: kill   => kill_bundle
end type carried_bundle

contains

    ! RUN CONTEXT

    !> the store of one add-on run: nrows and nrows_frozen are the row counts
    !! of the working and the frozen project, nstates is the length of the
    !! per-state frozen counts and nfrozen their sum
    subroutine new( self, run_id, nrows, nrows_frozen, nfrozen_state, carried )
        class(frozen_accum), intent(inout) :: self
        character(len=*),    intent(in)    :: run_id
        integer,             intent(in)    :: nrows, nrows_frozen
        integer,             intent(in)    :: nfrozen_state(:)
        logical, optional,   intent(in)    :: carried
        call self%kill
        if( present(carried) ) self%l_carried = carried
        if( len_trim(run_id) == 0 .or. index(trim(run_id), ' ') > 0 ) &
            &THROW_HARD('frozen context requires a blank-free run identifier')
        if( size(nfrozen_state) < 1 ) THROW_HARD('frozen context has no state layout')
        if( nrows < 1 .or. any(nfrozen_state < 0) ) THROW_HARD('frozen context particle counts are invalid')
        if( nrows_frozen < 1 ) THROW_HARD('frozen context row counts are invalid')
        self%run_id        = run_id
        self%nstates       = size(nfrozen_state)
        self%nrows         = nrows
        self%nrows_frozen  = nrows_frozen
        self%nfrozen_state = nfrozen_state
        self%nfrozen       = sum(nfrozen_state)
    end subroutine new

    !> atomically publish the run context (written to .tmp, then renamed)
    subroutine write( self, fname )
        class(frozen_accum), intent(in) :: self
        class(string),       intent(in) :: fname
        type(string) :: tmpname
        integer :: funit, io_stat, s
        if( self%nstates < 1 .or. .not. allocated(self%nfrozen_state) ) THROW_HARD('frozen context has no state layout')
        tmpname = fname//'.tmp'
        call del_file(tmpname)
        call fopen(funit, file=tmpname, status='REPLACE', action='WRITE', iostat=io_stat)
        call fileiochk('frozen_accum%write opening '//tmpname%to_char(), io_stat)
        write(funit,'(A,1X,I0)') FROZEN_CONTEXT_SCHEMA, FROZEN_SCHEMA_VERSION
        write(funit,'(A,1X,A)')  'run_id',  trim(self%run_id)
        write(funit,'(A,1X,A)')  'backend', FROZEN_BACKEND
        write(funit,'(A,1X,I0)') 'nstates', self%nstates
        write(funit,'(A,1X,I0)') 'nrows',   self%nrows
        write(funit,'(A,1X,I0)') 'nrows_frozen', self%nrows_frozen
        write(funit,'(A,1X,I0)') 'nfrozen', self%nfrozen
        do s = 1, self%nstates
            write(funit,'(A,1X,I0,1X,I0)') 'nfrozen_state', s, self%nfrozen_state(s)
        enddo
        write(funit,'(A)') 'end'
        call fclose(funit)
        call simple_rename(tmpname, fname, overwrite=.true.)
        call tmpname%kill
    end subroutine write

    !> parse a run context; status /= 0 with a message on any defect
    subroutine read( self, fname, status, msg )
        class(frozen_accum), intent(inout) :: self
        class(string),       intent(in)    :: fname
        integer,             intent(out)   :: status
        character(len=*),    intent(out)   :: msg
        character(len=XLONGSTRLEN) :: line
        character(len=64) :: key, schema, backend
        integer :: funit, io_stat, version, s, cnt, nseen
        logical :: l_end, l_seen(6)
        call self%kill
        status = 1
        msg    = ''
        l_end  = .false.
        l_seen = .false.
        nseen  = 0
        backend = ''
        if( .not. file_exists(fname) )then
            msg = 'frozen context file is missing: '//fname%to_char()
            return
        endif
        call fopen(funit, file=fname, status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 )then
            msg = 'frozen context file is unreadable'
            return
        endif
        read(funit,'(A)',iostat=io_stat) line
        if( io_stat == 0 ) read(line,*,iostat=io_stat) schema, version
        if( io_stat /= 0 .or. trim(schema) /= FROZEN_CONTEXT_SCHEMA )then
            msg = 'not a frozen context file'
            call fclose(funit)
            return
        endif
        if( version /= FROZEN_SCHEMA_VERSION )then
            msg = 'unsupported frozen context schema version'
            call fclose(funit)
            return
        endif
        do
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            read(line,*,iostat=io_stat) key
            if( io_stat /= 0 ) exit
            select case(trim(key))
                case('run_id')
                    read(line,*,iostat=io_stat) key, self%run_id
                    l_seen(1) = .true.
                case('backend')
                    read(line,*,iostat=io_stat) key, backend
                    l_seen(2) = .true.
                case('nstates')
                    read(line,*,iostat=io_stat) key, self%nstates
                    if( io_stat == 0 .and. self%nstates >= 1 .and. self%nstates <= MAXS )then
                        allocate(self%nfrozen_state(self%nstates), source=-1)
                        l_seen(3) = .true.
                    else
                        io_stat = 1
                    endif
                case('nrows')
                    read(line,*,iostat=io_stat) key, self%nrows
                    l_seen(4) = .true.
                case('nrows_frozen')
                    read(line,*,iostat=io_stat) key, self%nrows_frozen
                    l_seen(6) = .true.
                case('nfrozen')
                    read(line,*,iostat=io_stat) key, self%nfrozen
                    l_seen(5) = .true.
                case('nfrozen_state')
                    read(line,*,iostat=io_stat) key, s, cnt
                    if( io_stat == 0 )then
                        if( .not. allocated(self%nfrozen_state) )then
                            io_stat = 1
                        else if( s < 1 .or. s > self%nstates )then
                            io_stat = 1
                        else
                            self%nfrozen_state(s) = cnt
                            nseen = nseen + 1
                        endif
                    endif
                case('end')
                    l_end = .true.
                    exit
                case DEFAULT
                    msg = 'unknown frozen context field: '//trim(key)
                    call fclose(funit)
                    call self%kill
                    return
            end select
            if( io_stat /= 0 ) exit
        enddo
        call fclose(funit)
        if( io_stat /= 0 )then
            msg = 'corrupt frozen context record'
        else if( .not. l_end )then
            msg = 'truncated frozen context file'
        else if( .not. all(l_seen) )then
            msg = 'incomplete frozen context file'
        else if( trim(backend) /= FROZEN_BACKEND )then
            msg = 'frozen context backend '//trim(backend)//' is not gridding: add-on chains reconstruct by gridding only'
        else if( nseen /= self%nstates .or. any(self%nfrozen_state < 0) )then
            msg = 'frozen context state counts are incomplete'
        else if( self%nrows < 1 .or. self%nfrozen < 1 .or. sum(self%nfrozen_state) /= self%nfrozen )then
            msg = 'frozen context particle counts are inconsistent'
        else if( self%nrows_frozen < 1 )then
            msg = 'frozen context row counts are inconsistent'
        else
            status = 0
        endif
        if( status /= 0 ) call self%kill
    end subroutine read

    !> the context must describe the loading run: state layout and particle
    !! index space, which is the frozen project's for a producer (a
    !! reconstruct3D on the frozen project) and the working project's for a
    !! consumer (an add-on reconstruction)
    subroutine validate( self, nstates, nrows, producer, status, msg )
        class(frozen_accum), intent(in)  :: self
        integer,             intent(in)  :: nstates, nrows
        logical,             intent(in)  :: producer
        integer,             intent(out) :: status
        character(len=*),    intent(out) :: msg
        status = 1
        msg    = ''
        if( self%nstates /= nstates )then
            msg = 'frozen context state layout differs from the reconstruction'
        else if( producer .and. self%nrows_frozen /= nrows )then
            msg = 'frozen context particle index space differs from the frozen project'
        else if( .not. producer .and. self%nrows /= nrows )then
            msg = 'frozen context particle index space differs from the reconstruction project'
        else
            status = 0
        endif
    end subroutine validate

    !> Read and validate the context named by a handshake; any defect is fatal:
    !! a consumer never falls back to a reconstruction without its frozen term.
    !! cc_objfun is the loading run's reconstruction weighting; producer is
    !! true for the frozen_seed handshake, false for frozen_rec.
    subroutine load( self, fname, nstates, nrows, cc_objfun, producer, wset_id, carried )
        class(frozen_accum), intent(inout) :: self
        class(string),       intent(in)    :: fname
        integer,             intent(in)    :: nstates, nrows, cc_objfun
        logical,             intent(in)    :: producer
        integer(kind=8),     intent(in)    :: wset_id(2)
        logical, optional,   intent(in)    :: carried
        character(len=STDLEN) :: msg
        integer :: status
        call self%read(fname, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call self%validate(nstates, nrows, producer, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call self%set_weighting(cc_objfun, wset_id)
        if( present(carried) ) self%l_carried = carried
    end subroutine load

    !> the reconstruction weighting and state weight set identity the sets are bound to
    subroutine set_weighting( self, cc_objfun, wset_id )
        class(frozen_accum), intent(inout) :: self
        integer,             intent(in)    :: cc_objfun
        integer(kind=8),     intent(in)    :: wset_id(2)
        self%weighting = trim(merge('euclid', 'cc    ', cc_objfun == OBJFUN_EUCLID))
        self%wset_id   = wset_id
    end subroutine set_weighting

    !> the directory holding the sets (a carried bundle of an earlier run)
    subroutine set_dir( self, dir )
        class(frozen_accum), intent(inout) :: self
        class(string),       intent(in)    :: dir
        self%dir = dir
    end subroutine set_dir

    !> frozen particles of one inherited state
    integer function get_nfrozen_state( self, state ) result( n )
        class(frozen_accum), intent(in) :: self
        integer,             intent(in) :: state
        if( state < 1 .or. state > self%nstates ) THROW_HARD('state is outside the frozen context layout')
        n = self%nfrozen_state(state)
    end function get_nfrozen_state

    subroutine kill( self )
        class(frozen_accum), intent(inout) :: self
        self%run_id  = ''
        self%nstates = 0
        self%nrows   = 0
        self%nrows_frozen = 0
        self%nfrozen = 0
        if( allocated(self%nfrozen_state) ) deallocate(self%nfrozen_state)
        self%weighting = ''
        self%wset_id   = 0_8
        self%l_carried = .false.
        call self%dir%kill
    end subroutine kill

    ! GRIDDING SETS

    !> one of the four components of the set of state at box, in the store's directory
    function component( self, state, box, ifile ) result( fname )
        class(frozen_accum), intent(in) :: self
        integer,             intent(in) :: state, box, ifile
        type(string) :: fname, fbody
        if( self%l_carried )then
            fbody = refine3D_carried_rec_fbody(state, box)
        else
            fbody = refine3D_frozen_rec_fbody(state, box)
        endif
        select case(ifile)
            case(1)
                fname = fbody//'_even'//MRC_EXT
            case(2)
                fname = string('rho_')//fbody//'_even'//MRC_EXT
            case(3)
                fname = fbody//'_odd'//MRC_EXT
            case(4)
                fname = string('rho_')//fbody//'_odd'//MRC_EXT
            case DEFAULT
                THROW_HARD('invalid frozen gridding component index')
        end select
        if( self%dir%is_allocated() ) fname = self%dir//'/'//fname
        call fbody%kill
    end function component

    !> the manifest of the set of state at box, written last
    function set_manifest( self, state, box ) result( fname )
        class(frozen_accum), intent(in) :: self
        integer,             intent(in) :: state, box
        type(string) :: fname
        if( self%l_carried )then
            fname = refine3D_carried_rec_fbody(state, box)//TXT_EXT
        else
            fname = refine3D_frozen_rec_fbody(state, box)//TXT_EXT
        endif
        if( self%dir%is_allocated() ) fname = self%dir//'/'//fname
    end function set_manifest

    integer(kind=8) function file_size( fname ) result( sz )
        class(string), intent(in) :: fname
        sz = -1
        if( file_exists(fname) ) inquire(file=fname%to_char(), size=sz)
    end function file_size

    !> Publish the frozen even/odd accumulators of one state at the accumulators'
    !! own box as one artifact set: the manifest is removed first and written
    !! last with every component's byte size, so an interrupted write leaves no
    !! valid set.
    subroutine write_gridding_set( self, state, even_rec, odd_rec )
        class(frozen_accum),  intent(in)    :: self
        integer,              intent(in)    :: state
        class(reconstructor), intent(inout) :: even_rec, odd_rec
        type(string)    :: manifest
        integer(kind=8) :: sizes(4)
        integer :: ldim(3), ifile
        real    :: smpd
        if( state < 1 .or. state > self%nstates ) THROW_HARD('frozen gridding set state is outside the context layout')
        ldim = even_rec%get_ldim()
        smpd = even_rec%get_smpd()
        if( any(odd_rec%get_ldim() /= ldim) ) THROW_HARD('frozen gridding pair dimensions do not match')
        manifest = self%set_manifest(state, ldim(1))
        call del_file(manifest)
        call even_rec%write_raw_accum(self%component(state, ldim(1), 1), self%component(state, ldim(1), 2))
        call odd_rec%write_raw_accum( self%component(state, ldim(1), 3), self%component(state, ldim(1), 4))
        do ifile = 1, 4
            sizes(ifile) = file_size(self%component(state, ldim(1), ifile))
        enddo
        call self%write_set_manifest(state, ldim(1), smpd, sizes)
        if( self%l_carried )then
            write(logfhandle,'(A,I0,A,I0)') '>>> CARRIED GRIDDING SET WRITTEN, STATE ', state, ', BOX ', ldim(1)
        else
            write(logfhandle,'(A,I0,A,I0)') '>>> FROZEN GRIDDING SET WRITTEN, STATE ', state, ', BOX ', ldim(1)
        endif
        call manifest%kill
    end subroutine write_gridding_set

    !> the manifest of the set of state at box, bound to this store, through a temporary name
    subroutine write_set_manifest( self, state, box, smpd, sizes )
        class(frozen_accum), intent(in) :: self
        integer,             intent(in) :: state, box
        real,                intent(in) :: smpd
        integer(kind=8),     intent(in) :: sizes(4)
        type(string) :: manifest, tmpname
        integer :: funit, io_stat
        manifest = self%set_manifest(state, box)
        tmpname  = manifest//'.tmp'
        call del_file(tmpname)
        call fopen(funit, file=tmpname, status='REPLACE', action='WRITE', iostat=io_stat)
        call fileiochk('frozen_accum%write_set_manifest opening '//tmpname%to_char(), io_stat)
        write(funit,'(A,1X,I0)')     FROZEN_SET_SCHEMA, FROZEN_SCHEMA_VERSION
        write(funit,'(A,1X,A)')      'run_id',  trim(self%run_id)
        write(funit,'(A,1X,A)')      'backend', FROZEN_BACKEND
        write(funit,'(A,1X,I0)')     'state',   state
        write(funit,'(A,1X,I0)')     'nstates', self%nstates
        write(funit,'(A,1X,I0)')     'box',     box
        write(funit,'(A,1X,ES16.9)') 'smpd',    smpd
        write(funit,'(A,1X,I0)')     'nrows',   self%nrows
        write(funit,'(A,1X,I0)')     'nfrozen_state', self%nfrozen_state(state)
        write(funit,'(A,1X,A)')      'weighting', trim(merge(self%weighting, 'none    ', len_trim(self%weighting) > 0))
        write(funit,'(A,2(1X,I0))')  'wset',    self%wset_id
        write(funit,'(A,4(1X,I0))')  'sizes',   sizes
        write(funit,'(A)')           'end'
        call fclose(funit)
        call simple_rename(tmpname, manifest, overwrite=.true.)
        call manifest%kill
        call tmpname%kill
    end subroutine write_set_manifest

    !> The set of state at (box, smpd) of another store (a carried bundle), validated against that
    !! store, becomes this store's set: its components are copied under this store's names and its
    !! manifest is bound to this store (run, rows, counts, weighting). No component is decoded.
    subroutine import_gridding_set( self, src, state, box, smpd )
        class(frozen_accum), intent(in) :: self
        class(frozen_accum), intent(in) :: src
        integer,             intent(in) :: state, box
        real,                intent(in) :: smpd
        character(len=STDLEN) :: msg
        integer(kind=8) :: sizes(4)
        integer :: status, ifile
        if( state < 1 .or. state > self%nstates ) THROW_HARD('imported gridding set state is outside the context layout')
        if( self%nfrozen_state(state) /= src%nfrozen_state(state) ) &
            &THROW_HARD('an imported gridding set must keep its particle count')
        if( trim(self%weighting) /= trim(src%weighting) .or. any(self%wset_id /= src%wset_id) ) &
            &THROW_HARD('an imported gridding set must keep its weighting and state weight set')
        call src%gridding_set_status(state, box, smpd, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call del_file(self%set_manifest(state, box))
        do ifile = 1, 4
            call simple_copy_file(src%component(state, box, ifile), self%component(state, box, ifile))
            sizes(ifile) = file_size(self%component(state, box, ifile))
        enddo
        call self%write_set_manifest(state, box, smpd, sizes)
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> CARRIED GRIDDING SET ADOPTED, STATE ', state, ', BOX ', box, ' (', &
            &self%nfrozen_state(state), ' PARTICLES, NO PARTICLE READ)'
    end subroutine import_gridding_set

    !> Validate the gridding set a consumer at (box, smpd) would read for state:
    !! the manifest must be complete and bound to the context, the grid must be
    !! exactly the consumer's and every component must have its recorded size.
    subroutine gridding_set_status( self, state, box, smpd, status, msg )
        class(frozen_accum), intent(in)  :: self
        integer,             intent(in)  :: state, box
        real,                intent(in)  :: smpd
        integer,             intent(out) :: status
        character(len=*),    intent(out) :: msg
        character(len=XLONGSTRLEN) :: line
        character(len=64)  :: key, schema, run_id, backend, weighting
        type(string)       :: manifest
        integer(kind=8)    :: sizes(4), wset_set(2)
        integer :: funit, io_stat, version, state_set, nstates_set, box_set, nrows_set, nfrozen_set, ifile
        real    :: smpd_set
        logical :: l_end
        status = 1
        msg    = ''
        l_end  = .false.
        run_id = ''; backend = ''; weighting = ''
        state_set = -1; nstates_set = -1; box_set = -1; nrows_set = -1; nfrozen_set = -1; smpd_set = -1.; sizes = -2
        wset_set = -1_8
        manifest = self%set_manifest(state, box)
        if( .not. file_exists(manifest) )then
            msg = 'frozen gridding set is missing: '//manifest%to_char()
            call manifest%kill
            return
        endif
        call fopen(funit, file=manifest, status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 )then
            msg = 'frozen gridding manifest is unreadable'
            call manifest%kill
            return
        endif
        read(funit,'(A)',iostat=io_stat) line
        if( io_stat == 0 ) read(line,*,iostat=io_stat) schema, version
        if( io_stat /= 0 .or. trim(schema) /= FROZEN_SET_SCHEMA .or. version /= FROZEN_SCHEMA_VERSION )then
            msg = 'not a frozen gridding set manifest of the supported schema'
            call fclose(funit)
            call manifest%kill
            return
        endif
        do
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            read(line,*,iostat=io_stat) key
            if( io_stat /= 0 ) exit
            select case(trim(key))
                case('run_id');        read(line,*,iostat=io_stat) key, run_id
                case('backend');       read(line,*,iostat=io_stat) key, backend
                case('state');         read(line,*,iostat=io_stat) key, state_set
                case('nstates');       read(line,*,iostat=io_stat) key, nstates_set
                case('box');           read(line,*,iostat=io_stat) key, box_set
                case('smpd');          read(line,*,iostat=io_stat) key, smpd_set
                case('nrows');         read(line,*,iostat=io_stat) key, nrows_set
                case('nfrozen_state'); read(line,*,iostat=io_stat) key, nfrozen_set
                case('weighting');     read(line,*,iostat=io_stat) key, weighting
                case('wset');          read(line,*,iostat=io_stat) key, wset_set
                case('sizes');         read(line,*,iostat=io_stat) key, sizes
                case('end')
                    l_end = .true.
                    exit
                case DEFAULT
                    io_stat = 1
            end select
            if( io_stat /= 0 ) exit
        enddo
        call fclose(funit)
        call manifest%kill
        if( io_stat /= 0 .or. .not. l_end )then
            msg = 'corrupt or truncated frozen gridding manifest'
        else if( trim(run_id) /= trim(self%run_id) )then
            msg = 'frozen gridding set belongs to another add-on run'
        else if( trim(backend) /= FROZEN_BACKEND )then
            msg = 'frozen set backend is not gridding'
        else if( state_set /= state .or. nstates_set /= self%nstates )then
            msg = 'frozen gridding set state layout mismatch'
        else if( box_set /= box )then
            msg = 'frozen gridding set box differs from the consuming grid'
        else if( abs(smpd_set - smpd) > SMPD_RELTOL * smpd )then
            msg = 'frozen gridding set sampling differs from the consuming grid'
        else if( nrows_set /= self%nrows )then
            msg = 'frozen gridding set particle index space mismatch'
        else if( nfrozen_set /= self%nfrozen_state(state) )then
            msg = 'frozen gridding set frozen count differs from the run context'
        else if( trim(weighting) /= trim(merge(self%weighting, 'none    ', len_trim(self%weighting) > 0)) )then
            msg = 'frozen gridding set was accumulated with another objective-function weighting'
        else if( any(wset_set /= self%wset_id) )then
            msg = 'frozen gridding set was accumulated under another state weight set'
        else
            status = 0
            do ifile = 1, 4
                if( file_size(self%component(state, box, ifile)) /= sizes(ifile) )then
                    msg    = 'frozen gridding set component size mismatch (incomplete or mixed set)'
                    status = 1
                    exit
                endif
            enddo
        endif
    end subroutine gridding_set_status

    !> Add the validated frozen set of state into the consumer's even/odd
    !! accumulators with coefficient one; read_even/read_odd are scratch
    !! accumulators on the same grid. A missing or invalid set is fatal.
    subroutine add_gridding_set( self, state, even_rec, odd_rec, read_even, read_odd )
        use simple_imgfile, only: imgfile
        class(frozen_accum),  intent(in)    :: self
        integer,              intent(in)    :: state
        class(reconstructor), intent(inout) :: even_rec, odd_rec, read_even, read_odd
        character(len=STDLEN) :: msg
        type(imgfile) :: ioimg
        type(string)  :: fname
        integer :: ldim(3), status, funit, ierr
        real    :: smpd
        ldim = even_rec%get_ldim()
        smpd = even_rec%get_smpd()
        if( any(odd_rec%get_ldim() /= ldim) .or. any(read_even%get_ldim() /= ldim) .or. &
            &any(read_odd%get_ldim() /= ldim) ) THROW_HARD('frozen add requires accumulators on one grid')
        call self%gridding_set_status(state, ldim(1), smpd, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call read_one(read_even, 1)
        call read_one(read_odd,  3)
        call even_rec%sum_reduce(read_even)
        call odd_rec%sum_reduce(read_odd)
        call fname%kill

    contains

        subroutine read_one( rec, ifile )
            class(reconstructor), intent(inout) :: rec
            integer,              intent(in)    :: ifile
            call rec%reset
            fname = self%component(state, ldim(1), ifile)
            call ioimg%open(fname, ldim, smpd, formatchar='M', readhead=.false., rwaction='read')
            call rec%read_raw_mrc(ioimg)
            call ioimg%close
            fname = self%component(state, ldim(1), ifile+1)
            call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
            call fileiochk('frozen_accum%add_gridding_set opening '//fname%to_char(), ierr)
            call rec%read_raw_rho(funit)
            call fclose(funit)
        end subroutine read_one

    end subroutine add_gridding_set

    ! CARRIED BUNDLE

    !> publish the bundle manifest atomically: the sets and the membership are complete before it
    subroutine write_bundle( self, fname )
        class(carried_bundle), intent(in) :: self
        class(string),         intent(in) :: fname
        type(string) :: tmpname
        integer :: funit, io_stat, s
        if( self%nstates < 1 .or. .not. allocated(self%ncarried_state) .or. .not. allocated(self%boxes) ) &
            &THROW_HARD('carried bundle has no layout')
        tmpname = fname//'.tmp'
        call del_file(tmpname)
        call fopen(funit, file=tmpname, status='REPLACE', action='WRITE', iostat=io_stat)
        call fileiochk('carried_bundle%write opening '//tmpname%to_char(), io_stat)
        write(funit,'(A,1X,I0)')     BUNDLE_SCHEMA, BUNDLE_SCHEMA_VERSION
        write(funit,'(A,1X,A)')      'run_id',     trim(self%run_id)
        write(funit,'(A,1X,A)')      'parent_run', trim(self%parent_run)
        write(funit,'(A,1X,A)')      'route',      trim(self%route)
        write(funit,'(A,1X,A)')      'backend',    FROZEN_BACKEND
        write(funit,'(A,1X,I0)')     'nstates',    self%nstates
        write(funit,'(A,1X,I0)')     'nrows',      self%nrows
        write(funit,'(A,1X,A)')      'weighting',  trim(self%weighting)
        write(funit,'(A,2(1X,I0))')  'wset',       self%wset_id
        write(funit,'(A,1X,I0)')     'sigma2_generation', self%sigma2_generation
        write(funit,'(A,*(1X,I0))')  'boxes',      size(self%boxes), self%boxes
        write(funit,'(A,*(1X,ES16.9))') 'smpds',   self%smpds
        do s = 1, self%nstates
            write(funit,'(A,1X,I0,1X,I0)') 'ncarried_state', s, self%ncarried_state(s)
        enddo
        write(funit,'(A)') 'end'
        call fclose(funit)
        call simple_rename(tmpname, fname, overwrite=.true.)
        call tmpname%kill
    end subroutine write_bundle

    !> parse a bundle manifest; status /= 0 with a message on any defect
    subroutine read_bundle( self, fname, status, msg )
        class(carried_bundle), intent(inout) :: self
        class(string),         intent(in)    :: fname
        integer,               intent(out)   :: status
        character(len=*),      intent(out)   :: msg
        character(len=XLONGSTRLEN) :: line
        character(len=64) :: key, schema, backend
        integer :: funit, io_stat, version, s, cnt, nboxes, nseen
        logical :: l_end, l_states, l_boxes, l_smpds
        call self%kill
        l_states = .false.
        l_boxes  = .false.
        l_smpds  = .false.
        status  = 1
        msg     = ''
        l_end   = .false.
        backend = ''
        nseen   = 0
        if( .not. file_exists(fname) )then
            msg = 'carried bundle manifest is missing: '//fname%to_char()
            return
        endif
        call fopen(funit, file=fname, status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 )then
            msg = 'carried bundle manifest is unreadable'
            return
        endif
        read(funit,'(A)',iostat=io_stat) line
        if( io_stat == 0 ) read(line,*,iostat=io_stat) schema, version
        if( io_stat /= 0 .or. trim(schema) /= BUNDLE_SCHEMA .or. version /= BUNDLE_SCHEMA_VERSION )then
            msg = 'not a carried bundle manifest of the supported schema'
            call fclose(funit)
            return
        endif
        do
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            read(line,*,iostat=io_stat) key
            if( io_stat /= 0 ) exit
            select case(trim(key))
                case('run_id');            read(line,*,iostat=io_stat) key, self%run_id
                case('parent_run');        read(line,*,iostat=io_stat) key, self%parent_run
                case('route');             read(line,*,iostat=io_stat) key, self%route
                case('backend');           read(line,*,iostat=io_stat) key, backend
                case('smpds')
                    if( l_boxes .and. .not. l_smpds )then
                        l_smpds = .true.
                        allocate(self%smpds(size(self%boxes)))
                        read(line,*,iostat=io_stat) key, self%smpds
                    else
                        io_stat = 1
                    endif
                case('weighting');         read(line,*,iostat=io_stat) key, self%weighting
                case('wset');              read(line,*,iostat=io_stat) key, self%wset_id
                case('sigma2_generation'); read(line,*,iostat=io_stat) key, self%sigma2_generation
                case('nrows');             read(line,*,iostat=io_stat) key, self%nrows
                case('nstates')
                    if( l_states ) io_stat = 1
                    if( io_stat == 0 ) read(line,*,iostat=io_stat) key, self%nstates
                    l_states = .true.
                    if( io_stat == 0 .and. self%nstates >= 1 .and. self%nstates <= MAXS )then
                        allocate(self%ncarried_state(self%nstates), source=-1)
                    else
                        io_stat = 1
                    endif
                case('boxes')
                    if( l_boxes ) io_stat = 1
                    if( io_stat == 0 ) read(line,*,iostat=io_stat) key, nboxes
                    l_boxes = .true.
                    if( io_stat == 0 .and. nboxes >= 1 .and. nboxes <= 2 )then
                        allocate(self%boxes(nboxes))
                        read(line,*,iostat=io_stat) key, nboxes, self%boxes
                    else
                        io_stat = 1
                    endif
                case('ncarried_state')
                    read(line,*,iostat=io_stat) key, s, cnt
                    if( io_stat == 0 )then
                        if( .not. allocated(self%ncarried_state) )then
                            io_stat = 1
                        else if( s < 1 .or. s > self%nstates )then
                            io_stat = 1
                        else
                            self%ncarried_state(s) = cnt
                            nseen = nseen + 1
                        endif
                    endif
                case('end')
                    l_end = .true.
                    exit
                case DEFAULT
                    io_stat = 1
            end select
            if( io_stat /= 0 ) exit
        enddo
        call fclose(funit)
        if( io_stat /= 0 .or. .not. l_end )then
            msg = 'corrupt or truncated carried bundle manifest'
        else if( trim(backend) /= FROZEN_BACKEND )then
            msg = 'carried bundle backend '//trim(backend)//' is not gridding'
        else if( len_trim(self%run_id) == 0 .or. self%nrows < 1 .or. .not. l_states .or. .not. l_boxes .or. &
            &.not. l_smpds .or. self%sigma2_generation < 1_8 )then
            msg = 'incomplete carried bundle manifest'
        else if( nseen /= self%nstates )then
            msg = 'carried bundle state counts are incomplete'
        else if( trim(self%route) /= 'first' .and. trim(self%route) /= 'carried' )then
            msg = 'carried bundle route is unknown'
        else if( any(self%smpds <= 0.) )then
            msg = 'carried bundle sampling is invalid'
        else if( any(self%ncarried_state < 0) .or. sum(self%ncarried_state) < 1 )then
            msg = 'carried bundle particle counts are invalid'
        else
            status = 0
        endif
        if( status /= 0 ) call self%kill
    end subroutine read_bundle

    !> the carried store the bundle describes, its sets in dir
    function bundle_store( self, dir ) result( store )
        class(carried_bundle), intent(in) :: self
        class(string),         intent(in) :: dir
        type(frozen_accum) :: store
        call store%new(self%run_id, self%nrows, self%nrows, self%ncarried_state, carried=.true.)
        store%weighting = self%weighting
        store%wset_id   = self%wset_id
        call store%set_dir(dir)
    end function bundle_store

    !> every set (state x box) of the bundle in dir is complete and bound to it, and the
    !! membership file has the bundle's row count: headers and sizes only
    subroutine validate_bundle( self, dir, status, msg )
        class(carried_bundle), intent(in)  :: self
        class(string),         intent(in)  :: dir
        integer,               intent(out) :: status
        character(len=*),      intent(out) :: msg
        type(frozen_accum) :: store
        integer, allocatable :: states(:)
        integer :: s, b
        store = self%store(dir)
        do s = 1, self%nstates
            do b = 1, size(self%boxes)
                call store%gridding_set_status(s, self%boxes(b), self%smpds(b), status, msg)
                if( status /= 0 )then
                    msg = 'carried bundle: '//trim(msg)
                    call store%kill
                    return
                endif
            enddo
        enddo
        call store%kill
        call read_carried_membership(dir//'/'//refine3D_carried_membership_fname(), self%nrows, states, status)
        if( status /= 0 )then
            msg = 'carried bundle membership is missing or does not have the bundle''s rows'
            return
        endif
        do s = 1, self%nstates
            if( count(states == s) /= self%ncarried_state(s) )then
                status = 1
                msg    = 'carried bundle membership disagrees with its particle counts'
                return
            endif
        enddo
        status = 0
        msg    = ''
    end subroutine validate_bundle

    subroutine kill_bundle( self )
        class(carried_bundle), intent(inout) :: self
        self%run_id     = ''
        self%parent_run = ''
        self%route      = ''
        self%nstates    = 0
        self%nrows      = 0
        if( allocated(self%boxes) )          deallocate(self%boxes)
        if( allocated(self%smpds) )          deallocate(self%smpds)
        if( allocated(self%ncarried_state) ) deallocate(self%ncarried_state)
        self%weighting = ''
        self%wset_id   = 0_8
        self%sigma2_generation = 0_8
    end subroutine kill_bundle

    ! MEMBERSHIP

    !> the rows inside a carried bundle's sums: one int32 state label per row (0: not inside),
    !! written through a temporary name
    subroutine write_carried_membership( fname, states )
        class(string), intent(in) :: fname
        integer,       intent(in) :: states(:)
        type(string) :: tmpname
        integer :: funit, io_stat
        tmpname = fname//'.tmp'
        call del_file(tmpname)
        call fopen(funit, file=tmpname, status='REPLACE', action='WRITE', access='STREAM', iostat=io_stat)
        call fileiochk('write_carried_membership opening '//tmpname%to_char(), io_stat)
        write(funit) int(states, kind=4)
        call fclose(funit)
        call simple_rename(tmpname, fname, overwrite=.true.)
        call tmpname%kill
    end subroutine write_carried_membership

    !> the membership of nrows rows; status /= 0 when the file is missing or of another size
    subroutine read_carried_membership( fname, nrows, states, status )
        class(string),        intent(in)  :: fname
        integer,              intent(in)  :: nrows
        integer, allocatable, intent(out) :: states(:)
        integer,              intent(out) :: status
        integer(kind=4), allocatable :: buf(:)
        integer :: funit, io_stat
        status = 1
        if( nrows < 1 ) return
        if( file_size(fname) /= 4_8 * int(nrows,8) ) return
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=io_stat)
        if( io_stat /= 0 ) return
        allocate(buf(nrows))
        read(funit, iostat=io_stat) buf
        call fclose(funit)
        if( io_stat /= 0 ) return
        allocate(states(nrows))
        states = int(buf)
        if( any(states < 0) ) return
        status = 0
    end subroutine read_carried_membership

end module simple_frozen_accum

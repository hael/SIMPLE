!@descr: frozen accumulator store of one solve3D_addon run: run context, per-box sets, writers and validated adds into a reduction
! Raw frozen-particle accumulators (gridding even/odd S/rho, PCG B/D per half), one set per state
! and box, added with coefficient one before any restoration or prior; nothing here scales or restores.
! Sets are bound to the run context (run id, backend, states, row counts, weighting) and
! validated against it and the consumer grid before any payload is read; defects are fatal.
! Contract: doc/policies/3D/solve3D_addon_policy.md sec. 5 and 7.
module simple_frozen_accum
use simple_core_module_api
use simple_reconstructor,     only: reconstructor
use simple_reconstructor_pcg, only: reconstructor_pcg, read_pcg_raw_accum_header
use simple_refine3D_fnames,   only: refine3D_frozen_rec_fbody, refine3D_frozen_manifest_fname, &
    &refine3D_frozen_pcg_fname
implicit none

public :: frozen_accum
private
#include "simple_local_flags.inc"

character(len=*), parameter :: FROZEN_CONTEXT_SCHEMA = 'solve3D_addon_frozen_context'
character(len=*), parameter :: FROZEN_SET_SCHEMA     = 'solve3D_addon_frozen_set'
integer,          parameter :: FROZEN_SCHEMA_VERSION = 2
real,             parameter :: SMPD_RELTOL           = 1.e-5

!> one add-on run's frozen contribution: its identity (the run context) and
!! the per-box accumulator sets bound to it
type :: frozen_accum
    private
    character(len=64) :: run_id  = ''
    character(len=16) :: backend = ''
    integer :: nstates      = 0 !< inherited state count
    integer :: nrows        = 0 !< rows of the add-on's working project (the consumers' index space)
    integer :: nrows_frozen = 0 !< rows of the frozen project (the producers' index space)
    integer :: nfrozen      = 0 !< frozen particles (state > 0 and updatecnt > 0 in the frozen project)
    integer, allocatable :: nfrozen_state(:)
    character(len=8)  :: weighting = '' !< the loading run's reconstruction weighting (euclid|cc)
  contains
    ! the run context
    procedure          :: new
    procedure          :: write
    procedure          :: read
    procedure          :: validate
    procedure          :: load
    procedure          :: get_nfrozen_state
    ! gridding sets
    procedure          :: write_gridding_set
    procedure          :: gridding_set_status
    procedure          :: add_gridding_set
    ! PCG halves
    procedure, private :: pcg_provenance
    procedure          :: write_pcg_half
    procedure          :: pcg_half_status
    procedure          :: add_pcg_half
    procedure          :: kill
end type frozen_accum

contains

    ! RUN CONTEXT

    !> the store of one add-on run: nrows and nrows_frozen are the row counts
    !! of the working and the frozen project, nstates is the length of the
    !! per-state frozen counts and nfrozen their sum
    subroutine new( self, run_id, backend, nrows, nrows_frozen, nfrozen_state )
        class(frozen_accum), intent(inout) :: self
        character(len=*),    intent(in)    :: run_id, backend
        integer,             intent(in)    :: nrows, nrows_frozen
        integer,             intent(in)    :: nfrozen_state(:)
        call self%kill
        if( len_trim(run_id) == 0 .or. index(trim(run_id), ' ') > 0 ) &
            &THROW_HARD('frozen context requires a blank-free run identifier')
        if( size(nfrozen_state) < 1 ) THROW_HARD('frozen context has no state layout')
        if( nrows < 1 .or. any(nfrozen_state < 0) ) THROW_HARD('frozen context particle counts are invalid')
        if( nrows_frozen < 1 ) THROW_HARD('frozen context row counts are invalid')
        self%run_id        = run_id
        self%backend       = backend
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
        write(funit,'(A,1X,A)')  'backend', trim(self%backend)
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
        character(len=64) :: key, schema
        integer :: funit, io_stat, version, s, cnt, nseen
        logical :: l_end, l_seen(6)
        call self%kill
        status = 1
        msg    = ''
        l_end  = .false.
        l_seen = .false.
        nseen  = 0
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
                    read(line,*,iostat=io_stat) key, self%backend
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

    !> the context must describe the loading run: backend, state layout and
    !! particle index space, which is the frozen project's for a producer
    !! (a reconstruct3D on the frozen project) and the working project's for a
    !! consumer (an add-on reconstruction)
    subroutine validate( self, backend, nstates, nrows, producer, status, msg )
        class(frozen_accum), intent(in)  :: self
        character(len=*),    intent(in)  :: backend
        integer,             intent(in)  :: nstates, nrows
        logical,             intent(in)  :: producer
        integer,             intent(out) :: status
        character(len=*),    intent(out) :: msg
        status = 1
        msg    = ''
        if( trim(self%backend) /= trim(backend) )then
            msg = 'frozen context backend '//trim(self%backend)//' differs from the reconstruction backend '//trim(backend)
        else if( self%nstates /= nstates )then
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
    subroutine load( self, fname, backend, nstates, nrows, cc_objfun, producer )
        class(frozen_accum), intent(inout) :: self
        class(string),       intent(in)    :: fname
        character(len=*),    intent(in)    :: backend
        integer,             intent(in)    :: nstates, nrows, cc_objfun
        logical,             intent(in)    :: producer
        character(len=STDLEN) :: msg
        integer :: status
        call self%read(fname, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        call self%validate(backend, nstates, nrows, producer, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        self%weighting = trim(merge('euclid', 'cc    ', cc_objfun == OBJFUN_EUCLID))
    end subroutine load

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
        self%backend = ''
        self%nstates = 0
        self%nrows   = 0
        self%nrows_frozen = 0
        self%nfrozen = 0
        if( allocated(self%nfrozen_state) ) deallocate(self%nfrozen_state)
        self%weighting = ''
    end subroutine kill

    ! GRIDDING SETS

    function gridding_component( state, box, ifile ) result( fname )
        integer, intent(in) :: state, box, ifile
        type(string) :: fname, fbody
        fbody = refine3D_frozen_rec_fbody(state, box)
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
        call fbody%kill
    end function gridding_component

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
        type(string)    :: manifest, tmpname
        integer(kind=8) :: sizes(4)
        integer :: ldim(3), funit, io_stat, ifile
        real    :: smpd
        if( state < 1 .or. state > self%nstates ) THROW_HARD('frozen gridding set state is outside the context layout')
        ldim = even_rec%get_ldim()
        smpd = even_rec%get_smpd()
        if( any(odd_rec%get_ldim() /= ldim) ) THROW_HARD('frozen gridding pair dimensions do not match')
        manifest = refine3D_frozen_manifest_fname(state, ldim(1))
        call del_file(manifest)
        call even_rec%write_raw_accum(gridding_component(state, ldim(1), 1), gridding_component(state, ldim(1), 2))
        call odd_rec%write_raw_accum( gridding_component(state, ldim(1), 3), gridding_component(state, ldim(1), 4))
        do ifile = 1, 4
            sizes(ifile) = file_size(gridding_component(state, ldim(1), ifile))
        enddo
        tmpname = manifest//'.tmp'
        call del_file(tmpname)
        call fopen(funit, file=tmpname, status='REPLACE', action='WRITE', iostat=io_stat)
        call fileiochk('frozen_accum%write_gridding_set opening '//tmpname%to_char(), io_stat)
        write(funit,'(A,1X,I0)')     FROZEN_SET_SCHEMA, FROZEN_SCHEMA_VERSION
        write(funit,'(A,1X,A)')      'run_id',  trim(self%run_id)
        write(funit,'(A,1X,A)')      'backend', 'gridding'
        write(funit,'(A,1X,I0)')     'state',   state
        write(funit,'(A,1X,I0)')     'nstates', self%nstates
        write(funit,'(A,1X,I0)')     'box',     ldim(1)
        write(funit,'(A,1X,ES16.9)') 'smpd',    smpd
        write(funit,'(A,1X,I0)')     'nrows',   self%nrows
        write(funit,'(A,1X,I0)')     'nfrozen_state', self%nfrozen_state(state)
        write(funit,'(A,1X,A)')      'weighting', trim(merge(self%weighting, 'none    ', len_trim(self%weighting) > 0))
        write(funit,'(A,4(1X,I0))')  'sizes',   sizes
        write(funit,'(A)')           'end'
        call fclose(funit)
        call simple_rename(tmpname, manifest, overwrite=.true.)
        write(logfhandle,'(A,I0,A,I0)') '>>> FROZEN GRIDDING SET WRITTEN, STATE ', state, ', BOX ', ldim(1)
        call manifest%kill
        call tmpname%kill
    end subroutine write_gridding_set

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
        integer(kind=8)    :: sizes(4)
        integer :: funit, io_stat, version, state_set, nstates_set, box_set, nrows_set, nfrozen_set, ifile
        real    :: smpd_set
        logical :: l_end
        status = 1
        msg    = ''
        l_end  = .false.
        run_id = ''; backend = ''; weighting = ''
        state_set = -1; nstates_set = -1; box_set = -1; nrows_set = -1; nfrozen_set = -1; smpd_set = -1.; sizes = -2
        manifest = refine3D_frozen_manifest_fname(state, box)
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
        else if( trim(backend) /= 'gridding' )then
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
        else
            status = 0
            do ifile = 1, 4
                if( file_size(gridding_component(state, box, ifile)) /= sizes(ifile) )then
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
            fname = gridding_component(state, ldim(1), ifile)
            call ioimg%open(fname, ldim, smpd, formatchar='M', readhead=.false., rwaction='read')
            call rec%read_raw_mrc(ioimg)
            call ioimg%close
            fname = gridding_component(state, ldim(1), ifile+1)
            call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
            call fileiochk('frozen_accum%add_gridding_set opening '//fname%to_char(), ierr)
            call rec%read_raw_rho(funit)
            call fclose(funit)
        end subroutine read_one

    end subroutine add_gridding_set

    ! PCG HALVES

    !> provenance tag of a PCG frozen raw half: binds it to the add-on run, its
    !! particle index space and the reconstruction weighting it was accumulated with
    function pcg_provenance( self ) result( provenance )
        class(frozen_accum), intent(in) :: self
        character(len=256) :: provenance
        write(provenance,'(A,A,A,I0,A,I0,A,A)') 'frozen run=', trim(self%run_id), ' rows=', self%nrows, &
            &' nstates=', self%nstates, ' weighting=', trim(merge(self%weighting, 'none    ', len_trim(self%weighting) > 0))
    end function pcg_provenance

    !> Publish the open raw accumulator of one frozen (state,half) at box; the
    !! raw writer is atomic (.tmp then rename)
    subroutine write_pcg_half( self, state, eo, box, pcgop, nptcls )
        class(frozen_accum),      intent(in) :: self
        integer,                  intent(in) :: state, eo, box, nptcls
        class(reconstructor_pcg), intent(in) :: pcgop
        type(string) :: fname
        if( state < 1 .or. state > self%nstates ) THROW_HARD('frozen PCG half state is outside the context layout')
        fname = refine3D_frozen_pcg_fname(state, box, merge('odd ', 'even', eo == 1))
        call pcgop%write_raw_accum(fname, state, eo, 1, 1, nptcls, trim(self%pcg_provenance()))
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FROZEN PCG HALF WRITTEN, STATE ', state, ', HALF ', eo, &
            &', BOX ', box
        call fname%kill
    end subroutine write_pcg_half

    !> Validate the PCG frozen half a consumer at (box, smpd) would read: exact
    !! box and sampling, identity and the context's provenance tag
    subroutine pcg_half_status( self, state, eo, box, smpd, status, msg )
        class(frozen_accum), intent(in)  :: self
        integer,             intent(in)  :: state, eo, box
        real,                intent(in)  :: smpd
        integer,             intent(out) :: status
        character(len=*),    intent(out) :: msg
        character(len=256) :: provenance
        type(string) :: fname
        integer :: state_file, eo_file, part_file, nparts_file, nptcls_file, box_file, io_stat
        real    :: smpd_file
        status = 1
        msg    = ''
        fname  = refine3D_frozen_pcg_fname(state, box, merge('odd ', 'even', eo == 1))
        if( .not. file_exists(fname) )then
            msg = 'frozen PCG half is missing: '//fname%to_char()
            call fname%kill
            return
        endif
        call read_pcg_raw_accum_header(fname, state_file, eo_file, part_file, nparts_file, nptcls_file, &
            &box_file, smpd_file, provenance, io_stat)
        call fname%kill
        if( io_stat /= 0 )then
            msg = 'frozen PCG half is unreadable or not in the raw accumulator format'
        else if( trim(provenance) /= trim(self%pcg_provenance()) )then
            msg = 'frozen PCG half belongs to another add-on run or weighting'
        else if( state_file /= state .or. eo_file /= eo .or. part_file /= 1 .or. nparts_file /= 1 )then
            msg = 'frozen PCG half identity mismatch'
        else if( box_file /= box )then
            msg = 'frozen PCG half box differs from the consuming grid'
        else if( abs(smpd_file - smpd) > SMPD_RELTOL * smpd )then
            msg = 'frozen PCG half sampling differs from the consuming grid'
        else if( nptcls_file < 0 )then
            msg = 'frozen PCG half has an invalid particle count'
        else
            status = 0
        endif
    end subroutine pcg_half_status

    !> Add the validated frozen half into an open raw reduction with weight
    !! one; nptcls returns the frozen particle count of the half. A missing or
    !! invalid half is fatal.
    subroutine add_pcg_half( self, state, eo, box, smpd, pcgop, nptcls )
        class(frozen_accum),      intent(in)    :: self
        integer,                  intent(in)    :: state, eo, box
        real,                     intent(in)    :: smpd
        class(reconstructor_pcg), intent(inout) :: pcgop
        integer,                  intent(out)   :: nptcls
        character(len=STDLEN) :: msg
        type(string) :: fname
        integer :: status
        call self%pcg_half_status(state, eo, box, smpd, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        fname = refine3D_frozen_pcg_fname(state, box, merge('odd ', 'even', eo == 1))
        call pcgop%add_raw_accum_weighted(fname, state, eo, 1, 1, trim(self%pcg_provenance()), 1.0, nptcls)
        call fname%kill
    end subroutine add_pcg_half

end module simple_frozen_accum

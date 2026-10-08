!@descr: state_weight_set: a transactional set of per-state particle weights registered in the project
! One file per state and a manifest published last (simple_state_weights_file); the project's out
! segment points at the manifest (imgkind state_weights). `new` validates the whole set once against
! the current project; afterwards one state column is held at a time. Applied weights are the raw
! weights with every value at or below the consumer's threshold set to zero.
module simple_state_weight_set
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int8, int64, real32, real64
use simple_error,              only: simple_exception
use simple_fileio,             only: get_fpath
use simple_oris,               only: oris
use simple_ptcl_layout,        only: ptcl_layout_digest
use simple_sp_project,         only: sp_project
use simple_string,             only: string
use simple_string_utils,       only: int2str
use simple_syslib,             only: file_exists, del_file, simple_abspath
use simple_state_weights_file, only: state_weights_manifest, state_weights_file_header, state_weights_fname, &
    &state_weights_write_state_file, state_weights_read_state_file, state_weights_file_checksum, &
    &state_weights_write_manifest, state_weights_read_manifest, STATE_WEIGHTS_MANIFEST_FNAME, &
    &STATE_WEIGHTS_KIND_PARTITION, STATE_WEIGHTS_KIND_KERNEL
implicit none
private
#include "simple_local_flags.inc"

public :: state_weight_set, discard_state_weight_set, STATE_WEIGHTS_KIND_PARTITION, STATE_WEIGHTS_KIND_KERNEL

!> a weighted PARTITION row sums to one across the states within this
real(real64), parameter :: ROW_SUM_TOL    = 1.0e-3_real64
!> relative agreement of the populations recomputed from the files with the manifest
real(real64), parameter :: POPULATION_TOL = 1.0e-9_real64

type :: state_weight_set
    private
    type(state_weights_manifest) :: manifest
    type(string)                 :: manifest_path
    type(string)                 :: dir
    integer                      :: loaded = 0           !< the state whose column is held
    real(real32),  allocatable   :: column(:)
    integer(int8), allocatable   :: flags(:)
    logical                      :: exists = .false.
  contains
    procedure :: new
    procedure :: publish
    procedure :: publish_work_state
    procedure :: get_nstates
    procedure :: get_nptcls
    procedure :: get_kind
    procedure :: get_producer
    procedure :: get_generation
    procedure :: get_layout_digest
    procedure :: get_manifest_path
    procedure :: get_parent_generation
    procedure :: get_parent_layout_digest
    procedure :: get_parent_state
    procedure :: get_mass
    procedure :: get_ess
    procedure :: get_pop
    procedure :: get_weights
    procedure :: get_applied_weights
    procedure :: get_labels
    procedure :: get_selection
    procedure :: applied_mass
    procedure :: applied_ess
    procedure :: get_members
    procedure :: is_fractional
    procedure :: kill
    procedure, private :: load_state
end type state_weight_set

contains

    !> Open the set the project registers and validate it completely against the project and its
    !! particle field: manifest, file sizes and checksums, file headers, layout digest, populations.
    !! Without status, any failure stops the program.
    subroutine new( self, project, particles, status, message )
        class(state_weight_set),    intent(inout) :: self
        class(sp_project),          intent(in)    :: project
        class(oris),                intent(in)    :: particles
        integer,          optional, intent(out)   :: status
        character(len=*), optional, intent(out)   :: message
        type(string)       :: path
        character(len=512) :: msg
        integer            :: stat
        logical            :: found
        call self%kill
        call project%get_state_weights(path, found)
        if( found )then
            call open_set(self, path, project, particles, stat, msg)
        else
            stat = 1
            msg  = 'the project registers no state weight set'
        endif
        if( stat /= 0 ) call self%kill
        if( present(status) ) status = stat
        if( present(message) ) message = msg
        if( stat /= 0 .and. .not. present(status) ) THROW_HARD(trim(msg))
        call path%kill
    end subroutine new

    !> Publish a new generation in the working directory and register it in the project:
    !! state files, then the manifest (temporary name and rename), then the project pointer
    !! (written to projfile), then the deletion of the older generation's files in this directory.
    !! weights_sel(nsel,nstates) and labels_sel(nsel) belong to the project rows pinds(nsel); every other
    !! row has zero weight. The set is validated before anything is published; self then holds it.
    !! The kind is inferred from the rows unless given; a derived set records its parent's identity.
    subroutine publish( self, project, particles, projfile, pinds, weights_sel, labels_sel, producer, kind, &
        &parent_generation, parent_layout_digest, parent_state )
        class(state_weight_set),  intent(inout) :: self
        class(sp_project),        intent(inout) :: project
        class(oris),              intent(in)    :: particles
        class(string),            intent(in)    :: projfile
        integer,                  intent(in)    :: pinds(:), labels_sel(:)
        real,                     intent(in)    :: weights_sel(:,:)
        character(len=*),         intent(in)    :: producer
        integer,        optional, intent(in)    :: kind
        integer(int64), optional, intent(in)    :: parent_generation, parent_layout_digest
        integer,        optional, intent(in)    :: parent_state
        type(state_weights_manifest)    :: manifest, previous
        type(state_weights_file_header) :: header
        type(string)                    :: prev_path, prev_dir, fname, manifest_path
        real(real32),  allocatable      :: column(:)
        integer(int8), allocatable      :: flags(:)
        real(real64),  allocatable      :: rowsum(:)
        logical,       allocatable      :: seen(:)
        character(len=512) :: msg
        integer(int64)     :: digest
        integer            :: nptcls, nsel, nstates, i, s, stat
        logical            :: l_prev, l_partition
        call self%kill
        nptcls  = particles%get_noris()
        nsel    = size(pinds)
        nstates = size(weights_sel,2)
        if( nptcls < 1 .or. nsel < 1 .or. nstates < 1 ) THROW_HARD('empty state weight set; publish')
        ! the manifest stores the producer as one list-directed token
        if( len_trim(producer) == 0 .or. len_trim(producer) > 64 .or. scan(trim(producer), ' ,/') > 0 )then
            THROW_HARD('the producer must be one word of at most 64 characters; publish')
        endif
        if( size(weights_sel,1) /= nsel .or. size(labels_sel) /= nsel ) THROW_HARD('state weight table dimensions disagree; publish')
        digest = ptcl_layout_digest(project, particles)
        if( digest == 0_int64 ) THROW_HARD('the particle layout digest is undefined for this project; publish')
        ! rows: inside the project, distinct, finite in [0,1], labels pointing at positive weights
        allocate(seen(nptcls), source=.false.)
        allocate(rowsum(nptcls), source=0._real64)
        do i = 1, nsel
            if( pinds(i) < 1 .or. pinds(i) > nptcls ) THROW_HARD('state weight row outside the project; publish')
            if( seen(pinds(i)) ) THROW_HARD('duplicate project row in the state weight set; publish')
            seen(pinds(i)) = .true.
            if( .not. all(ieee_is_finite(weights_sel(i,:))) ) THROW_HARD('non-finite state weight; publish')
            if( any(weights_sel(i,:) < 0.) .or. any(weights_sel(i,:) > 1.) ) THROW_HARD('state weight outside [0,1]; publish')
            if( labels_sel(i) < 0 .or. labels_sel(i) > nstates ) THROW_HARD('state label outside the state range; publish')
            if( labels_sel(i) > 0 )then
                if( weights_sel(i,labels_sel(i)) <= 0. ) THROW_HARD('a state label points at a zero weight; publish')
            endif
            rowsum(pinds(i)) = sum(real(real(weights_sel(i,:),real32),real64))
        enddo
        l_partition = .true.
        do i = 1, nptcls
            if( rowsum(i) <= 0._real64 ) cycle
            if( abs(rowsum(i) - 1._real64) > ROW_SUM_TOL )then
                l_partition = .false.
                exit
            endif
        enddo
        deallocate(seen, rowsum)
        ! the generation follows the one the project registers and any manifest already in this directory
        manifest%generation = 1_int64
        call project%get_state_weights(prev_path, l_prev)
        if( l_prev )then
            call state_weights_read_manifest(prev_path, previous, stat, msg)
            if( stat == 0 ) manifest%generation = previous%generation + 1_int64
        endif
        manifest_path = simple_abspath(STATE_WEIGHTS_MANIFEST_FNAME, check_exists=.false.)
        if( file_exists(STATE_WEIGHTS_MANIFEST_FNAME) )then
            call state_weights_read_manifest(string(STATE_WEIGHTS_MANIFEST_FNAME), previous, stat, msg)
            if( stat == 0 )then
                manifest%generation = max(manifest%generation, previous%generation + 1_int64)
                ! the files this publication supersedes in this directory
                prev_path = manifest_path
                l_prev    = .true.
            endif
        endif
        if( l_prev )then
            call state_weights_read_manifest(prev_path, previous, stat, msg)
            l_prev = stat == 0 .and. prev_path == manifest_path
        endif
        manifest%kind          = merge(STATE_WEIGHTS_KIND_PARTITION, STATE_WEIGHTS_KIND_KERNEL, l_partition)
        if( present(kind) )then
            if( kind /= STATE_WEIGHTS_KIND_PARTITION .and. kind /= STATE_WEIGHTS_KIND_KERNEL ) &
                &THROW_HARD('unknown state weight kind; publish')
            manifest%kind = kind
        endif
        if( present(parent_generation)    ) manifest%parent_generation    = parent_generation
        if( present(parent_layout_digest) ) manifest%parent_layout_digest = parent_layout_digest
        if( present(parent_state)         ) manifest%parent_state         = parent_state
        manifest%producer      = producer
        manifest%nstates       = nstates
        manifest%nptcls        = nptcls
        manifest%layout_digest = digest
        allocate(manifest%fnames(nstates), manifest%nbytes(nstates), manifest%checksums(nstates), &
            &manifest%mass(nstates), manifest%ess(nstates), manifest%pop(nstates))
        allocate(column(nptcls), flags(nptcls))
        do s = 1, nstates
            column = 0._real32
            flags  = 0_int8
            do i = 1, nsel
                column(pinds(i)) = real(weights_sel(i,s), real32)
                if( labels_sel(i) == s ) flags(pinds(i)) = 1_int8
            enddo
            fname = state_weights_fname(manifest%generation, s)
            header%nptcls        = int(nptcls, int64)
            header%nstates       = int(nstates, int64)
            header%state         = int(s, int64)
            header%generation    = manifest%generation
            header%layout_digest = digest
            header%kind          = int(manifest%kind, int64)
            call state_weights_write_state_file(fname, header, column, flags, stat, msg)
            if( stat /= 0 ) THROW_HARD(trim(msg))
            manifest%fnames(s)    = fname
            manifest%nbytes(s)    = header%file_bytes
            manifest%checksums(s) = state_weights_file_checksum(fname)
            call column_populations(column, flags, manifest%mass(s), manifest%ess(s), manifest%pop(s))
        enddo
        deallocate(column, flags)
        call state_weights_write_manifest(string(STATE_WEIGHTS_MANIFEST_FNAME), manifest, stat, msg)
        if( stat /= 0 ) THROW_HARD(trim(msg))
        call project%add_state_weights2os_out(manifest_path)
        call project%write_segment_inside('out', projfile)
        ! only now may the superseded generation of this directory go
        if( l_prev )then
            prev_dir = get_fpath(prev_path)
            do s = 1, previous%nstates
                if( previous%generation == manifest%generation ) exit
                fname = prev_dir//previous%fnames(s)
                call del_file(fname)
            enddo
        endif
        call open_set(self, manifest_path, project, particles, stat, msg)
        if( stat /= 0 ) THROW_HARD('the published state weight set does not validate: '//trim(msg))
        call fname%kill
        call prev_path%kill
        call prev_dir%kill
        call manifest_path%kill
    end subroutine publish

    !> The single-state set of a work project derived from one state of parent (as refine3D_auto
    !! state=X builds it): work row i is parent row parent_rows(i); its weight and hard-label flag are
    !! the parent's for that state. The set keeps the parent's producer and kind and records the
    !! parent's identity and state. It is published in the working directory and registered in
    !! projfile, which must not be where the parent's own manifest lives.
    subroutine publish_work_state( self, parent, state, parent_rows, project, particles, projfile )
        class(state_weight_set), intent(inout) :: self
        class(state_weight_set), intent(inout) :: parent
        integer,                 intent(in)    :: state, parent_rows(:)
        class(sp_project),       intent(inout) :: project
        class(oris),             intent(in)    :: particles
        class(string),           intent(in)    :: projfile
        type(string)         :: target_path
        real,    allocatable :: weights(:,:)
        integer, allocatable :: labels(:), pinds(:)
        integer :: nwork, i
        call check_state(parent, state)
        nwork = size(parent_rows)
        if( nwork < 1 .or. nwork /= particles%get_noris() ) THROW_HARD('work rows and work particles disagree; publish_work_state')
        if( any(parent_rows < 1) .or. any(parent_rows > parent%manifest%nptcls) ) &
            &THROW_HARD('work row mapped outside the parent set; publish_work_state')
        ! the work set's manifest would replace the parent's
        target_path = simple_abspath(STATE_WEIGHTS_MANIFEST_FNAME, check_exists=.false.)
        if( target_path == parent%manifest_path ) &
            &THROW_HARD('a work project''s state weight set cannot be published beside its parent''s; publish_work_state')
        call parent%load_state(state)
        allocate(weights(nwork,1), labels(nwork), pinds(nwork))
        do i = 1, nwork
            pinds(i)     = i
            weights(i,1) = parent%column(parent_rows(i))
            labels(i)    = merge(1, 0, parent%flags(parent_rows(i)) == 1_int8)
        enddo
        call self%publish(project, particles, projfile, pinds, weights, labels, parent%get_producer(), &
            &kind=parent%manifest%kind, parent_generation=parent%manifest%generation, &
            &parent_layout_digest=parent%manifest%layout_digest, parent_state=state)
        deallocate(weights, labels, pinds)
        call target_path%kill
    end subroutine publish_work_state

    ! ---- getters ----

    integer function get_nstates( self )
        class(state_weight_set), intent(in) :: self
        get_nstates = self%manifest%nstates
    end function get_nstates

    integer function get_nptcls( self )
        class(state_weight_set), intent(in) :: self
        get_nptcls = self%manifest%nptcls
    end function get_nptcls

    !> STATE_WEIGHTS_KIND_PARTITION or STATE_WEIGHTS_KIND_KERNEL
    integer function get_kind( self )
        class(state_weight_set), intent(in) :: self
        get_kind = self%manifest%kind
    end function get_kind

    function get_producer( self ) result( producer )
        class(state_weight_set), intent(in) :: self
        character(len=:), allocatable :: producer
        producer = trim(self%manifest%producer)
    end function get_producer

    !> with the layout digest, the identity of the set
    integer(int64) function get_generation( self )
        class(state_weight_set), intent(in) :: self
        get_generation = self%manifest%generation
    end function get_generation

    integer(int64) function get_layout_digest( self )
        class(state_weight_set), intent(in) :: self
        get_layout_digest = self%manifest%layout_digest
    end function get_layout_digest

    function get_manifest_path( self ) result( path )
        class(state_weight_set), intent(in) :: self
        type(string) :: path
        path = self%manifest_path
    end function get_manifest_path

    !> generation of the set this one was derived from; 0 for a set with no parent
    integer(int64) function get_parent_generation( self )
        class(state_weight_set), intent(in) :: self
        get_parent_generation = self%manifest%parent_generation
    end function get_parent_generation

    integer(int64) function get_parent_layout_digest( self )
        class(state_weight_set), intent(in) :: self
        get_parent_layout_digest = self%manifest%parent_layout_digest
    end function get_parent_layout_digest

    !> the parent's state this set's single column is; 0 for a set with no parent
    integer function get_parent_state( self )
        class(state_weight_set), intent(in) :: self
        get_parent_state = self%manifest%parent_state
    end function get_parent_state

    !> applied mass of the raw weights (sum over all rows)
    real(real64) function get_mass( self, state )
        class(state_weight_set), intent(in) :: self
        integer,                 intent(in) :: state
        call check_state(self, state)
        get_mass = self%manifest%mass(state)
    end function get_mass

    !> effective sample size of the raw weights, mass**2 / sum of squared weights
    real(real64) function get_ess( self, state )
        class(state_weight_set), intent(in) :: self
        integer,                 intent(in) :: state
        call check_state(self, state)
        get_ess = self%manifest%ess(state)
    end function get_ess

    !> hard population: rows whose label is this state
    integer function get_pop( self, state )
        class(state_weight_set), intent(in) :: self
        integer,                 intent(in) :: state
        call check_state(self, state)
        get_pop = self%manifest%pop(state)
    end function get_pop

    !> raw weights of one state over every project row
    subroutine get_weights( self, state, weights )
        class(state_weight_set),   intent(inout) :: self
        integer,                   intent(in)    :: state
        real(real32), allocatable, intent(out)   :: weights(:)
        call self%load_state(state)
        weights = self%column
    end subroutine get_weights

    !> weights of one state with every value at or below threshold set to zero
    subroutine get_applied_weights( self, state, threshold, weights )
        class(state_weight_set),   intent(inout) :: self
        integer,                   intent(in)    :: state
        real,                      intent(in)    :: threshold
        real(real32), allocatable, intent(out)   :: weights(:)
        call self%load_state(state)
        weights = merge(self%column, 0._real32, self%column > threshold)
    end subroutine get_applied_weights

    !> hard labels over every project row: the state whose flag is set, 0 for none
    subroutine get_labels( self, labels )
        class(state_weight_set), intent(inout) :: self
        integer, allocatable,    intent(out)   :: labels(:)
        integer :: s
        allocate(labels(self%manifest%nptcls), source=0)
        do s = 1, self%manifest%nstates
            call self%load_state(s)
            where( self%flags == 1_int8 ) labels = s
        enddo
    end subroutine get_labels

    !> selected rows: those whose applied weights sum above zero
    subroutine get_selection( self, threshold, selected )
        class(state_weight_set), intent(inout) :: self
        real,                    intent(in)    :: threshold
        logical, allocatable,    intent(out)   :: selected(:)
        integer :: s
        allocate(selected(self%manifest%nptcls), source=.false.)
        do s = 1, self%manifest%nstates
            call self%load_state(s)
            selected = selected .or. self%column > threshold
        enddo
    end subroutine get_selection

    !> applied mass of one state over the given project rows
    real(real64) function applied_mass( self, state, threshold, rows )
        class(state_weight_set), intent(inout) :: self
        integer,                 intent(in)    :: state
        real,                    intent(in)    :: threshold
        integer,                 intent(in)    :: rows(:)
        real(real64) :: ess
        integer      :: pop
        call self%load_state(state)
        call row_populations(self%column, threshold, rows, applied_mass, ess, pop)
    end function applied_mass

    !> effective sample size of one state's applied weights over the given project rows
    real(real64) function applied_ess( self, state, threshold, rows )
        class(state_weight_set), intent(inout) :: self
        integer,                 intent(in)    :: state
        real,                    intent(in)    :: threshold
        integer,                 intent(in)    :: rows(:)
        real(real64) :: mass
        integer      :: pop
        call self%load_state(state)
        call row_populations(self%column, threshold, rows, mass, applied_ess, pop)
    end function applied_ess

    !> The rows of a subset whose weight for state is above threshold, in subset order, and their weights
    subroutine get_members( self, state, threshold, rows, members, weights )
        class(state_weight_set), intent(inout) :: self
        integer,                 intent(in)    :: state
        real,                    intent(in)    :: threshold
        integer,                 intent(in)    :: rows(:)
        integer, allocatable,    intent(out)   :: members(:)
        real,    allocatable,    intent(out)   :: weights(:)
        logical, allocatable :: l_member(:)
        integer :: i
        call self%load_state(state)
        allocate(l_member(size(rows)))
        do i = 1, size(rows)
            if( rows(i) < 1 .or. rows(i) > self%manifest%nptcls ) THROW_HARD('row outside the state weight set')
            l_member(i) = self%column(rows(i)) > threshold
        enddo
        members = pack(rows, l_member)
        allocate(weights(size(members)))
        do i = 1, size(members)
            weights(i) = real(self%column(members(i)))
        enddo
        deallocate(l_member)
    end subroutine get_members

    !> Whether some weight of the state lies strictly between zero and one (else the column is 0/1)
    logical function is_fractional( self, state )
        class(state_weight_set), intent(inout) :: self
        integer,                 intent(in)    :: state
        call self%load_state(state)
        is_fractional = any(self%column > 0._real32 .and. self%column < 1._real32)
    end function is_fractional

    subroutine kill( self )
        class(state_weight_set), intent(inout) :: self
        if( allocated(self%manifest%fnames)    ) deallocate(self%manifest%fnames)
        if( allocated(self%manifest%nbytes)    ) deallocate(self%manifest%nbytes)
        if( allocated(self%manifest%checksums) ) deallocate(self%manifest%checksums)
        if( allocated(self%manifest%mass)      ) deallocate(self%manifest%mass)
        if( allocated(self%manifest%ess)       ) deallocate(self%manifest%ess)
        if( allocated(self%manifest%pop)       ) deallocate(self%manifest%pop)
        self%manifest = state_weights_manifest()
        if( allocated(self%column) ) deallocate(self%column)
        if( allocated(self%flags)  ) deallocate(self%flags)
        call self%manifest_path%kill
        call self%dir%kill
        self%loaded = 0
        self%exists = .false.
    end subroutine kill

    !> Withdraw the set the project registers, for a consumer that turns the weights into hard labels and
    !! keeps no use for them: the out-segment entry goes first (and the segment is written), then the state
    !! files and the manifest. A project without a set is left alone.
    subroutine discard_state_weight_set( project, projfile )
        class(sp_project), intent(inout) :: project
        class(string),     intent(in)    :: projfile
        type(state_weights_manifest) :: manifest
        type(string)       :: path, dir
        character(len=512) :: msg
        integer            :: s, stat
        logical            :: found
        call project%get_state_weights(path, found)
        if( .not. found ) return
        call project%remove_state_weights_from_osout
        if( project%os_out%get_noris() > 0 )then
            call project%write_segment_inside('out', projfile)
        else
            ! the set was the segment's only entry: empty the segment on disk
            call project%clear_segment_inside('out', projfile)
        endif
        call state_weights_read_manifest(path, manifest, stat, msg)
        if( stat == 0 )then
            dir = get_fpath(path)
            do s = 1, manifest%nstates
                call del_file(dir//manifest%fnames(s))
            enddo
            call dir%kill
        endif
        if( file_exists(path) ) call del_file(path)
        call path%kill
    end subroutine discard_state_weight_set

    ! ---- private ----

    !> hold the column of one state (read on demand; the set was validated by new)
    subroutine load_state( self, state )
        class(state_weight_set), intent(inout) :: self
        integer,                 intent(in)    :: state
        type(state_weights_file_header) :: header
        character(len=512) :: msg
        integer :: stat
        call check_state(self, state)
        if( self%loaded == state ) return
        call state_weights_read_state_file(self%dir//self%manifest%fnames(state), header, self%column, &
            &self%flags, stat, msg)
        if( stat /= 0 ) THROW_HARD(trim(msg))
        if( header%generation /= self%manifest%generation ) THROW_HARD('a state weight file changed after validation')
        self%loaded = state
    end subroutine load_state

    subroutine check_state( self, state )
        class(state_weight_set), intent(in) :: self
        integer,                 intent(in) :: state
        if( .not. self%exists ) THROW_HARD('state weight set not opened')
        if( state < 1 .or. state > self%manifest%nstates ) THROW_HARD('state out of range of the state weight set')
    end subroutine check_state

    !> read the manifest at path and validate every file of it against the project
    subroutine open_set( self, path, project, particles, status, message )
        class(state_weight_set), intent(inout) :: self
        class(string),           intent(in)    :: path
        class(sp_project),       intent(in)    :: project
        class(oris),             intent(in)    :: particles
        integer,                 intent(out)   :: status
        character(len=*),        intent(out)   :: message
        type(state_weights_file_header) :: header
        type(string)               :: fpath
        real(real32),  allocatable :: column(:)
        integer(int8), allocatable :: flags(:)
        real(real64),  allocatable :: rowsum(:)
        integer,       allocatable :: nflag(:)
        real(real64)   :: mass, ess
        integer(int64) :: digest, fsize
        integer        :: s, pop, i
        call self%kill
        call state_weights_read_manifest(path, self%manifest, status, message)
        if( status /= 0 ) return
        status = 1
        if( self%manifest%nptcls /= particles%get_noris() )then
            message = 'the state weight set holds '//int2str(self%manifest%nptcls)//' rows, the project '//&
                &int2str(particles%get_noris()); return
        endif
        digest = ptcl_layout_digest(project, particles)
        if( digest == 0_int64 .or. digest /= self%manifest%layout_digest )then
            message = 'the state weight set belongs to another particle layout'; return
        endif
        self%dir = get_fpath(path)
        allocate(rowsum(self%manifest%nptcls), source=0._real64)
        allocate(nflag(self%manifest%nptcls),  source=0)
        do s = 1, self%manifest%nstates
            fpath = self%dir//self%manifest%fnames(s)
            if( .not. file_exists(fpath) )then
                message = 'missing state weight file '//fpath%to_char(); return
            endif
            inquire(file=fpath%to_char(), size=fsize)
            if( fsize /= self%manifest%nbytes(s) )then
                message = 'state weight file has the wrong size: '//fpath%to_char(); return
            endif
            if( state_weights_file_checksum(fpath) /= self%manifest%checksums(s) )then
                message = 'state weight file checksum mismatch: '//fpath%to_char(); return
            endif
            call state_weights_read_state_file(fpath, header, column, flags, status, message)
            if( status /= 0 ) return
            status = 1
            if( header%generation /= self%manifest%generation .or. header%state /= s .or. &
                &header%nstates /= self%manifest%nstates .or. header%nptcls /= self%manifest%nptcls .or. &
                &header%layout_digest /= self%manifest%layout_digest .or. header%kind /= self%manifest%kind )then
                message = 'state weight file belongs to another delivery: '//fpath%to_char(); return
            endif
            if( .not. all(ieee_is_finite(column)) .or. any(column < 0.) .or. any(column > 1.) )then
                message = 'state weight file holds values outside [0,1]: '//fpath%to_char(); return
            endif
            if( any(flags == 1_int8 .and. column <= 0.) .or. any(flags /= 0_int8 .and. flags /= 1_int8) )then
                message = 'state weight file has invalid label flags: '//fpath%to_char(); return
            endif
            call column_populations(column, flags, mass, ess, pop)
            if( .not. close_enough(mass, self%manifest%mass(s)) .or. .not. close_enough(ess, self%manifest%ess(s)) &
                &.or. pop /= self%manifest%pop(s) )then
                message = 'state weight populations disagree with the manifest: '//fpath%to_char(); return
            endif
            rowsum = rowsum + real(column, real64)
            nflag  = nflag  + int(flags)
        enddo
        if( any(nflag > 1) )then
            message = 'a particle carries the label of more than one state'; return
        endif
        ! a derived set holds one column of its parent: its rows sum to one across the parent's states
        if( self%manifest%kind == STATE_WEIGHTS_KIND_PARTITION .and. self%manifest%parent_generation == 0_int64 )then
            do i = 1, self%manifest%nptcls
                if( rowsum(i) <= 0._real64 ) cycle
                if( abs(rowsum(i) - 1._real64) > ROW_SUM_TOL )then
                    message = 'a PARTITION row does not sum to one across the states'; return
                endif
            enddo
        endif
        self%manifest_path = path
        self%loaded        = self%manifest%nstates
        call move_alloc(column, self%column)
        call move_alloc(flags,  self%flags)
        self%exists        = .true.
        status             = 0
        message            = ''
        call fpath%kill
    end subroutine open_set

    !> mass, effective sample size and hard population of one full column
    subroutine column_populations( column, flags, mass, ess, pop )
        real(real32),  intent(in)  :: column(:)
        integer(int8), intent(in)  :: flags(:)
        real(real64),  intent(out) :: mass, ess
        integer,       intent(out) :: pop
        real(real64) :: sumsq
        mass  = sum(real(column, real64))
        sumsq = sum(real(column, real64)**2)
        ess   = 0._real64
        if( sumsq > 0._real64 ) ess = mass**2 / sumsq
        pop   = count(flags == 1_int8)
    end subroutine column_populations

    !> applied mass, effective sample size and contributor count of one column over a row subset
    subroutine row_populations( column, threshold, rows, mass, ess, n )
        real(real32), intent(in)  :: column(:)
        real,         intent(in)  :: threshold
        integer,      intent(in)  :: rows(:)
        real(real64), intent(out) :: mass, ess
        integer,      intent(out) :: n
        real(real64) :: sumsq, w
        integer      :: i
        mass = 0._real64; sumsq = 0._real64; n = 0
        do i = 1, size(rows)
            if( rows(i) < 1 .or. rows(i) > size(column) ) THROW_HARD('row outside the state weight set')
            if( column(rows(i)) <= threshold ) cycle
            w     = real(column(rows(i)), real64)
            mass  = mass  + w
            sumsq = sumsq + w * w
            n     = n + 1
        enddo
        ess = 0._real64
        if( sumsq > 0._real64 ) ess = mass**2 / sumsq
    end subroutine row_populations

    pure logical function close_enough( a, b )
        real(real64), intent(in) :: a, b
        close_enough = abs(a - b) <= POPULATION_TOL * max(1._real64, abs(a), abs(b))
    end function close_enough

end module simple_state_weight_set

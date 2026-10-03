!@descr: manifest of a gridding trailing-reconstruction chain: provenance, generation, component sizes and represented population
! The four chain components (even/odd Fourier sums and densities of one state) are one artifact
! set; its manifest is deleted first and written last. It records the crop box and sampling,
! the project row count, the state layout, a generation counter, the byte size of each
! component and M, the population the chain represents (population rule:
! population_blend_weights in simple_oris).
! A manifest of any other version is unreadable, so the chain is discarded and re-seeded.
! Validation policy stays with the owner (volassemble).
module simple_trail_chain_manifest
use simple_core_module_api
implicit none

public :: trail_chain_manifest
public :: TRAIL_MANIFEST_OK, TRAIL_MANIFEST_MISSING, TRAIL_MANIFEST_UNREADABLE
private
#include "simple_local_flags.inc"

integer, parameter :: TRAIL_MANIFEST_OK         = 0
integer, parameter :: TRAIL_MANIFEST_MISSING    = 1
integer, parameter :: TRAIL_MANIFEST_UNREADABLE = 2
integer, parameter :: MANIFEST_VERSION          = 2   ! 2: records the represented population

type :: trail_chain_manifest
    private
    integer         :: box      = 0
    real            :: smpd     = 0.
    integer         :: nptcls   = 0
    integer         :: nstates  = 0
    integer         :: state    = 0
    integer         :: gen      = 0
    integer(kind=8) :: sizes(4) = 0_8
    real            :: mrep     = 0.
    logical         :: exists   = .false.
  contains
    procedure :: new
    procedure :: write
    procedure :: read
    procedure :: get_box
    procedure :: get_smpd
    procedure :: get_nptcls
    procedure :: get_nstates
    procedure :: get_state
    procedure :: get_gen
    procedure :: get_size
    procedure :: get_mrep
    procedure :: kill
end type trail_chain_manifest

contains

    subroutine new( self, box, smpd, nptcls, nstates, state, gen, sizes, mrep )
        class(trail_chain_manifest), intent(inout) :: self
        integer,                     intent(in)    :: box, nptcls, nstates, state, gen
        real,                        intent(in)    :: smpd, mrep
        integer(kind=8),             intent(in)    :: sizes(4)
        call self%kill
        if( mrep < 0. ) THROW_HARD('negative represented population; trail_chain_manifest%new')
        self%box     = box
        self%smpd    = smpd
        self%nptcls  = nptcls
        self%nstates = nstates
        self%state   = state
        self%gen     = gen
        self%sizes   = sizes
        self%mrep    = mrep
        self%exists  = .true.
    end subroutine new

    !> write the manifest, status /= 0 on failure (the caller warns and re-seeds next time)
    subroutine write( self, fname, status )
        class(trail_chain_manifest), intent(in)  :: self
        class(string),               intent(in)  :: fname
        integer,                     intent(out) :: status
        integer :: funit
        if( .not. self%exists ) THROW_HARD('empty manifest; trail_chain_manifest%write')
        call fopen(funit, file=fname, status='REPLACE', action='WRITE', iostat=status)
        if( status /= 0 ) return
        write(funit,*,iostat=status) MANIFEST_VERSION, self%box, self%smpd, self%nptcls, self%nstates, &
            &self%state, self%gen, self%sizes, self%mrep
        call fclose(funit)
    end subroutine write

    !> read a manifest: TRAIL_MANIFEST_OK, _MISSING or _UNREADABLE (corrupt, or a version
    !! other than MANIFEST_VERSION); empty unless OK
    subroutine read( self, fname, status )
        class(trail_chain_manifest), intent(inout) :: self
        class(string),               intent(in)    :: fname
        integer,                     intent(out)   :: status
        character(len=1024) :: line
        integer :: funit, io_stat, version
        call self%kill
        status = TRAIL_MANIFEST_MISSING
        if( .not. file_exists(fname) ) return
        status = TRAIL_MANIFEST_UNREADABLE
        call fopen(funit, file=fname, status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 ) return
        read(funit,'(A)',iostat=io_stat) line
        call fclose(funit)
        if( io_stat /= 0 ) return
        read(line,*,iostat=io_stat) version, self%box, self%smpd, self%nptcls, self%nstates, &
            &self%state, self%gen, self%sizes, self%mrep
        if( io_stat /= 0 .or. version /= MANIFEST_VERSION )then
            call self%kill
            return
        endif
        if( self%mrep < 0. )then
            call self%kill
            return
        endif
        self%exists = .true.
        status      = TRAIL_MANIFEST_OK
    end subroutine read

    integer function get_box( self )
        class(trail_chain_manifest), intent(in) :: self
        get_box = self%box
    end function get_box

    real function get_smpd( self )
        class(trail_chain_manifest), intent(in) :: self
        get_smpd = self%smpd
    end function get_smpd

    integer function get_nptcls( self )
        class(trail_chain_manifest), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    integer function get_nstates( self )
        class(trail_chain_manifest), intent(in) :: self
        get_nstates = self%nstates
    end function get_nstates

    integer function get_state( self )
        class(trail_chain_manifest), intent(in) :: self
        get_state = self%state
    end function get_state

    integer function get_gen( self )
        class(trail_chain_manifest), intent(in) :: self
        get_gen = self%gen
    end function get_gen

    integer(kind=8) function get_size( self, i )
        class(trail_chain_manifest), intent(in) :: self
        integer,                     intent(in) :: i
        if( i < 1 .or. i > 4 ) THROW_HARD('chain component index out of range; get_size')
        get_size = self%sizes(i)
    end function get_size

    real function get_mrep( self )
        class(trail_chain_manifest), intent(in) :: self
        get_mrep = self%mrep
    end function get_mrep

    subroutine kill( self )
        class(trail_chain_manifest), intent(inout) :: self
        self%box     = 0
        self%smpd    = 0.
        self%nptcls  = 0
        self%nstates = 0
        self%state   = 0
        self%gen     = 0
        self%sizes   = 0_8
        self%mrep    = 0.
        self%exists  = .false.
    end subroutine kill

end module simple_trail_chain_manifest

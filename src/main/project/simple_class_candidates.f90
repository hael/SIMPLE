!@descr: class_candidates: the class student's candidate classes per particle, a side file registered in the project's out segment (imgkind class_candidates) and tied to the particle rows by the layout digest
! predict_classes publishes it, or its master from the parts' range files; the sparse 2D table (refine=prob_neigh)
! opens it and scores a seeded particle's candidates alone. A file whose digest, row count or class count differs from the consumer's is refused, and the
! consumer runs as without it (doc/algorithms/class_student.md).
module simple_class_candidates
use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
use simple_error,       only: simple_exception
use simple_oris,        only: oris
use simple_ptcl_layout, only: ptcl_layout_digest
use simple_sp_project,  only: sp_project
use simple_string,      only: string
use simple_syslib,      only: file_exists, simple_abspath
implicit none
private
#include "simple_local_flags.inc"

public :: class_candidates, run_all_class_candidates_tests
public :: CANDS_OK, CANDS_NONE, CANDS_UNREADABLE, CANDS_STALE

character(len=*), parameter :: IMGKIND = 'class_candidates'
character(len=8), parameter :: MAGIC   = 'CLSCAND1'
character(len=8), parameter :: MAGIC_PART = 'CLSCANDP' !< a distributed part's rows fromp..top, merged by publish_parts
integer,          parameter :: VERSION = 1
integer, parameter :: CANDS_OK         = 0 !< opened
integer, parameter :: CANDS_NONE       = 1 !< the project registers no artifact
integer, parameter :: CANDS_UNREADABLE = 2 !< registered, but missing or not this format
integer, parameter :: CANDS_STALE      = 3 !< another particle layout or another class count

type :: class_candidates
    private
    integer        :: nrows = 0, ncls = 0, kmax = 0
    integer(int64) :: digest = 0_int64
    integer(int8),  allocatable :: route(:)   !< (nrows) 1 seeded, 2 dense, 0 not predicted
    integer(int8),  allocatable :: ncand(:)   !< (nrows)
    integer(int16), allocatable :: cands(:,:) !< (kmax, nrows)
    logical :: exists = .false.
  contains
    procedure :: publish
    procedure :: write_part
    procedure :: publish_parts
    procedure :: new
    procedure :: is_seeded
    procedure :: get_cands
    procedure :: get_nseeded
    procedure :: kill
    procedure, private :: fill_rows
    procedure, private :: write_and_register
end type class_candidates

interface
    module subroutine run_all_class_candidates_tests()
    end subroutine run_all_class_candidates_tests
end interface

contains

    !> The candidates of the predicted rows pinds of particles (route 1 seeded, 2 dense; the other rows 0) written
    !! as fname under the layout digest of the project's particles and registered in its out segment, which is
    !! written to projfile
    subroutine publish( self, project, particles, projfile, fname, pinds, routes, ncands, cands, ncls )
        class(class_candidates), intent(inout) :: self
        class(sp_project),       intent(inout) :: project
        class(oris),             intent(in)    :: particles
        class(string),           intent(in)    :: projfile, fname
        integer,                 intent(in)    :: pinds(:), routes(:), ncands(:), cands(:,:), ncls
        call self%fill_rows(particles%get_noris(consider_state=.false.), 1, particles%get_noris(consider_state=.false.), &
            &pinds, routes, ncands, cands, ncls)
        call self%write_and_register(project, particles, projfile, fname)
    end subroutine publish

    !> A distributed part's candidates: rows fromp..top of nrows (the predicted pinds among them, the others 0) written
    !! as fname for the master's publish_parts; the project is not touched
    subroutine write_part( self, nrows, fromp, top, fname, pinds, routes, ncands, cands, ncls )
        class(class_candidates), intent(inout) :: self
        integer,                 intent(in)    :: nrows, fromp, top
        class(string),           intent(in)    :: fname
        integer,                 intent(in)    :: pinds(:), routes(:), ncands(:), cands(:,:), ncls
        integer :: u, ios
        if( fromp < 1 .or. top > nrows .or. fromp > top ) THROW_HARD('a part''s range outside the rows; class_candidates%write_part')
        call self%fill_rows(nrows, fromp, top, pinds, routes, ncands, cands, ncls)
        open(newunit=u, file=fname%to_char(), access='stream', form='unformatted', status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('could not open the candidates part file; class_candidates%write_part')
        write(u) MAGIC_PART, int(VERSION,int32), int(nrows,int32), int(ncls,int32), int(self%kmax,int32), &
            &int(fromp,int32), int(top,int32)
        write(u) self%route(fromp:top), self%ncand(fromp:top), self%cands(:,fromp:top)
        close(u)
        call self%kill
    end subroutine write_part

    !> The parts' files (write_part) read in order into one artifact over the particles' rows, which they must cover
    !! without overlap, then written and registered as publish does
    subroutine publish_parts( self, project, particles, projfile, fname, part_fnames, ncls )
        class(class_candidates), intent(inout) :: self
        class(sp_project),       intent(inout) :: project
        class(oris),             intent(in)    :: particles
        class(string),           intent(in)    :: projfile, fname, part_fnames(:)
        integer,                 intent(in)    :: ncls
        logical, allocatable :: covered(:)
        character(len=8) :: magic_read
        integer(int32)   :: version_read, nrows_read, ncls_read, kmax_read, fromp, top
        integer :: ipart, u, ios
        call self%kill
        self%nrows = particles%get_noris(consider_state=.false.)
        self%ncls  = ncls
        allocate(covered(self%nrows), source=.false.)
        do ipart = 1,size(part_fnames)
            open(newunit=u, file=part_fnames(ipart)%to_char(), access='stream', form='unformatted', status='old', &
                &action='read', iostat=ios)
            if( ios /= 0 ) THROW_HARD('a candidates part file is missing: '//part_fnames(ipart)%to_char())
            read(u, iostat=ios) magic_read, version_read, nrows_read, ncls_read, kmax_read, fromp, top
            if( ios /= 0 .or. magic_read /= MAGIC_PART .or. version_read /= VERSION ) &
                &THROW_HARD('not a candidates part file: '//part_fnames(ipart)%to_char())
            if( nrows_read /= self%nrows .or. ncls_read /= ncls ) &
                &THROW_HARD('a candidates part of other rows or classes: '//part_fnames(ipart)%to_char())
            if( ipart == 1 )then
                self%kmax = kmax_read
                allocate(self%route(self%nrows), self%ncand(self%nrows), source=0_int8)
                allocate(self%cands(self%kmax,self%nrows), source=0_int16)
            else if( kmax_read /= self%kmax )then
                THROW_HARD('candidates parts of different caps: '//part_fnames(ipart)%to_char())
            endif
            if( fromp < 1 .or. top > self%nrows .or. fromp > top ) &
                &THROW_HARD('a candidates part''s range outside the rows: '//part_fnames(ipart)%to_char())
            if( any(covered(fromp:top)) ) THROW_HARD('candidates parts overlap: '//part_fnames(ipart)%to_char())
            covered(fromp:top) = .true.
            read(u, iostat=ios) self%route(fromp:top), self%ncand(fromp:top), self%cands(:,fromp:top)
            close(u)
            if( ios /= 0 ) THROW_HARD('a truncated candidates part file: '//part_fnames(ipart)%to_char())
        enddo
        if( .not. all(covered) ) THROW_HARD('the candidates parts leave rows uncovered; class_candidates%publish_parts')
        call self%write_and_register(project, particles, projfile, fname)
    end subroutine publish_parts

    !> the rows' arrays over nrows from the predicted rows pinds, which lie in fromp..top
    subroutine fill_rows( self, nrows, fromp, top, pinds, routes, ncands, cands, ncls )
        class(class_candidates), intent(inout) :: self
        integer,                 intent(in)    :: nrows, fromp, top
        integer,                 intent(in)    :: pinds(:), routes(:), ncands(:), cands(:,:), ncls
        integer :: i, iptcl
        call self%kill
        self%nrows = nrows
        self%ncls  = ncls
        self%kmax  = size(cands,1)
        if( ncls < 1 .or. ncls > huge(1_int16) ) THROW_HARD('class count out of range; class_candidates')
        if( self%kmax < 1 .or. self%kmax > huge(1_int8) ) THROW_HARD('candidate cap out of range; class_candidates')
        allocate(self%route(self%nrows), self%ncand(self%nrows), source=0_int8)
        allocate(self%cands(self%kmax,self%nrows), source=0_int16)
        do i = 1,size(pinds)
            iptcl = pinds(i)
            if( iptcl < fromp .or. iptcl > top ) THROW_HARD('a row outside the range; class_candidates')
            if( ncands(i) < 0 .or. ncands(i) > self%kmax ) THROW_HARD('candidate count out of range; class_candidates')
            if( any(cands(1:ncands(i),i) < 1) .or. any(cands(1:ncands(i),i) > ncls) ) &
                &THROW_HARD('a candidate outside the classes; class_candidates')
            self%route(iptcl)              = int(routes(i), int8)
            self%ncand(iptcl)              = int(ncands(i), int8)
            self%cands(1:ncands(i),iptcl) = int(cands(1:ncands(i),i), int16)
        enddo
    end subroutine fill_rows

    !> the filled rows written as fname under the layout digest of the project's particles and registered in its out
    !! segment, which is written to projfile
    subroutine write_and_register( self, project, particles, projfile, fname )
        class(class_candidates), intent(inout) :: self
        class(sp_project),       intent(inout) :: project
        class(oris),             intent(in)    :: particles
        class(string),           intent(in)    :: projfile, fname
        type(string) :: path
        integer :: u, ios, ind
        self%digest = ptcl_layout_digest(project, particles)
        if( self%digest == 0_int64 ) THROW_HARD('the particle layout digest is undefined for this project; class_candidates')
        open(newunit=u, file=fname%to_char(), access='stream', form='unformatted', status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('could not open the candidates file; class_candidates')
        write(u) MAGIC, int(VERSION,int32), self%digest, int(self%nrows,int32), int(self%ncls,int32), int(self%kmax,int32)
        write(u) self%route, self%ncand, self%cands
        close(u)
        self%exists = .true.
        path = simple_abspath(fname)
        call project%add_entry2os_out(IMGKIND, ind)
        call project%os_out%set(ind, 'imgkind', IMGKIND)
        call project%os_out%set(ind, IMGKIND,   path%to_char())
        call project%write_segment_inside('out', projfile)
        call path%kill
    end subroutine write_and_register

    !> The project's registered candidates for its particles and a table of ncls classes: status CANDS_OK when the
    !! file is this format, its digest the particles' layout, its rows and classes the consumer's; otherwise the
    !! status says why (message), and the consumer runs without it
    subroutine new( self, project, particles, ncls, status, message )
        class(class_candidates),       intent(inout) :: self
        class(sp_project),             intent(in)    :: project
        class(oris),                   intent(in)    :: particles
        integer,                       intent(in)    :: ncls
        integer,                       intent(out)   :: status
        character(len=:), allocatable, intent(out)   :: message
        type(string)     :: path, kind
        character(len=8) :: magic_read
        integer(int32)   :: version_read, nrows_read, ncls_read, kmax_read
        integer(int64)   :: digest_read
        integer :: i, u, ios
        logical :: found
        call self%kill
        status  = CANDS_NONE
        message = 'no class candidates registered'
        found   = .false.
        do i = 1,project%os_out%get_noris()
            if( .not. project%os_out%isthere(i, 'imgkind') ) cycle
            call project%os_out%getter(i, 'imgkind', kind)
            if( kind%to_char() /= IMGKIND ) cycle
            if( .not. project%os_out%isthere(i, IMGKIND) ) cycle
            call project%os_out%getter(i, IMGKIND, path)
            found = .true.
        enddo
        if( .not. found ) return
        status  = CANDS_UNREADABLE
        message = 'the class candidates file is missing or unreadable: '//path%to_char()
        if( .not. file_exists(path) ) return
        open(newunit=u, file=path%to_char(), access='stream', form='unformatted', status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(u, iostat=ios) magic_read, version_read, digest_read, nrows_read, ncls_read, kmax_read
        if( ios /= 0 .or. magic_read /= MAGIC .or. version_read /= VERSION .or. kmax_read < 1 .or. nrows_read < 1 )then
            close(u)
            return
        endif
        status  = CANDS_STALE
        message = 'the class candidates belong to another particle layout or class count: '//path%to_char()
        if( digest_read /= ptcl_layout_digest(project, particles) .or. &
            &nrows_read /= particles%get_noris(consider_state=.false.) .or. ncls_read /= ncls )then
            close(u)
            return
        endif
        self%nrows  = nrows_read
        self%ncls   = ncls_read
        self%kmax   = kmax_read
        self%digest = digest_read
        allocate(self%route(self%nrows), self%ncand(self%nrows), self%cands(self%kmax,self%nrows))
        read(u, iostat=ios) self%route, self%ncand, self%cands
        close(u)
        if( ios /= 0 )then
            call self%kill
            status  = CANDS_UNREADABLE
            message = 'the class candidates file is truncated: '//path%to_char()
            return
        endif
        self%exists = .true.
        status  = CANDS_OK
        message = ''
    end subroutine new

    !> whether row iptcl is seeded: predicted with its candidates within the cap and confident
    pure logical function is_seeded( self, iptcl )
        class(class_candidates), intent(in) :: self
        integer,                 intent(in) :: iptcl
        is_seeded = .false.
        if( .not. self%exists ) return
        if( iptcl < 1 .or. iptcl > self%nrows ) return
        is_seeded = self%route(iptcl) == 1_int8 .and. self%ncand(iptcl) > 0_int8
    end function is_seeded

    !> row iptcl's candidates, the first n entries of cands (in descending mass)
    subroutine get_cands( self, iptcl, n, cands )
        class(class_candidates), intent(in)    :: self
        integer,                 intent(in)    :: iptcl
        integer,                 intent(out)   :: n
        integer,                 intent(inout) :: cands(:)
        n = 0
        if( .not. self%exists ) return
        if( iptcl < 1 .or. iptcl > self%nrows ) return
        n = min(int(self%ncand(iptcl)), size(cands))
        cands(1:n) = int(self%cands(1:n,iptcl))
    end subroutine get_cands

    !> how many of the rows pinds are seeded
    integer function get_nseeded( self, pinds )
        class(class_candidates), intent(in) :: self
        integer,                 intent(in) :: pinds(:)
        integer :: i
        get_nseeded = 0
        do i = 1,size(pinds)
            if( self%is_seeded(pinds(i)) ) get_nseeded = get_nseeded + 1
        enddo
    end function get_nseeded

    subroutine kill( self )
        class(class_candidates), intent(inout) :: self
        if( allocated(self%route) ) deallocate(self%route)
        if( allocated(self%ncand) ) deallocate(self%ncand)
        if( allocated(self%cands) ) deallocate(self%cands)
        self%nrows  = 0
        self%ncls   = 0
        self%kmax   = 0
        self%digest = 0_int64
        self%exists = .false.
    end subroutine kill

end module simple_class_candidates

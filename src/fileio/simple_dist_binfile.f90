module simple_dist_binfile
include 'simple_lib.f08'
implicit none

public :: dist_binfile
private
#include "simple_local_flags.inc"

type dist_binfile
    private
    character(len=:), allocatable :: fname
    integer :: file_header(3) = 0
    integer :: fromp          = 0
    integer :: top            = 0
    integer :: nrefs          = 0
    integer :: headsz         = 0
    integer :: datasz         = 0
    logical :: exists         = .false.
contains
    ! constructor
    procedure          :: new
    procedure          :: new_from_file
    ! I/O
    procedure          :: read
    procedure          :: read_to_glob
    procedure          :: read_from_glob
    procedure          :: write
    procedure          :: write_info
    procedure, private :: create_empty
    procedure, private :: open_and_check_header
    procedure, private :: open_only
    procedure, private :: read_header
    ! destructor
    procedure          :: kill
end type dist_binfile

contains

    subroutine new( self, fname, fromp, top, nrefs )
        class(dist_binfile), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        integer,             intent(in)    :: fromp, top, nrefs
        real(sp) :: r
        call self%kill
        self%fname          = trim(fname)
        self%fromp          = fromp
        self%top            = top
        self%nrefs          = nrefs
        self%file_header(1) = fromp
        self%file_header(2) = top
        self%file_header(3) = nrefs
        self%headsz         = sizeof(self%file_header)
        self%datasz         = sizeof(r) * self%nrefs
        self%exists         = .true.
    end subroutine new

    subroutine new_from_file( self, fname )
        class(dist_binfile), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        real(sp) :: r
        call self%kill
        if (.not. file_exists(fname)) then
            THROW_HARD('dist_binfile: new_from_file; file ' // trim(fname) // ' does not exist')
        end if
        self%fname          = trim(fname)
        call self%read_header
        self%file_header(1) = self%fromp
        self%file_header(2) = self%top
        self%file_header(3) = self%nrefs
        self%headsz         = sizeof(self%file_header)
        self%datasz         = sizeof(r) * self%nrefs
        self%exists         = .true.
    end subroutine new_from_file

    ! I/O

    ! read in all corrs value from file, corrs shape matches self%fromp, self%to
    subroutine read( self, corrs )
        class(dist_binfile), intent(inout) :: self
        real,                intent(inout) :: corrs(self%nrefs,self%fromp:self%top,2)
        call self%read_to_glob(self%fromp, self%top, corrs)
    end subroutine read

    ! read in all corrs value from file, corrs shape is smaller than self%fromp, self%to
    subroutine read_from_glob( self, fromp, top, corrs )
        class(dist_binfile), intent(inout) :: self
        integer,             intent(in)    :: fromp, top
        real,                intent(inout) :: corrs(self%nrefs,fromp:top,2)
        integer :: funit
        logical :: success
        integer :: iptcl, addr
        success = self%open_only( funit, .true. )
        if( .not. success ) return
        ! read corr
        addr = self%headsz + 1 + (fromp - self%fromp) * self%datasz
        do iptcl = fromp, top
            read(unit=funit,pos=addr) corrs(:,iptcl,1)
            addr = addr + self%datasz
        end do
        addr = self%headsz + 1 + (fromp - self%fromp) * self%datasz + (self%top - self%fromp + 1) * self%datasz
        ! read loc
        do iptcl = fromp, top
            read(unit=funit,pos=addr) corrs(:,iptcl,2)
            addr = addr + self%datasz
        end do
        call fclose(funit)
    end subroutine read_from_glob

    ! read in all corrs value from file, corrs shape is larger than self%fromp, self%to
    subroutine read_to_glob( self, fromp, top, corrs )
        class(dist_binfile), intent(inout) :: self
        integer,             intent(in)    :: fromp, top
        real,                intent(inout) :: corrs(self%nrefs,fromp:top,2)
        integer :: funit
        logical :: success
        integer :: iptcl, addr
        success = self%open_only( funit, .true. )
        if( .not. success ) return
        ! read corr
        addr = self%headsz + 1
        do iptcl = self%fromp, self%top
            read(unit=funit,pos=addr) corrs(:,iptcl,1)
            addr = addr + self%datasz
        end do
        ! read loc
        do iptcl = self%fromp, self%top
            read(unit=funit,pos=addr) corrs(:,iptcl,2)
            addr = addr + self%datasz
        end do
        call fclose(funit)
    end subroutine read_to_glob

    subroutine write( self, corrs )
        class(dist_binfile), intent(inout) :: self
        real,                intent(in)    :: corrs(self%nrefs,self%fromp:self%top,2)
        integer :: funit
        logical :: success
        integer :: addr, iptcl
        if( .not. file_exists(self%fname) )then
            call self%create_empty( funit )
        else
            success = self%open_and_check_header( funit, .false. )
            if( .not. success ) return
        end if
        ! write corr
        addr = self%headsz + 1
        do iptcl = self%fromp,self%top
            write(funit,pos=addr) corrs(:,iptcl,1)
            addr = addr + self%datasz
        end do
        ! write loc
        do iptcl = self%fromp,self%top
            write(funit,pos=addr) corrs(:,iptcl,2)
            addr = addr + self%datasz
        end do
        call fclose(funit)
    end subroutine write

    subroutine write_info( self )
        class(dist_binfile), intent(in) :: self
        write(logfhandle,*) 'fromp:  ',self%fromp
        write(logfhandle,*) 'top:    ',self%top
        write(logfhandle,*) 'nrefs:  ',self%nrefs
        write(logfhandle,*) 'datasz: ',self%datasz
    end subroutine write_info

    function open_and_check_header( self, funit, readonly ) result ( success )
        class(dist_binfile), intent(inout) :: self
        integer,             intent(out)   :: funit
        logical,             intent(in)    :: readonly
        logical :: success
        integer :: io_stat
        integer :: fromp_here, top_here, nrefs
        if( .not. file_exists(trim(self%fname)) )then
            success = .false.
            return
        end if
        if( readonly )then
            call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        else
            call fopen(funit,trim(self%fname),access='STREAM',status='OLD', iostat=io_stat)
        end if
        call fileiochk('dist_binfile; open_and_check_header; file: '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        fromp_here = self%file_header(1)
        top_here   = self%file_header(2)
        nrefs      = self%file_header(3)
        if ((fromp_here.ne.self%fromp) .or. (top_here.ne.self%top) .or. (nrefs.ne.self%nrefs)) then
            THROW_WARN( 'dimensions in corr file do not match')
            write (*,*) 'self%fromp: ', self%fromp, ' ; in corr file: ', fromp_here
            write (*,*) 'self%top: ',   self%top,   ' ; in corr file: ', top_here
            write (*,*) 'self%nrefs: ', self%nrefs, ' ; in corr file: ', nrefs
            THROW_HARD( 'exiting')
            call fclose(funit)
            success = .false.
        else
            success = .true.
        end if
    end function open_and_check_header

    function open_only( self, funit, readonly ) result ( success )
        class(dist_binfile), intent(inout) :: self
        integer,             intent(out)   :: funit
        logical,             intent(in)    :: readonly
        logical :: success
        integer :: io_stat
        success = .true.
        if( .not. file_exists(trim(self%fname)) )then
            success = .false.
            return
        end if
        if( readonly )then
            call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        else
            call fopen(funit,trim(self%fname),access='STREAM',status='OLD', iostat=io_stat)
        end if
        call fileiochk('dist_binfile; open_and_check_header; file: '//trim(self%fname), io_stat)
    end function open_only

    subroutine read_header( self )
        class(dist_binfile), intent(inout) :: self
        integer :: fromp_here, top_here, nrefs_here
        integer :: funit, io_stat
        integer :: file_header(3)
        call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('dist_binfile; read_header; file: '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) file_header
        fromp_here = file_header(1)
        top_here   = file_header(2)
        nrefs_here = file_header(3)
        if( top_here < fromp_here )then
            THROW_WARN('dist_binfile; read_header; header dimensions not making sense')
            write (*,*) 'fromp: ', fromp_here, ' ; top: ', top_here
            THROW_HARD( 'exiting')
        end if
        call fclose(funit)
        self%fromp = fromp_here
        self%top   = top_here
        self%nrefs = nrefs_here
    end subroutine read_header

    subroutine create_empty( self, funit )
        class(dist_binfile), intent(in)  :: self
        integer,             intent(out) :: funit
        integer  :: io_stat
        real(sp) :: corr_empty(self%nrefs, self%fromp:self%top,2)
        corr_empty = 0.
        call fopen(funit,trim(self%fname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) corr_empty
    end subroutine create_empty

    ! destructor

    subroutine kill( self )
        class(dist_binfile), intent(inout) :: self
        self%file_header  = 0
        self%headsz       = 0
        self%datasz       = 0
        self%fromp        = 0
        self%top          = 0
        self%nrefs        = 0
        self%exists       = .false.
    end subroutine kill

end module simple_dist_binfile

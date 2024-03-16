module simple_dist_binfile
include 'simple_lib.f08'
implicit none

public :: dist_binfile
private
#include "simple_local_flags.inc"

type dist_binfile
    private
    character(len=:), allocatable :: fname
    integer :: file_header(2) = 0
    integer :: nrefs          = 0
    integer :: nptcls         = 0
    integer :: headsz         = 0
    logical :: exists         = .false.
contains
    ! constructor
    procedure          :: new
    procedure          :: new_from_file
    ! I/O
    procedure          :: read_to_glob
    procedure          :: write
    procedure          :: write_info
    procedure, private :: open_and_check_header
    procedure, private :: open_only
    procedure, private :: read_header
    ! destructor
    procedure          :: kill
end type dist_binfile

contains

    subroutine new( self, fname, nrefs, nptcls )
        class(dist_binfile), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        integer,             intent(in)    :: nrefs, nptcls
        call self%kill
        self%fname          = trim(fname)
        self%nrefs          = nrefs
        self%nptcls         = nptcls
        self%file_header(1) = nrefs
        self%file_header(2) = nptcls
        self%headsz         = sizeof(self%file_header)
        self%exists         = .true.
    end subroutine new

    subroutine new_from_file( self, fname )
        class(dist_binfile), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        call self%kill
        if (.not. file_exists(fname)) then
            THROW_HARD('dist_binfile: new_from_file; file ' // trim(fname) // ' does not exist')
        end if
        self%fname  = trim(fname)
        call self%read_header
        self%headsz = sizeof(self%file_header)
        self%exists = .true.
    end subroutine new_from_file

    ! I/O

    ! read in all dists value from file, dists shape is larger than self%fromp, self%to
    subroutine read_to_glob( self, nptcls_glob, mat_glob )
        class(dist_binfile), intent(inout) :: self
        integer,             intent(in)    :: nptcls_glob
        type(ptcl_ref),      intent(inout) :: mat_glob(self%nrefs,nptcls_glob)
        type(ptcl_ref) :: mat_loc(self%nrefs,self%nptcls)
        integer :: funit, i_loc, addr, i_glob
        logical :: success
        success = self%open_only( funit, .true. )
        if( .not. success ) return
        addr = self%headsz + 1
        ! read partition information
        read(unit=funit,pos=addr) mat_loc
        call fclose(funit)
        !$omp parallel do collapse(2) default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( mat_loc(1,i_loc)%pind == mat_glob(1,i_glob)%pind ) mat_glob(:,i_glob) = mat_loc(:,i_loc)
            end do
        end do
        !$omp end parallel do
    end subroutine read_to_glob

    subroutine write( self, mat )
        class(dist_binfile), intent(inout) :: self
        type(ptcl_ref),      intent(in)    :: mat(self%nrefs,self%nptcls)
        integer :: funit, addr
        logical :: success
        call fopen(funit,trim(self%fname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) self%file_header
        addr = self%headsz + 1
        write(funit,pos=addr) mat
        call fclose(funit)
    end subroutine write

    subroutine write_info( self )
        class(dist_binfile), intent(in) :: self
        write(logfhandle,*) 'nrefs:  ',self%nrefs
        write(logfhandle,*) 'nptcls: ',self%nptcls
    end subroutine write_info

    function open_and_check_header( self, funit, readonly ) result ( success )
        class(dist_binfile), intent(inout) :: self
        integer,             intent(out)   :: funit
        logical,             intent(in)    :: readonly
        logical :: success
        integer :: io_stat, nrefs, nptcls
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
        nrefs  = self%file_header(1)
        nptcls = self%file_header(2)
        if( (nrefs.ne.self%nrefs) .or. (nptcls.ne.self%nptcls) )then
            THROW_WARN( 'dimensions in dist file do not match')
            write (*,*) 'self%nrefs: ',  self%nrefs,  ' ; in dist file: ', nrefs
            write (*,*) 'self%nptcls: ', self%nptcls, ' ; in dist file: ', nptcls  
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
        integer :: funit, io_stat
        call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('dist_binfile; read_header; file: '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        call fclose(funit)
        self%nrefs  = self%file_header(1)
        self%nptcls = self%file_header(2)
    end subroutine read_header

    ! destructor

    subroutine kill( self )
        class(dist_binfile), intent(inout) :: self
        self%file_header = 0
        self%headsz      = 0
        self%nrefs       = 0
        self%nptcls      = 0
        self%exists      = .false.
    end subroutine kill

end module simple_dist_binfile

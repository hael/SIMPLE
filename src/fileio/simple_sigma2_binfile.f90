module simple_sigma2_binfile
include 'simple_lib.f08'
implicit none

public :: sigma2_binfile
private
#include "simple_local_flags.inc"

type sigma2_binfile
    private
    character(len=:), allocatable :: fname
    integer                       :: file_header(4) = 0
    integer                       :: fromp          = 0
    integer                       :: top            = 0
    integer                       :: kfromto(2)     = 0
    integer                       :: headsz         = 0
    integer                       :: sigmassz       = 0
    logical                       :: exists         = .false.

contains
    ! constructor
    procedure                     :: new
    procedure                     :: new_from_file
    ! I/O
    procedure                     :: read
    procedure                     :: write
    procedure                     :: write_info
    procedure, private            :: create_empty
    procedure, private            :: open_and_check_header
    procedure, private            :: read_header
    ! getters / setters
    procedure                     :: kill
end type sigma2_binfile

contains

    subroutine new( self, fname, fromp, top, kfromto )
        class(sigma2_binfile),         intent(inout) :: self
        character(len=:), allocatable, intent(in)    :: fname
        integer,                       intent(in)    :: fromp, top, kfromto(2)
        real(sp) :: r
        call self%kill
        self%fname   = trim(fname)
        self%fromp   = fromp
        self%top     = top
        self%kfromto = kfromto
        self%file_header(1)   = fromp
        self%file_header(2)   = top
        self%file_header(3:4) = kfromto
        self%headsz           = sizeof(self%file_header)
        self%sigmassz         = sizeof(r)*(self%kfromto(2)-self%kfromto(1)+1)
        self%exists = .true.
    end subroutine new

    subroutine new_from_file( self, fname )
        class(sigma2_binfile),         intent(inout) :: self
        character(len=:), allocatable, intent(in)    :: fname
        real(sp) :: r
        call self%kill
        if (.not. file_exists(fname)) then
            THROW_HARD('sigma2_binfile: new_from_file; file ' // trim(fname) // ' does not exist')
        end if
        self%fname   = trim(fname)
        call self%read_header
        self%file_header(1)   = self%fromp
        self%file_header(2)   = self%top
        self%file_header(3:4) = self%kfromto(:)
        self%headsz           = sizeof(self%file_header)
        self%sigmassz         = sizeof(r)*(self%kfromto(2)-self%kfromto(1)+1)
        self%exists = .true.
    end subroutine new_from_file

    ! utils

    subroutine write_info( self )
        class(sigma2_binfile), intent(in) :: self
        write(logfhandle,*) 'fromto:   ',self%kfromto
        write(logfhandle,*) 'fromp:    ',self%fromp
        write(logfhandle,*) 'top:      ',self%top
        write(logfhandle,*) 'sigmassz: ',self%sigmassz
    end subroutine write_info

    ! I/O

    subroutine read( self, sigma2 )! read in all sigmas from file
        class(sigma2_binfile), intent(inout) :: self
        real, allocatable,     intent(out)   :: sigma2(:,:)
        integer :: funit
        logical :: success
        integer :: iptcl, addr
        real    :: sigma2_noise_n(self%kfromto(1):self%kfromto(2))
        allocate(sigma2(self%kfromto(1):self%kfromto(2),self%fromp:self%top),source=0.0)
        success = self%open_and_check_header( funit, .true. )
        if( .not. success ) return
        do iptcl = self%fromp, self%top
            addr = self%headsz + (iptcl - self%fromp) * self%sigmassz + 1
            read(unit=funit,pos=addr) sigma2_noise_n
            sigma2(:,iptcl) = sigma2_noise_n
        end do
        call fclose(funit)
    end subroutine read

    subroutine write( self, sigma2 )
        class(sigma2_binfile),  intent(inout) :: self
        real, allocatable,      intent(in)    :: sigma2(:,:)
        integer :: funit
        logical :: success
        integer :: addr, iptcl
        ! make sure the dimensions of inout matrix agree
        if( (lbound(sigma2,2).ne.self%fromp)     .or.(ubound(sigma2,2).ne.self%top).or.&
            (lbound(sigma2,1).ne.self%kfromto(1)).or.(ubound(sigma2,1).ne.self%kfromto(2))) then
            THROW_WARN('sigma2_binfile: read; dimensions of sigma2 dont agree')
            write (*,*) 'self%fromp: ',   self%fromp,   ' ; lbound(sigma2,2): ', lbound(sigma2,2)
            write (*,*) 'self%top: ',     self%top,     ' ; ubound(sigma2,2): ', ubound(sigma2,2)
            write (*,*) 'self%kfromto: ', self%kfromto, ' ; size(sigma2,1): ', size(sigma2,1)
            THROW_HARD( 'exiting')
        end if
        if( .not. file_exists(self%fname) )then
            call self%create_empty( funit )
        else
            success = self%open_and_check_header( funit, .false. )
            if( .not. success ) return
        end if
        do iptcl = self%fromp,self%top
            addr = self%headsz + (iptcl - self%fromp) * self%sigmassz + 1
            write(funit,pos=addr) sigma2(:,iptcl)
        end do
        call fclose(funit)
    end subroutine

    function open_and_check_header( self, funit, readonly ) result ( success )
        class(sigma2_binfile), intent(inout) :: self
        integer,               intent(out)   :: funit
        logical,               intent(in)    :: readonly
        logical :: success
        integer :: io_stat
        integer :: fromp_here, top_here, kfromto_here(2)
        if( .not. file_exists(trim(self%fname)) )then
            success = .false.
            return
        end if
        if( readonly )then
            call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        else
            call fopen(funit,trim(self%fname),access='STREAM',status='OLD', iostat=io_stat)
        end if
        call fileiochk('sigma2_binfile; open_and_check_header; file: '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) self%file_header
        fromp_here      = self%file_header(1)
        top_here        = self%file_header(2)
        kfromto_here    = self%file_header(3:4)
        if ((fromp_here     .ne.self%fromp     ) .or. (top_here.ne.self%top              )  .or. &
            (kfromto_here(1).ne.self%kfromto(1)) .or. (kfromto_here(2).ne.self%kfromto(2))) then
            THROW_WARN( 'dimensions in sigmas file do not match')
            write (*,*) 'self%fromp: ',   self%fromp,   ' ; in sigmas file: ', fromp_here
            write (*,*) 'self%top: ',     self%top,     ' ; in sigmas file: ', top_here
            write (*,*) 'self%kfromto: ', self%kfromto, ' ; in sigmas file: ', kfromto_here
            THROW_HARD( 'exiting')
            call fclose(funit)
            success = .false.
        else
            success = .true.
        end if
    end function open_and_check_header

    subroutine read_header( self )
        class(sigma2_binfile), intent(inout) :: self
        integer :: fromp_here, top_here, kfromto_here(2)
        integer :: funit, io_stat
        integer :: file_header(4)
        call fopen(funit,trim(self%fname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        call fileiochk('sigma2_binfile; read_header; file: '//trim(self%fname), io_stat)
        read(unit=funit,pos=1) file_header
        fromp_here      = file_header(1)
        top_here        = file_header(2)
        kfromto_here    = file_header(3:4)
        if(( kfromto_here(2) < kfromto_here(1) ).or.(top_here < fromp_here))then
            THROW_WARN('sigma2_binfile; read_header; header dimensions not making sense')
            write (*,*) 'fromp: ',      fromp_here,      ' ; top: ',        top_here
            write (*,*) 'kfromto(1): ', kfromto_here(1), ' ; kfromto(2): ', kfromto_here(2)
            THROW_HARD( 'exiting')
        end if
        call fclose(funit)
        self%fromp   = fromp_here
        self%top     = top_here
        self%kfromto = kfromto_here
    end subroutine read_header

    subroutine create_empty( self, funit )
        class(sigma2_binfile), intent(in)  :: self
        integer,               intent(out) :: funit
        integer  :: io_stat
        real(sp) :: sigmas_empty(self%kfromto(1):self%kfromto(2), self%fromp:self%top)
        sigmas_empty = 0.
        call fopen(funit,trim(self%fname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) self%file_header
        write(unit=funit,pos=self%headsz + 1) sigmas_empty
    end subroutine create_empty

    ! destructor

    subroutine kill( self )
        class(sigma2_binfile), intent(inout) :: self
        self%file_header  = 0
        self%headsz       = 0
        self%sigmassz     = 0
        self%kfromto      = 0
        self%fromp        = 0
        self%top          = 0
        self%exists       = .false.
    end subroutine kill


end module simple_sigma2_binfile

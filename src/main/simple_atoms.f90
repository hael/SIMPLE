! atomic structures and pdb parser
module simple_atoms
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error, only: allocchk
use simple_strings
use simple_fileio
implicit none

public :: atoms
private
#include "simple_local_flags.inc"

character(len=74) :: pdbfmt      = "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)" ! custom 3.3
character(len=74) :: pdbfmt_long = "(A5,I6,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)" ! custom 3.3
character(len=74) :: pdbfmt_read = "(A11,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)"   ! custom 3.3

!>  \brief type for dealing with atomic structures
type :: atoms
    private
    integer                       :: n = 0
    character(len=4), allocatable :: name(:)
    character(len=1), allocatable :: altloc(:)
    character(len=1), allocatable :: chain(:)
    character(len=1), allocatable :: icode(:)
    character(len=3), allocatable :: resname(:)
    character(len=3), allocatable :: element(:)
    real,             allocatable :: charge(:)
    real,             allocatable :: xyz(:,:)
    real,             allocatable :: mw(:)
    real,             allocatable :: occupancy(:)
    real,             allocatable :: beta(:)
    integer,          allocatable :: num(:)
    integer,          allocatable :: resnum(:)
    integer,          allocatable :: Z(:)
    logical,          allocatable :: het(:)
    logical                       :: exists = .false.
  contains
    ! CONSTRUCTORS
    procedure, private :: new_instance
    procedure, private :: new_from_pdb
    generic            :: new => new_from_pdb, new_instance
    generic            :: assignment(=) => copy
    procedure, private :: copy
    ! GETTERS/SETTERS
    procedure          :: does_exist
    procedure          :: get_n
    procedure          :: get_name
    procedure          :: get_coord
    procedure          :: get_num
    procedure          :: set_coord
    procedure          :: set_chain
    procedure          :: set_num
    procedure          :: set_beta
    procedure          :: set_resnum
    ! I/O
    procedure          :: print_atom
    procedure          :: writepdb
    ! CALCULATORS
    procedure, private :: Z_from_name
    procedure          :: get_geom_center
    ! MODIFIERS
    procedure          :: translate
    ! DESTRUCTOR
    procedure          :: kill
end type atoms

contains

    ! CONSTRUCTORS

    subroutine new_from_pdb( self, fname )
        class(atoms),     intent(inout) :: self
        character(len=*), intent(in)    :: fname
        character(len=STDLEN) :: line
        character(len=11)     :: elevenfirst
        character(len=6)     :: atom_field
        integer               :: i, l, nl, filnum, io_stat, n, num
        call self%kill
        nl = nlines(trim(fname))
        if(nl == 0 .or. .not.file_exists(fname))then
            THROW_HARD('I/O, file: '//trim(fname)//'; new_from_pdb')
        endif
        call fopen(filnum, status='OLD', action='READ', file=fname, iostat=io_stat)
        call fileiochk('new_from_pdb; simple_atoms opening '//trim(fname), io_stat)
        ! first pass
        n = 0
        do i = 1, nl
            read(filnum,'(A11)')elevenfirst
            if( .not.is_valid_entry(elevenfirst(1:6)) )cycle
            ! support for over 100000 entries
            call str2int(elevenfirst(7:11), io_stat, num )
            if( io_stat .ne. 0 )call str2int(elevenfirst(6:11), io_stat, num )
            if( io_stat .ne. 0 )call str2int(elevenfirst(5:11), io_stat, num )
            if( io_stat .ne. 0 )cycle
            n = n + 1
        enddo
        if( n == 0 )then
            THROW_HARD('no atoms found; new_from_pdb')
        endif
        ! second pass
        call self%new_instance(n)
        rewind(filnum)
        i = 0
        do l = 1, nl
            read(filnum,'(A)')line
            if( .not.is_valid_entry(line(1:6)) )cycle
            call str2int(line(7:11), io_stat, num )
            if( io_stat .ne. 0 )call str2int(line(6:11), io_stat, num )
            if( io_stat .ne. 0 )call str2int(line(5:11), io_stat, num )
            if( io_stat .ne. 0 )cycle
            i = i + 1
            read(line,pdbfmt_read, iostat=io_stat)elevenfirst, self%name(i), self%altloc(i),&
            &self%resname(i), self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
            &self%occupancy(i), self%beta(i)
            self%num(i) = num
            call fileiochk('new_from_pdb; simple_atoms error reading line '//trim(fname), io_stat)
            self%het(i) = atom_field == 'HETATM'
        enddo
        ! done
        call fclose(filnum, errmsg='new_from_pdb; simple_atoms closing '//trim(fname))
        contains

            elemental logical function is_valid_entry( str )
                character(len=6), intent(in) :: str
                is_valid_entry = .false.
                if( str(1:6).eq.'HETATM' ) is_valid_entry = .true.
                if( str(1:4).eq.'ATOM' )   is_valid_entry = .true.
            end function
    end subroutine new_from_pdb

    subroutine new_instance( self, n, dummy )
        class(atoms),      intent(inout) :: self
        integer,           intent(inout) :: n
        logical, optional, intent(in)    :: dummy
        integer :: i,alloc_stat
        logical :: ddummy
        ddummy = .false.
        if(present(dummy)) ddummy = dummy
        call self%kill
        allocate(self%name(n), self%chain(n), self%resname(n), self%xyz(n,3), self%mw(n),&
            self%occupancy(n), self%beta(n), self%num(n), self%Z(n), self%het(n), self%icode(n),&
            self%altloc(n), self%resnum(n), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('new_instance :: simple_atoms', alloc_stat)
        self%name(:)    = '    '
        self%resname(:) = '   '
        self%chain(:)   = ' '
        self%altloc(:)  = ' '
        self%icode(:)   = ' '
        self%mw        = 0.
        self%xyz       = 0.
        self%beta      = 0.
        self%occupancy = 0.
        self%num    = 0
        self%resnum = 0
        self%Z      = 0
        self%n      = n
        self%het    = .false.
        if(ddummy)then
            do i=1,self%n
                self%name(i)    = ' X  '
                self%resname(:) = ' X '
                self%chain(:)   = 'A'
                self%beta      = 1.
                self%occupancy = 1.
                self%num    = i
                self%resnum = 1
            enddo
        endif
        self%exists = .true.
    end subroutine new_instance

    subroutine copy( self, self_in )
        class(atoms), intent(inout) :: self
        class(atoms), intent(in)    :: self_in
        integer :: n
        if(.not.self_in%exists) THROW_HARD('Unallocated input instance; copy')
        n = self_in%n
        call self%new_instance(n)
        self%name(:)    = self_in%name
        self%resname(:) = self_in%resname
        self%chain(:)   = self_in%chain
        self%altloc(:)  = self_in%altloc
        self%icode(:)   = self_in%icode
        self%mw         = self_in%mw
        self%xyz        = self_in%xyz
        self%beta       = self_in%beta
        self%occupancy  = self_in%occupancy
        self%num        = self_in%num
        self%resnum     = self_in%resnum
        self%Z          = self_in%Z
        self%n          = n
        self%het        = self_in%het
    end subroutine copy

    ! GETTERS / SETTERS

    logical function does_exist( self )
        class(atoms), intent(inout) :: self
        does_exist = self%exists
    end function does_exist

    integer function get_n( self )
        class(atoms), intent(in) :: self
        get_n = self%n
    end function get_n

    function get_coord( self, i )result( xyz )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        real :: xyz(3)
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_coord')
        xyz = self%xyz(i,:)
    end function get_coord

    integer function get_num( self, i )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_coord')
        get_num = self%num(i)
    end function get_num

    character(len=4) function get_name( self, i )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        get_name = self%name(i)
    end function get_name

    subroutine set_coord( self, i, xyz )
        class(atoms), intent(inout) :: self
        integer,      intent(in)    :: i
        real,         intent(in)    :: xyz(3)
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_coord')
        self%xyz(i,:) = xyz
    end subroutine set_coord

    subroutine set_chain( self, i, chain )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=1), intent(in)    :: chain
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_coord')
        self%chain(i) = chain
    end subroutine set_chain

    subroutine set_num( self, i, num )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        integer,          intent(in)    :: num
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_num')
        self%num(i) = num
    end subroutine set_num

    subroutine set_beta( self, i, bfac )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        real,             intent(in)    :: bfac
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_beta')
        self%beta(i) = bfac
    end subroutine set_beta

    subroutine set_resnum( self, i, resnum )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        integer,          intent(in)    :: resnum
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_resnum')
        self%resnum(i) = resnum
    end subroutine set_resnum

    subroutine print_atom( self, i )
        class(atoms), intent(inout) :: self
        integer,      intent(in)    :: i
        if(self%het(i))then
            write(logfhandle,'(A6)',advance='no')'HETERO-ATOM '
        else
            write(logfhandle,'(A6)',advance='no')'ATOM        '
        endif
        write(logfhandle,'(I6,1X)',advance='no')self%num(i)
        write(logfhandle,'(A4,1X)',advance='no')self%name(i)
        write(logfhandle,'(A1,1X)',advance='no')self%altloc(i)
        write(logfhandle,'(A3,1X)',advance='no')self%resname(i)
        write(logfhandle,'(A1,1X)',advance='no')self%chain(i)
        write(logfhandle,'(I4,1X)',advance='no')self%resnum(i)
        write(logfhandle,'(A1,1X)',advance='no')self%icode(i)
        write(logfhandle,'(3F8.3,1X)',advance='no')self%xyz(i,:)
        write(logfhandle,'(2F6.2,1X)',advance='yes')self%occupancy(i), self%beta(i)
    end subroutine print_atom

    ! I/O

    subroutine writepdb( self, fbody )
        class(atoms),     intent(in) :: self
        character(len=*), intent(in) :: fbody
        character(len=STDLEN) :: fname
        character(len=76)     :: line
        integer               :: i, funit, io_stat
        logical               :: long
        fname = trim(adjustl(fbody)) // '.pdb'
        long  = self%n >= 99999
        if(.not.self%exists) THROW_HARD('Cannot write non existent atoms type; writePDB')
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
        call fileiochk('writepdb; simple_atoms opening '//trim(fname), io_stat)
        do i = 1, self%n
            write(funit,'(A76)')pdbstr(i)
        enddo
        call fclose(funit, errmsg='writepdb; simple_atoms closing '//trim(fname))
        contains

            character(len=76) function pdbstr( ind )
                integer,           intent(in)  :: ind
                character(len=6) :: atom_field
                if(self%het(ind))then
                    atom_field(1:6) = 'HETATM'
                    write(pdbstr,pdbfmt)atom_field,self%num(ind),self%name(ind),self%altloc(ind),&
                        self%resname(ind),self%chain(ind), self%resnum(ind), self%icode(ind), self%xyz(ind,:),&
                        self%occupancy(ind), self%beta(ind)
                else
                    if( long )then
                        atom_field(1:5) = 'ATOM '
                        write(pdbstr,pdbfmt_long)atom_field,self%num(ind),self%name(ind),self%altloc(ind),&
                            self%resname(ind),self%chain(ind), self%resnum(ind), self%icode(ind), self%xyz(ind,:),&
                            self%occupancy(ind), self%beta(ind)

                    else
                        atom_field(1:6) = 'ATOM  '
                        write(pdbstr,pdbfmt)atom_field,self%num(ind),self%name(ind),self%altloc(ind),&
                            self%resname(ind),self%chain(ind), self%resnum(ind), self%icode(ind), self%xyz(ind,:),&
                            self%occupancy(ind), self%beta(ind)
                    endif
                endif
            end function pdbstr
    end subroutine writepdb

    ! CALCULATORS

    integer function Z_from_name( self, name )
        class(atoms),     intent(inout) :: self
        character(len=4), intent(in)    :: name
        character(len=4) :: uppercase_name
        uppercase_name = upperCase(name)
        select case(trim(adjustl(uppercase_name)))
        ! organic
        case('H')
            Z_from_name = 1
        case('C')
            Z_from_name = 6
        case('N')
            Z_from_name = 7
        case('O')
            Z_from_name = 8
        case('P')
            Z_from_name = 15
        case('S')
            Z_from_name = 16
        ! metals
        case('FE')
            Z_from_name = 26
        case('PD')
            Z_from_name = 46
        case('PT')
            Z_from_name = 78
        end select
    end function Z_from_name

    function get_geom_center(self) result(center)
        class(atoms), intent(in) :: self
        real :: center(3)
        center(1) = sum(self%xyz(:,1)) / real(self%n)
        center(2) = sum(self%xyz(:,2)) / real(self%n)
        center(3) = sum(self%xyz(:,3)) / real(self%n)
    end function get_geom_center

    ! MODIFIERS

    subroutine translate( self, shift )
        class(atoms), intent(inout) :: self
        real,         intent(in)    :: shift(3)
        self%xyz(:,1) = self%xyz(:,1) + shift(1)
        self%xyz(:,2) = self%xyz(:,2) + shift(2)
        self%xyz(:,3) = self%xyz(:,3) + shift(3)
    end subroutine translate

    ! DESTRUCTOR
    subroutine kill( self )
        class(atoms), intent(inout) :: self
        if( allocated(self%name) )deallocate(self%name)
        if( allocated(self%altloc) )deallocate(self%altloc)
        if( allocated(self%chain) )deallocate(self%chain)
        if( allocated(self%icode) )deallocate(self%icode)
        if( allocated(self%resname) )deallocate(self%resname)
        if( allocated(self%xyz) )deallocate(self%xyz)
        if( allocated(self%mw) )deallocate(self%mw)
        if( allocated(self%occupancy) )deallocate(self%occupancy)
        if( allocated(self%beta) )deallocate(self%beta)
        if( allocated(self%num) )deallocate(self%num)
        if( allocated(self%resnum) )deallocate(self%resnum)
        if( allocated(self%Z) )deallocate(self%Z)
        if( allocated(self%het) )deallocate(self%het)
        self%n      = 0
        self%exists = .false.
    end subroutine kill

end module

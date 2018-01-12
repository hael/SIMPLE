! atomic structures and pdb parser
module simple_atoms
!$ use omp_lib
!$ use omp_lib_kinds
#include "simple_lib.f08"
implicit none

public :: atoms
private
#include "simple_local_flags.inc"

!character(len=74) :: pdbfmt = "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)" ! v3.3
character(len=74) :: pdbfmt = "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)" ! custom 3.3

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
    procedure          :: set_resnum
    ! I/O
    procedure          :: print_atom
    procedure          :: writepdb
    ! CALCULATORS
    procedure, private :: Z_from_name
    procedure          :: get_geom_center
    ! MODIFIERS
    ! DESTRUCTOR
    procedure          :: kill
end type atoms

contains

    ! CONSTRUCTORS

    subroutine new_from_pdb( self, fname )
        class(atoms),     intent(inout) :: self
        character(len=*), intent(in)    :: fname
        character(len=STDLEN) :: line
        character(len=6)      :: atom_field
        integer               :: i, l, nl, filnum, io_stat, n
        call self%kill
        nl = nlines(trim(fname))
        if(nl == 0 .or. .not.file_exists(fname))then
            print *, 'IO problem with file:', trim(fname)
            stop 'simple_atoms :: new_from_pdb'
        endif
        call fopen(filnum, status='OLD', action='READ', file=fname, iostat=io_stat)
        call fileio_errmsg('new_from_pdb; simple_atoms opening '//trim(fname), io_stat)
        ! first pass
        n = 0
        do i = 1, nl
            read(filnum,'(A6)')atom_field
            if(is_valid_entry(atom_field))n = n + 1
        enddo
        if( n == 0 )then
            print *,'No atoms found in:'
            stop 'simple_atoms :: new_from_pdb'
        endif
        ! instance
        call self%new_instance(n)
        ! second pass
        rewind(filnum)
        i = 0
        do l = 1, nl
            read(filnum,'(A)')line
            if(.not.is_valid_entry(line(1:6)))cycle
            i = i + 1
            read(line,pdbfmt, iostat=io_stat)atom_field, self%num(i), self%name(i), self%altloc(i),&
            &self%resname(i), self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
            &self%occupancy(i), self%beta(i)
            call fileio_errmsg('new_from_pdb; simple_atoms error reading line '//trim(fname), io_stat)
            self%het(i) = atom_field == 'HETATM'
        enddo
        ! done
        call fclose(filnum, errmsg='new_from_pdb; simple_atoms closing '//trim(fname))
        contains

            logical function is_valid_entry( str )
                character(len=6), intent(in) :: str
                select case(str)
                case( 'ATOM ', 'HETATM' )
                    is_valid_entry = .true.
                case DEFAULT
                    is_valid_entry = .false.
                end select
            end function
    end subroutine new_from_pdb

    subroutine new_instance( self, n )
        class(atoms), intent(inout) :: self
        integer,      intent(inout) :: n
        integer :: i, alloc_stat
        call self%kill
        allocate(self%name(n), self%chain(n), self%resname(n), self%xyz(n,3), self%mw(n),&
            self%occupancy(n), self%beta(n), self%num(n), self%Z(n), self%het(n), self%icode(n),&
            self%altloc(n), self%resnum(n), stat=alloc_stat)
        allocchk('new_instance :: simple_atoms')
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
        self%exists = .true.
    end subroutine new_instance

    subroutine copy( self, self_in )
        class(atoms), intent(inout) :: self
        class(atoms), intent(in)    :: self_in
        integer :: n
        if(.not.self_in%exists)stop 'Unallocated input instance; simple_atoms%copy'
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
        if(i.lt.1 .or. i.gt.self%n)stop 'index out of range; simple_atoms%get_coord'
        xyz = self%xyz(i,:)
    end function get_coord

    integer function get_num( self, i )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        if(i.lt.1 .or. i.gt.self%n)stop 'index out of range; simple_atoms%get_coord'
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
        if(i.lt.1 .or. i.gt.self%n)stop 'index out of range; simple_atoms%set_coord'
        self%xyz(i,:) = xyz
    end subroutine set_coord

    subroutine set_chain( self, i, chain )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=1), intent(in)    :: chain
        if(i.lt.1 .or. i.gt.self%n)stop 'index out of range; simple_atoms%set_coord'
        self%chain(i) = chain
    end subroutine set_chain

    subroutine set_num( self, i, num )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        integer,          intent(in)    :: num
        if(i.lt.1 .or. i.gt.self%n)stop 'index out of range; simple_atoms%set_num'
        self%num(i) = num
    end subroutine set_num

    subroutine set_resnum( self, i, resnum )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        integer,          intent(in)    :: resnum
        if(i.lt.1 .or. i.gt.self%n)stop 'index out of range; simple_atoms%set_resnum'
        self%resnum(i) = resnum
    end subroutine set_resnum

    subroutine print_atom( self, i )
        class(atoms), intent(inout) :: self
        integer,      intent(in)    :: i
        if(self%het(i))then
            write(*,'(A6)',advance='no')'HETERO-ATOM '
        else
            write(*,'(A6)',advance='no')'ATOM        '
        endif
        write(*,'(I6,1X)',advance='no')self%num(i)
        write(*,'(A4,1X)',advance='no')self%name(i)
        write(*,'(A1,1X)',advance='no')self%altloc(i)
        write(*,'(A3,1X)',advance='no')self%resname(i)
        write(*,'(A1,1X)',advance='no')self%chain(i)
        write(*,'(I4,1X)',advance='no')self%resnum(i)
        write(*,'(A1,1X)',advance='no')self%icode(i)
        write(*,'(3F8.3,1X)',advance='no')self%xyz(i,:)
        write(*,'(2F6.2,1X)',advance='yes')self%occupancy(i), self%beta(i)
    end subroutine print_atom

    ! I/O

    subroutine writepdb( self, fbody )
        class(atoms),     intent(inout) :: self
        character(len=*), intent(in)    :: fbody
        character(len=STDLEN) :: fname
        character(len=76)     :: line
        integer               :: i, funit, io_stat
        fname = trim(adjustl(fbody)) // '.pdb'
        if(.not.self%exists)stop 'Cannot write non existent atoms type; simple_atoms :: writePDB'
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
        call fileio_errmsg('writepdb; simple_atoms opening '//trim(fname), io_stat)
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
                else
                    atom_field(1:6) = 'ATOM  '
                endif
                write(pdbstr,pdbfmt)atom_field,self%num(ind),self%name(ind),self%altloc(ind),&
                    self%resname(ind),self%chain(ind), self%resnum(ind), self%icode(ind), self%xyz(ind,:),&
                    self%occupancy(ind), self%beta(ind)
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

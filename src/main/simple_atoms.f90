! atomic structures and pdb parser
module simple_atoms
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_strings
use simple_fileio
use simple_defs_atoms
implicit none

public :: atoms
private
#include "simple_local_flags.inc"

character(len=78), parameter :: pdbfmt          = "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10x,A2)"  ! custom 3.3
character(len=78), parameter :: pdbfmt_long     = "(A5,I6,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10x,A2)"  ! custom 3.3
character(len=74), parameter :: pdbfmt_read     = "(A11,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)"           ! custom 3.3
character(len=78), parameter :: pdbfmt_longread = "(A11,1X,A4,A1,A3,1X,A1,I1,A4,3X,3F8.3,2F6.2,10x,A2)"    ! custom 3.3

!>  \brief type for dealing with atomic structures
type :: atoms
    private
    integer                       :: n = 0
    character(len=4), allocatable :: name(:)
    character(len=1), allocatable :: altloc(:)
    character(len=1), allocatable :: chain(:)
    character(len=1), allocatable :: icode(:)
    character(len=3), allocatable :: resname(:)
    character(len=2), allocatable :: element(:)
    real,             allocatable :: charge(:)
    real,             allocatable :: xyz(:,:) ! dim1 -> 1:#atms, dime2 -> 1:3 Cartesian dims
    real,             allocatable :: mw(:)
    real,             allocatable :: occupancy(:)
    real,             allocatable :: beta(:)
    real,             allocatable :: radius(:)
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
    procedure          :: copy
    ! CHECKER
    procedure          :: element_exists
    ! GETTERS/SETTERS
    procedure          :: does_exist
    procedure          :: get_beta
    procedure          :: get_n
    procedure          :: get_name
    procedure          :: get_element
    procedure          :: get_coord
    procedure          :: get_num
    procedure          :: get_atomicnumber
    procedure          :: get_radius
    procedure          :: set_coord
    procedure          :: set_chain
    procedure          :: set_name
    procedure          :: set_element
    procedure          :: set_num
    procedure          :: set_beta
    procedure          :: set_resnum
    procedure          :: set_occupancy
    ! I/O
    procedure          :: print_atom
    procedure          :: writepdb
    ! CALCULATORS
    procedure          :: guess_element
    procedure, private :: guess_an_element
    procedure, private :: Z_and_radius_from_name
    procedure          :: get_geom_center
    procedure          :: convolve
    procedure          :: geometry_analysis_pdb
    procedure          :: find_masscen
    ! MODIFIERS
    procedure          :: translate
    procedure          :: rotate
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
        character(len=6)      :: atom_field
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
            if( len_trim(line) < 68 )then
                read(line,pdbfmt_read, iostat=io_stat)elevenfirst, self%name(i), self%altloc(i),&
                &self%resname(i), self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
                &self%occupancy(i), self%beta(i)
            else
                read(line,pdbfmt_longread, iostat=io_stat)elevenfirst, self%name(i), self%altloc(i),&
                &self%resname(i), self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
                &self%occupancy(i), self%beta(i), self%element(i)
            endif
            self%num(i) = num
            call fileiochk('new_from_pdb; simple_atoms error reading line '//trim(fname), io_stat)
            self%het(i) = atom_field == 'HETATM'
        enddo
        call self%guess_element
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
        integer,           intent(in)    :: n
        logical, optional, intent(in)    :: dummy
        integer :: i
        logical :: ddummy
        ddummy = .false.
        if(present(dummy)) ddummy = dummy
        call self%kill
        allocate(self%name(n), self%chain(n), self%resname(n), self%xyz(n,3), self%mw(n),&
            self%occupancy(n), self%beta(n), self%num(n), self%Z(n), self%het(n), self%icode(n),&
            self%altloc(n), self%resnum(n), self%element(n), self%radius(n))
        self%name(:)    = '    '
        self%resname(:) = '   '
        self%chain(:)   = ' '
        self%altloc(:)  = ' '
        self%icode(:)   = ' '
        self%element(:) = '  '
        self%mw        = 0.
        self%xyz       = 0.
        self%beta      = 0.
        self%occupancy = 0.
        self%radius    = 0.
        self%num    = 0
        self%resnum = 0
        self%Z      = 0
        self%n      = n
        self%het    = .false.
        if(ddummy)then
            do i=1,self%n
                self%name(i)      = ' X  '
                self%resname(i)   = ' X '
                self%chain(i)     = 'A'
                self%beta(i)      = 1.
                self%occupancy(i) = 1.
                self%radius(i)    = 1.
                self%num(i)       = i
                self%resnum(i)    = 1
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
        self%element    = self_in%element
        self%radius     = self_in%radius
    end subroutine copy

    ! CHECKERS

    logical function element_exists( self, element)
        class(atoms),     intent(inout) :: self
        character(len=*), intent(in)    :: element
        select case(len_trim(element))
        case(1)
            element_exists = exists(element(1:1)//' ')
        case(2)
            element_exists = exists(element(1:2))
        case(3)
            element_exists = exists(element(1:1)//' ').and.exists(element(2:3))
            if(.not.element_exists) element_exists = exists(element(2:3)).and.exists(element(3:3)//' ')
        case(4)
            element_exists = exists(element(1:2)).and.exists(element(3:4))
        case DEFAULT
            element_exists = .false.
            THROW_WARN('Non complying format; atoms%element_exists : '//trim(element))
        end select
        if( .not.element_exists) THROW_WARN('Unknown element: '//trim(element))
        contains

            logical function exists(el)
                character(len=2), intent(in) :: el
                real    :: r
                integer :: Z
                call self%Z_and_radius_from_name(el//'  ', Z, r)
                exists = Z /= 0
            end function exists

    end function element_exists

    ! GETTERS / SETTERS

    logical function does_exist( self )
        class(atoms), intent(inout) :: self
        does_exist = self%exists
    end function does_exist

    function get_beta( self, i) result(beta)
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        real :: beta
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_beta')
        beta = self%beta(i)
    end function get_beta

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
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_name')
        get_name = self%name(i)
    end function get_name

    character(len=2) function get_element( self, i )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_element')
        get_element = self%element(i)
    end function get_element

    integer function get_atomicnumber( self, i )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_atomicnumber')
        get_atomicnumber = self%Z(i)
    end function get_atomicnumber

    real function get_radius( self, i )
        class(atoms), intent(in) :: self
        integer,      intent(in) :: i
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; get_atomicnumber')
        get_radius = self%radius(i)
    end function get_radius

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
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_chain')
        self%chain(i) = upperCase(chain)
    end subroutine set_chain

    subroutine set_name( self, i, name )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=4), intent(in)    :: name
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_name')
        self%name(i) = upperCase(name)
        call self%guess_an_element(i)
    end subroutine set_name

    subroutine set_element( self, i, element )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=2), intent(in)    :: element
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_element')
        self%element(i) = upperCase(element(1:2))
        call self%guess_an_element(i)
    end subroutine set_element

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

    subroutine set_occupancy( self, i, occupancy )
        class(atoms),     intent(inout) :: self
        integer,          intent(in)    :: i
        real,             intent(in)    :: occupancy
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; set_occupancy')
        self%occupancy(i) = occupancy
    end subroutine set_occupancy

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
        write(logfhandle,'(2F6.2,1X)',advance='no')self%occupancy(i), self%beta(i)
        write(logfhandle,'(A2)',advance='yes')self%element(i)
    end subroutine print_atom

    ! I/O

    subroutine writepdb( self, fbody )
        class(atoms),     intent(in) :: self
        character(len=*), intent(in) :: fbody
        character(len=LONGSTRLEN) :: fname
        integer          :: i, funit, io_stat
        logical          :: long
        fname = trim(adjustl(fbody)) // '.pdb'
        long  = self%n >= 99999
        if(.not.self%exists) THROW_HARD('Cannot write non existent atoms type; writePDB')
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
        call fileiochk('writepdb; simple_atoms opening '//trim(fname), io_stat)
        do i = 1, self%n
            if(self%het(i))then
                write(funit,pdbfmt)'HETATM',self%num(i),self%name(i),self%altloc(i),&
                    self%resname(i),self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
                    self%occupancy(i), self%beta(i),' '
            else
                if( long )then
                    write(funit,pdbfmt_long)'ATOM ',self%num(i),self%name(i),self%altloc(i),&
                        self%resname(i),self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
                        self%occupancy(i), self%beta(i), self%element(i)
                else
                    write(funit,pdbfmt)'ATOM  ',self%num(i),self%name(i),self%altloc(i),&
                        self%resname(i),self%chain(i), self%resnum(i), self%icode(i), self%xyz(i,:),&
                        self%occupancy(i), self%beta(i), self%element(i)
                endif
            endif

        enddo
        call fclose(funit, errmsg='writepdb; simple_atoms closing '//trim(fname))
    end subroutine writepdb

    ! CALCULATORS

    subroutine guess_element( self )
        class(atoms), intent(inout) :: self
        integer :: i
        if( .not.self%exists ) THROW_HARD('Non-existent atoms object! guess_atomic_number')
        do i=1,self%n
            call self%guess_an_element(i)
        enddo
    end subroutine guess_element

    subroutine guess_an_element( self, i )
        class(atoms), intent(inout) :: self
        integer,      intent(in)    :: i
        real    :: r
        integer :: z
        if( .not.self%exists ) THROW_HARD('Non-existent atoms object! guess_atomic_number')
        if(i.lt.1 .or. i.gt.self%n) THROW_HARD('index out of range; guess_an_element')
        if( len_trim(self%element(i)) == 0 )then
            call self%Z_and_radius_from_name(self%name(i),z,r,self%element(i))
        else
            call self%Z_and_radius_from_name(self%element(i)//'  ',z,r)
        endif
        if( z == 0 )then
            THROW_WARN('Unknown atom '//int2str(i)//' : '//trim(self%name(i))//' - '//self%element(i))
        else
            self%radius(i) = r
            self%Z(i) = z
        endif
    end subroutine guess_an_element

    subroutine Z_and_radius_from_name( self, name, Z, r, el )
        class(atoms),     intent(inout) :: self
        character(len=4), intent(in)    :: name
        integer,          intent(out)   :: Z
        real,             intent(out)   :: r
        character(len=4)                :: uppercase_name
        character(len=2)                :: uppercase_el
        character(len=2), optional, intent(inout) :: el
        uppercase_name = upperCase(name)
        uppercase_el   = trim(adjustl(uppercase_name))
        call get_element_Z_and_radius(uppercase_el,Z,r)
        if( Z == 0 ) THROW_HARD('Unknown element: '//uppercase_el)
        if( present(el) ) el = uppercase_el
    end subroutine Z_and_radius_from_name

    function get_geom_center(self) result(center)
        class(atoms), intent(in) :: self
        real :: center(3)
        center(1) = sum(self%xyz(:,1)) / real(self%n)
        center(2) = sum(self%xyz(:,2)) / real(self%n)
        center(3) = sum(self%xyz(:,3)) / real(self%n)
    end function get_geom_center

    ! simulate electrostatic potential from elastic scattering
    ! Using 5-gaussian atomic scattering factors from Rullgard et al, J of Microscopy, 2011
    ! and parametrization from Peng, Acta Cryst, 1996, A52, Table 1 (also in ITC)
    subroutine convolve( self, vol, cutoff, lp )
        use simple_image, only: image
        class(atoms),   intent(in)    :: self
        class(image),   intent(inout) :: vol
        real,           intent(in)    :: cutoff
        real, optional, intent(in)    :: lp
        real, parameter   :: C = 2132.79 ! eq B.6, conversion to eV
        real, parameter   :: fourpisq = 4.*PI*PI
        real, allocatable :: rmat(:,:,:)
        real    :: lp_here, a(5),b(5),aterm(5), xyz(3), smpd,r2,bfac,rjk2,cutoffsq,D,E
        integer :: bbox(3,2),ldim(3),pos(3),i,j,k,l,jj,kk,z,icutoff
        if( .not.vol%is_3d() .or. vol%is_ft() ) THROW_HARD('Only for real-space volumes')
        smpd    = vol%get_smpd()
        ldim    = vol%get_ldim()
        lp_here = 2.*smpd
        if( present(lp) ) lp_here = max(lp,lp_here)
        bfac = (4.*lp_here)**2.
        D    = sqrt(TWOPI) * 0.425 * lp_here
        E    = 0.5 * lp_here*lp_here
        allocate(rmat(ldim(1),ldim(2),ldim(3)),source=0.)
        icutoff  = ceiling(cutoff/smpd)
        cutoffsq = cutoff*cutoff
        !$omp parallel do default(shared) private(i,z,a,b,aterm,xyz,pos,bbox,j,k,l,r2,rjk2,jj,kk)&
        !$omp proc_bind(close) reduction(+:rmat)
        do i = 1,self%n
            z = self%Z(i)
            select case(z)
            case(1) ! H
                a = [0.0349, 0.1201, 0.1970, 0.0573, 0.1195]
                b = [0.5347, 3.5867, 12.3471, 18.9525, 38.6269]
            case(2) ! HE
                a = [0.0317, 0.0838, 0.1526, 0.1334, 0.0164]
                b = [0.2507, 1.4751, 4.4938, 12.6646, 31.1653]
            case(3) ! LI
                a = [0.0750, 0.2249, 0.5548, 1.4954, 0.9354]
                b = [0.3864, 2.9383, 15.3829, 53.5545, 138.7337]
            case(4) ! BE
                a = [0.0780, 0.2210, 0.6740, 1.3867, 0.6925]
                b = [0.3131, 2.2381, 10.1517, 30.9061, 78.3273]
            case(5) ! B
                a = [0.0909, 0.2551, 0.7738, 1.2136, 0.4606]
                b = [0.2995, 2.1155, 8.3816, 24.1292, 63.1314]
            case(6) ! C
                a = [0.0893, 0.2563, 0.7570, 1.0487, 0.3575]
                b = [0.2465, 1.7100, 6.4094,18.6113,50.2523]
            case(7) ! N
                a = [0.1022, 0.3219, 0.7982, 0.8197, 0.1715]
                b = [0.2451, 1.7481, 6.1925,17.3894,48.1431]
            case(8) ! O
                a = [0.0974, 0.2921, 0.6910,  0.6990, 0.2039]
                b = [0.2067, 1.3815, 4.6943, 12.7105,32.4726]
            case(9) ! F
                a = [0.1083, 0.3175, 0.6487, 0.5846, 0.1421]
                b = [0.2057, 1.3439, 4.2788, 11.3932, 28.7881]
            case(10) ! NE
                a = [0.1269, 0.3535, 0.5582, 0.4674, 0.1460]
                b = [0.2200, 1.3779, 4.0203, 9.4934, 23.1278]
            case(11) ! NA
                a = [0.2142, 0.6853, 0.7692, 1.6589, 1.4482]
                b = [0.3334, 2.3446, 10.0830, 48.3037, 138.2700]
            case(12) ! MG
                a = [0.2314, 0.6866, 0.9677, 2.1882, 1.1339]
                b = [0.3278, 2.2720, 10.9241, 39.2898, 101.9748]
            case(13) ! AL
                a = [0.2390, 0.6573, 1.2011, 2.5586, 1.2312]
                b = [0.3138, 2.1063, 10.4163, 34.4552, 98.5344]
            case(14) ! SI
                a = [0.2519, 0.6372, 1.3795, 2.5082, 1.0500]
                b = [0.3075, 2.0174, 9.6746, 29.3744, 80.4732]
            case(15) ! P
                a = [0.2548, 0.6106, 1.4541, 2.3204, 0.8477]
                b = [0.2908, 1.8740, 8.5176, 24.3434, 63.2996]
            case(16) ! S
                a = [0.2497, 0.5628, 1.3899, 2.1865, 0.7715]
                b = [0.2681, 1.6711, 7.0267, 19.5377, 50.3888]
            case(17) ! CL
                a = [0.2443, 0.5397, 1.3919, 2.0197, 0.6621]
                b = [0.2468, 1.5242, 6.1537, 16.6687, 42.3086]
            case(18) ! AR
                a = [0.2385, 0.5017, 1.3428, 1.8899, 0.6079]
                b = [0.2289, 1.3694, 5.2561, 14.0928, 35.5361]
            case(19) ! K
                a = [0.4115, 1.4031, 2.2784, 2.6742, 2.2162]
                b = [0.3703, 3.3874, 13.1029, 68.9592, 194.4329]
            case(20) ! CA
                a = [0.4054, 1.3880, 2.1602, 3.7532, 2.2063]
                b = [0.3499, 3.0991, 11.9608, 53.9353, 142.3892]
            case(21) ! SC
                a = [0.3787, 1.2181, 2.0594, 3.2618, 2.3870]
                b = [0.3133, 2.5856, 9.5813, 41.7688, 116.7282]
            case(22) ! TI
                a = [0.3825, 1.2598, 2.0008, 3.0617, 2.0694]
                b = [0.3040, 2.4863, 9.2783, 39.0751, 109.4583]
            case(23) ! V
                a = [0.3876, 1.2750, 1.9109, 2.8314, 1.8979]
                b = [0.2967, 2.3780, 8.7981, 35.9528, 101.7201]
            case(24) ! CR
                a = [0.4046, 1.3696, 1.8941, 2.0800, 1.2196]
                b = [0.2986, 2.3958, 9.1406, 37.4701, 113.7121]
            case(25) ! MN
                a = [0.3796, 1.2094, 1.7815, 2.5420, 1.5937]
                b = [0.2699, 2.0455, 7.4726, 31.0604, 91.5622]
            case(26) ! Fe
                a = [0.3946, 1.2725, 1.7031,  2.3140,  1.4795]
                b = [0.2717, 2.0443, 7.6007, 29.9714, 86.2265]
            case(27) ! CO
                a = [0.4118, 1.3161, 1.6493, 2.1930, 1.2830]
                b = [0.2742, 2.0372, 7.7205, 29.9680, 84.9383]
            case(28) ! NI
                a = [0.3860, 1.1765, 1.5451, 2.0730, 1.3814]
                b = [0.2478, 1.7660, 6.3107, 25.2204, 74.3146]
            case(29) ! CU
                a = [0.4314, 1.3208, 1.5236, 1.4671, 0.8562]
                b = [0.2694, 1.9223, 7.3474, 28.9892, 90.6246]
            case(30) ! CU
                a = [0.4288, 1.2646, 1.4472, 1.8294, 1.0934]
                b = [0.2593, 1.7998, 6.7500, 25.5860, 73.5284]
            case(31) ! GA
                a = [0.4818, 1.4032, 1.6561, 2.4605, 1.1054]
                b = [0.2825, 1.9785, 8.7546, 32.5238, 98.5523]
            case(32) ! GE
                a = [0.4655, 1.3014, 1.6088, 2.6998, 1.3003]
                b = [0.2647, 1.7926, 7.6071, 26.5541, 77.5238]
            case(33) ! AS
                a = [0.4517, 1.2229, 1.5852, 2.7958, 1.2638]
                b = [0.2493, 1.6436, 6.8154, 22.3681, 62.0390]
            case(34) ! Se
                a = [0.4477, 1.1678, 1.5843, 2.8087, 1.1956]
                b = [0.2405, 1.5442, 6.3231, 19.4610, 52.0233]
            case(35) ! BR
                a = [0.4798, 1.1948, 1.8695, 2.6953, 0.8203]
                b = [0.2504, 1.5963, 6.9653, 19.8492, 50.3233]
            case(36) ! KR
                a = [0.4546, 1.0993, 1.7696, 2.7068, 0.8672]
                b = [0.2309, 1.4279, 5.9449, 16.6752, 42.2243]
            case(37) ! RB
                a = [1.0160, 2.8528, 3.5466, -7.7804, 12.1148]
                b = [0.4853, 5.0925, 25.7851, 130.4515, 138.6775]
            case(38) ! SR
                a = [0.6703, 1.4926, 3.3368, 4.4600, 3.1501]
                b = [0.3190, 2.2287, 10.3504, 52.3291, 151.2216]
            case(39) ! Y
                a = [0.6894, 1.5474, 3.2450, 4.2126, 2.9764]
                b = [0.3189, 2.2904, 10.0062, 44.0771, 125.0120]
            case(40) ! ZR
                a = [0.6719, 1.4684, 3.1668, 3.9557, 2.8920]
                b = [0.3036, 2.1249, 8.9236, 36.8458, 108.2049]
            case(41) ! NB
                a = [0.6123, 1.2677, 3.0348, 3.3841, 2.3683]
                b = [0.2709, 1.7683, 7.2489, 27.9465, 98.5624]
            case(42) ! MO
                a = [0.6773, 1.4798, 3.1788, 3.0824, 1.8384]
                b = [0.2920, 2.0606, 8.1129, 30.5336, 100.0658]
            case(43) ! TC
                a = [0.7082, 1.6392, 3.1993, 3.4327, 1.8711]
                b = [0.2976, 2.2106, 8.5246, 33.1456, 96.6377]
            case(44) ! RU
                a = [0.6735, 1.4934, 3.0966, 2.7254, 1.5597]
                b = [0.2773, 1.9716, 7.3249, 26.6891, 90.5581]
            case(45) ! RH
                a = [0.6413, 1.3690, 2.9854, 2.6952, 1.5433]
                b = [0.2580, 1.7721, 6.3854, 23.2549, 85.1517]
            case(46) ! PD
                a = [0.5904, 1.1775, 2.6519, 2.2875, 0.8689]
                b = [0.2324, 1.5019, 5.1591, 15.5428, 46.8213]
            case(47) ! AG
                a = [0.6377, 1.3790, 2.8294, 2.3631, 1.4553]
                b = [0.2466, 1.674, 5.7656, 20.0943, 76.7372]
            case(48) ! CD
                a = [0.6364, 1.4247, 2.7802, 2.5973, 1.7886]
                b = [0.2407, 1.6823, 5.6588, 20.7219, 69.1409]
            case(49) ! IN
                a = [0.6768, 1.6589, 2.7740, 3.1835, 2.1326]
                b = [0.2522, 1.8545, 6.2936, 25.1457, 84.5448]
            case(50) ! SN
                a = [0.7224, 1.9610, 2.7161, 3.5603, 1.8972]
                b = [0.2651, 2.0604, 7.3011, 27.5493, 81.3349]
            case(51) ! SB
                a = [0.7106, 1.9247, 2.6149, 3.8322, 1.8899]
                b = [0.2562, 1.9646, 6.8852, 24.7648, 68.9168]
            case(52) ! TE
                a = [0.6947, 1.8690, 2.5356, 4.0013, 1.8955]
                b = [0.2459, 1.8542, 6.4411, 22.1730, 59.2206]
            case(53) ! I
                a = [0.7047, 1.9484, 2.5940, 4.1526, 1.5057]
                b = [0.2455, 1.8638, 6.7639, 21.8007, 56.4395]
            case(54) ! XE
                a = [0.6737, 1.7908, 2.4129, 4.2100, 1.7058]
                b = [0.2305, 1.6890, 5.8218, 18.3928, 47.2496]
            case(55) ! CS
                a = [1.2704, 3.8018, 5.6618, 0.9205, 4.8105]
                b = [0.4356, 4.2058, 23.4342, 136.778393, 171.7561]
            case(56) ! BA
                a = [0.9049, 2.6076, 4.8498, 5.1603, 4.7388]
                b = [0.3066, 2.4363, 12.1821, 54.6135, 161.9978]
            case(57) ! LA
                a = [0.8405, 2.3863, 4.6139, 5.1514, 4.7949]
                b = [0.2791, 2.1410, 10.3400, 41.9148, 132.0204]
            case(58) ! CE
                a = [0.8551, 2.3915, 4.5772, 5.0278, 4.5118]
                b = [0.2805, 2.1200, 10.1808, 42.0633, 130.9893]
            case(59) ! PR
                a = [0.9096, 2.5313, 4.5266, 4.6376, 4.3690]
                b = [0.2939, 2.2471, 10.8266, 48.8842, 147.6020]
            case(60) ! ND
                a = [0.8807, 2.4183, 4.4448, 4.6858, 4.1725]
                b = [0.2802, 2.0836, 10.0357, 47.4506, 146.9976]
            case(61) ! PM
                a = [0.9471, 2.5463, 4.3523, 4.4789, 3.9080]
                b = [0.2977, 2.2276, 10.5762, 49.3619, 145.3580]
            case(62) ! SM
                a = [0.9699, 2.5837, 4.2778, 4.4575, 3.5985]
                b = [0.3003, 2.2447, 10.6487, 50.7994, 146.4179]
            case(63) ! EU
                a = [0.8694, 2.2413, 3.9196, 3.9694, 4.5498]
                b = [0.2653, 1.8590, 8.3998, 36.7397, 125.7089]
            case(64) ! GD
                a = [0.9673, 2.4702, 4.1148, 4.4972, 3.2099]
                b = [0.2909, 2.1014, 9.7067, 43.4270, 125.9474]
            case(65) ! TB
                a = [0.9325, 2.3673, 3.8791, 3.9674, 3.7996]
                b = [0.2761, 1.9511, 8.9296, 41.5937, 131.0122]
            case(66) ! DY
                a = [0.9505, 2.3705, 3.8218, 4.0471, 3.4451]
                b = [0.2773, 1.9469, 8.8862, 43.0938, 133.1396]
            case(67) ! HO
                a = [0.9248, 2.2428, 3.6182, 3.7910, 3.7912]
                b = [0.2660, 1.8183, 7.9655, 33.1129, 101.8139]
            case(68) ! ER
                a = [1.0373, 2.4824, 3.6558, 3.8925, 3.0056]
                b = [0.2944, 2.0797, 9.4156, 45.8056, 132.7720]
            case(69) ! TM
                a = [1.0075, 2.3787, 3.5440, 3.6932, 3.1759]
                b = [0.2816, 1.9486, 8.7162, 41.8420, 125.0320]
            case(70) ! YB
                a = [1.0347, 2.3911, 3.4619, 3.6556, 3.0052]
                b = [0.2855, 1.9679, 8.7619, 42.3304, 125.6499]
            case(71) ! LU
                a = [0.9927, 2.2436, 3.3554, 3.7813, 3.0994]
                b = [0.2701, 1.8073, 7.8112, 34.4849, 103.3526]
            case(72) ! HF
                a = [1.0295, 2.2911, 3.4110, 3.9497, 2.4925]
                b = [0.2761, 1.8625, 8.0961, 34.2712, 98.5295]
            case(73) ! TA
                a = [1.0190, 2.2291, 3.4097, 3.9252, 2.2679]
                b = [0.2694, 1.7962, 7.6944, 31.0942, 91.1089]
            case(74) ! W
                a = [0.9853, 2.1167, 3.3570, 3.7981, 2.2798]
                b = [0.2569, 1.6745, 7.0098, 26.9234, 81.3910]
            case(75) ! RE
                a = [0.99141, 2.0858, 3.4531, 3.38812, 1.8526]
                b = [0.2548, 1.6518, 6.8845, 26.7234, 81.7215]
            case(76) ! OS
                a = [0.9813, 2.0322, 3.3665, 3.6235, 1.9741]
                b = [0.2487, 1.5973, 6.4737, 23.2817, 70.9254]
            case(77) ! IR
                a = [1.0194, 2.0645, 3.4425, 3.4914, 1.6976]
                b = [0.2554, 1.6475, 6.5966, 23.2269, 70.0272]
            case(78) ! Pt
                a = [0.9148, 1.8096, 3.2134,  3.2953,  1.5754]
                b = [0.2263, 1.3813, 5.3243, 17.5987, 60.0171]
            case(79) ! Au
                a = [0.9674, 1.8916, 3.3993,  3.0524,  1.2607]
                b = [0.2358, 1.4712, 5.6758, 18.7119, 61.5286]
            case(80) ! HG
                a = [1.0033, 1.9469, 3.4396,  3.1548,  1.4180]
                b = [0.2413, 1.5298, 5.8009, 19.4520, 60.5753]
            case(81) ! TL
                a = [1.0689, 2.1038, 3.6039,  3.4927,  1.8283]
                b = [0.2540, 1.6715, 6.3509, 23.1531, 78.7099]
            case(82) ! PB
                a = [1.0891, 2.1867, 3.6160, 3.8031, 1.8994]
                b = [0.2552, 1.7174, 6.5131, 23.9170, 74.7039]
            case(83) ! BI
                a = [1.1007, 2.2306, 3.5689, 4.1549, 2.0382]
                b = [0.2546, 1.7351, 6.4948, 23.6464, 70.3780]
            case(84) ! PO
                a = [1.1568, 2.4353, 3.6459, 4.4064, 1.7179]
                b = [0.2648, 1.8786, 7.1749, 25.1766, 69.2821]
            case(85) ! AT
                a = [1.0909, 2.1976, 3.3831, 4.6700, 2.1277]
                b = [0.2466, 1.6707, 6.0197, 20.7657, 57.2663]
            case(86) ! RN
                a = [1.0756, 2.1630, 3.3178, 4.8852, 2.0489]
                b = [0.2402, 1.6169, 5.7644, 19.4568, 52.5009]
            case(87) ! FR
                a = [1.4282, 3.5081, 5.6767, 4.1964, 3.8946]
                b = [0.3183, 2.6889, 13.4816, 54.3866, 200.8321]
            case(88) ! RA
                a = [1.3127, 3.1243, 5.2988, 5.3891, 5.4133]
                b = [0.2887, 2.2897, 10.8276, 43.5389, 145.6109]
            case(89) ! AC
                a = [1.3128, 3.1021, 5.3385, 5.9611, 4.7562]
                b = [0.2861, 2.2509, 10.5287, 41.7796, 128.2973]
            case(90) ! TH
                a = [1.2553, 2.9178, 5.0862, 6.1206, 4.7122]
                b = [0.2701, 2.0636, 9.3051, 34.5977, 107.9200]
            case(91) ! PA
                a = [1.3218, 3.1444, 5.4371, 5.6444, 4.0107]
                b = [0.2827, 2.2250, 10.2454, 41.1162, 124.4449]
            case(92) ! U
                a = [1.3382, 3.2043, 5.4558, 5.4839, 3.6342]
                b = [0.2838, 2.2452, 10.2519, 41.7251, 124.9023]
            case(93) ! NP
                a = [1.5193, 4.0053, 6.5327, -0.1402, 6.7489]
                b = [0.3213, 2.8206, 14.8878, 68.9103, 81.7257]
            case(94) ! PU
                a = [1.3517, 3.2937, 5.3213, 4.6466, 3.5714]
                b = [0.2813, 2.2418, 9.9952, 42.7939, 132.1739]
            case(95) ! AM
                a = [1.2135, 2.7962, 4.7545, 4.5731, 4.4786]
                b = [0.2483, 1.8437, 7.5421, 29.3841, 112.4579]
            case(96) ! CM
                a = [1.2937, 3.1100, 5.0393, 4.7546, 3.5031]
                b = [0.2638, 2.0341, 8.7101, 35.2992, 109.4972]
            case DEFAULT
                cycle
            end select
            b = b + bfac        ! eq B.6
            aterm = a/b**1.5    ! eq B.6
            xyz   = self%xyz(i,:)/smpd
            pos   = floor(xyz)
            bbox(:,1) = pos   - icutoff
            bbox(:,2) = pos+1 + icutoff
            if( any(bbox(:,2) < 1) )      cycle
            if( any(bbox(:,1) > ldim(1)) )cycle
            where( bbox < 1 ) bbox = 1
            where( bbox > ldim(1) ) bbox = ldim(1)
            do j = bbox(1,1),bbox(1,2)
                jj = j-1
                do k = bbox(2,1),bbox(2,2)
                    kk = k-1
                    rjk2 = sum((smpd*(xyz(1:2)-real([jj,kk])))**2.)
                    if(rjk2 > cutoffsq) cycle
                    do l = bbox(3,1),bbox(3,2)
                        r2 = rjk2 + (smpd*(xyz(3)-real(l-1)))**2.
                        rmat(j,k,l) = rmat(j,k,l) + epot(r2,aterm,b)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        call vol%set_rmat(rmat,.false.)
        deallocate(rmat)
    contains
        ! potential assuming static atoms (eq B.6)
        real function epot(r2,aterm,b)
            real, intent(in) :: r2, aterm(5),b(5)
            epot = C * sum( aterm * exp(-fourpisq*r2/b) )
        end function epot

        ! single gaussian convolution, unused
        elemental real function egau(r2,zi)
            real,    intent(in) :: r2
            integer, intent(in) :: zi
            egau = real(zi) / D * exp(-r2*E)
        end function egau

    end subroutine convolve

    subroutine geometry_analysis_pdb(self, pdbfile, thresh)
      use simple_math, only : plane_from_points
      class(atoms),     intent(inout) :: self
      character(len=*), intent(in)    :: pdbfile   ! all the atomic positions
      real, optional,   intent(in)    :: thresh    ! for belonging
      character(len=2)     :: element
      type(atoms)          :: init_atoms, final_atoms
      real,    allocatable :: radii(:),line(:,:), plane(:,:,:),points(:,:), distances_totheplane(:), distances_totheline(:)
      real,    allocatable :: w(:),v(:,:),d(:),pointsTrans(:,:)
      logical, allocatable :: flag(:) ! flags the atoms belonging to the plane/column
      integer, parameter   :: N_DISCRET = 500
      integer :: i, n, n_tot, t, s, filnum, io_stat, cnt_intersect, cnt
      real    :: atom1(3), atom2(3), atom3(3), dir_1(3), dir_2(3), vec(3), m(3), dist_plane, dist_line
      real    :: t_vec(N_DISCRET), s_vec(N_DISCRET), denominator, centroid(3), prod(3), tthresh
      call init_atoms%new(pdbfile)
      n = init_atoms%get_n()
      if(present(thresh)) then
          tthresh = thresh
      else
          tthresh = 1.2*sum(self%radius)/real(self%get_n())  ! avg of the radii*1.2
      endif
      if(n < 2 .or. n > 3 ) THROW_HARD('Inputted pdb file contains the wrong number of atoms!; geometry_analysis_pdb')
      do i = 1, N_DISCRET/2
          t_vec(i) = -real(i)/10.
      enddo
      t_vec(N_DISCRET/2+1:N_DISCRET) = -t_vec(1:N_DISCRET/2)
      s_vec(:) = t_vec(:)
      n_tot = self%n
      ! fetch thoretical radius
      element = self%element(1) ! pick the first atom (should be heterogeneous)
      allocate(flag(n_tot), source = .false.)
      if(n == 2) then
        write(logfhandle,*)'COLUMN IDENTIFICATION, INITIATION'
        allocate(line(3, N_DISCRET), source = 0.)
        atom1(:) = init_atoms%get_coord(1)
        atom2(:) = init_atoms%get_coord(2)
        dir_1 = atom1-atom2
        do t = 1, N_DISCRET
          line(1,t) = atom1(1) + t_vec(t)* dir_1(1)
          line(2,t) = atom1(2) + t_vec(t)* dir_1(2)
          line(3,t) = atom1(3) + t_vec(t)* dir_1(3)
        enddo
        ! calculate how many atoms does the line intersect and flag them
        do i = 1, n_tot
            do t = 1, N_DISCRET
              dist_line = euclid(self%xyz(i,:3),line(:3,t))
              if(dist_line <= 0.6*tthresh) then ! it intersects atoms
                   flag(i) = .true. !flags also itself
               endif
            enddo
        enddo
       ! generate pdb file for visualisation
       cnt_intersect = 0
       call final_atoms%new(count(flag), dummy=.true.)
       do i = 1, n_tot
           if(flag(i)) then
               cnt_intersect = cnt_intersect + 1
               call final_atoms%set_name(cnt_intersect,self%name(i))
               call final_atoms%set_element(cnt_intersect,self%element(i))
               call final_atoms%set_coord(cnt_intersect,(self%xyz(i,:3)))
               call final_atoms%set_occupancy(cnt_intersect,self%occupancy(i))
               call final_atoms%set_beta(cnt_intersect,self%beta(i))
           endif
       enddo
       call final_atoms%writePDB('AtomColumn')
       call final_atoms%kill
       ! Find the line that best fits the atoms
       allocate(points(3,count(flag)), source = 0.)
       cnt = 0
       do i = 1,n_tot
         if(flag(i)) then
           cnt = cnt + 1
           points(:3,cnt) = self%xyz(i,:3)
         endif
       enddo
       ! calculate centroid of the points
       centroid = sum(points(:,:), dim = 2)/real(count(flag))
       ! svd fit
       allocate(pointsTrans(count(flag),3), source = 0.) ! because svdcmp modifies its input
       ! translate
       do i = 1, count(flag)
         pointsTrans(i,:3) = points(:3,i) - centroid(:3)
       enddo
       allocate(w(3), v(3,3), source = 0.)
       allocate(d(3), source = 0.)
       call svdcmp(pointsTrans,w,v)
       d = v(:,1)
       write(logfhandle, *) 'Directional vector of the line', d
       ! line
       ! line(1,t) = centroid(1) + t_vec(t)* d(1)
       ! line(2,t) = centroid(2) + t_vec(t)* d(2)
       ! line(3,t) = centroid(3) + t_vec(t)* d(3)
       ! calculate the distance to the points from the identified line
       allocate(distances_totheline(cnt), source = 0.)
       allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
       denominator = sqrt(d(1)**2+d(2)**2+d(3)**2)
       m = sum(self%xyz(:,:), dim=1) /real(self%n)
       cnt = count(flag)
       do i = 1, cnt
         vec  = centroid(:3)-points(:3,i)
         prod = cross(vec,d)
         distances_totheline(i) = sqrt(prod(1)**2+prod(2)**2+prod(3)**2)/denominator
         radii(i) = euclid(points(:,i), m)
       enddo
       ! it's already in A
       call fopen(filnum, file='Radii.csv', iostat=io_stat)
       write (filnum,*) 'r'
       do i = 1, cnt
         write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
       end do
       call fclose(filnum)
       call fopen(filnum, file='DistancesToTheLine.csv',iostat=io_stat)
       write (filnum,*) 'd'
       do i = 1, cnt
         write (filnum,'(A)', advance='yes') trim(real2str(distances_totheline(i)))
       end do
       call fclose(filnum)
      elseif(n == 3) then
        write(logfhandle,*)'PLANE IDENTIFICATION, INITIATION'
        atom1(:) = init_atoms%get_coord(1)
        atom2(:) = init_atoms%get_coord(2)
        atom3(:) = init_atoms%get_coord(3)
        dir_1 = atom1-atom2
        dir_2 = atom1-atom3
        allocate(plane(3, N_DISCRET, N_DISCRET), source = 0.)
        do t = 1, N_DISCRET
            do s = 1, N_DISCRET
                plane(1,t,s) = atom1(1) + t_vec(t)* dir_1(1) + s_vec(s)* dir_2(1)
                plane(2,t,s) = atom1(2) + t_vec(t)* dir_1(2) + s_vec(s)* dir_2(2)
                plane(3,t,s) = atom1(3) + t_vec(t)* dir_1(3) + s_vec(s)* dir_2(3)
            enddo
        enddo
        ! calculate how many atoms does the plane intersect and flag them
        do i = 1, n_tot
            do t = 1, N_DISCRET
              do s = 1, N_DISCRET
                dist_plane = euclid(self%xyz(i,:3),plane(:3,t,s))
                if(dist_plane <= 0.6*tthresh) then ! it intersects atoms i
                     flag(i) = .true. !flags also itself
                 endif
              enddo
            enddo
        enddo
        ! generate pdb for visualisation
        cnt_intersect = 0
        call final_atoms%new(count(flag), dummy=.true.)
        do i = 1, n_tot
            if(flag(i)) then
                cnt_intersect = cnt_intersect + 1
                call final_atoms%set_name(cnt_intersect,self%name(i))
                call final_atoms%set_element(cnt_intersect,self%element(i))
                call final_atoms%set_coord(cnt_intersect,(self%xyz(i,:3)))
                call final_atoms%set_occupancy(cnt_intersect,self%occupancy(i))
                call final_atoms%set_beta(cnt_intersect,self%beta(i))
            endif
        enddo
        call final_atoms%writePDB('AtomPlane')
        call final_atoms%kill
        allocate(points(3, count(flag)), source = 0.)
        ! calculate center of mass of the points
        m = sum(self%xyz(:,:), dim=1) /real(self%n)
        cnt = 0
        do i = 1, n_tot
            if(flag(i)) then
                cnt = cnt + 1
                points(:3,cnt) = self%xyz(i,:3)-m(:)
            endif
        enddo
        vec = plane_from_points(points)
        allocate(distances_totheplane(cnt), source = 0.)
        allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
        cnt = 0
        denominator = sqrt(vec(1)**2+vec(2)**2+1.)
        write(logfhandle,*) 'Normal vector: [', vec(1), ',', vec(2), ',', -1., ']'
        do i = 1, n_tot
            if(flag(i)) then
                cnt = cnt + 1
                ! formula for distance of a point to a plane
                distances_totheplane(cnt) = abs(vec(1)*points(1,cnt)+vec(2)*points(2,cnt)-points(3,cnt)+vec(3))/denominator
                radii(cnt) = euclid(self%xyz(i,:3), m)
            endif
        enddo
        ! it's already in A
        call fopen(filnum, file='Radii.csv', iostat=io_stat)
        write (filnum,*) 'r'
        do i = 1, cnt
          write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
        end do
        call fclose(filnum)
        call fopen(filnum, file='DistancesToThePlane.csv',iostat=io_stat)
        write (filnum,*) 'd'
        do i = 1, cnt
          write (filnum,'(A)', advance='yes') trim(real2str(distances_totheplane(i)))
        end do
        call fclose(filnum)
      endif
      call init_atoms%kill
      if(allocated(line))  deallocate(line)
      if(allocated(plane)) deallocate(plane)
    contains

      ! Compute the cross product of 2 3D real vectors
      function cross(a, b) result(c)
          real, intent(in) :: a(3),b(3)
          real :: c(3)
          c(1) = a(2) * b(3) - a(3) * b(2)
          c(2) = a(3) * b(1) - a(1) * b(3)
          c(3) = a(1) * b(2) - a(2) * b(1)
      end function cross
    end subroutine geometry_analysis_pdb

    function find_masscen( self ) result( m )
        class(atoms), intent(in) :: self
        real    :: m(3) ! mass center vector
        integer :: i
        m = 0.
        do i = 1, self%n ! #atms set in the constructor
            m = m + self%xyz(i,:)
        end do
        m = m / real(self%n)
    end function find_masscen

    ! subroutine shift2masscen( self, m )
    !     class(atoms), intent(in) :: self
    !     real, intent(in)         :: m(3)
    !     call translate(m)
    !
    ! end subroutine shift2masscen

    ! MODIFIERS

    subroutine translate( self, shift )
        class(atoms), intent(inout) :: self
        real,         intent(in)    :: shift(3)
        self%xyz(:,1) = self%xyz(:,1) + shift(1)
        self%xyz(:,2) = self%xyz(:,2) + shift(2)
        self%xyz(:,3) = self%xyz(:,3) + shift(3)
    end subroutine translate

    subroutine rotate( self, mat )
        class(atoms), intent(inout) :: self
        real,         intent(in)    :: mat(3,3)
        self%xyz = matmul(self%xyz,transpose(mat))
    end subroutine rotate

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
        if( allocated(self%element) )deallocate(self%element)
        if( allocated(self%radius) )deallocate(self%radius)
        self%n      = 0
        self%exists = .false.
    end subroutine kill

end module

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
    real,             allocatable :: xyz(:,:)
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
    procedure, private :: copy
    ! GETTERS/SETTERS
    procedure          :: does_exist
    procedure          :: get_beta
    procedure          :: get_n
    procedure          :: get_name
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
        integer :: i,alloc_stat
        logical :: ddummy
        ddummy = .false.
        if(present(dummy)) ddummy = dummy
        call self%kill
        allocate(self%name(n), self%chain(n), self%resname(n), self%xyz(n,3), self%mw(n),&
            self%occupancy(n), self%beta(n), self%num(n), self%Z(n), self%het(n), self%icode(n),&
            self%altloc(n), self%resnum(n), self%element(n), self%radius(n), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('new_instance :: simple_atoms', alloc_stat)
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
        get_name = self%name(i)
    end function get_name

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
        self%element(i) = upperCase(element)
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
        character(len=STDLEN) :: fname
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

            character(len=78) function pdbstr( ind )
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
                            self%occupancy(ind), self%beta(ind),self%element(ind)
                    else
                        atom_field(1:6) = 'ATOM  '
                        write(pdbstr,pdbfmt_long)atom_field,self%num(ind),self%name(ind),self%altloc(ind),&
                            self%resname(ind),self%chain(ind), self%resnum(ind), self%icode(ind), self%xyz(ind,:),&
                            self%occupancy(ind), self%beta(ind),self%element(ind)
                    endif
                endif
            end function pdbstr
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

    ! single covalent radii from Cordero, et al., 2008, "Covalent radii revisited"
    ! Dalton Trans. (21): 2832–2838. doi:10.1039/b801115j
    subroutine Z_and_radius_from_name( self, name, Z, r, el )
        class(atoms),     intent(inout) :: self
        character(len=4), intent(in)    :: name
        integer,          intent(out)   :: Z
        real,             intent(out)   :: r
        character(len=2), optional, intent(inout) :: el
        character(len=4) :: uppercase_name
        Z = 0
        r = 1.
        uppercase_name = upperCase(name)
        select case(trim(adjustl(uppercase_name)))
        ! organic
        case('H')
            Z = 1 ; r = 0.31
        case('C')
            Z = 6;  r = 0.76 ! sp3
        case('N')
            Z = 7;  r = 0.71
        case('O')
            Z = 8;  r = 0.66
        case('P')
            Z = 15; r = 1.07
        case('S')
            Z = 16; r = 1.05
        ! metals
        case('FE')
            Z = 26; r = 1.02
        case('PD')
            Z = 46; r = 1.12
        case('PT')
            Z = 78; r = 1.10
        case('AU')
            Z = 79; r = 1.23
        end select
        if( present(el) ) el = trim(adjustl(uppercase_name))
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
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_image, only: image
        class(atoms),   intent(in)    :: self
        class(image),   intent(inout) :: vol
        real,           intent(in)    :: cutoff
        real, optional, intent(in)    :: lp
        real, parameter   :: C = 2132.79 ! eq B.6, conversion to eV
        real, parameter   :: fourpisq = 4.*PI*PI
        real, allocatable :: rmat(:,:,:)
        real    :: a(5),b(5),aterm(5), xyz(3), smpd,r2,bfac,rjk2,cutoffsq
        integer :: bbox(3,2),ldim(3),pos(3),i,j,k,l,jj,kk,z,icutoff
        if( .not.vol%is_3d() .or. vol%is_ft() ) THROW_HARD('Only for real-space volumes')
        smpd = vol%get_smpd()
        ldim = vol%get_ldim()
        bfac = (4.*2.*smpd)**2.
        if( present(lp) ) bfac = max(bfac,(4.*lp)**2.)
        allocate(rmat(ldim(1),ldim(2),ldim(3)),source=0.)
        icutoff  = ceiling(cutoff/smpd)
        cutoffsq = cutoff*cutoff
        !$omp parallel do default(shared) private(i,z,a,b,aterm,xyz,pos,bbox,j,k,l,r2,rjk2,jj,kk)&
        !$omp proc_bind(close) reduction(+:rmat)
        do i = 1,self%n
            z = self%Z(i)
            select case(z)
            case(6) ! C
                a = [0.0893, 0.2563, 0.7570, 1.0487, 0.3575]
                b = [0.2465, 1.7100, 6.4094,18.6113,50.2523]
            case(7) ! N
                a = [0.1022, 0.3219, 0.7982, 0.8197, 0.1715]
                b = [0.2451, 1.7481, 6.1925,17.3894,48.1431]
            case(8) ! O
                a = [0.0974, 0.2921, 0.6910,  0.6990, 0.2039]
                b = [0.2067, 1.3815, 4.6943, 12.7105,32.4726]
                ! case('P')
                !     Z = 15; r = 1.07
                ! case('S')
                !     Z = 16; r = 1.05
                ! case('PD')
                !     Z = 46; r = 1.2
            case(26) ! Fe
                a = [0.3946, 1.2725, 1.7031,  2.3140,  1.4795]
                b = [0.2717, 2.0443, 7.6007, 29.9714, 86.2265]
            case(78) ! Pt
                a = [0.9148, 1.8096, 3.2134,  3.2953,  1.5754]
                b = [0.2263, 1.3813, 5.3243, 17.5987, 60.0171]
            case(79) ! Au
                a = [0.9674, 1.8916, 3.3993,  3.0524,  1.2607]
                b = [0.2358, 1.4712, 5.6758, 18.7119, 61.5286]
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
                        if(r2 > cutoffsq) cycle
                        rmat(j,k,l) = rmat(j,k,l) + epot(r2,aterm,b)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        call vol%set_rmat(rmat)
        deallocate(rmat)
    contains
        ! potential assuming static atoms (eq B.6)
        real function epot(r2,aterm,b)
            real, intent(in) :: r2, aterm(5),b(5)
            epot = C * sum( aterm * exp(-fourpisq*r2/b) )
        end function epot

    end subroutine convolve

    subroutine geometry_analysis_pdb(self, pdbfile)
      use simple_math, only : plane_from_points
      class(atoms),     intent(inout) :: self
      character(len=*), intent(in)    :: pdbfile   ! all the atomic positions
      character(len=2)     :: element
      type(atoms)          :: init_atoms, final_atoms
      real,    allocatable :: radii(:),line(:,:), plane(:,:,:),points(:,:), distances_totheplane(:), distances_totheline(:)
      real,    allocatable :: w(:),v(:,:),d(:),pointsTrans(:,:)
      logical, allocatable :: flag(:) ! flags the atoms belonging to the plane/column
      integer, parameter   :: N_DISCRET = 500
      integer :: i, j, n, n_tot, t, s, filnum, io_stat, cnt_intersect, cnt, n_cc
      real    :: atom1(3), atom2(3), atom3(3), dir_1(3), dir_2(3), vec(3), m(3), dist_plane, dist_line
      real    :: t_vec(N_DISCRET), s_vec(N_DISCRET), denominator, centroid(3), prod(3)
      call init_atoms%new(pdbfile)
      n = init_atoms%get_n()
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
              if(dist_line <= 0.6*self%radius(i)) then ! it intersects atoms
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
       print *, 'Directional vector of the line', d
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
                if(dist_plane <= 0.6*self%radius(i)) then ! it intersects atoms i
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
            endif
        enddo
        call final_atoms%writePDB('AtomPlane')
        call final_atoms%kill
        stop
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
        distances_totheplane = (distances_totheplane)
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

!@descr: two PDB entries the tests use (1JYX beta-galactosidase, 6VXX SARS-CoV-2 spike), read from testdata/pdb
! The atoms live in testdata/pdb/1jyx.pdb and 6vxx.pdb of the source tree (until 2026-11-02 they were array
! constructors here, 56 000 lines compiled into the library). The files are read by fixed PDB columns, so
! every field comes back exactly as written.
module simple_molecule_data
use simple_defs,  only: SIMPLE_SOURCE_PATH
use simple_error, only: simple_exception
implicit none
private
public :: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
#include "simple_local_flags.inc"

character(len=*), parameter :: PDB_DIR = SIMPLE_SOURCE_PATH//'/testdata/pdb/'

type :: molecule_data
    integer :: n = 0
    character(len=4), allocatable :: name(:)
    character(len=3), allocatable :: resname(:)
    character(len=1), allocatable :: chain(:)
    character(len=2), allocatable :: element(:)
    integer,          allocatable :: num(:)
    integer,          allocatable :: resnum(:)
    real,             allocatable :: xyz(:,:)   ! (n,3)
    real,             allocatable :: occupancy(:)
    real,             allocatable :: beta(:)
    logical,          allocatable :: het(:)
end type molecule_data

contains

    function betagal_1jyx() result( mol )
        type(molecule_data) :: mol
        mol = read_molecule_pdb(PDB_DIR//'1jyx.pdb')
    end function betagal_1jyx

    function sars_cov2_spkgp_6vxx() result( mol )
        type(molecule_data) :: mol
        mol = read_molecule_pdb(PDB_DIR//'6vxx.pdb')
    end function sars_cov2_spkgp_6vxx

    !> every ATOM and HETATM record of fname, by the fixed columns of the wwPDB coordinate section
    function read_molecule_pdb( fname ) result( mol )
        character(len=*), intent(in) :: fname
        type(molecule_data) :: mol
        character(len=80)   :: line
        integer :: funit, ios, n, i
        open(newunit=funit, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cannot open the test molecule '//fname)
        n = 0
        do
            read(funit, '(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( is_atom(line) ) n = n + 1
        end do
        if( n == 0 ) THROW_HARD('no ATOM or HETATM records in '//fname)
        mol%n = n
        allocate(mol%name(n), mol%resname(n), mol%chain(n), mol%element(n), mol%num(n), mol%resnum(n), &
            &mol%xyz(n,3), mol%occupancy(n), mol%beta(n), mol%het(n))
        rewind(funit)
        i = 0
        do while( i < n )
            read(funit, '(A)', iostat=ios) line
            if( ios /= 0 ) THROW_HARD('cannot read the records of '//fname)
            if( .not. is_atom(line) ) cycle
            i = i + 1
            mol%het(i)     = line(1:6) == 'HETATM'
            mol%name(i)    = line(13:16)
            mol%resname(i) = line(18:20)
            mol%chain(i)   = line(22:22)
            mol%element(i) = adjustl(line(77:78))
            read(line, '(6X,I5,11X,I4,4X,3F8.3,2F6.2)', iostat=ios) mol%num(i), mol%resnum(i), mol%xyz(i,1:3), &
                &mol%occupancy(i), mol%beta(i)
            if( ios /= 0 ) THROW_HARD('malformed record '//trim(line)//' in '//fname)
        end do
        close(funit)
    end function read_molecule_pdb

    logical function is_atom( line )
        character(len=*), intent(in) :: line
        is_atom = line(1:6) == 'ATOM  ' .or. line(1:6) == 'HETATM'
    end function is_atom

end module simple_molecule_data

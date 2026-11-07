!@descr: the atoms of a PDB entry as plain arrays (molecule_data), the input of atoms%new_from_molecule and pdb2mrc
module simple_molecule_data
implicit none
private
public :: molecule_data

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

end module simple_molecule_data

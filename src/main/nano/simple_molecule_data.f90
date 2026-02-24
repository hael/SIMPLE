!@descr: example of molecule data used for simple testing
module simple_molecule_data
implicit none
private 
public :: molecule_data, betagal_1jyx_molecule

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

    function betagal_1jyx_molecule() result( mol )
        type(molecule_data) :: mol
        integer, parameter  :: n = 20
        mol%n = n
        allocate(mol%name(n), mol%resname(n), mol%chain(n), mol%element(n))
        allocate(mol%num(n),  mol%resnum(n))
        allocate(mol%xyz(n,3))
        allocate(mol%occupancy(n), mol%beta(n))
        allocate(mol%het(n))
        ! First 20 ATOM records from 1JYX.pdb
        mol%name    = [character(len=4) :: &
            ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NE ', ' CZ ', ' NH1', &
            ' NH2', ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NE ', ' CZ '  ]
        mol%resname = [character(len=3) :: &
            'ARG','ARG','ARG','ARG','ARG','ARG','ARG','ARG','ARG','ARG', &
            'ARG','ARG','ARG','ARG','ARG','ARG','ARG','ARG','ARG','ARG' ]
        mol%chain   = [character(len=1) :: &
            'A','A','A','A','A','A','A','A','A','A', &
            'A','A','A','A','A','A','A','A','A','A' ]
        mol%element = [character(len=2) :: &
            'N ','C ','C ','O ','C ','C ','C ','N ','C ','N ', &
            'N ','N ','C ','C ','O ','C ','C ','C ','N ','C ' ]
        mol%num     = [ &
            1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20 ]
        mol%resnum  = [ &
            13,13,13,13,13,13,13,13,13,13, &
            13,14,14,14,14,14,14,14,14,14 ]
        mol%occupancy = [ &
            1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00, &
            1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00 ]
        mol%beta = [ &
            40.72,25.65,30.42,40.19,41.05,69.42,50.40,38.70,37.48,30.35, &
            39.88,28.98,24.45,26.58,27.63,34.79,32.88,37.80,40.72,36.85 ]
        mol%het = [ &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false. ]
        ! Coordinates: xyz(:,1)=x, xyz(:,2)=y, xyz(:,3)=z  (PDB Å)
        mol%xyz = reshape( [ &
            ! x (1..20)
            -4.688, -4.288, -4.636, -4.021, -2.826, -2.483, -2.993, -2.605, -1.789, -1.196, &
            -1.583, -5.644, -5.938, -5.519, -6.332, -6.642, -7.722, -9.144, -9.052, -10.042, &
            ! y (1..20)
            -57.052, -55.651, -54.762, -53.679, -55.643, -55.404, -55.413, -56.491, -56.754, -56.091, &
            -57.686, -54.708, -53.360, -52.928, -51.802, -53.515, -54.580, -54.486, -55.681, -55.613, &
            ! z (1..20)
            -5.121, -5.262, -3.874, -3.784, -5.684, -7.155, -7.674, -7.042, -7.580, -8.718, &
            -6.821, -2.588, -2.323, -3.379, -3.314, -1.015, -0.754, -1.824, -1.835, -10.937  &
        ], [n,3] )
    end function betagal_1jyx_molecule

end module simple_molecule_data





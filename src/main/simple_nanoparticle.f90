module simple_nanoparticle
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'

use simple_image,      only: image
use simple_binimage,   only: binimage
use simple_atoms,      only: atoms
use simple_parameters, only: params_glob
!use simple_qr_solve
use simple_nanoparticle_utils
implicit none

public :: nanoparticle
private
#include "simple_local_flags.inc"

! module global constants
real,             parameter :: CORR_THRES_SIGMA    = -2.0    ! sigma for valid_corr thresholding
integer,          parameter :: NBIN_THRESH         = 15      ! number of thresholds for binarization
integer,          parameter :: CN_THRESH_XTAL      = 5       ! cn-threshold highly crystalline NPs
integer,          parameter :: NVOX_THRESH         = 3       ! min # voxels per atom is 3
logical,          parameter :: DEBUG               = .false. ! for debugging purposes
logical,          parameter :: WRITE_OUTPUT        = .false. ! for figures generation
logical,          parameter :: ATOMS_STATS_OMIT    = .false. ! omit = shorter atoms stats output
integer,          parameter :: SOFT_EDGE           = 6
integer,          parameter :: N_DISCRET           = 1000
integer,          parameter :: CNMIN               = 3
integer,          parameter :: CNMAX               = 13
integer,          parameter :: NSTRAIN_COMPS       = 7
character(len=*), parameter :: ATOMS_STATS_FILE    = 'atoms_stats.csv'
character(len=*), parameter :: NP_STATS_FILE       = 'nanoparticle_stats.csv'
character(len=*), parameter :: CN_STATS_FILE       = 'cn_dependent_stats.csv'
character(len=*), parameter :: ATOM_VAR_CORRS_FILE = 'atom_param_corrs.txt'

character(len=*), parameter :: ATOM_STATS_HEAD = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//&
&'CN_STD'//CSV_DELIM//'NN_BONDL'//CSV_DELIM//'CN_GEN'//CSV_DELIM//'DIAM'//CSV_DELIM//&
&'ADJ_CN13'//CSV_DELIM//'AVG_INT'//CSV_DELIM//'MAX_INT'//CSV_DELIM//'CENDIST'//CSV_DELIM//'VALID_CORR'&
&//CSV_DELIM//'U_ISO'//CSV_DELIM//'U_MAJ'//CSV_DELIM//'U_MED'//CSV_DELIM//'U_MIN'//CSV_DELIM//'AZIMUTH'//&
&CSV_DELIM//'POLAR'//CSV_DELIM//'DOI'//CSV_DELIM//'DOI_MIN'//CSV_DELIM//'ISO_CORR'//CSV_DELIM//'ANISO_CORR'//CSV_DELIM//'X'//CSV_DELIM&
&//'Y'//CSV_DELIM//'Z'//CSV_DELIM//'EXX_STRAIN'//CSV_DELIM//'EYY_STRAIN'//CSV_DELIM//'EZZ_STRAIN'//&
&CSV_DELIM//'EXY_STRAIN'//CSV_DELIM//'EYZ_STRAIN'//CSV_DELIM//'EXZ_STRAIN'//CSV_DELIM//'RADIAL_STRAIN'

character(len=*), parameter :: ATOM_STATS_HEAD_OMIT = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//&
&'CN_STD'//CSV_DELIM//'NN_BONDL'//CSV_DELIM//'CN_GEN'//CSV_DELIM//'DIAM'//CSV_DELIM//'ADJ_CN13'//&
&CSV_DELIM//'AVG_INT'//CSV_DELIM//'MAX_INT'//CSV_DELIM//'CENDIST'//CSV_DELIM//'VALID_CORR'//CSV_DELIM//&
&'U_ISO'//CSV_DELIM//'U_MAJ'//CSV_DELIM//'U_MED'//CSV_DELIM//&
&'U_MIN'//CSV_DELIM//'AZIMUTH'//CSV_DELIM//'POLAR'//CSV_DELIM//'DOI'//CSV_DELIM//'DOI_MIN'//CSV_DELIM//'ISO_CORR'&
&//CSV_DELIM//'ANISO_CORR'//CSV_DELIM//'RADIAL_STRAIN'

character(len=*), parameter :: NP_STATS_HEAD = 'NATOMS'//CSV_DELIM//'NANISO'//CSV_DELIM//'DIAM'//&
&CSV_DELIM//'AVG_NVOX'//CSV_DELIM//'MED_NVOX'//CSV_DELIM//'SDEV_NVOX'//&
&CSV_DELIM//'AVG_CN_STD'//CSV_DELIM//'MED_CN_STD'//CSV_DELIM//'SDEV_CN_STD'//&
&CSV_DELIM//'AVG_NN_BONDL'//CSV_DELIM//'MED_NN_BONDL'//CSV_DELIM//'SDEV_NN_BONDL'//&
&CSV_DELIM//'AVG_CN_GEN'//CSV_DELIM//'MED_CN_GEN'//CSV_DELIM//'SDEV_CN_GEN'//&
&CSV_DELIM//'AVG_DIAM'//CSV_DELIM//'MED_DIAM'//CSV_DELIM//'SDEV_DIAM'//&
&CSV_DELIM//'AVG_AVG_INT'//CSV_DELIM//'MED_AVG_INT'//CSV_DELIM//'SDEV_AVG_INT'//&
&CSV_DELIM//'AVG_MAX_INT'//CSV_DELIM//'MED_MAX_INT'//CSV_DELIM//'SDEV_MAX_INT'//&
&CSV_DELIM//'AVG_VALID_CORR'//CSV_DELIM//'MED_VALID_CORR'//CSV_DELIM//'SDEV_VALID_CORR'//&
&CSV_DELIM//'AVG_U_ISO'//CSV_DELIM//'MED_U_ISO'//CSV_DELIM//'SDEV_U_ISO'//&
&CSV_DELIM//'AVG_U_MAJ'//CSV_DELIM//'MED_U_MAJ'//CSV_DELIM//'SDEV_U_MAJ'//&
&CSV_DELIM//'AVG_U_MED'//CSV_DELIM//'MED_U_MED'//CSV_DELIM//'SDEV_U_MED'//&
&CSV_DELIM//'AVG_U_MIN'//CSV_DELIM//'MED_U_MIN'//CSV_DELIM//'SDEV_U_MIN'//&
&CSV_DELIM//'AVG_AZIMUTH'//CSV_DELIM//'MED_AZIMUTH'//CSV_DELIM//'SDEV_AZIMUTH'//&
&CSV_DELIM//'AVG_POLAR'//CSV_DELIM//'MED_POLAR'//CSV_DELIM//'SDEV_POLAR'//&
&CSV_DELIM//'AVG_DOI'//CSV_DELIM//'MED_DOI'//CSV_DELIM//'SDEV_DOI'//&
&CSV_DELIM//'AVG_DOI_MIN'//CSV_DELIM//'MED_DOI_MIN'//CSV_DELIM//'SDEV_DOI_MIN'//&
&CSV_DELIM//'AVG_ISO_CORR'//CSV_DELIM//'MED_ISO_CORR'//CSV_DELIM//'SDEV_ISO_CORR'//&
&CSV_DELIM//'AVG_ANISO_CORR'//CSV_DELIM//'MED_ANISO_CORR'//CSV_DELIM//'SDEV_ANISO_CORR'//&
&CSV_DELIM//'AVG_RADIAL_STRAIN'//CSV_DELIM//'MED_RADIAL_STRAIN'//CSV_DELIM//'SDEV_RADIAL_STRAIN'//&
&CSV_DELIM//'MIN_RADIAL_STRAIN'//CSV_DELIM//'MAX_RADIAL_STRAIN'

character(len=*), parameter :: CN_STATS_HEAD = 'CN_STD'//CSV_DELIM//'NATOMS'//CSV_DELIM//'NANISO'//&
&CSV_DELIM//'AVG_NVOX'//CSV_DELIM//'MED_NVOX'//CSV_DELIM//'SDEV_NVOX'//&
&CSV_DELIM//'AVG_NN_BONDL'//CSV_DELIM//'MED_NN_BONDL'//CSV_DELIM//'SDEV_NN_BONDL'//&
&CSV_DELIM//'AVG_CN_GEN'//CSV_DELIM//'MED_CN_GEN'//CSV_DELIM//'SDEV_CN_GEN'//&
&CSV_DELIM//'AVG_DIAM'//CSV_DELIM//'MED_DIAM'//CSV_DELIM//'SDEV_DIAM'//&
&CSV_DELIM//'AVG_AVG_INT'//CSV_DELIM//'MED_AVG_INT'//CSV_DELIM//'SDEV_AVG_INT'//&
&CSV_DELIM//'AVG_MAX_INT'//CSV_DELIM//'MED_MAX_INT'//CSV_DELIM//'SDEV_MAX_INT'//&
&CSV_DELIM//'AVG_VALID_CORR'//CSV_DELIM//'MED_VALID_CORR'//CSV_DELIM//'SDEV_VALID_CORR'//&
&CSV_DELIM//'AVG_U_ISO'//CSV_DELIM//'MED_U_ISO'//CSV_DELIM//'SDEV_U_ISO'//&
&CSV_DELIM//'AVG_U_MAJ'//CSV_DELIM//'MED_U_MAJ'//CSV_DELIM//'SDEV_U_MAJ'//&
&CSV_DELIM//'AVG_U_MED'//CSV_DELIM//'MED_U_MED'//CSV_DELIM//'SDEV_U_MED'//&
&CSV_DELIM//'AVG_U_MIN'//CSV_DELIM//'MED_U_MIN'//CSV_DELIM//'SDEV_U_MIN'//&
&CSV_DELIM//'AVG_AZIMUTH'//CSV_DELIM//'MED_AZIMUTH'//CSV_DELIM//'SDEV_AZIMUTH'//&
&CSV_DELIM//'AVG_POLAR'//CSV_DELIM//'MED_POLAR'//CSV_DELIM//'SDEV_POLAR'//&
&CSV_DELIM//'AVG_DOI'//CSV_DELIM//'MED_DOI'//CSV_DELIM//'SDEV_DOI'//&
&CSV_DELIM//'AVG_DOI_MIN'//CSV_DELIM//'MED_DOI_MIN'//CSV_DELIM//'SDEV_DOI_MIN'//&
&CSV_DELIM//'AVG_ISO_CORR'//CSV_DELIM//'MED_ISO_CORR'//CSV_DELIM//'SDEV_ISO_CORR'//&
&CSV_DELIM//'AVG_ANISO_CORR'//CSV_DELIM//'MED_ANISO_CORR'//CSV_DELIM//'SDEV_ANISO_CORR'//&
&CSV_DELIM//'AVG_RADIAL_STRAIN'//CSV_DELIM//'MED_RADIAL_STRAIN'//CSV_DELIM//'SDEV_RADIAL_STRAIN'//&
&CSV_DELIM//'MIN_RADIAL_STRAIN'//CSV_DELIM//'MAX_RADIAL_STRAIN'

! container for per-atom statistics
type :: atom_stats
    ! various per-atom parameters                                                                   ! csv file labels
    integer :: cc_ind            = 0  ! index of the connected component                            INDEX
    integer :: size              = 0  ! number of voxels in connected component                     NVOX
    integer :: cn_std            = 0  ! standard coordination number                                CN_STD
    integer :: adjacent_cn13     = 0  ! number of neighboring atoms w/ cn > 12                      ADJ_CN13
    real    :: bondl             = 0. ! nearest neighbour bond lenght in A                          NN_BONDL
    real    :: cn_gen            = 0. ! generalized coordination number                             CN_GEN
    real    :: diam              = 0. ! atom diameter                                               DIAM
    real    :: avg_int           = 0. ! average grey level intensity across the connected component AVG_INT
    real    :: max_int           = 0. ! maximum            -"-                                      MAX_INT
    real    :: cendist           = 0. ! distance from the centre of mass of the nanoparticle        CENDIST
    real    :: valid_corr        = 0. ! per-atom correlation with the simulated map                 VALID_CORR
    real    :: u_iso             = 0. ! isotropic displacement parameters (b-factor)                U_ISO
    real    :: u_evals(3)        = 0. ! eigenvalues of aniso displacement matrix (maj,med,min)      U_(MAJ,MED,MIN)
    real    :: azimuth           = 0. ! Azimuthal angle of major eigenvector [0,Pi)                 AZIMUTH
    real    :: polar             = 0. ! Polar angle of major eigenvector [0, Pi)                    POLAR
    real    :: doi               = 0. ! degree of isotropy                                          DOI
    real    :: doi_min           = 0. ! degree of isotropy w/r to average eigenvalues of aniso      DOI_MIN
    real    :: isocorr           = 0. ! Correlation of isotropic B-Factor fit to input map          ISO_CORR
    real    :: anisocorr         = 0. ! Correlation of anisotropic B-Factor fit to input map        ANISO_CORR
    real    :: center(3)         = 0. ! atom center                                                 X Y Z
    
    ! strain
    real    :: exx_strain        = 0. ! tensile strain in %                                         EXX_STRAIN
    real    :: eyy_strain        = 0. ! -"-                                                         EYY_STRAIN
    real    :: ezz_strain        = 0. ! -"-                                                         EZZ_STRAIN
    real    :: exy_strain        = 0. ! -"-                                                         EXY_STRAIN
    real    :: eyz_strain        = 0. ! -"-                                                         EYZ_STRAIN
    real    :: exz_strain        = 0. ! -"-                                                         EXZ_STRAIN
    real    :: radial_strain     = 0. ! -"-                                                         RADIAL_STRAIN

    ! Auxillary (non-output)
    real    :: aniso(3,3)        = 0. ! ADP Matrix for ANISOU PDB file                              N/A
    logical :: tossADP           = .false. ! True if atom inadequate for ADP calculations           N/A

end type atom_stats

type :: nanoparticle
    private
    type(image)           :: img, img_raw
    type(binimage)        :: img_bin, img_cc         ! binary and connected component images
    integer               :: ldim(3)            = 0  ! logical dimension of image
    integer               :: n_cc               = 0  ! number of atoms (connected components)                NATOMS
    integer               :: n_aniso            = 0  ! number of atoms with aniso calculations               NANISO
    integer               :: n4stats            = 0  ! number of atoms in subset used for stats calc
    real                  :: smpd               = 0. ! sampling distance
    real                  :: NPcen(3)           = 0. ! coordinates of the center of mass of the nanoparticle
    real                  :: NPdiam             = 0. ! diameter of the nanoparticle                          DIAM
    real                  :: theoretical_radius = 0. ! theoretical atom radius in A
    ! GLOBAL NP STATS
    type(stats_struct)    :: map_stats
    ! -- the rest
    type(stats_struct)    :: size_stats
    type(stats_struct)    :: cn_std_stats
    type(stats_struct)    :: bondl_stats
    type(stats_struct)    :: cn_gen_stats
    type(stats_struct)    :: diam_stats
    type(stats_struct)    :: avg_int_stats
    type(stats_struct)    :: max_int_stats
    type(stats_struct)    :: valid_corr_stats
    type(stats_struct)    :: u_iso_stats
    type(stats_struct)    :: u_maj_stats
    type(stats_struct)    :: u_med_stats
    type(stats_struct)    :: u_min_stats
    type(stats_struct)    :: azimuth_stats
    type(stats_struct)    :: polar_stats
    type(stats_struct)    :: isocorr_stats
    type(stats_struct)    :: anisocorr_stats
    type(stats_struct)    :: doi_stats
    type(stats_struct)    :: doi_min_stats
    type(stats_struct)    :: radial_strain_stats
    ! CN-DEPENDENT STATS
    ! -- # atoms
    real                  :: natoms_cns(CNMIN:CNMAX) = 0.       ! # of atoms per cn_std                            NATOMS
    real                  :: natoms_aniso_cns(CNMIN:CNMAX) = 0. ! # of atoms with aniso calculated per cn_std      NATOMS
    ! -- the rest
    type(stats_struct)    :: size_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: bondl_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: cn_gen_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: diam_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: avg_int_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: max_int_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: valid_corr_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: u_iso_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: u_maj_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: u_med_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: u_min_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: azimuth_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: polar_stats_cns(CNMIN:CNMAX)  
    type(stats_struct)    :: doi_stats_cns(CNMIN:CNMAX) 
    type(stats_struct)    :: doi_min_stats_cns(CNMIN:CNMAX) 
    type(stats_struct)    :: isocorr_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: anisocorr_stats_cns(CNMIN:CNMAX) 
    type(stats_struct)    :: radial_strain_stats_cns(CNMIN:CNMAX)
    ! PER-ATOM STATISTICS
    type(atom_stats), allocatable :: atominfo(:)
    real,             allocatable :: coords4stats(:,:)
    ! OTHER
    character(len=2)      :: element   = ' '
    character(len=4)      :: atom_name = '    '
    character(len=STDLEN) :: npname    = '' ! fname
    character(len=STDLEN) :: fbody     = '' ! fbody
  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_ldim
    procedure          :: get_natoms
    procedure          :: get_valid_corrs
    procedure          :: set_img
    procedure          :: set_atomic_coords
    procedure          :: set_coords4stats
    procedure, private :: pack_instance4stats
    ! utils
    procedure, private :: atominfo2centers
    procedure, private :: atominfo2centers_A
    procedure, private :: center_on_atom
    procedure          :: update_ncc
    ! atomic position determination
    procedure          :: identify_lattice_params
    procedure          :: identify_atomic_pos
    procedure, private :: binarize_and_find_centers
    procedure, private :: find_centers
    procedure, private :: discard_atoms_with_low_contact_score
    procedure, private :: discard_lowly_coordinated
    procedure, private :: discard_low_valid_corr_atoms
    procedure, private :: split_atoms
    procedure          :: validate_atoms
    ! calc stats
    procedure          :: fillin_atominfo
    procedure, private :: masscen
    procedure, private :: calc_longest_atm_dist
    procedure, private :: check_neighbors_cn
    procedure, private :: calc_isotropic_disp
    procedure, private :: calc_anisotropic_disp
    ! visualization and output
    procedure          :: simulate_atoms
    procedure, private :: write_centers_1
    procedure, private :: write_centers_2
    generic            :: write_centers => write_centers_1, write_centers_2
    procedure          :: write_centers_aniso
    procedure          :: write_individual_atoms
    procedure          :: write_csv_files
    procedure, private :: write_atominfo
    procedure, private :: write_np_stats
    procedure, private :: write_cn_stats
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    subroutine new_nanoparticle( self, fname, msk )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        real, optional,      intent(in)    :: msk
        character(len=2) :: el_ucase
        integer :: nptcls
        integer :: Z ! atomic number
        real    :: smpd
        call self%kill
        self%npname    = fname
        self%fbody     = get_fbody(trim(basename(fname)), trim(fname2ext(fname)))
        self%smpd      = params_glob%smpd
        self%atom_name = ' '//params_glob%element
        self%element   = params_glob%element
        el_ucase       = upperCase(trim(adjustl(params_glob%element)))
        call get_element_Z_and_radius(el_ucase, Z, self%theoretical_radius)
        if( Z == 0 ) THROW_HARD('Unknown element: '//el_ucase)
        call find_ldim_nptcls(self%npname, self%ldim, nptcls, smpd)
        call self%img%new(self%ldim, self%smpd)
        call self%img_bin%new_bimg(self%ldim, self%smpd)
        call self%img_bin%new(self%ldim, self%smpd)
        call self%img%read(fname)
        if( present(msk) ) call self%img%mask(msk, 'soft')
        call self%img_raw%copy(self%img)
        call self%img_raw%stats(self%map_stats%avg, self%map_stats%sdev, self%map_stats%maxv, self%map_stats%minv)
    end subroutine new_nanoparticle

    ! getters/setters

    subroutine get_ldim(self,ldim)
        class(nanoparticle), intent(in)  :: self
        integer,             intent(out) :: ldim(3)
        ldim = self%img%get_ldim()
    end subroutine get_ldim

    function get_natoms(self) result(n)
       class(nanoparticle), intent(inout)  :: self
       integer :: n
       call self%img_cc%get_nccs(n)
    end function get_natoms

    function get_valid_corrs( self ) result( corrs )
        class(nanoparticle), intent(in) :: self
        real, allocatable :: corrs(:)
        if( allocated(self%atominfo) )then
            allocate(corrs(size(self%atominfo)), source=self%atominfo(:)%valid_corr)
        endif
    end function get_valid_corrs

    ! set one of the images of the nanoparticle type
    subroutine set_img( self, imgfile, which )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: imgfile
        character(len=*),    intent(in)    :: which
        select case(which)
            case('img')
                call self%img%new(self%ldim, self%smpd)
                call self%img%read(imgfile)
            case('img_bin')
                call self%img_bin%new_bimg(self%ldim, self%smpd)
                call self%img_bin%read_bimg(imgfile)
            case('img_cc')
                call self%img_cc%new_bimg(self%ldim, self%smpd)
                call self%img_cc%read_bimg(imgfile)
            case('img_raw')
                call self%img%new(self%ldim, self%smpd)
                call self%img%read(imgfile)
            case DEFAULT
                THROW_HARD('Wrong input parameter img type (which); set_img')
        end select
    end subroutine set_img

    ! sets the atom positions to be the ones in the inputted PDB file.
    subroutine set_atomic_coords( self, pdb_file )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: pdb_file
        type(atoms) :: a
        integer     :: i, N
        if( fname2ext(pdb_file) .ne. 'pdb' ) THROW_HARD('Inputted filename has to have pdb extension; set_atomic_coords')
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        call a%new(pdb_file)
        N = a%get_n() ! number of atoms
        allocate(self%atominfo(N))
        do i = 1, N
            self%atominfo(i)%center(:) = a%get_coord(i)/self%smpd + 1.
        enddo
        self%n_cc = N
        call a%kill
    end subroutine set_atomic_coords

    subroutine set_coords4stats( self, pdb_file )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: pdb_file
        call read_pdb2matrix(pdb_file, self%coords4stats)
        self%n4stats = size(self%coords4stats, dim=2)
    end subroutine set_coords4stats

    subroutine pack_instance4stats( self, strain_array )
        class(nanoparticle), intent(inout) :: self
        real, allocatable,   intent(inout) :: strain_array(:,:)
        real,                allocatable   :: centers_A(:,:), strain_array_new(:,:)
        logical,             allocatable   :: mask(:)
        integer,             allocatable   :: imat_cc(:,:,:), imat_cc_new(:,:,:), imat_bin_new(:,:,:)
        type(atom_stats),    allocatable   :: atominfo_new(:)
        integer :: n_cc_orig, cc, cnt, nx, ny, nz
        if( .not. allocated(self%coords4stats) ) return
        centers_A = self%atominfo2centers_A()
        n_cc_orig = size(centers_A, dim=2)
        allocate(mask(n_cc_orig), source=.false.)
        call find_atoms_subset(self%coords4stats, centers_A, mask)
        ! remove atoms not in mask
        if( n_cc_orig /= self%n_cc ) THROW_HARD('incongruent # cc:s')
        ! (1) update img_cc & img_bin
        call self%img_cc%get_imat(imat_cc)
        nx = size(imat_cc, dim=1)
        ny = size(imat_cc, dim=2)
        nz = size(imat_cc, dim=3)
        allocate(imat_cc_new(nx,ny,nz), imat_bin_new(nx,ny,nz), source=0)
        cnt = 0
        do cc = 1, self%n_cc
            if( mask(cc) )then
                cnt = cnt + 1
                where(imat_cc == cc) imat_cc_new  = cnt
                where(imat_cc == cc) imat_bin_new = 1
            endif
        end do
        call self%img_cc%set_imat(imat_cc_new)
        call self%img_bin%set_imat(imat_bin_new)
        ! (2) update atominfo & strain_array
        allocate(atominfo_new(cnt), strain_array_new(cnt,NSTRAIN_COMPS))
        cnt = 0
        do cc = 1, self%n_cc
            if( mask(cc) )then
                cnt = cnt + 1
                atominfo_new(cnt)       = self%atominfo(cc)
                strain_array_new(cnt,:) = strain_array(cc,:)
            endif
        end do
        deallocate(self%atominfo, strain_array)
        allocate(self%atominfo(cnt), source=atominfo_new)
        allocate(strain_array(cnt,NSTRAIN_COMPS), source=strain_array_new)
        deallocate(centers_A, mask, imat_cc, imat_cc_new, imat_bin_new, atominfo_new)
        ! (3) update number of connected components
        self%n_cc = cnt
    end subroutine pack_instance4stats

    ! utils

    function atominfo2centers( self, mask ) result( centers )
        class(nanoparticle), intent(in) :: self
        logical, optional,   intent(in) :: mask(size(self%atominfo))
        real, allocatable :: centers(:,:)
        logical :: mask_present
        integer :: sz, i, cnt
        sz = size(self%atominfo)
        mask_present = .false.
        if( present(mask) ) mask_present = .true.
        if( mask_present )then
            cnt = count(mask)
            allocate(centers(3,cnt), source=0.)
            cnt = 0
            do i = 1, sz
                if( mask(i) )then
                    cnt = cnt + 1
                    centers(:,cnt) = self%atominfo(i)%center(:)
                endif
            end do
        else
            allocate(centers(3,sz), source=0.)
            do i = 1, sz
                centers(:,i) = self%atominfo(i)%center(:)
            end do
        endif
    end function atominfo2centers

    function atominfo2centers_A( self, mask ) result( centers_A )
        class(nanoparticle), intent(in) :: self
        logical, optional,   intent(in) :: mask(size(self%atominfo))
        real, allocatable :: centers_A(:,:)
        logical :: mask_present
        integer :: sz, i, cnt
        sz = size(self%atominfo)
        mask_present = .false.
        if( present(mask) ) mask_present = .true.
        if( mask_present )then
            cnt = count(mask)
            allocate(centers_A(3,cnt), source=0.)
            cnt = 0
            do i = 1, sz
                if( mask(i) )then
                    cnt = cnt + 1
                    centers_A(:,cnt) = (self%atominfo(i)%center(:) - 1.) * self%smpd
                endif
            end do
        else
            allocate(centers_A(3,sz), source=0.)
            do i = 1, sz
                centers_A(:,i) = (self%atominfo(i)%center(:) - 1.) * self%smpd
            end do
        endif
    end function atominfo2centers_A

    ! Translate the identified atomic positions so that the center of mass
    ! of the nanoparticle coincides with its closest atom
    subroutine center_on_atom( self, pdbfile_in, pdbfile_out )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: pdbfile_in
        character(len=*),    intent(inout) :: pdbfile_out
        type(atoms) :: atom_centers
        real        :: m(3), vec(3), d, d_before
        integer     :: i
        call atom_centers%new(pdbfile_in)
        m(:)     = self%masscen()
        d_before = huge(d_before)
        do i = 1, self%n_cc
            d = euclid(m,self%atominfo(i)%center(:))
            if(d < d_before) then
                vec(:)   = m(:) - self%atominfo(i)%center(:)
                d_before = d
            endif
        enddo
        do i = 1, self%n_cc
            self%atominfo(i)%center(:) = self%atominfo(i)%center(:) + vec
            call atom_centers%set_coord(i,(self%atominfo(i)%center(:)-1.)*self%smpd)
        enddo
        call atom_centers%writePDB(pdbfile_out)
        call atom_centers%kill
    end subroutine center_on_atom

    subroutine update_ncc( self, img_cc )
        class(nanoparticle),      intent(inout) :: self
        type(binimage), optional, intent(inout) :: img_cc
        if( present(img_cc) )then
            call img_cc%get_nccs(self%n_cc)
        else
            call self%img_cc%get_nccs(self%n_cc)
        endif
    end subroutine update_ncc

    ! atomic position determination

    subroutine identify_lattice_params( self, a, use_auto_corr_thres )
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: a(3)                ! lattice parameters
        logical,             intent(in)    :: use_auto_corr_thres ! true -> use automatic corr thres
        real, allocatable :: centers_A(:,:)                       ! coordinates of the atoms in ANGSTROMS
        type(image)       :: simatms
        integer           :: n_discard
        ! MODEL BUILDING
        ! Phase correlation approach
        call phasecorr_one_atom(self%img, self%img, self%element)
        ! Nanoparticle binarization
        call self%binarize_and_find_centers()
        ! atom splitting by correlation map validation
        call self%split_atoms()
        ! OUTLIERS DISCARDING
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms)
        call self%validate_atoms(simatms)
        ! discard atoms with low valid_corr
        call self%discard_low_valid_corr_atoms(use_auto_corr_thres, n_discard)
        ! fit lattice
        centers_A = self%atominfo2centers_A()
        call fit_lattice(self%element, centers_A, a)
        deallocate(centers_A)
        call simatms%kill
    end subroutine identify_lattice_params

    subroutine identify_atomic_pos( self, a, l_fit_lattice, use_cs_thres, use_auto_corr_thres, cs_thres, split_fname )
        class(nanoparticle),        intent(inout) :: self
        real,                       intent(inout) :: a(3)                ! lattice parameters
        logical,                    intent(in)    :: l_fit_lattice       ! fit lattice or use inputted
        logical,                    intent(inout) :: use_cs_thres        ! use or not contact score thres
        logical,                    intent(in)    :: use_auto_corr_thres ! true -> use automatic corr thres
        integer,          optional, intent(in)    :: cs_thres
        character(len=*), optional, intent(in)    :: split_fname
        logical     :: use_cn_thresh, fixed_cs_thres
        type(image) :: simatms, img_cos
        integer     :: n_discard
        ! MODEL BUILDING
        ! Phase correlation approach
        call phasecorr_one_atom(self%img, self%img, self%element)
        if( WRITE_OUTPUT ) call self%img%write('denoised.mrc')
        ! Nanoparticle binarization
        call self%binarize_and_find_centers()
        ! atom splitting by correlation map validation
        call self%split_atoms(split_fname)
        ! OUTLIERS DISCARDING
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms)
        if (WRITE_OUTPUT) call simatms%write(trim(self%fbody)//'_SIM_pre_validation.mrc')
        call self%validate_atoms(simatms)
        if ( WRITE_OUTPUT )call self%write_centers('valid_corr_in_bfac_field_pre_discard', 'valid_corr')
        ! discard atoms with low valid_corr
        call self%discard_low_valid_corr_atoms(use_auto_corr_thres, n_discard)
        ! discard lowly coordinated atoms
        fixed_cs_thres = present(cs_thres)
        if( fixed_cs_thres )then
            call self%discard_atoms_with_low_contact_score(use_cn_thresh, cs_thres)
        else if( use_cs_thres )then
            call self%discard_atoms_with_low_contact_score(use_cn_thresh)
            if( use_cn_thresh ) call self%discard_lowly_coordinated(CN_THRESH_XTAL, a, l_fit_lattice)
        endif
        ! re-calculate valid_corr:s (since they are otherwise lost from the B-factor field due to reallocations of atominfo)
        call self%simulate_atoms(simatms)
        call self%validate_atoms(simatms)
        ! WRITE OUTPUT
        call self%img_bin%write_bimg(trim(self%fbody)//'_BIN.mrc')
        write(logfhandle,'(A)') 'output, binarized map:            '//trim(self%fbody)//'_BIN.mrc'
        call self%img_bin%grow_bins(1)
        call self%img_bin%cos_edge(SOFT_EDGE, img_cos)
        call img_cos%write(trim(self%fbody)//'_MSK.mrc')
        write(logfhandle,'(A)') 'output, envelope mask map:        '//trim(self%fbody)//'_MSK.mrc'
        call self%img_cc%write_bimg(trim(self%fbody)//'_CC.mrc')
        write(logfhandle,'(A)') 'output, connected components map: '//trim(self%fbody)//'_CC.mrc'
        call self%write_centers
        call simatms%write(trim(self%fbody)//'_SIM.mrc')
        write(logfhandle,'(A)') 'output, simulated atomic density: '//trim(self%fbody)//'_SIM.mrc'
        ! destruct
        call img_cos%kill
        call simatms%kill
    end subroutine identify_atomic_pos

    ! This subrotuine takes in input a nanoparticle and
    ! binarizes it by thresholding. The gray level histogram is split
    ! in 20 parts, which corrispond to 20 possible threshold.
    ! Among those threshold, the selected one is the for which
    ! tha correlation between the raw map and a simulated distribution
    ! obtained with that threshold reaches the maximum value.
    subroutine binarize_and_find_centers( self )
        class(nanoparticle), intent(inout) :: self
        type(binimage)       :: img_bin_t
        type(binimage)       :: img_ccs_t
        type(atoms)          :: atom
        type(image)          :: simulated_distrib
        integer, allocatable :: imat_t(:,:,:)
        real,    allocatable :: x_mat(:)  ! vectorization of the volume
        real,    allocatable :: coords(:,:)
        real,    allocatable :: rmat(:,:,:)
        integer :: i, fnr
        real    :: otsu_thresh, corr, step, step_refine, max_corr, thresh, thresh_opt, lbt, rbt
        logical, parameter      :: L_BENCH = .false.
        real(timer_int_kind)    :: rt_find_ccs, rt_find_centers, rt_gen_sim, rt_real_corr, rt_tot
        integer(timer_int_kind) ::  t_find_ccs,  t_find_centers,  t_gen_sim,  t_real_corr,  t_tot
        call otsu_nano(self%img,otsu_thresh) ! find initial threshold
        write(logfhandle,'(A)') '>>> BINARIZATION'
        rmat = self%img%get_rmat()
        x_mat = pack(rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                    &rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) >= 0.)
        allocate(imat_t(self%ldim(1), self%ldim(2), self%ldim(3)), source = 0)
        step = (maxval(x_mat)-otsu_thresh )/real(NBIN_THRESH)
        step_refine = step / 6.
        deallocate(x_mat)
        call simulated_distrib%new(self%ldim,self%smpd)
        rt_find_ccs     =  0.
        rt_find_centers =  0.
        rt_gen_sim      =  0.
        rt_real_corr    =  0.
        t_tot           =  tic()
        max_corr        = -1.
        ! discrete search
        do i = 1, NBIN_THRESH
            if( i == 1 )then
                thresh = otsu_thresh
            else
                thresh = thresh + step
            endif
            corr = t2c( thresh )
            write(logfhandle,*) 'threshold: ', thresh , 'corr: ', corr
            if( corr > max_corr )then
                max_corr  = corr
                thresh_opt = thresh
            endif
            if( corr < max_corr ) exit ! convex goal function
        enddo
        ! refinement
        lbt    = thresh_opt - step + step_refine
        rbt    = thresh_opt + step - step_refine
        thresh = lbt
        i      = NBIN_THRESH
        do while( thresh <= rbt )
            i = i + 1
            corr = t2c( thresh )
            write(logfhandle,*) 'threshold: ', thresh, 'corr: ', corr
            if( corr > max_corr )then
                max_corr  = corr
                thresh_opt = thresh
            endif
            thresh = thresh + step_refine
        end do
        rt_tot = toc(t_tot)
        if( L_BENCH )then
            call fopen(fnr, FILE='BINARIZE_AND_FIND_CENTERS_BENCH.txt', STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'find_ccs       : ', rt_find_ccs
            write(fnr,'(a,1x,f9.2)') 'find_centers   : ', rt_find_centers
            write(fnr,'(a,1x,f9.2)') 'gen_sim        : ', rt_gen_sim
            write(fnr,'(a,1x,f9.2)') 'real_corr      : ', rt_real_corr
            write(fnr,'(a,1x,f9.2)') 'total time     : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'find_ccs       : ', (rt_find_ccs/rt_tot)     * 100.
            write(fnr,'(a,1x,f9.2)') 'find_centers   : ', (rt_find_centers/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'gen_sim        : ', (rt_gen_sim/rt_tot)      * 100.
            write(fnr,'(a,1x,f9.2)') 'real_corr      : ', (rt_real_corr/rt_tot)    * 100.
            write(fnr,'(a,1x,f9.2)') 'total time     : ', rt_tot
            write(fnr,'(a,1x,f9.2)') '% accounted for: ',&
            &((rt_find_ccs+rt_find_centers+rt_gen_sim+rt_real_corr)/rt_tot)     * 100.
            call fclose(fnr)
        endif
        write(logfhandle,*) 'optimal threshold: ', thresh_opt, 'max_corr: ', max_corr
        ! Update img_bin and img_cc
        corr = t2c( thresh_opt )
        call self%img_bin%copy_bimg(img_bin_t)
        call self%img_cc%copy_bimg(img_ccs_t)
        call self%update_ncc()
        call self%find_centers()
        call img_bin_t%kill_bimg
        call img_ccs_t%kill_bimg
        ! deallocate and kill
        if(allocated(rmat))   deallocate(rmat)
        if(allocated(imat_t)) deallocate(imat_t)
        if(allocated(coords)) deallocate(coords)
        call simulated_distrib%kill
        write(logfhandle,'(A)') '>>> BINARIZATION, COMPLETED'

    contains

        ! Otsu binarization for nanoparticle maps
        ! It considers the grey level value only in the positive range.
        ! It doesn't threshold the map. It just returns the ideal threshold.
        ! This is based on the implementation of 1D otsu
        subroutine otsu_nano(img, scaled_thresh)
            type(image),    intent(inout) :: img
            real,           intent(out)   :: scaled_thresh ! returns the threshold in the correct range
            real, pointer     :: rmat(:,:,:)
            real, allocatable :: x(:)
            call img%get_rmat_ptr(rmat)
            x = pack(rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                    &rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.)
            call otsu(size(x), x, scaled_thresh)
        end subroutine otsu_nano

        real function t2c( thres )
            real, intent(in) :: thres
            where(rmat > thres)
                imat_t = 1
            elsewhere
                imat_t = 0
            endwhere
            ! Generate binary image and cc image
            call img_bin_t%new_bimg(self%ldim, self%smpd)
            call img_bin_t%set_imat(imat_t)
            t_find_ccs = tic()
            call img_ccs_t%new_bimg(self%ldim, self%smpd)
            call img_bin_t%find_ccs(img_ccs_t)
            rt_find_ccs = rt_find_ccs + toc(t_find_ccs)
            ! Find atom centers in the generated distributions
            call self%update_ncc(img_ccs_t) ! self%n_cc is needed in find_centers
            t_find_centers = tic()
            call self%find_centers(img_bin_t, img_ccs_t, coords)
            rt_find_centers = rt_find_centers + toc(t_find_centers)
            ! Generate a simulated distribution based on those center
            t_gen_sim = tic()
            call self%write_centers('centers_'//trim(int2str(i))//'_iteration', coords)
            call atom%new          ('centers_'//trim(int2str(i))//'_iteration.pdb')
            call atom%convolve(simulated_distrib, cutoff = 8.*self%smpd)
            call del_file('centers_'//trim(int2str(i))//'_iteration.pdb')
            call atom%kill
            rt_gen_sim = rt_gen_sim + toc(t_gen_sim)
            ! correlate volumes
            t_real_corr = tic()
            t2c = self%img%real_corr(simulated_distrib)
            rt_real_corr = rt_real_corr + toc(t_real_corr)
            if( WRITE_OUTPUT ) call simulated_distrib%write('simvol_thres'//trim(real2str(thres))//'_corr'//trim(real2str(t2c))//'.mrc')
        end function t2c

    end subroutine binarize_and_find_centers

    ! Find the centers coordinates of the atoms in the particle
    ! and save it in the global variable centers.
    ! If coords is present, it saves it also in coords.
    subroutine find_centers( self, img_bin, img_cc, coords )
        class(nanoparticle),         intent(inout) :: self
        type(binimage), optional,    intent(inout) :: img_bin, img_cc
        real, optional, allocatable, intent(out)   :: coords(:,:)
        real,        pointer :: rmat_raw(:,:,:)
        integer, allocatable :: imat_cc_in(:,:,:)
        logical, allocatable :: mask(:,:,:)
        integer :: i, ii, jj, kk
        real    :: m(3), sum_mass
        ! sanity check
        if( present(img_bin) .and. .not. present(img_cc) ) THROW_HARD('img_bin and img_cc have to be both present')
        ! global variables allocation
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        allocate( self%atominfo(self%n_cc) )
        if( present(img_cc) )then
            call img_cc%get_imat(imat_cc_in)
        else
            call self%img_cc%get_imat(imat_cc_in)
        endif
        call self%img_raw%get_rmat_ptr(rmat_raw)
        allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
        !$omp parallel do default(shared) private(i,ii,jj,kk,mask,m,sum_mass) schedule(static) proc_bind(close)
        do i=1,self%n_cc
            mask     = .true.
            where( imat_cc_in /= i ) mask = .false.
            m        = 0.
            sum_mass = 0.
            do ii = 1, self%ldim(1)
                do jj = 1, self%ldim(2)
                    do kk = 1, self%ldim(3)
                        if( mask(ii,jj,kk) )then
                            m = m + real([ii,jj,kk]) * rmat_raw(ii,jj,kk)
                            sum_mass = sum_mass + rmat_raw(ii,jj,kk)
                        endif
                    enddo
                enddo
            enddo
            self%atominfo(i)%center(:) = m / sum_mass
        enddo
        !$omp end parallel do
        ! saving centers coordinates, optional
        if( present(coords) )then
            allocate(coords(3,self%n_cc))
            do i=1,self%n_cc
                coords(:,i) = self%atominfo(i)%center(:)
            end do
        endif
    end subroutine find_centers

    subroutine split_atoms( self, fname )
        class(nanoparticle),        intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        type(binimage)       :: img_split_ccs
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:)
        integer, allocatable :: imat(:,:,:), imat_cc(:,:,:), imat_bin(:,:,:), imat_split_ccs(:,:,:)
        integer, parameter   :: RANK_THRESH = 4
        integer :: icc, cnt, cnt_split
        integer :: rank, m(1)
        real    :: new_centers(3,2*self%n_cc) ! will pack it afterwards if it has too many elements
        real    :: pc, radius
        write(logfhandle, '(A)') '>>> SPLITTING CONNECTED ATOMS'
        call self%img%get_rmat_ptr(rmat_pc) ! rmat_pc contains the phase correlation
        call self%img_cc%get_imat(imat_cc)  ! to pass to the subroutine split_atoms
        allocate(imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),           source = imat_cc)
        allocate(imat_split_ccs(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = 0)
        call img_split_ccs%new_bimg(self%ldim, self%smpd)
        call img_split_ccs%new(self%ldim, self%smpd)
        cnt       = 0
        cnt_split = 0
        do icc = 1, self%n_cc ! for each cc check if the center corresponds with the local max of the phase corr
            pc = rmat_pc(nint(self%atominfo(icc)%center(1)),nint(self%atominfo(icc)%center(2)),nint(self%atominfo(icc)%center(3)))
            ! calculate the rank
            x = pack(rmat_pc(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=imat == icc)
            call hpsort(x)
            m(:) = minloc(abs(x - pc))
            rank = size(x) - m(1)
            deallocate(x)
            ! calculate radius
            call self%calc_longest_atm_dist(icc, radius)
            ! split
            if( rank > RANK_THRESH .or. radius > 1.5 * self%theoretical_radius )then
                where(imat == icc)
                    imat_split_ccs = 1
                end where
                cnt_split = cnt_split + 1
                call split_atom(new_centers,cnt)
            else
                cnt = cnt + 1 ! new number of centers derived from splitting
                new_centers(:,cnt) = self%atominfo(icc)%center(:)
            endif
        enddo
        write(logfhandle,*) '# atoms split: ', cnt_split
        deallocate(self%atominfo)
        self%n_cc = cnt ! update
        allocate(self%atominfo(cnt))
        ! update centers
        do icc = 1, cnt
            self%atominfo(icc)%center(:) = new_centers(:,icc)
        enddo
        call self%img_bin%get_imat(imat_bin)
        ! update binary image
        where( imat_cc > 0 )
            imat_bin = 1
        elsewhere
            imat_bin = 0
        endwhere
        ! update relevant data fields
        call img_split_ccs%set_imat(imat_split_ccs)
        if( present(fname) )then
            call img_split_ccs%write(trim(fname))
        else
            call img_split_ccs%write('split_ccs.mrc')
        endif
        call img_split_ccs%kill
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%update_img_rmat()
        call self%img_bin%find_ccs(self%img_cc)
        call self%update_ncc(self%img_cc)
        call self%find_centers()
        write(logfhandle, '(A)') '>>> SPLITTING CONNECTED ATOMS, COMPLETED'

    contains

        subroutine split_atom(new_centers,cnt)
            real,    intent(inout) :: new_centers(:,:) ! updated coordinates of the centers
            integer, intent(inout) :: cnt              ! atom counter, to update the center coords
            integer :: new_center1(3), new_center2(3), new_center3(3)
            integer :: i, j, k
            logical :: found3d_cen
            logical :: mask(self%ldim(1),self%ldim(2),self%ldim(3)) ! false in the layer of connection of the atom to be split
            mask = .false. ! initialization
            ! Identify first new center
            new_center1 = maxloc(rmat_pc(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=imat == icc)
            cnt = cnt + 1
            new_centers(:,cnt) = real(new_center1)
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if( imat(i,j,k) == icc )then
                            if(((real(i - new_center1(1)))**2 + (real(j - new_center1(2)))**2 + &
                            &   (real(k - new_center1(3)))**2) * self%smpd  <=  (0.9 * self%theoretical_radius)**2) then
                                mask(i,j,k) = .true.
                            endif
                        endif
                    enddo
                enddo
            enddo
            ! Second likely center.
            new_center2 = maxloc(rmat_pc(:self%ldim(1),:self%ldim(2),:self%ldim(3)), (imat == icc) .and. .not. mask)
            if( any(new_center2 > 0) )then ! if anything was found
                ! Validate second center (check if it's 2 merged atoms, or one pointy one)
                if( sum(real(new_center2 - new_center1)**2.) * self%smpd <= (0.9 * self%theoretical_radius)**2) then
                    ! the new_center2 is within the diameter of the atom position at new_center1
                    ! therefore, it is not another atom and should be removed
                    where( imat_cc == icc .and. (.not.mask) ) imat_cc = 0
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center2)
                    ! In the case of two merged atoms, build the second atom
                    do i = 1, self%ldim(1)
                        do j = 1, self%ldim(2)
                            do k = 1, self%ldim(3)
                                if( imat(i,j,k) == icc )then
                                    if(((real(i - new_center2(1)))**2 + (real(j - new_center2(2)))**2 +&
                                    &   (real(k - new_center2(3)))**2) * self%smpd <= (0.9 * self%theoretical_radius)**2 )then
                                        mask(i,j,k) = .true.
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif
            ! Third likely center.
            new_center3 = maxloc(rmat_pc(:self%ldim(1),:self%ldim(2),:self%ldim(3)), (imat == icc) .and. .not. mask)
            if( any(new_center3 > 0) )then ! if anything was found
                ! Validate third center
                if(sum(real(new_center3 - new_center1)**2.) * self%smpd <= (0.9 * self%theoretical_radius)**2 .or. &
                &  sum(real(new_center3 - new_center2)**2.) * self%smpd <= (0.9 * self%theoretical_radius)**2 )then
                    ! the new_center3 is within the diameter of the atom position at new_center1 or new_center2
                    ! therefore, it is not another atom and should be removed
                    where( imat_cc == icc .and. (.not.mask) ) imat_cc = 0
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center3)
                    found3d_cen = .false.
                    ! In the case of two merged atoms, build the second atom
                    do i = 1, self%ldim(1)
                        do j = 1, self%ldim(2)
                            do k = 1, self%ldim(3)
                                if( imat(i,j,k) == icc )then
                                    if( ((real(i - new_center3(1)))**2 + (real(j - new_center3(2)))**2 + &
                                    &    (real(k - new_center3(3)))**2) * self%smpd <= (0.9 * self%theoretical_radius)**2 )then
                                         found3d_cen = .not.mask(i,j,k)
                                         mask(i,j,k) = .true.
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif
            ! Set the merged cc back to 0
            where(imat_cc == icc .and. (.not.mask) ) imat_cc = 0
            call self%img_cc%set_imat(imat_cc)
        end subroutine split_atom

    end subroutine split_atoms

    subroutine validate_atoms( self, simatms )
        class(nanoparticle), intent(inout) :: self
        class(image),        intent(in)    :: simatms
        real, allocatable :: centers(:,:)           ! coordinates of the atoms in PIXELS
        real, allocatable :: pixels1(:), pixels2(:) ! pixels extracted around the center
        real    :: maxrad
        integer :: ijk(3), npix_in, npix_out1, npix_out2, i, winsz
        type(stats_struct) :: corr_stats
        maxrad  = (self%theoretical_radius * 1.5) / self%smpd ! in pixels
        winsz   = ceiling(maxrad)
        npix_in = (2 * winsz + 1)**3 ! cubic window size (-winsz:winsz in each dim)
        centers = self%atominfo2centers()
        allocate(pixels1(npix_in), pixels2(npix_in), source=0.)
        ! calculate per-atom correlations
        do i = 1, self%n_cc
            ijk = nint(centers(:,i))
            call self%img_raw%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, npix_out1, pixels1)
            call simatms%win2arr_rad(     ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, npix_out2, pixels2)
            self%atominfo(i)%valid_corr = pearsn_serial(pixels1(:npix_out1),pixels2(:npix_out2))
        end do
        call calc_stats(self%atominfo(:)%valid_corr, corr_stats)
        write(logfhandle,'(A)') '>>> VALID_CORR (PER-ATOM CORRELATION WITH SIMULATED DENSITY) STATS BELOW'
        write(logfhandle,'(A,F8.4)') 'Average: ', corr_stats%avg
        write(logfhandle,'(A,F8.4)') 'Median : ', corr_stats%med
        write(logfhandle,'(A,F8.4)') 'Sigma  : ', corr_stats%sdev
        write(logfhandle,'(A,F8.4)') 'Max    : ', corr_stats%maxv
        write(logfhandle,'(A,F8.4)') 'Min    : ', corr_stats%minv
    end subroutine validate_atoms

    subroutine discard_low_valid_corr_atoms( self, use_auto_corr_thres, n_discard )
        use simple_stat, only: robust_sigma_thres
        class(nanoparticle), intent(inout) :: self
        logical,             intent(in)    :: use_auto_corr_thres
        integer,             intent(inout) :: n_discard
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        integer :: cc
        real    :: corr_thres
        write(logfhandle, '(A)') '>>> DISCARDING ATOMS WITH LOW VALID_CORR'
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        if( use_auto_corr_thres )then
            corr_thres = min(max(robust_sigma_thres(self%atominfo(:)%valid_corr, CORR_THRES_SIGMA), 0.3), 0.7)
            write(logfhandle, *) 'Valid_corr threshold calculated: ', corr_thres
        else
            corr_thres = params_glob%corr_thres
        endif
        write(logfhandle, *) 'Valid_corr threshold applied:    ', corr_thres
        n_discard = 0
        do cc = 1, self%n_cc
            if( self%atominfo(cc)%valid_corr < corr_thres )then
                where(imat_cc == cc) imat_bin = 0
                n_discard = n_discard + 1
            endif
        end do
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        deallocate(imat_bin, imat_cc)
        write(logfhandle, *) 'Numbers of atoms discarded because of low valid_corr ', n_discard
        write(logfhandle, *) 'Total number of atoms after discarding atoms with low valid_corr ', self%n_cc
        write(logfhandle, '(A)') '>>> DISCARDING ATOMS WITH LOW VALID_CORR, COMPLETED'
    end subroutine discard_low_valid_corr_atoms

    subroutine discard_atoms_with_low_contact_score( self, use_cn_thresh, cs_thres )
        class(nanoparticle), intent(inout) :: self
        logical,             intent(inout) :: use_cn_thresh
        integer, optional,   intent(in)    :: cs_thres
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        real, allocatable    :: centers_A(:,:) ! coordinates of the atoms in ANGSTROMS
        integer :: cscores(self%n_cc), cc, n_discard, cthresh, new_cthresh
        type(stats_struct)   :: cscore_stats
        centers_A = self%atominfo2centers_A()
        call calc_contact_scores(self%element,centers_A,cscores)
        call calc_stats(real(cscores), cscore_stats)
        write(logfhandle,'(A)') '>>> CONTACT SCORE STATS BELOW'
        write(logfhandle,'(A,F8.4)') 'Average: ', cscore_stats%avg
        write(logfhandle,'(A,F8.4)') 'Median : ', cscore_stats%med
        write(logfhandle,'(A,F8.4)') 'Sigma  : ', cscore_stats%sdev
        write(logfhandle,'(A,F8.4)') 'Max    : ', cscore_stats%maxv
        write(logfhandle,'(A,F8.4)') 'Min    : ', cscore_stats%minv
        if( .not. present(cs_thres) )then
            cthresh = min(5,max(3,nint(cscore_stats%avg - cscore_stats%sdev)))
            write(logfhandle,'(A,I3)') 'CONTACT SCORE THRESHOLD: ', cthresh
            use_cn_thresh = .false.
            if( cthresh == 5 )then ! highly crystalline
                use_cn_thresh = .true.
                return
            endif
        else
            cthresh = cs_thres
        endif
        write(logfhandle, '(A)') '>>> DISCARDING OUTLIERS BASED ON CONTACT SCORE'
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        ! Removing outliers from the binary image and the connected components image
        ! remove atoms with < NVOX_THRESH voxels
        n_discard = 0
        do cc = 1, self%n_cc
            if( count(imat_cc == cc) < NVOX_THRESH )then
                where(imat_cc == cc) imat_bin = 0
                n_discard = n_discard + 1
            endif
        end do
        ! Removing outliers based on coordination number
        call remove_lowly_contacted( cthresh )
        ! don't leave behind any atoms with cscore < cthresh - 2
        new_cthresh = cthresh - 2
        if( new_cthresh < 1 )then
            ! we're done
        else
            call remove_lowly_contacted( new_cthresh )
        endif
        deallocate(imat_bin, imat_cc, centers_A)
        write(logfhandle, *) 'Numbers of atoms discarded because of low cscore ', n_discard
        write(logfhandle, *) 'Total number of atoms after discarding outliers based on cscore  ', self%n_cc
        write(logfhandle, '(A)') '>>> DISCARDING OUTLIERS BASED ON CONTACT SCORE, COMPLETED'

        contains

            subroutine remove_lowly_contacted( cthresh )
                integer, intent(in) :: cthresh
                do cc = 1, self%n_cc
                    if( cscores(cc) < cthresh )then
                        where(imat_cc == cc) imat_bin = 0
                        n_discard = n_discard + 1
                    endif
                enddo
                call self%img_bin%set_imat(imat_bin)
                call self%img_bin%find_ccs(self%img_cc)
                ! update number of connected components
                call self%img_cc%get_nccs(self%n_cc)
                call self%find_centers()
                if( allocated(centers_A) ) deallocate(centers_A)
                centers_A = self%atominfo2centers_A()
                call calc_contact_scores(self%element,centers_A,cscores)
            end subroutine remove_lowly_contacted

    end subroutine discard_atoms_with_low_contact_score

    ! This subroutine discards outliers that resisted binarization
    ! It calculates the standard coordination number (cn) of each atom and discards
    ! the atoms with cn_std < cn_thresh
    ! It modifies the img_bin and img_cc instances deleting the identified outliers.
    subroutine discard_lowly_coordinated( self, cn_thresh, a, l_fit_lattice )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: cn_thresh     ! threshold for discarding outliers based on coordination number
        real,                intent(inout) :: a(3)          ! lattice parameter
        logical,             intent(in)    :: l_fit_lattice ! fit lattice or use inputted
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        real, allocatable    :: centers_A(:,:) ! coordinates of the atoms in ANGSTROMS
        real    :: cn_gen(self%n_cc)
        integer :: cn(self%n_cc), cc, n_discard, new_cn_thresh
        write(logfhandle, '(A)') '>>> DISCARDING OUTLIERS BASED ON CN'
        centers_A = self%atominfo2centers_A()
        if( l_fit_lattice ) call fit_lattice(self%element, centers_A, a) ! else use inputted lattice params
        call run_cn_analysis(self%element,centers_A,a,cn,cn_gen)
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        ! Removing outliers from the binary image and the connected components image
        ! remove atoms with < NVOX_THRESH voxels
        n_discard = 0
        do cc = 1, self%n_cc
            if( count(imat_cc == cc) < NVOX_THRESH )then
                where(imat_cc == cc) imat_bin = 0
                n_discard = n_discard + 1
            endif
        end do
        ! Removing outliers based on coordination number
        n_discard = 0
        call remove_lowly_coordinated( cn_thresh )
        ! don't leave behind any atoms with cn_std < cn_thresh - 2
        new_cn_thresh = cn_thresh - 2
        if( new_cn_thresh < 1 )then
            ! we're done
        else
            call remove_lowly_coordinated( new_cn_thresh )
        endif
        deallocate(imat_bin, imat_cc, centers_A)
        write(logfhandle, *) 'Numbers of atoms discarded because of low cn ', n_discard
        write(logfhandle, *) 'Total number of atoms after discarding outliers based on cn      ', self%n_cc
        write(logfhandle, '(A)') '>>> DISCARDING OUTLIERS BASED ON CN, COMPLETED'

        contains

            subroutine remove_lowly_coordinated( cn_thresh )
                integer, intent(in) :: cn_thresh
                do cc = 1, self%n_cc
                    if( cn(cc) < cn_thresh )then
                        where(imat_cc == cc) imat_bin = 0
                        n_discard = n_discard + 1
                    endif
                enddo
                call self%img_bin%set_imat(imat_bin)
                call self%img_bin%find_ccs(self%img_cc)
                ! update number of connected components
                call self%img_cc%get_nccs(self%n_cc)
                call self%find_centers()
                if( allocated(centers_A) ) deallocate(centers_A)
                centers_A = self%atominfo2centers_A()
                call run_cn_analysis(self%element,centers_A,a,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
            end subroutine remove_lowly_coordinated

    end subroutine discard_lowly_coordinated

    ! calc stats
    subroutine fillin_atominfo( self, a0 )
        class(nanoparticle),        intent(inout) :: self
        real,             optional, intent(in)    :: a0(3) ! lattice parameters
        type(image)          :: simatms, fit_isotropic, fit_anisotropic
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:), strain_array(:,:)
        real,    pointer     :: rmat_raw(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        logical, allocatable :: cc_mask(:)
        real                 :: tmp_diam, a(3)
        integer              :: i, cc, cn, max_size
        character(*), parameter :: fn_fit_isotropic="fit_isotropic.mrc", fn_fit_anisotropic="fit_anisotropic.mrc"

        write(logfhandle, '(A)') '>>> EXTRACTING ATOM STATISTICS'
        write(logfhandle, '(A)') '---Dev Note: ADP and Max Neighboring Displacements Under Testing---'
        ! calc cn and cn_gen
        centers_A = self%atominfo2centers_A()
        if( present(a0) )then
            a = a0
        else
            call fit_lattice(self%element, centers_A, a)
        endif
        call run_cn_analysis(self%element,centers_A,a,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        ! calc strain and lattice displacements for all atoms
        allocate(strain_array(self%n_cc,NSTRAIN_COMPS), source=0.)
        call strain_analysis(self%element, centers_A, a, strain_array)
        if( allocated(self%coords4stats) ) call self%pack_instance4stats(strain_array)
        allocate(cc_mask(self%n_cc), source=.true.) ! because self%n_cc might change after pack_instance4stats
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms)
        call self%validate_atoms(simatms)
        ! calc NPdiam & NPcen
        tmpcens     = self%atominfo2centers()
        self%NPdiam = 0.
        do i = 1, self%n_cc
            tmp_diam = pixels_dist(self%atominfo(i)%center(:), tmpcens, 'max', cc_mask)
            if( tmp_diam > self%NPdiam ) self%NPdiam = tmp_diam
        enddo
        cc_mask     = .true. ! restore
        self%NPdiam = self%NPdiam * self%smpd ! in A
        write(logfhandle,*) 'nanoparticle diameter (A): ', self%NPdiam
        self%NPcen  = self%masscen()
        write(logfhandle,*) 'nanoparticle mass center: ', self%NPcen
        ! CALCULATE PER-ATOM PARAMETERS
        ! extract atominfo
        allocate(mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = .false.)
        call self%img_raw%get_rmat_ptr(rmat_raw)
        call self%img_cc%get_imat(imat_cc)

        call fit_isotropic%new(self%img_raw%get_ldim(), self%img_raw%get_smpd())
        call fit_anisotropic%new(self%img_raw%get_ldim(), self%img_raw%get_smpd())
        max_size = 0
        do cc = 1, self%n_cc
            call progress(cc, self%n_cc)
            ! index of the connected component
            self%atominfo(cc)%cc_ind = cc
            ! number of voxels in connected component
            where( imat_cc == cc ) mask = .true.
            self%atominfo(cc)%size = count(mask)
            ! distance from the centre of mass of the nanoparticle
            self%atominfo(cc)%cendist = euclid(self%atominfo(cc)%center(:), self%NPcen) * self%smpd
            ! atom diameter
            call self%calc_longest_atm_dist(cc, self%atominfo(cc)%diam)
            self%atominfo(cc)%diam = 2.*self%atominfo(cc)%diam ! radius --> diameter in A
            ! whether atom is neighbors with an atom of CN > 12
            if (maxval(self%atominfo(:)%cn_std) > 12) then
                call self%check_neighbors_cn(cc, a)
            end if
            ! maximum grey level intensity across the connected component
            self%atominfo(cc)%max_int = maxval(rmat_raw(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            ! average grey level intensity across the connected component
            self%atominfo(cc)%avg_int = sum(rmat_raw(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            self%atominfo(cc)%avg_int = self%atominfo(cc)%avg_int / real(count(mask))
            ! bond length of nearest neighbour...
            self%atominfo(cc)%bondl = pixels_dist(self%atominfo(cc)%center(:), tmpcens, 'min', mask=cc_mask) ! Use all the atoms
            self%atominfo(cc)%bondl = self%atominfo(cc)%bondl * self%smpd ! convert to A
            ! Atomic displacement
            call self%calc_isotropic_disp(cc, a, rmat_raw, fit_isotropic)
            call self%calc_anisotropic_disp(cc, a, rmat_raw, fit_anisotropic)
            ! set strain values
            self%atominfo(cc)%exx_strain    = strain_array(cc,1)
            self%atominfo(cc)%eyy_strain    = strain_array(cc,2)
            self%atominfo(cc)%ezz_strain    = strain_array(cc,3)
            self%atominfo(cc)%exy_strain    = strain_array(cc,4)
            self%atominfo(cc)%eyz_strain    = strain_array(cc,5)
            self%atominfo(cc)%exz_strain    = strain_array(cc,6)
            self%atominfo(cc)%radial_strain = strain_array(cc,7)
            ! reset masks
            mask    = .false.
            cc_mask = .true.
        end do
        write(logfhandle,'(a,i5)') "ADP Tossed: ", count(self%atominfo(:)%tossADP)
        write(logfhandle, '(A)') '>>> WRITING OUTPUT'
        self%n_aniso = self%n_cc - count(self%atominfo(:)%tossADP)
        call fit_isotropic%write(fn_fit_isotropic)
        call fit_anisotropic%write(fn_fit_anisotropic)
        if (WRITE_OUTPUT) call write_2D_slice()

        ! CALCULATE GLOBAL NP PARAMETERS
        call calc_stats(  real(self%atominfo(:)%size),    self%size_stats,         mask=self%atominfo(:)%size >= NVOX_THRESH )
        call calc_stats(  real(self%atominfo(:)%cn_std),  self%cn_std_stats        )
        call calc_stats(  self%atominfo(:)%bondl,         self%bondl_stats         )
        call calc_stats(  self%atominfo(:)%cn_gen,        self%cn_gen_stats        )
        call calc_stats(  self%atominfo(:)%diam,          self%diam_stats,         mask=self%atominfo(:)%size >= NVOX_THRESH )
        call calc_zscore( self%atominfo(:)%avg_int ) ! to get comparable intensities between different particles
        call calc_zscore( self%atominfo(:)%max_int ) ! -"-
        call calc_stats(  self%atominfo(:)%avg_int,       self%avg_int_stats       )
        call calc_stats(  self%atominfo(:)%max_int,       self%max_int_stats       )
        call calc_stats(  self%atominfo(:)%valid_corr,    self%valid_corr_stats    )
        call calc_stats(  self%atominfo(:)%u_iso,         self%u_iso_stats         )
        call calc_stats(  self%atominfo(:)%u_evals(1),    self%u_maj_stats,         mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%u_evals(2),    self%u_med_stats,         mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%u_evals(3),    self%u_min_stats,         mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%azimuth,       self%azimuth_stats,       mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%polar,         self%polar_stats,         mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%doi,           self%doi_stats,           mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%doi_min,       self%doi_min_stats,       mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%isocorr,       self%isocorr_stats       )
        call calc_stats(  self%atominfo(:)%anisocorr,     self%anisocorr_stats,     mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%radial_strain, self%radial_strain_stats )
        ! CALCULATE CN-DEPENDENT STATS & WRITE CN-ATOMS
        do cn = CNMIN, CNMAX
            call calc_cn_stats( cn )
            call write_cn_atoms( cn )
        end do
        ! write pdf files with valid_corr and max_int in the B-factor field (for validation/visualisation)
        call self%write_centers('valid_corr_in_bfac_field', 'valid_corr')
        call self%write_centers('max_int_in_bfac_field',    'max_int')
        call self%write_centers('cn_std_in_bfac_field',     'cn_std')
        call self%write_centers('u_iso_in_bfac_field',      'u_iso')
        call self%write_centers('doi_in_bfac_field',        'doi')
        call self%write_centers('doi_min_in_bfac_field',    'doi_min')
        call self%write_centers_aniso('aniso_bfac_field')
        ! Write a pdb file containing the ideal lattice positions (with no missing atoms) for visualization
        ! destruct
        deallocate(mask, cc_mask, imat_cc, tmpcens, strain_array, centers_A)
        call fit_isotropic%kill
        call fit_anisotropic%kill
        call simatms%kill
        write(logfhandle, '(A)') '>>> EXTRACTING ATOM STATISTICS, COMPLETED'

        contains

            subroutine calc_zscore( arr )
                real, intent(inout) :: arr(:)
                arr = (arr - self%map_stats%avg) / self%map_stats%sdev
            end subroutine calc_zscore

            subroutine calc_cn_stats( cn )
                integer, intent(in)  :: cn ! calculate stats for given std cn
                integer :: n
                logical :: cn_mask(self%n_cc), size_mask(self%n_cc), adp_mask(self%n_cc)
                ! Generate masks
                cn_mask   = self%atominfo(:)%cn_std == cn
                size_mask = self%atominfo(:)%size >= NVOX_THRESH .and. cn_mask
                adp_mask = (.not. self%atominfo(:)%tossADP) .and. cn_mask
                n         = count(cn_mask)
                if( n == 0 ) return
                ! -- # atoms
                self%natoms_cns(cn) = real(n)
                self%natoms_aniso_cns(cn) = count(adp_mask)
                if( n < 2 ) return
                ! -- the rest
                call calc_stats( real(self%atominfo(:)%size),    self%size_stats_cns(cn),          mask=size_mask )
                call calc_stats( self%atominfo(:)%bondl,         self%bondl_stats_cns(cn),         mask=cn_mask   )
                call calc_stats( self%atominfo(:)%cn_gen,        self%cn_gen_stats_cns(cn),        mask=cn_mask   )
                call calc_stats( self%atominfo(:)%diam,          self%diam_stats_cns(cn),          mask=size_mask )
                call calc_stats( self%atominfo(:)%avg_int,       self%avg_int_stats_cns(cn),       mask=cn_mask   )
                call calc_stats( self%atominfo(:)%max_int,       self%max_int_stats_cns(cn),       mask=cn_mask   )
                call calc_stats( self%atominfo(:)%valid_corr,    self%valid_corr_stats_cns(cn),    mask=cn_mask   )
                call calc_stats( self%atominfo(:)%u_iso,         self%u_iso_stats_cns(cn),         mask=cn_mask   )
                call calc_stats( self%atominfo(:)%u_evals(1),    self%u_maj_stats_cns(cn),         mask=adp_mask  )
                call calc_stats( self%atominfo(:)%u_evals(2),    self%u_med_stats_cns(cn),         mask=adp_mask  )
                call calc_stats( self%atominfo(:)%u_evals(3),    self%u_min_stats_cns(cn),         mask=adp_mask  )
                call calc_stats( self%atominfo(:)%azimuth,       self%azimuth_stats_cns(cn),       mask=adp_mask  )
                call calc_stats( self%atominfo(:)%polar,         self%polar_stats_cns(cn),         mask=adp_mask  )
                call calc_stats( self%atominfo(:)%doi,           self%doi_stats_cns(cn),           mask=adp_mask  )
                call calc_stats( self%atominfo(:)%doi_min,       self%doi_min_stats_cns(cn),       mask=adp_mask  )
                call calc_stats( self%atominfo(:)%isocorr,       self%isocorr_stats_cns(cn),       mask=cn_mask   )
                call calc_stats( self%atominfo(:)%anisocorr,     self%anisocorr_stats_cns(cn),     mask=adp_mask  )
                call calc_stats( self%atominfo(:)%radial_strain, self%radial_strain_stats_cns(cn), mask=cn_mask   )
            end subroutine calc_cn_stats

            subroutine write_cn_atoms( cn_std )
                integer, intent(in)  :: cn_std
                type(binimage)       :: img_atom
                type(image)          :: simatms
                type(atoms)          :: atoms_obj
                integer, allocatable :: imat(:,:,:), imat_atom(:,:,:)
                logical :: cn_mask(self%n_cc)
                integer :: i
                ! make cn mask
                cn_mask = self%atominfo(:)%cn_std == cn_std
                call self%simulate_atoms(simatms, mask=cn_mask, atoms_obj=atoms_obj)
                call atoms_obj%writepdb('atoms_cn'//int2str_pad(cn_std,2))
                call simatms%write('simvol_cn'//int2str_pad(cn_std,2)//'.mrc')
                ! make binary image of atoms with given cn_std
                call img_atom%copy_bimg(self%img_cc)
                allocate(imat_atom(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
                call img_atom%get_imat(imat)
                do i = 1, self%n_cc
                    if( cn_mask(i) )then
                        where( imat == i ) imat_atom = 1
                    endif
                enddo
                call img_atom%set_imat(imat_atom)
                call img_atom%write_bimg('binvol_cn'//int2str_pad(cn_std,2)//'.mrc')
                deallocate(imat,imat_atom)
                call img_atom%kill_bimg
                call simatms%kill
                call atoms_obj%kill
            end subroutine write_cn_atoms

            ! Outputs the intensities of a 2D slice containing the center of the largest CC 
            ! as a CSV file for visualization in Python/Matlab. Output positions are the voxel
            ! positions with respect to the CC center.
            subroutine write_2D_slice()
                integer :: cc, cc_largest, max_size, center(3), i, j, funit
                character(*), parameter :: fn_slice='cc_2D_slice.csv'
                character(*), parameter :: header='X'//CSV_DELIM//'Y'//CSV_DELIM//'INTENSITY'
                ! Find largest CC for best visualization
                max_size = 0
                cc_largest = 0
                do cc=1, self%n_cc
                    if (self%atominfo(cc)%size > max_size) then
                        cc_largest = cc
                        max_size = self%atominfo(cc)%size
                    end if
                end do
                center(:) = self%atominfo(cc_largest)%center(:)
                ! Write output as CSV file with the cc center at the origin for easy analysis
                call fopen(funit, file=fn_slice, status='replace')
                write(funit, '(A)') 'X'//CSV_DELIM//'Y'//CSV_DELIM//'INTENSITY'
                601 format(F10.6,A2)
                602 format(F10.6)
                do i=1,self%ldim(1)
                    do j=1, self%ldim(2)
                        if ( imat_cc(i, j, center(3)) == cc_largest ) then
                            write(funit,601,advance='no') real(i-center(1)),    CSV_DELIM ! X
                            write(funit,601,advance='no') real(j-center(2)),    CSV_DELIM ! Y
                            write(funit,602) rmat_raw(i,j,center(3))                      ! Intensity
                        end if
                    end do
                end do
                call fclose(funit)
            end subroutine write_2D_slice

    end subroutine fillin_atominfo

    ! calc the avg of the centers coords
    function masscen( self ) result( m )
       class(nanoparticle), intent(inout) :: self
       real    :: m(3) ! mass center coords
       integer :: i
       m = 0.
       do i = 1, self%n_cc
           m = m + self%atominfo(i)%center(:)
       end do
       m = m / real(self%n_cc)
    end function masscen

    subroutine calc_longest_atm_dist( self, label, longest_dist )
       class(nanoparticle), intent(inout) :: self
       integer,             intent(in)    :: label
       real,                intent(out)   :: longest_dist
       integer, allocatable :: pos(:,:)
       integer, allocatable :: imat_cc(:,:,:)
       logical, allocatable :: mask_dist(:) ! for min and max dist calculation
       integer :: location(1)               ! location of vxls of the atom farthest from its center
       call self%img_cc%get_imat(imat_cc)
       where(imat_cc.eq.label)
           imat_cc = 1
       elsewhere
           imat_cc = 0
       endwhere
       call get_pixel_pos( imat_cc, pos ) ! pxls positions of the shell
       allocate(mask_dist(size(pos, dim=2)), source = .true.)
       if( size(pos,2) == 1 ) then ! if the connected component has size 1 (just 1 vxl)
           longest_dist  = self%smpd
           return
       else
           longest_dist  = pixels_dist(self%atominfo(label)%center(:), real(pos),'max', mask_dist, location) * self%smpd
       endif
       deallocate(imat_cc, pos, mask_dist)
    end subroutine calc_longest_atm_dist

    ! For a given cc in an NP with lattice params a, finds the number
    ! of neighboring atoms with CN > 12.
    subroutine check_neighbors_cn(self, cc, a)
        class(nanoparticle), intent(inout) :: self
        integer,          intent(in)    :: cc
        real,             intent(in)    :: a(3)
        character(len=2) :: el_ucase
        character(len=8) :: crystal_system
        real             :: a0, foo, d
        integer          :: i
        a0 = sum(a)/real(size(a)) ! alrithmetic mean of fitted lattice parameters
        ! Identify nearest neighbors
        el_ucase = uppercase(trim(adjustl(self%element)))
        call get_lattice_params(el_ucase, crystal_system, foo)
        select case(trim(adjustl(crystal_system)))
            case('rocksalt')
                d = a0 * ((1. / 2. + 1. / sqrt(2.)) / 2.)
            case('bcc')
                d = a0 * ((1. + sqrt(3.) / 2.) / 2.)
            case DEFAULT ! FCC by default
                d = a0 * ((1. + 1. / sqrt(2.)) / 2.)
        end select
        do i=1, self%n_cc
            if (i/=cc .and. self%atominfo(i)%cn_std > 12 .and.&
                    &euclid(self%atominfo(cc)%center(:3),self%atominfo(i)%center(:3))*self%smpd < d) then
                self%atominfo(cc)%adjacent_cn13 = self%atominfo(cc)%adjacent_cn13 + 1
            end if
        end do
    end subroutine check_neighbors_cn


    ! Calculates the isotropic displacement parameter U_ISO of a given cc by fitting
    ! the real space intensity distribution of the 3D reconstructed volume stored
    ! in rmat with a 3D Gaussian with isotropic variance in a sphere of radius 3*a/(8sqrt(2))
    ! about the cc center. The best fit solution is output into the 3D map fit, and
    ! the correlation between fit and rmat is stored in atominfo.
    subroutine calc_isotropic_disp(self, cc, a, rmat, fit)
        class(nanoparticle), intent(inout)  :: self
        class(image), intent(inout)         :: fit
        integer, intent(in)                 :: cc
        real, intent(in)                    :: a(3) ! Lattice params
        real, pointer, intent(in)           :: rmat(:,:,:)
        real        :: output_rad, fit_rad, center(3), maxrad, r, var, amp, int, max_int_out
        real        :: XTWX(2,2), XTWX_inv(2,2), XTWY(2,1), B(2,1)
        integer     :: ilo, ihi, jlo, jhi, klo, khi, i, j, k, errflg
        logical     :: fit_mask(self%ldim(1),self%ldim(2),self%ldim(3))

        output_rad = (sum(a)/3)/(2.*sqrt(2.))/self%smpd  ! 1/2 FCC nearest neighbor dist
        fit_rad = 0.75 * output_rad ! 0.75 prevents fitting tails of other atoms

        ! Create search window containing sphere of fit rad to speed up loops.
        center = self%atominfo(cc)%center(:)
        maxrad  = 0.5*(sum(a)/3) / self%smpd
        ilo = max(nint(center(1) - maxrad), 1)
        ihi = min(nint(center(1) + maxrad), self%ldim(1))
        jlo = max(nint(center(2) - maxrad), 1)
        jhi = min(nint(center(2) + maxrad), self%ldim(2))
        klo = max(nint(center(3) - maxrad), 1)
        khi = min(nint(center(3) + maxrad), self%ldim(3))

        ! Linear least squares to calculate the best fit params B.
        ! Note the sum of the square error of the logs in minimized instead
        ! of the sum of the square error to make the problem linear.
        ! Solution: B = ((X^T)WX)^-1((X^T)WY)
        XTWX = 0.
        XTWX_inv = 0.
        XTWY = 0.
        do k=klo, khi
            do j=jlo, jhi
                do i=ilo, ihi
                    r = euclid(1.*(/i, j, k/), 1.*center)
                    if (r < fit_rad) then
                        int = rmat(i,j,k)
                        if (int > 0) then
                            XTWX(1,1) = XTWX(1,1) + int
                            XTWX(1,2) = XTWX(1,2) + r**2 * int
                            XTWX(2,2) = XTWX(2,2) + r**4 * int
                            XTWY(1,1) = XTWY(1,1) + int * log(int)
                            XTWY(2,1) = XTWY(2,1) + r**2 * int * log(int)
                        end if
                    end if
                end do
            end do
        end do
        XTWX(2,1) = XTWX(1,2)
        call matinv(XTWX, XTWX_inv, 2, errflg)
        B = matmul(XTWX_inv, XTWY)
        amp = exp(B(1,1)) ! Best fit peak
        var = -0.5 / B(2,1) ! Best fit variance

        ! Sample fit for goodness of fit and visualization
        max_int_out = 0.
        fit_mask = .false.
        do k=klo, khi 
            do j=jlo, jhi
                do i=ilo, ihi
                    r = euclid(1.*(/i, j, k/), 1.*center) 
                    if (r < output_rad) then
                        fit_mask(i,j,k) = .true.
                        int =  amp * exp(-0.5 * r**2 / var)
                        call fit%set_rmat_at(i, j, k, int)
                        if (int > max_int_out) then
                            max_int_out = int
                        end if
                    end if
                end do
            end do
        end do
        self%atominfo(cc)%isocorr = fit%real_corr(self%img_raw, mask=fit_mask)
        self%atominfo(cc)%u_iso = var * self%smpd**2
    end subroutine calc_isotropic_disp

    ! Calculates the anisotropic displacement parameter matrix of a given cc by fitting
    ! the real space intensity distribution of the 3D reconstructed volume stored
    ! in rmat with a 3D multivariate Gaussian in a sphere of radius 3*a/(8sqrt(2))
    ! about the cc center. The best fit solution is output into the 3D map fit, and
    ! the correlation between fit and rmat is stored in atominfo. The variances
    ! along the 3 principal axes and the angular orientation of the major axis are also
    ! calculated.
    subroutine calc_anisotropic_disp(self, cc, a, rmat, fit)
        class(nanoparticle), intent(inout)  :: self
        class(image), intent(inout)         :: fit
        integer, intent(in)                 :: cc
        real, intent(in)                    :: a(3) ! Lattice params
        real, pointer, intent(in)           :: rmat(:,:,:)
        real(kind=8), allocatable           :: X(:,:), XTW(:,:), Y(:,:)
        real(kind=8)    :: XTWX(7,7), XTWX_inv(7,7), XTWY(7,1), B(7,1), cov(3,3), cov_inv(3,3)
        real(kind=8)    :: eigenvals(3), eigenvecs(3,3)
        real(kind=8)    :: majvector(3), rvec(3,1), beta(1,1), r(3)
        real    :: center(3), output_rad, fit_rad, maxrad, int, amp, max_int_out, corr
        integer :: ilo, ihi, jlo, jhi, klo, khi, i, j, k, nvoxels, n, errflg, nrot
        logical :: fit_mask(self%ldim(1),self%ldim(2),self%ldim(3))

        output_rad = (sum(a)/3)/(2.*sqrt(2.))  ! 1/2 FCC nearest neighbor dist
        fit_rad = 0.75 * output_rad ! 0.75 prevents fitting tails of other atoms

        ! Create search window containing sphere of fit rad to speed up loops.
        center = self%atominfo(cc)%center(:)
        maxrad  = 0.5*(sum(a)/3) / self%smpd
        ilo = max(nint(center(1) - maxrad), 1)
        ihi = min(nint(center(1) + maxrad), self%ldim(1))
        jlo = max(nint(center(2) - maxrad), 1)
        jhi = min(nint(center(2) + maxrad), self%ldim(2))
        klo = max(nint(center(3) - maxrad), 1)
        khi = min(nint(center(3) + maxrad), self%ldim(3))

        ! Get nvoxels within sphere of fitting
        nvoxels = 0
        do k=klo, khi
            do j=jlo, jhi
                do i=ilo, ihi
                    if (euclid(1.*(/i, j, k/), 1.*center)*self%smpd < fit_rad) then
                        nvoxels = nvoxels + 1
                    end if
                end do
            end do
        end do
        allocate(X(nvoxels, 7), XTW(7, nvoxels), Y(nvoxels,1), source = 0._dp)

        ! Linear least squares to calculate the best fit params B
        ! Conduct fit in units of Angstroms (easier on matrix operations)
        ! Note the sum of the square error of the logs is minimized instead
        ! of the sum of the square error to make the problem linear.
        ! Solution: B = ((X^T)WX)^-1((X^T)WY)
        XTWX = 0._dp
        XTWX_inv = 0._dp
        XTWY = 0._dp
        n = 1
        X(:nvoxels,1) = 1.
        do k=klo, khi
            do j=jlo, jhi
                do i=ilo, ihi
                    r = (1.*(/i, j, k/) - center) * self%smpd
                    if (norm_2(r) < fit_rad) then
                        int = rmat(i,j,k)
                        if (int > 0) then
                            X(n,2:4)  = r(1:3)**2
                            X(n,5)    = r(1)*r(2)
                            X(n,6)    = r(1)*r(3)
                            X(n,7)    = r(2)*r(3)
                            XTW(:7,n) = int*X(n,:7)
                            Y(n,1)    = log(int)
                            n = n + 1
                        end if
                    end if
                end do
            end do
        end do
        XTWX = matmul(XTW, X)
        XTWY = matmul(XTW, Y)
        deallocate(X,XTW,Y)
        call matinv(XTWX, XTWX_inv, 7, errflg)
        B = matmul(XTWX_inv, XTWY)
        amp = exp(B(1,1)) ! Best fit peak
        ! Generate covariance matrix
        cov_inv(1,1) = -2. * B(2,1)
        cov_inv(2,2) = -2. * B(3,1)
        cov_inv(3,3) = -2. * B(4,1)
        cov_inv(1,2) = -1. * B(5,1)
        cov_inv(1,3) = -1. * B(6,1)
        cov_inv(2,3) = -1. * B(7,1)
        cov_inv(2,1) = cov_inv(1,2)
        cov_inv(3,1) = cov_inv(1,3)
        cov_inv(3,2) = cov_inv(2,3)
        call matinv(cov_inv, cov, 3, errflg)
        self%atominfo(cc)%aniso = cov

        ! Sample fit for goodness of fit and visualization
        max_int_out = 0.
        fit_mask = .false.
        do k=klo, khi 
            do j=jlo, jhi   
                do i=ilo, ihi
                    rvec(:3,1) = (1.*(/i, j, k/) - center) * self%smpd
                    if (norm_2(rvec(:3,1)) < output_rad) then
                        fit_mask(i,j,k) = .true.
                        beta = matmul(matmul(transpose(rvec),cov_inv),rvec)
                        int =  amp * exp(-0.5 * beta(1,1))
                        call fit%set_rmat_at(i, j, k, int)
                        if (int > max_int_out) then
                            max_int_out = int
                        end if
                    end if
                end do
            end do
        end do
        corr = fit%real_corr(self%img_raw, mask=fit_mask)
        self%atominfo(cc)%anisocorr = corr 

        ! Calculate eigenvalues and orientation of major eigenvector
        call jacobi(cov, 3, 3, eigenvals, eigenvecs, nrot)
        call eigsrt(eigenvals, eigenvecs, 3, 3)
        ! Any negative eigenvalue or large eigenvalue implies intensity 
        ! distribution can't be approximated as a multivariate Gaussian
        if (any(eigenvals <= 0._dp) .or. any(eigenvals > (2*fit_rad)**2)) then
            self%atominfo(cc)%tossADP = .true.
            self%atominfo(cc)%u_evals = 0.
            self%atominfo(cc)%anisocorr = 0.
            self%atominfo(cc)%azimuth = 0.
            self%atominfo(cc)%polar = 0.
            self%atominfo(cc)%doi = 0.
            self%atominfo(cc)%doi_min = 0.
            self%atominfo(cc)%aniso = 0.
            return
        end if
        self%atominfo(cc)%u_evals(:3) = eigenvals(:3)
        self%atominfo(cc)%doi = eigenvals(3)/eigenvals(1)
        self%atominfo(cc)%doi_min = eigenvals(2)/eigenvals(1)
        ! Find azimuthal and polar angles of the major eigenvector
        majvector = eigenvecs(:,1)
        ! Sign of majvector is arbitrary. Use convention that y-coordinate must be >= 0
        if (majvector(2) < 0.) then
            majvector = -1. * majvector
        end if
        write(*, '(i4, 3f8.2)') cc, majvector(:3)
        self%atominfo(cc)%azimuth = atan(majvector(2)/majvector(1))
        if (majvector(1) > 0.) then
            self%atominfo(cc)%azimuth = atan(majvector(2)/majvector(1))
        else if (majvector(1) < 0.) then
            self%atominfo(cc)%azimuth = atan(majvector(2)/majvector(1)) + PI
        else
            self%atominfo(cc)%azimuth = PI / 2
        end if
        if (majvector(3) > 0.) then
            self%atominfo(cc)%polar = atan(sqrt(majvector(1)**2 + majvector(2)**2) / majvector(3))
        else if (majvector(3) < 0.) then
            self%atominfo(cc)%polar = atan(sqrt(majvector(1)**2 + majvector(2)**2) / majvector(3)) + PI
        else
            self%atominfo(cc)%polar = PI / 2
        end if
    end subroutine calc_anisotropic_disp

    ! visualization and output

    subroutine simulate_atoms( self, simatms, betas, mask, atoms_obj)
        class(nanoparticle),           intent(inout) :: self
        class(image),                  intent(inout) :: simatms
        real,        optional,         intent(in)    :: betas(self%n_cc) ! in pdb file b-factor
        logical,     optional,         intent(in)    :: mask(self%n_cc)
        type(atoms), optional, target, intent(inout) :: atoms_obj
        type(atoms), target  :: atms_here
        type(atoms), pointer :: atms_ptr => null()
        logical :: betas_present, mask_present, atoms_obj_present
        integer :: i, cnt
        betas_present     = present(betas)
        mask_present      = present(mask)
        atoms_obj_present = present(atoms_obj)
        if( atoms_obj_present )then
            atms_ptr => atoms_obj
        else
            atms_ptr => atms_here
        endif
        ! generate atoms object
        if( mask_present )then
            cnt = count(mask)
            call atms_ptr%new(cnt)
        else
            call atms_ptr%new(self%n_cc)
        endif
        if( mask_present )then
            cnt = 0
            do i = 1, self%n_cc
                if( mask(i) )then
                    cnt = cnt + 1
                    call set_atom(cnt, i)
                endif
            end do
        else
            do i = 1, self%n_cc
                call set_atom(i, i)
            enddo
        endif
        call simatms%new(self%ldim,self%smpd)
        call atms_ptr%convolve(simatms, cutoff = 8.*self%smpd)
        if( .not. atoms_obj_present ) call atms_ptr%kill

        contains

            subroutine set_atom( atms_obj_ind, ainfo_ind )
                integer, intent(in) :: atms_obj_ind, ainfo_ind
                call atms_ptr%set_name(     atms_obj_ind, self%atom_name)
                call atms_ptr%set_element(  atms_obj_ind, self%element)
                call atms_ptr%set_coord(    atms_obj_ind, (self%atominfo(i)%center(:)-1.)*self%smpd)
                call atms_ptr%set_num(      atms_obj_ind, atms_obj_ind)
                call atms_ptr%set_resnum(   atms_obj_ind, atms_obj_ind)
                call atms_ptr%set_chain(    atms_obj_ind, 'A')
                call atms_ptr%set_occupancy(atms_obj_ind, 1.)
                if( betas_present )then
                    call atms_ptr%set_beta(atms_obj_ind, betas(ainfo_ind))
                else
                    call atms_ptr%set_beta(atms_obj_ind, self%atominfo(ainfo_ind)%cn_gen) ! use generalised coordination number
                endif
            end subroutine set_atom

    end subroutine simulate_atoms

    subroutine write_centers_1( self, fname, coords )
        class(nanoparticle),        intent(inout) :: self
        character(len=*), optional, intent(in)    :: fname
        real,             optional, intent(in)    :: coords(:,:)
        type(atoms) :: centers_pdb
        integer     :: cc
        if( present(coords) )then
            call centers_pdb%new(size(coords, dim = 2), dummy=.true.)
            do cc=1,size(coords, dim = 2)
                call set_atm_info
            enddo
        else
            call centers_pdb%new(self%n_cc, dummy=.true.)
            do cc=1,self%n_cc
                call set_atm_info
            enddo
        endif
        if( present(fname) ) then
            call centers_pdb%writepdb(fname)
        else
           call centers_pdb%writepdb(trim(self%fbody)//'_ATMS')
           write(logfhandle,'(A)') 'output, atomic coordinates:       '//trim(self%fbody)//'_ATMS'
        endif

        contains

            subroutine set_atm_info
                call centers_pdb%set_name(cc,self%atom_name)
                call centers_pdb%set_element(cc,self%element)
                call centers_pdb%set_coord(cc,(self%atominfo(cc)%center(:)-1.)*self%smpd)
                call centers_pdb%set_beta(cc,self%atominfo(cc)%valid_corr) ! use per atom valid corr
                call centers_pdb%set_resnum(cc,cc)
            end subroutine set_atm_info

    end subroutine write_centers_1

    subroutine write_centers_2( self, fname, which )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=*),    intent(in)    :: which ! parameter in the B-factor field of the pdb file
        type(atoms) :: centers_pdb
        integer     :: cc
        call centers_pdb%new(self%n_cc, dummy=.true.)
        do cc=1,self%n_cc
            call centers_pdb%set_name(cc,self%atom_name)
            call centers_pdb%set_element(cc,self%element)
            call centers_pdb%set_coord(cc,(self%atominfo(cc)%center(:)-1.)*self%smpd)
            select case(which)
                case('cn_gen')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%cn_gen)      ! use generalised coordination number
                case('cn_std')
                    call centers_pdb%set_beta(cc,real(self%atominfo(cc)%cn_std))! use standard coordination number
                case('max_int')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%max_int)     ! use z-score of maximum intensity
                case('doi')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%doi)     ! use isotropic b-factor
                case('doi_min')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%doi_min)     ! use isotropic b-factor
                case('u_iso')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%u_iso)     ! use isotropic b-factor
                case DEFAULT
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%valid_corr)  ! use per-atom validation correlation
            end select
            call centers_pdb%set_resnum(cc,cc)
        enddo
        call centers_pdb%writepdb(fname)
    end subroutine write_centers_2

    subroutine write_centers_aniso( self, fname)
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        type(atoms) :: centers_pdb
        real        :: aniso(3, 3, self%n_cc)
        integer     :: cc
        call centers_pdb%new(self%n_cc, dummy=.true.)
        do cc=1,self%n_cc
            call centers_pdb%set_name(cc,self%atom_name)
            call centers_pdb%set_element(cc,self%element)
            call centers_pdb%set_coord(cc,(self%atominfo(cc)%center(:)-1.)*self%smpd)
            call centers_pdb%set_resnum(cc,cc)
            aniso(:,:,cc) = self%atominfo(cc)%aniso(:,:) ! in Angstroms
        enddo
        call centers_pdb%writepdb_aniso(fname, aniso)
    end subroutine write_centers_aniso

    subroutine write_individual_atoms( self )
        class(nanoparticle), intent(inout) :: self
        type(binimage)       :: img_atom
        integer, allocatable :: imat(:,:,:), imat_atom(:,:,:)
        integer :: i
        call img_atom%copy_bimg(self%img_cc)
        allocate(imat_atom(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
        call img_atom%get_imat(imat)
        do i = 1, maxval(imat)
            where(imat == i)
                imat_atom = 1
            elsewhere
                imat_atom = 0
            endwhere
            call img_atom%set_imat(imat_atom)
            call img_atom%write_bimg('Atom'//trim(int2str(i))//'.mrc')
        enddo
        deallocate(imat,imat_atom)
        call img_atom%kill_bimg
    end subroutine write_individual_atoms

    subroutine write_csv_files( self )
        class(nanoparticle), intent(in) :: self
        integer            :: ios, funit, cc, cn
        character(len=256) :: io_msg
        ! NANOPARTICLE STATS
        call fopen(funit, file=NP_STATS_FILE, iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//NP_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        write(funit,'(a)') NP_STATS_HEAD
        ! write record
        call self%write_np_stats(funit)
        call fclose(funit)
        ! CN-DEPENDENT STATS
        call fopen(funit, file=CN_STATS_FILE, iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//CN_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        write(funit,'(a)') CN_STATS_HEAD
        ! write records
        do cn = CNMIN, CNMAX
            call self%write_cn_stats(cn, funit)
        enddo
        call fclose(funit)
        ! PER-ATOM STATS
        call fopen(funit, file=ATOMS_STATS_FILE, iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//ATOMS_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        if (ATOMS_STATS_OMIT) then
            write(funit,'(a)') ATOM_STATS_HEAD_OMIT
        else
            write(funit,'(a)') ATOM_STATS_HEAD
        endif
        ! write records
        do cc = 1, size(self%atominfo)
            call self%write_atominfo(cc, funit, omit=ATOMS_STATS_OMIT)
        enddo
        call fclose(funit)
    end subroutine write_csv_files

    subroutine write_atominfo( self, cc, funit, omit )
        class(nanoparticle), intent(in) :: self
        integer,             intent(in) :: cc, funit
        logical, optional,   intent(in) :: omit
        logical :: omit_here
        if( self%atominfo(cc)%size < NVOX_THRESH ) return
        omit_here = .false.
        if( present(omit) ) omit_here = omit
        601 format(F8.4,A2)
        602 format(F8.4)
        ! various per-atom parameters
        write(funit,601,advance='no') real(self%atominfo(cc)%cc_ind),           CSV_DELIM ! INDEX          (1)
        write(funit,601,advance='no') real(self%atominfo(cc)%size),             CSV_DELIM ! NVOX           (2)
        write(funit,601,advance='no') real(self%atominfo(cc)%cn_std),           CSV_DELIM ! CN_STD         (3)
        write(funit,601,advance='no') self%atominfo(cc)%bondl,                  CSV_DELIM ! NN_BONDL       (4)
        write(funit,601,advance='no') self%atominfo(cc)%cn_gen,                 CSV_DELIM ! CN_GEN         (5)
        write(funit,601,advance='no') self%atominfo(cc)%diam,                   CSV_DELIM ! DIAM           (6)
        write(funit,601,advance='no') real(self%atominfo(cc)%adjacent_cn13),    CSV_DELIM ! ADJ_CN13       (7)
        write(funit,601,advance='no') self%atominfo(cc)%avg_int,                CSV_DELIM ! AVG_INT        (8)
        write(funit,601,advance='no') self%atominfo(cc)%max_int,                CSV_DELIM ! MAX_INT        (9)
        write(funit,601,advance='no') self%atominfo(cc)%cendist,                CSV_DELIM ! CENDIST        (10)
        write(funit,601,advance='no') self%atominfo(cc)%valid_corr,             CSV_DELIM ! VALID_CORR     (11)
        write(funit,601,advance='no') self%atominfo(cc)%u_iso,                  CSV_DELIM ! BFAC           (12)
        ! anisotropic displacement
        write(funit,601,advance='no') self%atominfo(cc)%u_evals(1),             CSV_DELIM ! U_MAJ
        write(funit,601,advance='no') self%atominfo(cc)%u_evals(2),             CSV_DELIM ! U_MED
        write(funit,601,advance='no') self%atominfo(cc)%u_evals(3),             CSV_DELIM ! U_MIN
        write(funit,601,advance='no') self%atominfo(cc)%azimuth,                CSV_DELIM ! AZIMUTH
        write(funit,601,advance='no') self%atominfo(cc)%polar,                  CSV_DELIM ! POLAR
        write(funit,601,advance='no') self%atominfo(cc)%doi,                    CSV_DELIM ! DOI
        write(funit,601,advance='no') self%atominfo(cc)%doi_min,                CSV_DELIM ! DOI_MIN
        ! atomic displacement fit correlations
        write(funit,601,advance='no') self%atominfo(cc)%isocorr,                CSV_DELIM ! ISO_CORR
        write(funit,601,advance='no') self%atominfo(cc)%anisocorr,              CSV_DELIM ! ANISO_CORR
        if( .not. omit_here )then
        write(funit,601,advance='no') self%atominfo(cc)%center(1)*self%smpd,    CSV_DELIM ! X
        write(funit,601,advance='no') self%atominfo(cc)%center(2)*self%smpd,    CSV_DELIM ! Y
        write(funit,601,advance='no') self%atominfo(cc)%center(3)*self%smpd,    CSV_DELIM ! Z
        ! strain
        write(funit,601,advance='no') self%atominfo(cc)%exx_strain,             CSV_DELIM ! EXX_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%eyy_strain,             CSV_DELIM ! EYY_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%ezz_strain,             CSV_DELIM ! EZZ_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%exy_strain,             CSV_DELIM ! EXY_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%eyz_strain,             CSV_DELIM ! EYZ_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%exz_strain,             CSV_DELIM ! EXZ_STRAIN
        endif
        if( .not. omit_here )then
        write(funit,602)              self%atominfo(cc)%radial_strain                     ! RADIAL_STRAIN
        else
        write(funit,602)              self%atominfo(cc)%radial_strain                     ! RADIAL_STRAIN
        endif
    end subroutine write_atominfo

    subroutine write_np_stats( self, funit )
        class(nanoparticle), intent(in) :: self
        integer,             intent(in) :: funit
        601 format(F8.4,A2)
        602 format(F8.4)
        ! -- # atoms
        write(funit,601,advance='no') real(self%n_cc),               CSV_DELIM ! NATOMS               (1)
        write(funit,601,advance='no') real(self%n_aniso),            CSV_DELIM ! NANISO               (2)
        ! -- NP diameter
        write(funit,601,advance='no') self%NPdiam,                   CSV_DELIM ! DIAM                 (3)
        ! -- atom size
        write(funit,601,advance='no') self%size_stats%avg,           CSV_DELIM ! AVG_NVOX             (4)
        write(funit,601,advance='no') self%size_stats%med,           CSV_DELIM ! MED_NVOX             (5)
        write(funit,601,advance='no') self%size_stats%sdev,          CSV_DELIM ! SDEV_NVOX            (6)
        ! -- standard coordination number
        write(funit,601,advance='no') self%cn_std_stats%avg,         CSV_DELIM ! AVG_CN_STD           (7)
        write(funit,601,advance='no') self%cn_std_stats%med,         CSV_DELIM ! MED_CN_STD           (8)
        write(funit,601,advance='no') self%cn_std_stats%sdev,        CSV_DELIM ! SDEV_CN_STD          (9)
        ! -- bond length 
        write(funit,601,advance='no') self%bondl_stats%avg,          CSV_DELIM ! AVG_NN_BONDL         (10)
        write(funit,601,advance='no') self%bondl_stats%med,          CSV_DELIM ! MED_NN_BONDL         (11)
        write(funit,601,advance='no') self%bondl_stats%sdev,         CSV_DELIM ! SDEV_NN_BONDL        (12)
        ! -- generalized coordination number
        write(funit,601,advance='no') self%cn_gen_stats%avg,         CSV_DELIM ! AVG_CN_GEN           (13)
        write(funit,601,advance='no') self%cn_gen_stats%med,         CSV_DELIM ! MED_CN_GEN           (14)
        write(funit,601,advance='no') self%cn_gen_stats%sdev,        CSV_DELIM ! SDEV_CN_GEN          (15)
        ! -- atom diameter
        write(funit,601,advance='no') self%diam_stats%avg,           CSV_DELIM ! AVG_DIAM             (16)
        write(funit,601,advance='no') self%diam_stats%med,           CSV_DELIM ! MED_DIAM             (17)
        write(funit,601,advance='no') self%diam_stats%sdev,          CSV_DELIM ! SDEV_DIAM            (18)
        ! -- average intensity
        write(funit,601,advance='no') self%avg_int_stats%avg,        CSV_DELIM ! AVG_AVG_INT          (19)
        write(funit,601,advance='no') self%avg_int_stats%med,        CSV_DELIM ! MED_AVG_INT          (20)
        write(funit,601,advance='no') self%avg_int_stats%sdev,       CSV_DELIM ! SDEV_AVG_INT         (21)
        ! -- maximum intensity
        write(funit,601,advance='no') self%max_int_stats%avg,        CSV_DELIM ! AVG_MAX_INT          (22)
        write(funit,601,advance='no') self%max_int_stats%med,        CSV_DELIM ! MED_MAX_INT          (23)
        write(funit,601,advance='no') self%max_int_stats%sdev,       CSV_DELIM ! SDEV_MAX_INT         (24)
        ! -- maximum correlation
        write(funit,601,advance='no') self%valid_corr_stats%avg,     CSV_DELIM ! AVG_VALID_CORR       (25)
        write(funit,601,advance='no') self%valid_corr_stats%med,     CSV_DELIM ! MED_VALID_CORR       (26)
        write(funit,601,advance='no') self%valid_corr_stats%sdev,    CSV_DELIM ! SDEV_VALID_CORR      (27)
        ! -- Isotropic displacement parameter
        write(funit,601,advance='no') self%u_iso_stats%avg,          CSV_DELIM ! AVG_U_ISO            (28)
        write(funit,601,advance='no') self%u_iso_stats%med,          CSV_DELIM ! MED_U_ISO            (29)
        write(funit,601,advance='no') self%u_iso_stats%sdev,         CSV_DELIM ! SDEV_U_ISO           (30)
        ! -- anisotropic displacement parameters major eigenvalue
        write(funit,601,advance='no') self%u_maj_stats%avg,          CSV_DELIM ! AVG_U_MAJ            (31)
        write(funit,601,advance='no') self%u_maj_stats%med,          CSV_DELIM ! MED_U_MAJ            (32)
        write(funit,601,advance='no') self%u_maj_stats%sdev,         CSV_DELIM ! SDEV_U_MAJ           (33)
        ! -- anisotropic displacement parameters medium eigenvalue 
        write(funit,601,advance='no') self%u_med_stats%avg,          CSV_DELIM ! AVG_U_MED            (34)
        write(funit,601,advance='no') self%u_med_stats%med,          CSV_DELIM ! MED_U_MED            (35)
        write(funit,601,advance='no') self%u_med_stats%sdev,         CSV_DELIM ! SDEV_U_MED           (36)
        ! -- anisotropic displacement parameters minor eigenvalue
        write(funit,601,advance='no') self%u_min_stats%avg,          CSV_DELIM ! AVG_U_MIN            (37)
        write(funit,601,advance='no') self%u_min_stats%med,          CSV_DELIM ! MED_U_MIN            (38)
        write(funit,601,advance='no') self%u_min_stats%sdev,         CSV_DELIM ! SDEV_U_MIN           (39)
        ! -- Azimuthal angle of major eigenvector 
        write(funit,601,advance='no') self%azimuth_stats%avg,        CSV_DELIM ! AVG_AZIMUTH          (40)
        write(funit,601,advance='no') self%azimuth_stats%med,        CSV_DELIM ! MED_AZIMUTH          (41)
        write(funit,601,advance='no') self%azimuth_stats%sdev,       CSV_DELIM ! SDEV_AZIMUTH         (42)
        ! -- Polar angle of major eigenvector
        write(funit,601,advance='no') self%polar_stats%avg,          CSV_DELIM ! AVG_POLAR            (43)
        write(funit,601,advance='no') self%polar_stats%med,          CSV_DELIM ! MED_POLAR            (44)
        write(funit,601,advance='no') self%polar_stats%sdev,         CSV_DELIM ! SDEV_POLAR           (45)
        ! -- degree of isotropy
        write(funit,601,advance='no') self%doi_stats%avg,            CSV_DELIM ! AVG_DOI              (46)
        write(funit,601,advance='no') self%doi_stats%med,            CSV_DELIM ! MED_DOI              (47)
        write(funit,601,advance='no') self%doi_stats%sdev,           CSV_DELIM ! SDEV_DOI             (48)
        ! -- degree of isotropy minimum w/r to anisotropy average
        write(funit,601,advance='no') self%doi_min_stats%avg,        CSV_DELIM ! AVG_DOI_MIN          (49)
        write(funit,601,advance='no') self%doi_min_stats%med,        CSV_DELIM ! MED_DOI_MIN          (50)
        write(funit,601,advance='no') self%doi_min_stats%sdev,       CSV_DELIM ! SDEV_DOI_MIN         (51)
        ! -- Isotropic displacement fit correlation
        write(funit,601,advance='no') self%isocorr_stats%avg,        CSV_DELIM ! AVG_ISO_CORR         (52)
        write(funit,601,advance='no') self%isocorr_stats%med,        CSV_DELIM ! MED_ISO_CORR         (53)
        write(funit,601,advance='no') self%isocorr_stats%sdev,       CSV_DELIM ! SDEV_ISO_CORR        (54)
        ! -- Anisotropic displacement fit correlation
        write(funit,601,advance='no') self%anisocorr_stats%avg,      CSV_DELIM ! AVG_ANISO_CORR       (55)
        write(funit,601,advance='no') self%anisocorr_stats%med,      CSV_DELIM ! MED_ANISO_CORR       (56)
        write(funit,601,advance='no') self%anisocorr_stats%sdev,     CSV_DELIM ! SDEV_ANISO_CORR      (57)
        ! -- radial strain
        write(funit,601,advance='no') self%radial_strain_stats%avg,  CSV_DELIM ! AVG_RADIAL_STRAIN    (58)
        write(funit,601,advance='no') self%radial_strain_stats%med,  CSV_DELIM ! MED_RADIAL_STRAIN    (59)
        write(funit,601,advance='no') self%radial_strain_stats%sdev, CSV_DELIM ! SDEV_RADIAL_STRAIN   (60)
        write(funit,601,advance='no') self%radial_strain_stats%minv, CSV_DELIM ! MIN_RADIAL_STRAIN    (61)
        write(funit,602)              self%radial_strain_stats%maxv            ! MAX_RADIAL_STRAIN    (62)
    end subroutine write_np_stats

    subroutine write_cn_stats( self, cn, funit )
        class(nanoparticle), intent(in) :: self
        integer,             intent(in) :: cn, funit
        601 format(F8.4,A2)
        602 format(F8.4)
        if( count(self%atominfo(:)%cn_std == cn) < 2 ) return
        ! -- coordination number
        write(funit,601,advance='no') real(cn),                              CSV_DELIM ! CN_STD               (1)
        ! -- # atoms per cn
        write(funit,601,advance='no') self%natoms_cns(cn),                   CSV_DELIM ! NATOMS               (2)
        write(funit,601,advance='no') self%natoms_aniso_cns(cn),             CSV_DELIM ! NANISO               (3)
        ! -- atom size
        write(funit,601,advance='no') self%size_stats_cns(cn)%avg,           CSV_DELIM ! AVG_NVOX             (4)
        write(funit,601,advance='no') self%size_stats_cns(cn)%med,           CSV_DELIM ! MED_NVOX             (5)
        write(funit,601,advance='no') self%size_stats_cns(cn)%sdev,          CSV_DELIM ! SDEV_NVOX            (6)
        ! -- bond length
        write(funit,601,advance='no') self%bondl_stats_cns(cn)%avg,          CSV_DELIM ! AVG_NN_BONDL         (7)
        write(funit,601,advance='no') self%bondl_stats_cns(cn)%med,          CSV_DELIM ! MED_NN_BONDL         (8)
        write(funit,601,advance='no') self%bondl_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_NN_BONDL        (9)
        ! -- generalized coordination number
        write(funit,601,advance='no') self%cn_gen_stats_cns(cn)%avg,         CSV_DELIM ! AVG_CN_GEN           (10)
        write(funit,601,advance='no') self%cn_gen_stats_cns(cn)%med,         CSV_DELIM ! MED_CN_GEN           (11)
        write(funit,601,advance='no') self%cn_gen_stats_cns(cn)%sdev,        CSV_DELIM ! SDEV_CN_GEN          (12)
        ! -- atom diameter
        write(funit,601,advance='no') self%diam_stats_cns(cn)%avg,           CSV_DELIM ! AVG_DIAM             (13)
        write(funit,601,advance='no') self%diam_stats_cns(cn)%med,           CSV_DELIM ! MED_DIAM             (14)
        write(funit,601,advance='no') self%diam_stats_cns(cn)%sdev,          CSV_DELIM ! SDEV_DIAM            (15)
        ! -- average intensity
        write(funit,601,advance='no') self%avg_int_stats_cns(cn)%avg,        CSV_DELIM ! AVG_AVG_INT          (16)
        write(funit,601,advance='no') self%avg_int_stats_cns(cn)%med,        CSV_DELIM ! MED_AVG_INT          (17)
        write(funit,601,advance='no') self%avg_int_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_AVG_INT         (18)
        ! -- maximum intensity
        write(funit,601,advance='no') self%max_int_stats_cns(cn)%avg,        CSV_DELIM ! AVG_MAX_INT          (19)
        write(funit,601,advance='no') self%max_int_stats_cns(cn)%med,        CSV_DELIM ! MED_MAX_INT          (20)
        write(funit,601,advance='no') self%max_int_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_MAX_INT         (21)
        ! -- maximum correlation
        write(funit,601,advance='no') self%valid_corr_stats_cns(cn)%avg,     CSV_DELIM ! AVG_VALID_CORR       (22)
        write(funit,601,advance='no') self%valid_corr_stats_cns(cn)%med,     CSV_DELIM ! MED_VALID_CORR       (23)
        write(funit,601,advance='no') self%valid_corr_stats_cns(cn)%sdev,    CSV_DELIM ! SDEV_VALID_CORR      (24)
        ! -- Isotropic displacement parameter
        write(funit,601,advance='no') self%u_iso_stats_cns(cn)%avg,          CSV_DELIM ! AVG_U_ISO            (25)
        write(funit,601,advance='no') self%u_iso_stats_cns(cn)%med,          CSV_DELIM ! MED_U_ISO            (26)
        write(funit,601,advance='no') self%u_iso_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_U_ISO           (27)
        ! -- anisotropic displacement parameters major eigenvalue
        write(funit,601,advance='no') self%u_maj_stats_cns(cn)%avg,          CSV_DELIM ! AVG_U_MAJ
        write(funit,601,advance='no') self%u_maj_stats_cns(cn)%med,          CSV_DELIM ! MED_U_MAJ
        write(funit,601,advance='no') self%u_maj_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_U_MAJ
        ! -- anisotropic displacement parameters medium eigenvalue
        write(funit,601,advance='no') self%u_med_stats_cns(cn)%avg,          CSV_DELIM ! AVG_U_MED
        write(funit,601,advance='no') self%u_med_stats_cns(cn)%med,          CSV_DELIM ! MED_U_MED
        write(funit,601,advance='no') self%u_med_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_U_MED
        ! -- anisotropic displacement parameters minor eigenvalue
        write(funit,601,advance='no') self%u_min_stats_cns(cn)%avg,          CSV_DELIM ! AVG_U_MIN
        write(funit,601,advance='no') self%u_min_stats_cns(cn)%med,          CSV_DELIM ! MED_U_MIN
        write(funit,601,advance='no') self%u_min_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_U_MIN
        ! -- Azimuthal angle of major eigenvector
        write(funit,601,advance='no') self%azimuth_stats_cns(cn)%avg,        CSV_DELIM ! AVG_AZIMUTH
        write(funit,601,advance='no') self%azimuth_stats_cns(cn)%med,        CSV_DELIM ! MED_AZIMUTH
        write(funit,601,advance='no') self%azimuth_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_AZIMUTH
        ! -- Polar angle of major eigenvector
        write(funit,601,advance='no') self%polar_stats_cns(cn)%avg,          CSV_DELIM ! AVG_POLAR
        write(funit,601,advance='no') self%polar_stats_cns(cn)%med,          CSV_DELIM ! MED_POLAR
        write(funit,601,advance='no') self%polar_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_POLAR
        ! -- degree of isotropy
        write(funit,601,advance='no') self%doi_stats_cns(cn)%avg,            CSV_DELIM ! AVG_DOI
        write(funit,601,advance='no') self%doi_stats_cns(cn)%med,            CSV_DELIM ! MED_DOI
        write(funit,601,advance='no') self%doi_stats_cns(cn)%sdev,           CSV_DELIM ! SDEV_DOI   
        ! -- degree of isotropy min w/r to average anisotropy
        write(funit,601,advance='no') self%doi_min_stats_cns(cn)%avg,        CSV_DELIM ! AVG_DOI_MIN
        write(funit,601,advance='no') self%doi_min_stats_cns(cn)%med,        CSV_DELIM ! MED_DOI_MIN
        write(funit,601,advance='no') self%doi_min_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_DOI_MIN
        ! -- isotropic displacement fit correlation
        write(funit,601,advance='no') self%isocorr_stats_cns(cn)%avg,        CSV_DELIM ! AVG_ISO_CORR
        write(funit,601,advance='no') self%isocorr_stats_cns(cn)%med,        CSV_DELIM ! MED_ISO_CORR
        write(funit,601,advance='no') self%isocorr_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_ISO_CORR
        ! -- anisotropic displacement fit correlation
        write(funit,601,advance='no') self%anisocorr_stats_cns(cn)%avg,      CSV_DELIM ! AVG_ANISO_CORR
        write(funit,601,advance='no') self%anisocorr_stats_cns(cn)%med,      CSV_DELIM ! MED_ANISO_CORR
        write(funit,601,advance='no') self%anisocorr_stats_cns(cn)%sdev,     CSV_DELIM ! SDEV_ANISO_CORR
        ! -- radial strain
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%avg,  CSV_DELIM ! AVG_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%med,  CSV_DELIM ! MED_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%sdev, CSV_DELIM ! SDEV_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%minv, CSV_DELIM ! MIN_RADIAL_STRAIN
        write(funit,602)              self%radial_strain_stats_cns(cn)%maxv            ! MAX_RADIAL_STRAIN
    end subroutine write_cn_stats

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        call self%img%kill()
        call self%img_raw%kill
        call self%img_bin%kill_bimg()
        call self%img_cc%kill_bimg()
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
    end subroutine kill_nanoparticle

end module simple_nanoparticle

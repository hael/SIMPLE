
!@descr: the nanoparticle abstract data type, used for automated atomic model building in SINGLE
module simple_nanoparticle
use simple_core_module_api
use simple_defs_atoms
use simple_image,      only: image
use simple_image_msk,  only: image_msk
use simple_image_bin,  only: image_bin
use simple_atoms,      only: atoms
use simple_parameters, only: parameters
!use simple_linalg
use simple_nanoparticle_utils
implicit none

public :: nanoparticle
private
#include "simple_local_flags.inc"

logical :: DEBUG = .false.

! module global constants
integer,          parameter :: NBIN_THRESH         = 20      ! number of thresholds for binarization
integer,          parameter :: NVOX_THRESH         = 3       ! min # voxels per atom is 3
logical,          parameter :: WRITE_OUTPUT        = .false. ! for figures generation
logical,          parameter :: ATOMS_STATS_OMIT    = .false. ! omit = shorter atoms stats output
integer,          parameter :: SOFT_EDGE           = 6
integer,          parameter :: CNMIN               = 3
integer,          parameter :: CNMAX               = 13
integer,          parameter :: NSTRAIN_COMPS       = 7
character(len=*), parameter :: ATOMS_STATS_FILE    = 'atoms_stats.csv'
character(len=*), parameter :: NP_STATS_FILE       = 'nanoparticle_stats.csv'
character(len=*), parameter :: CN_STATS_FILE       = 'cn_dependent_stats.csv'
! species-free detection (no element): pseudo-atom template and length scales from the nearest-neighbour distance d_NN
real,             parameter :: B_REF_PSEUDO        = 13.9    ! B factor of the pseudo-atom template, A**2 (sigma 0.42 A)
real,             parameter :: RADIUS_NN_FRAC      = 0.4     ! atom radius in units of d_NN
real,             parameter :: SPLIT_EXCL_NN_FRAC  = 0.21    ! exclusion radius of atom splitting in units of d_NN
integer,          parameter :: CSCORE_CEIL_NOEL    = 12      ! contact-score ceiling of close packing

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
    real    :: bondl             = 0. ! nearest neighbor bond length in A                           NN_BONDL
    real    :: cn_gen            = 0. ! generalized coordination number                             CN_GEN
    real    :: diam              = 0. ! atom diameter                                               DIAM
    real    :: avg_int           = 0. ! average grey level intensity across the connected component AVG_INT
    real    :: max_int           = 0. ! maximum            -"-                                      MAX_INT
    real    :: cendist           = 0. ! distance from the center of mass of the nanoparticle        CENDIST
    real    :: valid_corr        = 0. ! per-atom correlation with the simulated map                 VALID_CORR
    real    :: u_iso             = 0. ! isotropic displacement parameters (b-factor)                U_ISO
    real    :: u_evals(3)        = 0. ! eigenvalues of aniso displacement matrix (maj,med,min)      U_(MAJ,MED,MIN)
    real    :: azimuth           = 0. ! azimuthal angle of major eigenvector [0,Pi)                 AZIMUTH
    real    :: polar             = 0. ! polar angle of major eigenvector [0, Pi)                    POLAR
    real    :: doi               = 0. ! degree of isotropy                                          DOI
    real    :: doi_min           = 0. ! degree of isotropy ratio of the two smallest eigenvalues    DOI_MIN
    real    :: isocorr           = 0. ! correlation of isotropic B-Factor fit to input map          ISO_CORR
    real    :: anisocorr         = 0. ! correlation of anisotropic B-Factor fit to input map        ANISO_CORR
    real    :: center(3)         = 0. ! atom center                                                 X Y Z
    ! strain
    real    :: exx_strain        = 0. ! tensile strain in %                                         EXX_STRAIN
    real    :: eyy_strain        = 0. ! -"-                                                         EYY_STRAIN
    real    :: ezz_strain        = 0. ! -"-                                                         EZZ_STRAIN
    real    :: exy_strain        = 0. ! -"-                                                         EXY_STRAIN
    real    :: eyz_strain        = 0. ! -"-                                                         EYZ_STRAIN
    real    :: exz_strain        = 0. ! -"-                                                         EXZ_STRAIN
    real    :: radial_strain     = 0. ! -"-                                                         RADIAL_STRAIN
    ! Auxiliary (non-output)
    real    :: aniso(3,3)        = 0. ! ADP Matrix for ANISOU PDB file                              N/A
    logical :: tossADP           = .false. ! true if atom inadequate for ADP calculations           N/A
    ! species discovery (discover_species=yes), written to the _species files only
    real    :: amp               = 0. ! fitted Gaussian amplitude
    real    :: bfac              = 0. ! B factor of the class-tied fit (stage 2), A**2
    real    :: bfac_stage1       = 0. ! B factor of the free-amplitude fit (stage 1), A**2
    real    :: aper_int          = 0. ! aperture intensity
    real    :: aper_snr          = 0. ! aperture intensity over its noise
    real    :: det_z             = 0. ! detection z of a recovered atom
    integer :: det_stage         = 0  ! 0 level 1, 1 residual stage A, 2 residual stage B
    integer :: det_level         = 0  ! residual level that accepted the atom
    logical :: gate_used         = .false. ! accepted through the neighbour gate
    integer :: species           = 0  ! intensity class, numbered by decreasing intensity
end type atom_stats

type :: nanoparticle
    private
    type(image)           :: img, img_raw
    type(image_bin)       :: img_bin, img_cc         ! binary and connected component images
    integer               :: ldim(3)            = 0  ! logical dimension of image
    integer               :: n_cc               = 0  ! number of atoms (connected components)                NATOMS
    integer               :: n_aniso            = 0  ! number of atoms with aniso calculations               NANISO
    integer               :: n4stats            = 0  ! number of atoms in subset used for stats calc
    real                  :: smpd               = 0. ! sampling distance
    real                  :: NPcen(3)           = 0. ! coordinates of the center of mass of the nanoparticle
    real                  :: NPdiam             = 0. ! diameter of the nanoparticle                          DIAM
    real                  :: theoretical_radius = 0. ! theoretical atom radius in A
    real                  :: d_nn               = 0. ! median nearest-neighbour distance in A (species-free path)
    ! the pruning policy of discard_atoms, kept for the recovered atoms of discovery: centre of mass (voxels), the
    ! radius beyond which atoms are subject to deletion (A) and the contact-score threshold
    real                  :: prune_cen(3)       = 0.
    real                  :: prune_rad_thres    = 0.
    integer               :: prune_cs_thres     = 0
    logical               :: l_species_free     = .false. ! no element given: pseudo-atom template, d_NN length scales
    ! species discovery (discover_species=yes)
    logical               :: l_discover_species = .false.
    integer               :: min_nbrs           = 3       ! recovered atoms need this many found neighbours, 0 = no gate
    integer               :: nspecies           = 0  ! 0: number of species from the data
    real                  :: msk_rad            = 0. ! radius of the soft spherical mask in voxels
    type(string)          :: vol_even, vol_odd       ! optional half maps, the noise reference
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
    real                  :: natoms_cns(CNMIN:CNMAX)       = 0. ! # of atoms per cn_std                            NATOMS
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
    type(atom_stats), allocatable :: recovered(:)    ! atoms recovered from the residual map, diagnostics only
    real,             allocatable :: species_post(:,:) ! class posteriors of the level-1 then the recovered atoms
    real,             allocatable :: coords4stats(:,:)
    ! OTHER
    character(len=2)      :: element     = '  '      ! symbol of the first species: radii, templates, atom names
    character(len=5)      :: element_key = '     '   ! element as given (a compound selector): lattice lookups
    character(len=4)      :: atom_name   = '    '
    character(len=2), allocatable :: species(:)      ! species named by a list or compound selector, brightest first
    type(string)          :: npname
    type(string)          :: fbody
  contains
    ! constructor
    procedure          :: new
    ! getters/setters
    procedure          :: get_ldim
    procedure          :: get_natoms
    procedure          :: get_valid_corrs
    procedure          :: get_img
    procedure          :: get_img_raw
    procedure          :: set_ncc
    procedure          :: set_half_maps
    procedure          :: set_img
    procedure, private :: set_atomic_coords_from_pdb
    procedure, private :: set_atomic_coords_from_xyz
    generic            :: set_atomic_coords => set_atomic_coords_from_pdb, set_atomic_coords_from_xyz
    procedure          :: set_coords4stats
    procedure, private :: pack_instance4stats
    ! utils
    procedure, private :: atominfo2centers
    procedure, private :: atominfo2centers_A
    procedure, private :: center_on_atom
    procedure          :: update_ncc
    ! atomic position determination
    procedure          :: conv_denoise 
    procedure          :: identify_lattice_params
    procedure          :: identify_atomic_pos
    procedure, private :: binarize_and_find_centers
    procedure          :: find_centers
    procedure, private :: discard_small_ccs
    procedure, private :: discard_atoms
    procedure, private :: set_nn_length_scales
    procedure, private :: contact_scores
    procedure, private :: discover_species
    procedure, private :: split_atoms
    procedure          :: validate_atoms
    ! calc stats
    procedure          :: fillin_atominfo
    procedure, private :: masscen
    procedure, private :: calc_longest_atm_dist
    procedure, private :: check_neighbors_cn
    procedure, private :: binary_lattice => np_binary_lattice
    procedure, private :: lattice_bond   => np_lattice_bond
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
    procedure          :: kill
end type nanoparticle

contains

    subroutine new( self, params, fname, msk )
        class(nanoparticle), intent(inout) :: self
        class(parameters),   intent(in)    :: params
        class(string),       intent(in)    :: fname
        real, optional,      intent(in)    :: msk
        type(image_msk)  :: mskvol
        character(len=2) :: el_ucase
        integer :: nptcls
        integer :: Z ! atomic number
        real    :: msk_in_pix
        call self%kill
        self%npname    = fname
        self%fbody     = get_fbody(basename(fname), fname2ext(fname))
        self%smpd      = params%smpd
        self%l_species_free = len_trim(params%element) == 0
        if( self%l_species_free )then
            ! the radius follows from d_NN once the first binarisation has found centres
            self%atom_name   = 'X1  '
            self%element     = 'X1'
            self%element_key = 'X1'
            self%theoretical_radius = 0.
        else
            self%element_key = params%element(1:len(self%element_key))
            self%element     = params%element(1:len(self%element))
            self%atom_name   = self%element
            el_ucase         = upperCase(self%element)
            call get_element_Z_and_radius(el_ucase, Z, self%theoretical_radius)
            if( Z == 0 ) THROW_HARD('Unknown element: '//el_ucase)
        endif
        call find_ldim_nptcls(self%npname, self%ldim, nptcls)
        call self%img%new(self%ldim, self%smpd)
        call self%img_bin%new_bimg(self%ldim, self%smpd)
        call self%img%read(fname)
        if( present(msk) )then
            call self%img%mask3D_soft(msk)
            self%msk_rad = msk
        else
            call mskvol%estimate_spher_mask_diam(params, self%img, AMSKLP_NANO, msk_in_pix)
            write(logfhandle,*) 'mask diameter in A: ', 2. * msk_in_pix * self%smpd
            call self%img%mask3D_soft(msk_in_pix)
            call mskvol%kill_bimg
            self%msk_rad = msk_in_pix
        endif
        self%l_discover_species = params%l_discover_species
        self%min_nbrs           = params%min_nbrs
        self%nspecies           = params%nspecies
        if( allocated(params%species) ) self%species = params%species
        if( DEBUG ) call self%img%write(string('masked_input_vol.mrc'))
        call self%img_raw%copy(self%img)
        call self%img_raw%stats(self%map_stats%avg, self%map_stats%sdev, self%map_stats%maxv, self%map_stats%minv)
    end subroutine new

    ! getters/setters

    subroutine set_half_maps( self, vol_even, vol_odd )
        class(nanoparticle), intent(inout) :: self
        class(string),       intent(in)    :: vol_even, vol_odd
        self%vol_even = vol_even
        self%vol_odd  = vol_odd
    end subroutine set_half_maps

    subroutine get_ldim( self, ldim )
        class(nanoparticle), intent(in)  :: self
        integer,             intent(out) :: ldim(3)
        ldim = self%img%get_ldim()
    end subroutine get_ldim

    function get_natoms(self) result(n)
       class(nanoparticle), intent(inout)  :: self
       integer :: n
       call self%img_cc%get_nccs(n)
    end function get_natoms

    ! returning center positions alongside valid corr to 
    function get_valid_corrs( self ) result( corrs )
        class(nanoparticle), intent(in) :: self
        real, allocatable :: corrs(:)
        if( allocated(self%atominfo) )then
            allocate(corrs(size(self%atominfo)), source = self%atominfo(:)%valid_corr)
        endif
    end function get_valid_corrs

    subroutine get_img( self, img )
        class(nanoparticle), intent(in)  :: self
        type(image),         intent(out) :: img
        img = self%img
    end subroutine get_img

    subroutine get_img_raw( self, raw_img )
        class(nanoparticle), intent(in)  :: self
        type(image),         intent(out) :: raw_img
        raw_img = self%img_raw
    end subroutine get_img_raw

    subroutine set_ncc( self, ncc )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: ncc
        self%n_cc = ncc
    end subroutine set_ncc

    ! set one of the images of the nanoparticle type
    subroutine set_img( self, imgfile, which )
        class(nanoparticle), intent(inout) :: self
        class(string),       intent(in)    :: imgfile
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

    subroutine set_atomic_coords_from_xyz( self, xyz )
        class(nanoparticle), intent(inout) :: self
        real,                intent(in)    :: xyz(:,:)
        type(atoms) :: a
        integer     :: N, i
        if( size(xyz,dim=2) /= 3 ) THROW_HARD("Error! Non-conforming dimensions of xyz; set_atomic_coords_from_xyz")
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        N = size(xyz, dim=1)
        allocate(self%atominfo(N))
        call a%new(N)
        do i = 1, N
            call a%set_coord(i,xyz(i,:))
            self%atominfo(i)%center(:) = a%get_coord(i) 
        enddo
        self%n_cc = N
        call a%kill
    end subroutine set_atomic_coords_from_xyz

    ! sets the atom positions to be the ones in the inputted PDB file.
    subroutine set_atomic_coords_from_pdb( self, pdb_file )
        class(nanoparticle),     intent(inout) :: self
        class(string),           intent(in)    :: pdb_file
        type(atoms) :: a
        integer     :: i, N
        if( fname2ext(pdb_file) .ne. 'pdb' ) THROW_HARD('Inputted filename has to have pdb extension; set_atomic_coords_from_pdb')
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        call a%new(pdb_file)
        N = a%get_n() ! number of atoms
        allocate(self%atominfo(N))
        do i = 1, N
            self%atominfo(i)%center(:) = a%get_coord(i)/self%smpd + 1.
        enddo
        self%n_cc = N
        call a%kill
    end subroutine set_atomic_coords_from_pdb

    subroutine set_coords4stats( self, pdb_file )
        class(nanoparticle), intent(inout) :: self
        class(string),       intent(in)    :: pdb_file
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
        sz           = size(self%atominfo)
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
            enddo
        else
            allocate(centers(3,sz), source=0.)
            do i = 1, sz
                centers(:,i) = self%atominfo(i)%center(:)
            enddo
        endif
    end function atominfo2centers

    function atominfo2centers_A( self, mask ) result( centers_A )
        class(nanoparticle), intent(in) :: self
        logical, optional,   intent(in) :: mask(size(self%atominfo))
        real, allocatable :: centers_A(:,:)
        logical :: mask_present
        integer :: sz, i, cnt
        sz           = size(self%atominfo)
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
            enddo
        else
            allocate(centers_A(3,sz), source=0.)
            do i = 1, sz
                centers_A(:,i) = (self%atominfo(i)%center(:) - 1.) * self%smpd
            enddo
        endif
    end function atominfo2centers_A

    ! Translate the identified atomic positions so that the center of mass
    ! of the nanoparticle coincides with its closest atom
    subroutine center_on_atom( self, pdbfile_in, pdbfile_out )
        class(nanoparticle), intent(inout) :: self
        class(string),       intent(in)    :: pdbfile_in
        class(string),       intent(inout) :: pdbfile_out
        type(atoms) :: atom_centers
        real        :: m(3), vec(3), d, d_before
        integer     :: i
        call atom_centers%new(pdbfile_in)
        m(:)     = self%masscen()
        d_before = huge(d_before)
        vec      = 0.
        do i = 1, self%n_cc
            d = euclid(m,self%atominfo(i)%center(:))
            if( d < d_before )then
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
        type(image_bin), optional, intent(inout) :: img_cc
        if( present(img_cc) )then
            call img_cc%get_nccs(self%n_cc)
        else
            call self%img_cc%get_nccs(self%n_cc)
        endif
    end subroutine update_ncc

    ! atomic position determination

    subroutine conv_denoise( self, fname )
        class(nanoparticle), intent(inout) :: self
        class(string),       intent(in)    :: fname 
        if( self%l_species_free ) THROW_HARD('conv_denoise needs an element')
        call phasecorr_one_atom(self%img, self%element)
        call self%img%write(fname)
    end subroutine conv_denoise

    subroutine identify_lattice_params( self, a )
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: a(3) ! lattice parameters
        real, allocatable :: centers_A(:,:)        ! coordinates of the atoms in ANGSTROMS
        type(image)       :: simatms
        if( self%l_species_free ) THROW_HARD('identify_lattice_params needs an element')
        ! MODEL BUILDING
        ! phase correlation approach
        call phasecorr_one_atom(self%img, self%element)
        ! nanoparticle binarization
        call self%binarize_and_find_centers()
        ! discard small connected components
        call self%discard_small_ccs
        ! atom splitting by correlation map validation
        call self%split_atoms()
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms)
        call self%validate_atoms(simatms, l_print=.true.)
        ! discard atoms
        call self%discard_atoms
        ! fit lattice
        centers_A = self%atominfo2centers_A()
        call fit_lattice(self%element_key, centers_A, a)
        deallocate(centers_A)
        call simatms%kill
    end subroutine identify_lattice_params

    ! 1. take different position along the mrc volume and compute bfactors for the peaks
    ! 2. do a clutering in two classes, expecte for Cd and Se
    ! 3. simulate atomic model based on those positions and b-factors
    ! 4. compute correlation of each atom with the simulated density
    ! 5. discard atoms with low correlation and split atoms with high correlation but large volumes (i.e. likely to be two merged atoms)

    subroutine identify_atomic_pos( self, a, l_atom_thres, split_fname, l_print )
        class(nanoparticle),     intent(inout) :: self
        real,                    intent(inout) :: a(3)                ! lattice parameters
        logical,                 intent(in)    :: l_atom_thres        ! do atomic thresholding or not
        class(string), optional, intent(in)    :: split_fname
        logical,       optional, intent(in)    :: l_print
        type(image)  :: simatms, img_cos
        type(string) :: errmsg
        logical      :: ll_print
        ll_print = .true.
        if( present( l_print) ) ll_print = l_print
        ! MODEL BUILDING
        ! Phase correlation approach
        if( self%l_species_free )then
            call phasecorr_one_atom(self%img, self%element, bfac_ref=B_REF_PSEUDO)
        else
            call phasecorr_one_atom(self%img, self%element)
        endif
        if( DEBUG ) call self%img%write(string('after_phasecorr.mrc'))
        ! Nanoparticle binarization
        call self%binarize_and_find_centers(l_print=ll_print)
        if( self%l_discover_species .and. self%n_cc < 2 )then
            errmsg = 'discover_species needs two or more level-1 atoms to measure d_NN; level 1 found '//int2str(self%n_cc)
            THROW_HARD(errmsg%to_char())
        endif
        if( self%l_species_free .or. self%l_discover_species ) call self%set_nn_length_scales
        ! discard small connected components
        call self%discard_small_ccs
        ! atom splitting by correlation map validation
        call self%split_atoms(split_fname)
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms)
        if( WRITE_OUTPUT ) call simatms%write(self%fbody//'_SIM_pre_validate.mrc')
        call self%validate_atoms(simatms, l_print=.false.)
        if( WRITE_OUTPUT ) call self%write_centers(string('valid_corr_in_bfac_field_pre_discard.pdb'), 'valid_corr')
        if( l_atom_thres ) call self%discard_atoms(l_print=ll_print)
        ! re-calculate valid_corr:s (since they are otherwise lost from the B-factor field due to reallocations of atominfo)
        call self%simulate_atoms(simatms)
        call self%validate_atoms(simatms, l_print=ll_print)
        ! WRITE OUTPUT
        call self%img_bin%write_bimg(self%fbody//'_BIN.mrc')
        if( ll_print ) write(logfhandle,'(A)') 'output, binarized map:            '//self%fbody%to_char()//'_BIN.mrc'
        call self%img_bin%grow_bins(1)
        call self%img_bin%cos_edge(SOFT_EDGE, img_cos)
        call img_cos%write(self%fbody//'_MSK.mrc')
        if( ll_print ) write(logfhandle,'(A)') 'output, envelope mask map:        '//self%fbody%to_char()//'_MSK.mrc'
        call self%img_cc%write_bimg(self%fbody//'_CC.mrc')
        if( ll_print ) write(logfhandle,'(A)') 'output, connected components map: '//self%fbody%to_char()//'_CC.mrc'
        call self%write_centers
        call simatms%write(self%fbody//'_SIM.mrc')
        if( ll_print ) write(logfhandle,'(A)') 'output, simulated atomic density: '//self%fbody%to_char()//'_SIM.mrc'
        ! the present products are written; discovery reads the level-1 atoms and writes its own files only
        if( self%l_discover_species ) call self%discover_species
        ! destruct
        call img_cos%kill
        call simatms%kill
    end subroutine identify_atomic_pos

    ! This subrotuine takes in input a nanoparticle and
    ! binarizes it by thresholding. The gray level histogram is split
    ! in 20 parts, which corresponds to 20 possible thresholds
    ! Among those thresholds, the selected one is the for which
    ! that correlation between the raw map and a simulated distribution
    ! obtained with that threshold reaches the maximum value.
    subroutine binarize_and_find_centers( self, l_print )
        class(nanoparticle), intent(inout) :: self
        logical, optional,   intent(in)    :: l_print
        type(image_bin)      :: img_bin_t
        type(image_bin)      :: img_ccs_t
        type(atoms)          :: atom
        type(image)          :: simulated_distrib
        integer, allocatable :: imat_t(:,:,:)
        real,    allocatable :: coords(:,:)
        real,    allocatable :: rmat(:,:,:)
        logical, parameter   :: L_BENCH = .false.
        logical :: ll_print
        real    :: ts(NBIN_THRESH)
        integer :: fnr, low, high, mid
        real    :: corr, max_corr, thresh_opt
        real(timer_int_kind)    :: rt_find_ccs, rt_find_centers, rt_gen_sim, rt_real_corr, rt_tot
        integer(timer_int_kind) ::  t_find_ccs,  t_find_centers,  t_gen_sim,  t_real_corr,  t_tot
        ll_print = .true.
        if( present(l_print) ) ll_print = l_print
        if( ll_print ) write(logfhandle,'(A)') '>>> BINARIZATION'
        rmat = self%img%get_rmat()
        allocate(imat_t(self%ldim(1), self%ldim(2), self%ldim(3)), source = 0)
        call simulated_distrib%new(self%ldim,self%smpd)
        rt_find_ccs     =  0.
        rt_find_centers =  0.
        rt_gen_sim      =  0.
        rt_real_corr    =  0.
        t_tot           =  tic()
        call thres_detect_conv_atom_denoised(self%img, NBIN_THRESH, ts)
        thresh_opt = ts(1)
        max_corr   = t2c(thresh_opt)
        ! if( DEBUG )then
        !     ! exhaustive evaluation of all thresholds for debugging
        !     do ithres = 2, NBIN_THRESH
        !         corr = t2c(ts(ithres))
        !         if( ll_print ) write(logfhandle,*) 'threshold: ', ts(ithres), 'corr: ', corr
        !         if( corr > max_corr )then
        !             max_corr   = corr
        !             thresh_opt = ts(ithres)
        !         endif
        !     enddo
        ! else
            low  = 1
            high = NBIN_THRESH
            do while( low <= high )
                mid  = (low + high) / 2
                corr = t2c(ts(mid))
                if( ll_print ) write(logfhandle,*) 'threshold: ', ts(mid), 'corr: ', corr
                if( corr > max_corr )then
                    max_corr   = corr
                    thresh_opt = ts(mid)
                    low        = mid + 1
                else
                    high       = mid - 1
                endif
            enddo
        ! endif
        rt_tot = toc(t_tot)
        if( L_BENCH )then
            call fopen(fnr, FILE=string('BINARIZE_AND_FIND_CENTERS_BENCH.txt'), STATUS='REPLACE', action='WRITE')
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
        if( ll_print ) write(logfhandle,*) 'optimal threshold: ', thresh_opt, 'max_corr: ', max_corr
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
        if( ll_print ) write(logfhandle,'(A)') '>>> BINARIZATION, COMPLETED'

    contains

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
            call self%find_centers(img_ccs_t, coords)
            rt_find_centers = rt_find_centers + toc(t_find_centers)
            ! Generate a simulated distribution based on those center
            t_gen_sim = tic()
            call self%write_centers(string('centers_iteration.pdb'), coords)
            call atom%new(string('centers_iteration.pdb'))
            if( self%l_species_free ) call set_pseudo_bfacs(atom)
            call atom%convolve(simulated_distrib, cutoff = 8.*self%smpd)
            call del_file('centers_iteration.pdb')
            call atom%kill
            rt_gen_sim = rt_gen_sim + toc(t_gen_sim)
            ! correlate volumes
            t_real_corr = tic()
            t2c = self%img%real_corr(simulated_distrib)
            rt_real_corr = rt_real_corr + toc(t_real_corr)
            if( WRITE_OUTPUT ) call simulated_distrib%write(string('simvol_thres'//trim(real2str(thres))//'_corr'//trim(real2str(t2c))//'.mrc'))
        end function t2c

        ! the B column of write_centers holds valid_corr; the pseudo-atom template width is rendered instead
        subroutine set_pseudo_bfacs( atms )
            type(atoms), intent(inout) :: atms
            integer :: iatm
            do iatm = 1,atms%get_n()
                call atms%set_beta(iatm, B_REF_PSEUDO)
            enddo
        end subroutine set_pseudo_bfacs

    end subroutine binarize_and_find_centers

    subroutine find_centers( self, img_cc, coords, imat )
        class(nanoparticle),            intent(inout) :: self
        type(image_bin),       optional, intent(inout) :: img_cc
        integer,              optional, intent(in)    :: imat(:,:,:)
        real,    allocatable, optional, intent(out)   :: coords(:,:)
        integer, allocatable :: imat_cc_in(:,:,:)
        integer  :: i, ii, jj, kk
        real(dp) :: m(3,self%n_cc), sum_mass(self%n_cc), val
        ! global variables allocation
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        allocate( self%atominfo(self%n_cc) )
        if( present(img_cc) )then
            call img_cc%get_imat(imat_cc_in)
        else if (present(imat)) then
            imat_cc_in = imat
        else
            call self%img_cc%get_imat(imat_cc_in)
        endif
        m        = 0._dp
        sum_mass = 0._dp
        ! voxels of one component are spread over threads, so the per-component sums are reductions
        !$omp parallel do collapse(3) default(shared) private(i,ii,jj,kk,val) reduction(+:m,sum_mass)&
        !$omp schedule(static) proc_bind(close)
        do kk = 1, self%ldim(3)
            do jj = 1, self%ldim(2)
                do ii = 1, self%ldim(1)
                    i = imat_cc_in(ii,jj,kk)
                    if( i >= 1 .and. i <= self%n_cc )then
                        val         = real(self%img_raw%get([ii,jj,kk]),dp)
                        m(:,i)      = m(:,i)      + val * real([ii,jj,kk],dp)
                        sum_mass(i) = sum_mass(i) + val
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
        do i = 1, self%n_cc
            self%atominfo(i)%center = real([self%ldim(1),self%ldim(2),self%ldim(3)]) / 2.
            if( sum_mass(i) > DTINY ) self%atominfo(i)%center = real(m(:,i) / sum_mass(i))
        enddo
        ! saving centers coordinates, optional
        if( present(coords) )then
            allocate(coords(3,self%n_cc))
            do i=1,self%n_cc
                coords(:,i) = self%atominfo(i)%center
            enddo
        endif
    end subroutine find_centers

    subroutine discard_small_ccs( self )
        class(nanoparticle), intent(inout) :: self
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        integer :: cc
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        do cc = 1, self%n_cc
            if( count(imat_cc == cc) < NVOX_THRESH )then ! removes artificial small densities
                where(imat_cc == cc) imat_bin = 0
            endif
        enddo
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        deallocate(imat_bin, imat_cc)
    end subroutine discard_small_ccs

    subroutine split_atoms( self, fname )
        class(nanoparticle),     intent(inout) :: self
        class(string), optional, intent(in)    :: fname
        type(image_bin)       :: img_split_ccs
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:)
        integer, allocatable :: imat(:,:,:), imat_cc(:,:,:), imat_bin(:,:,:), imat_split_ccs(:,:,:)
        integer, parameter   :: RANK_THRESH = 4
        integer :: icc, cnt, cnt_split
        integer :: rank, m(1)
        real    :: new_centers(3,3*self%n_cc) ! will pack it afterwards if it has too many elements
        real    :: pc, radius, split_excl
        write(logfhandle, '(A)') '>>> SPLITTING CONNECTED ATOMS'
        ! squared exclusion radius as the distance tests below use it: voxels**2 * smpd
        if( self%l_species_free )then
            split_excl = (SPLIT_EXCL_NN_FRAC * self%d_nn)**2 / self%smpd
        else
            split_excl = (0.9 * self%theoretical_radius)**2
        endif
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
        write(logfhandle,*) '# atoms split:    ', cnt_split
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
            call img_split_ccs%write(fname)
        else
            call img_split_ccs%write(string('split_ccs.mrc'))
        endif
        call img_split_ccs%kill_bimg
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%update_img_rmat()
        call self%img_bin%find_ccs(self%img_cc)
        call self%update_ncc(self%img_cc)
        call self%find_centers()
        write(logfhandle,*) '# atoms detected: ', self%n_cc
        write(logfhandle, '(A)') '>>> SPLITTING CONNECTED ATOMS, COMPLETED'

    contains

        subroutine split_atom( new_centers, cnt )
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
                            &   (real(k - new_center1(3)))**2) * self%smpd  <=  split_excl) then
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
                if( sum(real(new_center2 - new_center1)**2.) * self%smpd <= split_excl) then
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
                                    &   (real(k - new_center2(3)))**2) * self%smpd <= split_excl )then
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
                if(sum(real(new_center3 - new_center1)**2.) * self%smpd <= split_excl .or. &
                &  sum(real(new_center3 - new_center2)**2.) * self%smpd <= split_excl )then
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
                                    &    (real(k - new_center3(3)))**2) * self%smpd <= split_excl )then
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

    subroutine validate_atoms( self, simatms, l_print )
        class(nanoparticle), intent(inout) :: self
        class(image),        intent(in)    :: simatms
        logical,             intent(in)    :: l_print
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
        enddo
        call calc_stats(self%atominfo(:)%valid_corr, corr_stats)
        if( l_print )then
            write(logfhandle,'(A)') '>>> VALID_CORR (PER-ATOM CORRELATION WITH SIMULATED DENSITY) STATS BELOW'
            write(logfhandle,'(A,F8.4)') 'VALID_CORR Average: ', corr_stats%avg
            write(logfhandle,'(A,F8.4)') 'VALID_CORR Median : ', corr_stats%med
            write(logfhandle,'(A,F8.4)') 'VALID_CORR Sigma  : ', corr_stats%sdev
            write(logfhandle,'(A,F8.4)') 'VALID_CORR Max    : ', corr_stats%maxv
            write(logfhandle,'(A,F8.4)') 'VALID_CORR Min    : ', corr_stats%minv
        endif
    end subroutine validate_atoms

    subroutine discard_atoms( self, l_print )
        class(nanoparticle), intent(inout) :: self
        logical, optional,   intent(in)    :: l_print
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:), cscores(:)
        real,    allocatable :: centers_A(:,:), cendists(:), cendists_sorted(:)
        logical, allocatable :: atom_del_mask(:)
        integer              :: cscore_thres
        type(stats_struct)   :: cscore_stats
        real    :: percen, cendist_thres, foo(3)
        integer :: cc, cn, n_discard, cnt_discard, it_contact_score
        logical :: ll_print
        character(len=5)     :: el_ucase
        character(len=10)    :: crystal_system
        ll_print = .true.
        if( self%l_species_free )then
            it_contact_score = CSCORE_CEIL_NOEL
        else
            el_ucase = uppercase(trim(adjustl(self%element_key)))
            call get_lattice_params(el_ucase, crystal_system, foo)
            select case( crystal_system )
                case('wurtzite', 'zincblende')
                    it_contact_score = 4
                case default !fcc bcc rocksalt
                    it_contact_score = 12
            end select
        endif
        if( present(l_print) ) ll_print = l_print
        if( ll_print ) write(logfhandle, '(A)') '>>> DISCARDING ATOMS'
        ! calculate contact scores
        centers_A = self%atominfo2centers_A()
        allocate(cscores(self%n_cc), source=0)
        call self%contact_scores(centers_A,cscores)
        ! calculate atomic distances from the center of mass of the nanoparticle
        self%NPcen = self%masscen()
        allocate(cendists(self%n_cc), cendists_sorted(self%n_cc), atom_del_mask(self%n_cc))
        do cc = 1, self%n_cc
            cendists(cc)        = euclid(self%atominfo(cc)%center(:), self%NPcen) * self%smpd
            cendists_sorted(cc) = cendists(cc)
        end do
        call hpsort(cendists_sorted)
        ! only subject 15% of the atoms farthest away from the center of mass to deletion
        cendist_thres = cendists_sorted(nint(0.85 * real(self%n_cc)))
        where( cendists > cendist_thres )
            atom_del_mask = .true.
        elsewhere
            atom_del_mask = .false.
        endwhere
        cscore_thres = it_contact_score ! clamped to it_contact_score/2 below if no cn qualifies
        do cn = 1,it_contact_score
            percen = (real(count(cscores >= cn)) / real(self%n_cc)) * 100.
            if( ll_print ) write(logfhandle,*) 'percen atoms with contact score > '//int2str(cn)//':', percen
            if( percen <= 95. )then
                cscore_thres = cn
                exit
            endif
        end do
        if( cscore_thres > it_contact_score/2 ) cscore_thres = it_contact_score/2
        if( ll_print ) write(logfhandle,*) 'contact score threshold: ', cscore_thres
        self%prune_cen       = self%NPcen
        self%prune_rad_thres = cendist_thres
        self%prune_cs_thres  = cscore_thres
        ! get connected components and binary matrices
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        ! discard atoms
        n_discard = 0
        call remove_lowly_contacted(cscore_thres - 1) ! remove these without consideration to size
        if( ll_print ) write(logfhandle, *) '# atoms, discarded based on cs ', n_discard
        cnt_discard = 1
        do while( cnt_discard > 0 )
            call remove_small_and_lowly_contacted(cscore_thres)
            if( ll_print ) write(logfhandle, *) '# atoms, discarded based on sz ', cnt_discard
        end do
        call calc_stats(real(cscores), cscore_stats)
        deallocate(imat_bin, imat_cc)
        if( ll_print )then
            write(logfhandle,'(A)') '>>> CONTACT SCORE STATS BELOW'
            write(logfhandle,'(A,F8.4)') 'Average: ', cscore_stats%avg
            write(logfhandle,'(A,F8.4)') 'Median : ', cscore_stats%med
            write(logfhandle,'(A,F8.4)') 'Sigma  : ', cscore_stats%sdev
            write(logfhandle,'(A,F8.4)') 'Max    : ', cscore_stats%maxv
            write(logfhandle,'(A,F8.4)') 'Min    : ', cscore_stats%minv
            write(logfhandle, *) '# atoms, discarded in total ', n_discard
            write(logfhandle, *) '# atoms, final              ', self%n_cc
            write(logfhandle, '(A)') '>>> DISCARDING ATOMS, COMPLETED'
        endif

    contains

        subroutine remove_lowly_contacted( cthresh )
            integer, intent(in) :: cthresh
            ! discard
            cnt_discard  = 0
            do cc = 1, self%n_cc
                if( cscores(cc) < cthresh .and. atom_del_mask(cc) )then
                    where(imat_cc == cc) imat_bin = 0
                    n_discard   = n_discard   + 1
                    cnt_discard = cnt_discard + 1
                endif
            enddo
            if( cnt_discard > 0 )then
                ! update atoms and contact scores
                call self%img_bin%set_imat(imat_bin)
                call self%img_bin%find_ccs(self%img_cc)
                call self%img_cc%get_nccs(self%n_cc)
                call self%find_centers()
                if( allocated(centers_A) ) deallocate(centers_A)
                if( allocated(cscores)   ) deallocate(cscores)
                allocate(cscores(self%n_cc), source=0)
                centers_A = self%atominfo2centers_A()
                call self%contact_scores(centers_A,cscores)
            endif
        end subroutine remove_lowly_contacted

        subroutine remove_small_and_lowly_contacted( cthresh )
            integer, intent(in) :: cthresh
            real, allocatable   :: radii(:)
            real :: radius_thres
            ! calculate atomic radii
            allocate(radii(self%n_cc), source=0.)
            do cc = 1, self%n_cc
                call self%calc_longest_atm_dist(cc, radii(cc), imat=imat_cc)
            end do
            radius_thres = self%theoretical_radius/2.
            ! discard
            cnt_discard  = 0
            do cc = 1, self%n_cc
                if( atom_del_mask(cc) )then
                    if( radii(cc) < radius_thres .and. cscores(cc) < cthresh )then
                        where(imat_cc == cc) imat_bin = 0
                        n_discard   = n_discard   + 1
                        cnt_discard = cnt_discard + 1
                    endif
                endif
            enddo
            if( cnt_discard > 0 )then
                ! update atoms and contact scores
                call self%img_bin%set_imat(imat_bin)
                call self%img_bin%find_ccs(self%img_cc)
                call self%img_cc%get_nccs(self%n_cc)
                call self%find_centers()
                if( allocated(centers_A) ) deallocate(centers_A)
                if( allocated(cscores)   ) deallocate(cscores)
                allocate(cscores(self%n_cc), source=0)
                centers_A = self%atominfo2centers_A()
                call self%contact_scores(centers_A,cscores)
            endif
        end subroutine remove_small_and_lowly_contacted

    end subroutine discard_atoms

    ! d_NN from the centres of the first binarisation; on the species-free path the atom radius follows from it
    subroutine set_nn_length_scales( self )
        class(nanoparticle), intent(inout) :: self
        real, allocatable :: centers_A(:,:)
        centers_A = self%atominfo2centers_A()
        self%d_nn = est_nn_dist(centers_A)
        write(logfhandle,'(A,F8.4)') 'nearest-neighbour distance d_NN (A): ', self%d_nn
        if( self%l_species_free )then
            self%theoretical_radius = RADIUS_NN_FRAC * self%d_nn
            write(logfhandle,'(A,F8.4)') 'atom radius 0.4 d_NN (A):            ', self%theoretical_radius
        endif
        deallocate(centers_A)
    end subroutine set_nn_length_scales

    ! contact scores of discard_atoms; the neighbour cutoff comes from d_NN on the species-free path
    subroutine contact_scores( self, centers_A, cscores )
        class(nanoparticle), intent(in)    :: self
        real, allocatable,   intent(in)    :: centers_A(:,:)
        integer,             intent(inout) :: cscores(:)
        if( self%l_species_free )then
            call calc_contact_scores(self%element_key, centers_A, cscores, d_nn=self%d_nn)
        else
            call calc_contact_scores(self%element_key, centers_A, cscores)
        endif
    end subroutine contact_scores

    ! Residual recovery and species call on the level-1 atoms (doc/implementation_notes/planned/species_discovery.md,
    ! sections 3.2 to 3.9). Reads the map and the level-1 atoms, fills recovered(:) and the discovery fields of
    ! atominfo(:), and writes the three _species files; the present products are not touched.
    subroutine discover_species( self )
        use simple_nano_species, only: fit_species_mixture, class_separation, enclosed_fraction, calibrate_threshold,&
            &expected_false_count, fit_gauss_width, gauss_filter3D, local_maxima, robust_spread
        class(nanoparticle), intent(inout) :: self
        real,    parameter :: TARGET_FALSE   = 0.2    ! expected noise maxima above k_A in the stage A search volume
        real,    parameter :: CAL_LEVELS(3)  = [2.5, 3.0, 3.5]
        real,    parameter :: REGION_NN_FRAC = 1.5    ! noise region: farther than this from every atom, in d_NN
        real,    parameter :: NMS_NN_FRAC    = 0.35   ! non-maximum suppression and centroid radius, in d_NN
        real,    parameter :: EXCL_A_NN_FRAC = 0.7    ! stage A (and ungated stage B): no atom closer, in d_NN
        real,    parameter :: EXCL_B_NN_FRAC = 0.85   ! gate: no atom closer, in d_NN
        real,    parameter :: GATE_NN_FRAC   = 1.15   ! gate: at least min_nbrs found atoms within, in d_NN
        integer, parameter :: MAX_LEVELS     = 8
        integer, parameter :: NSWEEPS_STAGE1 = 3
        integer, parameter :: NROUNDS_STAGE2 = 2
        integer, parameter :: NSWEEPS_ROUND  = 2
        integer, parameter :: NPHANTOM       = 1000   ! phantom sites for s_A and s_I
        integer, parameter :: NREGION_MIN    = 20000  ! smallest noise region, voxels
        integer, parameter :: NSHELL_MAX     = 5
        integer, parameter :: NSHELL_ATOMS   = 20     ! a radial shell holds at least this many atoms
        real,    parameter :: STRONG_SNR     = 10.    ! strong atoms: amplitude above STRONG_SNR * s_A
        integer, parameter :: NPRIOR_MIN     = 5      ! strong atoms needed for a shell or global width prior
        real,    parameter :: PRIOR_SD_MIN   = 0.1    ! floor of the spread of ln B in the width prior
        real,    parameter :: PRIOR_SD_START = 0.5    ! spread of ln B around B_ref before any fit
        real,    parameter :: TAU_CLASS_MIN  = 0.05   ! floor of the intrinsic class spread, fraction of I_k
        real,    parameter :: LOW_SNR        = 6.5    ! shells below this predicted signal-to-noise are flagged
        real,    parameter :: INNER_FRAC     = 0.3    ! interior atoms for the coordination deficit
        real,    parameter :: CN_CLOSE_PACKED = 12.
        type(atom_stats), allocatable :: tab(:)
        real,    allocatable :: map(:,:,:), hdiff(:,:,:), wmsk(:,:,:), model(:,:,:), backg(:,:,:), resid(:,:,:)
        real,    allocatable :: cmap(:,:,:), zmap(:,:,:), cnoise(:,:,:), emap(:,:,:), omap(:,:,:)
        logical, allocatable :: region(:,:,:), sphere(:,:,:), smask(:,:,:)
        real,    allocatable :: post(:,:), mu(:), var(:), prior_m(:), prior_s(:), rad(:), x(:)
        integer, allocatable :: labels(:), coord(:)
        type(image)  :: img_tmp, simimg
        real         :: s, d, dv, sig_ref, sig2n, s_a, s_i, rob_ratio, k_a, k_b, c_search_a, c_region, z_mu, z_sd
        real,    allocatable :: bic(:)
        logical, allocatable :: adm(:)
        real         :: efalse_a, efalse_b, half_agree, half_corr
        integer      :: ldim(3), n, n1, nreg, nsearch_a, nsearch_b, nadded(2,MAX_LEVELS), npruned, K, stage, level
        integer      :: nnew, i, iround, isweep
        logical      :: l_halves
        if( self%d_nn <= 0. ) THROW_HARD('d_NN was not measured; discover_species')
        write(logfhandle,'(A)') '>>> SPECIES DISCOVERY (DIAGNOSTICS ONLY)'
        s       = self%smpd
        d       = self%d_nn
        dv      = d / s
        sig_ref = sqrt(B_REF_PSEUDO / (8. * PI**2))
        ldim    = self%ldim
        map     = self%img_raw%get_rmat()
        ! mask weights of the soft sphere new applied to the map
        call img_tmp%new(ldim, s)
        allocate(wmsk(ldim(1),ldim(2),ldim(3)), source=1.)
        call img_tmp%set_rmat(wmsk, .false.)
        call img_tmp%mask3D_soft(self%msk_rad, backgr=0.)
        wmsk   = img_tmp%get_rmat()
        sphere = wmsk > 0.
        nsearch_a = count(sphere)
        ! half maps, masked as the map
        l_halves = self%vol_even%is_allocated()
        if( l_halves )then
            call read_half(self%vol_even, emap)
            call read_half(self%vol_odd,  omap)
            hdiff = 0.5 * (emap - omap)
            allocate(cnoise(ldim(1),ldim(2),ldim(3)))
            call gauss_filter3D(hdiff, sig_ref / s, cnoise)
        endif
        allocate(model(ldim(1),ldim(2),ldim(3)), backg(ldim(1),ldim(2),ldim(3)), resid(ldim(1),ldim(2),ldim(3)),&
            &cmap(ldim(1),ldim(2),ldim(3)), zmap(ldim(1),ldim(2),ldim(3)), source=0.)
        allocate(region(ldim(1),ldim(2),ldim(3)), smask(ldim(1),ldim(2),ldim(3)))
        ! level-1 atoms
        n1  = self%n_cc
        n   = n1
        tab = self%atominfo(1:n1)
        do i = 1,n
            tab(i)%amp       = 0.
            tab(i)%bfac      = B_REF_PSEUDO
            tab(i)%det_stage = 0
            tab(i)%det_level = 0
        enddo
        ! provisional stage-1 fit of level 1
        call build_region
        call noise_stats(.true.)
        call fit_atoms_joint(NSWEEPS_STAGE1)
        call build_region
        call noise_stats(.true.)
        ! calibrated thresholds and the residual levels
        call compute_zmap
        call calibrate
        nadded = 0
        do stage = 1,2
            do level = 1,MAX_LEVELS
                call detect_residual_level(stage, level, nnew)
                nadded(stage,level) = nnew
                if( nnew == 0 ) exit
            enddo
        enddo
        call prune_recovered
        ! both fit stages and the species call on the merged, pruned set
        model = 0.
        do i = 1,n
            call render(i, 1.)
        enddo
        call fit_atoms_joint(NSWEEPS_STAGE1)
        do i = 1,n
            tab(i)%bfac_stage1 = tab(i)%bfac
        enddo
        call build_region
        call noise_stats(.true.)
        call calc_aperture_int
        call assign_species
        do iround = 1,NROUNDS_STAGE2
            do isweep = 1,NSWEEPS_ROUND
                call tied_sweep
            enddo
            call calc_aperture_int
            call assign_species
        enddo
        call recovered_valid_corr
        call geometry
        if( l_halves ) call halfmap_agreement
        call write_species_report
        call write_radial_profiles
        call write_species_table
        call write_species_pdb
        ! tables: discovery fields of the level-1 atoms, the recovered atoms on their own
        self%atominfo(1:n1) = tab(1:n1)
        if( allocated(self%recovered) ) deallocate(self%recovered)
        self%recovered    = tab(n1+1:n)
        self%species_post = post
        call img_tmp%kill
        call simimg%kill
        write(logfhandle,'(A,I0,A,I0,A,I0)') 'level-1 atoms: ', n1, ', recovered atoms: ', n - n1, ', classes: ', K
        write(logfhandle,'(A)') '>>> SPECIES DISCOVERY, COMPLETED'

    contains

        subroutine read_half( fname, arr )
            class(string),     intent(in)  :: fname
            real, allocatable, intent(out) :: arr(:,:,:)
            integer :: ldim_h(3), nptcls
            call find_ldim_nptcls(fname, ldim_h, nptcls)
            if( any(ldim_h /= ldim) ) THROW_HARD('half maps and map differ in size; discover_species')
            call img_tmp%new(ldim, s)
            call img_tmp%read(fname)
            call img_tmp%mask3D_soft(self%msk_rad)
            arr = img_tmp%get_rmat()
        end subroutine read_half

        ! voxel window of a ball of radius rv voxels around c (voxel coordinates)
        subroutine ball_window( c, rv, lo, hi )
            real,    intent(in)  :: c(3), rv
            integer, intent(out) :: lo(3), hi(3)
            lo = max(1, floor(c - rv))
            hi = min(ldim, ceiling(c + rv))
        end subroutine ball_window

        ! set arr to val within rv voxels of c
        subroutine mark_ball( c, rv, arr, val )
            real,    intent(in)    :: c(3), rv
            logical, intent(inout) :: arr(:,:,:)
            logical, intent(in)    :: val
            integer :: lo(3), hi(3), ii, jj, kk
            call ball_window(c, rv, lo, hi)
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        if( sum((real([ii,jj,kk]) - c)**2) <= rv * rv ) arr(ii,jj,kk) = val
                    enddo
                enddo
            enddo
        end subroutine mark_ball

        ! render cutoff of atom i in A: four standard deviations, and always beyond the fit sphere
        real function render_cutoff( i )
            integer, intent(in) :: i
            render_cutoff = max(4. * sqrt(tab(i)%bfac / (8. * PI**2)), 0.5 * d + s)
        end function render_cutoff

        ! add sgn times the fitted density of atom i to the model map
        subroutine render( i, sgn )
            integer, intent(in) :: i
            real,    intent(in) :: sgn
            integer :: lo(3), hi(3), ii, jj, kk
            real    :: rc, r2
            if( tab(i)%amp <= 0. ) return
            rc = render_cutoff(i)
            call ball_window(tab(i)%center, rc / s, lo, hi)
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        r2 = sum((real([ii,jj,kk]) - tab(i)%center)**2) * s * s
                        if( r2 > rc * rc ) cycle
                        model(ii,jj,kk) = model(ii,jj,kk) + sgn * tab(i)%amp * exp(-4. * PI**2 * r2 / tab(i)%bfac)
                    enddo
                enddo
            enddo
        end subroutine render

        ! samples of src - background - model within rad A of c, with atom iself's own density added back
        subroutine gather( src, c, rad_a, iself, y, r2 )
            real,              intent(in)  :: src(:,:,:), c(3), rad_a
            integer,           intent(in)  :: iself
            real, allocatable, intent(out) :: y(:), r2(:)
            real    :: ybuf((2 * (ceiling(rad_a / s) + 1) + 1)**3), rbuf(size(ybuf)), dd
            integer :: lo(3), hi(3), ii, jj, kk, m
            call ball_window(c, rad_a / s, lo, hi)
            m = 0
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        dd = sum((real([ii,jj,kk]) - c)**2) * s * s
                        if( dd > rad_a * rad_a ) cycle
                        m = m + 1
                        rbuf(m) = dd
                        ybuf(m) = src(ii,jj,kk) - backg(ii,jj,kk) - model(ii,jj,kk)
                        if( iself > 0 )then
                            if( tab(iself)%amp > 0. ) ybuf(m) = ybuf(m) + tab(iself)%amp * exp(-4. * PI**2 * dd / tab(iself)%bfac)
                        endif
                    enddo
                enddo
            enddo
            y  = ybuf(:m)
            r2 = rbuf(:m)
        end subroutine gather

        ! voxels of full mask weight farther than REGION_NN_FRAC d_NN from every atom
        subroutine build_region
            character(len=:), allocatable :: msg
            region = wmsk >= 0.9999
            do i = 1,n
                call mark_ball(tab(i)%center, REGION_NN_FRAC * dv, region, .false.)
            enddo
            nreg = count(region)
            if( nreg < NREGION_MIN )then
                msg = 'noise region of '//int2str(nreg)//' voxels is too small; give a larger mskdiam'
                THROW_HARD(msg)
            endif
        end subroutine build_region

        ! sigma_n from the half-map difference or the map over the region; with l_full also the spread ratio,
        ! and s_A and s_I from phantom sites of the region
        subroutine noise_stats( l_full )
            logical, intent(in) :: l_full
            real, allocatable :: vals(:), amps(:), ints(:), y(:), r2(:), g(:)
            integer :: stride, cnt, nsite, ii, jj, kk
            real    :: c(3)
            if( l_halves )then
                vals = pack(hdiff, region)
            else
                vals = pack(map, region)
            endif
            sig2n = sdev(vals)**2
            if( sig2n <= 0. )then
                THROW_HARD('no noise in the noise region: discover_species needs a noisy map')
            endif
            if( .not. l_full ) return
            rob_ratio = robust_spread(vals) / sqrt(sig2n)
            stride = max(1, nreg / NPHANTOM)
            allocate(amps(nreg / stride + 1), ints(nreg / stride + 1))
            cnt   = 0
            nsite = 0
            do kk = 1,ldim(3)
                do jj = 1,ldim(2)
                    do ii = 1,ldim(1)
                        if( .not. region(ii,jj,kk) ) cycle
                        cnt = cnt + 1
                        if( mod(cnt - 1, stride) /= 0 ) cycle
                        c = real([ii,jj,kk])
                        call gather(map, c, 0.5 * d, 0, y, r2)
                        g = exp(-r2 / (2. * sig_ref**2))
                        nsite = nsite + 1
                        amps(nsite) = sum(y * g) / sum(g * g)
                        ints(nsite) = s**3 * sum(y) / enclosed_fraction(0.5 * d / sig_ref)
                    enddo
                enddo
            enddo
            s_a = sdev(amps(:nsite))
            s_i = sdev(ints(:nsite))
        end subroutine noise_stats

        real function sdev( v )
            real, intent(in) :: v(:)
            real :: m
            m    = sum(v) / real(size(v))
            sdev = sqrt(sum((v - m)**2) / real(max(size(v) - 1, 1)))
        end function sdev

        ! width priors of stage 1: ln B of the strong atoms within d_NN in radius, else of all strong atoms
        subroutine width_priors
            real, allocatable :: lnb(:), r(:)
            logical, allocatable :: strong(:), sel(:)
            real    :: cen(3)
            integer :: j
            if( allocated(prior_m) ) deallocate(prior_m, prior_s)
            allocate(prior_m(n), source=log(B_REF_PSEUDO))
            allocate(prior_s(n), source=PRIOR_SD_START)
            cen = 0.
            do j = 1,n
                cen = cen + tab(j)%center
            enddo
            cen = cen / real(n)
            allocate(r(n), lnb(n), strong(n), sel(n))
            do j = 1,n
                r(j)      = sqrt(sum((tab(j)%center - cen)**2)) * s
                lnb(j)    = log(tab(j)%bfac)
                strong(j) = tab(j)%amp > STRONG_SNR * s_a
            enddo
            if( count(strong) < NPRIOR_MIN ) return
            do j = 1,n
                sel = strong .and. abs(r - r(j)) <= d
                if( count(sel) < NPRIOR_MIN ) sel = strong
                prior_m(j) = median(pack(lnb, sel))
                prior_s(j) = max(robust_spread(pack(lnb, sel)), PRIOR_SD_MIN)
            enddo
        end subroutine width_priors

        ! refit atom i on its own residual within d_NN / 2, free (stage 1) or tied to its class (stage 2)
        subroutine fit_one( i, i_class, tau_class )
            integer,        intent(in) :: i
            real, optional, intent(in) :: i_class, tau_class
            real, allocatable :: y(:), r2(:)
            real :: a, b
            call gather(map, tab(i)%center, 0.5 * d, i, y, r2)
            if( present(i_class) )then
                call fit_gauss_width(y, r2, sig2n, 0., 1., a, b, i_class=i_class, tau_class=tau_class)
            else
                call fit_gauss_width(y, r2, sig2n, prior_m(i), prior_s(i), a, b)
            endif
            call render(i, -1.)
            tab(i)%amp  = a
            tab(i)%bfac = b
            call render(i, 1.)
        end subroutine fit_one

        ! background: the residual smoothed with a Gaussian of standard deviation d_NN
        subroutine update_background
            resid = map - model
            call gauss_filter3D(resid, dv, backg)
        end subroutine update_background

        ! stage 1: Gauss-Seidel sweeps of free-amplitude fits, the background updated after each
        subroutine fit_atoms_joint( nsweeps )
            integer, intent(in) :: nsweeps
            integer :: isw, j
            do isw = 1,nsweeps
                call width_priors
                do j = 1,n
                    call fit_one(j)
                enddo
                call update_background
            enddo
        end subroutine fit_atoms_joint

        ! stage 2: every atom tied to the intensity of its class
        subroutine tied_sweep
            real    :: tau
            integer :: j, kc
            do j = 1,n
                kc = tab(j)%species
                if( mu(kc) <= 0. ) cycle
                tau = max(sqrt(max(var(kc) - s_i**2, 0.)), TAU_CLASS_MIN * mu(kc))
                call fit_one(j, i_class=mu(kc), tau_class=tau)
            enddo
            call update_background
        end subroutine tied_sweep

        ! standardised matched-filter map of the residual
        subroutine compute_zmap
            real, allocatable :: vals(:)
            resid = map - backg - model
            call gauss_filter3D(resid, sig_ref / s, cmap)
            if( l_halves )then
                vals = pack(cnoise, region)
            else
                vals = pack(cmap, region)
            endif
            z_mu = sum(vals) / real(size(vals))
            z_sd = sdev(vals)
            zmap = (cmap - z_mu) / z_sd
        end subroutine compute_zmap

        ! thresholds from the maxima of the noise field in the region (section 3.2)
        subroutine calibrate
            real,    allocatable :: vals(:), nz(:,:,:)
            integer, allocatable :: ijk(:,:)
            real    :: counts(size(CAL_LEVELS)), ratio
            integer :: l
            if( l_halves )then
                nz = (cnoise - z_mu) / z_sd
                call local_maxima(nz, region, CAL_LEVELS(1), 0, ijk, vals)
            else
                call local_maxima(zmap, region, CAL_LEVELS(1), 0, ijk, vals)
            endif
            do l = 1,size(CAL_LEVELS)
                counts(l) = real(count(vals > CAL_LEVELS(l)))
            enddo
            ratio = real(nsearch_a) / real(nreg)
            call calibrate_threshold(CAL_LEVELS, counts, ratio, TARGET_FALSE, k_a, c_search_a)
            c_region = c_search_a / ratio
            k_b      = k_a - 1.
            call stage_b_mask
            nsearch_b = count(smask)
            efalse_a  = expected_false_count(c_search_a, k_a)
            efalse_b  = expected_false_count(c_region * real(nsearch_b) / real(nreg), k_b)
            write(logfhandle,'(A,3F8.1)')  'region maxima above 2.5, 3.0, 3.5: ', counts
            write(logfhandle,'(A,2F8.3)')  'calibrated thresholds k_A, k_B:     ', k_a, k_b
        end subroutine calibrate

        ! stage B searches within GATE_NN_FRAC d_NN of the atoms when gated, else the whole sphere
        subroutine stage_b_mask
            integer :: j
            if( self%min_nbrs > 0 )then
                smask = .false.
                do j = 1,n
                    call mark_ball(tab(j)%center, GATE_NN_FRAC * dv, smask, .true.)
                enddo
                smask = smask .and. sphere
            else
                smask = sphere
            endif
        end subroutine stage_b_mask

        ! one residual level: candidates are local maxima of the z map above the stage threshold, positioned at
        ! the centroid of the positive residual, accepted by the exclusion distance or the neighbour gate (3.3)
        subroutine detect_residual_level( stage, level, nnew )
            integer, intent(in)  :: stage, level
            integer, intent(out) :: nnew
            real,    allocatable :: vals(:), considered(:,:), dists(:)
            integer, allocatable :: ijk(:,:)
            type(atom_stats) :: atm
            real    :: thres, pos(3)
            integer :: ic, nc, nbefore, j
            logical :: l_gated, l_ok
            call build_region
            call noise_stats(.false.)
            call compute_zmap
            l_gated = stage == 2 .and. self%min_nbrs > 0
            if( stage == 1 )then
                thres = k_a
                smask = sphere
            else
                thres = k_b
                call stage_b_mask
            endif
            ! NVOX_THRESH of the 27 voxels above the threshold, the maximum itself included
            call local_maxima(zmap, smask, thres, NVOX_THRESH - 1, ijk, vals)
            allocate(considered(3,max(size(vals),1)))
            nc      = 0
            nbefore = n
            do ic = 1,size(vals)
                ! non-maximum suppression between the maxima of this level
                if( nc > 0 )then
                    if( any(sum((considered(:,:nc) - spread(real(ijk(:,ic)), 2, nc))**2, dim=1) < (NMS_NN_FRAC * dv)**2) ) cycle
                endif
                nc = nc + 1
                considered(:,nc) = real(ijk(:,ic))
                pos = centroid(ijk(:,ic))
                allocate(dists(n))
                do j = 1,n
                    dists(j) = sqrt(sum((tab(j)%center - pos)**2)) * s
                enddo
                if( l_gated )then
                    ! the gate counts the atoms accepted before this level: one shell per level
                    l_ok = minval(dists) > EXCL_B_NN_FRAC * d .and. count(dists(:nbefore) <= GATE_NN_FRAC * d) >= self%min_nbrs
                else
                    l_ok = minval(dists) > EXCL_A_NN_FRAC * d
                endif
                deallocate(dists)
                if( .not. l_ok ) cycle
                atm = atom_stats()
                atm%center    = pos
                atm%amp       = 0.
                atm%bfac      = B_REF_PSEUDO
                atm%det_z     = vals(ic)
                atm%det_stage = stage
                atm%det_level = level
                atm%gate_used = l_gated
                tab = [tab(:n), atm]
                n   = n + 1
            enddo
            nnew = n - nbefore
            if( nnew > 0 )then
                call width_priors
                do j = nbefore+1,n
                    call fit_one(j)
                enddo
            endif
            write(logfhandle,'(A,I1,A,I1,A,I0)') 'residual stage ', stage, ', level ', level, ': atoms added ', nnew
        end subroutine detect_residual_level

        ! centroid of the positive residual within NMS_NN_FRAC d_NN of a voxel
        function centroid( ijk ) result( pos )
            integer, intent(in) :: ijk(3)
            real    :: pos(3), w, wsum, c(3)
            integer :: lo(3), hi(3), ii, jj, kk
            c    = real(ijk)
            pos  = 0.
            wsum = 0.
            call ball_window(c, NMS_NN_FRAC * dv, lo, hi)
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        if( sum((real([ii,jj,kk]) - c)**2) > (NMS_NN_FRAC * dv)**2 ) cycle
                        w = max(resid(ii,jj,kk), 0.)
                        pos  = pos + w * real([ii,jj,kk])
                        wsum = wsum + w
                    enddo
                enddo
            enddo
            if( wsum > 0. )then
                pos = pos / wsum
            else
                pos = c
            endif
        end function centroid

        ! the pruning of discard_atoms applied to the recovered atoms (ruling of 2026-10-09, second): in the outer zone
        ! of the level-1 set, a recovered atom with fewer level-1 atoms within the contact cutoff than the threshold
        ! discard_atoms derived goes; recovered atoms never vouch for each other, atoms inside the zone stay
        subroutine prune_recovered
            logical, allocatable :: keep(:)
            real    :: rmax_c, cendist
            integer :: j, l, ncont
            npruned = 0
            if( n == n1 ) return
            if( self%l_species_free )then
                rmax_c = find_rMax(self%element_key, self%d_nn)
            else
                rmax_c = find_rMax(self%element_key)
            endif
            allocate(keep(n), source=.true.)
            do j = n1+1,n
                cendist = euclid(tab(j)%center, self%prune_cen) * s
                if( cendist <= self%prune_rad_thres ) cycle
                ncont = 0
                do l = 1,n1
                    if( euclid(tab(j)%center, tab(l)%center) * s < rmax_c ) ncont = ncont + 1
                enddo
                if( ncont < self%prune_cs_thres ) keep(j) = .false.
            enddo
            npruned = count(.not. keep)
            if( npruned > 0 )then
                tab = pack(tab(:n), keep)
                n   = size(tab)
            endif
            write(logfhandle,'(A,I0)') 'recovered atoms pruned by contact score: ', npruned
        end subroutine prune_recovered

        subroutine merged_centers_a( centers_a )
            real, allocatable, intent(out) :: centers_a(:,:)
            integer :: j
            allocate(centers_a(3,n))
            do j = 1,n
                centers_a(:,j) = (tab(j)%center - 1.) * s
            enddo
        end subroutine merged_centers_a

        ! aperture intensities (3.5) in the units of the continuous integral
        subroutine calc_aperture_int
            real, allocatable :: y(:), r2(:)
            integer :: j
            do j = 1,n
                call gather(map, tab(j)%center, 0.5 * d, j, y, r2)
                tab(j)%aper_int = s**3 * sum(y) / enclosed_fraction(0.5 * d / sqrt(tab(j)%bfac / (8. * PI**2)))
                tab(j)%aper_snr = tab(j)%aper_int / s_i
            enddo
        end subroutine calc_aperture_int

        ! mixture of the aperture intensities (3.6)
        subroutine assign_species
            integer :: j
            if( allocated(labels) ) deallocate(labels)
            allocate(labels(n))
            x = tab(:n)%aper_int
            call fit_species_mixture(x, s_i, self%nspecies, K, labels, post, mu, var, bic, adm)
            do j = 1,n
                tab(j)%species = labels(j)
            enddo
        end subroutine assign_species

        ! valid_corr of the recovered atoms against the simulation with class intensities and own widths
        subroutine recovered_valid_corr
            type(atoms)       :: atms
            real, allocatable :: pix1(:), pix2(:)
            real    :: maxrad, cutoff
            integer :: j, winsz, npix_in, nout1, nout2, ijk(3)
            if( n == n1 ) return
            call atms%new(n, dummy=.true.)
            cutoff = 8. * s
            do j = 1,n
                call atms%set_element(j, 'X1')
                call atms%set_coord(j, (tab(j)%center - 1.) * s)
                call atms%set_occupancy(j, max(mu(tab(j)%species), 0.))
                call atms%set_beta(j, tab(j)%bfac)
                cutoff = max(cutoff, 4. * sqrt(tab(j)%bfac / (8. * PI**2)))
            enddo
            call simimg%new(ldim, s)
            call atms%convolve(simimg, cutoff)
            call atms%kill
            maxrad  = 1.5 * RADIUS_NN_FRAC * dv
            winsz   = ceiling(maxrad)
            npix_in = (2 * winsz + 1)**3
            allocate(pix1(npix_in), pix2(npix_in), source=0.)
            do j = n1+1,n
                ijk = nint(tab(j)%center)
                call self%img_raw%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, nout1, pix1)
                call simimg%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, nout2, pix2)
                tab(j)%valid_corr = pearsn_serial(pix1(:nout1), pix2(:nout2))
            enddo
        end subroutine recovered_valid_corr

        ! radius from the centre of the atom positions and coordination within the neighbour cutoff
        subroutine geometry
            real, allocatable :: centers_a(:,:)
            real    :: cen(3)
            integer :: j
            call merged_centers_a(centers_a)
            if( allocated(coord) ) deallocate(coord)
            allocate(coord(n), source=0)
            call self%contact_scores(centers_a, coord)
            cen = sum(centers_a, dim=2) / real(n)
            rad = sqrt(sum((centers_a - spread(cen, 2, n))**2, dim=1))
            do j = 1,n
                tab(j)%cendist = rad(j)
            enddo
        end subroutine geometry

        ! aperture intensities in each half map, classified with the final mixture
        subroutine halfmap_agreement
            real, allocatable :: y(:), r2(:), ie(:), io(:)
            integer :: j, nagree
            real    :: frac
            allocate(ie(n), io(n))
            do j = 1,n
                frac = enclosed_fraction(0.5 * d / sqrt(tab(j)%bfac / (8. * PI**2)))
                call gather(emap, tab(j)%center, 0.5 * d, j, y, r2)
                ie(j) = s**3 * sum(y) / frac
                call gather(omap, tab(j)%center, 0.5 * d, j, y, r2)
                io(j) = s**3 * sum(y) / frac
            enddo
            nagree = 0
            do j = 1,n
                if( classify(ie(j)) == classify(io(j)) ) nagree = nagree + 1
            enddo
            half_agree = real(nagree) / real(n)
            half_corr  = pearsn_serial(ie, io)
        end subroutine halfmap_agreement

        integer function classify( xi )
            real, intent(in) :: xi
            real    :: lp(K), w
            integer :: kc
            do kc = 1,K
                w = real(count(labels == kc)) / real(n)
                if( w > 0. )then
                    lp(kc) = log(w) - 0.5 * log(var(kc)) - 0.5 * (xi - mu(kc))**2 / var(kc)
                else
                    lp(kc) = -huge(1.)
                endif
            enddo
            classify = maxloc(lp, dim=1)
        end function classify

        ! best possible detection signal-to-noise of an atom of intensity ik and width sigma (section 7): its amplitude
        ! over s_A, scaled from the template width to its own as (sigma / sig_ref)**1.5, so it goes as ik sigma**(-1.5)
        real function pred_snr( ik, sigma )
            real, intent(in) :: ik, sigma
            pred_snr = ik / (2. * PI * sigma**2)**1.5 * (sigma / sig_ref)**1.5 / s_a
        end function pred_snr

        ! equal-count radial shells of class kc: first and last index into the radius order
        subroutine class_shells( kc, order, nsh, bnd )
            integer,              intent(in)  :: kc
            integer, allocatable, intent(out) :: order(:), bnd(:)
            integer,              intent(out) :: nsh
            real,    allocatable :: rk(:)
            integer :: j, nk
            order = pack([(j, j=1,n)], tab(:n)%species == kc)
            nk    = size(order)
            nsh   = max(1, min(NSHELL_MAX, nk / NSHELL_ATOMS))
            if( nk > 0 )then
                rk = rad(order)
                call hpsort(rk, order)
            endif
            allocate(bnd(nsh+1))
            do j = 1,nsh+1
                bnd(j) = ((j-1) * nk) / nsh
            enddo
        end subroutine class_shells

        subroutine write_radial_profiles
            integer, allocatable :: order(:), bnd(:), sel(:)
            real,    allocatable :: sig1(:), sig2(:), ints(:)
            integer :: funit, kc, ish, nsh, m
            real    :: rms2, rms1
            call fopen(funit, file=self%fbody//'_species_radial.csv', status='REPLACE', action='WRITE')
            write(funit,'(A)') 'CLASS,SHELL,NATOMS,MEAN_RADIUS,RMS_SIGMA,SE_SIGMA,B,RMS_SIGMA_STAGE1,MEAN_INT,SE_INT,'//&
                &'MEAN_COORD,PRED_SNR,LOW_SNR'
            do kc = 1,K
                call class_shells(kc, order, nsh, bnd)
                do ish = 1,nsh
                    m = bnd(ish+1) - bnd(ish)
                    if( m < 1 ) cycle
                    sel  = order(bnd(ish)+1:bnd(ish+1))
                    sig2 = sqrt(tab(sel)%bfac / (8. * PI**2))
                    sig1 = sqrt(tab(sel)%bfac_stage1 / (8. * PI**2))
                    ints = tab(sel)%aper_int
                    rms2 = sqrt(sum(sig2**2) / real(m))
                    rms1 = sqrt(sum(sig1**2) / real(m))
                    write(funit,'(I0,",",I0,",",I0,9(",",ES16.8),",",I0)') kc, ish, m, sum(rad(sel)) / real(m), rms2,&
                        &sdev(sig2) / sqrt(real(m)), 8. * PI**2 * rms2**2, rms1, sum(ints) / real(m), sdev(ints) / sqrt(real(m)),&
                        &real(sum(coord(sel))) / real(m), pred_snr(mu(kc), rms2), merge(1, 0, pred_snr(mu(kc), rms2) < LOW_SNR)
                enddo
            enddo
            call fclose(funit)
        end subroutine write_radial_profiles

        ! one posterior column per class, POST1 to POSTK
        subroutine write_species_table
            type(string) :: header, fmt
            integer      :: funit, j, origin, kc
            header = 'INDEX,RECOVERED,X,Y,Z,RADIUS,COORD,STAGE,LEVEL,DET_Z,AMP,B_STAGE1,B_STAGE2,APER_INT,APER_SNR,CLASS'
            do kc = 1,K
                header = header//',POST'//int2str(kc)
            enddo
            header = header//',VALID_CORR,GATE_USED'
            fmt    = '(I0,",",I0,5(",",ES16.8),2(",",I0),6(",",ES16.8),",",I0,'//int2str(K+1)//'(",",ES16.8),",",I0)'
            call fopen(funit, file=self%fbody//'_species.csv', status='REPLACE', action='WRITE')
            write(funit,'(A)') header%to_char()
            do j = 1,n
                origin = merge(1, 0, j > n1)
                write(funit,fmt%to_char()) j, origin,&
                    &(tab(j)%center - 1.) * s, rad(j), real(coord(j)), tab(j)%det_stage, tab(j)%det_level, tab(j)%det_z,&
                    &tab(j)%amp, tab(j)%bfac_stage1, tab(j)%bfac, tab(j)%aper_int, tab(j)%aper_snr, tab(j)%species,&
                    &post(j,:), tab(j)%valid_corr, merge(1, 0, tab(j)%gate_used)
            enddo
            call fclose(funit)
        end subroutine write_species_table

        ! every atom with its class: element column and chain by class (X1, or the first named species, and chain A
        ! the strongest), occupancy the class intensity over the strongest class's, B column the stage-2 B factor;
        ! a viewer selects a class by chain
        subroutine write_species_pdb
            type(atoms)      :: atms
            character(len=2) :: el
            integer          :: j, kc
            real             :: occ
            if( K > 36 ) THROW_HARD('more classes than PDB chains A-Z, 0-9; write_species_pdb')
            call atms%new(n, dummy=.true.)
            do j = 1,n
                kc  = tab(j)%species
                occ = 1.
                if( mu(1) > TINY ) occ = max(mu(kc), 0.) / mu(1)
                if( allocated(self%species) )then
                    el = self%species(kc)
                else
                    el = 'X'//int2str(kc)
                endif
                call atms%set_name(     j, el//'  ')
                call atms%set_element(  j, el)
                if( kc <= 26 )then
                    call atms%set_chain(j, achar(iachar('A') + kc - 1))
                else
                    call atms%set_chain(j, achar(iachar('0') + kc - 27))
                endif
                call atms%set_coord(    j, (tab(j)%center - 1.) * s)
                call atms%set_num(      j, j)
                call atms%set_resnum(   j, j)
                call atms%set_occupancy(j, occ)
                call atms%set_beta(     j, tab(j)%bfac)
            enddo
            call atms%writepdb(self%fbody//'_species.pdb')
            call atms%kill
            write(logfhandle,'(A)') 'output, atoms with their class:   '//self%fbody%to_char()//'_species.pdb'
        end subroutine write_species_pdb

        subroutine write_species_report
            real, allocatable :: lnb(:), r(:)
            logical, allocatable :: strong(:)
            integer, allocatable :: order(:), bnd(:)
            real, parameter :: SURF_CONFINED = 0.8   ! a class with more of its atoms in the outer zone is flagged
            integer :: funit, kc, l, ninner, nsh, ish, nlow, nin
            real    :: bmed, cn_inner, misclass, rms2, cn_out, frac_out
            call fopen(funit, file=self%fbody//'_species.txt', status='REPLACE', action='WRITE')
            write(funit,'(A)') '# species discovery of detect_atoms (discover_species=yes): diagnostics, not used by any product'
            call write_kv(funit, 'K', int2str(K))
            call write_kv(funit, 'nspecies_fixed', int2str(self%nspecies))
            do kc = 1,K
                call write_kv(funit, 'class_intensity_'//int2str(kc), rstr(mu(kc)))
                call write_kv(funit, 'class_ratio_'//int2str(kc),     rstr(mu(kc) / mu(1)))
                call write_kv(funit, 'class_fraction_'//int2str(kc),  rstr(real(count(labels == kc)) / real(n)))
                call write_kv(funit, 'class_atoms_'//int2str(kc),     int2str(count(labels == kc)))
                call write_kv(funit, 'class_sd_'//int2str(kc),        rstr(sqrt(var(kc))))
            enddo
            do kc = 1,K-1
                call write_kv(funit, 'separation_D_'//int2str(kc)//'_'//int2str(kc+1), rstr(class_separation(mu(kc), var(kc), mu(kc+1), var(kc+1))))
            enddo
            do l = 1,size(bic)
                call write_kv(funit, 'bic_K'//int2str(l),        rstr(bic(l)))
                call write_kv(funit, 'admissible_K'//int2str(l), merge('yes', 'no ', adm(l)))
            enddo
            misclass = sum(1. - maxval(post, dim=2)) / real(n)
            call write_kv(funit, 'expected_misclassification', rstr(misclass))
            call write_kv(funit, 'd_NN', rstr(d))
            call write_kv(funit, 'B_ref', rstr(B_REF_PSEUDO))
            allocate(strong(n))
            strong = tab(:n)%amp > STRONG_SNR * s_a
            if( count(strong) > 0 )then
                lnb  = pack(tab(:n)%bfac_stage1, strong)
                bmed = median(lnb)
                call write_kv(funit, 'B_stage1_median_strong', rstr(bmed))
                call write_kv(funit, 'B_stage1_over_B_ref', rstr(bmed / B_REF_PSEUDO))
            endif
            call write_kv(funit, 'noise_source', merge('half_maps', 'region   ', l_halves))
            call write_kv(funit, 'sigma_n', rstr(sqrt(sig2n)))
            call write_kv(funit, 's_A', rstr(s_a))
            call write_kv(funit, 's_I', rstr(s_i))
            call write_kv(funit, 'robust_over_plain_spread', rstr(rob_ratio))
            call write_kv(funit, 'region_voxels', int2str(nreg))
            call write_kv(funit, 'k_A', rstr(k_a))
            call write_kv(funit, 'k_B', rstr(k_b))
            call write_kv(funit, 'search_voxels_A', int2str(nsearch_a))
            call write_kv(funit, 'search_voxels_B', int2str(nsearch_b))
            call write_kv(funit, 'expected_false_A', rstr(efalse_a))
            call write_kv(funit, 'expected_false_B', rstr(efalse_b))
            call write_kv(funit, 'min_nbrs', int2str(self%min_nbrs))
            call write_kv(funit, 'level1_atoms', int2str(n1))
            do l = 1,MAX_LEVELS
                call write_kv(funit, 'added_A_level_'//int2str(l), int2str(nadded(1,l)))
            enddo
            do l = 1,MAX_LEVELS
                call write_kv(funit, 'added_B_level_'//int2str(l), int2str(nadded(2,l)))
            enddo
            call write_kv(funit, 'prune_zone_radius', rstr(self%prune_rad_thres))
            call write_kv(funit, 'prune_contact_threshold', int2str(self%prune_cs_thres))
            call write_kv(funit, 'pruned', int2str(npruned))
            call write_kv(funit, 'recovered_atoms', int2str(n - n1))
            call write_kv(funit, 'total_atoms', int2str(n))
            ! coordination of the inner atoms: a deficit flags undetected sites (section 7)
            r = rad
            order = [(l, l=1,n)]
            call hpsort(r, order)
            ninner   = max(1, nint(INNER_FRAC * real(n)))
            cn_inner = real(sum(coord(order(:ninner)))) / real(ninner)
            call write_kv(funit, 'interior_coordination', rstr(cn_inner))
            call write_kv(funit, 'interior_coordination_deficit', rstr(CN_CLOSE_PACKED - cn_inner))
            nlow = 0
            do kc = 1,K
                call class_shells(kc, order, nsh, bnd)
                do ish = 1,nsh
                    if( bnd(ish+1) <= bnd(ish) ) cycle
                    rms2 = sqrt(sum(tab(order(bnd(ish)+1:bnd(ish+1)))%bfac) / real(bnd(ish+1) - bnd(ish)) / (8. * PI**2))
                    if( pred_snr(mu(kc), rms2) < LOW_SNR ) nlow = nlow + 1
                enddo
            enddo
            call write_kv(funit, 'low_snr_shells', int2str(nlow))
            ! a class confined to the outer zone of the pruning may be partial occupancy rather than a species
            do kc = 1,K
                nin = 0
                cn_out = 0.
                do l = 1,n
                    if( tab(l)%species /= kc ) cycle
                    if( euclid(tab(l)%center, self%prune_cen) * s > self%prune_rad_thres )then
                        nin    = nin + 1
                        cn_out = cn_out + real(coord(l))
                    endif
                enddo
                frac_out = real(nin) / real(max(count(labels == kc), 1))
                call write_kv(funit, 'class_outer_zone_fraction_'//int2str(kc), rstr(frac_out))
                call write_kv(funit, 'class_outer_zone_mean_coord_'//int2str(kc), rstr(cn_out / real(max(nin, 1))))
                call write_kv(funit, 'class_surface_confined_'//int2str(kc), merge('yes', 'no ', frac_out > SURF_CONFINED))
            enddo
            if( l_halves )then
                call write_kv(funit, 'halfmap_label_agreement', rstr(half_agree))
                call write_kv(funit, 'halfmap_intensity_corr', rstr(half_corr))
            endif
            call fclose(funit)
        end subroutine write_species_report

        subroutine write_kv( funit, key, val )
            integer,          intent(in) :: funit
            character(len=*), intent(in) :: key, val
            write(funit,'(A,T34,A)') key, '= '//trim(val)
        end subroutine write_kv

        function rstr( v ) result( str )
            real, intent(in) :: v
            character(len=16) :: str
            write(str,'(ES16.8)') v
            str = adjustl(str)
        end function rstr

    end subroutine discover_species

    subroutine fillin_atominfo( self, a0, imat )
        class(nanoparticle),        intent(inout) :: self
        real,             optional, intent(in)    :: a0(3) ! lattice parameters
        integer,          optional, intent(in)    :: imat(:,:,:)
        type(image)          :: simatms, fit_isotropic, fit_anisotropic
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:), strain_array(:,:)
        real,    pointer     :: rmat_raw(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        logical, allocatable :: cc_mask(:)
        real                 :: tmp_diam, a(3)
        integer              :: i, cc, cn, max_size
        character(len=*), parameter :: fn_fit_isotropic="fit_isotropic.mrc", fn_fit_anisotropic="fit_anisotropic.mrc"
        if( self%l_species_free ) THROW_HARD('fillin_atominfo needs an element')
        write(logfhandle, '(A)') '>>> EXTRACTING ATOM STATISTICS'
        write(logfhandle, '(A)') '---Dev Note: ADP and Max Neighboring Displacements Under Testing---'
        ! calc cn and cn_gen
        centers_A = self%atominfo2centers_A()
        if( present(a0) )then
            a = a0
        else
            call fit_lattice(self%element_key, centers_A, a)
        endif
        call run_cn_analysis(self%element_key,centers_A,a,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        ! calc strain and lattice displacements for all atoms
        allocate(strain_array(self%n_cc,NSTRAIN_COMPS), source=0.)
        call strain_analysis(self%element_key, centers_A, a, strain_array)
        if( allocated(self%coords4stats) ) call self%pack_instance4stats(strain_array)
        allocate(cc_mask(self%n_cc), source=.true.) ! because self%n_cc might change after pack_instance4stats
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms)
        call self%validate_atoms(simatms, l_print=.true.)
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
        if( present(imat) )then
            imat_cc = imat
            call self%img_cc%new_bimg(self%ldim, self%smpd)
            call self%img_cc%set_imat(imat_cc)
        else
            call self%img_cc%get_imat(imat_cc)
        endif
        call fit_isotropic%new(self%img_raw%get_ldim(), self%img_raw%get_smpd())
        call fit_anisotropic%new(self%img_raw%get_ldim(), self%img_raw%get_smpd())
        max_size = 0
        do cc = 1, self%n_cc
            call progress(cc, self%n_cc)
            ! index of the connected component
            self%atominfo(cc)%cc_ind = cc
            ! number of voxels in connected component
            where( imat_cc == cc ) mask = .true.
            self%atominfo(cc)%size    = count(mask)
            ! distance from the centre of mass of the nanoparticle
            self%atominfo(cc)%cendist = euclid(self%atominfo(cc)%center(:), self%NPcen) * self%smpd
            ! atom diameter
            call self%calc_longest_atm_dist(cc, self%atominfo(cc)%diam, imat=imat)
            self%atominfo(cc)%diam = 2.*self%atominfo(cc)%diam ! radius --> diameter in A
            ! whether atom is neighbors with an atom of CN > 12
            if( maxval(self%atominfo(:)%cn_std) > 12 )then
                call self%check_neighbors_cn(cc, a)
            endif
            ! maximum grey level intensity across the connected component
            self%atominfo(cc)%max_int = maxval(rmat_raw(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            ! average grey level intensity across the connected component
            self%atominfo(cc)%avg_int = sum(rmat_raw(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            self%atominfo(cc)%avg_int = self%atominfo(cc)%avg_int / real(count(mask))
            ! bond length of nearest neighbour...
            self%atominfo(cc)%bondl   = pixels_dist(self%atominfo(cc)%center(:), tmpcens, 'min', mask=cc_mask) ! Use all the atoms
            self%atominfo(cc)%bondl   = self%atominfo(cc)%bondl * self%smpd ! convert to A
            ! atomic displacement
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
        enddo
        write(logfhandle,'(a,i5)') "ADP Tossed: ", count(self%atominfo(:)%tossADP)
        write(logfhandle, '(A)') '>>> WRITING OUTPUT'
        self%n_aniso = self%n_cc - count(self%atominfo(:)%tossADP)
        call fit_isotropic%write(string(fn_fit_isotropic))
        call fit_anisotropic%write(string(fn_fit_anisotropic))
        if( WRITE_OUTPUT ) call write_2D_slice()
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
        enddo
        ! write pdf files with valid_corr and max_int in the B-factor field (for validation/visualisation)
        call self%write_centers(string('valid_corr_in_bfac_field.pdb'), 'valid_corr')
        call self%write_centers(string('max_int_in_bfac_field.pdb'),    'max_int')
        call self%write_centers(string('cn_std_in_bfac_field.pdb'),     'cn_std')
        call self%write_centers(string('u_iso_in_bfac_field.pdb'),      'u_iso')
        call self%write_centers(string('doi_in_bfac_field.pdb'),        'doi')
        call self%write_centers(string('doi_min_in_bfac_field.pdb'),    'doi_min')
        call self%write_centers_aniso(string('aniso_bfac_field.pdb'))
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
                ! generate masks
                cn_mask   = self%atominfo(:)%cn_std == cn
                size_mask = self%atominfo(:)%size >= NVOX_THRESH .and. cn_mask
                adp_mask  = (.not. self%atominfo(:)%tossADP) .and. cn_mask
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

            ! work around with imat
            subroutine write_cn_atoms( cn_std )
                integer, intent(in)  :: cn_std
                type(image_bin)       :: img_atom
                type(image)          :: simatms
                type(atoms)          :: atoms_obj
                integer, allocatable :: imat(:,:,:), imat_atom(:,:,:)
                logical :: cn_mask(self%n_cc)
                integer :: i
                ! make cn mask
                cn_mask = self%atominfo(:)%cn_std == cn_std
                call self%simulate_atoms(simatms, mask=cn_mask, atoms_obj=atoms_obj)
                call atoms_obj%writepdb(string('atoms_cn'//int2str_pad(cn_std,2)//'.pdb'))
                call simatms%write(string('simvol_cn'//int2str_pad(cn_std,2)//'.mrc'))
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
                call img_atom%write_bimg(string('binvol_cn'//int2str_pad(cn_std,2)//'.mrc'))
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
                character(len=*), parameter :: fn_slice='cc_2D_slice.csv'
                character(len=*), parameter :: header='X'//CSV_DELIM//'Y'//CSV_DELIM//'INTENSITY'
                ! Find largest CC for best visualization
                max_size   = 0
                cc_largest = 0
                do cc = 1, self%n_cc
                    if( self%atominfo(cc)%size > max_size )then
                        cc_largest = cc
                        max_size   = self%atominfo(cc)%size
                    endif
                enddo
                center(:) = int(self%atominfo(cc_largest)%center(:))
                ! Write output as CSV file with the cc center at the origin for easy analysis
                call fopen(funit, file=string(fn_slice), status='replace')
                write(funit, '(A)') 'X'//CSV_DELIM//'Y'//CSV_DELIM//'INTENSITY'
                601 format(F10.6,A2)
                602 format(F10.6)
                do i = 1, self%ldim(1)
                    do j = 1, self%ldim(2)
                        if( imat_cc(i, j, center(3)) == cc_largest )then
                            write(funit,601,advance='no') real(i-center(1)),    CSV_DELIM ! X
                            write(funit,601,advance='no') real(j-center(2)),    CSV_DELIM ! Y
                            write(funit,602) rmat_raw(i,j,center(3))                      ! Intensity
                        endif
                    enddo
                enddo
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
        enddo
        m = m / real(self%n_cc)
    end function masscen

    subroutine calc_longest_atm_dist( self, label, longest_dist, imat )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: label
        real,                intent(out)   :: longest_dist
        integer, optional,   intent(in)    :: imat(:,:,:)
        integer, allocatable :: pos(:,:)
        integer, allocatable :: imat_cc(:,:,:)
        logical, allocatable :: mask_dist(:) ! for min and max dist calculation
        integer :: location(1)               ! location of vxls of the atom farthest from its center
        if( present(imat) )then
            imat_cc = imat
        else
            call self%img_cc%get_imat(imat_cc)
        endif
        where( imat_cc .eq. label )
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
    subroutine check_neighbors_cn( self, cc, a )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: cc
        real,                intent(in)    :: a(3)
        character(len=5)  :: el_ucase
        character(len=10) :: crystal_system
        real              :: a0, foo(3), d
        integer           :: i
        a0 = sum(a)/real(size(a)) ! alrithmetic mean of fitted lattice parameters
        ! identify nearest neighbors within the first-shell cutoff, as run_cn_analysis
        el_ucase = uppercase(trim(adjustl(self%element_key)))
        call get_lattice_params(el_ucase, crystal_system, foo)
        if( trim(crystal_system) == 'wurtzite' )then
            d = lattice_cutoff(crystal_system, a)
        else
            d = lattice_cutoff(crystal_system, [a0, a0, a0])
        endif
        do i = 1, self%n_cc
            if( i/=cc .and. self%atominfo(i)%cn_std > 12 .and. &
                    & euclid(self%atominfo(cc)%center(:3),self%atominfo(i)%center(:3))*self%smpd < d )then
                self%atominfo(cc)%adjacent_cn13 = self%atominfo(cc)%adjacent_cn13 + 1
            endif
        enddo
    end subroutine check_neighbors_cn

    ! the crystal system of the element is one of the binary crystals (rocksalt, zincblende, wurtzite)
    logical function np_binary_lattice( self )
        class(nanoparticle), intent(in) :: self
        character(len=10) :: crystal_system
        character(len=5)  :: el_ucase
        real              :: foo(3)
        el_ucase = uppercase(trim(adjustl(self%element_key)))
        call get_lattice_params(el_ucase, crystal_system, foo)
        np_binary_lattice = binary_lattice(crystal_system)
    end function np_binary_lattice

    ! bond length of the element's crystal system for the fitted lattice parameters a
    real function np_lattice_bond( self, a )
        class(nanoparticle), intent(in) :: self
        real,                intent(in) :: a(3)
        character(len=10) :: crystal_system
        character(len=5)  :: el_ucase
        real              :: foo(3)
        el_ucase = uppercase(trim(adjustl(self%element_key)))
        call get_lattice_params(el_ucase, crystal_system, foo)
        np_lattice_bond = lattice_bond(crystal_system, a)
    end function np_lattice_bond

    ! Calculates the isotropic displacement parameter U_ISO of a given cc by fitting
    ! the real space intensity distribution of the 3D reconstructed volume stored
    ! in rmat with a 3D Gaussian with isotropic variance in a sphere of radius 3*a/(8sqrt(2))
    ! about the cc center. The best fit solution is output into the 3D map fit, and
    ! the correlation between fit and rmat is stored in atominfo.
    subroutine calc_isotropic_disp( self, cc, a, rmat, fit )
        class(nanoparticle), intent(inout) :: self
        class(image), intent(inout)        :: fit
        integer, intent(in)                :: cc
        real, intent(in)                   :: a(3) ! lattice params
        real, pointer, intent(in)          :: rmat(:,:,:)
        real    :: output_rad, fit_rad, center(3), maxrad, r, var, amp, int, max_int_out
        real    :: XTWX(2,2), XTWX_inv(2,2), XTWY(2,1), B(2,1)
        integer :: ilo, ihi, jlo, jhi, klo, khi, i, j, k, errflg
        logical :: fit_mask(self%ldim(1),self%ldim(2),self%ldim(3))

        if( self%binary_lattice() )then
            ! half the bond of the binary crystal, in the same proportions as the fcc radii below
            output_rad = 0.5 * self%lattice_bond(a) / self%smpd
            fit_rad    = 0.75 * output_rad
            center     = self%atominfo(cc)%center(:)
            maxrad     = sqrt(2.) * output_rad
        else
            output_rad = (sum(a)/3)/(2.*sqrt(2.))/self%smpd  ! 1/2 FCC nearest neighbor dist
            fit_rad    = 0.75 * output_rad ! 0.75 prevents fitting tails of other atoms
            ! Create search window containing sphere of fit rad to speed up loops.
            center     = self%atominfo(cc)%center(:)
            maxrad     = 0.5*(sum(a)/3) / self%smpd
        endif
        ilo        = max(nint(center(1) - maxrad), 1)
        ihi        = min(nint(center(1) + maxrad), self%ldim(1))
        jlo        = max(nint(center(2) - maxrad), 1)
        jhi        = min(nint(center(2) + maxrad), self%ldim(2))
        klo        = max(nint(center(3) - maxrad), 1)
        khi        = min(nint(center(3) + maxrad), self%ldim(3))

        ! Linear least squares to calculate the best fit params B.
        ! Note the sum of the square error of the logs in minimized instead
        ! of the sum of the square error to make the problem linear.
        ! Solution: B = ((X^T)WX)^-1((X^T)WY)
        XTWX     = 0.
        XTWX_inv = 0.
        XTWY     = 0.
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi
                    r = euclid(1.*(/i, j, k/), 1.*center)
                    if( r < fit_rad )then
                        int = rmat(i,j,k)
                        if( int > 0 )then
                            XTWX(1,1) = XTWX(1,1) + int
                            XTWX(1,2) = XTWX(1,2) + r**2 * int
                            XTWX(2,2) = XTWX(2,2) + r**4 * int
                            XTWY(1,1) = XTWY(1,1) + int * log(int)
                            XTWY(2,1) = XTWY(2,1) + r**2 * int * log(int)
                        endif
                    endif
                enddo
            enddo
        enddo
        XTWX(2,1) = XTWX(1,2)
        call matinv(XTWX, XTWX_inv, 2, errflg)
        B   = matmul(XTWX_inv, XTWY)
        amp = exp(B(1,1))   ! Best fit peak
        var = -0.5 / B(2,1) ! Best fit variance

        ! Sample fit for goodness of fit and visualization
        max_int_out = 0.
        fit_mask    = .false.
        do k = klo, khi 
            do j = jlo, jhi
                do i = ilo, ihi
                    r = euclid(1.*(/i, j, k/), 1.*center) 
                    if( r < output_rad )then
                        fit_mask(i,j,k) = .true.
                        int =  amp * exp(-0.5 * r**2 / var)
                        call fit%set_rmat_at(i, j, k, int)
                        if( int > max_int_out )then
                            max_int_out = int
                        endif
                    endif
                enddo
            enddo
        enddo
        self%atominfo(cc)%isocorr = fit%real_corr(self%img_raw, mask=fit_mask)
        self%atominfo(cc)%u_iso   = var * self%smpd**2
    end subroutine calc_isotropic_disp

    ! Calculates the anisotropic displacement parameter matrix of a given cc by fitting
    ! the real space intensity distribution of the 3D reconstructed volume stored
    ! in rmat with a 3D multivariate Gaussian in a sphere of radius 3*a/(8sqrt(2))
    ! about the cc center. The best fit solution is output into the 3D map fit, and
    ! the correlation between fit and rmat is stored in atominfo. The variances
    ! along the 3 principal axes and the angular orientation of the major axis are also
    ! calculated.
    subroutine calc_anisotropic_disp( self, cc, a, rmat, fit )
        class(nanoparticle), intent(inout) :: self
        class(image),        intent(inout) :: fit
        integer,             intent(in)    :: cc
        real,                intent(in)    :: a(3) ! Lattice params
        real, pointer,       intent(in)    :: rmat(:,:,:)
        real(kind=8), allocatable :: X(:,:), XTW(:,:), Y(:,:)
        real(kind=8) :: XTWX(7,7), XTWX_inv(7,7), XTWY(7,1), B(7,1), cov(3,3), cov_inv(3,3)
        real(kind=8) :: eigenvals(3), eigenvecs(3,3)
        real(kind=8) :: majvector(3), rvec(3,1), beta(1,1), r(3)
        real         :: center(3), output_rad, fit_rad, maxrad, int, amp, max_int_out, corr
        integer      :: ilo, ihi, jlo, jhi, klo, khi, i, j, k, nvoxels, n, errflg, nrot
        logical      :: fit_mask(self%ldim(1), self%ldim(2), self%ldim(3))

        if( self%binary_lattice() )then
            ! half the bond of the binary crystal, in the same proportions as the fcc radii below
            output_rad = 0.5 * self%lattice_bond(a)
            fit_rad    = 0.75 * output_rad
            center     = self%atominfo(cc)%center(:)
            maxrad     = sqrt(2.) * output_rad / self%smpd
        else
            output_rad = ( sum(a) / 3 ) / ( 2. * sqrt(2.) )  ! 1/2 FCC nearest neighbor dist
            fit_rad    = 0.75 * output_rad ! 0.75 prevents fitting tails of other atoms
            ! Create search window containing sphere of fit rad to speed up loops.
            center  = self%atominfo(cc)%center(:)
            maxrad  = 0.5 * (sum(a)/3) / self%smpd
        endif
        ilo     = max(nint(center(1) - maxrad), 1)
        ihi     = min(nint(center(1) + maxrad), self%ldim(1))
        jlo     = max(nint(center(2) - maxrad), 1)
        jhi     = min(nint(center(2) + maxrad), self%ldim(2))
        klo     = max(nint(center(3) - maxrad), 1)
        khi     = min(nint(center(3) + maxrad), self%ldim(3))
        ! Get nvoxels within sphere of fitting
        nvoxels = 0
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi
                    if( euclid(1.*(/i, j, k/), 1.*center)*self%smpd < fit_rad )then
                        nvoxels = nvoxels + 1
                    endif
                enddo
            enddo
        enddo
        allocate(X(nvoxels, 7), XTW(7, nvoxels), Y(nvoxels,1), source = 0._dp)
        ! Linear least squares to calculate the best fit params B
        ! Conduct fit in units of Angstroms (easier on matrix operations)
        ! Note the sum of the square error of the logs is minimized instead
        ! of the sum of the square error to make the problem linear.
        ! Solution: B = ((X^T)WX)^-1((X^T)WY)
        XTWX          = 0._dp
        XTWX_inv      = 0._dp
        XTWY          = 0._dp
        n             = 1
        X(:nvoxels,1) = 1.
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi
                    r = (1.*(/i, j, k/) - center) * self%smpd
                    if( norm_2(r) < fit_rad )then
                        int = rmat(i,j,k)
                        if( int > 0 .and. n <= nvoxels )then
                            X(n,2:4)  = r(1:3)**2
                            X(n,5)    = r(1)*r(2)
                            X(n,6)    = r(1)*r(3)
                            X(n,7)    = r(2)*r(3)
                            XTW(:7,n) = int*X(n,:7)
                            Y(n,1)    = log(int)
                            n = n + 1
                        endif
                    endif
                enddo
            enddo
        enddo
        XTWX = matmul(XTW, X)
        XTWY = matmul(XTW, Y)
        deallocate(X,XTW,Y)
        call matinv(XTWX, XTWX_inv, 7, errflg)
        B   = matmul(XTWX_inv, XTWY)
        amp = real(exp(B(1,1)), kind=kind(amp)) ! Best fit peak
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
        self%atominfo(cc)%aniso = real(cov, kind=kind(self%atominfo(cc)%aniso))
        ! Sample fit for goodness of fit and visualization
        max_int_out = 0.
        fit_mask    = .false.
        do k = klo, khi 
            do j = jlo, jhi   
                do i = ilo, ihi
                    rvec(:3,1) = (1.*(/i, j, k/) - center) * self%smpd
                    if( norm_2(rvec(:3,1)) < output_rad )then
                        fit_mask(i,j,k) = .true.
                        beta            = matmul(matmul(transpose(rvec),cov_inv),rvec)
                        int             = real(amp * exp(-0.5 * beta(1,1)), kind=kind(int))
                        call fit%set_rmat_at(i, j, k, int)
                        if( int > max_int_out )then
                            max_int_out = int
                        endif
                    endif
                enddo
            enddo
        enddo
        corr = fit%real_corr(self%img_raw, mask=fit_mask)
        self%atominfo(cc)%anisocorr = corr 
        ! Calculate eigenvalues and orientation of major eigenvector
        call jacobi(cov, 3, 3, eigenvals, eigenvecs, nrot)
        call eigsrt(eigenvals, eigenvecs, 3, 3)
        ! Any negative eigenvalue or large eigenvalue implies intensity 
        ! distribution can't be approximated as a multivariate Gaussian
        if( any(eigenvals <= 0._dp) .or. any(eigenvals > (1.25*fit_rad)**2) )then
            self%atominfo(cc)%tossADP   = .true.
            self%atominfo(cc)%u_evals   = 0.
            self%atominfo(cc)%anisocorr = 0.
            self%atominfo(cc)%azimuth   = 0.
            self%atominfo(cc)%polar     = 0.
            self%atominfo(cc)%doi       = 1. 
            self%atominfo(cc)%doi_min   = 0.
            self%atominfo(cc)%aniso     = 0.
            return
        endif
        self%atominfo(cc)%u_evals(:3) = real(eigenvals(:3), kind=kind(self%atominfo(cc)%u_evals))
        self%atominfo(cc)%doi         = real(eigenvals(3) / eigenvals(1), kind=kind(self%atominfo(cc)%doi))
        self%atominfo(cc)%doi_min     = real(eigenvals(2) / eigenvals(1), kind=kind(self%atominfo(cc)%doi_min))
        ! Find azimuthal and polar angles of the major eigenvector
        majvector = eigenvecs(:,1)
        ! Sign of majvector is arbitrary. Use convention that y-coordinate must be >= 0
        if( majvector(2) < 0. )then
            majvector = -1. * majvector
        endif
        self%atominfo(cc)%azimuth = real(atan(majvector(2) / majvector(1)), kind=kind(self%atominfo(cc)%azimuth))
        if( majvector(1) > 0. )then
            self%atominfo(cc)%azimuth = real(atan(majvector(2) / majvector(1)), kind=kind(self%atominfo(cc)%azimuth))
        elseif( majvector(1) < 0. )then
            self%atominfo(cc)%azimuth = real(atan(majvector(2) / majvector(1)) + PI, kind=kind(self%atominfo(cc)%azimuth))
        else
            self%atominfo(cc)%azimuth = PI / 2
        endif
        if( majvector(3) > 0. )then
            self%atominfo(cc)%polar = real(atan(sqrt(majvector(1)**2 + majvector(2)**2) / majvector(3)), &
                &kind=kind(self%atominfo(cc)%polar))
        elseif( majvector(3) < 0. )then
            self%atominfo(cc)%polar = real(atan(sqrt(majvector(1)**2 + majvector(2)**2) / majvector(3)) + PI, &
                &kind=kind(self%atominfo(cc)%polar))
        else
            self%atominfo(cc)%polar = PI / 2
        endif
    end subroutine calc_anisotropic_disp

    ! visualization and output

    subroutine simulate_atoms( self, simatms, betas, mask, atoms_obj )
        class(nanoparticle),           intent(inout) :: self
        class(image),                  intent(out)   :: simatms
        real,        optional,         intent(in)    :: betas(self%n_cc) ! in pdb file b-factor
        logical,     optional,         intent(in)    :: mask(self%n_cc)
        type(atoms), optional, target, intent(inout) :: atoms_obj
        type(atoms), target  :: atms_here
        type(atoms), pointer :: atms_ptr => null()
        logical              :: betas_present, mask_present, atoms_obj_present
        integer              :: i, cnt
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
            enddo
        else
            do i = 1, self%n_cc
                call set_atom(i, i)
            enddo
        endif
        call simatms%new(self%ldim, self%smpd)
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
                ! a pseudo-atom is rendered at the template width
                if( self%l_species_free ) call atms_ptr%set_beta(atms_obj_ind, B_REF_PSEUDO)
            end subroutine set_atom

    end subroutine simulate_atoms

    subroutine write_centers_1( self, fname, coords )
        class(nanoparticle),     intent(inout) :: self
        class(string), optional, intent(in)    :: fname
        real,          optional, intent(in)    :: coords(:,:)
        type(atoms) :: centers_pdb
        integer     :: cc
        if( present(coords) )then
            call centers_pdb%new(size(coords, dim = 2), dummy=.true.)
            do cc = 1, size(coords, dim = 2)
                call set_atm_info
            enddo
        else
            call centers_pdb%new(self%n_cc, dummy=.true.)
            do cc = 1, self%n_cc
                call set_atm_info
            enddo
        endif
        if( present(fname) )then
            call centers_pdb%writepdb(fname)
        else
            call centers_pdb%writepdb(self%fbody//'_ATMS.pdb')
            write(logfhandle,'(A)') 'output, atomic coordinates:       '//self%fbody%to_char()//'_ATMS.pdb'
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
        class(string),       intent(in)    :: fname
        character(len=*),    intent(in)    :: which ! parameter in the B-factor field of the pdb file
        type(atoms) :: centers_pdb
        integer     :: cc
        call centers_pdb%new(self%n_cc, dummy=.true.)
        do cc = 1, self%n_cc
            call centers_pdb%set_name(cc,self%atom_name)
            call centers_pdb%set_element(cc,self%element)
            call centers_pdb%set_coord(cc,(self%atominfo(cc)%center(:)-1.)*self%smpd)
            select case(which)
                case('cn_gen')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%cn_gen)       ! use generalised coordination number
                case('cn_std')
                    call centers_pdb%set_beta(cc,real(self%atominfo(cc)%cn_std)) ! use standard coordination number
                case('max_int')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%max_int)      ! use z-score of maximum intensity
                case('doi')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%doi)          ! use isotropic b-factor
                case('doi_min')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%doi_min)      ! use isotropic b-factor
                case('u_iso')
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%u_iso)        ! use isotropic b-factor
                case DEFAULT
                    call centers_pdb%set_beta(cc,self%atominfo(cc)%valid_corr)   ! use per-atom validation correlation
            end select
            call centers_pdb%set_resnum(cc,cc)
        enddo
        call centers_pdb%writepdb(fname)
    end subroutine write_centers_2

    subroutine write_centers_aniso( self, fname )
        class(nanoparticle), intent(inout) :: self
        class(string),       intent(in)    :: fname
        type(atoms) :: centers_pdb
        real        :: aniso(3, 3, self%n_cc)
        integer     :: cc
        call centers_pdb%new(self%n_cc, dummy=.true.)
        do cc = 1, self%n_cc
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
        type(image_bin)      :: img_atom
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
            call img_atom%write_bimg(string('Atom'//int2str(i)//'.mrc'))
        enddo
        deallocate(imat,imat_atom)
        call img_atom%kill_bimg
    end subroutine write_individual_atoms

    subroutine write_csv_files( self )
        class(nanoparticle), intent(in) :: self
        integer               :: ios, funit, cc, cn
        character(len=STDLEN) :: io_msg
        ! NANOPARTICLE STATS
        call fopen(funit, file=string(NP_STATS_FILE), iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//NP_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        write(funit,'(a)') NP_STATS_HEAD
        ! write record
        call self%write_np_stats(funit)
        call fclose(funit)
        ! CN-DEPENDENT STATS
        call fopen(funit, file=string(CN_STATS_FILE), iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//CN_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        write(funit,'(a)') CN_STATS_HEAD
        ! write records
        do cn = CNMIN, CNMAX
            call self%write_cn_stats(cn, funit)
        enddo
        call fclose(funit)
        ! PER-ATOM STATS
        call fopen(funit, file=string(ATOMS_STATS_FILE), iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//ATOMS_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        if( ATOMS_STATS_OMIT )then
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
        ! -- isotropic displacement parameter
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
        ! -- azimuthal angle of major eigenvector 
        write(funit,601,advance='no') self%azimuth_stats%avg,        CSV_DELIM ! AVG_AZIMUTH          (40)
        write(funit,601,advance='no') self%azimuth_stats%med,        CSV_DELIM ! MED_AZIMUTH          (41)
        write(funit,601,advance='no') self%azimuth_stats%sdev,       CSV_DELIM ! SDEV_AZIMUTH         (42)
        ! -- polar angle of major eigenvector
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
        ! -- isotropic displacement fit correlation
        write(funit,601,advance='no') self%isocorr_stats%avg,        CSV_DELIM ! AVG_ISO_CORR         (52)
        write(funit,601,advance='no') self%isocorr_stats%med,        CSV_DELIM ! MED_ISO_CORR         (53)
        write(funit,601,advance='no') self%isocorr_stats%sdev,       CSV_DELIM ! SDEV_ISO_CORR        (54)
        ! -- anisotropic displacement fit correlation
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
        if( count(self%atominfo(:)%cn_std == cn) < 2 )return
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
        ! -- azimuthal angle of major eigenvector
        write(funit,601,advance='no') self%azimuth_stats_cns(cn)%avg,        CSV_DELIM ! AVG_AZIMUTH
        write(funit,601,advance='no') self%azimuth_stats_cns(cn)%med,        CSV_DELIM ! MED_AZIMUTH
        write(funit,601,advance='no') self%azimuth_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_AZIMUTH
        ! -- polar angle of major eigenvector
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

    subroutine kill( self )
        class(nanoparticle), intent(inout) :: self
        call self%img%kill()
        call self%img_raw%kill
        call self%img_bin%kill_bimg()
        call self%img_cc%kill_bimg()
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
        if( allocated(self%recovered) ) deallocate(self%recovered)
        if( allocated(self%species_post) ) deallocate(self%species_post)
        if( allocated(self%species) ) deallocate(self%species)
        self%l_species_free     = .false.
        self%d_nn               = 0.
        self%prune_cen          = 0.
        self%prune_rad_thres    = 0.
        self%prune_cs_thres     = 0
        self%l_discover_species = .false.
        self%min_nbrs           = 3
        self%nspecies           = 0
        self%msk_rad            = 0.
        call self%vol_even%kill
        call self%vol_odd%kill
    end subroutine kill

end module simple_nanoparticle

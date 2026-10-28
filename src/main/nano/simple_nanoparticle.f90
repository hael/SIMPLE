
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
real,             parameter :: MSKDIAM_OVER_DIAM   = 1.2     ! estimated mask diameter over the particle diameter
! detection is species-free for every project: pseudo-atom template, length scales from the nearest-neighbour
! distance d_NN and the contact-score ceiling from the data; an element serves the lattice analyses of atoms_stats only
real,             parameter :: B_REF_PSEUDO        = 13.9    ! B factor of the pseudo-atom template, A**2 (sigma 0.42 A)
real,             parameter :: RADIUS_NN_FRAC      = 0.4     ! atom radius in units of d_NN
real,             parameter :: SPLIT_EXCL_NN_FRAC  = 0.21    ! exclusion radius of atom splitting in units of d_NN
real,             parameter :: CSCORE_INNER_FRAC   = 0.3     ! the contact ceiling is the median score of this inner fraction

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
    ! the weak atoms of the residual recovery
    real    :: aper_int          = 0. ! aperture intensity, the cross-check of the fitted intensity
    real    :: det_z             = 0. ! detection z of a weak atom
    integer :: det_stage         = 0  ! 0 strong (level 1), 1 weak of residual stage A, 2 weak of stage B
    integer :: det_level         = 0  ! residual level that accepted the atom
    logical :: gate_used         = .false. ! accepted through the neighbour gate
    ! the atomic density model (atomic_density_model.md), written to _ATMS.pdb and _atoms.csv
    real    :: int_fit           = 0. ! integrated intensity I
    real    :: bfac_fit          = 0. ! width B, A**2
    real    :: sd_lnint          = 0. ! measurement spread of ln I
    real    :: sd_lnbfac         = 0. ! measurement spread of ln B
    real    :: center_fit(3)     = 0. ! fitted position, voxels
    integer :: cn_fit            = 0  ! coordination number on the fitted positions
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
    real                  :: d_nn               = 0. ! median nearest-neighbour distance in A
    real                  :: int_scale          = 0. ! intensity scale S of the density model, the occupancy unit
    real                  :: int_noise          = 0. ! measurement spread of an intensity at fixed width (s_int)
    logical               :: l_noise_resid      = .false. ! the model's noise is the residual between the atoms
    ! the pruning policy of discard_atoms, kept for the weak atoms: centre of mass (voxels), the radius beyond which
    ! atoms are subject to deletion (A) and the contact-score threshold
    real                  :: prune_cen(3)       = 0.
    real                  :: prune_rad_thres    = 0.
    integer               :: prune_cs_thres     = 0
    integer               :: cscore_ceil        = 0       ! contact-score ceiling discard_atoms derived from the data
    logical               :: l_species_free     = .false. ! no element given: atoms written as X1, no lattice analyses
    integer               :: min_nbrs           = 3       ! gated weak atoms need this many found neighbours, 0 = no gate
    real                  :: msk_rad            = 0. ! radius of the soft spherical mask in voxels
    type(string)          :: transfer_tab            ! a transfer table to use instead of measuring it (_transfer.tsv)
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
    real,             allocatable :: coords4stats(:,:)
    ! OTHER
    character(len=2)      :: element     = '  '      ! symbol of the first species: radii, templates, atom names
    character(len=5)      :: element_key = '     '   ! element as given (a compound selector): lattice lookups
    character(len=4)      :: atom_name   = '    '
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
    procedure, private :: fit_density_model
    procedure, private :: empty_zone
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
        real    :: msk_in_pix, diam
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
            call mskvol%estimate_spher_mask_diam(params, self%img, AMSKLP_NANO, diam)
            if( diam <= 0. ) THROW_HARD('the mask estimate found no particle; give mskdiam')
            msk_in_pix = MSKDIAM_OVER_DIAM * diam / (2. * self%smpd)
            write(logfhandle,*) 'particle diameter in A: ', diam
            write(logfhandle,*) 'mask diameter in A: ', 2. * msk_in_pix * self%smpd
            call self%img%mask3D_soft(msk_in_pix)
            call mskvol%kill_bimg
            self%msk_rad = msk_in_pix
        endif
        self%min_nbrs           = params%min_nbrs
        if( params%transfer_tab%is_allocated() )then
            if( .not. params%transfer_tab%is_blank() ) self%transfer_tab = params%transfer_tab
        endif
        if( DEBUG ) call self%img%write(string('masked_input_vol.mrc'))
        call self%img_raw%copy(self%img)
        call self%img_raw%stats(self%map_stats%avg, self%map_stats%sdev, self%map_stats%maxv, self%map_stats%minv)
    end subroutine new

    ! getters/setters

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
        call phasecorr_one_atom(self%img, B_REF_PSEUDO)
        call self%img%write(fname)
    end subroutine conv_denoise

    subroutine identify_lattice_params( self, a )
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: a(3) ! lattice parameters
        real, allocatable :: centers_A(:,:)        ! coordinates of the atoms in ANGSTROMS
        type(image)       :: simatms
        type(string)      :: errmsg
        if( self%l_species_free ) THROW_HARD('identify_lattice_params needs an element for the lattice fit')
        ! MODEL BUILDING, species-free as in detect_atoms; the element serves the lattice fit only
        ! phase correlation approach
        call phasecorr_one_atom(self%img, B_REF_PSEUDO)
        ! nanoparticle binarization
        call self%binarize_and_find_centers()
        if( self%n_cc < 2 )then
            errmsg = 'the lattice fit needs two or more level-1 atoms to measure d_NN; level 1 found '//int2str(self%n_cc)
            THROW_HARD(errmsg%to_char())
        endif
        call self%set_nn_length_scales
        ! discard small connected components
        call self%discard_small_ccs
        ! atom splitting by correlation map validation
        call self%split_atoms()
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms, l_pseudo=.true.)
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
        call phasecorr_one_atom(self%img, B_REF_PSEUDO)
        if( DEBUG ) call self%img%write(string('after_phasecorr.mrc'))
        ! Nanoparticle binarization
        call self%binarize_and_find_centers(l_print=ll_print)
        if( self%n_cc < 2 )then
            errmsg = 'detect_atoms needs two or more level-1 atoms to measure d_NN; level 1 found '//int2str(self%n_cc)
            THROW_HARD(errmsg%to_char())
        endif
        call self%set_nn_length_scales
        ! discard small connected components
        call self%discard_small_ccs
        ! atom splitting by correlation map validation
        call self%split_atoms(split_fname)
        ! validation through per-atom correlation with the simulated density
        call self%simulate_atoms(simatms, l_pseudo=.true.)
        if( WRITE_OUTPUT ) call simatms%write(self%fbody//'_SIM_pre_validate.mrc')
        call self%validate_atoms(simatms, l_print=.false.)
        if( WRITE_OUTPUT ) call self%write_centers(string('valid_corr_in_bfac_field_pre_discard.pdb'), 'valid_corr')
        if( l_atom_thres ) call self%discard_atoms(l_print=ll_print)
        ! re-calculate valid_corr:s (since they are otherwise lost from the B-factor field due to reallocations of atominfo)
        call self%simulate_atoms(simatms, l_pseudo=.true.)
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
        ! the density model of the atoms (_SIM, _FIT, _RES, _atoms.csv, _model.txt), whose columns _ATMS.pdb carries
        call self%fit_density_model
        call self%write_centers
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
        use simple_bench_io, only: open_bench_report
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
            call open_bench_report(fnr, 'BINARIZE_AND_FIND_CENTERS')
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
            call set_pseudo_bfacs(atom)
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

        ! pseudo-atoms at the template width (the B column of write_centers holds valid_corr, the element column the
        ! element of atoms_stats)
        subroutine set_pseudo_bfacs( atms )
            type(atoms), intent(inout) :: atms
            integer :: iatm
            do iatm = 1,atms%get_n()
                call atms%set_element(iatm, 'X1')
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
        integer, allocatable :: imat_bin(:,:,:), first(:), vox(:,:)
        integer :: cc, iv
        call self%img_cc%get_cc_voxels(self%n_cc, first, vox)
        call self%img_bin%get_imat(imat_bin)
        do cc = 1, self%n_cc
            if( first(cc+1) - first(cc) < NVOX_THRESH )then ! removes artificial small densities
                do iv = first(cc), first(cc+1) - 1
                    imat_bin(vox(1,iv),vox(2,iv),vox(3,iv)) = 0
                enddo
            endif
        enddo
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        deallocate(imat_bin, first, vox)
    end subroutine discard_small_ccs

    subroutine split_atoms( self, fname )
        class(nanoparticle),     intent(inout) :: self
        class(string), optional, intent(in)    :: fname
        type(image_bin)       :: img_split_ccs
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:)
        integer, allocatable :: imat_cc(:,:,:), imat_bin(:,:,:), imat_split_ccs(:,:,:), first(:), vox(:,:)
        logical, allocatable :: claimed(:) ! voxels of the component being split within split_excl of a centre
        integer, parameter   :: RANK_THRESH = 4
        integer :: icc, cnt, cnt_split, iv
        integer :: rank, m(1)
        real    :: new_centers(3,3*self%n_cc) ! will pack it afterwards if it has too many elements
        real    :: pc, radius, split_excl
        write(logfhandle, '(A)') '>>> SPLITTING CONNECTED ATOMS'
        ! squared exclusion radius as the distance tests below use it: voxels**2 * smpd
        split_excl = (SPLIT_EXCL_NN_FRAC * self%d_nn)**2 / self%smpd
        call self%img%get_rmat_ptr(rmat_pc) ! rmat_pc contains the phase correlation
        call self%img_cc%get_imat(imat_cc)
        ! the voxels of each component as found here; splitting edits imat_cc only inside the component it splits
        call self%img_cc%get_cc_voxels(self%n_cc, first, vox)
        allocate(imat_split_ccs(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = 0)
        call img_split_ccs%new_bimg(self%ldim, self%smpd)
        call img_split_ccs%new(self%ldim, self%smpd)
        cnt       = 0
        cnt_split = 0
        do icc = 1, self%n_cc ! for each cc check if the center corresponds with the local max of the phase corr
            pc = rmat_pc(nint(self%atominfo(icc)%center(1)),nint(self%atominfo(icc)%center(2)),nint(self%atominfo(icc)%center(3)))
            ! calculate the rank
            allocate(x(first(icc+1)-first(icc)))
            do iv = first(icc), first(icc+1) - 1
                x(iv-first(icc)+1) = rmat_pc(vox(1,iv),vox(2,iv),vox(3,iv))
            enddo
            call hpsort(x)
            m(:) = minloc(abs(x - pc))
            rank = size(x) - m(1)
            deallocate(x)
            ! calculate radius
            call self%calc_longest_atm_dist(icc, vox(:,first(icc):first(icc+1)-1), radius)
            ! split
            if( rank > RANK_THRESH .or. radius > 1.5 * self%theoretical_radius )then
                do iv = first(icc), first(icc+1) - 1
                    imat_split_ccs(vox(1,iv),vox(2,iv),vox(3,iv)) = 1
                enddo
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
            if( allocated(claimed) ) deallocate(claimed)
            allocate(claimed(first(icc):first(icc+1)-1), source=.false.)
            ! Identify first new center
            new_center1 = brightest()
            cnt = cnt + 1
            new_centers(:,cnt) = real(new_center1)
            call claim(new_center1)
            ! Second likely center.
            new_center2 = brightest()
            if( any(new_center2 > 0) )then ! if anything was found
                ! Validate second center (check if it's 2 merged atoms, or one pointy one)
                if( sum(real(new_center2 - new_center1)**2.) * self%smpd <= split_excl) then
                    ! the new_center2 is within the diameter of the atom position at new_center1
                    ! therefore, it is not another atom and should be removed
                    call release_unclaimed
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center2)
                    call claim(new_center2)
                endif
            endif
            ! Third likely center.
            new_center3 = brightest()
            if( any(new_center3 > 0) )then ! if anything was found
                ! Validate third center
                if(sum(real(new_center3 - new_center1)**2.) * self%smpd <= split_excl .or. &
                &  sum(real(new_center3 - new_center2)**2.) * self%smpd <= split_excl )then
                    ! the new_center3 is within the diameter of the atom position at new_center1 or new_center2
                    ! therefore, it is not another atom and should be removed
                    call release_unclaimed
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center3)
                    call claim(new_center3)
                endif
            endif
            ! Set the merged cc back to 0
            call release_unclaimed
        end subroutine split_atom

        ! the first maximum of the phase correlation over the unclaimed voxels of the component, as maxloc
        ! with a mask returns it; zeros when every voxel is claimed
        function brightest() result( pos )
            integer :: pos(3), iv
            logical :: found
            real    :: vmax
            pos   = 0
            found = .false.
            vmax  = 0.
            do iv = lbound(claimed,1), ubound(claimed,1)
                if( claimed(iv) ) cycle
                if( .not. found )then
                    found = .true.
                    vmax  = rmat_pc(vox(1,iv),vox(2,iv),vox(3,iv))
                    pos   = vox(:,iv)
                else if( rmat_pc(vox(1,iv),vox(2,iv),vox(3,iv)) > vmax )then
                    vmax  = rmat_pc(vox(1,iv),vox(2,iv),vox(3,iv))
                    pos   = vox(:,iv)
                endif
            enddo
        end function brightest

        subroutine claim( cen )
            integer, intent(in) :: cen(3)
            integer :: iv
            do iv = lbound(claimed,1), ubound(claimed,1)
                if( ((real(vox(1,iv) - cen(1)))**2 + (real(vox(2,iv) - cen(2)))**2 + &
                &    (real(vox(3,iv) - cen(3)))**2) * self%smpd <= split_excl ) claimed(iv) = .true.
            enddo
        end subroutine claim

        subroutine release_unclaimed
            integer :: iv
            do iv = lbound(claimed,1), ubound(claimed,1)
                if( claimed(iv) ) cycle
                if( imat_cc(vox(1,iv),vox(2,iv),vox(3,iv)) == icc ) imat_cc(vox(1,iv),vox(2,iv),vox(3,iv)) = 0
            enddo
        end subroutine release_unclaimed

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
        integer, allocatable :: imat_bin(:,:,:), cscores(:), first(:), vox(:,:)
        real,    allocatable :: centers_A(:,:), cendists(:), cendists_sorted(:)
        logical, allocatable :: atom_del_mask(:)
        integer              :: cscore_thres
        type(stats_struct)   :: cscore_stats
        real,    allocatable :: inner(:)
        real    :: percen, cendist_thres, inner_thres, med_inner
        integer :: cc, cn, n_discard, cnt_discard, it_contact_score, ninner
        logical :: ll_print
        ll_print = .true.
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
        ! the contact-score ceiling from the data: the median contact score of the inner CSCORE_INNER_FRAC of the atoms by
        ! distance from the centre of mass, rounded (12 for close packing, 4 for the tetrahedral lattices)
        ninner      = max(1, nint(CSCORE_INNER_FRAC * real(self%n_cc)))
        inner_thres = cendists_sorted(ninner)
        inner       = real(pack(cscores, cendists <= inner_thres))
        call hpsort(inner)
        if( mod(size(inner), 2) == 0 )then
            med_inner = 0.5 * (inner(size(inner)/2) + inner(size(inner)/2 + 1))
        else
            med_inner = inner(size(inner)/2 + 1)
        endif
        it_contact_score = max(1, nint(med_inner))
        self%cscore_ceil = it_contact_score
        if( ll_print ) write(logfhandle,'(A,I0)') 'contact-score ceiling from the inner atoms: ', it_contact_score
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
        ! component voxels and binary matrix
        call self%img_cc%get_cc_voxels(self%n_cc, first, vox)
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
        deallocate(imat_bin, first, vox)
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
                    call zero_component(cc)
                    n_discard   = n_discard   + 1
                    cnt_discard = cnt_discard + 1
                endif
            enddo
            if( cnt_discard > 0 ) call relabel
        end subroutine remove_lowly_contacted

        subroutine remove_small_and_lowly_contacted( cthresh )
            integer, intent(in) :: cthresh
            real, allocatable   :: radii(:)
            real :: radius_thres
            ! calculate atomic radii
            allocate(radii(self%n_cc), source=0.)
            do cc = 1, self%n_cc
                call self%calc_longest_atm_dist(cc, vox(:,first(cc):first(cc+1)-1), radii(cc))
            end do
            radius_thres = self%theoretical_radius/2.
            ! discard
            cnt_discard  = 0
            do cc = 1, self%n_cc
                if( atom_del_mask(cc) )then
                    if( radii(cc) < radius_thres .and. cscores(cc) < cthresh )then
                        call zero_component(cc)
                        n_discard   = n_discard   + 1
                        cnt_discard = cnt_discard + 1
                    endif
                endif
            enddo
            if( cnt_discard > 0 ) call relabel
        end subroutine remove_small_and_lowly_contacted

        ! after discards the components are relabelled: renew everything indexed by atom, keeping the zone
        subroutine relabel
            integer :: icc
            call self%img_bin%set_imat(imat_bin)
            call self%img_bin%find_ccs(self%img_cc)
            call self%img_cc%get_nccs(self%n_cc)
            call self%find_centers()
            if( allocated(centers_A) ) deallocate(centers_A)
            if( allocated(cscores)   ) deallocate(cscores)
            allocate(cscores(self%n_cc), source=0)
            centers_A = self%atominfo2centers_A()
            call self%contact_scores(centers_A,cscores)
            call self%img_cc%get_cc_voxels(self%n_cc, first, vox)
            deallocate(atom_del_mask)
            allocate(atom_del_mask(self%n_cc))
            do icc = 1, self%n_cc
                atom_del_mask(icc) = euclid(self%atominfo(icc)%center(:), self%NPcen) * self%smpd > cendist_thres
            end do
        end subroutine relabel

        subroutine zero_component( label )
            integer, intent(in) :: label
            integer :: iv
            do iv = first(label), first(label+1) - 1
                imat_bin(vox(1,iv),vox(2,iv),vox(3,iv)) = 0
            enddo
        end subroutine zero_component

    end subroutine discard_atoms

    ! d_NN from the centres of the first binarisation, and the atom radius of detection from it
    subroutine set_nn_length_scales( self )
        class(nanoparticle), intent(inout) :: self
        real, allocatable :: centers_A(:,:)
        centers_A = self%atominfo2centers_A()
        self%d_nn = est_nn_dist(centers_A)
        write(logfhandle,'(A,F8.4)') 'nearest-neighbour distance d_NN (A): ', self%d_nn
        self%theoretical_radius = RADIUS_NN_FRAC * self%d_nn
        write(logfhandle,'(A,F8.4)') 'atom radius 0.4 d_NN (A):            ', self%theoretical_radius
        deallocate(centers_A)
    end subroutine set_nn_length_scales

    ! contact scores of discard_atoms, the neighbour cutoff from d_NN
    subroutine contact_scores( self, centers_A, cscores )
        class(nanoparticle), intent(in)    :: self
        real, allocatable,   intent(in)    :: centers_A(:,:)
        integer,             intent(inout) :: cscores(:)
        call calc_contact_scores(centers_A, cscores, self%d_nn)
    end subroutine contact_scores

    ! voxels that can hold only noise once the atoms are excluded: inside the box's inscribed sphere and inside any
    ! mask the map already carries (its nearest exact zero), both less two soft-edge widths
    subroutine empty_zone( self, nmap, zone )
        class(nanoparticle),  intent(in)  :: self
        real,                 intent(in)  :: nmap(:,:,:)
        logical, allocatable, intent(out) :: zone(:,:,:)
        real    :: cen(3), rlim, r2
        integer :: ii, jj, kk
        logical :: l_zero
        cen  = real(self%ldim / 2 + 1)
        rlim = real(minval(self%ldim)) / 2.
        do kk = 1,self%ldim(3)
            do jj = 1,self%ldim(2)
                do ii = 1,self%ldim(1)
                    l_zero = .not. abs(nmap(ii,jj,kk)) > 0.
                    if( l_zero ) rlim = min(rlim, sqrt(sum((real([ii,jj,kk]) - cen)**2)))
                enddo
            enddo
        enddo
        rlim = rlim - 2. * COSMSKHALFWIDTH
        allocate(zone(self%ldim(1),self%ldim(2),self%ldim(3)))
        do kk = 1,self%ldim(3)
            do jj = 1,self%ldim(2)
                do ii = 1,self%ldim(1)
                    r2 = sum((real([ii,jj,kk]) - cen)**2)
                    zone(ii,jj,kk) = r2 <= rlim * rlim
                enddo
            enddo
        enddo
    end subroutine empty_zone

    ! The density model of the atoms (atomic_density_model.md, sections 2, 4, 5 and 12), fitted on the map as read less
    ! the mask's constant: the provisional fit of the level-1 (strong) atoms, the residual recovery of the weak atoms
    ! and their pruning, the full fit; the transfer is measured in the full fit, or read from transfer_tab and held; fills
    ! atominfo (grown by the weak atoms); writes _SIM, _FIT, _RES, _atoms.csv, _model.txt, _transfer.tsv and
    ! _resid_extrema.csv
    subroutine fit_density_model( self )
        use simple_nano_density, only: nano_density, render_atoms, gauss_filter3D
        use simple_nano_species, only: calibrate_threshold, expected_false_count, local_maxima
        class(nanoparticle), intent(inout) :: self
        integer, parameter :: NREGION_MIN      = 20000 ! smallest noise region, voxels
        real,    parameter :: REGION_NN_FRAC   = 1.5   ! noise region: farther than this from every atom, in d_NN
        real,    parameter :: ENVELOPE_NN_FRAC = 1.    ! envelope: within this of an atom, in d_NN
        real,    parameter :: NOISELESS_FRAC   = 1.e-3 ! a region spread below this of the envelope's is no noise
        real,    parameter :: TARGET_FALSE     = 0.2   ! expected noise maxima above k_A in the stage A search volume
        real,    parameter :: CAL_LEVELS(3)    = [2.5, 3.0, 3.5]
        real,    parameter :: NMS_NN_FRAC      = 0.35  ! non-maximum suppression and centroid radius, in d_NN
        real,    parameter :: EXCL_A_NN_FRAC   = 0.7   ! stage A (and ungated stage B): no atom closer, in d_NN
        real,    parameter :: EXCL_B_NN_FRAC   = 0.85  ! gate: no atom closer, in d_NN
        real,    parameter :: GATE_NN_FRAC     = 1.15  ! gate: at least min_nbrs found atoms within, in d_NN
        integer, parameter :: MAX_LEVELS       = 8
        real,    parameter :: BANDS(2,3)       = reshape([5., 3.3, 2.5, 2.0, 2.0, 1.6], [2,3]) ! correlation bands, A
        class(nano_density), allocatable :: dens
        type(image)          :: img_tmp, img_band
        type(atom_stats)     :: blank
        real,    allocatable :: map(:,:,:), arr(:,:,:), resid(:,:,:), cmap(:,:,:), zmap(:,:,:), xyz(:,:), ints(:), bfacs(:)
        real,    allocatable :: vi(:), vb(:), tau_i(:), wpos(:,:), wdetz(:), aper(:)
        real,    allocatable :: mband(:,:,:,:)
        logical, allocatable :: zone(:,:,:), region(:,:,:), envelope(:,:,:), cal(:,:,:), sphere(:,:,:), smask(:,:,:)
        logical, allocatable :: wgate(:)
        integer, allocatable :: cn(:), bin_lo(:), wstage(:), wlevel(:)
        real,    allocatable :: t_meas(:), t_sm(:), pow_t(:), t_rep(:)
        real,    allocatable :: round_fields(:,:), round_curves(:,:)
        logical, allocatable :: t_rel(:)
        real    :: dv, sig_n, s_int, f_corr, ln_s, beta, ln_bp, gamma, tau_b, resid_ratio, sig_ref, ds_t, t_change
        real    :: gauge_scale, gauge_bfac
        real    :: z_mu, z_sd, k_a, k_b, c_search_a, c_region, efalse_a, efalse_b, corrs(3,4)
        integer :: n, n1, i, nsearch_a, nsearch_b, nadded(2,MAX_LEVELS), npruned, stage, level, nnew
        integer :: npeaks, nwells, n_held, n_gauge, round_kept
        logical :: l_region, l_resid, l_prof, l_measure, l_calibrated, l_rounds_capped
        integer(timer_int_kind) :: t_stage
        real(timer_int_kind)    :: rt_prov, rt_rec, rt_full, rt_prod
        write(logfhandle,'(A)') '>>> ATOMIC DENSITY MODEL'
        t_stage = tic()
        n1      = self%n_cc
        dv      = self%d_nn / self%smpd
        sig_ref = sqrt(B_REF_PSEUDO / (8. * PI**2))
        call img_tmp%new(self%ldim, self%smpd)
        call img_tmp%read(self%npname)
        ! the zone from the exact zeros of the map as read, then the map less the constant the mask took off
        map = img_tmp%get_rmat()
        call self%empty_zone(map, zone)
        call img_tmp%zero_background
        map = img_tmp%get_rmat()
        ! the search volume: the soft sphere new applied to the map
        allocate(arr(self%ldim(1),self%ldim(2),self%ldim(3)), source=1.)
        call img_tmp%set_rmat(arr, .false.)
        call img_tmp%mask3D_soft(self%msk_rad, backgr=0.)
        arr       = img_tmp%get_rmat()
        sphere    = arr > 0.
        nsearch_a = count(sphere)
        allocate(xyz(3,n1))
        do i = 1,n1
            xyz(:,i) = self%atominfo(i)%center
        enddo
        allocate(dens)
        call dens%new(map, self%smpd, self%d_nn, xyz)
        ! the transfer: measured in each fit, or the table of transfer_tab held throughout
        l_measure = .not. self%transfer_tab%is_allocated()
        if( .not. l_measure )then
            call read_transfer_tab(ds_t, t_sm)
            call dens%set_transfer(ds_t, t_sm)
            write(logfhandle,'(A)') 'transfer read from '//self%transfer_tab%to_char()
        endif
        l_region = .false.
        call build_masks
        l_region = count(region) >= NREGION_MIN
        if( l_region ) l_region = sdev_of(pack(map, region)) > NOISELESS_FRAC * sdev_of(pack(map, envelope))
        call build_masks
        allocate(resid(self%ldim(1),self%ldim(2),self%ldim(3)), cmap(self%ldim(1),self%ldim(2),self%ldim(3)),&
            &zmap(self%ldim(1),self%ldim(2),self%ldim(3)))
        allocate(wpos(3,0), wdetz(0), wgate(0), wstage(0), wlevel(0))
        ! the provisional fit of the strong atoms
        ! the provisional fit through the plain Gaussian, or the given transfer: a transfer measured on the strong atoms
        ! alone takes up the density of the atoms still missing where they share the lattice's frequencies, and the
        ! residual the recovery searches loses them (the light particle of species_discovery, phase 6); the transfer is
        ! measured once, after the recovery, on the complete model
        call set_noise
        call dens%fit_model
        call set_noise
        rt_prov = toc(t_stage)
        t_stage = tic()
        ! the residual recovery of the weak atoms and their pruning
        nadded  = 0
        npruned = 0
        call compute_zmap
        call calibrate
        if( l_calibrated )then
            do stage = 1,2
                do level = 1,MAX_LEVELS
                    call detect_residual_level(stage, level, nnew)
                    nadded(stage,level) = nnew
                    if( nnew == 0 ) exit
                enddo
            enddo
            call prune_weak
        endif
        rt_rec  = toc(t_stage)
        t_stage = tic()
        ! the transfer in rounds on the complete model (strong and weak atoms): measured, held fixed through a full fit,
        ! until the curve changes by 0.01 or less from one round to the next; the last round is kept unless clearly
        ! worse by the residual over the noise than the best, and measured again for the report only (_transfer.tsv,
        ! transfer_change); with transfer_tab the one full fit
        call build_masks
        call set_noise
        if( l_measure )then
            call dens%fit_transfer_rounds(envelope)
        else
            call dens%fit_model
        endif
        call dens%get_rounds(round_fields, round_curves, round_kept, l_rounds_capped)
        call build_masks
        call set_noise
        rt_full = toc(t_stage)
        t_stage = tic()
        ! the residual after the fit against the calibrated threshold: peaks are atoms still missing, wells atoms too many
        npeaks = 0
        nwells = 0
        if( l_calibrated )then
            call compute_zmap
            call count_extrema
        endif
        ! the atom table grows by the weak atoms; the model fields of every atom
        n = dens%get_natoms()
        call dens%get_atoms(xyz, ints, bfacs, vi, vb)
        cn = dens%get_coordination()
        self%atominfo = [self%atominfo(1:n1), spread(blank, 1, n - n1)]
        do i = n1+1,n
            self%atominfo(i)%center    = wpos(:,i-n1)
            self%atominfo(i)%det_stage = wstage(i-n1)
            self%atominfo(i)%det_level = wlevel(i-n1)
            self%atominfo(i)%det_z     = wdetz(i-n1)
            self%atominfo(i)%gate_used = wgate(i-n1)
        enddo
        allocate(aper(n))
        call dens%aperture_intensities(map, aper)
        do i = 1,n
            self%atominfo(i)%int_fit    = ints(i)
            self%atominfo(i)%bfac_fit   = bfacs(i)
            ! -1: no measurement spread (no intensity, or a singular fit)
            self%atominfo(i)%sd_lnint   = -1.
            self%atominfo(i)%sd_lnbfac  = -1.
            if( vi(i) < huge(1.) ) self%atominfo(i)%sd_lnint  = sqrt(vi(i))
            if( vb(i) < huge(1.) ) self%atominfo(i)%sd_lnbfac = sqrt(vb(i))
            self%atominfo(i)%center_fit = xyz(:,i)
            self%atominfo(i)%cn_fit     = cn(i)
            self%atominfo(i)%aper_int   = aper(i)
        enddo
        call dens%get_transfer(ds_t, t_meas, t_sm, pow_t, t_rel, t_rep, n_held, t_change)
        call dens%get_gauge(gauge_scale, gauge_bfac, n_gauge)
        call dens%get_noise(sig_n, s_int, f_corr, l_resid)
        self%int_noise     = s_int
        self%l_noise_resid = l_resid
        l_prof = dens%has_profiles()
        if( l_prof )then
            call dens%get_intensity_profile(ln_s, beta, bin_lo, tau_i)
            call dens%get_width_profile(ln_bp, gamma, tau_b)
            self%int_scale = exp(ln_s)
        else
            THROW_WARN('too few strong atoms for the profiles; the intensity scale is the median intensity')
            self%int_scale = median(ints)
        endif
        if( .not. self%int_scale > 0. ) self%int_scale = 1.
        ! products: the atoms without the transfer (and the weak atoms' valid_corr against it), the full model, the residual;
        ! the map correlations of the identical-atom render, _SIM and _FIT inside the mask, overall and by band
        ! single-threaded transforms: the system FFTW's threaded ones are slow above 16 threads on some machines
        call img_band%new(self%ldim, self%smpd, wthreads=.false.)
        allocate(mband(self%ldim(1),self%ldim(2),self%ldim(3),3))
        do i = 1,3
            mband(:,:,:,i) = band_passed(map, i)
        enddo
        call render_atoms(arr, self%smpd, reshape([(self%atominfo(i)%center, i=1,n)], [3,n]), spread(1., 1, n),&
            &spread(B_REF_PSEUDO, 1, n), 0., 0.)
        call map_corrs(arr, corrs(:,1), corrs(1,4))
        call dens%render_sim(arr)
        call map_corrs(arr, corrs(:,2), corrs(2,4))
        call weak_valid_corr
        call write_vol('_SIM.mrc', 'simulated atomic density: ')
        call dens%render_fit(arr)
        call map_corrs(arr, corrs(:,3), corrs(3,4))
        call write_vol('_FIT.mrc', 'density model of the map: ')
        resid_ratio = sqrt(sum((map - arr)**2, mask=envelope) / real(count(envelope))) / sig_n
        arr = map - arr
        call write_vol('_RES.mrc', 'residual of the model:    ')
        call write_atoms_csv
        call write_model_txt
        call write_transfer_tsv
        write(logfhandle,'(A,I0,A,I0,A,I0)') 'strong atoms: ', n1, ', weak atoms: ', n - n1, ', weak atoms pruned: ', npruned
        write(logfhandle,'(A,I0,A,F8.4,A,F8.4)') 'transfer shells held: ', n_held, ', change on the final model: ', t_change,&
            &', residual over noise in the envelope: ', resid_ratio
        rt_prod = toc(t_stage)
        write(logfhandle,'(A,4(A,F7.2))') 'density model wall times (s):', ' provisional fit ', rt_prov, ', recovery ', rt_rec,&
            &', full fit ', rt_full, ', products ', rt_prod
        call dens%kill
        call img_tmp%kill
        call img_band%kill
        write(logfhandle,'(A)') '>>> ATOMIC DENSITY MODEL, COMPLETED'

    contains

        ! the noise region (the zone farther than REGION_NN_FRAC d_NN from every atom), the envelope (within
        ! ENVELOPE_NN_FRAC d_NN of an atom) and the calibration voxels, on the current atoms
        subroutine build_masks
            real, allocatable :: pos(:,:), d1(:), d2(:), d3(:), d4(:)
            integer :: j
            call dens%get_atoms(pos, d1, d2, d3, d4)
            region = zone
            if( allocated(envelope) ) deallocate(envelope)
            allocate(envelope(self%ldim(1),self%ldim(2),self%ldim(3)), source=.false.)
            do j = 1,size(pos,2)
                call mark_ball(pos(:,j), REGION_NN_FRAC * dv, region, .false.)
                call mark_ball(pos(:,j), ENVELOPE_NN_FRAC * dv, envelope, .true.)
            enddo
            ! without an empty region the noise reference is the residual between the atoms (section 4, item 5)
            cal = envelope
            do j = 1,size(pos,2)
                call mark_ball(pos(:,j), 0.5 * dv, cal, .false.)
            enddo
            if( l_region ) cal = region
        end subroutine build_masks

        ! the noise from the empty region, or from the residual between the atoms when there is none
        subroutine set_noise
            if( l_region )then
                call dens%noise_from_region(region)
            else
                call dens%noise_from_residual(envelope)
            endif
        end subroutine set_noise

        ! the standardised matched-filter map of the residual: the residual filtered by the detection template,
        ! standardised over the calibration voxels
        subroutine compute_zmap
            real, allocatable :: vals(:)
            call dens%get_residual(resid)
            call gauss_filter3D(resid, sig_ref / self%smpd, cmap)
            vals = pack(cmap, cal)
            z_mu = sum(vals) / real(size(vals))
            z_sd = sdev_of(vals)
            zmap = 0.
            if( z_sd > 0. ) zmap = (cmap - z_mu) / z_sd
        end subroutine compute_zmap

        ! the thresholds from the maxima of the noise field over the calibration voxels (species section 3.2): k_A at an
        ! expected TARGET_FALSE false atoms in the search volume, k_B one unit lower
        subroutine calibrate
            real,    allocatable :: vals(:)
            integer, allocatable :: ijk(:,:)
            real    :: counts(size(CAL_LEVELS)), ratio
            integer :: l
            l_calibrated = .false.
            k_a = 0.
            k_b = 0.
            efalse_a = 0.
            efalse_b = 0.
            nsearch_b = 0
            if( .not. z_sd > 0. ) return
            call local_maxima(zmap, cal, CAL_LEVELS(1), 0, ijk, vals)
            do l = 1,size(CAL_LEVELS)
                counts(l) = real(count(vals > CAL_LEVELS(l)))
            enddo
            if( .not. sum(counts) > 0. )then
                THROW_WARN('no noise maxima above the calibration levels: no weak atoms are recovered')
                return
            endif
            ratio = real(nsearch_a) / real(count(cal))
            call calibrate_threshold(CAL_LEVELS, counts, ratio, TARGET_FALSE, k_a, c_search_a)
            c_region  = c_search_a / ratio
            k_b       = k_a - 1.
            call stage_b_mask
            nsearch_b = count(smask)
            efalse_a  = expected_false_count(c_search_a, k_a)
            efalse_b  = expected_false_count(c_region * real(nsearch_b) / real(count(cal)), k_b)
            l_calibrated = .true.
            write(logfhandle,'(A,3F8.1)') 'noise maxima above 2.5, 3.0, 3.5: ', counts
            write(logfhandle,'(A,2F8.3)') 'calibrated thresholds k_A, k_B:     ', k_a, k_b
        end subroutine calibrate

        ! stage B searches within GATE_NN_FRAC d_NN of the atoms when gated, else the whole search volume
        subroutine stage_b_mask
            real, allocatable :: pos(:,:), d1(:), d2(:), d3(:), d4(:)
            integer :: j
            if( self%min_nbrs > 0 )then
                call dens%get_atoms(pos, d1, d2, d3, d4)
                if( allocated(smask) ) deallocate(smask)
                allocate(smask(self%ldim(1),self%ldim(2),self%ldim(3)), source=.false.)
                do j = 1,size(pos,2)
                    call mark_ball(pos(:,j), GATE_NN_FRAC * dv, smask, .true.)
                enddo
                smask = smask .and. sphere
            else
                smask = sphere
            endif
        end subroutine stage_b_mask

        ! one residual level: local maxima of the z map above the stage threshold, at the centroid of the positive
        ! residual, accepted by the exclusion distance or the neighbour gate (species section 3.3); the accepted atoms
        ! join the model, each fitted on its residual
        subroutine detect_residual_level( stage, level, nnew )
            integer, intent(in)  :: stage, level
            integer, intent(out) :: nnew
            real,    allocatable :: vals(:), considered(:,:), dists(:), pos(:,:), newpos(:,:), d1(:), d2(:), d3(:), d4(:)
            integer, allocatable :: ijk(:,:)
            real    :: thres, p(3)
            integer :: ic, nc, nbefore
            logical :: l_gated, l_ok
            call build_masks
            call set_noise
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
            call dens%get_atoms(pos, d1, d2, d3, d4)
            nbefore = size(pos,2)
            allocate(considered(3,max(size(vals),1)), newpos(3,max(size(vals),1)))
            nc   = 0
            nnew = 0
            do ic = 1,size(vals)
                ! non-maximum suppression between the maxima of this level
                if( nc > 0 )then
                    if( any(sum((considered(:,:nc) - spread(real(ijk(:,ic)), 2, nc))**2, dim=1) < (NMS_NN_FRAC * dv)**2) ) cycle
                endif
                nc = nc + 1
                considered(:,nc) = real(ijk(:,ic))
                p = centroid(ijk(:,ic))
                dists = sqrt(sum((pos - spread(p, 2, size(pos,2)))**2, dim=1)) * self%smpd
                if( nnew > 0 ) dists = [dists, sqrt(sum((newpos(:,:nnew) - spread(p, 2, nnew))**2, dim=1)) * self%smpd]
                if( l_gated )then
                    ! the gate counts the atoms accepted before this level: one shell per level
                    l_ok = minval(dists) > EXCL_B_NN_FRAC * self%d_nn
                    if( l_ok ) l_ok = count(dists(:nbefore) <= GATE_NN_FRAC * self%d_nn) >= self%min_nbrs
                else
                    l_ok = minval(dists) > EXCL_A_NN_FRAC * self%d_nn
                endif
                if( .not. l_ok ) cycle
                nnew = nnew + 1
                newpos(:,nnew) = p
                wpos   = reshape([wpos, p], [3,size(wstage)+1])
                wdetz  = [wdetz, vals(ic)]
                wgate  = [wgate, l_gated]
                wstage = [wstage, stage]
                wlevel = [wlevel, level]
            enddo
            if( nnew > 0 ) call dens%add_atoms(newpos(:,:nnew))
            write(logfhandle,'(A,I1,A,I1,A,I0)') 'residual stage ', stage, ', level ', level, ': weak atoms added ', nnew
        end subroutine detect_residual_level

        ! centroid of the positive residual within NMS_NN_FRAC d_NN of a voxel
        function centroid( ijk ) result( pos )
            integer, intent(in) :: ijk(3)
            real    :: pos(3), w, wsum, c(3)
            integer :: lo(3), hi(3), ii, jj, kk
            c    = real(ijk)
            pos  = 0.
            wsum = 0.
            lo   = max(1, floor(c - NMS_NN_FRAC * dv))
            hi   = min(self%ldim, ceiling(c + NMS_NN_FRAC * dv))
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        if( sum((real([ii,jj,kk]) - c)**2) > (NMS_NN_FRAC * dv)**2 ) cycle
                        w    = max(resid(ii,jj,kk), 0.)
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

        ! the pruning of discard_atoms applied to the weak atoms (ruling of 2026-10-09 of the species work): in the outer
        ! zone of the level-1 set, a weak atom with fewer strong atoms within the contact cutoff than the threshold
        ! discard_atoms derived goes; weak atoms never vouch for each other, atoms inside the zone stay
        subroutine prune_weak
            real,    allocatable :: pos(:,:), d1(:), d2(:), d3(:), d4(:)
            logical, allocatable :: keep(:)
            real    :: rmax_c
            integer :: j, l, ncont, nw
            call dens%get_atoms(pos, d1, d2, d3, d4)
            nw = size(pos,2) - n1
            if( nw < 1 ) return
            rmax_c = contact_cutoff(self%d_nn)
            allocate(keep(size(pos,2)), source=.true.)
            do j = n1+1,size(pos,2)
                if( euclid(pos(:,j), self%prune_cen) * self%smpd <= self%prune_rad_thres ) cycle
                ncont = 0
                do l = 1,n1
                    if( euclid(pos(:,j), pos(:,l)) * self%smpd < rmax_c ) ncont = ncont + 1
                enddo
                if( ncont < self%prune_cs_thres ) keep(j) = .false.
            enddo
            npruned = count(.not. keep)
            if( npruned > 0 )then
                call dens%keep_atoms(keep)
                wpos   = reshape(pack(wpos, spread(keep(n1+1:), 1, 3)), [3,nw-npruned])
                wdetz  = pack(wdetz,  keep(n1+1:))
                wgate  = pack(wgate,  keep(n1+1:))
                wstage = pack(wstage, keep(n1+1:))
                wlevel = pack(wlevel, keep(n1+1:))
            endif
            write(logfhandle,'(A,I0)') 'weak atoms pruned by contact score: ', npruned
        end subroutine prune_weak

        ! residual maxima above k_A (atoms still missing) and minima below -k_A (atoms too many) in the search volume,
        ! with the detection's neighbour condition; their positions (A) and z in _resid_extrema.csv
        subroutine count_extrema
            real,    allocatable :: vals(:), nz(:,:,:)
            integer, allocatable :: ijk(:,:)
            integer :: funit, ios, l
            call fopen(funit, file=self%fbody//'_resid_extrema.csv', iostat=ios, status='replace')
            call fileiochk('simple_nanoparticle :: fit_density_model; ERROR when opening the _resid_extrema.csv file', ios)
            write(funit,'(A)') 'KIND,X,Y,Z,Z_SCORE'
            call local_maxima(zmap, sphere, k_a, NVOX_THRESH - 1, ijk, vals)
            npeaks = size(vals)
            do l = 1,npeaks
                write(funit,'(A,4(A,ES13.6))') 'peak', ',', real(ijk(1,l) - 1) * self%smpd, ',', real(ijk(2,l) - 1) * self%smpd,&
                    &',', real(ijk(3,l) - 1) * self%smpd, ',', vals(l)
            enddo
            nz = -zmap
            call local_maxima(nz, sphere, k_a, NVOX_THRESH - 1, ijk, vals)
            nwells = size(vals)
            do l = 1,nwells
                write(funit,'(A,4(A,ES13.6))') 'well', ',', real(ijk(1,l) - 1) * self%smpd, ',', real(ijk(2,l) - 1) * self%smpd,&
                    &',', real(ijk(3,l) - 1) * self%smpd, ',', -vals(l)
            enddo
            call fclose(funit)
        end subroutine count_extrema

        ! valid_corr of the weak atoms against the model's atoms (arr holds _SIM), in the window of validate_atoms
        subroutine weak_valid_corr
            type(image)       :: simimg
            real, allocatable :: pix1(:), pix2(:)
            real    :: maxrad
            integer :: j, winsz, npix_in, nout1, nout2, ijk(3)
            if( n == n1 ) return
            call simimg%new(self%ldim, self%smpd)
            call simimg%set_rmat(arr, .false.)
            maxrad  = 1.5 * RADIUS_NN_FRAC * dv
            winsz   = ceiling(maxrad)
            npix_in = (2 * winsz + 1)**3
            allocate(pix1(npix_in), pix2(npix_in), source=0.)
            do j = n1+1,n
                ijk = nint(self%atominfo(j)%center)
                call self%img_raw%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, nout1, pix1)
                call simimg%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, nout2, pix2)
                self%atominfo(j)%valid_corr = pearsn_serial(pix1(:nout1), pix2(:nout2))
            enddo
            call simimg%kill
        end subroutine weak_valid_corr

        ! correlation of a render with the map inside the search volume, in the bands BANDS (band-passed both) and overall
        subroutine map_corrs( render, cband, coverall )
            real, intent(in)  :: render(:,:,:)
            real, intent(out) :: cband(3), coverall
            integer :: ib_
            coverall = corr_in(map, render)
            do ib_ = 1,3
                cband(ib_) = corr_in(mband(:,:,:,ib_), band_passed(render, ib_))
            enddo
        end subroutine map_corrs

        ! a3 band-passed to band ib_ of BANDS
        function band_passed( a3, ib_ ) result( b3 )
            real,    intent(in) :: a3(:,:,:)
            integer, intent(in) :: ib_
            real, allocatable   :: b3(:,:,:)
            call img_band%set_rmat(a3, .false.)
            call img_band%bp(BANDS(1,ib_), BANDS(2,ib_), width=1.)
            b3 = img_band%get_rmat()
        end function band_passed

        real function corr_in( a, b )
            real, intent(in) :: a(:,:,:), b(:,:,:)
            real(dp) :: ma, mb, saa, sbb, sab
            integer  :: m
            m   = count(sphere)
            ma  = sum(real(a, dp), mask=sphere) / real(m, dp)
            mb  = sum(real(b, dp), mask=sphere) / real(m, dp)
            saa = sum((real(a, dp) - ma)**2, mask=sphere)
            sbb = sum((real(b, dp) - mb)**2, mask=sphere)
            sab = sum((real(a, dp) - ma) * (real(b, dp) - mb), mask=sphere)
            corr_in = 0.
            if( saa > 0._dp .and. sbb > 0._dp ) corr_in = real(sab / sqrt(saa * sbb))
        end function corr_in

        real function sdev_of( v )
            real, intent(in) :: v(:)
            real :: m
            sdev_of = 0.
            if( size(v) < 2 ) return
            m       = sum(v) / real(size(v))
            sdev_of = sqrt(sum((v - m)**2) / real(size(v) - 1))
        end function sdev_of

        ! set arr3 to val within rv voxels of c
        subroutine mark_ball( c, rv, arr3, val )
            real,    intent(in)    :: c(3), rv
            logical, intent(inout) :: arr3(:,:,:)
            logical, intent(in)    :: val
            integer :: lo(3), hi(3), ii, jj, kk
            lo = max(1, floor(c - rv))
            hi = min(self%ldim, ceiling(c + rv))
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        if( sum((real([ii,jj,kk]) - c)**2) <= rv * rv ) arr3(ii,jj,kk) = val
                    enddo
                enddo
            enddo
        end subroutine mark_ball

        subroutine write_vol( suffix, what )
            character(len=*), intent(in) :: suffix, what
            call img_tmp%new(self%ldim, self%smpd)
            call img_tmp%set_rmat(arr, .false.)
            call img_tmp%write(self%fbody//suffix)
            write(logfhandle,'(A)') 'output, '//what//' '//self%fbody%to_char()//suffix
        end subroutine write_vol

        ! one row per atom in the order of _ATMS.pdb, which keeps the detection positions; here the fitted ones, in A
        subroutine write_atoms_csv
            integer :: funit, ios, j
            character(len=6) :: label
            call fopen(funit, file=self%fbody//'_atoms.csv', iostat=ios, status='replace')
            call fileiochk('simple_nanoparticle :: fit_density_model; ERROR when opening the _atoms.csv file', ios)
            write(funit,'(A)') 'INDEX,X,Y,Z,INT,BFAC,CN,SD_LNINT,SD_LNBFAC,OCC,VALID_CORR,LABEL,STAGE,LEVEL,DET_Z,GATE_USED,APER_INT'
            do j = 1,n
                label = 'strong'
                if( j > n1 ) label = 'weak'
                write(funit,'(I0,5(A,ES13.6),A,I0,4(A,ES13.6),A,A,2(A,I0),A,ES13.6,A,I0,A,ES13.6)') j,&
                    &',', (self%atominfo(j)%center_fit(1) - 1.) * self%smpd, ',', (self%atominfo(j)%center_fit(2) - 1.) * self%smpd,&
                    &',', (self%atominfo(j)%center_fit(3) - 1.) * self%smpd, ',', self%atominfo(j)%int_fit,&
                    &',', self%atominfo(j)%bfac_fit, ',', self%atominfo(j)%cn_fit, ',', self%atominfo(j)%sd_lnint,&
                    &',', self%atominfo(j)%sd_lnbfac, ',', self%atominfo(j)%int_fit / self%int_scale,&
                    &',', self%atominfo(j)%valid_corr, ',', trim(label), ',', self%atominfo(j)%det_stage,&
                    &',', self%atominfo(j)%det_level, ',', self%atominfo(j)%det_z, ',', merge(1, 0, self%atominfo(j)%gate_used),&
                    &',', self%atominfo(j)%aper_int
            enddo
            call fclose(funit)
            write(logfhandle,'(A)') 'output, per-atom model table:     '//self%fbody%to_char()//'_atoms.csv'
        end subroutine write_atoms_csv

        ! the particle's parameters, key = value
        subroutine write_model_txt
            character(len=*), parameter :: SETS(3) = ['equal', 'sim  ', 'fit  ']
            character(len=*), parameter :: BAND_TAGS(3) = ['5.0-3.3', '2.5-2.0', '2.0-1.6']
            integer :: funit, ios, k, l
            call fopen(funit, file=self%fbody//'_model.txt', iostat=ios, status='replace')
            call fileiochk('simple_nanoparticle :: fit_density_model; ERROR when opening the _model.txt file', ios)
            call kv(funit, 'S', self%int_scale)
            if( l_prof )then
                call kv(funit, 'beta', beta)
                call kv(funit, 'gamma', gamma)
                call kv(funit, 'B_p', exp(ln_bp))
                call kv(funit, 'tau_B', tau_b)
                write(funit,'(A,I0)') 'tau_I_bins = ', size(tau_i)
                do k = 1,size(tau_i)
                    write(funit,'(A,I0,A,I0)') 'tau_I_bin_', k, '_cn_min = ', bin_lo(k)
                    call kv(funit, 'tau_I_bin_'//int2str(k), tau_i(k))
                enddo
            endif
            write(funit,'(A,A)') 'profiles = ', merge('yes', 'no ', l_prof)
            write(funit,'(A,A)') 'transfer = ', merge('measured', 'file    ', l_measure)
            write(funit,'(A,I0)') 'transfer_shells_held = ', n_held
            call kv(funit, 'transfer_change', t_change)
            ! the rounds of the measured transfer, the one kept, and per round the residual over the noise (the noise as
            ! set before the rounds) and the curve's change from the previous round
            write(funit,'(A,I0)') 'transfer_rounds = ', size(round_fields, 2)
            write(funit,'(A,I0)') 'transfer_round = ', round_kept
            write(funit,'(A,A)') 'transfer_rounds_capped = ', merge('yes', 'no ', l_rounds_capped)
            do k = 1,size(round_fields, 2)
                call kv(funit, 'transfer_round_'//int2str(k)//'_resid_over_noise', round_fields(1,k))
                if( k > 1 ) call kv(funit, 'transfer_round_'//int2str(k)//'_curve_change', round_fields(6,k))
            enddo
            write(funit,'(A,I0)') 'transfer_gauge_shells = ', n_gauge
            call kv(funit, 'transfer_gauge_scale', gauge_scale)
            call kv(funit, 'transfer_gauge_bfac', gauge_bfac)
            call kv(funit, 'sigma_n', sig_n)
            write(funit,'(A,A)') 'noise_source = ', merge('residual', 'region  ', l_resid)
            call kv(funit, 's_int', s_int)
            call kv(funit, 'f_corr', f_corr)
            call kv(funit, 'resid_over_noise', resid_ratio)
            call kv(funit, 'd_NN', self%d_nn)
            write(funit,'(A,I0)') 'strong_atoms = ', n1
            write(funit,'(A,I0)') 'weak_atoms = ', n - n1
            write(funit,'(A,A)') 'recovery_calibrated = ', merge('yes', 'no ', l_calibrated)
            call kv(funit, 'k_A', k_a)
            call kv(funit, 'k_B', k_b)
            call kv(funit, 'expected_false_A', efalse_a)
            call kv(funit, 'expected_false_B', efalse_b)
            write(funit,'(A,I0)') 'search_voxels_A = ', nsearch_a
            write(funit,'(A,I0)') 'search_voxels_B = ', nsearch_b
            write(funit,'(A,I0)') 'min_nbrs = ', self%min_nbrs
            do l = 1,2
                do k = 1,MAX_LEVELS
                    write(funit,'(A,A,A,I0,A,I0)') 'added_', merge('A', 'B', l == 1), '_level_', k, ' = ', nadded(l,k)
                enddo
            enddo
            write(funit,'(A,I0)') 'weak_pruned = ', npruned
            call kv(funit, 'prune_zone_radius', self%prune_rad_thres)
            write(funit,'(A,I0)') 'prune_contact_threshold = ', self%prune_cs_thres
            write(funit,'(A,I0)') 'contact_ceiling = ', self%cscore_ceil
            write(funit,'(A,I0)') 'residual_peaks = ', npeaks
            write(funit,'(A,I0)') 'residual_wells = ', nwells
            do k = 1,3
                call kv(funit, 'corr_'//trim(SETS(k))//'_all', corrs(k,4))
                do l = 1,3
                    call kv(funit, 'corr_'//trim(SETS(k))//'_'//BAND_TAGS(l), corrs(l,k))
                enddo
            enddo
            call fclose(funit)
            write(logfhandle,'(A)') 'output, particle model parameters: '//self%fbody%to_char()//'_model.txt'
        end subroutine write_model_txt

        ! the transfer per Fourier shell: s range (1/A), resolution range (A), as measured (0 without a measurement),
        ! as used by the fit, the atoms' power, whether the shell is reliable, and the report's measurement on the final
        ! model (0 without); as used (smoothed and held), the atoms' power sum |F_S|**2 and whether the shell is reliable
        subroutine write_transfer_tsv
            integer :: funit, ios, k
            real    :: tm, pw, tr
            logical :: l_rel
            call fopen(funit, file=self%fbody//'_transfer.tsv', iostat=ios, status='replace')
            call fileiochk('simple_nanoparticle :: fit_density_model; ERROR when opening the _transfer.tsv file', ios)
            write(funit,'(A)') 's_lo'//char(9)//'s_hi'//char(9)//'res_lo_A'//char(9)//'res_hi_A'//char(9)//'T_measured'//&
                &char(9)//'T_smoothed'//char(9)//'power_S'//char(9)//'reliable'//char(9)//'T_final_model'
            do k = 1,size(t_sm)
                tm    = 0.
                pw    = 0.
                l_rel = .false.
                tr    = 0.
                if( k <= size(t_meas) )then
                    tm    = t_meas(k)
                    pw    = pow_t(k)
                    l_rel = t_rel(k)
                endif
                if( k <= size(t_rep) ) tr = t_rep(k)
                write(funit,'(F8.4,A,F8.4,A,F10.3,A,F10.3,A,ES12.4,A,ES12.4,A,ES12.4,A,A,A,ES12.4)') real(k - 1) * ds_t, char(9),&
                    &real(k) * ds_t, char(9), 1. / max(real(k - 1) * ds_t, 1.e-3), char(9), 1. / (real(k) * ds_t), char(9),&
                    &tm, char(9), t_sm(k), char(9), pw, char(9), merge('yes', 'no ', l_rel), char(9), tr
            enddo
            call fclose(funit)
            write(logfhandle,'(A)') 'output, transfer of the map:      '//self%fbody%to_char()//'_transfer.tsv'
        end subroutine write_transfer_tsv

        ! the smoothed transfer of a _transfer.tsv (columns s_lo, s_hi, ..., T_smoothed), its shell width
        subroutine read_transfer_tab( ds, t )
            real,              intent(out) :: ds
            real, allocatable, intent(out) :: t(:)
            character(len=512) :: line
            real    :: vals(6), tbuf(10000)
            integer :: funit, ios, m
            call fopen(funit, file=self%transfer_tab, status='old', action='read', iostat=ios)
            call fileiochk('simple_nanoparticle :: fit_density_model; ERROR when opening transfer_tab', ios)
            read(funit,'(A)', iostat=ios) line
            m  = 0
            ds = 0.
            do
                read(funit,'(A)', iostat=ios) line
                if( ios /= 0 ) exit
                if( len_trim(line) == 0 ) cycle
                read(line,*, iostat=ios) vals
                if( ios /= 0 ) THROW_HARD('transfer_tab: a row does not start with six numbers (s_lo s_hi res_lo res_hi T_measured T_smoothed)')
                m = m + 1
                if( m > size(tbuf) ) THROW_HARD('transfer_tab: too many shells')
                tbuf(m) = vals(6)
                if( m == 1 ) ds = vals(2) - vals(1)
            enddo
            call fclose(funit)
            if( m < 1 .or. .not. ds > 0. ) THROW_HARD('transfer_tab holds no shell')
            t = tbuf(:m)
        end subroutine read_transfer_tab

        subroutine kv( funit, key, val )
            integer,          intent(in) :: funit
            character(len=*), intent(in) :: key
            real,             intent(in) :: val
            write(funit,'(A,A,ES14.6)') key, ' = ', val
        end subroutine kv

    end subroutine fit_density_model

    subroutine fillin_atominfo( self, a0, imat )
        class(nanoparticle),        intent(inout) :: self
        real,             optional, intent(in)    :: a0(3) ! lattice parameters
        integer,          optional, intent(in)    :: imat(:,:,:)
        type(image)          :: simatms, fit_isotropic, fit_anisotropic
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:), strain_array(:,:)
        real,    pointer     :: rmat_raw(:,:,:)
        integer, allocatable :: imat_cc(:,:,:), first(:), vox(:,:)
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
        call self%img_cc%get_cc_voxels(self%n_cc, first, vox)
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
            call self%calc_longest_atm_dist(cc, vox(:,first(cc):first(cc+1)-1), self%atominfo(cc)%diam)
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

    ! largest distance from the centre of atom label to the voxels pos(3,:) of its connected component
    subroutine calc_longest_atm_dist( self, label, pos, longest_dist )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: label
        integer,             intent(in)    :: pos(:,:)
        real,                intent(out)   :: longest_dist
        logical, allocatable :: mask_dist(:) ! for min and max dist calculation
        integer :: location(1)               ! location of vxls of the atom farthest from its center
        allocate(mask_dist(size(pos, dim=2)), source = .true.)
        if( size(pos,2) == 1 ) then ! if the connected component has size 1 (just 1 vxl)
            longest_dist  = self%smpd
        else
            longest_dist  = pixels_dist(self%atominfo(label)%center(:), real(pos),'max', mask_dist, location) * self%smpd
        endif
        deallocate(mask_dist)
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

    subroutine simulate_atoms( self, simatms, betas, mask, atoms_obj, l_pseudo )
        class(nanoparticle),           intent(inout) :: self
        class(image),                  intent(out)   :: simatms
        real,        optional,         intent(in)    :: betas(self%n_cc) ! in pdb file b-factor
        logical,     optional,         intent(in)    :: mask(self%n_cc)
        type(atoms), optional, target, intent(inout) :: atoms_obj
        logical,     optional,         intent(in)    :: l_pseudo         ! detection: pseudo-atoms at the template width
        type(atoms), target  :: atms_here
        type(atoms), pointer :: atms_ptr => null()
        logical              :: betas_present, mask_present, atoms_obj_present, l_pseudo_here
        integer              :: i, cnt
        l_pseudo_here = self%l_species_free
        if( present(l_pseudo) ) l_pseudo_here = l_pseudo_here .or. l_pseudo
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
                if( l_pseudo_here )then
                    call atms_ptr%set_element(atms_obj_ind, 'X1')
                    call atms_ptr%set_beta(atms_obj_ind, B_REF_PSEUDO)
                endif
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
        else if( present(fname) )then
            call centers_pdb%new(self%n_cc, dummy=.true.)
            do cc = 1, self%n_cc
                call set_atm_info
            enddo
        else
            ! every atom of the density model, strong and weak
            call centers_pdb%new(size(self%atominfo), dummy=.true.)
            do cc = 1, size(self%atominfo)
                call set_atm_info
                call set_model_info
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

            ! _ATMS.pdb: occupancy I / S and B of the density model; chain A the strong atoms, B the weak ones
            subroutine set_model_info
                call centers_pdb%set_beta(cc, self%atominfo(cc)%bfac_fit)
                call centers_pdb%set_occupancy(cc, self%atominfo(cc)%int_fit / self%int_scale)
                if( self%atominfo(cc)%det_stage > 0 )then
                    call centers_pdb%set_chain(cc, 'B')
                else
                    call centers_pdb%set_chain(cc, 'A')
                endif
            end subroutine set_model_info

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
        self%l_species_free     = .false.
        self%d_nn               = 0.
        self%int_scale          = 0.
        self%int_noise          = 0.
        self%l_noise_resid      = .false.
        self%prune_cen          = 0.
        self%prune_rad_thres    = 0.
        self%prune_cs_thres     = 0
        self%cscore_ceil        = 0
        self%min_nbrs           = 3
        self%msk_rad            = 0.
        call self%transfer_tab%kill
    end subroutine kill

end module simple_nanoparticle

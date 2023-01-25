module simple_nanoparticle
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'

use simple_image,      only: image
use simple_binimage,   only: binimage
use simple_atoms,      only: atoms
use simple_neighs,     only: neigh_4_3D
use simple_parameters, only: params_glob
use simple_qr_solve
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
integer,          parameter :: NPARAMS_ADP         = 6       ! min # voxels for ADP calc. = # params in cov.
integer,          parameter :: SCALE_FAC           = 4       ! scale factor for ADP calc
logical,          parameter :: DEBUG               = .false. ! for debugging purposes
logical,          parameter :: WRITE_OUTPUT        = .false. ! for figures generation
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
&'CN_STD'//CSV_DELIM//'NN_BONDL'//CSV_DELIM//'CN_GEN'//CSV_DELIM//'DIAM'//CSV_DELIM//'AVG_INT'//&
&CSV_DELIM//'MAX_INT'//CSV_DELIM//'CENDIST'//CSV_DELIM//'VALID_CORR'//CSV_DELIM//CSV_DELIM//&
&'DISPL'//CSV_DELIM//'MAX_NDISPL'//CSV_DELIM//'SEMIAX_MAJ'//CSV_DELIM//'SEMIAX_MED'//CSV_DELIM//&
&'SEMIAX_MIN'//CSV_DELIM//'ANISO_X'//CSV_DELIM//'ANISO_Y'//CSV_DELIM//'ANISO_Z'//'X'//CSV_DELIM//&
&'Y'//CSV_DELIM//'Z'//CSV_DELIM//'EXX_STRAIN'//CSV_DELIM//'EYY_STRAIN'//CSV_DELIM//'EZZ_STRAIN'//&
&CSV_DELIM//'EXY_STRAIN'//CSV_DELIM//'EYZ_STRAIN'//CSV_DELIM//'EXZ_STRAIN'//CSV_DELIM//'RADIAL_STRAIN'

character(len=*), parameter :: ATOM_STATS_HEAD_OMIT = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//&
&'CN_STD'//CSV_DELIM//'NN_BONDL'//CSV_DELIM//'CN_GEN'//CSV_DELIM//'DIAM'//CSV_DELIM//'AVG_INT'//&
&CSV_DELIM//'MAX_INT'//CSV_DELIM//'CENDIST'//CSV_DELIM//'VALID_CORR'//CSV_DELIM//'DISPL'//CSV_DELIM//&
&'MAX_NDISPL'//CSV_DELIM//'SEMIAX_MAJ'//CSV_DELIM//'SEMIAX_MED'//CSV_DELIM//'SEMIAX_MIN'//CSV_DELIM//&
&'ANISO_X'//CSV_DELIM//'ANISO_Y'//CSV_DELIM//'ANISO_Z'//CSV_DELIM//'RADIAL_STRAIN'

character(len=*), parameter :: NP_STATS_HEAD = 'NATOMS'//CSV_DELIM//'NANISO'//CSV_DELIM//'DIAM'//&
&CSV_DELIM//'AVG_NVOX'//CSV_DELIM//'MED_NVOX'//CSV_DELIM//'SDEV_NVOX'//&
&CSV_DELIM//'AVG_CN_STD'//CSV_DELIM//'MED_CN_STD'//CSV_DELIM//'SDEV_CN_STD'//&
&CSV_DELIM//'AVG_NN_BONDL'//CSV_DELIM//'MED_NN_BONDL'//CSV_DELIM//'SDEV_NN_BONDL'//&
&CSV_DELIM//'AVG_CN_GEN'//CSV_DELIM//'MED_CN_GEN'//CSV_DELIM//'SDEV_CN_GEN'//&
&CSV_DELIM//'AVG_DIAM'//CSV_DELIM//'MED_DIAM'//CSV_DELIM//'SDEV_DIAM'//&
&CSV_DELIM//'AVG_AVG_INT'//CSV_DELIM//'MED_AVG_INT'//CSV_DELIM//'SDEV_AVG_INT'//&
&CSV_DELIM//'AVG_MAX_INT'//CSV_DELIM//'MED_MAX_INT'//CSV_DELIM//'SDEV_MAX_INT'//&
&CSV_DELIM//'AVG_VALID_CORR'//CSV_DELIM//'MED_VALID_CORR'//CSV_DELIM//'SDEV_VALID_CORR'//&
&CSV_DELIM//'AVG_DISPL'//CSV_DELIM//'MED_DISPL'//CSV_DELIM//'SDEV_DISPL'//&
&CSV_DELIM//'AVG_MAX_NDISPL'//CSV_DELIM//'MED_MAX_NDISPL'//CSV_DELIM//'SDEV_MAX_NDISPL'//&
&CSV_DELIM//'AVG_SEMIAX_MAJ'//CSV_DELIM//'MED_SEMIAX_MAJ'//CSV_DELIM//'SDEV_SEMIAX_MAJ'//&
&CSV_DELIM//'AVG_SEMIAX_MED'//CSV_DELIM//'MED_SEMIAX_MED'//CSV_DELIM//'SDEV_SEMIAX_MED'//&
&CSV_DELIM//'AVG_SEMIAX_MIN'//CSV_DELIM//'MED_SEMIAX_MIN'//CSV_DELIM//'SDEV_SEMIAX_MIN'//&
&CSV_DELIM//'AVG_ANISO_X'//CSV_DELIM//'MED_ANISO_X'//CSV_DELIM//'SDEV_ANISO_X'//&
&CSV_DELIM//'AVG_ANISO_Y'//CSV_DELIM//'MED_ANISO_Y'//CSV_DELIM//'SDEV_ANISO_Y'//&
&CSV_DELIM//'AVG_ANISO_Z'//CSV_DELIM//'MED_ANISO_Z'//CSV_DELIM//'SDEV_ANISO_Z'//&
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
&CSV_DELIM//'AVG_DISPL'//CSV_DELIM//'MED_DISPL'//CSV_DELIM//'SDEV_DISPL'//&
&CSV_DELIM//'AVG_MAX_NDISPL'//CSV_DELIM//'MED_MAX_NDISPL'//CSV_DELIM//'SDEV_MAX_NDISPL'//&
&CSV_DELIM//'AVG_SEMIAX_MAJ'//CSV_DELIM//'MED_SEMIAX_MAJ'//CSV_DELIM//'SDEV_SEMIAX_MAJ'//&
&CSV_DELIM//'AVG_SEMIAX_MED'//CSV_DELIM//'MED_SEMIAX_MED'//CSV_DELIM//'SDEV_SEMIAX_MED'//&
&CSV_DELIM//'AVG_SEMIAX_MIN'//CSV_DELIM//'MED_SEMIAX_MIN'//CSV_DELIM//'SDEV_SEMIAX_MIN'//&
&CSV_DELIM//'AVG_ANISO_X'//CSV_DELIM//'MED_ANISO_X'//CSV_DELIM//'SDEV_ANISO_X'//&
&CSV_DELIM//'AVG_ANISO_Y'//CSV_DELIM//'MED_ANISO_Y'//CSV_DELIM//'SDEV_ANISO_Y'//&
&CSV_DELIM//'AVG_ANISO_Z'//CSV_DELIM//'MED_ANISO_Z'//CSV_DELIM//'SDEV_ANISO_Z'//&
&CSV_DELIM//'AVG_RADIAL_STRAIN'//CSV_DELIM//'MED_RADIAL_STRAIN'//CSV_DELIM//'SDEV_RADIAL_STRAIN'//&
&CSV_DELIM//'MIN_RADIAL_STRAIN'//CSV_DELIM//'MAX_RADIAL_STRAIN'

! container for per-atom statistics
type :: atom_stats
    ! various per-atom parameters                                                                   ! csv file labels
    integer :: cc_ind            = 0  ! index of the connected component                            INDEX
    integer :: size              = 0  ! number of voxels in connected component                     NVOX
    integer :: cn_std            = 0  ! standard coordination number                                CN_STD
    integer :: niso              = 0  ! Number of voxels in isotropic displ calculations            NISO
    real    :: bondl             = 0. ! nearest neighbour bond lenght in A                          NN_BONDL
    real    :: cn_gen            = 0. ! generalized coordination number                             CN_GEN
    real    :: diam              = 0. ! atom diameter                                               DIAM
    real    :: avg_int           = 0. ! average grey level intensity across the connected component AVG_INT
    real    :: max_int           = 0. ! maximum            -"-                                      MAX_INT
    real    :: cendist           = 0. ! distance from the centre of mass of the nanoparticle        CENDIST
    real    :: valid_corr        = 0. ! per-atom correlation with the simulated map                 VALID_CORR
    real    :: semiaxes(3)       = 0. ! lengths of semiaxes of ellipsoidal fit                      SEMIAX_(MAJ,MED,MIN)
    real    :: aniso_xyz(3)      = 0. ! anisotropic displacement in x,y,z directions                ANISO_()
    real    :: displ             = 0. ! Lattice displacement                                        DISPL
    real    :: max_ndispl        = 0. ! Maximum lattice displacement of neighboring atoms           NDISPL
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
    real    :: doa               = 0. ! Legacy (will clean up)

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
    type(stats_struct)    :: displ_stats
    type(stats_struct)    :: max_ndispl_stats
    type(stats_struct)    :: semiaxis_maj_stats
    type(stats_struct)    :: semiaxis_med_stats
    type(stats_struct)    :: semiaxis_min_stats
    type(stats_struct)    :: aniso_x_stats
    type(stats_struct)    :: aniso_y_stats
    type(stats_struct)    :: aniso_z_stats    
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
    type(stats_struct)    :: semiaxis_maj_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: semiaxis_med_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: semiaxis_min_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: aniso_x_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: aniso_y_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: aniso_z_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: displ_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: max_ndispl_stats_cns(CNMIN:CNMAX)
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
    procedure, private :: calc_aniso
    procedure, private :: lattice_displ_analysis
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
        call self%validate_atoms(simatms)
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
        type(image)          :: pc
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
        real    :: maxrad, corrs(self%n_cc)
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
        do cc = 1, self%n_cc
            if( count(imat_cc == cc) < NVOX_THRESH )then
                where(imat_cc == cc) imat_bin = 0
            endif
        end do
        ! Removing outliers based on coordination number
        n_discard = 0
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
        do cc = 1, self%n_cc
            if( count(imat_cc == cc) < NVOX_THRESH )then
                where(imat_cc == cc) imat_bin = 0
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
        type(image)          :: simatms
        type(binimage)       :: img_cc_scaled
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:), strain_array(:,:), lattice_displ(:,:)
        real,    pointer     :: rmat_raw(:,:,:)
        integer, allocatable :: imat_cc(:,:,:), imat_cc_scaled(:,:,:), neigh_4_pixs(:)
        character(len=256)   :: io_msg
        logical, allocatable :: cc_mask(:), displ_neighbor(:), border(:,:,:)
        logical, parameter   :: test_fit = .true.
        real    :: tmp_diam, a(3), res_fsc05, res_fsc0143
        integer :: i, j, k, cc, cn, n, m, l, x, y, z, ios, nsz, funit
        character(*), parameter :: fn_img_cc_scaled="cc_map_scaled.mrc", fn_muA="adp_info.txt"
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
        allocate(strain_array(self%n_cc,NSTRAIN_COMPS), lattice_displ(self%n_cc, 3), source=0.)
        call strain_analysis(self%element, centers_A, a, strain_array, lattice_displ)
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

        ! Create scaled CC map from unscale CC map (i,j,k)->(x,y,z) for aniso calculations
        call img_cc_scaled%new_bimg(SCALE_FAC*self%ldim, self%smpd / SCALE_FAC)
        allocate(imat_cc_scaled(SCALE_FAC*self%ldim(1), SCALE_FAC*self%ldim(2), SCALE_FAC*self%ldim(3)), source=0)
        do k=1, self%ldim(3)
            do j=1, self%ldim(2)
                do i=1, self%ldim(1)
                    do n=0, SCALE_FAC-1
                        x = i*SCALE_FAC - n
                        do m=0, SCALE_FAC-1 
                            y = j*SCALE_FAC - m
                            do l=0, SCALE_FAC-1
                                z=k*SCALE_FAC - l
                                imat_cc_scaled(x,y,z) = imat_cc(i,j,k)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call img_cc_scaled%set_imat(imat_cc_scaled)
        write(logfhandle, '(A, i3)') "ADP CALCULATIONS: CC MAP SCALED BY ", SCALE_FAC**3

        ! Find the surfaces of each atom (We assume atomic surfaces don't overlap)
        allocate(border(SCALE_FAC*self%ldim(1), SCALE_FAC*self%ldim(2), SCALE_FAC*self%ldim(3)), source=.false.)
        allocate(neigh_4_pixs(6),  source=0) ! There are 6 neighbors since we don't count diagonals
        do k=1, SCALE_FAC*self%ldim(3)
            do j=1, SCALE_FAC*self%ldim(2)
                do i=1, SCALE_FAC*self%ldim(1)
                    if(imat_cc_scaled(i,j,k) /= 0) then
                        call neigh_4_3D(SCALE_FAC*self%ldim, imat_cc_scaled, [i,j,k], neigh_4_pixs, nsz)
                        if(any(neigh_4_pixs(:nsz)<1)) border(i,j,k) = .true.
                    endif
                enddo
            enddo
        enddo
        deallocate(neigh_4_pixs)
        write(logfhandle, '(A, i3)') "ADP CALCULATIONS: IDENTIFIED ATOMIC BORDERS"
        call fopen(funit, FILE=trim(fn_muA), STATUS='REPLACE', action='WRITE')

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
            ! maximum grey level intensity across the connected component
            self%atominfo(cc)%max_int = maxval(rmat_raw(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            ! average grey level intensity across the connected component
            self%atominfo(cc)%avg_int = sum(rmat_raw(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            self%atominfo(cc)%avg_int = self%atominfo(cc)%avg_int / real(count(mask))
            ! bond length of nearest neighbour...
            self%atominfo(cc)%bondl = pixels_dist(self%atominfo(cc)%center(:), tmpcens, 'min', mask=cc_mask) ! Use all the atoms
            self%atominfo(cc)%bondl = self%atominfo(cc)%bondl * self%smpd ! convert to A
            ! Lattice displacement magnitudes ( |center - expected center| ) and max neighboring lattice displ
            call self%lattice_displ_analysis(cc, centers_A, a, lattice_displ)
            ! calculate anisotropic displacement parameters.  
            call self%calc_aniso(cc, imat_cc_scaled, border)

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
        call img_cc_scaled%write_bimg(fn_img_cc_scaled)
        call fclose(funit)

        ! CALCULATE GLOBAL NP PARAMETERS
        call calc_stats(  real(self%atominfo(:)%size),    self%size_stats, mask=self%atominfo(:)%size >= NVOX_THRESH )
        call calc_stats(  real(self%atominfo(:)%cn_std),  self%cn_std_stats        )
        call calc_stats(  self%atominfo(:)%bondl,         self%bondl_stats         )
        call calc_stats(  self%atominfo(:)%cn_gen,        self%cn_gen_stats        )
        call calc_stats(  self%atominfo(:)%diam,          self%diam_stats, mask=self%atominfo(:)%size >= NVOX_THRESH )
        call calc_zscore( self%atominfo(:)%avg_int ) ! to get comparable intensities between different particles
        call calc_zscore( self%atominfo(:)%max_int ) ! -"-
        call calc_stats(  self%atominfo(:)%avg_int,       self%avg_int_stats       )
        call calc_stats(  self%atominfo(:)%max_int,       self%max_int_stats       )
        call calc_stats(  self%atominfo(:)%valid_corr,    self%valid_corr_stats    )
        call calc_stats(  self%atominfo(:)%displ,         self%displ_stats         )
        call calc_stats(  self%atominfo(:)%max_ndispl,    self%max_ndispl_stats    )
        call calc_stats(  self%atominfo(:)%semiaxes(1),   self%semiaxis_maj_stats, mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%semiaxes(2),   self%semiaxis_med_stats, mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%semiaxes(3),   self%semiaxis_min_stats, mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%aniso_xyz(1),  self%aniso_x_stats, mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%aniso_xyz(2),  self%aniso_y_stats, mask=.not.self%atominfo(:)%tossADP )
        call calc_stats(  self%atominfo(:)%aniso_xyz(3),  self%aniso_z_stats, mask=.not.self%atominfo(:)%tossADP )
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
        call self%write_centers('doa_in_bfac_field',        'doa')
        call self%write_centers_aniso('aniso_bfac_field')
        ! destruct
        deallocate(mask, cc_mask, imat_cc, imat_cc_scaled, tmpcens, strain_array, centers_A, border, &
            &lattice_displ)
        call img_cc_scaled%kill_bimg
        call simatms%kill
        write(logfhandle, '(A)') '>>> EXTRACTING ATOM STATISTICS, COMPLETED'

        contains

            subroutine calc_zscore( arr )
                real, intent(inout) :: arr(:)
                arr = (arr - self%map_stats%avg) / self%map_stats%sdev
            end subroutine calc_zscore

            subroutine calc_cn_stats( cn )
                integer, intent(in)  :: cn ! calculate stats for given std cn
                integer :: cc, n, n_size, n_diam
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
                call calc_stats( self%atominfo(:)%displ,         self%displ_stats_cns(cn),         mask=cn_mask   )
                call calc_stats( self%atominfo(:)%max_ndispl,    self%max_ndispl_stats_cns(cn),    mask=cn_mask   )
                call calc_stats( self%atominfo(:)%semiaxes(1),   self%semiaxis_maj_stats_cns(cn),  mask=adp_mask  )
                call calc_stats( self%atominfo(:)%semiaxes(2),   self%semiaxis_med_stats_cns(cn),  mask=adp_mask  )
                call calc_stats( self%atominfo(:)%semiaxes(3),   self%semiaxis_min_stats_cns(cn),  mask=adp_mask  )
                call calc_stats( self%atominfo(:)%aniso_xyz(1),  self%aniso_x_stats_cns(cn),       mask=adp_mask  )
                call calc_stats( self%atominfo(:)%aniso_xyz(2),  self%aniso_y_stats_cns(cn),       mask=adp_mask  )
                call calc_stats( self%atominfo(:)%aniso_xyz(3),  self%aniso_z_stats_cns(cn),       mask=adp_mask  )
                call calc_stats( self%atominfo(:)%radial_strain, self%radial_strain_stats_cns(cn), mask=adp_mask  )
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

    ! For a given connected component cc, calculates the Anisotropic Displacement Parameters
    ! from the scaled CC map (imat_scaled) and the mask of border voxels in the scaled CC map (border).  
    ! The border voxels of cc are fit with an ellipsoid using an algebraic least squares fit
    ! (i.e. finding the best fit parameters B1,...,B6 for the function B1x^2 + B2y^2 + B3z^2 + B4xy + B5xz + B6yz = 1).  
    ! Reports the best fit ellipsoid semi-axes and the Anisotropic Displacement in the x,y, and z directions in Angstroms.
    ! Also reports a 3x3 anisotropic displacement matrix for writing an ANISOU PDB file to visualize ellipsoids in Chimera.
    ! Author: Henry Wietfeldt
    subroutine calc_aniso( self, cc, imat_scaled, border)
        class(nanoparticle), intent(inout)  :: self
        integer, intent(in)                 :: cc
        integer, allocatable, intent(in)    :: imat_scaled(:,:,:)
        logical, allocatable, intent(in)    :: border(:,:,:)
        real(kind=8), allocatable           :: A(:,:), AT(:,:), ones(:)
        real(kind=8)        :: beta(6), ATA(6,6), ATA_inv(6,6), matavg, B(3,3), eigenvals(3), eigenvecs(3,3)
        real(kind=8)        :: eigenvecs_inv(3,3), tmpeigenvecs(3,3), aniso(3,3), aniso_sq(3,3)
        real                :: center_scaled(3), maxrad, com(3), u, v, w
        integer             :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, ifoo, nborder, errflg

        ! To avoid overfitting, must have more voxels in the unscaled CC than fitting params
        if (self%atominfo(cc)%size <= NPARAMS_ADP) then
            self%atominfo(cc)%tossADP = .true.
            return
        end if

        ! Create search window that definitely contains the cc to speed up the iterations
        ! by avoiding having to iterate over the entire scaled images for each connected component.
        center_scaled = self%atominfo(cc)%center(:)*SCALE_FAC - 0.5*(SCALE_FAC-1)
        maxrad  = (self%theoretical_radius * 3) / (self%smpd / SCALE_FAC)
        ilo = max(nint(center_scaled(1) - maxrad), 1)
        ihi = min(nint(center_scaled(1) + maxrad), SCALE_FAC*self%ldim(1))
        jlo = max(nint(center_scaled(2) - maxrad), 1)
        jhi = min(nint(center_scaled(2) + maxrad), SCALE_FAC*self%ldim(2))
        klo = max(nint(center_scaled(3) - maxrad), 1)
        khi = min(nint(center_scaled(3) + maxrad), SCALE_FAC*self%ldim(3))

        ! Calculate the unweighted center of mass of the scaled connected component.
        com = 0.
        do k=klo, khi
            do j=jlo, jhi
                do i=ilo, ihi
                    if (imat_scaled(i,j,k) == cc) then
                        com(1) = com(1) + i
                        com(2) = com(2) + j
                        com(3) = com(3) + k
                    end if
                end do
            end do
        end do
        com = com / (self%atominfo(cc)%size * (SCALE_FAC**3))

        ! Find the number of border voxels in the scaled cc
        nborder = 0
        do k=klo, khi
            do j=jlo, jhi
                do i=ilo, ihi
                    if (imat_scaled(i,j,k) == cc .and. border(i,j,k)) then
                        nborder = nborder + 1
                    end if
                end do
            end do
        end do

        ! Fill in matrices for the least squares fit, with the center of mass as the origin.
        allocate(A(nborder, 6), AT(6, nborder), source=0.0_dp)
        allocate(ones(nborder), source=1.0_dp)
        nborder = 1
        do k=klo, khi
            do j=jlo, jhi
                do i=ilo, ihi
                    if (imat_scaled(i,j,k) == cc .and. border(i,j,k)) then
                        u = (i-com(1)) * self%smpd / SCALE_FAC
                        v = (j-com(2)) * self%smpd / SCALE_FAC
                        w = (k-com(3)) * self%smpd / SCALE_FAC
                        A(nborder, 1) = u**2
                        A(nborder, 2) = v**2
                        A(nborder, 3) = w**2
                        A(nborder, 4) = 2*u*v
                        A(nborder, 5) = 2*u*w
                        A(nborder, 6) = 2*v*w
                        nborder = nborder + 1
                    end if
                end do
            end do
        end do
        ! Find the least squares solution to A*beta=ones where beta are the fitting params.
        ! This minimizes the squared error in B1x^2 + B2y^2 + B3z^2 + B4xy + B5xz + B6yz = 1
        AT = transpose(A)
        ATA = matmul(AT, A) 
        matavg = sum(ATA)/size(ATA) !Normalize ATA to make matinv easier on the computer
        ATA = ATA / matavg
        call matinv(ATA, ATA_inv, 6, errflg)
        beta = matmul(ATA_inv, matmul(AT, ones)/matavg)
        
        ! Find the principal axes of the ellipsoid
        B(1,1) = beta(1)
        B(2,2) = beta(2)
        B(3,3) = beta(3)
        B(1,2) = beta(4)
        B(1,3) = beta(5)
        B(2,3) = beta(6)
        B(2,1) = B(1,2)
        B(3,1) = B(1,3)
        B(3,2) = B(2,3)
        call jacobi(B, 3, 3, eigenvals, eigenvecs, ifoo)
        eigenvals = 1/sqrt(eigenvals) ! Convert eigenvalues to semi-axes
        call eigsrt(eigenvals, eigenvecs, 3, 3)

        ! Fill in the anisotropic displacement matrix in the principal basis.
        aniso = 0.
        do i=1,3
            aniso(i,i) = eigenvals(i)
            self%atominfo(cc)%semiaxes(i) = aniso(i,i) / self%theoretical_radius 
        end do
        call matinv(eigenvecs, eigenvecs_inv, 3, errflg)
        if (errflg /= 0) then
            ! The above call to matinv is prone to numerical errors if some matrix elements are very small. 
            ! Find the approximate inverse by setting all very small numbers to 0.
            tmpeigenvecs = eigenvecs
            do j=1,3
                do i=1,3
                    if (abs(tmpeigenvecs(i,j)) < 0.00001) tmpeigenvecs(i,j) = 0._dp
                end do
            end do
            call matinv(tmpeigenvecs, eigenvecs_inv, 3, errflg)
            if (errflg /= 0) then
                write(logfhandle, '(a,3f10.5,a,3f10.5,a,3f10.5)') "Numerical error taking inverse of eigenvectors: ", &
                    &eigenvecs(:,1), ", ", eigenvecs(:,2), ", ", eigenvecs(:,3)
                write(logfhandle, '(a, i5)') "Tossing ADP Calculations of cc: ", cc
                self%atominfo(cc)%tossADP = .true.
            end if
        end if
        aniso_sq = aniso**2   ! ANISOU format uses the squared matrix values
        aniso = matmul(matmul(eigenvecs, aniso), eigenvecs_inv) ! Principal basis -> x,y,z basis
        do i=1,3
            self%atominfo(cc)%aniso_xyz(i) = aniso(i,i) / self%theoretical_radius
        end do
        self%atominfo(cc)%aniso = matmul(matmul(eigenvecs, aniso_sq), eigenvecs_inv) 
    end subroutine calc_aniso

    ! For the input cc, model of atomic centers, and lattice displacements |center - expected center|,
    ! finds the maximum lattice displacement of the atoms neighboring the atom cc.
    subroutine lattice_displ_analysis(self, cc, model, a, lattice_displ )
        class(nanoparticle), intent(inout) :: self
        integer,          intent(in)    :: cc
        real,             intent(in)    :: model(:,:), a(3), lattice_displ(:, :) ! a = lattice params
        character(len=2) :: el_ucase
        character(len=8) :: crystal_system
        real             :: a0, foo, d, max, max_cc, displ
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
        ! Calculate displacement magnitude of atom cc if not yet calculated
        if (self%atominfo(cc)%displ == 0.) then
            ! displ magnitude is not yet calculated
            self%atominfo(cc)%displ = norm_2(lattice_displ(cc, :))
        end if
        ! Find all atoms withing sphere of radius d and calculate max lattice displacement
        max = 0
        max_cc = -1
        do i=1, self%n_cc
            if (i/=cc .and. euclid(model(:, i), model(:, cc)) < d) then
                if (self%atominfo(i)%displ == 0.) then
                    ! displ magnitude is not yet calculated
                    self%atominfo(i)%displ = norm_2(lattice_displ(i, :))
                end if
                displ = self%atominfo(i)%displ
                if (displ > max) then
                    max = displ
                    max_cc = i
                end if
            end if
        end do
        self%atominfo(cc)%max_ndispl = max
        !print *, cc, self%atominfo(cc)%displ, max, max_cc
    end subroutine lattice_displ_analysis

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
        write(funit,'(a)') ATOM_STATS_HEAD_OMIT
        ! write records
        do cc = 1, size(self%atominfo)
            call self%write_atominfo(cc, funit, omit=.true.)
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
        write(funit,601,advance='no') real(self%atominfo(cc)%cc_ind),           CSV_DELIM ! INDEX
        write(funit,601,advance='no') real(self%atominfo(cc)%size),             CSV_DELIM ! NVOX
        write(funit,601,advance='no') real(self%atominfo(cc)%cn_std),           CSV_DELIM ! CN_STD
        write(funit,601,advance='no') self%atominfo(cc)%bondl,                  CSV_DELIM ! NN_BONDL
        write(funit,601,advance='no') self%atominfo(cc)%cn_gen,                 CSV_DELIM ! CN_GEN
        write(funit,601,advance='no') self%atominfo(cc)%diam,                   CSV_DELIM ! DIAM
        write(funit,601,advance='no') self%atominfo(cc)%avg_int,                CSV_DELIM ! AVG_INT
        write(funit,601,advance='no') self%atominfo(cc)%max_int,                CSV_DELIM ! MAX_INT
        write(funit,601,advance='no') self%atominfo(cc)%cendist,                CSV_DELIM ! CENDIST
        write(funit,601,advance='no') self%atominfo(cc)%valid_corr,             CSV_DELIM ! VALID_CORR
        write(funit,601,advance='no') self%atominfo(cc)%displ,                  CSV_DELIM ! DISPL
        write(funit,601,advance='no') self%atominfo(cc)%max_ndispl,             CSV_DELIM ! MAX_NDISPL
        ! aniso
        write(funit,601,advance='no') self%atominfo(cc)%semiaxes(1),            CSV_DELIM ! SEMIAX_MAJ
        write(funit,601,advance='no') self%atominfo(cc)%semiaxes(2),            CSV_DELIM ! SEMIAX_MED
        write(funit,601,advance='no') self%atominfo(cc)%semiaxes(3),            CSV_DELIM ! SEMIAX_MIN
        write(funit,601,advance='no') self%atominfo(cc)%aniso_xyz(1),           CSV_DELIM ! ANISO_X
        write(funit,601,advance='no') self%atominfo(cc)%aniso_xyz(2),           CSV_DELIM ! ANISO_Y
        write(funit,601,advance='no') self%atominfo(cc)%aniso_xyz(3),           CSV_DELIM ! ANISO_Z
        if( .not. omit_here )then
        write(funit,601,advance='no') self%atominfo(cc)%center(1),              CSV_DELIM ! X
        write(funit,601,advance='no') self%atominfo(cc)%center(2),              CSV_DELIM ! Y
        write(funit,601,advance='no') self%atominfo(cc)%center(3),              CSV_DELIM ! Z
        ! strain
        write(funit,601,advance='no') self%atominfo(cc)%exx_strain,             CSV_DELIM ! EXX_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%eyy_strain,             CSV_DELIM ! EYY_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%ezz_strain,             CSV_DELIM ! EZZ_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%exy_strain,             CSV_DELIM ! EXY_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%eyz_strain,             CSV_DELIM ! EYZ_STRAIN
        write(funit,601,advance='no') self%atominfo(cc)%exz_strain,             CSV_DELIM ! EXZ_STRAIN
        endif
        if( .not. omit_here )then
        write(funit,601,advance='no') self%atominfo(cc)%radial_strain,          CSV_DELIM ! RADIAL_STRAIN
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
        write(funit,601,advance='no') real(self%n_cc),               CSV_DELIM ! NATOMS
        write(funit,601,advance='no') real(self%n_aniso),            CSV_DELIM ! NANISO
        ! -- NP diameter
        write(funit,601,advance='no') self%NPdiam,                   CSV_DELIM ! DIAM
        ! -- atom size
        write(funit,601,advance='no') self%size_stats%avg,           CSV_DELIM ! AVG_NVOX
        write(funit,601,advance='no') self%size_stats%med,           CSV_DELIM ! MED_NVOX
        write(funit,601,advance='no') self%size_stats%sdev,          CSV_DELIM ! SDEV_NVOX
        ! -- standard coordination number
        write(funit,601,advance='no') self%cn_std_stats%avg,         CSV_DELIM ! AVG_CN_STD
        write(funit,601,advance='no') self%cn_std_stats%med,         CSV_DELIM ! MED_CN_STD
        write(funit,601,advance='no') self%cn_std_stats%sdev,        CSV_DELIM ! SDEV_CN_STD
        ! -- bond length
        write(funit,601,advance='no') self%bondl_stats%avg,          CSV_DELIM ! AVG_NN_BONDL
        write(funit,601,advance='no') self%bondl_stats%med,          CSV_DELIM ! MED_NN_BONDL
        write(funit,601,advance='no') self%bondl_stats%sdev,         CSV_DELIM ! SDEV_NN_BONDL
        ! -- generalized coordination number
        write(funit,601,advance='no') self%cn_gen_stats%avg,         CSV_DELIM ! AVG_CN_GEN
        write(funit,601,advance='no') self%cn_gen_stats%med,         CSV_DELIM ! MED_CN_GEN
        write(funit,601,advance='no') self%cn_gen_stats%sdev,        CSV_DELIM ! SDEV_CN_GEN
        ! -- atom diameter
        write(funit,601,advance='no') self%diam_stats%avg,           CSV_DELIM ! AVG_DIAM
        write(funit,601,advance='no') self%diam_stats%med,           CSV_DELIM ! MED_DIAM
        write(funit,601,advance='no') self%diam_stats%sdev,          CSV_DELIM ! SDEV_DIAM
        ! -- average intensity
        write(funit,601,advance='no') self%avg_int_stats%avg,        CSV_DELIM ! AVG_AVG_INT
        write(funit,601,advance='no') self%avg_int_stats%med,        CSV_DELIM ! MED_AVG_INT
        write(funit,601,advance='no') self%avg_int_stats%sdev,       CSV_DELIM ! SDEV_AVG_INT
        ! -- maximum intensity
        write(funit,601,advance='no') self%max_int_stats%avg,        CSV_DELIM ! AVG_MAX_INT
        write(funit,601,advance='no') self%max_int_stats%med,        CSV_DELIM ! MED_MAX_INT
        write(funit,601,advance='no') self%max_int_stats%sdev,       CSV_DELIM ! SDEV_MAX_INT
        ! -- maximum correlation
        write(funit,601,advance='no') self%valid_corr_stats%avg,     CSV_DELIM ! AVG_VALID_CORR
        write(funit,601,advance='no') self%valid_corr_stats%med,     CSV_DELIM ! MED_VALID_CORR
        write(funit,601,advance='no') self%valid_corr_stats%sdev,    CSV_DELIM ! SDEV_VALID_CORR
        ! -- lattice displacement
        write(funit,601,advance='no') self%displ_stats%avg,          CSV_DELIM ! AVG_DOA
        write(funit,601,advance='no') self%displ_stats%med,          CSV_DELIM ! MED_DOA
        write(funit,601,advance='no') self%displ_stats%sdev,         CSV_DELIM ! SDEV_DOA
        ! -- maximum neighboring displacement
        write(funit,601,advance='no') self%max_ndispl_stats%avg,     CSV_DELIM ! AVG_DOA
        write(funit,601,advance='no') self%max_ndispl_stats%med,     CSV_DELIM ! MED_DOA
        write(funit,601,advance='no') self%max_ndispl_stats%sdev,    CSV_DELIM ! SDEV_DOA
        ! -- Major semi-axis lengths of elliptical fits of CC shell
        write(funit,601,advance='no') self%semiaxis_maj_stats%avg,   CSV_DELIM ! AVG_SEMIAX_MAJ
        write(funit,601,advance='no') self%semiaxis_maj_stats%med,   CSV_DELIM ! MED_SEMIAX_MAJ
        write(funit,601,advance='no') self%semiaxis_maj_stats%sdev,  CSV_DELIM ! SDEV_SEMIAX_MAJ
        ! -- Medium semi-axis lengths of elliptical fits of CC shell
        write(funit,601,advance='no') self%semiaxis_med_stats%avg,   CSV_DELIM ! AVG_SEMIAX_MED
        write(funit,601,advance='no') self%semiaxis_med_stats%med,   CSV_DELIM ! MED_SEMIAX_MED
        write(funit,601,advance='no') self%semiaxis_med_stats%sdev,  CSV_DELIM ! SDEV_SEMIAX_MED
        ! -- Minor semi-axis lengths of elliptical fits of CC shell
        write(funit,601,advance='no') self%semiaxis_min_stats%avg,   CSV_DELIM ! AVG_SEMIAX_MIN
        write(funit,601,advance='no') self%semiaxis_min_stats%med,   CSV_DELIM ! MED_SEMIAX_MIN
        write(funit,601,advance='no') self%semiaxis_min_stats%sdev,  CSV_DELIM ! SDEV_SEMIAX_MIN
        ! -- Anisotropic displacement along x-axis
        write(funit,601,advance='no') self%aniso_x_stats%avg,        CSV_DELIM ! AVG_ANISO_X
        write(funit,601,advance='no') self%aniso_x_stats%med,        CSV_DELIM ! MED_ANISO_X
        write(funit,601,advance='no') self%aniso_x_stats%sdev,      CSV_DELIM ! SDEV_ANISO_X
        ! -- Anisotropic displacement along y-axis
        write(funit,601,advance='no') self%aniso_y_stats%avg,        CSV_DELIM ! AVG_ANISO_Y
        write(funit,601,advance='no') self%aniso_y_stats%med,        CSV_DELIM ! MED_ANISO_Y
        write(funit,601,advance='no') self%aniso_y_stats%sdev,      CSV_DELIM ! SDEV_ANISO_Y
        ! -- Anisotropic displacement along z-axis
        write(funit,601,advance='no') self%aniso_z_stats%avg,        CSV_DELIM ! AVG_ANISO_Z
        write(funit,601,advance='no') self%aniso_z_stats%med,        CSV_DELIM ! MED_ANISO_Z
        write(funit,601,advance='no') self%aniso_z_stats%sdev,      CSV_DELIM ! SDEV_ANISO_Z
        ! -- radial strain
        write(funit,601,advance='no') self%radial_strain_stats%avg,  CSV_DELIM ! AVG_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats%med,  CSV_DELIM ! MED_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats%sdev, CSV_DELIM ! SDEV_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats%minv, CSV_DELIM ! MIN_RADIAL_STRAIN
        write(funit,602)              self%radial_strain_stats%maxv            ! MAX_RADIAL_STRAIN
    end subroutine write_np_stats

    subroutine write_cn_stats( self, cn, funit )
        class(nanoparticle), intent(in) :: self
        integer,             intent(in) :: cn, funit
        601 format(F8.4,A2)
        602 format(F8.4)
        if( count(self%atominfo(:)%cn_std == cn) < 2 ) return
        ! -- coordination number
        write(funit,601,advance='no') real(cn),                              CSV_DELIM ! CN_STD
        ! -- # atoms per cn
        write(funit,601,advance='no') self%natoms_cns(cn),                   CSV_DELIM ! NATOMS
        write(funit,601,advance='no') self%natoms_aniso_cns(cn),             CSV_DELIM ! NANISO
        ! -- atom size
        write(funit,601,advance='no') self%size_stats_cns(cn)%avg,           CSV_DELIM ! AVG_NVOX
        write(funit,601,advance='no') self%size_stats_cns(cn)%med,           CSV_DELIM ! MED_NVOX
        write(funit,601,advance='no') self%size_stats_cns(cn)%sdev,          CSV_DELIM ! SDEV_NVOX
        ! -- bond length
        write(funit,601,advance='no') self%bondl_stats_cns(cn)%avg,          CSV_DELIM ! AVG_NN_BONDL
        write(funit,601,advance='no') self%bondl_stats_cns(cn)%med,          CSV_DELIM ! MED_NN_BONDL
        write(funit,601,advance='no') self%bondl_stats_cns(cn)%sdev,         CSV_DELIM ! SDEV_NN_BONDL
        ! -- generalized coordination number
        write(funit,601,advance='no') self%cn_gen_stats_cns(cn)%avg,         CSV_DELIM ! AVG_CN_GEN
        write(funit,601,advance='no') self%cn_gen_stats_cns(cn)%med,         CSV_DELIM ! MED_CN_GEN
        write(funit,601,advance='no') self%cn_gen_stats_cns(cn)%sdev,        CSV_DELIM ! SDEV_CN_GEN
        ! -- atom diameter
        write(funit,601,advance='no') self%diam_stats_cns(cn)%avg,           CSV_DELIM ! AVG_DIAM
        write(funit,601,advance='no') self%diam_stats_cns(cn)%med,           CSV_DELIM ! MED_DIAM
        write(funit,601,advance='no') self%diam_stats_cns(cn)%sdev,          CSV_DELIM ! SDEV_DIAM
        ! -- average intensity
        write(funit,601,advance='no') self%avg_int_stats_cns(cn)%avg,        CSV_DELIM ! AVG_AVG_INT
        write(funit,601,advance='no') self%avg_int_stats_cns(cn)%med,        CSV_DELIM ! MED_AVG_INT
        write(funit,601,advance='no') self%avg_int_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_AVG_INT
        ! -- maximum intensity
        write(funit,601,advance='no') self%max_int_stats_cns(cn)%avg,        CSV_DELIM ! AVG_MAX_INT
        write(funit,601,advance='no') self%max_int_stats_cns(cn)%med,        CSV_DELIM ! MED_MAX_INT
        write(funit,601,advance='no') self%max_int_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_MAX_INT
        ! -- maximum correlation
        write(funit,601,advance='no') self%valid_corr_stats_cns(cn)%avg,     CSV_DELIM ! AVG_VALID_CORR
        write(funit,601,advance='no') self%valid_corr_stats_cns(cn)%med,     CSV_DELIM ! MED_VALID_CORR
        write(funit,601,advance='no') self%valid_corr_stats_cns(cn)%sdev,    CSV_DELIM ! SDEV_VALID_CORR
        ! -- lattice displacement
        write(funit,601,advance='no') self%displ_stats_cns(cn)%avg,          CSV_DELIM ! AVG_MAX_NDISPL
        write(funit,601,advance='no') self%displ_stats_cns(cn)%med,          CSV_DELIM ! MED_MAX_NDISPL
        write(funit,601,advance='no') self%displ_stats_cns(cn)%sdev,         CSV_DELIM! SDEV_MAX_NDISPL
        ! -- maximum displacement of neighboring atoms
        write(funit,601,advance='no') self%max_ndispl_stats_cns(cn)%avg,     CSV_DELIM ! AVG_MAX_NDISPL
        write(funit,601,advance='no') self%max_ndispl_stats_cns(cn)%med,     CSV_DELIM ! MED_MAX_NDISPL
        write(funit,601,advance='no') self%max_ndispl_stats_cns(cn)%sdev,    CSV_DELIM! SDEV_MAX_NDISPL
        ! -- Major semi-axis lengths of elliptical fits of CC shell
        write(funit,601,advance='no') self%semiaxis_maj_stats_cns(cn)%avg,   CSV_DELIM ! AVG_SEMIAX_MAJ
        write(funit,601,advance='no') self%semiaxis_maj_stats_cns(cn)%med,   CSV_DELIM ! MED_SEMIAX_MAJ
        write(funit,601,advance='no') self%semiaxis_maj_stats_cns(cn)%sdev,  CSV_DELIM ! SDEV_SEMIAX_MAJ
        ! -- Medium semi-axis lengths of elliptical fits of CC shell
        write(funit,601,advance='no') self%semiaxis_med_stats_cns(cn)%avg,   CSV_DELIM ! AVG_SEMIAX_MED
        write(funit,601,advance='no') self%semiaxis_med_stats_cns(cn)%med,   CSV_DELIM ! MED_SEMIAX_MED
        write(funit,601,advance='no') self%semiaxis_med_stats_cns(cn)%sdev,  CSV_DELIM ! SDEV_SEMIAX_MED
        ! -- Minor semi-axis lengths of elliptical fits of CC shell
        write(funit,601,advance='no') self%semiaxis_min_stats_cns(cn)%avg,   CSV_DELIM ! AVG_SEMIAX_MIN
        write(funit,601,advance='no') self%semiaxis_min_stats_cns(cn)%med,   CSV_DELIM ! MED_SEMIAX_MIN
        write(funit,601,advance='no') self%semiaxis_min_stats_cns(cn)%sdev,  CSV_DELIM ! SDEV_SEMIAX_MIN
        ! -- Anisotropic displacement along x-axis
        write(funit,601,advance='no') self%aniso_x_stats_cns(cn)%avg,        CSV_DELIM ! AVG_ANISO_X
        write(funit,601,advance='no') self%aniso_x_stats_cns(cn)%med,        CSV_DELIM ! MED_ANISO_X
        write(funit,601,advance='no') self%aniso_x_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_ANISO_X
        ! -- Anisotropic displacement along y-axis
        write(funit,601,advance='no') self%aniso_y_stats_cns(cn)%avg,        CSV_DELIM ! AVG_ANISO_Y
        write(funit,601,advance='no') self%aniso_y_stats_cns(cn)%med,        CSV_DELIM ! MED_ANISO_Y
        write(funit,601,advance='no') self%aniso_y_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_ANISO_Y
        ! -- Anisotropic displacement along z-axis
        write(funit,601,advance='no') self%aniso_z_stats_cns(cn)%avg,        CSV_DELIM ! AVG_ANISO_Z
        write(funit,601,advance='no') self%aniso_z_stats_cns(cn)%med,        CSV_DELIM ! MED_ANISO_Z
        write(funit,601,advance='no') self%aniso_z_stats_cns(cn)%sdev,       CSV_DELIM ! SDEV_ANISO_Z
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

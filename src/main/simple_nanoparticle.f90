module simple_nanoparticle
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,    only: image
use simple_binimage, only: binimage
use simple_atoms,    only: atoms, get_Z_and_radius_from_name
use simple_nanoparticle_utils
implicit none

public :: nanoparticle
private
#include "simple_local_flags.inc"

! module global constants
integer,          parameter :: NBIN_THRESH         = 20      ! number of thresholds for binarization
integer,          parameter :: NVOX_THRESH         = 3       ! min # voxels per atom is 3
logical,          parameter :: DEBUG               = .false. ! for debugging purposes
logical,          parameter :: GENERATE_FIGS       = .false. ! for figures generation
integer,          parameter :: SOFT_EDGE           = 6
integer,          parameter :: N_DISCRET           = 1000
integer,          parameter :: CNMIN               = 5
integer,          parameter :: CNMAX               = 12
integer,          parameter :: NSTRAIN_COMPS       = 7
character(len=3), parameter :: CSV_DELIM           = ', '
character(len=*), parameter :: ATOM_STATS_FILE     = 'atom_stats.csv'
character(len=*), parameter :: NP_STATS_FILE       = 'nanoparticle_stats.csv'
character(len=*), parameter :: CN_STATS_FILE       = 'cn_dependent_stats.csv'
character(len=*), parameter :: ATOM_VAR_CORRS_FILE = 'atom_param_corrs.txt'
character(len=*), parameter :: BIN_CLS_STATS_FILE  = 'binary_class_stats.csv'

character(len=*), parameter :: ATOM_STATS_HEAD = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//'CN_STD'//CSV_DELIM//'NN_BONDL'//&
&CSV_DELIM//'CN_GEN'//CSV_DELIM//'ASPECT_RATIO'//CSV_DELIM//'POLAR_ANGLE'//CSV_DELIM//'DIAM'//CSV_DELIM//'AVG_INT'//&
&CSV_DELIM//'MAX_INT'//CSV_DELIM//'SDEV_INT'//CSV_DELIM//'RADIAL_POS'//CSV_DELIM//'MAX_CORR'//CSV_DELIM//'X'//&
&CSV_DELIM//'Y'//CSV_DELIM//'Z'//CSV_DELIM//'EXX_STRAIN'//CSV_DELIM//'EYY_STRAIN'//&
&CSV_DELIM//'EZZ_STRAIN'//CSV_DELIM//'EXY_STRAIN'//CSV_DELIM//'EYZ_STRAIN'//CSV_DELIM//'EXZ_STRAIN'//CSV_DELIM//'RADIAL_STRAIN'//&
&CSV_DELIM//'NVOX_CLASS'//CSV_DELIM//'NN_BONDL_CLASS'//CSV_DELIM//'ASPECT_RATIO_CLASS'//CSV_DELIM//'POLAR_ANGLE_CLASS'//&
&CSV_DELIM//'DIAM_CLASS'//CSV_DELIM//'MAX_INT_CLASS'//CSV_DELIM//'MAX_CORR_CLASS'//CSV_DELIM//'RADIAL_STRAIN_CLASS'

character(len=*), parameter :: ATOM_STATS_HEAD_OMIT = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//'CN_STD'//CSV_DELIM//'NN_BONDL'//&
&CSV_DELIM//'CN_GEN'//CSV_DELIM//'ASPECT_RATIO'//CSV_DELIM//'POLAR_ANGLE'//CSV_DELIM//'DIAM'//CSV_DELIM//'AVG_INT'//&
&CSV_DELIM//'MAX_INT'//CSV_DELIM//'SDEV_INT'//CSV_DELIM//'RADIAL_POS'//CSV_DELIM//'MAX_CORR'//CSV_DELIM//'RADIAL_STRAIN'

character(len=*), parameter :: NP_STATS_HEAD = 'NATOMS'//CSV_DELIM//'DIAM'//&
&CSV_DELIM//'AVG_NVOX'//CSV_DELIM//'MED_NVOX'//CSV_DELIM//'SDEV_NVOX'//&
&CSV_DELIM//'AVG_CN_STD'//CSV_DELIM//'MED_CN_STD'//CSV_DELIM//'SDEV_CN_STD'//&
&CSV_DELIM//'AVG_NN_BONDL'//CSV_DELIM//'MED_NN_BONDL'//CSV_DELIM//'SDEV_NN_BONDL'//&
&CSV_DELIM//'AVG_CN_GEN'//CSV_DELIM//'MED_CN_GEN'//CSV_DELIM//'SDEV_CN_GEN'//&
&CSV_DELIM//'AVG_ASPECT_RATIO'//CSV_DELIM//'MED_ASPECT_RATIO'//CSV_DELIM//'SDEV_ASPECT_RATIO'//&
&CSV_DELIM//'AVG_POLAR_ANGLE'//CSV_DELIM//'MED_POLAR_ANGLE'//CSV_DELIM//'SDEV_POLAR_ANGLE'//&
&CSV_DELIM//'AVG_DIAM'//CSV_DELIM//'MED_DIAM'//CSV_DELIM//'SDEV_DIAM'//&
&CSV_DELIM//'AVG_AVG_INT'//CSV_DELIM//'MED_AVG_INT'//CSV_DELIM//'SDEV_AVG_INT'//&
&CSV_DELIM//'AVG_MAX_INT'//CSV_DELIM//'MED_MAX_INT'//CSV_DELIM//'SDEV_MAX_INT'//&
&CSV_DELIM//'AVG_MAX_CORR'//CSV_DELIM//'MED_MAX_CORR'//CSV_DELIM//'SDEV_MAX_CORR'//&
&CSV_DELIM//'AVG_RADIAL_STRAIN'//CSV_DELIM//'MED_RADIAL_STRAIN'//CSV_DELIM//'SDEV_RADIAL_STRAIN'//&
&CSV_DELIM//'MIN_RADIAL_STRAIN'//CSV_DELIM//'MAX_RADIAL_STRAIN'

character(len=*), parameter :: CN_STATS_HEAD = 'CN_STD'//CSV_DELIM//'NATOMS'//CSV_DELIM//'NATOMS_PER_VOL'//&
&CSV_DELIM//'AVG_NVOX'//CSV_DELIM//'MED_NVOX'//CSV_DELIM//'SDEV_NVOX'//&
&CSV_DELIM//'AVG_NN_BONDL'//CSV_DELIM//'MED_NN_BONDL'//CSV_DELIM//'SDEV_NN_BONDL'//&
&CSV_DELIM//'AVG_CN_GEN'//CSV_DELIM//'MED_CN_GEN'//CSV_DELIM//'SDEV_CN_GEN'//&
&CSV_DELIM//'AVG_ASPECT_RATIO'//CSV_DELIM//'MED_ASPECT_RATIO'//CSV_DELIM//'SDEV_ASPECT_RATIO'//&
&CSV_DELIM//'AVG_POLAR_ANGLE'//CSV_DELIM//'MED_POLAR_ANGLE'//CSV_DELIM//'SDEV_POLAR_ANGLE'//&
&CSV_DELIM//'AVG_DIAM'//CSV_DELIM//'MED_DIAM'//CSV_DELIM//'SDEV_DIAM'//&
&CSV_DELIM//'AVG_AVG_INT'//CSV_DELIM//'MED_AVG_INT'//CSV_DELIM//'SDEV_AVG_INT'//&
&CSV_DELIM//'AVG_MAX_INT'//CSV_DELIM//'MED_MAX_INT'//CSV_DELIM//'SDEV_MAX_INT'//&
&CSV_DELIM//'AVG_MAX_CORR'//CSV_DELIM//'MED_MAX_CORR'//CSV_DELIM//'SDEV_MAX_CORR'//&
&CSV_DELIM//'AVG_RADIAL_STRAIN'//CSV_DELIM//'MED_RADIAL_STRAIN'//CSV_DELIM//'SDEV_RADIAL_STRAIN'//&
&CSV_DELIM//'MIN_RADIAL_STRAIN'//CSV_DELIM//'MAX_RADIAL_STRAIN'

character(len=*), parameter :: BIN_CLS_STATS_HEAD = 'PARAMETER'//&
&CSV_DELIM//'AVG_1'//CSV_DELIM//'MED_1'//CSV_DELIM//'SDEV_1'//CSV_DELIM//'MIN_1'//CSV_DELIM//'MAX_1'//CSV_DELIM//'POP_1'//&
&CSV_DELIM//'AVG_2'//CSV_DELIM//'MED_2'//CSV_DELIM//'SDEV_2'//CSV_DELIM//'MIN_2'//CSV_DELIM//'MAX_2'//CSV_DELIM//'POP_2'

! container for per-atom statistics
type :: atom_stats
    ! various per-atom parameters
    integer :: cc_ind            = 0  ! index of the connected component                            INDEX
    integer :: size              = 0  ! number of voxels in connected component                     NVOX
    integer :: cn_std            = 0  ! standard coordination number                                CN_STD
    integer :: loc_ldist(3)      = 0  ! vxl that determins the longest dim of the atom
    real    :: bondl             = 0. ! nearest neighbour bond lenght in A                          NN_BONDL
    real    :: cn_gen            = 0. ! generalized coordination number                             CN_GEN
    real    :: aspect_ratio      = 0. !                                                             ASPECT_RATIO
    real    :: polar_vec(3)      = 0. ! polarization vector
    real    :: polar_angle       = 0. ! polarization angle                                          POLAR_ANGLE
    real    :: diam              = 0. ! atom diameter                                               DIAM
    real    :: avg_int           = 0. ! average grey level intensity across the connected component AVG_INT
    real    :: max_int           = 0. ! maximum            -"-                                      MAX_INT
    real    :: sdev_int          = 0. ! standard deviation -"-                                      SDEV_INT
    real    :: cendist           = 0. ! distance from the centre of mass of the nanoparticle        RADIAL_POS
    real    :: max_corr          = 0. ! maximum atom correlation within the connected component     MAX_CORR
    real    :: center(3)         = 0. ! atom center                                                 X Y Z
    ! strain
    real    :: exx_strain        = 0. ! tensile strain in %                                         EXX_STRAIN
    real    :: eyy_strain        = 0. ! -"-                                                         EYY_STRAIN
    real    :: ezz_strain        = 0. ! -"-                                                         EZZ_STRAIN
    real    :: exy_strain        = 0. ! -"-                                                         EXY_STRAIN
    real    :: eyz_strain        = 0. ! -"-                                                         EYZ_STRAIN
    real    :: exz_strain        = 0. ! -"-                                                         EXZ_STRAIN
    real    :: radial_strain     = 0. ! -"-                                                         RADIAL_STRAIN
    ! binary cluster assignments
    integer :: size_cls          = 0  !                                                             NVOX_CLASS
    integer :: bondl_cls         = 0  !                                                             NN_BONDL_CLASS
    integer :: aspect_ratio_cls  = 0  !                                                             ASPECT_RATIO_CLASS
    integer :: polar_angle_cls   = 0  !                                                             POLAR_ANGLE_CLASS
    integer :: diam_cls          = 0  !                                                             DIAM_CLASS
    integer :: max_int_cls       = 0  !                                                             MAX_INT_CLASS
    integer :: max_corr_cls      = 0  !                                                             MAX_CORR_CLASS
    integer :: radial_strain_cls = 0  !                                                             RADIAL_STRAIN_CLASS
end type atom_stats

type :: nanoparticle
    private
    type(atoms)           :: centers_pdb
    type(image)           :: img, img_raw
    type(binimage)        :: img_bin, img_cc
    integer               :: ldim(3)            = 0  ! logical dimension of image
    integer               :: n_cc               = 0  ! number of atoms (connected components)                NATOMS
    real                  :: smpd               = 0. ! sampling distance
    real                  :: NPcen(3)           = 0. ! coordinates of the center of mass of the nanoparticle
    real                  :: NPdiam             = 0. ! diameter of the nanoparticle                          DIAM
    real                  :: theoretical_radius = 0. ! theoretical atom radius in A
    ! GLOBAL NP STATS
    ! -- dipole stats
    real                  :: net_dipole(3)      = 0. ! sum of all the directions of polarization             X_POLAR Y_POLAR Z_POLAR
    real                  :: net_dipole_mag     = 0. ! net dipole magnitude                                  POLAR_MAG
    real                  :: net_dipole_ang     = 0. ! net dipole angle                                      POLAR_ANGLE
    ! -- the rest
    type(stats_struct)    :: size_stats
    type(stats_struct)    :: cn_std_stats
    type(stats_struct)    :: bondl_stats
    type(stats_struct)    :: cn_gen_stats
    type(stats_struct)    :: aspect_ratio_stats
    type(stats_struct)    :: polar_angle_stats
    type(stats_struct)    :: diam_stats
    type(stats_struct)    :: avg_int_stats
    type(stats_struct)    :: max_int_stats
    type(stats_struct)    :: max_corr_stats
    type(stats_struct)    :: radial_strain_stats
    ! CN-DEPENDENT STATS
    ! -- cendist
    real                  :: max_cendist_cns(CNMIN:CNMAX)     = 0. ! max distance from the centre of mass of the NP   MAX_RADIAL_POS
    real                  :: min_cendist_cns(CNMIN:CNMAX)     = 0. ! min -"-                                          MIN_RADIAL_POS
    ! -- # atoms
    real                  :: natoms_cns(CNMIN:CNMAX)          = 0. ! # of atoms per cn_std                            NATOMS
    ! -- # atoms per volume
    real                  :: atoms_per_vol_cns(CNMIN:CNMAX)   = 0. ! # atoms per volume (#/A**3)                      NATOMS_PER_VOL
    ! -- the rest
    type(stats_struct)    :: size_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: bondl_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: cn_gen_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: aspect_ratio_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: polar_angle_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: diam_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: avg_int_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: max_int_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: max_corr_stats_cns(CNMIN:CNMAX)
    type(stats_struct)    :: radial_strain_stats_cns(CNMIN:CNMAX)
    ! PER-ATOM STATISTICS
    type(atom_stats), allocatable :: atominfo(:)
    ! BINARY CLASS STATISTICS
    type(stats_struct)    :: size_cls_stats(2)
    type(stats_struct)    :: bondl_cls_stats(2)
    type(stats_struct)    :: aspect_ratio_cls_stats(2)
    type(stats_struct)    :: polar_angle_cls_stats(2)
    type(stats_struct)    :: diam_cls_stats(2)
    type(stats_struct)    :: max_int_cls_stats(2)
    type(stats_struct)    :: max_corr_cls_stats(2)
    type(stats_struct)    :: radial_strain_cls_stats(2)
    ! OTHER
    character(len=2)      :: element     = ' '
    character(len=4)      :: atom_name   = '    '
    character(len=STDLEN) :: partname    = '' ! fname
    character(len=STDLEN) :: fbody       = '' ! fbody
  contains
    ! constructor
    procedure          :: new => new_nanoparticle
    ! getters/setters
    procedure          :: get_ldim
    procedure          :: get_natoms
    procedure          :: set_img
    procedure          :: set_atomic_coords
    ! utils
    procedure, private :: atominfo2centers
    procedure, private :: atominfo2centers_A
    procedure, private :: center_on_atom
    procedure          :: mask
    procedure          :: norm
    procedure          :: make_soft_mask
    procedure          :: update_ncc
    ! atomic position determination
    procedure          :: identify_atomic_pos
    procedure, private :: phasecorr_one_atom
    procedure, private :: binarize_and_find_centers
    procedure, private :: find_centers
    procedure, private :: discard_outliers
    procedure, private :: split_atoms
    ! calc stats
    procedure          :: fillin_atominfo
    procedure, private :: masscen
    procedure, private :: calc_aspect_ratio
    ! visualization and output
    procedure          :: write_centers
    procedure          :: write_atoms
    procedure          :: write_csv_files
    procedure, private :: write_atominfo
    procedure, private :: write_np_stats
    procedure, private :: write_cn_stats
    ! clustering
    procedure, private :: id_corr_vars
    procedure, private :: bicluster_otsu
    procedure, private :: cluster_atom_features
    procedure, private :: cluster_atom_intensity
    procedure          :: cluster_atom_maxint
    procedure          :: cluster_atom_intint
    procedure          :: cluster_ang
    procedure          :: cluster_ar
    procedure          :: cluster_bondl
    ! others
    procedure          :: geometry_analysis
    ! kill
    procedure          :: kill => kill_nanoparticle
end type nanoparticle

contains

    subroutine new_nanoparticle(self, fname, cline_smpd, element)
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        real,                intent(in)    :: cline_smpd
        character(len=2),    intent(inout) :: element
        integer :: nptcls
        integer :: Z ! atomic number
        real    :: smpd
        call self%kill
        self%partname  = fname
        self%fbody     = get_fbody(trim(basename(fname)), trim(fname2ext(fname)))
        self%smpd      = cline_smpd
        self%atom_name = ' '//element
        self%element   = element
        call get_Z_and_radius_from_name(self%atom_name, Z, self%theoretical_radius)
        call find_ldim_nptcls(self%partname, self%ldim, nptcls, smpd)
        call self%img%new(self%ldim, self%smpd)
        call self%img_raw%new(self%ldim, self%smpd)
        call self%img_bin%new_bimg(self%ldim, self%smpd)
        call self%img_bin%new(self%ldim, self%smpd)
        call self%img%read(fname)
        call self%img_raw%read(fname)
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

    ! utils

    function atominfo2centers( self ) result( centers )
        class(nanoparticle), intent(in) :: self
        real, allocatable :: centers(:,:)
        integer :: sz, i
        sz = size(self%atominfo)
        allocate(centers(3,sz), source=0.)
        do i = 1, sz
            centers(:,i) = self%atominfo(i)%center(:)
        end do
    end function atominfo2centers

    function atominfo2centers_A( self ) result( centers_A )
        class(nanoparticle), intent(in) :: self
        real, allocatable :: centers_A(:,:)
        integer :: sz, i
        sz = size(self%atominfo)
        allocate(centers_A(3,sz), source=0.)
        do i = 1, sz
            centers_A(:,i) = (self%atominfo(i)%center(:) - 1.) * self%smpd
        end do
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

    subroutine mask( self, msk_rad )
        class(nanoparticle), intent(inout) :: self
        real,                intent(in)    :: msk_rad
        call self%img%mask(msk_rad, 'soft')
    end subroutine mask

    subroutine norm( self )
        class(nanoparticle), intent(inout) :: self
        call self%img%norm
    end subroutine norm

    subroutine make_soft_mask(self)
        class(nanoparticle), intent(inout) :: self
        type(binimage) :: simulated_density
        type(image)    :: img_cos
        type(atoms)    :: atomic_pos
        call simulated_density%new_bimg(self%ldim, self%smpd)
        call atomic_pos%new(trim(self%fbody)//'_atom_centers.pdb')
        call atomic_pos%convolve(simulated_density, cutoff=8.*self%smpd)
        call simulated_density%grow_bins(nint(0.5*self%theoretical_radius/self%smpd)+1)
        call simulated_density%cos_edge(SOFT_EDGE, img_cos)
        call img_cos%write(trim(self%fbody)//'SoftMask.mrc')
        call img_cos%kill
        call simulated_density%kill_bimg
        call atomic_pos%kill
    end subroutine make_soft_mask

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

    subroutine identify_atomic_pos( self, cn_thresh )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: cn_thresh
        ! Phase correlation approach
        call self%phasecorr_one_atom(self%img)
        ! Nanoparticle binarization
        call self%binarize_and_find_centers()
        ! atom splitting by correlation map validation
        call self%split_atoms()
        ! Outliers discarding
        call self%discard_outliers(cn_thresh)
        ! write output
        call self%img_bin%write_bimg(trim(self%fbody)//'BIN.mrc')
        write(logfhandle,*) 'output, binarized map:            ', trim(self%fbody)//'BIN.mrc'
        call self%img_cc%write_bimg(trim(self%fbody)//'CC.mrc')
        write(logfhandle,*) 'output, connected components map: ', trim(self%fbody)//'CC.mrc'
        call self%write_centers()
    end subroutine identify_atomic_pos

    ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
    subroutine phasecorr_one_atom( self, out_img )
        use simple_segmentation
        class(nanoparticle), intent(inout) :: self
        type(image),         intent(inout) :: out_img  ! where to save the correlation values
        type(image) :: one_atom, img_copy, phasecorr
        type(atoms) :: atom
        real :: cutoff
        call img_copy%new(self%ldim, self%smpd)
        call img_copy%copy(self%img)
        call phasecorr%new(self%ldim, self%smpd)
        call phasecorr%set_ft(.true.)
        call one_atom%new(self%ldim,self%smpd)
        cutoff = 8.*self%smpd
        call atom%new(1)
        call atom%set_element(1,self%element)
        call atom%set_coord(1,self%smpd*(real(self%ldim)/2.)) ! DO NOT NEED THE +1
        call atom%convolve(one_atom, cutoff)
        call one_atom%fft()
        call img_copy%fft()
        call img_copy%phase_corr(one_atom,phasecorr,1.)
        if( GENERATE_FIGS ) call phasecorr%write(trim(self%fbody)//'CorrFiltered.mrc')
        ! Save Output
        call out_img%kill() ! if it exists already
        call out_img%new(self%ldim, self%smpd)
        call out_img%fft()
        call out_img%copy(phasecorr)
        call one_atom%kill()
        call phasecorr%kill()
    end subroutine phasecorr_one_atom

    ! This subrotuine takes in input a nanoparticle and
    ! binarizes it by thresholding. The gray level histogram is split
    ! in 20 parts, which corrispond to 20 possible threshold.
    ! Among those threshold, the selected one is the for which
    ! tha correlation between the raw map and a simulated distribution
    ! obtained with that threshold reaches the maximum value.
    subroutine binarize_and_find_centers( self )
        class(nanoparticle), intent(inout) :: self
        type(binimage)       :: img_bin_thresh(NBIN_THRESH/2-1)
        type(binimage)       :: img_ccs_thresh(NBIN_THRESH/2-1)
        type(image)          :: pc
        type(atoms)          :: atom
        type(image)          :: simulated_distrib
        integer, allocatable :: imat_t(:,:,:)
        real,    allocatable :: x_mat(:)  ! vectorization of the volume
        real,    allocatable :: coords(:,:)
        real,    allocatable :: rmat(:,:,:)
        integer :: i, t
        real    :: otsu_thresh
        real    :: x_thresh(NBIN_THRESH/2-1)
        real    :: step, maximum, mm(2)
        call otsu_nano(self%img,otsu_thresh) ! find initial threshold
        write(logfhandle,'(A)') '>>> BINARIZATION'
        rmat = self%img%get_rmat()
        x_mat = pack(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) >= 0.)
        allocate(imat_t(self%ldim(1), self%ldim(2), self%ldim(3)), source = 0)
        step = (maxval(x_mat)-otsu_thresh )/real(NBIN_THRESH)
        deallocate(x_mat)
        call simulated_distrib%new(self%ldim,self%smpd)
        call pc%new(self%ldim,self%smpd)
        call self%img%fft ! for pc calculation
        t = 1
        maximum = 0.
        do i = 1, NBIN_THRESH/2-1
            call progress(i, NBIN_THRESH/2-1)
            call img_bin_thresh(i)%new_bimg(self%ldim, self%smpd)
            if( i == 1 )then
                x_thresh(i) = otsu_thresh
            else
                x_thresh(i) = x_thresh(i-1) + step
            endif
            where(rmat > x_thresh(i))
                imat_t = 1
            elsewhere
                imat_t = 0
            endwhere
            ! Generate binary image and cc image
            call img_bin_thresh(i)%set_imat(imat_t)
            call img_bin_thresh(i)%find_ccs(img_ccs_thresh(i))
            ! Find atom centers in the generated distributions
            call self%update_ncc(img_ccs_thresh(i)) ! self%n_cc is needed in find_centers
            call self%find_centers(img_bin_thresh(i), img_ccs_thresh(i), coords)
            ! Generate a simulated distribution based on those center
            call self%write_centers('centers_'//trim(int2str(i))//'_iteration', coords)
            call atom%new          ('centers_'//trim(int2str(i))//'_iteration.pdb')
            call atom%convolve(simulated_distrib, cutoff = 8.*self%smpd)
            call del_file('centers_'//trim(int2str(i))//'_iteration.pdb')
            call atom%kill
            ! Take care of Fourier status, for phase_corr calculation
            call simulated_distrib%fft
            call pc%zero_and_flag_ft
            ! Correlation volume generation
            call self%img%phase_corr(simulated_distrib, pc, lp=3.)
            ! Calculation and update of the maxval the correlation reaches
            mm = pc%minmax()
            if( mm(2) > maximum )then
                maximum = mm(2)
                t       = i
            endif
            call simulated_distrib%set_ft(.false.)
        enddo
        call self%img%ifft ! To remove
        write(logfhandle,*) 'Selected threshold: ', x_thresh(t)
        ! Update img_bin and img_cc
        call self%img_bin%copy_bimg(img_bin_thresh(t))
        if( GENERATE_FIGS ) call self%img_bin%write_bimg(trim(self%fbody)//'SelectedThreshold.mrc')
        call self%img_cc%copy_bimg(img_ccs_thresh(t))
        call self%update_ncc()
        call self%find_centers()
        do i = 1,  NBIN_THRESH/2-1
            call img_bin_thresh(i)%kill_bimg
            call img_ccs_thresh(i)%kill_bimg
        enddo
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
             use simple_math, only : otsu
             type(image),    intent(inout) :: img
             real,           intent(out)   :: scaled_thresh ! returns the threshold in the correct range
             real, pointer     :: rmat(:,:,:)
             real, allocatable :: x(:)
             call img%get_rmat_ptr(rmat)
             x = pack(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) > 0.)
             call otsu(x, scaled_thresh)
         end subroutine otsu_nano

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

    subroutine split_atoms( self )
        class(nanoparticle), intent(inout) :: self
        real,    allocatable :: x(:)
        real,    pointer     :: rmat_pc(:,:,:)
        integer, allocatable :: imat(:,:,:), imat_cc(:,:,:), imat_bin(:,:,:)
        integer, parameter   :: RANK_THRESH = 4
        integer :: icc, cnt, cnt_split
        integer :: rank, m(1)
        real    :: new_centers(3,2*self%n_cc) ! will pack it afterwards if it has too many elements
        real    :: pc, aspect_ratio, radius
        write(logfhandle, '(A)') '>>> VALIDATING ATOMIC POSITIONS'
        call self%img%get_rmat_ptr(rmat_pc) ! now img contains the phase correlation
        call self%img_cc%get_imat(imat_cc)  ! to pass to the subroutine split_atoms
        allocate(imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), source = imat_cc)
        cnt       = 0
        cnt_split = 0
        ! Remember to update the centers
        do icc =1, self%n_cc ! for each cc check if the center corresponds with the local max of the phase corr
            pc = rmat_pc(nint(self%atominfo(icc)%center(1)),nint(self%atominfo(icc)%center(2)),nint(self%atominfo(icc)%center(3)))
            ! calculate the rank
            x = pack(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) == icc)
            call hpsort(x)
            m(:) = minloc(abs(x - pc))
            rank = size(x) - m(1)
            deallocate(x)
            ! calculate radius
            call self%calc_aspect_ratio(icc, aspect_ratio, radius)
            ! split
            if( rank > RANK_THRESH .or. radius > 1.5 * self%theoretical_radius )then
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
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%update_img_rmat()
        call self%img_bin%find_ccs(self%img_cc)
        call self%update_ncc(self%img_cc)
        call self%find_centers()
        write(logfhandle, '(A)') '>>> VALIDATING ATOMIC POSITIONS, COMPLETED'

    contains

        subroutine split_atom(new_centers,cnt)
            real,    intent(inout) :: new_centers(:,:) ! updated coordinates of the centers
            integer, intent(inout) :: cnt              ! atom counter, to update the center coords
            integer :: new_center1(3), new_center2(3), new_center3(3)
            integer :: i, j, k
            logical :: mask(self%ldim(1),self%ldim(2),self%ldim(3)) ! false in the layer of connection of the atom to be split
            mask = .false. ! initialization
            ! Identify first new center
            new_center1 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), imat == icc)
            cnt = cnt + 1
            new_centers(:,cnt) = real(new_center1)
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if(((real(i - new_center1(1)))**2 + (real(j - new_center1(2)))**2 + &
                        &   (real(k - new_center1(3)))**2) * self%smpd  <=  (0.9 * self%theoretical_radius)**2) then
                            if( imat(i,j,k) == icc ) mask(i,j,k) = .true.
                        endif
                    enddo
                enddo
            enddo
            ! Second likely center.
            new_center2 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), (imat == icc) .and. .not. mask)
            if( any(new_center2 > 0) )then ! if anything was found
                ! Validate second center (check if it's 2 merged atoms, or one pointy one)
                if( sum(real(new_center2 - new_center1)**2.) * self%smpd <= 2. * self%theoretical_radius) then
                    ! Set the merged cc back to 0
                    where( imat_cc == icc .and. (.not.mask) ) imat_cc = 0
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center2)
                    ! In the case two merged atoms, build the second atom
                    do i = 1, self%ldim(1)
                        do j = 1, self%ldim(2)
                            do k = 1, self%ldim(3)
                                if(((real(i - new_center2(1)))**2 + (real(j - new_center2(2)))**2 +&
                                &(real(k - new_center2(3)))**2) * self%smpd < (0.9 * self%theoretical_radius)**2 )then
                                    if( imat(i,j,k) == icc ) mask(i,j,k) = .true.
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif
            ! Third likely center.
            new_center3 = maxloc(rmat_pc(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), (imat == icc) .and. .not. mask)
            if( any(new_center3 > 0) )then ! if anything was found
                ! Validate third center
                if(sum(real(new_center3 - new_center1)**2.) * self%smpd <= 2. * self%theoretical_radius .or. &
                &  sum(real(new_center3 - new_center2)**2.) * self%smpd <= 2. * self%theoretical_radius )then
                    ! Set the merged cc back to 0
                    where( imat_cc == icc .and. (.not.mask) ) imat_cc = 0
                    return
                else
                    cnt = cnt + 1
                    new_centers(:,cnt) = real(new_center3)
                    ! In the case two merged atoms, build the second atom
                    do i = 1, self%ldim(1)
                        do j = 1, self%ldim(2)
                            do k = 1, self%ldim(3)
                                ! a little smaller to be sure
                                if( ((real(i - new_center3(1)))**2 + (real(j - new_center3(2)))**2 + &
                                &(real(k - new_center3(3)))**2) * self%smpd < (0.9 * self%theoretical_radius)**2) then
                                    if( imat(i,j,k) == icc ) mask(i,j,k) = .true.
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

    ! This subroutine discards outliers that resisted binarization
    ! It calculates the standard coordination number (cn) of each atom and discards
    ! the atoms with cn_std < cn_thresh
    ! It modifies the img_bin and img_cc instances deleting the identified outliers.
    subroutine discard_outliers( self, cn_thresh )
        class(nanoparticle), intent(inout) :: self
        integer,             intent(in)    :: cn_thresh ! threshold for discarding outliers based on coordination number
        integer, allocatable :: imat_bin(:,:,:), imat_cc(:,:,:)
        logical, allocatable :: mask(:)
        real, allocatable    :: centers_A(:,:) ! coordinates of the atoms in ANGSTROMS
        real    :: radius  ! radius of the sphere to consider for cn calculation
        real    :: a(3)    ! lattice parameter
        real    :: cn_gen(self%n_cc)
        integer :: cn(self%n_cc), cc, n_discard, filnum, filnum1, filnum2, i, nvox
        write(logfhandle, '(A)') '>>> DISCARDING OUTLIERS'
        centers_A = self%atominfo2centers_A()
        call fit_lattice(centers_A,a)
        call find_cn_radius(a,radius)
        call run_coord_number_analysis(centers_A,radius,cn,cn_gen)
        allocate(mask(self%n_cc), source=.true.)
        where( cn < cn_thresh ) mask = .false. ! false where atom has to be discarded
        n_discard = count(cn < cn_thresh)
        write(logfhandle, *) 'Numbers of atoms discarded because of low cn ', n_discard
        call self%img_cc%get_imat(imat_cc)
        call self%img_bin%get_imat(imat_bin)
        ! Removing outliers from the binary image and the connected components image
        ! remove atoms with < NVOX_THRESH voxels
        do cc = 1, self%n_cc
            nvox = count(imat_cc == cc)
            if( nvox < NVOX_THRESH )then
                where(imat_cc == cc) imat_bin = 0
            endif
        end do
        ! Removing outliers based on coordination number
        do cc = 1, self%n_cc
            if( cn(cc) < cn_thresh )then
                where(imat_cc == cc) imat_bin = 0
            endif
        enddo
        call self%img_bin%set_imat(imat_bin)
        call self%img_bin%find_ccs(self%img_cc)
        ! update number of connected components
        call self%img_cc%get_nccs(self%n_cc)
        call self%find_centers()
        deallocate(centers_A)
        centers_A = self%atominfo2centers_A()
        call run_coord_number_analysis(centers_A,radius,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        ! ATTENTION: you will see low coord numbers because they are UPDATED, after elimination
        ! of the atoms with low cn. It is like this in order to be consistent with the figure.
        write(logfhandle, '(A)') '>>> DISCARDING OUTLIERS, COMPLETED'
    end subroutine discard_outliers

    ! calc stats

    subroutine fillin_atominfo( self )
        class(nanoparticle), intent(inout) :: self
        type(image)          :: phasecorr
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: centers_A(:,:), tmpcens(:,:), strain_array(:,:)
        real,    pointer     :: rmat(:,:,:), rmat_corr(:,:,:)
        integer, allocatable :: imat_cc(:,:,:)
        character(len=256)   :: io_msg
        real    :: tmp_diam, a(3), radius_cn
        logical :: cc_mask(self%n_cc)
        integer :: i, j, k, cc, cn, n, funit, ios
        write(logfhandle, '(A)') '>>> EXTRACTING ATOM STATISTICS'
        ! calc cn and cn_gen
        centers_A = self%atominfo2centers_A()
        call fit_lattice(centers_A, a)
        call find_cn_radius(a, radius_cn)
        call run_coord_number_analysis(centers_A,radius_cn,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
        ! calc strain
        allocate(strain_array(self%n_cc,NSTRAIN_COMPS), source=0.)
        call strain_analysis(centers_A, a, strain_array)
        ! calc NPdiam & NPcen
        tmpcens     = self%atominfo2centers()
        cc_mask     = .true.
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
        call self%phasecorr_one_atom(phasecorr)
        call phasecorr%get_rmat_ptr(rmat_corr)
        call self%img%get_rmat_ptr(rmat)
        call self%img_cc%get_imat(imat_cc)
        do cc = 1, self%n_cc
            call progress(cc, self%n_cc)
            ! index of the connected component
            self%atominfo(cc)%cc_ind = cc
            ! number of voxels in connected component
            where( imat_cc == cc ) mask = .true.
            self%atominfo(cc)%size = count(mask)
            ! distance from the centre of mass of the nanoparticle
            self%atominfo(cc)%cendist = euclid(self%atominfo(cc)%center(:), self%NPcen) * self%smpd
            ! atom aspect ratio and diameter
            call self%calc_aspect_ratio(cc, self%atominfo(cc)%aspect_ratio, self%atominfo(cc)%diam)
            self%atominfo(cc)%diam = 2.*self%atominfo(cc)%diam ! radius --> diameter in A
            ! maximum grey level intensity across the connected component
            self%atominfo(cc)%max_int = maxval(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            ! average grey level intensity across the connected component
            self%atominfo(cc)%avg_int = sum(rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)), mask)
            self%atominfo(cc)%avg_int = self%atominfo(cc)%avg_int / real(count(mask))
            ! standard deviation of the grey level intensity across the connected component
            self%atominfo(cc)%sdev_int = 0.
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    do k = 1, self%ldim(3)
                        if( mask(i,j,k) ) self%atominfo(cc)%sdev_int = &
                        &self%atominfo(cc)%sdev_int + (rmat(i,j,k) - self%atominfo(cc)%avg_int)**2.
                    enddo
                enddo
            enddo
            if( count(mask) > 1 )then
                self%atominfo(cc)%sdev_int = sqrt(self%atominfo(cc)%sdev_int / real(count(mask)-1))
            else ! atom composed by one voxel
                self%atominfo(cc)%sdev_int = 0.
            endif
            ! atom correlation maximum within the connected component
            self%atominfo(cc)%max_corr = maxval(rmat_corr(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),mask)
            ! vector and angle of polarization
            if( self%atominfo(cc)%size > 2 )then
                self%atominfo(cc)%polar_vec   = real(self%atominfo(cc)%loc_ldist) - self%atominfo(cc)%center
                self%atominfo(cc)%polar_angle = ang3D_zvec(self%atominfo(cc)%polar_vec)
            else
                self%atominfo(cc)%polar_vec   = [0.,0.,0.]
                self%atominfo(cc)%polar_angle = 0. ! If the cc is too small it doesn't make sense
            endif
            ! bond length of nearest neighbour...
            self%atominfo(cc)%bondl = pixels_dist(self%atominfo(cc)%center(:), tmpcens, 'min', mask=cc_mask) ! Use all the atoms
            self%atominfo(cc)%bondl = self%atominfo(cc)%bondl * self%smpd ! convert to A
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
        ! CALCULATE GLOBAL NP PARAMETERS
        ! --global dipole
        self%net_dipole        = calc_net_dipole()
        self%net_dipole_mag    = arg(self%net_dipole)
        self%net_dipole_ang    = ang3D_zvec(self%net_dipole)
        ! --the rest
        call calc_stats( real(self%atominfo(:)%size),    self%size_stats, mask=self%atominfo(:)%size >= NVOX_THRESH )
        call calc_stats( real(self%atominfo(:)%cn_std),  self%cn_std_stats        )
        call calc_stats( self%atominfo(:)%bondl,         self%bondl_stats         )
        call calc_stats( self%atominfo(:)%cn_gen,        self%cn_gen_stats        )
        call calc_stats( self%atominfo(:)%aspect_ratio,  self%aspect_ratio_stats  )
        call calc_stats( self%atominfo(:)%polar_angle,   self%polar_angle_stats   )
        call calc_stats( self%atominfo(:)%diam,          self%diam_stats, mask=self%atominfo(:)%size >= NVOX_THRESH )
        call calc_stats( self%atominfo(:)%avg_int,       self%avg_int_stats       )
        call calc_stats( self%atominfo(:)%max_int,       self%max_int_stats       )
        call calc_stats( self%atominfo(:)%max_corr,      self%max_corr_stats      )
        call calc_stats( self%atominfo(:)%radial_strain, self%radial_strain_stats )
        ! CALCULATE CN-DEPENDENT STATS & WRITE CN-ATOMS
        do cn = CNMIN, CNMAX
            call calc_cn_stats( cn )
            call write_cn_atoms( cn )
        end do
        ! BINARY CLUSTERING OF RELEVANT PARAMETERS
        call fopen(funit, file=BIN_CLS_STATS_FILE, iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: bicluster_otsu; ERROR when opening file "//BIN_CLS_STATS_FILE//'; '//trim(io_msg),ios)
        ! write header
        write(funit,'(a)') BIN_CLS_STATS_HEAD
        call self%bicluster_otsu('size',          funit)
        call self%bicluster_otsu('bondl',         funit)
        call self%bicluster_otsu('aspect_ratio',  funit)
        call self%bicluster_otsu('polar_angle',   funit)
        call self%bicluster_otsu('diam',          funit)
        call self%bicluster_otsu('max_int',       funit)
        call self%bicluster_otsu('max_corr',      funit)
        call self%bicluster_otsu('radial_strain', funit)
        call fclose(funit)
        ! identify correlated variables with Pearson's product moment correation coefficient
        call self%id_corr_vars
        ! destruct
        deallocate(mask, imat_cc, tmpcens, strain_array, centers_A)
        call phasecorr%kill
        write(logfhandle, '(A)') '>>> EXTRACTING ATOM STATISTICS, COMPLETED'

        contains

            function ang3D_zvec( vec ) result( ang )
                real, intent(in) :: vec(3)
                real, parameter  :: ZVEC(3) = [0.,0.,1.]
                real :: vec_here(3)
                real :: ang     ! output angle in degrees
                real :: ang_rad ! angle in radians
                real :: dot_prod, sqsum
                ! normalise
                sqsum = vec(1)**2. + vec(2)**2. + vec(3)**2.
                if( sqsum > TINY )then
                    vec_here = vec / sqrt(sqsum)
                else
                    ang = 0.
                    return
                endif
                ! dot product
                dot_prod = ZVEC(1) * vec_here(1) + ZVEC(2) * vec_here(2) + ZVEC(3) * vec_here(3)
                ! sanity check
                if( dot_prod > 1. .or. dot_prod < -1. )then
                    THROW_WARN('Out of the domain of definition of arccos; fillin_atominfo :: ang3D_zvec')
                    ang_rad = 0.
                else
                    ang_rad = acos(dot_prod)
                endif
                ang = rad2deg(ang_rad)
            end function ang3D_zvec

            function calc_net_dipole( ) result( net_dipole )
                real    :: net_dipole(3), m_adjusted(3)
                integer :: i, io_stat, filnum
                net_dipole = 0. ! inizialization
                m_adjusted = 0.
                do i = 1, self%n_cc
                    net_dipole = net_dipole + self%atominfo(i)%polar_vec
                    m_adjusted = m_adjusted + (tmpcens(:,i) - 1.) * self%smpd
                enddo
                m_adjusted = m_adjusted / real(self%n_cc)
                net_dipole = (net_dipole - 1.) * self%smpd + m_adjusted
                write(logfhandle, *) 'Generating output global dipole *.bild file'
                call fopen(filnum, file='NetDipoleGlobal.bild', iostat=io_stat)
                write (filnum,'(a6,6i4)') trim('.arrow '), nint(m_adjusted), nint(net_dipole)
                call fclose(filnum)
            end function calc_net_dipole

            subroutine calc_cn_stats( cn )
                integer, intent(in)  :: cn ! calculate stats for given std cn
                integer :: cc, n, n_size, n_diam, natoms
                logical :: cn_mask(self%n_cc), size_mask(self%n_cc)
                real    :: vshell
                ! Generate masks
                cn_mask   = self%atominfo(:)%cn_std == cn
                size_mask = self%atominfo(:)%size >= NVOX_THRESH .and. cn_mask
                n         = count(cn_mask)
                if( n == 0 ) return
                ! -- cendist
                self%max_cendist_cns(cn) = maxval(self%atominfo(:)%cendist,  mask=cn_mask)
                self%min_cendist_cns(cn) = minval(self%atominfo(:)%cendist,  mask=cn_mask)
                ! -- # atoms
                self%natoms_cns(cn)       = real(n)
                ! -- # atoms per volume
                vshell = (FOURPI * self%max_cendist_cns(cn)**3.) / 3. - (FOURPI * self%min_cendist_cns(cn)**3.) / 3.
                natoms = count(self%atominfo(:)%cendist <= self%max_cendist_cns(cn) .and.&
                              &self%atominfo(:)%cendist >= self%min_cendist_cns(cn))
                self%atoms_per_vol_cns(cn) = 0.
                if( natoms > 0 .and. vshell > TINY ) self%atoms_per_vol_cns(cn) = real(natoms) / vshell
                if( n < 2 ) return
                ! -- the rest
                call calc_stats( real(self%atominfo(:)%size),    self%size_stats_cns(cn),          mask=size_mask )
                call calc_stats( self%atominfo(:)%bondl,         self%bondl_stats_cns(cn),         mask=cn_mask   )
                call calc_stats( self%atominfo(:)%cn_gen,        self%cn_gen_stats_cns(cn),        mask=cn_mask   )
                call calc_stats( self%atominfo(:)%aspect_ratio,  self%aspect_ratio_stats_cns(cn),  mask=cn_mask   )
                call calc_stats( self%atominfo(:)%polar_angle,   self%polar_angle_stats_cns(cn),   mask=cn_mask   )
                call calc_stats( self%atominfo(:)%diam,          self%diam_stats_cns(cn),          mask=size_mask )
                call calc_stats( self%atominfo(:)%avg_int,       self%avg_int_stats_cns(cn),       mask=cn_mask   )
                call calc_stats( self%atominfo(:)%max_int,       self%max_int_stats_cns(cn),       mask=cn_mask   )
                call calc_stats( self%atominfo(:)%max_corr,      self%max_corr_stats_cns(cn),      mask=cn_mask   )
                call calc_stats( self%atominfo(:)%radial_strain, self%radial_strain_stats_cns(cn), mask=cn_mask   )
            end subroutine calc_cn_stats

            subroutine write_cn_atoms( cn_std )
                integer, intent(in)  :: cn_std
                type(binimage)       :: img_atom
                integer, allocatable :: imat(:,:,:), imat_atom(:,:,:)
                logical :: cn_mask(self%n_cc)
                integer :: i
                ! make cn mask
                cn_mask = self%atominfo(:)%cn_std == cn_std
                ! make binary matrix of atoms with given cn_std
                call img_atom%copy_bimg(self%img_cc)
                allocate(imat_atom(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
                call img_atom%get_imat(imat)
                do i = 1, self%n_cc
                    if( cn_mask(i) )then
                        where( imat == i ) imat_atom = 1
                    endif
                enddo
                call img_atom%set_imat(imat_atom)
                call img_atom%write_bimg('Atoms_with_cn_std'//trim(int2str(cn_std))//'.mrc')
                deallocate(imat,imat_atom)
                call img_atom%kill_bimg
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

    ! This subroutine takes in input a connected component (cc) image
    ! and the label of one of its ccs and calculates the aspect ratio of the cc,
    ! defined as the ratio of the width and the height.
    ! The idea behind this is that the center of the cc is calculated,
    ! then everything is deleted except the borders of the cc. Finally,
    ! in order to calculate the width and the height, the min/max
    ! distances between the center and the borders are calculated. The
    ! aspect ratio is the ratio of those 2 distances.
    subroutine calc_aspect_ratio( self, label, ratio, longest_dist )
       class(nanoparticle), intent(inout) :: self
       integer,             intent(in)    :: label
       real,                intent(out)   :: ratio
       real,                intent(out)   :: longest_dist
       integer, allocatable :: pos(:,:)
       integer, allocatable :: imat_cc(:,:,:)
       logical, allocatable :: border(:,:,:)
       logical, allocatable :: mask_dist(:) ! for min and max dist calculation
       integer :: location(1)               ! location of vxls of the atom farthest from its center
       real    :: shortest_dist
       call self%img_cc%get_imat(imat_cc)
       call self%img_cc%border_mask(border, label, .true.) ! use 4neigh instead of 8neigh
       where(border .eqv. .true.)
           imat_cc = 1
       elsewhere
           imat_cc = 0
       endwhere
       call get_pixel_pos( imat_cc, pos ) ! pxls positions of the shell
       allocate(mask_dist(size(pos, dim=2)), source = .true.)
       shortest_dist = pixels_dist(self%atominfo(label)%center(:), real(pos), 'min', mask_dist, location)
       if( size(pos,2) == 1 ) then ! if the connected component has size 1 (just 1 vxl)
           shortest_dist = 0.
           longest_dist  = shortest_dist
           ratio         = 1.
           self%atominfo(label)%loc_ldist(:) = nint(self%atominfo(label)%center(:))
           return
       else
           longest_dist  = pixels_dist(self%atominfo(label)%center(:), real(pos),'max', mask_dist, location)
           self%atominfo(label)%loc_ldist(:) =  pos(:,location(1))
       endif
       if( abs(longest_dist) > TINY .and. size(pos,2) > 1 )then
           ratio = shortest_dist/longest_dist
       else
           ratio = 0.
       endif
       longest_dist  = longest_dist  * self%smpd  ! in A
       shortest_dist = shortest_dist * self%smpd
       deallocate(imat_cc, border, pos, mask_dist)
    end subroutine calc_aspect_ratio

    ! visualization and output

    subroutine write_centers(self, fname, coords)
       class(nanoparticle),        intent(inout) :: self
       character(len=*), optional, intent(in)    :: fname
       real,             optional, intent(in)    :: coords(:,:)
       type(image) :: simulated_distrib
       type(atoms) :: atom
       integer     :: cc
       if( present(coords) )then
           call self%centers_pdb%new(size(coords, dim = 2), dummy=.true.)
           do cc=1,size(coords, dim = 2)
               call self%centers_pdb%set_name(cc,self%atom_name)
               call self%centers_pdb%set_element(cc,self%element)
               call self%centers_pdb%set_coord(cc,(coords(:,cc)-1.)*self%smpd)
           enddo
       else
           call self%centers_pdb%new(self%n_cc, dummy=.true.)
           do cc=1,self%n_cc
               call self%centers_pdb%set_name(cc,self%atom_name)
               call self%centers_pdb%set_element(cc,self%element)
               call self%centers_pdb%set_coord(cc,(self%atominfo(cc)%center(:)-1.)*self%smpd)
               call self%centers_pdb%set_beta(cc,self%atominfo(cc)%cn_gen) ! use generalised coordination number
               call self%centers_pdb%set_resnum(cc,cc)
           enddo
       endif
       if( present(fname) ) then
           call self%centers_pdb%writepdb(fname)
       else
           call self%centers_pdb%writepdb(trim(self%fbody)//'_atom_centers')
           call atom%new                 (trim(self%fbody)//'_atom_centers.pdb')
           call simulated_distrib%new(self%ldim,self%smpd)
           call atom%convolve(simulated_distrib, cutoff = 8.*self%smpd)
           call simulated_distrib%write(trim(self%fbody)//'SIM.mrc')
           write(logfhandle,*) 'output, atomic coordinates:       ', trim(self%fbody)//'_atom_centers.pdb'
           write(logfhandle,*) 'output, simulated atomic density: ', trim(self%fbody)//'SIM.mrc'
           call atom%kill
           call simulated_distrib%kill
       endif
    end subroutine write_centers

    subroutine write_atoms( self )
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
    end subroutine write_atoms

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
        call fopen(funit, file=ATOM_STATS_FILE, iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: write_csv_files; ERROR when opening file "//ATOM_STATS_FILE//'; '//trim(io_msg),ios)
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
        601 format(F8.4,A3)
        602 format(F8.4)
        ! various per-atom parameters
        write(funit,601,advance='no') real(self%atominfo(cc)%cc_ind),           CSV_DELIM ! INDEX
        write(funit,601,advance='no') real(self%atominfo(cc)%size),             CSV_DELIM ! NVOX
        write(funit,601,advance='no') real(self%atominfo(cc)%cn_std),           CSV_DELIM ! CN_STD
        write(funit,601,advance='no') self%atominfo(cc)%bondl,                  CSV_DELIM ! NN_BONDL
        write(funit,601,advance='no') self%atominfo(cc)%cn_gen,                 CSV_DELIM ! CN_GEN
        write(funit,601,advance='no') self%atominfo(cc)%aspect_ratio,           CSV_DELIM ! ASPECT_RATIO
        write(funit,601,advance='no') self%atominfo(cc)%polar_angle,            CSV_DELIM ! POLAR_ANGLE
        write(funit,601,advance='no') self%atominfo(cc)%diam,                   CSV_DELIM ! DIAM
        write(funit,601,advance='no') self%atominfo(cc)%avg_int,                CSV_DELIM ! AVG_INT
        write(funit,601,advance='no') self%atominfo(cc)%max_int,                CSV_DELIM ! MAX_INT
        write(funit,601,advance='no') self%atominfo(cc)%sdev_int,               CSV_DELIM ! SDEV_INT
        write(funit,601,advance='no') self%atominfo(cc)%cendist,                CSV_DELIM ! RADIAL_POS
        write(funit,601,advance='no') self%atominfo(cc)%max_corr,               CSV_DELIM ! MAX_CORR
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
        if( .not. omit_here )then
        ! cluster assignments
        write(funit,601,advance='no') real(self%atominfo(cc)%size_cls),         CSV_DELIM ! NVOX_CLASS
        write(funit,601,advance='no') real(self%atominfo(cc)%bondl_cls),        CSV_DELIM ! NN_BONDL_CLASS
        write(funit,601,advance='no') real(self%atominfo(cc)%aspect_ratio_cls), CSV_DELIM ! ASPECT_RATIO_CLASS
        write(funit,601,advance='no') real(self%atominfo(cc)%polar_angle_cls),  CSV_DELIM ! POLAR_ANGLE_CLASS
        write(funit,601,advance='no') real(self%atominfo(cc)%diam_cls),         CSV_DELIM ! DIAM_CLASS
        write(funit,601,advance='no') real(self%atominfo(cc)%max_int_cls),      CSV_DELIM ! MAX_INT_CLASS
        write(funit,601,advance='no') real(self%atominfo(cc)%max_corr_cls),     CSV_DELIM ! MAX_CORR_CLASS
        write(funit,602)              real(self%atominfo(cc)%radial_strain_cls)           ! RADIAL_STRAIN_CLASS
        endif
    end subroutine write_atominfo

    subroutine write_np_stats( self, funit )
        class(nanoparticle), intent(in) :: self
        integer,             intent(in) :: funit
        601 format(F8.4,A3)
        602 format(F8.4)
        ! -- # atoms
        write(funit,601,advance='no') real(self%n_cc),               CSV_DELIM ! NATOMS
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
        ! -- aspect ratio
        write(funit,601,advance='no') self%aspect_ratio_stats%avg,   CSV_DELIM ! AVG_ASPECT_RATIO
        write(funit,601,advance='no') self%aspect_ratio_stats%med,   CSV_DELIM ! MED_ASPECT_RATIO
        write(funit,601,advance='no') self%aspect_ratio_stats%sdev,  CSV_DELIM ! SDEV_ASPECT_RATIO
        ! -- polarization angle
        write(funit,601,advance='no') self%polar_angle_stats%avg,    CSV_DELIM ! AVG_POLAR_ANGLE
        write(funit,601,advance='no') self%polar_angle_stats%med,    CSV_DELIM ! MED_POLAR_ANGLE
        write(funit,601,advance='no') self%polar_angle_stats%sdev,   CSV_DELIM ! SDEV_POLAR_ANGLE
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
        write(funit,601,advance='no') self%max_corr_stats%avg,       CSV_DELIM ! AVG_MAX_CORR
        write(funit,601,advance='no') self%max_corr_stats%med,       CSV_DELIM ! MED_MAX_CORR
        write(funit,601,advance='no') self%max_corr_stats%sdev,      CSV_DELIM ! SDEV_MAX_CORR
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
        601 format(F8.4,A3)
        602 format(F8.4)
        if( count(self%atominfo(:)%cn_std == cn) < 2 ) return
        ! -- coordination number
        write(funit,601,advance='no') real(cn),                              CSV_DELIM ! CN_STD
        ! -- # atoms per cn
        write(funit,601,advance='no') self%natoms_cns(cn),                   CSV_DELIM ! NATOMS
        ! -- # atoms per volume (#/A**3)
        write(funit,601,advance='no') self%atoms_per_vol_cns(cn),            CSV_DELIM ! NATOMS_PER_VOL
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
        ! -- aspect ratio
        write(funit,601,advance='no') self%aspect_ratio_stats_cns(cn)%avg,   CSV_DELIM ! AVG_ASPECT_RATIO
        write(funit,601,advance='no') self%aspect_ratio_stats_cns(cn)%med,   CSV_DELIM ! MED_ASPECT_RATIO
        write(funit,601,advance='no') self%aspect_ratio_stats_cns(cn)%sdev,  CSV_DELIM ! SDEV_ASPECT_RATIO
        ! -- polarization angle
        write(funit,601,advance='no') self%polar_angle_stats_cns(cn)%avg,    CSV_DELIM ! AVG_POLAR_ANGLE
        write(funit,601,advance='no') self%polar_angle_stats_cns(cn)%med,    CSV_DELIM ! MED_POLAR_ANGLE
        write(funit,601,advance='no') self%polar_angle_stats_cns(cn)%sdev,   CSV_DELIM ! SDEV_POLAR_ANGLE
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
        write(funit,601,advance='no') self%max_corr_stats_cns(cn)%avg,       CSV_DELIM ! AVG_MAX_CORR
        write(funit,601,advance='no') self%max_corr_stats_cns(cn)%med,       CSV_DELIM ! MED_MAX_CORR
        write(funit,601,advance='no') self%max_corr_stats_cns(cn)%sdev,      CSV_DELIM ! SDEV_MAX_CORR
        ! -- radial strain
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%avg,  CSV_DELIM ! AVG_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%med,  CSV_DELIM ! MED_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%sdev, CSV_DELIM ! SDEV_RADIAL_STRAIN
        write(funit,601,advance='no') self%radial_strain_stats_cns(cn)%minv, CSV_DELIM ! MIN_RADIAL_STRAIN
        write(funit,602)              self%radial_strain_stats_cns(cn)%maxv            ! MAX_RADIAL_STRAIN
    end subroutine write_cn_stats

    ! identify correlated variables with Pearson's product moment correation coefficient
    subroutine id_corr_vars( self )
        class(nanoparticle), target, intent(in) :: self
        integer, parameter :: NFLAGS = 12
        integer, parameter :: NPAIRS = (NFLAGS * (NFLAGS - 1)) / 2
        character(len=13)  :: flags(NFLAGS), flags_i(NPAIRS), flags_j(NPAIRS)
        character(len=256) :: io_msg
        real, pointer      :: ptr1(:), ptr2(:)
        real, allocatable, target :: sizes(:)
        real    :: corrs(NPAIRS), corrs_copy(NPAIRS), corr
        integer :: i, j, inds(NPAIRS), cnt, funit, ios
        ! make sizes a real array
        allocate(sizes(self%n_cc), source=real(self%atominfo(:)%size))
        ! variables to correlate
        flags(1)  = 'NVOX'          ! size
        flags(2)  = 'NN_BONDL'      ! bondl
        flags(3)  = 'CN_GEN'        ! cn_gen
        flags(4)  = 'ASPECT_RATIO'  ! aspect_ratio
        flags(5)  = 'POLAR_ANGLE'   ! polar_angle
        flags(6)  = 'DIAM'          ! diam
        flags(7)  = 'AVG_INT'       ! avg_int
        flags(8)  = 'MAX_INT'       ! max_int
        flags(9)  = 'SDEV_INT'      ! sdev_int
        flags(10) = 'RADIAL_POS'    ! cendist
        flags(11) = 'MAX_CORR'      ! max_corr
        flags(12) = 'RADIAL_STRAIN' ! radial_strain
        ! calculate correlations
        cnt = 0
        do i = 1, NFLAGS - 1
            do j = i + 1, NFLAGS
                cnt          = cnt + 1
                inds(cnt)    = cnt
                flags_i(cnt) = flags(i)
                flags_j(cnt) = flags(j)
                call set_ptr(flags(i), ptr1)
                call set_ptr(flags(j), ptr2)
                corrs(cnt) = pearsn_serial(ptr1, ptr2)
            end do
        end do
        ! sort
        corrs_copy = corrs
        call hpsort(corrs_copy, inds)
        ! write output
        call fopen(funit, file=ATOM_VAR_CORRS_FILE, iostat=ios, status='replace', iomsg=io_msg)
        call fileiochk("simple_nanoparticle :: id_corr_vars; ERROR when opening file "//ATOM_VAR_CORRS_FILE//'; '//trim(io_msg),ios)
        do i = 1, NPAIRS
            write(funit,'(A,F7.4)')'PEARSONS CORRELATION BTW '//flags_i(inds(i))//' & '//flags_j(inds(i))//' IS ', corrs(inds(i))
        end do
        call fclose(funit)

        contains

            subroutine set_ptr( flag, ptr )
                character(len=*), intent(in)    :: flag
                real, pointer,    intent(inout) :: ptr(:)
                 select case(trim(flag))
                    case('NVOX')
                        ptr => sizes
                    case('NN_BONDL')
                        ptr => self%atominfo(:)%bondl
                    case('CN_GEN')
                        ptr => self%atominfo(:)%cn_gen
                    case('ASPECT_RATIO')
                        ptr => self%atominfo(:)%aspect_ratio
                    case('POLAR_ANGLE')
                        ptr => self%atominfo(:)%polar_angle
                    case('DIAM')
                        ptr => self%atominfo(:)%diam
                    case('AVG_INT')
                        ptr => self%atominfo(:)%avg_int
                    case('MAX_INT')
                        ptr => self%atominfo(:)%max_int
                    case('SDEV_INT')
                        ptr => self%atominfo(:)%sdev_int
                    case('RADIAL_POS')
                        ptr => self%atominfo(:)%cendist
                    case('MAX_CORR')
                        ptr => self%atominfo(:)%max_corr
                    case('RADIAL_STRAIN')
                        ptr => self%atominfo(:)%radial_strain
                end select
            end subroutine set_ptr

    end subroutine id_corr_vars

    subroutine bicluster_otsu( self, which, funit )
        class(nanoparticle), intent(inout) :: self
        character(len=*),    intent(in)    :: which
        integer,             intent(in)    :: funit
        integer :: n, bicls(self%n_cc)
        real    :: vals(self%n_cc), thresh
        logical :: cls_mask(self%n_cc), mask(self%n_cc)
        select case(which)
            case('size')
                vals = real(self%atominfo(:)%size)
                call bicluster_local(self%atominfo(:)%size >= NVOX_THRESH, which, 'NVOX')
                self%atominfo(:)%size_cls = bicls
            case('bondl')
                vals = self%atominfo(:)%bondl
                call bicluster_local(self%atominfo(:)%bondl > 0., which, 'NN_BONDL')
                self%atominfo(:)%bondl_cls = bicls
            case('aspect_ratio')
                vals = self%atominfo(:)%aspect_ratio
                call bicluster_local(self%atominfo(:)%aspect_ratio > 0., which, 'ASPECT_RATIO')
                self%atominfo(:)%aspect_ratio_cls = bicls
            case('polar_angle')
                vals = self%atominfo(:)%polar_angle
                mask = .true.
                call bicluster_local(mask, which, 'POLAR_ANGLE')
                self%atominfo(:)%polar_angle_cls = bicls
            case('diam')
                vals = self%atominfo(:)%diam
                call bicluster_local(self%atominfo(:)%diam > 0., which, 'DIAM')
                self%atominfo(:)%diam_cls = bicls
            case('max_int')
                vals = self%atominfo(:)%max_int
                mask = .true.
                call bicluster_local(mask, which, 'MAX_INT')
                self%atominfo(:)%max_int_cls = bicls
            case('max_corr')
                vals = self%atominfo(:)%max_corr
                mask = .true.
                call bicluster_local(mask, which, 'MAX_CORR')
                self%atominfo(:)%max_corr_cls = bicls
            case('radial_strain')
                vals = self%atominfo(:)%radial_strain
                mask = .true.
                call bicluster_local(mask, which, 'RADIAL_STRAIN')
                self%atominfo(:)%radial_strain_cls = bicls
        case DEFAULT
            THROW_HARD('unsupported parameter for bicluster_otsu')
        end select

        contains

            subroutine bicluster_local( mask, which, param )
                logical,          intent(in) :: mask(self%n_cc)
                character(len=*), intent(in) :: which, param
                real, allocatable  :: vals4otsu(:)
                type(stats_struct) :: stats(2)
                character(len=13)  :: param_str
                integer            :: pop1, pop2
                601 format(F8.4,A3)
                602 format(F8.4)
                param_str = trim(param)
                vals4otsu = pack(vals, mask=mask)
                call otsu(vals4otsu, thresh)
                deallocate(vals4otsu)
                where( vals > thresh )
                    bicls = 1 ! 1 is large
                elsewhere
                    bicls = 2 ! 2 is small
                endwhere
                cls_mask = bicls == 1 .and. mask
                pop1     = count(cls_mask)
                call calc_stats(vals, stats(1), mask=cls_mask)
                call write_cls_atoms(1, which)
                cls_mask = bicls == 2 .and. mask
                pop2     = count(cls_mask)
                call calc_stats(vals, stats(2), mask=cls_mask)
                call write_cls_atoms(2, which)
                write(logfhandle,'(A,A21,1X,F8.4,1X,F8.4)') param_str, ' class averages 1 & 2:', stats(1)%avg, stats(2)%avg
                write(funit,'(A,A3)',advance='no') param_str,  CSV_DELIM ! PARAMETER
                write(funit,601,  advance='no') stats(1)%avg,  CSV_DELIM ! AVG_1
                write(funit,601,  advance='no') stats(1)%med,  CSV_DELIM ! MED_1
                write(funit,601,  advance='no') stats(1)%sdev, CSV_DELIM ! SDEV_1
                write(funit,601,  advance='no') stats(1)%minv, CSV_DELIM ! MIN_1
                write(funit,601,  advance='no') stats(1)%maxv, CSV_DELIM ! MAX_1
                write(funit,601,  advance='no') real(pop1),    CSV_DELIM ! POP_1
                write(funit,601,  advance='no') stats(2)%avg,  CSV_DELIM ! AVG_2
                write(funit,601,  advance='no') stats(2)%med,  CSV_DELIM ! MED_2
                write(funit,601,  advance='no') stats(2)%sdev, CSV_DELIM ! SDEV_2
                write(funit,601,  advance='no') stats(2)%minv, CSV_DELIM ! MIN_2
                write(funit,601,  advance='no') stats(2)%maxv, CSV_DELIM ! MAX_2
                write(funit,602)                real(pop2)               ! POP_2
            end subroutine bicluster_local

            subroutine write_cls_atoms( cls, which )
                integer,          intent(in) :: cls
                character(len=*), intent(in) :: which
                type(binimage)       :: img_atom
                integer, allocatable :: imat(:,:,:), imat_atom(:,:,:)
                integer :: i
                call img_atom%copy_bimg(self%img_cc)
                allocate(imat_atom(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
                call img_atom%get_imat(imat)
                do i = 1, self%n_cc
                    if( cls_mask(i) )then
                        where( imat == i ) imat_atom = 1
                    endif
                enddo
                call img_atom%set_imat(imat_atom)
                if( cls == 1 )then
                    call img_atom%write_bimg('Atoms_'//trim(which)//'_large.mrc')
                else
                    call img_atom%write_bimg('Atoms_'//trim(which)//'_small.mrc')
                endif
                deallocate(imat,imat_atom)
                call img_atom%kill_bimg
            end subroutine write_cls_atoms

    end subroutine bicluster_otsu

    subroutine cluster_atom_features( self )
        use simple_aff_prop, only: aff_prop
        class(nanoparticle), intent(inout) :: self
        integer, allocatable :: centers(:), labels(:)
        type(aff_prop) :: ap
        real           :: smat(self%n_cc,self%n_cc), simsum, datvecs(self%n_cc,8)
        integer        :: i, j, cc
        do cc = 1, self%n_cc
            ! real-valued data
            datvecs(cc,1) = real(self%atominfo(cc)%size)
            datvecs(cc,2) =      self%atominfo(cc)%bondl
            datvecs(cc,3) =      self%atominfo(cc)%aspect_ratio
            datvecs(cc,4) =      self%atominfo(cc)%polar_angle
            datvecs(cc,5) =      self%atominfo(cc)%diam
            datvecs(cc,6) =      self%atominfo(cc)%max_int
            datvecs(cc,7) =      self%atominfo(cc)%max_corr
            datvecs(cc,8) =      self%atominfo(cc)%radial_strain
        end do
        ! rescaling to avoid the polarization angle to dominate the analysis
        call norm_minmax(datvecs(:,1))
        call norm_minmax(datvecs(:,2))
        call norm_minmax(datvecs(:,3))
        call norm_minmax(datvecs(:,4))
        call norm_minmax(datvecs(:,5))
        call norm_minmax(datvecs(:,6))
        call norm_minmax(datvecs(:,7))
        call norm_minmax(datvecs(:,8))
        ! calculate similarity matrix
        do i = 1, self%n_cc - 1
            do j = i + 1, self%n_cc
                ! smat(i,j) =   pearsn_serial(datvecs(i,:), datvecs(j,:))
                smat(i,j) = - euclid(datvecs(i,:), datvecs(j,:))
                smat(j,i) = smat(i,j)
            end do
        end do
        ! affinity propagation
        call ap%new(self%n_cc, smat)
        call ap%propagate(centers, labels, simsum)
        write(logfhandle,*) '# clusters found with affinity propagation', size(centers)

        contains

            subroutine norm_minmax( arr )
                real, intent(inout) :: arr(:)
                real                :: smin, smax, delta
                smin  = minval(arr)
                smax  = maxval(arr)
                delta = smax - smin
                arr = (arr - smin)/delta
            end subroutine norm_minmax

    end subroutine cluster_atom_features

    ! This subroutine clusters the atoms with respect to the maximum intensity
    ! or the integrated density (according to the values contained in feature)
    ! using kmeans algorithm with 2 classes. The initial guess fo the centers
    ! is intentionally biased. It supposes there are two distinguished classes
    ! with different avgs (proved with simulated data).
    subroutine cluster_atom_intensity( self, feature )
        use gnufor2
        use simple_nanoML, only : nanoML
        class(nanoparticle), intent(inout) :: self
        real,                intent(inout) :: feature(:)
        integer, parameter   :: MAX_IT = 50 ! maximum number of iterations for
        integer, allocatable :: imat_cc(:,:,:)
        real, pointer        :: rmat1(:,:,:), rmat2(:,:,:)
        type(image)          :: class1, class2
        type(nanoML)         :: emfit
        real    :: centers_kmeans(2) ! output of k-means
        real    :: avgs(2),  vars(2), gammas(self%n_cc,2) ! output of ML
        integer :: i, cnt1, cnt2, filnum, io_stat
        feature = feature/maxval(feature)*10. ! normalise with maxval*10
        call hist(feature, 20)
        write(logfhandle,*) '****clustering wrt maximum intensity, init'
        ! Report clusters on images in dedicated directory
        call class1%new(self%ldim, self%smpd)
        call class2%new(self%ldim, self%smpd)
        call class1%get_rmat_ptr(rmat1)
        call class2%get_rmat_ptr(rmat2)
        call self%img_cc%get_imat(imat_cc)
        call fopen(filnum, file='ClusterIntensities.txt', iostat=io_stat)
        ! kmeans
        call emfit%kmeans_biased2classes(feature, centers_kmeans)
        write(filnum,*) 'centers_kmeans', centers_kmeans
        ! ML, fit
        call emfit%new(self%n_cc,2)
        call emfit%set_data(feature)
        call emfit%fit(MAX_IT,centers_kmeans)
        avgs   = emfit%get_avgs()
        vars   = emfit%get_vars()
        gammas = emfit%get_gammas()
        write(filnum,*) 'AVG/VAR 1:', avgs(1), vars(1)
        write(filnum,*) 'AVG/VAR 2:', avgs(2), vars(2)
        cnt2 = 0
        do i = 1, self%n_cc
            if( (avgs(1) - feature(i))**2. < (avgs(2) - feature(i))**2. ) then
                cnt2 = cnt2 + 1
                write(filnum,*) 'connected component #', i, 'belongs to class 1 with probability', max(gammas(i,1),gammas(i,2))
              else
                write(filnum,*) 'connected component #', i, 'belongs to class 2 with probability', max(gammas(i,1),gammas(i,2))
            endif
        enddo
        cnt1 = count((avgs(1) - feature)**2. <  (avgs(2) - feature)**2. )
        cnt2 = count((avgs(1) - feature)**2. >= (avgs(2) - feature)**2. )
        do i = 1, self%n_cc
            if( (avgs(1) - feature(i))**2. < (avgs(2) - feature(i))**2. ) then
                where( imat_cc == i ) rmat1(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
            else
                where( imat_cc == i ) rmat2(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
            endif
        enddo
        call class1%write('Class1.mrc')
        call class2%write('Class2.mrc')
        write(filnum,*)  'Class1 contains ', cnt1 ,' atoms'
        write(filnum,*)  'Class2 contains ', cnt2 ,' atoms'
        write(filnum,*)  'Total ', cnt1+cnt2, ' atoms'
        call fclose(filnum)
        call class1%kill
        call class2%kill
        deallocate(imat_cc)
        ! come back to root directory
        write(logfhandle,*) '****clustering wrt maximum intensity, completed'

    contains

        subroutine initialise_centers(data,cen1,cen2)
            real, intent(inout) :: cen1,cen2
            real, intent(inout) :: data(:)
            !>   rheapsort from numerical recepies (largest last)
            call hpsort(data)
            cen1 = sum(data(1:self%n_cc/2))/real(self%n_cc/2)
            cen2 = sum(data(self%n_cc/2+1:size(data)))/real(self%n_cc/2)
        end subroutine initialise_centers

        subroutine update_centers(cen1,cen2,converged,val1,val2)
            real,    intent(inout) :: cen1,cen2
            logical, intent(inout) :: converged
            integer, intent(inout) :: val1, val2
            integer :: i
            integer :: cnt1, cnt2
            real    :: sum1, sum2
            real :: cen1_new, cen2_new
            sum1 = 0.
            cnt1 = 0
            do i=1,self%n_cc
                if( (cen1-feature(i))**2. < (cen2-feature(i))**2. )then
                    cnt1 = cnt1 + 1 ! number of elements in cluster 1
                    sum1 = sum1 + feature(i)
                endif
            end do
            cnt2 = self%n_cc - cnt1       ! number of elements in cluster 2
            sum2 = sum(feature)- sum1
            cen1_new = sum1 / real(cnt1)
            cen2_new = sum2 / real(cnt2)
            if(abs(cen1_new - cen1) < TINY .and. abs(cen2_new - cen2) < TINY) then
                converged = .true.
            else
                converged = .false.
            endif
            ! update
            cen1 = cen1_new
            cen2 = cen2_new
            ! assign values to the centers
            if( cen1 > cen2 )then
                val1           = 1
                val2           = 0
            else
                val1           = 0
                val2           = 1
            endif
        end subroutine update_centers

    end subroutine cluster_atom_intensity

    ! Cluster the atoms wrt the maximum intensity
    ! k-means 2 classes + ML
    subroutine cluster_atom_maxint( self )
        class(nanoparticle), intent(inout) :: self
        real    :: max_intensity(self%n_cc)
        integer :: i, io_stat, filnum
        call fopen(filnum, file='MaxIntensity.csv', action='readwrite', iostat=io_stat)
        if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file MaxIntensity.csv Did you run atoms_stats?; cluster_atom_maxint')
        endif
        read(filnum,*) ! first line is variable name
        do i = 1, self%n_cc
            read(filnum,*) max_intensity(i)
        enddo
        call fclose(filnum)
        call self%cluster_atom_intensity(max_intensity)
    end subroutine cluster_atom_maxint

    ! Cluster the atoms wrt the integrated density
    ! k-means 2 classes + ML
    subroutine cluster_atom_intint( self )
        class(nanoparticle), intent(inout) :: self
        real    :: int_intensity(self%n_cc)
        integer :: i, io_stat, filnum
        call fopen(filnum, file='IntIntensity.csv', action='readwrite', iostat=io_stat)
        if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file IntIntensity.csv Did you run atoms_stats?; cluster_atom_intint')
        endif
        read(filnum,*) ! first line is variable name
        do i = 1, self%n_cc
            read(filnum,*) int_intensity(i)
        enddo
        call fclose(filnum)
        call self%cluster_atom_intensity(int_intensity)
    end subroutine cluster_atom_intint

    ! Cluster wrt to the angle between the direction of the longest dimension and a fixed vector
    subroutine cluster_ang(self, thresh)
        class(nanoparticle), intent(inout) :: self
        real,    intent(in)  :: thresh ! threshold for class definition, user inputted
        real,    allocatable :: stdev_within(:), centroids(:)
        integer, allocatable :: labels(:), populations(:)
        integer              :: i, j, ncls, dim, filnum, io_stat
        real                 :: avg, stdev
        type(binimage)       :: img_1clss
        integer, allocatable :: imat_cc(:,:,:), imat_1clss(:,:,:)
        character(len=4)     :: str_thres
        ! Preparing for clustering
        ! Aspect ratios and loc_longest dist estimation
        ! call self%est_aspect_ratios(print_ar=.false.) ! it's needed to identify the dir of longest dim
        ! call self%search_polarization(output_files =.false.)
        dim = size(self%atominfo)
        allocate(labels(dim), source = 0)
        call hac_1d(self%atominfo(:)%polar_angle,thresh,labels,centroids,populations)
        ! Stats calculations
        ncls = maxval(labels)
        allocate(stdev_within(ncls), source = 0.)
        !stdev within the same class
        do i = 1, ncls
            do j = 1, dim
                if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%atominfo(j)%polar_angle-centroids(i))**2
            enddo
        enddo
        where (populations>1)
            stdev_within = sqrt(stdev_within/(real(populations)-1.))
        elsewhere
            stdev_within = 0.
        endwhere
        ! avg and stdev among different classes
        avg = 0.
        do i = 1, ncls
            avg = avg + centroids(i)
        enddo
        avg = avg/real(ncls)
        stdev = 0.
        if(ncls>1) then
          do i = 1, ncls
              stdev = stdev + (centroids(i) - avg)**2
          enddo
          stdev = sqrt(stdev/(real(ncls)-1.))
        endif
        ! Output on a file
        str_thres = trim(real2str(thresh))
        call fopen(filnum, file='ClusterAngThresh'//str_thres//'.txt', iostat=io_stat)
        write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
        write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
        do i = 1, dim
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3)') 'Atom #: ', i, '; data (deg): ', self%atominfo(i)%polar_angle, '; class: ', labels(i)
        enddo
        write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
        do i = 1, ncls
            write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), ' degrees; stdev within the class: ', stdev_within(i), ' degrees'
        enddo
        write(unit = filnum,fmt ='(a,f6.2,a,f6.2,a)') 'AVG among the classes: ', avg, ' degrees; STDEV among the classes: ', stdev, ' degrees'
        call fclose(filnum)
        call fopen(filnum, file='Ang.csv', iostat=io_stat)
        write(filnum,*) 'ang'
        do i  = 1, self%n_cc
          write(filnum,*) self%atominfo(i)%polar_angle
        enddo
        call fclose(filnum)
        !Generate one figure for each class
        if(GENERATE_FIGS) then
          call img_1clss%new_bimg(self%ldim, self%smpd)
          call self%img_cc%get_imat(imat_cc)
          allocate(imat_1clss(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
          do i = 1, ncls
              imat_1clss = 0
              do j = 1, dim
                  if(labels(j) == i) then
                      where(imat_cc == j) imat_1clss = i
                  endif
              enddo
              call img_1clss%set_imat(imat_1clss)
              call img_1clss%write_bimg('AngClass'//int2str(i)//'.mrc')
          enddo
          call img_1clss%kill_bimg
        endif
        if(allocated(stdev_within)) deallocate(stdev_within)
        deallocate(centroids, labels, populations)
    end subroutine cluster_ang

   ! Cluster the atoms wrt to the aspect ratio
   ! Additionally, calculate the median std cn
   ! (coordination number) in each class
   subroutine cluster_ar(self, thresh)
      class(nanoparticle), intent(inout) :: self
      real,    intent(in)  :: thresh ! threshold for class definition, user inputted
      real,    allocatable :: centroids(:)
      integer, allocatable :: labels(:), populations(:)
      real,    allocatable :: stdev_within(:), centers_A(:,:), cn_cls(:), median_cn(:)
      integer              :: i, j, ncls, dim, filnum, io_stat, cnt
      real                 :: avg, stdev, radius_cn, a(3)
      type(binimage)       :: img_1clss
      integer, allocatable :: imat_cc(:,:,:), imat_1clss(:,:,:)
      character(len=4)     :: str_thres
      ! Calculate cn and cn_gen
      centers_A = self%atominfo2centers_A()
      call fit_lattice(centers_A,a)
      call find_cn_radius(a,radius_cn)
      call run_coord_number_analysis(centers_A,radius_cn,self%atominfo(:)%cn_std,self%atominfo(:)%cn_gen)
      deallocate(centers_A)
      if(thresh > 1. .or. thresh < 0.) THROW_HARD('Invalid input threshold! AR is in [0,1]; cluster_ar')
      ! Preparing for clustering
      ! Aspect ratios calculations
      ! call self%est_aspect_ratios(print_ar=.false.) ! it's needed to identify the dir of longest dim
      dim = size(self%atominfo)
      allocate(labels(dim), source = 0)
      ! classify
      call hac_1d(self%atominfo(:)%aspect_ratio,thresh,labels,centroids,populations)
      ! Stats calculations
      ncls = maxval(labels)
      allocate(stdev_within(ncls), source = 0.)
      ! stdev within the same class
      do i = 1, ncls
          do j = 1, dim
              if(labels(j) == i) stdev_within(i) = stdev_within(i) + (self%atominfo(j)%aspect_ratio-centroids(i))**2
          enddo
      enddo
      where (populations>1)
          stdev_within = sqrt(stdev_within/(real(populations)-1.))
      elsewhere
          stdev_within = 0.
      endwhere
      ! avg and stdev among different classes
      avg = 0.
      do i = 1, ncls
          avg = avg + centroids(i)
      enddo
      avg = avg/real(ncls)
      stdev = 0.
      if(ncls>1) then
        do i = 1, ncls
            stdev = stdev + (centroids(i) - avg)**2
        enddo
        stdev = sqrt(stdev/(real(ncls)-1.))
      endif
      ! median std cn on each class
      allocate(median_cn(ncls), source = 0.)
      do i = 1, ncls
        if(allocated(cn_cls)) deallocate(cn_cls)
        allocate(cn_cls(populations(i)), source = 0.)
        cnt = 0
        do j = 1, dim
            if(labels(j) == i) then
              cnt = cnt + 1
              cn_cls(cnt) = self%atominfo(j)%cn_std
            endif
        enddo
        median_cn(i) = median(cn_cls)
      enddo
      ! Output on a file
      str_thres = trim(real2str(thresh))
      call fopen(filnum, file='ClusterARThresh'//str_thres//'.txt', iostat=io_stat)
      write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
      write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
      do i = 1, dim
        write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3,a,f6.2)') 'Atom #: ', i, '; data (adimensional): ', self%atominfo(i)%aspect_ratio, '; class: ', labels(i)
      enddo
      write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
      do i = 1, ncls
        write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a,f6.2)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), '; stdev within the class: ', stdev_within(i), '; median std cn: ', median_cn(i)
      enddo
      write(unit = filnum,fmt ='(a,f6.2,a,f6.2)') 'AVG among the classes: ', avg, '; STDEV among the classes: ', stdev
      call fclose(filnum)
      ! Generate file that contains the calculater AR, that can be read afterwards
      call fopen(filnum, file='AspectRatio.csv', iostat=io_stat)
      write(filnum,*) 'ar'
      do i  = 1, self%n_cc
        write(filnum,*) self%atominfo(i)%aspect_ratio
      enddo
      call fclose(filnum)
      ! Generate one figure for each class
      if(GENERATE_FIGS) then
        call img_1clss%new_bimg(self%ldim, self%smpd)
        call self%img_cc%get_imat(imat_cc)
        allocate(imat_1clss(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
        do i = 1, ncls
            imat_1clss = 0
            do j = 1, dim
                if(labels(j) == i) then
                    where(imat_cc == j) imat_1clss = i
                endif
            enddo
            call img_1clss%set_imat(imat_1clss)
            call img_1clss%write_bimg('ArClass'//int2str(i)//'.mrc')
        enddo
        call img_1clss%kill_bimg
      endif
      if(allocated(stdev_within)) deallocate(stdev_within)
      deallocate(centroids, labels, populations)
    end subroutine cluster_ar

    ! Cluster the atoms wrt to the interatomic distances
    subroutine cluster_bondl(self, thresh)
        class(nanoparticle), intent(inout) :: self
        real,    intent(in)  :: thresh ! threshold for class definition, user inputted
        real,    allocatable :: centroids(:)
        integer, allocatable :: labels(:), populations(:)
        real,    allocatable :: stdev_within(:), avg_dist_cog(:)
        integer              :: i, j, ncls, dim, filnum, io_stat
        real                 :: avg, stdev, cog(3)
        type(binimage)       :: img_1clss
        integer, allocatable :: imat_cc(:,:,:), imat_1clss(:,:,:)
        character(len=4)     :: str_thres
        ! Preparing for clustering
        dim = size(self%atominfo)
        allocate(labels(dim), source = 0)
        ! classify
        call hac_1d(self%atominfo(:)%bondl, thresh, labels, centroids, populations)
        ! Stats calculations
        ncls = maxval(labels)
        allocate(stdev_within(ncls), source = 0.)
        allocate(avg_dist_cog(ncls), source = 0.)
        cog = self%masscen()
        ! stdev within the same class
        ! avg dist to the center of gravity of each class
        do i = 1, ncls
            do j = 1, dim
                if(labels(j) == i) then
                   stdev_within(i) = stdev_within(i) + (self%atominfo(j)%bondl - centroids(i))**2.
                   avg_dist_cog(i) = avg_dist_cog(i) + euclid(cog,self%atominfo(j)%center(:))
                endif
            enddo
        enddo
        avg_dist_cog = (avg_dist_cog*self%smpd)/real(populations) ! in A
        where (populations>1)
            stdev_within = sqrt(stdev_within/(real(populations)-1.))
        elsewhere
            stdev_within = 0.
        endwhere
        ! avg and stdev among different classes
        avg = 0.
        do i = 1, ncls
            avg = avg + centroids(i)
        enddo
        avg = avg/real(ncls)
        stdev = 0.
        if(ncls>1) then
          do i = 1, ncls
              stdev = stdev + (centroids(i) - avg)**2
          enddo
          stdev = sqrt(stdev/(real(ncls)-1.))
        endif
        ! Output on a file
        str_thres = trim(real2str(thresh))
        call fopen(filnum, file='ClusterInterDistThresh'//str_thres//'.txt', iostat=io_stat)
        write(unit = filnum,fmt ='(a,i2,a,f6.2)') 'NR OF IDENTIFIED CLUSTERS:', ncls, ' SELECTED THRESHOLD: ',  thresh
        write(unit = filnum,fmt ='(a)') 'CLASSIFICATION '
        do i = 1, dim
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a,i3)') 'Atom #: ', i, '; data (A): ', self%atominfo(i)%bondl, '; class: ', labels(i)
        enddo
        write(unit = filnum,fmt ='(a)') 'CLASS STATISTICS '
        do i = 1, ncls
          write(unit = filnum,fmt ='(a,i3,a,i3,a,f6.2,a,f6.2,a)') 'class: ', i, '; cardinality: ', populations(i), '; centroid: ', centroids(i), ' A; stdev within the class: ', stdev_within(i), ' A'
        enddo
        do i = 1, ncls
          write(unit = filnum,fmt ='(a,i3,a,f6.2,a)') 'class: ', i, '; average distance to the center of gravity: ', avg_dist_cog(i), ' A'
        enddo
        write(unit = filnum,fmt ='(a,f6.2,a,f6.2,a)') 'AVG among the classes: ', avg, ' A; STDEV among the classes: ', stdev, ' A'
        call fclose(filnum)
        call fopen(filnum, file='Dist.csv', iostat=io_stat)
        write(filnum,*) 'dist'
        do i  = 1, self%n_cc
          write(filnum,*) self%atominfo(i)%bondl
        enddo
        call fclose(filnum)
        ! Generate one figure for each class
        if(GENERATE_FIGS) then
          call img_1clss%new_bimg(self%ldim, self%smpd)
          call self%img_cc%get_imat(imat_cc)
          allocate(imat_1clss(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
          do i = 1, ncls
              imat_1clss = 0
              do j = 1, dim
                  if(labels(j) == i) then
                      where(imat_cc == j) imat_1clss = i
                  endif
              enddo
              call img_1clss%set_imat(imat_1clss)
              call img_1clss%write_bimg('DistClass'//int2str(i)//'.mrc')
          enddo
          call img_1clss%kill_bimg
        endif
        if(allocated(stdev_within)) deallocate(stdev_within)
        deallocate(centroids, labels, populations)
    end subroutine cluster_bondl

    ! others

    subroutine geometry_analysis(self, pdbfile2, thresh)
        use simple_math, only : plane_from_points
        class(nanoparticle),        intent(inout) :: self
        character(len=*),           intent(in)    :: pdbfile2  ! atomic pos of the 2 (or 3) selected atoms
        real,             optional, intent(in)    :: thresh    ! for belonging/not belonging to the plane/column
        type(atoms)           :: init_atoms, final_atoms
        type(binimage)        :: img_out
        character(len=3)      :: aux_var
        integer, allocatable  :: imat_cc(:,:,:), imat(:,:,:)
        real,    allocatable  :: line(:,:), plane(:,:,:),points(:,:),distances_totheplane(:),pointsTrans(:,:)
        real,    allocatable  :: radii(:),max_intensity(:),w(:),v(:,:),d(:),distances_totheline(:)
        integer :: i, n, t, s, filnum, io_stat, cnt_intersect, cnt
        logical :: flag(self%n_cc)
        real    :: atom1(3), atom2(3), atom3(3), dir_1(3), dir_2(3), vec(3), m(3), dist_plane, dist_line
        real    :: t_vec(N_DISCRET), s_vec(N_DISCRET), denominator, prod(3), centroid(3), tthresh
        if(present(thresh)) then
            tthresh = thresh
        else
            tthresh = 0.9*self%theoretical_radius
        endif
        call init_atoms%new(pdbfile2)
        n = init_atoms%get_n()
        if(n < 2 .or. n > 3 ) THROW_HARD('Inputted pdb file contains the wrong number of atoms!; geometry_analysis')
        do i = 1, N_DISCRET/2
            t_vec(i) = -real(i)/10.
        enddo
        t_vec(N_DISCRET/2+1:N_DISCRET) = -t_vec(1:N_DISCRET/2)
        s_vec(:) = t_vec(:)
        call self%img_cc%get_imat(imat_cc)
        call img_out%new(self%ldim, self%smpd)
        allocate(imat(self%ldim(1),self%ldim(2),self%ldim(3)), source = 0)
        flag(:) = .false.  ! initialization
        if( n == 2 )then
          write(logfhandle,*)'COLUMN IDENTIFICATION, INITIATION'
          allocate(line(3, N_DISCRET), source = 0.)
          atom1(:) = init_atoms%get_coord(1)/self%smpd + 1.
          atom2(:) = init_atoms%get_coord(2)/self%smpd + 1.
          dir_1 = atom1-atom2
          do t = 1, N_DISCRET
              line(1,t) = atom1(1) + t_vec(t)* dir_1(1)
              line(2,t) = atom1(2) + t_vec(t)* dir_1(2)
              line(3,t) = atom1(3) + t_vec(t)* dir_1(3)
          enddo
          ! calculate how many atoms does the line intersect and flag them
          do i = 1, self%n_cc
            do t = 1, N_DISCRET
              dist_line = euclid(self%atominfo(i)%center(:),line(:3,t))
              if(dist_line*self%smpd <= tthresh) then ! it intersects atoms i
                flag(i) = .true. !flags also itself
              endif
            enddo
          enddo
          ! generate volume for visualisation
          imat          = 0
          cnt_intersect = 0
          call img_out%new_bimg(self%ldim, self%smpd)
          call final_atoms%new(count(flag), dummy=.true.)
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt_intersect = cnt_intersect + 1
              call final_atoms%set_name(cnt_intersect,self%atom_name)
              call final_atoms%set_element(cnt_intersect,self%element)
              call final_atoms%set_coord(cnt_intersect,(self%atominfo(i)%center(:)-1.)*self%smpd)
              where(imat_cc == i) imat = 1
            endif
          enddo
          call img_out%set_imat(imat)
          call img_out%write_bimg('ImageColumn.mrc')
          call final_atoms%writePDB('AtomColumn')
          call final_atoms%kill
          ! Find the line that best fits the atoms
          allocate(points(3,count(flag)), source = 0.)
          cnt = 0
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              points(:3,cnt) = self%atominfo(i)%center(:)
            endif
          enddo
          ! svd fit (https://au.mathworks.com/matlabcentral/answers/424591-3d-best-fit-line)
          allocate(pointsTrans(count(flag),3), source = 0.) ! because svdcmp modifies its input
          ! translate
          centroid = sum(points, dim = 2)/real(count(flag)) ! it belongs to the line
          do i = 1, count(flag)
            pointsTrans(i,:3) = points(:3,i) - centroid(:3)
          enddo
          allocate(w(3), v(3,3), source = 0.)
          allocate(d(3), source = 0.)
          call svdcmp(pointsTrans,w,v)
          d = v(:,1)
          write(logfhandle,*) 'Directional vector of the line', d
          ! line
          ! line(1,t) = centroid(1) + t_vec(t)* d(1)
          ! line(2,t) = centroid(2) + t_vec(t)* d(2)
          ! line(3,t) = centroid(3) + t_vec(t)* d(3)
          ! calculate the distance to the points from the identified line
          allocate(distances_totheline(cnt), source = 0.)
          allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
          denominator = sqrt(d(1)**2+d(2)**2+d(3)**2)
          m = self%masscen()
          cnt = count(flag)
          do i = 1, cnt
            vec  = centroid(:3)-points(:3,i)
            prod = cross(vec,d)
            distances_totheline(i) = sqrt(prod(1)**2+prod(2)**2+prod(3)**2)/denominator
            radii(i) = euclid(points(:,i), m)*self%smpd
          enddo
          distances_totheline = distances_totheline*self%smpd ! go to A
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
          ! Read ratios, ang_var, dists, max_intensity
          if(.not. allocated(self%atominfo)) allocate(self%atominfo(self%n_cc))
          call fopen(filnum, file='../'//'../'//'AspectRatio.csv', action='readwrite',iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file AspectRatio.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*) ! first line is variable name
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%aspect_ratio
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Ang.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Ang.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%polar_angle
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Dist.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Dist.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%bondl
          enddo
          call fclose(filnum)
          allocate(max_intensity(self%n_cc))
          call fopen(filnum, file='../'//'../'//'MaxIntensity.csv', action='readwrite', iostat=io_stat)
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) max_intensity(i)
          enddo
          call fclose(filnum)
          ! Output on Excel file all the stats on the atoms belonging to the plane
          cnt = 0
          call fopen(filnum, file='AtomColumnInfo.txt',iostat=io_stat)
          write (filnum,*) '        Atom #    ','   Ar        ','       Dist     ','         Ang     ','         MaxInt     '
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              write (filnum,*) i, '   ', self%atominfo(i)%aspect_ratio, self%atominfo(i)%bondl, self%atominfo(i)%polar_angle, max_intensity(i)
            endif
          end do
          call fclose(filnum)
        elseif(n == 3) then
          write(logfhandle,*)'PLANE IDENTIFICATION, INITIATION'
          atom1(:) = init_atoms%get_coord(1)/self%smpd + 1.
          atom2(:) = init_atoms%get_coord(2)/self%smpd + 1.
          atom3(:) = init_atoms%get_coord(3)/self%smpd + 1.
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
          do i = 1, self%n_cc
            do t = 1, N_DISCRET
              do s = 1, N_DISCRET
                dist_plane = euclid(self%atominfo(i)%center(:),plane(:3,t,s))
                if(dist_plane*self%smpd <= tthresh) then ! it intersects atoms i
                  flag(i) = .true. !flags also itself
                endif
              enddo
            enddo
          enddo
          ! generate volume for visualisation
          ! reset
          imat    = 0
          cnt_intersect = 0
          call img_out%new_bimg(self%ldim, self%smpd)
          call final_atoms%new(count(flag), dummy=.true.)
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt_intersect = cnt_intersect + 1
              call final_atoms%set_name(cnt_intersect,self%atom_name)
              call final_atoms%set_element(cnt_intersect,self%element)
              call final_atoms%set_coord(cnt_intersect,(self%atominfo(i)%center(:)-1.)*self%smpd)
              where(imat_cc == i) imat = 1
            endif
          enddo
          call img_out%set_imat(imat)
          call img_out%write_bimg('ImagePlane.mrc')
          call final_atoms%writePDB('AtomPlane')
          call final_atoms%kill
          allocate(points(3, count(flag)), source = 0.)
          m = self%masscen()
          cnt = 0
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              points(:3,cnt) = self%atominfo(i)%center(:)-m(:)
            endif
          enddo
          vec = plane_from_points(points)
          allocate(distances_totheplane(cnt), source = 0.)
          allocate(radii(cnt), source = 0.) ! which radius is the atom center belonging to
          cnt = 0
          denominator = sqrt(vec(1)**2+vec(2)**2+1.)
          write(logfhandle,*) 'Normal vector: [', vec(1), ',', vec(2), ',', -1., ']'
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              ! formula for distance of a point to a plane
              distances_totheplane(cnt) = abs(vec(1)*points(1,cnt)+vec(2)*points(2,cnt)-points(3,cnt)+vec(3))/denominator
              radii(cnt) = euclid(self%atominfo(i)%center(:), m)*self%smpd
            endif
          enddo
          distances_totheplane = (distances_totheplane)*self%smpd
          call fopen(filnum, file='Radii.csv', action='readwrite', iostat=io_stat)
          write (filnum,*) 'r'
          do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(radii(i)))
          end do
          call fclose(filnum)
          call fopen(filnum, file='DistancesToThePlane.csv', action='readwrite',iostat=io_stat)
          write (filnum,*) 'd'
          do i = 1, cnt
            write (filnum,'(A)', advance='yes') trim(real2str(distances_totheplane(i)))
          end do
          call fclose(filnum)
          ! Read ratios, ang_var, dists, max_intensity
          if(.not. allocated(self%atominfo)) allocate(self%atominfo(self%n_cc))
          call fopen(filnum, file='../'//'../'//'AspectRatio.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file AspectRatio.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)  ! first line is variable name
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%aspect_ratio
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Ang.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Ang.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%polar_angle
          enddo
          call fclose(filnum)
          call fopen(filnum, file='../'//'../'//'Dist.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file Dist.csv Did you run cluster_analysis?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) self%atominfo(i)%bondl
          enddo
          call fclose(filnum)
          allocate(max_intensity(self%n_cc))
          call fopen(filnum, file='../'//'../'//'MaxIntensity.csv', action='readwrite', iostat=io_stat)
          if( io_stat .ne. 0 ) then
            THROW_HARD('Unable to read file MaxIntensity.csv Did you run atoms_stats?; geometry_analysis')
          endif
          read(filnum,*)
          do i = 1, self%n_cc
            read(filnum,*) max_intensity(i)
          enddo
          call fclose(filnum)
          ! Output on txt file all the stats on the atoms belonging to the plane
          cnt = 0
          call fopen(filnum, file='AtomPlaneInfo.txt',action='write',iostat=io_stat)
          write (filnum,*) '        Atom #    ','   Ar        ','       Dist     ','         Ang     ','      MaxInt     '
          do i = 1, self%n_cc
            if(flag(i)) then
              cnt = cnt + 1
              write (filnum,*) i, '   ', self%atominfo(i)%aspect_ratio, self%atominfo(i)%bondl, self%atominfo(i)%polar_angle, max_intensity(i)
            endif
          end do
          call fclose(filnum)
        endif
        call img_out%kill_bimg
        call init_atoms%kill
        if(allocated(line))  deallocate(line)
        if(allocated(plane)) deallocate(plane)
        deallocate(imat_cc, imat)

    contains

      ! Compute the cross product of 2 3D real vectors
      function cross(a, b) result(c)
          real, intent(in) :: a(3),b(3)
          real :: c(3)
          c(1) = a(2) * b(3) - a(3) * b(2)
          c(2) = a(3) * b(1) - a(1) * b(3)
          c(3) = a(1) * b(2) - a(2) * b(1)
      end function cross

    end subroutine geometry_analysis

    subroutine kill_nanoparticle(self)
        class(nanoparticle), intent(inout) :: self
        call self%img%kill()
        call self%img_raw%kill
        call self%img_bin%kill_bimg()
        call self%img_cc%kill_bimg()
        call self%centers_pdb%kill
        if( allocated(self%atominfo) ) deallocate(self%atominfo)
    end subroutine kill_nanoparticle

end module simple_nanoparticle

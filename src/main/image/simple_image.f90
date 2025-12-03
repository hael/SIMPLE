! the abstract image data type and its methods. 2D/3D & FT/real all implemented by this class
! and Fourier transformations done in-place to reduce memory usage
module simple_image
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ftiter,  only: ftiter
use simple_imgfile, only: imgfile
use simple_winfuns, only: winfuns
use simple_neighs
use gnufor2
implicit none

public :: image, image_ptr, test_image, image_stack
private
#include "simple_local_flags.inc"

type image_ptr
    real(kind=c_float),            public, pointer :: rmat(:,:,:)
    complex(kind=c_float_complex), public, pointer :: cmat(:,:,:)
end type image_ptr

type :: image
    private
    logical                                :: ft=.false.           !< Fourier transformed or not
    integer                                :: ldim(3)=[1,1,1]      !< logical image dimensions
    integer                                :: nc                   !< number of F-comps
    real                                   :: smpd                 !< sampling distance
    type(ftiter)                           :: fit                  !< Fourier iterator object
    type(c_ptr)                            :: p                    !< c pointer for fftw allocation
    real(kind=c_float),            pointer :: rmat(:,:,:)=>null()  !< image pixels/voxels (in data)
    complex(kind=c_float_complex), pointer :: cmat(:,:,:)=>null()  !< Fourier components
    real                                   :: shconst(3)           !< shift constant
    type(c_ptr)                            :: plan_fwd             !< fftw plan for the image (fwd)
    type(c_ptr)                            :: plan_bwd             !< fftw plan for the image (bwd)
    integer                                :: array_shape(3)       !< shape of complex array
    logical                                :: wthreads  = .false.  !< with threads flag
    logical                                :: existence = .false.  !< indicates existence
contains
    ! CONSTRUCTORS / LIFECYCLE, file: simple_image_core.f90
    procedure          :: new
    procedure          :: set_wthreads
    procedure          :: construct_thread_safe_tmp_imgs
    procedure          :: kill
    procedure          :: kill_thread_safe_tmp_imgs
    procedure          :: copy
    procedure          :: copy_fast
    procedure, private :: disc_1, disc_2
    generic            :: disc => disc_1, disc_2
    procedure          :: ring, soft_ring
    procedure          :: disc_sideview
    procedure          :: cylinder
    ! ARITHMETICS, file, file: simple_image_arith.f90
    procedure, private :: assign
    procedure, private :: assign_r2img
    procedure, private :: assign_c2img
    generic :: assignment(=) => assign, assign_r2img, assign_c2img
    procedure, private :: subtraction
    generic :: operator(-) => subtraction
    procedure, private :: addition
    procedure, private :: addition_const_real
    generic :: operator(+) => addition, addition_const_real
    procedure, private :: multiplication
    procedure, private :: multiplication_const_real
    procedure, private :: multiplication_const_int
    generic :: operator(*) => multiplication, multiplication_const_real, multiplication_const_int
    procedure, private :: division
    generic :: operator(/) => division
    procedure, private :: add_1, add_2, add_3, add_4, add_5
    generic            :: add => add_1, add_2, add_3, add_4, add_5
    procedure          :: add_workshare
    procedure, private :: subtr_1, subtr_2, subtr_3, subtr_4
    generic            :: subtr => subtr_1, subtr_2, subtr_3, subtr_4
    procedure, private :: div_1, div_2, div_3, div_4
    generic            :: div => div_1, div_2, div_3, div_4
    procedure, private :: mul_1, mul_2, mul_3, mul_4, mul_5
    generic            :: mul => mul_1, mul_2, mul_3, mul_4, mul_5
    procedure, private :: conjugate
    generic            :: conjg => conjugate
    procedure, private :: mul_rmat_at_1, mul_rmat_at_2
    generic            :: mul_rmat_at => mul_rmat_at_1, mul_rmat_at_2
    procedure, private :: div_rmat_at_1, div_rmat_at_2
    generic            :: div_rmat_at => div_rmat_at_1, div_rmat_at_2
    procedure, private :: add_cmat_at_1, add_cmat_at_2
    generic            :: add_cmat_at => add_cmat_at_1, add_cmat_at_2
    procedure, private :: mul_cmat_at_1, mul_cmat_at_2, mul_cmat_at_3, mul_cmat_at_4
    generic            :: mul_cmat_at => mul_cmat_at_1, mul_cmat_at_2, mul_cmat_at_3, mul_cmat_at_4
    procedure, private :: mul_cmat_1, mul_cmat_2
    generic            :: mul_cmat => mul_cmat_1, mul_cmat_2
    procedure, private :: div_cmat_at_1, div_cmat_at_2, div_cmat_at_3, div_cmat_at_4
    generic            :: div_cmat_at => div_cmat_at_1, div_cmat_at_2, div_cmat_at_3, div_cmat_at_4
    procedure          :: sq_rt
    procedure          :: add_cmats_to_cmats
    ! FFT, file: simple_image_fft.f90
    procedure          :: fft  => fwd_ft
    procedure          :: ifft => bwd_ft
    procedure          :: fft_noshift
    procedure          :: img2ft
    procedure          :: ft2img
    procedure          :: pad_fft
    procedure          :: norm_noise_pad_fft
    procedure          :: expand_ft
    ! I/O
    procedure, private :: open
    procedure          :: read
    procedure          :: write
    procedure          :: update_header_stats
    procedure          :: write_jpg
    ! ACCESS, file: simple_image_access.f90
    ! Basic shape / metadata
    procedure          :: get_array_shape
    procedure          :: get_ldim, get_box
    procedure          :: get_smpd
    procedure          :: get_nyq
    procedure          :: get_filtsz
    procedure          :: get_shconst
    procedure          :: set_smpd
    ! Getters
    procedure          :: get
    procedure          :: get_rmat
    procedure          :: get_mat_ptrs
    procedure          :: get_rmat_ptr
    procedure          :: get_rmat_sub
    procedure, private :: get_rmat_at_1, get_rmat_at_2
    generic            :: get_rmat_at => get_rmat_at_1, get_rmat_at_2
    procedure          :: get_cmat
    procedure          :: get_cmat_ptr
    procedure          :: get_cmat_sub
    procedure, private :: get_cmat_at_1, get_cmat_at_2
    generic            :: get_cmat_at => get_cmat_at_1, get_cmat_at_2
    procedure          :: get_fcomp
    procedure          :: get_fcomp2D
    ! Setters
    procedure, private :: set_1, set_2
    generic            :: set => set_1, set_2
    procedure          :: set_rmat
    procedure          :: set_rmat_at
    procedure, private :: set_cmat_1, set_cmat_2, set_cmat_3, set_cmat_4
    generic            :: set_cmat => set_cmat_1, set_cmat_2, set_cmat_3, set_cmat_4
    procedure, private :: set_cmat_at_1, set_cmat_at_2
    generic            :: set_cmat_at => set_cmat_at_1, set_cmat_at_2
    procedure          :: set_fcomp
    procedure          :: set_within
    procedure          :: set_cmats_from_cmats
    ! Slices, sub-images, freq info
    procedure          :: get_slice
    procedure          :: set_slice
    procedure          :: get_subimg
    procedure          :: get_lfny
    procedure          :: get_lhp
    procedure          :: get_lp
    procedure          :: get_spat_freq
    procedure          :: get_find
    ! Flags
    procedure          :: rmat_associated
    procedure          :: cmat_associated
    procedure          :: is_wthreads
    procedure          :: set_ft
    ! Serialization / misc
    procedure, private :: serialize_1, serialize_2, serialize_3
    generic            :: serialize => serialize_1, serialize_2, serialize_3
    procedure          :: unserialize
    procedure          :: winserialize
    ! CHECKS, file: simple_image_checks.f90
    procedure          :: exists
    procedure          :: is_2d
    procedure          :: is_3d
    procedure          :: even_dims
    procedure          :: square_dims
    procedure, private :: same_dims_1
    generic            :: operator(.eqdims.) => same_dims_1
    procedure          :: same_dims
    procedure          :: same_smpd
    generic            :: operator(.eqsmpd.) => same_smpd
    procedure          :: is_ft
    procedure          :: is_empty
    ! FILTERS/DENOISE, file: simple_image_filt.f90
    procedure          :: bp, lp
    procedure          :: lp_background
    procedure          :: bpgau2D, bpgau3D
    procedure          :: tophat
    procedure          :: phase_rand
    procedure          :: ran_phases_below_noise_power
    procedure          :: whiten_noise_power
    procedure          :: real_space_filter
    procedure          :: hannw
    procedure          :: apply_bfac
    procedure, private :: apply_filter_1, apply_filter_2
    generic            :: apply_filter => apply_filter_1, apply_filter_2
    procedure          :: apply_filter_serial
    procedure          :: NLmean2D, NLmean2D_eo, NLmean3D, NLmean3D_eo
    procedure          :: ICM2D, ICM2D_eo, ICM3D, ICM3D_eo
    procedure          :: GLCM
    ! FREQUENCY ANALYSIS, file: simple_image_freq_anal.f90
    procedure          :: acf
    procedure          :: ccf
    procedure          :: spectrum
    procedure          :: power_spectrum
    procedure          :: guinier_bfac
    procedure          :: guinier
    procedure          :: fsc, fsc_scaled
    procedure          :: get_res
    procedure          :: frc_pspec
    procedure          :: fcomps_below_noise_power_stats
    procedure          :: resmsk
    procedure          :: img2spec
    procedure          :: mic2spec
    procedure          :: pspec_graphene_mask
    procedure          :: dampen_pspec_central_cross
    ! CALCULATORS, file: simple_image_calc.f90
    ! Basic stats / local stats
    procedure          :: get_sum_int
    procedure          :: minmax
    procedure          :: loc_sdev
    procedure          :: avg_loc_sdev
    procedure          :: loc_var, loc_var3D
    procedure          :: rmsd
    procedure, private :: stats_1, stats_2
    generic            :: stats => stats_1, stats_2
    procedure          :: variance
    procedure          :: skew, kurt
    procedure          :: noisesdev
    procedure          :: mean
    procedure          :: contains_nans
    procedure          :: checkimg4nans
    procedure          :: cure
    procedure          :: dead_hot_positions
    ! Gradients / geometry
    procedure          :: calc_gradient
    procedure          :: gradients_magnitude
    procedure          :: gradient
    procedure          :: calc_ice_score
    procedure          :: calc_principal_axes_rotmat
    ! Physical coords helpers
    procedure          :: loop_lims
    procedure, private :: comp_addr_phys1, comp_addr_phys2, comp_addr_phys3
    generic            :: comp_addr_phys =>  comp_addr_phys1, comp_addr_phys2, comp_addr_phys3
    !  Correlation / distances
    procedure          :: corr
    procedure          :: corr_shifted
    procedure, private :: real_corr_1, real_corr_2
    generic            :: real_corr => real_corr_1, real_corr_2
    procedure          :: euclid_dist_two_imgs
    procedure          :: phase_corr
    procedure          :: fcorr_shift, fcorr_shift3D
    procedure          :: prenorm4real_corr_1, prenorm4real_corr_2, prenorm4real_corr_3
    generic            :: prenorm4real_corr => prenorm4real_corr_1, prenorm4real_corr_2, prenorm4real_corr_3
    procedure, private :: real_corr_prenorm_1, real_corr_prenorm_2, real_corr_prenorm_3
    generic            :: real_corr_prenorm => real_corr_prenorm_1, real_corr_prenorm_2, real_corr_prenorm_3
    procedure          :: radial_cc
    procedure          :: sqeuclid
    procedure, private :: sqeuclid_matrix_1, sqeuclid_matrix_2
    generic            :: sqeuclid_matrix => sqeuclid_matrix_1, sqeuclid_matrix_2
    procedure          :: euclid_norm
    ! cost / shift
    procedure          :: opt_filter_costfun
    procedure          :: opt_filter_costfun_workshare
    procedure, private :: oshift_1, oshift_2
    generic            :: oshift => oshift_1, oshift_2
    procedure, private :: gen_argtransf_comp
    ! VISUALIZATION, file: simple_image_vis.f90
    ! print/vis helpers
    procedure          :: print_cmat
    procedure          :: print_rmat
    procedure          :: vis
    procedure          :: before_after
    procedure          :: scale_pspec4viz
    procedure          :: generate_orthogonal_reprojs
    procedure          :: collage
    procedure          :: tile
    ! SEGMENTATION, file: simple_image_seg.f90
    procedure          :: nforeground
    procedure          :: nbackground
    procedure, private :: binarize_1, binarize_2, binarize_3
    generic            :: binarize => binarize_1, binarize_2, binarize_3
    procedure          :: cendist
    procedure          :: bin_inv
    procedure          :: remove_edge
    procedure          :: zero2one
    procedure          :: one_at_edge
    procedure          :: bin2logical
    procedure          :: logical2bin
    procedure          :: density_inoutside
    procedure          :: calc_bin_thres
    procedure          :: mask
    procedure          :: taper_edges, taper_edges_hann
    ! OPERATIONS, file: simple_image_ops.f90
    ! ctf
    procedure          :: ctf_dens_correct
    procedure          :: ctf_dens_correct_wiener
    ! insertions
    procedure          :: insert
    procedure          :: insert_lowres
    procedure          :: insert_lowres_serial
    ! noise
    procedure          :: ran
    procedure          :: gauran
    procedure          :: add_gauran
    procedure          :: salt_n_pepper
    ! background
    procedure          :: div_w_instrfun
    procedure          :: estimate_background
    procedure          :: subtr_backgr
    procedure          :: subtr_backgr_ramp
    procedure          :: subtract_background
    procedure          :: upsample_square_background
    procedure          :: remove_neg
    ! arithmetics
    procedure          :: neg
    procedure          :: div_below
    procedure          :: inv
    ! zeroing
    procedure          :: zero
    procedure          :: zero_and_flag_ft
    procedure          :: zero_and_unflag_ft
    procedure          :: zero_background
    procedure          :: zero_below
    procedure          :: zero_edgeavg
    procedure          :: zero_env_background
    procedure          :: zero_neg
    ! GEOMETRICAL, file: simple_image_geom.f90
    ! windowing
    procedure          :: window
    procedure          :: window_slim
    procedure          :: window_center
    procedure          :: add_window
    procedure          :: win2arr
    procedure          :: win2arr_rad
    ! shapes
    procedure          :: corner
    procedure          :: ellipse
    procedure          :: corners
    procedure, private :: gauimg_1, gauimg_2
    generic            :: gauimg => gauimg_1, gauimg_2
    procedure          :: gauimg2D
    procedure          :: gauimg3D
    procedure          :: reshape2cube
    procedure          :: square
    ! pad/clip/crop
    procedure          :: pad
    procedure          :: pad_inplace
    procedure, private :: pad_mirr_1, pad_mirr_2
    generic            :: pad_mirr => pad_mirr_1, pad_mirr_2
    procedure          :: clip
    procedure          :: clip_inplace
    procedure          :: read_and_crop
    ! flip/mirror/rotate/shift
    procedure          :: flip
    procedure          :: mirror
    procedure          :: calc_shiftcen, calc_shiftcen_serial
    procedure          :: roavg
    procedure          :: rtsq
    procedure          :: rtsq_serial
    procedure          :: shift_phorig
    procedure          :: shift
    procedure, private :: shift2Dserial_1, shift2Dserial_2
    generic            :: shift2Dserial => shift2Dserial_1, shift2Dserial_2
    procedure          :: masscen
    procedure          :: masscen_adjusted
    ! NORMALIZE, file: simple_image_norm.f90
    procedure          :: scale_pixels
    procedure          :: norm
    procedure          :: norm_minmax
    procedure          :: norm4viz
    procedure          :: norm_ext
    procedure          :: norm_noise
    procedure          :: norm_within
    procedure          :: cure_outliers
    procedure          :: quantize_fwd
    procedure          :: quantize_bwd
end type image

interface

    ! ===== Constructors/lifecycle =====

    module subroutine new( self, ldim, smpd, wthreads )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: wthreads
    end subroutine new

    module subroutine set_wthreads( self, wthreads )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: wthreads
    end subroutine set_wthreads

    module subroutine construct_thread_safe_tmp_imgs( self, nthr )
        class(image), intent(in) :: self
        integer,      intent(in) :: nthr
    end subroutine construct_thread_safe_tmp_imgs

    module subroutine kill( self )
        class(image), intent(inout) :: self
    end subroutine kill

    module subroutine kill_thread_safe_tmp_imgs( self )
        class(image), intent(in) :: self
    end subroutine kill_thread_safe_tmp_imgs

    module subroutine copy( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
    end subroutine copy

    module subroutine copy_fast( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
    end subroutine copy_fast

    module subroutine disc_1( self, ldim, smpd, radius, npix )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, radius
        integer, optional, intent(inout) :: npix
    end subroutine disc_1

    module subroutine disc_2( self, ldim, smpd, radius, lmsk, npix )
        class(image),         intent(inout) :: self
        integer,              intent(in)    :: ldim(3)
        real,                 intent(in)    :: smpd, radius
        logical, allocatable, intent(out)   :: lmsk(:,:,:)
        integer, optional,    intent(inout) :: npix
    end subroutine disc_2

    module subroutine ring( self, ldim, smpd, outer_radius, inner_radius, npix )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, outer_radius, inner_radius
        integer, optional, intent(inout) :: npix
    end subroutine ring

    module subroutine soft_ring( self, ldim, smpd, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        real,              intent(in)    :: radius ! in pixels
    end subroutine soft_ring

    module subroutine disc_sideview( self, ldim, smpd, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        real,              intent(in)    :: radius ! in pixels
    end subroutine disc_sideview

    module subroutine cylinder( self, ldim, smpd, height, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, height, radius
    end subroutine cylinder

    ! ===== Arithmetic procedure interfaces =====

    !--- assignment(=) ---!
    module subroutine assign( selfout, selfin )
        class(image), intent(inout) :: selfout
        class(image), intent(in)    :: selfin
    end subroutine assign

    module subroutine assign_r2img( self, realin )
        class(image), intent(inout) :: self
        real,         intent(in)    :: realin
    end subroutine assign_r2img

    module subroutine assign_c2img( self, compin )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: compin
    end subroutine assign_c2img

    !--- operator(-) ---!
    module function subtraction( self_from, self_to ) result( self )
        class(image), intent(in) :: self_from, self_to
        type(image) :: self
    end function subtraction

    !--- operator(+) ---!
    module function addition( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
    end function addition

    module function addition_const_real( self1, rconst ) result( self )
        class(image), intent(in) :: self1
        real,         intent(in) :: rconst
        type(image) :: self
    end function addition_const_real

    !--- operator(*) ---!
    module function multiplication( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
    end function multiplication

    module function multiplication_const_real( self1, rconst ) result( self )
        class(image), intent(in) :: self1
        real,         intent(in) :: rconst
        type(image) :: self
    end function multiplication_const_real

    module function multiplication_const_int( self1, iconst ) result( self )
        class(image), intent(in) :: self1
        integer,      intent(in) :: iconst
        type(image) :: self
    end function multiplication_const_int

    !--- operator(/) ---!
    module function division( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
    end function division

    !--- add_* family (adapt args to your real code) ---!
    module subroutine add_1( self, self_to_add, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_add
        real, optional, intent(in)    :: w
    end subroutine add_1

    module subroutine add_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(in)    :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
    end subroutine add_2

    module subroutine add_3( self, rcomp, i, j, k )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rcomp
        integer,      intent(in)    :: i, j, k
    end subroutine add_3

    module subroutine add_4( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
    end subroutine add_4

    module subroutine add_5( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
    end subroutine add_5

    module subroutine add_workshare( self, self_to_add )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_add
    end subroutine add_workshare

    !--- subtr_* ---!
    module subroutine subtr_1( self, self_to_subtr, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_subtr
        real, optional, intent(in)    :: w
    end subroutine subtr_1

    module subroutine subtr_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(out)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
    end subroutine subtr_2

    module subroutine subtr_3( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
    end subroutine subtr_3

    module subroutine subtr_4( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
    end subroutine subtr_4

    !--- div_* ---!
    module subroutine div_1( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
    end subroutine div_1

    module subroutine div_2( self, logi, k, square )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: k(:,:,:)
        logical,      intent(in)    :: square
    end subroutine div_2

    module subroutine div_3( self, logi, k, phys_in )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        real,              intent(in)    :: k
        integer, optional, intent(in)    :: phys_in(3)
    end subroutine div_3

    module subroutine div_4( self, self2div )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2div
    end subroutine div_4

    !--- mul_* ---!
    module subroutine mul_1( self, logi, rc, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        real,              intent(in)    :: rc
        integer, optional, intent(in)    :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
    end subroutine mul_1

    module subroutine mul_2( self, rc )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rc
    end subroutine mul_2

    module subroutine mul_3( self, self2mul )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
    end subroutine mul_3

    module subroutine mul_4( self, self2mul, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
        real,         intent(in)    :: lp
    end subroutine mul_4

    module subroutine mul_5( self, logi, c)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: c
    end subroutine mul_5

    module function conjugate( self ) result ( self_out )
        class(image), intent(in) :: self
        type(image) :: self_out
    end function conjugate

    !--- index-based per-element ops ---!
    module subroutine mul_rmat_at_1( self, logi, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: rval
    end subroutine mul_rmat_at_1

    module elemental pure subroutine mul_rmat_at_2( self,i, j, k, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: rval
    end subroutine mul_rmat_at_2

    module subroutine div_rmat_at_1( self, logi, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: rval
    end subroutine div_rmat_at_1

    module subroutine div_rmat_at_2( self, i, j, k, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: rval
    end subroutine div_rmat_at_2

    module pure subroutine add_cmat_at_1( self , phys , comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: comp
    end subroutine add_cmat_at_1

    module pure subroutine add_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        complex,      intent(in) :: comp
    end subroutine add_cmat_at_2

    module pure subroutine mul_cmat_at_1( self, phys, rval)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        real,         intent(in)    :: rval
    end subroutine mul_cmat_at_1

    module pure subroutine mul_cmat_at_2( self, phys, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: cval
    end subroutine mul_cmat_at_2

    module pure subroutine mul_cmat_at_3( self, h,k,l,rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        real,         intent(in)    :: rval
    end subroutine mul_cmat_at_3

    module pure subroutine mul_cmat_at_4( self, h,k,l, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
    end subroutine mul_cmat_at_4

    module subroutine mul_cmat_1( self, rmat )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
    end subroutine mul_cmat_1

    module subroutine mul_cmat_2( self, rmat, resmsk )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        logical,      intent(in)    :: resmsk(self%array_shape(1),self%array_shape(2),self%array_shape(3))
    end subroutine mul_cmat_2

    module pure subroutine div_cmat_at_1( self, phys, rval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        real,              intent(in)    :: rval
    end subroutine div_cmat_at_1

    module pure subroutine div_cmat_at_2( self,h,k,l, rval)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        real,         intent(in) :: rval
    end subroutine div_cmat_at_2

    module pure subroutine div_cmat_at_3( self, phys, cval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        complex,           intent(in)    :: cval
    end subroutine div_cmat_at_3

    module pure subroutine div_cmat_at_4( self,h,k,l, cval)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
    end subroutine div_cmat_at_4

    module subroutine sq_rt( self )
        class(image), intent(inout) :: self
    end subroutine sq_rt

    module subroutine add_cmats_to_cmats( self1 , self2 , self3, self4, self2set1, self2set2, lims, expcmat3, expcmat4)
        class(image), intent(in)    :: self1, self2,self3,self4
        class(image), intent(inout) :: self2set1, self2set2
        integer,      intent(in)    :: lims(3,2)
        real,         intent(inout) :: expcmat3(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,         intent(inout) :: expcmat4(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
    end subroutine add_cmats_to_cmats

    ! ===== FFT procedure interfaces =====
    
    module subroutine fwd_ft(self)
        class(image), intent(inout) :: self
    end subroutine fwd_ft

    module subroutine bwd_ft( self )
        class(image), intent(inout) :: self
    end subroutine bwd_ft

    module subroutine fft_noshift( self )
        class(image), intent(inout) :: self
    end subroutine fft_noshift

    module subroutine img2ft( self, img )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: img
    end subroutine img2ft

    module subroutine ft2img( self, which, img )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: which
        class(image),     intent(inout) :: img
    end subroutine ft2img

    module subroutine pad_fft( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: self_out
    end subroutine pad_fft

    module subroutine norm_noise_pad_fft( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
    end subroutine norm_noise_pad_fft

    module function expand_ft( self ) result( fplane )
        class(image), intent(in) :: self
        complex, allocatable :: fplane(:,:)
    end function expand_ft

    ! ===== I/O procedure interfaces =====

    module subroutine open( self, fname, ioimg, formatchar, readhead, rwaction )
        class(image),               intent(inout) :: self
        class(string),              intent(in)    :: fname
        class(imgfile),             intent(inout) :: ioimg
        character(len=1), optional, intent(in)    :: formatchar
        logical,          optional, intent(in)    :: readhead
        character(len=*), optional, intent(in)    :: rwaction
    end subroutine open

    module subroutine read( self, fname, i, readhead )
        class(image),               intent(inout) :: self
        class(string),              intent(in)    :: fname
        integer,          optional, intent(in)    :: i
        logical,          optional, intent(in)    :: readhead
    end subroutine read

    module subroutine write( self, fname, i, del_if_exists)
        class(image),      intent(inout) :: self
        class(string),     intent(in)    :: fname
        integer, optional, intent(in)    :: i
        logical, optional, intent(in)    :: del_if_exists
    end subroutine write

    module subroutine update_header_stats( self, fname, stats)
        class(image),  intent(inout) :: self
        class(string), intent(in)    :: fname
        real,          intent(in)    :: stats(4) ! stats to update: min, max, mean, rms
    end subroutine update_header_stats

    module subroutine write_jpg( self, fname, quality, colorspec, norm )
        use simple_jpg, only: jpg_img
        class(image),               intent(inout) :: self
        class(string),              intent(in)    :: fname
        integer,          optional, intent(in)    :: quality, colorspec
        logical,          optional, intent(in)    :: norm
    end subroutine write_jpg

    ! ===== access procedure interfaces =====

    !--- Basic shape / metadata ---!

    module pure function get_array_shape( self ) result( shape)
        class(image), intent(in) :: self
        integer :: shape(3)
    end function get_array_shape

    module pure function get_ldim( self ) result( ldim )
        class(image), intent(in) :: self
        integer :: ldim(3)
    end function get_ldim

    module pure integer function get_box( self )
        class(image), intent(in) :: self
    end function get_box

    module pure function get_smpd( self ) result( smpd )
        class(image), intent(in) :: self
        real :: smpd
    end function get_smpd

    module pure function get_nyq( self ) result( nyq )
        class(image), intent(in) :: self
        integer :: nyq
    end function get_nyq

    module pure function get_filtsz( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
    end function get_filtsz

    module pure function get_shconst( self ) result( shconst )
        class(image), intent(in) :: self
        real :: shconst(3)
    end function get_shconst

    module subroutine set_smpd( self, smpd )
        class(image), intent(inout) :: self
        real,         intent(in)    :: smpd
    end subroutine set_smpd

    !--- Getters ---!

    module function get( self, logi ) result( val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real :: val
    end function get

    module pure function get_rmat( self ) result( rmat )
        class(image), intent(in) :: self
        real, allocatable :: rmat(:,:,:)
    end function get_rmat

    module subroutine get_mat_ptrs( self, mat_ptrs )
        class(image), target, intent(in)  :: self
        class(image_ptr),     intent(out) :: mat_ptrs
    end subroutine get_mat_ptrs

    module subroutine get_rmat_ptr( self, rmat_ptr )
        class(image), target,        intent(in)  :: self
        real(kind=c_float), pointer, intent(out) :: rmat_ptr(:,:,:)
    end subroutine get_rmat_ptr

    module pure subroutine get_rmat_sub( self, rmat )
        class(image), intent(in)  :: self
        real,         intent(out) :: rmat(self%ldim(1),self%ldim(2),self%ldim(3))
    end subroutine get_rmat_sub

    module pure function get_rmat_at_1( self, logi ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3)
        real :: val
    end function get_rmat_at_1

    module pure function get_rmat_at_2( self, i,j,k ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: i,j,k
        real :: val
    end function get_rmat_at_2

    module pure function get_cmat( self ) result( cmat )
        class(image), intent(in) :: self
        complex, allocatable :: cmat(:,:,:)
    end function get_cmat

    module subroutine get_cmat_ptr( self, cmat_ptr )
        class(image), target,                   intent(in)  :: self
        complex(kind=c_float_complex), pointer, intent(out) :: cmat_ptr(:,:,:)
    end subroutine get_cmat_ptr

    module pure subroutine get_cmat_sub( self, cmat )
        class(image), intent(in)  :: self
        complex,      intent(out) :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
    end subroutine get_cmat_sub

    module pure function get_cmat_at_1( self, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: phys(3)
        complex :: comp
    end function get_cmat_at_1

    module pure function get_cmat_at_2( self, h,k,l ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k,l
        complex :: comp
    end function get_cmat_at_2

    module pure function get_fcomp( self, logi, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3), phys(3)
        complex :: comp
    end function get_fcomp

    module elemental complex function get_fcomp2D(self, h, k)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k
    end function get_fcomp2D

    !--- Setters ---!

    module subroutine set_1( self, logi, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: val
    end subroutine set_1

    module subroutine set_2( self, self2set )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2set
    end subroutine set_2

    module subroutine set_rmat( self, rmat, ft )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(:,:,:)
        logical,      intent(in)    :: ft
    end subroutine set_rmat

    module pure subroutine set_rmat_at( self, i,j,k, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: val
    end subroutine set_rmat_at

    module subroutine set_cmat( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(:,:,:)
    end subroutine set_cmat

    module pure subroutine set_cmat_1( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
    end subroutine set_cmat_1

    module pure subroutine set_cmat_2( self, cval )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cval
    end subroutine set_cmat_2

    module pure subroutine set_cmat_3( self, self2copy )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2copy
    end subroutine set_cmat_3

    module pure subroutine set_cmat_4( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2))
    end subroutine set_cmat_4

    module pure subroutine set_cmat_at_1( self, phys ,comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: comp
    end subroutine set_cmat_at_1

    module pure subroutine set_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer, intent(in) :: h,k,l
        complex, intent(in) :: comp
    end subroutine set_cmat_at_2

    module subroutine set_fcomp( self, logi, phys, comp )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3), phys(3)
        complex,      intent(in)    :: comp
    end subroutine set_fcomp

    module subroutine set_within( self, xyz, radius, val )
        class(image), intent(inout) :: self
        real,         intent(in)    :: xyz(3), radius, val
    end subroutine set_within

    module subroutine set_cmats_from_cmats( self1, self2 , self3, self4, self2set1, self2set2, lims, expcmat3, expcmat4)
        class(image), intent(in)    :: self1, self2, self3, self4
        class(image), intent(inout) :: self2set1, self2set2
        integer,      intent(in)    :: lims(3,2)
        real,         intent(inout) :: expcmat3(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,         intent(inout) :: expcmat4(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
    end subroutine set_cmats_from_cmats

    !--- Slices, sub-images, freq info ---!

    module subroutine get_slice( self3D, slice, self2D )
        class(image), intent(in)    :: self3D
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self2D
    end subroutine get_slice

    module subroutine set_slice( self3D, slice, self2D )
        class(image), intent(in)    :: self2D
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self3D
    end subroutine set_slice

    module subroutine get_subimg(self, binning, xoffset, yoffset, img)
        class(image), intent(in)  :: self
        integer,      intent(in)  :: binning, xoffset, yoffset
        class(image), intent(out) :: img
    end subroutine get_subimg

    module pure function get_lfny( self, which ) result( fnyl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: fnyl
    end function get_lfny

    module pure function get_lhp( self, which ) result( hpl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: hpl
    end function get_lhp

    module pure function get_lp( self, ind ) result( lp )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: lp
    end function get_lp

    module pure function get_spat_freq( self, ind ) result( spat_freq )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: spat_freq
    end function get_spat_freq

    module pure function get_find( self, res ) result( ind )
        class(image), intent(in) :: self
        real,         intent(in) :: res
        integer :: ind
    end function get_find

    !--- Flags ---!

    module function rmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
    end function rmat_associated

    module function cmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
    end function cmat_associated

    module function is_wthreads( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
    end function is_wthreads

    module subroutine set_ft( self, is )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: is
    end subroutine set_ft

    !--- Serialization / misc ---!

    module function serialize_1( self ) result( vec )
        class(image), intent(in) :: self
        real,    allocatable :: vec(:)
        complex, allocatable :: cvec(:)
    end function serialize_1

    module function serialize_2( self, l_msk )result( pcavec )
        class(image), intent(in) :: self
        logical,      intent(in) :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable :: pcavec(:)
    end function serialize_2

    module function serialize_3( self, thres ) result( vec )
        class(image), intent(in) :: self
        real,         intent(in) :: thres
        real, allocatable :: vec(:)
    end function serialize_3

    module subroutine unserialize( self, pcavec, l_msk )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: pcavec(:)
        logical, optional, intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
    end subroutine unserialize

    module subroutine winserialize( self, coord, winsz, pcavec )
        class(image),      intent(inout) :: self
        real, allocatable, intent(inout) :: pcavec(:)
        integer,           intent(in)    :: coord(:), winsz
    end subroutine winserialize

    ! ===== check procedure interfaces =====

    module pure function exists( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
    end function exists

    module pure logical function is_2d(self)
        class(image), intent(in)  ::  self
    end function is_2d

    module pure logical function is_3d(self)
        class(image), intent(in)  ::  self
    end function is_3d

    module pure function even_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep, test(2)
    end function even_dims

    module pure function square_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep
    end function square_dims

    module pure logical function same_dims_1( self1, self2 )
        class(image), intent(in) :: self1, self2
    end function same_dims_1

    module pure logical function same_dims( self, ldim )
        class(image), intent(in) :: self
        integer,      intent(in) :: ldim(3) !< dimensions
    end function same_dims

    module logical pure function same_smpd( self1, self2 )
        class(image), intent(in) :: self1, self2
    end function same_smpd

    module pure function is_ft( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
    end function is_ft

    module function is_empty( self ) result( is )
        class(image), intent(in) :: self
        real    :: minmax(2)
        logical :: is
    end function is_empty

    ! ===== filter procedure interfaces =====

    module subroutine bp( self, hplim, lplim, width )
        class(image), intent(inout) :: self
        real, intent(in)            :: hplim, lplim
        real, intent(in), optional  :: width
    end subroutine bp

    module subroutine lp( self, find, width )
        class(image),   intent(inout) :: self
        integer,        intent(in)    :: find
        real, optional, intent(in)    :: width
    end subroutine lp

    module subroutine lp_background( self, mskvol, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: mskvol
        real,         intent(in)    :: lp
    end subroutine lp_background

    module subroutine bpgau2D( self, hp, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
    end subroutine bpgau2D

    module subroutine bpgau3D( self, hp, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
    end subroutine bpgau3D

    module subroutine tophat( self, shell, halfwidth )
        class(image),   intent(inout) :: self
        integer,        intent(in)    :: shell
        real, optional, intent(in)    :: halfwidth
    end subroutine tophat

    module subroutine phase_rand( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
    end subroutine phase_rand

    module subroutine ran_phases_below_noise_power( self_even, self_odd )
        class(image), intent(inout) :: self_even, self_odd
    end subroutine ran_phases_below_noise_power

    module subroutine whiten_noise_power( self_even, self_odd, is_ptcl )
        class(image), intent(inout) :: self_even, self_odd
        logical,      intent(in)    :: is_ptcl
    end subroutine whiten_noise_power

    module subroutine real_space_filter( self, winsz, which )
        class(image),     intent(inout) :: self
        integer,          intent(in)    :: winsz
        character(len=*), intent(in)    :: which
    end subroutine real_space_filter

    module function hannw( self, oshoot_in ) result( w )
        class(image), intent(inout) :: self
        real, intent(in), optional  :: oshoot_in
        real, allocatable           :: w(:)
    end function hannw

    module subroutine apply_bfac( self, b )
        class(image), intent(inout) :: self
        real,         intent(in)    :: b
    end subroutine apply_bfac

    module subroutine apply_filter_1( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
    end subroutine apply_filter_1

    module subroutine apply_filter_2( self, filter )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: filter(:)
    end subroutine apply_filter_2

    module subroutine apply_filter_serial( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
    end subroutine apply_filter_serial

    module subroutine NLmean2D( self, msk, sdev_noise )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: msk
        real, optional, intent(in)    :: sdev_noise
    end subroutine NLmean2D

    module subroutine NLmean2D_eo( even, odd, avg )
        class(image), intent(inout) :: even, odd, avg
    end subroutine NLmean2D_eo

    module subroutine NLmean3D( self, msk, sdev_noise )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: msk
        real, optional, intent(in)    :: sdev_noise
    end subroutine NLmean3D

    module subroutine NLmean3D_eo( even, odd, avg )
        class(image), intent(inout) :: even, odd, avg
    end subroutine NLmean3D_eo

    module subroutine ICM2D( self, lambda, verbose )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: verbose
    end subroutine ICM2D

    module subroutine ICM2D_eo( even, odd, lambda, verbose )
        class(image),      intent(inout) :: even, odd
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: verbose
    end subroutine ICM2D_eo

    module subroutine ICM3D( self, lambda )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lambda
    end subroutine ICM3D

    module subroutine ICM3D_eo( even, odd, lambda, l_msk )
        class(image),      intent(inout) :: even, odd
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: l_msk(even%ldim(1),even%ldim(2),even%ldim(3))
    end subroutine ICM3D_eo

    module subroutine GLCM( self, nquanta, pmat )
        class(image), intent(in)    :: self
        integer,      intent(in)    :: nquanta
        real,         intent(inout) :: pmat(nquanta,nquanta)
    end subroutine GLCM

    ! ===== frequency analysis procedure interfaces =====

    module subroutine acf( self )
        class(image), intent(inout) :: self
    end subroutine acf

    module function ccf( self1, self2 ) result( cc )
        class(image), intent(inout) :: self1, self2
        type(image) :: cc
    end function ccf

    module subroutine spectrum( self, which, spec, norm )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real, allocatable, intent(inout) :: spec(:)
        logical, optional, intent(in)    :: norm
    end subroutine spectrum

    module subroutine power_spectrum( self, spec )
        class(image), intent(in)    :: self
        real,         intent(inout) :: spec(fdim(self%ldim(1)) - 1)
    end subroutine power_spectrum

    module function guinier_bfac( self, hp, lp ) result( bfac )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
        real :: bfac
    end function guinier_bfac

    module function guinier( self, verbose ) result( plot )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: verbose
        real, allocatable :: plot(:,:)
    end function guinier

    module subroutine fsc( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
    end subroutine fsc

    module subroutine fsc_scaled( self1, self2, sz, corrs )
        class(image), intent(inout) :: self1, self2
        integer,      intent(in)    :: sz
        real,         intent(out)   :: corrs(sz)
    end subroutine fsc_scaled

    module function get_res( self ) result( res )
        class(image), intent(in) :: self
        real, allocatable        :: res(:)
    end function get_res

    module subroutine frc_pspec( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
    end subroutine frc_pspec

    module subroutine fcomps_below_noise_power_stats( self, noise_vol )
        class(image), intent(inout) :: self, noise_vol
    end subroutine fcomps_below_noise_power_stats

    module subroutine resmsk( self, hplim, lplim )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hplim, lplim
    end subroutine resmsk

    module subroutine img2spec( self, speckind, lp_backgr_subtr, img_out, postproc )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        logical, optional, intent(in)    :: postproc
    end subroutine img2spec

    module subroutine mic2spec( self, box, speckind, lp_backgr_subtr, img_out, postproc )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: box
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        logical, optional, intent(in)    :: postproc
    end subroutine mic2spec

    module subroutine pspec_graphene_mask( self, ldim, smpd )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        real,         intent(in)    :: smpd
    end subroutine pspec_graphene_mask

    module subroutine dampen_pspec_central_cross( self )
        class(image), intent(inout) :: self
    end subroutine dampen_pspec_central_cross

    ! ===== calculator procedure interfaces =====

    !--- Basic stats / local stats ---!

    module function get_sum_int(self) result(sum_int)
        class(image), intent(in) :: self
        real :: sum_int
    end function get_sum_int

    module function minmax( self, radius )result( mm )
        class(image),   intent(in) :: self
        real, optional, intent(in) :: radius
        real    :: mm(2)
    end function minmax

    module subroutine loc_sdev( self, winsz, sdevimg, asdev )
        class(image),   intent(in)    :: self
        integer,        intent(in)    :: winsz
        class(image),   intent(inout) :: sdevimg
        real, optional, intent(inout) :: asdev
    end subroutine loc_sdev

    module function avg_loc_sdev( self, winsz ) result( asdev )
        class(image), intent(in) :: self
        integer,      intent(in) :: winsz
        real :: asdev
    end function avg_loc_sdev

    module subroutine loc_var( self, varimg, avar )
        class(image),   intent(in)    :: self
        class(image),   intent(inout) :: varimg
        real, optional, intent(inout) :: avar
    end subroutine loc_var

    module subroutine loc_var3D( self, varimg, avar )
        class(image),   intent(in)    :: self
        class(image),   intent(inout) :: varimg
        real, optional, intent(inout) :: avar
    end subroutine loc_var3D

    module subroutine rmsd( self, dev, mean )
        class(image),   intent(inout) :: self
        real,           intent(out)   :: dev
        real, optional, intent(out)   :: mean
    end subroutine rmsd

    module subroutine stats_1( self, which, ave, sdev, maxv, minv, msk, med, errout )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: ave, sdev, maxv, minv
        real,    optional, intent(in)    :: msk
        real,    optional, intent(out)   :: med
        logical, optional, intent(out)   :: errout
    end subroutine stats_1

    module subroutine stats_2( self, ave, sdev, maxv, minv, mskimg, med, errout )
        class(image),           intent(inout) :: self
        real,                   intent(out)   :: ave, sdev, maxv, minv
        class(image), optional, intent(in)    :: mskimg
        real,         optional, intent(out)   :: med
        logical,      optional, intent(out)   :: errout
    end subroutine stats_2

    module function variance( self ) result( var )
        class(image), intent(in) :: self
        real    :: ave, var, ep, rmat_subtr_avg(self%ldim(1),self%ldim(2),self%ldim(3))
    end function variance

    module real function skew( self, mask )
        class(image),      intent(in)  :: self
        logical, optional, intent(in)  :: mask(self%ldim(1),self%ldim(2))
    end function skew

    module real function kurt( self, mask )
        class(image),      intent(in)  :: self
        logical, optional, intent(in)  :: mask(self%ldim(1),self%ldim(2))
    end function kurt

    module function noisesdev( self, msk ) result( sdev )
        use simple_online_var, only: online_var
        class(image), intent(inout) :: self
        real,         intent(in)    :: msk
        real :: sdev
    end function noisesdev

    module function mean( self ) result( avg )
        class(image), intent(inout) :: self
        real :: avg
    end function mean

    module logical function contains_nans( self )
        class(image), intent(in) :: self
    end function contains_nans

    module subroutine checkimg4nans( self )
        class(image), intent(in) :: self
    end subroutine checkimg4nans

    module subroutine cure( self, maxv, minv, ave, sdev, n_nans )
        class(image), intent(inout) :: self
        real,         intent(out)   :: maxv, minv, ave, sdev
        integer,      intent(out)   :: n_nans
    end subroutine cure

    module function dead_hot_positions( self, frac ) result( pos )
        class(image), intent(in) :: self
        real,         intent(in) :: frac
        logical, allocatable     :: pos(:,:)
    end function dead_hot_positions

    !--- Gradients / geometry ---!

    module subroutine calc_gradient(self, grad, Dc, Dr)
        class(image),   intent(inout) :: self
        real,           intent(out)   :: grad(self%ldim(1), self%ldim(2), self%ldim(3)) ! gradient matrix
        real, optional, intent(out)   :: Dc(self%ldim(1), self%ldim(2), self%ldim(3)), Dr(self%ldim(1), self%ldim(2), self%ldim(3)) ! derivates column and row matrices
    end subroutine calc_gradient

    module subroutine gradients_magnitude( self, self_out )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: self_out
    end subroutine gradients_magnitude

    module subroutine gradient(self, Dc, Dr, Dz, grad)
        class(image),   intent(inout) :: self
        real, optional, intent(out)   :: Dc(self%ldim(1), self%ldim(2), self%ldim(3)), & ! derivates column matrix
                                         Dr(self%ldim(1), self%ldim(2), self%ldim(3)), & ! derivates row matrix
                                         Dz(self%ldim(1), self%ldim(2), self%ldim(3))    ! derivates z matrix
        real, optional, intent(out)   :: grad(self%ldim(1), self%ldim(2), self%ldim(3))  ! gradient matrix
    end subroutine gradient

    module subroutine calc_ice_score( self, score )
        class(image), intent(in)  :: self
        real,         intent(out) :: score
    end subroutine calc_ice_score

    module subroutine calc_principal_axes_rotmat( self, radius, R )
        class(image), intent(in)  :: self
        real,         intent(in)  :: radius
        real,         intent(out) :: R(3,3)
    end subroutine calc_principal_axes_rotmat

    !--- Physical coords helpers ---!

    module function loop_lims( self, mode, lp_dyn ) result( lims )
        class(image), intent(in)   :: self
        integer, intent(in)        :: mode
        real, intent(in), optional :: lp_dyn
        integer                    :: lims(3,2)
    end function loop_lims

    module pure function comp_addr_phys1(self,logi) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3) !<  Logical address
        integer                  :: phys(3) !<  Physical address
    end function comp_addr_phys1

    module pure function comp_addr_phys2(self,h,k,m) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k,m   !<  Logical address
        integer                  :: phys(3) !<  Physical address
    end function comp_addr_phys2

    module pure function comp_addr_phys3(self,h,k) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k     !<  Logical address
        integer                  :: phys(2) !<  Physical address
    end function comp_addr_phys3

    !--- Correlation / distances ---!

    module function corr( self1, self2, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self1, self2
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real :: r
    end function corr

    module function corr_shifted( self_ref, self_ptcl, shvec, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self_ref, self_ptcl
        real,           intent(in)    :: shvec(3)
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real :: r
    end function corr_shifted

    module function real_corr_1( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real :: r
    end function real_corr_1

    module function real_corr_2( self1, self2, mask ) result( r )
        class(image), intent(inout) :: self1, self2
        logical,      intent(in)    :: mask(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: r
    end function real_corr_2

    module function euclid_dist_two_imgs(self1, self2, mask1) result(dist)
        class(image),      intent(inout) :: self1, self2
        logical, optional, intent(in)    :: mask1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: dist
    end function euclid_dist_two_imgs

    module subroutine phase_corr(self1, self2, pc, lp )
        class(image),      intent(inout) :: self1, self2, pc
        real,              intent(in)    :: lp
    end subroutine phase_corr

    module subroutine fcorr_shift( self1, self2, trs, shift, peak_interp )
        class(image),      intent(inout) :: self1, self2
        real,              intent(in)    :: trs
        real,              intent(inout) :: shift(2)
        logical, optional, intent(in)    :: peak_interp
    end subroutine fcorr_shift

    module subroutine fcorr_shift3D( self1, self2, trs, shift, peak_interp )
        class(image),      intent(inout) :: self1, self2
        real,              intent(in)    :: trs
        real,              intent(inout) :: shift(3)
        logical, optional, intent(in)    :: peak_interp
    end subroutine fcorr_shift3D

    module subroutine prenorm4real_corr_1( self, sxx )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
    end subroutine prenorm4real_corr_1

    module subroutine prenorm4real_corr_2( self, sxx, mask )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        logical,      intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))
    end subroutine prenorm4real_corr_2

    module subroutine prenorm4real_corr_3( self, err )
        class(image), intent(inout) :: self
        logical,      intent(out)   :: err
    end subroutine prenorm4real_corr_3

    module function real_corr_prenorm_1( self_ref, self_ptcl, sxx_ref ) result( r )
        class(image), intent(inout) :: self_ref, self_ptcl
        real,         intent(in)    :: sxx_ref
        real :: r
    end function real_corr_prenorm_1

    module function real_corr_prenorm_2( self_ref, self_ptcl, sxx_ref, mask ) result( r )
        class(image), intent(inout) :: self_ref, self_ptcl
        real,         intent(in)    :: sxx_ref
        logical,      intent(in)    :: mask(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: r
    end function real_corr_prenorm_2

    module real function real_corr_prenorm_3( self_ref, self_ptcl )
        class(image), intent(inout) :: self_ref, self_ptcl
    end function real_corr_prenorm_3

    module subroutine radial_cc( self1, self2, self_w, smpd, rad_corrs, rad_dists )
        class(image), intent(inout):: self1, self2, self_w
        real,         intent(in)   :: smpd
        real,         intent(out)  :: rad_corrs(int(self1%ldim(1)/2.)), rad_dists(int(self1%ldim(1)/2.))
    end subroutine radial_cc

    module function sqeuclid( self1, self2, mask ) result( r )
        class(image), intent(inout) :: self1, self2
        logical,      intent(in)    :: mask(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: r
    end function sqeuclid

    module subroutine sqeuclid_matrix_1( self1, self2, sqdiff )
        class(image), intent(in)    :: self1, self2
        real,         intent(inout) :: sqdiff(self1%ldim(1),self1%ldim(2),self1%ldim(3))
    end subroutine sqeuclid_matrix_1

    module subroutine sqeuclid_matrix_2( self1, self2, sqdiff_img )
        class(image), intent(in)    :: self1, self2
        class(image), intent(inout) :: sqdiff_img
    end subroutine sqeuclid_matrix_2

    module function euclid_norm( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real :: r 
    end function euclid_norm

    !--- cost / shift ---!

    module subroutine opt_filter_costfun( even_filt, odd_raw, odd_filt, even_raw, sqdiff_img )
        class(image), intent(in)    :: even_filt, odd_raw, odd_filt, even_raw
        class(image), intent(inout) :: sqdiff_img
    end subroutine opt_filter_costfun

    module subroutine opt_filter_costfun_workshare( even_filt, odd_raw, odd_filt, even_raw, sqdiff_img )
        class(image), intent(in)    :: even_filt, odd_raw, odd_filt, even_raw
        class(image), intent(inout) :: sqdiff_img
    end subroutine opt_filter_costfun_workshare

    module pure function oshift_1( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        real,         intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
    end function oshift_1

    module pure function oshift_2( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
    end function oshift_2

    module function gen_argtransf_comp( self, logi, ldim ) result( arg )
        class(image), intent(in)      :: self
        real, intent(in)              :: logi(3)
        integer, intent(in), optional :: ldim
        real                          :: arg(3)
    end function gen_argtransf_comp

    ! ===== visualization procedure interfaces =====

    module subroutine print_cmat( self )
        class(image), intent(in) :: self
    end subroutine print_cmat

    module subroutine print_rmat( self )
        class(image), intent(in) :: self
    end subroutine print_rmat

    module subroutine vis( self, sect, geomorsphr )
        class(image),      intent(in) :: self
        integer, optional, intent(in) :: sect
        logical, optional, intent(in) :: geomorsphr !< geometrical or spherical complex format
    end subroutine vis

    module subroutine before_after( left, right, ba, mask )
        class(image),      intent(in)    :: left, right
        type(image),       intent(inout) :: ba
        logical, optional, intent(in)    :: mask(left%ldim(1),left%ldim(2),left%ldim(3))
    end subroutine before_after

    module subroutine scale_pspec4viz( self, rsmpd4viz )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: rsmpd4viz
    end subroutine scale_pspec4viz

    module subroutine generate_orthogonal_reprojs( self, reprojs )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: reprojs
    end subroutine generate_orthogonal_reprojs

    module subroutine collage( self1, self2, img_out )
        class(image), intent(inout) :: self1, self2, img_out
    end subroutine collage

    module subroutine tile( self, stkimg, x, y)
        class(image), intent(inout) :: self
        class(image), intent(inout) :: stkimg
        integer,      intent(in)    :: x, y
    end subroutine tile

    ! ===== segmentation procedure interfaces =====

    module function nforeground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
    end function nforeground

    module function nbackground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
    end function nbackground

    module subroutine binarize_1( self_in, thres, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: thres
        class(image), optional, intent(inout) :: self_out
    end subroutine binarize_1

    module subroutine binarize_2( self, npix )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: npix
    end subroutine binarize_2

    module subroutine binarize_3( self, thres, mask )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        logical,      intent(inout) :: mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
    end subroutine binarize_3

    module subroutine cendist( self, c_point )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: c_point(3)
    end subroutine cendist

    module subroutine bin_inv( self )
        class(image), intent(inout) :: self
    end subroutine bin_inv

    module subroutine remove_edge( self )
        class(image), intent(inout) :: self
    end subroutine remove_edge

    module subroutine zero2one( self )
        class(image), intent(inout) :: self
    end subroutine zero2one

    module subroutine one_at_edge( self )
        class(image), intent(inout) :: self
    end subroutine one_at_edge

    module function bin2logical( self ) result( mask )
        class(image), intent(in)  :: self
        logical,      allocatable :: mask(:,:,:)
    end function bin2logical

    module subroutine logical2bin( self, mask )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
    end subroutine logical2bin

    module subroutine density_inoutside( self, msk, nin, nout, nmsk )
        class(image), intent(in)  :: self
        real,         intent(in)  :: msk
        integer,      intent(out) :: nin, nout, nmsk
    end subroutine density_inoutside

    module subroutine calc_bin_thres( self, frac_fg_target, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: frac_fg_target
        real,         intent(out)   :: thres 
    end subroutine calc_bin_thres

    module subroutine mask( self, mskrad, which, inner, width, backgr )
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        character(len=*), intent(in)    :: which
        real, optional,   intent(in)    :: inner, width, backgr
    end subroutine mask

    module subroutine taper_edges( self )
        class(image), intent(inout) :: self
    end subroutine taper_edges

    module subroutine taper_edges_hann( self, borders )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: borders(2)
    end subroutine taper_edges_hann

    ! ===== image operations procedure interfaces =====

    !--- CTF ---!

    module subroutine ctf_dens_correct( self_sum, self_rho )
        class(image), intent(inout) :: self_sum
        class(image), intent(inout) :: self_rho
    end subroutine ctf_dens_correct

    module subroutine ctf_dens_correct_wiener( self_sum, self_rho, ssnr )
        class(image), intent(inout) :: self_sum
        class(image), intent(in)    :: self_rho
        real,         intent(in)    :: ssnr(:)
    end subroutine ctf_dens_correct_wiener

    !--- Insertions ---!

    module subroutine insert(self_in, coord, self_out )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: coord(2)
        type(image),  intent(inout) :: self_out
    end subroutine insert

    module subroutine insert_lowres( self, self2insert, find )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2insert
        integer,      intent(in)    :: find
    end subroutine insert_lowres

    module subroutine insert_lowres_serial( self, self2insert, find )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2insert
        integer,      intent(in)    :: find
    end subroutine insert_lowres_serial

    !--- Noise ---!

    module subroutine ran( self, b )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: b
    end subroutine ran

    module subroutine gauran( self, mean, sdev )
        class(image), intent(inout) :: self
        real,         intent(in)    :: mean, sdev
    end subroutine gauran

    module subroutine add_gauran( self, snr )
        class(image), intent(inout) :: self
        real,         intent(in)    :: snr
    end subroutine add_gauran

    module subroutine salt_n_pepper( self, pos )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: pos(:,:)
    end subroutine salt_n_pepper

    !--- Background ---!

    module subroutine div_w_instrfun( self, interpfun, alpha, padded_dim )
        class(image),           intent(inout) :: self
        character(len =*),      intent(in)    :: interpfun
        real,         optional, intent(in)    :: alpha
        integer,      optional, intent(in)    :: padded_dim
    end subroutine div_w_instrfun

    module subroutine estimate_background( self, freq, backgr, mode )
        class(image),     intent(in)    :: self
        real,             intent(in)    :: freq
        class(image),     intent(inout) :: backgr
        character(len=*), intent(in)    :: mode
    end subroutine estimate_background

    module subroutine subtr_backgr( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
    end subroutine subtr_backgr

    module subroutine subtr_backgr_ramp( self, lmsk )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
    end subroutine subtr_backgr_ramp

    module subroutine subtract_background( self, freq, mode )
        class(image),               intent(inout) :: self
        real,                       intent(in)    :: freq
        character(len=*), optional, intent(in)    :: mode
    end subroutine subtract_background

    module subroutine upsample_square_background( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(2)
    end subroutine upsample_square_background

    !--- Arithmetics ---!

    module subroutine remove_neg( self )
        class(image), intent(inout) :: self
    end subroutine remove_neg

    module subroutine neg( self )
        class(image), intent(inout) :: self
    end subroutine neg

    module subroutine div_below( self, thres, val )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres, val
    end subroutine div_below

    module subroutine inv( self )
        class(image), intent(inout) :: self
    end subroutine inv

    !--- Zeroing ---!

    module subroutine zero(self)
        class(image), intent(inout) :: self
    end subroutine zero

    module subroutine zero_and_flag_ft(self)
        class(image), intent(inout) :: self
    end subroutine zero_and_flag_ft

    module subroutine zero_and_unflag_ft(self)
        class(image), intent(inout) :: self
    end subroutine zero_and_unflag_ft

    module subroutine zero_background( self )
        class(image), intent(inout) :: self
    end subroutine zero_background

    module subroutine zero_below( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
    end subroutine zero_below

    module subroutine zero_edgeavg( self )
        class(image), intent(inout) :: self
    end subroutine zero_edgeavg

    module subroutine zero_env_background( self, volmsk )
        class(image), intent(inout) :: self, volmsk
    end subroutine zero_env_background

    module subroutine zero_neg( self )
        class(image), intent(inout) :: self
    end subroutine zero_neg

    ! ===== geometric procedure interfaces =====

    !--- windowing ---!

    module subroutine window( self_in, coord, box, self_out, noutside )
        class(image),      intent(in)    :: self_in
        integer,           intent(in)    :: coord(2), box
        class(image),      intent(inout) :: self_out
        integer, optional, intent(inout) :: noutside
    end subroutine window

    module subroutine window_slim( self_in, coord, box, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: coord(:), box
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
    end subroutine window_slim

    module subroutine window_center( self_in, center, rad, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: center(:), rad
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
    end subroutine window_center

    module subroutine add_window( self, imgwin, coord, offset )
        class(image),      intent(inout) :: self
        class(image),      intent(in)    :: imgwin
        integer,           intent(in)    :: coord(2)
        integer, optional, intent(in)    :: offset 
    end subroutine add_window

    module function win2arr( self, i, j, k, winsz ) result( pixels )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i, j, k, winsz
        real, allocatable :: pixels(:)
    end function win2arr

    module subroutine win2arr_rad( self, i, j, k, winsz, npix_in, maxrad, npix_out, pixels )
        class(image), intent(in)    :: self
        integer,      intent(in)    :: i, j, k, winsz, npix_in
        real,         intent(in)    :: maxrad ! in pixels
        integer,      intent(inout) :: npix_out
        real,         intent(inout) :: pixels(npix_in)
    end subroutine win2arr_rad

    !--- shapes ---!

    module subroutine corner( self_in, box, self_out )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: box
        type(image),  intent(inout) :: self_out
    end subroutine corner

    module subroutine corners( self, sqrad )
        class(image), intent(inout) :: self
        integer, intent(in)         :: sqrad
    end subroutine corners

    module subroutine gauimg_1( self, wsz , alpha )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: wsz
        real, intent(in), optional  :: alpha
    end subroutine gauimg_1

    module subroutine gauimg_2( self, wsz, offx,offy)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: wsz, offx, offy
    end subroutine gauimg_2

    module subroutine gauimg2D( self, xsigma, ysigma, cutoff )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: xsigma, ysigma
        real, optional, intent(in)    :: cutoff
    end subroutine gauimg2D

    module subroutine gauimg3D( self, xsigma, ysigma, zsigma, cutoff )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: xsigma, ysigma, zsigma
        real, optional, intent(in)    :: cutoff
    end subroutine gauimg3D
    
    module subroutine reshape2cube( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(out)   :: self_out
    end subroutine reshape2cube

    module subroutine square( self, sqrad )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: sqrad
    end subroutine square

    !--- pad/clip/crop ---!

    module subroutine pad( self_in, self_out, backgr, antialiasing )
        class(image),      intent(inout) :: self_in, self_out
        real,    optional, intent(in)    :: backgr
        logical, optional, intent(in)    :: antialiasing
    end subroutine pad

    module subroutine pad_inplace( self, ldim, backgr, antialiasing )
        class(image),   intent(inout) :: self
        integer,           intent(in) :: ldim(3)
        real,    optional, intent(in) :: backgr
        logical, optional, intent(in) :: antialiasing
    end subroutine pad_inplace

    module subroutine pad_mirr_1( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
    end subroutine pad_mirr_1

    module subroutine pad_mirr_2( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
    end subroutine pad_mirr_2

    module subroutine clip( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
    end subroutine clip

    module subroutine clip_inplace( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
    end subroutine clip_inplace

    module subroutine read_and_crop( self, volfname, smpd, box_crop, smpd_crop )
        class(image),   intent(inout) :: self
        class(string),  intent(in)    :: volfname
        integer,        intent(in)    :: box_crop
        real,           intent(in)    :: smpd, smpd_crop
    end subroutine read_and_crop

    !--- flip/mirror/rotate/shift ---!

    module recursive subroutine flip( self, mode )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: mode
    end subroutine flip

    module subroutine mirror( self, md, fourier )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: md
        logical, optional, intent(in)    :: fourier
    end subroutine mirror

    module function calc_shiftcen( self, lp, msk, hp ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real, optional, intent(in)    :: msk, hp
        real :: xyz(3)
    end function calc_shiftcen

    module function calc_shiftcen_serial( self, lp, msk, hp ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real,           intent(in)    :: msk
        real, optional, intent(in)    :: hp
        real :: xyz(3)
    end function calc_shiftcen_serial

    module subroutine roavg( self, angstep, avg, ang_stop )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: angstep
        class(image),      intent(inout) :: avg
        integer, optional, intent(in)    :: ang_stop
    end subroutine roavg

    module subroutine rtsq( self_in, ang, shxi, shyi, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: ang,shxi,shyi
        class(image), optional, intent(inout) :: self_out
    end subroutine rtsq

    module subroutine rtsq_serial( self_in, ang, shxi, shyi, rmat_out )
        class(image), intent(inout) :: self_in
        real,         intent(in)    :: ang,shxi,shyi
        real,         intent(inout) :: rmat_out(self_in%ldim(1),self_in%ldim(2),1)
    end subroutine rtsq_serial

    module subroutine shift_phorig( self )
        class(image), intent(inout) :: self
    end subroutine shift_phorig

    module subroutine shift( self, shvec )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(3)
    end subroutine shift

    module subroutine shift2Dserial_1( self, shvec  )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(2)
    end subroutine shift2Dserial_1

    module subroutine shift2Dserial_2( self, shvec, self_out )
        class(image), intent(inout) :: self, self_out
        real,         intent(in)    :: shvec(2)
    end subroutine shift2Dserial_2

    module subroutine masscen( self, xyz, mask_in )
        class(image),      intent(inout) :: self
        real        ,      intent(out)   :: xyz(3)
        logical, optional, intent(in)    :: mask_in(:,:,:)
    end subroutine masscen

    module subroutine masscen_adjusted( self, xyz, mask_in )
        class(image),      intent(inout) :: self
        real        ,      intent(out)   :: xyz(3)
        logical, optional, intent(in)    :: mask_in(:,:,:)
    end subroutine masscen_adjusted

    ! ===== normalization procedure interfaces =====

    module subroutine scale_pixels(self, new_range, ssc, oold_range)
        class(image),   intent(inout) :: self
        real,           intent(in)    :: new_range(2)
        real, optional, intent(out)   :: oold_range(2), ssc
    end subroutine scale_pixels

    module subroutine norm( self, a_s )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: a_s(2)
    end subroutine norm

    module subroutine norm_minmax( self  )
        class(image), intent(inout) :: self
    end subroutine norm_minmax

    module subroutine norm4viz( self, brightness, maxmin)
        class(image),      intent(inout) :: self
        real,    optional, intent(in)    :: brightness
        logical, optional, intent(in)    :: maxmin
    end subroutine norm4viz

    module subroutine norm_ext( self, avg, sdev )
        class(image), intent(inout) :: self
        real,         intent(in)    :: avg, sdev
    end subroutine norm_ext

    module subroutine norm_noise( self, lmsk, sdev_noise )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3)) ! foreground must be true
        real,         intent(inout) :: sdev_noise
    end subroutine norm_noise

    module subroutine norm_within( self, mask )
        class(image),      intent(inout) :: self
        logical,           intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))    ! foreground must be true
    end subroutine norm_within

    module subroutine cure_outliers( self, ncured, nsigma, deadhot, outliers )
        class(image),                   intent(inout) :: self
        integer,                        intent(inout) :: ncured
        real,                           intent(in)    :: nsigma
        integer,                        intent(out)   :: deadhot(2)
        logical, allocatable, optional, intent(out)   :: outliers(:,:)
    end subroutine cure_outliers

    module subroutine quantize_fwd( self, nquanta, transl_tab, l_msk )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: nquanta
        real,              intent(inout) :: transl_tab(nquanta)
        logical, optional, intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3)) 
    end subroutine quantize_fwd

    module subroutine quantize_bwd( self, nquanta, transl_tab )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nquanta
        real,         intent(in)    :: transl_tab(nquanta)
    end subroutine quantize_bwd

end interface

interface image
    module procedure constructor
end interface image

type :: image_stack
    type(image), allocatable :: stack(:)
end type image_stack

! CLASS PARAMETERS/VARIABLES
logical,     parameter   :: shift_to_phase_origin=.true.
type(image), allocatable :: thread_safe_tmp_imgs(:)

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor( ldim, smpd, wthreads ) result( self ) !(FAILS W PRESENT GFORTRAN)
        integer,           intent(in) :: ldim(:)
        real,              intent(in) :: smpd
        logical, optional, intent(in) :: wthreads
        type(image) :: self
        call self%new( ldim, smpd, wthreads )
    end function constructor

    subroutine ring_stats( self, stats )
        class(image),                    intent(inout) :: self
        type(stats_struct), allocatable, intent(inout) :: stats(:)
        integer     :: iring, npix, nrings
        type(image) :: maskimg
        nrings = fdim(self%ldim(1))-3
        if( allocated(stats) ) deallocate(stats)
        allocate(stats(nrings))
        do iring = 1, nrings
            call maskimg%ring(self%ldim, self%smpd, real(iring+2), real(iring-2), npix)
            call self%stats_2(stats(iring)%avg, stats(iring)%sdev, stats(iring)%maxv, stats(iring)%minv, maskimg)
        end do
    end subroutine ring_stats

    !>  \brief ellipse constructs an ellipse of given axes.
    !    optional parameter 'hole' (yes|no) allows the user to choose
    !    between the full ellipse or just its borders. Default: full.
    !    It has the optiono f rotating the ellipse. The ellipse is built 'on top' of the image,
    !    it doesn't cancel what is already present
    subroutine ellipse(self, center, axes, hole)
        class(image),               intent(inout) :: self
        real,                       intent(in)    :: axes(2)
        integer,                    intent(in)    :: center(2)
        character(len=*), optional, intent(in)    :: hole
        integer :: i, j
        real, allocatable :: rmat_t(:,:,:)
        if(.not. self%existence) THROW_HARD('Image has to be created before; ellipse')
        if(self%ldim(3) /= 1) THROW_HARD('For 2D images only; ellipse')
        rmat_t = self%get_rmat()
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                if((real(i)-real(center(1)))**2/(axes(1)**2) + (real(j)-real(center(2)))**2/(axes(2)**2) - 1 < TINY) then
                    if( maxval(self%rmat(:,:,:)) - minval(self%rmat(:,:,:)) > TINY) then
                      self%rmat(i,j,1) = maxval(self%rmat(:,:,:))
                    else
                      self%rmat(i,j,1) = 1.
                    endif
                endif
            enddo
        enddo
        if(present(hole)) then
            if(hole .eq. 'yes') then
              do i = 1, self%ldim(1)
                  do j = 1, self%ldim(2)
                      if((real(i)-real(center(1)))**2/(axes(1)-1)**2 + (real(j)-real(center(2)))**2/(axes(2)-1)**2 - 1 < TINY) then
                          if( maxval(self%rmat(:,:,:)) - minval(self%rmat(:,:,:)) > TINY) then
                              self%rmat(i,j,1) = rmat_t(i,j,1) !minval(self%rmat(:,:,:))
                          else
                              self%rmat(i,j,1) = rmat_t(i,j,1)
                          endif
                        endif
                  enddo
              enddo
            elseif(hole .ne. 'no') then
              THROW_HARD('Input error for hole parameter; ellipse')
            endif
          endif
    end subroutine ellipse

    !>  \brief  is the image class unit test
    subroutine test_image( doplot )
        logical, intent(in)  :: doplot
        write(logfhandle,'(a)') '**info(simple_image_unit_test): testing square dimensions'
        call test_image_local( 100, 100, 100, doplot )
        !        write(logfhandle,'(a)') '**info(simple_image_unit_test): testing non-square dimensions'
        !        call test_image_local( 120, 90, 80, doplot )
        write(logfhandle,'(a)') 'SIMPLE_IMAGE_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'

    contains

        subroutine test_image_local( ld1, ld2, ld3, doplot )
            integer, intent(in)  :: ld1, ld2, ld3
            logical, intent(in)  :: doplot
            type(image)          :: img, img_2, img_3, img_4, img3d
            type(image)          :: imgs(20)
            integer              :: i, j, k, cnt, ldim(3)
            real                 :: input, ave, sdev, med
            real                 :: corr, corr_lp, maxv, minv
            real                 :: smpd=2.
            logical              :: passed, test(6)

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 1): testing basal constructors'
            call img%new([ld1,ld2,1], 1.)
            call img_3%new([ld1,ld2,1], 1.)
            call img3d%new([ld1,ld2,ld3], 1.)
            if( .not. img%exists() )   THROW_HARD('ERROR, in constructor or in exists function, 1')
            if( .not. img3d%exists() ) THROW_HARD('ERROR, in constructor or in exists function, 2')

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 2): testing getters/setters'
            passed = .true.
            cnt = 1
            do i=1,ld1
                do j=1,ld2
                    input = real(cnt)
                    call img%set([i,j,1], input)
                    if( img%get([i,j,1]) /= input) passed = .false.
                    do k=1,ld3
                        input = real(cnt)
                        call img3d%set([i,j,k],input)
                        if( img3d%get([i,j,k]) /= input) passed = .false.
                        cnt = cnt+1
                    end do
                end do
            end do
            if( .not. passed )  THROW_HARD('getters/setters test failed')

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 4): testing checkups'
            img_2 = img
            test(1) = img%even_dims()
            if( ld1 == ld2 )then
                test(2) = img%square_dims()
            else
                test(2) = .not. img%square_dims()
            endif
            test(3) = img.eqdims.img_2
            test(4) = img.eqsmpd.img_2
            test(5) = img%is_2d()
            test(6) = .not. img%is_3d()
            passed = all(test)
            if( .not. passed ) then
                write(logfhandle,*) ""
                write(logfhandle,*) ' checkups ', test
                write(logfhandle,*) ""
                THROW_HARD('checkups test failed')
            endif
            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 6): testing stats'
            passed = .false.
            call img%gauran( 5., 15. )
            call img%stats( 'foreground', ave, sdev, maxv, minv, 40., med )
            if( ave >= 4. .and. ave <= 6. .and. sdev >= 14. .and.&
                sdev <= 16. .and. med >= 4. .and. med <= 6. ) passed = .true.
            if( .not. passed )  THROW_HARD('stats test failed')

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 9): testing lowpass filter'
            call img%square( 10 )
            if( doplot ) call img%vis
            call img%bp(0., 5.)
            if( doplot ) call img%vis
            call img%bp(0., 10.)
            if( doplot ) call img%vis
            call img%bp(0., 20.)
            if( doplot ) call img%vis
            call img%bp(0., 30.)
            if( doplot ) call img%vis

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 10): testing spherical mask'
            call img%ran
            if( doplot ) call img%vis
            call img%mask(35.,'hard')
            if( doplot ) call img%vis
            call img%ran
            call img%mask(35.,'soft')
            if( doplot ) call img%vis

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 13): testing bicubic rots'
            cnt = 0
            call img_3%square(20)
            if( ld1 == ld2 )then
                call img_4%new([ld1,ld2,1], 1.)
                do i=0,360,30
                    call img_3%rtsq(real(i), 0., 0., img_4)
                    cnt = cnt+1
                    if( doplot ) call img_4%vis
                end do
            endif

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 14): testing binary imgproc routines'
            passed = .false.
            call img%gauimg(20)
            call img%norm_minmax
            if( doplot ) call img%vis
            call img%binarize(0.5)
            if( doplot ) call img%vis
            call img%gauimg(20)
            call img%binarize(500)
            if( doplot ) call img%vis

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 15): testing auto correlation function'
            call img%square( 10 )
            if( doplot ) call img%vis
            call img%acf
            if( doplot ) call img%vis
            call img%square( 10 )
            call img%shift([5.,-5.,0.])
            if( doplot ) call img%vis
            call img%acf
            if( doplot ) call img%vis

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 16): testing correlation functions'
            passed = .false.
            ldim = [100,100,1]
            call img%new(ldim, smpd)
            call img_2%new(ldim, smpd)
            call img%gauimg(10)
            call img%fft()
            call img_2%gauimg(13)
            call img_2%fft()
            corr = img%corr(img_2)
            corr_lp = img%corr(img_2,20.)
            if( corr > 0.96 .and. corr < 0.98 .and. corr_lp > 0.96 .and. corr_lp < 0.98 ) passed = .true.
            if( .not. passed ) THROW_HARD('corr test failed')

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 17): testing downscaling'
            if( ld1 == ld2 )then
                call img%gauimg(20)
                if( doplot )  call img%vis
                if( doplot ) call img_2%vis
            endif

            if( img%square_dims() .and. nthr_glob > 2 )then
                write(logfhandle,'(a)') '**info(simple_image_unit_test, part 19): testing rotational averager'
                call img%square( 10 )
                if( doplot ) call img%vis
                call img%roavg(5,img_2)
                if( doplot ) call img_2%vis
            endif

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 20): testing the read/write capabilities'
            ! create a square
            ldim = [120,120,1]
            call img%new(ldim, smpd)
            call img%square(20)
            ! write stacks of 5 squares
            do i=1,5
                call img%write(string('squares_spider.spi'),i)
                call img%write(string('squares_mrc.mrc'),i)
            end do
            ! convert the squares from SPIDER to MRC & vice versa
            do i=1,5
                call img%read(string('squares_spider.spi'),i)
                call img%write(string('squares_spider_converted.mrc'),i)
                call img%read(string('squares_mrc.mrc'),i)
                call img%write(string('squares_mrc_converted.spi'),i)
            end do
            ! test SPIDER vs. MRC & converted vs. nonconverted
            do i=1,20
                call imgs(i)%new(ldim, smpd)
            end do
            cnt = 0
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read(string('squares_spider.spi'),i)
            end do
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read(string('squares_spider_converted.mrc'),i)
            end do
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read(string('squares_mrc.mrc'),i)
            end do
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read(string('squares_mrc_converted.spi'),i)
            end do
            do i=1,19
                do j=i+1,20
                    corr = imgs(i)%corr(imgs(j))
                    if( corr < 0.99999 )then
                        THROW_HARD('SPIDER vs. MRC & converted vs. nonconverted test failed')
                    endif
                end do
            end do
            ! create a cube
            ldim = [120,120,120]
            call img%new(ldim, smpd)
            call img%square(20)
            ! write volume files
            do i=1,5
                call img%write(string('cube_spider.spi'))
                call img%write(string('cube_mrc.mrc'))
            end do
            ! convert the cubes from SPIDER to MRC & vice versa
            do i=1,5
                call img%read(string('cube_spider.spi'))
                call img%write(string('cube_spider_converted.mrc'))
                call img%read(string('cube_mrc.mrc'))
                call img%write(string('cube_mrc_converted.spi'))
            end do
            ! test SPIDER vs. MRC & converted vs. nonconverted
            do i=1,4
                call imgs(i)%new(ldim, smpd)
                call imgs(i)%read(string('cube_spider.spi'))
                call imgs(i)%read(string('cube_spider_converted.mrc'))
                call imgs(i)%read(string('cube_mrc.mrc'))
                call imgs(i)%read(string('cube_mrc_converted.spi'))
            end do
            do i=1,3
                do j=i+1,4
                    corr = imgs(i)%corr(imgs(j))
                    if( corr < 0.99999 )then
                        THROW_HARD('SPIDER vs. MRC & converted vs. nonconverted test failed')
                    endif
                end do
            end do

            write(logfhandle,'(a)') '**info(simple_image_unit_test, part 21): testing destructor'
            passed = .false.
            call img%kill()
            call img3d%kill()
            test(1) = .not. img%exists()
            test(2) = .not. img3d%exists()
            passed = all(test)
            if( .not. passed )  THROW_HARD('destructor test failed')
        end subroutine test_image_local

    end subroutine test_image

end module simple_image

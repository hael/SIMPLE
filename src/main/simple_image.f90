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
    ! CONSTRUCTORS
    procedure          :: new
    procedure          :: set_wthreads
    procedure          :: construct_thread_safe_tmp_imgs
    procedure, private :: disc_1, disc_2
    generic            :: disc => disc_1, disc_2
    procedure          :: ring, soft_ring
    procedure          :: disc_sideview
    procedure          :: cylinder
    procedure          :: copy
    procedure          :: copy_fast
    procedure          :: img2spec
    procedure          :: mic2spec
    procedure          :: pspec_graphene_mask
    procedure          :: dampen_pspec_central_cross
    procedure          :: scale_pspec4viz
    procedure          :: window
    procedure          :: window_slim
    procedure          :: window_center
    procedure          :: add_window
    procedure          :: win2arr
    procedure          :: win2arr_rad
    procedure          :: corner
    ! I/O
    procedure, private :: open
    procedure          :: read
    procedure          :: write
    procedure          :: update_header_stats
    procedure          :: write_jpg
    ! GETTERS/SETTERS
    procedure          :: get_array_shape
    procedure          :: get_ldim, get_box
    procedure          :: get_smpd
    procedure          :: get_nyq
    procedure          :: get_filtsz
    procedure          :: get_shconst
    procedure          :: get
    procedure          :: get_rmat
    procedure          :: get_mat_ptrs
    procedure          :: get_rmat_ptr
    procedure          :: get_rmat_sub
    procedure          :: get_cmat
    procedure          :: get_cmat_ptr
    procedure          :: get_cmat_sub
    procedure, private :: get_cmat_at_1, get_cmat_at_2
    generic            :: get_cmat_at => get_cmat_at_1, get_cmat_at_2
    procedure, private :: get_rmat_at_1, get_rmat_at_2
    generic            :: get_rmat_at => get_rmat_at_1, get_rmat_at_2
    procedure          :: get_avg_int
    procedure          :: get_sum_int
    procedure          :: set_rmat_at
    procedure, private :: set_1, set_2
    generic            :: set => set_1, set_2
    procedure          :: set_rmat
    procedure, private :: set_cmat_1, set_cmat_2, set_cmat_3, set_cmat_4
    generic            :: set_cmat => set_cmat_1, set_cmat_2, set_cmat_3, set_cmat_4
    procedure, private :: set_cmat_at_1, set_cmat_at_2
    generic            :: set_cmat_at => set_cmat_at_1, set_cmat_at_2
    procedure          :: set_cmats_from_cmats
    procedure          :: add_cmats_to_cmats
    procedure          :: print_cmat
    procedure          :: print_rmat
    procedure          :: expand_ft
    procedure          :: set_smpd
    procedure          :: get_slice
    procedure          :: set_slice
    procedure          :: get_subimg
    procedure          :: get_lfny
    procedure          :: get_lhp
    procedure          :: get_lp
    procedure          :: get_spat_freq
    procedure          :: get_find
    procedure          :: rmat_associated
    procedure          :: cmat_associated
    procedure          :: is_wthreads
    procedure, private :: serialize_1, serialize_2, serialize_3
    generic            :: serialize => serialize_1, serialize_2, serialize_3
    procedure          :: unserialize
    procedure          :: winserialize
    procedure          :: zero2one
    procedure          :: get_fcomp
    procedure          :: get_fcomp2D
    procedure          :: set_fcomp
    procedure          :: vis
    procedure          :: set_ft
    ! CHECKUPS
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
    ! ARITHMETICS
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
    procedure          :: ctf_dens_correct
    procedure          :: ctf_dens_correct_wiener
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
    ! BINARY IMAGE METHODS (heavy binary lifters deferred to binimage that extends this class)
    procedure          :: nforeground
    procedure          :: nbackground
    procedure, private :: binarize_1, binarize_2, binarize_3
    generic            :: binarize => binarize_1, binarize_2, binarize_3
    procedure          :: cendist
    procedure          :: masscen
    procedure          :: masscen_adjusted
    procedure          :: box_cen_arg
    procedure          :: calc_shiftcen, calc_shiftcen_serial
    procedure          :: bin_inv
    procedure          :: remove_edge
    procedure          :: one_at_edge
    procedure          :: bin2logical
    procedure          :: logical2bin
    procedure          :: density_inoutside
    procedure          :: collage
    procedure          :: tile
    procedure          :: generate_orthogonal_reprojs
    ! FILTERS
    procedure          :: acf
    procedure          :: ccf
    procedure          :: guinier_bfac
    procedure          :: guinier
    procedure          :: spectrum
    procedure          :: power_spectrum
    procedure          :: ran_phases_below_noise_power
    procedure          :: whiten_noise_power
    procedure          :: fcomps_below_noise_power_stats
    procedure          :: apply_bfac
    procedure          :: bp
    procedure          :: lp
    procedure          :: bpgau2D
    procedure          :: tophat
    procedure, private :: apply_filter_1, apply_filter_2
    generic            :: apply_filter => apply_filter_1, apply_filter_2
    procedure          :: apply_filter_serial
    procedure, private :: imfilter1, imfilter2, imfilter3
    generic            :: imfilter => imfilter1, imfilter2, imfilter3
    procedure          :: phase_rand
    procedure          :: hannw
    procedure          :: real_space_filter
    procedure          :: NLmean2D, NLmean2D_eo, NLmean3D, NLmean3D_eo
    procedure          :: ICM2D, ICM2D_eo, ICM3D, ICM3D_eo
    procedure          :: GLCM
    ! CALCULATORS
    procedure          :: minmax
    procedure          :: loc_sdev
    procedure          :: avg_loc_sdev
    procedure          :: loc_var, loc_var3D
    procedure          :: rmsd
    procedure, private :: stats_1, stats_2
    generic            :: stats => stats_1, stats_2
    procedure          :: skew, kurt
    procedure          :: noisesdev
    procedure          :: mean
    procedure          :: contains_nans
    procedure          :: checkimg4nans
    procedure          :: cure
    procedure          :: loop_lims
    procedure          :: calc_gradient
    procedure          :: gradients_magnitude
    procedure          :: calc_ice_score
    procedure          :: gradient
    procedure, private :: comp_addr_phys1, comp_addr_phys2, comp_addr_phys3
    generic            :: comp_addr_phys =>  comp_addr_phys1, comp_addr_phys2, comp_addr_phys3
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
    procedure          :: sqeuclid
    procedure, private :: sqeuclid_matrix_1, sqeuclid_matrix_2
    generic            :: sqeuclid_matrix => sqeuclid_matrix_1, sqeuclid_matrix_2
    procedure          :: euclid_norm
    procedure          :: opt_filter_costfun
    procedure          :: opt_filter_costfun_workshare
    procedure          :: fsc, fsc_scaled
    procedure          :: get_res
    procedure, private :: oshift_1, oshift_2
    generic            :: oshift => oshift_1, oshift_2
    procedure, private :: gen_argtransf_comp
    procedure          :: radial_cc
    procedure          :: calc_principal_axes_rotmat
    ! MODIFIERS
    procedure          :: lp_background
    procedure          :: combine_fgbg_filt
    procedure          :: insert
    procedure          :: insert_lowres
    procedure          :: insert_lowres_serial
    procedure          :: inv
    procedure          :: zero_neg
    procedure          :: remove_neg
    procedure          :: ran
    procedure          :: gauran
    procedure          :: add_gauran
    procedure          :: dead_hot_positions
    procedure          :: taper_edges, taper_edges_hann
    procedure          :: subtr_backgr_ramp
    procedure          :: subtract_background, estimate_background, upsample_square_background
    procedure          :: zero
    procedure          :: zero_and_unflag_ft
    procedure          :: zero_and_flag_ft
    procedure          :: zero_background
    procedure          :: zero_env_background
    procedure          :: pad_fft
    procedure          :: norm_noise_pad_fft
    procedure          :: div_w_instrfun
    procedure          :: salt_n_pepper
    procedure          :: square
    procedure          :: corners
    procedure          :: before_after
    procedure, private :: gauimg_1, gauimg_2
    generic            :: gauimg => gauimg_1, gauimg_2
    procedure          :: gauimg2D
    procedure          :: gauimg3D
    procedure          :: subtr_backgr
    procedure          :: resmsk
    procedure          :: frc_pspec
    procedure          :: mask
    procedure          :: neg
    procedure          :: pad
    procedure          :: pad_inplace
    procedure, private :: pad_mirr_1, pad_mirr_2
    generic            :: pad_mirr => pad_mirr_1, pad_mirr_2
    procedure          :: clip
    procedure          :: clip_inplace
    procedure          :: read_and_crop
    procedure          :: scale_pixels
    procedure          :: mirror
    procedure          :: norm
    procedure          :: variance
    procedure          :: norm_minmax
    procedure          :: norm4viz
    procedure          :: norm_ext
    procedure          :: norm_noise
    procedure          :: norm_within
    procedure          :: calc_bin_thres
    procedure          :: quantize_fwd
    procedure          :: quantize_bwd
    procedure          :: zero_edgeavg
    procedure          :: roavg
    procedure          :: rtsq
    procedure          :: rtsq_serial
    procedure          :: shift_phorig
    procedure          :: shift
    procedure, private :: shift2Dserial_1, shift2Dserial_2
    generic            :: shift2Dserial => shift2Dserial_1, shift2Dserial_2
    procedure          :: set_within
    procedure          :: ft2img
    procedure          :: img2ft
    procedure          :: cure_outliers
    procedure          :: zero_below
    procedure          :: div_below
    procedure          :: ellipse
    procedure          :: flip
    procedure          :: reshape2cube
    ! FFTs
    procedure          :: fft  => fwd_ft
    procedure          :: ifft => bwd_ft
    procedure          :: fft_noshift
    ! DESTRUCTOR
    procedure :: kill
    procedure :: kill_thread_safe_tmp_imgs
end type image

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

    !>  \brief  Constructor for simple_image class
    subroutine new( self, ldim, smpd, wthreads )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: wthreads
        integer(kind=c_int) :: rc
        integer             :: i
        logical             :: do_allocate
        ! we need to be clever about allocation (because it is costly)
        if( self%existence )then
            if( any(self%ldim /= ldim) )then
                do_allocate = .true.
                call self%kill()
            else
                do_allocate = .false.
                !$omp critical
                call fftwf_destroy_plan(self%plan_fwd)
                call fftwf_destroy_plan(self%plan_bwd)
                !$omp end critical
            endif
        else
            do_allocate = .true.
        endif
        self%wthreads = .true.
        if( present(wthreads) ) self%wthreads = wthreads
        self%wthreads = self%wthreads .and. nthr_glob > 1
        self%ldim = ldim
        self%smpd = smpd
        ! Make Fourier iterator
        call self%fit%new(ldim, smpd)
        ! Work out dimensions of the complex array
        self%array_shape(1)   = fdim(self%ldim(1))
        self%array_shape(2:3) = self%ldim(2:3)
        self%nc = int(product(self%array_shape)) ! # components
        if( do_allocate )then
            ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
            self%p = fftwf_alloc_complex(int(product(self%array_shape),c_size_t))
            ! Set up the complex array which will point at the allocated memory
            call c_f_pointer(self%p,self%cmat,self%array_shape)
            ! Work out the shape of the real array
            self%array_shape(1) = 2*(self%array_shape(1))
            ! Set up the real array
            call c_f_pointer(self%p,self%rmat,self%array_shape)
        endif
        ! put back the shape of the complex array
        self%array_shape(1) = fdim(self%ldim(1))
        ! init
        self%rmat = 0.
        self%ft   = .false.
        !$omp critical
        ! make fftw plans
        if( self%wthreads .and. (any(ldim >= 200) .or. ldim(3) >= 100) )then
            rc = fftwf_init_threads()
            call fftwf_plan_with_nthreads(nthr_glob)
        endif
        if(self%ldim(3) > 1)then
            self%plan_fwd = fftwf_plan_dft_r2c_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        else
            self%plan_fwd = fftwf_plan_dft_r2c_2d(self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_2d(self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        endif
        if( self%wthreads .and. (any(ldim >= 200) .or. ldim(3) >= 100) )then
            ! disable threads for subsequent plans
            call fftwf_plan_with_nthreads(1)
        endif
        !$omp end critical
        ! set shift constant (shconst)
        do i=1,3
            if( self%ldim(i) == 1 )then
                self%shconst(i) = 0.
                cycle
            endif
            if( is_even(self%ldim(i)) )then
                self%shconst(i) = PI/real(self%ldim(i)/2.)
            else
                self%shconst(i) = PI/real((self%ldim(i)-1)/2.)
            endif
        end do
        self%existence = .true.
    end subroutine new

    subroutine set_wthreads( self, wthreads )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: wthreads
        integer(kind=c_int) :: rc
        if( self%wthreads .eqv. wthreads ) return
        self%wthreads = wthreads
        self%wthreads = self%wthreads .and. nthr_glob > 1
        !$omp critical
        call fftwf_destroy_plan(self%plan_fwd)
        call fftwf_destroy_plan(self%plan_bwd)
        ! make fftw plans
        if( self%wthreads )then
            rc = fftwf_init_threads()
            call fftwf_plan_with_nthreads(nthr_glob)
        endif
        if( self%ldim(3) > 1 )then
            self%plan_fwd = fftwf_plan_dft_r2c_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        else
            self%plan_fwd = fftwf_plan_dft_r2c_2d(self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_2d(self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        endif
        if( self%wthreads )then
            ! disable threads for subsequent plans
            call fftwf_plan_with_nthreads(1)
        endif
        !$omp end critical
    end subroutine set_wthreads

    subroutine construct_thread_safe_tmp_imgs( self, nthr )
        class(image), intent(in) :: self
        integer,      intent(in) :: nthr
        integer :: i, sz, ldim(3)
        logical :: do_allocate
        if( allocated(thread_safe_tmp_imgs) )then
            ldim = thread_safe_tmp_imgs(1)%get_ldim()
            sz   = size(thread_safe_tmp_imgs)
            if( any(self%ldim /= ldim) .or. sz /= nthr )then
                do i=1,size(thread_safe_tmp_imgs)
                    call thread_safe_tmp_imgs(i)%kill
                end do
                deallocate(thread_safe_tmp_imgs)
                do_allocate = .true.
            else
                do_allocate = .false.
            endif
        else
            do_allocate = .true.
        endif
        if( do_allocate )then
            allocate( thread_safe_tmp_imgs(nthr) )
            do i=1,nthr
                call thread_safe_tmp_imgs(i)%new(self%ldim, self%smpd, .false.)
            end do
        endif
    end subroutine construct_thread_safe_tmp_imgs

    !>  \brief disc constructs a binary disc of given radius and returns the number of 1:s
    subroutine disc_1( self, ldim, smpd, radius, npix )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, radius
        integer, optional, intent(inout) :: npix
        call self%new(ldim, smpd)
        call self%cendist
        where(self%rmat <= radius)
            self%rmat = 1.
        else where
            self%rmat = 0.
        end where
        if( present(npix) ) npix = count(self%rmat>0.5)
    end subroutine disc_1

    subroutine disc_2( self, ldim, smpd, radius, lmsk, npix )
        class(image),         intent(inout) :: self
        integer,              intent(in)    :: ldim(3)
        real,                 intent(in)    :: smpd, radius
        logical, allocatable, intent(out)   :: lmsk(:,:,:)
        integer, optional,    intent(inout) :: npix
        call self%new(ldim, smpd)
        call self%cendist
        if( allocated(lmsk) ) deallocate(lmsk)
        allocate(lmsk(ldim(1),ldim(2),ldim(3)))
        where(self%rmat(:ldim(1),:ldim(2),:ldim(3)) <= radius)
            self%rmat(:ldim(1),:ldim(2),:ldim(3)) = 1.
            lmsk = .true.
        else where
            self%rmat(:ldim(1),:ldim(2),:ldim(3)) = 0.
            lmsk = .false.
        end where
        if( present(npix) ) npix = count(lmsk)
    end subroutine disc_2

    !>  \brief ring constructs a binary ring and returns the number of 1:s
    subroutine ring( self, ldim, smpd, outer_radius, inner_radius, npix )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, outer_radius, inner_radius
        integer, optional, intent(inout) :: npix
        call self%new(ldim, smpd)
        call self%cendist
        where(self%rmat <= outer_radius .and. self%rmat >= inner_radius )
            self%rmat = 1.
        else where
            self%rmat = 0.
        end where
        if( present(npix) )npix = count(self%rmat>0.5)
    end subroutine ring

    !>  \brief soft ring based on gamma distribution
    subroutine soft_ring( self, ldim, smpd, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        real,              intent(in)    :: radius ! in pixels
        real, parameter :: K     = 2.0
        real, parameter :: THETA = 2.0
        real    :: d,km1,mode,val,scale
        integer :: c(3),i,j,l
        call self%new(ldim, smpd)
        c     = nint(real(self%ldim)/2.)+1
        km1   = K-1.
        mode  = km1 * THETA
        val   = mode**km1 * exp(-mode/THETA)
        scale = 1./val
        do l=1,self%ldim(3)
        do j=1,self%ldim(2)
        do i=1,self%ldim(1)
            d = hyp(i-c(1),j-c(2),l-c(3))
            if( d > radius )then
                val = 0.
            else
                d   = 20. * (radius - d) / radius
                val = max(0.,min(1.0,scale * (d**km1 * exp(-d/THETA))))
            endif
            self%rmat(i,j,l) = val
        enddo
        enddo
        enddo
    end subroutine soft_ring

    ! is an empirical reprojection of the side view of a thick disc
    subroutine disc_sideview( self, ldim, smpd, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        real,              intent(in)    :: radius ! in pixels
        integer, parameter :: B = 3
        real    :: t,d
        integer :: c(3),i,j,h,f
        if( ldim(3) > 1 ) THROW_HARD('2D only; membrane')
        call self%new(ldim, smpd)
        c = nint(real(self%ldim)/2.)+1
        ! central line
        j = c(2)
        do i = 1,self%ldim(1)
            h = i-c(1)
            d = real(abs(h))
            if( d >= radius )then
                if( d > radius+real(B) )then
                    self%rmat(i,j,1) = 0.
                else
                    ! horizontal soft cosine edge
                    d = (d-radius)/real(B)
                    self%rmat(i,j,1) = 2.*sqrt(radius-0.25)*max(0.,cos(d*PIO2))
                endif
            else
                self%rmat(i,j,1) = self%rmat(i,j,1) + 2.*sqrt(radius**2-real(h**2))
            endif
        enddo
        ! thickness taken as ~40Angs, soft edge excluded
        t = 40./self%smpd
        f = floor(t/2.)
        do j = c(2)-f-B,c(2)+f+B
            if( j == c(2) ) cycle
            if( (j >= c(2)-f) .and. (j <= c(2)+f) )then
                self%rmat(:,j,1) = self%rmat(:,c(2),1)
            else
                ! vertical soft edge
                d = (real( abs(j-c(2)) - f)) / real(B)
                self%rmat(:,j,1) = self%rmat(:,c(2),1) * max(0.,cos(d*PIO2))
            endif
        enddo
        self%rmat = self%rmat / maxval(self%rmat(1:self%ldim(1),c(2),1))
    end subroutine disc_sideview

    !>  \brief constructs a binary cylinder along the z-axis
    subroutine cylinder( self, ldim, smpd, height, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, height, radius
        real    :: rsq,radiussq
        integer :: icenter(3),i,j,k
        call self%new(ldim, smpd)
        icenter  = nint(real(self%ldim/2.))+1
        radiussq = radius**2
        do k = 1,self%ldim(3)
            if( real(abs(k-icenter(3))) > height/2. ) cycle
            do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                rsq = (i-icenter(1))**2 + (j-icenter(2))**2
                if( rsq < radiussq ) self%rmat(i,j,k) = 1
            enddo
            enddo
        enddo
    end subroutine cylinder

    subroutine copy( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
        call self%new(self_in%ldim, self_in%smpd, self%wthreads)
        self%rmat = self_in%rmat
        self%ft   = self_in%ft
    end subroutine copy

    subroutine copy_fast( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
        if( self%wthreads )then
            !$omp parallel workshare
            self%rmat = self_in%rmat
            !$omp end parallel workshare
        else
            self%rmat = self_in%rmat
        endif
        self%ft = self_in%ft
    end subroutine copy_fast

    !> img2spec calculates the powerspectrum of the input image
    !!          the resulting spectrum has dampened central cross and subtracted background
    subroutine img2spec( self, speckind, lp_backgr_subtr, img_out, postproc )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        logical, optional, intent(in)    :: postproc
        type(image) :: tmp, tmp2
        logical     :: didft, l_postproc
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images')
        l_postproc = .true.
        if( present(postproc) ) l_postproc = postproc
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        call img_out%new(self%ldim, self%smpd)
        call tmp%copy(self)
        call tmp%norm()
        call tmp%zero_edgeavg
        call tmp%fft()
        call tmp%ft2img(speckind, img_out)
        if( l_postproc )then
            call img_out%dampen_pspec_central_cross
            call img_out%subtr_backgr(lp_backgr_subtr)
        endif
        if( didft ) call self%fft()
        call tmp%kill
        call tmp2%kill
    end subroutine img2spec

    !> mic2spec calculates the average powerspectrum over a micrograph
    !!          the resulting spectrum has dampened central cross and subtracted background
    subroutine mic2spec( self, box, speckind, lp_backgr_subtr, img_out, postproc )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: box
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        logical, optional, intent(in)    :: postproc
        type(image) :: tmp, tmp2
        integer     :: xind, yind, cnt
        logical     :: didft, outside, l_postproc
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            THROW_HARD('cannot use a box larger than the image')
        endif
        l_postproc = .true.
        if( present(postproc) ) l_postproc = postproc
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        call img_out%new([box,box,1], self%smpd)
        call tmp%new([box,box,1], self%smpd)
        call tmp2%new([box,box,1], self%smpd)
        cnt = 0
        do xind=0,self%ldim(1)-box,box/2
            do yind=0,self%ldim(2)-box,box/2
                call self%window_slim([xind,yind],box,tmp,outside)
                call tmp%norm()
                call tmp%zero_edgeavg
                call tmp%fft()
                call tmp%ft2img(speckind, tmp2)
                call img_out%add(tmp2)
                cnt = cnt+1
                call tmp%zero_and_unflag_ft
                call tmp2%zero_and_unflag_ft
            end do
        end do
        call img_out%div(real(cnt))
        if( l_postproc )then
            call img_out%dampen_pspec_central_cross
            call img_out%subtr_backgr(lp_backgr_subtr)
        endif
        if( didft ) call self%fft()
        call tmp%kill
        call tmp2%kill
    end subroutine mic2spec

    subroutine pspec_graphene_mask( self, ldim, smpd )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        real,         intent(in)    :: smpd
        logical, allocatable :: graphene_mask(:)
        type(image) :: tmp
        integer     :: h, k, l, lims(3,2), phys(3), sh, lfny
        call self%new(ldim, smpd)
        self%ft = .true.
        call tmp%new(ldim, smpd)
        graphene_mask = calc_graphene_mask(ldim(1), self%smpd)
        lims = self%fit%loop_lims(2)
        lfny = self%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > lfny ) cycle
                    phys = self%fit%comp_addr_phys([h,k,l])
                    if( graphene_mask(sh) )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(1.,0.)
                    else
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
        call self%ft2img('real', tmp)
        call self%copy(tmp)
        call tmp%kill
    end subroutine pspec_graphene_mask

    !> \brief dampens the central cross of a powerspectrum by mean filtering
    subroutine dampen_pspec_central_cross( self )
        class(image), intent(inout) :: self
        integer, parameter :: HDAMPWINSZ=2
        integer, parameter :: DAMPWINSZ =5
        real :: medi(self%ldim(1)),medj(self%ldim(2)),vals(DAMPWINSZ**2)
        integer :: h,mh,k,mk,lims(3,2),i,j,l,r,u,d,n
        if( self%ft )          THROW_HARD('not intended for FTs; dampen_pspec_central_cross')
        if( self%ldim(3) > 1 ) THROW_HARD('not intended for 3D imgs; dampen_pspec_central_cross')
        lims = self%loop_lims(3)
        mh   = abs(lims(1,1))
        mk   = abs(lims(2,1))
        n    = DAMPWINSZ*DAMPWINSZ
        ! along h
        h = 0
        i = min(max(1,h+mh+1),self%ldim(1))
        l = min(max(1,h-HDAMPWINSZ+mh+1),self%ldim(1))
        r = min(max(1,h+HDAMPWINSZ+mh+1),self%ldim(1))
        do j=1,self%ldim(2)
            d = max(1,j-HDAMPWINSZ)
            u = min(d+DAMPWINSZ-1,self%ldim(2))
            d = u-DAMPWINSZ+1
            vals = reshape(self%rmat(l:r,d:u,1),(/n/))
            medj(j) = median_nocopy(vals)
        enddo
        ! along k
        k = 0
        j = min(max(1,k+mk+1),self%ldim(2))
        d = min(max(1,k-HDAMPWINSZ+mk+1),self%ldim(2))
        u = min(max(1,k+HDAMPWINSZ+mk+1),self%ldim(2))
        do i=1,self%ldim(1)
            l = max(1,i-HDAMPWINSZ)
            r = min(l+DAMPWINSZ-1,self%ldim(1))
            l = r-DAMPWINSZ+1
            vals = reshape(self%rmat(l:r,d:u,1),(/n/))
            medi(i) = median_nocopy(vals)
        enddo
        ! replace
        h = 0
        i = min(max(1,h+mh+1),self%ldim(1))
        self%rmat(i,1:self%ldim(2),1) = medj
        k = 0
        j = min(max(1,k+mk+1),self%ldim(2))
        self%rmat(1:self%ldim(1),j,1) = medi
    end subroutine dampen_pspec_central_cross

    subroutine scale_pspec4viz( self, rsmpd4viz )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: rsmpd4viz
        type(image) :: tmp
        real        :: scale4viz, smpd4viz_here
        integer     :: box4clip
        if( self%ft )          THROW_HARD('pspec input assumed to be in real-space; scale_pspec4viz')
        if( self%ldim(3) > 1 ) THROW_HARD('pspec input assumed to be 2D; scale_pspec4viz')
        smpd4viz_here = SMPD4VIZ
        if( present(rsmpd4viz) ) smpd4viz_here = rsmpd4viz
        scale4viz = min(self%smpd / smpd4viz_here, 1.)
        if( scale4viz < 1. )then
            box4clip = round2even(scale4viz * real(self%ldim(1)))
        else
            return
        endif
        call tmp%new([box4clip,box4clip,1], smpd4viz_here)
        call self%clip(tmp)
        call tmp%zero_edgeavg
        call tmp%fft
        call tmp%pad(self)
        call self%ifft
        call tmp%kill
    end subroutine scale_pspec4viz

    !>  \brief window extracts a particle image from a box as defined by EMAN 1.9
    subroutine window( self_in, coord, box, self_out, noutside )
        class(image),      intent(in)    :: self_in
        integer,           intent(in)    :: coord(2), box
        class(image),      intent(inout) :: self_out
        integer, optional, intent(inout) :: noutside
        integer :: i,j, fromc(2), toc(2), xoshoot, yoshoot, xushoot, yushoot, xboxrange(2), yboxrange(2)
        if( self_in%ldim(3) > 1 ) THROW_HARD('only 4 2D images; window')
        if( self_in%is_ft() )     THROW_HARD('only 4 real images; window')
        if( self_out%exists() )then
            if( self_out%is_ft() ) THROW_HARD('only 4 real images; window')
            if( self_out%ldim(1) == box .and. self_out%ldim(2) == box .and. self_out%ldim(3) == 1 )then
                ! go ahead
            else
                call self_out%new([box,box,1], self_in%smpd)
            endif
        else
            call self_out%new([box,box,1], self_in%smpd)
        endif
        fromc = coord+1       ! compensate for the c-range that starts at 0
        toc   = fromc+(box-1) ! the lower left corner is 1,1
        if( any(fromc < 1) .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) )then
            if( present(noutside) )then
                noutside = noutside+1
            else
                THROW_WARN('box extends outside micrograph; window')
            endif
        endif
        xoshoot = 0
        yoshoot = 0
        xushoot = 0
        yushoot = 0
        if( toc(1)   > self_in%ldim(1) ) xoshoot =  toc(1)   - self_in%ldim(1)
        if( toc(2)   > self_in%ldim(2) ) yoshoot =  toc(2)   - self_in%ldim(2)
        if( fromc(1) < 1               ) xushoot = -fromc(1) + 1
        if( fromc(2) < 1               ) yushoot = -fromc(2) + 1
        toc(1)        = toc(1)   - xoshoot
        toc(2)        = toc(2)   - yoshoot
        fromc(1)      = fromc(1) + xushoot
        fromc(2)      = fromc(2) + yushoot
        xboxrange(1)  = xushoot  + 1
        xboxrange(2)  = box      - xoshoot
        yboxrange(1)  = yushoot  + 1
        yboxrange(2)  = box      - yoshoot
        self_out%rmat = 0.
        self_out%rmat(xboxrange(1):xboxrange(2),yboxrange(1):yboxrange(2),1) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),1)
        ! Rather than stretching the last inside line that creates stripes
        ! it is scrambled to mitigate Gibbs phenomenon and preserve image stats
        if( xboxrange(1) > 1 )then
            do i = 1,xboxrange(1)-1
                do j = 1,self_out%ldim(2)
                    self_out%rmat(i,j,1) = self_out%rmat(xboxrange(1),irnd_uni(self_out%ldim(2)),1)
                enddo
            enddo
        endif
        if( xboxrange(2) < self_out%ldim(1) )then
            do i = xboxrange(2)+1,self_out%ldim(1)
                do j = 1,self_out%ldim(2)
                    self_out%rmat(i,j,1) = self_out%rmat(xboxrange(2),irnd_uni(self_out%ldim(2)),1)
                enddo
            enddo
        endif
        if( yboxrange(1) > 1 )then
            do j = 1,yboxrange(1)-1
                do i = 1,self_out%ldim(1)
                    self_out%rmat(i,j,1) = self_out%rmat(irnd_uni(self_out%ldim(1)),yboxrange(1),1)
                enddo
            enddo
        endif
        if( yboxrange(2) < self_out%ldim(2) )then
            do j = yboxrange(2)+1,self_out%ldim(2)
                do i = 1,self_out%ldim(1)
                    self_out%rmat(i,j,1) = self_out%rmat(irnd_uni(self_out%ldim(1)),yboxrange(2),1)
                enddo
            enddo
        endif
    end subroutine window

    !>  window_slim  extracts a particle image from a box as defined by EMAN 1.9
    subroutine window_slim( self_in, coord, box, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: coord(:), box
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
        logical                     :: isvol
        integer, allocatable        :: fromc(:), toc(:)
        isvol         = self_in%is_3d()
        if( isvol )then
            allocate(fromc(3), toc(3))
        else
            allocate(fromc(2), toc(2))
        endif
        fromc         = coord + 1         ! compensate for the c-range that starts at 0
        toc           = fromc + (box - 1) ! the lower left corner is 1,1
        self_out%rmat = 0.
        self_out%ft   = .false.
        outside       = .false.
        if( isvol )then
            if( size(coord) /= 3 ) THROW_HARD("Error! expecting 3D coordinates; window_slim")
            if( fromc(1) < 1 .or. fromc(2) < 1 .or. fromc(3) < 1 .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) .or. toc(3) > self_in%ldim(3) )then
                outside = .true.
            else
                self_out%rmat(1:box,1:box,1:box) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),fromc(3):toc(3))
            endif
        else
            if( size(coord) /= 2 ) THROW_HARD("Error! expecting 2D coordinates; window_slim")
            if( fromc(1) < 1 .or. fromc(2) < 1 .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) )then
                outside = .true.
            else
                self_out%rmat(1:box,1:box,1) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),1)
            endif
        endif
    end subroutine window_slim

    subroutine window_center( self_in, center, rad, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: center(:), rad
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
        integer :: coord(3), box
        box = rad*2
        if( self_in%is_3d() )then
            coord = center - rad
            if( any(coord < 1) )      THROW_HARD('Error! window is out of the image')
            call window_slim(self_in, coord, box, self_out, outside)
        else
            coord(1:2) = center - rad
            if( any(coord(1:2) < 1) ) THROW_HARD('Error! window is out of the image')
            call window_slim(self_in, coord(1:2), box, self_out, outside)
        endif
    end subroutine window_center

    ! for re-generation of micrograph after convolutional PPCA
    subroutine add_window( self, imgwin, coord, offset )
        class(image), intent(inout)   :: self
        class(image), intent(in)      :: imgwin
        integer,      intent(in)      :: coord(2)
        integer, optional, intent(in) :: offset  !if offset is present it doesn't sum windows, but take just the inner part
        integer ::  fromc(2), toc(2), ld(3), box
        ld    = imgwin%get_ldim()
        box   = ld(1)
        fromc = coord + 1         ! compensate for the c-range that starts at 0
        toc   = fromc + (box - 1) ! the lower left corner is 1,1
        if( fromc(1) < 1 .or. fromc(2) < 1 .or. toc(1) > self%ldim(1) .or. toc(2) > self%ldim(2) )then
          return
        endif
        if(present(offset)) then  !no sovrapposition
            if(offset == int(box/2)) THROW_HARD("invalid offset choice; add_window") ! box is supposet to be even
            if (offset > int(box/2)) then
                self%rmat(fromc(1)+box-offset-1:offset,fromc(2)+box-offset-1:offset,1) = &
                    &imgwin%rmat(box-offset:offset,box-offset:offset,1)
            else
                self%rmat(fromc(1)+offset-1:box-offset,fromc(2)+offset-1:box-offset,1) = &
                    &imgwin%rmat(offset:box-offset,offset:box-offset,1)
            endif
        else
            self%rmat(fromc(1):toc(1),fromc(2):toc(2),1) = self%rmat(fromc(1):toc(1),fromc(2):toc(2),1) &
                &+ imgwin%rmat(1:box,1:box,1) !add everything
        endif
    end subroutine add_window

    !>  \brief win2arr extracts a small window into an array (circular indexing)
    function win2arr( self, i, j, k, winsz ) result( pixels )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i, j, k, winsz
        real, allocatable :: pixels(:)
        integer :: s, ss, t, tt, u, uu, cnt, npix
        if( self%is_ft() ) THROW_HARD('only 4 real images; win2arr')
        if( self%is_3d() )then
            npix = (2*winsz+1)**3
        else
            npix = (2*winsz+1)**2
        endif
        allocate(pixels(npix))
        cnt = 1
        do s=i-winsz,i+winsz
            ss = cyci_1d_static(self%ldim(1), s)
            do t=j-winsz,j+winsz
                tt = cyci_1d_static(self%ldim(2), t)
                if( self%ldim(3) > 1 )then
                    do u=k-winsz,k+winsz
                        uu          = cyci_1d_static(self%ldim(3), u)
                        pixels(cnt) = self%rmat(ss,tt,uu)
                        cnt         = cnt+1
                    end do
                else
                    pixels(cnt) = self%rmat(ss,tt,1)
                    cnt         = cnt+1
                endif
            end do
        end do
    end function win2arr

    !>  \brief win2arr_rad extracts a small window into an array (circular indexing) with radial thres
    subroutine win2arr_rad( self, i, j, k, winsz, npix_in, maxrad, npix_out, pixels )
        class(image), intent(in)    :: self
        integer,      intent(in)    :: i, j, k, winsz, npix_in
        real,         intent(in)    :: maxrad ! in pixels
        integer,      intent(inout) :: npix_out
        real,         intent(inout) :: pixels(npix_in)
        integer :: s, ss, t, tt, u, uu
        real    :: vec_cen(3), vec_displ(3), maxradsq
        if( self%is_ft() ) THROW_HARD('only 4 real images; win2arr_rad')
        if( self%is_3d() )then
            vec_cen = real([i,j,k])
        else
            vec_cen = real([i,j,1])
        endif
        maxradsq = maxrad**2.
        npix_out = 1
        do s = i - winsz,i + winsz
            ss = cyci_1d_static(self%ldim(1), s)
            do t = j - winsz, j + winsz
                tt = cyci_1d_static(self%ldim(2), t)
                if( self%ldim(3) > 1 )then
                    do u = k - winsz, k + winsz
                        vec_displ = real([s,t,u]) - vec_cen
                        if( sum(vec_displ**2.) <= maxradsq )then
                            uu               = cyci_1d_static(self%ldim(3), u)
                            pixels(npix_out) = self%rmat(ss,tt,uu)
                            npix_out         = npix_out + 1
                        endif
                    end do
                else
                    vec_displ = real([s,t,1]) - vec_cen
                    if( sum(vec_displ(1:2))**2. <= maxradsq )then
                        pixels(npix_out) = self%rmat(ss,tt,1)
                        npix_out         = npix_out + 1
                    endif
                endif
            end do
        end do
    end subroutine win2arr_rad

    !>  \brief corner extracts a corner of a volume with size box
    subroutine corner( self_in, box, self_out )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: box
        type(image),  intent(inout) :: self_out
        if( self_in%ldim(3) <= 1 ) THROW_HARD('only 4 3D images; corner')
        if( self_in%is_ft() )      THROW_HARD('only 4 real images; corner')
        call self_out%new([box,box,box], self_in%smpd)
        self_out%rmat(:box,:box,:box) = self_in%rmat(:box,:box,:box)
    end subroutine corner

    ! I/O

    !>  \brief  open: for reading 2D images from stack or volumes from volume files
    subroutine open( self, fname, ioimg, formatchar, readhead, rwaction )
        class(image),               intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        class(imgfile),             intent(inout) :: ioimg
        character(len=1), optional, intent(in)    :: formatchar
        logical,          optional, intent(in)    :: readhead
        character(len=*), optional, intent(in)    :: rwaction
        character(len=1) :: form
        integer          :: mode
        if( self%existence )then
            if( .not. file_exists(fname) )then
                write(logfhandle,*) 'file: ', trim(fname)
                THROW_HARD('The file you are trying to open does not exists; open')
            endif
            if( present(formatchar) )then
                form = formatchar
            else
                form = fname2format(fname)
            endif
            self%ft = .false.
            select case(form)
                case('M')
                    call ioimg%open(fname, self%ldim, self%smpd, formatchar=formatchar, readhead=readhead, rwaction=rwaction)
                    ! data type: 0 image: signed 8-bit bytes range -128 to 127
                    !            1 image: 16-bit halfwords
                    !            2 image: 32-bit reals (DEFAULT MODE)
                    !            3 transform: complex 16-bit integers
                    !            4 transform: complex 32-bit reals (THIS WOULD BE THE DEFAULT FT MODE)
                    mode = ioimg%getMode()
                    if( mode == 3 .or. mode == 4 ) self%ft = .true.
                case('F','S','J','L')
                    call ioimg%open(fname, self%ldim, self%smpd, formatchar=formatchar, readhead=readhead, rwaction=rwaction)
            end select
        else
            THROW_HARD('image need to be constructed before read/write; open')
        endif
    end subroutine open

    !>  \brief read: for reading 2D images from stack or volumes from volume files
    subroutine read( self, fname, i, readhead )
        class(image),               intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        integer,          optional, intent(in)    :: i
        logical,          optional, intent(in)    :: readhead
        type(imgfile)    :: ioimg
        character(len=1) :: form
        integer          :: ldim(3), first_slice
        integer          :: last_slice, ii
        real             :: smpd
        logical          :: isvol
        ldim  = self%ldim
        smpd  = self%smpd
        isvol = .true. ! assume volume by default
        ii    = 1      ! default location
        if( present(i) )then
            ! we are reading from a stack & in SIMPLE volumes are not allowed
            ! to be stacked so the image object must be 2D
            isvol = .false.
            ii = i ! replace default location
        endif
        form = fname2format(fname)
        select case(form)
            case('M', 'F', 'S', 'J', 'L')
                call self%open(fname, ioimg, form, readhead, rwaction='READ')
            case DEFAULT
                write(logfhandle,*) 'Trying to read from file: ', trim(fname)
                THROW_HARD('unsupported file format; read')
        end select
        if( isvol )then
            first_slice = 1
            last_slice = ldim(3)
        else
            first_slice = ii
            last_slice = ii
        endif
        call ioimg%rSlices(first_slice,last_slice,self%rmat,is_mrc=form.eq.'M')
        call ioimg%close
    end subroutine read

    !>  \brief  for writing any kind of images to stack or volumes to volume files
    subroutine write( self, fname, i, del_if_exists)
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: fname
        integer, optional, intent(in)    :: i
        logical, optional, intent(in)    :: del_if_exists
        real             :: dev, mean
        type(imgfile)    :: ioimg
        character(len=1) :: form
        integer          :: first_slice, last_slice, iform, ii
        logical          :: isvol, die
        isvol = .false.
        die   = .false.
        isvol = self%is_3d()
        if( present(del_if_exists) ) die = del_if_exists
        ii = 1 ! default location
        if( present(i) )then
            ! we are writing to a stack & in SIMPLE volumes are not allowed
            ! to be stacked so the image object must be 2D
            if( isvol ) THROW_HARD('trying to write 3D image to stack ; write')
            ii = i ! replace default location
        endif
        ! work out the slice range
        if( isvol )then
            first_slice = 1
            last_slice = self%ldim(3)
        else
            first_slice = ii
            last_slice = ii
        endif
        form = fname2format(fname)
        select case(form)
            case('M','F')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=die, formatchar=form, readhead=.false.)
                if( self%ft )then
                    call ioimg%setMode(4)
                else
                    call ioimg%setMode(2)
                endif
                call self%rmsd(dev, mean=mean)
                call ioimg%setRMSD(dev)
                call ioimg%setMean(mean)
                ! write slice(s) to disk & close
                call ioimg%wmrcSlices(first_slice,last_slice,self%rmat,self%ldim,self%ft)
                call ioimg%close
            case('S')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=die, formatchar=form, readhead=.false.)
                ! iform file type specifier:
                !   1 = 2D image
                !   3 = 3D volume
                ! -11 = 2D Fourier odd
                ! -12 = 2D Fourier even
                ! -21 = 3D Fourier odd
                ! -22 = 3D Fourier even
                if( self%is_2d() )then
                    if( self%ft )then
                        if( self%even_dims() )then
                            iform = -12
                        else
                            iform = -11
                        endif
                    else
                        iform = 1
                    endif
                else
                    if( self%ft )then
                        if( self%even_dims() )then
                            iform = -22
                        else
                            iform = -21
                        endif
                    else
                        iform = 3
                    endif
                endif
                call ioimg%setIform(iform)
                ! write slice(s) to disk & close
                call ioimg%wSlices(first_slice,last_slice,self%rmat,self%ldim,self%ft,self%smpd)
                call ioimg%close
            case DEFAULT
                write(logfhandle,*) 'format descriptor: ', form
                THROW_HARD('unsupported file format; write')
        end select
    end subroutine write

    !>  \brief  for updating header stats in a real space MRC image file only
    subroutine update_header_stats( self, fname, stats)
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: fname
        real,             intent(in)    :: stats(4) ! stats to update: min, max, mean, rms
        type(imgfile)    :: ioimg
        character(len=1) :: form
        if( self%existence )then
            if( self%ft )return
            form = fname2format(fname)
            select case(form)
            case('M','F')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=.false., formatchar=form, readhead=.true.)
                call ioimg%setMode(2) ! 32-bit reals (DEFAULT MODE)
                !  updates header
                call ioimg%update_MRC_stats(stats)
                ! writes header & close unit
                call ioimg%close
            case('S')
                ! spider stacks have one header each so we do nothing
                return
            case DEFAULT
                write(logfhandle,*) 'format descriptor: ', form
                THROW_HARD('unsupported file format; update_stats_header')
            end select
        else
            THROW_HARD('nonexisting image cannot be updated; update_header_stats')
        endif
    end subroutine update_header_stats

    subroutine write_jpg( self, fname, quality, colorspec, norm )
        use simple_jpg, only: jpg_img
        class(image),               intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        integer,          optional, intent(in)    :: quality, colorspec
        logical,          optional, intent(in)    :: norm
        type(jpg_img)     :: jpg
        integer           :: status
        logical           :: norm_here
        norm_here = .false.
        if(present(norm))norm_here = norm
        if( norm_here )call self%norm4viz
        if( self%is_2d() )then
            status = jpg%writejpg(trim(fname), self%rmat(:self%ldim(1),:self%ldim(2),1),&
                &quality=quality, colorspec=colorspec)
        else
            status = jpg%writejpg(trim(fname), self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                &quality=quality, colorspec=colorspec)
        endif
    end subroutine write_jpg

    ! GETTERS/SETTERS

    !> \brief get_array_shape  is a getter
    pure function get_array_shape( self ) result( shape)
        class(image), intent(in) :: self
        integer :: shape(3)
        shape = self%array_shape
    end function get_array_shape

    !> \brief get_ldim  is a getter
    pure function get_ldim( self ) result( ldim )
        class(image), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    !> \brief get_ldim  is a getter
    pure integer function get_box( self )
        class(image), intent(in) :: self
        get_box = self%ldim(1)
    end function get_box

    !> \brief get_smpd  is a getter
    pure function get_smpd( self ) result( smpd )
        class(image), intent(in) :: self
        real :: smpd
        smpd = self%smpd
    end function get_smpd

    !>  \brief get_nyq get the Nyquist Fourier index
    pure function get_nyq( self ) result( nyq )
        class(image), intent(in) :: self
        integer :: nyq
        nyq = fdim(self%ldim(1)) - 1
    end function get_nyq

    !> \brief get_filtsz  to get the size of the filters
    pure function get_filtsz( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = fdim(self%ldim(1)) - 1
    end function get_filtsz

    pure function get_shconst( self ) result( shconst )
        class(image), intent(in) :: self
        real :: shconst(3)
        shconst = self%shconst
    end function get_shconst

    function get( self, logi ) result( val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real :: val
        if( logi(1) > self%ldim(1) .or. logi(1) < 1 )then
            val = 0.
            return
        endif
        if( logi(2) > self%ldim(2) .or. logi(2) < 1 )then
            val = 0.
            return
        endif
        if( logi(3) > self%ldim(3) .or. logi(3) < 1 )then
            val = 0.
            return
        endif
        val = self%rmat(logi(1),logi(2),logi(3))
    end function get

    pure function get_rmat( self ) result( rmat )
        class(image), intent(in) :: self
        real, allocatable :: rmat(:,:,:)
        integer :: ldim(3)
        ldim = self%ldim
        allocate(rmat(ldim(1),ldim(2),ldim(3)), source=self%rmat(:ldim(1),:ldim(2),:ldim(3)))
    end function get_rmat

    subroutine get_mat_ptrs( self, mat_ptrs )
        class(image),      target, intent(in)  :: self
        class(image_ptr),          intent(out) :: mat_ptrs
        mat_ptrs%cmat => self%cmat
        mat_ptrs%rmat => self%rmat
    end subroutine get_mat_ptrs

    subroutine get_rmat_ptr( self, rmat_ptr )
        class(image), target,        intent(in)  :: self
        real(kind=c_float), pointer, intent(out) :: rmat_ptr(:,:,:)
        rmat_ptr => self%rmat
    end subroutine get_rmat_ptr

    pure subroutine get_rmat_sub( self, rmat )
        class(image), intent(in)  :: self
        real,         intent(out) :: rmat(self%ldim(1),self%ldim(2),self%ldim(3))
        rmat = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
    end subroutine get_rmat_sub

    pure function get_rmat_at_1( self, logi ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3)
        real :: val
        val = self%rmat(logi(1),logi(2),logi(3))
    end function get_rmat_at_1

    pure function get_rmat_at_2( self, i,j,k ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: i,j,k
        real :: val
        val = self%rmat(i,j,k)
    end function get_rmat_at_2

    function get_avg_int(self) result(avg_int)
        class(image), intent(in) :: self
        real    :: avg_int
        integer :: num_pixels
        if (self%is_3d()) then
            num_pixels = self%ldim(1) * self%ldim(2) * self%ldim(3)
            avg_int = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / num_pixels
        else
            num_pixels = self%ldim(1) * self%ldim(2)
            avg_int = sum(self%rmat(:self%ldim(1),:self%ldim(2),1)) / num_pixels
        end if
    end function get_avg_int

    function get_sum_int(self) result(sum_int)
        class(image), intent(in) :: self
        real :: sum_int
        sum_int = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
    end function get_sum_int

    pure subroutine set_rmat_at( self, i,j,k, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: val
        self%rmat(i,j,k) = val
    end subroutine set_rmat_at

    !>  \brief   get_cmat get the image object's complex matrix
    pure function get_cmat( self ) result( cmat )
        class(image), intent(in) :: self
        complex, allocatable :: cmat(:,:,:)
        allocate(cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3)), source=self%cmat)
    end function get_cmat

    subroutine get_cmat_ptr( self, cmat_ptr )
        class(image), target,                   intent(in)  :: self
        complex(kind=c_float_complex), pointer, intent(out) :: cmat_ptr(:,:,:)
        cmat_ptr => self%cmat
    end subroutine get_cmat_ptr

    !>  \brief   get_cmat get the image object's complex matrix
    pure subroutine get_cmat_sub( self, cmat )
        class(image), intent(in)  :: self
        complex,      intent(out) :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        cmat=self%cmat
    end subroutine get_cmat_sub

    pure function get_cmat_at_1( self, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: phys(3)
        complex :: comp
        comp = self%cmat(phys(1),phys(2),phys(3))
    end function get_cmat_at_1

    pure function get_cmat_at_2( self, h,k,l ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k,l
        complex :: comp
        comp = self%cmat(h,k,l)
    end function get_cmat_at_2

    pure subroutine set_cmat_1( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        self%ft   = .true.
        self%cmat = cmat
    end subroutine set_cmat_1

    pure subroutine set_cmat_2( self, cval )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cval
        self%ft   = .true.
        self%cmat = cval
    end subroutine set_cmat_2

    pure subroutine set_cmat_3( self, self2copy )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2copy
        self%ft   = .true.
        self%cmat = self2copy%cmat
    end subroutine set_cmat_3

    pure subroutine set_cmat_4( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2))
        self%ft          = .true.
        self%cmat(:,:,1) = cmat
    end subroutine set_cmat_4

    ! set comp to cmat at index phys
    pure subroutine set_cmat_at_1( self, phys ,comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: comp
        self%cmat(phys(1),phys(2),phys(3)) = comp
    end subroutine set_cmat_at_1

    !> add comp to cmat at index phys
    pure subroutine set_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer, intent(in) :: h,k,l
        complex, intent(in) :: comp
        self%cmat(h,k,l) =  comp
    end subroutine set_cmat_at_2

    !> add comp to cmat at index phys
    pure subroutine add_cmat_at_1( self , phys , comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: comp
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) + comp
    end subroutine add_cmat_at_1

    !> add comp to cmat at index (h,k,l)
    pure subroutine add_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        complex,      intent(in) :: comp
        self%cmat(h,k,l) = self%cmat(h,k,l) + comp
    end subroutine add_cmat_at_2

    !>  set complex matrices from images & arrays. Specialized routine for simple_classaverager
    subroutine set_cmats_from_cmats( self1 , self2 , self3, self4, self2set1, self2set2, lims, expcmat3, expcmat4)
        class(image), intent(in)    :: self1, self2,self3,self4
        class(image), intent(inout) :: self2set1, self2set2
        integer,      intent(in)    :: lims(3,2)
        real,         intent(inout) :: expcmat3(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,         intent(inout) :: expcmat4(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        integer :: h, k, logi(3), phys(3)
        !$omp parallel default(shared) private(h,k,logi,phys) proc_bind(close)
        !$omp workshare
        self1%cmat = self2set1%cmat
        self2%cmat = self2set2%cmat
        !$omp end workshare
        !$omp do collapse(2) schedule(static)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                logi = [h,k,0]
                phys = self1%comp_addr_phys(logi)
                self3%cmat(phys(1),phys(2),phys(3)) = cmplx(expcmat3(h,k),0.)
                self4%cmat(phys(1),phys(2),phys(3)) = cmplx(expcmat4(h,k),0.)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine set_cmats_from_cmats

    !>  adds complex matrices from images & arrays. Specialized routine for simple_classaverager
    subroutine add_cmats_to_cmats( self1 , self2 , self3, self4, self2set1, self2set2, lims, expcmat3, expcmat4)
        class(image), intent(in)    :: self1, self2,self3,self4
        class(image), intent(inout) :: self2set1, self2set2
        integer,      intent(in)    :: lims(3,2)
        real,         intent(inout) :: expcmat3(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,         intent(inout) :: expcmat4(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        integer :: h, k, logi(3), phys(3)
        !$omp parallel default(shared) private(h,k,logi,phys) proc_bind(close)
        !$omp workshare
        self1%cmat = self1%cmat + self2set1%cmat
        self2%cmat = self2%cmat + self2set2%cmat
        !$omp end workshare
        !$omp do collapse(2) schedule(static)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                logi = [h,k,0]
                phys = self1%comp_addr_phys(logi)
                self3%cmat(phys(1),phys(2),phys(3)) = self3%cmat(phys(1),phys(2),phys(3)) + cmplx(expcmat3(h,k),0.)
                self4%cmat(phys(1),phys(2),phys(3)) = self4%cmat(phys(1),phys(2),phys(3)) + cmplx(expcmat4(h,k),0.)
            enddo
        enddo
        !$omp end do nowait
        !$omp end parallel
    end subroutine add_cmats_to_cmats

    pure subroutine div_cmat_at_1( self, phys, rval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        real,              intent(in)    :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) / rval
        else
            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
        endif
    end subroutine div_cmat_at_1

    pure subroutine div_cmat_at_2( self,h,k,l, rval)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        real,         intent(in) :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / rval
        else
            self%cmat(h,k,l) =cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_2

    pure subroutine div_cmat_at_3( self, phys, cval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        complex,           intent(in)    :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) / cval
        else
            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
        endif
    end subroutine div_cmat_at_3

    pure subroutine div_cmat_at_4( self,h,k,l, cval)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / cval
        else
            self%cmat(h,k,l) = cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_4

    pure subroutine mul_cmat_at_1( self, phys, rval)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        real,         intent(in)    :: rval
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * rval
    end subroutine mul_cmat_at_1

    pure subroutine mul_cmat_at_2( self, phys, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: cval
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * cval
    end subroutine mul_cmat_at_2

    pure subroutine mul_cmat_at_3( self, h,k,l,rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        real,         intent(in)    :: rval
        self%cmat(h,k,l) = self%cmat(h,k,l) * rval
    end subroutine mul_cmat_at_3

    pure subroutine mul_cmat_at_4( self, h,k,l, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
        self%cmat(h,k,l) = self%cmat(h,k,l) * cval
    end subroutine mul_cmat_at_4

    subroutine mul_cmat_1( self, rmat )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        self%cmat = self%cmat * rmat
    end subroutine mul_cmat_1

    subroutine mul_cmat_2( self, rmat, resmsk )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        logical,      intent(in)    :: resmsk(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        where( resmsk )
            self%cmat = self%cmat * rmat
        end where
    end subroutine mul_cmat_2

    subroutine print_cmat( self )
        class(image), intent(in) :: self
        write(logfhandle,*) self%cmat
    end subroutine print_cmat

    subroutine print_rmat( self )
        class(image), intent(in) :: self
        write(logfhandle,*) self%rmat
    end subroutine print_rmat

    !>  \brief  expand_ft is for getting a Fourier plane using the old SIMPLE logics
    function expand_ft( self ) result( fplane )
        class(image), intent(in) :: self
        complex, allocatable :: fplane(:,:)
        integer :: xdim, ydim, h, k, phys(3)
        if(is_even(self%ldim(1)))then
            xdim = self%ldim(1)/2
            ydim = self%ldim(2)/2
        else
            xdim = (self%ldim(1)-1)/2
            ydim = (self%ldim(2)-1)/2
        endif
        allocate(fplane(-xdim:xdim,-ydim:ydim))
        fplane = cmplx(0.,0.)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                fplane(h,k) = self%get_fcomp([h,k,0],phys)
            end do
        end do
    end function expand_ft

    !>  \brief  set image value at position x,y,z
    subroutine set_1( self, logi, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: val
        if( logi(1) <= self%ldim(1) .and. logi(1) >= 1 .and. logi(2) <= self%ldim(2)&
            .and. logi(2) >= 1 .and. logi(3) <= self%ldim(3) .and. logi(3) >= 1 )then
            self%rmat(logi(1),logi(2),logi(3)) = val
        endif
    end subroutine set_1

    subroutine set_2( self, self2set )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2set
        self%ft       = self2set%ft
        self%rmat     = self2set%rmat
        self%wthreads = self2set%wthreads
    end subroutine set_2

    subroutine set_rmat( self, rmat, ft )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(:,:,:)
        logical,      intent(in)    :: ft
        self%rmat = 0.
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
        self%ft   = ft
    end subroutine set_rmat

    !>  \brief  set (replace) image data with new 3D data
    subroutine set_cmat( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(:,:,:)
        integer :: cdim(3)
        cdim(1) = size(cmat,1)
        cdim(2) = size(cmat,2)
        cdim(3) = size(cmat,3)
        if( all(self%array_shape .eq. cdim) )then
            self%ft   = .true.
            self%cmat = cmplx(0.,0.)
            self%cmat(:cdim(1),:cdim(2),:cdim(3)) = cmat
        else
            write(logfhandle,*) 'dim(cmat): ', cdim
            write(logfhandle,*) 'dim(self%cmat): ', self%array_shape
            THROW_HARD('nonconforming dims; set_cmat')
        endif
    end subroutine set_cmat

    !>  \brief set_smpd for setting smpd
    subroutine set_smpd( self, smpd )
        class(image), intent(inout) :: self
        real,         intent(in)    :: smpd
        self%smpd = smpd
        call self%fit%new(self%ldim, self%smpd)
    end subroutine set_smpd

    !> \brief get_slice is for getting a slice from a volume
    subroutine get_slice( self3D, slice, self2D )
        class(image), intent(in)    :: self3D
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self2D
        self2D%rmat(:,:,1) = self3D%rmat(:,:,slice)
    end subroutine get_slice

    !>  \brief set_slice is for putting a slice into a volume
    subroutine set_slice( self3D, slice, self2D )
        class(image), intent(in)    :: self2D
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self3D
        self3D%rmat(:,:,slice) = self2D%rmat(:,:,1)
    end subroutine set_slice

    pure function get_lfny( self, which ) result( fnyl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: fnyl
        fnyl = self%fit%get_lfny(which)
    end function get_lfny

    pure function get_lhp( self, which ) result( hpl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: hpl
        hpl = self%fit%get_lhp(which)
    end function get_lhp

    subroutine get_subimg(self, binning, xoffset, yoffset, img)
        class(image), intent(in)  :: self
        integer,      intent(in)  :: binning, xoffset, yoffset
        class(image), intent(out) :: img
        integer :: i,j,ii,jj,ldim(3)
        if( .not.is_even(binning) .and. binning > 0 )then
            THROW_HARD('Binning must be even; get_sub_img')
        endif
        ldim(1:2) = self%ldim(1:2) / binning
        ldim(3)   = 1
        call img%new(ldim, real(binning)*self%smpd)
        do j=1,ldim(2)
            jj = binning*(j-1) + yoffset + 1
            do i=1,ldim(1)
                ii = binning*(i-1) + xoffset + 1
                img%rmat(i,j,1) = self%rmat(ii,jj,1)
            enddo
        enddo
    end subroutine get_subimg

    pure function get_lp( self, ind ) result( lp )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: lp
        lp = self%fit%get_lp(1, ind)
    end function get_lp

    pure function get_spat_freq( self, ind ) result( spat_freq )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: spat_freq
        spat_freq = self%fit%get_spat_freq(1, ind)
    end function get_spat_freq

    pure function get_find( self, res ) result( ind )
        class(image), intent(in) :: self
        real,         intent(in) :: res
        integer :: ind
        ind = self%fit%get_find(1, res)
    end function get_find

    function rmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
        assoc = associated(self%rmat)
    end function rmat_associated

    function cmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
        assoc = associated(self%cmat)
    end function cmat_associated

    function is_wthreads( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%wthreads
    end function is_wthreads

    function serialize_1( self ) result( vec )
        class(image), intent(in) :: self
        real,    allocatable :: vec(:)
        complex, allocatable :: cvec(:)
        if( self%is_ft() )then
            cvec = pack(self%cmat(:self%array_shape(1),:self%array_shape(2),:self%array_shape(3)), mask=.true.)
            vec  = sqrt(real(cvec * conjg(cvec)))
        else
            vec = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
        endif
    end function serialize_1

    function serialize_2( self, l_msk )result( pcavec )
        class(image), intent(in) :: self
        logical,      intent(in) :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable :: pcavec(:)
        integer :: sz, cnt, i, j, k
        sz = count(l_msk)
        allocate(pcavec(sz))
        cnt = 0
        do k=1,self%ldim(3)
            do j=1,self%ldim(2)
                do i=1,self%ldim(1)
                    if( l_msk(i,j,k) )then
                        cnt         = cnt + 1
                        pcavec(cnt) = self%rmat(i,j,k)
                    endif
                end do
            end do
        end do
    end function serialize_2

    function serialize_3( self, thres ) result( vec )
        class(image), intent(in) :: self
        real,         intent(in) :: thres
        real, allocatable :: vec(:)
        vec = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > thres)
    end function serialize_3

    subroutine unserialize( self, pcavec, l_msk )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: pcavec(:)
        logical, optional, intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        integer :: sz, sz_msk, i, j, k, cnt
        if( present(l_msk) )then
            sz     = size(pcavec)
            sz_msk = count(l_msk)
            if( sz /= sz_msk )then
                write(logfhandle,*) 'ERROR! Nonconforming sizes'
                write(logfhandle,*) 'sizeof(pcavec): ', sz
                write(logfhandle,*) 'sizeof(l_msk) : ', sz_msk
                THROW_HARD('unserialize')
            endif
        endif
        if( self%ft ) self%ft = .false.
        self%rmat = 0.
        cnt = 0
        if( present(l_msk) )then
            do k=1,self%ldim(3)
                do j=1,self%ldim(2)
                    do i=1,self%ldim(1)
                        if( l_msk(i,j,k) )then
                            cnt = cnt + 1
                            self%rmat(i,j,k) =  pcavec(cnt)
                        endif
                    end do
                end do
            end do
        else
            do k=1,self%ldim(3)
                do j=1,self%ldim(2)
                    do i=1,self%ldim(1)
                        cnt = cnt + 1
                        self%rmat(i,j,k) =  pcavec(cnt)
                    end do
                end do
            end do
        endif
    end subroutine unserialize

    !>  \brief winserialize is for packing/unpacking a serialized image vector for convolutional pca analysis
    subroutine winserialize( self, coord, winsz, pcavec )
        class(image),      intent(inout) :: self
        real, allocatable, intent(inout) :: pcavec(:)
        integer,           intent(in)    :: coord(:), winsz
        integer :: i, j, k, cnt, npix
        logical :: pack
        if( self%ft ) THROW_HARD('winserialization not yet implemented for Fourier')
        if( self%is_2d() )then
            npix = winsz**2
            call set_action
            cnt = 0
            do i=coord(1),coord(1)+winsz-1
                do j=coord(2),coord(2)+winsz-1
                    cnt = cnt+1
                    if( pack )then
                        if( i > self%ldim(1) .or. j > self%ldim(2) )then
                            pcavec(cnt) = 0.
                        else
                            pcavec(cnt) = self%rmat(i,j,1)
                        endif
                    else
                        if( i > self%ldim(1) .or. j > self%ldim(2) )then
                        else
                            self%rmat(i,j,1) = self%rmat(i,j,1)+pcavec(cnt)
                        endif
                    endif
                end do
            end do
        else
            if( size(coord) < 3 ) THROW_HARD('need a 3D coordinate for a 3D image; winserialize')
            npix = winsz**3
            call set_action
            cnt = 0
            do i=coord(1),coord(1)+winsz-1
                do j=coord(2),coord(2)+winsz-1
                    do k=coord(3),coord(3)+winsz-1
                        cnt = cnt+1
                        if( pack )then
                            if( i > self%ldim(1) .or. j > self%ldim(2) .or. k > self%ldim(3) )then
                                pcavec(cnt) = 0.
                            else
                                pcavec(cnt) = self%rmat(i,j,k)
                            endif
                        else
                            if( i > self%ldim(1) .or. j > self%ldim(2) .or. k > self%ldim(3) )then
                            else
                                self%rmat(i,j,k) = self%rmat(i,j,k)+pcavec(cnt)
                            endif
                        endif
                    end do
                end do
            end do
        endif

        contains

            subroutine set_action
                if( allocated(pcavec) )then
                    if( size(pcavec) /= npix ) THROW_HARD('size mismatch mask/npix; winserialize')
                    pack = .false.
                else
                    pack = .true.
                    allocate( pcavec(npix) )
                    pcavec = 0.
                endif
            end subroutine set_action

    end subroutine winserialize

    !>  \brief  for swapping all zeroes in image with ones
    subroutine zero2one( self )
        class(image), intent(inout) :: self
        integer :: i, j, k
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if( is_zero(self%rmat(i,j,k)) ) self%rmat(i,j,k) = 1.
                end do
            end do
        end do
    end subroutine zero2one

    !>  \brief get_fcomp for getting a Fourier component from the compact representation
    pure function get_fcomp( self, logi, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3), phys(3)
        complex :: comp
        comp = merge(conjg(self%cmat(phys(1),phys(2),phys(3))),&
            &self%cmat(phys(1),phys(2),phys(3)), logi(1)<0)
    end function get_fcomp

    !>  \brief get_fcomp for getting a 2D only Fourier component from the compact representation
    elemental complex function get_fcomp2D(self, h, k)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k
        integer :: phys1, phys2
        if (h .ge. 0) then
            phys1 = h + 1
            phys2 = k + 1 + merge(self%ldim(2),0, k<0)
            get_fcomp2D = self%cmat(phys1,phys2,1)
        else
            phys1 = -h + 1
            phys2 = -k + 1 + merge(self%ldim(2),0, -k<0)
            get_fcomp2D = conjg(self%cmat(phys1,phys2,1))
        endif
    end function get_fcomp2D

    !> \brief set_fcomp  for setting a Fourier component in the compact representation
    subroutine set_fcomp( self, logi, phys, comp )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3), phys(3)
        complex,      intent(in)    :: comp
        complex :: comp_here
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = comp_here
    end subroutine set_fcomp

    !>  \brief vis is for plotting an image
    subroutine vis( self, sect, geomorsphr )
        class(image),      intent(in) :: self
        integer, optional, intent(in) :: sect
        logical, optional, intent(in) :: geomorsphr !< geometrical or spherical complex format
        complex, allocatable :: fplane(:,:)
        integer              :: sect_here
        logical              :: geomorsphr_here
        sect_here = 1
        if( present(sect) ) sect_here = sect
        geomorsphr_here=.true.
        if (present(geomorsphr))geomorsphr_here=geomorsphr
        if( self%ft )then
            if( self%ldim(3) == 1 ) sect_here = 0
            fplane = self%expand_ft()
            if(geomorsphr_here)then
                call gnufor_image(real(fplane), palette='gray')
                call gnufor_image(aimag(fplane), palette='gray')
            else
                call gnufor_image(cabs(fplane), palette='gray')
                call gnufor_image(atan2(real(fplane),aimag(fplane)), palette='gray')
            endif
            deallocate(fplane)
        else
            if( self%ldim(3) == 1 ) sect_here = 1
            call gnufor_image(self%rmat(:self%ldim(1),:self%ldim(2),sect_here), palette='gray')
        endif
    end subroutine vis

    !>  \brief  set_ft sets image ft state
    subroutine set_ft( self, is )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: is
        self%ft = is
    end subroutine set_ft

    ! CHECKUPS

    !>  \brief  Checks for existence
    pure function exists( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%existence
    end function exists

    !>  \brief  Checks whether the image is 2D
    pure logical function is_2d(self)
        class(image), intent(in)  ::  self
        is_2d = count(self%ldim .eq. 1) .eq. 1
    end function is_2d

    !>  \brief  Checks whether the image is 3D
    pure logical function is_3d(self)
        class(image), intent(in)  ::  self
        is_3d = .not. any(self%ldim .eq. 1)
    end function is_3d

    !>  \brief  checks for even dimensions
    pure function even_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep, test(2)
        test = .false.
        test(1) = is_even(self%ldim(1))
        test(2) = is_even(self%ldim(2))
        yep = all(test)
    end function even_dims

    !>  \brief  checks for square dimensions
    pure function square_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep
        yep = self%ldim(1) == self%ldim(2)
        if( self%ldim(3) == 1 .and. yep )then
        else
            yep = self%ldim(3) == self%ldim(1)
        endif
    end function square_dims

    !>  \brief  checks for same dimensions, overloaded as (.eqdims.)
    pure logical function same_dims_1( self1, self2 )
        class(image), intent(in) :: self1, self2
        same_dims_1 = all(self1%ldim == self2%ldim)
    end function same_dims_1

    !>  \brief  checks for same dimensions
    pure logical function same_dims( self, ldim )
        class(image), intent(in) :: self
        integer,      intent(in) :: ldim(3) !< dimensions
        same_dims = all(self%ldim == ldim)
    end function same_dims

    !>  \brief  checks for same sampling distance, overloaded as (.eqsmpd.)
    logical pure function same_smpd( self1, self2 )
        class(image), intent(in) :: self1, self2
        same_smpd = abs(self1%smpd-self2%smpd) < 0.0001
    end function same_smpd

    !>  \brief  checks if image is ft
    pure function is_ft( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%ft
    end function is_ft

    function is_empty( self ) result( is )
        class(image), intent(in) :: self
        real    :: minmax(2)
        logical :: is
        minmax = self%minmax()
        is     = is_equal(minmax(2)-minmax(1),0.) ! empty image
    end function is_empty

    ! ARITHMETICS

    !>  \brief  assign, polymorphic assignment (=)
    subroutine assign( selfout, selfin )
        class(image), intent(inout) :: selfout
        class(image), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    !>  \brief assign_r2img real constant to image assignment(=) operation
    subroutine assign_r2img( self, realin )
        class(image), intent(inout) :: self
        real,         intent(in)    :: realin
        self%rmat = realin
        self%ft = .false.
    end subroutine assign_r2img

    !>  \brief  assign_c2img  complex constant to image assignment(=) operation
    subroutine assign_c2img( self, compin )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: compin
        self%cmat = compin
        self%ft = .true.
    end subroutine assign_c2img

    !>  \brief  is for image addition(+) addition
    function addition( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .neqv. self2%ft )then
                THROW_HARD('cannot add images of different FT state; addition(+)')
            endif
            self%rmat = self1%rmat+self2%rmat
        else
            THROW_HARD('cannot add images of different dims; addition(+)')
        endif
        self%ft = self1%ft
    end function addition

    !>  \brief  is for image addition(+) addition
    !! \param self1 image object 1
    !! \param self2  image object 2
    !! \return lhs, copy of added images
    !!
    function addition_const_real( self1, rconst ) result( self )
        class(image), intent(in) :: self1
        real, intent(in) :: rconst
        type(image) :: self
        call self%new(self1%ldim, self1%smpd)
        self%rmat = self1%rmat+rconst
        self%ft = self1%ft
    end function addition_const_real

    !>  \brief  square root of image, real space only
    subroutine sq_rt( self )
        class(image), intent(inout) :: self
        if( self%ft ) THROW_HARD('Real space only; sq_rt')
        where( self%rmat > 0. )
            self%rmat = sqrt(self%rmat)
        else where
            self%rmat = 0.
        end where
    end subroutine sq_rt

    !>  \brief add_1 is for image summation, not overloaded
    subroutine add_1( self, self_to_add, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_add
        real, optional, intent(in)    :: w
        real :: ww
        ww = 1.
        if( present(w) ) ww = w
        if( self%ft )then
            self%cmat = self%cmat+ww*self_to_add%cmat
        else
            self%rmat = self%rmat+ww*self_to_add%rmat
        endif
    end subroutine add_1

    !>  \brief add_2 is for componentwise summation, not overloaded
    subroutine add_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(in)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot add complex number to real image; add_2')
        if( present(phys_in) )then
            phys = phys_in
        else
            phys = self%fit%comp_addr_phys(logi)
        endif
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))+comp_here
        if( present(phys_out) ) phys_out = phys
    end subroutine add_2

    !> \brief add_3  is for componentwise summation, not overloaded
    subroutine add_3( self, rcomp, i, j, k )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rcomp
        integer,      intent(in)    :: i, j, k
        if(  self%ft ) THROW_HARD('cannot add real number to transform; add_3')
        self%rmat(i,j,k) = self%rmat(i,j,k)+rcomp
    end subroutine add_3

    subroutine add_4( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot add complex number to real image; add_2')
        phys = self%fit%comp_addr_phys(logi)
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        if( abs(k(phys(1),phys(2),phys(3))) > 1e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))+(comp_here/k(phys(1),phys(2),phys(3)))*w
        endif
    end subroutine add_4

    !>  \brief  is for adding a constant
    subroutine add_5( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        if( self%ft )then
            self%cmat = self%cmat+cmplx(c,0.)
        else
            self%rmat = self%rmat+c
        endif
    end subroutine add_5

    subroutine add_workshare( self, self_to_add )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_add
        if( self%ft )then
            !$omp parallel workshare proc_bind(close)
            self%cmat = self%cmat + self_to_add%cmat
            !$omp end parallel workshare
        else
            !$omp parallel workshare proc_bind(close)
            self%rmat = self%rmat + self_to_add%rmat
            !$omp end parallel workshare
        endif
    end subroutine add_workshare

    !>  \brief subtraction is for image subtraction(-)
    function subtraction( self_from, self_to ) result( self )
        class(image), intent(in) :: self_from, self_to
        type(image) :: self
        if( self_from.eqdims.self_to )then
            call self%new(self_from%ldim, self_from%smpd)
            if( self_from%ft .neqv. self_to%ft )then
                THROW_HARD('cannot subtract images of different FT state; subtraction(+)')
            endif
            self%rmat = self_from%rmat-self_to%rmat
        else
            THROW_HARD('cannot subtract images of different dims; subtraction(-)')
        endif
    end function subtraction

    !>  \brief subtr_1 is for image subtraction, not overloaded
    subroutine subtr_1( self, self_to_subtr, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_subtr
        real, optional, intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( self%ft )then
            self%cmat = self%cmat-ww*self_to_subtr%cmat
        else
            self%rmat = self%rmat-ww*self_to_subtr%rmat
        endif
    end subroutine subtr_1

    !>  \brief subtr_2 is for componentwise subtraction, not overloaded
    subroutine subtr_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(out)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot subtract complex number from real image; subtr_2')
        if( present(phys_in) )then
            phys = phys_in
        else
            phys = self%fit%comp_addr_phys(logi)
        endif
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))-comp_here
        if( present(phys_out) ) phys_out = phys
    end subroutine subtr_2

    subroutine subtr_3( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot subtract complex number from real image; subtr_3')
        phys = self%fit%comp_addr_phys(logi)
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        if( .not. is_zero(k(phys(1),phys(2),phys(3))) )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) -&
                &(comp_here/k(phys(1),phys(2),phys(3)))*w
        endif
    end subroutine subtr_3

    !>  \brief subtr_4 is for subtracting a constant from a real image, not overloaded
    subroutine subtr_4( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        self%rmat = self%rmat-c
    end subroutine subtr_4

    !>  \brief multiplication is for image multiplication(*)
    function multiplication( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .and. self2%ft )then
                self%cmat = self1%cmat*self2%cmat
                self%ft = .true.
            else if( self1%ft .eqv. self2%ft )then
                self%rmat = self1%rmat*self2%rmat
                self%ft = .false.
            else if(self1%ft)then
                self%cmat = self1%cmat*self2%rmat
                self%ft = .true.
            else
                self%cmat = self1%rmat*self2%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot multiply images of different dims; multiplication(*)')
        endif
    end function multiplication

    function multiplication_const_real( self1, rconst ) result( self )
        class(image), intent(in) :: self1
        real,         intent(in) :: rconst
        type(image) :: self
        call self%new(self1%ldim, self1%smpd)
        self%ft = self1%ft
        if(self1%ft)then
            self%cmat = self1%cmat*rconst
        else
            self%rmat = self1%rmat*rconst
        endif
    end function multiplication_const_real

    function multiplication_const_int( self1, iconst ) result( self )
        class(image), intent(in) :: self1
        integer,      intent(in) :: iconst
        type(image) :: self
        call self%new(self1%ldim, self1%smpd)
        self%ft = self1%ft
        if(self1%ft)then
            self%cmat = self1%cmat*iconst
        else
            self%rmat = self1%rmat*iconst
        endif
    end function multiplication_const_int

    ! elementwise multiplication in real-space
    subroutine mul_rmat_at_1( self, logi, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: rval
        self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))*rval
    end subroutine mul_rmat_at_1

    ! elementwise multiplication in real-space
    elemental pure subroutine mul_rmat_at_2( self,i, j, k, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: rval
        self%rmat(i,j,k) = self%rmat(i,j,k)*rval
    end subroutine mul_rmat_at_2

    !>  \brief mul_1 is for component-wise multiplication of an image with a real constant
    subroutine mul_1( self, logi, rc, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        real,              intent(in)    :: rc
        integer, optional, intent(in)    :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        if( self%is_ft() )then
            if( present(phys_in) )then
                phys = phys_in
            else
                phys = self%fit%comp_addr_phys(logi)
            endif
            if( present(phys_out) ) phys_out = phys
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*rc
        else
            self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))*rc
        endif
    end subroutine mul_1

    !>  \brief mul_2 is for  multiplication of an image with a real constant
    subroutine mul_2( self, rc )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rc
        if( self%ft )then
            self%cmat = self%cmat*rc
        else
            self%rmat = self%rmat*rc
        endif
    end subroutine mul_2

    !>  \brief mul_3 is for multiplication of images
    subroutine mul_3( self, self2mul )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
        if( self.eqdims.self2mul )then
            if( self%ft .and. self2mul%ft )then
                self%cmat = self%cmat*self2mul%cmat
            else if( self%ft .eqv. self2mul%ft )then
                self%rmat = self%rmat*self2mul%rmat
                self%ft = .false.
            else if(self%ft)then
                self%cmat = self%cmat*self2mul%rmat
            else
                self%cmat = self%rmat*self2mul%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot multiply images of different dims; mul_3')
        endif
    end subroutine mul_3

    !>  \brief mul_4 is for low-pass limited multiplication of images
    subroutine mul_4( self, self2mul, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
        real,         intent(in)    :: lp
        integer                     :: lims(3,2),sqlim,h,k,l,phys(3)
        if( .not. self%is_ft() )     THROW_HARD('low-pass limited multiplication requires self to be FT')
        if( .not. self2mul%is_ft() ) THROW_HARD('low-pass limited multiplication requires self2mul to be FT')
        if( self.eqdims.self2mul )then
            lims = self%fit%loop_lims(1,lp)
            sqlim = (maxval(lims(:,2)))**2
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        if( h * h + k * k + l * l <= sqlim )then
                            phys = self%fit%comp_addr_phys([h,k,l])
                            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
                                self2mul%cmat(phys(1),phys(2),phys(3))
                        endif
                    end do
                end do
            end do
        else
            THROW_HARD('cannot multiply images of different dims; mul_3')
        endif
    end subroutine mul_4

    !>  \brief mul_1 is for component-wise multiplication of an image with a complex constant
    subroutine mul_5( self, logi, c)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: c
        integer :: phys(3)
        phys = self%fit%comp_addr_phys(logi)
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * c
    end subroutine mul_5

    !>  \brief division is for image division(/)
    function division( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .and. self2%ft )then
                self%cmat = self1%cmat/self2%cmat
                self%ft = .true.
            else if( self1%ft .eqv. self2%ft )then
                self%rmat = self1%rmat/self2%rmat
                self%ft = .false.
            else if(self1%ft)then
                self%cmat = self1%cmat/self2%rmat
                self%ft = .true.
            else
                self%cmat = self1%rmat/self2%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot divide images of different dims; division(/)')
        endif
    end function division

    ! component-wise division of an image with a real number
    subroutine div_rmat_at_1( self, logi, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: rval
        if( abs(rval) > 1e-6 )then
            self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))/rval
        endif
    end subroutine div_rmat_at_1

    ! component-wise division of an image with a real number
    subroutine div_rmat_at_2( self, i, j, k, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: rval
        if( abs(rval) > 1e-6 )then
            self%rmat(i,j,k) = self%rmat(i,j,k)/rval
        endif
    end subroutine div_rmat_at_2

    !>  \brief div_1 is for dividing image with real constant, not overloaded
    subroutine div_1( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        if( abs(c) < 1e-6 )then
            THROW_HARD('division with zero; div_1')
        else
            if( self%ft )then
                self%cmat = self%cmat/c
            else
                self%rmat = self%rmat/c
            endif
        endif
    end subroutine div_1

    !>  \brief div_2 is for component-wise matrix division of a Fourier transform with a real matrix, k
    subroutine div_2( self, logi, k, square )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: k(:,:,:)
        logical,      intent(in)    :: square
        integer :: phys(3)
        if( self%ft )then
            phys = self%fit%comp_addr_phys(logi)
            if( abs(k(phys(1),phys(2),phys(3))) > 1e-6 )then
                self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/k(phys(1),phys(2),phys(3))
            else
                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
            endif
        else
            THROW_HARD('Image need to be Fourier transformed; div_2')
        endif
    end subroutine div_2

    !>  \brief div_3 is for component-wise division of an image with a real number
    subroutine div_3( self, logi, k, phys_in )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        real,              intent(in)    :: k
        integer, optional, intent(in)    :: phys_in(3)
        integer :: phys(3)
        if( self%ft )then
            if( present(phys_in) )then
                phys = phys_in
            else
                phys = self%fit%comp_addr_phys(logi)
            endif
            if( abs(k) > 1e-6 )then
                self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/k
            else
                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
            endif
        else
            if( abs(k) > 1e-6 )then
                self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))/k
            endif
        endif
    end subroutine div_3

    !>  \brief div_4 is for division of images
    subroutine div_4( self, self2div )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2div
        if( self.eqdims.self2div )then
            if( self%ft .and. self2div%ft )then
                self%cmat = self%cmat/self2div%cmat
            else if( self%ft .eqv. self2div%ft )then
                where(abs(self2div%rmat) > 1.e-6) self%rmat = self%rmat/self2div%rmat
                self%ft = .false.
            else if(self%ft)then
              where(abs(self2div%rmat) > 1.e-6) self%cmat = self%cmat/self2div%rmat
            else
                self%cmat = self%rmat/self2div%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot divide images of different dims; div_4')
        endif
    end subroutine div_4

    !> \brief ctf_dens_correct for sampling density compensation & Wiener normalization
    !! \param self_sum sum image
    !! \param self_rho density image
    !! KEEP SERIAL
    subroutine ctf_dens_correct( self_sum, self_rho )
        class(image), intent(inout) :: self_sum
        class(image), intent(inout) :: self_rho
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh    = nint(hyp(h,k,l))
                    phys  = self_sum%comp_addr_phys([h,k,l])
                    denom = real(self_rho%cmat(phys(1),phys(2),phys(3)))
                    if(sh <= nyq .and. abs(denom) > 1.e-10 )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3)) / denom
                    else
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct

    subroutine ctf_dens_correct_wiener( self_sum, self_rho, ssnr )
        class(image), intent(inout) :: self_sum
        class(image), intent(in)    :: self_rho
        real,         intent(in)    :: ssnr(:)
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(h,k,l))
                    phys = self_sum%comp_addr_phys([h,k,l])
                    if(sh > nyq )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh==0 )then
                            denom = real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        else
                            denom = ssnr(sh)*real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = ssnr(sh)*self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        endif
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct_wiener

    !>  \brief conjugate is for complex conjugation of a FT
    function conjugate( self ) result ( self_out )
        class(image), intent(in) :: self
        type(image) :: self_out
        if( self%ft )then
            call self_out%copy(self)
            self_out%cmat = conjg(self%cmat)
        else
            THROW_WARN('cannot conjugate real image')
        endif
    end function conjugate

    ! BINARY IMAGE METHODS

    !>  \brief nforeground counts the number of foreground (white) pixels in a binary image
    function nforeground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = count(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.5)
    end function nforeground

    !>  \brief nbackground counts the number of background (black) pixels in a binary image
    function nbackground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = product(self%ldim)-self%nforeground()
    end function nbackground

    subroutine binarize_1( self_in, thres, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: thres
        class(image), optional, intent(inout) :: self_out
        integer :: n_foreground
        n_foreground = count(self_in%rmat > thres)
        if( n_foreground < 1 ) THROW_HARD('Binarization produces empty image!')
        if( self_in%ft ) THROW_HARD('only for real images; bin_1')
        if( present(self_out) )then
            if( any(self_in%ldim /= self_out%ldim)) THROW_HARD('Images dimensions are not compatible; binarize_1')
            where( self_in%rmat >= thres )
                self_out%rmat  = 1.
            elsewhere
                self_out%rmat  = 0.
            end where
        else
            where( self_in%rmat >= thres )
                self_in%rmat = 1.
            elsewhere
                self_in%rmat = 0.
            end where
        endif
    end subroutine binarize_1

    subroutine binarize_2( self, npix )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: npix
        real, allocatable :: forsort(:)
        real    :: thres
        integer :: npixtot
        if( self%ft ) THROW_HARD('only for real images')
        npixtot = product(self%ldim)
        forsort = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), .true.)
        call hpsort(forsort)
        thres = forsort(npixtot-npix-1) ! everyting above this value 1 else 0
        call self%binarize_1( thres )
        deallocate( forsort )
    end subroutine binarize_2

    subroutine binarize_3( self, thres, mask )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        logical,      intent(inout) :: mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        if( self%ft ) THROW_HARD('only for real images')
        where( self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) >= thres )
            mask = .true.
        elsewhere
            mask = .false.
        end where
    end subroutine binarize_3

    !>  \brief cendist produces an image holding the distances from the centre of the image or from an arbitrary input point
    subroutine cendist( self, c_point )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: c_point(3)
        real    :: centre(3)
        integer :: i
        if( self%ft ) THROW_HARD('real space only; cendist')
        ! Builds square distance image
        self   = 0.
        if( present(c_point) )then
            centre = c_point
        else
            centre = real(self%ldim)/2.+1.
        endif
        if( self%is_2d() )then
            ! 2D
            do i=1,self%ldim(1)
                self%rmat(i,:,1) = self%rmat(i,:,1) + (real(i)-centre(1))**2.
            enddo
            do i=1,self%ldim(2)
                self%rmat(:,i,1) = self%rmat(:,i,1) + (real(i)-centre(2))**2.
            enddo
        else
            ! 3D
            do i=1,self%ldim(1)
                self%rmat(i,:,:) = self%rmat(i,:,:) + (real(i)-centre(1))**2.
            enddo
            do i=1,self%ldim(2)
                self%rmat(:,i,:) = self%rmat(:,i,:) + (real(i)-centre(2))**2.
            enddo
            do i=1,self%ldim(3)
                self%rmat(:,:,i) = self%rmat(:,:,i) + (real(i)-centre(3))**2.
            enddo
        endif
        self%rmat = sqrt(self%rmat)
    end subroutine cendist

    !> to find centre of gravity
    subroutine masscen( self, xyz, mask_in )
        class(image),      intent(inout) :: self
        real        ,      intent(out)   :: xyz(3)
        logical, optional, intent(in)    :: mask_in(:,:,:)
        real                 ::  spix, ci, cj, ck
        integer              :: i, j, k
        logical, allocatable :: mask_here(:,:,:)
        allocate(mask_here(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
        if (present(mask_in)) then 
            if (.not.(size(mask_in, dim=1) .eq. self%ldim(1))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=2) .eq. self%ldim(2))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=3) .eq. self%ldim(3))) THROW_HARD('mask_in dimension must match dimension of image')
            mask_here = mask_in   
        end if
        if( self%is_ft() ) THROW_HARD('masscen not implemented for FTs; masscen')
        spix = 0.
        xyz  = 0.
        ci   = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if (mask_here(i,j,k)) then
                        xyz  = xyz  + self%rmat(i,j,k) * [ci, cj, ck]
                        spix = spix + self%rmat(i,j,k) 
                    end if
                    ck = ck + 1.
                end do
                cj = cj + 1.
            end do
            ci = ci + 1.
        end do
        if( is_equal(spix,0.) ) return
        xyz = xyz / spix
        if( self%ldim(3) == 1 ) xyz(3) = 0.
    end subroutine masscen

    !> to find center of gravity in absolute, not relative, terms
    subroutine masscen_adjusted( self, xyz, mask_in )
        class(image),      intent(inout) :: self
        real        ,      intent(out)   :: xyz(3)
        logical, optional, intent(in)    :: mask_in(:,:,:)
        real                 ::  spix
        integer              :: i, j, k
        logical, allocatable :: mask_here(:,:,:)
        allocate(mask_here(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
        if (present(mask_in)) then 
            if (.not.(size(mask_in, dim=1) .eq. self%ldim(1))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=2) .eq. self%ldim(2))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=3) .eq. self%ldim(3))) THROW_HARD('mask_in dimension must match dimension of image')
            mask_here = mask_in   
        end if
        if( self%is_ft() ) THROW_HARD('masscen not implemented for FTs; masscen')
        spix = 0.
        xyz  = 0.
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if (mask_here(i,j,k)) then
                        xyz  = xyz  + self%rmat(i,j,k) * [i, j, k]
                        spix = spix + self%rmat(i,j,k) 
                    end if
                end do
            end do
        end do
        if( is_equal(spix,0.) ) return
        xyz = xyz / spix
        if( self%ldim(3) == 1 ) xyz(3) = 0.
    end subroutine masscen_adjusted

    function box_cen_arg( self, tmp ) result( a )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: tmp
        real :: a, xyz(3)
        call tmp%copy_fast(self)
        call tmp%norm_minmax
        where(tmp%rmat < TINY) tmp%rmat=0.
        call tmp%mask(real(self%ldim(1))/2., 'hard')
        call tmp%masscen(xyz)
        if (self%ldim(3) == 1) then
            a = arg(xyz(:2))
        else
            a = arg(xyz)
        end if
    end function box_cen_arg

    !>  \brief is for estimating the center of an image based on center of mass
    function calc_shiftcen( self, lp, msk, hp ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real, optional, intent(in)    :: msk, hp
        type(image) :: tmp
        real        :: xyz(3), rmsk
        if( present(msk) )then
            rmsk = msk
        else
            rmsk = real( self%ldim(1) )/2. - 5. ! 5 pixels outer width
        endif
        call tmp%copy(self)
        if( present(hp) )then
            call tmp%bp(hp, lp)
        else
            call tmp%bp(0., lp)
        endif
        call tmp%ifft()
        call tmp%mask(rmsk, 'hard')
        ! such that norm_minmax will neglect everything < 0. and preserve zero
        where(tmp%rmat < TINY) tmp%rmat=0.
        call tmp%norm_minmax
        call tmp%masscen(xyz)
        call tmp%kill
    end function calc_shiftcen

    function calc_shiftcen_serial( self, lp, msk, hp ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real,           intent(in)    :: msk
        real, optional, intent(in)    :: hp
        real    :: xyz(3)
        integer :: ithr
        ! get thread index
        ithr = omp_get_thread_num() + 1
        if( all(self%ldim == thread_safe_tmp_imgs(ithr)%ldim) )then
            ! copy rmat
            thread_safe_tmp_imgs(ithr)%rmat = self%rmat
            thread_safe_tmp_imgs(ithr)%ft   = .false.
            call thread_safe_tmp_imgs(ithr)%fft()
            if( present(hp) )then
                call thread_safe_tmp_imgs(ithr)%bp(hp, lp)
            else
                call thread_safe_tmp_imgs(ithr)%bp(0., lp)
            endif
            call thread_safe_tmp_imgs(ithr)%ifft()
            call thread_safe_tmp_imgs(ithr)%mask(msk, 'hard')
            where(thread_safe_tmp_imgs(ithr)%rmat < TINY) thread_safe_tmp_imgs(ithr)%rmat = 0.
            call thread_safe_tmp_imgs(ithr)%norm_minmax
            call thread_safe_tmp_imgs(ithr)%masscen(xyz)
        else
            THROW_HARD('Incompatible dimensions bwetween self and thread_safe_tmp_imgs; calc_shiftcen_serial')
        endif
    end function calc_shiftcen_serial

    !>  \brief bin_inv inverts a binary image
    subroutine bin_inv( self )
        class(image), intent(inout) :: self
        self%rmat = -1.*(self%rmat-1.)
    end subroutine bin_inv

    !>  \brief  remove edge from binary image
    subroutine remove_edge( self )
        class(image), intent(inout) :: self
        if( self%ft ) THROW_HARD('only for real binary images (not FTed ones); remove_edge')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input to remove edge not binary; remove_edge')
        where( self%rmat < 0.999 ) self%rmat = 0.
    end subroutine remove_edge

    !>  \brief  set the edge pixels to one
    subroutine one_at_edge( self )
        class(image), intent(inout) :: self
        if( self%ft ) THROW_HARD('only for real binary images (not FTed ones); one_at_edge')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input to one_at_edge not binary')
        where( self%rmat < 0.999 .and. self%rmat > TINY ) self%rmat = 1.
    end subroutine one_at_edge

    !>  \brief  generates a logical mask from a binary one
    function bin2logical( self ) result( mask )
        class(image), intent(in)  :: self
        logical,      allocatable :: mask(:,:,:)
        allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)),source=.false.)
        where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > TINY )
            mask = .true.
        end where
    end function bin2logical

    !>  \brief  generates a binary image from a logical mask
    subroutine logical2bin( self, mask )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        self%rmat = 0.
        where(mask) self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
    end subroutine logical2bin

    subroutine density_inoutside( self, msk, nin, nout, nmsk )
        class(image), intent(in)  :: self
        real,         intent(in)  :: msk
        integer,      intent(out) :: nin, nout, nmsk
        integer :: i
        real    :: centre(3), cendists(self%ldim(1),self%ldim(2)), msksq
        ! calculate distances from the center
        centre = real(self%ldim)/2.+1.
        cendists = 0.
        do i=1,self%ldim(1)
            cendists(i,:) = cendists(i,:) + (real(i)-centre(1))**2.
        enddo
        do i=1,self%ldim(2)
            cendists(:,i) = cendists(:,i) + (real(i)-centre(2))**2.
        enddo
        ! pixels forming the mask
        msksq = msk**2
        nmsk  = count(cendists <= msksq)
        ! pixels inside mask that are foreground
        nin = count(cendists <= msksq .and. self%rmat(:self%ldim(1),:self%ldim(2),1) > 0.5)
        ! pixels outside mask that are foreground
        nout = count(cendists > msksq .and. self%rmat(:self%ldim(1),:self%ldim(2),1) > 0.5)
    end subroutine density_inoutside

    ! performs left/right collage of input images, un-modified on output
    subroutine collage( self1, self2, img_out )
        class(image), intent(inout) :: self1, self2, img_out
        real, parameter :: background = 128. ! taken as centre of [0.255] for jpegs
        type(image)     :: img_pad
        integer         :: ldim(3), ldim_col(3), border
        if( .not.self1%is_2d() ) THROW_HARD('2D only; collage')
        if( self1%is_ft() )      THROW_HARD('Real space only; collage')
        if( .not.self2%is_2d() ) THROW_HARD('2D only; collage')
        if( self2%is_ft() )      THROW_HARD('Real space only; collage')
        border   = 1
        ldim(1)  = max(self1%ldim(1),self2%ldim(1))
        ldim(2)  = max(self1%ldim(2),self2%ldim(2))
        ldim(1)  = max(ldim(1), ldim(2))
        ldim(2)  = ldim(1)
        ldim(3)  = 1
        ldim_col = [2*ldim(1)+border, ldim(2), 1]
        call img_out%new(ldim_col,1.)
        img_out%rmat = background
        ! pad & copy left image
        call img_pad%new(ldim,self1%get_smpd())
        img_pad%rmat = background
        call self1%norm4viz
        call self1%pad(img_pad, backgr=background)
        img_out%rmat(:ldim(1),:ldim(2),1) = img_pad%rmat(:ldim(1),:ldim(2),1)
        ! pad & copy right image
        img_pad%rmat = background
        call self2%norm4viz
        call img_pad%set_smpd(self2%get_smpd())
        call self2%pad(img_pad, backgr=background)
        img_out%rmat(ldim(1)+border+1:ldim_col(1),:ldim_col(2),1) = img_pad%rmat(:ldim(1),:ldim(2),1)
        call img_pad%kill()
    end subroutine collage

    ! place stack i in tile position x,y
    subroutine tile( self, stkimg, x, y)
        class(image), intent(inout) :: self
        class(image), intent(inout) :: stkimg
        integer,      intent(in)    :: x, y
        integer                     :: stkimg_ldim(3), x_start, y_start, x_end, y_end
        stkimg_ldim = stkimg%get_ldim()
        x_start = (x - 1) * stkimg_ldim(1) + 1
        y_start = (y - 1) * stkimg_ldim(2) + 1
        x_end = x * stkimg_ldim(1)
        y_end = y * stkimg_ldim(2)
        if(x_start .lt. 1 .or. y_start .lt. 1) THROW_HARD('tile: out of bounds')
        if(x_end .gt. self%ldim(1) .or. y_end .gt. self%ldim(2)) THROW_HARD('tile: out of bounds')
        call stkimg%norm4viz(brightness=80.0, maxmin=.true.)
        self%rmat(x_start:x_end, y_start:y_end, 1) = stkimg%rmat(:stkimg_ldim(1), :stkimg_ldim(2), 1)
    end subroutine tile

    ! Calculates the rotation matrix that aligns the inertia tensor of object to xyz cartesian axes
    subroutine calc_principal_axes_rotmat( self, radius, R )
        class(image), intent(in)  :: self
        real,         intent(in)  :: radius
        real,         intent(out) :: R(3,3)
        real(dp) :: coord(3), ixx, iyy, izz, ixz, ixy, iyz, m
        real(dp) :: inertia(3,3), eigvals(3), eigvecs(3,3)
        real     :: radiussq
        integer  :: icenter(3),i,j,k
        if( self%is_ft() )      THROW_HARD('Real space only; calc_principal_axes_rotmat')
        if( .not.self%is_3d() ) THROW_HARD('Volumes only; calc_principal_axes_rotmat')
        icenter  = nint(real(self%ldim)/2.)+1
        radiussq = radius**2
        ! Inertia Tensor
        ixx = 0.d0; iyy = 0.d0; izz = 0.d0
        ixy = 0.d0; ixz = 0.d0; iyz = 0.d0
        do k =1,self%ldim(3)
        do j =1,self%ldim(2)
        do i =1,self%ldim(1)
            if( (real(sum(([i,j,k]-icenter)**2)) < radiussq) .and. (self%rmat(i,j,k)>0.0) )then
                coord = real([i,j,k]-icenter, dp)
                m     = real(self%rmat(i,j,k), dp)
                ixx   = ixx + m * (coord(2)**2 + coord(3)**2)
                iyy   = iyy + m * (coord(1)**2 + coord(3)**2)
                izz   = izz + m * (coord(1)**2 + coord(2)**2)
                ixy   = ixy + m * coord(1) * coord(2)
                ixz   = ixz + m * coord(1) * coord(3)
                iyz   = iyz + m * coord(2) * coord(3)
            endif
        enddo
        enddo
        enddo
        inertia(1,:) = [ ixx, -ixy, -ixz]
        inertia(2,:) = [-ixy,  iyy, -iyz]
        inertia(3,:) = [-ixz, -iyz,  izz]
        ! Spectral analysis
        call svdcmp(inertia, eigvals, eigvecs)
        call eigsrt(eigvals, eigvecs, 3, 3)
        ! double checking
        ! identity = matmul(eigvecs, transpose(eigvecs))
        ! inertia  = matmul(eigvecs, matmul(eye(3)*eigvals, transpose(eigvecs)))
        ! Reverse rotation matrix
        R = real(transpose(eigvecs))
    end subroutine calc_principal_axes_rotmat

    ! generate the 3 orthogonal reprojections from a volume into a single image
    subroutine generate_orthogonal_reprojs( self, reprojs )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: reprojs
        integer, parameter :: b=3
        type(image) :: reproj
        integer     :: ldim_reproj(3),ldim_reprojs(3)
        if( self%is_ft() )      THROW_HARD('Real space only; generate_orthogonal_reprojs')
        if( .not.self%is_3d() ) THROW_HARD('Volumes only; generate_orthogonal_reprojs')
        ldim_reproj(1:2) = self%ldim(1:2)
        ldim_reproj(3)   = 1
        ldim_reprojs(1)  = 3*ldim_reproj(1) + 4*b
        ldim_reprojs(2)  = ldim_reproj(2) + 2*b
        ldim_reprojs(3)  = 1
        call reproj%new(ldim_reproj,   self%smpd)
        call reprojs%new(ldim_reprojs, self%smpd)
        !$omp parallel workshare
        reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1) = sum(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),dim=3)
        !$omp end parallel workshare
        call reproj%norm
        reprojs%rmat(b+1:b+ldim_reproj(1),b+1:b+ldim_reproj(2),1) = reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1)
        !$omp parallel workshare
        reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1) = sum(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),dim=2)
        !$omp end parallel workshare
        call reproj%norm
        reprojs%rmat(ldim_reproj(1)+2*b+1:2*(b+ldim_reproj(1)),b+1:b+ldim_reproj(2),1) = reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1)
        !$omp parallel workshare
        reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1) = sum(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),dim=1)
        !$omp end parallel workshare
        call reproj%norm
        reprojs%rmat(2*ldim_reproj(1)+3*b+1:3*(b+ldim_reproj(1)),b+1:b+ldim_reproj(2),1) = reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1)
        call reproj%kill
    end subroutine generate_orthogonal_reprojs

    ! FILTERS

    !>  \brief  acf calculates the autocorrelation function of an image
    subroutine acf( self )
        class(image), intent(inout) :: self
        if( .not. self%ft )then
            call self%fft()
        endif
        self%cmat = self%cmat*conjg(self%cmat) / sqrt(sum(csq(self%cmat)))
        call self%ifft()
    end subroutine acf

    !>  \brief ccf calculates the cross-correlation function between two images
    function ccf( self1, self2 ) result( cc )
        class(image), intent(inout) :: self1, self2
        type(image) :: cc
        if( .not. self1%ft )then
            call self1%fft()
        endif
        if( .not. self2%ft )then
            call self2%fft()
        endif
        cc      = self1
        cc%cmat = cc%cmat*conjg(self2%cmat)
        call cc%ifft()
    end function ccf

    !>  \brief guinier_bfac  generates the bfactor from the Guinier plot of the unfiltered volume
    function guinier_bfac( self, hp, lp ) result( bfac )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
        real, allocatable :: plot(:,:)
        integer :: fromk, tok, nk
        real    :: slope, intercept, corr, bfac
        plot  = self%guinier(.false.)
        fromk = self%get_find(hp)
        tok   = self%get_find(lp)
        nk    = tok-fromk+1
        call fit_straight_line(nk, plot(fromk:tok,:), slope, intercept, corr)
        bfac  = 4. * slope
        deallocate(plot)
    end function guinier_bfac

    !>  \brief guinier generates the Guinier plot for a volume, which should be unfiltered
    function guinier( self, verbose ) result( plot )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: verbose
        real, allocatable :: spec(:), plot(:,:)
        integer           :: lfny, k
        call self%spectrum('absreal',spec=spec)
        lfny = self%get_lfny(1)
        allocate( plot(lfny,2) )
        do k=1,lfny
            plot(k,1) = 1./(self%get_lp(k)**2.)
            plot(k,2) = log(spec(k))
            if( verbose ) write(logfhandle,'(A,1X,F8.4,1X,A,1X,F7.3)') '>>> RECIPROCAL SQUARE RES:', plot(k,1), '>>> LOG(ABS(REAL(F))):', plot(k,2)
        end do
        deallocate( spec )
    end function guinier

    !>  \brief spectrum generates the rotationally averaged spectrum of an image
    !>  keep serial
    subroutine spectrum( self, which, spec, norm )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real, allocatable, intent(inout) :: spec(:)
        logical, optional, intent(in)    :: norm
        real, allocatable :: counts(:)
        integer :: lims(3,2), phys(3), sh, lfny, h, k, l
        logical :: didft, nnorm
        nnorm = .true.
        if( present(norm) ) nnorm = norm
        didft = .false.
        if( which .ne. 'count' )then
            if( .not. self%ft )then
                call self%fft()
                didft = .true.
            endif
        endif
        lfny = self%get_lfny(1)
        if(allocated(spec))deallocate(spec)
        allocate( spec(lfny), counts(lfny) )
        spec   = 0.
        counts = 0.
        lims   = self%fit%loop_lims(2)
        select case(which)
        case('real')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + real(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('power')
            do l=lims(3,1),lims(3,2)
                do k=lims(2,1),lims(2,2)
                    do h=lims(1,1),lims(1,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + csq(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('sqrt')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + sqrt(csq(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('log')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + log(csq(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('absreal')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + abs(real(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('absimag')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + abs(aimag(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('abs')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + cabs(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('phase')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + phase_angle(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('count')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + 1.
                        counts(sh) = counts(sh)+1.
                    end do
                end do
            end do
        case DEFAULT
            write(logfhandle,*) 'Spectrum kind: ', trim(which)
            THROW_HARD('Unsupported spectrum kind; spectrum')
        end select
        if( which .ne. 'count' .and. nnorm )then
            where(counts > 0.)
                spec = spec/counts
            end where
        endif
        if( didft ) call self%ifft()
    end subroutine spectrum

    subroutine power_spectrum( self, spec )
        class(image), intent(in)    :: self
        real,         intent(inout) :: spec(fdim(self%ldim(1)) - 1)
        integer  :: counts(fdim(self%ldim(1)) - 1)
        real(dp) :: dspec( fdim(self%ldim(1)) - 1)
        integer  :: lims(3,2), phys(3), filtsz, h, k, l, sh
        if( .not.self%is_ft() ) THROW_HARD('Only for FTed images! power_spectrum')
        filtsz = fdim(self%ldim(1)) - 1
        dspec  = 0.d0
        counts = 0
        lims   = self%fit%loop_lims(2)
        if( self%is_2d() )then
            do k=lims(2,1),lims(2,2)
                do h=lims(1,1),lims(1,2)
                    sh = nint(hyp(h,k))
                    if( sh == 0 .or. sh > filtsz ) cycle
                    phys(1:2)  = self%fit%comp_addr_phys(h,k)
                    dspec(sh)  = dspec(sh)  + real(csq_fast(self%cmat(phys(1),phys(2),1)),dp)
                    counts(sh) = counts(sh) + 1
                end do
            end do
        else
            do l=lims(3,1),lims(3,2)
                do k=lims(2,1),lims(2,2)
                    do h=lims(1,1),lims(1,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > filtsz ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        dspec(sh)  = dspec(sh)  + real(csq_fast(self%cmat(phys(1),phys(2),phys(3))),dp)
                        counts(sh) = counts(sh) + 1
                    end do
                end do
            end do
        endif
        where(counts > 0)
            dspec = dspec / real(counts,dp)
        end where
        spec = real(dspec,kind=sp)
    end subroutine power_spectrum

    !> \brief apply_bfac  is for applying bfactor to an image
    subroutine apply_bfac( self, b )
        class(image), intent(inout) :: self
        real, intent(in)            :: b
        integer                     :: i,j,k,phys(3),lims(3,2)
        real                        :: wght, res
        logical                     :: didft
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lims = self%fit%loop_lims(2)
        if( self%wthreads )then
            !$omp parallel do collapse(3) default(shared) proc_bind(close)&
            !$omp private(k,j,i,res,phys,wght) schedule(static)
            do k=lims(3,1),lims(3,2)
                do j=lims(2,1),lims(2,2)
                    do i=lims(1,1),lims(1,2)
                        res = sqrt(real(k*k+j*j+i*i))/(real(self%ldim(1))*self%smpd) ! assuming square dimensions
                        phys = self%fit%comp_addr_phys([i,j,k])
                        wght = max(0.,exp(-(b/4.)*res*res))
                        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*wght
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do k=lims(3,1),lims(3,2)
                do j=lims(2,1),lims(2,2)
                    do i=lims(1,1),lims(1,2)
                        res = sqrt(real(k*k+j*j+i*i))/(real(self%ldim(1))*self%smpd) ! assuming square dimensions
                        phys = self%fit%comp_addr_phys([i,j,k])
                        wght = max(0.,exp(-(b/4.)*res*res))
                        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*wght
                    end do
                end do
            end do
        endif
        if( didft ) call self%ifft()
    end subroutine apply_bfac

    !> \brief bp  is for band-pass filtering an image with a cosine filter
    subroutine bp( self, hplim, lplim, width )
        class(image), intent(inout) :: self
        real, intent(in)            :: hplim, lplim
        real, intent(in), optional  :: width
        integer :: h, k, l, lims(3,2), phys(3)
        logical :: didft, dohp, dolp
        real    :: freq, hplim_freq, lplim_freq, wwidth, w
        wwidth =10.
        if( present(width) ) wwidth = width
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        dohp = abs(hplim) > TINY
        dolp = abs(lplim) > TINY
        hplim_freq = self%fit%get_find(1,hplim)
        lplim_freq = self%fit%get_find(1,lplim)
        lims = self%fit%loop_lims(2)
        if( self%wthreads )then
            !$omp parallel do private(h,k,l,freq,phys,w) default(shared)&
            !$omp collapse(3) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if( dohp )then
                            if(freq .lt. hplim_freq) then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .le. hplim_freq + wwidth) then
                                w = (1.-cos(((freq-hplim_freq)/wwidth)*pi))/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                        if( dolp )then
                            if(freq .gt. lplim_freq)then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .ge. lplim_freq - wwidth)then
                                w = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if( dohp )then
                            if(freq .lt. hplim_freq) then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .le. hplim_freq + wwidth) then
                                w = (1.-cos(((freq-hplim_freq)/wwidth)*pi))/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                        if( dolp )then
                            if(freq .gt. lplim_freq)then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .ge. lplim_freq - wwidth)then
                                w = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                    end do
                end do
            end do
        endif
        if( didft ) call self%ifft()
    end subroutine bp

     !> \brief lp  is for low-pass filtering an image with a cosine filter
    subroutine lp( self, find, width )
        class(image),   intent(inout) :: self
        integer,        intent(in)    :: find
        real, optional, intent(in)    :: width
        integer :: h, k, l, lims(3,2), phys(3)
        real    :: freq, wwidth, w
        wwidth =10.
        if( present(width) ) wwidth = width
        lims = self%fit%loop_lims(2)
        !$omp parallel do private(h,k,l,freq,phys,w) default(shared)&
        !$omp collapse(3) proc_bind(close)
        do l=lims(3,1),lims(3,2)
            do k=lims(2,1),lims(2,2)
                do h=lims(1,1),lims(1,2)
                    freq = hyp(h,k,l)
                    phys = self%comp_addr_phys([h,k,l])
                    if(freq .gt. find)then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else if(freq .ge. find - wwidth)then
                        w = (cos(((freq-(find-wwidth))/wwidth)*pi)+1.)/2.
                        self%cmat(phys(1),phys(2),phys(3)) = &
                            &self%cmat(phys(1),phys(2),phys(3)) * w
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine lp

    ! Band-pass gaussian filter, 2D images only
    subroutine bpgau2D( self, hp, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
        real    :: hp_fwhm(2), hp_halfinvsigsq(2), hpa
        real    :: lp_fwhm(2), lp_halfinvsigsq(2), lpa
        integer :: phys(2), lims(3,2), h,k
        logical :: l_hp, l_lp
        if(.not.self%ft) THROW_HARD('Input image must be in the reciprocal domain')
        if(.not.self%is_2d()) THROW_HARD('Input image must be two-dimensional')
        lims = self%fit%loop_lims(2)
        l_hp = .false.
        if( hp > TINY )then
            l_hp = .true.
            hp_fwhm         = hp / self%smpd / real(self%ldim(1:2))
            hp_halfinvsigsq = 0.5 * (PI * 2.0 * hp_fwhm / 2.35482)**2
        endif
        l_lp = .false.
        if( lp > TINY )then
            l_lp = .true.
            lp_fwhm         = lp / self%smpd / real(self%ldim(1:2))
            lp_halfinvsigsq = 0.5 * (PI * 2.0 * lp_fwhm / 2.35482)**2
        endif
        !$omp parallel do collapse(2) schedule(static) default(shared) proc_bind(close)&
        !$omp private(h,k,hpa,lpa,phys)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                phys = self%comp_addr_phys(h,k)
                if( l_hp )then
                    hpa  = real(h*h) * hp_halfinvsigsq(1) + real(k*k) * hp_halfinvsigsq(2)
                    self%cmat(phys(1),phys(2),1) = self%cmat(phys(1),phys(2),1) * (1.0-exp(-hpa))
                endif
                if( l_lp )then
                    lpa  = real(h*h) * lp_halfinvsigsq(1) + real(k*k) * lp_halfinvsigsq(2)
                    self%cmat(phys(1),phys(2),1) = self%cmat(phys(1),phys(2),1) * exp(-lpa)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine bpgau2D

    !> \brief bp  is for tophat band-pass filtering an image
    subroutine tophat( self, shell, halfwidth )
        class(image),   intent(inout) :: self
        integer,        intent(in)    :: shell
        real, optional, intent(in)    :: halfwidth
        integer :: h, k, l, lims(3,2), phys(3)
        logical :: didft
        real    :: freq, hplim_freq, lplim_freq, hwdth
        hwdth = 0.5
        if( present(halfwidth) ) hwdth = halfwidth
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        hplim_freq = max(0., real(shell) - hwdth)
        lplim_freq = real(shell) + hwdth
        lims = self%fit%loop_lims(2)
        if( self%wthreads )then
            !$omp parallel do private(h,k,l,freq,phys) default(shared)&
            !$omp collapse(3) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if(freq .lt. hplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        else if(freq .gt. lplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if(freq .lt. hplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        else if(freq .gt. lplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        endif
                    end do
                end do
            end do
        endif
        if( didft ) call self%ifft()
    end subroutine tophat

    !> \brief apply_filter_1  is for application of an arbitrary 1D filter function
    subroutine apply_filter_1( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        integer :: nyq, sh, h, k, l, lims(3,2)
        logical :: didft
        real    :: fwght, wzero
        nyq = size(filter)
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        wzero = maxval(filter)
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    ! set filter weight
                    if( sh > nyq )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    ! multiply with the weight
                    call self%mul([h,k,l], fwght)
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine apply_filter_1

    !> \brief apply_filter_2  is for application of an arbitrary 1D complex-valued filter function
    subroutine apply_filter_2( self, filter )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: filter(:)
        integer :: nyq, sh, h, k, l, lims(3,2)
        logical :: didft
        complex :: fwght, wzero, wnyq
        nyq = size(filter)
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        wzero = filter(1)
        wnyq  = cmplx(0.,0.)
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    ! set filter weight
                    if( sh > nyq )then
                        fwght = wnyq
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    ! multiply with the weight
                    call self%mul([h,k,l], fwght)
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine apply_filter_2

    !> \brief apply_filter_1  is for application of an arbitrary 1D filter function
    subroutine apply_filter_serial( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        integer :: nyq, sh, h, k, l, lims(3,2)
        real    :: fwght, wzero
        nyq   = size(filter)
        wzero = maxval(filter)
        lims  = self%fit%loop_lims(2)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    ! set filter weight
                    if( sh > nyq )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    ! multiply with the weight
                    call self%mul([h,k,l], fwght)
                end do
            end do
        end do
    end subroutine apply_filter_serial

    subroutine fcomps_below_noise_power_stats( self, noise_vol )
        class(image), intent(inout) :: self, noise_vol
        real, allocatable :: res(:)
        integer :: h, k, l, lims(3,2), phys(3), cnt, counts(self%get_filtsz())
        integer :: sh, counts_all(self%get_filtsz()), filtsz
        real    :: sig_pow, noise_pow
        lims  = self%fit%loop_lims(2)
        if( .not.self%is_ft() ) THROW_HARD('Instance need to be FTed')
        if( .not.self%is_ft() ) THROW_HARD('noise_vol need to be FTed')
        cnt        = 0
        counts     = 0
        counts_all = 0
        filtsz     = self%get_filtsz()
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh             = nint(hyp(h,k,l))
                    if( sh == 0 ) cycle
                    phys           = self%comp_addr_phys([h,k,l])
                    sig_pow        = csq(self%cmat(phys(1),phys(2),phys(3)))
                    noise_pow      = csq(noise_vol%cmat(phys(1),phys(2),phys(3)))
                    if( sh <= filtsz ) counts_all(sh) = counts_all(sh) + 1
                    if( noise_pow > sig_pow )then
                        call self%mul([h,k,l], 0.)
                        cnt = cnt + 1
                    else
                        if( sh <= filtsz ) counts(sh) = counts(sh) + 1
                    endif
                enddo
            enddo
        enddo
        res = self%get_res()
        do k=1,size(res)
            write(logfhandle,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:',&
            &res(k), '>>> % pow(FCOMPS) > pow(NOISE):', 100 * (real(counts(k)) / real(counts_all(k)))
        end do
        deallocate(res)
        write(logfhandle,'(a,f4.1)')&
        'fraction (%) of Fourier components with power below noise power zeroed: ',&
        &100 * (real(cnt) / product(self%ldim))
    end subroutine fcomps_below_noise_power_stats

    subroutine ran_phases_below_noise_power( self_even, self_odd )
        class(image), intent(inout) :: self_even, self_odd
        integer :: h, k, l, lims(3,2), phys(3)
        real    :: noise_pow, even_pow, odd_pow, phase
        complex :: diff
        lims  = self_even%fit%loop_lims(2)
        if( .not.self_even%is_ft() ) THROW_HARD('even image needs to be FTed')
        if( .not.self_odd%is_ft()  ) THROW_HARD('odd  image needs to be FTed')
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,diff,noise_pow,even_pow,odd_pow,phase)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    phys      = self_even%comp_addr_phys([h,k,l])
                    diff      = self_even%cmat(phys(1),phys(2),phys(3)) -&
                               &self_odd%cmat(phys(1),phys(2),phys(3))
                    noise_pow = csq_fast(diff)
                    even_pow  = csq_fast(self_even%cmat(phys(1),phys(2),phys(3)))
                    odd_pow   = csq_fast(self_odd%cmat(phys(1),phys(2),phys(3)))
                    if( noise_pow > even_pow .or. noise_pow > odd_pow )then
                        phase = ran3() * TWOPI
                        self_even%cmat(phys(1),phys(2),phys(3)) = mycabs(self_even%cmat(phys(1),phys(2),phys(3))) * cmplx(cos(phase), sin(phase))
                        phase = ran3() * TWOPI
                        self_odd%cmat(phys(1),phys(2),phys(3))  = mycabs(self_odd%cmat(phys(1),phys(2),phys(3)))  * cmplx(cos(phase), sin(phase))
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine ran_phases_below_noise_power

    ! self_even can be particle if is_ptcl==.true.
    ! then self_odd is CTF modulated reprojection
    subroutine whiten_noise_power( self_even, self_odd, is_ptcl )
        class(image), intent(inout) :: self_even, self_odd
        logical,      intent(in)    :: is_ptcl
        real(dp) :: counts(fdim(self_even%ldim(1)) - 1)
        real(dp) ::  dspec(fdim(self_even%ldim(1)) - 1)
        complex  :: diff
        integer  :: filtsz, h, k, l, sh, lims(3,2), phys(3)
        filtsz = fdim(self_even%ldim(1)) - 1
        dspec  = 0.d0
        counts = 0.d0
        lims   = self_even%fit%loop_lims(2)
        do l = lims(3,1),lims(3,2)
            do k = lims(2,1),lims(2,2)
                do h = lims(1,1),lims(1,2)
                    phys = self_even%fit%comp_addr_phys([h,k,l])
                    sh   = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > filtsz ) cycle
                    diff      = self_even%cmat(phys(1),phys(2),phys(3)) -&
                                &self_odd%cmat(phys(1),phys(2),phys(3))
                    dspec(sh)  = dspec(sh)  + csq_fast(dcmplx(diff))
                    counts(sh) = counts(sh) + 1.d0
                end do
            end do
        end do        
        if( is_ptcl )then
            where(counts > DTINY) dspec =         dspec / counts
        else
            where(counts > DTINY) dspec = 0.5d0 * dspec / counts ! 0.5 because of the e/o split
        endif
        do l = lims(3,1),lims(3,2)
            do k = lims(2,1),lims(2,2)
                do h = lims(1,1),lims(1,2)
                    phys = self_even%fit%comp_addr_phys([h,k,l])
                    sh   = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > filtsz ) cycle
                    if( dspec(sh) > DTINY )then
                                            self_even%cmat(phys(1),phys(2),phys(3)) = &
                                           &self_even%cmat(phys(1),phys(2),phys(3)) / real(dsqrt(dspec(sh)),kind=sp)
                        if( .not. is_ptcl ) self_odd%cmat( phys(1),phys(2),phys(3)) = &
                                           &self_odd%cmat( phys(1),phys(2),phys(3)) / real(dsqrt(dspec(sh)),kind=sp)
                    endif
                end do
            end do
        end do
    end subroutine whiten_noise_power

    ! This function performs image filtering by convolution
    ! with the 1D kernel filt. REAL SPACE.
    subroutine imfilter1(img,filt)
        class(image), intent(inout) :: img
        type(image) :: img_p
        real, allocatable   ::  rmat(:,:,:)
        real, intent(in)    :: filt(:)
        real, allocatable   :: shifted_filt(:)
        real, allocatable   :: rmat_t(:,:,:)
        integer :: ldim(3), sz_f(1), L1
        integer :: i, j, m
        ldim = img%get_ldim()
        sz_f = shape(filt)
        L1 = sz_f(1)
        allocate(shifted_filt(-(L1-1)/2:(L1-1)/2), source = filt)
        allocate(rmat_t(ldim(1),ldim(2),ldim(3)), source = 0.)
        call img_p%new([ldim(1)+L1-1,ldim(2)+L1-1,1],1.)
        call img%pad(img_p)
        rmat = img_p%get_rmat()
        !$omp parallel do collapse(2) default(shared) private(i,j,m)&
        !$omp schedule(static) proc_bind(close)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do m = -(L1-1)/2,(L1-1)/2
                    rmat_t(i,j,1) = rmat_t(i,j,1)+rmat(i+m+1,j+1,1)*shifted_filt(m)
                enddo
            end do
        end do
        !$omp end parallel do
        call img%set_rmat(rmat_t,.false.)
        deallocate(rmat, rmat_t, shifted_filt)
    end subroutine imfilter1

    ! This function performs image filtering by convolution
    ! with the 2D kernel filt. REAL SPACE.
    subroutine imfilter2(img,filt)
        class(image), intent(inout) :: img
        type(image) :: img_p
        real, allocatable   ::  rmat(:,:,:)
        real, intent(in)    :: filt(:,:)
        real, allocatable   :: shifted_filt(:,:)
        real, allocatable   :: rmat_t(:,:,:)
        integer :: ldim(3), sz_f(2), L1, L2
        integer :: i, j, m, n
        ldim = img%get_ldim()
        sz_f = shape(filt)
        L1 = sz_f(1)
        L2 = sz_f(2)
        allocate(shifted_filt(-(L1-1)/2:(L1-1)/2, -(L2-1)/2:(L2-1)/2), source = filt)
        allocate(rmat_t(ldim(1),ldim(2),ldim(3)), source = 0.)
        call img_p%new([ldim(1)+L1-1,ldim(2)+L2-1,1],1.)
        call img%pad(img_p)
        rmat = img_p%get_rmat()
        !$omp parallel do collapse(2) default(shared) private(i,j,m,n)&
        !$omp schedule(static) proc_bind(close)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do m = -(L1-1)/2,(L1-1)/2
                    do n = -(L2-1)/2,(L2-1)/2
                        rmat_t(i,j,1) = rmat_t(i,j,1)+rmat(i+m+1,j+n+1,1)*shifted_filt(m,n)
                    enddo
                enddo
            end do
        end do
        !$omp end parallel do
        call img%set_rmat(rmat_t,.false.)
        deallocate(rmat, rmat_t, shifted_filt)
    end subroutine imfilter2

    ! This function performs image filtering by convolution
    ! with the 3D kernel filt. REAL SPACE.
    subroutine imfilter3(img,filt)
        class(image), intent(inout) :: img
        type(image) :: img_p
        real, allocatable   ::  rmat(:,:,:)
        real, intent(in)    :: filt(:,:,:)
        real, allocatable   :: shifted_filt(:,:,:)
        real, allocatable   :: rmat_t(:,:,:)
        integer :: ldim(3), sz_f(3), L1, L2, L3
        integer :: i, j, k, m, n, o
        ldim = img%get_ldim()
        sz_f = shape(filt)
        L1 = sz_f(1)
        L2 = sz_f(2)
        L3 = sz_f(3)
        allocate(shifted_filt(-(L1-1)/2:(L1-1)/2,-(L2-1)/2:(L2-1)/2,-(L3-1)/2:(L3-1)/2), source = filt)
        allocate(rmat_t(ldim(1),ldim(2),ldim(3)), source = 0.)
        call img_p%new([ldim(1)+L1-1,ldim(2)+L2-1,ldim(3)+L3-1],1.)
        call img%pad(img_p)
        rmat = img_p%get_rmat()
        !$omp parallel do collapse(2) default(shared) private(i,j,k,m,n,o)&
        !$omp schedule(static) proc_bind(close)
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do k = 1, ldim(3)
                    do m = -(L1-1)/2,(L1-1)/2
                        do n = -(L2-1)/2,(L2-1)/2
                            do o = -(L3-1)/2,(L3-1)/2
                                rmat_t(i,j,k) = rmat_t(i,j,k)+rmat(i+m+1,j+n+1,k+o+1)*shifted_filt(m,n,o)
                            enddo
                        enddo
                    enddo
                enddo
            end do
        end do
        !$omp end parallel do
        call img%set_rmat(rmat_t,.false.)
        deallocate(rmat, rmat_t, shifted_filt)
    end subroutine imfilter3

    !> \brief phase_rand  is for randomzing the phases of the FT of an image from lp and out
    subroutine phase_rand( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
        integer                     :: h,k,l,phys(3),lims(3,2)
        logical                     :: didft
        real                        :: freq,lp_freq, amp,phase
        real, parameter             :: errfrac=0.5
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lp_freq = real(self%fit%get_find(1,lp)) ! assuming square 4 now
        lims    = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,freq,phys,amp,phase)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    freq = hyp(h,k,l)
                    if(freq .gt. lp_freq)then
                        phys  = self%fit%comp_addr_phys([h,k,l])
                        amp   = mycabs(self%cmat(phys(1),phys(2),phys(3)))
                        phase = ran3() * TWOPI
                        self%cmat(phys(1),phys(2),phys(3)) = amp * cmplx(cos(phase), sin(phase))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine phase_rand

    !> \brief hannw a constructor that constructs an antialiasing Hanning window
    !! \param oshoot_in overshoot
    !! \return  w Hanning window
    !!
    function hannw( self, oshoot_in ) result( w )
        class(image), intent(inout) :: self
        real, intent(in), optional  :: oshoot_in
        integer                     :: lims(3,2), k, kmax, maxl
        type(winfuns)               :: wfuns
        character(len=STDLEN)       :: wstr
        real, allocatable           :: w(:)
        real                        :: oshoot
        oshoot = 0.3
        if( present(oshoot_in) ) oshoot = oshoot_in
        lims = self%loop_lims(2)
        maxl = maxval(lims)
        kmax = maxl+int(oshoot*real(maxl))
        allocate( w(kmax) )
        wstr = 'hann'
        wfuns = winfuns(wstr, real(kmax), 2.)
        do k=1,kmax
            w(k) = wfuns%eval_apod(real(k))
        end do
    end function hannw

    !>  \brief average and median filtering in real-space
    subroutine real_space_filter( self, winsz, which )
        class(image),     intent(inout) :: self
        integer,          intent(in)    :: winsz
        character(len=*), intent(in)    :: which
        real, allocatable     :: pixels(:), wfvals(:)
        integer               :: n, i, j, k, cnt, npix
        real                  :: rn, wfun(-winsz:winsz), norm, avg, sdev
        type(winfuns)         :: fwin
        character(len=STDLEN) :: wstr
        type(image)           :: img_filt
        real                  :: k2(3,3), k3(3,3,3) !laplacian kernels (2D-3D)
        ! check the number of pixels in window
        if( self%is_3d() )then
            npix = (2*winsz+1)**3
        else
            npix = (2*winsz+1)**2
        endif
        pixels = self%win2arr(1, 1, 1, winsz)
        n = size(pixels)
        rn = real(n)
        allocate(wfvals(n))
        ! make the window function
        wstr = 'bman'
        fwin = winfuns(wstr, real(WINSZ), 1.0)
        ! sample the window function
        do i=-winsz,winsz
            wfun(i) = fwin%eval_apod(real(i))
        end do
        ! memoize wfun vals & normalisation constant
        norm = 0.
        cnt  = 0
        if( self%ldim(3) == 1 )then
            do i=-winsz,winsz
                do j=-winsz,winsz
                    cnt = cnt + 1
                    wfvals(cnt) = wfun(i) * wfun(j)
                    norm = norm + wfvals(cnt)
                end do
            end do
        else
            if (which == 'bman')then
                do i=-winsz,winsz
                    do j=-winsz,winsz
                        do k=-winsz,winsz
                            cnt = cnt + 1
                            wfvals(cnt) = wfun(i) * wfun(j) * wfun(k)
                            norm = norm + wfvals(cnt)
                        end do
                    end do
                end do
            else
                !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)&
                !$omp reduction(+:norm)
                do i=-winsz,winsz
                    do j=-winsz,winsz
                        do k=-winsz,winsz
                            norm = norm + wfun(i) * wfun(j) * wfun(k)
                        end do
                    end do
                end do
                !$omp end parallel do
            endif
        endif
        ! make the output image
        call img_filt%new(self%ldim, self%smpd)
        ! filter
        if( self%ldim(3) == 1 )then
            select case(which)
            case('median')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = median_nocopy(pixels)
                    end do
                end do
                !$omp end parallel do
            case('average')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = sum(pixels)/rn
                    end do
                end do
                !$omp end parallel do
            case('stdev')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels,avg,sdev) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        call avg_sdev(pixels, avg, sdev)
                        img_filt%rmat(i,j,1) = sdev
                    end do
                end do
                !$omp end parallel do
            case('bman')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = sum(pixels * wfvals) / norm
                    end do
                end do
                !$omp end parallel do
            case('NLmean')
                call self%NLmean2D()
                img_filt%rmat = self%rmat
            case('laplacian')
                k2 = (1./8.)*reshape([0.,1.,0.,1.,-4., 1., 0., 1., 0.], [3,3])
                call self%imfilter(k2)
                img_filt = self
            case DEFAULT
                THROW_HARD('unknown filter type; real_space_filter')
            end select
        else !3D
            select case(which)
            case('median')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = median_nocopy(pixels)
                        end do
                    end do
                end do
                !$omp end parallel do
            case('average')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) !schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = sum(pixels)/rn
                        end do
                    end do
                end do
                !$omp end parallel do
            case('stdev')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels,avg,sdev) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            call avg_sdev(pixels, avg, sdev)
                            img_filt%rmat(i,j,k) = sdev
                        end do
                    end do
                end do
                !$omp end parallel do
            case('bman')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = sum(pixels * wfvals) / norm
                        end do
                    end do
                end do
                !$omp end parallel do
            case('laplacian')
                k3 = (1./12.)*reshape([0.,0.,0., 0.,1.,0., 0.,0.,0.,&
                &                     0.,1.,0., 1.,-6.,1., 0.,1.,0.,0.,0.,0., 0.,1.,0., 0.,0.,0.], [3,3,3])
                call self%imfilter(k3)
                img_filt = self
            case DEFAULT
                THROW_HARD('unknown filter type; real_space_filter')
            end select
        endif
        call self%copy(img_filt)
        call img_filt%kill()
    end subroutine real_space_filter

    !>  \brief Non-local mean filter, don't touch default parameters here. This routine is being actively used
    subroutine NLmean2D( self, msk, sdev_noise )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: msk
        real, optional, intent(in)    :: sdev_noise
        real,  allocatable :: rmat_pad(:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 3   ! good if use Euclidean distance
        integer, parameter :: CFR_BOX = 10  ! as suggested in the paper, use a box 21x21
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW)
        real    :: z, sigma, h, h_sq, avg, mmsk
        integer :: i, j, m, n, pad, ithr
        if( self%is_3d() ) THROW_HARD('2D images only; NLmean2D')
        if( self%ft )      THROW_HARD('Real space only;NLmean2D')
        mmsk = real(self%ldim(1)) / 2. - real(DIM_SW)
        if( present(msk) ) mmsk = msk
        if( present(sdev_noise) )then
            if( sdev_noise > SMALL )then
                sigma = sdev_noise
            else
                sigma = self%noisesdev(mmsk) ! estimation of noise
            endif
        else
            sigma = self%noisesdev(mmsk) ! estimation of noise
        endif
        pad   = CFR_BOX + 2
        h     = 4.*sigma
        h_sq  = h**2.
        avg   = sum(self%rmat(:self%ldim(1),:self%ldim(2),1)) / real(product(self%ldim))
        allocate(rmat_threads(nthr_glob,self%ldim(1),self%ldim(2),1), source=0.)
        allocate(rmat_pad(-pad:self%ldim(1)+pad,-pad:self%ldim(2)+pad), source=avg)
        rmat_pad(1:self%ldim(1),1:self%ldim(2)) = self%rmat(1:self%ldim(1),1:self%ldim(2),1)
        !$omp parallel do schedule(static) default(shared) private(m,n,ithr,sw_px,exponentials,i,j,z)&
        !$omp proc_bind(close) firstprivate(rmat_pad) collapse(2)
        do n = 1,self%ldim(2)
            do m = 1,self%ldim(1)
                ithr  = omp_get_thread_num() + 1
                sw_px = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1)
                exponentials = 0.
                do j = -CFR_BOX,CFR_BOX
                    do i = -CFR_BOX,CFR_BOX
                      exponentials(i,j) = &
                      & exp( -sum( (sw_px - rmat_pad(m+i:m+i+DIM_SW-1, n+j:n+j+DIM_SW-1))**2. )/h_sq) ! Euclidean norm
                  enddo
                enddo
                z = sum(exponentials)
                if( z < 0.0000001 ) cycle
                rmat_threads(ithr,m,n,1) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX)) / z
            enddo
        enddo
        !$omp end parallel do
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = sum(rmat_threads, dim=1)
    end subroutine NLmean2D

    subroutine NLmean2D_eo( even, odd, avg )
        class(image), intent(inout) :: even, odd, avg
        type(image)        :: noise, noise_var
        real,  allocatable :: rmat_pad(:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 1
        integer, parameter :: CFR_BOX = 3
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW)
        real    :: z, pix_avg, sigma2t2
        integer :: i, j, m, n, pad, ithr
        if( even%is_3d() ) THROW_HARD('2D images only; NLmean2D_eo')
        if( even%ft )      THROW_HARD('Real space only;NLmean2D_eo')
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var(noise_var)
        call noise_var%norm([5.,2.])
        call avg%copy(even)
        call avg%add(odd)
        call avg%mul(0.5)
        pad     = CFR_BOX + 2
        pix_avg = sum(avg%rmat(:avg%ldim(1),:avg%ldim(2),1)) / real(product(avg%ldim))
        allocate(rmat_threads(nthr_glob,avg%ldim(1),avg%ldim(2),1),   source=0.)
        allocate(rmat_pad(-pad:avg%ldim(1)+pad,-pad:avg%ldim(2)+pad), source=pix_avg)
        rmat_pad(1:avg%ldim(1),1:avg%ldim(2)) = avg%rmat(1:avg%ldim(1),1:avg%ldim(2),1)
        do n = 1,avg%ldim(2)
            do m = 1,avg%ldim(1)
                ithr         = omp_get_thread_num() + 1
                sw_px        = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1)
                exponentials = 0.
                sigma2t2     = 2. * noise_var%rmat(m,n,1)
                do j = -CFR_BOX,CFR_BOX
                    do i = -CFR_BOX,CFR_BOX
                      exponentials(i,j) = &
                      & exp( -sum( (sw_px - rmat_pad(m+i:m+i+DIM_SW-1, n+j:n+j+DIM_SW-1))**2. )/sigma2t2) ! Euclidean norm
                  enddo
                enddo
                z = sum(exponentials)
                if( z < 0.0000001 ) cycle
                rmat_threads(ithr,m,n,1) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX)) / z
            enddo
        enddo
        avg%rmat(:avg%ldim(1),:avg%ldim(2),:avg%ldim(3)) = sum(rmat_threads, dim=1)
        call noise%kill
        call noise_var%kill
    end subroutine NLmean2D_eo

    !>  \brief Non-local mean filter
    subroutine NLmean3D( self, msk, sdev_noise )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: msk
        real, optional, intent(in)    :: sdev_noise
        real,  allocatable :: rmat_pad(:,:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 3
        integer, parameter :: CFR_BOX = 10 ! as suggested in the paper, use a box 21x21
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW,DIM_SW)
        real    :: z, sigma, h, h_sq, mmsk, avg
        integer :: i, j, k, m, n, o, pad, ithr
        if( self%is_2d() ) THROW_HARD('3D images only; NLmean3D')
        if( self%ft )      THROW_HARD('Real space only; 3DNLmean3D')
        mmsk = real(self%ldim(1)) / 2. - real(DIM_SW)
        if( present(msk) ) mmsk = msk
        if( present(sdev_noise) )then
            if( sdev_noise > SMALL )then
                sigma = sdev_noise
            else
                sigma = self%noisesdev(mmsk) ! estimation of noise
            endif
        else
            sigma = self%noisesdev(mmsk) ! estimation of noise
        endif
        pad  = CFR_BOX + 2
        h    = 4.*sigma
        h_sq = h**2.
        avg  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        allocate(rmat_threads(nthr_glob,self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(rmat_pad(-pad:self%ldim(1)+pad,-pad:self%ldim(2)+pad,-pad:self%ldim(3)+pad), source=avg)
        rmat_pad(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        !$omp parallel do schedule(static) default(shared) private(m,n,o,ithr,sw_px,exponentials,i,j,k,z)&
        !$omp proc_bind(close) firstprivate(rmat_pad) collapse(3)
        do o = 1, self%ldim(3)
            do n = 1, self%ldim(2)
                do m = 1, self%ldim(1)
                    ithr  = omp_get_thread_num() + 1
                    sw_px = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1,o:o+DIM_SW-1)
                    exponentials = 0.
                    do k = -CFR_BOX,CFR_BOX
                        do j = -CFR_BOX,CFR_BOX
                            do i = -CFR_BOX,CFR_BOX
                                exponentials(i,j,k) = &
                                & exp(-sum((sw_px - rmat_pad(m+i:m+i+DIM_SW-1,n+j:n+j+DIM_SW-1,o+k:o+k+DIM_SW-1))**2.)/h_sq) ! euclidean norm
                            enddo
                        enddo
                    enddo
                    z = sum(exponentials)
                    if( z < 0.0000001 ) cycle
                    rmat_threads(ithr,m,n,o) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX,o-CFR_BOX:o+CFR_BOX)) / z
                enddo
            enddo
        enddo
        !$omp end parallel do
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = sum(rmat_threads, dim=1)
    end subroutine NLmean3D

    !>  \brief Non-local mean filter
    subroutine NLmean3D_eo( even, odd, avg )
        class(image), intent(inout) :: even, odd, avg
        real,  allocatable :: rmat_pad(:,:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 1
        integer, parameter :: CFR_BOX = 3
        type(image) :: noise, noise_var
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW,DIM_SW)
        real    :: z, sigma2t2, pixavg
        integer :: i, j, k, m, n, o, pad, ithr
        if( even%is_2d() ) THROW_HARD('3D images only; NLmean3D_eo')
        if( even%ft )      THROW_HARD('Real space only; 3DNLmean3D_eo')
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var3D(noise_var)
        call noise_var%norm([5.,2.])
        call avg%copy(even)
        call avg%add(odd)
        call avg%mul(0.5)
        pad     = CFR_BOX + 2
        pixavg  = sum(avg%rmat(:avg%ldim(1),:avg%ldim(2),:avg%ldim(3))) / real(product(avg%ldim))
        allocate(rmat_threads(nthr_glob,avg%ldim(1),avg%ldim(2),avg%ldim(3)), source=0.)
        allocate(rmat_pad(-pad:avg%ldim(1)+pad,-pad:avg%ldim(2)+pad,-pad:avg%ldim(3)+pad), source=pixavg)
        rmat_pad(1:avg%ldim(1),1:avg%ldim(2),1:avg%ldim(3)) = avg%rmat(1:avg%ldim(1),1:avg%ldim(2),1:avg%ldim(3))
        !$omp parallel do schedule(static) default(shared) private(m,n,o,ithr,sw_px,exponentials,sigma2t2,i,j,k,z)&
        !$omp proc_bind(close) firstprivate(rmat_pad) collapse(3)
        do o = 1, avg%ldim(3)
            do n = 1, avg%ldim(2)
                do m = 1, avg%ldim(1)
                    ithr         = omp_get_thread_num() + 1
                    sw_px        = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1,o:o+DIM_SW-1)
                    exponentials = 0.
                    sigma2t2     = 2. * noise_var%rmat(m,n,o)
                    do k = -CFR_BOX,CFR_BOX
                        do j = -CFR_BOX,CFR_BOX
                            do i = -CFR_BOX,CFR_BOX
                                exponentials(i,j,k) = &
                                & exp(-sum((sw_px - rmat_pad(m+i:m+i+DIM_SW-1,n+j:n+j+DIM_SW-1,o+k:o+k+DIM_SW-1))**2.)/sigma2t2) ! euclidean norm
                            enddo
                        enddo
                    enddo
                    z = sum(exponentials)
                    if( z < 0.0000001 ) cycle
                    rmat_threads(ithr,m,n,o) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX,o-CFR_BOX:o+CFR_BOX)) / z
                enddo
            enddo
        enddo
        !$omp end parallel do
        avg%rmat(:avg%ldim(1),:avg%ldim(2),:avg%ldim(3)) = sum(rmat_threads, dim=1)
    end subroutine NLmean3D_eo

    ! uniform noise variance (sigma2 constant, set to 5)
    subroutine ICM2D( self, lambda, verbose )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: verbose
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        integer     :: n_8(3,8), nsz, i, j, k, m, n
        real        :: pot_term, pix, min, proba, sigma2t2, x, xmin, transl_tab(NQUANTA), eucl, y, sy, syy, diff, rnsz
        type(image) :: self_prev
        logical     :: l_verbose
        if( self%is_3d() ) THROW_HARD('2D images only; ICM')
        if( self%ft )      THROW_HARD('Real space only; ICM')
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        call self%quantize_fwd(NQUANTA, transl_tab)
        call self_prev%copy(self)
        sigma2t2 = 10.
        do i = 1, MAXITS
            do m = 1,self%ldim(2)
                do n = 1,self%ldim(1)
                    pix  = self_prev%rmat(n,m,1)
                    call neigh_8(self%ldim, [n,m,1], n_8, nsz)
                    rnsz = real(nsz)
                    sy   = 0.
                    syy  = 0.
                    do j = 1, nsz
                        y   = self_prev%rmat(n_8(1,j),n_8(2,j),1)
                        sy  = sy  + y
                        syy = syy + y*y
                    end do
                    pot_term = syy
                    min      = (pix * pix) / sigma2t2 + lambda * pot_term
                    xmin     = 0.
                    ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                    do k = 1,NQUANTA - 1
                        x        = real(k)
                        pot_term = syy + rnsz*x*x - 2.0*sy*x
                        diff     = pix - x
                        proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                        if( min > proba )then
                            min  = proba
                            xmin = x
                        endif
                    end do
                    self%rmat(n,m,1) = xmin
                end do
            end do
            if( l_verbose )then
                eucl = self%euclid_norm(self_prev)
                ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl
            endif
            if( i < MAXITS ) self_prev%rmat = self%rmat
        end do
        call self%quantize_bwd(NQUANTA, transl_tab)
        call self_prev%kill
    end subroutine ICM2D

    ! nonuniform
    subroutine ICM2D_eo( even, odd, lambda, verbose )
        class(image),      intent(inout) :: even, odd
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: verbose
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        type(image) :: even_prev, odd_prev, noise, noise_var
        integer     :: n_8(3,8), nsz, i, j, k, m, n
        real        :: transl_tab_even(NQUANTA), transl_tab_odd(NQUANTA)
        real        :: pot_term(2), pix(2), minv(2), proba(2), sigma2t2
        real        :: x, xmin(2), eucl, y(2), sy(2), syy(2), diff(2), rnsz
        logical     :: l_verbose
        if( even%is_3d() ) THROW_HARD('2D images only; ICM2D_eo')
        if( even%ft )      THROW_HARD('Real space only; ICM2D_eo')
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var(noise_var)
        call noise_var%norm([5.,2.])
        call even%quantize_fwd(NQUANTA, transl_tab_even)
        call odd%quantize_fwd(NQUANTA, transl_tab_odd)
        call even_prev%copy(even)
        call odd_prev%copy(odd)
        do i = 1, MAXITS
            do m = 1,even%ldim(2)
                do n = 1,even%ldim(1)
                    pix(1)   = odd_prev%rmat(n,m,1)
                    pix(2)   = even_prev%rmat(n,m,1)
                    sigma2t2 = 2. * noise_var%rmat(n,m,1)
                    call neigh_8(even%ldim, [n,m,1], n_8, nsz)
                    rnsz = real(nsz)
                    sy   = 0.
                    syy  = 0.
                    do j = 1, nsz
                        y(1) = odd_prev%rmat(n_8(1,j),n_8(2,j),1)
                        y(2) = even_prev%rmat(n_8(1,j),n_8(2,j),1)
                        sy   = sy  + y
                        syy  = syy + y*y
                    end do
                    xmin     = 0.
                    pot_term = syy ! x=0.
                    minv     = (pix * pix) / sigma2t2 + lambda * pot_term
                    ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                    do k = 1,NQUANTA - 1
                        x        = real(k)
                        pot_term = syy + rnsz*x*x - 2.0*sy*x
                        diff     = pix - x
                        proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                        if( minv(1) > proba(1) )then
                            minv(1) = proba(1)
                            xmin(1) = x
                        endif
                        if( minv(2) > proba(2) )then
                            minv(2) = proba(2)
                            xmin(2) = x
                        endif
                    end do
                    odd%rmat(n,m,1)  = xmin(1)
                    even%rmat(n,m,1) = xmin(2)
                end do
            end do
            if( l_verbose )then
                eucl = (even%euclid_norm(even_prev) + odd%euclid_norm(odd_prev)) / 2.
                ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl
            endif
            if( i < MAXITS)then
                even_prev%rmat = even%rmat
                odd_prev%rmat  = odd%rmat
            endif
        end do
        call even%quantize_bwd(NQUANTA, transl_tab_even)
        call odd%quantize_bwd(NQUANTA, transl_tab_odd)
        call even_prev%kill
        call odd_prev%kill
        call noise%kill
        call noise_var%kill
    end subroutine ICM2D_eo

    ! uniform noise variance (sigma2 constant, set to 5)
    subroutine ICM3D( self, lambda )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lambda
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        type(image) :: self_prev
        integer     :: n_4(3,6), nsz, i, j, k, m, n, l
        real        :: pot_term, pix, min, proba, sigma2t2, x, xmin, transl_tab(NQUANTA), eucl,y, sy, syy, diff, rnsz
        if( self%is_2d() ) THROW_HARD('3D images only; ICM')
        if( self%ft )      THROW_HARD('Real space only; ICM')
        call self%quantize_fwd(NQUANTA, transl_tab)
        call self_prev%copy(self)
        sigma2t2 = 10.
        do i = 1, MAXITS
            !$omp parallel do schedule(static) default(shared) private(n,m,l,pix,n_4,nsz,rnsz,pot_term,diff,j,xmin,min,k,x,proba,y,sy,syy)&
            !$omp proc_bind(close) collapse(3)
            do l = 1,self%ldim(3)
                do m = 1,self%ldim(2)
                    do n = 1,self%ldim(1)
                        pix  = self_prev%rmat(n,m,l)
                        call neigh_4_3D(self%ldim, [n,m,l], n_4, nsz)
                        rnsz = real(nsz)
                        sy   = 0.
                        syy  = 0.
                        do j = 1, nsz
                            y   = self_prev%rmat(n_4(1,j),n_4(2,j),n_4(3,j))
                            sy  = sy  + y
                            syy = syy + y*y
                        end do
                        pot_term = syy
                        min      = (pix * pix) / sigma2t2 + lambda * pot_term
                        xmin     = 0.
                        ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                        do k = 1,NQUANTA - 1
                            x        = real(k)
                            pot_term = syy + rnsz*x*x - 2.0*sy*x
                            diff     = pix - x
                            proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                            if( min > proba )then
                                min  = proba
                                xmin = x
                            endif
                        end do
                        self%rmat(n,m,l) = xmin
                    end do
                end do
            end do
            !$omp end parallel do
            eucl = self%euclid_norm(self_prev)
            ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl 
            call self_prev%copy(self)
        end do
        call self%quantize_bwd(NQUANTA, transl_tab)
        call self_prev%kill
    end subroutine ICM3D

    ! nonuniform
    subroutine ICM3D_eo( even, odd, lambda, l_msk )
        class(image),      intent(inout) :: even, odd
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: l_msk(even%ldim(1),even%ldim(2),even%ldim(3))
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        type(image) :: even_prev, odd_prev, noise, noise_var
        integer     :: n_4(3,6), nsz, i, j, k, m, n, l
        real        :: transl_tab_even(NQUANTA), transl_tab_odd(NQUANTA)
        real        :: sy(2), syy(2), y(2), pot_term(2), pix(2), minv(2)
        real        :: proba(2), sigma2t2, x, xmin(2), eucl, diff(2), rnsz
        if( even%is_2d() ) THROW_HARD('3D images only; ICM3D_eo')
        if( even%ft      ) THROW_HARD('Real space only; ICM3D_eo')
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var3D(noise_var)
        call noise_var%norm([5.,2.])
        if( present(l_msk) )then
            call even%quantize_fwd(NQUANTA, transl_tab_even, l_msk)
            call odd%quantize_fwd(NQUANTA, transl_tab_odd, l_msk)
        else
            call even%quantize_fwd(NQUANTA, transl_tab_even)
            call odd%quantize_fwd(NQUANTA, transl_tab_odd)
        endif
        call even_prev%copy(even)
        call odd_prev%copy(odd)
        do i = 1, MAXITS
            !$omp parallel do private(n,m,l,pix,sigma2t2,n_4,nsz,rnsz,pot_term,diff,j,xmin,minv,k,x,proba,y,syy,sy)&
            !$omp proc_bind(close) collapse(3) schedule(static) default(shared)
            do l = 1,even%ldim(3)
                do m = 1,even%ldim(2)
                    do n = 1,even%ldim(1)
                        pix(1)   = odd_prev%rmat(n,m,l)
                        pix(2)   = even_prev%rmat(n,m,l)
                        sigma2t2 = 2. * noise_var%rmat(n,m,l)
                        call neigh_4_3D(even%ldim, [n,m,l], n_4, nsz)
                        rnsz = real(nsz)
                        ! x: central pixel/candidate value
                        ! y: nsz neighbours, constants
                        ! pot_term = SUMi((x-yi)**2) = SUMi(yi**2) + nsz.x**2 - 2.x.SUMi(yi)
                        sy  = 0.
                        syy = 0.
                        do j = 1, nsz
                            y(1) = odd_prev%rmat( n_4(1,j),n_4(2,j),n_4(3,j))
                            y(2) = even_prev%rmat(n_4(1,j),n_4(2,j),n_4(3,j))
                            sy   = sy  + y
                            syy  = syy + y * y
                        end do
                        xmin     = 0.
                        pot_term = syy ! x=0
                        minv     = (pix * pix) / sigma2t2 + lambda * pot_term
                        ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                        do k = 1,NQUANTA - 1
                            x        = real(k)
                            pot_term = syy + rnsz*x*x - 2.0*sy*x
                            diff     = pix - x
                            proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                            if( minv(1) > proba(1) )then
                                minv(1) = proba(1)
                                xmin(1) = x
                            endif
                            if( minv(2) > proba(2) )then
                                minv(2) = proba(2)
                                xmin(2) = x
                            endif
                        end do
                        odd%rmat(n,m,l)  = xmin(1)
                        even%rmat(n,m,l) = xmin(2)
                    end do
                end do
            end do
            !$omp end parallel do
            eucl = (even%euclid_norm(even_prev) + odd%euclid_norm(odd_prev)) / 2.
            ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl
            call even_prev%copy(even)
            call odd_prev%copy(odd)
        end do
        call even%quantize_bwd(NQUANTA, transl_tab_even)
        call odd%quantize_bwd(NQUANTA, transl_tab_odd)
        call even_prev%kill
        call odd_prev%kill
        call noise%kill
        call noise_var%kill
    end subroutine ICM3D_eo

    ! generate gray level co-occurence matrix
    subroutine GLCM( self, nquanta, pmat )
        class(image), intent(in)    :: self
        integer,      intent(in)    :: nquanta
        real,         intent(inout) :: pmat(nquanta,nquanta)
        type(image) :: img_q
        real        :: transl_tab(nquanta)
        integer     :: m, n, j, n_8(3,8), nsz, ipix, jpix
        call img_q%copy(self)
        call img_q%quantize_fwd(nquanta, transl_tab)
        do m = 1,img_q%ldim(2)
            do n = 1,img_q%ldim(1)
                ipix = nint(img_q%rmat(n,m,1) + 1.)
                call neigh_8(img_q%ldim, [n,m,1], n_8, nsz)
                do j = 1, nsz
                    jpix = nint(img_q%rmat(n_8(1,j),n_8(2,j),1) + 1.)
                    if( jpix == ipix ) pmat(ipix,jpix) = pmat(ipix,jpix) + 1.
                end do
            end do
        end do
        call normalize_minmax(pmat)
        call img_q%kill
    end subroutine GLCM

    ! CALCULATORS

    !> \brief stats  is for providing foreground/background statistics
    subroutine stats_1( self, which, ave, sdev, maxv, minv, msk, med, errout )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: ave, sdev, maxv, minv
        real,    optional, intent(in)    :: msk
        real,    optional, intent(out)   :: med
        logical, optional, intent(out)   :: errout
        integer           :: i, j, k, npix, minlen
        real              :: ci, cj, ck, mskrad, e, var
        logical           :: err, didft, background
        real, allocatable :: pixels(:)
        ! FT
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ! 2d/3d
        if( self%ldim(3) > 1 )then
            minlen = minval(self%ldim)
        else
            minlen = minval(self%ldim(1:2))
        endif
        ! mask
        if( present(msk) )then
            mskrad = msk
        else
            mskrad = real(minlen)/2.
        endif
        ! back/foreground
        if( which.eq.'background' )then
            background = .true.
        else if( which.eq.'foreground' )then
            background = .false.
        else
            THROW_HARD('unrecognized parameter: which; stats_1')
        endif
        allocate( pixels(product(self%ldim)) )
        pixels = 0.
        npix   = 0
        if( self%ldim(3) > 1 )then
            ! 3d
            ci = -real(self%ldim(1))/2.
            do i=1,self%ldim(1)
                cj = -real(self%ldim(2))/2.
                do j=1,self%ldim(2)
                    ck = -real(self%ldim(3))/2.
                    do k=1,self%ldim(3)
                        e = hardedge(ci,cj,ck,mskrad)
                        if( background )then
                            if( e < 0.5 )then
                                npix = npix+1
                                pixels(npix) = self%rmat(i,j,k)
                            endif
                        else
                            if( e > 0.5 )then
                                npix = npix+1
                                pixels(npix) = self%rmat(i,j,k)
                            endif
                        endif
                        ck = ck + 1.
                    end do
                    cj = cj + 1.
                end do
                ci = ci + 1.
            end do
        else
            ! 2d
            ci = -real(self%ldim(1))/2.
            do i=1,self%ldim(1)
                cj = -real(self%ldim(2))/2.
                do j=1,self%ldim(2)
                    e = hardedge(ci,cj,mskrad)
                    if( background )then
                        if( e < 0.5 )then
                            npix = npix+1
                            pixels(npix) = self%rmat(i,j,1)
                        endif
                    else
                        if( e > 0.5 )then
                            npix = npix+1
                            pixels(npix) = self%rmat(i,j,1)
                        endif
                    endif
                    cj = cj + 1.
                end do
                ci = ci + 1.
            end do
        endif
        maxv = maxval(pixels(:npix))
        minv = minval(pixels(:npix))
        if(npix>1)then
            call moment( pixels(:npix), ave, sdev, var, err )
            if( present(med) ) med  = median_nocopy(pixels(:npix))
        end if
        deallocate( pixels )
        if( present(errout) )then
            errout = err
        else
            if( err ) THROW_WARN('variance zero; stats_1')
        endif
        if( didft ) call self%fft()
    end subroutine stats_1

    !> \brief stats  is for providing within mask statistics
    subroutine stats_2( self, ave, sdev, maxv, minv, mskimg, med, errout )
        class(image),           intent(inout) :: self
        real,                   intent(out)   :: ave, sdev, maxv, minv
        class(image), optional, intent(in)    :: mskimg
        real,         optional, intent(out)   :: med
        logical,      optional, intent(out)   :: errout
        real              :: var
        logical           :: err
        real, allocatable :: pixels(:)
        ! FT
        if( self%ft ) THROW_HARD('not for FTed imgs; stats_2')
        if( present(mskimg) )then
            pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                &mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.95 )
        else
          pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), .true.)
        endif
        maxv = maxval(pixels)
        minv = minval(pixels)
        call moment( pixels, ave, sdev, var, err )
        if( present(med) ) med  = median_nocopy(pixels)
        deallocate( pixels )
        if( present(errout) )then
            errout = err
        else
            if( err ) THROW_WARN('variance zero; stats_2')
        endif
    end subroutine stats_2

    !>  \brief minmax to get the minimum and maximum values in an image
    !! \return  mm 2D element (minimum , maximum)
    function minmax( self, radius )result( mm )
        class(image),   intent(in) :: self
        real, optional, intent(in) :: radius
        real    :: mm(2), radsq
        integer :: c(3),i,j,k,dksq,djsq,dsq
        if( present(radius) )then
            radsq = radius**2
            c     = self%ldim/2+1
            if( self%ldim(3)==1 ) c(3)=1
            mm    = [huge(0.),-huge(0.)]
            do k = 1,self%ldim(3)
                dksq = (k-c(3))**2
                do j = 1,self%ldim(2)
                    djsq = dksq + (j-c(2))**2
                    do i = 1,self%ldim(1)
                        dsq = djsq + (i-c(1))**2
                        if( real(dsq) > radsq ) cycle
                        mm(1) = min(mm(1),self%rmat(i,j,k))
                        mm(2) = max(mm(2),self%rmat(i,j,k))
                    enddo
                enddo
            enddo
        else
            mm(1) = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
            mm(2) = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        endif
    end function minmax

    subroutine loc_sdev( self, winsz, sdevimg, asdev )
        class(image),   intent(in)    :: self
        integer,        intent(in)    :: winsz
        class(image),   intent(inout) :: sdevimg
        real, optional, intent(inout) :: asdev
        real                 :: avg
        integer              :: i, j, k, ir(2), jr(2), kr(2), isz, jsz, ksz, npix
        logical              :: isvol
        isvol = .false.
        isvol = self%is_3d()
        if( isvol )then
            call sdevimg%new(self%ldim, self%smpd)
            !$omp parallel do private(i,j,k,ir,jr,kr,isz,jsz,ksz,npix,avg) default(shared) proc_bind(close)
            do i = 1, self%ldim(1)
                ir(1) = max(1,            i - winsz)
                ir(2) = min(self%ldim(1), i + winsz)
                isz   = ir(2) - ir(1) + 1
                do j = 1, self%ldim(2)
                    jr(1) = max(1,            j - winsz)
                    jr(2) = min(self%ldim(2), j + winsz)
                    jsz   = jr(2) - jr(1) + 1
                    do k = 1, self%ldim(3)
                        kr(1)               = max(1,            k - winsz)
                        kr(2)               = min(self%ldim(3), k + winsz)
                        ksz                 = kr(2) - kr(1) + 1
                        npix                = isz * jsz * ksz
                        avg                 = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))) / real(npix)
                        sdevimg%rmat(i,j,1) = sqrt(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2)) - avg)**2.0) / real(npix - 1)) 
                    enddo    
                enddo
            enddo
            !$omp end parallel do
            if( present(asdev) )then
                asdev = sum(sdevimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(self%ldim(1) * self%ldim(2) * self%ldim(3))
            endif
        else
            call sdevimg%new(self%ldim, self%smpd)
            !$omp parallel do private(i,j,ir,jr,isz,jsz,npix,avg) default(shared) proc_bind(close)
            do i = 1,self%ldim(1)
               ir(1) = max(1,            i - winsz)
               ir(2) = min(self%ldim(1), i + winsz)
               isz   = ir(2) - ir(1) + 1
               do j = 1,self%ldim(2)
                   jr(1)               = max(1,            j - winsz)
                   jr(2)               = min(self%ldim(2), j + winsz)
                   jsz                 = jr(2) - jr(1) + 1
                   npix                = isz * jsz
                   avg                 = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),1)) / real(npix)
                   sdevimg%rmat(i,j,1) = sqrt(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),1) - avg)**2.0) / real(npix - 1)) 
                enddo
            enddo
            !$omp end parallel do
            if( present(asdev) )then
                asdev = sum(sdevimg%rmat(:self%ldim(1),:self%ldim(2),1)) / real(self%ldim(1) * self%ldim(2))
            endif
        endif
    end subroutine loc_sdev

    function avg_loc_sdev( self, winsz ) result( asdev )
        class(image), intent(in) :: self
        integer,      intent(in) :: winsz
        real(dp)             :: sum_sdevs
        real                 :: avg, asdev
        integer              :: i, j, k, ir(2), jr(2), kr(2), isz, jsz, ksz, npix
        logical              :: isvol
        isvol     = .false.
        isvol     = self%is_3d()
        sum_sdevs = 0.d0
        if( isvol )then
            do i = 1, self%ldim(1)
                ir(1) = max(1,            i - winsz)
                ir(2) = min(self%ldim(1), i + winsz)
                isz   = ir(2) - ir(1) + 1
                do j = 1, self%ldim(2)
                    jr(1)     = max(1,            j - winsz)
                    jr(2)     = min(self%ldim(2), j + winsz)
                    jsz       = jr(2) - jr(1) + 1
                    do k = 1, self%ldim(3)
                        kr(1) = max(1,            k - winsz)
                        kr(2) = min(self%ldim(3), k + winsz)
                        ksz       = kr(2) - kr(1) + 1
                        npix      = isz * jsz * ksz
                        avg       = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))) / real(npix)
                        sum_sdevs = sum_sdevs + sqrt(real(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2)) - avg)**2.0),dp) / real(npix-1,dp))
                    enddo
                enddo
            enddo
            asdev = real(sum_sdevs / real(self%ldim(1) * self%ldim(2) * self%ldim(3),dp))
        else
            do i = 1, self%ldim(1)
                ir(1) = max(1,            i - winsz)
                ir(2) = min(self%ldim(1), i + winsz)
                isz   = ir(2) - ir(1) + 1
                do j = 1, self%ldim(2)
                    jr(1)     = max(1,            j - winsz)
                    jr(2)     = min(self%ldim(2), j + winsz)
                    jsz       = jr(2) - jr(1) + 1
                    npix      = isz * jsz
                    avg       = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),1)) / real(npix)
                    sum_sdevs = sum_sdevs + sqrt(real(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),1) - avg)**2.0),dp) / real(npix-1,dp))
                enddo
            enddo
            asdev = real(sum_sdevs / real(self%ldim(1) * self%ldim(2),dp))
        endif
    end function avg_loc_sdev

    subroutine loc_var( self, varimg, avar )
        class(image),   intent(in)    :: self
        class(image),   intent(inout) :: varimg
        real, optional, intent(inout) :: avar
        real    :: avg, ep, val, var
        integer :: i, j, l, nsz, n_4(3,8)
        if( self%ldim(3) /= 1 ) THROW_HARD('not for 3d')
        call varimg%new(self%ldim, self%smpd)
        do i = 1,self%ldim(1)
            do j = 1,self%ldim(2)
                call neigh_8(self%ldim, [i,j,1], n_4, nsz)
                avg = 0.
                do l = 1,nsz
                    avg = avg + self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                end do
                avg = avg / real(nsz)
                ep  = 0.
                var = 0.
                do l = 1,nsz
                    val = self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                    ep  = ep  + val
                    var = var + val * val
                end do
                var = (var-ep**2./real(nsz))/(real(nsz)-1.) ! corrected two-pass formula
                varimg%rmat(i,j,1) = var
            end do
        end do
        if( present(avar) )then
            avar = sum(varimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        endif
    end subroutine loc_var

    subroutine loc_var3D( self, varimg, avar )
        class(image),   intent(in)    :: self
        class(image),   intent(inout) :: varimg
        real, optional, intent(inout) :: avar
        real    :: avg, ep, val, var
        integer :: i, j, k, l, nsz, n_4(3,6)
        if( self%ldim(3) == 1 ) THROW_HARD('not for 2d')
        call varimg%new(self%ldim, self%smpd)
        !$omp parallel do private(i,j,k,l,n_4,nsz,avg,ep,val,var) default(shared) proc_bind(close)
        do i = 1,self%ldim(1)
            do j = 1,self%ldim(2)
                do k = 1,self%ldim(3)
                    call neigh_4_3D(self%ldim, [i,j,k], n_4, nsz)
                    avg = 0.
                    do l = 1,nsz
                        avg = avg + self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                    end do
                    avg = avg / real(nsz)
                    ep  = 0.
                    var = 0.
                    do l = 1,nsz
                        val = self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                        ep  = ep  + val
                        var = var + val * val
                    end do
                    var = (var-ep**2./real(nsz))/(real(nsz)-1.) ! corrected two-pass formula
                    varimg%rmat(i,j,k) = var 
                end do
            end do
        end do
        !$omp end parallel do
        if( present(avar) )then
            avar = sum(varimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        endif
    end subroutine loc_var3D

    !>  \brief rmsd for calculating the RMSD of a map
    !! \return  dev root mean squared deviation
    subroutine rmsd( self, dev, mean )
        class(image),   intent(inout) :: self
        real,           intent(out)   :: dev
        real, optional, intent(out)   :: mean
        real :: avg
        if( self%ft )then
            dev = 0.
            if( present(mean) ) mean = 0.
        else
            avg = self%mean()
            if(present(mean)) mean = avg
            dev = sum((self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - avg)**2.0)&
                  &/ real(product(self%ldim))
            if( dev > 0. )then
                dev = sqrt(dev)
            else
                dev = 0.
            endif
        endif
    end subroutine rmsd

    !> \brief noisesdev is for estimating the noise variance of an image
    !!          by online estimation of the variance of the background pixels
    function noisesdev( self, msk ) result( sdev )
        use simple_online_var, only: online_var
        class(image), intent(inout) :: self
        real, intent(in)            :: msk
        type(online_var)            :: ovar
        integer                     :: i, j, k
        real                        :: ci, cj, ck, e, sdev, mv(2)
        logical                     :: didft
        ovar = online_var( )
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ci = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        e = hardedge(ci,cj,ck,msk)
                    else
                        e = hardedge(ci,cj,msk)
                    endif
                    if( e < 0.5 )then
                        call ovar%add(self%rmat(i,j,k))
                    endif
                    ck = ck+1
                end do
                cj = cj+1.
            end do
            ci = ci+1.
        end do
        mv(1) = ovar%get_mean()
        mv(2) = ovar%get_var()
        sdev = 0.
        if( mv(2) > 0. ) sdev = sqrt(mv(2))
        if( didft ) call self%fft()
    end function noisesdev

    !>  \brief  is for calculating the mean of an image
    function mean( self ) result( avg )
        class(image), intent(inout) :: self
        real    :: avg
        logical :: didft
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        avg = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))/real(product(self%ldim))
        if( didft ) call self%ifft()
    end function mean

    !>  \brief  is for checking the numerical soundness of an image
    logical function contains_nans( self )
        class(image), intent(in) :: self
        integer :: i, j, k
        contains_nans = .false.
        do i=1,size(self%rmat,1)
            do j=1,size(self%rmat,2)
                do k=1,size(self%rmat,3)
                    if( .not. is_a_number(self%rmat(i,j,k)) )then
                        contains_nans = .true.
                        return
                    endif
                end do
            end do
        end do
    end function contains_nans

    !> \brief checkimg4nans  is for checking the numerical soundness of an image
    subroutine checkimg4nans( self )
        class(image), intent(in) :: self
        if( self%ft )then
            call check4nans3D(self%cmat)
        else
            call check4nans3D(self%rmat)
        endif
    end subroutine checkimg4nans

    !> \brief cure_2  is for checking the numerical soundness of an image and curing it if necessary
    subroutine cure( self, maxv, minv, ave, sdev, n_nans )
        class(image), intent(inout) :: self
        real,         intent(out)   :: maxv, minv, ave, sdev
        integer,      intent(out)   :: n_nans
        integer                     :: i, j, k, npix
        real                        :: var, ep, dev
        if( self%ft )then
            THROW_WARN('cannot cure FTs; cure')
            return
        endif
        npix   = product(self%ldim)
        n_nans = 0
        ave    = 0.
        !$omp parallel do default(shared) private(i,j,k) schedule(static)&
        !$omp collapse(3) proc_bind(close) reduction(+:n_nans,ave)
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if( .not. is_a_number(self%rmat(i,j,k)) )then
                        n_nans = n_nans + 1
                    else
                        ave = ave + self%rmat(i,j,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( n_nans > 0 )then
            write(logfhandle,*) 'found NaNs in simple_image; cure:', n_nans
        endif
        ave       = ave/real(npix)
        maxv      = maxval( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) )
        minv      = minval( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) )
        self%rmat = self%rmat - ave
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,j,k,dev) schedule(static)&
        !$omp collapse(3) proc_bind(close) reduction(+:ep,var)
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    dev = self%rmat(i,j,k)
                    ep  = ep + dev
                    var = var + dev * dev
                end do
            end do
        end do
        !$omp end parallel do
        var  = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        sdev = sqrt(var)
        if( sdev > 0. ) self%rmat = self%rmat/sdev
    end subroutine cure

    !>  \brief loop_lims is for determining loop limits for transforms
    function loop_lims( self, mode, lp_dyn ) result( lims )
        class(image), intent(in)   :: self
        integer, intent(in)        :: mode
        real, intent(in), optional :: lp_dyn
        integer                    :: lims(3,2)
        if( present(lp_dyn) )then
            lims = self%fit%loop_lims(mode, lp_dyn)
        else
            lims = self%fit%loop_lims(mode)
        endif
    end function loop_lims

    ! This function returns a the gradient matrix of the input image.
    ! It is also possible to have derivates row and column
    ! as output (optional).
    ! It uses masks found in http://www.holoborodko.com/pavel/image-processing/edge-detection/
    ! which is better than Sobel masks because of:
    !                     1) isotropic noise suppression
    !                     2) the estimation of the gradient is still precise
    subroutine calc_gradient(self, grad, Dc, Dr)
        class(image),   intent(inout) :: self
        real,           intent(out)   :: grad(self%ldim(1), self%ldim(2), self%ldim(3)) ! gradient matrix
        real, optional, intent(out)   :: Dc(self%ldim(1), self%ldim(2), self%ldim(3)), Dr(self%ldim(1), self%ldim(2), self%ldim(3)) ! derivates column and row matrices
        type(image)        :: img_p                         ! padded image
        real, allocatable  :: wc(:,:), wr(:,:)              ! row and column Sobel masks
        integer, parameter :: L1 = 5 , L2 = 3               ! dimension of the masks
        integer            :: ldim(3)                       ! dimension of the image, save just for comfort
        integer            :: i,j,m,n                       ! loop indeces
        real :: Ddc(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: Ddr(self%ldim(1),self%ldim(2),self%ldim(3)) ! column and row derivates
        ldim = self%ldim
        allocate(wc((-(L1-1)/2):((L1-1)/2),(-(L2-1)/2):((L2-1)/2)),&
        &        wr(-(L2-1)/2:(L2-1)/2,-(L1-1)/2:(L1-1)/2), source = 0.)
        wc = (1./32.)*reshape([-1,-2,0,2,1,-2,-4,0,4,2,-1,-2,0,2,1], [L1,L2])
        wr = (1./32.)*reshape([-1,-2,-1,-2,-4,-2,0,0,0,2,4,2,1,2,1], [L2,L1])
        Ddc  = 0. ! initialisation
        Ddr  = 0.
        grad = 0.
        call img_p%new([ldim(1)+L1-1,ldim(2)+L1-1,1],1.) ! pad with the biggest among L1 and L2
        call self%pad(img_p)
        do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L1-1)/2,(L1-1)/2
                  do n = -(L2-1)/2,(L2-1)/2
                      Ddc(i,j,1) = Ddc(i,j,1)+img_p%rmat(i+m+2,j+n+2,1)*wc(m,n)
                  end do
              end do
          end do
        end do
        do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L2-1)/2,(L2-1)/2
                  do n = -(L1-1)/2,(L1-1)/2
                      Ddr(i,j,1) = Ddr(i,j,1)+img_p%rmat(i+m+2,j+n+2,1)*wr(m,n)
                  end do
              end do
          end do
        end do
        deallocate(wc, wr)
        grad = sqrt(Ddc**2 + Ddr**2)
        if(present(Dc)) Dc = Ddc
        if(present(Dr)) Dr = Ddr
        call img_p%kill
    end subroutine calc_gradient

    subroutine gradients_magnitude( self, self_out )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: self_out
        integer :: i,j,ni,nj
        if( self%is_ft() ) THROW_HARD('Image input must be in the spatial domain!')
        if( .not.self%is_2d() ) THROW_HARD('Image input must be in 2D!')
        if( .not.self_out%exists() ) call self_out%copy(self)
        ni = self%ldim(1)
        nj = self%ldim(2)
        !$omp parallel private(i,j) proc_bind(close) default(shared)
        !$omp do
        do i = 1,ni
            if( i == 1 )then
                self_out%rmat(i,:nj,1) = (self%rmat(2,:nj,1) - self%rmat(1,:nj,1))**2
            else if( i == ni )then
                self_out%rmat(i,:nj,1) = (self%rmat(i,:nj,1) - self%rmat(i-1,:nj,1))**2
            else
                self_out%rmat(i,:nj,1) = (self%rmat(i+1,:nj,1) - self%rmat(i-1,:nj,1))**2
            endif
        enddo
        !$omp end do
        !$omp do
        do j = 1,nj
            if( j == 1 )then
                self_out%rmat(:ni,j,1) = self_out%rmat(:ni,j,1) + (self%rmat(:ni,2,1) - self%rmat(:ni,1,1))**2
            else if( j == nj )then
                self_out%rmat(:ni,j,1) = self_out%rmat(:ni,j,1) + (self%rmat(:ni,j,1) - self%rmat(:ni,j-1,1))**2
            else
                self_out%rmat(:ni,j,1) = self_out%rmat(:ni,j,1) + (self%rmat(:ni,j+1,1) - self%rmat(:ni,j-1,1))**2
            endif
        enddo
        !$omp end do
        !$omp workshare
        self_out%rmat(1:ni,1:nj,1) = self_out%rmat(1:ni,1:nj,1)/4.
        self_out%rmat(1:ni,1:nj,1) = merge(sqrt(self_out%rmat(1:ni,1:nj,1)), 0., self_out%rmat(1:ni,1:nj,1)>0.)
        !$omp end workshare
        !$omp end parallel
    end subroutine gradients_magnitude

    subroutine calc_ice_score( self, score )
        class(image), intent(in)  :: self
        real,         intent(out) :: score
        real,   parameter :: START_FREQ = 15.
        real,   parameter :: END_FREQ   = 6.
        real, allocatable :: res(:), tmp(:)
        real    :: powspec(fdim(self%ldim(1)) - 1)
        real    :: g, gs, ge, mag, mag_max, band_avg, ice_avg
        integer :: lims(3,2), ice_maxind, start_find, end_find
        integer :: nbands, s, e, h, k, hmax, kmax, sh, cnt
        score = 0.
        if( self%smpd > (ICE_BAND1/2.) ) return
        if( .not.self%is_ft() ) THROW_HARD('Image input must be in the Fourier domain!; calc_ice_score')
        lims = self%loop_lims(2)
        res  = get_resarr(self%ldim(1), self%smpd)
        ice_maxind = get_find_at_res(res, ICE_BAND1)
        start_find = get_find_at_res(res, START_FREQ)
        end_find   = get_find_at_res(res, END_FREQ)
        call self%power_spectrum(powspec)
        nbands = end_find-start_find+1
        tmp = powspec(start_find:end_find)
        call hpsort(tmp)
        e = nbands
        s = nint(0.5 *real(nbands))
        band_avg = sum(tmp(s:e)) / real(e-s+1)
        ! location of maximum in ice band
        mag_max = -1.
        gs = real(max(            1,   ice_maxind-3)) / real(self%ldim(1))
        ge = real(min(size(powspec)-1, ice_maxind+3)) / real(self%ldim(1))
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                sh = nint(hyp(h,k))
                g  = real(sh) / real(self%ldim(1))
                if( g > gs .and. g < ge )then
                    mag = csq_fast(self%get_fcomp2D(h,k))
                    if( mag > mag_max )then
                        hmax = h
                        kmax = k
                        mag_max = mag
                    endif
                endif
            end do
        end do
        ! ice peak
        ice_avg = 0.
        cnt     = 0
        do k = kmax-1,kmax+1
            do h = hmax-1,hmax+1
                sh = nint(hyp(h,k))
                g  = real(sh) / real(self%ldim(1))
                if( g < 0.5 )then
                    ice_avg = ice_avg + csq_fast(self%get_fcomp2D(h,k))
                    cnt     = cnt+1
                endif
            enddo
        enddo
        ice_avg = ice_avg / real(cnt)
        score   = ice_avg / (band_avg + TINY)
    end subroutine calc_ice_score

    ! This function returns the derivates row, column, and z as outputs,
    ! together with the optional gradient matrix/volume.
    ! It uses standard central difference scheme
    subroutine gradient(self, Dc, Dr, Dz, grad)
        class(image),   intent(inout) :: self
        real, optional, intent(out)   :: Dc(self%ldim(1), self%ldim(2), self%ldim(3)), & ! derivates column matrix
                                         Dr(self%ldim(1), self%ldim(2), self%ldim(3)), & ! derivates row matrix
                                         Dz(self%ldim(1), self%ldim(2), self%ldim(3))    ! derivates z matrix
        real, optional, intent(out)   :: grad(self%ldim(1), self%ldim(2), self%ldim(3))  ! gradient matrix
        integer :: ldim(3), k
        real    :: Ddc(self%ldim(1),self%ldim(2),self%ldim(3))
        real    :: Ddr(self%ldim(1),self%ldim(2),self%ldim(3))
        real    :: Ddz(self%ldim(1),self%ldim(2),self%ldim(3))
        ldim = self%ldim
        Ddc  = 0. ! initialisation
        Ddr  = 0.
        Ddz  = 0.
        do k = 2, ldim(1)-1
            Ddc(k,1:ldim(2),1:ldim(3)) = 0.5*(self%rmat(k+1,1:ldim(2),1:ldim(3)) - self%rmat(k-1,1:ldim(2),1:ldim(3)))
        enddo
        do k = 2, ldim(2)-1
            Ddr(1:ldim(1),k,1:ldim(3)) = 0.5*(self%rmat(1:ldim(1),k+1,1:ldim(3)) - self%rmat(1:ldim(1),k-1,1:ldim(3)))
        enddo
        Ddc(1        ,1:ldim(2),1:ldim(3)) = self%rmat(2        ,1:ldim(2),1:ldim(3)) - self%rmat(1          ,1:ldim(2)  ,1:ldim(3))
        Ddc(  ldim(1),1:ldim(2),1:ldim(3)) = self%rmat(  ldim(1),1:ldim(2),1:ldim(3)) - self%rmat(  ldim(1)-1,1:ldim(2)  ,1:ldim(3))
        Ddr(1:ldim(1),1        ,1:ldim(3)) = self%rmat(1:ldim(1),2        ,1:ldim(3)) - self%rmat(1:ldim(1)  ,1          ,1:ldim(3))
        Ddr(1:ldim(1),  ldim(2),1:ldim(3)) = self%rmat(1:ldim(1),  ldim(2),1:ldim(3)) - self%rmat(1:ldim(1)  ,  ldim(2)-1,1:ldim(3))
        if( ldim(3) > 1 )then
            do k = 2, ldim(3)-1
                Ddz(1:ldim(1),1:ldim(2),k) = 0.5*(self%rmat(1:ldim(1),1:ldim(2),k+1) - self%rmat(1:ldim(1),1:ldim(2),k-1))
            enddo
            Ddz(1:ldim(1),1:ldim(2),1      ) = self%rmat(1:ldim(1),1:ldim(2),2)       - self%rmat(1:ldim(1),1:ldim(2),1)
            Ddz(1:ldim(1),1:ldim(2),ldim(3)) = self%rmat(1:ldim(1),1:ldim(2),ldim(3)) - self%rmat(1:ldim(1),1:ldim(2),ldim(3)-1)
        endif
        if(present(Dc))   Dc   = Ddc
        if(present(Dr))   Dr   = Ddr
        if(present(Dz))   Dz   = Ddz
        if(present(grad)) grad = sqrt(Ddc**2 + Ddr**2 + Ddz**2)
    end subroutine gradient

    !>  \brief  Convert logical address to physical address. Complex image.
    pure function comp_addr_phys1(self,logi) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3) !<  Logical address
        integer                  :: phys(3) !<  Physical address
        phys = self%fit%comp_addr_phys(logi)
    end function comp_addr_phys1

    !>  \brief  Convert logical address to physical address. Complex image.
    pure function comp_addr_phys2(self,h,k,m) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k,m !<  Logical address
        integer                  :: phys(3) !<  Physical address
        phys = self%fit%comp_addr_phys(h,k,m)
    end function comp_addr_phys2

    !>  \brief  Convert 2D logical address to physical address. Complex image.
    pure function comp_addr_phys3(self,h,k) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k     !<  Logical address
        integer                  :: phys(2) !<  Physical address
        phys = self%fit%comp_addr_phys(h,k)
    end function comp_addr_phys3

    !>  \brief corr is for correlating two images
    function corr( self1, self2, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self1, self2
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real    :: r, sumasq, sumbsq, eps
        integer :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        logical :: didft1, didft2
        r = 0.
        sumasq = 0.
        sumbsq = 0.
        eps    = epsilon(sumasq)
        if( self1.eqdims.self2 )then
            didft1 = .false.
            if( .not. self1%ft )then
                call self1%fft()
                didft1 = .true.
            endif
            didft2 = .false.
            if( .not. self2%ft )then
                call self2%fft()
                didft2 = .true.
            endif
            if( present(lp_dyn) )then
                lims = self1%fit%loop_lims(1,lp_dyn)
            else
                lims = self1%fit%loop_lims(2) ! Nyqvist default low-pass limit
            endif
            sqlp = (maxval(lims(:,2)))**2
            if( present(hp_dyn) )then
                sqhp = max(2,self1%get_find(hp_dyn))**2
            else
                sqhp = 2 ! index 2 default high-pass limit
            endif
            !$omp parallel do collapse(3) default(shared) private(h,k,l,sqarg,phys)&
            !$omp reduction(+:r,sumasq,sumbsq) schedule(static) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sqarg = h*h + k*k + l*l
                        if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                            phys = self1%fit%comp_addr_phys([h,k,l])
                            ! real part of the complex mult btw 1 and 2*
                            r = r + real(self1%cmat(phys(1),phys(2),phys(3))*conjg(self2%cmat(phys(1),phys(2),phys(3))))
                            sumasq = sumasq + csq(self2%cmat(phys(1),phys(2),phys(3)))
                            sumbsq = sumbsq + csq(self1%cmat(phys(1),phys(2),phys(3)))
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
            if( r < eps .and. sumasq < eps .and. sumbsq < eps )then
                r = 1.
            elseif( sqrt(sumasq * sumbsq) < eps )then
                r = 0.
            else
                r = r / sqrt(sumasq * sumbsq)
            endif
            if( didft1 ) call self1%ifft()
            if( didft2 ) call self2%ifft()
        else
            write(logfhandle,*) 'self1%ldim:', self1%ldim
            write(logfhandle,*) 'self2%ldim:', self2%ldim
            THROW_HARD('images to be correlated need to have same dimensions; corr')
        endif
    end function corr

    function corr_shifted( self_ref, self_ptcl, shvec, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self_ref, self_ptcl
        real,           intent(in)    :: shvec(3)
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real                          :: r, sumasq, sumbsq
        complex                       :: shcomp
        integer                       :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        ! this is for highly optimised code, so we assume that images are always Fourier transformed beforehand
        if( .not. self_ref%ft  ) THROW_HARD('self_ref not FTed;  corr_shifted')
        if( .not. self_ptcl%ft ) THROW_HARD('self_ptcl not FTed; corr_shifted')
        r = 0.
        sumasq = 0.
        sumbsq = 0.
        if( present(lp_dyn) )then
            lims = self_ref%fit%loop_lims(1,lp_dyn)
        else
            lims = self_ref%fit%loop_lims(2) ! Nyqvist default low-pass limit
        endif
        sqlp = (maxval(lims(:,2)))**2
        if( present(hp_dyn) )then
            sqhp = max(2,self_ref%get_find(hp_dyn))**2
        else
            sqhp = 2 ! index 2 default high-pass limit
        endif
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sqarg,phys,shcomp)&
        !$omp reduction(+:r,sumasq,sumbsq) schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sqarg = h*h + k*k + l*l
                    if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                        phys = self_ref%fit%comp_addr_phys(h,k,l)
                        ! shift particle
                        shcomp = self_ptcl%cmat(phys(1),phys(2),phys(3))*&
                            &self_ptcl%oshift([h,k,l], shvec)
                        ! real part of the complex mult btw 1 and 2*
                        r = r + real(self_ref%cmat(phys(1),phys(2),phys(3))*conjg(shcomp))
                        sumasq = sumasq + csq(shcomp)
                        sumbsq = sumbsq + csq(self_ref%cmat(phys(1),phys(2),phys(3)))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( sumasq > 0. .and. sumbsq > 0. )then
            r = r / sqrt(sumasq * sumbsq)
        else
            r = 0.
        endif
    end function corr_shifted

    !>  \brief is for calculating a real-space correlation coefficient between images
    !! \param self1,self2 image objects
    !! \return  r correlation coefficient
    function real_corr_1( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real    :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3)), diff1sc
        real    :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3)), diff2sc
        real    :: r, ax, ay, sxx, syy, sxy, npix
        integer :: i,j,k
        if( self1%wthreads .and. self2%wthreads )then
            npix = real(product(self1%ldim))
            ax   = self1%mean()
            ay   = self2%mean()
            sxx = 0.
            syy = 0.
            sxy = 0.
            !$omp parallel do default(shared) private(i,j,k,diff1sc,diff2sc) collapse(3) proc_bind(close) schedule(static) reduction(+:sxx,syy,sxy)
            do i=1,self1%ldim(1)
                do j=1,self1%ldim(2)
                    do k=1,self1%ldim(3)
                        diff1sc = self1%rmat(i,j,k) -ax
                        diff2sc = self2%rmat(i,j,k) -ay
                        sxx     = sxx + diff1sc * diff1sc
                        syy     = syy + diff2sc * diff2sc
                        sxy     = sxy + diff1sc * diff2sc
                    end do
                end do
            end do
            !$omp end parallel do
            if( sxx > TINY .and. syy > TINY )then
                r = sxy / sqrt(sxx * syy)
            else
                r = 0.
            endif
        else
            diff1 = 0.
            diff2 = 0.
            npix  = real(product(self1%ldim))
            ax    = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3))) / npix
            ay    = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3))) / npix
            diff1 = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
            diff2 = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
            sxx   = sum(diff1 * diff1)
            syy   = sum(diff2 * diff2)
            sxy   = sum(diff1 * diff2)
            if( sxx > TINY .and. syy > TINY )then
                r = sxy / sqrt(sxx * syy)
            else
                r = 0.
            endif
        endif
    end function real_corr_1

    !>  \brief real_corr_2 is for calculating a real-space correlation coefficient between images within a mask
    function real_corr_2( self1, self2, mask ) result( r )
        class(image), intent(inout) :: self1, self2
        logical,      intent(in)    :: mask(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3))
        real :: r, sxx, syy, sxy, npix, ax, ay
        diff1 = 0.
        diff2 = 0.
        npix  = real(count(mask))
        ax    = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)), mask=mask) / npix
        ay    = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)), mask=mask) / npix
        diff1 = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
        diff2 = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
        sxx   = sum(diff1 * diff1, mask=mask)
        syy   = sum(diff2 * diff2, mask=mask)
        sxy   = sum(diff1 * diff2, mask=mask)
        if( sxx > TINY .and. syy > TINY )then
            r = sxy / sqrt(sxx * syy)
        else
            r = 0.
        endif
    end function real_corr_2

    function euclid_dist_two_imgs(self1, self2, mask1) result(dist)
        use simple_linalg, only: euclid
        class(image),      intent(inout) :: self1, self2
        logical, optional, intent(in)    :: mask1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real              :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real              :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3))
        real              :: ax, ay, npix, dist
        real, allocatable :: diff1_flat(:), diff2_flat(:)
        logical           :: mask_here(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        if (present(mask1)) then
            mask_here = mask1
        else 
            mask_here = .true.
        end if
        npix       = real(count(mask_here))
        ax         = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)), mask=mask_here) / npix
        ay         = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)), mask=mask_here) / npix
        diff1      = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
        diff2      = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
        diff1_flat = pack(diff1, mask=.true.)
        diff2_flat = pack(diff2, mask=.true.)
        dist       = euclid(diff1_flat, diff2_flat)
    end function euclid_dist_two_imgs

    !> /brief phase_corr calculates the phase correlation between
    !> self1 and self2 and returns it in pc.
    !> if border is present, then it calculates it discarding
    !> the external frame in self1 self2 with width border.
    ! FORMULA: phasecorr = ifft2(fft2(self1).*conj(fft2(self2)));
    ! Preserves Fourier state
    ! lp: frequencies below lp are set to 0 in a
    ! smooth way (cos edge).
    subroutine phase_corr(self1, self2, pc, lp )
        class(image),      intent(inout) :: self1, self2, pc
        real,              intent(in)    :: lp
        real, parameter :: width = 3.
        complex     :: c1,c2
        real        :: w,rw,rlplim,rsh,normsq
        real(dp)    :: sqsum1,sqsum2
        integer     :: nrflims(3,2),phys(3)
        integer     :: h,k,l,shsq, lplim, lplimsq, bplplimsq
        if( .not. all([self1%is_ft(),self2%is_ft(),pc%is_ft()]) )then
            THROW_HARD('All inputted images must be FTed')
        endif
        sqsum1    = 0.d0
        sqsum2    = 0.d0
        nrflims   = self1%loop_lims(2)
        lplim     = calc_fourier_index(lp,minval(self1%ldim(1:2)),self1%smpd)
        rlplim    = real(lplim)
        lplimsq   = lplim*lplim
        bplplimsq = min(minval(nrflims(1:2,2)),lplim-nint(WIDTH))**2
        call pc%zero_and_flag_ft
        if( self1%wthreads .or. self2%wthreads )then
            !$omp parallel do default(shared) private(h,k,l,w,rw,shsq,phys,rsh,c1,c2)&
            !$omp proc_bind(close) schedule(static) reduction(+:sqsum1,sqsum2)
            do h = nrflims(1,1),nrflims(1,2)
                rw = merge(1., 2., h==0)
                do k = nrflims(2,1),nrflims(2,2)
                    do l = nrflims(3,1),nrflims(3,2)
                        shsq = h*h+k*k+l*l
                        w    = 1.
                        if( shsq > lplimsq )then
                            cycle
                        else if( shsq == 0 )then
                            cycle
                        else if( shsq > bplplimsq )then
                            rsh = sqrt(real(shsq))
                            w   = 0.5*(1.+cos(PI*(rsh-(rlplim-width))/width))
                        endif
                        phys = pc%comp_addr_phys(h,k,l)
                        c1 = w*self1%cmat(phys(1),phys(2),phys(3))
                        c2 = w*self2%cmat(phys(1),phys(2),phys(3))
                        sqsum1 = sqsum1 + real(rw*csq(c1),dp)
                        sqsum2 = sqsum2 + real(rw*csq(c2),dp)
                        pc%cmat(phys(1),phys(2),phys(3)) = c1 * conjg(c2)
                    enddo
                enddo
            enddo
            !$omp end parallel do
        else
            do h = nrflims(1,1),nrflims(1,2)
                rw = merge(1., 2., h==0)
                do k = nrflims(2,1),nrflims(2,2)
                    do l = nrflims(3,1),nrflims(3,2)
                        shsq = h*h+k*k+l*l
                        w    = 1.
                        if( shsq > lplimsq )then
                            cycle
                        else if( shsq == 0 )then
                            cycle
                        else if( shsq > bplplimsq )then
                            rsh = sqrt(real(shsq))
                            w   = 0.5*(1.+cos(PI*(rsh-(rlplim-width))/width))
                        endif
                        phys = pc%comp_addr_phys(h,k,l)
                        c1 = w*self1%cmat(phys(1),phys(2),phys(3))
                        c2 = w*self2%cmat(phys(1),phys(2),phys(3))
                        sqsum1 = sqsum1 + real(rw*csq(c1),dp)
                        sqsum2 = sqsum2 + real(rw*csq(c2),dp)
                        pc%cmat(phys(1),phys(2),phys(3)) = c1 * conjg(c2)
                    enddo
                enddo
            enddo
        endif
        call pc%ifft()
        normsq = real(sqsum1*sqsum2)
        if( is_a_number(normsq) )then
            if( normsq > 1.0e-12 ) pc%rmat = pc%rmat / sqrt(normsq)
        endif
    end subroutine phase_corr

    ! returns the discrete shift that registers self2 to self1. self1 is the unnormalized correlation image on output
    subroutine fcorr_shift( self1, self2, trs, shift, peak_interp )
        class(image),      intent(inout) :: self1, self2
        real,              intent(in)    :: trs
        real,              intent(inout) :: shift(2)
        logical, optional, intent(in)    :: peak_interp
        real    :: alpha, beta, gamma, denom
        integer :: center(2), pos(2), itrs
        logical :: l_interp
        if( self1%is_3d() .or. self2%is_3d() ) THROW_HARD('2d only supported')
        if( .not.(self1%is_ft() .and. self2%is_ft()) ) THROW_HARD('FTed only supported')
        if( .not.(self1.eqdims.self2) ) THROW_HARD('Inconsistent dimensions in fcorr_shift')
        l_interp = .false.
        if(present(peak_interp)) l_interp = peak_interp
        ! dimensions
        center = self1%ldim(1:2)/2+1
        itrs   = min(floor(trs),minval(center)-1)
        ! Correlation image
        self1%cmat = self1%cmat * conjg(self2%cmat)
        call self1%ifft
        ! maximum correlation & offset
        pos   = maxloc(self1%rmat(center(1)-itrs:center(1)+itrs, center(2)-itrs:center(2)+itrs, 1)) -itrs-1
        ! peak interpolation
        if( l_interp )then
            shift = real(pos)
            beta  = self1%rmat(center(1)+pos(1),center(2)+pos(2), 1)
            ! along x
            if( abs(pos(1)) < itrs )then ! within limits
                alpha = self1%rmat(center(1)+pos(1)-1,center(2)+pos(2), 1)
                gamma = self1%rmat(center(1)+pos(1)+1,center(2)+pos(2), 1)
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(1) = shift(1) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! along y
            if( abs(pos(2)) < itrs )then
                alpha = self1%rmat(center(1)+pos(1),center(2)+pos(2)-1, 1)
                gamma = self1%rmat(center(1)+pos(1),center(2)+pos(2)+1, 1)
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(2) = shift(2) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! convention
            shift = -shift
        else
            shift = -real(pos)
        endif
    end subroutine fcorr_shift

    ! returns the discrete shift that registers self2 to self1. self1 is the unnormalized correlation image on output
    subroutine fcorr_shift3D( self1, self2, trs, shift, peak_interp )
        class(image),      intent(inout) :: self1, self2
        real,              intent(in)    :: trs
        real,              intent(inout) :: shift(3)
        logical, optional, intent(in)    :: peak_interp
        real    :: alpha, beta, gamma, denom
        integer :: cen(3), pos(3), itrs
        logical :: l_interp
        if( self1%is_2d() .or. self2%is_2d() ) THROW_HARD('3d only supported')
        if( .not.(self1%is_ft() .and. self2%is_ft()) ) THROW_HARD('FTed only supported')
        if( .not.(self1.eqdims.self2) ) THROW_HARD('Inconsistent dimensions in fcorr_shift')
        l_interp = .false.
        if(present(peak_interp)) l_interp = peak_interp
        ! dimensions
        cen  = self1%ldim/2+1
        itrs = min(floor(trs),minval(cen)-1)
        ! Correlation image
        self1%cmat = self1%cmat * conjg(self2%cmat)
        call self1%ifft
        ! maximum correlation & offset
        pos = maxloc(self1%rmat(cen(1)-itrs:cen(1)+itrs, cen(2)-itrs:cen(2)+itrs, cen(3)-itrs:cen(3)+itrs)) -itrs-1
        ! peak interpolation
        if( l_interp )then
            shift = real(pos)
            ! pos   = pos+itrs+1
            beta  = self1%rmat(cen(1)+pos(1),cen(2)+pos(2), cen(3)+pos(3))
            ! along x
            if( abs(pos(1)) < itrs )then ! within limits
                alpha = self1%rmat(cen(1)+pos(1)-1, cen(2)+pos(2), cen(3)+pos(3))
                gamma = self1%rmat(cen(1)+pos(1)+1, cen(2)+pos(2), cen(3)+pos(3))
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(1) = shift(1) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! along y
            if( abs(pos(2)) < itrs )then
                alpha = self1%rmat(cen(1)+pos(1), cen(2)+pos(2)-1, cen(3)+pos(3))
                gamma = self1%rmat(cen(1)+pos(1), cen(2)+pos(2)+1, cen(3)+pos(3))
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(2) = shift(2) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! along z
            if( abs(pos(3)) < itrs )then
                alpha = self1%rmat(cen(1)+pos(1), cen(2)+pos(2), cen(3)+pos(3)-1)
                gamma = self1%rmat(cen(1)+pos(1), cen(2)+pos(2), cen(3)+pos(3)+1)
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(3) = shift(3) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! convention
            shift = -shift
        else
            shift = -real(pos)
        endif
    end subroutine fcorr_shift3D

    !> \brief prenorm4real_corr pre-normalises the reference in preparation for real_corr_prenorm
    subroutine prenorm4real_corr_1( self, sxx )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        real :: npix, ax
        npix = real(product(self%ldim))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / npix
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ax
        sxx  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
    end subroutine prenorm4real_corr_1

    !> \brief prenorm4real_corr pre-normalises the reference in preparation for real_corr_prenorm
    subroutine prenorm4real_corr_2( self, sxx, mask )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        logical,      intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: npix, ax
        npix = real(count(mask))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask) / npix
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ax
        sxx  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask)
    end subroutine prenorm4real_corr_2

    !> \brief prenorm4real_corr_3 pre-normalises the reference in preparation for real_corr_prenorm_3
    !>  The image is centered and standardized
    subroutine prenorm4real_corr_3( self, err )
        class(image), intent(inout) :: self
        logical,      intent(out)   :: err
        real :: npix, ave, var
        npix      = real(product(self%ldim))
        ave       = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        self%rmat = self%rmat - ave
        var       = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        if( var > TINY )then
            self%rmat = self%rmat / sqrt(var)
            err = .false.
        else
            err = .true.
        endif
    end subroutine prenorm4real_corr_3

    real function skew( self, mask )
        class(image),      intent(in)  :: self
        logical, optional, intent(in)  :: mask(self%ldim(1),self%ldim(2))
        if( self%is_3d() ) THROW_HARD('2D only!')
        if( self%is_ft() ) THROW_HARD('Real space only!')
        skew = skewness(self%rmat(1:self%ldim(1),1:self%ldim(2),1), mask)
    end function skew

    real function kurt( self, mask )
        class(image),      intent(in)  :: self
        logical, optional, intent(in)  :: mask(self%ldim(1),self%ldim(2))
        if( self%is_3d() ) THROW_HARD('2D only!')
        if( self%is_ft() ) THROW_HARD('Real space only!')
        kurt = kurtosis(self%rmat(1:self%ldim(1),1:self%ldim(2),1), mask)
    end function kurt

    !>  \brief real_corr_prenorm is for calculating a real-space correlation coefficient between images (reference is pre-normalised)
    function real_corr_prenorm_1( self_ref, self_ptcl, sxx_ref ) result( r )
        class(image), intent(inout) :: self_ref, self_ptcl
        real,         intent(in)    :: sxx_ref
        real :: diff(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: r, ay, syy, sxy, npix
        npix = real(product(self_ptcl%ldim))
        ay   = sum(self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3))) / npix
        diff = self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3)) - ay
        syy  = sum(diff * diff)
        sxy  = sum(self_ref%rmat(:self_ref%ldim(1),:self_ref%ldim(2),:self_ref%ldim(3)) * diff)
        if( sxx_ref > 0. .or. syy > 0. )then
            r = sxy / sqrt(sxx_ref * syy)
        else
            r = 0.
        endif
    end function real_corr_prenorm_1

    !>  \brief real_corr_prenorm is for calculating a real-space correlation coefficient between images (reference is pre-normalised)
    function real_corr_prenorm_2( self_ref, self_ptcl, sxx_ref, mask ) result( r )
        class(image), intent(inout) :: self_ref, self_ptcl
        real,         intent(in)    :: sxx_ref
        logical,      intent(in)    :: mask(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: diff(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: r, ay, syy, sxy, npix
        npix = real(count(mask))
        ay   = sum(self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3)), mask=mask) / npix
        where( mask ) diff = self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3)) - ay
        syy  = sum(diff * diff, mask=mask)
        sxy  = sum(self_ref%rmat(:self_ref%ldim(1),:self_ref%ldim(2),:self_ref%ldim(3)) * diff, mask=mask)
        if( sxx_ref > 0. .or. syy > 0. )then
            r = sxy / sqrt(sxx_ref * syy)
        else
            r = 0.
        endif
    end function real_corr_prenorm_2

    !>  \brief real_corr_prenorm_3 is for calculating a real-space correlation coefficient between images
    !>  both inputs are assumed centered & standardized
    real function real_corr_prenorm_3( self_ref, self_ptcl )
        class(image), intent(inout) :: self_ref, self_ptcl
        real_corr_prenorm_3 = sum(self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3))&
                                & * self_ref%rmat(:self_ref%ldim(1),:self_ref%ldim(2),:self_ref%ldim(3)))
        real_corr_prenorm_3 = real_corr_prenorm_3 / product(self_ref%ldim)
    end function real_corr_prenorm_3

    function sqeuclid( self1, self2, mask ) result( r )
        class(image), intent(inout) :: self1, self2
        logical,      intent(in)    :: mask(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: r
        if( self1%wthreads )then
            !$omp parallel workshare
            r = sum((self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
            &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0, mask=mask)
            !$omp end parallel workshare
        else
            r = sum((self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
            &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0, mask=mask)
        endif
    end function sqeuclid

    subroutine sqeuclid_matrix_1( self1, self2, sqdiff )
        class(image), intent(in)    :: self1, self2
        real,         intent(inout) :: sqdiff(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        sqdiff = (self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
        &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0
    end subroutine sqeuclid_matrix_1

    subroutine sqeuclid_matrix_2( self1, self2, sqdiff_img )
        class(image), intent(in)    :: self1, self2
        class(image), intent(inout) :: sqdiff_img
        sqdiff_img%rmat = (self1%rmat - self2%rmat)**2.0
    end subroutine sqeuclid_matrix_2

    function euclid_norm( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real :: r
        r = sqrt(sum((self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
        &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0)) / real(product(self1%ldim))
    end function euclid_norm

    subroutine opt_filter_costfun( even_filt, odd_raw, odd_filt, even_raw, sqdiff_img )
        class(image), intent(in)    :: even_filt, odd_raw, odd_filt, even_raw
        class(image), intent(inout) :: sqdiff_img
        sqdiff_img%rmat = abs(even_filt%rmat - odd_raw%rmat + odd_filt%rmat - even_raw%rmat)
    end subroutine opt_filter_costfun

    subroutine opt_filter_costfun_workshare( even_filt, odd_raw, odd_filt, even_raw, sqdiff_img )
        class(image), intent(in)    :: even_filt, odd_raw, odd_filt, even_raw
        class(image), intent(inout) :: sqdiff_img
        !$omp parallel workshare
        sqdiff_img%rmat = abs(even_filt%rmat - odd_raw%rmat + odd_filt%rmat - even_raw%rmat)
        !$omp end parallel workshare
    end subroutine opt_filter_costfun_workshare

    !> \brief fsc is for calculation of Fourier ring/shell correlation in double precision
    subroutine fsc( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
        real(dp)    :: corrs_8(fdim(self1%ldim(1))-1)
        real(dp)    :: sumasq(fdim(self1%ldim(1))-1), sumbsq(fdim(self1%ldim(1))-1)
        complex(dp) :: comp1, comp2
        integer     :: n, lims(3,2), phys(3), sh, h, k, l
        corrs_8 = 0.d0
        sumasq  = 0.d0
        sumbsq  = 0.d0
        lims    = self1%fit%loop_lims(2)
        n       = self1%get_filtsz()
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,sh,comp1,comp2)&
        !$omp schedule(static) reduction(+:corrs_8,sumasq,sumbsq) proc_bind(close)
        do k=lims(2,1),lims(2,2)
            do h=lims(1,1),lims(1,2)
                do l=lims(3,1),lims(3,2)
                    ! compute physical address
                    phys = self1%fit%comp_addr_phys(h,k,l)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > n ) cycle
                    ! real part of the complex mult btw self1 and targ*
                    comp1 = self1%cmat(phys(1),phys(2),phys(3))
                    comp2 = self2%cmat(phys(1),phys(2),phys(3))
                    corrs_8(sh) = corrs_8(sh)+ real(comp1 * conjg(comp2), kind=dp)
                    sumasq(sh) = sumasq(sh) + csq(comp1)
                    sumbsq(sh) = sumbsq(sh) + csq(comp2)
                end do
            end do
        end do
        !$omp end parallel do
        ! normalize correlations and compute resolutions
        do k=1,n
            if( sumasq(k) > 0.d0 .and. sumbsq(k) > 0.d0 )then
                corrs_8(k) = corrs_8(k)/sqrt(sumasq(k) * sumbsq(k))
            else
                corrs_8(k) = 0.
            endif
        end do
        ! return single-precision corrs
        corrs = real(corrs_8)
    end subroutine fsc

    !> \brief fsc is for calculation of Fourier ring/shell correlation in double precision
    subroutine fsc_scaled( self1, self2, sz, corrs )
        class(image), intent(inout) :: self1, self2
        integer,      intent(in)    :: sz
        real,         intent(out)   :: corrs(sz)
        real(dp)    :: corrs_8(sz), sumasq(sz), sumbsq(sz)
        real        :: scale, rsh
        complex(dp) :: comp1, comp2
        integer     :: lims(3,2), phys(3), sh, sh_sc, h, k, l, n
        scale   = real(2*sz) / real(self1%ldim(1))
        corrs_8 = 0.d0
        sumasq  = 0.d0
        sumbsq  = 0.d0
        lims    = self1%fit%loop_lims(2)
        n       = self1%get_filtsz()
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,sh,sh_sc,rsh,comp1,comp2)&
        !$omp schedule(static) reduction(+:corrs_8,sumasq,sumbsq) proc_bind(close)
        do k=lims(2,1),lims(2,2)
            do h=lims(1,1),lims(1,2)
                do l=lims(3,1),lims(3,2)
                    ! shell
                    rsh = sqrt(real(h*h) + real(k*k) + real(l*l))
                    sh  = nint(rsh)
                    if( sh == 0 .or. sh > n ) cycle
                    sh_sc = nint(scale*rsh)
                    if( sh_sc == 0 .or. sh_sc > sz ) cycle
                    ! compute physical address
                    phys = self1%fit%comp_addr_phys(h,k,l)
                    ! real part of the complex mult btw self1 and targ*
                    comp1          = self1%cmat(phys(1),phys(2),phys(3))
                    comp2          = self2%cmat(phys(1),phys(2),phys(3))
                    corrs_8(sh_sc) = corrs_8(sh_sc)+ real(comp1 * conjg(comp2), kind=dp)
                    sumasq(sh_sc)  = sumasq(sh_sc) + csq_fast(comp1)
                    sumbsq(sh_sc)  = sumbsq(sh_sc) + csq_fast(comp2)
                end do
            end do
        end do
        !$omp end parallel do
        ! normalize correlations and compute resolutions
        where( sumasq>DTINY.and. sumbsq>DTINY )
            corrs = real(corrs_8/dsqrt(sumasq * sumbsq))
        else where
            corrs = 0.
        end where
    end subroutine fsc_scaled

    !>  \brief get array of resolution steps
    function get_res( self ) result( res )
        class(image), intent(in) :: self
        real, allocatable        :: res(:)
        integer                  :: n, k
        n = self%get_filtsz()
        allocate( res(n) )
        do k=1,n
            res(k) = self%fit%get_lp(1,k)
        end do
    end function get_res

    !>  \brief  returns the real and imaginary parts of the phase shift at point
    !!          logi in a Fourier transform caused by the origin shift in shvec
    pure function oshift_1( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        real,         intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
        real    :: arg
        integer :: ldim
        if( self%ldim(3) == 1 )then
            ldim = 2
        else
            ldim = 3
        endif
        arg  = sum(logi(:ldim)*shvec(:ldim)*self%shconst(:ldim))
        comp = cmplx(cos(arg),sin(arg))
    end function oshift_1

    !>  \brief  returns the real and imaginary parts of the phase shift at point
    !!          logi in a Fourier transform caused by the origin shift in shvec
    pure function oshift_2( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
        comp = self%oshift_1(real(logi), shvec)
    end function oshift_2

    !>  \brief  returns the real argument transfer matrix components at point logi in a Fourier transform
    function gen_argtransf_comp( self, logi, ldim ) result( arg )
        class(image), intent(in)      :: self
        real, intent(in)              :: logi(3)
        integer, intent(in), optional :: ldim
        real                          :: arg(3)
        integer                       :: lstop, i
        lstop = 2
        if( self%ldim(3) > 1 ) lstop = 3
        if( present(ldim) )    lstop = ldim
        arg = 0.
        do i=1,lstop
            if( self%ldim(i) == 1 )then
                cycle
            else
                if( is_even(self%ldim(i)) )then
                    arg = arg+logi(i)*(PI/real(self%ldim(i)/2.))
                else
                    arg = arg+logi(i)*(PI/real((self%ldim(i)-1)/2.))
                endif
            endif
        end do
    end function gen_argtransf_comp

    ! MODIFIERS

    subroutine lp_background( self, mskvol, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: mskvol
        real,         intent(in)    :: lp
        type(image) :: weights, self_filt
        if( self%is_ft() ) THROW_HARD('only 4 real images; lp_background')
        if( self%ldim(3) == 1 ) THROW_HARD('only 4 volumes; lp_background')
        if( any((self%ldim-mskvol%ldim)/=0) ) THROW_HARD('inconsistent image/msk dimensions; lp_background')
        call self%zero_background
        call self_filt%new(self%ldim,self%smpd)
        call self_filt%copy(self)
        call weights%new(self%ldim,1.)
        call weights%copy(mskvol)
        ! self
        call self%mul(weights)
        ! self low-pass
        call self_filt%bp(0., lp)
        weights%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
            &1. - weights%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        call self_filt%mul(weights)
        ! addition
        call self%add(self_filt)
        ! clean
        call weights%kill
        call self_filt%kill
    end subroutine lp_background

    subroutine combine_fgbg_filt( self_fg, self_bg, mask )
        class(image), intent(inout) :: self_fg
        class(image), intent(in)    :: self_bg, mask
        self_fg%rmat = mask%rmat * self_fg%rmat + (1. - mask%rmat) * self_bg%rmat
    end subroutine combine_fgbg_filt

    !> \brief insert  inserts a box*box particle image into a micrograph
    !! \param self_in input image
    !! \param coord box coordinates
    !! \param self_out output image
    !!
    subroutine insert(self_in, coord, self_out )
        class(image), intent(in)   :: self_in
        integer, intent(in)        :: coord(2)
        type(image), intent(inout) :: self_out
        integer :: xllim, xulim, yllim, yulim
        if( self_in%ldim(3) > 1 )       THROW_HARD('only 4 2D images; insert')
        if( self_in%is_ft() )           THROW_HARD('only 4 real images; insert')
        if( .not. self_in%even_dims() ) THROW_HARD('only 4 even particle dims; insert')
        if( self_out%exists() )then
            if( self_out%ldim(3) > 1 )  THROW_HARD('only 4 2D images; insert')
            if( self_out%is_ft() )      THROW_HARD('only 4 real images; insert')
            if( self_out%ldim(1) > self_in%ldim(1) .and. self_out%ldim(2) > self_in%ldim(2) .and. self_out%ldim(3) == 1 )then
                if( (coord(1) < self_in%ldim(1)/2+1 .or. coord(1) > self_out%ldim(1)-self_in%ldim(1)/2-1) .or.&
                    (coord(2) < self_in%ldim(2)/2+1 .or. coord(2) > self_out%ldim(2)-self_in%ldim(2)/2-1) )then
                    THROW_HARD('particle outside micrograph area; insert')
                endif
            else
                THROW_HARD('micrograph needs to have dimensions larger than the particle; insert')
            endif
        else
            THROW_HARD('micrograph (self_out) does not exist; insert')
        endif
        ! set range
        xllim = coord(1)-self_in%ldim(1)/2
        xulim = coord(1)+self_in%ldim(1)/2-1
        yllim = coord(2)-self_in%ldim(2)/2
        yulim = coord(2)+self_in%ldim(2)/2-1
        ! insert particle image matrix into micrograph image matrix
        self_out%rmat(xllim:xulim,yllim:yulim,1) = self_in%rmat(1:self_in%ldim(1),1:self_in%ldim(2),1)
    end subroutine insert

    ! inserts the low-resolution information from one image into another
    subroutine insert_lowres( self, self2insert, find )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2insert
        integer,      intent(in)    :: find
        integer :: lims(3,2), phys(3), h, k, l, sh
        complex :: comp
        if( .not. self%ft        ) THROW_HARD('image to be modified assumed to be FTed; insert_lowres')
        if( .not. self2insert%ft ) THROW_HARD('image to insert assumed to be FTed; insert_lowres')
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys,comp)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    if( sh <= find )then
                        ! insert component
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self2insert%get_fcomp([h,k,l],phys)
                        call self%set_fcomp([h,k,l],phys,comp)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine insert_lowres

    ! inserts the low-resolution information from one image into another
    subroutine insert_lowres_serial( self, self2insert, find )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2insert
        integer,      intent(in)    :: find
        integer :: lims(3,2), phys(3), h, k, l, sh
        complex :: comp
        if( .not. self%ft        ) THROW_HARD('image to be modified assumed to be FTed; insert_lowres')
        if( .not. self2insert%ft ) THROW_HARD('image to insert assumed to be FTed; insert_lowres')
        lims = self%fit%loop_lims(2)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    if( sh <= find )then
                        ! insert component
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self2insert%get_fcomp([h,k,l],phys)
                        call self%set_fcomp([h,k,l],phys,comp)
                    endif
                end do
            end do
        end do
    end subroutine insert_lowres_serial

    !>  \brief  is for inverting an image
    subroutine inv( self )
        class(image), intent(inout) :: self
        self%rmat = -1.*self%rmat
    end subroutine inv

    subroutine zero_neg( self )
        class(image), intent(inout) :: self
        where( self%rmat < TINY ) self%rmat = 0.
    end subroutine zero_neg

    subroutine remove_neg( self )
        class(image), intent(inout) :: self
        real :: minv
        minv = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        if( minv < 0. )self%rmat = self%rmat + abs(minv)
    end subroutine remove_neg

    !>  \brief  is for making a random image (0,b)
    subroutine ran( self, b )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: b      ! support upper bound, defaults to 1
        integer :: i, j, k
        do k=1,self%ldim(3)
            do j=1,self%ldim(2)
                do i=1,self%ldim(1)
                    self%rmat(i,j,k) = ran3()
                end do
            end do
        end do
        self%ft = .false.
        if( present(b) ) self%rmat = self%rmat * b
    end subroutine ran

    !> \brief gauran  is for making a Gaussian random image (0,1)
    !! \param mean Mean of noise
    !! \param sdev Standard deviation of noise
    !!
    subroutine gauran( self, mean, sdev )
        class(image), intent(inout) :: self
        real, intent(in) :: mean, sdev
        integer :: i, j, k
        do k=1,self%ldim(3)
            do j=1,self%ldim(2)
                do i=1,self%ldim(1)
                    self%rmat(i,j,k) = gasdev( mean, sdev )
                end do
            end do
        end do
        self%ft = .false.
    end subroutine gauran

    !> \brief add_gauran  is for adding Gaussian noise to an image
    !! \param snr signal-to-noise ratio
    !! \param noiseimg output image
    !!
    subroutine add_gauran( self, snr )
        class(image), intent(inout) :: self
        real,         intent(in)    :: snr
        real    :: sdev_noise, var
        integer :: i, j, k
        var = self%variance()
        if( var > TINY )then
            sdev_noise = sqrt(var/snr)
        else
            THROW_HARD('variance of image is zero')
        endif
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    self%rmat(i,j,k) = self%rmat(i,j,k) + gasdev(0., sdev_noise)
                end do
            end do
        end do
    end subroutine add_gauran

    !>  \brief dead_hot_positions is for generating dead/hot pixel positions in an image
    !! \param frac fraction of ON/OFF pixels
    !! \return  pos binary 2D map
    !!
    function dead_hot_positions( self, frac ) result( pos )
        class(image), intent(in) :: self
        real, intent(in)         :: frac
        logical, allocatable     :: pos(:,:)
        integer :: ipix, jpix, cnt
        allocate(pos(self%ldim(1),self%ldim(2)))
        pos = .false.
        cnt = 0
        do ipix=1,self%ldim(1)
            do jpix=1,self%ldim(2)
                if( ran3() <= frac )then
                    pos(ipix,jpix) = .true.
                    cnt = cnt+1
                endif
            end do
        end do
    end function dead_hot_positions

    !>  \brief zero image
    subroutine zero(self)
        class(image), intent(inout) :: self
        if( self%ft )then
            self%cmat = cmplx(0.,0.)
        else
            self%rmat = 0.
        endif
    end subroutine zero

    !>  \brief zero image
    subroutine zero_and_unflag_ft(self)
        class(image), intent(inout) :: self
        self%rmat = 0.
        self%ft   = .false.
    end subroutine zero_and_unflag_ft

    !>  \brief zero image
    subroutine zero_and_flag_ft(self)
        class(image), intent(inout) :: self
        self%cmat = cmplx(0.,0.)
        self%ft   = .true.
    end subroutine zero_and_flag_ft

    ! estimates median of background along edges of box and subtracts it to flatten the background
    ! modifies pixels along the edges of the box, which ought to be safe as we are masking
    subroutine zero_background( self )
        class(image), intent(inout) :: self
        integer :: k
        real    :: med, val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12
        k = self%ldim(1)/2
        if( self%ldim(3) == 1 )then
            val1  = selec(k,self%ldim(1),self%rmat( 1           , :self%ldim(2),1))
            val2  = selec(k,self%ldim(1),self%rmat( self%ldim(1), :self%ldim(2),1))
            val3  = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  1,           1))
            val4  = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  self%ldim(2),1))
            med   = (val1+val2+val3+val4) / 4.
        else
            val1  = selec(k,self%ldim(1),self%rmat( 1           ,  1            , :self%ldim(3)))
            val2  = selec(k,self%ldim(1),self%rmat( 1           ,  self%ldim(2) , :self%ldim(3)))
            val3  = selec(k,self%ldim(1),self%rmat( self%ldim(1),  1            , :self%ldim(3)))
            val4  = selec(k,self%ldim(1),self%rmat( self%ldim(1),  self%ldim(2) , :self%ldim(3)))
            val5  = selec(k,self%ldim(1),self%rmat( 1           , :self%ldim(2) ,  1           ))
            val6  = selec(k,self%ldim(1),self%rmat( 1           , :self%ldim(2) ,  self%ldim(3)))
            val7  = selec(k,self%ldim(1),self%rmat( self%ldim(1), :self%ldim(2) ,  1           ))
            val8  = selec(k,self%ldim(1),self%rmat( self%ldim(1), :self%ldim(2) ,  self%ldim(3)))
            val9  = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  1            ,  1           ))
            val10 = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  1            ,  self%ldim(3)))
            val11 = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  self%ldim(2) ,  1           ))
            val12 = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  self%ldim(2) ,  self%ldim(3)))
            med   = (val1+val2+val3+val4+val5+val6+val7+val8+val9+val10+val11+val12) / 12.
        endif
        if(abs(med) > TINY) self%rmat = self%rmat - med
    end subroutine zero_background

    ! substracts median of background defined by an enveloppe mask (eg 0< <1)
    subroutine zero_env_background( self, volmsk )
        class(image), intent(inout) :: self, volmsk
        real, allocatable :: vals(:)
        real              :: med
        integer           :: cnt, npix, npix_env, i,j,k
        logical           :: env_mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        npix = product(self%ldim)
        where( volmsk%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) > 0.0001&
        &.and. volmsk%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) < 0.9999 )
            env_mask = .true.
        else where
            env_mask = .false.
        end where
        npix_env = count(env_mask)
        if( npix_env==0 .or. npix_env==npix )then
            THROW_HARD('Volume mask is not a volume mask; simple_image :: zero_env_background')
        endif
        allocate(vals(npix_env))
        cnt = 0
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if(env_mask(i,j,k))then
                        cnt = cnt + 1
                        vals(cnt) = self%rmat(i,j,k)
                    endif
                enddo
            enddo
        enddo
        med = median_nocopy(vals)
        if(abs(med) > TINY) self%rmat = self%rmat - med
        deallocate(vals)
    end subroutine zero_env_background

    !>  \brief generates instrument function divided image
    subroutine div_w_instrfun( self, interpfun, alpha, padded_dim )
        class(image),           intent(inout) :: self
        character(len =*),      intent(in)    :: interpfun
        real,         optional, intent(in)    :: alpha
        integer,      optional, intent(in)    :: padded_dim
        type(kbinterpol)  :: kbwin
        real, allocatable :: w(:)
        real    :: arg
        integer :: center(3), i,j,k, dim, iarg
        if( any(self%ldim==0) .or. self%is_ft() .or. .not.self%square_dims() )then
            THROW_HARD('Erroneous image in div_w_instrfun')
        endif
        center = self%ldim/2+1
        dim    = self%ldim(1)
        if( present(padded_dim) ) dim = padded_dim
        select case(trim(interpfun))
        case('kb')
            ! kaiser-bessel window
            if(.not.present(alpha)) THROW_HARD('alpha must be given for KB interpolator')
            kbwin = kbinterpol(KBWINSZ,alpha)
            allocate(w(self%ldim(1)),source=1.)
            do i = 1,self%ldim(1)
                arg  = real(i-center(1))/real(dim)
                w(i) = kbwin%instr(arg)
            end do
            if( self%is_2d() )then
                !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close) schedule(static)
                do i = 1,self%ldim(1)
                    do j = 1,self%ldim(2)
                        self%rmat(i,j,1) = self%rmat(i,j,1) / (w(i)*w(j))
                    enddo
                enddo
                !$omp end parallel do
            else
                !$omp parallel do collapse(3) private(i,j,k) default(shared) proc_bind(close) schedule(static)
                do i = 1,self%ldim(1)
                    do j = 1,self%ldim(2)
                        do k = 1,self%ldim(3)
                            self%rmat(i,j,k) = self%rmat(i,j,k) / (w(i)*w(j)*w(k))
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            endif
        case('linear')
            ! Tri-linear interpolation
            !$omp parallel do collapse(3) private(i,j,k,iarg,arg) default(shared) proc_bind(close) schedule(static)
            do i = 1,self%ldim(1)
                do j = 1,self%ldim(2)
                    do k = 1,self%ldim(3)
                        iarg = sum(([i,j,k]-center)**2)
                        if( iarg == 0 )cycle
                        arg = PI * sqrt(real(iarg)) / real(dim)
                        arg = sin(arg) / arg
                        arg = arg*arg ! normalized sinc^2
                        self%rmat(i,j,k) = self%rmat(i,j,k) / arg
                    enddo
                enddo
            enddo
            !$omp end parallel do
        case('nn')
            ! Nearest-neighbour interpolation
            !$omp parallel do collapse(3) private(i,j,k,iarg,arg) default(shared) proc_bind(close) schedule(static)
            do i = 1,self%ldim(1)
                do j = 1,self%ldim(2)
                    do k = 1,self%ldim(3)
                        iarg = sum(([i,j,k]-center)**2)
                        if( iarg == 0 )cycle
                        arg = PI * sqrt(real(iarg)) / real(dim)
                        arg = sin(arg) / arg ! normalized sinc
                        self%rmat(i,j,k) = self%rmat(i,j,k) / arg
                    enddo
                enddo
            enddo
            !$omp end parallel do
        case DEFAULT
            THROW_HARD('Unsupported interpolation method')
        end select
    end subroutine div_w_instrfun

    subroutine pad_fft( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine pad_fft

    subroutine norm_noise_pad_fft( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        real    :: sdev_noise
        call self%norm_noise(lmsk, sdev_noise)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine norm_noise_pad_fft

    !>  \brief  Taper edges of image so that there are no sharp discontinuities in real space
    !!          This is a re-implementation of the MRC program taperedgek.for (Richard Henderson, 1987)
    !!          I stole it from CTFFIND4 (thanks Alexis for the beautiful re-implementation)
    subroutine taper_edges( self )
        class(image), intent(inout) :: self
        real, allocatable  :: avg_curr_edge_start(:,:)
        real, allocatable  :: avg_curr_edge_stop(:,:)
        real, allocatable  :: avg_curr_edge_avg(:,:)
        real, allocatable  :: smooth_avg_curr_edge_start(:,:)
        real, allocatable  :: smooth_avg_curr_edge_stop(:,:)
        integer            :: curr_dim, ndims
        integer            :: dim2, dim3
        integer            :: i,j,k
        integer            :: j_shift,k_shift
        integer            :: jj,kk
        integer            :: nvals_runnavg
        integer, parameter :: avg_strip_width(3)   = 100
        integer, parameter :: taper_strip_width(3) = 500
        integer, parameter :: smooth_half_width(3) = 1
        ndims = 2
        ! initialise vars
        dim2 = 2;  dim3 = 3; nvals_runnavg = 0
        if (self%is_3d()) ndims = 3
        do curr_dim=1,ndims
            ! take care of dimensions
            select case (curr_dim)
            case (1)
                dim2 = 2
                dim3 = 3
            case (2)
                dim2 = 1
                dim3 = 3
            case (3)
                dim2 = 1
                dim3 = 2
            end select
            ! take care of allocation & initialisation
            if(allocated(avg_curr_edge_start))        deallocate(avg_curr_edge_start)
            if(allocated(avg_curr_edge_stop))         deallocate(avg_curr_edge_stop)
            if(allocated(avg_curr_edge_avg))          deallocate(avg_curr_edge_avg)
            if(allocated(smooth_avg_curr_edge_start)) deallocate(smooth_avg_curr_edge_start)
            if(allocated(smooth_avg_curr_edge_stop))  deallocate(smooth_avg_curr_edge_stop)
            allocate( avg_curr_edge_start(self%ldim(dim2),self%ldim(dim3)),&
                avg_curr_edge_stop (self%ldim(dim2),self%ldim(dim3)),&
                avg_curr_edge_avg(self%ldim(dim2),self%ldim(dim3)),&
                smooth_avg_curr_edge_start(self%ldim(dim2),self%ldim(dim3)),&
                smooth_avg_curr_edge_stop (self%ldim(dim2),self%ldim(dim3)))
            avg_curr_edge_start        = 0.0e0
            avg_curr_edge_stop         = 0.0e0
            avg_curr_edge_avg          = 0.0e0
            smooth_avg_curr_edge_start = 0.0e0
            smooth_avg_curr_edge_stop  = 0.0e0
            ! Deal with X=0 and X=self%ldim(1) edges
            i=1
            do k=1,self%ldim(dim3)
                do j=1,self%ldim(dim2)
                    select case (curr_dim)
                    case (1)
                        avg_curr_edge_start(j,k) =&
                            sum(self%rmat(1:avg_strip_width(curr_dim),j,k))&
                            /avg_strip_width(curr_dim)
                        avg_curr_edge_stop(j,k)  =&
                            sum(self%rmat(self%ldim(curr_dim)-avg_strip_width(1)+1:self%ldim(curr_dim),j,k))&
                            /avg_strip_width(curr_dim)
                    case (2)
                        avg_curr_edge_start(j,k) =&
                            sum(self%rmat(j,1:avg_strip_width(curr_dim),k))&
                            /avg_strip_width(curr_dim)
                        avg_curr_edge_stop(j,k) =&
                            sum(self%rmat(j,self%ldim(curr_dim)-avg_strip_width(1)+1:self%ldim(curr_dim),k))&
                            /avg_strip_width(curr_dim)
                    case (3)
                        avg_curr_edge_start(j,k) =&
                            sum(self%rmat(j,k,1:avg_strip_width(curr_dim)))&
                            /avg_strip_width(curr_dim)
                        avg_curr_edge_stop(j,k) =&
                            sum(self%rmat(j,k,self%ldim(curr_dim)-avg_strip_width(1)+1:self%ldim(curr_dim)))&
                            /avg_strip_width(curr_dim)
                    end select
                enddo
            enddo
            avg_curr_edge_avg   = 0.5e0*(avg_curr_edge_stop + avg_curr_edge_start)
            avg_curr_edge_start = avg_curr_edge_start - avg_curr_edge_avg
            avg_curr_edge_stop  = avg_curr_edge_stop - avg_curr_edge_avg
            ! Apply smoothing parallel to edge in the form of a running average
            do k=1,self%ldim(dim3)
                do j=1,self%ldim(dim2)
                    nvals_runnavg = 0
                    ! Loop over neighbourhood of non-smooth arrays
                    do k_shift=-smooth_half_width(dim3),smooth_half_width(dim3)
                        kk = k+k_shift
                        if (kk .lt. 1 .or. kk .gt. self%ldim(dim3)) cycle
                        do j_shift=-smooth_half_width(dim2),smooth_half_width(dim2)
                            jj = j+j_shift
                            if (jj .lt. 1 .or. jj .gt. self%ldim(dim2)) cycle
                            nvals_runnavg = nvals_runnavg + 1
                            smooth_avg_curr_edge_start (j,k) =&
                                smooth_avg_curr_edge_start(j,k)+avg_curr_edge_start(jj,kk)
                            smooth_avg_curr_edge_stop(j,k)   =&
                                smooth_avg_curr_edge_stop(j,k)+avg_curr_edge_stop(jj,kk)
                        enddo
                    enddo
                    ! Now we can compute the average
                    smooth_avg_curr_edge_start(j,k) = smooth_avg_curr_edge_start(j,k)/nvals_runnavg
                    smooth_avg_curr_edge_stop(j,k)   = smooth_avg_curr_edge_stop(j,k)/nvals_runnavg
                enddo
            enddo
            ! Taper the image
            do i=1,self%ldim(curr_dim)
                if (i .le. taper_strip_width(curr_dim)) then
                    select case (curr_dim)
                    case (1)
                        self%rmat(i,:,:) = self%rmat(i,:,:)&
                            - smooth_avg_curr_edge_start (:,:)&
                            * (taper_strip_width(curr_dim)-i+1)&
                            / taper_strip_width(curr_dim)
                    case (2)
                        self%rmat(1:self%ldim(1),i,:) = self%rmat(1:self%ldim(1),i,:)&
                            - smooth_avg_curr_edge_start(:,:)&
                            * (taper_strip_width(curr_dim)-i+1)&
                            / taper_strip_width(curr_dim)
                    case (3)
                        self%rmat(1:self%ldim(1),:,i) = self%rmat(1:self%ldim(1),:,i)&
                            - smooth_avg_curr_edge_start (:,:)&
                            * (taper_strip_width(curr_dim)-i+1)&
                            / taper_strip_width(curr_dim)
                    end select
                else if (i .ge. self%ldim(curr_dim)-taper_strip_width(curr_dim)+1) then
                    select case (curr_dim)
                    case (1)
                        self%rmat(i,:,:) = self%rmat(i,:,:)&
                            - smooth_avg_curr_edge_stop(:,:)&
                            * (taper_strip_width(curr_dim)+i&
                            - self%ldim(curr_dim))&
                            / taper_strip_width(curr_dim)
                    case (2)
                        self%rmat(1:self%ldim(1),i,:) = self%rmat(1:self%ldim(1),i,:)&
                            - smooth_avg_curr_edge_stop(:,:)&
                            * (taper_strip_width(curr_dim)+i&
                            - self%ldim(curr_dim))&
                            / taper_strip_width(curr_dim)
                    case (3)
                        self%rmat(1:self%ldim(1),:,i) = self%rmat(1:self%ldim(1),:,i)&
                            - smooth_avg_curr_edge_stop(:,:)&
                            * (taper_strip_width(curr_dim)+i&
                            - self%ldim(curr_dim))&
                            / taper_strip_width(curr_dim)
                    end select
                endif
            enddo
        enddo
        if(allocated(avg_curr_edge_start))        deallocate(avg_curr_edge_start)
        if(allocated(avg_curr_edge_stop))         deallocate(avg_curr_edge_stop)
        if(allocated(avg_curr_edge_avg))          deallocate(avg_curr_edge_avg)
        if(allocated(smooth_avg_curr_edge_start)) deallocate(smooth_avg_curr_edge_start)
        if(allocated(smooth_avg_curr_edge_stop))  deallocate(smooth_avg_curr_edge_stop)
    end subroutine taper_edges

    subroutine taper_edges_hann( self, borders )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: borders(2)
        real    :: w
        integer :: i,j,n
        if( self%is_ft() )THROW_HARD('Real space only!')
        if( self%is_3d() )THROW_HARD('2D images only!')
        n = 0
        do j = 1,borders(2)
            n = n+1
            w = 1. - cos(PIO2* real(j-1)/real(borders(2)) )
            w = min(1.,max(0.,w))
            self%rmat(1:self%ldim(1),j,1) = w * self%rmat(1:self%ldim(1),j,1)
        enddo
        n = 0
        do j = self%ldim(2)-borders(2)+1,self%ldim(2)
            n = n+1
            w = cos(PIO2* real(j-self%ldim(2)+borders(2))/real(borders(2)))
            w = min(1.,max(0.,w))
            self%rmat(1:self%ldim(1),j,1) = w * self%rmat(1:self%ldim(1),j,1)
        enddo
        do i = 1,borders(1)
            w = 1. - cos(PIO2* real(i-1)/real(borders(1)) )
            w = min(1.,max(0.,w))
            self%rmat(i,1:self%ldim(2),1) = w * self%rmat(i,1:self%ldim(2),1)
        enddo
        do i = self%ldim(1)-borders(1)+1,self%ldim(1)
            w = cos(PIO2* real(i-self%ldim(1)+borders(1))/real(borders(1)))
            w = min(1.,max(0.,w))
            self%rmat(i,1:self%ldim(2),1) = w * self%rmat(i,1:self%ldim(2),1)
        enddo
    end subroutine taper_edges_hann

    !>  subtracts background linear ramp including mean
    subroutine subtr_backgr_ramp( self, lmsk )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable :: xyz(:,:)
        real    :: A,B,C,D
        integer :: cen(2),i,j,npix
        logical :: err
        if( self%ft )      THROW_HARD('Real space only!, subtr_backr_ramp')
        if( self%is_3d() ) THROW_HARD('2D images only!, subtr_backr_ramp')
        npix = product(self%ldim) - count(lmsk)
        if( npix < 2 ) return
        allocate(xyz(npix,3))
        cen  = self%ldim(1:2)/2 + 1
        npix = 0
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                if( lmsk(i,j,1) ) cycle
                npix = npix + 1
                xyz(npix,:) = [real(i-cen(1)), real(j-cen(2)), self%rmat(i,j,1)]
            enddo
        enddo
        call fit_lsq_plane(npix, xyz, A,B,C, err)
        if( err ) return
        do j = 1,self%ldim(2)
            D = B*real(j-cen(2)) + C
            do i = 1,self%ldim(1)
                self%rmat(i,j,1) = self%rmat(i,j,1) - (D + A*real(i-cen(1)))
            enddo
        enddo
    end subroutine subtr_backgr_ramp

    subroutine upsample_square_background( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(2)
        type(image) :: upsampled
        integer :: i,j,fx,fy,fpx,fpy
        real :: scales(2), scale, x,y, dx,dy,b
        if( ldim(1) /= ldim(2) ) THROW_HARD('Unsupported dimensions!')
        scales = real(self%ldim(1:2)) / real(ldim(1:2))
        scale  = product(scales)
        call upsampled%new([ldim(1), ldim(2), 1], self%smpd/scales(1))
        !$omp parallel do default(shared) collapse(2) proc_bind(close)&
        !$omp private(i,j,x,y,fx,fy,fpx,fpy,dx,dy,b)
        do j = 1,ldim(1)
            do i = 1,ldim(2)
                x   = (real(i-1) * scales(1)) + 1.
                y   = (real(j-1) * scales(2)) + 1.
                fx  = floor(x)
                fy  = floor(y)
                if( fx == self%ldim(1) )then
                    fpx = fx-1
                else
                    fpx = fx+1
                endif
                if( fy == self%ldim(2) )then
                    fpy = fy-1
                else
                    fpy = fy+1
                endif
                dx = x - real(fx)
                dy = y - real(fy)
                b =     self%rmat(fx, fy, 1) * (1.-dx) * (1.-dy)
                b = b + self%rmat(fpx,fy, 1) *     dx  * (1.-dy)
                b = b + self%rmat(fx, fpy,1) * (1.-dx) * dy
                b = b + self%rmat(fpx,fpy,1) *     dx  * dy
                upsampled%rmat(i,j,1) = b * scale
            enddo
        enddo
        !$omp end parallel do
        call self%copy(upsampled)
        call upsampled%kill
    end subroutine upsample_square_background

    !>  Subtracts background, for micrographs
    subroutine subtract_background( self, freq, mode )
        class(image),               intent(inout) :: self
        real,                       intent(in)    :: freq
        character(len=*), optional, intent(in)    :: mode
        integer,    parameter :: CS_DIM  = 1024
        type(image)           :: tmpimg
        character(len=STDLEN) :: cmode
        if( self%ft )      THROW_HARD('Real space only!, subtract_background')
        if( self%is_3d() ) THROW_HARD('2D images only!, subtract_background')
        cmode = trim(NIL)
        if( present(mode) ) cmode = mode
        call self%estimate_background(freq, tmpimg, cmode)
        select case(trim(mode))
        case('cryosparc','cs')
            call tmpimg%upsample_square_background(self%ldim(1:2))
        case DEFAULT
            ! all done
        end select
        !$omp parallel workshare proc_bind(close)
        self%rmat = self%rmat - tmpimg%rmat
        !$omp end parallel workshare
        call tmpimg%kill
    end subroutine subtract_background

    !>  Estimates background from gaussian filtered image, for micrographs
    subroutine estimate_background( self, freq, backgr, mode )
        class(image),     intent(in)    :: self
        real,             intent(in)    :: freq
        class(image),     intent(inout) :: backgr
        character(len=*), intent(in)    :: mode
        real,    parameter :: PADDING = sqrt(2.)
        integer, parameter :: CS_DIM  = 1024
        type(image)     :: img_pad, msk
        real            :: smpd_bin
        integer         :: ldim_pd(3), ldim_bin(3), ldim_cs(3), bin_dim
        integer         :: i,j,is,js,ie,je, binning
        if( self%ft )      THROW_HARD('Real space only!, estimate_background')
        if( self%is_3d() ) THROW_HARD('2D images only!, estimate_background')
        select case(trim(mode))
        case('cryosparc','cs')
            ! produces backround of shape CS_DIM x CS_DIM
            ! bin micrograph
            ldim_cs = [CS_DIM, CS_DIM,1]
            bin_dim = 2**ceiling( log(real(maxval(self%ldim))) / log(2.0) )
            binning = ceiling(real(bin_dim) / real(CS_DIM))
            ldim_bin(1:2) = nint(self%ldim(1:2) / real(binning))
            ldim_bin(3)   = 1
            smpd_bin      = self%smpd * real(binning)
            call backgr%new(ldim_bin, smpd_bin)
            !$omp parallel do private(i,j,is,ie,js,je) proc_bind(close) collapse(2) default(shared)
            do j = 1,ldim_bin(2)
                do i = 1,ldim_bin(1)
                    is = (i-1)*binning + 1
                    ie = i*binning
                    js = (j-1)*binning + 1
                    je = j*binning
                    backgr%rmat(i,j,1) = sum(self%rmat(is:ie,js:je,1))
                enddo
            enddo
            !$omp end parallel do
            ! low-pass padded mic & crop
            ldim_pd = [2*CS_DIM, 2*CS_DIM, 1]
            call img_pad%new(ldim_pd, smpd_bin)
            call backgr%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0., 2.*freq)
            call img_pad%ifft
            call img_pad%clip(backgr)
            ! low pass padded mask & crop
            call msk%new(ldim_cs,smpd_bin)
            msk = 1.
            call msk%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0., 2.*freq)
            call img_pad%ifft
            call img_pad%clip(msk)
            ! correct for padding
            !$omp parallel workshare proc_bind(close)
            where( msk%rmat(1:ldim_cs(1),1:ldim_cs(2),1) > TINY )
                backgr%rmat(1:ldim_cs(1),1:ldim_cs(2),1) = &
                    &backgr%rmat(1:ldim_cs(1),1:ldim_cs(2),1) / msk%rmat(1:ldim_cs(1),1:ldim_cs(2),1)
            else where
                backgr%rmat(1:ldim_cs(1),1:ldim_cs(2),1) = 0.
            end where
            !$omp end parallel workshare
        case DEFAULT
            ! produces backround of shape self%ldim
            ! padded dimensions
            ldim_pd(1:2) = nint(PADDING*real(self%ldim(1:2)))
            ldim_pd(1:2) = find_larger_magic_box(ldim_pd(1:2))
            ldim_pd(3)   = 1
            ! low-pass padded mic & crop
            call backgr%copy(self)
            call img_pad%new(ldim_pd, self%smpd)
            call backgr%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0.,freq)
            call img_pad%ifft
            call img_pad%clip(backgr)
            ! low pass padded mask & crop
            call msk%new(self%ldim, self%smpd)
            msk = 1.
            call msk%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0.,freq)
            call img_pad%ifft
            call img_pad%clip(msk)
            ! correct for padding
            !$omp parallel workshare proc_bind(close)
            where( msk%rmat(1:self%ldim(1),1:self%ldim(2),1) > TINY )
                backgr%rmat(1:self%ldim(1),1:self%ldim(2),1) = &
                    &backgr%rmat(1:self%ldim(1),1:self%ldim(2),1) / msk%rmat(1:self%ldim(1),1:self%ldim(2),1)
            else where
                backgr%rmat(1:self%ldim(1),1:self%ldim(2),1) = 0.
            end where
            !$omp end parallel workshare
        end select
        call img_pad%kill
        call msk%kill
    end subroutine estimate_background

    !> \brief salt_n_pepper  is for adding salt and pepper noise to an image
    !! \param pos 2D mask
    subroutine salt_n_pepper( self, pos )
        class(image), intent(inout) :: self
        logical, intent(in)         :: pos(:,:)
        integer :: ipix, jpix
        if( .not. self%is_2d() ) THROW_HARD('only for 2D images; salt_n_pepper')
        call self%norm_minmax
        do ipix=1,self%ldim(1)
            do jpix=1,self%ldim(2)
                if( pos(ipix,jpix) )then
                    if( ran3() < 0.5 )then
                        self%rmat(ipix,jpix,1) = 0.
                    else
                        self%rmat(ipix,jpix,1) = 1.
                    endif
                endif
            end do
        end do
    end subroutine salt_n_pepper

    !>  \brief square just a binary square for testing purposes
    !! \param sqrad half width of square
    subroutine square( self, sqrad )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: sqrad
        integer :: i, j, k
        self%rmat = 0.
        if( all(self%ldim(1:2) .gt. sqrad) .and. self%ldim(3) == 1 ) then
            do i=self%ldim(1)/2-sqrad+1,self%ldim(1)/2+sqrad
                do j=self%ldim(2)/2-sqrad+1,self%ldim(2)/2+sqrad
                    self%rmat(i,j,1) = 1.
                end do
            end do
        else if( all(self%ldim .gt. sqrad) .and. self%ldim(3) > 1 )then
            do i=self%ldim(1)/2-sqrad+1,self%ldim(1)/2+sqrad
                do j=self%ldim(2)/2-sqrad+1,self%ldim(2)/2+sqrad
                    do k=self%ldim(3)/2-sqrad+1,self%ldim(3)/2+sqrad
                        self%rmat(i,j,k) = 1.
                    end do
                end do
            end do
        else
            THROW_HARD('image is to small to fit the square; square')
        endif
        self%ft = .false.
    end subroutine square

    !>  \brief  just a corner filling fun for testing purposes
    !! \param sqrad half width of square
    subroutine corners( self, sqrad )
        class(image), intent(inout) :: self
        integer, intent(in)         :: sqrad
        integer :: i, j
        self%rmat = 0.
        do i=self%ldim(1)-sqrad+1,self%ldim(1)
            do j=self%ldim(2)-sqrad+1,self%ldim(2)
                self%rmat(i,j,1) = 1.
            end do
        end do
        do i=1,sqrad
            do j=1,sqrad
                self%rmat(i,j,1) = 1.
            end do
        end do
        do i=self%ldim(1)-sqrad+1,self%ldim(1)
            do j=1,sqrad
                self%rmat(i,j,1) = 1.
            end do
        end do
        do i=1,sqrad
            do j=self%ldim(2)-sqrad+1,self%ldim(2)
                self%rmat(i,j,1) = 1.
            end do
        end do
        self%ft = .false.
    end subroutine corners

    !> before_after to generate a before (left) and after (right) image
    !! \param left,right input images
    !! \return ba output montage
    subroutine before_after( left, right, ba, mask )
        class(image),      intent(in)    :: left, right
        type(image),       intent(inout) :: ba
        logical, optional, intent(in)    :: mask(left%ldim(1),left%ldim(2),left%ldim(3))
        integer     :: ldim(3), i, j
        if( left.eqdims.right )then
            if( left.eqsmpd.right )then
                if( left%ft .or. right%ft ) THROW_HARD('not for FTs; before_after')
                if( left%is_3d() .or. right%is_3d() ) THROW_HARD('not for 3D imgs; before_after')
                ldim = left%ldim
                ba = left
                ba%rmat(:ldim(1)/2,:ldim(2),1)   = left%rmat(:ldim(1)/2,:ldim(2),1)
                ba%rmat(ldim(1)/2+1:,:ldim(2),1) = 0.
                if( present(mask) )then
                    do i=ldim(1)/2+1,ldim(1)
                        do j=1,ldim(2)
                            if( mask(i,j,1) )then
                                ba%rmat(i,j,1) = right%rmat(i,j,1)
                            endif
                        end do
                    end do
                else
                    ba%rmat(ldim(1)/2+1:,:ldim(2),1) = right%rmat(ldim(1)/2+1:,:ldim(2),1)
                endif
                ba%rmat(ldim(1)/2:ldim(1)/2+1,:,1) = 0.
            else
                THROW_HARD('before (left) and after (right) not of same smpd; before_after')
            endif
        else
            THROW_HARD('before (left) and after (right) not of same dim; before_after')
        endif
    end subroutine before_after

    !> \brief gauimg  just a Gaussian fun for testing purposes
    !! \param wsz window size
    subroutine gauimg_1( self, wsz , alpha )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: wsz
        real, intent(in), optional  :: alpha
        real    :: x, y, z, xw, yw, zw, a
        integer :: i, j, k
        a = 0.5
        if(present(alpha)) a= alpha
        x = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            xw = gauwfun(x, a*real(wsz))
            y = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                yw = gauwfun(y, a*real(wsz))
                z = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        zw = gauwfun(z, a*real(wsz))
                    else
                        zw = 1.
                    endif
                    self%rmat(i,j,k) = xw*yw*zw
                    z = z+1.
                end do
                y = y+1.
            end do
            x = x+1.
        end do
        self%ft = .false.
    end subroutine gauimg_1

    !> \brief gauimg  just a Gaussian fun for testing purposes
    !! \param wsz window size
    subroutine gauimg_2( self, wsz, offx,offy)
        class(image), intent(inout) :: self
        integer, intent(in) :: wsz, offx, offy
        real    :: x, y, z, xw, yw, zw
        integer :: i, j, k
        x = -real(self%ldim(1))/2. -real(offx)
        do i=1,self%ldim(1)
            xw = gauwfun(x, 0.5*real(wsz))
            y = -real(self%ldim(2))/2. - real(offy)
            do j=1,self%ldim(2)
                yw = gauwfun(y, 0.5*real(wsz))
                z = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        zw = gauwfun(z, 0.5*real(wsz))
                    else
                        zw = 1.
                    endif
                    self%rmat(i,j,k) =  self%rmat(i,j,k) + xw*yw*zw
                    z = z+1.
                end do
                y = y+1.
            end do
            x = x+1.
        end do
        self%ft = .false.
    end subroutine gauimg_2

    subroutine gauimg2D( self, xsigma, ysigma, cutoff )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: xsigma, ysigma
        real, optional, intent(in)    :: cutoff
        real    :: xx, x, y, cutoffsq_here, center(2)
        integer :: i, j
        cutoffsq_here = huge(cutoffsq_here)
        if(present(cutoff)) cutoffsq_here = cutoff*cutoff
        ! Center of the img is assumed self%ldim/2 + 1
        center = real(self%ldim(1:2))/2.+1.
        do i=1,self%ldim(1)
            x = real(i)-center(1)
            xx = x*x
            do j=1,self%ldim(2)
                y = real(j)-center(2)
                if(xx+y*y > cutoffsq_here) cycle
                self%rmat(i,j,1) = gaussian2D( [0.,0.], x, y, xsigma, ysigma )
            enddo
        enddo
        self%ft = .false.
    end subroutine gauimg2D

    subroutine gauimg3D( self, xsigma, ysigma, zsigma, cutoff )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: xsigma, ysigma, zsigma
        real, optional, intent(in)    :: cutoff
        real    :: xx, x, yy, y, z, cutoffsq_here, center(3)
        integer :: i, j, k
        if(self%ldim(3) == 1) then
            call self%gauimg2D(xsigma, ysigma, cutoff)
            return
        endif
        cutoffsq_here = huge(cutoffsq_here)
        if(present(cutoff)) cutoffsq_here = cutoff*cutoff
        ! Center of the img is assumed self%ldim/2 + 1
        center = real(self%ldim(1:3))/2.+1.
        do i=1,self%ldim(1)
            x = real(i)-center(1)
            xx = x*x
            do j=1,self%ldim(2)
                y = real(j)-center(2)
                yy = y*y
                do k=1,self%ldim(3)
                    z = real(k)-center(3)
                    if(xx+yy+z*z > cutoffsq_here) cycle
                    self%rmat(i,j,k) = gaussian3D( [0.,0.,0.], x, y, z, xsigma, ysigma, zsigma )
                enddo
            enddo
        enddo
        self%ft = .false.
    end subroutine gauimg3D

    subroutine fwd_ft(self)
        class(image), intent(inout) :: self
        if( self%ft ) return
        if( shift_to_phase_origin ) call self%shift_phorig
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        ! now scale the values so that a ifft() of the output yields the
        ! original image back following FFTW
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fwd_ft

    subroutine bwd_ft( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            call fftwf_execute_dft_c2r(self%plan_bwd,self%cmat,self%rmat)
            self%ft = .false.
            if( shift_to_phase_origin ) call self%shift_phorig
        endif
    end subroutine bwd_ft

    subroutine fft_noshift( self )
        class(image), intent(inout) :: self
        if( self%ft ) return
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fft_noshift

    !> \brief ft2img  generates images for visualization of a Fourier transform
    subroutine ft2img( self, which, img )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: which
        class(image),     intent(inout) :: img
        integer :: h,mh,k,mk,l,ml,lims(3,2),inds(3),phys(3)
        integer :: which_flag
        logical :: didft
        complex :: comp
        if( .not.(self.eqdims.img) )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            write(logfhandle,*) 'img%ldim:  ', img%ldim
            THROW_HARD('non-equal dims; ft2img')
        endif
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        select case(which)
            case ('real')
                which_flag = 0
            case('power')
                which_flag = 1
            case('sqrt')
                which_flag = 2
            case ('log')
                which_flag = 3
            case('phase')
                which_flag = 4
            case DEFAULT
                THROW_HARD('unsupported mode: '//trim(which)//'; ft2img')
        end select
        call img%zero_and_unflag_ft
        lims = self%loop_lims(3)
        mh   = abs(lims(1,1))
        mk   = abs(lims(2,1))
        ml   = abs(lims(3,1))
        if( .not.self%wthreads .and. self%is_2d() )then
            do k=lims(2,1),lims(2,2)
                inds(2) = min(max(1,k+mk+1),self%ldim(2))
                do h=lims(1,1),lims(1,2)
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    comp    = self%get_fcomp2D(h,k)
                    select case(which_flag)
                        case(0)
                            img%rmat(inds(1),inds(2),1) = real(comp)
                        case(1)
                            img%rmat(inds(1),inds(2),1) = csq(comp)
                        case(2)
                            img%rmat(inds(1),inds(2),1) = sqrt(csq(comp))
                        case(3)
                            img%rmat(inds(1),inds(2),1) = log(csq(comp))
                        case(4)
                            img%rmat(inds(1),inds(2),1) = phase_angle(comp)
                    end select
                end do
            end do
        else
            !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,comp,inds)&
            !$omp schedule(static) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self%get_fcomp([h,k,l],phys)
                        inds(1) = min(max(1,h+mh+1),self%ldim(1))
                        inds(2) = min(max(1,k+mk+1),self%ldim(2))
                        inds(3) = min(max(1,l+ml+1),self%ldim(3))
                        select case(which_flag)
                        case (0)
                            call img%set(inds,real(comp))
                        case(1)
                            call img%set(inds,csq(comp))
                        case(2)
                            call img%set(inds,sqrt(csq(comp)))
                        case (3)
                            call img%set(inds,log(csq(comp)))
                        case(4)
                            call img%set(inds,phase_angle(comp))
                        end select
                    end do
                end do
            end do
            !$omp end parallel do
        endif
        if( didft ) call self%ifft()
    end subroutine ft2img

    subroutine img2ft( self, img )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: img
        integer :: h,k,l,lims(3,2),logi(3),phys(3)
        integer :: xcnt,ycnt,zcnt
        if( .not.(self.eqdims.img) )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            write(logfhandle,*) 'img%ldim:  ', img%ldim
            THROW_HARD('non-equal dims; img2ft')
        endif
        call img%zero_and_flag_ft
        xcnt = 0
        ycnt = 0
        zcnt = 0
        lims = self%loop_lims(3)
        do h=lims(1,1),lims(1,2)
            xcnt = xcnt + 1
            if( xcnt > self%ldim(1) ) cycle
            ycnt = 0
            do k=lims(2,1),lims(2,2)
                ycnt = ycnt + 1
                if( ycnt > self%ldim(2) ) cycle
                zcnt = 0
                do l=lims(3,1),lims(3,2)
                    zcnt = zcnt + 1
                    if( zcnt > self%ldim(3) ) cycle
                    logi = [h,k,l]
                    phys = self%comp_addr_phys(logi)
                    call img%set_fcomp(logi, phys, cmplx(self%rmat(xcnt,ycnt,zcnt),0.))
                end do
            end do
        end do
    end subroutine img2ft

    !> \brief subtracts the background of an image by subtracting a low-pass filtered
    !!        version of itself
    subroutine subtr_backgr( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
        type(image) :: tmp
        integer     :: winsz
        call tmp%copy(self)
        winsz = nint(real(self%ldim(1)/2)*self%smpd / lp / sqrt(2.))
        call tmp%real_space_filter(winsz, 'average')
        self%rmat = self%rmat - tmp%rmat
        call tmp%kill()
    end subroutine subtr_backgr

    !> \brief generates a real-space resolution mask for matching power-spectra
    subroutine resmsk( self, hplim, lplim )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hplim, lplim
        integer :: h, k, lims(3,2), mh, mk, inds(3)
        real    :: freq, hplim_freq, lplim_freq
        hplim_freq = self%fit%get_find(1,hplim)
        lplim_freq = self%fit%get_find(1,lplim)
        lims = self%loop_lims(3)
        mh = abs(lims(1,1))
        mk = abs(lims(2,1))
        inds = 1
        self%rmat = 0.0
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                freq = hyp(h,k)
                if(freq >= hplim_freq .and. freq <= lplim_freq )then
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    inds(2) = min(max(1,k+mk+1),self%ldim(2))
                    call self%set(inds, 1.)
                endif
            end do
        end do
    end subroutine resmsk

    subroutine frc_pspec( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
        integer     :: k, npix
        type(image) :: maskimg
        logical, allocatable :: l_mask(:,:,:)
        corrs = 0.
        do k = 1, fdim(self1%ldim(1))-3
            call maskimg%ring(self1%ldim, self1%smpd, real(k+2), real(k-2), npix )
            l_mask = bin2logical(maskimg)
            corrs(k) = self1%real_corr(self2, l_mask)
        end do
        call maskimg%kill
        deallocate(l_mask)
    end subroutine frc_pspec

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

    !>  \brief  an image shifter to prepare for Fourier transformation
    subroutine shift_phorig( self )
        class(image), intent(inout) :: self
        complex(sp) :: cswap
        real(sp)    :: rswap
        integer     :: i, j, k, h, ok, ii,jj, kfrom, kto, koffset
        if( .not.self%even_dims() )then
            write(logfhandle,*) 'ldim: ', self%ldim
            THROW_HARD('even dimensions assumed; shift_phorig')
        endif
        if( self%ft )then
            if( self%ldim(3) == 1 )then
                ! serial for now
                koffset = self%array_shape(2) / 2
                do k = 1,koffset
                    ok = k + koffset
                    do h = 1,self%array_shape(1)
                        cswap             = self%cmat(h, k,1)
                        self%cmat(h, k,1) = self%cmat(h,ok,1)
                        self%cmat(h,ok,1) = cswap
                    end do
                end do
            else
                THROW_HARD('3D FT not supported yet; shift_phorig')
            endif
        else
            if( self%ldim(3) == 1 )then
                if( self%wthreads )then
                    !$omp parallel do collapse(2) default(shared) private(rswap,i,j)&
                    !$omp schedule(static) proc_bind(close)
                    do i=1,self%ldim(1)/2
                        do j=1,self%ldim(2)/2
                            !(1)
                            rswap = self%rmat(i,j,1)
                            self%rmat(i,j,1) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,1)
                            self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,1) = rswap
                            !(2)
                            rswap = self%rmat(i,self%ldim(2)/2+j,1)
                            self%rmat(i,self%ldim(2)/2+j,1) = self%rmat(self%ldim(1)/2+i,j,1)
                            self%rmat(self%ldim(1)/2+i,j,1) = rswap
                        end do
                    end do
                    !$omp end parallel do
                else
                    do j=1,self%ldim(2)/2
                        jj = self%ldim(2)/2+j
                        do i=1,self%ldim(1)/2
                            ii = self%ldim(1)/2+i
                            !(1)
                            rswap = self%rmat(i,j,1)
                            self%rmat(i,j,1) = self%rmat(ii,jj,1)
                            self%rmat(ii,jj,1) = rswap
                            !(2)
                            rswap = self%rmat(i,jj,1)
                            self%rmat(i,jj,1) = self%rmat(ii,j,1)
                            self%rmat(ii,j,1) = rswap
                        end do
                    end do
                endif
            else
                kfrom = 1
                kto   = self%ldim(3)/2
                if( self%wthreads )then
                    !$omp parallel do collapse(3) default(shared) private(rswap,i,j,k)&
                    !$omp schedule(static) proc_bind(close)
                    do i=1,self%ldim(1)/2
                        do j=1,self%ldim(2)/2
                            do k=1,kto
                                !(1)
                                rswap = self%rmat(i,j,k)
                                self%rmat(i,j,k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k) = rswap
                                !(2)
                                rswap = self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,j,k)
                                self%rmat(self%ldim(1)/2+i,j,k) = rswap
                                !(3)
                                rswap = self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k) = self%rmat(i,self%ldim(2)/2+j,k)
                                self%rmat(i,self%ldim(2)/2+j,k) = rswap
                                !(4)
                                rswap = self%rmat(i,j,self%ldim(3)/2+k)
                                self%rmat(i,j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k) = rswap
                            end do
                        end do
                    end do
                    !$omp end parallel do
                else
                    do i=1,self%ldim(1)/2
                        do j=1,self%ldim(2)/2
                            do k=1,kto
                                !(1)
                                rswap = self%rmat(i,j,k)
                                self%rmat(i,j,k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k) = rswap
                                !(2)
                                rswap = self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,j,k)
                                self%rmat(self%ldim(1)/2+i,j,k) = rswap
                                !(3)
                                rswap = self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k) = self%rmat(i,self%ldim(2)/2+j,k)
                                self%rmat(i,self%ldim(2)/2+j,k) = rswap
                                !(4)
                                rswap = self%rmat(i,j,self%ldim(3)/2+k)
                                self%rmat(i,j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k) = rswap
                            end do
                        end do
                    end do
                endif
            endif
        endif
    end subroutine shift_phorig

    !> \brief shift  is for origin shifting an image
    !! \param x position in axis 0
    !! \param y position in axis 1
    !! \param z position in axis 2
    subroutine shift( self, shvec )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(3)
        integer :: h, k, l, lims(3,2), phys(3)
        real    :: shvec_here(3)
        logical :: didft
        shvec_here = shvec
        if( self%ldim(3) == 1 ) shvec_here(3) = 0.0
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(phys,h,k,l)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    phys = self%fit%comp_addr_phys([h,k,l])
                    self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
                    self%oshift([h,k,l], shvec_here)
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine shift

    subroutine shift2Dserial_1( self, shvec  )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(2)
        real(dp),allocatable :: hcos(:), hsin(:)
        real(dp) :: sh(2), arg, ck, sk
        integer  :: h,k, hphys,kphys, lims(3,2)
        lims = self%fit%loop_lims(2)
        sh   = real(shvec * self%shconst(1:2),dp)
        allocate(hcos(lims(1,1):lims(1,2)),hsin(lims(1,1):lims(1,2)))
        do h=lims(1,1),lims(1,2)
            arg = real(h,dp)*sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        do k=lims(2,1),lims(2,2)
            kphys = k + 1 + merge(self%ldim(2),0,k<0)
            arg = real(k,dp)*sh(2)
            ck  = dcos(arg)
            sk  = dsin(arg)
            do h=lims(1,1),lims(1,2)
                hphys = h + 1
                self%cmat(hphys,kphys,1) = self%cmat(hphys,kphys,1)&
                    &* cmplx(ck*hcos(h)-sk*hsin(h), ck*hsin(h)+sk*hcos(h),sp)
            end do
        end do
    end subroutine shift2Dserial_1

    subroutine shift2Dserial_2( self, shvec, self_out )
        class(image), intent(inout) :: self, self_out
        real,         intent(in)    :: shvec(2)
        real(dp),allocatable :: hcos(:), hsin(:)
        real(dp) :: sh(2), arg, ck, sk
        integer  :: h,k, hphys,kphys, lims(3,2)
        lims = self%fit%loop_lims(2)
        sh   = real(shvec * self%shconst(1:2),dp)
        allocate(hcos(lims(1,1):lims(1,2)),hsin(lims(1,1):lims(1,2)))
        do h=lims(1,1),lims(1,2)
            arg = real(h,dp)*sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        do k=lims(2,1),lims(2,2)
            kphys = k + 1 + merge(self%ldim(2),0,k<0)
            arg = real(k,dp)*sh(2)
            ck  = dcos(arg)
            sk  = dsin(arg)
            do h=lims(1,1),lims(1,2)
                hphys = h + 1
                self_out%cmat(hphys,kphys,1) = self%cmat(hphys,kphys,1)&
                    &* cmplx(ck*hcos(h)-sk*hsin(h), ck*hsin(h)+sk*hcos(h),sp)
            end do
        end do
    end subroutine shift2Dserial_2

    !> \brief mask  is for spherical masking
    !! \param mskrad mask radius in pixels
    !! \param which mask type
    !! \param inner include cosine edge material
    !! \param width width of inner patch
    subroutine mask( self, mskrad, which, inner, width, backgr )
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        character(len=*), intent(in)    :: which
        real, optional,   intent(in)    :: inner, width, backgr
        real(dp) :: sumv
        real     :: e, wwidth, d_sq, rad_sq, ave
        real     :: cis(self%ldim(1)), cjs(self%ldim(2)), cks(self%ldim(3))
        integer  :: i, j, k, minlen, ir, jr, kr, npix
        logical  :: didft, doinner, soft, avg_backgr
        ! width
        wwidth = 10.
        if( present(width) ) wwidth = width
        ! inner
        doinner = .false.
        if( present(inner) ) doinner = .true.
        ! FT
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ! minlen
        if( self%is_3d() )then
            minlen = minval(self%ldim)
        else
            minlen = minval(self%ldim(1:2))
        endif
        ! soft mask width limited to +/- COSMSKHALFWIDTH pixels
        minlen = min(nint(2.*(mskrad+COSMSKHALFWIDTH)), minlen)
        ! soft/hard
        soft       = .true.
        avg_backgr = .false.
        select case(trim(which))
            case('soft')
                soft  = .true.
                if(present(backgr))then
                    self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
                        &self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - backgr
                else
                    call self%zero_background
                endif
            case('softavg')
                ! background is set to its average value
                if( doinner )THROW_HARD('Inner masking not supported with softavg; simple_image :: mask')
                soft       = .true.
                avg_backgr = .true.
                rad_sq     = mskrad*mskrad
                sumv       = 0.d0
                npix       = 0
            case('hard')
                soft  = .false.
                if(present(backgr))THROW_HARD('no backround subtraction with hard masking; simple_image :: mask')
            case DEFAULT
                THROW_HARD('undefined which parameter; simple_image :: mask')
        end select
        ! init center as origin
        forall(i=1:self%ldim(1)) cis(i) = -real(self%ldim(1))/2. + real(i-1)
        forall(i=1:self%ldim(2)) cjs(i) = -real(self%ldim(2))/2. + real(i-1)
        if(self%is_3d())forall(i=1:self%ldim(3)) cks(i) = -real(self%ldim(3))/2. + real(i-1)
        ! MASKING
        if( soft )then
            ! Soft masking
            if( self%is_3d() )then
                if( avg_backgr )then
                    do k=1,self%ldim(3)
                        do j=1,self%ldim(2)
                            do i=1,self%ldim(1)
                                d_sq = cis(i)*cis(i) + cjs(j)*cjs(j) + cks(k)*cks(k)
                                if( d_sq > rad_sq )then
                                    npix = npix + 1
                                    sumv = sumv + real(self%rmat(i,j,k),dp)
                                endif
                            enddo
                        enddo
                    enddo
                    if( npix > 0 )then
                        ave = real(sumv/real(npix,dp))
                        do i=1,self%ldim(1)
                            do j=1,self%ldim(2)
                                do k=1,self%ldim(3)
                                    e = cosedge(cis(i),cjs(j),cks(k),minlen,mskrad)
                                    if( e < 0.0001 )then
                                        self%rmat(i,j,k) = ave
                                    else
                                        if( e < 0.9999 )then
                                            self%rmat(i,j,k) = e*self%rmat(i,j,k) + (1.-e)*ave
                                        endif
                                    endif
                                enddo
                            enddo
                        enddo
                    endif
                else
                    ! 3d
                    do i=1,self%ldim(1)/2
                        ir = self%ldim(1)+1-i
                        do j=1,self%ldim(2)/2
                            jr = self%ldim(2)+1-j
                            do k=1,self%ldim(3)/2
                                e = cosedge(cis(i),cjs(j),cks(k),minlen,mskrad)
                                if( doinner )e = e * cosedge_inner(cis(i),cjs(j),cks(k),wwidth,inner)
                                if(e > 0.9999) cycle
                                kr = self%ldim(3)+1-k
                                self%rmat(i,j,k)    = e * self%rmat(i,j,k)
                                self%rmat(i,j,kr)   = e * self%rmat(i,j,kr)
                                self%rmat(i,jr,k)   = e * self%rmat(i,jr,k)
                                self%rmat(i,jr,kr)  = e * self%rmat(i,jr,kr)
                                self%rmat(ir,j,k)   = e * self%rmat(ir,j,k)
                                self%rmat(ir,j,kr)  = e * self%rmat(ir,j,kr)
                                self%rmat(ir,jr,k)  = e * self%rmat(ir,jr,k)
                                self%rmat(ir,jr,kr) = e * self%rmat(ir,jr,kr)
                            enddo
                        enddo
                    enddo
                endif
            else
                ! 2d
                if( avg_backgr )then
                    do j=1,self%ldim(2)
                        do i=1,self%ldim(1)
                            d_sq = cis(i)*cis(i) + cjs(j)*cjs(j)
                            if( d_sq > rad_sq )then
                                npix = npix + 1
                                sumv = sumv + real(self%rmat(i,j,1),dp)
                            endif
                        enddo
                    enddo
                    if( npix > 0 )then
                        ave  = real(sumv/real(npix,dp))
                        do i=1,self%ldim(1)
                            do j=1,self%ldim(2)
                                e = cosedge(cis(i),cjs(j),minlen,mskrad)
                                if( e < 0.0001 )then
                                    self%rmat(i,j,1) = ave
                                else
                                    if( e < 0.9999 )then
                                        self%rmat(i,j,1) = e*self%rmat(i,j,1) + (1.-e)*ave
                                    endif
                                endif
                            enddo
                        enddo
                    endif
                else
                    do i=1,self%ldim(1)/2
                        ir = self%ldim(1)+1-i
                        do j=1,self%ldim(2)/2
                            e = cosedge(cis(i),cjs(j),minlen,mskrad)
                            if( doinner )e = e * cosedge_inner(cis(i),cjs(j),wwidth,inner)
                            if(e > 0.9999)cycle
                            jr = self%ldim(2)+1-j
                            self%rmat(i,j,1)   = e * self%rmat(i,j,1)
                            self%rmat(i,jr,1)  = e * self%rmat(i,jr,1)
                            self%rmat(ir,j,1)  = e * self%rmat(ir,j,1)
                            self%rmat(ir,jr,1) = e * self%rmat(ir,jr,1)
                        enddo
                    enddo
                endif
            endif
        else
            ! Hard masking
            if( self%is_3d() )then
                ! 3d
                do i=1,self%ldim(1)/2
                    ir = self%ldim(1)+1-i
                    do j=1,self%ldim(2)/2
                        jr = self%ldim(2)+1-j
                        do k=1,self%ldim(3)/2
                            e = hardedge(cis(i),cjs(j),cks(k),mskrad)
                            if( doinner )e = e * hardedge_inner(cis(i),cjs(j),cks(k),inner)
                            kr = self%ldim(3)+1-k
                            self%rmat(i,j,k)    = e * self%rmat(i,j,k)
                            self%rmat(i,j,kr)   = e * self%rmat(i,j,kr)
                            self%rmat(i,jr,k)   = e * self%rmat(i,jr,k)
                            self%rmat(i,jr,kr)  = e * self%rmat(i,jr,kr)
                            self%rmat(ir,j,k)   = e * self%rmat(ir,j,k)
                            self%rmat(ir,j,kr)  = e * self%rmat(ir,j,kr)
                            self%rmat(ir,jr,k)  = e * self%rmat(ir,jr,k)
                            self%rmat(ir,jr,kr) = e * self%rmat(ir,jr,kr)
                        enddo
                    enddo
                enddo
            else
                ! 2d
                do i=1,self%ldim(1)/2
                    ir = self%ldim(1)+1-i
                    do j=1,self%ldim(2)/2
                        jr = self%ldim(2)+1-j
                        e = hardedge(cis(i),cjs(j),mskrad)
                        if( doinner )e = e * hardedge_inner(ir,jr,inner)
                        self%rmat(i,j,1)   = e * self%rmat(i,j,1)
                        self%rmat(i,jr,1)  = e * self%rmat(i,jr,1)
                        self%rmat(ir,j,1)  = e * self%rmat(ir,j,1)
                        self%rmat(ir,jr,1) = e * self%rmat(ir,jr,1)
                    enddo
                enddo
            endif
        endif
        if( didft ) call self%fft()
    end subroutine mask

    !> \brief neg  is for inverting the contrast
    subroutine neg( self )
        class(image), intent(inout) :: self
        logical :: didft
        didft = .false.
        if( self%ft )then
        else
            call self%fft()
            didft = .true.
        endif
        call self%mul(-1.)
        if( didft ) call self%ifft()
    end subroutine neg

    !> \brief pad is a constructor that pads the input image to input ldim
    !! \param self_in image object
    !! \param self_out image object
    !! \param backgr
    subroutine pad( self_in, self_out, backgr, antialiasing )
        class(image),      intent(inout) :: self_in, self_out
        real,    optional, intent(in)    :: backgr
        logical, optional, intent(in)    :: antialiasing
        real, allocatable :: antialw(:)
        real              :: w, ratio
        integer           :: starts(3), stops(3), lims(3,2)
        integer           :: h, k, l, phys_in(3), phys_out(3)
        logical           :: l_antialiasing
        if( self_in.eqdims.self_out )then
            call self_out%copy(self_in)
            return
        endif
        l_antialiasing = .true.
        if( present(antialiasing) ) l_antialiasing = antialiasing
        if( self_out%ldim(1) >= self_in%ldim(1) .and. self_out%ldim(2) >= self_in%ldim(2)&
        .and. self_out%ldim(3) >= self_in%ldim(3) )then
            if( self_in%ft )then
                self_out = cmplx(0.,0.)
                lims = self_in%fit%loop_lims(2)
                if( l_antialiasing )then
                    antialw = self_in%hannw()
                    !$omp parallel do collapse(3) schedule(static) default(shared)&
                    !$omp private(h,k,l,w,phys_out,phys_in) proc_bind(close)
                    do h=lims(1,1),lims(1,2)
                        do k=lims(2,1),lims(2,2)
                            do l=lims(3,1),lims(3,2)
                                w = antialw(max(1,abs(h)))*antialw(max(1,abs(k)))*antialw(max(1,abs(l)))
                                phys_out = self_out%fit%comp_addr_phys(h,k,l)
                                phys_in  = self_in%fit%comp_addr_phys(h,k,l)
                                self_out%cmat(phys_out(1),phys_out(2),phys_out(3))=&
                                self_in%cmat(phys_in(1),phys_in(2),phys_in(3))*w
                            end do
                        end do
                    end do
                    !$omp end parallel do
                    deallocate(antialw)
                else
                    !$omp parallel do collapse(3) schedule(static) default(shared)&
                    !$omp private(h,k,l,phys_out,phys_in) proc_bind(close)
                    do h=lims(1,1),lims(1,2)
                        do k=lims(2,1),lims(2,2)
                            do l=lims(3,1),lims(3,2)
                                phys_out = self_out%fit%comp_addr_phys(h,k,l)
                                phys_in  = self_in%fit%comp_addr_phys(h,k,l)
                                self_out%cmat(phys_out(1),phys_out(2),phys_out(3))= self_in%cmat(phys_in(1),phys_in(2),phys_in(3))
                            end do
                        end do
                    end do
                    !$omp end parallel do
                endif
                ratio = real(self_in%ldim(1))/real(self_out%ldim(1))
                self_out%smpd = self_in%smpd*ratio ! padding Fourier transform, so sampling is finer
                self_out%ft = .true.
            else
                starts = (self_out%ldim-self_in%ldim)/2+1
                stops  = self_out%ldim-starts+1
                if( self_in%ldim(3) == 1 )then
                    starts(3) = 1
                    stops(3)  = 1
                endif
                if( present(backgr) )then
                    self_out%rmat = backgr
                else
                    self_out%rmat = 0.
                endif
                self_out%rmat(starts(1):stops(1),starts(2):stops(2),starts(3):stops(3)) =&
                self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),:self_in%ldim(3))
                self_out%ft = .false.
            endif
        endif
    end subroutine pad

    subroutine pad_inplace( self, ldim, backgr, antialiasing )
        class(image),   intent(inout) :: self
        integer,           intent(in) :: ldim(3)
        real,    optional, intent(in) :: backgr
        logical, optional, intent(in) :: antialiasing
        type(image) :: tmp
        call tmp%new(ldim, self%smpd, wthreads=self%wthreads)
        call self%pad(tmp, backgr=backgr, antialiasing=antialiasing)
        call self%copy(tmp)
        call tmp%kill()
    end subroutine pad_inplace

    !> \brief pad_mirr is a constructor that pads the input image to input ldim in real space using mirroring
    !! \param self_in image object
    !! \param self_out image object
    subroutine pad_mirr_1( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
        integer :: starts(3), stops(3)
        integer :: i,j, i_in, j_in
        if( self_in.eqdims.self_out )then
            call self_out%copy(self_in)
            return
        endif
        if(self_in%is_3d())THROW_HARD('2D images only; pad_mirr')
        if(self_in%ft)THROW_HARD('real space 2D images only; pad_mirr')
        if( self_out%ldim(1) >= self_in%ldim(1) .and. self_out%ldim(2) >= self_in%ldim(2))then
            self_out%rmat = 0.
            starts  = (self_out%ldim-self_in%ldim)/2+1
            stops   = self_out%ldim-starts+1
            ! actual image
            self_out%rmat(starts(1):stops(1),starts(2):stops(2),1) =&
                &self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1)
            ! left border
            i_in = 0
            do i = starts(1)-1,1,-1
                i_in = i_in + 1
                if(i_in > self_in%ldim(1))exit
                self_out%rmat(i,starts(2):stops(2),1) = self_in%rmat(i_in,:self_in%ldim(2),1)
            enddo
            ! right border
            i_in = self_in%ldim(1)+1
            do i=stops(1)+1,self_out%ldim(1)
                i_in = i_in - 1
                if(i_in < 1)exit
                self_out%rmat(i,starts(2):stops(2),1) = self_in%rmat(i_in,:self_in%ldim(2),1)
            enddo
            ! upper border & corners
            j_in = starts(2)
            do j = starts(2)-1,1,-1
                j_in = j_in + 1
                if(i_in > self_in%ldim(1))exit
                self_out%rmat(:self_out%ldim(1),j,1) = self_out%rmat(:self_out%ldim(1),j_in,1)
            enddo
            ! lower border & corners
            j_in = stops(2)+1
            do j = stops(2)+1, self_out%ldim(2)
                j_in = j_in - 1
                if(j_in < 1)exit
                self_out%rmat(:self_out%ldim(1),j,1) = self_out%rmat(:self_out%ldim(1),j_in,1)
            enddo
            self_out%ft = .false.
        else
            THROW_HARD('inconsistent dimensions; pad_mirr')
        endif
    end subroutine pad_mirr_1

    subroutine pad_mirr_2( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        type(image) :: tmp
        call tmp%new(ldim , self%smpd)
        call self%pad_mirr_1(tmp)
        call self%copy(tmp)
        call tmp%kill
    end subroutine pad_mirr_2

    !> \brief clip is a constructor that clips the input image to input ldim
    !!             DO NOT PARALLELISE
    !! \param self_in image object
    !! \param self_out image object
    subroutine clip( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
        real                        :: ratio
        integer                     :: starts(3), stops(3), lims(3,2)
        integer                     :: phys_out(3), phys_in(3), h,k,l,kpi,kpo,hp
        if( self_out%ldim(1) <= self_in%ldim(1) .and. self_out%ldim(2) <= self_in%ldim(2)&
        .and. self_out%ldim(3) <= self_in%ldim(3) )then
            if( self_in%ft )then
                lims = self_out%fit%loop_lims(2)
                if( self_in%is_2d() )then
                    do k=lims(2,1),lims(2,2)
                        kpi = k + 1 + merge(self_in%ldim(2) ,0,k<0)
                        kpo = k + 1 + merge(self_out%ldim(2),0,k<0)
                        do h=lims(1,1),lims(1,2)
                            hp = h+1
                            self_out%cmat(hp,kpo,1) = self_in%cmat(hp,kpi,1)
                        end do
                    end do
                else
                    do l=lims(3,1),lims(3,2)
                        do k=lims(2,1),lims(2,2)
                            do h=lims(1,1),lims(1,2)
                                phys_out = self_out%fit%comp_addr_phys(h,k,l)
                                phys_in = self_in%fit%comp_addr_phys(h,k,l)
                                self_out%cmat(phys_out(1),phys_out(2),phys_out(3)) =&
                                self_in%cmat(phys_in(1),phys_in(2),phys_in(3))
                            end do
                        end do
                    end do
                endif
                ratio = real(self_in%ldim(1))/real(self_out%ldim(1))
                call self_out%set_smpd(self_in%smpd*ratio) ! clipping Fourier transform, so sampling is coarser
                self_out%ft = .true.
            else
                starts = (self_in%ldim-self_out%ldim)/2+1
                stops  = self_in%ldim-starts+1
                if( self_in%ldim(3) == 1 )then
                    starts(3) = 1
                    stops(3)  = 1
                endif
                self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),:self_out%ldim(3))&
                = self_in%rmat(starts(1):stops(1),starts(2):stops(2),starts(3):stops(3))
                self_out%ft = .false.
            endif
        endif
    end subroutine clip

    subroutine clip_inplace( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        type(image) :: tmp
        call tmp%new(ldim, self%smpd, wthreads=self%wthreads)
        call self%clip(tmp)
        call self%copy(tmp)
        call tmp%kill()
    end subroutine clip_inplace

    subroutine read_and_crop( self, volfname, smpd, box_crop, smpd_crop )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: volfname
        integer,           intent(in)    :: box_crop
        real,              intent(in)    :: smpd, smpd_crop
        integer :: ldim(3), ifoo, box
        call find_ldim_nptcls(volfname, ldim, ifoo)
        ! HE, I would not trust the smpd from the header
        if( ldim(3) /= ldim(1) ) THROW_HARD('Only for volumes')
        box = ldim(1)
        call self%new(ldim, smpd)
        call self%read(volfname)
        if( box < box_crop )then
            ! pad
            call self%fft
            call self%pad_inplace([box_crop,box_crop,box_crop], antialiasing=.false.)
            call self%ifft
        else if( box > box_crop )then
            ! clip
            call self%fft
            call self%clip_inplace([box_crop,box_crop,box_crop])
            call self%ifft
        endif
        call self%set_smpd(smpd_crop) ! safety
    end subroutine read_and_crop

    ! This subroutine rescales the pixel intensities to a new input range.
    subroutine scale_pixels(self, new_range, ssc, oold_range)
          class(image),   intent(inout) :: self
          real,           intent(in)    :: new_range(2)
          real, optional, intent(out)   :: oold_range(2), ssc
          real :: old_range(2), sc
          old_range(1) = minval(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))
          old_range(2) = maxval(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))
          sc = (new_range(2) - new_range(1))/(old_range(2) - old_range(1))
          self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = sc*self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))+new_range(1)-sc*old_range(1)
          if(present(ssc)) ssc = sc
          if(present(oold_range)) oold_range = old_range
    end subroutine scale_pixels

    !>  \brief  is for mirroring an image
    !!          mirror('x') corresponds to mirror2d
    subroutine mirror( self, md, fourier )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: md
        logical, optional, intent(in)    :: fourier
        integer :: i, j
        logical :: didft, l_fourier
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        l_fourier = .false.
        if( present(fourier) ) l_fourier = fourier
        if( md == 'x' )then
            do i=1,self%ldim(2)
                do j=1,self%ldim(3)
                    if( l_fourier )then
                        call reverse_f(self%rmat(1:self%ldim(1),i,j))
                    else
                        call reverse(self%rmat(1:self%ldim(1),i,j))
                    endif
                end do
            end do
        else if( md == 'y' )then
            do i=1,self%ldim(1)
                do j=1,self%ldim(3)
                    if( l_fourier )then
                        call reverse_f(self%rmat(i,1:self%ldim(2),j))
                    else
                        call reverse(self%rmat(i,1:self%ldim(2),j))
                    endif
                end do
            end do
        else if( md == 'z' )then
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    if( l_fourier )then
                        call reverse_f(self%rmat(i,j,1:self%ldim(3)))
                    else
                        call reverse(self%rmat(i,j,1:self%ldim(3)))
                    endif
                end do
            end do
        else
            write(logfhandle,'(a)') 'Mode needs to be either x, y or z; mirror; simple_image'
        endif
        if( didft ) call self%fft()
    end subroutine mirror

    subroutine norm( self, a_s )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: a_s(2)
        integer :: npix
        real    :: ave, var, ep
        npix   = product(self%ldim)
        ave    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        self%rmat = self%rmat - ave
        ep     = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        var    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0)
        var    = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        if( is_a_number(var) )then
            if( var > 0. ) self%rmat = self%rmat / (sqrt(var))
        endif
        if( present(a_s) )then
            self%rmat = self%rmat * a_s(2)
            self%rmat = self%rmat + a_s(1)
        endif
    end subroutine norm

    function variance( self ) result( var )
        class(image), intent(in) :: self
        real    :: ave, var, ep, rmat_subtr_avg(self%ldim(1),self%ldim(2),self%ldim(3))
        integer :: npix
        npix           = product(self%ldim)
        ave            = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        rmat_subtr_avg = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ave
        ep             = sum(rmat_subtr_avg)
        var            = sum(rmat_subtr_avg**2.0)
        var            = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
    end function variance

    subroutine norm_minmax( self  )
        class(image), intent(inout) :: self
        real    :: smin, smax, delta
        smin  = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        smax  = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        delta = smax - smin
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
        &(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - smin)/delta
    end subroutine norm_minmax

    subroutine norm4viz( self, brightness, maxmin)
        class(image),      intent(inout) :: self
        real,    optional, intent(in)    :: brightness
        logical, optional, intent(in)    :: maxmin
        real    :: brightness_l
        logical :: maxmin_l
        brightness_l = 128.0
        maxmin_l     = .false.
        if(present(brightness)) brightness_l = brightness
        if(present(maxmin))     maxmin_l     = maxmin
        if(self%is_ft())THROW_HARD('real space only; norm4viz')
        if(maxmin_l) then
            call self%norm_minmax
            brightness_l = brightness_l - 128.0
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = brightness_l + 256 *&
            &self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
        else
            call self%norm
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = brightness_l + 10.5 *& ! magic numbers from Joe
            &self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
        endif
        ! thresholding
        where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 255. )
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 255.
        else where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0. )
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 0.
        end where
    end subroutine norm4viz

    !> \brief norm_ext  is for normalization of an image using inputted average and standard deviation
    !! \param avg Average
    !! \param sdev Standard deviation
    subroutine norm_ext( self, avg, sdev )
        class(image), intent(inout) :: self
        real, intent(in)            :: avg, sdev
        if( self%ft )then
            THROW_WARN('cannot normalize FTs; norm_ext')
            return
        endif
        if( abs(avg) > TINY ) self%rmat = self%rmat - avg
        if( sdev     > 0.   ) self%rmat = self%rmat / sdev
    end subroutine norm_ext

    !> \brief normalize whole image to standardize image background (outside of logical mask)
    subroutine norm_noise( self, lmsk, sdev_noise )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3)) ! foreground must be true
        real,         intent(inout) :: sdev_noise
        integer :: npix
        real    :: ave, var, ep
        npix = product(self%ldim) - count(lmsk) ! # background pixels
        ave  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.not. lmsk) / real(npix) ! background average
        if( abs(ave) > TINY ) self%rmat = self%rmat - ave
        ep         = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),      mask=.not. lmsk)
        var        = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0, mask=.not. lmsk)
        var        = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        sdev_noise = 0.
        if( is_a_number(var) )then
            sdev_noise = sqrt(var)
            if( var > 0. ) self%rmat = self%rmat / sdev_noise
        endif
    end subroutine norm_noise

    !> \brief normalize whole image to standardize image foreground (within logical mask)
    subroutine norm_within( self, mask )
        class(image),      intent(inout) :: self
        logical,           intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))    ! foreground must be true
        real :: npix, ax, sxx
        npix = real(count(mask))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask) / npix
        self%rmat = self%rmat - ax
        sxx  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask)
        sxx  = sxx/real(npix)
        if( sxx > TINY ) self%rmat = self%rmat / sqrt(sxx)
    end subroutine norm_within

    subroutine calc_bin_thres( self, frac_fg_target, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: frac_fg_target
        real,         intent(out)   :: thres 
        integer, parameter   :: NQUANTA = 10
        real,    allocatable :: pixvals(:), means(:), pix_ts(:), frac_fgs(:)
        integer, allocatable :: labels(:)
        integer :: iq, n_fg, npix, loc(1), cnt
        allocate( means(NQUANTA), labels(NQUANTA), frac_fgs(NQUANTA), pix_ts(NQUANTA) )
        means   = 0.
        labels  = 0
        pixvals = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
        call sortmeans(pixvals, NQUANTA, means, labels)
        npix    = product(self%ldim)
        cnt     = 0
        do iq = 1, NQUANTA
            if( count(labels == iq) == 0 )cycle
            cnt         = cnt + 1
            pix_ts(cnt) = minval(pixvals, mask=labels == iq )
            n_fg        = count(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) >= pix_ts(cnt))
            if( n_fg == npix )then
                frac_fgs(cnt) = 1.
            else if( n_fg == 0 )then
                frac_fgs(cnt) = 0.
            else
                frac_fgs(cnt) = real(n_fg) / real(npix)
            endif
        end do
        loc      = minloc(abs(frac_fgs(:cnt) - frac_fg_target))
        thres    = pix_ts(loc(1))
    end subroutine calc_bin_thres

    subroutine quantize_fwd( self, nquanta, transl_tab, l_msk )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: nquanta
        real,              intent(inout) :: transl_tab(nquanta)
        logical, optional, intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3)) 
        real, allocatable :: pixvals(:)
        integer :: i, j, k, ind
        real    :: dist
        if( present(l_msk) )then
            pixvals = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=l_msk)
        else
            pixvals = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
        endif
        if( self%ldim(3) > 1 )then
            call quantize_vec(pixvals, nquanta, transl_tab)
            !$omp parallel do schedule(static) default(shared) private(k,j,i,ind,dist) collapse(3) proc_bind(close)
            do k = 1,self%ldim(3)
                do j = 1,self%ldim(2)
                    do i = 1,self%ldim(1)
                        call find(transl_tab, nquanta, self%rmat(i,j,k), ind, dist)
                        self%rmat(i,j,k) = real(ind - 1)
                    end do
                end do
            end do
            !$omp end parallel do
        else
            call quantize_vec_serial(pixvals, nquanta, transl_tab)
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    call find(transl_tab, nquanta, self%rmat(i,j,1), ind, dist)
                    self%rmat(i,j,1) = real(ind - 1)
                end do
            end do
        endif
    end subroutine quantize_fwd

    subroutine quantize_bwd( self, nquanta, transl_tab )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nquanta
        real,         intent(in)    :: transl_tab(nquanta)
        real    :: pixvals(nquanta)
        integer :: i, j, k, ind
        real    :: dist
        pixvals = real((/(k,k=1,nquanta)/)) - 1.0
        if( self%ldim(3) > 1 )then
            !$omp parallel do schedule(static) default(shared) private(k,j,i,ind,dist) collapse(3) proc_bind(close)
            do k = 1,self%ldim(3)
                do j = 1,self%ldim(2)
                    do i = 1,self%ldim(1)
                        call find(pixvals, nquanta, self%rmat(i,j,k), ind, dist)
                        self%rmat(i,j,k) = transl_tab(ind)
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    call find(pixvals, nquanta, self%rmat(i,j,1), ind, dist)
                    self%rmat(i,j,1) = transl_tab(ind)
                end do
            end do
        endif
    end subroutine quantize_bwd

    !>  \brief  putting the edge around the image to zero (necessary for avoiding FT artefacts)
    subroutine zero_edgeavg( self )
        class(image), intent(inout) :: self
        real :: edges_sum, edges_ave
        if( self%ft )           THROW_HARD('not for Fted images; zero_edgeavg')
        if( .not.self%is_2d() ) THROW_HARD('only for 2d images; zero_edgeavg')
        edges_sum = sum(self%rmat(1:self%ldim(1),1,1))
        edges_sum = edges_sum + sum(self%rmat(1:self%ldim(1),self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(1,1:self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(self%ldim(1),1:self%ldim(2),1))
        edges_ave = edges_sum / real( 2*(self%ldim(1)+self%ldim(2)) )
        if( abs(edges_ave) > TINY )self%rmat = self%rmat - edges_ave
    end subroutine zero_edgeavg

    !> \brief roavg  is for creating a rotation average of self
    !! \param angstep angular step
    !! \param avg output image rotation average
    subroutine roavg( self, angstep, avg, ang_stop )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: angstep
        class(image),      intent(inout) :: avg
        integer, optional, intent(in)    :: ang_stop
        real(dp):: avgs_rmat(nthr_glob,self%ldim(1),self%ldim(2),1)
        real    :: rotated(nthr_glob,self%ldim(1),self%ldim(2),1)
        integer :: irot, ithr, aang_stop
        aang_stop = 359
        if( present(ang_stop) ) aang_stop = ang_stop
        avgs_rmat = 0._dp
        rotated   = 0.
        !$omp parallel do schedule(static) default(shared) private(irot,ithr) proc_bind(close)
        do irot = 0 + angstep,aang_stop,angstep
            ! get thread index
            ithr = omp_get_thread_num() + 1
            ! rotate & sum
            call self%rtsq_serial(real(irot), 0., 0., rotated(ithr,:,:,:))
            avgs_rmat(ithr,:,:,:) = avgs_rmat(ithr,:,:,:) + real(rotated(ithr,:,:,:), dp)
        end do
        !$omp end parallel do
        ! add in the zero rotation
        avgs_rmat(1,:,:,:) = avgs_rmat(1,:,:,:) + real(self%rmat(:self%ldim(1),:self%ldim(2),:), dp)
        ! normalize and set output image object
        call avg%new(self%ldim, self%smpd)
        call avg%set_rmat(real(sum(avgs_rmat, dim=1)/real(360/angstep,dp)),.false.)
    end subroutine roavg

    !> \brief rtsq  rotation of image by quadratic interpolation (from spider)
    !! \param self_in image object
    !! \param ang angle of rotation
    !! \param shxi shift in x axis
    !! \param shyi shift in y axis
    !! \param self_out optional copy of processed result
    subroutine rtsq( self_in, ang, shxi, shyi, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: ang,shxi,shyi
        class(image), optional, intent(inout) :: self_out
        real    :: shx,shy,ry1,rx1,ry2,rx2,cod,sid,xi,fixcenmshx,fiycenmshy
        real    :: rye2,rye1,rxe2,rxe1,yi,ycod,ysid,yold,xold
        integer :: iycen,ixcen,ix,iy
        real    :: mat_in(self_in%ldim(1),self_in%ldim(2))
        real    :: mat_out(self_in%ldim(1),self_in%ldim(2))
        logical :: didft
        if( self_in%ldim(3) > 1 )         THROW_HARD('only for 2D images; rtsq')
        if( .not. self_in%square_dims() ) THROW_HARD('only for square dims; rtsq;')
        didft = .false.
        if( self_in%ft )then
            call self_in%ifft()
            didft = .true.
        endif
        mat_out = 0. ! this is necessary, because it bugs out if I try to use the 3D matrix
        mat_in = self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1)
        ! shift within image boundary
        shx = amod(shxi,float(self_in%ldim(1)))
        shy = amod(shyi,float(self_in%ldim(2)))
        ! spider image center
        ixcen = self_in%ldim(1)/2+1
        iycen = self_in%ldim(2)/2+1
        ! image dimensions around origin
        rx1 = -self_in%ldim(1)/2
        rx2 =  self_in%ldim(1)/2
        ry1 = -self_in%ldim(2)/2
        ry2 =  self_in%ldim(2)/2
        rye1 = 0.
        rye2 = 0.
        if(mod(self_in%ldim(1),2) == 0)then
            rx2  =  rx2-1.0
            rxe1 = -self_in%ldim(1)
            rxe2 =  self_in%ldim(1)
        else
            rxe1 = -self_in%ldim(1)-1
            rxe2 =  self_in%ldim(1)+1
        endif
        if(mod(self_in%ldim(2),2) == 0)then
            ry2  =  ry2-1.0
            rye1 = -self_in%ldim(2)
            rye2 =  self_in%ldim(2)
        else
            ry2  = -self_in%ldim(2)-1
            rye2 =  self_in%ldim(2)+1
        endif
        ! create transformation matrix
        cod = cos(deg2rad(ang))
        sid = sin(deg2rad(ang))
        !-(center plus shift)
        fixcenmshx = -ixcen-shx
        fiycenmshy = -iycen-shy
        !$omp parallel do default(shared) private(iy,yi,ycod,ysid,ix,xi,xold,yold)&
        !$omp schedule(static) proc_bind(close)
        do iy=1,self_in%ldim(2)
            yi = iy+fiycenmshy
            if(yi < ry1) yi = min(yi+rye2, ry2)
            if(yi > ry2) yi = max(yi+rye1, ry1)
            ycod =  yi*cod+iycen
            ysid = -yi*sid+ixcen
            do ix=1,self_in%ldim(1)
                xi = ix+fixcenmshx
                if(xi < rx1) xi = min(xi+rxe2, rx2)
                if(xi > rx2) xi = max(xi+rxe1, rx1)
                yold = xi*sid+ycod
                xold = xi*cod+ysid
                mat_out(ix,iy) = quadri(xold,yold,mat_in,self_in%ldim(1),self_in%ldim(2))
            enddo
        enddo
        !$omp end parallel do
        if( present(self_out) )then
            if( (self_in.eqdims.self_out) .and. (self_in.eqsmpd.self_out) )then
                ! no need to update dimensions
            else
                call self_out%copy(self_in)
            endif
            self_out%rmat(:self_in%ldim(1),:self_in%ldim(2),1) = mat_out
            self_out%ft = .false.
        else
            self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1) = mat_out
            self_in%ft = .false.
        endif
        if( didft )then
            call self_in%fft()
        endif
    end subroutine rtsq

    !> \brief rtsq  rotation of image by quadratic interpolation (from spider)
    !! \param self_in image object
    !! \param ang angle of rotation
    !! \param shxi shift in x axis
    !! \param shyi shift in y axis
    !! \param rmat_out is the rotated pixels
    subroutine rtsq_serial( self_in, ang, shxi, shyi, rmat_out )
        class(image), intent(inout) :: self_in
        real,         intent(in)    :: ang,shxi,shyi
        real,         intent(inout) :: rmat_out(self_in%ldim(1),self_in%ldim(2),1)
        real    :: shx,shy,ry1,rx1,ry2,rx2,cod,sid,xi,fixcenmshx,fiycenmshy
        real    :: rye2,rye1,rxe2,rxe1,yi,ycod,ysid,yold,xold
        integer :: iycen,ixcen,ix,iy
        real    :: mat_in(self_in%ldim(1),self_in%ldim(2))
        mat_in = self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1)
        ! shift within image boundary
        shx = amod(shxi,float(self_in%ldim(1)))
        shy = amod(shyi,float(self_in%ldim(2)))
        ! spider image center
        ixcen = self_in%ldim(1)/2+1
        iycen = self_in%ldim(2)/2+1
        ! image dimensions around origin
        rx1 = -self_in%ldim(1)/2
        rx2 =  self_in%ldim(1)/2
        ry1 = -self_in%ldim(2)/2
        ry2 =  self_in%ldim(2)/2
        rye1 = 0.
        rye2 = 0.
        if(mod(self_in%ldim(1),2) == 0)then
            rx2  =  rx2-1.0
            rxe1 = -self_in%ldim(1)
            rxe2 =  self_in%ldim(1)
        else
            rxe1 = -self_in%ldim(1)-1
            rxe2 =  self_in%ldim(1)+1
        endif
        if(mod(self_in%ldim(2),2) == 0)then
            ry2  =  ry2-1.0
            rye1 = -self_in%ldim(2)
            rye2 =  self_in%ldim(2)
        else
            ry2  = -self_in%ldim(2)-1
            rye2 =  self_in%ldim(2)+1
        endif
        ! create transformation matrix
        cod = cos(deg2rad(ang))
        sid = sin(deg2rad(ang))
        !-(center plus shift)
        fixcenmshx = -ixcen-shx
        fiycenmshy = -iycen-shy
        do iy=1,self_in%ldim(2)
            yi = iy+fiycenmshy
            if(yi < ry1) yi = min(yi+rye2, ry2)
            if(yi > ry2) yi = max(yi+rye1, ry1)
            ycod =  yi*cod+iycen
            ysid = -yi*sid+ixcen
            do ix=1,self_in%ldim(1)
                xi = ix+fixcenmshx
                if(xi < rx1) xi = min(xi+rxe2, rx2)
                if(xi > rx2) xi = max(xi+rxe1, rx1)
                yold = xi*sid+ycod
                xold = xi*cod+ysid
                rmat_out(ix,iy,1) = quadri(xold,yold,mat_in,self_in%ldim(1),self_in%ldim(2))
            enddo
        enddo
    end subroutine rtsq_serial

    !>  \brief  set pixels to value within a sphere
    subroutine set_within( self, xyz, radius, val )
        class(image), intent(inout) :: self
        real,         intent(in)    :: xyz(3), radius, val
        real    :: rpos(3), vec(3), dist_sq, radius_sq
        integer :: i,j,k, win(3,2)
        radius_sq = radius**2.
        rpos      = xyz / self%smpd
        win(:,1)  = 1 + floor(rpos - radius)
        where( win(:,1) < 1 )win(:,1) = 1
        win(:,2)  = 1 + ceiling(rpos + radius)
        if(win(1,2) > self%ldim(1)) win(1,2) = self%ldim(1)
        if(win(2,2) > self%ldim(2)) win(2,2) = self%ldim(2)
        if(win(3,2) > self%ldim(3)) win(3,2) = self%ldim(3)
        do i = win(1,1),win(1,2)
            do j = win(2,1),win(2,2)
                do k = win(3,1),win(3,2)
                    vec     = real([i,j,k] - 1) * self%smpd - xyz
                    dist_sq = dot_product(vec,vec)
                    if(dist_sq <= radius_sq)self%rmat(i,j,k) = val
                end do
            end do
      end do
    end subroutine set_within

    !>  \brief  cure_outliers for replacing extreme outliers with median of a 13x13 neighbourhood window
    !!          only done on negative values, assuming white ptcls on black bkgr
    !! \param ncured number of corrected points
    !! \param nsigma number of std. dev. to set upper and lower limits
    !! \param deadhot output index of corrected pixels
    !! \param outliers
    subroutine cure_outliers( self, ncured, nsigma, deadhot, outliers )
        class(image),                   intent(inout) :: self
        integer,                        intent(inout) :: ncured
        real,                           intent(in)    :: nsigma
        integer,                        intent(out)   :: deadhot(2)
        logical, allocatable, optional, intent(out)   :: outliers(:,:)
        real, allocatable :: win(:,:), rmat_pad(:,:)
        real    :: ave, sdev, var, lthresh, uthresh
        integer :: i, j, hwinsz, winsz
        logical :: was_fted, err, present_outliers
        if( self%ldim(3)>1 )THROW_HARD('for 2D images only; cure_outliers')
        was_fted = self%is_ft()
        if( was_fted )THROW_HARD('for real space images only; cure_outliers')
        present_outliers = present(outliers)
        ncured   = 0
        hwinsz   = 6
        deadhot  = 0
        if( present_outliers )then
            if( allocated(outliers) ) deallocate(outliers)
            allocate( outliers(self%ldim(1),self%ldim(2)),source =.false. )
        endif
        call moment( self%rmat(1:self%ldim(1),1:self%ldim(2),1), ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh = ave - real(nsigma) * sdev
        uthresh = ave + real(nsigma) * sdev
        if( any(self%rmat(1:self%ldim(1),1:self%ldim(2),1)<lthresh).or.&
            &any(self%rmat(1:self%ldim(1),1:self%ldim(2),1)>uthresh) )then
            winsz = 2*hwinsz+1
            allocate(rmat_pad(1-hwinsz:self%ldim(1)+hwinsz,1-hwinsz:self%ldim(2)+hwinsz), win(winsz,winsz))
            rmat_pad(:,:) = median( reshape(self%rmat(1:self%ldim(1),1:self%ldim(2),1), (/(self%ldim(1)*self%ldim(2))/)) )
            rmat_pad(1:self%ldim(1), 1:self%ldim(2)) = self%rmat(1:self%ldim(1),1:self%ldim(2),1)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,win)&
            !$omp reduction(+:ncured,deadhot) proc_bind(close)
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    if( self%rmat(i,j,1)<lthresh .or. self%rmat(i,j,1)>uthresh )then
                        if( present_outliers )then
                            outliers(i,j)=.true.
                            if (self%rmat(i,j,1)<lthresh) deadhot(1) = deadhot(1) + 1
                            if (self%rmat(i,j,1)>uthresh) deadhot(2) = deadhot(2) + 1
                        else
                            win = rmat_pad( i-hwinsz:i+hwinsz, j-hwinsz:j+hwinsz )
                            self%rmat(i,j,1) = median( reshape(win,(/winsz**2/)) )
                            ncured = ncured + 1
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do
            deallocate( win, rmat_pad )
        endif
    end subroutine cure_outliers

    !>  \brief  zero pixels below threshold
    subroutine zero_below( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        where( self%rmat < thres ) self%rmat = 0.
    end subroutine zero_below

    !>  \brief  divides below threshold
    subroutine div_below( self, thres, val )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres, val
        if( is_equal(val,0.) ) return
        where( self%rmat < thres ) self%rmat = self%rmat / val
    end subroutine div_below

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

    recursive subroutine flip( self, mode )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: mode
        integer :: i,j,ii,jj
        real    :: val
        if( self%ldim(3) > 1 )THROW_HARD('for 2D images only; flip')
        select case(upperCase(trim(mode)))
        case('X')
            !$omp parallel do private(i,j,ii,val) proc_bind(close) default(shared) schedule(static)
            do i = 1,self%ldim(1)/2
                ii = self%ldim(1) - i - 1
                do j = 1,self%ldim(2)
                    val               = self%rmat(i,j,1)
                    self%rmat(i,j,1)  = self%rmat(ii,j,1)
                    self%rmat(ii,j,1) = val
                enddo
            enddo
            !$omp end parallel do
        case('Y')
            !$omp parallel do private(i,j,jj,val) proc_bind(close) default(shared) schedule(static)
            do j = 1,self%ldim(2)/2
                jj = self%ldim(2) - j - 1
                do i = 1,self%ldim(1)
                    val               = self%rmat(i,j,1)
                    self%rmat(i,j,1)  = self%rmat(i,jj,1)
                    self%rmat(i,jj,1) = val
                enddo
            enddo
            !$omp end parallel do
        case('XY','YX')
            call self%flip('X')
            call self%flip('Y')
        case DEFAULT
            THROW_HARD('Unsupported flip mode')
        end select
    end subroutine flip  

    !> \brief rad_cc calculates the radial correlation function between two images/volumes and weight the intensities of the original image/volume
    subroutine radial_cc( self1, self2, self_w, smpd, rad_corrs, rad_dists )
        class(image), intent(inout):: self1, self2, self_w
        real,         intent(in)   :: smpd
        real,         intent(out)  :: rad_corrs(int(self1%ldim(1)/2.)), rad_dists(int(self1%ldim(1)/2.))
        real                 :: rad_weights(int(self1%ldim(1)/2.))
        type(image)          :: distimg
        logical, allocatable :: mask(:,:,:), shell_mask(:,:,:)
        real,    parameter   :: shell_size_pix = 1
        integer :: ldim3, n, n_shells
        real    :: dist_lbound, dist_ubound
        if( .not. (self1.eqdims.self2) ) THROW_HARD('Nonconforming dimensions in image; radial_cc')
        call distimg%new(self1%ldim,smpd)
        n_shells    = int(self1%ldim(1) / 2.)
        if( self1%is_3d() )then
            ! 3D
            ldim3 = self1%ldim(3)
        else
            ! 2D
            ldim3 = 1
        endif
        allocate(mask(self1%ldim(1), self1%ldim(2), ldim3),&
        &  shell_mask(self1%ldim(1), self1%ldim(2), ldim3), source=.true.)
        call distimg%cendist
        do n = 0, n_shells-1
            dist_lbound = real(n) * shell_size_pix
            dist_ubound = dist_lbound + shell_size_pix
            where( (distimg%rmat(:distimg%ldim(1),:distimg%ldim(2),:ldim3) > dist_lbound) .and. &
                  &(distimg%rmat(:distimg%ldim(1),:distimg%ldim(2),:ldim3) < dist_ubound) .and. &
                  &(mask(:distimg%ldim(1),:distimg%ldim(2),:ldim3) ) )
                shell_mask = .true.
            else where
                shell_mask = .false.
            end where
            if( count(shell_mask) < 3 )then
                rad_corrs(n+1) = 0.
            else
                rad_corrs(n+1) = self1%real_corr(self2, shell_mask)
            endif
            rad_dists(n+1) = ( ( dist_lbound * smpd + dist_ubound * smpd ) / 2. )
            if( rad_corrs(n+1)   > 0.      ) rad_weights(n+1) = 2 * rad_corrs(n+1) / (rad_corrs(n+1) +1)
            if( rad_weights(n+1) > 0.99999 ) rad_weights(n+1) = 0.99999
            rad_dists(n+1) = ( ( dist_lbound * smpd + dist_ubound * smpd ) / 2. )
            where( shell_mask(:,:,:) .eqv. .true. )
                self_w%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) = self_w%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) + rad_weights(n+1)
                self1%rmat (:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) = self1%rmat (:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) * rad_weights(n+1)
            end where
        enddo
    end subroutine radial_cc

    subroutine reshape2cube( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(out)   :: self_out
        logical :: isvol
        integer :: ldim(3), ldim_max
        real    :: smpd
        smpd  = self%get_smpd()
        isvol = self%is_3d()
        if(.not.isvol) THROW_HARD('this is only for volumes; reshape2cube')
        ldim  = self%ldim
        if( ldim(1) == ldim(2) .and. ldim(2) == ldim(3) ) return
        ldim_max      = max(ldim(1),ldim(2),ldim(3))
        call self_out%new([ldim_max,ldim_max,ldim_max], self%smpd, wthreads=self%wthreads)
        self_out%rmat = 0.
        self_out%rmat = self%rmat
        self_out%ft   = .false.
        call self_out%set_smpd(smpd)
    end subroutine reshape2cube

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
                call img%write('squares_spider.spi',i)
                call img%write('squares_mrc.mrc',i)
            end do
            ! convert the squares from SPIDER to MRC & vice versa
            do i=1,5
                call img%read('squares_spider.spi',i)
                call img%write('squares_spider_converted.mrc',i)
                call img%read('squares_mrc.mrc',i)
                call img%write('squares_mrc_converted.spi',i)
            end do
            ! test SPIDER vs. MRC & converted vs. nonconverted
            do i=1,20
                call imgs(i)%new(ldim, smpd)
            end do
            cnt = 0
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read('squares_spider.spi',i)
            end do
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read('squares_spider_converted.mrc',i)
            end do
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read('squares_mrc.mrc',i)
            end do
            do i=1,5
                cnt = cnt+1
                call imgs(cnt)%read('squares_mrc_converted.spi',i)
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
                call img%write('cube_spider.spi')
                call img%write('cube_mrc.mrc')
            end do
            ! convert the cubes from SPIDER to MRC & vice versa
            do i=1,5
                call img%read('cube_spider.spi')
                call img%write('cube_spider_converted.mrc')
                call img%read('cube_mrc.mrc')
                call img%write('cube_mrc_converted.spi')
            end do
            ! test SPIDER vs. MRC & converted vs. nonconverted
            do i=1,4
                call imgs(i)%new(ldim, smpd)
                call imgs(i)%read('cube_spider.spi')
                call imgs(i)%read('cube_spider_converted.mrc')
                call imgs(i)%read('cube_mrc.mrc')
                call imgs(i)%read('cube_mrc_converted.spi')
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

    !>  \brief  is a destructor
    subroutine kill( self )
        class(image), intent(inout) :: self
        if( self%existence )then
            call fftwf_free(self%p)
            self%rmat=>null()
            self%cmat=>null()
            !$omp critical
            call fftwf_destroy_plan(self%plan_fwd)
            call fftwf_destroy_plan(self%plan_bwd)
            !$omp end critical
            self%existence = .false.
        endif
    end subroutine kill

    subroutine kill_thread_safe_tmp_imgs( self )
        class(image), intent(in) :: self
        integer :: i
        if( allocated(thread_safe_tmp_imgs) )then
            do i=1,size(thread_safe_tmp_imgs)
                call thread_safe_tmp_imgs(i)%kill
            end do
            deallocate(thread_safe_tmp_imgs)
        endif
    end subroutine kill_thread_safe_tmp_imgs

end module simple_image

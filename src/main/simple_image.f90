! the abstract image data type and its methods. 2D/3D & FT/real all implemented by this class
! and Fourier transformations done in-place to reduce memory usage
module simple_image
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
!!import classes
use simple_ftiter,  only: ftiter
use simple_imgfile, only: imgfile
!!import functions
use simple_winfuns, only: winfuns
use simple_fftw3
use gnufor2
#ifdef PGI
use simple_cufft
#endif

implicit none

public :: image, test_image
#ifdef PGI
public :: test_image_pgi_cuda
#endif

private
#include "simple_local_flags.inc"

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
    procedure          :: construct_thread_safe_tmp_imgs
    procedure          :: disc
    procedure          :: ring
    procedure          :: copy
    procedure          :: mic2spec
    procedure          :: mic2eospecs
    procedure          :: boxconv
    procedure          :: window
    procedure          :: window_slim
    procedure          :: win2arr
    procedure          :: extr_pixels
    procedure          :: corner
    ! I/O
    procedure          :: open
    procedure          :: read
    procedure          :: write
    procedure          :: write_jpg
    ! GETTERS/SETTERS
    procedure          :: get_array_shape
    procedure          :: get_ldim
    procedure          :: get_smpd
    procedure          :: get_nyq
    procedure          :: get_filtsz
    procedure          :: get_shconst
    procedure          :: cyci
    procedure          :: get
    procedure          :: get_rmat
    procedure          :: get_cmat
    procedure, private :: get_cmat_at_1
    procedure, private :: get_cmat_at_2
    generic            :: get_cmat_at => get_cmat_at_1, get_cmat_at_2
    procedure, private :: get_rmat_at_1
    procedure, private :: get_rmat_at_2
    generic            :: get_rmat_at => get_rmat_at_1, get_rmat_at_2
    procedure          :: get_cmatfull
    procedure          :: set
    procedure          :: set_rmat
    procedure, private :: set_cmat_1
    procedure, private :: set_cmat_2
    generic            :: set_cmat => set_cmat_1, set_cmat_2
    procedure          :: set_cmat_at_1
    procedure          :: set_cmat_at_2
    generic            :: set_cmat_at => set_cmat_at_1, set_cmat_at_2
    procedure          :: set_cmats_from_cmats
    procedure          :: add_cmats_to_cmats
    procedure          :: print_cmat
    procedure          :: print_rmat
    procedure          :: expand_ft
    procedure          :: set_ldim
    procedure          :: set_smpd
    procedure          :: get_slice
    procedure          :: set_slice
    procedure          :: get_npix
    procedure          :: get_lfny
    procedure          :: get_lhp
    procedure          :: get_lp
    procedure          :: get_spat_freq
    procedure          :: get_find
    procedure          :: get_clin_lims
    procedure          :: rmat_associated
    procedure          :: cmat_associated
    procedure          :: serialize
    procedure          :: winserialize
    procedure          :: zero2one
    procedure          :: get_fcomp
    procedure          :: set_fcomp
    procedure          :: add_fcomp
    procedure          :: subtr_fcomp
    procedure          :: vis
    procedure          :: set_ft
    procedure          :: extr_fcomp
    procedure          :: packer
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
    ! ARITHMETICS
    procedure, private :: assign
    procedure, private :: assign_r2img
    procedure, private :: assign_c2img
    generic :: assignment(=) => assign, assign_r2img, assign_c2img
    procedure, private :: subtraction
    generic :: operator(-) => subtraction
    procedure, private :: addition
    generic :: operator(+) => addition
    procedure, private :: multiplication
    generic :: operator(*) => multiplication
    procedure, private :: division
    generic :: operator(/) => division
    procedure, private :: l1norm_1
    procedure, private :: l1norm_2
    generic :: operator(.lone.) => l1norm_1
    generic :: operator(.lonesum.) => l1norm_2
    procedure          :: l1weights
    procedure, private :: add_1
    procedure, private :: add_2
    procedure, private :: add_3
    procedure, private :: add_4
    procedure, private :: add_5
    generic            :: add => add_1, add_2, add_3, add_4, add_5
    procedure, private :: add_workshare_1
    procedure, private :: add_workshare_2
    generic            :: add_workshare => add_workshare_1, add_workshare_2
    procedure, private :: subtr_1
    procedure, private :: subtr_2
    procedure, private :: subtr_3
    procedure, private :: subtr_4
    generic            :: subtr => subtr_1, subtr_2, subtr_3, subtr_4
    procedure, private :: div_1
    procedure, private :: div_2
    procedure, private :: div_3
    procedure, private :: div_4
    generic            :: div => div_1, div_2, div_3, div_4
    procedure          :: ctf_dens_correct
    procedure, private :: mul_1
    procedure, private :: mul_2
    procedure, private :: mul_3
    procedure, private :: mul_4
    generic            :: mul => mul_1, mul_2, mul_3, mul_4
    procedure, private :: conjugate
    generic            :: conjg => conjugate
    procedure, private :: mul_rmat_at_1
    procedure, private :: mul_rmat_at_2
    generic            :: mul_rmat_at => mul_rmat_at_1, mul_rmat_at_2
    procedure, private :: div_rmat_at_1
    procedure, private :: div_rmat_at_2
    generic            :: div_rmat_at => div_rmat_at_1, div_rmat_at_2
    procedure, private :: add_cmat_at_1
    procedure, private :: add_cmat_at_2
    generic            :: add_cmat_at => add_cmat_at_1, add_cmat_at_2
    procedure, private :: mul_cmat_at_1
    procedure, private :: mul_cmat_at_2
    procedure, private :: mul_cmat_at_3
    procedure, private :: mul_cmat_at_4
    generic            :: mul_cmat_at => mul_cmat_at_1, mul_cmat_at_2, mul_cmat_at_3, mul_cmat_at_4
    procedure, private :: div_cmat_at_1
    procedure, private :: div_cmat_at_2
    procedure, private :: div_cmat_at_3
    procedure, private :: div_cmat_at_4
    generic            :: div_cmat_at => div_cmat_at_1, div_cmat_at_2, div_cmat_at_3, div_cmat_at_4
    procedure          :: sqpow
    procedure          :: signswap_aimag
    procedure          :: signswap_real
    ! BINARY IMAGE METHODS
    procedure          :: nforeground
    procedure          :: nbackground
    procedure, private :: bin_1
    procedure, private :: bin_2
    generic            :: bin => bin_1, bin_2
    procedure          :: bin_kmeans
    procedure          :: bin_filament
    procedure          :: bin_cylinder
    procedure          :: cendist
    procedure          :: masscen
    procedure          :: center
    procedure          :: center_serial
    procedure          :: bin_inv
    procedure          :: grow_bin
    procedure          :: grow_bins
    procedure          :: shrink_bin
    procedure          :: shrink_bins
    procedure          :: binary_erosion
    procedure          :: binary_dilation
    procedure          :: binary_opening
    procedure          :: binary_closing
    procedure          :: cos_edge
    procedure          :: remove_edge
    procedure          :: increment
    procedure          :: bin2logical
    procedure          :: collage
    ! FILTERS
    procedure          :: acf
    procedure          :: ccf
    procedure          :: guinier_bfac
    procedure          :: guinier
    procedure          :: spectrum
    procedure          :: shellnorm
    procedure, private :: shellnorm_and_apply_filter_1
    procedure, private :: shellnorm_and_apply_filter_2
    generic            :: shellnorm_and_apply_filter =>&
        &shellnorm_and_apply_filter_1, shellnorm_and_apply_filter_2
    procedure, private :: shellnorm_and_apply_filter_serial_1
    procedure, private :: shellnorm_and_apply_filter_serial_2
    generic            :: shellnorm_and_apply_filter_serial =>&
        &shellnorm_and_apply_filter_serial_1, shellnorm_and_apply_filter_serial_2
    procedure          :: apply_bfac
    procedure          :: bp
    procedure          :: gen_lpfilt
    procedure, private :: apply_filter_1
    procedure, private :: apply_filter_2
    generic            :: apply_filter => apply_filter_1, apply_filter_2
    procedure          :: apply_filter_serial
    procedure          :: phase_rand
    procedure          :: hannw
    procedure          :: real_space_filter
    ! CALCULATORS
    procedure          :: square_root
    procedure          :: maxcoord
    procedure          :: ccpeak_offset
    procedure          :: minmax
    procedure          :: rmsd
    procedure, private :: stats_1
    procedure, private :: stats_2
    generic            :: stats => stats_1, stats_2
    procedure          :: noisesdev
    procedure          :: mean
    procedure          :: contains_nans
    procedure          :: checkimg4nans
    procedure          :: cure
    procedure          :: loop_lims
    procedure          :: comp_addr_phys
    procedure          :: get_2Dphys_ind_mapping
    procedure          :: corr
    procedure          :: corr_shifted
    procedure, private :: real_corr_1
    procedure, private :: real_corr_2
    generic            :: real_corr => real_corr_1, real_corr_2
    procedure          :: prenorm4real_corr_1
    procedure          :: prenorm4real_corr_2
    generic            :: prenorm4real_corr => prenorm4real_corr_1, prenorm4real_corr_2
    procedure, private :: real_corr_prenorm_1
    procedure, private :: real_corr_prenorm_2
    generic            :: real_corr_prenorm => real_corr_prenorm_1, real_corr_prenorm_2
    procedure          :: fsc
    procedure          :: get_nvoxshell
    procedure          :: get_res
    procedure, private :: oshift_1
    procedure, private :: oshift_2
    generic            :: oshift => oshift_1, oshift_2
    procedure, private :: gen_argtransf_comp
    procedure          :: gen_argtransf_mats
    ! MODIFIERS
    procedure          :: insert
    procedure          :: insert_lowres
    procedure          :: insert_lowres_serial
    procedure          :: inv
    procedure          :: ran
    procedure          :: gauran
    procedure          :: add_gauran
    procedure          :: dead_hot_positions
    procedure          :: taper_edges
    procedure          :: zero_and_unflag_ft
    procedure          :: zero_and_flag_ft
    procedure          :: zero_background
    procedure          :: subtr_backgr_pad_divwinstr_fft
    procedure          :: subtr_backgr_pad_divwinstr_fft_serial
    procedure          :: subtr_backgr_pad_divwinstr_serial
    procedure          :: salt_n_pepper
    procedure          :: square
    procedure          :: corners
    procedure          :: before_after
    procedure          :: gauimg
    procedure          :: gauimg2
    procedure          :: gauimg2D
    procedure          :: dampen_central_cross
    procedure          :: subtr_backgr
    procedure          :: subtr_avg_and_square
    procedure          :: resmsk
    procedure          :: frc_pspec
    procedure          :: mask
    procedure          :: neg
    procedure          :: pad
    procedure          :: pad_mirr
    procedure          :: clip
    procedure          :: clip_inplace
    procedure          :: mirror
    procedure          :: norm
    procedure          :: norm_ext
    procedure          :: radius_norm
    procedure          :: noise_norm
    procedure          :: edges_norm
    procedure          :: norm_bin
    procedure          :: roavg
    procedure          :: rtsq
    procedure          :: rtsq_serial
    procedure          :: shift_phorig
    procedure          :: shift
    procedure          :: shift2Dserial
    procedure          :: set_within
    procedure          :: ft2img
    procedure          :: img2ft
    ! FFTs and signal processing operations
    procedure          :: fwd_ft    ! FFTW forward transform
    procedure          :: bwd_ft    ! FFTW inverse transform
    procedure          :: fft_pgi_cuda  ! CUFFT forward transform
    !        procedure         :: fft_pgi_stream
    procedure          :: ifft_pgi_cuda ! CUFFT inverse transform
#ifdef PGI
    generic            :: fft => fft_pgi_cuda
    generic            :: ifft => ifft_pgi_cuda
#else
    generic            :: fft => fwd_ft
    generic            :: ifft => bwd_ft
#endif
    ! DENOISING FUNCTIONS
    procedure          :: cure_outliers
    procedure          :: denoise_NLM
    procedure          :: zero_below
    ! DESTRUCTOR
    procedure :: kill
end type image

interface image
    module procedure constructor
end interface image

! CLASS PARAMETERS/VARIABLES
logical,     parameter   :: shift_to_phase_origin=.true.
type(image), allocatable :: thread_safe_tmp_imgs(:)

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    !! \param ldim image dimesions
    !! \param smpd sampling distance
    !! \param backgr  constant initial background
    !! \return  self new image object
    !!
    function constructor( ldim, smpd, wthreads ) result( self ) !(FAILS W PRESENT GFORTRAN)
        integer,           intent(in) :: ldim(:)
        real,              intent(in) :: smpd
        logical, optional, intent(in) :: wthreads
        type(image) :: self
        call self%new( ldim, smpd, wthreads )
    end function constructor

    !>  \brief  Constructor for simple_image class
    !!
    !!\param self this image object
    !!\param ldim 3D dimensions
    !!\param smpd sampling distance
    !!\param backgr constant initial background
    !!
    !!\return new image object
    subroutine new( self, ldim, smpd, wthreads )
        !! have to have a type-bound constructor here because we get a sigbus error with the function construct
        !! "program received signal sigbus: access to an undefined portion of a memory object."
        !! this seems to be related to how the cstyle-allocated matrix is referenced by the gfortran compiler
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
                call self%kill
            else
                do_allocate = .false.
                call fftwf_destroy_plan(self%plan_fwd)
                call fftwf_destroy_plan(self%plan_bwd)
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

    subroutine construct_thread_safe_tmp_imgs( self, nthr )
        class(image), intent(in) :: self
        integer,      intent(in) :: nthr
        integer :: i, ldim(3), sz
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
                call thread_safe_tmp_imgs(i)%new(self%ldim, self%smpd, wthreads=.false.)
            end do
        endif
    end subroutine construct_thread_safe_tmp_imgs

    !>  \brief disc constructs a binary disc of given radius and returns the number of 1:s
    !>
    !! \param ldim  image dimensions
    !! \param smpd sampling distance
    !! \param radius  radius of disc
    !! \param npix  num of ON bits in mask
    !!
    subroutine disc( self, ldim, smpd, radius, npix )
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
        if( present(npix) )npix = count(self%rmat>0.5)
    end subroutine disc

    !>  \brief disc constructs a binary ring and returns the number of 1:s
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

    !>  \brief copy is a constructor that copies the input object
    !! \param self image object
    !! \param self_in rhs object
    !!
    subroutine copy( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
        call self%new(self_in%ldim, self_in%smpd, wthreads=self%wthreads)
        self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
            &self_in%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        self%ft = self_in%ft
    end subroutine copy

    !> mic2spec calculates the average powerspectrum over a micrograph
    !! \param self image object
    !! \param box boxwidth filter size
    !! \param speckind
    !! \return img_out processed image
    !!
    function mic2spec( self, box, speckind ) result( img_out )
        class(image),     intent(inout) :: self
        integer,          intent(in)    :: box
        character(len=*), intent(in)    :: speckind
        type(image) :: img_out, tmp, tmp2
        integer     :: xind, yind, cnt
        logical     :: didft
        if( self%ldim(3) /= 1 ) call simple_stop('only for 2D images; mic2spec; simple_image')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            call simple_stop('cannot use a box larger than the image; mic2spec; simple_image')
        endif
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
                call self%window([xind,yind],box,tmp)
                call tmp%norm()
                call tmp%edges_norm
                call tmp%fft()
                call tmp%ft2img(speckind, tmp2)
                call img_out%add(tmp2)
                cnt = cnt+1
                call tmp%zero_and_unflag_ft
                call tmp2%zero_and_unflag_ft
            end do
        end do
        call img_out%div(real(cnt))
        if( didft ) call self%fft()
    end function mic2spec

    subroutine mic2eospecs( self, box, speckind, pspec_lower, pspec_upper, pspec_all )
        class(image),     intent(inout) :: self
        integer,          intent(in)    :: box
        character(len=*), intent(in)    :: speckind
        class(image),     intent(inout) :: pspec_lower, pspec_upper, pspec_all
        type(image) :: tmp, tmp2
        integer     :: xind, yind, cnt, neven_lim, nodd_lim, neven, nodd
        logical     :: didft, outside
        if( self%ldim(3) /= 1 ) stop 'only for 2D images; mic2eoimgs; simple_image'
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            stop 'cannot use a box larger than the image; mic2eoimgs; simple_image'
        endif
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ! count # images
        cnt = 0
        do xind=0,self%ldim(1)-box,box/2
            do yind=0,self%ldim(2)-box,box/2
                cnt = cnt + 1
            end do
        end do
        ! divide eo
        neven_lim = cnt / 2
        nodd_lim  = cnt - neven_lim
        ! prepate lower/upper/all spectra
        cnt   = 0
        neven = 0
        nodd  = 0
        call tmp%new([box,box,1], self%smpd)
        call tmp2%new([box,box,1], self%smpd)
        do xind=0,self%ldim(1)-box,box/2
            do yind=0,self%ldim(2)-box,box/2
                cnt = cnt + 1
                call self%window_slim([xind,yind],box,tmp,outside)
                call tmp%norm()
                call tmp%edges_norm
                call tmp%fft()
                if( cnt <= neven_lim )then
                    neven = neven + 1
                    call tmp%ft2img(speckind, tmp2)
                    call pspec_lower%add(tmp2)
                else
                    nodd = nodd + 1
                    call tmp%ft2img(speckind, tmp2)
                    call pspec_upper%add(tmp2)
                endif
                call tmp%zero_and_unflag_ft
                call tmp2%zero_and_unflag_ft
            end do
        end do
        call pspec_all%add(pspec_lower)
        call pspec_all%add(pspec_upper)
        call tmp%kill
        call tmp2%kill
        if( didft ) call self%fft()
    end subroutine mic2eospecs

    function boxconv( self, box ) result( img_out )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: box
        type(image) :: img_out, tmp
        integer     :: xind, yind, cnt
        logical     :: didft
        if( self%ldim(3) /= 1 ) call simple_stop('only for 2D images; boxconvolute; simple_image')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            call simple_stop('cannot boxconvolute using a box larger than the image; boxconvolute; simple_image')
        endif
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        call img_out%new([box,box,1], self%smpd)
        cnt = 0
        do xind=0,self%ldim(1)-box,box/2
            do yind=0,self%ldim(2)-box,box/2
                call self%window([xind,yind],box,tmp)
                call img_out%add(tmp)
                cnt = cnt+1
                call tmp%kill()
            end do
        end do
        call img_out%div(real(cnt))
        if( didft ) call self%fft()
    end function boxconv

    !>  \brief window extracts a particle image from a box as defined by EMAN 1.9
    !! \param self_in input image object
    !! \param coord window coordinates
    !! \param box boxwidth filter size
    !! \param self_out return image object
    !! \param noutside num window pixels outside image
    !!
    subroutine window( self_in, coord, box, self_out, noutside )
        class(image),      intent(in)    :: self_in
        integer,           intent(in)    :: coord(2), box
        class(image),      intent(inout) :: self_out
        integer, optional, intent(inout) :: noutside
        integer :: fromc(2), toc(2), xoshoot, yoshoot, xushoot, yushoot, xboxrange(2), yboxrange(2)
        if( self_in%ldim(3) > 1 ) call simple_stop('only 4 2D images; window; simple_image')
        if( self_in%is_ft() )     call simple_stop('only 4 real images; window; simple_image')
        if( self_out%exists() )then
            if( self_out%is_ft() ) call simple_stop('only 4 real images; window; simple_image')
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
                write(*,*) 'WARNING! Box extends outside micrograph; window; simple_image'
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
    end subroutine window

    !>  window_slim  extracts a particle image from a box as defined by EMAN 1.9
    !! \param self_in image object
    !! \param coord x,y coordinates
    !! \param box  boxwidth filter size
    !! \param self_out output image object
    !! \param noutside  num window pixels outside image
    !!
    subroutine window_slim( self_in, coord, box, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: coord(2), box !< boxwidth filter size
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
        integer :: fromc(2), toc(2)
        fromc = coord + 1         ! compensate for the c-range that starts at 0
        toc   = fromc + (box - 1) ! the lower left corner is 1,1
        self_out%rmat = 0.
        outside = .false.
        if( fromc(1) < 1 .or. fromc(2) < 1 .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) )then
            outside = .true.
        else
            self_out%rmat(1:box,1:box,1) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),1)
        endif
    end subroutine window_slim

    !>  \brief win2arr extracts a small window into an array (circular indexing)
    !! \param i,j,k window coords
    !! \param winsz window half-width size (odd)
    !! \return  pixels index array to pixels in window
    !!
    function win2arr( self, i, j, k, winsz ) result( pixels )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i, j, k, winsz
        real, allocatable :: pixels(:)
        integer :: s, ss, t, tt, u, uu, cnt, npix
        if( self%is_ft() ) call simple_stop('only 4 real images; win2arr; simple_image')
        if( self%is_3d() )then
            npix = (2*winsz+1)**3
        else
            npix = (2*winsz+1)**2
        endif
        allocate(pixels(npix), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: win2arr; simple_image',alloc_stat)
        cnt = 1
        do s=i-winsz,i+winsz
            ss = cyci_1d([1,self%ldim(1)], s)
            do t=j-winsz,j+winsz
                tt = cyci_1d([1,self%ldim(2)], t)
                if( self%ldim(3) > 1 )then
                    do u=k-winsz,k+winsz
                        uu          = cyci_1d([1,self%ldim(3)], u)
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

    !>  \brief extr_pixels extracts the pixels under the mask
    !! \param mskimg
    !! \return  pixels index array to pixels in mask
    !!
    function extr_pixels( self, mskimg ) result( pixels )
        class(image), intent(in) :: self
        class(image), intent(in) :: mskimg
        real, allocatable :: pixels(:)   !< 1D pixel array: self(mskimg==1)
        if( self%is_ft() ) call simple_stop('only 4 real images; extr_pixels; simple_image')
        if( self.eqdims.mskimg )then
            ! pixels = self%packer(mskimg) ! Intel hickup
            !$ allocate(pixels(self%ldim(1)*self%ldim(2)*self%ldim(3)), stat=alloc_stat)
            !$ if(alloc_stat.ne.0)call allocchk("In simple_image::extr_pixels PGI prealloc ",alloc_stat)
            pixels = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                &mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))>0.5 )
        else
            call simple_stop('mask and image of different dims; extr_pixels; simple_image')
        endif
    end function extr_pixels

    !>  \brief corner extracts a corner of a volume with size box
    !! \param self_in input image volume
    !! \param box  size of box
    !! \param self_out extracted corner of a volume
    !!
    subroutine corner( self_in, box, self_out )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: box
        type(image),  intent(inout) :: self_out
        if( self_in%ldim(3) <= 1 ) call simple_stop('only 4 3D images; corner; simple_image')
        if( self_in%is_ft() )      call simple_stop('only 4 real images; corner; simple_image')
        call self_out%new([box,box,box], self_in%smpd)
        self_out%rmat(:box,:box,:box) = self_in%rmat(:box,:box,:box)
    end subroutine corner

    ! I/O

    !>  \brief  open: for reading 2D images from stack or volumes from volume files
    !! \param fname   filename of image
    !! \param ioimg   IO file object
    !! \param formatchar  image type format (M,F,S)
    !! \param readhead  get header flag
    !! \param rwaction  read/write flag
    !!
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
                print *, 'file: ', trim(fname)
                call simple_stop('The file you are trying to open does not exists; open; simple_image')
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
                ! data type: 0 image: signed 8-bit bytes rante -128 to 127
                !            1 image: 16-bit halfwords
                !            2 image: 32-bit reals (DEFAULT MODE)
                !            3 transform: complex 16-bit integers
                !            4 transform: complex 32-bit reals (THIS WOULD BE THE DEFAULT FT MODE)
                mode = ioimg%getMode()
                if( mode == 3 .or. mode == 4 ) self%ft = .true.
            case('F','S')
                call ioimg%open(fname, self%ldim, self%smpd, formatchar=formatchar, readhead=readhead, rwaction=rwaction)
            end select
        else
            call simple_stop('ERROR, image need to be constructed before read/write; open; simple_image')
        endif
    end subroutine open

    !>  \brief read: for reading 2D images from stack or volumes from volume files
    !! \param fname            filename of image
    !! \param i                file index in stack
    !! \param ioimg            image IO object
    !! \param formatchar       image type (M,F,S)
    !! \param readhead         get header info flag
    !! \param rwaction         read mode flag
    !!
    subroutine read( self, fname, i, formatchar, readhead )
        class(image),               intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        integer,          optional, intent(in)    :: i
        character(len=1), optional, intent(in)    :: formatchar
        logical,          optional, intent(in)    :: readhead
        type(imgfile)         :: ioimg
        character(len=1)      :: form
        integer               :: ldim(3), iform, first_slice
        integer               :: last_slice, ii
        real                  :: smpd
        logical               :: isvol
        ldim          = self%ldim
        smpd          = self%smpd
        isvol = .true. ! assume volume by default
        ii    = 1      ! default location
        if( present(i) )then
            ! we are reading from a stack & in SIMPLE volumes are not allowed
            ! to be stacked so the image object must be 2D
            isvol = .false.
            ii = i ! replace default location
        endif
        if( present(formatchar) )then
            form = formatchar
        else
            form = fname2format(fname)
        endif
        select case(form)
        case('M', 'F', 'S')
            call self%open(fname, ioimg, formatchar, readhead, rwaction='READ')
        case DEFAULT
            write(*,*) 'Trying to read from file: ', trim(fname)
            stop 'ERROR, unsupported file format; read; simple_image'
        end select
        call exception_handler(ioimg)
        call read_local(ioimg)

    contains

        !> read_local
        !! \param ioimg Image file object
        !!
        subroutine read_local( ioimg )
            class(imgfile) :: ioimg
            ! work out the slice range
            if( isvol )then
                if( ii .gt. 1 ) stop 'ERROR, stacks of volumes not supported; read; simple_image'
                first_slice = 1
                last_slice = ldim(3)
            else
                first_slice = ii
                last_slice = ii
            endif
            call ioimg%rSlices(first_slice,last_slice,self%rmat)
            call ioimg%close
        end subroutine read_local

        !> exception_handler
        !! \param ioimg Image file object
        !!
        subroutine exception_handler( ioimg )
            class(imgfile) :: ioimg
            if( form .eq. 'S' ) call spider_exception_handler(ioimg)
            if( form .ne. 'F' )then
                ! make sure that the logical image dimensions of self are consistent with the overall header
                ldim = ioimg%getDims()
                if( .not. all(ldim(1:2) == self%ldim(1:2)) )then
                    write(*,*) 'ldim of image object: ', self%ldim
                    write(*,*) 'ldim in ioimg (fhandle) object: ', ldim
                    call simple_stop('ERROR, logical dims of overall header & image object do not match; read; simple_image')
                endif
            endif
        end subroutine exception_handler

        !> spider_exception_handler
        !! \param ioimg Image IO object to get Iform
        !! iform file type specifier:
        !!   1 = 2D image
        !!   3 = 3D volume
        !! -11 = 2D Fourier odd
        !! -12 = 2D Fourier even
        !! -21 = 3D Fourier odd
        !! -22 = 3D Fourier even
        subroutine spider_exception_handler(ioimg)
            class(imgfile) :: ioimg
            iform = ioimg%getIform()
            select case(iform)
            case(1,-11,-12)
                ! we are processing a stack of 2D images (single 2D images not allowed in SIMPLE)
                if( present(i) )then
                    ! all good
                else
                    call simple_stop('ERROR, optional argument i required for reading from stack; read; simple_image')
                endif
                if( self%ldim(3) == 1 )then
                    ! all good
                else if( self%ldim(3) > 1 )then
                    call simple_stop('ERROR, trying to read from a stack into a volume; read; simple_image')
                else
                    call simple_stop('ERROR, nonconforming logical dimension of image; read; simple_image')
                endif
                if( iform == -11 .or. iform == -12 ) self%ft = .true.
            case(3,-21,-22)
                ! we are processing a 3D image (stacks of 3D volumes not allowed in SIMPLE)
                if( present(i) )then
                    call simple_stop('ERROR, stacks of 3D volumes not allowed in SIMPLE; read; simple_image')
                endif
                if( self%ldim (3) > 1 )then
                    ! all good
                else if( self%ldim(3) == 1)then
                    call simple_stop('ERROR, trying to read from a volume into a 2D image; read; simple_image')
                else
                    call simple_stop('ERROR, nonconforming logical dimension of image; read; simple_image')
                endif
                if( iform == -21 .or. iform == -22 ) self%ft = .true.
            case DEFAULT
                write(*,*) 'iform = ', iform
                call simple_stop('Unsupported iform flag; simple_image :: read')
            end select
        end subroutine spider_exception_handler

    end subroutine read

    !>  \brief  for writing any kind of images to stack or volumes to volume files
    !! \param fname filename of image
    !! \param i  file index in stack/part
    !! \param del_if_exists overwrite if present
    !! \param formatchar
    !! \todo fix optional args
    subroutine write( self, fname, i, del_if_exists, formatchar )
        class(image),               intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        integer,          optional, intent(in)    :: i
        logical,          optional, intent(in)    :: del_if_exists
        character(len=1), optional, intent(in)    :: formatchar
        real             :: dev
        type(imgfile)    :: ioimg
        character(len=1) :: form
        integer          :: first_slice, last_slice, iform, ii
        logical          :: isvol, die
        if( self%existence )then
            dev = self%rmsd()
            die = .false.
            if( present(del_if_exists) ) die = del_if_exists
            if( self%is_2d() )then
                isvol = .false.
            else if( self%is_3d() )then
                isvol = .true.
            endif
            ii = 1 ! default location
            if( present(i) )then
                ! we are writing to a stack & in SIMPLE volumes are not allowed
                ! to be stacked so the image object must be 2D
                if( isvol )then
                    call simple_stop('ERROR, trying to write 3D image to stack ; write; simple_image')
                endif
                ii = i ! replace default location
            endif
            ! find format
            if( present(formatchar) )then
                form = formatchar
            else
                form = fname2format(fname)
            endif
            select case(form)
            case('M','F')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=die,&
                    formatchar=formatchar, readhead=.false.)
                ! data type: 0 image: signed 8-bit bytes rante -128 to 127
                !            1 image: 16-bit halfwords
                !            2 image: 32-bit reals (DEFAULT MODE)
                !            3 transform: complex 16-bit integers
                !            4 transform: complex 32-bit reals (THIS WOULD BE THE DEFAULT FT MODE)
                if( self%ft )then
                    call ioimg%setMode(4)
                else
                    call ioimg%setMode(2)
                endif
                call ioimg%setRMSD(dev)
            case('S')
                ! pixel size of object overrides pixel size in header
                call ioimg%open(fname, self%ldim, self%smpd, del_if_exists=die,&
                    formatchar=formatchar, readhead=.false.)
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
            case DEFAULT
                write(*,*) 'format descriptor: ', form
                call simple_stop('ERROR, unsupported file format; write; simple_image')
            end select
            ! work out the slice range
            if( isvol )then
                if( ii .gt. 1 ) call simple_stop('ERROR, stacks of volumes not supported; write; simple_image')
                first_slice = 1
                last_slice = self%ldim(3)
            else
                first_slice = ii
                last_slice = ii
            endif
            ! write slice(s) to disk
            call ioimg%wSlices(first_slice,last_slice,self%rmat,self%ldim,self%ft,self%smpd)
            call ioimg%close
        else
            call simple_stop('ERROR, nonexisting image cannot be written to disk; write; simple_image')
        endif
    end subroutine write

    subroutine write_jpg( self, fname, quality, colorspec, norm )
        use simple_jpg, only: jpg_img
        class(image),               intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        integer,          optional, intent(in)    :: quality, colorspec
        logical,          optional, intent(in)    :: norm
        type(jpg_img)     :: jpg
        real, allocatable :: pixels(:)
        real              :: ave, sdev, var
        integer           :: status
        logical           :: norm_here, err
        norm_here = .true.
        if(present(norm))norm_here = norm
        if( norm_here )then
            pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
            call moment(pixels, ave, sdev, var, err)
            deallocate(pixels)
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 128. + (self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ave) * (10. / sdev)
            where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 255. )
                self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 255.
            else where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0. )
                self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 0.
            end where
        endif
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
    !! \return  shape array dimensions
    !!
    pure function get_array_shape( self ) result( shape)
        class(image), intent(in) :: self
        integer :: shape(3)
        shape = self%array_shape
    end function get_array_shape

    !> \brief get_ldim  is a getter
    !! \return  ldim
    !!
    pure function get_ldim( self ) result( ldim )
        class(image), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim

    !> \brief get_smpd  is a getter
    !! \return smpd
    !!
    pure function get_smpd( self ) result( smpd )
        class(image), intent(in) :: self
        real :: smpd
        smpd = self%smpd
    end function get_smpd

    !>  \brief get_nyq get the Nyquist Fourier index
    !! \return nyq Nyquist Fourier index
    pure function get_nyq( self ) result( nyq )
        class(image), intent(in) :: self
        integer :: nyq
        nyq = fdim(self%ldim(1)) - 1
    end function get_nyq

    !> \brief get_filtsz  to get the size of the filters
    !! \return  n size of the filter
    !!
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

    !> \brief cyci  cyclic index generation
    !! \param logi
    !! \return  inds
    !!
    function cyci( self, logi ) result( inds )
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3)
        integer                  :: inds(3), lims(3,2)
        if( self%is_ft() )then
            lims = self%loop_lims(3)
        else
            lims = 1
            lims(1,2) = self%ldim(1)
            lims(2,2) = self%ldim(2)
            lims(3,2) = self%ldim(3)
        endif
        inds(1) = cyci_1d(lims(1,:), logi(1))
        inds(2) = cyci_1d(lims(2,:), logi(2))
        inds(3) = cyci_1d(lims(3,:), logi(3))
    end function cyci

    !> \brief get  is a getter
    !! \param logi
    !! \return  val
    !!
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

    !> \brief get_rmat  is a getter
    !! \return  rmat
    !!
    function get_rmat( self ) result( rmat )
        class(image), intent(in) :: self
        real, allocatable :: rmat(:,:,:)
        integer :: ldim(3)
        ldim = self%ldim
        allocate(rmat(ldim(1),ldim(2),ldim(3)), source=self%rmat(:ldim(1),:ldim(2),:ldim(3)), stat=alloc_stat)
        if (alloc_stat /= 0)call allocchk("simple_image::get_rmat ",alloc_stat)
    end function get_rmat

    !! get rmat value at index logi
    function get_rmat_at_1( self, logi ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3)
        real :: val
        val = self%rmat(logi(1),logi(2),logi(3))
    end function get_rmat_at_1

    !! get cmat value at index phys
    elemental pure function get_rmat_at_2( self, i,j,k ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: i,j,k
        real :: val
        val = self%rmat(i,j,k)
    end function get_rmat_at_2

    !>  \brief   get_cmat get the image object's complex matrix
    !! \return cmat a copy of this image object's cmat
    !!
    pure function get_cmat( self ) result( cmat )
        class(image), intent(in) :: self
        !integer :: array_shape(3)
        complex, allocatable :: cmat(:,:,:)
        allocate(cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3)), source=self%cmat)
        !! pure function cannot call allocchk
        !! if (alloc_stat /= 0)call allocchk("simple_image::get_cmat ")
    end function get_cmat

    !>  \brief   get_cmatfull get the image object's complex matrix, including index limits
    !! \return cmat a copy of this image object's cmat
    !!
    function get_cmatfull( self , llims) result( cmat )
        class(image), intent(in) :: self
        integer, intent(inout), optional :: llims(3,2)
        complex, allocatable :: cmat(:,:,:)
        integer :: llimits(3,2)
        if (present(llims))then
            llimits=llims
        else
            llimits=1
            llimits(:,2) = self%array_shape
        end if
        allocate(cmat(llimits(1,1):llimits(1,2),&
            &llimits(2,1):llimits(2,2),&
            &llimits(3,1):llimits(3,2)), source=self%cmat, stat=alloc_stat)
        if (alloc_stat /= 0)call allocchk("simple_image::get_cmat ",alloc_stat)
    end function get_cmatfull

    !! get cmat value at index phys
    pure function get_cmat_at_1( self, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  ::  phys(3)
        complex :: comp
        comp = self%cmat(phys(1),phys(2),phys(3))
    end function get_cmat_at_1

    !! get cmat value at index phys
    elemental pure function get_cmat_at_2( self, h,k,l ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: h,k,l
        complex :: comp
        comp = self%cmat(h,k,l)
    end function get_cmat_at_2

    pure subroutine set_cmat_1( self, cmat )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        self%ft = .true.
        self%cmat = cmat
    end subroutine set_cmat_1

    pure subroutine set_cmat_2( self, cval )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: cval
        self%ft = .true.
        self%cmat = cval
    end subroutine set_cmat_2

    !! set comp to cmat at index phys
    subroutine set_cmat_at_1( self , phys , comp)
        class(image), intent(in) :: self
        integer, intent(in) :: phys(3)
        complex, intent(in) :: comp
        self%cmat(phys(1),phys(2),phys(3)) = comp
    end subroutine set_cmat_at_1

    !> add comp to cmat at index phys
    subroutine set_cmat_at_2( self, h, k, l, comp)
        class(image), intent(in) :: self
        integer, intent(in) :: h,k,l
        complex, intent(in) :: comp
        self%cmat(h,k,l) =  comp
    end subroutine set_cmat_at_2

    !> add comp to cmat at index phys
    subroutine add_cmat_at_1( self , phys , comp)
        class(image), intent(in) :: self
        integer, intent(in) :: phys(3)
        complex, intent(in) :: comp
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) + comp
    end subroutine add_cmat_at_1

    !> add comp to cmat at index (h,k,l)
    subroutine add_cmat_at_2( self, h, k, l, comp)
        class(image), intent(in) :: self
        integer, intent(in) :: h,k,l
        complex, intent(in) :: comp
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
        !$omp end workshare nowait
        !$omp do collapse(2) schedule(static)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                logi = [h,k,0]
                phys = self1%comp_addr_phys(logi)
                self3%cmat(phys(1),phys(2),phys(3)) = cmplx(expcmat3(h,k),0.)
                self4%cmat(phys(1),phys(2),phys(3)) = cmplx(expcmat4(h,k),0.)
            enddo
        enddo
        !$omp end do nowait
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
        !$omp end workshare nowait
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

    ! divide comp by cmat at index phys
    subroutine div_cmat_at_1( self, phys, rval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        real,              intent(in)    :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) / cmplx(rval,0.)
        else
            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
        endif
    end subroutine div_cmat_at_1

    !> div comp to cmat at index h,k,l
    subroutine div_cmat_at_2( self,h,k,l, rval)
        class(image), intent(in) :: self
        integer, intent(in) :: h,k,l
        real, intent(in) :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / cmplx(rval,0.)
        else
            self%cmat(h,k,l) =cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_2
    ! divide cmat at index phys by cval
    subroutine div_cmat_at_3( self, phys, cval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        complex,           intent(in)    :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) / cval
        else
            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
        endif
    end subroutine div_cmat_at_3

    !> div comp to cmat at index h,k,l
    subroutine div_cmat_at_4( self,h,k,l, cval)
        class(image), intent(in) :: self
        integer, intent(in)      :: h,k,l
        complex, intent(in)      :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / cval
        else
            self%cmat(h,k,l) =cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_4

    !! multiply cmat with real k at index phys
    pure subroutine mul_cmat_at_1( self, phys, rval)
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        real,              intent(in)    :: rval
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * cmplx(rval,0.)
    end subroutine mul_cmat_at_1

    !! multiply cmat with complex k at index phys
    pure subroutine mul_cmat_at_2( self, phys, cval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        complex,           intent(in)    :: cval
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * cval
    end subroutine mul_cmat_at_2

    !! multiply cmat at index h,k,l by real
    subroutine mul_cmat_at_3( self, h,k,l,rval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: h,k,l
        real,              intent(in)    :: rval
        self%cmat(h,k,l) = self%cmat(h,k,l) * cmplx(rval,0.)
    end subroutine mul_cmat_at_3

    !! multiply cmat at index h,k,l by complex
    subroutine mul_cmat_at_4( self, h,k,l,cval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: h,k,l
        complex,           intent(in)    :: cval
        self%cmat(h,k,l) = self%cmat(h,k,l) * cval
    end subroutine mul_cmat_at_4

    !> print_cmat
    !!
    subroutine print_cmat( self )
        class(image), intent(in) :: self
        print *, self%cmat
    end subroutine print_cmat

    !> print_rmat
    !!
    subroutine print_rmat( self )
        class(image), intent(in) :: self
        print *, self%rmat
    end subroutine print_rmat

    !>  \brief  expand_ft is for getting a Fourier plane using the old SIMPLE logics
    !! \return fplane a copy of this image object's fplane
    !!
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
    !! \param logi coordinates
    !! \param val new value
    !!
    subroutine set( self, logi, val )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: val
        if( logi(1) <= self%ldim(1) .and. logi(1) >= 1 .and. logi(2) <= self%ldim(2)&
            .and. logi(2) >= 1 .and. logi(3) <= self%ldim(3) .and. logi(3) >= 1 )then
            self%rmat(logi(1),logi(2),logi(3)) = val
        endif
    end subroutine set

    !>  \brief  set (replace) image data with new 3D data
    !! \param rmat new 3D data
    !!
    subroutine set_rmat( self, rmat )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(:,:,:)
        integer :: ldim(3)
        ldim(1) = size(rmat,1)
        ldim(2) = size(rmat,2)
        ldim(3) = size(rmat,3)
        if( all(self%ldim .eq. ldim) )then
            self%ft   = .false.
            self%rmat = 0.
            self%rmat(:ldim(1),:ldim(2),:ldim(3)) = rmat
        else
            write(*,*) 'ldim(rmat): ', ldim
            write(*,*) 'ldim(img): ', self%ldim
            call simple_stop('nonconforming dims; simple_image :: set_rmat')
        endif
    end subroutine set_rmat

    !>  \brief  set (replace) image data with new 3D data
    !! \param cmat new 3D complex data
    !!
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
            write(*,*) 'dim(cmat): ', cdim
            write(*,*) 'dim(self%cmat): ', self%array_shape
            call simple_stop('nonconforming dims; simple_image :: set_cmat')
        endif
    end subroutine set_cmat

    !> \brief  set_ldim replace image dimensions new 3D size
    !! \param ldim new 3D dimensions
    !!
    subroutine set_ldim( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        self%ldim = ldim
    end subroutine set_ldim

    !>  \brief set_smpd for setting smpd
    !! \param smpd  sampling distance
    !!
    subroutine set_smpd( self, smpd )
        class(image), intent(inout) :: self
        real,         intent(in)    :: smpd
        self%smpd = smpd
    end subroutine set_smpd

    !> \brief get_slice is for getting a slice from a volume
    !! \param self3d this image object
    !! \param slice index of slice in image
    !! \return self2d copy of slice as 2D image
    !!
    function get_slice( self3d, slice ) result( self2d )
        class(image), intent(in) :: self3d
        integer,      intent(in) :: slice
        type(image)              :: self2d
        call self2d%new([self3d%ldim(1),self3d%ldim(2),1],self3d%smpd)
        self2d%rmat(:,:,1) = self3d%rmat(:,:,slice)
    end function get_slice

    !>  \brief set_slice is for putting a slice into a volume
    !! \param self3d  this image object new slice as 2D image
    !! \param slice  index of slice in image
    !! \param self2d new slice as 2D image
    !!
    subroutine set_slice( self3d, slice, self2d )
        class(image), intent(in)    :: self2d
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self3d
        self3d%rmat(:,:,slice) = self2d%rmat(:,:,1)
    end subroutine set_slice

    !>  \brief get_npix is for getting the number of pixels for serialization
    !! \param mskrad mask radius
    !! \return npix num of pixels
    !!
    function get_npix( self, mskrad ) result( npix )
        class(image), intent(in) :: self
        real,         intent(in) :: mskrad
        real                     :: ci, cj, ck, e
        integer                  :: npix, i, j, k
        npix = 0
        ci = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        e = hardedge(ci,cj,ck,mskrad)
                    else
                        e = hardedge(ci,cj,mskrad)
                    endif
                    if( e > 0.5 )then
                        npix = npix+1
                    endif
                    ck = ck+1
                end do
                cj = cj+1.
            end do
            ci = ci+1.
        end do
    end function get_npix

    !>  \brief   get_lfny
    !! \param which
    !! \return fnyl
    !!
    pure function get_lfny( self, which ) result( fnyl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: fnyl
        fnyl = self%fit%get_lfny(which)
    end function get_lfny

    !>  \brief   get_lhp
    !! \param which
    !! \return lhp
    !!
    pure function get_lhp( self, which ) result( hpl )
        class(image), intent(in) :: self
        integer,      intent(in) :: which
        integer :: hpl
        hpl = self%fit%get_lhp(which)
    end function get_lhp

    !>  \brief   get_lp
    !! \param ind
    !! \return lp
    !!
    pure function get_lp( self, ind ) result( lp )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: lp
        lp = self%fit%get_lp(1, ind)
    end function get_lp

    !>  \brief   get_spat_freq
    !! \param ind
    !! \return spat_freq
    !!
    pure function get_spat_freq( self, ind ) result( spat_freq )
        class(image), intent(in) :: self
        integer,      intent(in) :: ind
        real                     :: spat_freq
        spat_freq = self%fit%get_spat_freq(1, ind)
    end function get_spat_freq

    !>  \brief  get_find
    !! \param res
    !! \return  ind
    !!
    pure function get_find( self, res ) result( ind )
        class(image), intent(in) :: self
        real,         intent(in) :: res
        integer :: ind
        ind = self%fit%get_find(1, res)
    end function get_find

    !>  \brief   get_clin_lims
    !! \param lp_dyn
    !! \return lims
    !!
    function get_clin_lims( self, lp_dyn ) result( lims )
        class(image), intent(in) :: self
        real,         intent(in) :: lp_dyn
        integer                  :: lims(2)
        lims = self%fit%get_clin_lims(lp_dyn)
    end function get_clin_lims

    !>  \brief  rmat_associated check rmat association
    !! \return  assoc
    !!
    function rmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
        assoc = associated(self%rmat)
    end function rmat_associated

    !>  \brief cmat_associated check cmat association
    !! \return  assoc
    !!
    function cmat_associated( self ) result( assoc )
        class(image), intent(in) :: self
        logical :: assoc
        assoc = associated(self%cmat)
    end function cmat_associated

    !>  \brief serialize is for packing/unpacking a serialized image vector for pca analysis
    !! \param pcavec analysis vector
    !! \param mskrad mask radius
    !!
    subroutine serialize( self, pcavec, mskrad )
        class(image),      intent(inout) :: self
        real, allocatable, intent(inout) :: pcavec(:)
        real, optional,    intent(in)    :: mskrad
        integer                          :: i, j, k, npix
        real                             :: ci, cj, ck, e
        logical                          :: pack, usemsk = .false.
        if( present(mskrad) )usemsk=.true.
        if( self%ft ) call simple_stop('ERROR, serialization not yet implemented for FTs; serialize; simple_image')
        if( usemsk )then
            npix = self%get_npix(mskrad)
        else
            npix = self%ldim(1)*self%ldim(2)*self%ldim(3)
        endif
        if( allocated(pcavec) )then
            if( size(pcavec) /= npix ) call simple_stop('size mismatch mask/npix; serialize; simple_image')
            pack = .false.
        else
            pack = .true.
            allocate( pcavec(npix), stat=alloc_stat )
            if(alloc_stat/=0)call allocchk('serialize; simple_image')
            pcavec = 0.
        endif
        npix = 0
        ci = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( usemsk )then
                        if( self%ldim(3) > 1 )then
                            e = hardedge(ci,cj,ck,mskrad)
                        else
                            e = hardedge(ci,cj,mskrad)
                        endif
                    endif
                    if( (e>0.5).or.(.not.usemsk) )then
                        npix = npix+1
                        if( pack )then
                            pcavec(npix) = self%rmat(i,j,k)
                        else
                            self%rmat(i,j,k) = pcavec(npix)
                        endif
                    endif
                    ck = ck+1
                end do
                cj = cj+1.
            end do
            ci = ci+1.
        end do
    end subroutine serialize

    !>  \brief winserialize is for packing/unpacking a serialized image vector for convolutional pca analysis
    !! \param coord coordinate offset
    !! \param winsz window size
    !! \param pcavec analysis vector
    !!
    subroutine winserialize( self, coord, winsz, pcavec )
        class(image),      intent(inout) :: self
        real, allocatable, intent(inout) :: pcavec(:)
        integer,           intent(in)    :: coord(:), winsz
        integer :: i, j, k, cnt, npix
        logical :: pack
        if( self%ft ) call simple_stop('ERROR, winserialization not yet implemented for Fourier; winserialize; simple_image')
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
            if( size(coord) < 3 ) call simple_stop('need a 3D coordinate for a 3D image; winserialize; simple_image')
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
                if( size(pcavec) /= npix ) call simple_stop('size mismatch mask/npix; winserialize; simple_image')
                pack = .false.
            else
                pack = .true.
                allocate( pcavec(npix), stat=alloc_stat )
                if(alloc_stat.ne.0)call allocchk('winserialize; simple_image',alloc_stat)
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
                    if( self%rmat(i,j,k) == 0. ) self%rmat(i,j,k) = 1.
                end do
            end do
        end do
    end subroutine zero2one

    !>  \brief get_fcomp for getting a Fourier component from the compact representation
    !! \param logi
    !! \param phys
    !! \return  comp
    !!
    pure function get_fcomp( self, logi, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3), phys(3)
        complex :: comp
        comp = self%cmat(phys(1),phys(2),phys(3))
        if( logi(1) < 0 ) comp = conjg(comp)
    end function get_fcomp

    !> \brief set_fcomp  for setting a Fourier component in the compact representation
    !! \param logi
    !! \param phys
    !! \param comp
    !!
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

    !> \brief add_fcomp  is for componentwise summation
    !! \param logi
    !! \param phys
    !! \param comp
    !!
    subroutine add_fcomp( self, logi, phys, comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3), phys(3)
        complex,      intent(in)    :: comp
        complex :: comp_here
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) + comp_here
    end subroutine add_fcomp

    !> \brief subtr_fcomp  is for componentwise summation
    !! \param logi
    !! \param phys
    !! \param comp
    !!
    subroutine subtr_fcomp( self, logi, phys, comp )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3), phys(3)
        complex,      intent(in)    :: comp
        complex :: comp_here
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) - comp_here
    end subroutine subtr_fcomp

    !>  \brief vis is for plotting an image
    !! \param sect
    !!
    subroutine vis( self, sect , geomorsphr)
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
            if(geomorsphr)then
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
    !! \param is
    !!
    subroutine set_ft( self, is )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: is
        self%ft = is
    end subroutine set_ft

    !>  \brief extr_fcomp is for extracting a Fourier component at arbitrary
    !>          position in a 2D transform using windowed sinc interpolation
    !! \param h,k Fourier coordinates
    !! \param x,y Image coordinates
    !! \return  comp
    !!
    function extr_fcomp( self, h, k, x, y ) result( comp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: h, k, x, y
        complex, allocatable :: comps(:,:)
        complex :: comp
        integer :: win(2,2), i, j, phys(3)
        if( self%ldim(3) > 1 )         call simple_stop('only 4 2D images; extr_fcomp; simple_image')
        if( .not. self%ft )            call simple_stop('image need to be FTed; extr_fcomp; simple_image')
        ! evenness and squareness are checked in the comlin class
        call sqwin_2d(h, k, 1., win)
        allocate( comps(win(1,1):win(1,2),win(2,1):win(2,2)) )
        do i=win(1,1),win(1,2)
            do j=win(2,1),win(2,2)
                phys       = self%comp_addr_phys([i,j,0])
                comps(i,j) = self%get_fcomp([i,j,0], phys)
            end do
            comps(i,:) = comps(i,:) * sinc(h-real(i))
        end do
        do i = win(2,1), win(2,2)
            comps(:,i) = comps(:,i) * sinc(k-real(i))
        enddo
        comp = sum(comps)
        deallocate(comps)
        ! origin shift
        if( (abs(x) > TINY) .and. (abs(y) > TINY) )then
            comp = comp*oshift_here(self%ldim(1)/2, h, k, x, y)
        endif

    contains

        !> oshift_here
        !! \param xdim
        !! \param x,y Image coords
        !! \param dx,dy pixel width
        !! \return  comp
        !!
        pure function oshift_here( xdim, x, y, dx, dy ) result( comp )
            integer, intent(in)  :: xdim
            real, intent(in)     :: x, y, dx, dy
            complex              :: comp
            real                 :: arg
            arg = (pi/real(xdim)) * dot_product([x,y], [dx,dy])
            comp = cmplx(cos(arg),sin(arg))
        end function oshift_here

    end function extr_fcomp

    !> \brief packer  replaces the pack intrinsic because the Intel compiler bugs out
    !! \param mskimg
    !! \return  pixels
    !!
    function packer( self, mskimg ) result( pixels )
        class(image),           intent(in) :: self
        class(image), optional, intent(in) :: mskimg
        real, allocatable :: pixels(:)
        integer :: nsel,i,j,k,cnt
        if( present(mskimg) )then
            nsel = count(mskimg%rmat<0.5)
            allocate(pixels(nsel))
            cnt = 0
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    do k=1,self%ldim(3)
                        if( mskimg%rmat(i,j,k) > 0.5 )then
                            cnt = cnt + 1
                            pixels(cnt) = self%rmat(i,j,k)
                        endif
                    end do
                end do
            end do
        else
            nsel = product(self%ldim)
            allocate(pixels(nsel))
            cnt = 0
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    do k=1,self%ldim(3)
                        cnt = cnt + 1
                        pixels(cnt) = self%rmat(i,j,k)
                    end do
                end do
            end do
        endif
    end function packer

    ! CHECKUPS

    !>  \brief  Checks for existence
    !! \return logical flag if image object exists
    pure function exists( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%existence
    end function exists

    !>  \brief  Checks whether the image is 2D
    !! \return logical flag if image object is 2D
    pure logical function is_2d(self)
        class(image), intent(in)  ::  self
        is_2d = count(self%ldim .eq. 1) .eq. 1
    end function is_2d

    !>  \brief  Checks whether the image is 3D
    !! \return logical flag if image object is 3D
    pure logical function is_3d(self)
        class(image), intent(in)  ::  self
        is_3d = .not. any(self%ldim .eq. 1)
    end function is_3d

    !>  \brief  checks for even dimensions
    !! \return logical flag if image object has even dimensions
    pure function even_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep, test(2)
        test = .false.
        test(1) = is_even(self%ldim(1))
        test(2) = is_even(self%ldim(2))
        yep = all(test)
    end function even_dims

    !>  \brief  checks for square dimensions
    !! \return logical flag if image object has square dimensions
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
    !!
    !! \param self1 image object
    !! \param self2 image object
    !! \return logical flag if two image objects have same dimensions
    pure function same_dims_1( self1, self2 ) result( yep )
        class(image), intent(in) :: self1, self2
        logical :: yep, test(3)
        test = .false.
        test(1) = self1%ldim(1) == self2%ldim(1)
        test(2) = self1%ldim(2) == self2%ldim(2)
        test(3) = self1%ldim(3) == self2%ldim(3)
        yep = all(test)
    end function same_dims_1

    !>  \brief  checks for same dimensions
    !! \return logical flag if image object has same dimensions as ldim
    pure function same_dims( self1, ldim ) result( yep )
        class(image), intent(in) :: self1
        integer,      intent(in) :: ldim(3) !< dimensions
        logical :: yep, test(3)
        test = .false.
        test(1) = self1%ldim(1) == ldim(1)
        test(2) = self1%ldim(2) == ldim(2)
        test(3) = self1%ldim(3) == ldim(3)
        yep = all(test)
    end function same_dims

    !>  \brief  checks for same sampling distance, overloaded as (.eqsmpd.)
    !!
    !! \param self1 image object
    !! \param self2 image object
    !! \return logical flag if image objects have same sampling distance
    pure  function same_smpd( self1, self2 ) result( yep )
        class(image), intent(in) :: self1, self2
        logical :: yep
        if( abs(self1%smpd-self2%smpd) < 0.0001 )then
            yep = .true.
        else
            yep = .false.
        endif
    end function same_smpd

    !>  \brief  checks if image is ft
    !! \return logical flag if image objects have same kind
    pure function is_ft( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%ft
    end function is_ft

    ! ARITHMETICS

    !>  \brief  assign, polymorphic assignment (=)
    !! \param selfout rhs
    !! \param selfin lhs
    !!
    subroutine assign( selfout, selfin )
        class(image), intent(inout) :: selfout
        class(image), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    !>  \brief assign_r2img real constant to image assignment(=) operation
    !! \param realin fixed value for image
    !!
    subroutine assign_r2img( self, realin )
        class(image), intent(inout) :: self
        real,         intent(in)    :: realin
        self%rmat = realin
        self%ft = .false.
    end subroutine assign_r2img

    !>  \brief  assign_c2img  complex constant to image assignment(=) operation
    !! \param compin fixed value for image
    !!
    subroutine assign_c2img( self, compin )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: compin
        self%cmat = compin
        self%ft = .true.
    end subroutine assign_c2img

    !>  \brief  is for image addition(+) addition
    !! \param self1 image object 1
    !! \param self2  image object 2
    !! \return lhs, copy of added images
    !!
    function addition( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .neqv. self2%ft )then
                call simple_stop('cannot add images of different FT state; addition(+); simple_image')
            endif
            self%rmat = self1%rmat+self2%rmat
        else
            call simple_stop('cannot add images of different dims; addition(+); simple_image')
        endif
        self%ft = self1%ft
    end function addition

    !>  \brief  l1norm_1 is for l1 norm calculation
    !! \param self1 image object 1
    !! \param self2 image object 2
    !! \return  lhs, copy of l1 normed images
    !!
    function l1norm_1( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .neqv. self2%ft )then
                call simple_stop('cannot process images of different FT state; l1norm_1; simple_image')
            endif
            if( self1%ft )then
                self%cmat = cabs(self1%cmat-self2%cmat)
            else
                self%rmat = abs(self1%rmat-self2%rmat)
            endif
        else
            call simple_stop('cannot process images of different dims; l1norm_1; simple_image')
        endif
    end function l1norm_1

    !>  \brief l1norm_2 is for l1 norm calculation
    !! \param self1 image object 1
    !! \param self2 image object 2
    !! \return  lhs, copy of l1 normed images
    !!
    function l1norm_2( self1, self2 ) result( l1 )
        class(image), intent(in) :: self1, self2
        real :: l1
        if( self1%ft .neqv. self2%ft ) call simple_stop('images in different FT states; l1norm_2; simple_image')
        if(.not.(self1.eqdims.self2 ))call simple_stop('cannot process images of different dims; l1norm_2; simple_image')
        if( self1%ft .and. self2%ft )then
            l1 = sum(cabs(self1%cmat-self2%cmat))
        else
            l1 = sum(abs(self1%rmat-self2%rmat))
        endif
    end function l1norm_2

    !>  \brief l1weights is for l1 norm weight generation
    !! \param self1 image object 1
    !! \param self2 image object 2
    !! \param nvar normalisation variable
    !! \return  lhs, copy of l1 normed images
    !!
    function l1weights( self1, self2, nvar ) result( self )
        class(image),   intent(in) :: self1, self2
        real, optional, intent(in) :: nvar
        type(image) :: self
        integer :: i, j, k
        real :: sumw, nnvar
        if( self1%ft .or. self2%ft ) call simple_stop('not impemented for FTs; l1weights; simple_image')
        nnvar = 1.
        if( present(nvar) ) nnvar = nvar
        self = self1.lone.self2
        sumw = 0.
        do i=1,self1%ldim(1)
            do j=1,self1%ldim(2)
                do k=1,self1%ldim(3)
                    self%rmat(i,j,k) = exp(-self%rmat(i,j,k)/nvar)
                    sumw = sumw+self%rmat(i,j,k)
                end do
            end do
        end do
        call self%div(sumw)
    end function l1weights

    !>  \brief add_1 is for image summation, not overloaded
    !! \param self_to_add image object
    !! \param w
    !!
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
    !! \param logi image dimensions
    !! \param comp
    !! \param phys_in
    !! \param phys_out
    !!
    subroutine add_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(in)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) call simple_stop('cannot add complex number to real image; add_2; simple_image')
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
    !! \param rcomp
    !! \param i,j,k index
    !!
    subroutine add_3( self, rcomp, i, j, k )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rcomp
        integer,      intent(in)    :: i, j, k
        if(  self%ft ) call simple_stop('cannot add real number to transform; add_3; simple_image')
        self%rmat(i,j,k) = self%rmat(i,j,k)+rcomp
    end subroutine add_3

    !>  \brief add_4 is for componentwise weighted summation with kernel division, not overloaded
    !! \param logi coordinates
    !! \param comp complex additive input (denominator)
    !! \param w additive input (numerator multiplier)
    !! \param k componentwise additive input (numerator)
    !!
    subroutine add_4( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) call simple_stop('cannot add complex number to real image; add_2; simple_image')
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
    !! \param c complex additive input
    subroutine add_5( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        if( self%ft )then
            self%cmat = self%cmat+cmplx(c,0.)
        else
            self%rmat = self%rmat+c
        endif
    end subroutine add_5

    subroutine add_workshare_1( self1, self1_to_add, self2, self2_to_add, self3, self3_to_add, self4, self4_to_add )
        class(image),   intent(inout) :: self1, self2, self3, self4
        class(image),   intent(in)    :: self1_to_add, self2_to_add, self3_to_add, self4_to_add
        !$omp parallel workshare proc_bind(close)
        self1%cmat = self1%cmat + self1_to_add%cmat
        self2%cmat = self2%cmat + self2_to_add%cmat
        self3%cmat = self3%cmat + self3_to_add%cmat
        self4%cmat = self4%cmat + self4_to_add%cmat
        !$omp end parallel workshare
    end subroutine add_workshare_1

    subroutine add_workshare_2( self, self_to_add, rho, rho_to_add )
        class(image),       intent(inout) :: self
        class(image),       intent(in)    :: self_to_add
        real(kind=c_float), intent(inout) :: rho(:,:,:)
        real(kind=c_float), intent(in)    :: rho_to_add(:,:,:)
        if( self%ft )then
            !$omp parallel workshare proc_bind(close)
            self%cmat = self%cmat+self_to_add%cmat
            rho       = rho + rho_to_add
            !$omp end parallel workshare
        else
            !$omp parallel workshare proc_bind(close)
            self%rmat = self%rmat+self_to_add%rmat
            rho       = rho + rho_to_add
            !$omp end parallel workshare
        endif
    end subroutine add_workshare_2

    !>  \brief subtraction is for image subtraction(-)
    !! \param self_from lhs
    !! \param self_to  lhs subtractor
    !! \return copy of self
    !!
    function subtraction( self_from, self_to ) result( self )
        class(image), intent(in) :: self_from, self_to
        type(image) :: self
        if( self_from.eqdims.self_to )then
            call self%new(self_from%ldim, self_from%smpd)
            if( self_from%ft .neqv. self_to%ft )then
                call simple_stop('cannot subtract images of different FT state; subtraction(+); simple_image')
            endif
            self%rmat = self_from%rmat-self_to%rmat
        else
            call simple_stop('cannot subtract images of different dims; subtraction(-); simple_image')
        endif
    end function subtraction

    !>  \brief subtr_1 is for image subtraction,  not overloaded
    !! \param self_to_subtr image object
    !! \param w
    !!
    subroutine subtr_1( self, self_to_subtr, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_subtr
        real, optional, intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( self.eqdims.self_to_subtr )then
            if( self%ft .eqv. self_to_subtr%ft )then
                if( self%ft )then
                    self%cmat = self%cmat-ww*self_to_subtr%cmat
                else
                    self%rmat = self%rmat-ww*self_to_subtr%rmat
                endif
            else
                call simple_stop('cannot subtract images with different FT status; subtr_1; simple_image')
            endif
        else
            call simple_stop('cannot subtract images of different dims; subtr_1; simple_image')
        endif
    end subroutine subtr_1

    !>  \brief subtr_2 is for componentwise subtraction, not overloaded
    !! \param logi
    !! \param comp
    !! \param phys_in
    !! \param phys_out
    !!
    subroutine subtr_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(out)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) call simple_stop('cannot subtract complex number from real image; subtr_2; simple_image')
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

    !>  \brief subtr_3 is for componentwise weighted subtraction with kernel division, not overloaded
    !! \param logi   coordinates
    !! \param comp   complex subtraction input (denominator)
    !! \param w      subtraction input (numerator multiplier)
    !! \param k      componentwise subtraction input (numerator)
    !!
    subroutine subtr_3( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) call simple_stop('cannot subtract complex number from real image; subtr_3; simple_image')
        phys = self%fit%comp_addr_phys(logi)
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        if( abs(k(phys(1),phys(2),phys(3))) > TINY )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) -&
                &(comp_here/k(phys(1),phys(2),phys(3)))*w
        endif
    end subroutine subtr_3

    !>  \brief subtr_4 is for subtracting a constant from a real image, not overloaded
    !! \param c constant
    !!
    subroutine subtr_4( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        self%rmat = self%rmat-c
    end subroutine subtr_4

    !>  \brief multiplication is for image multiplication(*)
    !! \param self1 image object
    !! \param self2 image object
    !! \return  self
    !!
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
            call simple_stop('cannot multiply images of different dims; multiplication(*); simple_image')
        endif
    end function multiplication

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
    !! \param logi
    !! \param rc
    !! \param phys_in
    !! \param phys_out
    !!
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
    !! \param rc multiplier
    !!
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
    !! \param self2mul 3D multiplier
    !!
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
            call simple_stop('cannot multiply images of different dims; mul_3; simple_image')
        endif
    end subroutine mul_3

    !>  \brief mul_4 is for low-pass limited multiplication of images
    !! \param self2mul 3D multiplier
    !! \param lp cut off filter frequency
    !!
    subroutine mul_4( self, self2mul, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
        real,         intent(in)    :: lp
        integer                     :: lims(3,2),sqlim,h,k,l,phys(3)
        if( .not. self%is_ft() )     call simple_stop('low-pass limited multiplication requires self to be FT')
        if( .not. self2mul%is_ft() ) call simple_stop('low-pass limited multiplication requires self2mul to be FT')
        if( self.eqdims.self2mul )then
            lims = self%fit%loop_lims(1,lp)
            sqlim = (maxval(lims(:,2)))**2
            !$omp parallel do collapse(3) default(shared)&
            !$omp private(h,k,l,phys) schedule(static) proc_bind(close)
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
            !$omp end parallel do
        else
            call simple_stop('cannot multiply images of different dims; mul_3; simple_image')
        endif
    end subroutine mul_4

    !>  \brief division is for image division(/)
    !! \param self1 lhs numerator
    !! \param self2 lhs denominator
    !! \return rhs image copy
    !!
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
            call simple_stop('cannot divide images of different dims; division(/); simple_image')
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
    !! \param c divisor
    !!
    subroutine div_1( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        if( abs(c) < 1e-6 )then
            call simple_stop('division with zero; div; simple_image')
        else
            if( self%ft )then
                self%cmat = self%cmat/c
            else
                self%rmat = self%rmat/c
            endif
        endif
    end subroutine div_1

    !>  \brief div_2 is for component-wise matrix division of a Fourier transform with a real matrix, k
    !! \param logi coordinates
    !! \param k 3D divisor
    !! \param square logical flag if return image should be squared
    !!
    subroutine div_2( self, logi, k, square )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: k(:,:,:)
        logical,      intent(in)    :: square
        integer :: phys(3)
        if( self%ft )then
            phys = self%fit%comp_addr_phys(logi)
            if( abs(k(phys(1),phys(2),phys(3))) > 1e-6 )then
                if( square )then
                    self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/(k(phys(1),phys(2),phys(3))**2.)
                else
                    self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/k(phys(1),phys(2),phys(3))
                endif
            else
                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
            endif
        else
            call simple_stop('Image need to be Fourier transformed; simple_image::div_2')
        endif
    end subroutine div_2

    !>  \brief div_3 is for component-wise division of an image with a real number
    !! \param logi coordinates
    !! \param k 3D divisor
    !! \param phys_in
    !!
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
    !! \param self2div image object divisor
    !!
    subroutine div_4( self, self2div )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2div
        if( self.eqdims.self2div )then
            if( self%ft .and. self2div%ft )then
                self%cmat = self%cmat/self2div%cmat
            else if( self%ft .eqv. self2div%ft )then
                self%rmat = self%rmat/self2div%rmat
                self%ft = .false.
            else if(self%ft)then
                self%cmat = self%cmat/self2div%rmat
            else
                self%cmat = self%rmat/self2div%cmat
                self%ft = .true.
            endif
        else
            call simple_stop('cannot divide images of different dims; div_4; simple_image')
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
                    sh    = nint(hyp(real(h),real(k),real(l)))
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

    !>  \brief conjugate is for complex conjugation of a FT
    !! \return self_out
    !!
    function conjugate( self ) result ( self_out )
        class(image), intent(in) :: self
        type(image) :: self_out
        if( self%ft )then
            call self_out%copy(self)
            self%cmat = conjg(self%cmat)
        else
            write(*,'(a)') "WARNING! Cannot conjugate real image"
        endif
    end function conjugate

    !>  \brief sqpow is for calculating the square power of an image
    !!
    subroutine sqpow( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            self%cmat = (self%cmat*conjg(self%cmat))**2.
        else
            self%rmat = self%rmat*self%rmat
        endif
    end subroutine sqpow

    !>  \brief signswap_aimag is changing the sign of the imaginary part of the Fourier transform
    subroutine signswap_aimag( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            self%cmat = cmplx(real(self%cmat),-aimag(self%cmat))
        else
            call self%fft()
            self%cmat = cmplx(real(self%cmat),-aimag(self%cmat))
            call self%ifft()
        endif
    end subroutine signswap_aimag

    !>  \brief  signswap_real is changing the sign of the real part of the Fourier transform
    subroutine signswap_real( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            self%cmat = cmplx(-real(self%cmat),aimag(self%cmat))
        else
            call self%fft()
            self%cmat = cmplx(-real(self%cmat),aimag(self%cmat))
            call self%ifft()
        endif
    end subroutine signswap_real

    ! BINARY IMAGE METHODS

    !>  \brief nforeground counts the number of foreground (white) pixels in a binary image
    !! \return num of ON pixels
    !!
    function nforeground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = count(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.5)
    end function nforeground

    !>  \brief nbackground counts the number of background (black) pixels in a binary image
    !! \return num of OFF pixels
    !!
    function nbackground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = product(self%ldim)-self%nforeground()
    end function nbackground

    !>  \brief  is for binarizing an image with given threshold value
    !!          binary normalization (norm_bin) assumed!> bin_1
    !! \param thres threshold value
    !!
    subroutine bin_1( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        if( self%ft ) call simple_stop('only for real images; bin_1; simple image')
        where( self%rmat >= thres )
            self%rmat = 1.
        elsewhere
            self%rmat = 0.
        end where
    end subroutine bin_1

    !>  \brief  bin_2 is for binarizing an image using nr of pixels/voxels threshold
    !! \param npix
    !!
    subroutine bin_2( self, npix )
        class(image), intent(inout) :: self
        integer, intent(in)         :: npix
        real, allocatable           :: forsort(:)
        real                        :: thres
        integer                     :: npixtot
        if( self%ft ) call simple_stop('only for real images; bin_2; simple image')
        npixtot = product(self%ldim)
        ! forsort = self%packer() ! Intel hickup
        !$ allocate(forsort(self%ldim(1)*self%ldim(2)*self%ldim(3)),stat=alloc_stat)
        !$ if(alloc_stat.ne.0)call allocchk("simple_image::bin_2 packing rmat for sorting - PGI prealloc", alloc_stat)
        forsort = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), .true.)
        call hpsort(forsort)
        thres = forsort(npixtot-npix-1) ! everyting above this value 1 else 0
        call self%bin( thres )
        deallocate( forsort )
    end subroutine bin_2

    !>  \brief  is for binarizing an image using k-means to identify the background/
    !!          foreground distributions for the image
    subroutine bin_kmeans( self )
        class(image), intent(inout) :: self
        logical :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        real    :: cen1, cen2, sum1, sum2, val1, val2, sumval
        real    :: foreground_cen, background_cen
        integer :: cnt1, cnt2,  l, npix
        integer, parameter :: MAXITS=100
        if( self%ft ) stop 'only for real images; bin_kmeans; simple image'
        ! estimate background value around the edges of the box
        cen1 = 0.
        if( self%ldim(3) == 1 )then
            cen1  = cen1 + sum(self%rmat( 1           , :self%ldim(2),1))
            cen1  = cen1 + sum(self%rmat( self%ldim(1), :self%ldim(2),1))
            cen1  = cen1 + sum(self%rmat(:self%ldim(1),  1,           1))
            cen1  = cen1 + sum(self%rmat(:self%ldim(1),  self%ldim(2),1))
            cen1  = cen1 / real(4 * self%ldim(1))
        else
            cen1  = cen1 + sum(self%rmat( 1           ,  1            , :self%ldim(3)))
            cen1  = cen1 + sum(self%rmat( 1           ,  self%ldim(2) , :self%ldim(3)))
            cen1  = cen1 + sum(self%rmat( self%ldim(1),  1            , :self%ldim(3)))
            cen1  = cen1 + sum(self%rmat( self%ldim(1),  self%ldim(2) , :self%ldim(3)))
            cen1  = cen1 + sum(self%rmat( 1           , :self%ldim(2) ,  1           ))
            cen1  = cen1 + sum(self%rmat( 1           , :self%ldim(2) ,  self%ldim(3)))
            cen1  = cen1 + sum(self%rmat( self%ldim(1), :self%ldim(2) ,  1           ))
            cen1  = cen1 + sum(self%rmat( self%ldim(1), :self%ldim(2) ,  self%ldim(3)))
            cen1  = cen1 + sum(self%rmat(:self%ldim(1),  1            ,  1           ))
            cen1  = cen1 + sum(self%rmat(:self%ldim(1),  1            ,  self%ldim(3)))
            cen1  = cen1 + sum(self%rmat(:self%ldim(1),  self%ldim(2) ,  1           ))
            cen1  = cen1 + sum(self%rmat(:self%ldim(1),  self%ldim(2) ,  self%ldim(3)))
            cen1  = cen1 / real(12 * self%ldim(1))
        endif
        sumval = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        npix   = product(self%ldim)
        ! estimate foreground value as global average
        cen2   = sumval / real(npix)
        ! foreground/background identification with k-means
        do l=1,MAXITS
            where( (cen1 - self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))**2.0 <&
                &(cen2 - self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))**2.0 )
                l_msk = .true.
            else where
                l_msk = .false.
            end where
            cnt1 = count(l_msk)
            cnt2 = npix - cnt1
            sum1 = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), l_msk)
            sum2 = sumval - sum1
            cen1 = sum1 / real(cnt1)
            cen2 = sum2 / real(cnt2)
        end do
        ! assign values to the centers
        if( cen1 > cen2 )then
            val1           = 1.
            val2           = 0.
            foreground_cen = cen1
            background_cen = cen2
        else
            val1           = 0.
            val2           = 1.
            foreground_cen = cen2
            background_cen = cen1
        endif
        ! binarize the image
        where( (background_cen - self%rmat)**2. < (foreground_cen - self%rmat)**2. )
            self%rmat = val1
        elsewhere
            self%rmat = val2
        end where
    end subroutine bin_kmeans

    !>  \brief bin_filament is for creating a binary filament
    !! \param width_A physical width
    !!
    subroutine bin_filament( self, width_A )
        class(image), intent(inout) :: self
        real, intent(in)            :: width_A
        real    :: halfwidth_pix, cen_xdim
        integer :: xstart, xstop
        if( self%ldim(3) > 1 ) call simple_stop('only for 2D images; simple_imge :: bin_filament ')
        halfwidth_pix = (width_A/self%smpd)/2.
        cen_xdim      = real(self%ldim(1))/2.
        xstart        = floor(cen_xdim-halfwidth_pix)
        xstop         = ceiling(cen_xdim+halfwidth_pix)
        if( self%ft ) self%ft = .false.
        self%rmat = 0.
        self%rmat(xstart:xstop,:,1) = 1.
    end subroutine bin_filament

    !>  \brief bin_cylinder  is for creating a binary cyclinder along z-axis
    !! \param rad radius
    !! \param height height of cylinder
    !!
    subroutine bin_cylinder( self, rad, height )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rad, height
        type(image)       :: mask2d
        real, allocatable :: plane(:,:,:)
        real        :: centre(3)
        integer     :: k
        if( self%ldim(3) == 1 ) call simple_stop('only for 3D images; simple_imge :: bin_cylinder ')
        centre = 1. + real(self%ldim-1)/2.
        if( self%ft ) self%ft = .false.
        call mask2d%new([self%ldim(1), self%ldim(2), 1], self%smpd)
        mask2d = 1.
        call mask2d%mask(rad, 'hard')
        !$ allocate(plane(self%ldim(1),self%ldim(2),self%ldim(3)),stat=alloc_stat)
        !$ if(alloc_stat.ne.0)call allocchk("In simple_image::bin_cylinder PGI prealloc plane",alloc_stat)
        plane = mask2d%get_rmat()
        self%rmat = 0.
        do k = 1,self%ldim(3)
            if( abs(real(k)-centre(3)) < height/2. )then
                self%rmat(:self%ldim(1),:self%ldim(2),k) = plane(:,:,1)
            endif
        enddo
        call mask2d%kill()
        deallocate(plane)
    end subroutine bin_cylinder

    !>  \brief cendist produces an image with square distance from the centre of the image
    !!
    subroutine cendist( self )
        class(image), intent(inout) :: self
        real    :: centre(3)
        integer :: i
        if( self%ft ) call simple_stop('real space only; simple_image%cendist')
        ! Builds square distance image
        self   = 0.
        centre = real(self%ldim-1)/2.
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

    !>  \brief masscen is for determining the center of mass of binarised image
    !!          only use this function for integer pixels shifting
    !! \return  xyz
    !!
    function masscen( self ) result( xyz )
        class(image), intent(inout) :: self
        real    :: xyz(3), spix, ci, cj, ck
        integer :: i, j, k
        if( self%ft ) call simple_stop('masscen not implemented for FTs; masscen; simple_image')
        spix = 0.
        xyz  = 0.
        ci   = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    xyz  = xyz  + self%rmat(i,j,k) * [ci, cj, ck]
                    spix = spix + self%rmat(i,j,k)
                    ck   = ck + 1.
                end do
                cj = cj + 1.
            end do
            ci = ci + 1.
        end do
        xyz = xyz / spix
        if( self%ldim(3) == 1 ) xyz(3) = 0.
    end function masscen

    !> masscen as a subroutine
    subroutine masscens( self, xyz )
        class(image), intent(inout) :: self
        real        , intent(out)   :: xyz(3)
        real    ::  spix, ci, cj, ck
        integer :: i, j, k
        if( self%is_ft() ) call simple_stop('masscen not implemented for FTs; masscen; simple_image')
        spix = 0.
        xyz  = 0.
        ci   = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    xyz  = xyz  + self%rmat(i,j,k) * [ci, cj, ck]
                    spix = spix + self%rmat(i,j,k)
                    ck   = ck + 1.
                end do
                cj = cj + 1.
            end do
            ci = ci + 1.
        end do
        xyz = xyz / spix
        if( self%ldim(3) == 1 ) xyz(3) = 0.
    end subroutine masscens

    !>  \brief center is for centering an image based on center of mass
    !! \param lp low-pass cut-off freq
    !! \param msk mask
    !! \return xyz shift
    !!
    function center( self, lp, msk ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real, optional, intent(in)    :: msk
        type(image) :: tmp
        real        :: xyz(3), rmsk
        call tmp%copy(self)
        call tmp%bp(0., lp)
        call tmp%ifft()
        if( present(msk) )then
            rmsk = msk
        else
            rmsk = real( self%ldim(1) )/2. - 5. ! 5 pixels outer width
        endif
        call tmp%mask(rmsk, 'soft')
        ! where( tmp%rmat < 0. )
        !     tmp%rmat = 0.
        ! end where
        call tmp%bin_kmeans
        xyz = tmp%masscen()
        call tmp%kill
    end function center

    function center_serial( self, lp, msk ) result( xyz )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
        real,         intent(in)    :: msk
        real    :: xyz(3)
        integer :: ithr
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! copy rmat
        thread_safe_tmp_imgs(ithr)%rmat = self%rmat
        thread_safe_tmp_imgs(ithr)%ft   = .false.
        call thread_safe_tmp_imgs(ithr)%fft()
        call thread_safe_tmp_imgs(ithr)%bp(0., lp)
        call thread_safe_tmp_imgs(ithr)%ifft()
        call thread_safe_tmp_imgs(ithr)%mask(msk, 'soft')
        call thread_safe_tmp_imgs(ithr)%bin_kmeans
        xyz = thread_safe_tmp_imgs(ithr)%masscen()
    end function center_serial

    !>  \brief bin_inv inverts a binary image
    subroutine bin_inv( self )
        class(image), intent(inout) :: self
        self%rmat = -1.*(self%rmat-1.)
    end subroutine bin_inv

    !>  \brief grow_bin adds one layer of pixels bordering the background in a binary image
    !! Classical dilation of binary image
    subroutine grow_bin( self )
        class(image), intent(inout) :: self
        integer                     :: i,j,k
        integer                     :: il,ir,jl,jr,kl,kr
        logical, allocatable        :: add_pixels(:,:,:)
        if( self%ft ) call simple_stop('only for real images; grow_bin; simple image')
        allocate( add_pixels(self%ldim(1),self%ldim(2),self%ldim(3)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('grow_bin; simple_image')
        ! Figure out which pixels to add
        add_pixels = .false.
        if( self%ldim(3) == 1 )then
            do i=1,self%ldim(1)
                il = max(1,i-1)
                ir = min(self%ldim(1),i+1)
                do j=1,self%ldim(2)
                    if (self%rmat(i,j,1) < TINY) then
                        jl = max(1,j-1)
                        jr = min(self%ldim(2),j+1)
                        if( any(abs(self%rmat(il:ir,jl:jr,1)-1.) < TINY) )add_pixels(i,j,1) = .true.
                    end if
                end do
            end do
            ! add
            forall( i=1:self%ldim(1), j=1:self%ldim(2), add_pixels(i,j,1) )self%rmat(i,j,1) = 1.
        else
            do i=1,self%ldim(1)
                il = max(1,i-1)
                ir = min(self%ldim(1),i+1)
                do j=1,self%ldim(2)
                    jl = max(1,j-1)
                    jr = min(self%ldim(2),j+1)
                    do k=1,self%ldim(3)
                        if (abs(self%rmat(i,j,k)) < TINY) then
                            kl = max(1,k-1)
                            kr = min(self%ldim(3),k+1)
                            if( any(abs(self%rmat(il:ir,jl:jr,kl:kr)-1.) < TINY )) add_pixels(i,j,k) = .true.
                        end if
                    end do
                end do
            end do
            ! add
            forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), add_pixels(i,j,k) ) &
                & self%rmat(i,j,k) = 1.
        endif
        deallocate( add_pixels )
    end subroutine grow_bin

    !>  \brief shrink_bin removes one layer of pixels bordering the background in a binary image
    !! Classical erosion of binary image
    subroutine shrink_bin( self )
        class(image), intent(inout) :: self
        integer                     :: i,j,k
        integer                     :: il,ir,jl,jr,kl,kr
        logical, allocatable        :: sub_pixels(:,:,:)
        if( self%ft ) call simple_stop('only for real images; shrink_bin; simple image')
        allocate( sub_pixels(self%ldim(1),self%ldim(2),self%ldim(3)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('shrink_bin; simple_image')
        ! Figure out which pixels to remove
        sub_pixels = .false.
        if( self%ldim(3) == 1 )then
            do i=1,self%ldim(1)
                il = max(1,i-1)
                ir = min(self%ldim(1),i+1)
                do j=1,self%ldim(2)
                    if (self%rmat(i,j,1) == 0.) then
                        jl = max(1,j-1)
                        jr = min(self%ldim(2),j+1)
                        if( any(self%rmat(il:ir,jl:jr,1) == 1.) ) &
                            sub_pixels(i,j,1) = .true.
                    end if
                end do
            end do
            ! remove
            forall( i=1:self%ldim(1), j=1:self%ldim(2), sub_pixels(i,j,1) )self%rmat(i,j,1) = 0.
        else
            do i=1,self%ldim(1)
                il = max(1,i-1)
                ir = min(self%ldim(1),i+1)
                do j=1,self%ldim(2)
                    jl = max(1,j-1)
                    jr = min(self%ldim(2),j+1)
                    do k=1,self%ldim(3)
                        if (self%rmat(i,j,k) == 0.) then
                            kl = max(1,k-1)
                            kr = min(self%ldim(3),k+1)
                            if( any(self%rmat(il:ir,jl:jr,kl:kr) == 1.) ) &
                                sub_pixels(i,j,k) = .true.
                        end if
                    end do
                end do
            end do
            ! remove
            forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), sub_pixels(i,j,k) ) &
                & self%rmat(i,j,k) = 0.
        endif
        deallocate( sub_pixels )
    end subroutine shrink_bin

    !> \brief grow_bins adds one layer of pixels bordering the background in a binary image
    !! \param nlayers
    !! Classical iterative dilation of binary image
    subroutine grow_bins( self, nlayers )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nlayers
        integer                     :: i,j,k, tsz(3,2), win(3,2), pdsz(3,2)
        logical, allocatable        :: add_pixels(:,:,:), template(:,:,:)
        if( self%ft ) call simple_stop('only for real images; grow_bin; simple image')
        tsz(:,1) = -nlayers
        tsz(:,2) = nlayers
        if(self%is_2d())tsz(3,:) = 1
        allocate( template(tsz(1,1):tsz(1,2), tsz(2,1):tsz(2,2), tsz(3,1):tsz(3,2)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('grow_bins; simple_image 2')
        pdsz(:,1) = 1 - nlayers
        pdsz(:,2) = self%ldim + nlayers
        if(self%is_2d())pdsz(3,:) = 1
        allocate( add_pixels(pdsz(1,1):pdsz(1,2), pdsz(2,1):pdsz(2,2),&
            &pdsz(3,1):pdsz(3,2)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('grow_bins; simple_image 1')
        ! template matrix
        template = .true.
        do i = tsz(1,1), tsz(1,2)
            do j = tsz(2,1), tsz(2,2)
                if(self%is_2d())then
                    if(dot_product([i,j], [i,j]) > nlayers**2) template(i,j,1) = .false.
                else
                    do k = tsz(3,1), tsz(3,2)
                        if(dot_product([i,j,k],[i,j,k]) > nlayers**2) template(i,j,k) = .false.
                    enddo
                endif
            enddo
        enddo
        ! init paddedd logical array
        add_pixels = .false.
        forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), self%rmat(i,j,k)==1. )&
            & add_pixels(i,j,k) = .true.
        ! cycle
        if( self%is_3d() )then
            do i = 1, self%ldim(1)
                if( .not.any(self%rmat(i,:,:) > 0.5) )cycle
                do j = 1, self%ldim(2)
                    if( .not.any(self%rmat(i,j,:) > 0.5) )cycle
                    win(1:2,1) = [i, j] - nlayers
                    win(1:2,2) = [i, j] + nlayers
                    do k = 1, self%ldim(3)
                        if (self%rmat(i,j,k) <= 0.5)cycle
                        win(3,1) = k - nlayers
                        win(3,2) = k + nlayers
                        add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
                            &add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2))&
                            &.or.template
                    enddo
                enddo
            enddo
        else
            do i=1,self%ldim(1)
                if( .not.any(self%rmat(i,:,1) > 0.5) )cycle
                do j=1,self%ldim(2)
                    win(1:2,1) = [i, j] - nlayers
                    win(1:2,2) = [i, j] + nlayers
                    if (self%rmat(i,j,1) <= 0.5)cycle
                    add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1) =&
                        &add_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1).or.template(:,:,1)
                enddo
            enddo
        endif
        ! finalize
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 0.
        forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), add_pixels(i,j,k) ) &
            & self%rmat(i,j,k) = 1.
        deallocate( template, add_pixels )
    end subroutine grow_bins

    !> \brief shrink_bins removes n layers of pixels bordering the background in a binary image
    !! \param nlayers
    !! Classical iterative erosion of binary image
    subroutine shrink_bins( self, nlayers )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nlayers
        integer                     :: i,j,k, tsz(3,2), win(3,2), pdsz(3,2)
        logical, allocatable        :: sub_pixels(:,:,:), template(:,:,:)
        if( self%ft ) call simple_stop('only for real images; shrink_bin; simple image')
        tsz(:,1) = -nlayers
        tsz(:,2) = nlayers
        if(self%is_2d())tsz(3,:) = 1
        allocate( template(tsz(1,1):tsz(1,2), tsz(2,1):tsz(2,2), tsz(3,1):tsz(3,2)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('shrink_bins; simple_image 2')
        pdsz(:,1) = 1 - nlayers
        pdsz(:,2) = self%ldim + nlayers
        if(self%is_2d())pdsz(3,:) = 1
        allocate( sub_pixels(pdsz(1,1):pdsz(1,2), pdsz(2,1):pdsz(2,2),&
            &pdsz(3,1):pdsz(3,2)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('shrink_bins; simple_image 1')
        ! template matrix
        template = .true.
        do i = tsz(1,1), tsz(1,2)
            do j = tsz(2,1), tsz(2,2)
                if(self%is_2d())then
                    if(dot_product([i,j], [i,j]) > nlayers**2) template(i,j,1) = .false.
                else
                    do k = tsz(3,1), tsz(3,2)
                        if(dot_product([i,j,k],[i,j,k]) > nlayers**2) template(i,j,k) = .false.
                    enddo
                endif
            enddo
        enddo
        ! init paddedd logical array
        sub_pixels = .false.

        forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), self%rmat(i,j,k)==1. )&
            & sub_pixels(i,j,k) = .true.
        ! cycle
        if( self%is_3d() )then
            do i = 1, self%ldim(1)
                if( .not.any(self%rmat(i,:,:) > 0.5) )cycle
                do j = 1, self%ldim(2)
                    if( .not.any(self%rmat(i,j,:) > 0.5) )cycle
                    win(1:2,1) = [i, j] - nlayers
                    win(1:2,2) = [i, j] + nlayers
                    do k = 1, self%ldim(3)
                        if (self%rmat(i,j,k) <= 0.5)cycle
                        win(3,1) = k - nlayers
                        win(3,2) = k + nlayers
                        sub_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
                            &sub_pixels(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2))&
                            &.or.template
                    enddo
                enddo
            enddo
        else
            do i=1,self%ldim(1)
                if( .not.any(self%rmat(i,:,1) > 0.5) )cycle
                do j=1,self%ldim(2)
                    win(1:2,1) = [i, j] - nlayers
                    win(1:2,2) = [i, j] + nlayers
                    if (self%rmat(i,j,1) <= 0.5)cycle
                    sub_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1) =&
                        &sub_pixels(win(1,1):win(1,2), win(2,1):win(2,2), 1).or.template(:,:,1)
                enddo
            enddo
        endif
        ! finalize
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 0.
        forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), sub_pixels(i,j,k) ) &
            & self%rmat(i,j,k) = 0.
        deallocate( template, sub_pixels )
    end subroutine shrink_bins

    !> binary_dilation wrapper for grow_bin(s)
    !! \param nlayers number of layers
    !!
    subroutine binary_dilation(self,nlayers)
        class(image), intent(inout) :: self
        integer, intent(inout),optional :: nlayers
        if(.not.present(nlayers))then
            call   self%grow_bin()
        else
            call   self%grow_bins(nlayers)
        end if
    end subroutine binary_dilation

    !> binary_erosion wrapper for shrink_bin(s)
    !! \param nlayers number of layers
    !!
    subroutine binary_erosion(self,nlayers)
        class(image), intent(inout) :: self
        integer, intent(inout),optional :: nlayers
        if(.not.present(nlayers))then
            call   self%shrink_bin()
        else
            call   self%shrink_bins(nlayers)
        end if
    end subroutine binary_erosion

    subroutine binary_opening (self,nlayers)
        class(image), intent(inout) :: self
        integer, intent(inout),optional :: nlayers
        integer:: i,n
        n=1; if(present(nlayers))n=nlayers
        do i=1,nlayers
            call  self%binary_erosion()
            call  self%binary_dilation()
        end do
    end subroutine binary_opening

    subroutine binary_closing (self,nlayers)
        class(image), intent(inout) :: self
        integer, intent(inout),optional  :: nlayers
        integer:: i,n
        n=1; if(present(nlayers))n=nlayers
        do i=1,n
            call  self%binary_dilation()
            call  self%binary_erosion()
        end do
    end subroutine binary_closing

    !>  \brief cos_edge applies cosine squared edge to a binary image
    !! \param falloff
    !!
    subroutine cos_edge( self, falloff )
        class(image), intent(inout) :: self
        integer, intent(in)         :: falloff
        real, allocatable           :: rmat(:,:,:)
        real                        :: rfalloff, scalefactor
        integer                     :: i, j, k, is, js, ks, ie, je, ke
        integer                     :: il, ir, jl, jr, kl, kr, falloff_sq
        if( falloff<=0 ) call simple_stop('stictly positive values for edge fall-off allowed; simple_image::cos_edge')
        if( self%ft )    call simple_stop('not intended for FTs; simple_image :: cos_edge')
        self%rmat   = self%rmat/maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        rfalloff    = real( falloff )
        falloff_sq  = falloff**2
        scalefactor = PI / rfalloff
        allocate( rmat(self%ldim(1),self%ldim(2),self%ldim(3)),stat=alloc_stat )
        if(alloc_stat /= 0)call allocchk("In simple_image::cos_edge")
        rmat = self%rmat(1:self%ldim(1),:,:)
        do i=1,self%ldim(1)
            is = max(1,i-1)                  ! left neighbour
            ie = min(i+1,self%ldim(1))       ! right neighbour
            il = max(1,i-falloff)            ! left bounding box limit
            ir = min(i+falloff,self%ldim(1)) ! right bounding box limit
            if( any(rmat(i,:,:)==1.) )cycle ! no values equal to one
            do j=1,self%ldim(2)
                js = max(1,j-1)
                je = min(j+1,self%ldim(2))
                jl = max(1,j-falloff)
                jr = min(j+falloff,self%ldim(2))
                if( self%ldim(3)==1 )then
                    ! 2d
                    if( rmat(i,j,1) == 1. ) cycle ! cycle if not equal to one
                    ! within mask region
                    ! update if has a masked neighbour
                    if( any( rmat(is:ie,js:je,1) < 1.) )call update_mask_2d
                else
                    ! 3d
                    if(.not. any(rmat(i,j,:) == 1.))cycle ! cycle if equal to one
                    do k=1,self%ldim(3)

                        if( rmat(i,j,k) /= 1. )cycle
                        ! within mask region
                        ks = max(1,k-1)
                        ke = min(k+1,self%ldim(3))
                        if( any( rmat(is:ie,js:je,ks:ke) < 1.) )then
                            ! update since has a masked neighbour
                            kl = max(1,k-falloff)
                            kr = min(k+falloff,self%ldim(3))
                            call update_mask_3d
                        endif
                    end do
                endif
            end do
        end do
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = rmat
        deallocate(rmat)
    contains

        !> updates neighbours with cosine weight
        subroutine update_mask_2d
            integer :: ii, jj, di_sq, dist_sq
            do ii=il,ir
                di_sq = (ii-i)**2                 ! 1D squared distance in x dim
                do jj=jl,jr
                    dist_sq = di_sq + (jj-j)**2   ! 2D squared distance in x & y dim
                    if(dist_sq > falloff_sq)cycle
                    ! masked neighbour
                    if( rmat(ii,jj,1)<1. )&
                        &rmat(ii,jj,1) = max(local_versine(real(dist_sq)), rmat(ii,jj,1))
                enddo
            enddo
        end subroutine update_mask_2d

        !> updates neighbours with cosine weight
        subroutine update_mask_3d
            integer :: ii, jj, kk, di_sq, dij_sq, dist_sq
            do ii=il,ir
                di_sq = (ii-i)**2
                do jj=jl,jr
                    dij_sq = di_sq+(jj-j)**2
                    do kk=kl,kr
                        dist_sq = dij_sq + (kk-k)**2
                        if(dist_sq > falloff_sq)cycle
                        if( rmat(ii,jj,kk)<1. )&
                            &rmat(ii,jj,kk) = max(local_versine(real(dist_sq)), rmat(ii,jj,kk))
                    enddo
                enddo
            enddo
        end subroutine update_mask_3d

        !> Local elemental cosine edge function
        !> this is not a replacement of math%cosedge, which is not applicable here
        elemental real function local_versine( r_sq )result( c )
            real, intent(in) :: r_sq
            c = 0.5 * (1. - cos(scalefactor*(sqrt(r_sq)-rfalloff)) )
        end function local_versine

    end subroutine cos_edge

    !>  \brief  remove edge from binary image
    subroutine remove_edge( self )
        class(image), intent(inout) :: self
        if( self%ft ) call simple_stop('only for real binary images (not FTed ones); simple_image :: remove_edge')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            call simple_stop('input to remove edge not binary; simple_image :: remove_edge')
        where( self%rmat < 0.999 ) self%rmat = 0.
    end subroutine remove_edge


    !>  \brief  increments the logi pixel value with incr
    !! \param logi coordinates
    !! \param incr increment
    !!
    subroutine increment( self, logi, incr )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: incr
        self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))+incr
    end subroutine increment

    !>  \brief  generates a logical mask from a binary one
    function bin2logical( self ) result( mask )
        class(image), intent(in) :: self
        logical, allocatable :: mask(:,:,:)
        allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)))
        where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.0 )
            mask = .true.
        else where
            mask = .false.
        end where
    end function bin2logical

    ! performs left/right collage of input images
    subroutine collage( self1, self2, img_out )
        class(image), intent(inout) :: self1, self2, img_out
        type(image) :: img_pad
        integer     :: ldim(3), ldim_col(3), border
        if( .not.self1%is_2d() )stop '2D only; simple_mage::collage'
        if( self1%is_ft() )stop 'Real space only; simple_mage::collage'
        if( .not.self2%is_2d() )stop '2D only; simple_mage::collage'
        if( self2%is_ft() )stop 'Real space only; simple_mage::collage'
        border   = 2
        ldim(1)  = max(self1%ldim(1),self2%ldim(1))
        ldim(2)  = max(self1%ldim(2),self2%ldim(2))
        ldim(1)  = max(ldim(1), ldim(2))
        ldim(2)  = ldim(1)
        ldim(3)  = 1
        ldim_col = [2*ldim(1)+border, ldim(2), 1]
        call img_out%new(ldim_col,1.)
        ! pad & copy left image
        call img_pad%new(ldim,self1%get_smpd())
        call self1%pad(img_pad)
        call img_pad%norm
        img_out%rmat(:ldim(1),:ldim(2),1) = img_pad%rmat(:ldim(1),:ldim(2),1)
        ! pad & copy right image
        img_pad%rmat = 0.
        call img_pad%set_smpd(self2%get_smpd())
        call self2%pad(img_pad)
        call img_pad%norm
        img_out%rmat(ldim(1)+border+1:ldim_col(1),:ldim_col(2),1) = img_pad%rmat(:ldim(1),:ldim(2),1)
        call img_pad%kill
    end subroutine collage

    ! FILTERS

    !>  \brief  acf calculates the autocorrelation function of an image
    !!
    subroutine acf( self )
        class(image), intent(inout) :: self
        if( .not. self%ft )then
            call self%fft()
        endif
        self%cmat = self%cmat*conjg(self%cmat)
        call self%ifft()
    end subroutine acf

    !>  \brief ccf calculates the cross-correlation function between two images
    !! \param self1 image object
    !! \param self2 image object
    !! \return  cc
    !!  calculates thecross-correlation function between two images
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
    !! \param hp high-pass
    !! \param lp low-pass
    !! \return  bfac
    !!
    function guinier_bfac( self, hp, lp ) result( bfac )
        class(image), intent(inout) :: self
        real, intent(in)            :: hp, lp
        real, allocatable           :: plot(:,:)
        integer                     :: fromk, tok, nk
        real                        :: slope, intercept, corr, bfac
        plot  = self%guinier()
        fromk = self%get_find(hp)
        tok   = self%get_find(lp)
        nk    = tok-fromk+1
        call fit_straight_line(nk, plot(fromk:tok,:), slope, intercept, corr)
        bfac=4.*slope
        deallocate(plot)
    end function guinier_bfac

    !>  \brief guinier generates the Guinier plot for a volume, which should be unfiltered
    !! \return  plot
    !!
    function guinier( self ) result( plot )
        class(image), intent(inout) :: self
        real, allocatable :: spec(:), plot(:,:)
        integer           :: lfny, k
        if( .not. self%is_3d() ) call simple_stop('Only for 3D images; guinier; simple_image')
        call self%spectrum('absreal',spec=spec)
        lfny = self%get_lfny(1)
        allocate( plot(lfny,2), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk("In: guinier; simple_image")
        do k=1,lfny
            plot(k,1) = 1./(self%get_lp(k)**2.)
            plot(k,2) = log(spec(k))
            write(*,'(A,1X,F8.4,1X,A,1X,F7.3)') '>>> RECIPROCAL SQUARE RES:', plot(k,1), '>>> LOG(ABS(REAL(F))):', plot(k,2)
        end do
        deallocate(spec)
    end function guinier

    !>  \brief spectrum generates the rotationally averaged spectrum of an image
    !! \param which accepts 'real''power''absreal''absimag''abs''phase''count'
    !! \param norm normalise result
    !! \return spec Power spectrum array
    !!
    subroutine spectrum( self, which, spec, norm)
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real, allocatable, intent(inout) :: spec(:)
        logical, optional, intent(in)    :: norm
        real, allocatable :: counts(:)
        integer :: lfny, h, k, l
        integer :: sh, lims(3,2), phys(3)
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
        DebugPrint "In simple_image::spectrum lfny", lfny
        if(allocated(spec))deallocate(spec)
        allocate( spec(lfny), counts(lfny), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('spectrum; simple_image',alloc_stat)
        DebugPrint "In simple_image::spectrum spec size", size(spec)
        spec   = 0.
        counts = 0.
        lims   = self%fit%loop_lims(2)
        select case(which)
        case('real')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + real(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('power')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + csq(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('sqrt')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + sqrt(csq(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('log')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + log(csq(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('absreal')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + abs(real(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('absimag')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + abs(aimag(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('abs')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + cabs(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('phase')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + phase_angle(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('count')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sh = nint(hyp(real(h),real(k),real(l)))
                        if( sh == 0 .or. sh > lfny ) cycle
                        spec(sh) = spec(sh) + 1.
                        counts(sh) = counts(sh)+1.
                    end do
                end do
            end do
        case DEFAULT
            write(*,*) 'Spectrum kind: ', trim(which)
            call simple_stop('Unsupported spectrum kind; simple_image; spectrum', __FILENAME__,__LINE__)
        end select
        if( which .ne. 'count' .and. nnorm )then
            where(counts > 0.)
                spec = spec/counts
            end where
        endif
        if( didft ) call self%ifft()
        DebugPrint "In simple_image::spectrum done"
    end subroutine spectrum

    !> \brief shellnorm for normalising each shell to uniform (=1) power
    !!
    subroutine shellnorm( self )
        class(image), intent(inout) :: self
        real, allocatable  :: expec_pow(:)
        logical            :: didbwdft
        integer            :: sh, h, k, l, phys(3), lfny, lims(3,2)
        real               :: icomp, avg
        ! subtract average in real space
        didbwdft = .false.
        if( self%ft )then
            call self%ifft()
            didbwdft = .true.
        endif
        avg = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))/real(product(self%ldim))
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))-avg
        call self%fft()
        lfny  = self%get_lfny(1)
        lims  = self%fit%loop_lims(2)
        ! calculate the expectation value of the signal power in each shell
        !$ allocate(expec_pow(1))
        call self%spectrum('power',expec_pow)
        ! normalise
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh = nint(hyp(real(h),real(k),real(l)))
                    phys = self%fit%comp_addr_phys([h,k,l])
                    if( sh > lfny )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh == 0 ) cycle
                        if( expec_pow(sh) > 0. )then
                            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/sqrt(expec_pow(sh))
                        endif
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        ! take care of the central spot
        phys  = self%fit%comp_addr_phys([0,0,0])
        icomp = aimag(self%cmat(phys(1),phys(2),phys(3)))
        self%cmat(phys(1),phys(2),phys(3)) = cmplx(1.,icomp)
        ! Fourier plan upon return
        if( didbwdft )then
            ! return in Fourier space
        else
            ! return in real space
            call self%ifft()
        endif
    end subroutine shellnorm

    !> \brief  for normalising each shell to uniform (=1) power (assuming average has been
    !!         subtracted in real-space) and applying a filter function
    subroutine shellnorm_and_apply_filter_serial_1( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        real, allocatable  :: expec_pow(:)
        integer            :: sh, h, k, l, phys(3), lfny, lims(3,2), filtsz
        real               :: icomp, wzero, fwght
        filtsz    = size(filter)
        wzero     = maxval(filter)
        lfny      = self%get_lfny(1)
        lims      = self%fit%loop_lims(2)
        ! calculate the expectation value of the signal power in each shell
        call self%spectrum('power',expec_pow )
        ! normalise
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(real(h),real(k),real(l)))
                    phys = self%fit%comp_addr_phys([h,k,l])
                    ! set filter weight
                    if( sh > filtsz )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    if( sh > lfny )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh == 0 ) cycle
                        if( expec_pow(sh) > 0. )then
                            self%cmat(phys(1),phys(2),phys(3)) =&
                                &fwght * (self%cmat(phys(1),phys(2),phys(3)) / sqrt(expec_pow(sh)))
                        endif
                    endif
                end do
            end do
        end do
        ! take care of the central spot
        phys  = self%fit%comp_addr_phys([0,0,0])
        icomp = aimag(self%cmat(phys(1),phys(2),phys(3)))
        self%cmat(phys(1),phys(2),phys(3)) = cmplx(wzero,icomp)
    end subroutine shellnorm_and_apply_filter_serial_1

    !> \brief  for normalising each shell to uniform (=1) power (assuming average has been
    !!         subtracted in real-space) and applying a filter function
    subroutine shellnorm_and_apply_filter_serial_2( self, filter )
        class(image), intent(inout) :: self, filter
        real, allocatable  :: expec_pow(:)
        complex            :: comp
        integer            :: sh, h, k, l, phys(3), lfny, lims(3,2)
        real               :: icomp, fwght
        lfny      = self%get_lfny(1)
        lims      = self%fit%loop_lims(2)
        ! calculate the expectation value of the signal power in each shell
        call self%spectrum('poer',expec_pow)
        ! normalise
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(real(h),real(k),real(l)))
                    phys = self%fit%comp_addr_phys([h,k,l])
                    comp  = filter%get_fcomp([h,k,l],phys)
                    fwght = real(comp)
                    if( sh > lfny )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh == 0 ) cycle
                        if( expec_pow(sh) > 0. )then
                            self%cmat(phys(1),phys(2),phys(3)) =&
                                &fwght * (self%cmat(phys(1),phys(2),phys(3)) / sqrt(expec_pow(sh)))
                        endif
                    endif
                end do
            end do
        end do
        ! take care of the central spot
        phys  = self%fit%comp_addr_phys([0,0,0])
        comp  = filter%get_fcomp([0,0,0],phys)
        icomp = aimag(self%cmat(phys(1),phys(2),phys(3)))
        self%cmat(phys(1),phys(2),phys(3)) = cmplx(real(comp),icomp)
    end subroutine shellnorm_and_apply_filter_serial_2

    !> \brief  for normalising each shell to uniform (=1) power (assuming average has been
    !!         subtracted in real-space) and applying a filter function
    subroutine shellnorm_and_apply_filter_1( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        real, allocatable  :: expec_pow(:)
        integer            :: sh, h, k, l, phys(3), lfny, lims(3,2), filtsz
        real               :: icomp, wzero, fwght
        filtsz    = size(filter)
        wzero     = maxval(filter)
        lfny      = self%get_lfny(1)
        lims      = self%fit%loop_lims(2)
        ! calculate the expectation value of the signal power in each shell
        call self%spectrum('power', expec_pow )
        ! normalise
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(real(h),real(k),real(l)))
                    phys = self%fit%comp_addr_phys([h,k,l])
                    ! set filter weight
                    if( sh > filtsz )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    if( sh > lfny )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh == 0 ) cycle
                        if( expec_pow(sh) > 0. )then
                            self%cmat(phys(1),phys(2),phys(3)) =&
                                &fwght * (self%cmat(phys(1),phys(2),phys(3)) / sqrt(expec_pow(sh)))
                        endif
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        ! take care of the central spot
        phys  = self%fit%comp_addr_phys([0,0,0])
        icomp = aimag(self%cmat(phys(1),phys(2),phys(3)))
        self%cmat(phys(1),phys(2),phys(3)) = cmplx(wzero,icomp)
    end subroutine shellnorm_and_apply_filter_1

    !> \brief  for normalising each shell to uniform (=1) power (assuming average has been
    !!         subtracted in real-space) and applying a filter function
    subroutine shellnorm_and_apply_filter_2( self, filter )
        class(image), intent(inout) :: self, filter
        real, allocatable  :: expec_pow(:)
        complex            :: comp
        integer            :: sh, h, k, l, phys(3), lfny, lims(3,2)
        real               :: icomp, fwght
        lfny      = self%get_lfny(1)
        lims      = self%fit%loop_lims(2)
        ! calculate the expectation value of the signal power in each shell
        call self%spectrum('power',expec_pow )
        ! normalise
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys,comp,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(real(h),real(k),real(l)))
                    phys = self%fit%comp_addr_phys([h,k,l])
                    comp  = filter%get_fcomp([h,k,l],phys)
                    fwght = real(comp)
                    if( sh > lfny )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh == 0 ) cycle
                        if( expec_pow(sh) > 0. )then
                            self%cmat(phys(1),phys(2),phys(3)) =&
                                &fwght * (self%cmat(phys(1),phys(2),phys(3)) / sqrt(expec_pow(sh)))
                        endif
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        ! take care of the central spot
        phys  = self%fit%comp_addr_phys([0,0,0])
        comp  = filter%get_fcomp([0,0,0],phys)
        icomp = aimag(self%cmat(phys(1),phys(2),phys(3)))
        self%cmat(phys(1),phys(2),phys(3)) = cmplx(real(comp),icomp)
    end subroutine shellnorm_and_apply_filter_2

    !> \brief apply_bfac  is for applying bfactor to an image
    !! \param b
    !!
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
        if( didft ) call self%ifft()
    end subroutine apply_bfac

    !> \brief bp  is for band-pass filtering an image
    !! \param hplim
    !! \param lplim
    !! \param width
    !!
    subroutine bp( self, hplim, lplim, width )
        class(image), intent(inout) :: self
        real, intent(in)            :: hplim, lplim
        real, intent(in), optional  :: width
        integer                     :: h, k, l, lims(3,2), phys(3)
        logical                     :: didft, dohp, dolp
        real                        :: freq, hplim_freq, lplim_freq, wwidth, w
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
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    freq = hyp(real(h),real(k),real(l))
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
        if( didft ) call self%ifft()
    end subroutine bp

    !>  \brief gen_lpfilt is for generating low-pass filter weights
    !! \param lplim
    !! \param width
    !! \return  filter array
    !!
    function gen_lpfilt( self, lplim, width ) result( filter )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lplim
        real, optional, intent(in)    :: width
        integer                       :: nyq, k
        real                          :: wwidth, lplim_freq, freq
        real, allocatable             :: filter(:)
        wwidth = 5.
        if( present(width) ) wwidth = width
        nyq = self%get_nyq()
        allocate( filter(nyq), stat=alloc_stat )
        lplim_freq = self%fit%get_find(1,lplim)
        filter = 1.
        do k=1,nyq
            freq = real(k)
            if(freq .gt. lplim_freq)then
                filter(k) = 0.
            else if(freq .ge. lplim_freq-wwidth)then
                filter(k) = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
            endif
        end do
    end function gen_lpfilt

    !> \brief apply_filter_1  is for application of an arbitrary 1D filter function
    !! \param filter
    subroutine apply_filter_1( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        integer                     :: nyq, sh, h, k, l, lims(3,2)
        logical                     :: didft
        real                        :: fwght, wzero
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
                    sh = nint(hyp(real(h),real(k),real(l)))
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

    !> \brief apply_filter_2  is for application of an arbitrary filter function
    !! \param filter
    !!
    subroutine apply_filter_2( self, filter )
        class(image), intent(inout) :: self, filter
        real    :: fwght
        integer :: phys(3), lims(3,2), h, k, l
        complex :: comp
        if( self.eqdims.filter )then
            if( filter%ft )then
                if( .not. self%ft )then
                    call simple_stop('ERROR, image to be filtered is not in the Fourier domain; apply_filter_2; simple_image')
                endif
                lims = self%fit%loop_lims(2)
                !$omp parallel do collapse(3) default(shared) private(h,k,l,comp,fwght,phys)&
                !$omp schedule(static) proc_bind(close)
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        do l=lims(3,1),lims(3,2)
                            phys  = self%comp_addr_phys([h,k,l])
                            comp  = filter%get_fcomp([h,k,l],phys)
                            fwght = real(comp)
                            call self%mul([h,k,l],fwght,phys_in=phys)
                        end do
                    end do
                end do
                !$omp end parallel do
            else
                call simple_stop('assumed that the inputted filter is in the Fourier domain; apply_filter_2; simple_image')
            endif
        else
            call simple_stop('equal dims assumed; apply_filter_2; simple_image')
        endif
    end subroutine apply_filter_2

    !> \brief apply_filter_1  is for application of an arbitrary 1D filter function
    !! \param filter
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
                    sh = nint(hyp(real(h),real(k),real(l)))
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

    !> \brief phase_rand  is for randomzing the phases of the FT of an image from lp and out
    !! \param lp
    !!
    subroutine phase_rand( self, lp )
        class(image), intent(inout) :: self
        real, intent(in)            :: lp
        integer                     :: h,k,l,phys(3),lims(3,2)
        logical                     :: didft
        real                        :: freq,lp_freq,sgn1,sgn2,sgn3
        real, parameter             :: errfrac=0.5
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lp_freq = self%fit%get_find(1,lp) ! assuming square 4 now
        lims    = self%fit%loop_lims(2)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    freq = hyp(real(h),real(k),real(l))
                    if(freq .gt. lp_freq)then
                        phys = self%fit%comp_addr_phys([h,k,l])
                        sgn1 = 1.
                        sgn2 = 1.
                        sgn3 = 1.
                        if( ran3() > 0.5 ) sgn1 = -1.
                        if( ran3() > 0.5 ) sgn2 = -1.
                        if( ran3() > 0.5 ) sgn3 = -1.
                        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
                            self%oshift([h,k,l],[sgn1*ran3()*errfrac*self%ldim(1),&
                            sgn2*ran3()*errfrac*self%ldim(2),sgn3*ran3()*errfrac*self%ldim(3)])
                    endif
                end do
            end do
        end do
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
        allocate( w(kmax), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk("In: hannw; simple_image")
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
        integer               :: n, i, j, k, cnt
        real                  :: rn, wfun(-winsz:winsz), norm
        type(winfuns)         :: fwin
        character(len=STDLEN) :: wstr
        type(image)           :: img_filt
        ! check the number of pixels in window
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
            do i=-winsz,winsz
                do j=-winsz,winsz
                    do k=-winsz,winsz
                        cnt = cnt + 1
                        wfvals(cnt) = wfun(i) * wfun(j) * wfun(k)
                        norm = norm + wfvals(cnt)
                    end do
                end do
            end do
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
            case('bman')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = sum(pixels * wfvals) / norm
                    end do
                end do
                !$omp end parallel do
            case DEFAULT
                call simple_stop('unknown filter type; simple_image :: real_space_filter')
            end select
        else
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
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = sum(pixels)/rn
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
            case DEFAULT
                call simple_stop('unknown filter type; simple_image :: real_space_filter')
            end select
        endif
        call self%copy(img_filt)
        call img_filt%kill()
    end subroutine real_space_filter

    ! CALCULATORS

    !> \brief square_root  is for calculating the square root of an image
    !!
    subroutine square_root( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            where(real(self%cmat) > 0. )
                self%cmat = sqrt(real(self%cmat))
            end where
        else
            where(self%rmat > 0. )
                self%rmat = sqrt(self%rmat)
            end where
        endif
    end subroutine square_root

    !>  \brief maxcoord is for providing location of the maximum pixel value
    function maxcoord(self) result(loc)
        class(image), intent(inout) :: self
        integer :: loc(3)
        if( self%ft )then
            call simple_stop('maxloc not implemented 4 FTs! simple_image')
        else
            loc = maxloc(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        endif
    end function maxcoord

    !> \brief for finding the offset of the cross-correlation peak
    function ccpeak_offset( self ) result( xyz )
        class(image), intent(inout) :: self
        integer :: loc(3)
        real    :: xyz(3)
        if( self%ft ) call simple_stop('not implemented 4 FTs! simple_image :: ccpeak_offset')
        loc    = maxloc(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        xyz(1) = real(loc(1) - self%ldim(1)/2)
        xyz(2) = real(loc(2) - self%ldim(2)/2)
        if( self%ldim(3) > 1 )then
            xyz(3) = real(loc(3) - self%ldim(3)/2)
        else
            xyz(3) = 0.0
        endif
    end function ccpeak_offset

    !> \brief stats  is for providing foreground/background statistics
    !! \param which foreground or background
    !! \param ave Geometric mean
    !! \param sdev Standard Deviation
    !! \param var Variance
    !! \param msk optional input mask
    !! \param med median
    !! \param errout error flag
    !!
    subroutine stats_1( self, which, ave, sdev, maxv, minv, msk, med, errout )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: ave, sdev, maxv, minv
        real,    optional, intent(in)    :: msk
        real,    optional, intent(out)   :: med
        logical, optional, intent(out)   :: errout
        integer           :: i, j, k, npix, alloc_stat, minlen
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
            call simple_stop('unrecognized parameter: which; stats_1; simple_image', __FILENAME__,__LINE__)
        endif
        allocate( pixels(product(self%ldim)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('backgr; simple_image')
        pixels = 0.
        npix = 0
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
            if( err ) write(*,'(a)') 'WARNING: variance zero; stats_1; simple_image'
        endif
        if( didft ) call self%fft()
    end subroutine stats_1

    !> \brief stats  is for providing within mask statistics
    !! \param mskimg binary mask image
    !! \param ave Geometric mean of pixels within mask
    !! \param sdev Standard Deviation of pixels within mask
    !! \param maxv maximum value within mask
    !! \param minv minimum value within mask
    !! \param m optional input mask
    !! \param med median
    !! \param errout error flag
    !!
    subroutine stats_2( self, mskimg, ave, sdev, maxv, minv, med, errout )
        class(image),      intent(inout) :: self
        class(image),      intent(in)    :: mskimg
        real,              intent(out)   :: ave, sdev, maxv, minv
        real,    optional, intent(out)   :: med
        logical, optional, intent(out)   :: errout
        real              :: var
        logical           :: err
        real, allocatable :: pixels(:)
        ! FT
        if( self%ft ) stop 'not for FTed imgs; simple_image :: stats_2'
        pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
            &mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.95 )
        maxv = maxval(pixels)
        minv = minval(pixels)
        call moment( pixels, ave, sdev, var, err )
        if( present(med) ) med  = median_nocopy(pixels)
        deallocate( pixels )
        if( present(errout) )then
            errout = err
        else
            if( err ) write(*,'(a)') 'WARNING: variance zero; stats_2; simple_image'
        endif
    end subroutine stats_2

    !>  \brief minmax to get the minimum and maximum values in an image
    !! \return  mm 2D element (minimum , maximum)
    function minmax( self )result( mm )
        class(image), intent(in) :: self
        real :: mm(2)
        mm(1) = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        mm(2) = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
    end function minmax

    !>  \brief rmsd for calculating the RMSD of a map
    !! \return  dev root mean squared deviation
    !!
    function rmsd( self ) result( dev )
        class(image), intent(inout) :: self
        real :: devmat(self%ldim(1),self%ldim(2),self%ldim(3)), dev, avg
        if( self%ft )then
            dev = 0.
        else
            avg    = self%mean()
            devmat = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - avg
            dev    = sum(devmat**2.0)/real(product(self%ldim))
            if( dev > 0. )then
                dev = sqrt(dev)
            else
                dev = 0.
            endif
        endif
    end function rmsd

    !> \brief noisesdev is for estimating the noise variance of an image
    !!          by online estimation of the variance of the background pixels
    !>
    !! \param msk mask threshold
    !! \return  sdev
    !!
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
        call ovar%finalize
        mv = ovar%get_mean_var()
        sdev = 0.
        if( mv(2) > 0. ) sdev = sqrt(mv(2))
        if( didft ) call self%fft()
    end function noisesdev

    !>  \brief  is for calculating the mean of an image
    !> mean
    !! \return  avg
    !!
    function mean( self ) result( avg )
        class(image), intent(inout) :: self
        real :: avg
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
    !!
    subroutine checkimg4nans( self )
        class(image), intent(in) :: self
        if( self%ft )then
            call check4nans3D(self%cmat)
        else
            call check4nans3D(self%rmat)
        endif
    end subroutine checkimg4nans

    !> \brief cure_2  is for checking the numerical soundness of an image and curing it if necessary
    !! \param maxv
    !! \param minv
    !! \param ave
    !! \param sdev
    !! \param n_nans
    !!
    subroutine cure( self, maxv, minv, ave, sdev, n_nans )
        class(image), intent(inout) :: self
        real,         intent(out)   :: maxv, minv, ave, sdev
        integer,      intent(out)   :: n_nans
        integer                     :: i, j, k, npix
        real                        :: var, ep, dev
        if( self%ft )then
            write(*,*) 'WARNING: Cannot cure FTs; cure; simple_image'
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
            write(*,*) 'found NaNs in simple_image; cure:', n_nans
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
    !! \param mode
    !! \param lp_dyn
    !! \return  lims
    !!
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

    !>  \brief  Convert logical address to physical address. Complex image.
    pure function comp_addr_phys(self,logi) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3) !<  Logical address
        integer                  :: phys(3) !<  Physical address
        phys = self%fit%comp_addr_phys(logi)
    end function comp_addr_phys

    !>  \brief  generate physical index mapping array
    subroutine get_2Dphys_ind_mapping( self, lims, ind_map )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: lims(2,2)
        integer,      intent(out) :: ind_map(lims(1,1):lims(1,2),lims(2,1):lims(2,2),3)
        integer :: h, k
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                ind_map(h,k,:) = self%fit%comp_addr_phys([h,k,0])
            end do
        end do
    end subroutine get_2Dphys_ind_mapping

    !>  \brief corr is for correlating two images
    !! \param self1 input image 1
    !! \param self2 input image 2
    !! \param lp_dyn low-pass cut-off freq
    !! \param hp_dyn high-pass cut-off freq
    !! \return  r Correlation coefficient
    !!
    function corr( self1, self2, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self1, self2
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real    :: r, sumasq, sumbsq
        integer :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        logical :: didft1, didft2
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
            r = 0.
            sumasq = 0.
            sumbsq = 0.
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
            if( sumasq < TINY .or. sumbsq < TINY )then
                r = 0.
            else
                r = r / sqrt(sumasq * sumbsq)
            endif
            if( didft1 ) call self1%ifft()
            if( didft2 ) call self2%ifft()
        else
            write(*,*) 'self1%ldim:', self1%ldim
            write(*,*) 'self2%ldim:', self2%ldim
            call simple_stop('images to be correlated need to have same dimensions; corr; simple_image')
        endif
    end function corr

    !> \brief is for highly optimized correlation between 2D images, particle is
    !> shifted by shvec so remember to take care of this properly in the calling
    !> module corr_shifted
    !! \param self_ref reference image
    !! \param self_ptcl particle image object
    !! \param shvec shift vector
    !! \param lp_dyn  low-pass
    !! \param hp_dyn  high-pass
    !! \return  r correlation coefficient
    !!
    function corr_shifted( self_ref, self_ptcl, shvec, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self_ref, self_ptcl
        real,           intent(in)    :: shvec(3)
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real                          :: r, sumasq, sumbsq
        complex                       :: shcomp
        integer                       :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        ! this is for highly optimised code, so we assume that images are always Fourier transformed beforehand
        if( .not. self_ref%ft  ) call simple_stop('self_ref not FTed;  corr_shifted; simple_image')
        if( .not. self_ptcl%ft ) call simple_stop('self_ptcl not FTed; corr_shifted; simple_image')
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
                        phys = self_ref%fit%comp_addr_phys([h,k,l])
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
    !!
    function real_corr_1( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3))
        real :: r, ax, ay, sxx, syy, sxy, npix
        npix  = real(product(self1%ldim))
        ax    = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3))) / npix
        ay    = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3))) / npix
        diff1 = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
        diff2 = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
        sxx   = sum(diff1 * diff1)
        syy   = sum(diff2 * diff2)
        sxy   = sum(diff1 * diff2)
        if( sxx > 0. .and. syy > 0. )then
            r = sxy / sqrt(sxx * syy)
        else
            r = 0.
        endif
    end function real_corr_1

    !>  \brief real_corr_2 is for calculating a real-space correlation coefficient between images within a mask
    !>  Input mask is assumed binarized
    !! \param self1 image object
    !! \param self2 image object
    !! \param maskimg
    !! \return  r
    !!
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
        if( sxx > 0. .and. syy > 0. )then
            r = sxy / sqrt(sxx * syy)
        else
            r = 0.
        endif
    end function real_corr_2

    !> \brief prenorm4real_corr pre-normalises the reference in preparation for real_corr_prenorm
    !! \param sxx
    !!
    subroutine prenorm4real_corr_1( self, sxx )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        real :: diff(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: npix, ax
        npix = real(product(self%ldim))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / npix
        diff = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ax
        sxx  = sum(diff * diff)
    end subroutine prenorm4real_corr_1

    !> \brief prenorm4real_corr pre-normalises the reference in preparation for real_corr_prenorm
    !! \param sxx
    !!
    subroutine prenorm4real_corr_2( self, sxx, mask )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        logical,      intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: diff(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: npix, ax
        npix = real(count(mask))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask) / npix
        where( mask ) diff = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ax
        sxx  = sum(diff * diff, mask=mask)
    end subroutine prenorm4real_corr_2

    !>  \brief real_corr_prenorm is for calculating a real-space correlation coefficient between images (reference is pre-normalised)
    !! \param self_ref image object
    !! \param self_ptcl image object
    !! \param sxx_ref
    !! \return  r
    !!
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
    !! \param self_ref image object
    !! \param self_ptcl image object
    !! \param sxx_ref
    !! \return  r
    !!
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

    !> \brief fsc is for calculation of Fourier ring/shell correlation
    !! \param self1 image object
    !! \param self2 image object
    !! \param res
    !! \param corrs
    !!
    subroutine fsc( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
        real    :: sumasq(fdim(self1%ldim(1))-1), sumbsq(fdim(self1%ldim(1))-1)
        integer :: n, lims(3,2), phys(3), sh, h, k, l
        corrs  = 0.
        sumasq = 0.
        sumbsq = 0.
        lims   = self1%fit%loop_lims(2)
        n      = self1%get_filtsz()
        do k=lims(2,1),lims(2,2)
            do h=lims(1,1),lims(1,2)
                do l=lims(3,1),lims(3,2)
                    ! compute physical address
                    phys = self1%fit%comp_addr_phys([h,k,l])
                    ! find shell
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if( sh == 0 .or. sh > n ) cycle
                    ! real part of the complex mult btw self1 and targ*
                    corrs(sh) = corrs(sh)+&
                        real(self1%cmat(phys(1),phys(2),phys(3))*conjg(self2%cmat(phys(1),phys(2),phys(3))))
                    sumasq(sh) = sumasq(sh) + csq(self2%cmat(phys(1),phys(2),phys(3)))
                    sumbsq(sh) = sumbsq(sh) + csq(self1%cmat(phys(1),phys(2),phys(3)))
                end do
            end do
        end do
        ! normalize correlations and compute resolutions
        do k=1,n
            if( sumasq(k) > 0. .and. sumbsq(k) > 0. )then
                corrs(k) = corrs(k)/sqrt(sumasq(k) * sumbsq(k))
            else
                corrs(k) = 0.
            endif
        end do
    end subroutine fsc

    !> \brief get_nvoxshell is for calculation of voxels per Fourier shell
    !! \param voxs
    !!
    subroutine get_nvoxshell( self, voxs )
        class(image)     , intent(inout) :: self
        real, allocatable, intent(inout) :: voxs(:)
        integer                          :: n, lims(3,2), sh, h, k, l
        logical                          :: didft
        if( .not. square_dims(self) ) call simple_stop('square dimensions only! fsc; simple_image')
        didft = .false.
        if( .not. self%is_ft() )then
            call self%fft()
            didft = .true.
        endif
        n = self%get_filtsz()
        if( allocated(voxs) )deallocate(voxs)
        allocate( voxs(n), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('In: get_nvoxshell, module: simple_image')
        voxs = 0.
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if( sh == 0 .or. sh > n ) cycle
                    voxs( sh ) = voxs( sh )+1.
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine get_nvoxshell

    !>  \brief get array of resolution steps
    !> get_res
    !! \return  res
    !!
    function get_res( self ) result( res )
        class(image), intent(in) :: self
        real, allocatable        :: res(:)
        integer                  :: n, k
        n = self%get_filtsz()
        allocate( res(n), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('In: get_res, module: simple_image')
        do k=1,n
            res(k) = self%fit%get_lp(1,k)
        end do
    end function get_res

    !>  \brief  returns the real and imaginary parts of the phase shift at point
    !!          logi in a Fourier transform caused by the origin shift in shvec
    !> oshift_1
    !! \param logi index in Fourier domain
    !! \param shvec origin shift
    !! \return  comp
    !!
    function oshift_1( self, logi, shvec ) result( comp )
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
    !> oshift_2
    !! \param logi index in Fourier domain
    !! \param shvec origin shift
    !! \return  comp
    !!
    function oshift_2( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
        comp = self%oshift_1(real(logi), shvec)
    end function oshift_2

    !>  \brief  returns the real argument transfer matrix components at point logi in a Fourier transform
    !> gen_argtransf_comp
    !! \param logi index in Fourier domain
    !! \param ldim image dimensions
    !! \return  arg
    !!
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

    !> \brief gen_argtransf_mats  is for generating the argument transfer matrix for fast shifting of a FT
    !! \param transfmats transfer matrix
    !!
    subroutine gen_argtransf_mats( self, transfmats )
        class(image), intent(inout) :: self, transfmats(3)
        integer                     :: h, k, l, lims(3,2), phys(3)
        real                        :: arg(3)
        call transfmats(1)%new(self%ldim,self%smpd)
        call transfmats(2)%new(self%ldim,self%smpd)
        call transfmats(3)%new(self%ldim,self%smpd)
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(arg,phys,h,k,l)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    arg  = self%gen_argtransf_comp(real([h,k,l]))
                    phys = self%fit%comp_addr_phys([h,k,l])
                    transfmats(1)%cmat(phys(1),phys(2),phys(3)) = cmplx(arg(1),0.)
                    transfmats(2)%cmat(phys(1),phys(2),phys(3)) = cmplx(arg(2),0.)
                    transfmats(3)%cmat(phys(1),phys(2),phys(3)) = cmplx(arg(3),0.)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gen_argtransf_mats

    ! MODIFIERS

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
        if( self_in%ldim(3) > 1 )       call simple_stop('only 4 2D images; insert; simple_image')
        if( self_in%is_ft() )           call simple_stop('only 4 real images; insert; simple_image')
        if( .not. self_in%even_dims() ) call simple_stop('only 4 even particle dims; insert; simple_image')
        if( self_out%exists() )then
            if( self_out%ldim(3) > 1 )  call simple_stop('only 4 2D images; insert; simple_image')
            if( self_out%is_ft() )      call simple_stop('only 4 real images; insert; simple_image')
            if( self_out%ldim(1) > self_in%ldim(1) .and. self_out%ldim(2) > self_in%ldim(2) .and. self_out%ldim(3) == 1 )then
                if( (coord(1) < self_in%ldim(1)/2+1 .or. coord(1) > self_out%ldim(1)-self_in%ldim(1)/2-1) .or.&
                    (coord(2) < self_in%ldim(2)/2+1 .or. coord(2) > self_out%ldim(2)-self_in%ldim(2)/2-1) )then
                    call simple_stop('particle outside micrograph area; insert; simple_image')
                endif
            else
                call simple_stop('micrograph needs to have dimensions larger than the particle; insert; simple_image')
            endif
        else
            call simple_stop('micrograph (self_out) does not exist; insert; simple_image')
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
        if( .not. self%ft        ) stop 'image to be modified assumed to be FTed; image :: insert_lowres'
        if( .not. self2insert%ft ) stop 'image to insert assumed to be FTed; image :: insert_lowres'
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys,comp)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(real(h),real(k),real(l)))
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
        if( .not. self%ft        ) stop 'image to be modified assumed to be FTed; image :: insert_lowres'
        if( .not. self2insert%ft ) stop 'image to insert assumed to be FTed; image :: insert_lowres'
        lims = self%fit%loop_lims(2)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(real(h),real(k),real(l)))
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

    !>  \brief  is for making a random image (0,1)
    subroutine ran( self )
        class(image), intent(inout) :: self
        integer :: i, j, k
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    self%rmat(i,j,k) = ran3()
                end do
            end do
        end do
        self%ft = .false.
    end subroutine ran

    !> \brief gauran  is for making a Gaussian random image (0,1)
    !! \param mean Mean of noise
    !! \param sdev Standard deviation of noise
    !!
    subroutine gauran( self, mean, sdev )
        class(image), intent(inout) :: self
        real, intent(in) :: mean, sdev
        integer :: i, j, k
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
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
    subroutine add_gauran( self, snr, noiseimg )
        class(image), intent(inout)        :: self
        real, intent(in)                   :: snr
        type(image), intent(out), optional :: noiseimg
        real    :: noisesdev, ran
        integer :: i, j, k
        logical :: noiseimg_present
        call self%norm
        noiseimg_present = present(noiseimg)
        if( noiseimg_present ) call noiseimg%new(self%ldim, self%smpd)
        noisesdev = sqrt(1/snr)
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    ran = gasdev(0., noisesdev)
                    self%rmat(i,j,k) = self%rmat(i,j,k)+ran
                    if( noiseimg_present ) call noiseimg%set([i,j,k], ran)
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
        ! for the time where we care about modifying the image
        ! class(image), intent(inout) :: self
        ! integer :: k, npixels, start
        ! real, allocatable :: pixels(:)
        ! real    :: med
        ! if( self%ldim(3) == 1 )then
        !     npixels = 2*sum(self%ldim(1:2))
        !     allocate(pixels(npixels))
        !     pixels(1:self%ldim(2))                               = self%rmat(1           , :self%ldim(2),1)
        !     pixels(self%ldim(2)+1:2*self%ldim(2))                = self%rmat(self%ldim(1), :self%ldim(2),1)
        !     pixels(2*self%ldim(2)+1:2*self%ldim(2)+self%ldim(1)) = self%rmat(:self%ldim(1), 1, 1)
        !     pixels(2*self%ldim(2)+self%ldim(1)+1:)               = self%rmat(:self%ldim(1), self%ldim(2),1)
        ! else
        !     npixels = 4*sum(self%ldim(1:3))
        !     allocate(pixels(npixels))
        !     pixels(1:self%ldim(2))                                              = self%rmat(1           , :self%ldim(2),1)
        !     pixels(self%ldim(2)+1:2*self%ldim(2))                               = self%rmat(self%ldim(1), :self%ldim(2),1)
        !     pixels(2*self%ldim(2)+1:2*self%ldim(2)+self%ldim(1))                = self%rmat(:self%ldim(1), 1, 1)
        !     pixels(2*self%ldim(2)+self%ldim(1)+1:2*(self%ldim(2)+self%ldim(1))) = self%rmat(:self%ldim(1), self%ldim(2),1)
        !     start = 2*(self%ldim(2)+self%ldim(1))
        !     pixels(start:start+self%ldim(2))                                 = self%rmat(1           , :self%ldim(2),self%ldim(3))
        !     pixels(start+self%ldim(2)+1:start+2*self%ldim(2))                = self%rmat(self%ldim(1), :self%ldim(2),self%ldim(3))
        !     pixels(start+2*self%ldim(2)+1:start+2*self%ldim(2)+self%ldim(1)) = self%rmat(:self%ldim(1), 1,self%ldim(3))
        !     pixels(start+2*self%ldim(2)+self%ldim(1)+1:2*start)              = self%rmat(:self%ldim(1), self%ldim(2),self%ldim(3))
        !     start = 4*(self%ldim(2)+self%ldim(1))
        !     pixels(start:start+self%ldim(3))                    = self%rmat(1           ,1,           :self%ldim(3))
        !     pixels(start+self%ldim(3)+1:start+2*self%ldim(3))   = self%rmat(1           ,self%ldim(2),:self%ldim(3))
        !     pixels(start+2*self%ldim(3)+1:start+3*self%ldim(3)) = self%rmat(self%ldim(1),1,           :self%ldim(3))
        !     pixels(start+3*self%ldim(3)+1:)                     = self%rmat(self%ldim(1),self%ldim(2),:self%ldim(3))
        ! endif
        ! k = nint(real(npixels)/2.)
        ! med = selec(k,npixels,pixels)
        ! if(abs(med) > TINY) self%rmat = self%rmat - med
        ! deallocate(pixels)
    end subroutine zero_background

    subroutine subtr_backgr_pad_divwinstr_fft( self, mskimg, instr_fun, self_out )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: mskimg
        real,         intent(in)    :: instr_fun(:,:,:)
        class(image), intent(inout) :: self_out
        real, allocatable :: pixels(:)
        integer :: starts(3), stops(3)
        real    :: med
        pixels   = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
            &mask=mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0.5 )
        med = median_nocopy(pixels)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        !$omp parallel workshare proc_bind(close)
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - med
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),1) = &
            &self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),1) / instr_fun(:self_out%ldim(1),:self_out%ldim(2),1)
        !$omp end parallel workshare
        call self_out%fft()
    end subroutine subtr_backgr_pad_divwinstr_fft

    subroutine subtr_backgr_pad_divwinstr_fft_serial( self, mskimg, instr_fun, self_out )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: mskimg
        real,         intent(in)    :: instr_fun(:,:,:)
        class(image), intent(inout) :: self_out
        real, allocatable :: pixels(:)
        integer :: starts(3), stops(3)
        real    :: med
        pixels   = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
            &mask=mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0.5 )
        med = median_nocopy(pixels)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - med
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),1) = &
            &self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),1) / instr_fun(:self_out%ldim(1),:self_out%ldim(2),1)
        call self_out%fft()
    end subroutine subtr_backgr_pad_divwinstr_fft_serial

    subroutine subtr_backgr_pad_divwinstr_serial( self, mskimg, instr_fun, self_out )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: mskimg
        real,         intent(in)    :: instr_fun(:,:,:)
        class(image), intent(inout) :: self_out
        real, allocatable :: pixels(:)
        integer :: starts(3), stops(3)
        real    :: med
        pixels   = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
            &mask=mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0.5 )
        med = median_nocopy(pixels)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - med
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),1) = &
            &self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),1) / instr_fun(:self_out%ldim(1),:self_out%ldim(2),1)
    end subroutine subtr_backgr_pad_divwinstr_serial

    !>  \brief  Taper edges of image so that there are no sharp discontinuities in real space
    !!          This is a re-implementation of the MRC program taperedgek.for (Richard Henderson, 1987)
    !!          I stole it from CTFFIND4 (thanks Alexis for the beautiful re-implementation)
    subroutine taper_edges(self)
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
                smooth_avg_curr_edge_stop (self%ldim(dim2),self%ldim(dim3)),&
                stat=alloc_stat)
            if(alloc_stat/=0)call allocchk("In simple_image::taper_edges avg_curr etc.")
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

    !> \brief salt_n_pepper  is for adding salt and pepper noise to an image
    !! \param pos 2D mask
    !!
    subroutine salt_n_pepper( self, pos )
        class(image), intent(inout) :: self
        logical, intent(in)         :: pos(:,:)
        integer :: ipix, jpix
        if( .not. self%is_2d() ) stop 'only for 2D images; salt_n_pepper; simple_image'
        call self%norm_bin
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
    !!
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
            call simple_stop('image is to small to fit the square; square; simple_image')
        endif
        self%ft = .false.
    end subroutine square

    !>  \brief  just a corner filling fun for testing purposes
    !! \param sqrad half width of square
    !!
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
    !!
    function before_after( left, right, mask ) result( ba )
        class(image),      intent(in) :: left, right
        logical, optional, intent(in) :: mask(left%ldim(1),left%ldim(2),left%ldim(3))
        integer     :: ldim(3), i, j
        type(image) :: ba
        if( left.eqdims.right )then
            if( left.eqsmpd.right )then
                if( left%ft .or. right%ft ) call simple_stop('not for FTs; before_after; simple_image')
                if( left%is_3d() .or. right%is_3d() ) call simple_stop('not for 3D imgs; before_after; simple_image')
                ldim = left%ldim
                ba = left
                ba%rmat(:ldim(1)/2,:ldim(2),1) = left%rmat(:ldim(1)/2,:ldim(2),1)
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
            else
                call simple_stop('before (left) and after (right) not of same smpd; before_after; simple_image')
            endif
        else
            call simple_stop('before (left) and after (right) not of same dim; before_after; simple_image')
        endif
    end function before_after

    !> \brief gauimg  just a Gaussian fun for testing purposes
    !! \param wsz window size
    !!
    subroutine gauimg( self, wsz , alpha )
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
    end subroutine gauimg

    !> \brief gauimg2  just a Gaussian fun for testing purposes
    !! \param wsz window size
    subroutine gauimg2( self, wsz, offx,offy)
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
    end subroutine gauimg2


    subroutine gauimg2D( self, xsigma, ysigma )
        class(image), intent(inout) :: self
        real,         intent(in)    :: xsigma, ysigma
        real    :: x, y
        integer :: i, j
        x = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            y = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                self%rmat(i,j,1) = gaussian2D( [0.,0.], x, y, xsigma, ysigma )
                y = y+1.
            end do
            x = x+1.
        end do
        self%ft = .false.
    end subroutine gauimg2D

    subroutine fwd_ft(self)
        class(image), intent(inout) :: self
        if( self%ft ) return
        if( shift_to_phase_origin ) call self%shift_phorig
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        ! now scale the values so that a ifft() of the output yields the
        ! original image back, rather than a scaled version
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fwd_ft

    !> \brief bwd_ft  backward Fourier transform FFTW3
    !!
    subroutine bwd_ft( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            call fftwf_execute_dft_c2r(self%plan_bwd,self%cmat,self%rmat)
            self%ft = .false.
            if( shift_to_phase_origin ) call self%shift_phorig
        endif
    end subroutine bwd_ft

    !> \brief ft2img  generates images for visualization of a Fourier transform
    subroutine ft2img( self, which, img )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: which
        class(image),     intent(inout) :: img
        integer :: h,mh,k,mk,l,ml,lims(3,2),inds(3),phys(3)
        logical :: didft
        complex :: comp
        if( .not.(self.eqdims.img) )then
            print *, 'self%ldim: ', self%ldim
            print *, 'img%ldim:  ', img%ldim
            stop 'non-equal dims; simple_image :: ft2img'
        endif
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        call img%zero_and_unflag_ft
        lims = self%loop_lims(3)
        mh = maxval(lims(1,:))
        mk = maxval(lims(2,:))
        ml = maxval(lims(3,:))
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    phys = self%comp_addr_phys([h,k,l])
                    comp = self%get_fcomp([h,k,l],phys)
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    inds(2) = min(max(1,k+mk+1),self%ldim(2))
                    inds(3) = min(max(1,l+ml+1),self%ldim(3))
                    select case(which)
                    case ('real')
                        call img%set(inds,real(comp))
                    case('power')
                        call img%set(inds,csq(comp))
                    case('sqrt')
                        call img%set(inds,sqrt(csq(comp)))
                    case ('log')
                        call img%set(inds,log(csq(comp)))
                    case('phase')
                        call img%set(inds,phase_angle(comp))
                    case DEFAULT
                        write(*,*) 'Usupported mode: ', trim(which)
                        stop 'simple_image :: ft2img'
                    end select
                end do
            end do
        end do
        if( didft ) call self%ifft()
    end subroutine ft2img

    subroutine img2ft( self, img )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: img
        integer :: h,k,l,lims(3,2),logi(3),phys(3)
        integer :: xcnt,ycnt,zcnt
        if( .not.(self.eqdims.img) )then
            print *, 'self%ldim: ', self%ldim
            print *, 'img%ldim:  ', img%ldim
            stop 'non-equal dims; simple_image :: img2ft'
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


    subroutine fft_pgi_cuda( self )
#ifndef PGI
        class(image), intent(inout)      :: self
        print *,"simple_image::fft_pgi_cuda cannot be called from this compiler"
#else
        !! IN-PLACE 3D FFT using CUDA's cufft library
        use simple_cufft
#ifdef _OPENACC
        use openacc
#endif
        use cudafor
        class(image), intent(inout)      :: self
        complex, allocatable, device     :: coutput_d(:,:,:)
        complex, allocatable             :: coutput(:,:,:)
        real,    allocatable             :: rinput(:,:,:)
        integer      :: plan, planType
        integer      :: i,j,k,h,l,n,istat
        integer      :: nerrors, cdim(3),ldim(3),lims(3,2),phys(3)
        real         :: nscale
        integer(timer_int_kind) :: t1
        if( self%ft ) return
        ldim = self%ldim
        verbose=.true.
        if(debug)t1=tic()
        ! synchronise
        !        istat=cudaDeviceSynchronize()
        !        if(istat.ne.0)then
        !            call simple_cuda_stop("In simple_image::fft sync ",__FILENAME__,__LINE__)
        !        endif

        ! Allocate arrays on the host
        allocate(coutput(ldim(1),ldim(2),ldim(3)),stat=istat)
        if(alloc_stat.ne.0)call allocchk("In simple_image::fft coutput",istat)
        !$ allocate( rinput(ldim(1), ldim(2), ldim(3)),source=0.,stat=alloc_stat)
        !$ if(alloc_stat.ne.0)call allocchk("In simple_image::pgi_fft rinput",alloc_stat)
        rinput = self%get_rmat()

        ! Allocate arrays on the device
        allocate(coutput_d(ldim(1),ldim(2),ldim(3)))
        ! Set planType to either single or double precision
        planType = CUFFT_C2C    !! (fp_kind == singlePrecision)

        ! Copy arrays to device
        coutput=cmplx(rinput,0.)
        coutput_d = coutput

        ! Initialize the plan for real to complex transform
        call cufftPlan3d(plan,ldim(1),ldim(2),ldim(3),planType)
        DebugPrint "In Image::fft  PGI ", toc()

        ! Execute  Forward transform in place
        call cufftExec(plan,planType,coutput_d,coutput_d,CUFFT_FORWARD )

        ! Scale
        !$acc kernels
        coutput_d = coutput_d / nscale
        !$acc end kernels

        ! Copy results back to host
        coutput = coutput_d


        self%ft =.true.
        self%cmat(:cdim(1),:cdim(2),:cdim(3))=coutput(:cdim(1),:cdim(2),:cdim(3))


        ! Release memory on the host and on the device
        deallocate(coutput_d,coutput,rinput,stat=istat)
        if(alloc_stat /= 0)call allocchk("In simple_image::fft dealloc host coutput",istat)
        ! Destroy the plans
        call cufftDestroy(plan)
        DebugPrint "In simple_cufft::fft_pgi  PGI destroyed ", toc()

        ! DebugPrint "In simple_cufft::fft_pgi  deallocating host arrays", toc()
        ! ! Remove host storage
        ! if(allocated(coutput))then
        !     deallocate(coutput,stat=alloc_stat)
        !     if(alloc_stat /= 0)call allocchk("In simple_image::fft dealloc host coutput",alloc_stat)
        ! end if
        ! if(allocated(rinput))then
        !     deallocate(rinput,stat=alloc_stat)
        !     if(alloc_stat /= 0)call allocchk("In simple_image::fft dealloc host rinput",alloc_stat)
        ! end if
        ! DebugPrint "In simple_cufft::fft_pgi deallocating done "
        ! istat=cudaThreadSynchronize()
        !   VerbosePrint "In simple_cufft::fft_pgi  PGI scaling ", size(coutput), LBOUND(coutput), UBOUND(coutput)

        !istat=istat+cudaDeviceSynchronize()
        !if(istat.ne.0)then
        !    call simple_cuda_stop("In simple_image::fft post fft sync ",__FILENAME__,__LINE__)
        !endif
        DebugPrint "In simple_cufft::fft_pgi finished ", toc(t1)
#endif
    end subroutine fft_pgi_cuda

    subroutine ifft_pgi_cuda( self )
#ifndef PGI
        class(image), intent(inout)      :: self
        print *,"ERROR>>> simple_image::ifft_pgi_cuda cannot be called with compiler ", FC_COMPILER
#else
        !! IN-PLACE 3D FFT using CUDA's cufft
#ifdef _OPENACC
        use openacc
#endif
        use simple_cufft
        use gnufor2
        !        complex,intent(in) ::cin(:,:,:)
        !        integer, intent(in) :: ldim(3)
        class(image), intent(inout)       :: self
        complex, device, allocatable :: inplace_d(:,:,:)
        complex, allocatable         :: cinout(:,:,:)
        ! real, allocatable         :: routput(:,:,:)

        integer      :: i,j,k,f,g,h, xdim(2),ydim(2),zdim(2),ii,jj,kk
        integer      :: nerrors,  ldim(3), cdim(3),lims(3,2),phys(3), cmat_shape(3)
        real         :: nscale, rswap
        complex      :: comp, cswap
        integer      :: planType, plan,istat


        ldim = self%ldim

        cmat_shape = self%array_shape
        cmat_shape(1) = size(self%cmat,1)
        self%ft = .false.

        ! allocate arrays on the host
        !  if (allocated(routput))deallocate(routput)
        if (allocated(cinout))deallocate(cinout)
        !allocate(cin(cmat_shape(1),cmat_shape(2),cmat_shape(3)))
        !cin = self%get_cmatfull()
        !  allocate(routput(ldim(2),ldim(1),1))
        allocate(cinout(ldim(1),ldim(2),ldim(3)),stat=istat)
        if(alloc_stat.ne.0)call allocchk("In simple_image::ifft coutput",istat)
        ! allocate arrays on the device
        allocate(inplace_d(ldim(1),ldim(2),ldim(3)))

        ! if(is_even(n))then
        ! lowerbounds(1) = lbound(cin,1);lowerbounds(2) = lbound(cin,2);
        ! upperbounds(1) = ubound(cin,1); upperbounds(2) = ubound(cin,2)
        !     xbounds = (\ lbound(cin,1), ubound(cin,1) \) ! n/2
        !     ybounds = (\ lbound(cin,2), ubound(cin,2) \) !  m/2
        ! else
        xdim = [ -(ldim(1)-1)/2 , ((ldim(1)-1)/2) - 1  ]
        ydim = [ -(ldim(2)/2)   , (ldim(2)/2) - 1      ]
        zdim = [ -(ldim(3)/2)   , (ldim(3)/2) - 1      ]
        ! endif

        DebugPrint "simple_image::ifft_pgi_cuda CIN    ", xdim, ydim, zdim, cmat_shape
        DebugPrint "simple_image::ifft_pgi_cuda CINPUT ", shape(self%cmat)
        cinout = cmplx(0.,0.)
        lims = self%fit%loop_lims(2)
        DebugPrint "simple_image::ifft_pgi_cuda CMAT LIMS ", lims
        phys=1

        !$omp parallel do collapse(2) default(shared) private(i,j,k,f,g,h,phys,comp)
        do j=1,ldim(2)
            do k=1,ldim(3)
                g=ydim(1)+(j-1)
                f=zdim(1)+(k-1)
                phys=1
                i=1
                do h= xdim(1),-1
                    phys(1) = -h + 1
                    phys(2) = -g + 1 + MERGE(ldim(2),0, -g < 0)
                    phys(3) = -f + 1 + MERGE(ldim(3),0, -f < 0)
#ifdef _DEBUG
                    if (phys(1) < 1 .or. phys(1) > cmat_shape(1) ) print *," CIN ind 1 err ", i, h, phys(1)
                    if (phys(2)< 1 .or. phys(2) > cmat_shape(2)) print *," CIN ind 2 err ", j, g, phys(2)
                    if (phys(3)< 1 .or. phys(3) > cmat_shape(3)) print *," CIN ind 3 err ", k,f, phys(3)

                    if (i < 1 .or. i > ldim(1)) print *," CINOUT ind 1 err ", i, h
                    if (j < 1 .or. j > ldim(2)) print *," CINOUT ind 2 err ", j, g
                    if (k < 1 .or. k > ldim(3)) print *," CINOUT ind 3 err ", k, f
#endif

                    cinout(i,j,k) = conjg(self%cmat(phys(1),phys(2),phys(3)))
                    i=i+1
                end do
                do h=0,xdim(2)
                    phys(1) = h + 1
                    phys(2) = g + 1 + MERGE(ldim(2),0, g < 0)
                    phys(3) = f + 1 + MERGE(ldim(3),0, f < 0)
#ifdef _DEBUG
                    if (phys(1) < 1 .or. phys(1) > cmat_shape(1) ) print *," CIN ind 1 err ", i, h, phys(1)
                    if (phys(2)< 1 .or. phys(2) > cmat_shape(2)) print *," CIN ind 2 err ", j, g, phys(2)
                    if (phys(3)< 1 .or. phys(3) > cmat_shape(3)) print *," CIN ind 3 err ", k,f, phys(3)
                    if (i < 1 .or. i > ldim(1)) print *," CINOUT ind 1 err ", i, h
                    if (j < 1 .or. j > ldim(2)) print *," CINOUT ind 2 err ", j, g
                    if (k < 1 .or. k > ldim(3)) print *," CINOUT ind 3 err ", k, f
#endif
                    cinout(i,j,k) = self%cmat(phys(1),phys(2),phys(3))
                    i=i+1
                enddo
            end do
        end do
        !$omp end parallel do

        ! copy input to device
        inplace_d = cinout
        ! Initialize the plan for complex to complex transform
        call cufftPlan3D(plan,ldim(1),ldim(2),ldim(3),CUFFT_C2C)
        ! Execute  Backward transform in place
        call cufftExec(plan,CUFFT_C2C,inplace_d,inplace_d,CUFFT_INVERSE)
        ! Copy results back to host
        cinout = inplace_d

        ! Save to image obj
        !$omp parallel do collapse(3) default(shared) private(i,j,k,ii,jj,kk)
        do k=1,ldim(3)
            do j=1, ldim(2)
                do i=1, ldim(1)
                    jj = mod(j+(ldim(2)+1)/2,ldim(2))
                    if(jj==0) jj = ldim(2)
                    kk = mod(k+(ldim(3)+1)/2,ldim(3))
                    if(kk==0) kk = ldim(3)
                    ii = mod(i+(ldim(1)+1)/2,ldim(1))
                    if(ii==0) ii = ldim(1)
                    self%rmat(ii,jj,kk) = cabs(cinout(i,j,k))
                end do
            end do
        end do
        !$omp end parallel do

        ! Release memory on the host and on the device
        deallocate( cinout, inplace_d ,stat=alloc_stat)

        ! Destroy the plans
        call cufftDestroy(plan)
        DebugPrint "In simple_image::ifft_pgi done "
#endif
    end subroutine ifft_pgi_cuda

    !> \brief dampens the central cross of a powerspectrum by median filtering
    subroutine dampen_central_cross( self )
        class(image), intent(inout) :: self
        integer            :: h,mh,k,mk,lims(3,2),inds(3)
        integer, parameter :: XDAMPWINSZ=2
        real, allocatable  :: pixels(:)
        if( self%ft )          stop 'not intended for FTs; simple_image :: dampen_central_cross'
        if( self%ldim(3) > 1 ) stop 'not intended for 3D imgs; simple_image :: dampen_central_cross'
        lims = self%loop_lims(3)
        mh = maxval(lims(1,:))
        mk = maxval(lims(2,:))
        inds = 1
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                if( h == 0 .or. k == 0 )then
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    inds(2) = min(max(1,k+mk+1),self%ldim(2))
                    pixels = self%win2arr(inds(1), inds(2), 1, XDAMPWINSZ)
                    call self%set(inds, median_nocopy(pixels))
                endif
            end do
        end do
    end subroutine dampen_central_cross

    !> \brief subtracts the background of an image by subtracting a low-pass filtered
    !!        version of itself
    subroutine subtr_backgr( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
        type(image) :: tmp
        integer     :: winsz
        call tmp%copy(self)
        winsz = nint((self%ldim(1) * self%smpd) / lp)
        call tmp%real_space_filter(winsz, 'average')
        self%rmat = self%rmat - tmp%rmat
        call tmp%kill
    end subroutine subtr_backgr

    subroutine subtr_avg_and_square( self )
        class(image), intent(inout) :: self
        real :: avg
        avg = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        self%rmat = self%rmat - avg
        self%rmat = self%rmat * self%rmat
    end subroutine subtr_avg_and_square

    !> \brief generates a real-space resolution mask for matching power-spectra
    subroutine resmsk( self, hplim, lplim )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hplim, lplim
        integer :: h, k, lims(3,2), mh, mk, inds(3)
        real    :: freq, hplim_freq, lplim_freq
        hplim_freq = self%fit%get_find(1,hplim)
        lplim_freq = self%fit%get_find(1,lplim)
        lims = self%loop_lims(3)
        mh = maxval(lims(1,:))
        mk = maxval(lims(2,:))
        inds = 1
        self%rmat = 1.0
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                inds(1) = min(max(1,h+mh+1),self%ldim(1))
                inds(2) = min(max(1,k+mk+1),self%ldim(2))
                freq = hyp(real(h),real(k))
                if(freq .lt. hplim_freq .or. freq .gt. lplim_freq )then
                    call self%set(inds, 0.)
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
        do k=1,fdim(self1%ldim(1))-3
            call maskimg%ring(self1%ldim, self1%smpd, real(k+2), real(k-2), npix )
            l_mask = bin2logical(maskimg)
            corrs(k) = self1%real_corr(self2, l_mask)
        end do
    end subroutine frc_pspec

    !>  \brief  an image shifter to prepare for Fourier transformation
    subroutine shift_phorig( self )
        class(image), intent(inout) :: self
        integer :: i, j, k
        real    :: rswap
        integer :: kfrom,kto
        if( self%ft ) call simple_stop('ERROR, this method is intended for real images; shift_phorig; simple_image')
        if( self%even_dims() )then
            if( self%ldim(3) == 1 )then
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
                kfrom = 1
                kto   = self%ldim(3)/2
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
            endif
        else
            write(*,*) 'ldim: ', self%ldim
            call simple_stop('even dimensions assumed; shift_phorig; simple_image')
        endif
    end subroutine shift_phorig

    !> \brief shift  is for origin shifting an image
    !! \param x position in axis 0
    !! \param y position in axis 1
    !! \param z position in axis 2
    !! \param lp_dyn low-pass cut-off freq
    !! \param imgout processed image
    !!
    ! subroutine shift( self, shvec, lp_dyn, imgout )
    !     class(image),           intent(inout) :: self
    !     real,                   intent(in)    :: shvec(3)
    !     real,         optional, intent(in)    :: lp_dyn
    !     class(image), optional, intent(inout) :: imgout
    !     integer :: h, k, l, lims(3,2), phys(3)
    !     real    :: shvec_here(3)
    !     logical :: didft, imgout_present, lp_dyn_present
    !     imgout_present = present(imgout)
    !     lp_dyn_present = present(lp_dyn)
    !     if( all(abs(shvec) < TINY) )then
    !         if( imgout_present ) call imgout%copy(self)
    !         return
    !     endif
    !     shvec_here = shvec
    !     if( self%ldim(3) == 1 ) shvec_here(3) = 0.0
    !     didft = .false.
    !     if( .not. self%ft )then
    !         call self%fft()
    !         didft = .true.
    !     endif
    !     if( lp_dyn_present )then
    !         lims = self%fit%loop_lims(1,lp_dyn)
    !     else
    !         lims = self%fit%loop_lims(2)
    !     endif
    !     if( imgout_present )then
    !         imgout%ft = .true.
    !         !$omp parallel do collapse(3) default(shared) private(phys,h,k,l)&
    !         !$omp schedule(static) proc_bind(close)
    !         do h=lims(1,1),lims(1,2)
    !             do k=lims(2,1),lims(2,2)
    !                 do l=lims(3,1),lims(3,2)
    !                     phys = self%fit%comp_addr_phys([h,k,l])
    !                     imgout%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
    !                     self%oshift([h,k,l], shvec_here)
    !                 end do
    !             end do
    !         end do
    !         !$omp end parallel do
    !     else
    !         !$omp parallel do collapse(3) default(shared) private(phys,h,k,l)&
    !         !$omp schedule(static) proc_bind(close)
    !         do h=lims(1,1),lims(1,2)
    !             do k=lims(2,1),lims(2,2)
    !                 do l=lims(3,1),lims(3,2)
    !                     phys = self%fit%comp_addr_phys([h,k,l])
    !                     self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
    !                     self%oshift([h,k,l], shvec_here)
    !                 end do
    !             end do
    !         end do
    !         !$omp end parallel do
    !     endif
    !     if( didft )then
    !         call self%ifft()
    !         if( imgout_present ) call imgout%ifft()
    !     endif
    ! end subroutine shift

    !> \brief shift  is for origin shifting an image
    !! \param x position in axis 0
    !! \param y position in axis 1
    !! \param z position in axis 2
    !!
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

    subroutine shift2Dserial( self, shvec  )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(2)
        integer :: h, k, lims(3,2), phys(3)
        real    :: shvec_here(3)
        lims            = self%fit%loop_lims(2)
        shvec_here(1:2) = shvec
        shvec_here(3)   = 0.
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                phys = self%fit%comp_addr_phys([h,k,0])
                self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) *&
                    &self%oshift([h,k,0], shvec_here)
            end do
        end do
    end subroutine shift2Dserial

    !> \brief mask  is for spherical masking
    !! \param mskrad mask radius
    !! \param which image type
    !! \param inner include cosine edge material
    !! \param width width of inner patch
    !! \param msksum masking sum
    !!
    subroutine mask( self, mskrad, which, inner, width )
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        character(len=*), intent(in)    :: which
        real, optional,   intent(in)    :: inner, width
        real    :: e, wwidth
        real    :: cis(self%ldim(1)), cjs(self%ldim(2)), cks(self%ldim(3))
        integer :: i, j, k, minlen, ir, jr, kr
        logical :: didft, doinner, soft
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
        select case(trim(which))
        case('soft')
            soft  = .true.
            call self%zero_background
        case('hard')
            soft  = .false.
        case DEFAULT
            call simple_stop('undefined which parameter; mask; simple_image')
        end select
        ! init center as origin
        forall(i=1:self%ldim(1)) cis(i) = -real(self%ldim(1))/2. + real(i-1)
        forall(i=1:self%ldim(2)) cjs(i) = -real(self%ldim(2))/2. + real(i-1)
        if(self%is_3d())forall(i=1:self%ldim(3)) cks(i) = -real(self%ldim(3))/2. + real(i-1)
        ! MASKING
        if( soft )then
            ! Soft masking
            if( self%is_3d() )then
                ! 3d
                do i=1,self%ldim(1)/2
                    ir = self%ldim(1)+1-i
                    do j=1,self%ldim(2)/2
                        jr = self%ldim(2)+1-j
                        do k=1,self%ldim(3)/2
                            kr = self%ldim(3)+1-k
                            e = cosedge(cis(i),cjs(j),cks(k),minlen,mskrad)
                            if( doinner )e = e * cosedge_inner(cis(i),cjs(j),cks(k),wwidth,inner)
                            if(e > 0.9999) cycle
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
                        e = cosedge(cis(i),cjs(j),minlen,mskrad)
                        if( doinner )e = e * cosedge_inner(cis(i),cjs(j),wwidth,inner)
                        if(e > 0.9999)cycle
                        self%rmat(i,j,1)   = e * self%rmat(i,j,1)
                        self%rmat(i,jr,1)  = e * self%rmat(i,jr,1)
                        self%rmat(ir,j,1)  = e * self%rmat(ir,j,1)
                        self%rmat(ir,jr,1) = e * self%rmat(ir,jr,1)
                    enddo
                enddo
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
                            kr = self%ldim(3)+1-k
                            e = hardedge(cis(i),cjs(j),cks(k),mskrad)
                            if( doinner )e = e * hardedge_inner(cis(i),cjs(j),cks(k),inner)
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
                        if( doinner )e = e * hardedge_inner(cis(i),cjs(j),inner)
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
    !!
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
    !!
    subroutine pad( self_in, self_out, backgr )
        class(image),  intent(inout) :: self_in, self_out
        real, optional, intent(in)   :: backgr
        real              :: w, ratio
        integer           :: starts(3), stops(3), lims(3,2)
        integer           :: h, k, l, phys_in(3), phys_out(3)
        real, allocatable :: antialw(:)
        if( self_in.eqdims.self_out )then
            call self_out%copy(self_in)
            return
        endif
        if( self_out%ldim(1) >= self_in%ldim(1) .and. self_out%ldim(2) >= self_in%ldim(2)&
        .and. self_out%ldim(3) >= self_in%ldim(3) )then
            if( self_in%ft )then
                self_out = cmplx(0.,0.)
                !$ allocate(antialw(1), stat=alloc_stat)
                !$ if(alloc_stat.ne.0)call allocchk("In simple_image::pad PGI prealloc antialw",alloc_stat)
                antialw = self_in%hannw()
                lims = self_in%fit%loop_lims(2)
                !$omp parallel do collapse(3) schedule(static) default(shared)&
                !$omp private(h,k,l,w,phys_out,phys_in) proc_bind(close)
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        do l=lims(3,1),lims(3,2)
                            w = antialw(max(1,abs(h)))*antialw(max(1,abs(k)))*antialw(max(1,abs(l)))
                            phys_out = self_out%fit%comp_addr_phys([h,k,l])
                            phys_in  = self_in%fit%comp_addr_phys([h,k,l])
                            self_out%cmat(phys_out(1),phys_out(2),phys_out(3))=&
                            self_in%cmat(phys_in(1),phys_in(2),phys_in(3))*w
                        end do
                    end do
                end do
                !$omp end parallel do
                deallocate(antialw)
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

    !> \brief pad_mirr is a constructor that pads the input image to input ldim in real space using mirroring
    !! \param self_in image object
    !! \param self_out image object
    !!
    subroutine pad_mirr( self_in, self_out )
        !use simple_winfuns, only: winfuns
        class(image),   intent(inout) :: self_in, self_out
        integer :: starts(3), stops(3)
        integer :: i,j, i_in, j_in
        if( self_in.eqdims.self_out )then
            call self_out%copy(self_in)
            return
        endif
        if(self_in%is_3d())call simple_stop('2D images only; simple_image::pad_mirr')
        if(self_in%ft)call simple_stop('real space 2D images only; simple_image::pad_mirr')
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
            call simple_stop('Inconsistent dimensions; simple_image::pad_mirr')
        endif
    end subroutine pad_mirr

    !> \brief clip is a constructor that clips the input image to input ldim
    !!             DO NOT PARALLELISE
    !! \param self_in image object
    !! \param self_out image object
    !!
    subroutine clip( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
        real                        :: ratio
        integer                     :: starts(3), stops(3), lims(3,2)
        integer                     :: phys_out(3), phys_in(3), h, k, l
        if( self_out%ldim(1) <= self_in%ldim(1) .and. self_out%ldim(2) <= self_in%ldim(2)&
        .and. self_out%ldim(3) <= self_in%ldim(3) )then
            if( self_in%ft )then
                lims = self_out%fit%loop_lims(2)
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        do l=lims(3,1),lims(3,2)
                            phys_out = self_out%fit%comp_addr_phys([h,k,l])
                            phys_in = self_in%fit%comp_addr_phys([h,k,l])
                            self_out%cmat(phys_out(1),phys_out(2),phys_out(3)) =&
                            self_in%cmat(phys_in(1),phys_in(2),phys_in(3))
                        end do
                    end do
                end do
                ratio = real(self_in%ldim(1))/real(self_out%ldim(1))
                self_out%smpd = self_in%smpd*ratio ! clipping Fourier transform, so sampling is coarser
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

    !> \brief clip_inplace is a constructor that clips the input image to input ldim
    !! \param ldim
    !!
    subroutine clip_inplace( self, ldim )
        class(image), intent(inout) :: self
        integer, intent(in)         :: ldim(3)
        type(image)                 :: tmp
        call tmp%new(ldim, self%smpd)
        call self%clip(tmp)
        call self%copy(tmp)
        call tmp%kill()
    end subroutine clip_inplace

    !>  \brief  is for mirroring an image
    !!          mirror('x') corresponds to mirror2d
    subroutine mirror( self, md )
        class(image), intent(inout) :: self
        character(len=*), intent(in) :: md
        integer :: i, j
        logical :: didft
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        if( md == 'x' )then
            do i=1,self%ldim(2)
                do j=1,self%ldim(3)
                    call reverse(self%rmat(1:self%ldim(1),i,j))
                end do
            end do
        else if( md == 'y' )then
            do i=1,self%ldim(1)
                do j=1,self%ldim(3)
                    call reverse(self%rmat(i,1:self%ldim(2),j))
                end do
            end do
        else if( md == 'z' )then
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    call reverse(self%rmat(i,j,1:self%ldim(3)))
                end do
            end do
        else
            write(*,'(a)') 'Mode needs to be either x, y or z; mirror; simple_image'
        endif
        if( didft ) call self%fft()
    end subroutine mirror

    subroutine norm( self  )
        class(image), intent(inout) :: self
        integer :: npix
        real    :: ave, var, ep
        npix   = product(self%ldim)
        ave    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        if( abs(ave) > TINY ) self%rmat = self%rmat - ave
        ep     = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        var    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0)
        var    = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        if( is_a_number(var) )then
            if( var > 0. ) self%rmat = self%rmat / sqrt(var)
        endif
    end subroutine norm

    !> \brief norm_ext  is for normalization of an image using inputted average and standard deviation
    !! \param avg Average
    !! \param sdev Standard deviation
    !!
    subroutine norm_ext( self, avg, sdev )
        class(image), intent(inout) :: self
        real, intent(in)            :: avg, sdev
        if( self%ft )then
            write(*,*) 'WARNING: Cannot normalize FTs; norm_ext; simple_image'
            return
        endif
        if( abs(avg) > TINY ) self%rmat = self%rmat - avg
        if( sdev     > 0.   ) self%rmat = self%rmat / sdev
    end subroutine norm_ext

    !> \brief noise_norm  normalizes the image according to the background noise
    !! \param msk Mask value
    !! \param errout error flag
    !!
    subroutine noise_norm( self, msk, errout )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: msk
        logical, optional, intent(out)   :: errout
        type(image)       :: maskimg
        integer           :: npix,nbackgr,npix_tot
        real              :: med,ave,sdev,var
        real, allocatable :: pixels(:)
        logical           :: err
        if( self%ft )then
            write(*,*) 'WARNING: Cannot normalize FTs; noise_norm; simple_image'
            return
        endif
        call maskimg%disc(self%ldim, self%smpd, msk, npix)
        npix_tot = product(self%ldim)
        nbackgr = npix_tot-npix
        ! possible realloc
        pixels = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
            &maskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0.5 )
        med = median_nocopy(pixels)
        call moment(pixels, ave, sdev, var, err)
        deallocate(pixels)
        if( err )then
            call self%norm()
        else
            ! we subtract the pixel values with this ruboust estimate of the background average
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            (self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))-med)
            if( present(errout) )then
                errout = err
            else
                if( err ) write(*,'(a)') 'WARNING: variance zero; noise_norm; simple_image'
            endif
            if( sdev > 1e-6 )then
                self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
                    self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))/sdev
            endif
        endif
    end subroutine noise_norm

    !>  \brief  putting the edge around the image to zero (necessary for avoiding FT artefacts)
    subroutine edges_norm( self )
        class(image), intent(inout) :: self
        real :: edges_sum, edges_ave
        if( self%ft ) stop 'not for Fted images; simple_image :: edge_norm'
        if( .not.self%is_2d() ) stop 'only for 2d images; simple_image :: edge_norm'
        edges_sum = sum(self%rmat(1:self%ldim(1),1,1))
        edges_sum = edges_sum + sum(self%rmat(1:self%ldim(1),self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(1,1:self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(self%ldim(1),1:self%ldim(2),1))
        edges_ave = edges_sum / real( 2*(self%ldim(1)+self%ldim(2)) )
        if( abs(edges_ave) > TINY )self%rmat = self%rmat - edges_ave
    end subroutine edges_norm

    !> \brief radius_norm  normalizes the image based on a central sphere of input radius
    !! \param radius Radius of sphere
    !! \param errout error flag
    !!
    subroutine radius_norm( self, radius, errout )
        class(image),      intent(inout) :: self
        real,    optional, intent(in)    :: radius
        logical, optional, intent(out)   :: errout
        type(image)       :: maskimg
        integer           :: npix
        real              :: ave,sdev,var,irad
        real, allocatable :: pixels(:)
        logical           :: err
        if( self%ft )then
            write(*,*) 'WARNING: Cannot normalize FTs; noise_norm; simple_image'
            return
        endif
        if( .not.present(radius) )then
            irad = real(self%ldim(1)+1)/2.
        else
            irad = radius
        endif
        if( irad<=0. )call simple_stop('Invalid radius value in rad_norm')
        call maskimg%disc(self%ldim, self%smpd, irad, npix)
        ! pixels = self%packer(maskimg) ! Intel hickup
        pixels = pack( self%rmat, maskimg%rmat > 0.5 )
        call moment( pixels, ave, sdev, var, err )
        deallocate(pixels)
        if( err )then
            call self%norm()
        else
            self%rmat = self%rmat - ave
            if( present(errout) )then
                errout = err
            else
                if( err ) write(*,'(a)') 'WARNING: variance zero; rad_norm; simple_image'
            endif
            if( sdev > 1e-6 )self%rmat = self%rmat / sdev
        endif
    end subroutine radius_norm

    !>  \brief  is for [0,1] interval normalization of an image
    subroutine norm_bin( self )
        class(image), intent(inout) :: self
        real                        :: smin, smax
        logical                     :: didft
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ! find minmax
        smin  = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        smax  = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        ! create [0,1]-normalized image
        self%rmat = (self%rmat - smin)  / (smax-smin)
        self%rmat = (exp(self%rmat)-1.) / (exp(1.)-1.)
        if( didft ) call self%fft()
    end subroutine norm_bin

    !> \brief roavg  is for creating a rotation average of self
    !! \param angstep angular step
    !! \param avg output image rotation average
    !!
    subroutine roavg( self, angstep, avg )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: angstep
        class(image), intent(inout) :: avg
        real    :: avgs_rmat(nthr_glob,self%ldim(1),self%ldim(2),1)
        real    :: rotated(nthr_glob,self%ldim(1),self%ldim(2),1)
        integer :: irot, ithr
        avgs_rmat = 0.
        rotated   = 0.
        !$omp parallel do schedule(static) default(shared) private(irot) proc_bind(close)
        do irot = 0 + angstep,359,angstep
            ! get thread index
            ithr = omp_get_thread_num() + 1
            ! rotate & sum
            call self%rtsq_serial(real(irot), 0., 0., rotated(ithr,:,:,:))
            avgs_rmat(ithr,:,:,:) = avgs_rmat(ithr,:,:,:) + rotated(ithr,:,:,:)
        end do
        !$omp end parallel do
        ! add in the zero rotation
        avgs_rmat(1,:,:,:) = avgs_rmat(1,:,:,:) + self%rmat(:self%ldim(1),:self%ldim(2),:)
        ! set output image object
        call avg%new(self%ldim, self%smpd)
        call avg%set_rmat(sum(avgs_rmat, dim=1))
        ! normalise
        call avg%div(real(360/angstep))
    end subroutine roavg

    !> \brief rtsq  rotation of image by quadratic interpolation (from spider)
    !! \param self_in image object
    !! \param ang angle of rotation
    !! \param shxi shift in x axis
    !! \param shyi shift in y axis
    !! \param self_out optional copy of processed result
    !!
    subroutine rtsq( self_in, ang, shxi, shyi, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: ang,shxi,shyi
        class(image), optional, intent(inout) :: self_out
        type(image) :: self_here
        real    :: shx,shy,ry1,rx1,ry2,rx2,cod,sid,xi,fixcenmshx,fiycenmshy
        real    :: rye2,rye1,rxe2,rxe1,yi,ycod,ysid,yold,xold
        integer :: iycen,ixcen,ix,iy
        real    :: mat_in(self_in%ldim(1),self_in%ldim(2))
        real    :: mat_out(self_in%ldim(1),self_in%ldim(2))
        logical :: didft
        if( self_in%ldim(3) > 1 )         stop 'only for 2D images; rtsq; simple_image'
        if( .not. self_in%square_dims() ) stop 'only for square dims (need to sort shifts out); rtsq; simple_image'
        call self_here%new(self_in%ldim, self_in%smpd)
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
        iycen = self_in%ldim(1)/2+1
        ixcen = self_in%ldim(2)/2+1
        ! image dimensions around origin
        rx1 = -self_in%ldim(1)/2
        rx2 =  self_in%ldim(1)/2
        ry1 = -self_in%ldim(2)/2
        ry2 =  self_in%ldim(2)/2
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
        self_here%rmat(:self_here%ldim(1),:self_here%ldim(2),1) = mat_out
        self_here%ft = .false.
        if( present(self_out) )then
            call self_out%copy(self_here)
        else
            call self_in%copy(self_here)
        endif
        call self_here%kill()
        if( didft )then
            call self_in%ifft()
        endif
    end subroutine rtsq

    !> \brief rtsq  rotation of image by quadratic interpolation (from spider)
    !! \param self_in image object
    !! \param ang angle of rotation
    !! \param shxi shift in x axis
    !! \param shyi shift in y axis
    !! \param rmat_out is the rotated pixels
    !!
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
        iycen = self_in%ldim(1)/2+1
        ixcen = self_in%ldim(2)/2+1
        ! image dimensions around origin
        rx1 = -self_in%ldim(1)/2
        rx2 =  self_in%ldim(1)/2
        ry1 = -self_in%ldim(2)/2
        ry2 =  self_in%ldim(2)/2
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
    !! \param outliers -
    !!
    subroutine cure_outliers( self, ncured, nsigma, deadhot, outliers )
        class(image),      intent(inout) :: self
        integer,           intent(inout) :: ncured
        real,              intent(in)    :: nsigma
        integer,           intent(out)   :: deadhot(2)
        logical, optional, allocatable   :: outliers(:,:)
        real, allocatable :: win(:,:), rmat_pad(:,:)
        real    :: ave, sdev, var, lthresh, uthresh
        integer :: i, j, hwinsz, winsz
        logical :: was_fted, err, present_outliers
        if( self%ldim(3)>1 )stop 'for images only; simple_image::cure_outliers'
        was_fted = self%is_ft()
        if( was_fted )stop 'for real space images only; simple_image::cure_outliers'
        present_outliers = present(outliers)
        ncured   = 0
        hwinsz   = 6
        if( allocated(outliers) ) deallocate(outliers)
        allocate( outliers(self%ldim(1),self%ldim(2)), stat=alloc_stat)
        if(alloc_stat/=0)call allocchk("In simple_image::cure_outliers ")
        outliers = .false.
        call moment( self%rmat, ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh = ave - nsigma * sdev
        uthresh = ave + nsigma * sdev
        if( any(self%rmat<=lthresh) .or. any(self%rmat>=uthresh) )then
            winsz = 2*hwinsz+1
            deadhot = 0
            allocate(rmat_pad(1-hwinsz:self%ldim(1)+hwinsz,1-hwinsz:self%ldim(2)+hwinsz),&
                &win(winsz,winsz), stat=alloc_stat)
            if(alloc_stat/=0)call allocchk('In: cure_outliers; simple_image 1')
            rmat_pad(:,:) = median( reshape(self%rmat(:,:,1), (/(self%ldim(1)*self%ldim(2))/)) )
            rmat_pad(1:self%ldim(1), 1:self%ldim(2)) = &
                &self%rmat(1:self%ldim(1),1:self%ldim(2),1)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,win)&
            !$omp reduction(+:ncured) proc_bind(close)
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

    !>  \brief  denoise_NLM for denoising image with non-local means algorithm
    !! \param Hsigma power of noise cancelling
    !! \param searchRad search radius from pixel origin to patch origin
    !! \param patchSz patch width
    !! \param deadhot number of rejections for upper and lower limits
    !! \param outliers mask of rejection points
    !! \param l1normdiff return the L1-norm difference
    subroutine denoise_NLM( self, Hsigma, patchSz,searchRad, deadhot, outliers, l1normdiff)
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: patchSz
        real,              intent(in)    :: Hsigma,searchRad
        integer,           intent(out)   :: deadhot(2)
        logical, optional, allocatable   :: outliers(:,:)
        real, optional,    intent(inout) :: l1normdiff
        type(image) :: selfcopy
        real, allocatable :: patch(:,:), padded_image(:,:)
        real    :: ave, sdev, var, lthresh, uthresh, nsigma
        integer :: i, j, hwinsz, winsz,ncured
        logical :: was_fted, err, present_outliers, retl1norm
        if( self%ldim(3)>1 )stop 'for images only; simple_image:: denoise_NLM'
        was_fted = self%is_ft()
        if( was_fted )stop 'for real space images only; simple_image::denoise_NLM'
        present_outliers = present(outliers)
        retl1norm=.false.
        if(present(l1normdiff))then
            retl1norm=.true.
            selfcopy = self    ! create copy for comparision
        end if
        hwinsz   = 6
        nsigma = 1.0/ sqrt(1.0)
        if( allocated(outliers) ) deallocate(outliers)
        allocate( outliers(self%ldim(1),self%ldim(2)) )
        outliers = .false.
        call moment( self%rmat, ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh = ave - nsigma * sdev
        uthresh = ave + nsigma * sdev
        if( any(self%rmat<=lthresh) .or. any(self%rmat>=uthresh) )then
            winsz = 2*hwinsz+1
            deadhot = 0
            allocate(padded_image(1-hwinsz:self%ldim(1)+hwinsz,1-hwinsz:self%ldim(2)+hwinsz),&
                &patch(winsz,winsz), stat=alloc_stat)
            if(alloc_stat/=0)call allocchk('In: cure_outliers; simple_image 1')
            padded_image(:,:) = median( reshape(self%rmat(:,:,1), (/(self%ldim(1)*self%ldim(2))/)) )
            padded_image(1:self%ldim(1), 1:self%ldim(2)) = &
                &self%rmat(1:self%ldim(1),1:self%ldim(2),1)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,patch)&
            !$omp reduction(+:ncured) proc_bind(close)
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! if( self%rmat(i,j,1)<lthresh .or. self%rmat(i,j,1)>uthresh )then  !
                         if( present_outliers )then                                    !
                             outliers(i,j)=.true.                                      !
                    !         if (self%rmat(i,j,1)<lthresh) deadhot(1) = deadhot(1) + 1 !
                    !         if (self%rmat(i,j,1)>uthresh) deadhot(2) = deadhot(2) + 1 !
                         else                                                          !
                    !         patch = padded_image( i-hwinsz:i+hwinsz, j-hwinsz:j+hwinsz )    !
                    !         self%rmat(i,j,1) = median( reshape(patch,(/winsz**2/)) )    !
                            ncured = ncured + 1                                       !
                    !     endif                                                         !
                    endif                                                             !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                enddo
            enddo
            !$omp end parallel do
            deallocate( patch, padded_image )
        endif

        if(retl1norm)then
            l1normdiff = l1norm_2(selfcopy,self)
            call selfcopy%kill()
        end if
    end subroutine denoise_NLM

    !>  \brief  zero pixels below thres
    subroutine zero_below( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        where( self%rmat < thres ) self%rmat = 0.
    end subroutine zero_below

    !>  \brief  is the image class unit test
    subroutine test_image( doplot )
        logical, intent(in)  :: doplot
        write(*,'(a)') '**info(simple_image_unit_test): testing square dimensions'
        call test_image_local( 100, 100, 100, doplot )
        !        write(*,'(a)') '**info(simple_image_unit_test): testing non-square dimensions'
        !        call test_image_local( 120, 90, 80, doplot )
        write(*,'(a)') 'SIMPLE_IMAGE_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'

    contains

        subroutine test_image_local( ld1, ld2, ld3, doplot )
            integer, intent(in)  :: ld1, ld2, ld3
            logical, intent(in)  :: doplot
            type(image)          :: img, img_2, img_3, img_4, img3d
            type(image)          :: imgs(20)
            integer              :: i, j, k, cnt, ldim(3)
            real                 :: input, msk, ave, sdev, med, xyz(3)
            real                 :: corr, corr_lp, maxv, minv
            real, allocatable    :: pcavec1(:), pcavec2(:)
            real                 :: smpd=2.
            logical              :: passed, test(6)

            write(*,'(a)') '**info(simple_image_unit_test, part 1): testing basal constructors'
            call img%new([ld1,ld2,1], 1.)
            call img_3%new([ld1,ld2,1], 1.)
            call img3d%new([ld1,ld2,ld3], 1.)
            if( .not. img%exists() ) stop 'ERROR, in constructor or in exists function, 1'
            if( .not. img3d%exists() ) stop 'ERROR, in constructor or in exists function, 2'

            write(*,'(a)') '**info(simple_image_unit_test, part 2): testing getters/setters'
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
            if( .not. passed )  stop 'getters/setters test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 3): testing serialization'
            passed = .false.
            msk = 50.
            img_2 = img
            call img%ran
            if( doplot ) call img%vis
            call img%serialize(pcavec1, msk)
            img = 0.
            call img%serialize(pcavec1, msk)
            if( doplot ) call img%vis
            call img_2%serialize(pcavec1, msk)
            call img_2%serialize(pcavec2, msk)
            if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'serialization test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 4): testing checkups'
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
            if( .not. passed ) stop 'checkups test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 5): testing arithmetics'
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            passed = .false.
            msk = 50.
            call img%ran
            call img%serialize(pcavec1, msk)
            img_2 = img
            call img_2%serialize(pcavec2, msk)
            if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'polymorphic assignment test 1 failed'
            passed = .false.
            img = 5.
            img_2 = 10.
            img_3 = img_2-img
            call img%serialize(pcavec1, msk)
            call img_3%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'overloaded subtraction test failed'
            passed = .false.
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            img = 5.
            img_2 = 10.
            img_3 = 15.
            img_4 = img + img_2
            call img_3%serialize(pcavec1, msk)
            call img_4%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'overloaded addition test failed'
            passed = .false.
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            img = 5.
            img_2 = 2.
            img_3 = 10.
            img_4 = img*img_2
            call img_3%serialize(pcavec1, msk)
            call img_4%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'overloaded multiplication test failed'
            passed = .false.
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            img_4 = img_3/img_2
            call img%serialize(pcavec1, msk)
            call img_4%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'overloaded division test failed'
            passed = .false.
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            img = 0.
            img_2 = 5.
            img_3 = 1.
            do i=1,5
                call img%add(img_3 )
            end do
            call img_2%serialize(pcavec1, msk)
            call img%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'summation test failed'
            passed = .false.
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            do i=1,5
                call img%subtr(img_3)
            end do
            img_2 = 0.
            call img_2%serialize(pcavec1, msk)
            call img%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'subtraction test failed'
            passed = .false.
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            img_2 = 5.
            img_3 = 1.
            call img_2%div(5.)
            call img_2%serialize(pcavec1, msk)
            call img_3%serialize(pcavec2, msk)
            if( euclid(pcavec1, pcavec2) < 0.0001 ) passed = .true.
            if( .not. passed ) stop 'constant division test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 6): testing stats'
            passed = .false.
            call img%gauran( 5., 15. )
            call img%stats( 'foreground', ave, sdev, maxv, minv, 40., med )
            if( ave >= 4. .and. ave <= 6. .and. sdev >= 14. .and.&
                sdev <= 16. .and. med >= 4. .and. med <= 6. ) passed = .true.
            if( .not. passed )  stop 'stats test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 7): testing origin shift'
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            passed = .false.
            msk=50
            call img%gauimg(10)
            if( doplot ) call img%vis
            call img%serialize(pcavec1, msk)
            call img%shift([-9.345,-5.786,0.])
            if( doplot ) call img%vis
            call img%shift([9.345,5.786,0.])
            call img%serialize(pcavec2, msk)
            if( doplot ) call img%vis
            if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
            if( .not. passed )  stop 'origin shift test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 8): testing masscen'
            passed = .false.
            msk = 50
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            call img%square( 10 )
            if( doplot ) call img%vis
            call img%serialize(pcavec1, msk)
            call img%shift([10.,5.,0.])
            if( doplot ) call img%vis
            xyz = img%masscen()
            call img%shift([real(int(xyz(1))),real(int(xyz(2))),0.])
            if( doplot ) call img%vis
            call img%serialize(pcavec2, msk)
            if( pearsn(pcavec1, pcavec2) > 0.9 ) passed = .true.
            print *,'determined shift:', xyz
            if( .not. passed ) stop 'masscen test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 9): testing lowpass filter'
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

            write(*,'(a)') '**info(simple_image_unit_test, part 10): testing spherical mask'
            call img%ran
            if( doplot ) call img%vis
            call img%mask(35.,'hard')
            if( doplot ) call img%vis
            call img%ran
            call img%mask(35.,'soft')
            if( doplot ) call img%vis

            write(*,'(a)') '**info(simple_image_unit_test, part 11): testing padding/clipping'
            passed = .false.
            msk = 50
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            call img%ran
            call img%serialize(pcavec1, msk)
            if( doplot ) call img%vis
            call img_2%new([2*ld1,2*ld2,1],1.)
            call img%pad(img_2)
            if( doplot ) call img_2%vis
            call img_3%new([ld1,ld2,1],1.)
            call img%clip(img_3)
            call img_3%serialize(pcavec2, msk)
            if( doplot ) call img_3%vis
            if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
            if( .not. passed ) stop 'padding/clipping test failed'
            call img%square(10)
            if( doplot ) call img%vis
            call img%fft()
            call img%pad(img_2)
            call img_2%ifft()
            if( doplot ) call img_2%vis
            call img_2%square(20)
            if( doplot ) call img_2%vis
            call img_2%fft()
            call img_2%clip(img)
            call img%ifft()
            if( doplot ) call img%vis

            write(*,'(a)') '**info(simple_image_unit_test, part 13): testing bicubic rots'
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

            write(*,'(a)') '**info(simple_image_unit_test, part 14): testing binary imgproc routines'
            passed = .false.
            call img%gauimg(20)
            call img%norm_bin
            if( doplot ) call img%vis
            call img%bin(0.5)
            if( doplot ) call img%vis
            call img%gauimg(20)
            call img%bin(500)
            if( doplot ) call img%vis
            do i=1,10
                call img%grow_bin()
            end do
            if( doplot ) call img%vis

            write(*,'(a)') '**info(simple_image_unit_test, part 15): testing auto correlation function'
            call img%square( 10 )
            if( doplot ) call img%vis
            call img%acf
            if( doplot ) call img%vis
            call img%square( 10 )
            call img%shift([5.,-5.,0.])
            if( doplot ) call img%vis
            call img%acf
            if( doplot ) call img%vis

            write(*,'(a)') '**info(simple_image_unit_test, part 16): testing correlation functions'
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
            if( .not. passed ) stop 'corr test failed'

            write(*,'(a)') '**info(simple_image_unit_test, part 17): testing downscaling'
            if( ld1 == ld2 )then
                call img%gauimg(20)
                if( doplot )  call img%vis
                if( doplot ) call img_2%vis
            endif

            if( img%square_dims() .and. nthr_glob > 2 )then
                write(*,'(a)') '**info(simple_image_unit_test, part 19): testing rotational averager'
                call img%square( 10 )
                if( doplot ) call img%vis
                call img%roavg(5,img_2)
                if( doplot ) call img_2%vis
            endif

            write(*,'(a)') '**info(simple_image_unit_test, part 20): testing the read/write capabilities'
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
                        call simple_stop('SPIDER vs. MRC & converted vs. nonconverted test failed',&
                            __FILENAME__,__LINE__)
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
                        call simple_stop('SPIDER vs. MRC & converted vs. nonconverted test failed', &
                            & __FILENAME__,__LINE__)
                    endif
                end do
            end do

            write(*,'(a)') '**info(simple_image_unit_test, part 21): testing destructor'
            passed = .false.
            call img%kill()
            call img3d%kill()
            test(1) = .not. img%exists()
            test(2) = .not. img3d%exists()
            passed = all(test)
            if( .not. passed )  stop 'destructor test failed'
        end subroutine test_image_local

    end subroutine test_image

    subroutine test_image_pgi_cuda( doplot )
        logical, intent(in)  :: doplot
        write(*,'(a)') '**info(simple_image PGI unit_test): testing square dimensions'
        call test_image_pgi( 100, 100, 100, doplot )
        write(*,'(a)') '**info(simple_image PGI unit_test): testing non-square dimensions'
        ! call test_image_pgi( 120, 90, 80, doplot )
        write(*,'(a)') 'SIMPLE_IMAGE_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'



    contains



        subroutine test_image_pgi( ld1, ld2, ld3, doplot )
#ifdef PGI
            use cudafor
#endif
            integer, intent(in)  :: ld1, ld2, ld3
            logical, intent(in)  :: doplot
            type(image)          :: img, img_2, img_3, img_4, img3d
            integer              :: i, j, k, h, ldim(3), phys(3),lims(3,2)
            real                 :: msk, ave
            real, allocatable    :: pcavec1(:), pcavec2(:)
            real                 :: smpd=2.
            logical              :: passed


            verbose=.true.
            debug=.true.
            write(*,'(a)') '**info(simple_image PGI unit_test, part 1): testing basal constructors'

            call img%new([ld1,ld2,1], 1.)
            write(*,'(a)') '**info(simple_image PGI unit_test, part 1):'
            call img_3%new([ld1,ld2,1], 1.)
            write(*,'(a)') '**info(simple_image PGI unit_test, part 1):'
            call img3d%new([ld1,ld2,1], 1.)
            call img_4%new([ld1,ld2,1], 1.)
            write(*,'(a)') '**info(simple_image PGI unit_test, part 1):'
            if( .not. img%exists() ) call simple_stop('ERROR, in constructor or in exists function, 1')
            if( .not. img3d%exists() ) call simple_stop('ERROR, in constructor or in exists function, 2')

            write(*,'(a)') '**info(simple_image_PGI unit_test, PGI CUFFT 2): testing CUFFT and FFTW3'
            passed = .false.
            ldim = [ld1,ld2,1]
            call img%new(ldim, smpd)
            call img_2%new(ldim, smpd)
            call img_3%new(ldim, smpd)
            call img_4%new(ldim, smpd)
            img%cmat=cmplx(0.,0.)
            img%rmat=1.0
            call img%gauimg(10)
            call gnufor_image(log10(img%rmat(:ldim(1),:ldim(2),1)), palette='gray')
            call gnufor_image(log10(img%rmat(:,:,1)), palette='gray')
            i=1;j=1
            lims = img%fit%loop_lims(2)
            print *, "Rmat ", img%ldim, img%array_shape
            print *, "Rmat size", size(img%rmat,1),  size(img%rmat,2),  size(img%rmat,3)
            print *, "Loop lims ", lims
            do h = lims(1,1),lims(1,2)
                j=1
                do k = lims(2,1),lims(2,2)
                    phys = img%comp_addr_phys([h,k,0])
                    img%cmat(phys(1),phys(2),phys(3)) = cmplx(img%rmat(i,j,1), real(i*j)/100.)
                    j=j+1
                end do
                i=i+1
            end do
            call gnufor_image(log10(real(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1))), palette='gray')                               ! Show img
            call gnufor_image(log10(aimag(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1))), palette='gray')
            call gnufor_image(real(img%cmat(:,:,1)), palette='gray')                               ! Show img
            call gnufor_image(aimag(img%cmat(:,:,1)), palette='gray')
            write(*,'(a)') '**info(simple_image PGI test): testing fftw3 '
            read *
            call img%fft()

            write(*,'(a)') '**info(simple_image PGI unit_test, part 2): PGI FFT completed '
            call gnufor_image(log10(real(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1))), palette='gray')                               ! Show img
            call gnufor_image(log10(aimag(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1))), palette='gray')
            call gnufor_image(log10(cabs(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1))), palette='gray')                               ! Show img
            call gnufor_image(log10(atan2(real(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1)),&
                &aimag(img%cmat(lims(1,1):lims(1,2),lims(2,1):lims(2,2),1)))), palette='gray')

            call gnufor_image(real(img%cmat(:,:,1)), palette='gray')                               ! Show img
            call gnufor_image(aimag(img%cmat(:,:,1)), palette='gray')                               ! Show img

            call img%gauimg(10)
            call img%add_gauran(1.)
            write(*,'(a)') '**info(simple_image PGI unit_test, part 2): Noisy image'
            !                call img_2%gauimg(10)
            !                call img_2%add_gauimg(15)
            !                if(doplot)call img_2%vis                               ! Show img_2
            img_2 = img
            ave = (img.lonesum.img_2)
            print *," L1 norm sum of img copy [must be zero]", ave

            !  call img_2%ifftshift3r()
            !                if(doplot)call img_2%vis()                             ! Show img_2
            !                if(doplot)call img%vis                                 ! Show img
            write(*,'(a)') '**info(simple_image PGI unit_test, part 2): ifftshift3r prior to fft [PGI]'
            print *,'Press enter to run PGI R=>C FFT'
            read(*,*)
            call img_2%fft_pgi_cuda()
            write(*,'(a)') '**info(simple_image PGI unit_test, part 2): PGI FFT completed '
            !img_4 = img_2
            call gnufor_image(real(img_2%cmat(:,:,1)), palette='gray')                               ! Show img_2
            call gnufor_image(aimag(img_2%cmat(:,:,1)), palette='gray')                               ! Show img_2
            print *,'Press enter to run FFTW3 FFT'
            read(*,*)


            print *,'Press enter to show PGI amplitude of cmat'
            read *
            call img_2%ft2img('amp',img_4)
            if(doplot)call img_4%vis(geomorsphr=.false.)


            !  call img_2%fftshift3c()
            print *,'Press enter to show 2nd FFTSHIFT complex images'
            read *

            call gnufor_image(real(img_2%cmat(:,:,1)), palette='gray')                               ! Show img_2
            call gnufor_image(aimag(img_2%cmat(:,:,1)), palette='gray')                               ! Show img_2

            ! call img_2%ft2img('amp',img_4)
            ! if(doplot)call img_4%vis(geomorsphr=.false.)


            !    call img%ft2img('amp',img_3)
            !   if(doplot)call img_3%vis(geomorsphr=.false.)
            print *," L1 norm FFTW and CUFFT ", img.lonesum.img_2


            print *," PGI ", real(img_2%cmat(:,1,1))
            print *," FFTW ", real(img%cmat(:,1,1))

            ! write(*,'(a)') '**info(simple_image PGI test): testing subtracting CUFFT and FFTW3 '
            ! read *
            ! call img%subtr(img_2)
            ! call img%ft2img('amp',img_3)
            ! if(doplot)call img_3%vis(geomorsphr=.false.)
            ! write(*,'(a)') '**info(simple_image PGI test): testing difference '

            ! call img%add(img_2)
            ! read(*,*)
            ! VerbosePrint "In Image::fft_pgi testing "
            ! do i=1,ldim(1)
            !     do j=1,ldim(2)
            !         do k=1,ldim(3)
            !             if (abs(cabs(img%cmat(ldim(1)-i+1,j,k)) - cabs(img_2%cmat(i,j,k))) >= TINY )then
            !                 print *,'wrong FFT i/j/k fftw cufft', i,j,k,real(img%cmat(i,j,k)) ,&
            !                     &aimag(img%cmat(i,j,k)) , real( img_2%cmat(i,j,k)), aimag( img_2%cmat(i,j,k))
            !                 call simple_stop("image FFT inverse PGI test  stopped")
            !             end if
            !         end do
            !     end do
            ! end do
            print *,'Press enter to execute FFTW inverse FFT'
            read *
            call img%bwd_ft()
            print *,'Press enter to execute CUFFT inverse FFT'
            read *
            call img_2%ifft_pgi_cuda()

            call gnufor_image(img%rmat(:,:,1), palette='gray')                               ! Show img_2
            call gnufor_image(img_2%rmat(:,:,1), palette='gray')
            print *," L1 norm FFTW and CUFFT ", img.lonesum.img_2
            print *,'Press enter to compare FFTW vs FFT'
            read *

            do i=1,ldim(1)
                VerbosePrint "In Image::ifft testing ",i
                do j=1,ldim(2),2
                    do k=1,ldim(3)
                        if (abs( img%rmat(i,j,k) - img_2%rmat(i,j,k)) >= 1.e-06 )then
                            print *,'wrong FFT i/j/k fftw /cufft', i,j,k,img%rmat(i,j,k) , img_2%rmat(i,j,k)
                            call simple_stop("image FFT PGI test stopped")
                        end if
                    end do

                end do
            end do
            VerbosePrint "In Image::ifft testing done"


            write(*,'(a)') '**info(simple_image_unit_test, part 7): testing fftshift'
            if( allocated(pcavec1) ) deallocate(pcavec1)
            if( allocated(pcavec2) ) deallocate(pcavec2)
            passed = .false.
            msk=50
            call img%gauimg(10)
            write(*,'(a)') '**info(simple_image PGI unit_test, part 2):'
            if( doplot ) call img%vis
            write(*,'(a)') '**info(simple_image PGI unit_test, part 2):'
            call img%serialize(pcavec1, msk)
            call img%shift([-9.345,-5.786,0.])
            if( doplot ) call img%vis
            call img%shift([9.345,5.786,0.])
            call img%serialize(pcavec2, msk)
            if( doplot ) call img%vis
            if( pearsn(pcavec1, pcavec2) > 0.99 ) passed = .true.
            if( .not. passed )  call simple_stop('origin shift test failed')

        end subroutine test_image_pgi

    end subroutine test_image_pgi_cuda

    !>  \brief  is a destructor
    subroutine kill( self )
        class(image), intent(inout) :: self
        if( self%existence )then
            call fftwf_free(self%p)
            self%rmat=>null()
            self%cmat=>null()
            call fftwf_destroy_plan(self%plan_fwd)
            call fftwf_destroy_plan(self%plan_bwd)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_image

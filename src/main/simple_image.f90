! the abstract image data type and its methods. 2D/3D & FT/real all implemented by this class
! and Fourier transformations done in-place to reduce memory usage
module simple_image
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ftiter,  only: ftiter
use simple_imgfile, only: imgfile
use simple_winfuns, only: winfuns
use simple_fftw3
use gnufor2
implicit none
private
public :: image, test_image
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
    procedure, private :: disc_1
    procedure, private :: disc_2
    generic            :: disc => disc_1, disc_2
    procedure          :: ring
    procedure          :: copy
    procedure          :: mic2spec
    procedure          :: mic2eospecs
    procedure          :: boxconv
    procedure          :: window
    procedure          :: window_slim
    procedure          :: add_window
    procedure          :: win2arr
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
    procedure          :: get
    procedure          :: get_rmat
    procedure          :: get_rmat_ptr
    procedure          :: get_rmat_sub
    procedure          :: get_cmat
    procedure          :: get_cmat_sub
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
    procedure          :: get_lfny
    procedure          :: get_lhp
    procedure          :: get_lp
    procedure          :: get_spat_freq
    procedure          :: get_find
    procedure          :: get_clin_lims
    procedure          :: rmat_associated
    procedure          :: cmat_associated
    procedure          :: serialize
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
    procedure, private :: add_1
    procedure, private :: add_2
    procedure, private :: add_3
    procedure, private :: add_4
    procedure, private :: add_5
    generic            :: add => add_1, add_2, add_3, add_4, add_5
    procedure, private :: add_workshare_1
    procedure, private :: add_workshare_2
    procedure, private :: add_workshare_3
    generic            :: add_workshare => add_workshare_1, add_workshare_2, add_workshare_3
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
    ! BINARY IMAGE METHODS
    procedure          :: nforeground
    procedure          :: nbackground
    procedure, private :: bin_1
    procedure, private :: bin_2
    generic            :: bin => bin_1, bin_2
    procedure          :: bin_kmeans
    procedure          :: cendist
    procedure          :: masscen
    procedure          :: calc_shiftcen
    procedure          :: calc_shiftcen_serial
    procedure          :: bin_inv
    procedure          :: grow_bin
    procedure          :: grow_bins
    procedure          :: shrink_bin
    procedure          :: shrink_bins
    procedure          :: cos_edge
    procedure          :: remove_edge
    procedure          :: increment
    procedure          :: bin2logical
    procedure          :: collage
    procedure, private :: enumerate_white_pixels
    procedure, private :: enumerate_connected_comps
    procedure          :: find_connected_comps
    procedure          :: size_connected_comps
    procedure          :: is_cc_closed                !!!!!!!!!ADDED BY CHIARA
    procedure          :: prepare_connected_comps
    procedure          :: elim_cc                     !!!!!!!!!ADDED BY CHIARA
    procedure          :: dilatation                  !!!!!!!!!ADDED BY CHIARA
    procedure          :: erosion                     !!!!!!!!!ADDED BY CHIARA
    procedure          :: morpho_closing              !!!!!!!!!ADDED BY CHIARA
    procedure          :: morpho_opening              !!!!!!!!!ADDED BY CHIARA
    procedure          :: border_mask                 !!!!!!!!!ADDED BY CHIARA
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
    procedure          :: tophat
    procedure, private :: apply_filter_1
    procedure, private :: apply_filter_2
    generic            :: apply_filter => apply_filter_1, apply_filter_2
    procedure          :: apply_filter_serial
    procedure          :: phase_rand
    procedure          :: hannw
    procedure          :: real_space_filter
    !procedure          :: NLmean              !!!!!!!!!ADDED BY CHIARA, still to work on
    !procedure          :: similarity_window   !!!!!!!!!ADDED BY CHIARA
    ! CALCULATORS
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
    procedure          :: calc_gradient
    procedure          :: calc_neigh_8_1
    procedure          :: calc_neigh_8_2     !!!!!!!!!!ADDED BY CHIARA
    generic            :: calc_neigh_8   =>  calc_neigh_8_1, calc_neigh_8_2 !!!!!!!!!!ADDED BY CHIARA
    procedure          :: comp_addr_phys1
    procedure          :: comp_addr_phys2
    generic            :: comp_addr_phys =>  comp_addr_phys1, comp_addr_phys2
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
    procedure          :: get_res
    procedure, private :: oshift_1
    procedure, private :: oshift_2
    generic            :: oshift => oshift_1, oshift_2
    procedure, private :: gen_argtransf_comp
    procedure          :: median_value !!!!!!!!!!!!!!ADDED BY CHIARA
    procedure          :: center_edge  !!!!!!!!!!!!!!ADDED BY CHIARA
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
    procedure          :: zero_env_background
    procedure          :: noise_norm_pad_fft
    procedure          :: noise_norm_pad
    procedure          :: salt_n_pepper
    procedure          :: square
    procedure          :: draw_picked  !!!!!ADDED BY CHIARA
    procedure          :: corners
    procedure          :: before_after
    procedure, private :: gauimg_1
    procedure, private :: gauimg_2
    generic            :: gauimg => gauimg_1, gauimg_2
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
    procedure          :: scale_pixels
    procedure          :: mirror
    procedure          :: norm
    procedure, private :: norm4viz
    procedure          :: norm_ext
    procedure          :: noise_norm
    procedure          :: zero_edgeavg
    procedure          :: norm_bin
    procedure          :: roavg
    procedure          :: rtsq
    procedure          :: rtsq_serial
    procedure, private :: shift_phorig
    procedure          :: shift
    procedure, private :: shift2Dserial_1
    procedure, private :: shift2Dserial_2
    generic            :: shift2Dserial => shift2Dserial_1, shift2Dserial_2
    procedure          :: set_within
    procedure          :: ft2img
    procedure          :: img2ft
    procedure          :: cure_outliers
    procedure          :: zero_below
    procedure          :: build_ellipse    !!!!!!!!!!!!!!ADDED BY CHIARA
    ! FFTs
    procedure          :: fft  => fwd_ft
    procedure          :: ifft => bwd_ft
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
            self%rmat = 1.
            lmsk      = .true.
        else where
            self%rmat = 0.
            lmsk      = .false.
        end where
        if( present(npix) ) npix = count(lmsk)
    end subroutine disc_2

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
    subroutine copy( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
        call self%new(self_in%ldim, self_in%smpd, self%wthreads)
        self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
            &self_in%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        self%ft = self_in%ft
    end subroutine copy

    !> mic2spec calculates the average powerspectrum over a micrograph
    function mic2spec( self, box, speckind ) result( img_out )
        class(image),     intent(inout) :: self
        integer,          intent(in)    :: box
        character(len=*), intent(in)    :: speckind
        type(image) :: img_out, tmp, tmp2
        integer     :: xind, yind, cnt
        logical     :: didft
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            THROW_HARD('cannot use a box larger than the image')
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
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images; mic2eoimgs')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            THROW_HARD('cannot use a box larger than the image; mic2eoimgs')
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
                call tmp%zero_edgeavg
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
        call tmp%kill()
        call tmp2%kill()
        if( didft ) call self%fft()
    end subroutine mic2eospecs

    function boxconv( self, box ) result( img_out )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: box
        type(image) :: img_out, tmp
        integer     :: xind, yind, cnt
        logical     :: didft
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images; boxconv')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            THROW_HARD('cannot boxconv using a box larger than the image')
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
    subroutine window( self_in, coord, box, self_out, noutside )
        class(image),      intent(in)    :: self_in
        integer,           intent(in)    :: coord(2), box
        class(image),      intent(inout) :: self_out
        integer, optional, intent(inout) :: noutside
        integer :: fromc(2), toc(2), xoshoot, yoshoot, xushoot, yushoot, xboxrange(2), yboxrange(2)
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
    end subroutine window

    !>  window_slim  extracts a particle image from a box as defined by EMAN 1.9
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
        allocate(pixels(npix), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: win2arr; simple_image',alloc_stat)
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
                print *, 'file: ', trim(fname)
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
            THROW_HARD('image need to be constructed before read/write; open')
        endif
    end subroutine open

    !>  \brief read: for reading 2D images from stack or volumes from volume files
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
            THROW_HARD('unsupported file format; read')
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
                if( ii .gt. 1 ) THROW_HARD('stacks of volumes not supported; read')
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
                    THROW_HARD('logical dims of overall header & image object do not match; read')
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
                    THROW_HARD('optional argument i required for reading from stack; read')
                endif
                if( self%ldim(3) == 1 )then
                    ! all good
                else if( self%ldim(3) > 1 )then
                    THROW_HARD('trying to read from a stack into a volume; read')
                else
                    THROW_HARD('nonconforming logical dimension of image; read')
                endif
                if( iform == -11 .or. iform == -12 ) self%ft = .true.
            case(3,-21,-22)
                ! we are processing a 3D image (stacks of 3D volumes not allowed in SIMPLE)
                if( present(i) )then
                    THROW_HARD('stacks of 3D volumes not allowed in SIMPLE; read')
                endif
                if( self%ldim (3) > 1 )then
                    ! all good
                else if( self%ldim(3) == 1)then
                    THROW_HARD('trying to read from a volume into a 2D image; read')
                else
                    THROW_HARD('nonconforming logical dimension of image; read')
                endif
                if( iform == -21 .or. iform == -22 ) self%ft = .true.
            case DEFAULT
                write(*,*) 'iform = ', iform
                THROW_HARD('unsupported iform flag; read')
            end select
        end subroutine spider_exception_handler

    end subroutine read

    !>  \brief  for writing any kind of images to stack or volumes to volume files
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
        isvol = .false.
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
                    THROW_HARD('trying to write 3D image to stack ; write')
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
                THROW_HARD('unsupported file format; write')
            end select
            ! work out the slice range
            if( isvol )then
                if( ii .gt. 1 ) THROW_HARD('stacks of volumes not supported; write')
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
            THROW_HARD('nonexisting image cannot be written to disk; write')
        endif
    end subroutine write

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

    function get_rmat( self ) result( rmat )
        class(image), intent(in) :: self
        real, allocatable :: rmat(:,:,:)
        integer :: ldim(3)
        ldim = self%ldim
        allocate(rmat(ldim(1),ldim(2),ldim(3)), source=self%rmat(:ldim(1),:ldim(2),:ldim(3)), stat=alloc_stat)
        if (alloc_stat /= 0)call allocchk("simple_image::get_rmat ",alloc_stat)
    end function get_rmat

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

    function get_rmat_at_1( self, logi ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: logi(3)
        real :: val
        val = self%rmat(logi(1),logi(2),logi(3))
    end function get_rmat_at_1

    elemental pure function get_rmat_at_2( self, i,j,k ) result( val )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: i,j,k
        real :: val
        val = self%rmat(i,j,k)
    end function get_rmat_at_2

    !>  \brief   get_cmat get the image object's complex matrix
    pure function get_cmat( self ) result( cmat )
        class(image), intent(in) :: self
        complex, allocatable :: cmat(:,:,:)
        allocate(cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3)), source=self%cmat)
    end function get_cmat

    !>  \brief   get_cmat get the image object's complex matrix
    pure subroutine get_cmat_sub( self, cmat )
        class(image), intent(in)  :: self
        complex,      intent(out) :: cmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        cmat=self%cmat
    end subroutine get_cmat_sub

    !>  \brief   get_cmatfull get the image object's complex matrix, including index limits
    function get_cmatfull( self , llims) result( cmat )
        class(image), intent(in) :: self
        integer, intent(inout), optional :: llims(3,2)
        complex, allocatable :: cmat(:,:,:)
        integer :: llimits(3,2)
        if (present(llims))then
            llimits=llims
        else
            llimits      = 1
            llimits(:,2) = self%array_shape
        end if
        allocate(cmat(llimits(1,1):llimits(1,2),&
            &llimits(2,1):llimits(2,2),&
            &llimits(3,1):llimits(3,2)), source=self%cmat, stat=alloc_stat)
        if (alloc_stat /= 0)call allocchk("simple_image::get_cmat ",alloc_stat)
    end function get_cmatfull

    pure function get_cmat_at_1( self, phys ) result( comp )
        class(image), intent(in)  :: self
        integer,      intent(in)  :: phys(3)
        complex :: comp
        comp = self%cmat(phys(1),phys(2),phys(3))
    end function get_cmat_at_1

    elemental pure function get_cmat_at_2( self, h,k,l ) result( comp )
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

    ! set comp to cmat at index phys
    subroutine set_cmat_at_1( self , phys , comp)
        class(image), intent(inout) :: self
        integer,      intent(in) :: phys(3)
        complex,      intent(in) :: comp
        self%cmat(phys(1),phys(2),phys(3)) = comp
    end subroutine set_cmat_at_1

    !> add comp to cmat at index phys
    elemental subroutine set_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer, intent(in) :: h,k,l
        complex, intent(in) :: comp
        self%cmat(h,k,l) =  comp
    end subroutine set_cmat_at_2

    !> add comp to cmat at index phys
    subroutine add_cmat_at_1( self , phys , comp)
        class(image), intent(inout) :: self
        integer,      intent(in) :: phys(3)
        complex,      intent(in) :: comp
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) + comp
    end subroutine add_cmat_at_1

    !> add comp to cmat at index (h,k,l)
    elemental subroutine add_cmat_at_2( self, h, k, l, comp)
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

    elemental subroutine div_cmat_at_2( self,h,k,l, rval)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        real,         intent(in) :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / cmplx(rval,0.)
        else
            self%cmat(h,k,l) =cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_2

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

    elemental subroutine div_cmat_at_4( self,h,k,l, cval)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        complex,      intent(in) :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / cval
        else
            self%cmat(h,k,l) =cmplx(0.,0.)
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

    elemental subroutine mul_cmat_at_3( self, h,k,l,rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        real,         intent(in)    :: rval
        self%cmat(h,k,l) = self%cmat(h,k,l) * rval
    end subroutine mul_cmat_at_3

    elemental subroutine mul_cmat_at_4( self, h,k,l, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
        self%cmat(h,k,l) = self%cmat(h,k,l) * cval
    end subroutine mul_cmat_at_4

    subroutine print_cmat( self )
        class(image), intent(in) :: self
        print *, self%cmat
    end subroutine print_cmat

    subroutine print_rmat( self )
        class(image), intent(in) :: self
        print *, self%rmat
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
            THROW_HARD('nonconforming dims; set_rmat')
        endif
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
            write(*,*) 'dim(cmat): ', cdim
            write(*,*) 'dim(self%cmat): ', self%array_shape
            THROW_HARD('nonconforming dims; set_cmat')
        endif
    end subroutine set_cmat

    !> \brief  set_ldim replace image dimensions new 3D size
    subroutine set_ldim( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        self%ldim = ldim
    end subroutine set_ldim

    !>  \brief set_smpd for setting smpd
    subroutine set_smpd( self, smpd )
        class(image), intent(inout) :: self
        real,         intent(in)    :: smpd
        self%smpd = smpd
    end subroutine set_smpd

    !> \brief get_slice is for getting a slice from a volume
    function get_slice( self3d, slice ) result( self2d )
        class(image), intent(in) :: self3d
        integer,      intent(in) :: slice
        type(image)              :: self2d
        call self2d%new([self3d%ldim(1),self3d%ldim(2),1],self3d%smpd)
        self2d%rmat(:,:,1) = self3d%rmat(:,:,slice)
    end function get_slice

    !>  \brief set_slice is for putting a slice into a volume
    subroutine set_slice( self3d, slice, self2d )
        class(image), intent(in)    :: self2d
        integer,      intent(in)    :: slice
        class(image), intent(inout) :: self3d
        self3d%rmat(:,:,slice) = self2d%rmat(:,:,1)
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

    function get_clin_lims( self, lp_dyn ) result( lims )
        class(image), intent(in) :: self
        real,         intent(in) :: lp_dyn
        integer                  :: lims(2)
        lims = self%fit%get_clin_lims(lp_dyn)
    end function get_clin_lims

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

    function serialize( self, l_msk )result( pcavec )
        class(image), intent(in) :: self
        logical,      intent(in) :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable :: pcavec(:)
        integer :: sz, cnt, i, j, k
        sz = count(l_msk)
        allocate(pcavec(sz))
        cnt = 0
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if( l_msk(i,j,k) )then
                        cnt         = cnt + 1
                        pcavec(cnt) = self%rmat(i,j,k)
                    endif
                end do
            end do
        end do
    end function serialize

    subroutine unserialize( self, l_msk, pcavec )
        class(image),      intent(inout) :: self
        logical,           intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable, intent(in)    :: pcavec(:)
        !real    :: field(self%ldim(1),self%ldim(2),self%ldim(3))
        integer :: sz, sz_msk, i, j, k, cnt
        if( allocated(pcavec) )then
            sz     = size(pcavec)
            sz_msk = count(l_msk)
            if( sz /= sz_msk )then
                write(*,*) 'ERROR! Nonconforming sizes'
                write(*,*) 'sizeof(pcavec): ', sz
                write(*,*) 'sizeof(l_msk) : ', sz_msk
                THROW_HARD('unserialize')
            endif
        else
            THROW_HARD('pcavec unallocated; unserialize')
        endif
        if( self%ft ) self%ft = .false.
        self%rmat = 0.
        cnt = 0
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                  if( l_msk(i,j,k) )then
                      cnt = cnt + 1
                      self%rmat(i,j,k) =  pcavec(cnt)
                  endif
                end do
            end do
        end do
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

    !>  \brief get_fcomp for getting a Fourier component from the compact representation
    elemental complex function get_fcomp2D(self, h, k)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k
        integer :: phys1, phys2
        if (h .ge. 0) then
            phys1 = h + 1
            phys2 = k + 1 + merge(self%ldim(2),0, k<0)
            ! phys3 = l + 1 + merge(self%ldim(3),0, l<0)
            get_fcomp2D = self%cmat(phys1,phys2,1)
        else
            phys1 = -h + 1
            phys2 = -k + 1 + merge(self%ldim(2),0, -k<0)
            ! phys3 = -l + 1 + merge(self%ldim(3),0, -l<0)
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
    pure function is_ft( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%ft
    end function is_ft

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

    !>  \brief  is for image addition(+) addition
    !! \param self1 image object 1
    !! \param self2  image object 2
    !! \return lhs, copy of added images
    !!
    function sq_root( self1 ) result( self )
        class(image), intent(in) :: self1
        type(image) :: self
        call  self%new(self1%ldim, self1%smpd)
        self%rmat = self1%rmat
        self%rmat = sqrt(self%rmat)
        self%ft = self1%ft
    end function sq_root

    !>  \brief  l1norm_1 is for l1 norm calculation
    function l1norm_1( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .neqv. self2%ft )then
                THROW_HARD('cannot process images of different FT state; l1norm_1')
            endif
            if( self1%ft )then
                self%cmat = cabs(self1%cmat-self2%cmat)
            else
                self%rmat = abs(self1%rmat-self2%rmat)
            endif
        else
            THROW_HARD('cannot process images of different dims; l1norm_1')
        endif
    end function l1norm_1

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

    !>  \brief add_4 is for componentwise weighted summation with kernel division, not overloaded
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

    subroutine add_workshare_3( self, self_to_add )
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
    end subroutine add_workshare_3

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

    !>  \brief subtr_3 is for componentwise weighted subtraction with kernel division, not overloaded
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
                if( square )then
                    self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/(k(phys(1),phys(2),phys(3))**2.)
                else
                    self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/k(phys(1),phys(2),phys(3))
                endif
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
    function conjugate( self ) result ( self_out )
        class(image), intent(in) :: self
        type(image) :: self_out
        if( self%ft )then
            call self_out%copy(self)
            self%cmat = conjg(self%cmat)
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

    !>  \brief  is for binarizing an image with given threshold value
    !!          binary normalization (norm_bin) assumed!> bin_1
    subroutine bin_1( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        if( self%ft ) THROW_HARD('only for real images; bin_1')
        where( self%rmat >= thres )
            self%rmat = 1.
        elsewhere
            self%rmat = 0.
        end where
    end subroutine bin_1

    !>  \brief  bin_2 is for binarizing an image using nr of pixels/voxels threshold
    subroutine bin_2( self, npix )
        class(image), intent(inout) :: self
        integer, intent(in)         :: npix
        real, allocatable           :: forsort(:)
        real                        :: thres
        integer                     :: npixtot
        if( self%ft ) THROW_HARD('only for real images; bin_2')
        npixtot = product(self%ldim)
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
        if( self%ft ) THROW_HARD('only for real images; bin_kmeans')
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

    !>  \brief cendist produces an image with square distance from the centre of the image
    subroutine cendist( self )
        class(image), intent(inout) :: self
        real    :: centre(3)
        integer :: i
        if( self%ft ) THROW_HARD('real space only; cendist')
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

    !> to find centre of gravity
    subroutine masscen( self, xyz )
        class(image), intent(inout) :: self
        real        , intent(out)   :: xyz(3)
        real    ::  spix, ci, cj, ck
        integer :: i, j, k
        if( self%is_ft() ) THROW_HARD('masscen not implemented for FTs; masscen')
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
    end subroutine masscen

    !>  \brief is for estimating the center of an image based on center of mass
    function calc_shiftcen( self, lp, msk ) result( xyz )
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
        call tmp%mask(rmsk, 'hard')
        ! such that norm_bin will neglect everything < 0. and preserve zero
        where(tmp%rmat < TINY) tmp%rmat=0.
        call tmp%norm_bin
        call tmp%masscen(xyz)
        call tmp%kill
    end function calc_shiftcen

    function calc_shiftcen_serial( self, lp, msk ) result( xyz )
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
        call thread_safe_tmp_imgs(ithr)%mask(msk, 'hard')
        where(thread_safe_tmp_imgs(ithr)%rmat < TINY) thread_safe_tmp_imgs(ithr)%rmat = 0.
        call thread_safe_tmp_imgs(ithr)%norm_bin
        call thread_safe_tmp_imgs(ithr)%masscen(xyz)
    end function calc_shiftcen_serial

    !>  \brief bin_inv inverts a binary image
    subroutine bin_inv( self )
        class(image), intent(inout) :: self
        self%rmat = -1.*(self%rmat-1.)
    end subroutine bin_inv

    !>  \brief grow_bin adds one layer of pixels bordering the background in a binary image
    subroutine grow_bin( self )
        class(image), intent(inout) :: self
        integer                     :: i,j,k
        integer                     :: il,ir,jl,jr,kl,kr
        logical, allocatable        :: add_pixels(:,:,:)
        if( self%ft ) THROW_HARD('only for real images; grow_bin')
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
    subroutine shrink_bin( self )
        class(image), intent(inout) :: self
        integer                     :: i,j,k
        integer                     :: il,ir,jl,jr,kl,kr
        logical, allocatable        :: sub_pixels(:,:,:)
        if( self%ft ) THROW_HARD('only for real images; shrink_bin')
        allocate( sub_pixels(self%ldim(1),self%ldim(2),self%ldim(3)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('shrink_bin; simple_image')
        ! Figure out which pixels to remove
        sub_pixels = .false.
        if( self%ldim(3) == 1 )then
            do i=1,self%ldim(1)
                il = max(1,i-1)
                ir = min(self%ldim(1),i+1)
                do j=1,self%ldim(2)
                    if ( is_zero(self%rmat(i,j,1)) ) then
                        jl = max(1,j-1)
                        jr = min(self%ldim(2),j+1)
                        if( any( is_equal(self%rmat(il:ir,jl:jr,1), 1.) ) ) &
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
                        if ( is_zero(self%rmat(i,j,k)) ) then
                            kl = max(1,k-1)
                            kr = min(self%ldim(3),k+1)
                            if( any(is_equal(self%rmat(il:ir,jl:jr,kl:kr), 1.)) ) &
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
    subroutine grow_bins( self, nlayers )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nlayers
        integer                     :: i,j,k, tsz(3,2), win(3,2), pdsz(3,2)
        logical, allocatable        :: add_pixels(:,:,:), template(:,:,:)
        if( self%ft ) THROW_HARD('only for real images; grow_bins')
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
        forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), is_equal(self%rmat(i,j,k),1.) )&
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
    subroutine shrink_bins( self, nlayers )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nlayers
        integer                     :: i,j,k, tsz(3,2), win(3,2), pdsz(3,2)
        logical, allocatable        :: sub_pixels(:,:,:), template(:,:,:)
        if( self%ft ) THROW_HARD('only for real images; shrink_bins')
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

        forall( i=1:self%ldim(1), j=1:self%ldim(2), k=1:self%ldim(3), is_equal(self%rmat(i,j,k),1.) )&
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

    !>  \brief cos_edge applies cosine squared edge to a binary image
    subroutine cos_edge( self, falloff )
        class(image), intent(inout) :: self
        integer, intent(in)         :: falloff
        real, allocatable           :: rmat(:,:,:)
        real                        :: rfalloff, scalefactor
        integer                     :: i, j, k, is, js, ks, ie, je, ke
        integer                     :: il, ir, jl, jr, kl, kr, falloff_sq
        if( falloff<=0 ) THROW_HARD('stictly positive values for edge fall-off allowed; cos_edge')
        if( self%ft )    THROW_HARD('not intended for FTs; cos_edge')
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
            if( .not. any(is_equal(rmat(i,:,:),1.)) )cycle ! no values equal to one
            do j=1,self%ldim(2)
                js = max(1,j-1)
                je = min(j+1,self%ldim(2))
                jl = max(1,j-falloff)
                jr = min(j+falloff,self%ldim(2))
                if( self%ldim(3)==1 )then
                    ! 2d
                    if( is_equal(rmat(i,j,1), 1.) ) cycle ! cycle if not equal to one
                    ! within mask region
                    ! update if has a masked neighbour
                    if( any( rmat(is:ie,js:je,1) < 1.) )call update_mask_2d
                else
                    ! 3d
                    if(.not. any(is_equal(rmat(i,j,:), 1.)) )cycle ! cycle if equal to one
                    do k=1,self%ldim(3)
                        if(.not. is_equal(rmat(i,j,k) , 1.) )cycle
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
        if( self%ft ) THROW_HARD('only for real binary images (not FTed ones); remove_edge')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input to remove edge not binary; remove_edge')
        where( self%rmat < 0.999 ) self%rmat = 0.
    end subroutine remove_edge

    !>  \brief  increments the logi pixel value with incr
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

    ! This subroutine creates the initial label map of a binary image in order
    ! to identify the connected components of the image
    subroutine enumerate_white_pixels(self, label_matrix)
        class(image),      intent(in)   :: self
        real, allocatable, intent(out)  :: label_matrix(:,:,:)
        integer           :: i, j, cnt
        if( allocated(label_matrix) ) deallocate(label_matrix)
        allocate(label_matrix(self%ldim(1),self%ldim(2),1), source = 0.)
        cnt = 0 ! # labels
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                if( self%rmat(i,j,1) > 0.5 )then
                    cnt = cnt + 1
                    label_matrix(i,j,1) = cnt
                endif
            enddo
        enddo
    end subroutine enumerate_white_pixels


    !!!!!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!
    !This subroutine enumerates (increasing order) the connected components
    !of the image. This connected components are selected by
    !'find_connected_comps' procedure.
    subroutine enumerate_connected_comps(self_in, self_out)
        class(image), intent(inout) :: self_in
        class(image), intent(out)   :: self_out
        integer :: i, j, cnt
        real    :: tmp
        if(self_out%existence) call self_out%kill()  !Reset if present
        call self_out%new(self_in%ldim,self_in%smpd)
        if(self_in%ldim(3) /= 1) THROW_HARD('enumerate_connected_comps. Only for 2D images!')
        cnt = 0   !it is going to be the new label of the connected component
        do i = 1, self_in%ldim(1)
            do j = 1, self_in%ldim(2)
                if( self_in%rmat(i,j,1) > 0.5 ) then  !rmat == 0  --> background
                    cnt = cnt + 1
                    tmp = self_in%rmat(i,j,1)
                    where(self_in%rmat == tmp)
                        self_out%rmat = cnt
                        self_in%rmat  = 0.            !Not to consider this cc again
                    endwhere
                endif
            enddo
        enddo
    end subroutine enumerate_connected_comps

    ! Img_in should be a binary image. Img_out is the connencted component image.
    subroutine find_connected_comps(img_in, img_out)
        class(image), intent(in)  :: img_in
        type(image),  intent(out) :: img_out
        type(image)       :: img_cc   !in this img will be stored the cc with no specific order
        real, allocatable :: label_matrix(:,:,:),mat4compare(:,:,:),neigh_8(:)
        integer           :: i, j, n_it, n_maxit
        logical           :: finished_job
        if(img_out%existence) call img_out%kill()
        finished_job = .false.
        call img_cc%new (img_in%ldim,img_in%smpd)
        allocate(mat4compare(img_in%ldim(1),img_in%ldim(2),1), source = 0.)
        call enumerate_white_pixels(img_in, label_matrix)
        call img_cc%set_rmat(label_matrix)
        n_maxit = maxval(label_matrix)
        do n_it = 1, n_maxit
            if( all(abs(mat4compare - img_cc%rmat(:img_cc%ldim(1),:img_cc%ldim(2),:img_cc%ldim(3))) < TINY) )then
                finished_job = .true.! the matrices are identical
            endif
            if(.not. finished_job) then
                mat4compare = img_cc%rmat
                do i = 1, img_in%ldim(1)
                    do j = 1, img_in%ldim(2)
                        if(img_in%rmat(i,j,1) > 0.5) then  !not background
                            call img_cc%calc_neigh_8([i,j,1], neigh_8)
                            img_cc%rmat(i,j,1) = minval(neigh_8, neigh_8 > 0.5)
                        endif
                    enddo
                enddo
            endif
        enddo
        call img_cc%enumerate_connected_comps(img_out)
        deallocate(label_matrix, mat4compare)
    end subroutine find_connected_comps

    ! The result of the function is the size(# of pixels) of each cc. This
    ! value is stored in the 2nd column of sz. In the first one is recorded
    ! the label of the cc.  (cc = connected component)
    function size_connected_comps(self) result(sz)
        class(image), intent(in) :: self
        integer, allocatable :: sz(:,:)
        integer :: n_cc
        if(allocated(sz)) deallocate(sz)
        allocate(sz(2,int(maxval(self%rmat))), source = 0)
        do n_cc = int(maxval(self%rmat)),1,-1
            sz(1, n_cc) = n_cc
            sz(2, n_cc) = count(self%rmat == n_cc)
        enddo
    end function size_connected_comps

    ! This function takes in input a connected component (cc) image and the label
    ! and the label of one of its cc. It returns true if the cc is closed, false
    ! otherwise.
    function is_cc_closed(self,label) result(yes_no)
        class(image), intent(inout) :: self
        integer, intent(in) :: label
        logical :: yes_no
        integer :: i, j, n
        real, allocatable :: neigh_8(:)
        yes_no = .true.
        if(any(abs(self%rmat-real(label)) < TINY)) then   !if the label is present
            do i = 1, self%ldim(1)
              do j = 1, self%ldim(2)
                      if(abs(self%get([i,j,1])-real(label)) < TINY) then        !just pixels in the selected cc
                          call self%calc_neigh_8([i,j,1],neigh_8)
                          n = count(abs(neigh_8(:size(neigh_8)-1)-real(label)) < TINY) !do not consider the pixel itself
                          if(n < 2) then
                              yes_no = .false.
                              return
                          endif
                      endif
              enddo
            enddo
        else
            write(*,*) 'label ', label, 'cc selected not present; is_cc_closed'
            yes_no = .false.
            return
        endif
    end function is_cc_closed

    ! This function takes in input a connected component image and modifies
    ! it in order to prepare the centering process.
    ! Notation: cc = connected component.
    subroutine prepare_connected_comps(self, discard, min_sz)
        class(image),      intent(inout) :: self     !image which contains connected components
        logical, optional, intent(out)   :: discard
        integer, optional, intent(in)    :: min_sz
        integer, allocatable :: sz(:,:), biggest_cc(:), biggest_val(:)
        if(present(discard) .and. .not. present(min_sz)) THROW_HARD('Need min_sz;  prepare_connected_comps')
        if(present(min_sz) .and. .not. present(discard)) THROW_HARD('Need discard; prepare_connected_comps')
        sz = self%size_connected_comps()
        biggest_cc  = maxloc(sz(2,:))
        biggest_val = real(sz(1, biggest_cc))
        if(present(discard)) then
            if( sz(2,biggest_cc(1)) < min_sz ) then
                discard = .true.  !if the biggest cc is smaller than min_sz, discard image
            else
                discard = .false.
            endif
        endif
        where( abs(self%rmat - real(biggest_val(1))) > TINY)  !keep just the biggest cc
            self%rmat = 0.
        elsewhere
            self%rmat = 1.   !self is now binary
        endwhere
        deallocate(sz,biggest_cc,biggest_val)
    end subroutine prepare_connected_comps

    !!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!
    ! This subroutine takes in input a connected component (cc) image
    ! and sets to 0 the cc which has size (# pixels) smaller than min_sz = range(1)
    ! or bigger than max_sz = range(2).
    ! It is created in order to prepare a micrograph for picking particles.
    subroutine elim_cc(self, range)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: range(2)
        integer, allocatable :: sz(:,:)
        integer              :: n_cc
        real                 :: label
        sz = self%size_connected_comps()
        do n_cc = 1, size(sz, dim = 2)  !for each cc
            if(sz(2,n_cc) < range(1) .or. sz(2,n_cc) > range(2)) then  !if the cc has size < min_sz or > max_sz
                label = real(sz(1,n_cc))  !label of the cc
                where(abs(self%rmat - label) < TINY)  !rmat == label
                    self%rmat = 0.        !set to 0
                endwhere
            endif
        enddo
    end subroutine elim_cc

     !!!!!!!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!
     ! This subroutine is ment for 2D binary images. It implements
     ! the morphological operation dilatation.
    subroutine dilatation(self)
        class(image), intent(inout) :: self
        integer, allocatable :: neigh_8(:,:,:)
        type(image) :: self_copy
        integer     :: i, j, k
        call self_copy%copy(self)
        if(self%ldim(3) /= 1) THROW_HARD('This subroutine is for 2D images!; dilatation')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input for dilatation not binary; dilatation')
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                if(self_copy%rmat(i,j,1) == 1.) then  !just for white pixels
                    call self%calc_neigh_8([i,j,1], neigh_8)
                    do k = 1, size(neigh_8, dim=2)
                        call self%set(neigh_8(1:3,k,1), 1.) !self%rmat(neigh_8) = 1.
                    enddo
                    deallocate(neigh_8)
                endif
            enddo
        enddo
        call self_copy%kill
    end subroutine dilatation

     !!!!!!!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!
     ! This subroutine is ment for 2D binary images. It implements
     ! the morphological operation erosion.
    subroutine erosion(self)
        class(image) :: self
        logical, allocatable :: border(:,:,:)
        if(self%ldim(3) /= 1) THROW_HARD('This subroutine is for 2D images!; erosion')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input for erosion not binary; erosion')
        call self%border_mask(border)
        where(border)
            self%rmat = 0.
        endwhere
    end subroutine erosion

     !!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!
     ! This subroutine implements morphological operation of closing
     ! on the image self.
    subroutine morpho_closing(self)
        class(image), intent(inout) :: self
        if(self%ldim(3) /= 1) THROW_HARD('This subroutine is for 2D images!; morpho_closing')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input for morpho_closing not binary; morpho_closing')
        call self%dilatation()
        call self%erosion()
    end subroutine morpho_closing

     !!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!
     ! This subroutine implements morphological operation of opening
     ! on the image self.
    subroutine morpho_opening(self)
        class(image), intent(inout) :: self
        if(self%ldim(3) /= 1) THROW_HARD('This subroutine is for 2D images!; morpho_opening')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input for morpho_opening not binary; morpho_opening')
        call self%erosion()
        call self%dilatation()
    end subroutine morpho_opening

     !!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!
     ! This subroutine builds the logical array 'border' of the same dims of
     ! the input image self. Border is true in corrispondence of
     ! the border pixels in self. Self is meant to be binary.
    subroutine border_mask(self, border)
        class(image), intent(in)            :: self
        logical, allocatable, intent(inout) :: border(:,:,:)
        real, allocatable :: neigh_8(:)
        integer :: i, j
        if(self%ldim(3)/=1) THROW_HARD('This subroutine is for 2D images!; calc_border')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input to calculate border not binary; calc_border')
        if(allocated(border)) deallocate(border)
        allocate(border(self%ldim(1), self%ldim(2), self%ldim(3)), source = .false.)
        do i = 1,self%ldim(1)
            do j = 1, self%ldim(2)
                if(abs(self%rmat(i,j,1)-1.) < TINY ) then !white pixels
                    call self%calc_neigh_8([i,j,1], neigh_8)
                    if(any(abs(neigh_8)<TINY)) border(i,j,1) = .true.
                    deallocate(neigh_8)
                endif
            enddo
        enddo
    end subroutine border_mask

    ! FILTERS

    !>  \brief  acf calculates the autocorrelation function of an image
    subroutine acf( self )
        class(image), intent(inout) :: self
        if( .not. self%ft )then
            call self%fft()
        endif
        self%cmat = self%cmat*conjg(self%cmat)
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
    function guinier( self ) result( plot )
        class(image), intent(inout) :: self
        real, allocatable :: spec(:), plot(:,:)
        integer           :: lfny, k
        if( .not. self%is_3d() ) THROW_HARD('Only for 3D images; guinier')
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
    !>  keep serial
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
        if(allocated(spec))deallocate(spec)
        allocate( spec(lfny), counts(lfny), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('spectrum; simple_image',alloc_stat)
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
            THROW_HARD('Unsupported spectrum kind; spectrum')
        end select
        if( which .ne. 'count' .and. nnorm )then
            where(counts > 0.)
                spec = spec/counts
            end where
        endif
        if( didft ) call self%ifft()
    end subroutine spectrum

    !> \brief shellnorm for normalising each shell to uniform (=1) power
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
        call self%spectrum('power',expec_pow)
        ! normalise
        if( self%wthreads )then
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
        else
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
        endif
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
        call self%spectrum('power',expec_pow)
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
            !$omp end parallel do
        else
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
        endif
        if( didft ) call self%ifft()
    end subroutine bp

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
                        freq = hyp(real(h),real(k),real(l))
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
                        freq = hyp(real(h),real(k),real(l))
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
    subroutine apply_filter_2( self, filter )
        class(image), intent(inout) :: self, filter
        real    :: fwght
        integer :: phys(3), lims(3,2), h, k, l
        complex :: comp
        if( self.eqdims.filter )then
            if( filter%ft )then
                if( .not. self%ft )then
                    THROW_HARD('image to be filtered is not in the Fourier domain; apply_filter_2')
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
                THROW_HARD('assumed that the inputted filter is in the Fourier domain; apply_filter_2')
            endif
        else
            THROW_HARD('equal dims assumed; apply_filter_2')
        endif
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
        !$omp parallel do collapse(3) default(shared) private(h,k,l,freq,phys,sgn1,sgn2,sgn3)&
        !$omp schedule(static) proc_bind(close)
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
        integer               :: n, i, j, k, cnt, npix
        real                  :: rn, wfun(-winsz:winsz), norm
        type(winfuns)         :: fwin
        character(len=STDLEN) :: wstr
        type(image)           :: img_filt
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
                !$omp parallel do collapse(2) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = stdev(pixels)
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
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = stdev(pixels)
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
                THROW_HARD('unknown filter type; real_space_filter')
            end select
        endif
        call self%copy(img_filt)
        call img_filt%kill()
    end subroutine real_space_filter

    ! !!!!!!!!!ADDED BY CHIARA, for NLmean.
    ! function similarity_window(self,px) result(sw)
    !     class(image), intent(in) :: self
    !     integer,      intent(in) :: px(3)
    !     integer, parameter  :: DIM_SW = 7
    !     real                :: sw_mat(DIM_SW,DIM_SW,1), sw(DIM_SW*DIM_SW)
    !     if(px(3) /= 1) THROW_HARD('similarity window. Image has to be 2D!')
    !     if(px(1) < 1 .or. px(1)+DIM_SW-1 > self%ldim(1) .or. &
    !     &  px(2) < 1 .or. px(2)+DIM_SW-1 > self%ldim(2)) THROW_HARD('similarity window. Padding error!')
    !     sw_mat(:DIM_SW,:DIM_SW,1) = self%rmat(px(1):px(1)+DIM_SW-1, px(2):px(2)+DIM_SW-1,1)  !self is meant to be pad
    !     sw = reshape(sw_mat,[DIM_SW*DIM_SW])
    ! end function similarity_window

    !!!!!!!!ADDED BY CHIARA, still to work on
    ! subroutine NLmean(self)
    !   class(image), intent(inout) :: self
    !   real, allocatable  :: NL(:,:,:), Z(:,:,:)
    !   integer, parameter :: DIM_SW = 7
    !   real               :: sw(DIM_SW*DIM_SW), sw_px(DIM_SW*DIM_SW)
    !   integer            :: i, j, m, n,  pad
    !   integer, parameter :: cfr_box = 10
    !   real, parameter    :: sigma = 5., h =4.*sigma !SIGMA, H..?
    !   type(image)        :: img_p
    !   real               :: w
    !   pad = (DIM_SW-1)/2
    !   call img_p%new([self%ldim(1)+2*pad,self%ldim(2)+2*pad,1],1.) !you can also insert a control so that you pad only if necessary
    !   call self%scale_pixels([1.,500.])
    !   call self%pad(img_p)
    !   allocate(Z(self%ldim(1),self%ldim(2),1), NL(self%ldim(1),self%ldim(2),1), source = 0.)
    !   do m = 1,self%ldim(1)             !fix pixel (m,n)
    !     do n = 1,self%ldim(2)
    !       do i = -cfr_box,cfr_box       !As suggested in the paper, consider a box 21x21
    !           do j = -cfr_box,cfr_box
    !             if(i/=0 .or. j/=0) then !Not to take pixel (m,n)
    !                 if(m+i > 0 .and. m+i <= self%ldim(1) .and. n+j > 0 .and. n+j <= self%ldim(2)) then
    !                     Z(m,n,1) = Z(m,n,1)+ &
    !                     & exp(-((norm2( img_p%similarity_window([m,n,1]) &
    !                                 & - img_p%similarity_window([m+i,n+j,1])))**2./h**2.))
    !                 endif
    !             endif
    !           enddo
    !       enddo
    !     enddo
    !   enddo
    !
    !   do m = 1,self%ldim(1)
    !     do n = 1,self%ldim(2)
    !       sw_px = img_p%similarity_window([m,n,1])
    !       do i = -cfr_box,cfr_box   !As suggested in the paper, I just consider a box 21x21
    !           do j = -cfr_box,cfr_box
    !             if(i/=0 .or. j/=0) then !Not to take pixel (m,n)
    !                 if(m+i > 0 .and. m+i <= self%ldim(1) .and. n+j > 0 .and. n+j <= self%ldim(2)) then
    !                     if(abs(Z(m,n,1)) > TINY) then !do not divide by 0
    !                         sw = img_p%similarity_window([m+i,n+j,1])
    !                         w = (exp(-((norm2(sw_px-sw))**2./h**2.)))/Z(m,n,1)
    !                         NL(m,n,1) = NL(m,n,1)+w*self%rmat(m+i,n+j,1)
    !                     endif
    !                 endif
    !             endif
    !           enddo
    !       enddo
    !     enddo
    !   enddo
    !   call self%set_rmat(NL)
    !   deallocate(Z,NL)
    ! end subroutine NLmean

    ! CALCULATORS

    !> \brief stats  is for providing foreground/background statistics
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
            THROW_HARD('unrecognized parameter: which; stats_1')
        endif
        allocate( pixels(product(self%ldim)), stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('backgr; simple_image')
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
    function minmax( self )result( mm )
        class(image), intent(in) :: self
        real :: mm(2)
        mm(1) = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        mm(2) = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
    end function minmax

    !>  \brief rmsd for calculating the RMSD of a map
    !! \return  dev root mean squared deviation
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
    subroutine calc_gradient(self, grad, Dc, Dr)
        class(image),                intent(inout) :: self
        real, allocatable,           intent(out)   :: grad(:,:,:)  !gradient matrix
        real, allocatable, optional, intent(out)   :: Dc(:,:,:), Dr(:,:,:) ! derivates column and row matrices
        type(image)        :: img_p                         !padded image
        real, allocatable  :: wc(:,:,:), wr(:,:,:)          !row and column Sobel masks
        integer, parameter :: L = 3                         !dimension of the masks
        integer            :: ldim(3)                       !dimension of the image, save just for comfort
        integer            :: i,j,m,n                       !loop indeces
        real, allocatable  :: Ddc(:,:,:),Ddr(:,:,:)         !column and row derivates
        ldim = self%ldim
        if(ldim(3) /= 1) then
            THROW_HARD('image has to be 2D! calc_gradient')
        endif
        allocate( Ddc(ldim(1),ldim(2),1), Ddr(ldim(1),ldim(2),1), grad(ldim(1),ldim(2),1), &
               & wc(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1),wr(-(L-1)/2:(L-1)/2,-(L-1)/2:(L-1)/2,1), source = 0.)
        wc = (1./8.)*reshape([-1,0,1,-2,0,2,-1,0,1],[3,3,1]) ! Sobel masks
        wr = (1./8.)*reshape([-1,-2,-1,0,0,0,1,2,1],[3,3,1])
        call img_p%new([ldim(1)+L-1,ldim(2)+L-1,1],1.)
        call self%pad(img_p)                                 ! padding
        !$omp parallel do collapse(2) default(shared) private(i,j,m,n)&
        !$omp schedule(static) proc_bind(close)
        do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L-1)/2,(L-1)/2
                  do n = -(L-1)/2,(L-1)/2
                      Ddc(i,j,1) = Ddc(i,j,1)+img_p%rmat(i+m+1,j+n+1,1)*wc(m,n,1)
                      Ddr(i,j,1) = Ddr(i,j,1)+img_p%rmat(i+m+1,j+n+1,1)*wr(m,n,1)
                  end do
              end do
          end do
        end do
        !omp end parallel do
        deallocate(wc,wr)
        grad = sqrt(Ddc**2. + Ddr**2.)
        if(present(Dc)) allocate(Dc(ldim(1),ldim(2),1), source = Ddc)
        if(present(Dr)) allocate(Dr(ldim(1),ldim(2),1), source = Ddr)
        deallocate(Ddc,Ddr)
    end subroutine calc_gradient

    ! Returns 8-neighborhoods of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! and the value of the pixel itself (the last one).
    subroutine calc_neigh_8_1(self, px, neigh_8)
        class(image),      intent(in)    :: self
        integer,           intent(in)    :: px(3)
        real, allocatable, intent(inout) :: neigh_8(:)
        integer :: i, j
        if( px(3) /= 1 .or. self%ldim(3) /= 1 ) then
            THROW_HARD('image has to be 2D; calc_neigh_8')
        endif
        i = px(1)
        j = px(2) ! Assumes to have a 2-dim matrix
        ! sanity check
        if( i < 1 .or. i > self%ldim(1) )then
            THROW_HARD('x-coordinate out of bound; calc_neigh_8')
        endif
        if( j < 1 .or. j > self%ldim(2) )then
            THROW_HARD('y-coordinate out of bound; calc_neigh_8')
        endif
        if(allocated(neigh_8)) deallocate(neigh_8)
        ! identify neighborhood
        if( i-1 < 1 .and. j-1 < 1 )then                            ! NW corner
            allocate(neigh_8(4), source = 0.)
            neigh_8(1) = self%rmat(i+1,j,1)
            neigh_8(2) = self%rmat(i+1,j+1,1)
            neigh_8(3) = self%rmat(i,j+1,1)
            neigh_8(4) = self%rmat(i,j,1)
        else if (j+1 > self%ldim(2) .and. i+1 > self%ldim(1)) then ! SE corner
            allocate(neigh_8(4), source = 0.)
            neigh_8(1) = self%rmat(i-1,j,1)
            neigh_8(2) = self%rmat(i-1,j-1,1)
            neigh_8(3) = self%rmat(i,j-1,1)
            neigh_8(4) = self%rmat(i,j,1)
        else if (j-1 < 1  .and. i+1 >self%ldim(1)) then            ! SW corner
            allocate(neigh_8(4), source = 0.)
            neigh_8(3) = self%rmat(i-1,j,1)
            neigh_8(2) = self%rmat(i-1,j+1,1)
            neigh_8(1) = self%rmat(i,j+1,1)
            neigh_8(4) = self%rmat(i,j,1)
        else if (j+1 > self%ldim(2) .and. i-1 < 1) then            ! NE corner
            allocate(neigh_8(4), source = 0.)
            neigh_8(1) = self%rmat(i,j-1,1)
            neigh_8(2) = self%rmat(i+1,j-1,1)
            neigh_8(3) = self%rmat(i+1,j,1)
            neigh_8(4) = self%rmat(i,j,1)
        else if( j-1 < 1 ) then                                    ! N border
            allocate(neigh_8(6), source = 0.)
            neigh_8(5) = self%rmat(i-1,j,1)
            neigh_8(4) = self%rmat(i-1,j+1,1)
            neigh_8(3) = self%rmat(i,j+1,1)
            neigh_8(2) = self%rmat(i+1,j+1,1)
            neigh_8(1) = self%rmat(i+1,j,1)
            neigh_8(6) = self%rmat(i,j,1)
        else if ( j+1 > self%ldim(2) ) then                        ! S border
            allocate(neigh_8(6), source = 0.)
            neigh_8(1) = self%rmat(i-1,j,1)
            neigh_8(2) = self%rmat(i-1,j-1,1)
            neigh_8(3) = self%rmat(i,j-1,1)
            neigh_8(4) = self%rmat(i+1,j-1,1)
            neigh_8(5) = self%rmat(i+1,j,1)
            neigh_8(6) = self%rmat(i,j,1)
        else if ( i-1 < 1 ) then                                   ! W border
            allocate(neigh_8(6), source = 0.)
            neigh_8(1) = self%rmat(i,j-1,1)
            neigh_8(2) = self%rmat(i+1,j-1,1)
            neigh_8(3) = self%rmat(i+1,j,1)
            neigh_8(4) = self%rmat(i+1,j+1,1)
            neigh_8(5) = self%rmat(i,j+1,1)
            neigh_8(6) = self%rmat(i,j,1)
        else if ( i+1 > self%ldim(1) ) then                       ! E border
            allocate(neigh_8(6), source = 0.)
            neigh_8(1) = self%rmat(i,j+1,1)
            neigh_8(2) = self%rmat(i-1,j+1,1)
            neigh_8(3) = self%rmat(i-1,j,1)
            neigh_8(4) = self%rmat(i-1,j-1,1)
            neigh_8(5) = self%rmat(i,j-1,1)
            neigh_8(6) = self%rmat(i,j,1)
        else                                                     ! DEFAULT
            allocate(neigh_8(9), source = 0.)
            neigh_8(1) = self%rmat(i-1,j-1,1)
            neigh_8(2) = self%rmat(i,j-1,1)
            neigh_8(3) = self%rmat(i+1,j-1,1)
            neigh_8(4) = self%rmat(i+1,j,1)
            neigh_8(5) = self%rmat(i+1,j+1,1)
            neigh_8(6) = self%rmat(i,j+1,1)
            neigh_8(7) = self%rmat(i-1,j+1,1)
            neigh_8(8) = self%rmat(i-1,j,1)
            neigh_8(9) = self%rmat(i,j,1)
        endif
    end subroutine calc_neigh_8_1

    !!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!!!
    ! Returns 8-neighborhoods of the pixel position px in self
    ! it returns the pixel INDECES of the 8-neigh in a CLOCKWISE order,
    ! starting from any 8-neigh. It doesn't consider the pixel itself.
    subroutine calc_neigh_8_2(self, px, neigh_8)
        class(image),         intent(in)   :: self
        integer,              intent(in)   :: px(3)
        integer, allocatable, intent(inout):: neigh_8(:,:,:)
        integer :: i, j
        if( px(3) /= 1 .or. self%ldim(3) /= 1 ) then
            THROW_HARD('The image has to be 2D!; calc_neigh_8_2')
        endif
        i = px(1)
        j = px(2)            !Assumes to have a 2-dim matrix
        if(allocated(neigh_8)) deallocate(neigh_8)
        if ( i-1 < 1 .and. j-1 < 1 ) then
            allocate(neigh_8(3,3,1), source = 0)
            neigh_8(1:3,1,1) = [i+1,j,1]
            neigh_8(1:3,2,1) = [i+1,j+1,1]
            neigh_8(1:3,3,1) = [i,j+1,1]
        else if (j+1 > self%ldim(2) .and. i+1 > self%ldim(1)) then
            allocate(neigh_8(3,3,1), source = 0)
            neigh_8(1:3,1,1) = [i-1,j,1]
            neigh_8(1:3,2,1) = [i-1,j-1,1]
            neigh_8(1:3,3,1) = [i,j-1,1]
        else if (j-1 < 1  .and. i+1 >self%ldim(1)) then
            allocate(neigh_8(3,3,1), source = 0)
            neigh_8(1:3,3,1) = [i-1,j,1]
            neigh_8(1:3,2,1) = [i-1,j+1,1]
            neigh_8(1:3,1,1) = [i,j+1,1]
        else if (j+1 > self%ldim(2) .and. i-1 < 1) then
            allocate(neigh_8(3,3,1), source = 0)
            neigh_8(1:3,1,1) = [i,j-1,1]
            neigh_8(1:3,2,1) = [i+1,j-1,1]
            neigh_8(1:3,3,1) = [i+1,j,1]
        else if( j-1 < 1 ) then
            allocate(neigh_8(3,5,1), source = 0)
            neigh_8(1:3,5,1) = [i-1,j,1]
            neigh_8(1:3,4,1) = [i-1,j+1,1]
            neigh_8(1:3,3,1) = [i,j+1,1]
            neigh_8(1:3,2,1) = [i+1,j+1,1]
            neigh_8(1:3,1,1) = [i+1,j,1]
        else if ( j+1 > self%ldim(2) ) then
            allocate(neigh_8(3,5,1), source = 0)
            neigh_8(1:3,1,1) = [i-1,j,1]
            neigh_8(1:3,2,1) = [i-1,j-1,1]
            neigh_8(1:3,3,1) = [i,j-1,1]
            neigh_8(1:3,4,1) = [i+1,j-1,1]
            neigh_8(1:3,5,1) = [i+1,j,1]
        else if ( i-1 < 1 ) then
            allocate(neigh_8(3,5,1), source = 0)
            neigh_8(1:3,1,1) = [i,j-1,1]
            neigh_8(1:3,2,1) = [i+1,j-1,1]
            neigh_8(1:3,3,1) = [i+1,j,1]
            neigh_8(1:3,4,1) = [i+1,j+1,1]
            neigh_8(1:3,5,1) = [i,j+1,1]
        else if ( i+1 > self%ldim(1) ) then
            allocate(neigh_8(3,5,1), source = 0)
            neigh_8(1:3,1,1) = [i,j+1,1]
            neigh_8(1:3,2,1) = [i-1,j+1,1]
            neigh_8(1:3,3,1) = [i-1,j,1]
            neigh_8(1:3,4,1) = [i-1,j-1,1]
            neigh_8(1:3,5,1) = [i,j-1,1]
        else
            allocate(neigh_8(3,8,1), source = 0)
            neigh_8(1:3,1,1) = [i-1,j-1,1]
            neigh_8(1:3,2,1) = [i,j-1,1]
            neigh_8(1:3,3,1) = [i+1,j-1,1]
            neigh_8(1:3,4,1) = [i+1,j,1]
            neigh_8(1:3,5,1) = [i+1,j+1,1]
            neigh_8(1:3,6,1) = [i,j+1,1]
            neigh_8(1:3,7,1) = [i-1,j+1,1]
            neigh_8(1:3,8,1) = [i-1,j,1]
        endif
    end subroutine calc_neigh_8_2

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
    function corr( self1, self2, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self1, self2
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real    :: r, sumasq, sumbsq
        integer :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        logical :: didft1, didft2
        r = 0.
        sumasq = 0.
        sumbsq = 0.
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
    !!
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
            if( sxx > 0. .and. syy > 0. )then
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
            if( sxx > 0. .and. syy > 0. )then
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
        if( sxx > 0. .and. syy > 0. )then
            r = sxy / sqrt(sxx * syy)
        else
            r = 0.
        endif
    end function real_corr_2

    !> \brief prenorm4real_corr pre-normalises the reference in preparation for real_corr_prenorm
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
                    sh = nint(hyp(real(h),real(k),real(l)))
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
        corrs = corrs_8
    end subroutine fsc

    !>  \brief get array of resolution steps
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
    function oshift_2( self, logi, shvec ) result( comp )
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

   !!!!!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!
   ! The result of this funcion is the median value of the intensities
   ! of the pixels in the image self.
    function median_value(self) result(m)
        class(image), intent(in) :: self
        real              :: m
        real, allocatable :: pixels(:), pixels_nodup(:)
        pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
            & mask = .true.)
        call elim_dup(pixels,pixels_nodup)
        m = median(pixels_nodup)  !or median_nocopy?
        deallocate(pixels,pixels_nodup)
    end function median_value

   !!!!!!!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!
   ! This function has as output the coordinates of the mass center
   ! of img_in, calculated on the basis of the biggest connected
   ! component (cc) in img_in. If this cc is too small, than
   ! discard = .true.
    function center_edge( self, discard ) result( xyz )
        class(image), intent(inout) :: self
        logical,        intent(out) :: discard
        type(image) :: imgcc
        real        :: xyz(3)
        integer     :: min_sz !minumum size of the biggest cc in img_in
        min_sz = 20 !too small cc have already been eliminated, this control
        ! is in order to discard too offcentered windows, in which it wouldn t
        ! be possible to properly center the particle
        ! find all the cc in img_in
        call self%find_connected_comps(imgcc)
        !find the biggest cc
        call imgcc%prepare_connected_comps(discard, min_sz)
        call imgcc%masscen(xyz)
    end function center_edge

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
        if( .not. self%ft        ) THROW_HARD('image to be modified assumed to be FTed; insert_lowres')
        if( .not. self2insert%ft ) THROW_HARD('image to insert assumed to be FTed; insert_lowres')
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
    end subroutine zero_background

    ! substracts median of background defined by an enveloppe mask (eg 0< <1)
    subroutine zero_env_background( self, volmsk )
        class(image), intent(inout) :: self, volmsk
        real, allocatable :: vals(:)
        real              :: med
        integer           :: cnt, npix, npix_env, i,j,k
        logical           :: env_mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        npix = product(self%ldim)
        where( volmsk%rmat>0.0001 .and. volmsk%rmat<0.9999 )
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

    subroutine noise_norm_pad_fft( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        call self%noise_norm(lmsk)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine noise_norm_pad_fft

    subroutine noise_norm_pad( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        call self%noise_norm(lmsk)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine noise_norm_pad

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
    subroutine salt_n_pepper( self, pos )
        class(image), intent(inout) :: self
        logical, intent(in)         :: pos(:,:)
        integer :: ipix, jpix
        if( .not. self%is_2d() ) THROW_HARD('only for 2D images; salt_n_pepper')
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

    !!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!!
    ! This functions draws a window in self centered in part_coords.
    ! Required inputs: coordinates of the particle, radius of the particle.
    ! This function is meant to mark picked particles (visual purpose).
    subroutine draw_picked(self, part_coords, part_radius, border)
      class(image),intent(inout)    :: self
      integer,     intent(in)       :: part_coords(2)  !coordinates of the picked particle
      integer,     intent(in)       :: part_radius
      integer, optional, intent(in) :: border          !width of the border of the drawn line
      integer :: bborder
      real    :: value             !intensity value of the window
      integer :: wide              !side width of the window
      integer :: length            !length of the drawn side
      bborder = 1
      if(present(border)) bborder = border
      value = minval(self%rmat(:,:,:)) !maxval(self%rmat(:,:,:))
      wide = 4*part_radius
      length = int(part_radius/2)
      if( .not. self%is_2d() ) THROW_HARD('only for 2D images; draw_picked')
      if(part_coords(1)-wide/2-int((bborder-1)/2) < 1 .or. part_coords(1)+wide/2+int((bborder-1)/2) > self%ldim(1) .or. &
      &  part_coords(2)-wide/2-int((bborder-1)/2) < 1 .or. part_coords(2)+wide/2+int((bborder-1)/2) > self%ldim(2) ) then
        print *, 'The window is out of the border of the image!'
        return    !Do not throw error, just do not draw
      endif
      ! Edges of the window
      ! self%rmat(part_coords(1)-wide/2:part_coords(1)-wide/2+length, &
      !    & part_coords(2)-wide/2-int((bborder-1)/2):part_coords(2)-wide/2+int((bborder-1)/2), 1) = value
      ! self%rmat(part_coords(1)-wide/2-int((bborder-1)/2):part_coords(1)-wide/2+int((bborder-1)/2),&
      !    & part_coords(2)-wide/2:part_coords(2)-wide/2+length, 1) = value
      ! self%rmat(part_coords(1)+wide/2-length:part_coords(1)+wide/2,&
      !    & part_coords(2)-wide/2-int((bborder-1)/2):part_coords(2)-wide/2+int((bborder-1)/2), 1) = value
      ! self%rmat(part_coords(1)+wide/2-int((bborder-1)/2):part_coords(1)+wide/2+int((bborder-1)/2),&
      !    & part_coords(2)-wide/2:part_coords(2)-wide/2+length, 1) = value
      ! self%rmat(part_coords(1)-wide/2:part_coords(1)-wide/2+length,&
      !    & part_coords(2)+wide/2-int((bborder-1)/2):part_coords(2)+wide/2+int((bborder-1)/2), 1) = value
      ! self%rmat(part_coords(1)-wide/2-int((bborder-1)/2):part_coords(1)-wide/2+int((bborder-1)/2),&
      !    & part_coords(2)+wide/2-length:part_coords(2)+wide/2, 1) = value
      ! self%rmat(part_coords(1)+wide/2-length:part_coords(1)+wide/2,&
      !    & part_coords(2)+wide/2-int((bborder-1)/2):part_coords(2)+wide/2+int((bborder-1)/2), 1) = value
      ! self%rmat(part_coords(1)+wide/2-int((bborder-1)/2):part_coords(1)+wide/2+int((bborder-1)/2),&
      ! & part_coords(2)+wide/2-length:part_coords(2)+wide/2, 1) = value
      ! Central cross
      self%rmat(part_coords(1)-length:part_coords(1)+length, &
             &  part_coords(2)-int((bborder-1)/2):part_coords(2)+int((bborder-1)/2), 1) = value
      self%rmat(part_coords(1)-int((bborder-1)/2):part_coords(1)+int((bborder-1)/2), &
             &  part_coords(2)-length:part_coords(2)+length, 1) = value
    end subroutine draw_picked

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
    function before_after( left, right, mask ) result( ba )
        class(image),      intent(in) :: left, right
        logical, optional, intent(in) :: mask(left%ldim(1),left%ldim(2),left%ldim(3))
        integer     :: ldim(3), i, j
        type(image) :: ba
        if( left.eqdims.right )then
            if( left.eqsmpd.right )then
                if( left%ft .or. right%ft ) THROW_HARD('not for FTs; before_after')
                if( left%is_3d() .or. right%is_3d() ) THROW_HARD('not for 3D imgs; before_after')
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
                THROW_HARD('before (left) and after (right) not of same smpd; before_after')
            endif
        else
            THROW_HARD('before (left) and after (right) not of same dim; before_after')
        endif
    end function before_after

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
        self%cmat = self%cmat/sqrt(real(product(self%ldim)))
        self%ft = .true.
    end subroutine fwd_ft

    subroutine bwd_ft( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            call fftwf_execute_dft_c2r(self%plan_bwd,self%cmat,self%rmat)
            ! now scale the values accordingly
            self%rmat = self%rmat/sqrt(real(product(self%ldim)))
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
            THROW_HARD('non-equal dims; ft2img')
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
                        THROW_HARD('unsupported mode: '//trim(which)//'; ft2img')
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

    !> \brief dampens the central cross of a powerspectrum by median filtering
    subroutine dampen_central_cross( self )
        class(image), intent(inout) :: self
        integer            :: h,mh,k,mk,lims(3,2),inds(3)
        integer, parameter :: XDAMPWINSZ=2
        real, allocatable  :: pixels(:)
        if( self%ft )          THROW_HARD('not intended for FTs; dampen_central_cross')
        if( self%ldim(3) > 1 ) THROW_HARD('not intended for 3D imgs; dampen_central_cross')
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
        call tmp%kill()
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
        call maskimg%kill
        deallocate(l_mask)
    end subroutine frc_pspec

    !>  \brief  an image shifter to prepare for Fourier transformation
    subroutine shift_phorig( self )
        class(image), intent(inout) :: self
        integer :: i, j, k
        real    :: rswap
        integer :: kfrom,kto
        if( self%ft ) THROW_HARD('this method is intended for real images; shift_phorig')
        if( self%even_dims() )then
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
        else
            write(*,*) 'ldim: ', self%ldim
            THROW_HARD('even dimensions assumed; shift_phorig')
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
    end subroutine shift2Dserial_1

    subroutine shift2Dserial_2( self, shvec, self_out )
        class(image), intent(inout) :: self, self_out
        real,         intent(in)    :: shvec(2)
        integer :: h, k, lims(3,2), phys(3)
        real    :: shvec_here(3)
        lims            = self%fit%loop_lims(2)
        shvec_here(1:2) = shvec
        shvec_here(3)   = 0.
        self_out%ft     = .true.
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                phys = self%fit%comp_addr_phys([h,k,0])
                self_out%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) *&
                    &self%oshift([h,k,0], shvec_here)
            end do
        end do
    end subroutine shift2Dserial_2

    !> \brief mask  is for spherical masking
    !! \param mskrad mask radius
    !! \param which mask type
    !! \param inner include cosine edge material
    !! \param width width of inner patch
    subroutine mask( self, mskrad, which, inner, width, backgr )
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        character(len=*), intent(in)    :: which
        real, optional,   intent(in)    :: inner, width, backgr
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
        soft  = .true.
        select case(trim(which))
            case('soft')
                soft  = .true.
                if(present(backgr))then
                    self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
                        &self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - backgr
                else
                    call self%zero_background
                endif
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
            else
                ! 2d
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
    subroutine pad( self_in, self_out, backgr )
        class(image),   intent(inout) :: self_in, self_out
        real, optional, intent(in)    :: backgr
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
                antialw = self_in%hannw()
                lims = self_in%fit%loop_lims(2)
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
    subroutine pad_mirr( self_in, self_out )
        class(image),   intent(inout) :: self_in, self_out
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
    end subroutine pad_mirr

    !> \brief clip is a constructor that clips the input image to input ldim
    !!             DO NOT PARALLELISE
    !! \param self_in image object
    !! \param self_out image object
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
                            phys_out = self_out%fit%comp_addr_phys(h,k,l)
                            phys_in = self_in%fit%comp_addr_phys(h,k,l)
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
    subroutine clip_inplace( self, ldim )
        class(image), intent(inout) :: self
        integer, intent(in)         :: ldim(3)
        type(image)                 :: tmp
        call tmp%new(ldim, self%smpd)
        call self%clip(tmp)
        call self%copy(tmp)
        call tmp%kill()
    end subroutine clip_inplace

    !!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!!
    ! This subroutine rescales the pixel intensities to a new input range.
    subroutine scale_pixels(self, new_range)
          class(image), intent(inout) :: self
          real,         intent(in)    :: new_range(2)
          real :: old_range(2), sc
          if( .not. self%is_2d() ) THROW_HARD('only for 2D images; scale_pixels')
          old_range(1) = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
          old_range(2) = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
          sc = (new_range(2) - new_range(1))/(old_range(2) - old_range(1))
          self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = sc*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))+new_range(1)-sc*old_range(1)
    end subroutine scale_pixels

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

    subroutine norm4viz( self  )
        class(image), intent(inout) :: self
        if(self%is_ft())THROW_HARD('real space only; norm4viz')
        call self%norm
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 128. + 10.5 *& ! magic numbers from Joe
            &self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
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

    subroutine noise_norm( self, lmsk )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        integer :: npix
        real    :: ave, var, ep
        npix = product(self%ldim) - count(lmsk) ! # background pixels
        ave  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.not. lmsk) / real(npix) ! background average
        if( abs(ave) > TINY ) self%rmat = self%rmat - ave
        ep     = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),      mask=.not. lmsk)
        var    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0, mask=.not. lmsk)
        var    = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        if( is_a_number(var) )then
            if( var > 0. ) self%rmat = self%rmat / sqrt(var)
        endif
    end subroutine noise_norm

    !>  \brief  putting the edge around the image to zero (necessary for avoiding FT artefacts)
    subroutine zero_edgeavg( self )
        class(image), intent(inout) :: self
        real :: edges_sum, edges_ave
        if( self%ft )           THROW_HARD('not for Fted images; edge_norm')
        if( .not.self%is_2d() ) THROW_HARD('only for 2d images; edge_norm')
        edges_sum = sum(self%rmat(1:self%ldim(1),1,1))
        edges_sum = edges_sum + sum(self%rmat(1:self%ldim(1),self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(1,1:self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(self%ldim(1),1:self%ldim(2),1))
        edges_ave = edges_sum / real( 2*(self%ldim(1)+self%ldim(2)) )
        if( abs(edges_ave) > TINY )self%rmat = self%rmat - edges_ave
    end subroutine zero_edgeavg

    !>  \brief  is for [0,1] interval normalization of an image
    subroutine norm_bin( self )
        class(image), intent(inout) :: self
        real                        :: smin, smax
        if( self%ft ) THROW_HARD('image assumed to be real not FTed; norm_bin')
        ! find minmax
        smin  = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        smax  = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        ! create [0,1]-normalized image
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            &(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - smin)  / (smax-smin)
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
            &(exp(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))-1.) / (exp(1.)-1.)
    end subroutine norm_bin

    !> \brief roavg  is for creating a rotation average of self
    !! \param angstep angular step
    !! \param avg output image rotation average
    subroutine roavg( self, angstep, avg, ang_stop )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: angstep
        class(image),      intent(inout) :: avg
        integer, optional, intent(in) :: ang_stop
        real    :: avgs_rmat(nthr_glob,self%ldim(1),self%ldim(2),1)
        real    :: rotated(nthr_glob,self%ldim(1),self%ldim(2),1)
        integer :: irot, ithr, aang_stop
        aang_stop = 359
        if( present(ang_stop) ) aang_stop = ang_stop
        avgs_rmat = 0.
        rotated   = 0.
        !$omp parallel do schedule(static) default(shared) private(irot,ithr) proc_bind(close)
        do irot = 0 + angstep,aang_stop,angstep
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
        if( self_in%ldim(3) > 1 )         THROW_HARD('only for 2D images; rtsq')
        if( .not. self_in%square_dims() ) THROW_HARD('only for square dims; rtsq;')
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
        if( present_outliers )then
            if( allocated(outliers) ) deallocate(outliers)
            allocate( outliers(self%ldim(1),self%ldim(2)), stat=alloc_stat)
            if(alloc_stat/=0)call allocchk("In simple_image::cure_outliers ")
            outliers = .false.
        endif
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

    !>  \brief  zero pixels below thres
    subroutine zero_below( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        where( self%rmat < thres ) self%rmat = 0.
    end subroutine zero_below

    !!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!!
    ! build_ellipse construts an ellipse centered in center
    ! with axis length equal to axis and rotation angle rot.
    ! This subroutine is not a constructor since
    ! it supposes self has already been created.
    subroutine build_ellipse(self, center, axis, rot)
        class(image), intent(inout) :: self
        real,         intent(in)    :: center(2), axis(2), rot
        real, allocatable :: theta(:)
        integer           :: i, j, k
        if(rot < 0. .or. rot > 360. ) THROW_HARD("please insert an angle in the range [0,360]")
        if(self%ldim(3) /= 1) THROW_HARD("the image has to be 2D!")
        theta = (/ (deg2rad(real(i)),i=1,360,2) /)
        do k = 1,size(theta)
            do i = 1, self%ldim(1)
                do j = 1,self%ldim(2)
                    if(abs(real(i) - center(1) - axis(1)*cos(theta(k))*cos(deg2rad(rot))&
                                             & + axis(2)*sin(theta(k))*sin(deg2rad(rot)))<1 .and. &
                    &  abs(real(j) - center(1) - axis(1)*cos(theta(k))*sin(deg2rad(rot))&
                                             & - axis(2)*sin(theta(k))*cos(deg2rad(rot)))<1) then
                        call self%set([i,j,1], 1.)
                        call self%set([i+1,j+1,1], 0.)
                        call self%set([i-1,j-1,1], 0.)
                        call self%set([i+2,j+2,1], 0.)
                        call self%set([i-2,j-2,1], 0.)
                    end if
                enddo
            enddo
        enddo
        deallocate(theta)
    end subroutine build_ellipse

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
            if( .not. img%exists() )   THROW_HARD('ERROR, in constructor or in exists function, 1')
            if( .not. img3d%exists() ) THROW_HARD('ERROR, in constructor or in exists function, 2')

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
            if( .not. passed )  THROW_HARD('getters/setters test failed')

            write(*,'(a)') '**info(simple_image_unit_test, part 4): testing checkups'
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
                write(*,*) ""
                print *, ' checkups ', test
                write(*,*) ""
                THROW_HARD('checkups test failed')
            endif
            write(*,'(a)') '**info(simple_image_unit_test, part 6): testing stats'
            passed = .false.
            call img%gauran( 5., 15. )
            call img%stats( 'foreground', ave, sdev, maxv, minv, 40., med )
            if( ave >= 4. .and. ave <= 6. .and. sdev >= 14. .and.&
                sdev <= 16. .and. med >= 4. .and. med <= 6. ) passed = .true.
            if( .not. passed )  THROW_HARD('stats test failed')

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
            if( .not. passed ) THROW_HARD('corr test failed')

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

            write(*,'(a)') '**info(simple_image_unit_test, part 21): testing destructor'
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
            call fftwf_destroy_plan(self%plan_fwd)
            call fftwf_destroy_plan(self%plan_bwd)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_image

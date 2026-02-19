!@descr: Types and interfaces for production of Cartesian class averages
module simple_new_classaverager
use simple_core_module_api
use simple_builder,           only: build_glob
use simple_ctf,               only: ctf
use simple_discrete_stack_io, only: dstack_io
use simple_euclid_sigma2,     only: euclid_sigma2, eucl_sigma2_glob
use simple_image,             only: image
use simple_parameters,        only: params_glob
use simple_memoize_ft_maps
use simple_fftw3
use simple_ftiter
use simple_imgfile
implicit none

! Module public data & routines are prefixed with cavger_new_
! Init & book-keeping
public :: cavger_new_new, cavger_new_transf_oridat, cavger_new_gen2Dclassdoc
public :: cavger_new_read_euclid_sigma2, cavger_new_kill
! Interpolation & restoration
public :: cavger_new_assemble_sums, cavger_new_restore_cavgs
! I/O & handling of distributed sums
public :: cavger_new_write_eo, cavger_new_write_all, cavger_new_write_merged, cavger_new_read_all
public :: cavger_new_readwrite_partial_sums, cavger_new_assemble_sums_from_parts
! Stacks used for alignment
public :: cavgs_even_new, cavgs_odd_new, cavgs_merged_new
! Separate public utility to rotate particles
! public :: transform_ptcls
private
#include "simple_local_flags.inc"

! Data & private  types

! Convenience type to handle real/reciprocal domain individual slice/cavg pointers
type :: matrix_ptrs
    real(kind=c_float),            pointer :: r(:,:) => null()
    complex(kind=c_float_complex), pointer :: c(:,:) => null()
    logical                                :: ft =.false.
end type matrix_ptrs

! To handle sets of cavgs Fourier images & CTF^2
type :: stack
    type(c_ptr)                            :: p           = c_null_ptr  ! allocation pointer
    type(c_ptr)                            :: plan_fwd    = c_null_ptr  ! forward plan
    type(c_ptr)                            :: plan_bwd    = c_null_ptr  ! backward plan
    type(ftiter)                           :: fit                       ! convenience fourier iterator
    real(kind=c_float),            pointer :: rmat(:,:,:) => null()     ! pointer to real cagvs array
    complex(kind=c_float_complex), pointer :: cmat(:,:,:) => null()     ! pointer to complex cavgs array
    type(matrix_ptrs),         allocatable :: slices(:)                 ! pointers to individual slices/cavgs
    real,                      allocatable :: ctfsq(:,:,:)              ! CTF^2 array
    real,                      allocatable :: soft_mask(:,:)            ! convenience soft mask
    logical,                   allocatable :: nyq_mask(:,:)             ! convenience resolution mask
    integer            :: ldim(2)    = 0                                ! slice dimension
    integer            :: rshape(2)  = 0                                ! slice real space physical dimensions
    integer            :: cshape(2)  = 0                                ! slice complex physical dimensions
    integer            :: flims(3,2) = 0                                ! slice non-redundant Fourier limits
    integer            :: nslices    = 0                                ! # of slices/cavgs
  contains
    procedure          :: new_stack
    procedure          :: zero, zero_slice
    procedure          :: write, read_cmat
    procedure          :: write_cmat_debug, write_ctfsq_debug
    procedure          :: write_ctfsq, read_ctfsq
    procedure, private :: calc_cavgs_stats
    procedure, private :: quadrant_swap
    procedure          :: fft, ifft
    procedure          :: frc
    procedure          :: ctf_dens_correct
    procedure          :: softmask
    procedure          :: insert_lowres_serial
    procedure          :: add_invnoisepower2rho
    procedure          :: kill_stack
end type stack

! To handle even, odd & merged cavgs
type cavgs_set
    type(stack)        :: even               !< even cavgs
    type(stack)        :: odd                !< odd cavgs
    type(stack)        :: merged             !< merged cavgs
    integer            :: ldim(2)    = 0     !< real dimensions
    integer            :: ncls       = 0     !< # of cavgs
  contains
    procedure          :: new_set
    procedure          :: zero_set
    procedure          :: copy_fast
    procedure          :: kill_set
end type cavgs_set

! To hande particle metadata
type ptcl_record
    type(ctf)          :: tfun               !< transfer function
    real               :: pw         = 0.0   !< particle weight
    real               :: dfx        = 0.0   !< defocus in x (microns)
    real               :: dfy        = 0.0   !< defocus in y (microns)
    real               :: angast     = 0.0   !< angle of astigmatism (in degrees)
    real               :: e3         = 0.0   !< in-plane rotations
    real               :: shift(2)   = 0.0   !< rotational origin shift
    integer            :: pind       = 0     !< particle index
    integer            :: eo         = -1    !< even is 0, odd is 1, default is -1
    integer            :: class      = 0     !< class assignment
    integer            :: ind_in_stk = 0     !< index in stack
end type ptcl_record

! Main module variables
type(ptcl_record),   allocatable :: precs(:)                  !< Particle records
type(image), target, allocatable :: cavgs_even_new(:)         !< Even class averages for reading
type(image), target, allocatable :: cavgs_odd_new(:)          !< Odd class averages for reading
type(image), target, allocatable :: cavgs_merged_new(:)       !< Merged class averages for reading
type(cavgs_set)                  :: cavgs                     !< Class averages
type(euclid_sigma2)              :: eucl_sigma                !< Noise power estimates
logical,             allocatable :: pptcl_mask(:)             !< selected particles
integer                          :: ctfflag                   !< ctf flag <yes=1|no=0|flip=2>
integer                          :: istart      = 0, iend = 0 !< particle index range in partition
integer                          :: partsz      = 0           !< size of partition
integer                          :: ncls        = 0           !< # classes
integer                          :: ldim(3)        = [0,0,0]  !< logical dimension of image
integer                          :: ldim_crop(3)   = [0,0,0]  !< logical dimension of cropped image
integer                          :: ldim_pd(3)     = [0,0,0]  !< logical dimension of image, padded
integer                          :: ldim_croppd(3) = [0,0,0]  !< logical dimension of cropped image, padded
real                             :: smpd       = 0.           !< sampling distance
real                             :: smpd_crop  = 0.           !< cropped sampling distance
logical                          :: l_alloc_read_cavgs=.true. !< whether to allocate sums and read partial sums

interface

    !
    ! Type cavgs_set: for managing a set even, odd and merged cavgs 
    !

    module subroutine new_set( self, ldim, ncls )
        class(cavgs_set), intent(inout) :: self
        integer,         intent(in)    :: ldim(2), ncls
    end subroutine new_set

    module subroutine zero_set( self, ft )
        class(cavgs_set), intent(inout) :: self
        logical,          intent(in)    :: ft
    end subroutine zero_set

    module subroutine copy_fast( self, self2copy, is, ft )
        class(cavgs_set), intent(inout) :: self
        class(cavgs_set), intent(in)    :: self2copy
        integer,         intent(in)    :: is
        logical,         intent(in)    :: ft
    end subroutine

    module subroutine kill_set( self )
        class(cavgs_set), intent(inout) :: self
    end subroutine kill_set

    !
    ! Type stack: for managing complex and CTF2 arrays
    !

    module subroutine new_stack( self, ldim, nslices, alloc_ctfsq )
        class(stack),      intent(inout) :: self
        integer,           intent(in)    :: ldim(2), nslices
        logical, optional, intent(in)    :: alloc_ctfsq
    end subroutine new_stack

    module subroutine zero( self, ft )
        class(stack), intent(inout) :: self
        logical,      intent(in)    :: ft
    end subroutine zero

    module subroutine zero_slice( self, is, ft )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
        logical,      intent(in)    :: ft
    end subroutine zero_slice

    module subroutine write( self, fname, ft )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        logical,       intent(in)    :: ft
    end subroutine

    module subroutine write_cmat_debug( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
    end subroutine write_cmat_debug

    module subroutine write_ctfsq_debug( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
    end subroutine write_ctfsq_debug

    module subroutine write_ctfsq( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
    end subroutine write_ctfsq

    module subroutine read_cmat( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
    end subroutine read_cmat

    module subroutine read_ctfsq( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
    end subroutine read_ctfsq

    module subroutine calc_cavgs_stats( self, stats )
        class(stack), intent(in)  :: self
        real,         intent(out) :: stats(4)
    end subroutine calc_cavgs_stats

    module subroutine quadrant_swap( self, is )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
    end subroutine quadrant_swap

    module subroutine fft( self, i )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: i
    end subroutine fft

    module subroutine ifft( self, i )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: i
    end subroutine ifft

    module subroutine frc( self1, self2, is, corrs )
        class(stack), intent(in)  :: self1, self2
        integer,      intent(in)  :: is
        real,         intent(out) :: corrs(fdim(self1%ldim(1))-1)
    end subroutine frc

    module subroutine ctf_dens_correct( self, is )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
    end subroutine ctf_dens_correct

    module subroutine softmask( self, is )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
    end subroutine softmask

    module subroutine insert_lowres_serial( self, self2insert, is, find )
        class(stack), intent(inout) :: self
        class(stack), intent(in)    :: self2insert
        integer,      intent(in)    :: is, find
    end subroutine insert_lowres_serial

    module subroutine add_invnoisepower2rho( self, is, frcsz, frc )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is, frcsz
        real,         intent(in)    :: frc(1:frcsz)
    end subroutine add_invnoisepower2rho

    module subroutine kill_stack( self )
        class(stack), intent(inout) :: self
    end subroutine kill_stack

    !
    ! Module public routines
    !

    module subroutine cavger_new_new( pinds, alloccavgs )
        integer, optional, intent(in) :: pinds(:)
        logical, optional, intent(in) :: alloccavgs 
    end subroutine cavger_new_new

    ! Book-keeping & metadata

    module subroutine cavger_new_transf_oridat( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
    end subroutine cavger_new_transf_oridat

    module subroutine cavger_new_read_euclid_sigma2
    end subroutine cavger_new_read_euclid_sigma2

    module subroutine cavger_new_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), target, intent(inout) :: spproj
    end subroutine cavger_new_gen2Dclassdoc

    ! Restoration

    module subroutine cavger_new_assemble_sums( do_frac_update )
        logical, intent(in)      :: do_frac_update
    end subroutine cavger_new_assemble_sums

    module subroutine cavger_new_assemble_sums_conv( do_frac_update )
        logical, intent(in)      :: do_frac_update
    end subroutine cavger_new_assemble_sums_conv

    module subroutine cavger_new_restore_cavgs( frcs_fname )
        use simple_gridding, only: prep2D_inv_instrfun4mul
        class(string), intent(in) :: frcs_fname
    end subroutine cavger_new_restore_cavgs

    ! I/O & assembly

    module subroutine cavger_new_write_eo( fname_e, fname_o )
        class(string), intent(in) :: fname_e, fname_o
    end subroutine cavger_new_write_eo

    module subroutine cavger_new_write_all( fname, fname_e, fname_o )
        class(string), intent(in) :: fname, fname_e, fname_o
    end subroutine cavger_new_write_all

    module subroutine cavger_new_write_merged( fname )
        class(string), intent(in) :: fname
    end subroutine cavger_new_write_merged

    module subroutine cavger_new_read_all()
    end subroutine cavger_new_read_all

    module subroutine cavger_new_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
    end subroutine cavger_new_readwrite_partial_sums

    module subroutine cavger_new_assemble_sums_from_parts
    end subroutine cavger_new_assemble_sums_from_parts

    ! Destructors

    module subroutine cavger_new_kill( dealloccavgs )
        logical, optional, intent(in) :: dealloccavgs
    end subroutine cavger_new_kill

    ! Public utility

    module subroutine transform_ptcls( spproj, oritype, icls, timgs, pinds, phflip, cavg, imgs_ori)
        use simple_sp_project,          only: sp_project
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimgbatch
        use simple_memoize_ft_maps
        class(sp_project),                  intent(inout) :: spproj
        character(len=*),                   intent(in)    :: oritype
        integer,                            intent(in)    :: icls
        type(image),           allocatable, intent(inout) :: timgs(:)
        integer,               allocatable, intent(inout) :: pinds(:)
        logical,     optional,              intent(in)    :: phflip
        type(image), optional,              intent(inout) :: cavg
        type(image), optional, allocatable, intent(inout) :: imgs_ori(:)
    end subroutine transform_ptcls

end interface

end module simple_new_classaverager

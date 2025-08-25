! for calculation of band-pass limited cross-correlation of polar Fourier transforms
module simple_polarft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
implicit none

public :: polarft_corrcalc, pftcc_glob
private
#include "simple_local_flags.inc"

type fftw_cvec
    type(c_ptr)                            :: p
    complex(kind=c_float_complex), pointer :: c(:) => null()
end type fftw_cvec

type fftw_rvec
    type(c_ptr)                 :: p
    real(kind=c_float), pointer :: r(:) => null()
end type fftw_rvec

type fftw_drvec
    type(c_ptr)                  :: p
    real(kind=c_double), pointer :: r(:) => null()
end type fftw_drvec

type heap_vars
    complex(sp), pointer :: pft_tmp(:,:)       => null()
    complex(sp), pointer :: pft_tmp2(:,:)      => null()
    complex(sp), pointer :: pft_ref(:,:)       => null()
    complex(sp), pointer :: pft_ref_tmp(:,:)   => null()
    real(dp),    pointer :: argvec(:)          => null()
    complex(sp), pointer :: shmat(:,:)         => null()
    real(dp),    pointer :: kcorrs(:)          => null()
    complex(dp), pointer :: pft_ref_8(:,:)     => null()
    complex(dp), pointer :: pft_ref_tmp_8(:,:) => null()
    complex(dp), pointer :: pft_ptcl_8(:,:)    => null()
    complex(dp), pointer :: pft_dref_8(:,:,:)  => null()
    complex(dp), pointer :: shvec(:)           => null()
    complex(dp), pointer :: shmat_8(:,:)       => null()
    real(dp),    pointer :: pft_r1_8(:,:)      => null()
    real(dp),    pointer :: pft_r2_8(:,:)      => null()
    real(sp),    pointer :: pft_r(:,:)         => null()
end type heap_vars

type :: polarft_corrcalc
    private
    integer                          :: nptcls     = 1              !< the total number of particles in partition (logically indexded [fromp,top])
    integer                          :: nrefs      = 1              !< the number of references (logically indexded [1,nrefs])
    integer                          :: nrots      = 0              !< number of in-plane rotations for one pft (determined by radius of molecule)
    integer                          :: pftsz      = 0              !< size of reference and particle pft (nrots/2)
    integer                          :: pfromto(2) = 0              !< particle index range
    integer                          :: ldim(3)    = 0              !< logical dimensions of original cartesian image
    integer                          :: kfromto(2)                  !< band-pass Fourier index limits
    integer                          :: nk                          !< number of shells used during alignment
    real                             :: dang                        !< angular increment
    integer,             allocatable :: pinds(:)                    !< index array (to reduce memory when frac_update < 1)
    logical,             allocatable :: has_cache(:,:)              !< whether the corr is cached
    logical,             allocatable :: do_cache(:)                 !< should do caching or not
    real,                allocatable :: cached_vals(:,:)            !< cached vals
    real,                allocatable :: npix_per_shell(:)           !< number of (cartesian) pixels per shell
    real(dp),            allocatable :: sqsums_ptcls(:)             !< memoized square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),            allocatable :: ksqsums_ptcls(:)            !< memoized k-weighted square sums for the correlation calculations (taken from kfromto(1):kfromto(2))
    real(dp),            allocatable :: wsqsums_ptcls(:)            !< memoized square sums weighted by k and  sigmas^2 (taken from kfromto(1):kfromto(2))
    real(sp),            allocatable :: angtab(:)                   !< table of in-plane angles (in degrees)
    real(dp),            allocatable :: argtransf(:,:)              !< argument transfer constants for shifting the references
    real(sp),            allocatable :: polar(:,:)                  !< table of polar coordinates (in Cartesian coordinates)
    real(sp),            allocatable :: ctfmats(:,:,:)              !< expand set of CTF matrices (for efficient parallel exec)
    real(dp),            allocatable :: argtransf_shellone(:)       !< one dimensional argument transfer constants (shell k=1) for shifting the references
    complex(sp),         allocatable :: pfts_refs_even(:,:,:)       !< 3D complex matrix of polar reference sections (pftsz,nk,nrefs), even
    complex(sp),         allocatable :: pfts_refs_odd(:,:,:)        !< -"-, odd
    complex(sp),         allocatable :: pfts_drefs_even(:,:,:,:)    !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp),         allocatable :: pfts_drefs_odd(:,:,:,:)     !< derivatives w.r.t. orientation angles of 3D complex matrices
    complex(sp),         allocatable :: pfts_ptcls(:,:,:)           !< 3D complex matrix of particle sections
    ! FFTW plans
    type(c_ptr)                      :: plan_fwd1, plan_bwd1
    type(c_ptr)                      :: plan_mem_r2c
    ! Memoized terms
    complex(kind=c_float_complex), allocatable :: ft_ptcl_ctf(:,:,:)    !< Fourier Transform of particle times CTF
    complex(kind=c_float_complex), allocatable :: ft_absptcl_ctf(:,:,:) !< Fourier Transform of (particle times CTF)**2
    complex(kind=c_float_complex), allocatable :: ft_ctf2(:,:,:)        !< Fourier Transform of CTF squared modulus
    complex(kind=c_float_complex), allocatable :: ft_ref_even(:,:,:),     ft_ref_odd(:,:,:)     !< Fourier Transform of even/odd references
    complex(kind=c_float_complex), allocatable :: ft_ref2_even(:,:,:),    ft_ref2_odd(:,:,:)    !< Fourier Transform of even/odd references squared modulus
    complex(kind=c_float_complex), allocatable :: ft_refs2_even(:,:,:,:), ft_refs2_odd(:,:,:,:) !< Fourier Transform of even/odd references of different pairs
    ! Convenience vectors, thread memoization
    type(heap_vars),     allocatable :: heap_vars(:)
    type(fftw_cvec),     allocatable :: cvec1(:), cvec2(:)
    type(fftw_rvec),     allocatable :: rvec1(:)
    type(fftw_drvec),    allocatable :: drvec(:)
    ! Others
    logical,             allocatable :: iseven(:)                   !< eo assignment for gold-standard FSC
    real,                pointer     :: sigma2_noise(:,:) => null() !< for euclidean distances
    logical                          :: with_ctf  = .false.         !< CTF flag
    logical                          :: existence = .false.         !< to indicate existence

    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: reallocate_ptcls
    procedure          :: set_ref_pft
    procedure          :: set_ptcl_pft
    procedure          :: set_ref_fcomp
    procedure          :: set_dref_fcomp
    procedure          :: set_ptcl_fcomp
    procedure          :: cp_even2odd_ref
    procedure          :: cp_odd2even_ref
    procedure          :: cp_even_ref2ptcl
    procedure          :: cp_refs
    procedure          :: swap_ptclsevenodd
    procedure          :: set_eo
    procedure          :: set_eos
    procedure          :: set_with_ctf
    procedure          :: assign_sigma2_noise
    ! GETTERS
    procedure          :: get_nrots
    procedure          :: get_pdim
    procedure          :: get_kfromto
    procedure          :: get_pftsz
    procedure          :: get_rot
    procedure          :: get_roind, get_roind_fast
    procedure          :: get_dang
    procedure          :: get_coord
    procedure, private :: get_polar_coord_1, get_polar_coord_2
    generic            :: get_polar_coord => get_polar_coord_1, get_polar_coord_2
    procedure          :: get_ref_pft
    procedure          :: get_nrefs
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: get_pinds
    procedure          :: get_npix
    procedure          :: get_work_pft_ptr, get_work_rpft_ptr
    procedure          :: get_ptcls_ptr
    procedure          :: get_ctfmats_ptr
    procedure          :: get_refs_ptr
    procedure          :: get_linstates_irefs
    procedure          :: get_linstates_prefs
    procedure          :: set_cache
    procedure          :: reset_cache
    procedure          :: is_with_ctf
    ! PRINTERS/VISUALISERS
    procedure          :: print
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    procedure          :: prefs_to_cartesian
    procedure, private :: polar2cartesian_1, polar2cartesian_2
    generic            :: polar2cartesian => polar2cartesian_1, polar2cartesian_2
    procedure          :: read_ptcls_ctfs, write_ptcls_ctfs
    ! MODIFIERS
    procedure          :: shift_ptcl
    procedure          :: shift_ref
    procedure          :: mirror_ref_pft
    ! MEMOIZER
    procedure          :: memoize_sqsum_ptcl
    procedure, private :: setup_npix_per_shell
    procedure          :: memoize_ptcls, memoize_refs
    procedure, private :: kill_memoized_ptcls, kill_memoized_refs
    procedure          :: allocate_ptcls_memoization, allocate_refs_memoization
    ! CALCULATORS
    procedure          :: create_polar_absctfmats, calc_polar_ctf
    procedure          :: gen_shmat
    procedure, private :: gen_shmat_8
    procedure          :: calc_shift
    procedure          :: calc_corr_rot_shift, calc_magcorr_rot
    procedure          :: bidirectional_shift_search
    procedure          :: gencorrs_mag, gencorrs_mag_cc
    procedure          :: gencorrs_weighted_cc,      gencorrs_shifted_weighted_cc
    procedure          :: gencorrs_cc,               gencorrs_shifted_cc
    procedure, private :: gencorrs_euclid_1,         gencorrs_euclid_2
    procedure, private :: gencorrs_shifted_euclid_1, gencorrs_shifted_euclid_2
    generic            :: gencorrs_euclid         => gencorrs_euclid_1, gencorrs_euclid_2
    generic            :: gencorrs_shifted_euclid => gencorrs_shifted_euclid_1, gencorrs_shifted_euclid_2
    procedure, private :: gencorrs_1, gencorrs_2, gencorrs_3, gencorrs_4
    generic            :: gencorrs => gencorrs_1, gencorrs_2, gencorrs_3, gencorrs_4
    procedure, private :: gencorr_for_rot_8_1, gencorr_for_rot_8_2, gencorr_for_rot_8_3, gencorr_for_rot_8_4
    generic            :: gencorr_for_rot_8 => gencorr_for_rot_8_1, gencorr_for_rot_8_2, gencorr_for_rot_8_3, gencorr_for_rot_8_4
    procedure, private :: gencorr_grad_for_rot_8_1,      gencorr_grad_for_rot_8_2
    procedure, private :: gencorr_grad_only_for_rot_8_1, gencorr_grad_only_for_rot_8_2
    generic            :: gencorr_grad_for_rot_8      => gencorr_grad_for_rot_8_1,      gencorr_grad_for_rot_8_2
    generic            :: gencorr_grad_only_for_rot_8 => gencorr_grad_only_for_rot_8_1, gencorr_grad_only_for_rot_8_2
    procedure          :: gencorr_cc_for_rot_8
    procedure          :: gencorr_cc_grad_for_rot_8
    procedure          :: gencorr_cc_grad_only_for_rot_8
    procedure          :: gencorr_euclid_for_rot_8
    procedure          :: gencorr_euclid_grad_for_rot_8
    procedure          :: gencorr_cont_grad_cc_for_rot_8
    procedure          :: gencorr_cont_cc_for_rot_8
    procedure          :: gencorr_cont_shift_grad_cc_for_rot_8
    procedure          :: gencorr_cont_grad_euclid_for_rot_8
    procedure          :: gencorr_cont_shift_grad_euclid_for_rot_8
    procedure          :: gencorr_sigma_contrib
    procedure, private :: calc_frc, linear_ft_refs
    generic            :: linear_ft_ref2s => linear_ft_ref2s_1, linear_ft_ref2s_2
    procedure, private :: linear_ft_ref2s_1, linear_ft_ref2s_2
    procedure, private :: rotate_pft_1, rotate_pft_2, rotate_pft_3, rotate_pft_4
    generic            :: rotate_pft => rotate_pft_1, rotate_pft_2, rotate_pft_3, rotate_pft_4
    procedure, private :: rotate_iref_1, rotate_iref_2
    generic            :: rotate_iref => rotate_iref_1, rotate_iref_2
    ! DESTRUCTOR
    procedure          :: kill
end type polarft_corrcalc

! CLASS PARAMETERS/VARIABLES
complex(sp), parameter           :: zero            = cmplx(0.,0.) !< just a complex zero
integer,     parameter           :: FFTW_USE_WISDOM = 16
class(polarft_corrcalc), pointer :: pftcc_glob => null()

contains

    ! CONSTRUCTORS

    subroutine new( self, nrefs, pfromto, kfromto, eoarr )
        class(polarft_corrcalc), target, intent(inout) :: self
        integer,                         intent(in)    :: nrefs
        integer,                         intent(in)    :: pfromto(2), kfromto(2)
        integer, optional,               intent(in)    :: eoarr(pfromto(1):pfromto(2))
        real(sp), allocatable :: polar_here(:)
        real(dp)              :: A(2)
        integer               :: irot, k, ithr, i
        logical               :: even_dims, test(2)
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
        ! set band-pass Fourier index limits
        self%kfromto = kfromto
        self%nk      = self%kfromto(2) - self%kfromto(1) + 1
        ! error check
        if( self%pfromto(2) - self%pfromto(1) + 1 < 1 )then
            write(logfhandle,*) 'pfromto: ', self%pfromto(1), self%pfromto(2)
            THROW_HARD ('nptcls (# of particles) must be > 0; new')
        endif
        if( nrefs < 1 )then
            write(logfhandle,*) 'nrefs: ', nrefs
            THROW_HARD ('nrefs (# of reference sections) must be > 0; new')
        endif
        self%ldim = [params_glob%box,params_glob%box,1] !< logical dimensions of original cartesian image
        test      = .false.
        test(1)   = is_even(self%ldim(1))
        test(2)   = is_even(self%ldim(2))
        even_dims = all(test)
        if( .not. even_dims )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            THROW_HARD ('only even logical dims supported; new')
        endif
        ! set constants
        self%nptcls = self%pfromto(2) - self%pfromto(1) + 1   !< the total number of particles in partition
        self%nrefs  = nrefs                                   !< the number of references (logically indexded [1,nrefs])
        self%pftsz  = magic_pftsz(nint(params_glob%msk_crop)) !< size of reference (number of vectors used for matching,determined by radius of molecule)
        self%nrots  = 2 * self%pftsz                          !< number of in-plane rotations for one pft  (pftsz*2)
        ! generate polar coordinates
        allocate( self%polar(2*self%nrots,self%kfromto(1):self%kfromto(2)),&
                    &self%angtab(self%nrots), self%iseven(1:self%nptcls), polar_here(2*self%nrots))
        self%dang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*self%dang
            ! cycling over non-redundant logical dimensions
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            =  sin(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = -cos(self%angtab(irot))*real(k) ! y-coordinate
            end do
            ! for k = 1
            polar_here(irot)            =  sin(real(self%angtab(irot)))
            polar_here(irot+self%nrots) = -cos(real(self%angtab(irot)))
            ! angle (in degrees) from now
            self%angtab(irot) = rad2deg(self%angtab(irot))
        end do
        self%dang = rad2deg(self%dang)
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0 )
        self%pinds = (/(i,i=1,self%nptcls)/)
        ! eo assignment
        if( present(eoarr) )then
            if( all(eoarr == - 1) )then
                self%iseven = .true.
            else
                do i=self%pfromto(1),self%pfromto(2)
                    if( self%pinds(i) > 0 )then
                        if( eoarr(i) == 0 )then
                            self%iseven(self%pinds(i)) = .true.
                        else
                            self%iseven(self%pinds(i)) = .false.
                        endif
                    endif
                end do
            endif
        else
            self%iseven = .true.
        endif
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%kfromto(1):self%kfromto(2)), self%argtransf_shellone(self%nrots) )
        A = DPI / real(self%ldim(1:2)/2,dp) ! argument transfer matrix normalization constant
        ! shell = 1
        self%argtransf_shellone(:self%pftsz  ) = real(polar_here(:self%pftsz),dp)                        * A(1) ! x-part
        self%argtransf_shellone(self%pftsz+1:) = real(polar_here(self%nrots+1:self%nrots+self%pftsz),dp) * A(2) ! y-part
        ! all shells in resolution range
        self%argtransf(:self%pftsz,:)     = real(self%polar(:self%pftsz,:),dp)                          * A(1)  ! x-part
        self%argtransf(self%pftsz + 1:,:) = real(self%polar(self%nrots + 1:self%nrots+self%pftsz,:),dp) * A(2)  ! y-part
        ! allocate others
        allocate(self%pfts_refs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%pfts_refs_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%pfts_drefs_even(self%pftsz,self%kfromto(1):self%kfromto(2),3,params_glob%nthr),&
                &self%pfts_drefs_odd (self%pftsz,self%kfromto(1):self%kfromto(2),3,params_glob%nthr),&
                &self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                &self%sqsums_ptcls(1:self%nptcls),self%ksqsums_ptcls(1:self%nptcls),self%wsqsums_ptcls(1:self%nptcls),&
                &self%heap_vars(params_glob%nthr))
        do ithr=1,params_glob%nthr
            allocate(self%heap_vars(ithr)%pft_ref(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_ref_tmp(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_tmp(self%pftsz+1,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_tmp2(self%pftsz+1,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%argvec(self%pftsz),&
                    &self%heap_vars(ithr)%shvec(self%pftsz),&
                    &self%heap_vars(ithr)%shmat(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%kcorrs(self%nrots),&
                    &self%heap_vars(ithr)%pft_ref_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_ref_tmp_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_ptcl_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_dref_8(self%pftsz,self%kfromto(1):self%kfromto(2),3),&
                    &self%heap_vars(ithr)%shmat_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_r1_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_r2_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_r(self%pftsz,self%kfromto(1):self%kfromto(2)))
        end do
        self%pfts_refs_even = zero
        self%pfts_refs_odd  = zero
        self%pfts_ptcls     = zero
        self%sqsums_ptcls   = 0.d0
        self%ksqsums_ptcls  = 0.d0
        self%wsqsums_ptcls  = 0.d0
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! setup npix_per_shell
        call self%setup_npix_per_shell
        ! allocation for memoization
        call self%allocate_ptcls_memoization
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        pftcc_glob => self
    end subroutine new

    ! SETTERS

    subroutine reallocate_ptcls( self, nptcls, pinds )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: nptcls
        integer,                 intent(in)    :: pinds(nptcls)
        integer :: i,iptcl
        self%pfromto(1) = minval(pinds)
        self%pfromto(2) = maxval(pinds)
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( self%nptcls == nptcls )then
            ! just need to update particles indexing
        else
            ! re-index & reallocate
            self%nptcls = nptcls
            if( allocated(self%sqsums_ptcls) ) deallocate(self%sqsums_ptcls)
            if( allocated(self%ksqsums_ptcls)) deallocate(self%ksqsums_ptcls)
            if( allocated(self%wsqsums_ptcls)) deallocate(self%wsqsums_ptcls)
            if( allocated(self%iseven) )       deallocate(self%iseven)
            if( allocated(self%pfts_ptcls) )   deallocate(self%pfts_ptcls)
            allocate( self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                     &self%sqsums_ptcls(1:self%nptcls),self%ksqsums_ptcls(1:self%nptcls),self%wsqsums_ptcls(1:self%nptcls),self%iseven(1:self%nptcls))
            call self%kill_memoized_ptcls
            call self%allocate_ptcls_memoization
        endif
        self%pfts_ptcls    = zero
        self%sqsums_ptcls  = 0.d0
        self%ksqsums_ptcls = 0.d0
        self%wsqsums_ptcls = 0.d0
        self%iseven        = .true.
        allocate(self%pinds(self%pfromto(1):self%pfromto(2)), source=0)
        do i = 1,self%nptcls
            iptcl = pinds(i)
            self%pinds( iptcl ) = i
        enddo
        if( .not. self%with_ctf )return
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls), source=1.)
    end subroutine reallocate_ptcls

    subroutine set_ref_pft( self, iref, pft, iseven )
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iref     !< reference index
        complex(sp),             intent(in)    :: pft(:,:) !< reference pft
        logical,                 intent(in)    :: iseven   !< logical eo-flag
        if( iseven )then
            self%pfts_refs_even(:,:,iref) = pft
        else
            self%pfts_refs_odd(:,:,iref)  = pft
        endif
    end subroutine set_ref_pft

    subroutine set_ptcl_pft( self, iptcl, pft )
        class(polarft_corrcalc), intent(inout) :: self     !< this object
        integer,                 intent(in)    :: iptcl    !< particle index
        complex(sp),             intent(in)    :: pft(:,:) !< particle's pft
        self%pfts_ptcls(:,:,self%pinds(iptcl)) = pft
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    subroutine set_ref_fcomp( self, iref, irot, k, comp, iseven )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: comp
        logical,                 intent(in)    :: iseven
        if( iseven )then
            self%pfts_refs_even(irot,k,iref) = comp
        else
            self%pfts_refs_odd(irot,k,iref)  = comp
        endif
    end subroutine set_ref_fcomp

    subroutine set_dref_fcomp( self, iref, irot, k, dcomp, iseven )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, irot, k
        complex(sp),             intent(in)    :: dcomp(3)
        logical,                 intent(in)    :: iseven
        if( iseven )then
            self%pfts_drefs_even(irot,k,:,iref) = dcomp
        else
            self%pfts_drefs_odd(irot,k,:,iref)  = dcomp
        endif
    end subroutine set_dref_fcomp

    subroutine set_ptcl_fcomp( self, iptcl, irot, k, comp )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, irot, k
        complex(sp),             intent(in)    :: comp
        self%pfts_ptcls(irot,k,self%pinds(iptcl)) = comp
    end subroutine set_ptcl_fcomp

    subroutine cp_even2odd_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_odd(:,:,iref) = self%pfts_refs_even(:,:,iref)
    end subroutine cp_even2odd_ref

    subroutine cp_odd2even_ref( self, iref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        self%pfts_refs_even(:,:,iref) = self%pfts_refs_odd(:,:,iref)
    end subroutine cp_odd2even_ref

    subroutine cp_refs( self, self2 )
        class(polarft_corrcalc), intent(inout) :: self, self2
        self%pfts_refs_odd  = self2%pfts_refs_odd
        self%pfts_refs_even = self2%pfts_refs_even
    end subroutine cp_refs

    subroutine cp_even_ref2ptcl( self, iref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        self%pfts_ptcls(:,:,self%pinds(iptcl)) = self%pfts_refs_even(:,:,iref)
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
    end subroutine cp_even_ref2ptcl

    subroutine swap_ptclsevenodd( self )
        class(polarft_corrcalc), intent(inout) :: self
        self%iseven = .not.self%iseven
    end subroutine swap_ptclsevenodd

    subroutine set_eo( self, iptcl, is_even )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        logical,                 intent(in)    :: is_even
        self%iseven(self%pinds(iptcl)) = is_even
    end subroutine set_eo

    subroutine set_eos( self, eoarr )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: eoarr(self%nptcls)
        integer :: i
        if( all(eoarr == - 1) )then
            self%iseven = .true.
        else
            do i=1,self%nptcls
                if( eoarr(i) == 0 )then
                    self%iseven(i) = .true.
                else
                    self%iseven(i) = .false.
                endif
            end do
        endif
    end subroutine set_eos

    subroutine set_with_ctf( self, l_wctf )
        class(polarft_corrcalc), intent(inout) :: self
        logical,                 intent(in)    :: l_wctf
        self%with_ctf = l_wctf
    end subroutine set_with_ctf

    subroutine assign_sigma2_noise( self, sigma2_noise )
        class(polarft_corrcalc),      intent(inout) :: self
        real,    allocatable, target, intent(inout) :: sigma2_noise(:,:)
        self%sigma2_noise => sigma2_noise
    end subroutine assign_sigma2_noise

    ! GETTERS

    !>  \brief  for getting the number of in-plane rotations
    pure function get_nrots( self ) result( nrots )
        class(polarft_corrcalc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots

    !>  \brief  for getting the dimensions of the reference polar FT
    pure function get_pdim( self ) result( pdim )
        class(polarft_corrcalc), intent(in) :: self
        integer :: pdim(3)
        pdim = [self%pftsz,self%kfromto(1),self%kfromto(2)]
    end function get_pdim

    pure function get_kfromto( self ) result( kfromto )
        class(polarft_corrcalc), intent(in) :: self
        integer :: kfromto(2)
        kfromto = [self%kfromto(1),self%kfromto(2)]
    end function get_kfromto

    ! !>  \brief  for getting the dimension of the reference polar FT
    pure integer function get_pftsz( self )
        class(polarft_corrcalc), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    !>  \brief is for getting the continuous in-plane rotation
    !!         corresponding to in-plane rotation index roind
    function get_rot( self, roind ) result( rot )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: roind !< in-plane rotation index
        real(sp) :: rot
        if( roind < 1 .or. roind > self%nrots )then
            write(logfhandle,*) 'roind: ', roind
            write(logfhandle,*) 'nrots: ', self%nrots
            THROW_HARD('roind is out of range; get_rot')
        endif
        rot = self%angtab(roind)
    end function get_rot

    !>  \brief is for getting the discrete in-plane rotational
    !!         index corresponding to continuous rotation rot
    function get_roind( self, rot ) result( ind )
        class(polarft_corrcalc), intent(in) :: self
        real(sp),                intent(in) :: rot !<  continuous rotation
        real(sp) :: dists(self%nrots)
        integer  :: ind
        dists = abs(self%angtab-rot)
        where(dists>180.)dists = 360.-dists
        ind = minloc(dists, dim=1)
    end function get_roind

    pure integer function get_roind_fast( self, psi )
        class(polarft_corrcalc), intent(in) :: self
        real,                    intent(in) :: psi
        get_roind_fast = nint(psi / self%dang) + 1
        if( get_roind_fast <=         0 ) get_roind_fast = get_roind_fast + self%nrots
        if( get_roind_fast > self%nrots ) get_roind_fast = get_roind_fast - self%nrots
    end function get_roind_fast

    pure real function get_dang( self )
        class(polarft_corrcalc), intent(in) :: self
        get_dang = self%dang
    end function get_dang

    !>  \brief returns polar coordinate for rotation rot
    !!         and Fourier index k
    pure function get_coord( self, rot, k ) result( xy )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: rot, k
        real(sp) :: xy(2)
        xy(1) = self%polar(rot,k)
        xy(2) = self%polar(self%nrots+rot,k)
    end function get_coord

    subroutine get_polar_coord_1( self, xy, rot, k )
        class(polarft_corrcalc), intent(in)  :: self
        real,                    intent(in)  :: xy(2)
        real,                    intent(out) :: rot, k
        real :: tan_inv
        k       = sqrt(sum(xy**2))
        tan_inv = atan(xy(2), xy(1)) * 2 + PI
        rot     = 1. + tan_inv * real(self%pftsz) / TWOPI
    end subroutine get_polar_coord_1

    subroutine get_polar_coord_2( self, xy, rot, k )
        class(polarft_corrcalc), intent(in)  :: self
        real,                    intent(in)  :: xy(2)
        integer,                 intent(out) :: rot, k
        real :: r_real, k_real
        call self%get_polar_coord_1(xy, r_real, k_real)
        k   = nint(k_real)
        rot = nint(r_real)
    end subroutine get_polar_coord_2

    !>  \brief  returns polar Fourier transform of reference iref
    subroutine get_ref_pft( self, iref, iseven, pft )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: iref
        logical,                 intent(in)    :: iseven
        complex(sp),             intent(inout) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        if( iseven )then
            pft = self%pfts_refs_even(:,:,iref)
        else
            pft = self%pfts_refs_odd(:,:,iref)
        endif
    end subroutine get_ref_pft

    pure integer function get_nrefs( self )
        class(polarft_corrcalc), intent(in) :: self
        get_nrefs = self%nrefs
    end function get_nrefs

    logical function exists( self )
        class(polarft_corrcalc), intent(in) :: self
        exists = self%existence
    end function exists

    logical function ptcl_iseven( self, iptcl )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        ptcl_iseven = self%iseven(self%pinds(iptcl))
    end function ptcl_iseven

    integer function get_nptcls( self )
        class(polarft_corrcalc), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    subroutine get_pinds( self, pinds )
        class(polarft_corrcalc), intent(inout) :: self
        integer, allocatable,    intent(out)   :: pinds(:)
        pinds = self%pinds
    end subroutine get_pinds

    integer function get_npix( self )
        class(polarft_corrcalc), intent(in) :: self
        get_npix = sum(nint(self%npix_per_shell(self%kfromto(1):self%kfromto(2))))
    end function get_npix

    ! returns pointer to temporary pft according to current thread
    subroutine get_work_pft_ptr( self, ptr )
        class(polarft_corrcalc), intent(in)  :: self
        complex(sp),   pointer,  intent(out) :: ptr(:,:)
        integer :: ithr
        ithr = omp_get_thread_num()+1
        ptr => self%heap_vars(ithr)%pft_ref_tmp
    end subroutine get_work_pft_ptr

    ! returns pointer to temporary pft-sized real matrix (eg CTF) according to current thread
    subroutine get_work_rpft_ptr( self, ptr )
        class(polarft_corrcalc), intent(in)  :: self
        real(sp),      pointer,  intent(out) :: ptr(:,:)
        integer :: ithr
        ithr = omp_get_thread_num()+1
        ptr => self%heap_vars(ithr)%pft_r
    end subroutine get_work_rpft_ptr

    subroutine get_ptcls_ptr( self, ptr )
        class(polarft_corrcalc), target, intent(in)  :: self
        complex(sp),            pointer, intent(out) :: ptr(:,:,:)
        ptr => self%pfts_ptcls
    end subroutine get_ptcls_ptr

    subroutine get_ctfmats_ptr( self, ptr )
        class(polarft_corrcalc), target, intent(in)  :: self
        real(sp),               pointer, intent(out) :: ptr(:,:,:)
        ptr => self%ctfmats
    end subroutine get_ctfmats_ptr

    subroutine get_refs_ptr( self, ptre, ptro )
        class(polarft_corrcalc), target, intent(in)  :: self
        complex(sp),            pointer, intent(out) :: ptre(:,:,:), ptro(:,:,:)
        ptre => self%pfts_refs_even
        ptro => self%pfts_refs_odd
    end subroutine get_refs_ptr

    ! return irefs for all states of the same reprojection direction
    subroutine get_linstates_irefs( self, iref, irefs )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: iref
        integer,                 intent(inout) :: irefs(params_glob%nstates)
        integer :: istate, iproj
        iproj = mod(iref-1, params_glob%nspace) + 1
        do istate = 1, params_glob%nstates
            irefs(istate) = (istate - 1) * params_glob%nspace + iproj
        enddo
    end subroutine get_linstates_irefs

    ! compute weights of irefs, given the current prefs
    subroutine get_linstates_prefs( self, iref, irefs, prefs )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: iref
        integer,                 intent(in)    :: irefs(params_glob%nstates)
        real,                    intent(inout) :: prefs(params_glob%nstates)
        real    :: prefs_ori( params_glob%nstates), A( params_glob%nstates, params_glob%nstates),&
                  &prefs_comp(params_glob%nstates), Ai(params_glob%nstates, params_glob%nstates)
        integer :: istate, j, errflg
        prefs_ori = prefs / sum(prefs)
        select case(trim(params_glob%linstates_mode))
            case('forprob')
                do istate = 1, params_glob%nstates
                    prefs(istate) = 1.
                    if( iref == irefs(istate) ) cycle
                    do j = 1, params_glob%nstates
                        if( j == istate ) cycle
                        prefs(istate) = prefs(istate) * (1. - prefs_ori(j))
                    enddo
                enddo
                prefs = prefs / sum(prefs)
            case('backprob')
                prefs_comp = 1.
                do istate = 1, params_glob%nstates
                    do j = 1, params_glob%nstates
                        if( j == istate ) cycle
                        prefs_comp(istate) = prefs_comp(istate) * (1. - prefs_ori(j))
                    enddo
                enddo
                A = 1.
                do istate = 1, params_glob%nstates
                    do j = 1, params_glob%nstates
                        if( j == istate ) cycle
                        A(istate, j) = prefs_comp(j)
                    enddo
                enddo
                call matinv(A, Ai, params_glob%nstates, errflg)
                do istate = 1, params_glob%nstates
                    if( iref == irefs(istate) )then
                        prefs = Ai(istate,:)
                        exit
                    endif
                enddo
        end select
    end subroutine get_linstates_prefs

    ! PRINTERS/VISUALISERS

    subroutine vis_ptcl( self, iptcl )
        use gnufor2, only: gnufor_image
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        call gnufor_image( real(self%pfts_ptcls(:,:,self%pinds(iptcl))), palette='gray')
        call gnufor_image(aimag(self%pfts_ptcls(:,:,self%pinds(iptcl))), palette='gray')
    end subroutine vis_ptcl

    subroutine vis_ref( self, iref, iseven )
        use gnufor2, only: gnufor_image
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iref
        logical,                 intent(in) :: iseven
        if( iseven )then
            call gnufor_image( real(self%pfts_refs_even(:,:,iref)), palette='gray')
            call gnufor_image(aimag(self%pfts_refs_even(:,:,iref)), palette='gray')
        else
            call gnufor_image( real(self%pfts_refs_odd(:,:,iref)), palette='gray')
            call gnufor_image(aimag(self%pfts_refs_odd(:,:,iref)), palette='gray')
        endif
    end subroutine vis_ref

    subroutine polar2cartesian_1( self, i, isref, cmat, box, box_in )
        class(polarft_corrcalc), intent(in)    :: self
        integer,                 intent(in)    :: i
        logical,                 intent(in)    :: isref
        complex,    allocatable, intent(inout) :: cmat(:,:)
        integer,                 intent(out)   :: box
        integer,    optional,    intent(in)    :: box_in
        if( isref )then
            call self%polar2cartesian_2(self%pfts_refs_even(:,:,i), cmat, box, box_in)
        else
            call self%polar2cartesian_2(self%pfts_ptcls(:,:,i),     cmat, box, box_in)
        endif
    end subroutine polar2cartesian_1

    subroutine polar2cartesian_2( self, cmat_in, cmat, box, box_in )
        class(polarft_corrcalc), intent(in)    :: self
        complex,                 intent(in)    :: cmat_in(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex,    allocatable, intent(inout) :: cmat(:,:)
        integer,                 intent(out)   :: box
        integer,    optional,    intent(in)    :: box_in
        integer,    allocatable :: norm(:,:)
        complex :: comp
        integer :: k,c,irot,physh,physk
        if( allocated(cmat) ) deallocate(cmat)
        box = 2*self%kfromto(2)
        if( present(box_in) ) box = box_in
        c   = box/2+1
        allocate(cmat(box/2+1,box),source=cmplx(0.0,0.0))
        allocate(norm(box/2+1,box),source=0)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%kfromto(2)
                ! Nearest-neighbour interpolation
                physh = nint(self%polar(irot,k)) + 1
                physk = nint(self%polar(irot+self%nrots,k)) + c
                if( physk > box ) cycle
                comp              = cmat_in(irot,k)
                cmat(physh,physk) = cmat(physh,physk) + comp
                norm(physh,physk) = norm(physh,physk) + 1
            end do
        end do
        ! normalization
        where(norm>0)
            cmat = cmat / real(norm)
        end where
        ! irot = self%pftsz+1, eg. angle=180.
        do k = 1,box/2-1
            cmat(1,k+c) = conjg(cmat(1,c-k))
        enddo
        ! arbitrary magnitude
        cmat(1,c) = (0.0,0.0)
    end subroutine polar2cartesian_2

    subroutine prefs_to_cartesian( self, cavgs, box_in, is_even )
        use simple_image
        class(polarft_corrcalc), intent(in)    :: self
        type(image),             intent(inout) :: cavgs(self%nrefs)
        integer, optional,       intent(out)   :: box_in
        logical, optional,       intent(in)    :: is_even
        integer, allocatable :: norm(:,:)
        complex, allocatable :: cmat(:,:)
        complex :: comp
        integer :: k,c,irot,physh,physk,box,iref
        logical :: l_even
        l_even = .true.
        if( present(is_even) ) l_even = is_even
        box = 2*self%kfromto(2)
        if( present(box_in) ) box = box_in
        c   = box/2+1
        allocate(cmat(box/2+1,box),source=zero)
        allocate(norm(box/2+1,box),source=0)
        do iref = 1, self%nrefs
            if( l_even )then
                if( all(real(self%pfts_refs_even) .eq. zero) )then
                    call cavgs(iref)%new([box,box,1],1.0)
                    call cavgs(iref)%zero
                endif
            else
                if( all(self%pfts_refs_odd .eq. zero) )then
                    call cavgs(iref)%new([box,box,1],1.0)
                    call cavgs(iref)%zero
                endif
            endif
            cmat = zero
            norm = 0
            do irot=1,self%pftsz
                do k=self%kfromto(1),self%kfromto(2)
                    ! Nearest-neighbour interpolation
                    physh = nint(self%polar(irot,k)) + 1
                    physk = nint(self%polar(irot+self%nrots,k)) + c
                    if( physk > box ) cycle
                    if( l_even )then
                        comp          = self%pfts_refs_even(irot,k,iref)
                    else
                        comp          = self%pfts_refs_odd( irot,k,iref)
                    endif
                    cmat(physh,physk) = cmat(physh,physk) + comp
                    norm(physh,physk) = norm(physh,physk) + 1
                end do
            end do
            ! normalization
            where(norm>0)
                cmat = cmat / real(norm)
            end where
            ! irot = self%pftsz+1, eg. angle=180.
            do k = 1,box/2-1
                cmat(1,k+c) = conjg(cmat(1,c-k))
            enddo
            ! arbitrary magnitude
            cmat(1,c) = (0.0,0.0)
            call cavgs(iref)%new([box,box,1],1.0)
            call cavgs(iref)%set_cmat(cmat)
            call cavgs(iref)%shift_phorig()
            call cavgs(iref)%ifft
        enddo
    end subroutine prefs_to_cartesian

    subroutine read_ptcls_ctfs( self, ibatch, success )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: ibatch
        logical,                 intent(inout) :: success
        integer :: funit,io_stat
        character(len=LONGSTRLEN) :: fname
        success = .true.
        fname   = trim(POLARIZED_PTCLS)//'_'//int2str(ibatch)//BIN_EXT
        if( .not.file_exists(trim(fname)) )then
            success = .false.
            return
        endif
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//trim(fname), io_stat)
        read(unit=funit,pos=1) self%pfts_ptcls
        call fclose(funit)
        fname = trim(POLARIZED_CTFS)//'_'//int2str(ibatch)//BIN_EXT
        if( .not.file_exists(trim(fname)) )then
            success = .false.
            return
        endif
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//trim(fname), io_stat)
        read(unit=funit,pos=1) self%ctfmats
        call fclose(funit)
    end subroutine read_ptcls_ctfs

    subroutine write_ptcls_ctfs( self, ibatch )
        class(polarft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: ibatch
        integer :: funit,io_stat
        character(len=LONGSTRLEN) :: fname
        fname = trim(POLARIZED_PTCLS)//'_'//int2str(ibatch)//BIN_EXT
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//trim(fname),io_stat)
        write(unit=funit,pos=1) self%pfts_ptcls
        call fclose(funit)
        fname = trim(POLARIZED_CTFS)//'_'//int2str(ibatch)//BIN_EXT
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("write_pft_array: "//trim(fname),io_stat)
        write(unit=funit,pos=1) self%ctfmats
        call fclose(funit)
    end subroutine write_ptcls_ctfs

    subroutine print( self )
        class(polarft_corrcalc), intent(in) :: self
        write(logfhandle,*) "total n particles in partition         (self%nptcls): ", self%nptcls
        write(logfhandle,*) "number of references                    (self%nrefs): ", self%nrefs
        write(logfhandle,*) "number of rotations                     (self%nrots): ", self%nrots
        write(logfhandle,*) "size of pft                             (self%pftsz): ", self%pftsz
        write(logfhandle,*) "logical dim. of original Cartesian image (self%ldim): ", self%ldim
    end subroutine print

    ! MODIFIERS

    subroutine shift_ptcl( self, iptcl, shvec)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        real(sp),                intent(in)    :: shvec(2)
        complex(sp), pointer :: shmat(:,:)
        integer :: ithr, i
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        shmat => self%heap_vars(ithr)%shmat
        call self%gen_shmat(ithr, shvec, shmat)
        self%pfts_ptcls(:,:,i) = self%pfts_ptcls(:,:,i) * shmat
    end subroutine shift_ptcl

    subroutine shift_ref( self, iref, shvec)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        real(sp),                intent(in)    :: shvec(2)
        complex(dp), pointer :: shmat(:,:)
        integer :: ithr
        ithr = omp_get_thread_num() + 1
        shmat => self%heap_vars(ithr)%shmat_8
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat)
        self%pfts_refs_even(:,:,iref) = cmplx(dcmplx(self%pfts_refs_even(:,:,iref)) * shmat)
        self%pfts_refs_odd( :,:,iref) = cmplx(dcmplx(self%pfts_refs_odd( :,:,iref)) * shmat)
    end subroutine shift_ref

    ! mirror pft about h (mirror about y of cartesian image)
    subroutine mirror_ref_pft( self, iref )
        class(polarft_corrcalc), target, intent(in) :: self
        integer,                         intent(in) :: iref
        integer :: i,j
        complex(sp), pointer :: pft(:,:) => null(), pftmirr(:,:) => null()
        pft     => self%pfts_refs_even(:,:,iref)
        pftmirr => self%pfts_refs_odd(:,:,iref)
        pftmirr(1,:) = conjg(pft(1,:))
        if( is_even(self%pftsz) )then
            do i = 2,self%pftsz/2
                j = self%pftsz-i+2
                pftmirr(i,:) = pft(j,:)
                pftmirr(j,:) = pft(i,:)
            enddo
            i = self%pftsz/2 + 1
            pftmirr(i,:) = pft(i,:)
        else
            do i = 2,(self%pftsz+1)/2
                j = self%pftsz-i+2
                pftmirr(i,:) = pft(j,:)
                pftmirr(j,:) = pft(i,:)
            enddo
        endif
    end subroutine mirror_ref_pft

    ! MEMOIZERS

    subroutine memoize_sqsum_ptcl( self, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        real(dp) :: sumsqk
        integer  :: i, ik
        logical  :: l_sigma
        i       = self%pinds(iptcl)
        l_sigma = associated(self%sigma2_noise)
        self%sqsums_ptcls(i)  = 0.d0
        self%ksqsums_ptcls(i) = 0.d0
        if( l_sigma ) self%wsqsums_ptcls(i) = 0.d0
        do ik = self%kfromto(1),self%kfromto(2)
            sumsqk                = sum(real(self%pfts_ptcls(:,ik,i)*conjg(self%pfts_ptcls(:,ik,i)),dp))
            self%sqsums_ptcls(i)  = self%sqsums_ptcls(i) + sumsqk
            sumsqk                = real(ik,dp) * sumsqk
            self%ksqsums_ptcls(i) = self%ksqsums_ptcls(i) + sumsqk
            if( l_sigma ) self%wsqsums_ptcls(i) = self%wsqsums_ptcls(i) + sumsqk / real(self%sigma2_noise(ik,iptcl),dp)
        enddo
    end subroutine memoize_sqsum_ptcl

    subroutine rotate_pft_1( self, pft, irot, pft_rot )
        class(polarft_corrcalc), intent(in)  :: self
        complex(dp),             intent(in)  :: pft(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)  :: irot
        complex(dp),             intent(out) :: pft_rot(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = conjg(pft(mid+1:self%pftsz,:))
            pft_rot(irot:self%pftsz,:) =       pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = conjg(pft)
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = conjg(pft(1    :mid,       :))
            pft_rot(1              :irot-self%pftsz-1,:) =       pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_pft_1

    subroutine rotate_pft_2( self, pft, irot, pft_rot )
        class(polarft_corrcalc), intent(in)  :: self
        complex(sp),             intent(in)  :: pft(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)  :: irot
        complex(sp),             intent(out) :: pft_rot(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = conjg(pft(mid+1:self%pftsz,:))
            pft_rot(irot:self%pftsz,:) =       pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = conjg(pft)
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = conjg(pft(1    :mid,       :))
            pft_rot(1              :irot-self%pftsz-1,:) =       pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_pft_2

    subroutine rotate_pft_3( self, pft, irot, pft_rot )
        class(polarft_corrcalc), intent(in)  :: self
        real(sp),                intent(in)  :: pft(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)  :: irot
        real(sp),                intent(out) :: pft_rot(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = pft(mid+1:self%pftsz,:)
            pft_rot(irot:self%pftsz,:) = pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = pft
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = pft(1    :mid,       :)
            pft_rot(1              :irot-self%pftsz-1,:) = pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_pft_3

    subroutine rotate_pft_4( self, pft, irot, pft_rot )
        class(polarft_corrcalc), intent(in)  :: self
        real(dp),                intent(in)  :: pft(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)  :: irot
        real(dp),                intent(out) :: pft_rot(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer :: mid
        if( irot == 1 )then
            pft_rot = pft
        elseif( irot >= 2 .and. irot <= self%pftsz )then
            mid = self%pftsz - irot + 1
            pft_rot(   1:irot-1,    :) = pft(mid+1:self%pftsz,:)
            pft_rot(irot:self%pftsz,:) = pft(    1:mid,       :)
        elseif( irot == self%pftsz + 1 )then
            pft_rot = pft
        else
            mid = self%nrots - irot + 1
            pft_rot(irot-self%pftsz:self%pftsz       ,:) = pft(1    :mid,       :)
            pft_rot(1              :irot-self%pftsz-1,:) = pft(mid+1:self%pftsz,:)
        endif
    end subroutine rotate_pft_4

    subroutine rotate_iref_1( self, iref, irot, sh )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        integer,                 intent(in)    :: irot
        real,                    intent(in)    :: sh(2)
        complex(dp), pointer :: pft_ref(:,:), pft_ref_tmp_8(:,:), shmat(:,:)
        integer :: ithr
        ithr = omp_get_thread_num() + 1
        pft_ref       => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat         => self%heap_vars(ithr)%shmat_8
        ! even
        pft_ref_tmp_8 = dcmplx(self%pfts_refs_even(:,:,iref))
        call self%gen_shmat_8(ithr, real(sh,dp), shmat)
        call self%rotate_pft(pft_ref_tmp_8, irot, pft_ref)
        pft_ref                       = pft_ref * shmat
        self%pfts_refs_even(:,:,iref) = cmplx(pft_ref)
        ! odd
        pft_ref_tmp_8                 = dcmplx(self%pfts_refs_odd(:,:,iref))
        call self%rotate_pft(pft_ref_tmp_8, irot, pft_ref)
        pft_ref                       = pft_ref * shmat
        self%pfts_refs_odd( :,:,iref) = cmplx(pft_ref)
    end subroutine rotate_iref_1

    subroutine rotate_iref_2( self, pft_ref, irot, sh, pft_ref_out )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(in)    :: pft_ref(:,:)
        integer,                 intent(in)    :: irot
        real,                    intent(in)    :: sh(2)
        complex(dp), pointer,    intent(out)   :: pft_ref_out(:,:)
        complex(dp), pointer :: shmat(:,:)
        integer :: ithr
        ithr  = omp_get_thread_num() + 1
        shmat => self%heap_vars(ithr)%shmat_8
        call self%rotate_pft(pft_ref, irot, pft_ref_out)
        call self%gen_shmat_8(ithr, real(sh,dp), shmat)
        pft_ref_out = pft_ref_out * shmat
    end subroutine rotate_iref_2

    subroutine calc_polar_ctf( self, iptcl, smpd, kv, cs, fraca, dfx, dfy, angast )
        use simple_ctf, only: ctf
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        real,                    intent(in)    :: smpd, kv, cs, fraca, dfx, dfy, angast
        type(ctf) :: tfun
        real(sp)  :: spaFreqSq_mat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)  :: ang_mat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)  :: inv_ldim(3),hinv,kinv
        integer   :: i,irot,k
        if( .not.allocated(self%ctfmats) )then
            allocate(self%ctfmats(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls), source=1.)
        endif
        if(.not. self%with_ctf ) return
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do default(shared) private(irot,k,hinv,kinv) schedule(static) proc_bind(close)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%kfromto(2)
                hinv = self%polar(irot,k) * inv_ldim(1)
                kinv = self%polar(irot+self%nrots,k) * inv_ldim(2)
                spaFreqSq_mat(irot,k) = hinv*hinv+kinv*kinv
                ang_mat(irot,k)       = atan2(self%polar(irot+self%nrots,k),self%polar(irot,k))
            end do
        end do
        !$omp end parallel do
        i = self%pinds(iptcl)
        if( i > 0 )then
            tfun   = ctf(smpd, kv, cs, fraca)
            call tfun%init(dfx, dfy, angast)
            self%ctfmats(:,:,i) = tfun%eval(spaFreqSq_mat(:,:), ang_mat(:,:), 0.0, .not.params_glob%l_wiener_part)
        endif
    end subroutine calc_polar_ctf

    subroutine setup_npix_per_shell( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: h,k,sh
        if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
        allocate(self%npix_per_shell(self%kfromto(1):self%kfromto(2)),source=0.0)
        do h = 0,self%kfromto(2)
            do k = -self%kfromto(2),self%kfromto(2)
                if( (h==0) .and. (k>0) ) cycle
                sh = nint(sqrt(real(h**2+k**2)))
                if( sh < self%kfromto(1) ) cycle
                if( sh > self%kfromto(2) ) cycle
                self%npix_per_shell(sh) = self%npix_per_shell(sh) + 1.0
            end do
        end do
    end subroutine setup_npix_per_shell

    subroutine memoize_ptcls( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ithr,i,k
        logical :: l_memoize_absptcl_ctf
        l_memoize_absptcl_ctf = trim(params_glob%sh_inv).eq.'yes'
        !$omp parallel do collapse(2) private(i,k,ithr) default(shared) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            do k = self%kfromto(1),self%kfromto(2)
                ithr = omp_get_thread_num() + 1
                ! FT(X.CTF)
                if( self%with_ctf )then
                    self%cvec2(ithr)%c(1:self%pftsz) = self%pfts_ptcls(:,k,i) * self%ctfmats(:,k,i)
                else
                    self%cvec2(ithr)%c(1:self%pftsz) = self%pfts_ptcls(:,k,i)
                endif
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ptcl_ctf(:,k,i) = self%cvec2(ithr)%c(1:self%pftsz+1)
                ! FT(CTF2)
                if( self%with_ctf )then
                    self%rvec1(ithr)%r(1:self%pftsz)            = self%ctfmats(:,k,i)*self%ctfmats(:,k,i)
                    self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                else
                    self%rvec1(ithr)%r = 1.0
                endif
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ctf2(:,k,i) = self%cvec1(ithr)%c(1:self%pftsz+1)
                if( l_memoize_absptcl_ctf )then
                    if( self%with_ctf )then
                        self%cvec2(ithr)%c(1:self%pftsz) = abs(self%pfts_ptcls(:,k,i)) * self%ctfmats(:,k,i)
                    else
                        self%cvec2(ithr)%c(1:self%pftsz) = abs(self%pfts_ptcls(:,k,i))
                    endif
                    self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                    call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                    self%ft_absptcl_ctf(:,k,i) = self%cvec2(ithr)%c(1:self%pftsz+1)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine memoize_ptcls

    subroutine memoize_refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: k, ithr, iref, istate, jref, irefs(params_glob%nstates)
        ! allocations
        call self%allocate_refs_memoization
        ! memoization
        !$omp parallel do collapse(2) private(iref,k,ithr) default(shared) proc_bind(close) schedule(static)
        do iref = 1,self%nrefs
            do k = self%kfromto(1),self%kfromto(2)
                ithr = omp_get_thread_num() + 1
                ! FT(REFeven)*
                self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_even(:,k,iref)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ref_even(:,k,iref) = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! FT(REFodd)*
                self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_odd(:,k,iref)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ref_odd(:,k,iref) = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! FT(REF2even)*
                self%rvec1(ithr)%r(           1:self%pftsz) = real(self%pfts_refs_even(:,k,iref)*conjg(self%pfts_refs_even(:,k,iref)))
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ref2_even(:,k,iref) = conjg(self%cvec1(ithr)%c(1:self%pftsz+1))
                ! FT(REF2odd)*
                self%rvec1(ithr)%r(           1:self%pftsz) = real(self%pfts_refs_odd(:,k,iref)*conjg(self%pfts_refs_odd(:,k,iref)))
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ref2_odd(:,k,iref) = conjg(self%cvec1(ithr)%c(1:self%pftsz+1))
            enddo
        enddo
        !$omp end parallel do
        if( params_glob%l_linstates )then
            !$omp parallel do private(iref,k,ithr,istate,irefs,jref) default(shared) proc_bind(close) schedule(static)
            do iref = 1,self%nrefs
                call self%get_linstates_irefs(iref,irefs)
                do istate = 1,params_glob%nstates
                    jref = irefs(istate)
                    do k = self%kfromto(1),self%kfromto(2)
                        ithr = omp_get_thread_num() + 1
                        ! FT(iREFeven x jREFeven)*
                        self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_even(:,k,iref)*conjg(self%pfts_refs_even(:,k,jref))
                        self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                        call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                        self%ft_refs2_even(:,k,iref,istate) = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                        ! FT(iREFodd x jREFodd)*
                        self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_odd(:,k,iref)*conjg(self%pfts_refs_odd(:,k,jref))
                        self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                        call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                        self%ft_refs2_odd(:,k,iref,istate) = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine memoize_refs

    subroutine kill_memoized_ptcls( self )
        class(polarft_corrcalc), intent(inout) :: self
        if( allocated(self%ft_ptcl_ctf) )    deallocate(self%ft_ptcl_ctf,self%ft_ctf2)
        if( allocated(self%ft_absptcl_ctf) ) deallocate(self%ft_absptcl_ctf)
    end subroutine kill_memoized_ptcls

    subroutine kill_memoized_refs( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: i
        if( allocated(self%ft_ref_even) )then
            do i = 1,size(self%cvec1,dim=1)
                call fftwf_free(self%cvec1(i)%p)
                call fftwf_free(self%cvec2(i)%p)
                call fftw_free( self%drvec(i)%p)
            enddo
            deallocate(self%ft_ref_even,self%ft_ref_odd,self%ft_ref2_even,self%ft_ref2_odd,&
            &self%rvec1,self%cvec1,self%cvec2,self%drvec)
            call fftwf_destroy_plan(self%plan_fwd1)
            call fftwf_destroy_plan(self%plan_bwd1)
            call fftwf_destroy_plan(self%plan_mem_r2c)
        endif
        if( allocated(self%ft_refs2_even) ) deallocate(self%ft_refs2_even,self%ft_refs2_odd)
    end subroutine kill_memoized_refs

    subroutine allocate_ptcls_memoization( self )
        class(polarft_corrcalc), intent(inout) :: self
        allocate(self%ft_ptcl_ctf(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls),&
                &self%ft_ctf2(    self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls))
        if( trim(params_glob%sh_inv).eq.'yes' )then
            allocate(self%ft_absptcl_ctf(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls))
        endif
    end subroutine allocate_ptcls_memoization

    subroutine allocate_refs_memoization( self )
        class(polarft_corrcalc), intent(inout) :: self
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname ! FFTW wisdoms (per part or suffer I/O lag)
        integer(kind=c_int) :: wsdm_ret
        integer             :: ithr
        if( allocated(self%ft_ref_even) ) call self%kill_memoized_refs
        allocate(self%ft_ref_even( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref_odd(  self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref2_even(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref2_odd( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%rvec1(nthr_glob), self%cvec1(nthr_glob),self%cvec2(nthr_glob),&
                &self%drvec(nthr_glob))
        if( params_glob%l_linstates )then
            allocate(self%ft_refs2_even(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs,params_glob%nstates),&
                    &self%ft_refs2_odd( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs,params_glob%nstates),&
                    &self%has_cache(self%nrefs,nthr_glob),self%cached_vals(self%nrefs,nthr_glob),self%do_cache(nthr_glob))
            self%do_cache = .false.
        endif
        ! convenience objects
        do ithr = 1,nthr_glob
            self%cvec1(ithr)%p = fftwf_alloc_complex(int(self%pftsz+1, c_size_t))
            self%cvec2(ithr)%p = fftwf_alloc_complex(int(self%nrots, c_size_t))
            call c_f_pointer(self%cvec1(ithr)%p, self%cvec1(ithr)%c, [self%pftsz+1])
            call c_f_pointer(self%cvec1(ithr)%p, self%rvec1(ithr)%r, [self%nrots+2])
            call c_f_pointer(self%cvec2(ithr)%p, self%cvec2(ithr)%c, [self%nrots])
            self%drvec(ithr)%p = fftw_alloc_real(int(self%nrots, c_size_t))
            call c_f_pointer(self%drvec(ithr)%p, self%drvec(ithr)%r, [self%nrots])
        enddo
        ! plans & FFTW3 wisdoms
        if( params_glob%l_distr_exec )then
            allocate(fft_wisdoms_fname, source='fft_wisdoms_part'//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'//c_null_char)
        else
            allocate(fft_wisdoms_fname, source='fft_wisdoms.dat'//c_null_char)
        endif
        wsdm_ret = fftw_import_wisdom_from_filename(fft_wisdoms_fname)
        self%plan_fwd1    = fftwf_plan_dft_1d(    self%nrots, self%cvec2(1)%c, self%cvec2(1)%c, FFTW_FORWARD, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_bwd1    = fftwf_plan_dft_c2r_1d(self%nrots, self%cvec1(1)%c, self%rvec1(1)%r,               ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_mem_r2c = fftwf_plan_dft_r2c_1d(self%nrots, self%rvec1(1)%r, self%cvec1(1)%c,               ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        wsdm_ret = fftw_export_wisdom_to_filename(fft_wisdoms_fname)
        deallocate(fft_wisdoms_fname)
        if (wsdm_ret == 0) then
            write (*, *) 'Error: could not write FFTW3 wisdom file! Check permissions.'
        end if
    end subroutine allocate_refs_memoization

    ! CALCULATORS

    subroutine create_polar_absctfmats( self, spproj, oritype, pfromto )
        use simple_ctf,        only: ctf
        use simple_sp_project, only: sp_project
        class(polarft_corrcalc),   intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in)    :: oritype
        integer, optional,         intent(in)    :: pfromto(2)
        type(ctfparams) :: ctfparms(nthr_glob)
        type(ctf)       :: tfuns(nthr_glob)
        real(sp)        :: spaFreqSq_mat(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(sp)        :: ang_mat(self%pftsz,self%kfromto(1):self%kfromto(2)), hinv,kinv
        integer         :: i,irot,k,iptcl,ithr,ppfromto(2),ctfmatind
        logical         :: present_pfromto
        present_pfromto = present(pfromto)
        ppfromto = self%pfromto
        if( present_pfromto ) ppfromto = pfromto
        if(.not. self%with_ctf ) return
        !$omp parallel do default(shared) private(irot,k,hinv,kinv) schedule(static) proc_bind(close)
        do irot=1,self%pftsz
            do k=self%kfromto(1),self%kfromto(2)
                hinv = self%polar(irot,k) / self%ldim(1)
                kinv = self%polar(irot+self%nrots,k) / self%ldim(2)
                spaFreqSq_mat(irot,k) = hinv*hinv+kinv*kinv
                ang_mat(irot,k)       = atan2(self%polar(irot+self%nrots,k),self%polar(irot,k))
            end do
        end do
        !$omp end parallel do
        if( params_glob%l_wiener_part )then
            ! taking into account CTF is intact before limit
            !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr) schedule(static) proc_bind(close)
            do i=ppfromto(1),ppfromto(2)
                if( .not. present_pfromto )then
                    iptcl     = i
                    ctfmatind = i
                else
                    iptcl     = i
                    ctfmatind = i - ppfromto(1) + 1
                endif
                if( self%pinds(iptcl) > 0 )then
                    ithr           = omp_get_thread_num() + 1
                    ctfparms(ithr) = spproj%get_ctfparams(trim(oritype), iptcl)
                    tfuns(ithr)    = ctf(ctfparms(ithr)%smpd, ctfparms(ithr)%kv, ctfparms(ithr)%cs, ctfparms(ithr)%fraca)
                    call tfuns(ithr)%init(ctfparms(ithr)%dfx, ctfparms(ithr)%dfy, ctfparms(ithr)%angast)
                    if( ctfparms(ithr)%l_phaseplate )then
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), ctfparms(ithr)%phshift, .false. ))
                    else
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), 0.0,                    .false.))
                    endif
                endif
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr) schedule(static) proc_bind(close)
            do i=ppfromto(1),ppfromto(2)
                if( .not. present_pfromto )then
                    iptcl     = i
                    ctfmatind = i
                else
                    iptcl     = i
                    ctfmatind = i - ppfromto(1) + 1
                endif
                if( self%pinds(iptcl) > 0 )then
                    ithr           = omp_get_thread_num() + 1
                    ctfparms(ithr) = spproj%get_ctfparams(trim(oritype), iptcl)
                    tfuns(ithr)    = ctf(ctfparms(ithr)%smpd, ctfparms(ithr)%kv, ctfparms(ithr)%cs, ctfparms(ithr)%fraca)
                    call tfuns(ithr)%init(ctfparms(ithr)%dfx, ctfparms(ithr)%dfy, ctfparms(ithr)%angast)
                    if( ctfparms(ithr)%l_phaseplate )then
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:), ctfparms(ithr)%phshift) )
                    else
                        self%ctfmats(:,:,self%pinds(ctfmatind)) = abs(tfuns(ithr)%eval(spaFreqSq_mat(:,:), ang_mat(:,:)))
                    endif
                endif
            end do
            !$omp end parallel do
        endif
    end subroutine create_polar_absctfmats

    !>  Generate polar shift matrix by means of de Moivre's formula, double precision
    subroutine gen_shmat_8( self, ithr, shift_8 , shmat_8 )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        real(dp),                intent(in)    :: shift_8(2)
        complex(dp),    pointer, intent(inout) :: shmat_8(:,:)
        integer :: k
        ! first shell, analytic
        self%heap_vars(ithr)%argvec = self%argtransf(:self%pftsz,  self%kfromto(1)) * shift_8(1) +&
                                    & self%argtransf(self%pftsz+1:,self%kfromto(1)) * shift_8(2)
        shmat_8(:,self%kfromto(1))  = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
        ! one shell to the next
        self%heap_vars(ithr)%argvec = self%argtransf_shellone(:self%pftsz)   * shift_8(1) +&
                                    & self%argtransf_shellone(self%pftsz+1:) * shift_8(2)
        self%heap_vars(ithr)%shvec  = dcmplx(dcos(self%heap_vars(ithr)%argvec), dsin(self%heap_vars(ithr)%argvec))
        ! remaining shells, cos(kx)+isin(kx) = (cos(x)+isin(x))**k-1 * (cos(x)+isin(x))
        do k = self%kfromto(1)+1,self%kfromto(2)
            shmat_8(:,k) = shmat_8(:,k-1) * self%heap_vars(ithr)%shvec
        enddo
        ! alternative to:
        ! argmat  => self%heap_vars(ithr)%argmat_8
        ! argmat  =  self%argtransf(:self%pftsz,:)*shvec(1) + self%argtransf(self%pftsz + 1:,:)*shvec(2)
        ! shmat   =  cmplx(cos(argmat),sin(argmat),dp)
    end subroutine gen_shmat_8

    subroutine calc_shift( self, iref, iptcl, sh, rot_in )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref
        integer,                 intent(in)    :: iptcl
        real,                    intent(inout) :: sh(2)
        integer,  optional,      intent(in)    :: rot_in
        complex,  pointer   :: pft_ref(:,:), pft_ref_tmp(:,:)
        real(dp), pointer   :: args1(:,:), args2(:,:)
        integer,  parameter :: NPI   = 5,&       ! maximum number of trigonometry periods
                              &NLINS = 2,&       ! maximum number of 2x2 linear systems
                              &NEQS  = 2*NLINS
        integer :: ind, ithr, i, j, k, irot, cand1_cnt, cand2_cnt, eq_cnt, k_cands(NEQS), r_cands(NEQs)
        logical :: rots(self%pftsz), ks(self%kfromto(1):self%kfromto(2))
        complex :: AB
        real    :: C_1st, T_1st, C_2nd, T_2nd, RHS_1st(NPI), RHS_2nd(NPI), x_comp, y_comp, trs, abspft,&
                &cand1_x(NPI*NPI), cand1_y(NPI*NPI), cand2_x(NPI*NPI), cand2_y(NPI*NPI), minval, minx, miny, val
        ind  = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref
        pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp
        args1       => self%heap_vars(ithr)%pft_r1_8
        args2       => self%heap_vars(ithr)%pft_r2_8
        if( present(rot_in) )then
            if( self%iseven(ind) )then
                pft_ref_tmp = self%pfts_refs_even(:,:,iref)
            else
                pft_ref_tmp = self%pfts_refs_odd(:,:,iref)
            endif
            call self%rotate_pft(pft_ref_tmp, rot_in, pft_ref)
            call self%rotate_pft(self%argtransf(           1:self%pftsz,:), rot_in, args1)
            call self%rotate_pft(self%argtransf(self%pftsz+1:          ,:), rot_in, args2)
        else
            if( self%iseven(ind) )then
                pft_ref = self%pfts_refs_even(:,:,iref)
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)
            endif
            args1 = self%argtransf(           1:self%pftsz,:)
            args2 = self%argtransf(self%pftsz+1:          ,:)
        endif
        if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(:,:,ind)
        trs    = params_glob%trs
        eq_cnt = 0
        ! generating candidates for equations
        rots = .false.
        ks   = .false.
        do k = self%kfromto(1), self%kfromto(2)
            if( ks(k) )cycle
            do irot = 1, self%pftsz
                if( rots(irot) .or. ks(k) )cycle
                if( eq_cnt > NEQS )exit
                abspft = real(self%pfts_ptcls(irot,k,ind) * conjg(self%pfts_ptcls(irot,k,ind)))
                if( abspft < TINY )cycle
                abspft = real(pft_ref(irot,k) * conjg(pft_ref(irot,k)))
                if( abspft < TINY )cycle
                eq_cnt          = eq_cnt + 1
                k_cands(eq_cnt) = k
                r_cands(eq_cnt) = irot
                rots(irot)      = .true.
                ks(k)           = .true.
            enddo
        enddo
        ! not enough candidates to calculate the shifts
        if( eq_cnt < NEQS )then
            sh = 0.
            return
        endif
        ! first pair
        k       = k_cands(1)
        irot    = r_cands(1)
        AB      = pft_ref(irot,k)/self%pfts_ptcls(irot,k,ind)
        RHS_1st = atan(aimag(AB)/real(AB))
        C_1st   = args1(irot,k)
        T_1st   = args2(irot,k)
        do j = 1, NPI
            RHS_1st(j) = RHS_1st(j) + (j-1-NPI/2)*PI
        enddo
        k       = k_cands(2)
        irot    = r_cands(2)
        AB      = pft_ref(irot,k)/self%pfts_ptcls(irot,k,ind)
        RHS_2nd = atan(aimag(AB)/real(AB))
        C_2nd   = args1(irot,k)
        T_2nd   = args2(irot,k)
        do j = 1, NPI
            RHS_2nd(j) = RHS_2nd(j) + (j-1-NPI/2)*PI
        enddo
        cand1_cnt = 0
        do i = 1, NPI
            do j = 1, NPI
                x_comp = -(RHS_1st(i)*T_2nd - RHS_2nd(j)*T_1st)/(T_1st*C_2nd - T_2nd*C_1st)
                y_comp =  (RHS_1st(i)*C_2nd - RHS_2nd(j)*C_1st)/(T_1st*C_2nd - T_2nd*C_1st)
                if( x_comp < -trs .or. x_comp > trs .or. y_comp < -trs .or. y_comp > trs )cycle
                cand1_cnt          = cand1_cnt + 1
                cand1_x(cand1_cnt) = x_comp
                cand1_y(cand1_cnt) = y_comp
            enddo
        enddo
        ! another pair
        k       = k_cands(3)
        irot    = r_cands(3)
        AB      = pft_ref(irot,k)/self%pfts_ptcls(irot,k,ind)
        RHS_1st = atan(aimag(AB)/real(AB))
        C_1st   = args1(irot,k)
        T_1st   = args2(irot,k)
        do j = 1, NPI
            RHS_1st(j) = RHS_1st(j) + (j-1-NPI/2)*PI
        enddo
        k       = k_cands(4)
        irot    = r_cands(4)
        AB      = pft_ref(irot,k)/self%pfts_ptcls(irot,k,ind)
        RHS_2nd = atan(aimag(AB)/real(AB))
        C_2nd   = args1(irot,k)
        T_2nd   = args2(irot,k)
        do j = 1, NPI
            RHS_2nd(j) = RHS_2nd(j) + (j-1-NPI/2)*PI
        enddo
        cand2_cnt = 0
        do i = 1, NPI
            do j = 1, NPI
                x_comp = -(RHS_1st(i)*T_2nd - RHS_2nd(j)*T_1st)/(T_1st*C_2nd - T_2nd*C_1st)
                y_comp =  (RHS_1st(i)*C_2nd - RHS_2nd(j)*C_1st)/(T_1st*C_2nd - T_2nd*C_1st)
                if( x_comp < -trs .or. x_comp > trs .or. y_comp < -trs .or. y_comp > trs )cycle
                cand2_cnt          = cand2_cnt + 1
                cand2_x(cand2_cnt) = x_comp
                cand2_y(cand2_cnt) = y_comp
            enddo
        enddo
        ! finding the min between cand1 and cand2
        minval = huge(minval)
        do i = 1, cand1_cnt
            do j = 1, cand2_cnt
                val = sqrt((cand1_x(i)-cand2_x(j))**2 + (cand1_y(i)-cand2_y(j))**2)
                if( val < minval )then
                    minval = val
                    minx   = (cand1_x(i)+cand2_x(j)) / 2.
                    miny   = (cand1_y(i)+cand2_y(j)) / 2.
                endif
            enddo
        enddo
        sh = -[minx, miny]
    end subroutine calc_shift

    !>  Generate shift matrix following de Moivre's formula, single precision
    subroutine gen_shmat( self, ithr, shift, shmat )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: ithr
        real(sp),                intent(in)    :: shift(2)
        complex(sp),    pointer, intent(inout) :: shmat(:,:)
        call self%gen_shmat_8(ithr, real(shift,dp), self%heap_vars(ithr)%shmat_8)
        shmat = cmplx(self%heap_vars(ithr)%shmat_8)
    end subroutine gen_shmat

    ! Benchmarck for correlation calculation
    ! Is not FFT-accelerated, does not rely on memoization, for reference only
    real function calc_corr_rot_shift( self, iref, iptcl, shvec, irot, kweight )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        logical,       optional, intent(in)    :: kweight
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_rot_ref(:,:)
        real(dp)    :: sqsumref, sqsumptcl, num
        integer     :: i, k, ithr
        logical     :: kw
        kw = .true.
        if( present(kweight) ) kw = kweight
        calc_corr_rot_shift = 0.
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, real(shvec,dp),shmat)
        pft_ref = pft_ref * shmat
        call self%rotate_pft(pft_ref, irot, pft_rot_ref)
        if( self%with_ctf ) pft_rot_ref = pft_rot_ref * self%ctfmats(:,:,i)
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            sqsumref  = 0.d0
            sqsumptcl = 0.d0
            num       = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                if( kw )then
                    sqsumptcl = sqsumptcl + real(k,dp) * real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
                    sqsumref  = sqsumref  + real(k,dp) * real(sum(pft_rot_ref(:,k) * conjg(pft_rot_ref(:,k))),dp)
                    num       = num       + real(k,dp) * real(sum(pft_rot_ref(:,k) * conjg(self%pfts_ptcls(:,k,i))),dp)
                else
                    sqsumptcl = sqsumptcl + real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
                    sqsumref  = sqsumref  + real(sum(pft_rot_ref(:,k) * conjg(pft_rot_ref(:,k))),dp)
                    num       = num       + real(sum(pft_rot_ref(:,k) * conjg(self%pfts_ptcls(:,k,i))),dp)
                endif
            enddo
            calc_corr_rot_shift = real(num/sqrt(sqsumref*sqsumptcl))
        case(OBJFUN_EUCLID)
            pft_rot_ref = pft_rot_ref - self%pfts_ptcls(:,:,i)
            sqsumptcl = 0.d0
            num       = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                if( kw )then
                    sqsumptcl = sqsumptcl + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    num       = num       + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(csq_fast(pft_rot_ref(:,k)))
                else
                    sqsumptcl = sqsumptcl + (1.d0 / self%sigma2_noise(k,iptcl)) * sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    num       = num       + (1.d0 / self%sigma2_noise(k,iptcl)) * sum(csq_fast(pft_rot_ref(:,k)))
                endif
            end do
            calc_corr_rot_shift = real(exp( -num / sqsumptcl ))
        end select
    end function calc_corr_rot_shift

    real function calc_magcorr_rot( self, iref, iptcl, irot, kweight )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        logical,       optional, intent(in)    :: kweight
        complex(dp), pointer :: pft_ref(:,:), pft_rot_ref(:,:)
        real(dp),    pointer :: mag_rot_ref(:,:)
        real(dp)    :: sqsumref, sqsumptcl, num
        integer     :: i, k, ithr
        logical     :: kw
        kw = .true.
        if( present(kweight) ) kw = kweight
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        calc_magcorr_rot = 0.
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        mag_rot_ref => self%heap_vars(ithr)%pft_r1_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%rotate_pft(pft_ref, irot, pft_rot_ref)
        if( self%with_ctf ) pft_rot_ref = pft_rot_ref * self%ctfmats(:,:,i)
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            mag_rot_ref = abs(pft_rot_ref)
            sqsumref    = 0.d0
            sqsumptcl   = 0.d0
            num         = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                if( kw )then
                    sqsumptcl = sqsumptcl + real(k,dp) * sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    sqsumref  = sqsumref  + real(k,dp) * sum(mag_rot_ref(:,k)**2)
                    num       = num       + real(k,dp) * sum(mag_rot_ref(:,k) * real(abs(self%pfts_ptcls(:,k,i)),dp))
                else
                    sqsumptcl = sqsumptcl + sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    sqsumref  = sqsumref  + sum(mag_rot_ref(:,k)**2)
                    num       = num       + sum(mag_rot_ref(:,k) * real(abs(self%pfts_ptcls(:,k,i)),dp))
                endif
            enddo
            calc_magcorr_rot = real(num/sqrt(sqsumref*sqsumptcl))
        case(OBJFUN_EUCLID)
            ! not implemented
        end select
    end function calc_magcorr_rot

    subroutine gencorrs_mag_cc( self, iref, iptcl, ccs, kweight )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real,                    intent(inout) :: ccs(self%pftsz)
        logical,       optional, intent(in)    :: kweight
        complex(dp), pointer :: pft_ref(:,:)
        real(dp),    pointer :: pft_mag_ptcl(:,:)
        real(dp)    :: sumsqptcl, sumsqref
        integer     :: i, k, ithr
        logical     :: kw
        kw = .true.
        if( present(kweight) ) kw = kweight
        i            = self%pinds(iptcl)
        ithr         = omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref_8
        pft_mag_ptcl => self%heap_vars(ithr)%pft_r1_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        pft_mag_ptcl(:,:)           = real(abs(self%pfts_ptcls(:,:,i)),dp)
        sumsqptcl                   = 0.d0
        self%heap_vars(ithr)%kcorrs = 0.d0
        if( self%with_ctf )then
            self%drvec(ithr)%r = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! |X|2
                if( kw )then
                    sumsqptcl = sumsqptcl + real(k,dp) * sum(pft_mag_ptcl(:,k)**2)
                else
                    sumsqptcl = sumsqptcl + sum(pft_mag_ptcl(:,k)**2)
                endif
                ! FT(CTF2)
                self%rvec1(ithr)%r(1:self%pftsz) = self%ctfmats(:,k,i)**2
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%cvec2(ithr)%c(1:self%pftsz+1) = self%cvec1(ithr)%c
                ! FT(|REF|2)
                self%rvec1(ithr)%r(1:self%pftsz)            = real(pft_ref(:,k)*conjg(pft_ref(:,k)),sp)
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                ! FT(CTF2) x FT(|REF|2)*
                self%cvec1(ithr)%c = self%cvec2(ithr)%c(1:self%pftsz+1) * conjg(self%cvec1(ithr)%c)
                ! IFFT( FT(CTF2) x FT(|REF|2)* )
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                if( kw )then
                    self%drvec(ithr)%r(1:self%pftsz) = self%drvec(ithr)%r(1:self%pftsz) + real(k,dp) * real(self%rvec1(ithr)%r(1:self%pftsz),dp)
                else
                    self%drvec(ithr)%r(1:self%pftsz) = self%drvec(ithr)%r(1:self%pftsz) + real(self%rvec1(ithr)%r(1:self%pftsz),dp)
                endif
                ! FT(|REF|)*
                self%rvec1(ithr)%r(1:self%pftsz) = real(abs(pft_ref(:,k)), kind=sp)
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                ! FT(|X|.CTF) x FT(|REF|)*
                self%cvec1(ithr)%c = self%ft_absptcl_ctf(:,k,i) * conjg(self%cvec1(ithr)%c)
                ! IFFT( FT(|X|.CTF) x FT(|REF|)* )
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                if( kw )then
                    self%heap_vars(ithr)%kcorrs(1:self%pftsz) = self%heap_vars(ithr)%kcorrs(1:self%pftsz) +&
                        &real(k,dp) * real(self%rvec1(ithr)%r(1:self%pftsz),dp)
                else
                    self%heap_vars(ithr)%kcorrs(1:self%pftsz) = self%heap_vars(ithr)%kcorrs(1:self%pftsz) +&
                        &real(self%rvec1(ithr)%r(1:self%pftsz),dp)
                endif
            end do
            self%drvec(ithr)%r(1:self%pftsz) = self%drvec(ithr)%r(1:self%pftsz) * (sumsqptcl * real(2*self%nrots,dp))
            ccs = real(self%heap_vars(ithr)%kcorrs(1:self%pftsz) / dsqrt(self%drvec(ithr)%r(1:self%pftsz)))
        else
            sumsqref = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! FT(|REF|)*
                self%rvec1(ithr)%r(1:self%pftsz) = real(abs(pft_ref(:,k)), kind=sp)
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                ! FT(|X|) x FT(|REF|)*
                self%cvec1(ithr)%c = self%ft_absptcl_ctf(:,k,i) * conjg(self%cvec1(ithr)%c)
                ! IFFT( FT(|X|) x FT(|REF|)* )
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                if( kw )then
                    ! |X|2 & |REF|2
                    sumsqptcl = sumsqptcl + real(k,dp) * sum(pft_mag_ptcl(:,k)**2)
                    sumsqref  = sumsqref  + real(k,dp) * sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
                    self%heap_vars(ithr)%kcorrs(1:self%pftsz) = self%heap_vars(ithr)%kcorrs(1:self%pftsz) +&
                        &real(k,dp) * real(self%rvec1(ithr)%r(1:self%pftsz),dp)
                else
                    sumsqptcl = sumsqptcl + sum(pft_mag_ptcl(:,k)**2)
                    sumsqref  = sumsqref  + sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
                    self%heap_vars(ithr)%kcorrs(1:self%pftsz) = self%heap_vars(ithr)%kcorrs(1:self%pftsz) +&
                        &real(self%rvec1(ithr)%r(1:self%pftsz),dp)
                endif
            end do
            ccs = real(self%heap_vars(ithr)%kcorrs(1:self%pftsz) / (real(2*self%nrots,dp)*dsqrt(sumsqptcl*sumsqref)))
        endif
    end subroutine gencorrs_mag_cc

    subroutine gencorrs_mag( self, iref, iptcl, ccs, kweight )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real,                    intent(inout) :: ccs(self%pftsz)
        logical,       optional, intent(in)    :: kweight
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            call self%gencorrs_mag_cc( iref, iptcl, ccs, kweight )
        case(OBJFUN_EUCLID)
            ! unsupported
            ccs = 0.
        end select
    end subroutine gencorrs_mag

    subroutine calc_frc( self, iref, iptcl, irot, shvec, frc )
        class(polarft_corrcalc),  intent(inout) :: self
        integer,                  intent(in)    :: iref, iptcl, irot
        real(sp),                 intent(in)    :: shvec(2)
        real(sp),                 intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_rot_ref(:,:)
        real(dp) :: sumsqref, sumsqptcl, denom, num
        integer  :: k, ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat)
        pft_ref = pft_ref * shmat
        call self%rotate_pft(pft_ref, irot, pft_rot_ref)
        if( self%with_ctf ) pft_rot_ref = pft_rot_ref * real(self%ctfmats(:,:,i),dp)
        do k = self%kfromto(1),self%kfromto(2)
            num       = real(sum(pft_rot_ref(:,k)       * conjg(self%pfts_ptcls(:,k,i))),dp)
            sumsqptcl = real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
            sumsqref  = real(sum(pft_rot_ref(:,k)       * conjg(pft_rot_ref(:,k))),dp)
            denom     = sumsqptcl * sumsqref
            if( denom < 1.d-16 )then
                frc(k) = 0.
            else
                frc(k) = real(num / sqrt(denom))
            endif
        end do
    end subroutine calc_frc

    subroutine set_cache( self, val )
        class(polarft_corrcalc), intent(inout) :: self
        logical,                 intent(in)    :: val
        integer :: ithr
        ithr = omp_get_thread_num() + 1
        self%do_cache(ithr) = val
    end subroutine set_cache

    subroutine reset_cache( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ithr
        if( .not. params_glob%l_linstates ) return
        ithr = omp_get_thread_num() + 1
        self%do_cache(ithr)      = .true.
        self%has_cache(:,  ithr) = .false.
        self%cached_vals(:,ithr) = 0.
    end subroutine reset_cache

    pure logical function is_with_ctf( self )
        class(polarft_corrcalc), intent(in) :: self
        is_with_ctf = self%with_ctf
    end function is_with_ctf

    subroutine linear_ft_refs( self, iptcl, irefs, prefs, ft_ref, pft_ref_or_ft_ref )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        complex(sp),    pointer, intent(inout) :: ft_ref(:,:)
        logical,       optional, intent(in)    :: pft_ref_or_ft_ref     ! true for pft_ref
        integer  :: i, ithr, j, jref
        logical  :: even, l_pft
        l_pft = .false.
        if( present(pft_ref_or_ft_ref) ) l_pft = pft_ref_or_ft_ref
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        even  = self%iseven(i)
        if( l_pft )then
            ft_ref => self%heap_vars(ithr)%pft_ref
        else
            ft_ref => self%heap_vars(ithr)%pft_tmp
        endif
        ft_ref = cmplx(0.)
        do j = 1, size(irefs)
            jref = irefs(j)
            if( even )then
                if( l_pft )then
                    ft_ref = ft_ref + prefs(j) * self%pfts_refs_even(:,:,jref)
                else
                    ft_ref = ft_ref + prefs(j) * self%ft_ref_even(:,:,jref)
                endif
            else
                if( l_pft )then
                    ft_ref = ft_ref + prefs(j) * self%pfts_refs_odd(:,:,jref)
                else
                    ft_ref = ft_ref + prefs(j) * self%ft_ref_odd(:,:,jref)
                endif
            endif
        enddo
    end subroutine linear_ft_refs

    subroutine linear_ft_ref2s_1( self, iptcl, irefs, prefs, k, ft_ref2 )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: k
        complex(sp),    pointer, intent(inout) :: ft_ref2(:,:)
        integer  :: i, ithr, j, r, jref, rref, istate
        logical  :: even
        real     :: pjr
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        even         = self%iseven(i)
        ft_ref2(:,k) = cmplx(0.)
        do j = 1, size(irefs)
            jref = irefs(j)
            do r = 1, size(irefs)
                rref   = irefs(r)
                pjr    = prefs(j) * prefs(r)
                istate = floor(real(rref-1) / real(params_glob%nspace)) + 1
                if( even )then
                    ! ft_refs2_even is cached when linstates is on
                    if( params_glob%l_linstates )then
                        ft_ref2(:,k) = ft_ref2(:,k) + pjr * self%ft_refs2_even(:,k,jref,istate)
                    else
                        ! FT(jREFeven x rREFeven)*
                        self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_even(:,k,jref)*conjg(self%pfts_refs_even(:,k,rref))
                        self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                        call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                        ft_ref2(:,k) = ft_ref2(:,k) + pjr * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                    endif
                else
                    ! ft_refs2_odd is cached when linstates is on
                    if( params_glob%l_linstates )then
                        ft_ref2(:,k) = ft_ref2(:,k) + pjr * self%ft_refs2_odd(:,k,jref,istate)
                    else
                        ! FT(jREFodd x rREFodd)*
                        self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_odd(:,k,jref)*conjg(self%pfts_refs_odd(:,k,rref))
                        self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                        call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                        ft_ref2(:,k) = ft_ref2(:,k) + pjr * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                    endif
                endif
            enddo
        enddo
    end subroutine linear_ft_ref2s_1

    subroutine linear_ft_ref2s_2( self, iptcl, irefs, prefs, ft_ref2 )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        complex(sp),    pointer, intent(inout) :: ft_ref2(:,:)
        integer  :: i, j, r, jref, rref, istate
        logical  :: even
        real     :: pjr
        i       = self%pinds(iptcl)
        even    = self%iseven(i)
        ft_ref2 = cmplx(0.)
        do r = 1, size(irefs)
            rref   = irefs(r)
            istate = floor(real(rref-1) / real(params_glob%nspace)) + 1
            do j = 1, size(irefs)
                jref = irefs(j)
                pjr  = prefs(j) * prefs(r)
                if( even )then
                    ft_ref2 = ft_ref2 + pjr * self%ft_refs2_even(:,:,jref,istate)
                else
                    ft_ref2 = ft_ref2 + pjr * self%ft_refs2_odd(:,:,jref,istate)
                endif
            enddo
        enddo
    end subroutine linear_ft_ref2s_2

    subroutine gencorrs_1( self, iref, iptcl, cc, kweight, onestate )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        logical,       optional, intent(in)    :: kweight
        logical,       optional, intent(in)    :: onestate
        logical :: kw, l_onestate
        integer :: istate, irefs(params_glob%nstates), ithr
        real    :: prefs(params_glob%nstates)
        l_onestate = .false.
        if( present(onestate) ) l_onestate = onestate
        if( params_glob%l_linstates .and. .not.(l_onestate) )then
            ithr = omp_get_thread_num() + 1
            call self%get_linstates_irefs(iref, irefs)
            do istate = 1, params_glob%nstates
                if( self%do_cache(ithr) )then
                    if( self%has_cache(irefs(istate), ithr) )then
                        prefs(istate) = self%cached_vals(irefs(istate), ithr)
                        cycle
                    endif
                endif
                call self%gencorrs_1(irefs(istate), iptcl, cc, kweight, onestate=.true.)
                prefs(istate) = maxval(cc)
                if( self%do_cache(ithr) )then
                    self%has_cache(  irefs(istate), ithr) = .true.
                    self%cached_vals(irefs(istate), ithr) = prefs(istate)
                endif
            enddo
            call self%get_linstates_prefs(iref, irefs, prefs)
            call self%gencorrs_3( irefs, prefs, iptcl, cc, kweight )
        else
            select case(params_glob%cc_objfun)
                case(OBJFUN_CC)
                    kw = params_glob%l_kweight
                    if( present(kweight)) kw = kweight ! overrides params_glob%l_kweight
                    if( kw )then
                        call self%gencorrs_weighted_cc(iptcl, iref, cc)
                    else
                        call self%gencorrs_cc(iptcl, iref, cc)
                    endif
                case(OBJFUN_EUCLID)
                    call self%gencorrs_euclid(iptcl, iref, cc)
            end select
        endif
    end subroutine gencorrs_1

    subroutine gencorrs_2( self, iref, iptcl, shift, cc, kweight, onestate )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shift(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        logical,       optional, intent(in)    :: kweight
        logical,       optional, intent(in)    :: onestate
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer :: i, ithr
        logical :: kw, l_onestate
        integer :: istate, irefs(params_glob%nstates)
        real    :: prefs(params_glob%nstates)
        ithr       = omp_get_thread_num() + 1
        l_onestate = .false.
        if( present(onestate) ) l_onestate = onestate
        if( params_glob%l_linstates .and. .not.(l_onestate) )then
            call self%get_linstates_irefs(iref, irefs)
            do istate = 1, params_glob%nstates
                if( self%do_cache(ithr) )then
                    if( self%has_cache(irefs(istate), ithr) )then
                        prefs(istate) = self%cached_vals(irefs(istate), ithr)
                        cycle
                    endif
                endif
                call self%gencorrs_2(irefs(istate), iptcl, shift, cc, kweight, onestate=.true.)
                prefs(istate) = maxval(cc)
                if( self%do_cache(ithr) )then
                    self%has_cache(  irefs(istate), ithr) = .true.
                    self%cached_vals(irefs(istate), ithr) = prefs(istate)
                endif
            enddo
            call self%get_linstates_prefs(iref, irefs, prefs)
            call self%gencorrs_4( irefs, prefs, iptcl, shift, cc, kweight )
        else
            i       = self%pinds(iptcl)
            shmat   => self%heap_vars(ithr)%shmat
            pft_ref => self%heap_vars(ithr)%pft_ref
            call self%gen_shmat(ithr, shift, shmat)
            if( self%iseven(i) )then
                pft_ref = shmat * self%pfts_refs_even(:,:,iref)
            else
                pft_ref = shmat * self%pfts_refs_odd(:,:,iref)
            endif
            select case(params_glob%cc_objfun)
                case(OBJFUN_CC)
                    kw = params_glob%l_kweight
                    if( present(kweight)) kw = kweight ! overrides params_glob%l_kweight
                    if( kw )then
                        call self%gencorrs_shifted_weighted_cc(pft_ref, iptcl, iref, cc)
                    else
                        call self%gencorrs_shifted_cc(pft_ref, iptcl, iref, cc)
                    endif
                case(OBJFUN_EUCLID)
                    call self%gencorrs_shifted_euclid(pft_ref, iptcl, iref, cc)
            end select
        endif
    end subroutine gencorrs_2

    subroutine gencorrs_3( self, irefs, prefs, iptcl, cc, kweight )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: iptcl
        real(sp),                intent(out)   :: cc(self%nrots)
        logical,       optional, intent(in)    :: kweight
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                THROW_HARD('multi-refs gencorrs are not supported for cc yet!')
            case(OBJFUN_EUCLID)
                call self%gencorrs_euclid(iptcl, irefs, prefs, cc)
        end select
    end subroutine gencorrs_3

    subroutine gencorrs_4( self, irefs, prefs, iptcl, shift, cc, kweight )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: iptcl
        real(sp),                intent(in)    :: shift(2)
        real(sp),                intent(out)   :: cc(self%nrots)
        logical,       optional, intent(in)    :: kweight
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer :: i, ithr
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        shmat => self%heap_vars(ithr)%shmat
        call self%linear_ft_refs(iptcl, irefs, prefs, pft_ref, .true.)
        call self%gen_shmat(ithr, shift, shmat)
        pft_ref = shmat * pft_ref
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                THROW_HARD('multi-refs gencorrs are not supported for cc yet!')
            case(OBJFUN_EUCLID)
                call self%gencorrs_shifted_euclid(pft_ref, iptcl, irefs, prefs, cc)
        end select
    end subroutine gencorrs_4

    subroutine gencorrs_cc( self, iptcl, iref, corrs)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: corrs(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(dp) :: sqsumref
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        if( self%with_ctf )then
            self%drvec(ithr)%r = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! FT(CTF2) x FT(REF2)*)
                if( even )then
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
                else
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
                endif
                ! IFFT(FT(CTF2) x FT(REF2)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%drvec(ithr)%r(1:self%nrots) = self%drvec(ithr)%r(1:self%nrots) + real(self%rvec1(ithr)%r(1:self%nrots),dp)
                ! FT(X.CTF) x FT(REF)*
                if( even )then
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_even(:,k,iref)
                else
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_odd(:,k,iref)
                endif
                ! IFFT( FT(X.CTF) x FT(REF)* )
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs(1:self%nrots) = self%heap_vars(ithr)%kcorrs(1:self%nrots) + real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            self%drvec(ithr)%r(1:self%nrots) = self%drvec(ithr)%r(1:self%nrots) * (self%sqsums_ptcls(i) * real(2*self%nrots,dp))
            corrs = real(self%heap_vars(ithr)%kcorrs(1:self%nrots) / dsqrt(self%drvec(ithr)%r(1:self%nrots)))
        else
            pft_ref => self%heap_vars(ithr)%pft_ref
            if( even )then
                pft_ref = self%pfts_refs_even(:,:,iref)
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)
            endif
            sqsumref = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! |REF|2
                sqsumref = sqsumref + sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
                ! FT(X.CTF) x FT(REF)*
                if( even )then
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_even(:,k,iref)
                else
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_odd(:,k,iref)
                endif
                ! IFFT( FT(X.CTF) x FT(REF)* )
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs(1:self%nrots) = self%heap_vars(ithr)%kcorrs(1:self%nrots) + real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            corrs = real(self%heap_vars(ithr)%kcorrs(1:self%nrots) / (dsqrt(self%sqsums_ptcls(i)*sqsumref) * real(2*self%nrots,dp)))
        endif
    end subroutine gencorrs_cc

    subroutine gencorrs_shifted_cc( self, pft_ref, iptcl, iref, corrs)
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: corrs(self%nrots)
        real(dp) :: sqsumref
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        if( self%with_ctf )then
            self%drvec(ithr)%r = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! FT(CTF2) x FT(REF2)), REF2 is shift invariant
                if( even )then
                    self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
                else
                    self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
                endif
                ! IFFT(FT(CTF2) x FT(REF2))
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%drvec(ithr)%r = self%drvec(ithr)%r + real(self%rvec1(ithr)%r(1:self%nrots),dp)
                ! FT(S.REF), shifted reference
                self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                ! FT(X.CTF) x FT(S.REF)*
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! IFFT(FT(X.CTF) x FT(S.REF)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            self%drvec(ithr)%r = self%drvec(ithr)%r * real(self%sqsums_ptcls(i) * real(2*self%nrots),dp)
            corrs = real(self%heap_vars(ithr)%kcorrs / dsqrt(self%drvec(ithr)%r))
        else
            sqsumref = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! |REF|2
                sqsumref = sqsumref + sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
                ! FT(S.REF), shifted reference
                self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                ! FT(X) x FT(S.REF)*
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! IFFT(FT(X) x FT(S.REF)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            corrs = real(self%heap_vars(ithr)%kcorrs / (dsqrt(self%sqsums_ptcls(i)*sqsumref) * real(2*self%nrots,dp)))
        endif
    end subroutine gencorrs_shifted_cc

    subroutine gencorrs_weighted_cc( self, iptcl, iref, corrs)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: corrs(self%nrots)
        complex(sp), pointer :: pft_ref(:,:)
        real(dp) :: sqsumref
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        if( self%with_ctf )then
            self%drvec(ithr)%r = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! FT(CTF2) x FT(REF2)*)
                if( even )then
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
                else
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
                endif
                ! IFFT(FT(CTF2) x FT(REF2)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%drvec(ithr)%r(1:self%nrots) = self%drvec(ithr)%r(1:self%nrots) + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
                ! FT(X.CTF) x FT(REF)*
                if( even )then
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_even(:,k,iref)
                else
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_odd(:,k,iref)
                endif
                ! IFFT( FT(X.CTF) x FT(REF)* )
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs(1:self%nrots) = self%heap_vars(ithr)%kcorrs(1:self%nrots) + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            self%drvec(ithr)%r(1:self%nrots) = self%drvec(ithr)%r(1:self%nrots) * (self%ksqsums_ptcls(i) * real(2*self%nrots,dp))
            corrs = real(self%heap_vars(ithr)%kcorrs(1:self%nrots) / dsqrt(self%drvec(ithr)%r(1:self%nrots)))
        else
            pft_ref => self%heap_vars(ithr)%pft_ref
            if( even )then
                pft_ref = self%pfts_refs_even(:,:,iref)
            else
                pft_ref = self%pfts_refs_odd(:,:,iref)
            endif
            sqsumref = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! |REF|2
                sqsumref = sqsumref + real(k,dp) * sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
                ! FT(X) x FT(REF)*
                if( even )then
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_even(:,k,iref)
                else
                    self%cvec1(ithr)%c(1:self%pftsz+1) = self%ft_ptcl_ctf(:,k,i) * self%ft_ref_odd(:,k,iref)
                endif
                ! IFFT(FT(X) x FT(REF)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            corrs = real(self%heap_vars(ithr)%kcorrs / (dsqrt(self%ksqsums_ptcls(i)*sqsumref) * real(2*self%nrots,dp)))
        endif
    end subroutine gencorrs_weighted_cc

    subroutine gencorrs_shifted_weighted_cc( self, pft_ref, iptcl, iref, corrs)
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp),             intent(in)    :: pft_ref(1:self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: corrs(self%nrots)
        real(dp) :: sqsumref
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        if( self%with_ctf )then
            self%drvec(ithr)%r = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! FT(CTF2) x FT(REF2)), REF2 is shift invariant
                if( even )then
                    self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
                else
                    self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
                endif
                ! IFFT(FT(CTF2) x FT(REF2))
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%drvec(ithr)%r = self%drvec(ithr)%r + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
                ! FT(S.REF), shifted reference
                self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                ! FT(X.CTF) x FT(S.REF)*
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! IFFT(FT(X.CTF) x FT(S.REF)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            self%drvec(ithr)%r = self%drvec(ithr)%r * real(self%ksqsums_ptcls(i) * real(2*self%nrots),dp)
            corrs = real(self%heap_vars(ithr)%kcorrs / dsqrt(self%drvec(ithr)%r))
        else
            sqsumref = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                ! |REF|2
                sqsumref = sqsumref + real(k,dp) * sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
                ! FT(S.REF), shifted reference
                self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                ! FT(X) x FT(S.REF)*
                self%cvec1(ithr)%c = self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! IFFT(FT(X) x FT(S.REF)*)
                call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
                self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(k,dp) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            end do
            corrs = real(self%heap_vars(ithr)%kcorrs / (dsqrt(self%ksqsums_ptcls(i)*sqsumref) * real(2*self%nrots,dp)))
        endif
    end subroutine gencorrs_shifted_weighted_cc

    subroutine gencorrs_euclid_1( self, iptcl, iref, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(self%pfts_ptcls(:,k,i)*conjg(self%pfts_ptcls(:,k,i)),dp))
            ! FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)*
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i)     * self%ft_ref2_even(:,k,iref) - &
                                &2.0*self%ft_ptcl_ctf(:,k,i) * self%ft_ref_even(:,k,iref)
            else
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i)     * self%ft_ref2_odd(:,k,iref) - &
                                &2.0*self%ft_ptcl_ctf(:,k,i) * self%ft_ref_odd(:,k,iref)
            endif
            ! X.CTF.REF = IFFT( FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2.X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2.X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gencorrs_euclid_1

    subroutine gencorrs_euclid_2( self, iptcl, irefs, prefs, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iptcl
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        real(sp),                intent(out)   :: euclids(self%nrots)
        complex(sp), pointer :: ft_ref(:,:), ft_ref2(:,:)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ft_ref2 => self%heap_vars(ithr)%pft_tmp2
        call self%linear_ft_refs(iptcl, irefs, prefs, ft_ref, .false.)
        if( params_glob%l_linstates ) call self%linear_ft_ref2s(iptcl, irefs, prefs, ft_ref2)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            if( .not. params_glob%l_linstates ) call self%linear_ft_ref2s(iptcl, irefs, prefs, k, ft_ref2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(self%pfts_ptcls(:,k,i)*conjg(self%pfts_ptcls(:,k,i)),dp))
            ! FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)*
            self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * ft_ref2(:,k) - 2.0*self%ft_ptcl_ctf(:,k,i) * ft_ref(:,k)
            ! X.CTF.REF = IFFT( FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2.X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2.X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gencorrs_euclid_2

    subroutine gencorrs_shifted_euclid_1( self, pft_ref, iptcl, iref, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(in)    :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl, iref
        real(sp),                intent(out)   :: euclids(self%nrots)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        logical  :: even
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(self%pfts_ptcls(:,k,i)*conjg(self%pfts_ptcls(:,k,i)),dp))
            ! FT(CTF2) x FT(REF2)*
            if( even )then
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
            else
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
            endif
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            =       pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)*
            self%cvec1(ithr)%c = self%cvec1(ithr)%c - 2.0 * self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! IFFT( FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gencorrs_shifted_euclid_1

    subroutine gencorrs_shifted_euclid_2( self, pft_ref, iptcl, irefs, prefs, euclids )
        class(polarft_corrcalc), intent(inout) :: self
        complex(sp), pointer,    intent(in)    :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        real(sp),                intent(out)   :: euclids(self%nrots)
        complex(sp), pointer :: ft_ref2(:,:)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        self%heap_vars(ithr)%kcorrs = 0.d0
        ft_ref2 => self%heap_vars(ithr)%pft_tmp2
        if( params_glob%l_linstates ) call self%linear_ft_ref2s(iptcl, irefs, prefs, ft_ref2)
        do k = self%kfromto(1),self%kfromto(2)
            if( .not. params_glob%l_linstates ) call self%linear_ft_ref2s(iptcl, irefs, prefs, k, ft_ref2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(self%pfts_ptcls(:,k,i)*conjg(self%pfts_ptcls(:,k,i)),dp))
            ! FT(CTF2) x FT(REF2)*
            self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * ft_ref2(:,k)
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            =       pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)*
            self%cvec1(ithr)%c = self%cvec1(ithr)%c - 2.0 * self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! IFFT( FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gencorrs_shifted_euclid_2

    subroutine bidirectional_shift_search( self, iref, iptcl, irot, hn, shifts, grid1, grid2 )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot, hn
        real,                    intent(in)    :: shifts(-hn:hn)
        real,                    intent(out)   :: grid1(-hn:hn,-hn:hn), grid2(-hn:hn,-hn:hn)
        complex(dp), pointer :: pft_ref_8(:,:), diff_8(:,:), pft_shref_8(:,:)
        complex(sp), pointer :: pft_ptcl(:,:)
        real(sp),    pointer :: rctf(:,:)
        real(dp)             :: shvec(2), score, sqsum_ref, denom
        integer              :: ithr, i, ix, iy, k, prot
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ptcl    => self%heap_vars(ithr)%pft_ref
        pft_ref_8   => self%heap_vars(ithr)%pft_ref_8
        diff_8      => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_shref_8 => self%heap_vars(ithr)%shmat_8
        rctf        => self%heap_vars(ithr)%pft_r
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! Rotate particle
        prot = self%nrots-irot+2
        if( prot > self%nrots ) prot = prot-self%nrots
        call self%rotate_pft(self%pfts_ptcls(:,:,i), prot, pft_ptcl(:,:))
        if( self%with_ctf )then
            ! Reference CTF modulation
            call self%rotate_pft(self%ctfmats(:,:,i), prot, rctf)
            pft_ref_8 = pft_ref_8 * real(rctf,dp)
        endif
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            do ix = -hn,hn
                do iy = -hn,hn
                    shvec = real([shifts(ix), shifts(iy)],dp)
                    call self%gen_shmat_8(ithr, shvec, pft_shref_8) ! shift matrix
                    pft_shref_8 = pft_ref_8 * pft_shref_8           ! shifted reference
                    ! first orientation
                    sqsum_ref = 0.d0
                    score     = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_shref_8(:,k) * conjg(pft_shref_8(:,k)),dp))
                        score     = score     + real(k,kind=dp) * sum(real(pft_shref_8(:,k) * conjg(pft_ptcl(:,k)),dp))
                    end do
                    denom        = dsqrt(sqsum_ref * self%ksqsums_ptcls(i))
                    grid1(ix,iy) = real(score / denom)
                    ! second orientation (first+pi)
                    score = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        score = score + real(k,kind=dp) * sum(real(pft_shref_8(:,k) * pft_ptcl(:,k),dp))
                    end do
                    grid2(ix,iy) = real(score / denom)
                enddo
            enddo
        case(OBJFUN_EUCLID)
            do ix = -hn,hn
                do iy = -hn,hn
                    shvec = real([shifts(ix), shifts(iy)],dp)
                    call self%gen_shmat_8(ithr, shvec, pft_shref_8) ! shift matrix
                    pft_shref_8 = pft_ref_8 * pft_shref_8           ! shifted reference
                    ! first orientation
                    diff_8 = pft_shref_8 - dcmplx(pft_ptcl)
                    score = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        score = score + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(diff_8(:,k)*conjg(diff_8(:,k)),dp))
                    end do
                    grid1(ix,iy) = real(exp( -score / self%wsqsums_ptcls(i) ))
                    ! second orientation (first+pi)
                    diff_8 = pft_shref_8 - dcmplx(conjg(pft_ptcl))
                    score = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        score = score + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(diff_8(:,k)*conjg(diff_8(:,k)),dp))
                    end do
                    grid2(ix,iy) = real(exp( -score / self%wsqsums_ptcls(i) ))
                enddo
            enddo
        end select
    end subroutine bidirectional_shift_search

    real(dp) function gencorr_for_rot_8_1( self, iref, iptcl, irot, onestate )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl, irot
        logical, optional,       intent(in)    :: onestate
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:)
        integer :: ithr, i, istate, irefs(params_glob%nstates)
        logical :: l_onestate
        real    :: prefs(params_glob%nstates)
        l_onestate = .false.
        if( present(onestate) ) l_onestate = onestate
        if( params_glob%l_linstates .and. .not.(l_onestate) )then
            call self%get_linstates_irefs(iref, irefs)
            do istate = 1, params_glob%nstates
                prefs(istate) = real(self%gencorr_for_rot_8_1(irefs(istate), iptcl, irot, onestate=.true.))
            enddo
            call self%get_linstates_prefs(iref, irefs, prefs)
            gencorr_for_rot_8_1 = self%gencorr_for_rot_8_3( irefs, prefs, iptcl, irot )
        else
            i    =  self%pinds(iptcl)
            ithr = omp_get_thread_num() + 1
            pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
            pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
            if( self%iseven(i) )then
                pft_ref_8 = self%pfts_refs_even(:,:,iref)
            else
                pft_ref_8 = self%pfts_refs_odd(:,:,iref)
            endif
            ! rotation
            call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
            ! ctf
            if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
            gencorr_for_rot_8_1 = 0.d0
            select case(params_glob%cc_objfun)
                case(OBJFUN_CC)
                    gencorr_for_rot_8_1 = self%gencorr_cc_for_rot_8(pft_ref_tmp_8, i)
                case(OBJFUN_EUCLID)
                    gencorr_for_rot_8_1 = self%gencorr_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
            end select
        endif
    end function gencorr_for_rot_8_1

    real(dp) function gencorr_for_rot_8_2( self, iref, iptcl, shvec, irot, onestate )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        logical, optional,       intent(in)    :: onestate
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        integer :: ithr, i, istate, irefs(params_glob%nstates)
        logical :: l_onestate
        real    :: prefs(params_glob%nstates)
        l_onestate = .false.
        if( present(onestate) ) l_onestate = onestate
        if( params_glob%l_linstates .and. .not.(l_onestate) )then
            call self%get_linstates_irefs(iref, irefs)
            do istate = 1, params_glob%nstates
                prefs(istate) = real(self%gencorr_for_rot_8_2(irefs(istate), iptcl, shvec, irot, onestate=.true.))
            enddo
            call self%get_linstates_prefs(iref, irefs, prefs)
            gencorr_for_rot_8_2 = self%gencorr_for_rot_8_4(irefs, prefs, iptcl, shvec, irot )
        else
            i    =  self%pinds(iptcl)
            ithr = omp_get_thread_num() + 1
            pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
            pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
            shmat_8       => self%heap_vars(ithr)%shmat_8
            if( self%iseven(i) )then
                pft_ref_8 = self%pfts_refs_even(:,:,iref)
            else
                pft_ref_8 = self%pfts_refs_odd(:,:,iref)
            endif
            ! shift
            call self%gen_shmat_8(ithr, shvec, shmat_8)
            pft_ref_8 = pft_ref_8 * shmat_8
            ! rotation
            call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
            ! ctf
            if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
            gencorr_for_rot_8_2 = 0.d0
            select case(params_glob%cc_objfun)
                case(OBJFUN_CC)
                    gencorr_for_rot_8_2 = self%gencorr_cc_for_rot_8(pft_ref_tmp_8, i)
                case(OBJFUN_EUCLID)
                    gencorr_for_rot_8_2 = self%gencorr_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
            end select
        endif
    end function gencorr_for_rot_8_2

    real(dp) function gencorr_for_rot_8_3( self, irefs, prefs, iptcl, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: iptcl, irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:)
        complex(sp), pointer :: pft_ref(:,:)
        integer              :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref       => self%heap_vars(ithr)%pft_ref
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        call self%linear_ft_refs(iptcl, irefs, prefs, pft_ref, pft_ref_or_ft_ref=.true.)
        pft_ref_8 = dcmplx(pft_ref)
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        gencorr_for_rot_8_3 = 0.d0
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gencorr_for_rot_8_3 = self%gencorr_cc_for_rot_8(pft_ref_tmp_8, i)
            case(OBJFUN_EUCLID)
                gencorr_for_rot_8_3 = self%gencorr_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
        end select
    end function gencorr_for_rot_8_3

    real(dp) function gencorr_for_rot_8_4( self, irefs, prefs, iptcl, shvec, irot )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        complex(sp), pointer :: pft_ref(:,:)
        integer              :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref       => self%heap_vars(ithr)%pft_ref
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        call self%linear_ft_refs(iptcl, irefs, prefs, pft_ref, pft_ref_or_ft_ref=.true.)
        pft_ref_8 = dcmplx(pft_ref)
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        gencorr_for_rot_8_4 = 0.d0
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gencorr_for_rot_8_4 = self%gencorr_cc_for_rot_8(pft_ref_tmp_8, i)
            case(OBJFUN_EUCLID)
                gencorr_for_rot_8_4 = self%gencorr_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
        end select
    end function gencorr_for_rot_8_4

    real(dp) function gencorr_cc_for_rot_8( self, pft_ref, i )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: i
        real(dp) :: sqsum_ref
        integer  :: k
        sqsum_ref            = 0.d0
        gencorr_cc_for_rot_8 = 0.d0
        if( params_glob%l_kweight_shift )then
            do k = self%kfromto(1),self%kfromto(2)
                sqsum_ref            = sqsum_ref +            real(k,kind=dp) * sum(real(pft_ref(:,k) * conjg(pft_ref(:,k)),dp))
                gencorr_cc_for_rot_8 = gencorr_cc_for_rot_8 + real(k,kind=dp) * sum(real(pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            end do
            gencorr_cc_for_rot_8 = gencorr_cc_for_rot_8 / dsqrt(sqsum_ref * self%ksqsums_ptcls(i))
        else
            do k = self%kfromto(1),self%kfromto(2)
                sqsum_ref            = sqsum_ref +            sum(real(pft_ref(:,k) * conjg(pft_ref(:,k)),dp))
                gencorr_cc_for_rot_8 = gencorr_cc_for_rot_8 + sum(real(pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            end do
            gencorr_cc_for_rot_8 = gencorr_cc_for_rot_8 / dsqrt(sqsum_ref * self%sqsums_ptcls(i))
        endif
    end function gencorr_cc_for_rot_8

    real(dp) function gencorr_euclid_for_rot_8( self, pft_ref, iptcl )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:)
        integer,                 intent(in)    :: iptcl
        integer :: i,k
        i       = self%pinds(iptcl)
        pft_ref = pft_ref - self%pfts_ptcls(:,:,i)
        gencorr_euclid_for_rot_8 = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            gencorr_euclid_for_rot_8 = gencorr_euclid_for_rot_8 +&
                &(real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
        end do
        gencorr_euclid_for_rot_8 = dexp( -gencorr_euclid_for_rot_8 / self%wsqsums_ptcls(i))
    end function gencorr_euclid_for_rot_8

    subroutine gencorr_grad_for_rot_8_1( self, iref, iptcl, shvec, irot, f, grad, onestate )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        logical, optional,       intent(in)    :: onestate
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer :: ithr, i, istate, irefs(params_glob%nstates)
        logical :: l_onestate
        real    :: prefs(params_glob%nstates)
        l_onestate = .false.
        if( present(onestate) ) l_onestate = onestate
        if( params_glob%l_linstates .and. .not.(l_onestate) )then
            call self%get_linstates_irefs(iref, irefs)
            do istate = 1, params_glob%nstates
                prefs(istate) = real(self%gencorr_for_rot_8_2(irefs(istate), iptcl, shvec, irot, onestate=.true.))
            enddo
            call self%get_linstates_prefs(iref, irefs, prefs)
            call self%gencorr_grad_for_rot_8_2(irefs, prefs, iptcl, shvec, irot, f, grad)
        else
            i    = self%pinds(iptcl)
            ithr = omp_get_thread_num() + 1
            pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
            pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
            shmat_8       => self%heap_vars(ithr)%shmat_8
            if( self%iseven(i) )then
                pft_ref_8 = self%pfts_refs_even(:,:,iref)
            else
                pft_ref_8 = self%pfts_refs_odd(:,:,iref)
            endif
            call self%gen_shmat_8(ithr, shvec, shmat_8)
            pft_ref_8 = pft_ref_8 * shmat_8
            select case(params_glob%cc_objfun)
                case(OBJFUN_CC)
                    call self%gencorr_cc_grad_for_rot_8(    pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
                case(OBJFUN_EUCLID)
                    call self%gencorr_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            end select
        endif
    end subroutine gencorr_grad_for_rot_8_1

    subroutine gencorr_grad_for_rot_8_2( self, irefs, prefs, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        complex(sp), pointer :: pft_ref(:,:)
        integer              :: ithr, i
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref       => self%heap_vars(ithr)%pft_ref
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        call self%linear_ft_refs(iptcl, irefs, prefs, pft_ref, pft_ref_or_ft_ref=.true.)
        pft_ref_8 = dcmplx(pft_ref)
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_for_rot_8(    pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gencorr_grad_for_rot_8_2

    subroutine gencorr_cc_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        real(dp) :: sqsum_ref, sqsum_ptcl, denom
        integer  :: k, i
        i         = self%pinds(iptcl)
        sqsum_ref = 0.d0
        f         = 0.d0
        grad      = 0.d0
        if( self%with_ctf )then
            if( params_glob%l_kweight_shift )then
                sqsum_ptcl = self%ksqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                    f         = f         + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)                     * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            else
                sqsum_ptcl = self%sqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                    f         = f         + sum(real(self%ctfmats(:,k,i)                     * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            endif
        else
            if( params_glob%l_kweight_shift )then
                sqsum_ptcl = self%ksqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                    f         = f         + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            else
                sqsum_ptcl = self%sqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + sum(real(pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                    f         = f         + sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            endif
        endif
        denom = dsqrt(sqsum_ref * sqsum_ptcl)
        f     = f    / denom
        grad  = grad / denom
    end subroutine gencorr_cc_grad_for_rot_8

    subroutine gencorr_euclid_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: iptcl, irot
        real(dp),                intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_diff(:,:)
        real(dp) :: denom, w
        integer  :: k, i, ithr
        ithr     = omp_get_thread_num() + 1
        i        = self%pinds(iptcl)
        f        = 0.d0
        grad     = 0.d0
        denom    = self%wsqsums_ptcls(i)
        pft_diff => self%heap_vars(ithr)%shmat_8
        call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        pft_diff = pft_ref_tmp - self%pfts_ptcls(:,:,i) ! Ref(shift + rotation + CTF) - Ptcl
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w       = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            f       = f + w * sum(real(pft_diff(:,k)*conjg(pft_diff(:,k)),dp))
            grad(1) = grad(1) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
        end do
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
        if( self%with_ctf ) pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w       = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            grad(2) = grad(2) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
        end do
        f    = dexp( -f / denom )
        grad = -f * 2.d0 * grad / denom
    end subroutine gencorr_euclid_grad_for_rot_8

    subroutine gencorr_grad_only_for_rot_8_1( self, iref, iptcl, shvec, irot, grad, onestate )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        logical, optional,       intent(in)    :: onestate
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        real(dp) :: f
        integer  :: ithr, i, istate, irefs(params_glob%nstates)
        logical  :: l_onestate
        real     :: prefs(params_glob%nstates)
        l_onestate = .false.
        if( present(onestate) ) l_onestate = onestate
        if( params_glob%l_linstates .and. .not.(l_onestate) )then
            call self%get_linstates_irefs(iref, irefs)
            do istate = 1, params_glob%nstates
                prefs(istate) = real(self%gencorr_for_rot_8_2(irefs(istate), iptcl, shvec, irot, onestate=.true.))
            enddo
            call self%get_linstates_prefs(iref, irefs, prefs)
            call self%gencorr_grad_only_for_rot_8_2(irefs, prefs, iptcl, shvec, irot, grad)
        else
            i    = self%pinds(iptcl)
            ithr = omp_get_thread_num() + 1
            pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
            pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
            shmat_8       => self%heap_vars(ithr)%shmat_8
            if( self%iseven(i) )then
                pft_ref_8 = self%pfts_refs_even(:,:,iref)
            else
                pft_ref_8 = self%pfts_refs_odd(:,:,iref)
            endif
            call self%gen_shmat_8(ithr, shvec, shmat_8)
            pft_ref_8 = pft_ref_8 * shmat_8
            select case(params_glob%cc_objfun)
                case(OBJFUN_CC)
                    call self%gencorr_cc_grad_only_for_rot_8(pft_ref_8, pft_ref_tmp_8, i, irot, grad)
                case(OBJFUN_EUCLID)
                    call self%gencorr_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            end select
        endif
    end subroutine gencorr_grad_only_for_rot_8_1

    subroutine gencorr_grad_only_for_rot_8_2( self, irefs, prefs, iptcl, shvec, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: irefs(:)
        real,                    intent(in)    :: prefs(:)
        integer,                 intent(in)    :: iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        complex(sp), pointer :: pft_ref(:,:)
        real(dp) :: f
        integer  :: ithr, i
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref       => self%heap_vars(ithr)%pft_ref
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        call self%linear_ft_refs(iptcl, irefs, prefs, pft_ref, pft_ref_or_ft_ref=.true.)
        pft_ref_8 = dcmplx(pft_ref)
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gencorr_cc_grad_only_for_rot_8(pft_ref_8, pft_ref_tmp_8, i, irot, grad)
            case(OBJFUN_EUCLID)
                call self%gencorr_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gencorr_grad_only_for_rot_8_2

    subroutine gencorr_cc_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, i, irot, grad )
        class(polarft_corrcalc), intent(inout) :: self
        complex(dp), pointer,    intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                 intent(in)    :: i, irot
        real(dp),                intent(out)   :: grad(2)
        real(dp) :: sqsum_ref, sqsum_ptcl
        integer  :: k
        sqsum_ref = 0.d0
        grad      = 0.d0
        if( self%with_ctf )then
            if( params_glob%l_kweight_shift )then
                sqsum_ptcl = self%ksqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            else
                sqsum_ptcl = self%sqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            endif
        else
            if( params_glob%l_kweight_shift )then
                sqsum_ptcl = self%ksqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + real(k,kind=dp) * sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            else
                sqsum_ptcl = self%sqsums_ptcls(i)
                call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    sqsum_ref = sqsum_ref + sum(real(pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(1) = grad(1) + sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                enddo
                call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
                do k = self%kfromto(1),self%kfromto(2)
                    grad(2) = grad(2) + sum(real(pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
                end do
            endif
        endif
        grad = grad / dsqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gencorr_cc_grad_only_for_rot_8

    function gencorr_cont_cc_for_rot_8( self, iref, iptcl, shvec, irot )result(cc)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        real(dp)             :: sqsum_ref, cc
        integer              :: ithr, i, k
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        ! correlation
        sqsum_ref = 0.0
        cc        = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(pft_ref_tmp_8(:,k)),dp))
            cc        = cc        + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        end do
        cc = cc / dsqrt(sqsum_ref * self%ksqsums_ptcls(i))
    end function gencorr_cont_cc_for_rot_8

    function gencorr_cont_grad_cc_for_rot_8( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:), pft_dref_8(:,:,:)
        real(dp) :: T1(3), T2(3), cc, sqsum_ref, denom, num
        integer  :: ithr, j, k, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_dref_8    => self%heap_vars(ithr)%pft_dref_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8  = self%pfts_refs_even(:,:,iref)
            pft_dref_8 = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref_8  = self%pfts_refs_odd(:,:,iref)
            pft_dref_8 = self%pfts_drefs_odd(:,:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        do j = 1,3
            pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * shmat_8
        end do
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        pft_ref_8 = pft_ref_tmp_8
        do j = 1,3
            call self%rotate_pft(pft_dref_8(:,:,j), irot, pft_ref_tmp_8)
            pft_dref_8(:,:,j) = pft_ref_tmp_8
        end do
        ! ctf
        if( self%with_ctf )then
            pft_ref_8 = pft_ref_8 * self%ctfmats(:,:,i)
            do j = 1,3
                pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * self%ctfmats(:,:,i)
            end do
        endif
        ! correlation & derivatives
        sqsum_ref = 0.d0
        num       = 0.d0
        T1        = 0.d0
        T2        = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(pft_ref_8(:,k)),dp))
            num        = num       + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            do j = 1,3
                T1(j) = T1(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(self%pfts_ptcls(:,k,i))),dp)
                T2(j) = T2(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(pft_ref_8(:,k))),dp)
            enddo
        enddo
        denom = sqrt(sqsum_ref * self%ksqsums_ptcls(i))
        cc    = num / denom
        dcc   = (T1 - num * T2 / sqsum_ref) / denom
    end function gencorr_cont_grad_cc_for_rot_8

    subroutine gencorr_cont_shift_grad_cc_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:), pft_dref_8(:,:,:)
        real(dp) :: T1(3), T2(3), sqsum_ref, denom, num
        integer  :: ithr, j, k, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_dref_8    => self%heap_vars(ithr)%pft_dref_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8  = self%pfts_refs_even(:,:,iref)
            pft_dref_8 = self%pfts_drefs_even(:,:,:,iref)
        else
            pft_ref_8  = self%pfts_refs_odd(:,:,iref)
            pft_dref_8 = self%pfts_drefs_odd(:,:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        do j = 1,3
            pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * shmat_8
        end do
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        pft_ref_8 = pft_ref_tmp_8
        do j = 1,3
            call self%rotate_pft(pft_dref_8(:,:,j), irot, pft_ref_tmp_8)
            pft_dref_8(:,:,j) = pft_ref_tmp_8
        end do
        ! ctf
        if( self%with_ctf )then
            pft_ref_8 = pft_ref_8 * self%ctfmats(:,:,i)
            do j = 1,3
                pft_dref_8(:,:,j) = pft_dref_8(:,:,j) * self%ctfmats(:,:,i)
            end do
        endif
        ! correlation & orientation derivatives
        sqsum_ref = 0.d0
        num       = 0.d0
        T1        = 0.d0
        T2        = 0.d0
        grad      = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref  = sqsum_ref + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(pft_ref_8(:,k)),dp))
            num        = num       + real(k,kind=dp) * sum(real(pft_ref_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
            do j = 1,3
                T1(j) = T1(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(self%pfts_ptcls(:,k,i))),dp)
                T2(j) = T2(j) + real(k,kind=dp) * real(sum(pft_dref_8(:,k,j) * conjg(pft_ref_8(:,k))),dp)
            enddo
        enddo
        denom     = sqrt(sqsum_ref * self%ksqsums_ptcls(i))
        f         = num / denom
        grad(1:3) = (T1 - num * T2 / sqsum_ref) / denom
        ! shift derivatives
        pft_ref_tmp_8 = pft_ref_8 * (0.d0, 1.d0) * self%argtransf(:self%pftsz,:)
        num = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            num = num + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        grad(4) = num / denom
        pft_ref_tmp_8 = pft_ref_8 * (0.d0, 1.d0) * self%argtransf(self%pftsz+1:,:)
        num = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            num = num + real(k,kind=dp) * sum(real(pft_ref_tmp_8(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        grad(5) = num / denom
    end subroutine gencorr_cont_shift_grad_cc_for_rot_8

    function gencorr_cont_grad_euclid_for_rot_8( self, iref, iptcl, shvec, irot, dcc ) result( cc )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: dcc(3)
        real(dp) :: cc
        cc = 0._dp
        ! TODO: implement me
    end function gencorr_cont_grad_euclid_for_rot_8

    subroutine gencorr_cont_shift_grad_euclid_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(dp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(dp),                intent(out)   :: f
        real(dp),                intent(out)   :: grad(5) ! 3 orientation angles, 2 shifts
        ! TODO: implement me
    end subroutine gencorr_cont_shift_grad_euclid_for_rot_8

    subroutine gencorr_sigma_contrib( self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_corrcalc), intent(inout) :: self
        integer,                 intent(in)    :: iref, iptcl
        real(sp),                intent(in)    :: shvec(2)
        integer,                 intent(in)    :: irot
        real(sp),      optional, intent(out)   :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer  :: i, ithr
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        if( self%with_ctf ) pft_ref_tmp_8 = pft_ref_tmp_8 * real(self%ctfmats(:,:,i),dp)
        ! difference
        pft_ref_tmp_8 = pft_ref_tmp_8 - self%pfts_ptcls(:,:,i)
        ! sigma2
        if( present(sigma_contrib) )then
            sigma_contrib = real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8),dp), dim=1) / (2.d0*real(self%pftsz,dp)))
        else
            self%sigma2_noise(self%kfromto(1):self%kfromto(2), iptcl) = real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8),dp), dim=1) / (2.d0*real(self%pftsz,dp)))
        endif
    end subroutine gencorr_sigma_contrib

    ! DESTRUCTOR

    subroutine kill( self )
        class(polarft_corrcalc), intent(inout) :: self
        integer :: ithr
        if( self%existence )then
            do ithr=1,params_glob%nthr
                deallocate(self%heap_vars(ithr)%pft_ref,self%heap_vars(ithr)%pft_ref_tmp,&
                    &self%heap_vars(ithr)%argvec, self%heap_vars(ithr)%shvec,&
                    &self%heap_vars(ithr)%shmat,self%heap_vars(ithr)%kcorrs,&
                    &self%heap_vars(ithr)%pft_tmp,self%heap_vars(ithr)%pft_tmp2,&
                    &self%heap_vars(ithr)%pft_ref_8,self%heap_vars(ithr)%pft_ref_tmp_8,&
                    &self%heap_vars(ithr)%pft_ptcl_8,&
                    &self%heap_vars(ithr)%pft_dref_8,self%heap_vars(ithr)%pft_r,&
                    &self%heap_vars(ithr)%shmat_8,self%heap_vars(ithr)%pft_r1_8,self%heap_vars(ithr)%pft_r2_8)
            end do
            if( allocated(self%ctfmats)        ) deallocate(self%ctfmats)
            if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
            if( allocated(self%has_cache)      ) deallocate(self%has_cache)
            if( allocated(self%do_cache)       ) deallocate(self%do_cache)
            if( allocated(self%cached_vals)    ) deallocate(self%cached_vals)
            deallocate(self%sqsums_ptcls, self%ksqsums_ptcls, self%wsqsums_ptcls, self%angtab, self%argtransf,self%pfts_ptcls,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_drefs_even, self%pfts_drefs_odd,&
                &self%iseven, self%pinds, self%heap_vars, self%argtransf_shellone)
            call self%kill_memoized_ptcls
            call self%kill_memoized_refs
            nullify(self%sigma2_noise, pftcc_glob)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_polarft_corrcalc

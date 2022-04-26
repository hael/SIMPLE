module simple_cartft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none

public :: cartft_corrcalc, cartftcc_glob
private
#include "simple_local_flags.inc"

type :: cartft_corrcalc
    private
    integer                  :: nptcls        = 1       !< the total number of particles in partition (logically indexded [fromp,top])
    integer                  :: nrefs         = 1       !< the number of references (logically indexded [1,nrefs])
    integer                  :: lnyq          = 0       !< Nyqvist limit
    integer                  :: pfromto(2)    = 0       !< particle index range
    integer                  :: ldim(3)       = 0       !< logical dimensions of original cartesian image
    integer                  :: lims(3,2)     = 0       !< resolution mask limits
    integer                  :: cmat_shape(3) = 0       !< shape of complex matrix (dictated by the FFTW library)
    integer,     allocatable :: pinds(:)                !< index array (to reduce memory when frac_update < 1)
    real,        allocatable :: pxls_p_shell(:)         !< number of (cartesian) pixels per shell
    real(sp),    allocatable :: sqsums_ptcls(:)         !< memoized square sums for the correlation calculations (taken from kfromto(1):kstop)
    real(sp),    allocatable :: ctfmats(:,:,:)          !< expand set of CTF matrices (for efficient parallel exec)
    real(sp),    allocatable :: ref_optlp(:,:)          !< references optimal filter
    type(image), allocatable :: refs_eo(:,:)            !< reference images even/odd, dim1=1:nrefs, dim2=1:2 (1 is even 2 is odd)
    type(image), allocatable :: particles(:)            !< particle images
    logical,     allocatable :: iseven(:)               !< e/o assignment for gold-standard FSC
    logical,     allocatable :: resmsk(:,:,:)           !< resolution mask for corr calc
    logical                  :: l_match_filt = .false.  !< matched filter flag
    logical                  :: with_ctf     = .false.  !< CTF flag
    logical                  :: existence    = .false.  !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new

    ! SETTERS

    ! GETTERS

    ! MODIFIERS

    ! MEMOIZERS
    procedure, private :: memoize_resmsk
    procedure, private :: setup_pxls_p_shell
    procedure, private :: memoize_ptcl_sqsum
    
    

    ! CALCULATORS

    ! DESTRUCTOR
    procedure          :: kill
end type cartft_corrcalc

! CLASS PARAMETERS/VARIABLES
class(cartft_corrcalc), pointer :: cartftcc_glob => null()

contains

    ! CONSTRUCTOR

    subroutine new( self, nrefs, pfromto, l_match_filt, ptcl_mask, eoarr )
        class(cartft_corrcalc), target, intent(inout) :: self
        integer,                        intent(in)    :: nrefs
        integer,                        intent(in)    :: pfromto(2)
        logical,                        intent(in)    :: l_match_filt
        logical, optional,              intent(in)    :: ptcl_mask(pfromto(1):pfromto(2))
        integer, optional,              intent(in)    :: eoarr(pfromto(1):pfromto(2))
        logical :: even_dims, test(2)
        integer :: i, cnt
        ! kill possibly pre-existing object
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
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
        self%l_match_filt = l_match_filt !< do shellnorm and filtering here (needs to be local because in 3D we do it on the reference volumes)
        if( present(ptcl_mask) )then
            self%nptcls  = count(ptcl_mask)                      !< the total number of particles in partition
        else
            self%nptcls  = self%pfromto(2) - self%pfromto(1) + 1 !< the total number of particles in partition
        endif
        self%nrefs = nrefs               !< the number of references (logically indexded [1,nrefs])
        ! allocate optimal low-pass filter & memoized sqsums
        allocate(self%ref_optlp(params_glob%kfromto(1):params_glob%kstop,self%nrefs), self%sqsums_ptcls(1:self%nptcls), source=1.)
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0 )
        if( present(ptcl_mask) )then
            cnt = 0
            do i=self%pfromto(1),self%pfromto(2)
                if( ptcl_mask(i) )then
                    cnt = cnt + 1
                    self%pinds(i) = cnt
                endif
            end do
        else
            self%pinds = (/(i,i=1,self%nptcls)/)
        endif
        ! eo assignment
        allocate( self%iseven(1:self%nptcls), source=.true. )
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
        endif
        ! allocate the rest (we worry about heap variables similar to the polarft_corrcalc class when we get there)
        allocate( self%refs_eo(self%nrefs,2), self%particles(self%nptcls) )
        do i = 1,self%nrefs
            call self%refs_eo(i,1)%new(self%ldim, params_glob%smpd)
            call self%refs_eo(i,2)%new(self%ldim, params_glob%smpd)
        end do
        do i = 1,self%nptcls
            call self%particles(i)%new(self%ldim, params_glob%smpd)
        end do
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! setup pxls_p_shell
        call self%setup_pxls_p_shell
        ! memoize resolution mask
        call self%memoize_resmsk
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        cartftcc_glob => self
    end subroutine new

    ! MEMOIZERS

    subroutine memoize_resmsk( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: h, k, l, sh, phys(3), cmat_shape(3)
        self%lims = self%particles(1)%loop_lims(2)
        if( allocated(self%resmsk) ) deallocate(self%resmsk)
        cmat_shape = self%particles(1)%get_array_shape()
        allocate(self%resmsk(cmat_shape(1),cmat_shape(2),cmat_shape(3)), source=.false.)
        do k=self%lims(2,1),self%lims(2,2)
            do h=self%lims(1,1),self%lims(1,2)
                do l=self%lims(3,1),self%lims(3,2)
                    ! compute physical address
                    phys = self%particles(1)%comp_addr_phys(h,k,l)
                    ! find shell
                    sh = nint(hyp(real(h),real(k),real(l)))
                    if( sh < params_glob%kfromto(1) .or. sh > params_glob%kstop ) cycle
                    ! update logical mask
                    self%resmsk(phys(1),phys(2),phys(3)) = .true.
                enddo
            enddo
        enddo
    end subroutine memoize_resmsk

    subroutine setup_pxls_p_shell( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: h,k,sh
        if( allocated(self%pxls_p_shell) ) deallocate(self%pxls_p_shell)
        allocate(self%pxls_p_shell(params_glob%kfromto(1):params_glob%kfromto(2)))
        self%pxls_p_shell = 0.
        do h = 1,params_glob%kfromto(2)
            do k = 1,params_glob%kfromto(2)
                sh = nint(hyp(real(h),real(k)))
                if( ( sh >= params_glob%kfromto(1)) .and. ( sh <= params_glob%kfromto(2)) ) then
                    self%pxls_p_shell(sh) = self%pxls_p_shell(sh) + 1.
                end if
            end do
        end do
    end subroutine setup_pxls_p_shell

    subroutine memoize_ptcl_sqsum( self, i )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: i
        self%sqsums_ptcls(i) = self%particles(i)%calc_sumsq(self%resmsk)
    end subroutine memoize_ptcl_sqsum

    ! DESTRUCTOR

    subroutine kill( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: i
        if( self%existence )then
            if( allocated(self%pinds)        ) deallocate(self%pinds)
            if( allocated(self%pxls_p_shell) ) deallocate(self%pxls_p_shell)
            if( allocated(self%sqsums_ptcls) ) deallocate(self%sqsums_ptcls)
            if( allocated(self%ctfmats)      ) deallocate(self%ctfmats)
            if( allocated(self%ref_optlp)    ) deallocate(self%ref_optlp)
            if( allocated(self%iseven)       ) deallocate(self%iseven)
            do i = 1,self%nrefs
                call self%refs_eo(i,1)%kill
                call self%refs_eo(i,2)%kill
            end do
            do i = 1,self%nptcls
                call self%particles(i)%kill
            end do
            deallocate(self%refs_eo, self%particles)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_cartft_corrcalc
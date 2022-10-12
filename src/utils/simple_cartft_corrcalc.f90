module simple_cartft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
use simple_projector,  only: projector
use simple_ori,        only: ori
implicit none

public :: cartft_corrcalc, cartftcc_glob
private
#include "simple_local_flags.inc"

type heap_vars
    type(image)   :: img_ref
    type(image)   :: img_ref_tmp
    real, pointer :: frc(:) => null()
end type heap_vars

type :: cartft_corrcalc
    private
    type(projector), pointer     :: vol_even => null(), vol_odd => null() ! prepared e/o vols
    integer                      :: nptcls        = 1       !< the total number of particles in partition (logically indexded [fromp,top])
    integer                      :: filtsz        = 0       !< Nyqvist limit
    integer                      :: pfromto(2)    = 0       !< particle index range
    integer                      :: ldim(3)       = 0       !< logical dimensions of original cartesian image
    integer                      :: lims(3,2)     = 0       !< resolution mask limits
    integer                      :: cmat_shape(3) = 0       !< shape of complex matrix (dictated by the FFTW library)
    integer,         allocatable :: pinds(:)                !< index array (to reduce memory when frac_update < 1)
    real,            allocatable :: pxls_p_shell(:)         !< number of (cartesian) pixels per shell
    real(sp),        allocatable :: sqsums_ptcls(:)         !< memoized square sums for the correlation calculations (taken from kfromto(1):kstop)
    real(sp),        allocatable :: ctfmats(:,:,:,:)        !< expand set of CTF matrices (for efficient parallel exec)
    real(sp),        allocatable :: optlp(:)                !< references optimal filter
    type(image),     allocatable :: particles(:)            !< particle images
    logical,         allocatable :: iseven(:)               !< e/o assignment for gold-standard FSC
    logical,         allocatable :: resmsk(:,:,:)           !< resolution mask for corr calc
    logical                      :: l_clsfrcs    = .false.  !< CLS2D/3DRefs flag
    logical                      :: l_match_filt = .false.  !< matched filter flag
    logical                      :: l_filt_set   = .false.  !< to indicate whether filter is set
    logical                      :: with_ctf     = .false.  !< CTF flag
    logical                      :: existence    = .false.  !< to indicate existence
    type(heap_vars), allocatable :: heap_vars(:)            !< allocated fields to save stack allocation in subroutines and functions
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: reallocate_ptcls
    procedure          :: set_ptcl
    procedure          :: set_optlp
    procedure          :: set_eo
    ! GETTERS
    procedure          :: get_box
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: assign_pinds
    ! MODIFIERS
    procedure, private :: shellnorm_and_filter_ref
    ! MEMOIZERS
    procedure, private :: memoize_resmsk
    procedure, private :: setup_pxls_p_shell
    procedure, private :: memoize_sqsum_ptcl
    ! CALCULATORS
    procedure          :: create_absctfmats
    procedure          :: calc_corr
    procedure          :: project_and_correlate
    ! DESTRUCTOR
    procedure          :: kill
end type cartft_corrcalc

! CLASS PARAMETERS/VARIABLES
class(cartft_corrcalc), pointer :: cartftcc_glob => null()

contains

    ! CONSTRUCTOR

    subroutine new( self, vol_even, vol_odd, pfromto, l_match_filt, ptcl_mask, eoarr )
        class(cartft_corrcalc), target, intent(inout) :: self
        type(projector),        target, intent(in)    :: vol_even, vol_odd
        integer,                        intent(in)    :: pfromto(2)
        logical,                        intent(in)    :: l_match_filt
        logical, optional,              intent(in)    :: ptcl_mask(pfromto(1):pfromto(2))
        integer, optional,              intent(in)    :: eoarr(pfromto(1):pfromto(2))
        logical :: even_dims, test(2)
        integer :: i, cnt, ithr
        ! kill possibly pre-existing object
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
        ! error check
        if( self%pfromto(2) - self%pfromto(1) + 1 < 1 )then
            write(logfhandle,*) 'pfromto: ', self%pfromto(1), self%pfromto(2)
            THROW_HARD ('nptcls (# of particles) must be > 0; new')
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
        self%filtsz = fdim(params_glob%box) - 1
        self%l_match_filt = l_match_filt !< do shellnorm and filtering here (needs to be local because in 3D we do it on the reference volumes)
        if( present(ptcl_mask) )then
            self%nptcls  = count(ptcl_mask)                      !< the total number of particles in partition
        else
            self%nptcls  = self%pfromto(2) - self%pfromto(1) + 1 !< the total number of particles in partition
        endif
        ! input vols are prepared in preprefvol_2
        ! if( .not. vol_even%is_expanded() ) THROW_HARD('input vol_even expected to be prepared for interpolation') 
        ! if( .not. vol_odd%is_expanded()  ) THROW_HARD('input vol_odd expected to be prepared for interpolation')
        ! set pointers to projectors
        self%vol_even => vol_even
        self%vol_odd  => vol_odd
        ! allocate optimal low-pass filter & memoized sqsums
        allocate(self%optlp(params_glob%kfromto(1):params_glob%kstop), source=0.)
        allocate(self%sqsums_ptcls(1:self%nptcls), source=1.)
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
        ! allocate the rest
        allocate(self%particles(self%nptcls), self%heap_vars(params_glob%nthr))
        do i = 1,self%nptcls
            call self%particles(i)%new(self%ldim, params_glob%smpd)
        end do
        do ithr = 1,params_glob%nthr
            call self%heap_vars(ithr)%img_ref%new(self%ldim, params_glob%smpd)
            call self%heap_vars(ithr)%img_ref_tmp%new(self%ldim, params_glob%smpd)
            allocate(self%heap_vars(ithr)%frc(self%filtsz), source = 0.)
        end do
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! 2Dclass/3Drefs mapping on/off
        self%l_clsfrcs = params_glob%clsfrcs.eq.'yes'
        ! setup pxls_p_shell
        call self%setup_pxls_p_shell
        ! memoize resolution mask
        call self%memoize_resmsk
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        cartftcc_glob => self
    end subroutine new

    ! SETTERS

    subroutine reallocate_ptcls( self, nptcls, pinds )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: nptcls
        integer,                intent(in)    :: pinds(nptcls)
        integer :: i,iptcl,ik
        self%pfromto(1) = minval(pinds)
        self%pfromto(2) = maxval(pinds)
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( self%nptcls == nptcls )then
            ! just need to update particles indexing
        else
            ! re-index & reallocate
            self%nptcls = nptcls
            if( allocated(self%sqsums_ptcls) ) deallocate(self%sqsums_ptcls)
            if( allocated(self%iseven) )       deallocate(self%iseven)
            if( allocated(self%particles) )then
                do i = 1, size(self%particles)
                    call self%particles(i)%kill
                end do
                deallocate(self%particles)
            endif
            allocate( self%particles(self%nptcls), self%sqsums_ptcls(1:self%nptcls),self%iseven(1:self%nptcls) )
            do i = 1,self%nptcls
                call self%particles(i)%new(self%ldim, params_glob%smpd)
            end do
        endif
        self%sqsums_ptcls = 1.
        self%iseven       = .true.
        allocate(self%pinds(self%pfromto(1):self%pfromto(2)), source=0)
        do i = 1,self%nptcls
            iptcl = pinds(i)
            self%pinds( iptcl ) = i
        enddo
    end subroutine reallocate_ptcls

    subroutine set_ptcl( self, iptcl, img )
        class(cartft_corrcalc), intent(inout) :: self   !< this object
        integer,                intent(in)    :: iptcl  !< particle index
        class(image),           intent(in)    :: img    !< particle image
        if( .not. img%is_ft() ) THROW_HARD('input image expected to be FTed')
        if( .not. (img.eqdims.self%particles(self%pinds(iptcl))) ) THROW_HARD('inconsistent image dimensions, input vs class internal') 
        call self%particles(self%pinds(iptcl))%copy_fast(img)
        ! calculate the square sum required for correlation calculation
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
    end subroutine set_ptcl

    subroutine set_optlp( self, optlp )
        class(cartft_corrcalc), intent(inout) :: self
        real,                   intent(in)    :: optlp(params_glob%kfromto(1):params_glob%kstop)
        self%optlp(:) = optlp(:)
        self%l_filt_set = .true.
    end subroutine set_optlp

    subroutine set_eo( self, iptcl, is_even )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        logical,                intent(in)    :: is_even
        self%iseven(self%pinds(iptcl)) = is_even
    end subroutine set_eo

    ! GETTERS

    pure function get_box( self ) result( box )
        class(cartft_corrcalc), intent(in) :: self
        integer :: box
        box = self%ldim(1)
    end function get_box

    logical function exists( self )
        class(cartft_corrcalc), intent(in) :: self
        exists = self%existence
    end function exists

    logical function ptcl_iseven( self, iptcl )
        class(cartft_corrcalc), intent(in) :: self
        integer,                 intent(in) :: iptcl
        ptcl_iseven = self%iseven(self%pinds(iptcl))
    end function ptcl_iseven

    integer function get_nptcls( self )
        class(cartft_corrcalc), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    subroutine assign_pinds( self, pinds )
        class(cartft_corrcalc), intent(inout) :: self
        integer, allocatable,    intent(out)   :: pinds(:)
        pinds = self%pinds
    end subroutine assign_pinds

    ! MODIFIERS

    subroutine shellnorm_and_filter_ref( self, ref_img )
        class(cartft_corrcalc), intent(in)    :: self
        class(image),           intent(inout) :: ref_img
        if( self%l_match_filt .and. self%l_filt_set )then
            call ref_img%shellnorm_and_apply_filter_serial(self%optlp(:))
        else if( .not. params_glob%l_lpset .and. self%l_filt_set )then
            call ref_img%apply_filter_serial(self%optlp(:))
        endif
    end subroutine shellnorm_and_filter_ref

    ! MEMOIZERS

    subroutine memoize_resmsk( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: h, k, l, sh, phys(3)
        self%lims = self%particles(1)%loop_lims(2)
        if( allocated(self%resmsk) ) deallocate(self%resmsk)
        self%cmat_shape = self%particles(1)%get_array_shape()
        allocate(self%resmsk(self%cmat_shape(1),self%cmat_shape(2),self%cmat_shape(3)), source=.false.)
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
        do h = 0,params_glob%kfromto(2)
            do k = -params_glob%kfromto(2),params_glob%kfromto(2)
                if( (h==0) .and. (k>0) ) cycle
                sh = nint(hyp(real(h),real(k)))
                if( ( sh >= params_glob%kfromto(1)) .and. ( sh <= params_glob%kfromto(2)) ) then
                    self%pxls_p_shell(sh) = self%pxls_p_shell(sh) + 1.
                end if
            end do
        end do
    end subroutine setup_pxls_p_shell

    subroutine memoize_sqsum_ptcl( self, i )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: i
        self%sqsums_ptcls(i) = self%particles(i)%calc_sumsq(self%resmsk)
    end subroutine memoize_sqsum_ptcl

    ! CALCULATORS

    subroutine create_absctfmats( self, spproj, oritype, pfromto )
        use simple_ctf,        only: ctf
        use simple_sp_project, only: sp_project
        class(cartft_corrcalc),    intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        character(len=*),          intent(in)    :: oritype
        integer, optional,         intent(in)    :: pfromto(2)
        type(ctfparams) :: ctfparms(nthr_glob)
        type(ctf)       :: tfuns(nthr_glob)
        real(sp)        :: inv_ldim(3),inv(2)
        integer         :: i,h,k,iptcl,ithr,ppfromto(2),ctfmatind,phys(3)
        logical         :: present_pfromto
        present_pfromto = present(pfromto)
        ppfromto = self%pfromto
        if( present_pfromto ) ppfromto = pfromto
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%cmat_shape(1),self%cmat_shape(2),self%cmat_shape(3),1:self%nptcls), source=0.)
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr,h,k,phys,inv) schedule(static) proc_bind(close)
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
                do h = self%lims(1,1),self%lims(1,2)
                    do k = self%lims(2,1),self%lims(2,2)
                        ! compute physical address
                        phys = self%particles(1)%comp_addr_phys(h,k,0)        
                        inv  = real([h,k]) * inv_ldim(1:2)
                        if( ctfparms(ithr)%l_phaseplate )then
                            self%ctfmats(phys(1),phys(2),phys(3),self%pinds(ctfmatind))= abs(tfuns(ithr)%eval(dot_product(inv,inv), atan2(real(k),real(h)),&
                            &ctfparms(ithr)%phshift, .not.params_glob%l_wiener_part))
                        else
                            self%ctfmats(phys(1),phys(2),phys(3),self%pinds(ctfmatind))= abs(tfuns(ithr)%eval(dot_product(inv,inv), atan2(real(k),real(h)),&
                            &0., .not.params_glob%l_wiener_part))
                        endif
                    end do
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine create_absctfmats

    function calc_corr( self, iptcl, shvec, grad ) result( cc )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        real,                   intent(in)    :: shvec(2)
        real,         optional, intent(inout) :: grad(2)
        real(sp) :: cc
        integer  :: i, ithr
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        ! prep ref
        call self%heap_vars(ithr)%img_ref%fft
        ! shell normalization and filtering
        call self%shellnorm_and_filter_ref(self%heap_vars(ithr)%img_ref)
        ! multiply with CTF
        if( self%with_ctf ) call self%heap_vars(ithr)%img_ref%mul_cmat(self%ctfmats(:,:,:,i), self%resmsk)
        ! calc corr
        if( present(grad) )then
            cc = real(self%heap_vars(ithr)%img_ref%corr_grad_ad(self%particles(i), self%sqsums_ptcls(i), self%resmsk, shvec, params_glob%cc_objfun, grad), kind=sp)
        else
            cc = real(self%heap_vars(ithr)%img_ref%corr_grad_ad(self%particles(i), self%sqsums_ptcls(i), self%resmsk, shvec, params_glob%cc_objfun), kind=sp)
        endif
    end function calc_corr

    subroutine project_and_correlate( self, iptcl, o, corr, grad )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        real,                   intent(inout) :: corr
        real, optional,         intent(inout) :: grad(2)
        type(projector), pointer :: vol_ptr => null()
        logical :: iseven, present_grad
        integer :: ithr
        present_grad = present(grad)
        iseven = self%ptcl_iseven(iptcl)
        if( iseven )then
            vol_ptr => self%vol_even
        else
            vol_ptr => self%vol_odd
        endif
        ithr = omp_get_thread_num() + 1
        if( present_grad )then
            call vol_ptr%fproject_serial(o, self%heap_vars(ithr)%img_ref, params_glob%kstop)
            corr = self%calc_corr(iptcl, o%get_2Dshift(), grad(:))
        else
            call vol_ptr%fproject_serial(o, self%heap_vars(ithr)%img_ref, params_glob%kstop)
            corr = self%calc_corr(iptcl, o%get_2Dshift())
        endif
    end subroutine project_and_correlate

    ! DESTRUCTOR

    subroutine kill( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: i, ithr
        if( self%existence )then
            self%ldim     = 0
            self%vol_even => null()
            self%vol_odd  => null()
            if( allocated(self%pinds)        ) deallocate(self%pinds)
            if( allocated(self%pxls_p_shell) ) deallocate(self%pxls_p_shell)
            if( allocated(self%sqsums_ptcls) ) deallocate(self%sqsums_ptcls)
            if( allocated(self%ctfmats)      ) deallocate(self%ctfmats)
            if( allocated(self%optlp)        ) deallocate(self%optlp)
            if( allocated(self%iseven)       ) deallocate(self%iseven)
            do i = 1,self%nptcls
                call self%particles(i)%kill
            end do
            do ithr = 1,params_glob%nthr
                call self%heap_vars(ithr)%img_ref%kill
                call self%heap_vars(ithr)%img_ref_tmp%kill
                deallocate(self%heap_vars(ithr)%frc)
                self%heap_vars(ithr)%frc => null()
            end do
            deallocate(self%particles, self%heap_vars)
            self%l_filt_set = .false.
            self%existence  = .false.
        endif
    end subroutine kill

end module simple_cartft_corrcalc

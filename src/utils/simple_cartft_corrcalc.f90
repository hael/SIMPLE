module simple_cartft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
use simple_projector,  only: projector
implicit none

public :: cartft_corrcalc, cartftcc_glob
private
#include "simple_local_flags.inc"

type :: cartft_corrcalc
    private
    type(projector),               pointer     :: vol_even => null(), vol_odd => null() ! prepared e/o vols
    integer                                    :: nptcls     = 1         !< # particles in partition (logically indexded [fromp,top])
    integer                                    :: filtsz     = 0         !< Nyqvist limit
    integer                                    :: pfromto(2) = 0         !< particle index range
    integer                                    :: ldim(3)    = 0         !< logical dimensions of original cartesian image
    integer                                    :: lims(2,2)  = 0         !< resolution mask limits
    integer,                       allocatable :: pinds(:)               !< index array (to reduce memory when frac_update < 1)
    real,                          allocatable :: pxls_p_shell(:)        !< number of (cartesian) pixels per shell
    real(sp),                      allocatable :: ctfmats(:,:,:)         !< expand set of CTF matrices (for efficient parallel exec)
    real(sp),                      allocatable :: optlp(:)               !< references optimal filter
    complex(kind=c_float_complex), allocatable :: particles(:,:,:)       !< particle Fourier components (h,k,iptcl)
    complex(kind=c_float_complex), allocatable :: ref_heap(:,:,:)        !< reference Fourier components (h,k,ithr)
    logical,                       allocatable :: resmsk(:,:)            !< resolution mask
    logical,                       allocatable :: iseven(:)              !< e/o assignment for gold-standard FSC
    logical                                    :: l_clsfrcs    = .false. !< CLS2D/3DRefs flag
    logical                                    :: l_match_filt = .false. !< matched filter flag
    logical                                    :: l_filt_set   = .false. !< to indicate whether filter is set
    logical                                    :: with_ctf     = .false. !< CTF flag
    logical                                    :: existence    = .false. !< to indicate existence
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
    ! MEMOIZERS
    procedure          :: prep4shift_srch
    procedure, private :: setup_resmsk_and_pxls_p_shell
    ! CALCULATORS
    procedure          :: create_absctfmats
    procedure          :: project_and_correlate
    procedure          :: corr_shifted
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
        ! set pointers to projectors
        self%vol_even => vol_even
        self%vol_odd  => vol_odd
        ! allocate optimal low-pass filter
        allocate(self%optlp(params_glob%kfromto(1):params_glob%kfromto(2)), source=0.)
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
        ! Fourier index limits
        self%lims(1,1) =  0
        self%lims(1,2) =  params_glob%kfromto(2)
        self%lims(2,1) = -params_glob%kfromto(2)
        self%lims(2,2) =  params_glob%kfromto(2)
        ! allocate the rest
        allocate(self%particles(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),self%nptcls),&
        &         self%ref_heap(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),nthr_glob), source=cmplx(0.,0.))
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! 2Dclass/3Drefs mapping on/off
        self%l_clsfrcs = params_glob%clsfrcs.eq.'yes'
        ! setup resolution mask pxls_p_shell
        call self%setup_resmsk_and_pxls_p_shell
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
            if( allocated(self%iseven)    ) deallocate(self%iseven)
            if( allocated(self%particles) ) deallocate(self%particles)
            allocate(self%iseven(1:self%nptcls))
            allocate(self%particles(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),self%nptcls), source=cmplx(0.,0.))
        endif
        self%iseven = .true.
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
        integer :: ldim(3), phys(3), h, k
        if( .not. img%is_ft() ) THROW_HARD('input image expected to be FTed')
        ldim = img%get_ldim()
        if( .not. all(self%ldim .eq. ldim) )then
            THROW_HARD('inconsistent image dimensions, input vs class internal')
        endif
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                phys = img%comp_addr_phys(h,k)
                self%particles(h,k,iptcl) = img%get_cmat_at(phys(1),phys(2),1)
            end do
        end do
    end subroutine set_ptcl

    subroutine set_optlp( self, optlp )
        class(cartft_corrcalc), intent(inout) :: self
        real,                   intent(in)    :: optlp(params_glob%kfromto(1):params_glob%kfromto(2))
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

    subroutine prep4shift_srch( self, iptcl, o )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        type(projector), pointer :: vol_ptr => null()
        integer :: ithr
        logical :: iseven
        iseven = self%ptcl_iseven(iptcl)
        if( iseven )then
            vol_ptr => self%vol_even
        else
            vol_ptr => self%vol_odd
        endif
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! put reference projection in the heap
        call vol_ptr%fproject_serial(o, self%lims, self%ref_heap(:,:,ithr), self%resmsk)
    end subroutine prep4shift_srch

    subroutine setup_resmsk_and_pxls_p_shell( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: h,k,sh
        if( allocated(self%pxls_p_shell) ) deallocate(self%pxls_p_shell)
        allocate(self%pxls_p_shell(params_glob%kfromto(1):params_glob%kfromto(2)), source=0.)
        if( allocated(self%resmsk) ) deallocate(self%resmsk)
        allocate(self%resmsk(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2)), source=.false.)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                if( (h==0) .and. (k>0) ) cycle
                sh = nint(hyp(real(h),real(k)))
                if( ( sh >= params_glob%kfromto(1)) .and. ( sh <= params_glob%kfromto(2)) ) then
                    self%pxls_p_shell(sh) = self%pxls_p_shell(sh) + 1.
                    self%resmsk(h,k)      = .true.
                end if
            end do
        end do
    end subroutine setup_resmsk_and_pxls_p_shell

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
        integer         :: i,h,k,iptcl,ithr,ppfromto(2),ctfmatind
        logical         :: present_pfromto
        present_pfromto = present(pfromto)
        ppfromto = self%pfromto
        if( present_pfromto ) ppfromto = pfromto
        if( allocated(self%ctfmats) ) deallocate(self%ctfmats)
        allocate(self%ctfmats(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1:self%nptcls), source=1.)
        if(.not. self%with_ctf ) return
        inv_ldim = 1./real(self%ldim)
        !$omp parallel do default(shared) private(i,iptcl,ctfmatind,ithr,h,k,inv) schedule(static) proc_bind(close)
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
                        inv  = real([h,k]) * inv_ldim(1:2)
                        if( ctfparms(ithr)%l_phaseplate )then
                            self%ctfmats(h,k,self%pinds(ctfmatind))= abs(tfuns(ithr)%eval(dot_product(inv,inv), atan2(real(k),real(h)),&
                            &ctfparms(ithr)%phshift, .not.params_glob%l_wiener_part))
                        else
                            self%ctfmats(h,k,self%pinds(ctfmatind))= abs(tfuns(ithr)%eval(dot_product(inv,inv), atan2(real(k),real(h)),&
                            &0., .not.params_glob%l_wiener_part))
                        endif
                    end do
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine create_absctfmats

    function project_and_correlate( self, iptcl, o ) result( corr )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        type(projector), pointer :: vol_ptr => null()
        real    :: corr
        logical :: iseven
        integer :: i
        i      = self%pinds(iptcl)
        iseven = self%ptcl_iseven(iptcl)
        if( iseven )then
            vol_ptr => self%vol_even
        else
            vol_ptr => self%vol_odd
        endif
        corr = vol_ptr%fproject_and_correlate_serial(o, self%lims, self%particles(:,:,i), self%ctfmats(:,:,i), self%resmsk)
    end function project_and_correlate

    function corr_shifted( self, iptcl, o, shvec ) result( corr )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        real,                   intent(in)    :: shvec(2)
        complex(kind=c_float_complex), pointer :: cmat_ref(:,:,:), cmat_ptcl(:,:,:)
        integer  :: i, h, k, hphys, kphys, ithr
        complex  :: ref_comp, ptcl_comp, sh_comp
        real(dp) :: sh(2), arg, ck, sk, hcos(self%lims(1,1):self%lims(1,2)), hsin(self%lims(1,1):self%lims(1,2))
        real(sp) :: cc(3), corr, shconst
        logical  :: iseven
        ! physical particle index
        i = self%pinds(iptcl)
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! calculate constant factor (assumes self%ldim(1) == self%ldim(2))
        if( is_even(self%ldim(1)) )then
            shconst = PI / real(self%ldim(1)/2.)
        else
            shconst = PI / real((self%ldim(1)-1)/2.)
        endif
        ! optimized shift caluclation following (shift2Dserial_1 in image class)
        sh = real(shvec * shconst,dp)
        do h = self%lims(1,1),self%lims(1,2)
            arg     = real(h,dp) * sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        cc(:) = 0.
        do k = self%lims(2,1),self%lims(2,2)
            kphys = k + 1 + merge(self%ldim(2),0,k<0)
            arg   = real(k,dp) * sh(2)
            ck    = dcos(arg)
            sk    = dsin(arg)
            do h = self%lims(1,1),self%lims(1,2)
                if( .not. self%resmsk(h,k) ) cycle
                hphys   = h + 1
                sh_comp = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                ! retrieve reference component
                if( h > 0 )then
                    ref_comp =       self%ref_heap(h,k,ithr)  * self%ctfmats(h,k,i)
                else
                    ref_comp = conjg(self%ref_heap(h,k,ithr)) * self%ctfmats(h,k,i)
                endif
                ! shift the particle Fourier component
                ptcl_comp = self%particles(h,k,i) * sh_comp
                ! update cross product
                cc(1) = cc(1) + real(ref_comp  * conjg(ptcl_comp))
                ! update normalization terms
                cc(2) = cc(2) + real(ref_comp  * conjg(ref_comp))
                cc(3) = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
            end do
        end do
        corr = norm_corr(cc(1),cc(2),cc(3))
    end function corr_shifted

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
            if( allocated(self%ctfmats)      ) deallocate(self%ctfmats)
            if( allocated(self%optlp)        ) deallocate(self%optlp)
            if( allocated(self%iseven)       ) deallocate(self%iseven)
            if( allocated(self%particles)    ) deallocate(self%particles)
            if( allocated(self%ref_heap)     ) deallocate(self%ref_heap)
            if( allocated(self%resmsk)       ) deallocate(self%resmsk)
            self%l_filt_set = .false.
            self%existence  = .false.
        endif
    end subroutine kill

end module simple_cartft_corrcalc

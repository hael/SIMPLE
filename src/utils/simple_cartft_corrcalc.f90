module simple_cartft_corrcalc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
use simple_projector,  only: projector
implicit none

public :: cartft_corrcalc, cftcc_glob
private
#include "simple_local_flags.inc"

type :: cartft_corrcalc
    private
    type(projector),               pointer     :: vol_even => null(), vol_odd => null() ! prepared e/o vols
    real,                          pointer     :: sigma2_noise(:,:) => null() !< for euclidean distances
    integer                                    :: nptcls     = 1         !< # particles in partition (logically indexded [fromp,top])
    integer                                    :: filtsz     = 0         !< Nyqvist limit
    integer                                    :: pfromto(2) = 0         !< particle index range
    integer                                    :: ldim(3)    = 0         !< logical dimensions of original cartesian image
    integer                                    :: lims(2,2)  = 0         !< resolution mask limits
    integer,                       allocatable :: pinds(:)               !< index array (to reduce memory when frac_update < 1)
    real,                          allocatable :: pxls_p_shell(:)        !< number of (cartesian) pixels per shell
    real(sp),                      allocatable :: ctfmats(:,:,:)         !< logically indexed CTF matrices (for efficient parallel exec)
    real(sp),                      allocatable :: ref_filt_w(:,:,:)      !< reference filter weights
    complex(kind=c_float_complex), allocatable :: particles(:,:,:)       !< particle Fourier components (h,k,iptcl)
    complex(kind=c_float_complex), allocatable :: references(:,:,:)      !< reference Fourier components (h,k,ithr) or (h,k,iptcl) when expanded
    logical,                       allocatable :: resmsk(:,:)            !< resolution mask
    logical,                       allocatable :: iseven(:)              !< e/o assignment for gold-standard FSC
    logical                                    :: l_clsfrcs    = .false. !< CLS2D/3DRefs flag
    logical                                    :: l_match_filt = .false. !< matched filter flag
    logical                                    :: with_ctf     = .false. !< CTF flag
    logical                                    :: existence    = .false. !< to indicate existence
    contains
    ! CONSTRUCTORS
    procedure          :: new
    procedure          :: new_dev
    ! ON-THE-FLY REALLOCATION
    procedure          :: reallocate_ptcls
    procedure          :: expand_refs
    ! SETTERS
    procedure          :: set_ptcl
    procedure          :: set_ref
    procedure, private :: set_ref_optlp_1
    procedure, private :: set_ref_optlp_2
    generic            :: set_ref_optlp => set_ref_optlp_1, set_ref_optlp_2
    procedure          :: set_eo
    ! GETTERS
    procedure          :: get_box
    procedure          :: exists
    procedure          :: ptcl_iseven
    procedure          :: get_nptcls
    procedure          :: assign_pinds
    ! MEMOIZERS
    procedure, private :: setup_resmsk_and_pxls_p_shell
    procedure          :: create_absctfmats
    procedure          :: prep4shift_srch
    ! procedure          :: prep4parallel_shift_srch
    ! CALCULATORS
    procedure, private :: project_and_correlate_1
    procedure, private :: project_and_correlate_2
    generic            :: project_and_correlate => project_and_correlate_1, project_and_correlate_2
    procedure          :: project_and_srch_shifts
    procedure, private :: corr_shifted_1
    procedure, private :: corr_shifted_2
    generic            :: corr_shifted => corr_shifted_1, corr_shifted_2
    procedure, private :: corr_shifted_cc_1
    procedure, private :: corr_shifted_cc_2
    procedure, private :: corr_shifted_euclid_1
    procedure, private :: corr_shifted_euclid_2
    procedure          :: corr_shifted_ad
    generic            :: ori_chance => ori_chance_1, ori_chance_2
    procedure, private :: ori_chance_1
    procedure, private :: ori_chance_2
    procedure          :: assign_sigma2_noise
    procedure          :: calc_sigma_contrib
    ! DESTRUCTOR
    procedure          :: kill
end type cartft_corrcalc

! CLASS PARAMETERS/VARIABLES
class(cartft_corrcalc), pointer :: cftcc_glob => null()

contains

    ! CONSTRUCTORS

    subroutine new( self, vol_even, vol_odd, pfromto, l_match_filt, ptcl_mask, eoarr )
        class(cartft_corrcalc), target, intent(inout) :: self
        type(projector),        target, intent(in)    :: vol_even, vol_odd
        integer,                        intent(in)    :: pfromto(2)
        logical,                        intent(in)    :: l_match_filt
        logical, optional,              intent(in)    :: ptcl_mask(pfromto(1):pfromto(2))
        integer, optional,              intent(in)    :: eoarr(pfromto(1):pfromto(2))
        logical  :: even_dims, test(2)
        integer  :: i, cnt, ithr, lims(3,2)
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
        self%filtsz       = fdim(params_glob%box) - 1
        self%l_match_filt = l_match_filt                         !< do shellnorm and filtering here
        if( present(ptcl_mask) )then
            self%nptcls  = count(ptcl_mask)                      !< the total number of particles in partition
        else
            self%nptcls  = self%pfromto(2) - self%pfromto(1) + 1 !< the total number of particles in partition
        endif
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0 )
        ! set pointers to projectors
        self%vol_even => vol_even
        self%vol_odd  => vol_odd
        ! generate particle indices
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
        ! 2Dclass/3Drefs mapping on/off
        self%l_clsfrcs = params_glob%clsfrcs.eq.'yes'
        if( self%l_clsfrcs )then
            allocate(self%references(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1:self%nptcls), source=cmplx(0.,0.))
            allocate(self%ref_filt_w(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1:self%nptcls), source=1.)
        else
            allocate(self%references(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1:nthr_glob), source=cmplx(0.,0.))
            allocate(self%ref_filt_w(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1),           source=1.)
        endif
        ! allocate the rest
        allocate(self%particles(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),self%nptcls), source=cmplx(0.,0.))
        allocate(self%ctfmats(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1:self%nptcls), source=1.)
        ! set CTF flag
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! setup resolution mask pxls_p_shell
        call self%setup_resmsk_and_pxls_p_shell
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        cftcc_glob => self
    end subroutine new

    subroutine new_dev( self, pfromto )
        class(cartft_corrcalc), target, intent(inout) :: self
        integer,                        intent(in)    :: pfromto(2)
        logical      :: even_dims, test(2)
        integer      :: i, cnt, ithr, lims(3,2)
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
        self%filtsz       = fdim(params_glob%box) - 1
        self%l_match_filt = .false.
        self%nptcls       = self%pfromto(2) - self%pfromto(1) + 1 !< the total number of particles in partition
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0 )
        self%pinds = (/(i,i=1,self%nptcls)/)
        ! eo assignment
        allocate( self%iseven(1:self%nptcls), source=.true. )
        ! Fourier index limits
        self%lims(1,1) =  0
        self%lims(1,2) =  params_glob%kfromto(2)
        self%lims(2,1) = -params_glob%kfromto(2)
        self%lims(2,2) =  params_glob%kfromto(2)
        ! allocate the rest
        allocate(self%particles(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),self%nptcls),&
        &       self%references(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),nthr_glob), source=cmplx(0.,0.))
        allocate(self%ref_filt_w(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1),&
        &           self%ctfmats(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),1:self%nptcls), source=1.)
        ! set CTF flag
        self%with_ctf = .false.
        ! 2Dclass/3Drefs mapping on/off
        self%l_clsfrcs = .false.
        ! setup resolution mask pxls_p_shell
        call self%setup_resmsk_and_pxls_p_shell
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        cftcc_glob => self
    end subroutine new_dev

    ! ON-THE-FLY REALLOCATION

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

    ! expand references such that we have one per particle (for efficient parallel shift search)
    subroutine expand_refs( self )
        class(cartft_corrcalc), intent(inout) :: self
        if( allocated(self%references) )then
            if( self%nptcls == size(self%references, dim=3) )then
                return
            else
                deallocate(self%references)
            endif
        endif
        allocate(self%references(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2),self%nptcls), source=cmplx(0.,0.))
    end subroutine expand_refs

    ! SETTERS

    subroutine set_ptcl( self, iptcl, img )
        class(cartft_corrcalc), intent(inout) :: self   !< this object
        integer,                intent(in)    :: iptcl  !< particle index
        class(image),           intent(in)    :: img    !< particle image
        integer :: ldim(3), h, k, phys1, phys2
        if( .not. img%is_ft() ) THROW_HARD('input image expected to be FTed')
        ldim = img%get_ldim()
        if( .not. all(self%ldim .eq. ldim) )then
            THROW_HARD('inconsistent image dimensions, input vs class internal')
        endif
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                if( h .ge. 0 )then
                    phys1 = h + 1
                    phys2 = k + 1 + merge(self%ldim(2), 0, k < 0)
                else
                    phys1 = -h + 1
                    phys2 = -k + 1 + MERGE(self%ldim(2),0, -k < 0)
                endif
                self%particles(h,k,self%pinds(iptcl)) = img%get_cmat_at(phys1, phys2, 1)
            end do
        end do
    end subroutine set_ptcl

    subroutine set_ref( self, img )
        class(cartft_corrcalc), intent(inout) :: self   !< this object
        class(image),           intent(in)    :: img    !< particle image
        integer :: ldim(3), h, k, ithr, phys1, phys2
        if( .not. img%is_ft() ) THROW_HARD('input image expected to be FTed')
        ldim = img%get_ldim()
        if( .not. all(self%ldim .eq. ldim) )then
            THROW_HARD('inconsistent image dimensions, input vs class internal')
        endif
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                if( h .ge. 0 )then
                    phys1 = h + 1
                    phys2 = k + 1 + merge(self%ldim(2), 0, k < 0)
                else
                    phys1 = -h + 1
                    phys2 = -k + 1 + MERGE(self%ldim(2),0, -k < 0)
                endif
                self%references(h,k,1) = img%get_cmat_at(phys1, phys2, 1)
            end do
        end do
        do ithr = 2,nthr_glob
            self%references(:,:,ithr) = self%references(:,:,1)
        end do
    end subroutine set_ref

    subroutine set_ref_optlp_1( self, ref_optlp )
        class(cartft_corrcalc), intent(inout) :: self
        real,                   intent(in)    :: ref_optlp(params_glob%kfromto(1):params_glob%kfromto(2))
        integer :: h,k,sh,sqlp,sqarg
        self%ref_filt_w(:,:,1) = 0.
        sqlp = params_glob%kfromto(2) * params_glob%kfromto(2)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                if( (h==0) .and. (k==0) ) cycle
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                sh = nint(hyp(real(h),real(k)))
                self%ref_filt_w(h,k,1) = ref_optlp(sh)
            end do
        end do
    end subroutine set_ref_optlp_1

    subroutine set_ref_optlp_2( self, iref, ref_optlp )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iref
        real,                   intent(in)    :: ref_optlp(params_glob%kfromto(1):params_glob%kfromto(2))
        integer :: h,k,sh,sqlp,sqarg
        if( iref < 1 .or. iref > self%nptcls ) THROW_HARD('iref '//int2str(iref)//' is out of bounds')
        self%ref_filt_w(:,:,iref) = 0.
        sqlp = params_glob%kfromto(2) * params_glob%kfromto(2)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                if( (h==0) .and. (k==0) ) cycle
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                sh = nint(hyp(real(h),real(k)))
                self%ref_filt_w(h,k,iref) = ref_optlp(sh)
            end do
        end do
    end subroutine set_ref_optlp_2

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
        integer, allocatable,    intent(out)  :: pinds(:)
        pinds = self%pinds
    end subroutine assign_pinds

    ! MEMOIZERS

    subroutine setup_resmsk_and_pxls_p_shell( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: h,k,sh,sqlp,sqarg
        if( allocated(self%pxls_p_shell) ) deallocate(self%pxls_p_shell)
        allocate(self%pxls_p_shell(params_glob%kfromto(1):params_glob%kfromto(2)), source=0.)
        if( allocated(self%resmsk) ) deallocate(self%resmsk)
        allocate(self%resmsk(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2)), source=.false.)
        sqlp = params_glob%kfromto(2) * params_glob%kfromto(2)
        do h = self%lims(1,1),self%lims(1,2)
            do k = self%lims(2,1),self%lims(2,2)
                if( (h==0) .and. (k==0) ) cycle
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                self%resmsk(h,k) = .true.
                sh = nint(hyp(real(h),real(k)))
                self%pxls_p_shell(sh) = self%pxls_p_shell(sh) + 1.
            end do
        end do
    end subroutine setup_resmsk_and_pxls_p_shell

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
        if( self%l_match_filt )then
            call vol_ptr%fproject_shellnorm_serial(o, self%lims, self%references(:,:,ithr), self%resmsk(:,:))
        else
            call vol_ptr%fproject_serial(o, self%lims, self%references(:,:,ithr), self%resmsk(:,:))
        endif
    end subroutine prep4shift_srch

    ! subroutine prep4parallel_shift_srch( self, iptcl, os )
    !     class(cartft_corrcalc), intent(inout) :: self
    !     integer,                intent(in)    :: iptcl
    !     class(oris),            intent(in)    :: os
    !     type(projector), pointer :: vol_ptr => null()
    !     logical :: iseven
    !     integer :: i
    !     i      = self%pinds(iptcl)
    !     iseven = self%ptcl_iseven(iptcl)
    !     if( iseven )then
    !         vol_ptr => self%vol_even
    !     else
    !         vol_ptr => self%vol_odd
    !     endif
    !     ! expand references such that we have one per particle (for efficient parallel shift search)
    !     call self%expand_refs
    !     ! make one reference projection per particle
    !     call vol_ptr%fproject4cartftcc(os, self%lims, self%nptcls, self%references(:,:,:), self%resmsk(:,:))
    ! end subroutine prep4parallel_shift_srch

    ! CALCULATORS

    function project_and_correlate_1( self, iptcl, o ) result( corr )
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
        if( self%l_match_filt )then
            if( self%l_clsfrcs )then
                corr = vol_ptr%fproject_shellnorm_correlate_serial(o, self%lims, self%particles(:,:,i),&
                &self%ctfmats(:,:,i), self%resmsk(:,:), self%ref_filt_w(:,:,i))
            else
                corr = vol_ptr%fproject_shellnorm_correlate_serial(o, self%lims, self%particles(:,:,i),&
                &self%ctfmats(:,:,i), self%resmsk(:,:), self%ref_filt_w(:,:,1))
            endif
        else
            if( self%l_clsfrcs )then
                corr = vol_ptr%fproject_correlate_serial(o, self%lims, self%particles(:,:,i),&
                &self%ctfmats(:,:,i), self%resmsk(:,:), self%ref_filt_w(:,:,i))
            else
                corr = vol_ptr%fproject_correlate_serial(o, self%lims, self%particles(:,:,i),&
                &self%ctfmats(:,:,i), self%resmsk(:,:), self%ref_filt_w(:,:,1))
            endif
        endif
    end function project_and_correlate_1

    function project_and_correlate_2( self, iptcl, o, shvec ) result( corr )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        real,                   intent(in)    :: shvec(2)
        type(projector),        pointer       :: vol_ptr => null()
        logical :: iseven
        integer :: ithr
        real    :: corr
        iseven = self%ptcl_iseven(iptcl)
        if( iseven )then
            vol_ptr => self%vol_even
        else
            vol_ptr => self%vol_odd
        endif
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! put reference projection in the heap
        if( self%l_match_filt )then
            call vol_ptr%fproject_shellnorm_serial(o, self%lims, self%references(:,:,ithr), self%resmsk(:,:))
        else
            call vol_ptr%fproject_serial(o, self%lims, self%references(:,:,ithr), self%resmsk(:,:))
        endif
        call self%corr_shifted(iptcl, shvec, corr)
    end function project_and_correlate_2

    subroutine project_and_srch_shifts( self, iptcl, o, nsample, trs, shvec, corr_best, nevals )
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        class(ori),             intent(in)    :: o
        integer,                intent(in)    :: nsample
        real,                   intent(in)    :: trs
        real,                   intent(inout) :: shvec(2), corr_best
        integer,                intent(inout) :: nevals
        type(projector),        pointer       :: vol_ptr => null()
        logical :: iseven
        integer :: ithr, isample
        real    :: sigma, xshift, yshift, xavg, yavg, corr
        iseven = self%ptcl_iseven(iptcl)
        if( iseven )then
            vol_ptr => self%vol_even
        else
            vol_ptr => self%vol_odd
        endif
        ! get thread index
        ithr = omp_get_thread_num() + 1
        ! put reference projection in the heap
        if( self%l_match_filt )then
            call vol_ptr%fproject_shellnorm_serial(o, self%lims, self%references(:,:,ithr), self%resmsk(:,:))
        else
            call vol_ptr%fproject_serial(o, self%lims, self%references(:,:,ithr), self%resmsk(:,:))
        endif
        sigma = trs / 2. ! 2 sigma (soft) criterion, fixed for now
        call self%corr_shifted(iptcl, shvec, corr_best)
        nevals = 0
        do isample = 1,nsample
            xshift = gasdev(shvec(1), sigma)
            yshift = gasdev(shvec(2), sigma)
            call self%corr_shifted(iptcl, [xshift,yshift], corr)
            nevals = nevals + 1
            if( corr > corr_best )then
                corr_best = corr
                shvec(1)  = xshift
                shvec(2)  = yshift
                exit ! SHC
            endif
        end do
    end subroutine project_and_srch_shifts

    subroutine corr_shifted_1( self, iptcl, shvec, corr )
        class(cartft_corrcalc), target, intent(inout)  :: self
        integer,                        intent(in)     :: iptcl
        real,                           intent(in)     :: shvec(2)
        real,                           intent(inout)  :: corr
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call corr_shifted_cc_1( self, iptcl, shvec, corr )
            case(OBJFUN_EUCLID)
                call corr_shifted_euclid_1( self, iptcl, shvec, corr )
        end select
    end subroutine corr_shifted_1

    subroutine corr_shifted_2( self, iptcl, shvec, corr, grad )
        class(cartft_corrcalc), target, intent(inout)  :: self
        integer,                        intent(in)     :: iptcl
        real,                           intent(in)     :: shvec(2)
        real,                           intent(inout)  :: corr
        real,                           intent(inout)  :: grad(2)
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call corr_shifted_cc_2( self, iptcl, shvec, corr, grad )
            case(OBJFUN_EUCLID)
                call corr_shifted_euclid_2( self, iptcl, shvec, corr, grad )
        end select
    end subroutine corr_shifted_2

    subroutine corr_shifted_cc_1( self, iptcl, shvec, corr )
        class(cartft_corrcalc), target, intent(inout)  :: self
        integer,                        intent(in)     :: iptcl
        real,                           intent(in)     :: shvec(2)
        real,                           intent(inout)  :: corr
        integer  :: i, h, k, ithr
        complex  :: ref_comp, sh_comp, ptcl_comp, ref_ptcl_sh
        real(dp) :: sh(2), arg, hcos(self%lims(1,1):self%lims(1,2))
        real(dp) :: hsin(self%lims(1,1):self%lims(1,2)), ck, sk
        real(sp) :: cc(3), shconst
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
        ! optimized shift calculation following (shift2Dserial_1 in image class)
        sh = real(shvec * shconst,dp)
        do h = self%lims(1,1),self%lims(1,2)
            arg     = real(h,dp) * sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        cc(:) = 0.
        if( self%l_clsfrcs )then
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp  = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,i)
                    ! shift the particle Fourier component
                    ptcl_comp = self%particles(h,k,i)
                    ! update cross product
                    cc(1) = cc(1) + real(ref_comp  * conjg(ptcl_comp * sh_comp))
                    ! update normalization terms
                    cc(2) = cc(2) + real(ref_comp  * conjg(ref_comp))
                    cc(3) = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                end do
            end do
        else
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp  = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,1)
                    ! shift the particle Fourier component
                    ptcl_comp = self%particles(h,k,i)
                    ! update cross product
                    cc(1) = cc(1) + real(ref_comp  * conjg(ptcl_comp * sh_comp))
                    ! update normalization terms
                    cc(2) = cc(2) + real(ref_comp  * conjg(ref_comp))
                    cc(3) = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                end do
            end do
        endif
        corr = norm_corr(cc(1),cc(2),cc(3))
    end subroutine corr_shifted_cc_1

    subroutine corr_shifted_cc_2( self, iptcl, shvec, corr, grad )
        class(cartft_corrcalc), target, intent(inout)  :: self
        integer,                        intent(in)     :: iptcl
        real,                           intent(in)     :: shvec(2)
        real,                           intent(inout)  :: corr
        real,                           intent(inout)  :: grad(2)
        integer  :: i, h, k, ithr
        complex  :: ref_comp, sh_comp, ptcl_comp, ref_ptcl_sh
        real(dp) :: sh(2), arg, hcos(self%lims(1,1):self%lims(1,2)), hsin(self%lims(1,1):self%lims(1,2)), ck, sk
        real(sp) :: cc(3), shconst
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
        ! optimized shift calculation following (shift2Dserial_1 in image class)
        sh = real(shvec * shconst,dp)
        do h = self%lims(1,1),self%lims(1,2)
            arg     = real(h,dp) * sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        cc(:)   = 0.
        grad(:) = 0.
        if( self%l_clsfrcs )then
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp     = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp    = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,i)
                    ! shift the particle Fourier component
                    ptcl_comp   = self%particles(h,k,i)
                    ! update cross product
                    ref_ptcl_sh = ref_comp  * conjg(ptcl_comp * sh_comp)
                    cc(1)       = cc(1) + realpart(ref_ptcl_sh)
                    ! update normalization terms
                    cc(2)       = cc(2) + real(ref_comp  * conjg(ref_comp))
                    cc(3)       = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                    ! update the gradient
                    ref_ptcl_sh = imagpart(ref_ptcl_sh) * shconst
                    grad(1)     = grad(1) + real(ref_ptcl_sh) * h
                    grad(2)     = grad(2) + real(ref_ptcl_sh) * k
                end do
            end do
        else
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp     = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp    = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,1)
                    ! shift the particle Fourier component
                    ptcl_comp   = self%particles(h,k,i)
                    ! update cross product
                    ref_ptcl_sh = ref_comp  * conjg(ptcl_comp * sh_comp)
                    cc(1)       = cc(1) + realpart(ref_ptcl_sh)
                    ! update normalization terms
                    cc(2)       = cc(2) + real(ref_comp  * conjg(ref_comp))
                    cc(3)       = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                    ! update the gradient
                    ref_ptcl_sh = imagpart(ref_ptcl_sh) * shconst
                    grad(1)     = grad(1) + real(ref_ptcl_sh) * h
                    grad(2)     = grad(2) + real(ref_ptcl_sh) * k
                end do
            end do
        endif
        grad(1) = norm_corr(grad(1),cc(2), cc(3))
        grad(2) = norm_corr(grad(2),cc(2), cc(3))
        corr    = norm_corr(cc(1),  cc(2), cc(3))
    end subroutine corr_shifted_cc_2

    subroutine corr_shifted_euclid_1( self, iptcl, shvec, corr )
        class(cartft_corrcalc), target, intent(inout)  :: self
        integer,                        intent(in)     :: iptcl
        real,                           intent(in)     :: shvec(2)
        real,                           intent(inout)  :: corr
        integer  :: i, h, k, ithr
        complex  :: ref_comp, sh_comp, ptcl_comp, ref_ptcl_sh
        real(dp) :: sh(2), arg, hcos(self%lims(1,1):self%lims(1,2))
        real(dp) :: hsin(self%lims(1,1):self%lims(1,2)), ck, sk
        real(sp) :: cc(3), shconst
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
        ! optimized shift calculation following (shift2Dserial_1 in image class)
        sh = real(shvec * shconst,dp)
        do h = self%lims(1,1),self%lims(1,2)
            arg     = real(h,dp) * sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        cc(:) = 0.
        if( self%l_clsfrcs )then
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp  = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,i)
                    ! shift the particle Fourier component
                    ptcl_comp = self%particles(h,k,i) * sh_comp
                    ! update euclidean difference
                    cc(1) = cc(1) + real((ref_comp - ptcl_comp) * conjg(ref_comp - ptcl_comp))
                    ! update normalization terms
                    cc(2) = cc(2) + real( ref_comp * conjg( ref_comp))
                    cc(3) = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                end do
            end do
        else
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp  = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,1)
                    ! shift the particle Fourier component
                    ptcl_comp = self%particles(h,k,i) * sh_comp
                    ! update euclidean difference
                    cc(1) = cc(1) + real((ref_comp - ptcl_comp) * conjg(ref_comp - ptcl_comp))
                    ! update normalization terms
                    cc(2) = cc(2) + real( ref_comp * conjg( ref_comp))
                    cc(3) = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                end do
            end do
        endif
        corr = 1 - cc(1)/(cc(2)+cc(3))
    end subroutine corr_shifted_euclid_1

    subroutine corr_shifted_euclid_2( self, iptcl, shvec, corr, grad )
        class(cartft_corrcalc), target, intent(inout)  :: self
        integer,                        intent(in)     :: iptcl
        real,                           intent(in)     :: shvec(2)
        real,                           intent(inout)  :: corr
        real,                           intent(inout)  :: grad(2)
        integer  :: i, h, k, ithr
        complex  :: ref_comp, sh_comp, ptcl_comp, diff_comp
        real(dp) :: sh(2), arg, hcos(self%lims(1,1):self%lims(1,2)), hsin(self%lims(1,1):self%lims(1,2)), ck, sk
        real(sp) :: cc(3), shconst
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
        ! optimized shift calculation following (shift2Dserial_1 in image class)
        sh = real(shvec * shconst,dp)
        do h = self%lims(1,1),self%lims(1,2)
            arg     = real(h,dp) * sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        cc(:)   = 0.
        grad(:) = 0.
        if( self%l_clsfrcs )then
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp   = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp  = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,i)
                    ! shift the particle Fourier component
                    ptcl_comp = self%particles(h,k,i) * sh_comp
                    diff_comp = ref_comp  - ptcl_comp
                    ! update euclidean difference
                    cc(1)     = cc(1) + real(diff_comp * conjg(diff_comp))
                    ! update normalization terms
                    cc(2)     = cc(2) + real( ref_comp * conjg( ref_comp))
                    cc(3)     = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                    ! update the gradient
                    ptcl_comp = 2 * shconst * imagpart(ptcl_comp * conjg(diff_comp))
                    grad(1)   = grad(1) + real(ptcl_comp)*h
                    grad(2)   = grad(2) + real(ptcl_comp)*k
                end do
            end do
        else
            do k = self%lims(2,1),self%lims(2,2)
                arg = real(k,dp) * sh(2)
                ck  = dcos(arg)
                sk  = dsin(arg)
                do h = self%lims(1,1),self%lims(1,2)
                    if( .not. self%resmsk(h,k) ) cycle
                    sh_comp   = cmplx(ck * hcos(h) - sk * hsin(h), ck * hsin(h) + sk * hcos(h),sp)
                    ! retrieve reference component
                    ref_comp  = self%references(h,k,ithr) * self%ctfmats(h,k,i) * self%ref_filt_w(h,k,1)
                    ! shift the particle Fourier component
                    ptcl_comp = self%particles(h,k,i) * sh_comp
                    diff_comp = ref_comp  - ptcl_comp
                    ! update euclidean difference
                    cc(1)     = cc(1) + real(diff_comp * conjg(diff_comp))
                    ! update normalization terms
                    cc(2)     = cc(2) + real( ref_comp * conjg( ref_comp))
                    cc(3)     = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
                    ! update the gradient
                    ptcl_comp = 2 * shconst * imagpart(ptcl_comp * conjg(diff_comp))
                    grad(1)   = grad(1) + real(ptcl_comp)*h
                    grad(2)   = grad(2) + real(ptcl_comp)*k
                end do
            end do
        endif
        corr = 1 - cc(1)/(cc(2)+cc(3))
        grad =   -  grad/(cc(2)+cc(3))
    end subroutine corr_shifted_euclid_2

    ! auto differentiation gives substandard performance to anaytical gradients
    function corr_shifted_ad( self, iptcl, shvec, grad ) result( corr )
        use ADLib_NumParameters_m
        use ADdnSVM_m
        class(cartft_corrcalc), intent(inout)  :: self
        integer,                intent(in)     :: iptcl
        real,                   intent(in)     :: shvec(2)
        real,                   intent(inout)  :: grad(2)
        integer          :: i, h, k, ithr
        complex          :: ref_ptcl, ref_comp, ptcl_comp
        real(sp)         :: cc(3), corr, d0(2)
        logical          :: iseven
        type(dnS_t)      :: sh(2), arg, arg_k, numer_ad, arg_h(self%lims(1,1):self%lims(1,2))
        real(kind=Rkind) :: shconst
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
        sh(1)    = Variable( Val=real(shvec(1), kind=Rkind), nVar=2, iVar=1, nderiv=1 )
        sh(2)    = Variable( Val=real(shvec(2), kind=Rkind), nVar=2, iVar=2, nderiv=1 )
        sh(1)    = sh(1)*shconst
        sh(2)    = sh(2)*shconst
        cc(:)    = 0.
        numer_ad = ZERO
        do h = self%lims(1,1),self%lims(1,2)
            arg_h(h) = real(h, kind=Rkind) * sh(1)
        enddo
        do k = self%lims(2,1),self%lims(2,2)
            arg_k = real(k, kind=Rkind) * sh(2)
            do h = self%lims(1,1),self%lims(1,2)
                if( .not. self%resmsk(h,k) ) cycle
                ! shift component
                arg       = arg_k + arg_h(h)
                ! retrieve reference component
                ref_comp  = self%references(h,k,ithr) * self%ctfmats(h,k,i)
                ! shift the particle Fourier component
                ptcl_comp = self%particles(h,k,i)
                ref_ptcl  = ref_comp*conjg(ptcl_comp)
                ! update cross product
                numer_ad  = numer_ad + cos(arg)*real(realpart(ref_ptcl), kind=Rkind) &
                                    &+ sin(arg)*real(imagpart(ref_ptcl), kind=Rkind)
                ! update normalization terms
                cc(2)     = cc(2) + real(ref_comp  * conjg(ref_comp))
                cc(3)     = cc(3) + real(ptcl_comp * conjg(ptcl_comp))
            end do
        end do
        d0      = get_d0(numer_ad)
        cc(1)   = d0(1)
        grad    = get_d1(numer_ad)
        corr    = norm_corr(  cc(1),cc(2),cc(3))
        grad(1) = norm_corr(grad(1),cc(2),cc(3))
        grad(2) = norm_corr(grad(2),cc(2),cc(3))
    end function corr_shifted_ad

    ! computing the chance/probability that the particle's orientation is correct (prev_corrs are given)
    function ori_chance_1( self, iptcl, o, prev_oris, prev_corrs, R, n ) result(prob)
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: n               ! current iteration number
        integer,                intent(in)    :: iptcl
        type(ori),              intent(in)    :: o               ! current orientation
        type(ori),              intent(in)    :: prev_oris(n)    ! previous iterations' oris
        real,                   intent(in)    :: prev_corrs(n)   ! previous iterations' corrs, corresponding to the previous references
        real,                   intent(in)    :: R               ! a large number
        real                                  :: prob            ! chance/probability
        integer :: i
        real    :: o_euls(3), diff_euls(3)
        prob   = 0.
        o_euls = o%get_euler()
        do i = 1, n
            diff_euls = o_euls - prev_oris(i)%get_euler()
            prob      = prob + product(sin(R*diff_euls)/diff_euls)*prev_corrs(i)
        enddo
        prob = prob/(n*pi**3)
    end function ori_chance_1

    ! computing the chance/probability that the particle's orientation is correct (prev_corrs are recomputed with the new reference)
    function ori_chance_2( self, iptcl, o, prev_oris, R, n ) result(prob)
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: n               ! current iteration number
        integer,                intent(in)    :: iptcl
        type(ori),              intent(in)    :: o               ! current orientation
        type(ori),              intent(in)    :: prev_oris(n)    ! previous iterations' oris
        real,                   intent(in)    :: R               ! a large number
        real                                  :: prob            ! chance/probability
        integer :: i
        real    :: o_euls(3), diff_euls(3)
        prob   = 0.
        o_euls = o%get_euler()
        do i = 1, n
            diff_euls = o_euls - prev_oris(i)%get_euler()
            prob      = prob + product(sin(R*diff_euls)/diff_euls)*self%project_and_correlate(iptcl, prev_oris(i))
        enddo
        prob = prob/(n*pi**3)
    end function ori_chance_2

    subroutine assign_sigma2_noise( self, sigma2_noise )
        class(cartft_corrcalc),       intent(inout) :: self
        real,    allocatable, target, intent(inout) :: sigma2_noise(:,:)
        self%sigma2_noise => sigma2_noise
    end subroutine assign_sigma2_noise

    subroutine calc_sigma_contrib( self, iptcl, o, shvec, sigma_contrib)
        class(cartft_corrcalc), intent(inout) :: self
        integer,                intent(in)    :: iptcl
        type(ori),              intent(in)    :: o
        real(sp),               intent(in)    :: shvec(2)
        real(sp),               intent(out)   :: sigma_contrib(params_glob%kfromto(1):params_glob%kfromto(2))
        ! TODO
        sigma_contrib = 1.0
    end subroutine calc_sigma_contrib

    ! DESTRUCTOR

    subroutine kill( self )
        class(cartft_corrcalc), intent(inout) :: self
        integer :: i, ithr
        if( self%existence )then
            nullify(self%vol_even,self%vol_odd,self%sigma2_noise,cftcc_glob)
            if( allocated(self%pinds)        ) deallocate(self%pinds)
            if( allocated(self%pxls_p_shell) ) deallocate(self%pxls_p_shell)
            if( allocated(self%ctfmats)      ) deallocate(self%ctfmats)
            if( allocated(self%ref_filt_w)   ) deallocate(self%ref_filt_w)
            if( allocated(self%particles)    ) deallocate(self%particles)
            if( allocated(self%references)   ) deallocate(self%references)
            if( allocated(self%resmsk)       ) deallocate(self%resmsk)
            if( allocated(self%iseven)       ) deallocate(self%iseven)
            self%existence  = .false.
        endif
    end subroutine kill

end module simple_cartft_corrcalc

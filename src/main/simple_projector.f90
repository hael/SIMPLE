! projection of 3D volumes in the Fourier domain by convolution interpolation
! to generate band-pass limited Cartesian and polar 2D Fourier transforms
module simple_projector
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'

use simple_kbinterpol, only: kbinterpol
use simple_image,      only: image
use simple_ori,        only: ori
use simple_ori_light,  only: ori_light
use simple_oris,       only: oris
implicit none

public :: projector
private
#include "simple_local_flags.inc"

logical, parameter :: MEMOIZEKB     = .false.
integer, parameter :: NKBPOINTS     = 10
logical, parameter :: USE_WEISZFELD = .false.
logical, parameter :: USE_HYPERB    = .false.


type, extends(image) :: projector
    private
    type(kbinterpol)      :: kbwin                   !< window function object
    integer               :: ldim_exp(3,2) = 0       !< expanded FT matrix limits
    complex, allocatable  :: cmat_exp(:,:,:)         !< expanded FT matrix
    real,    allocatable  :: polweights_mat(:,:,:)   !< polar weights matrix for the image to polar transformer
    integer, allocatable  :: polcyc1_mat(:,:,:)      !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)      !< image cyclic adresses for the image to polar transformer
    logical, allocatable  :: is_in_mask(:,:,:)       !< neighbour matrix for the shape mask projector
    integer               :: wdim      = 0              !< dimension of K-B window
    integer               :: wdim_cubd = 0              !< dimension of K-B window
    integer               :: iwinsz    = 0              !< integer half-window size
    logical               :: expanded_exists=.false. !< indicates FT matrix existence
  contains
    ! CONSTRUCTORS
    procedure          :: expand_cmat
    ! SETTERS
    procedure          :: reset_expanded
    ! FOURIER PROJECTORS
    procedure          :: fproject
    procedure          :: fproject_serial
    procedure          :: fproject_polar
    procedure, private :: fproject_polar_memo
    procedure, private :: fproject_polar_weiszfeld
    procedure          :: dfproject_polar
    procedure          :: fdf_project_polar
    procedure, private :: fdf_project_polar_memo
    procedure          :: interp_fcomp
    procedure          :: interp_fcomp_hyperb
    procedure          :: interp_fcomp_weiszfeld
    procedure          :: interp_fcomp_weiszfeld_dp
    procedure          :: interp_fcomp_trilinear
    procedure, private :: interp_fcomp_memo
    procedure, private :: dinterp_fcomp
    procedure, private :: fdf_interp_fcomp
    procedure, private :: fdf_interp_fcomp_memo
    ! DESTRUCTOR
    procedure          :: kill_expanded
end type projector

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor of the expanded Fourier matrix
    subroutine expand_cmat( self, alpha )
        class(projector), intent(inout) :: self
        real,             intent(in)    :: alpha !< oversampling factor
        integer, allocatable :: cyck(:), cycm(:), cych(:)
        integer :: h, k, m, phys(3), logi(3), lims(3,2), ldim(3)
        !real    :: winsz
        call self%kill_expanded
        ldim = self%get_ldim()
        if( .not.self%is_ft() ) THROW_HARD('volume needs to be FTed before call; expand_cmat')
        if( ldim(3) == 1      ) THROW_HARD('only for volumes; expand_cmat')
        self%kbwin         = kbinterpol(KBWINSZ, alpha)
        self%iwinsz        = ceiling(self%kbwin%get_winsz() - 0.5)
        if( MEMOIZEKB ) call self%kbwin%memoize(NKBPOINTS)
        self%wdim          = self%kbwin%get_wdim()
        self%wdim_cubd     = self%wdim**3
        lims               = self%loop_lims(3)
        self%ldim_exp(:,2) = maxval(abs(lims)) + ceiling(self%kbwin%get_winsz())
        self%ldim_exp(:,1) = -self%ldim_exp(:,2)
        allocate( self%cmat_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),&
                                &self%ldim_exp(2,1):self%ldim_exp(2,2),&
                                &self%ldim_exp(3,1):self%ldim_exp(3,2)),&
                                &cych(self%ldim_exp(1,1):self%ldim_exp(1,2)),&
                                &cyck(self%ldim_exp(2,1):self%ldim_exp(2,2)),&
                                &cycm(self%ldim_exp(3,1):self%ldim_exp(3,2)), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: expand_cmat; simple_projector",alloc_stat)
        ! pre-compute addresses in 2nd and 3rd dimension
        do h = self%ldim_exp(1,1),self%ldim_exp(1,2)
            cych(h) = h
            if( h < lims(1,1) .or. h > lims(1,2) ) cych(h) = cyci_1d( lims(1,:),h )
        enddo
        do k = self%ldim_exp(2,1),self%ldim_exp(2,2)
            cyck(k) = k
            if( k < lims(2,1) .or. k > lims(2,2) ) cyck(k) = cyci_1d( lims(2,:),k )
        enddo
        do m = self%ldim_exp(3,1),self%ldim_exp(3,2)
            cycm(m) = m
            if( m < lims(3,1) .or. m > lims(3,2) ) cycm(m) = cyci_1d( lims(3,:),m )
        enddo
        ! build expanded fourier components matrix
        self%cmat_exp = cmplx(0.,0.)
        !$omp parallel do collapse(3) schedule(static) default(shared)&
        !$omp private(h,k,m,logi,phys) proc_bind(close)
        do h = self%ldim_exp(1,1),self%ldim_exp(1,2)
            do k = self%ldim_exp(2,1),self%ldim_exp(2,2)
                do m = self%ldim_exp(3,1),self%ldim_exp(3,2)
                    logi = [cych(h),cyck(k),cycm(m)]
                    phys = self%comp_addr_phys(logi)
                    self%cmat_exp(h,k,m) = self%get_fcomp(logi, phys)
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate( cych,cyck,cycm )
        self%expanded_exists = .true.
    end subroutine expand_cmat

    ! SETTERS

    !>  \brief  sets the expanded Fourier matrix to zero
    subroutine reset_expanded( self )
        class(projector), intent(inout) :: self
        if( .not.self%expanded_exists )&
            &THROW_HARD('expanded Fourier matrix does not exist; reset_expanded')
        self%cmat_exp = cmplx(0.,0.)
    end subroutine reset_expanded

    ! FOURIER PROJECTORS

    !> \brief  extracts a Fourier plane from the expanded FT matrix of a volume (self)
    subroutine fproject( self, e, fplane )
        class(projector), intent(inout) :: self
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real        :: loc(3), e_rotmat(3,3)
        integer     :: h, k, sqarg, sqlp, lims(3,2), phys(3), ldim(3)
        complex(sp) :: comp
        lims   = self%loop_lims(2)
        sqlp   = (maxval(lims(:,2)))**2
        ldim   = fplane%get_ldim()
        call fplane%zero_and_flag_ft()
        e_rotmat = e%get_mat()
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(h,k,sqarg,loc,phys,comp) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                loc  = matmul(real([h,k,0]), e_rotmat)
                if (USE_WEISZFELD) then
                    comp = self%interp_fcomp_weiszfeld(loc)
                else if (USE_HYPERB) then
                    comp = self%interp_fcomp_hyperb(loc)
                else
                    comp = self%interp_fcomp(loc)
                end if
                if (h > 0) then
                    phys(1) = h + 1
                    phys(2) = k + 1 + MERGE(ldim(2),0,k < 0)
                    phys(3) = 1
                    call fplane%set_cmat_at(phys, comp)
                else
                    phys(1) = -h + 1
                    phys(2) = -k + 1 + MERGE(ldim(2),0,-k < 0)
                    phys(3) = 1
                    call fplane%set_cmat_at(phys, conjg(comp))
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine fproject

    !> \brief  extracts a Fourier plane from the expanded FT matrix of a volume (self)
    subroutine fproject_serial( self, e, fplane )
        class(projector), intent(inout) :: self
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real        :: loc(3), e_rotmat(3,3)
        integer     :: h, k, sqarg, sqlp, lims(3,2), phys(3), ldim(3)
        complex(sp) :: comp
        lims   = self%loop_lims(2)
        sqlp   = (maxval(lims(:,2)))**2
        ldim   = fplane%get_ldim()
        call fplane%zero_and_flag_ft()
        e_rotmat = e%get_mat()
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                loc  = matmul(real([h,k,0]), e_rotmat)
                if (USE_WEISZFELD) then
                    comp = self%interp_fcomp_weiszfeld(loc)
                else if (USE_HYPERB) then
                    comp = self%interp_fcomp_hyperb(loc)
                else
                    comp = self%interp_fcomp(loc)
                end if
                if (h > 0) then
                    phys(1) = h + 1
                    phys(2) = k + 1 + MERGE(ldim(2),0,k < 0)
                    phys(3) = 1
                    call fplane%set_cmat_at(phys, comp)
                else
                    phys(1) = -h + 1
                    phys(2) = -k + 1 + MERGE(ldim(2),0,-k < 0)
                    phys(3) = 1
                    call fplane%set_cmat_at(phys, conjg(comp))
                endif
            end do
        end do
    end subroutine fproject_serial

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fproject_polar( self, iref, e, pftcc, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self   !< projector object
        integer,                 intent(in)    :: iref   !< which reference
        class(ori),              intent(in)    :: e      !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc  !< object that holds the polar image
        logical,                 intent(in)    :: iseven !< eo flag
        integer :: irot, k, pdim(3)
        real    :: vec(3), loc(3), e_rotmat(3,3)
        if( MEMOIZEKB )then
            call self%fproject_polar_memo(iref, e, pftcc, iseven)
            return
        endif
        if ( USE_WEISZFELD )then
            call self%fproject_polar_weiszfeld(iref, e, pftcc, iseven)
            return
        end if
        pdim = pftcc%get_pdim()
        e_rotmat = e%get_mat()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = pftcc%get_coord(irot,k)
                vec(3)  = 0.
                loc     = matmul(vec,e_rotmat)
                call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp(loc), iseven)
            end do
        end do
    end subroutine fproject_polar

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fproject_polar_memo( self, iref, e, pftcc, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self   !< projector object
        integer,                 intent(in)    :: iref   !< which reference
        class(ori),              intent(in)    :: e      !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc  !< object that holds the polar image
        logical,                 intent(in)    :: iseven !< eo flag
        integer :: irot, k, pdim(3)
        real    :: vec(3), loc(3), e_rotmat(3,3)
        pdim = pftcc%get_pdim()
        e_rotmat = e%get_mat()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = pftcc%get_coord(irot,k)
                vec(3)  = 0.
                loc     = matmul(vec,e_rotmat)
                call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp_memo(loc), iseven)
            end do
        end do
    end subroutine fproject_polar_memo

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fproject_polar_weiszfeld( self, iref, e, pftcc, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self   !< projector object
        integer,                 intent(in)    :: iref   !< which reference
        class(ori),              intent(in)    :: e      !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc  !< object that holds the polar image
        logical,                 intent(in)    :: iseven !< eo flag
        integer :: irot, k, pdim(3)
        real    :: vec(3), loc(3), e_rotmat(3,3)
        pdim = pftcc%get_pdim()
        e_rotmat = e%get_mat()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = pftcc%get_coord(irot,k)
                vec(3)  = 0.
                loc     = matmul(vec,e_rotmat)
                call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp_weiszfeld(loc), iseven)
            end do
        end do
    end subroutine fproject_polar_weiszfeld

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine dfproject_polar( self, iref, euls, pftcc, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self    !< projector object
        integer,                 intent(in)    :: iref    !< which reference
        real(dp),                intent(inout) :: euls(3) !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc   !< object that holds the polar image
        logical,                 intent(in)    :: iseven  !< eo flag
        integer         :: irot, k, pdim(3)
        real(dp)        :: vec(3), loc(3)
        complex         :: dfcomp(3)
        type(ori_light) :: or
        pdim = pftcc%get_pdim()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2)            = real(pftcc%get_coord(irot,k),kind=dp)
                vec(3)             = 0._dp
                loc                = matmul(vec,or%euler2m(euls))
                dfcomp             = self%dinterp_fcomp(euls, loc, vec )
                call pftcc%set_dref_fcomp( iref, irot, k, dfcomp, iseven )
            end do
        end do
    end subroutine dfproject_polar

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fdf_project_polar( self, iref, euls, pftcc, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self    !< projector object
        integer,                 intent(in)    :: iref    !< which reference
        real(dp),                intent(in)    :: euls(3) !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc   !< object that holds the polar image
        logical,                 intent(in)    :: iseven  !< eo flag
        integer         :: irot, k, pdim(3)
        real(dp)        :: vec(3), loc(3)
        complex         :: fcomp, dfcomp(3)
        type(ori_light) :: or
        if( MEMOIZEKB )then
            call self%fdf_project_polar_memo(iref, euls, pftcc, iseven)
            return
        endif
        pdim = pftcc%get_pdim()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = real(pftcc%get_coord(irot,k),kind=dp)
                vec(3)  = 0._dp
                loc     = matmul(vec,or%euler2m(euls))
                dfcomp  = self%fdf_interp_fcomp(euls, loc, vec, fcomp)
                call pftcc%set_ref_fcomp( iref, irot, k, fcomp, iseven )
                call pftcc%set_dref_fcomp( iref, irot, k, dfcomp, iseven )
            end do
        end do
    end subroutine fdf_project_polar

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fdf_project_polar_memo( self, iref, euls, pftcc, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self    !< projector object
        integer,                 intent(in)    :: iref    !< which reference
        real(dp),                intent(in)    :: euls(3) !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc   !< object that holds the polar image
        logical,                 intent(in)    :: iseven  !< eo flag
        integer         :: irot, k, pdim(3)
        real(dp)        :: vec(3), loc(3)
        complex         :: fcomp, dfcomp(3)
        type(ori_light) :: or
        pdim = pftcc%get_pdim()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = real(pftcc%get_coord(irot,k),kind=dp)
                vec(3)  = 0._dp
                loc     = matmul(vec,or%euler2m(euls))
                dfcomp  = self%fdf_interp_fcomp_memo(euls, loc, vec, fcomp)
                call pftcc%set_ref_fcomp( iref, irot, k, fcomp, iseven )
                call pftcc%set_dref_fcomp( iref, irot, k, dfcomp, iseven )
            end do
        end do
    end subroutine fdf_project_polar_memo

    !>  \brief is to interpolate from the expanded complex matrix
    pure function interp_fcomp( self, loc )result( comp )
        class(projector), intent(in) :: self
        real,             intent(in) :: loc(3)
        complex :: comp
        real    :: w(1:self%wdim,1:self%wdim,1:self%wdim)
        integer :: i, win(2,3) ! window boundary array in fortran contiguous format
        ! interpolation kernel window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        ! interpolation kernel matrix
        w = 1.
        do i=1,self%wdim
            w(i,:,:) = w(i,:,:) * self%kbwin%apod( real(win(1,1)+i-1)-loc(1) )
            w(:,i,:) = w(:,i,:) * self%kbwin%apod( real(win(1,2)+i-1)-loc(2) )
            w(:,:,i) = w(:,:,i) * self%kbwin%apod( real(win(1,3)+i-1)-loc(3) )
        end do
        ! SUM( kernel x components )
        comp = sum( w * self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2),win(1,3):win(2,3)) ) / sum( w )
    end function interp_fcomp

    function interp_fcomp_hyperb( self, loc )result( comp )
      class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        complex :: comp
        integer, parameter :: maxits     = 30
        real,    parameter :: relTol     = 1e-5
        real,    parameter :: min_dists  = 1e-8
        real,    parameter :: min_denom  = 5e-8
        real,    parameter :: min_denom2 = 1e-11
        integer :: i, j, k, win(2,3)            ! window boundary array in fortran contiguous format
        complex :: pts(self%wdim_cubd)          ! points contrib to geometric median
        real    :: dists(self%wdim_cubd)        ! distances
        real    :: rec_dists(self%wdim_cubd)    ! reciprocal distances
        real    :: denom(self%wdim_cubd)        ! denominator
        complex :: zi,zii                       ! iterate
        logical :: it_condition
        real    :: eps
        integer :: idx, cnt
        complex :: R
        real    :: gamma, eta, rr
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        ! setup problem
        do i=1,self%wdim
            do j = 1,self%wdim
                do k = 1,self%wdim
                    idx            = (i-1)*self%wdim**2 + (j-1)*self%wdim + k
                    pts(idx)       = self%cmat_exp(win(1,1)+i-1, win(1,2)+i-1, win(1,3)+i-1)
                    dists(idx)     = ( real(win(1,1)+i-1)-loc(1) )**2 + &
                        &( real(win(1,2)+i-1)-loc(2) )**2 + ( real(win(1,3)+i-1)-loc(3) )**2
                    rec_dists(idx) = 1. / sqrt( dists(idx) )
                end do
            end do
        end do
        ! if any distance(s) extremely small (weight approx. 1), then use that value(s)
        if (any(dists(:) .le. min_dists)) then
            comp = sum(pts, mask=(dists(:) .le. min_dists)) / real(count(dists(:) .le. min_dists))
            return
        end if
        rec_dists = rec_dists / sum(rec_dists(:))
        ! initialize with geometric mean
        comp = sum(rec_dists(:)*pts(:))
    end function interp_fcomp_hyperb

    function interp_fcomp_weiszfeld( self, loc )result( comp )
      class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        complex :: comp
        integer, parameter :: maxits     = 30
        real,    parameter :: relTol     = 1e-5
        real,    parameter :: min_dists  = 1e-8
        real,    parameter :: min_denom  = 5e-8
        real,    parameter :: min_denom2 = 1e-11
        integer :: i, j, k, win(2,3)            ! window boundary array in fortran contiguous format
        complex :: pts(self%wdim_cubd)          ! points contrib to geometric median
        real    :: dists(self%wdim_cubd)        ! distances
        real    :: rec_dists(self%wdim_cubd)    ! reciprocal distances
        real    :: denom(self%wdim_cubd)        ! denominator
        complex :: zi,zii                       ! iterate
        logical :: it_condition
        real    :: eps
        integer :: idx, cnt
        complex :: R
        real    :: gamma, eta, rr
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        ! setup problem
        do i=1,self%wdim
            do j = 1,self%wdim
                do k = 1,self%wdim
                    idx            = (i-1)*self%wdim**2 + (j-1)*self%wdim + k
                    pts(idx)       = self%cmat_exp(win(1,1)+i-1, win(1,2)+i-1, win(1,3)+i-1)
                    dists(idx)     = ( real(win(1,1)+i-1)-loc(1) )**2 + &
                        &( real(win(1,2)+i-1)-loc(2) )**2 + ( real(win(1,3)+i-1)-loc(3) )**2
                    rec_dists(idx) = 1. / sqrt( dists(idx) )
                end do
            end do
        end do
        ! if any distance(s) extremely small (weight approx. 1), then use that value(s)
        if (any(dists(:) .le. min_dists)) then
            comp = sum(pts, mask=(dists(:) .le. min_dists)) / real(count(dists(:) .le. min_dists))
            return
        end if
        rec_dists = rec_dists / sum(rec_dists(:))
        ! initialize with geometric mean
        zi = sum(rec_dists(:)*pts(:))
        it_condition = .true.
        cnt = 0
        do while (it_condition)
            denom = abs(pts(:)-zi)
            if (any(denom(:) .le. min_denom)) then
                ! if denominator too small (see Vardi, Zhang)
                eta   = sum(rec_dists(:), mask=(denom(:) .le. min_denom))
                zii   = sum(rec_dists(:)*pts(:)/denom(:), mask=(denom(:) .gt. min_denom))/&
                    sum(rec_dists(:)/denom(:), mask=(denom(:) .gt. min_denom))
                R     = sum(rec_dists(:)*(pts(:)-zi)/denom(:),mask=(denom(:) .gt. min_denom))
                rr    = abs(R)
                if ((eta .le. min_denom2).and.(rr .lt. min_denom2)) then
                    eta = 0.
                    rr  = 1.
                end if
                gamma = min(1., eta/rr)
                zii   = (1. - eta/rr) * zii + gamma * zi
            else
                zii = sum(rec_dists(:)*pts(:)/denom(:))/sum(rec_dists(:)/denom(:))
            end if
            eps = abs(zii-zi)
            zi  = zii
            cnt = cnt + 1
            it_condition = (eps > relTol).and.(cnt .lt. maxits)
        end do
        if ((cnt .gt. maxits).or.(is_nan(real(zi))).or.(is_nan(aimag(zi)))) then
            ! if convergence has not been achieved, try double precision algorithm
            comp = self%interp_fcomp_weiszfeld_dp( loc )
            if ((is_nan(real(comp))).or.(is_nan(aimag(comp)))) then
                ! if still not converged (due to NaN), then use KB kernel
                comp = self%interp_fcomp( loc )
            end if
        else
            comp = zi
        end if
    end function interp_fcomp_weiszfeld

    function interp_fcomp_weiszfeld_dp( self, loc ) result( comp )
      class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        complex :: comp
        integer,  parameter :: maxits     = 30
        real(dp), parameter :: relTol     = 1e-5_dp
        real(dp), parameter :: min_dists  = 1e-8_dp
        real(dp), parameter :: min_denom  = 5e-8_dp
        real(dp), parameter :: min_denom2 = 1e-11_dp
        integer     :: i, j, k, win(2,3)            ! window boundary array in fortran contiguous format
        complex     :: pts(self%wdim_cubd)          ! points contrib to geometric median
        real(dp)    :: dists(self%wdim_cubd)        ! distances
        real(dp)    :: rec_dists(self%wdim_cubd)    ! reciprocal distances
        real(dp)    :: denom(self%wdim_cubd)        ! denominator
        complex(dp) :: zi,zii                       ! iterate
        logical     :: it_condition
        real(dp)    :: eps
        integer     :: idx, cnt
        complex(dp) :: R
        real(dp)    :: gamma, eta, rr
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        ! setup problem
        do i=1,self%wdim
            do j = 1,self%wdim
                do k = 1,self%wdim
                    idx            = (i-1)*self%wdim**2 + (j-1)*self%wdim + k
                    pts(idx)       = self%cmat_exp(win(1,1)+i-1, win(1,2)+i-1, win(1,3)+i-1)
                    dists(idx)     = ( real(win(1,1)+i-1,dp)-real(loc(1),dp) )**2 + &
                        &( real(win(1,2)+i-1,dp)-real(loc(2),dp) )**2 + ( real(win(1,3)+i-1,dp)-real(loc(3),dp) )**2
                    rec_dists(idx) = 1._dp / sqrt( dists(idx) )
                end do
            end do
        end do
        ! if any distance(s) extremely small (weight approx. 1), then use that value(s)
        if (any(dists(:) .le. min_dists)) then
            comp = sum(pts, mask=(dists(:) .le. min_dists)) / real(count(dists(:) .le. min_dists))
            return
        end if
        rec_dists = rec_dists / sum(rec_dists(:))
        ! initialize with geometric mean
        zi = sum(rec_dists(:)*pts(:))
        it_condition = .true.
        cnt = 0
        do while (it_condition)
            denom = abs(pts(:)-zi)
            if (any(denom(:) .le. min_denom)) then
                ! if denominator too small (see Vardi, Zhang)
                eta   = sum(rec_dists(:), mask=(denom(:) .le. min_denom))
                zii   = sum(rec_dists(:)*pts(:)/denom(:), mask=(denom(:) .gt. min_denom))/&
                    sum(rec_dists(:)/denom(:), mask=(denom(:) .gt. min_denom))
                R     = sum(rec_dists(:)*(pts(:)-zi)/denom(:),mask=(denom(:) .gt. min_denom))
                rr    = abs(R)
                if ((eta .le. min_denom2).and.(rr .lt. min_denom2)) then
                    eta = 0._dp
                    rr  = 1._dp
                end if
                gamma = min(1._dp, eta/rr)
                zii   = (1._dp - eta/rr) * zii + gamma * zi
            else
                zii = sum(rec_dists(:)*pts(:)/denom(:))/sum(rec_dists(:)/denom(:))
            end if
            eps = abs(zii-zi)
            zi  = zii
            cnt     = cnt + 1
            it_condition = (eps > relTol).and.(cnt .lt. maxits)
        end do
        if ((cnt .gt. maxits).or.(is_nan(real(zi))).or.(is_nan(aimag(zi)))) then
            ! if not converging or results in NaN, fall back to kernel interpolation
            comp = self%interp_fcomp(loc)
        else
            comp = cmplx(zi,kind=sp)
        end if
    end function interp_fcomp_weiszfeld_dp

    !>  \brief is for tri-linear interpolation from the expanded complex matrix
    pure function interp_fcomp_trilinear( self, loc )result( comp )
        class(projector), intent(in) :: self
        real,             intent(in) :: loc(3)
        complex :: comp
        real    :: w(2,2,2), d(3), dp(3)
        integer :: lb(3)
        lb = floor(loc)
        d  = loc - real(lb)
        dp = 1. - d
        w(1,1,1) = product(dp)
        w(2,1,1) =  d(1) * dp(2) * dp(3)
        w(1,2,1) = dp(1) *  d(2) * dp(3)
        w(1,1,2) = dp(1) * dp(2) *  d(3)
        w(2,1,2) =  d(1) * dp(2) *  d(3)
        w(1,2,2) = dp(1) *  d(2) *  d(3)
        w(2,2,1) =  d(1) *  d(2) * dp(3)
        w(2,2,2) = product(d)
        comp = sum(w * self%cmat_exp(lb(1):lb(1)+1,lb(2):lb(2)+1,lb(3):lb(3)+1))
    end function interp_fcomp_trilinear

    !>  \brief is to interpolate from the expanded complex matrix
    function interp_fcomp_memo( self, loc )result( comp )
        class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        complex :: comp
        real    :: w(1:self%wdim,1:self%wdim,1:self%wdim)
        integer :: i, win(2,3) ! window boundary array in fortran contiguous format
        ! interpolation kernel window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        ! interpolation kernel matrix
        w = 1.
        do i=1,self%wdim
            w(i,:,:) = w(i,:,:) * self%kbwin%apod_memo( real(win(1,1)+i-1)-loc(1) )
            w(:,i,:) = w(:,i,:) * self%kbwin%apod_memo( real(win(1,2)+i-1)-loc(2) )
            w(:,:,i) = w(:,:,i) * self%kbwin%apod_memo( real(win(1,3)+i-1)-loc(3) )
        end do
        ! SUM( kernel x components )
        comp = sum( w * self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2),win(1,3):win(2,3)) ) / sum( w )
    end function interp_fcomp_memo

    !>  \brief is to compute the derivative of the interpolate from the expanded complex matrix
    !! \param loc 3-dimensional location in the volume
    !! \param q 2-dimensional location on the plane (h,k,0)
    !! \param e orientation class
    function dinterp_fcomp( self, euls, loc, q ) result(res)
        class(projector), intent(inout)                          :: self
        real(dp),         intent(in)                             :: euls(3)
        real(dp),         intent(in)                             :: loc(3)
        real(dp),         intent(in)                             :: q(3)
        complex(sp)                                              :: res(3)
        real(dp)                                                 :: drotmat(3,3,3)
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: w1, w2, w3, w
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_t, dw2_t, dw3_t    !theta
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_p, dw2_p, dw3_p    !psi
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_ph, dw2_ph, dw3_ph !phi
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: wt, wp, wph
        real(dp)                                                 :: dRdangle(3,3)
        real(dp)                                                 :: dapod_tmp(3)
        complex(dp)                                              :: N     !numerator
        real(dp)                                                 :: D, D2 !denominator
        integer                                                  :: i, win(2,3) ! window boundary array in fortran contiguous format
        type(ori_light)                                          :: or
        ! interpolation kernel window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        call or%euler2dm(euls, drotmat)
        dRdangle(:,1) = matmul(q, drotmat(:,:,1))
        dRdangle(:,2) = matmul(q, drotmat(:,:,2))
        dRdangle(:,3) = matmul(q, drotmat(:,:,3))
        w1  = 1.0_dp ; w2  = 1.0_dp ; w3  = 1.0_dp
        ! theta
        dw1_t = 1.0_dp ; dw2_t = 1.0_dp ; dw3_t = 1.0_dp
        ! psi
        dw1_p = 1.0_dp ; dw2_p = 1.0_dp ; dw3_p = 1.0_dp
        ! phi
        dw1_ph = 1.0_dp ; dw2_ph = 1.0_dp ; dw3_ph = 1.0_dp
        do i=1,self%wdim
            w1(i,:,:) = self%kbwin%apod_dp( real( win(1,1)+i-1,kind=dp) - loc(1) )
            w2(:,i,:) = self%kbwin%apod_dp( real( win(1,2)+i-1,kind=dp) - loc(2) )
            w3(:,:,i) = self%kbwin%apod_dp( real( win(1,3)+i-1,kind=dp) - loc(3) )
            dapod_tmp(1) = self%kbwin%dapod( real( win(1,1)+i-1,kind=dp) - loc(1) )
            dapod_tmp(2) = self%kbwin%dapod( real( win(1,2)+i-1,kind=dp) - loc(2) )
            dapod_tmp(3) = self%kbwin%dapod( real( win(1,3)+i-1,kind=dp) - loc(3) )
            dw1_t (i,:,:) = - dapod_tmp(1) * dRdangle(1,1)
            dw2_t (:,i,:) = - dapod_tmp(2) * dRdangle(2,1)
            dw3_t (:,:,i) = - dapod_tmp(3) * dRdangle(3,1)
            dw1_p (i,:,:) = - dapod_tmp(1) * dRdangle(1,2)
            dw2_p (:,i,:) = - dapod_tmp(2) * dRdangle(2,2)
            dw3_p (:,:,i) = - dapod_tmp(3) * dRdangle(3,2)
            dw1_ph(i,:,:) = - dapod_tmp(1) * dRdangle(1,3)
            dw2_ph(:,i,:) = - dapod_tmp(2) * dRdangle(2,3)
            dw3_ph(:,:,i) = - dapod_tmp(3) * dRdangle(3,3)
        end do
        wt  = dw1_t  * w2 * w3  + w1 * dw2_t  * w3  + w1 * w2 * dw3_t
        wp  = dw1_p  * w2 * w3  + w1 * dw2_p  * w3  + w1 * w2 * dw3_p
        wph = dw1_ph * w2 * w3  + w1 * dw2_ph * w3  + w1 * w2 * dw3_ph
        w   = w1 * w2 * w3
        N   = sum( w  *  self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) )
        D   = sum( w )
        D2   = D**2
        res(1) = cmplx( ( sum( wt * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum(wt) ) / D2 )
        res(2) = cmplx( ( sum( wp * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum(wp) ) / D2 )
        res(3) = cmplx( ( sum( wph * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N *sum(wph) ) / D2 )
    end function dinterp_fcomp

    !>  \brief is to compute the derivative of the interpolate from the expanded complex matrix
    !! \param loc 3-dimensional location in the volume
    !! \param q 2-dimensional location on the plane (h,k,0)
    !! \param e orientation class
    function fdf_interp_fcomp( self, euls, loc, q, fcomp ) result(res)
        class(projector), intent(inout)                          :: self
        real(dp),         intent(in)                             :: euls(3)
        real(dp),         intent(in)                             :: loc(3)
        real(dp),         intent(in)                             :: q(3)
        complex(sp),      intent(out)                            :: fcomp
        complex(sp)                                              :: res(3)
        real(dp)                                                 :: drotmat(3,3,3)
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: w1, w2, w3
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_t, dw2_t, dw3_t    !theta
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_p, dw2_p, dw3_p    !psi
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_ph, dw2_ph, dw3_ph !phi
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: wt, wp, wph, w
        real(dp)                                                 :: dRdangle(3,3)
        real(dp)                                                 :: dapod_tmp(3)
        complex(dp)                                              :: N     !numerator
        real(dp)                                                 :: D, D2 !denominator
        integer                                                  :: i, win(2,3) ! window boundary array in fortran contiguous format
        type(ori_light)                                          :: or
        ! interpolation kernel window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        call or%euler2dm(euls, drotmat)
        dRdangle(:,1) = matmul(q, drotmat(:,:,1))
        dRdangle(:,2) = matmul(q, drotmat(:,:,2))
        dRdangle(:,3) = matmul(q, drotmat(:,:,3))
        w1  = 1.0_dp ; w2  = 1.0_dp ; w3  = 1.0_dp
        ! theta
        dw1_t = 1.0_dp ; dw2_t = 1.0_dp ; dw3_t = 1.0_dp
        ! psi
        dw1_p = 1.0_dp ; dw2_p = 1.0_dp ; dw3_p = 1.0_dp
        ! phi
        dw1_ph = 1.0_dp ; dw2_ph = 1.0_dp ; dw3_ph = 1.0_dp
        do i=1,self%wdim
            w1(i,:,:) = self%kbwin%apod_dp( real( win(1,1)+i-1,kind=dp) - loc(1) )
            w2(:,i,:) = self%kbwin%apod_dp( real( win(1,2)+i-1,kind=dp) - loc(2) )
            w3(:,:,i) = self%kbwin%apod_dp( real( win(1,3)+i-1,kind=dp) - loc(3) )
            dapod_tmp(1) = self%kbwin%dapod( real( win(1,1)+i-1,kind=dp) - loc(1) )
            dapod_tmp(2) = self%kbwin%dapod( real( win(1,2)+i-1,kind=dp) - loc(2) )
            dapod_tmp(3) = self%kbwin%dapod( real( win(1,3)+i-1,kind=dp) - loc(3) )
            dw1_t (i,:,:) = - dapod_tmp(1) * dRdangle(1,1)
            dw2_t (:,i,:) = - dapod_tmp(2) * dRdangle(2,1)
            dw3_t (:,:,i) = - dapod_tmp(3) * dRdangle(3,1)
            dw1_p (i,:,:) = - dapod_tmp(1) * dRdangle(1,2)
            dw2_p (:,i,:) = - dapod_tmp(2) * dRdangle(2,2)
            dw3_p (:,:,i) = - dapod_tmp(3) * dRdangle(3,2)
            dw1_ph(i,:,:) = - dapod_tmp(1) * dRdangle(1,3)
            dw2_ph(:,i,:) = - dapod_tmp(2) * dRdangle(2,3)
            dw3_ph(:,:,i) = - dapod_tmp(3) * dRdangle(3,3)
        end do
        wt  = dw1_t  * w2 * w3  + w1 * dw2_t  * w3  + w1 * w2 * dw3_t
        wp  = dw1_p  * w2 * w3  + w1 * dw2_p  * w3  + w1 * w2 * dw3_p
        wph = dw1_ph * w2 * w3  + w1 * dw2_ph * w3  + w1 * w2 * dw3_ph
        w   = w1 * w2 * w3
        N   = sum( w  *  self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) )
        D   = sum( w )
        D2   = D**2
        !! FIXME D2 uninitialized
        res(1) = cmplx( ( sum( wt  * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wt  ) ) / D2 )
        res(2) = cmplx( ( sum( wp  * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wp  ) ) / D2 )
        res(3) = cmplx( ( sum( wph * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wph ) ) / D2 )
        fcomp  = cmplx( N / D )
    end function fdf_interp_fcomp

    !>  \brief is to compute the derivative of the interpolate from the expanded complex matrix
    !! \param loc 3-dimensional location in the volume
    !! \param q 2-dimensional location on the plane (h,k,0)
    !! \param e orientation class
    function fdf_interp_fcomp_memo( self, euls, loc, q, fcomp ) result(res)
        class(projector), intent(inout)                          :: self
        real(dp),         intent(in)                             :: euls(3)
        real(dp),         intent(in)                             :: loc(3)
        real(dp),         intent(in)                             :: q(3)
        complex(sp),      intent(out)                            :: fcomp
        complex(sp)                                              :: res(3)
        real(dp)                                                 :: drotmat(3,3,3)
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: w1, w2, w3
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_t, dw2_t, dw3_t    !theta
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_p, dw2_p, dw3_p    !psi
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: dw1_ph, dw2_ph, dw3_ph !phi
        real(dp), dimension(1:self%wdim,1:self%wdim,1:self%wdim) :: wt, wp, wph, w
        real(dp)                                                 :: dRdangle(3,3)
        real(dp)                                                 :: dapod_tmp(3)
        complex(dp)                                              :: N     !numerator
        real(dp)                                                 :: D, D2 !denominator
        integer                                                  :: i, win(2,3) ! window boundary array in fortran contiguous format
        type(ori_light)                                          :: or
        ! interpolation kernel window
        win(1,:) = nint(loc)
        win(2,:) = win(1,:) + self%iwinsz
        win(1,:) = win(1,:) - self%iwinsz
        call or%euler2dm(euls, drotmat)
        dRdangle(:,1) = matmul(q, drotmat(:,:,1))
        dRdangle(:,2) = matmul(q, drotmat(:,:,2))
        dRdangle(:,3) = matmul(q, drotmat(:,:,3))
        w1  = 1.0_dp ; w2  = 1.0_dp ; w3  = 1.0_dp
        ! theta
        dw1_t = 1.0_dp ; dw2_t = 1.0_dp ; dw3_t = 1.0_dp
        ! psi
        dw1_p = 1.0_dp ; dw2_p = 1.0_dp ; dw3_p = 1.0_dp
        ! phi
        dw1_ph = 1.0_dp ; dw2_ph = 1.0_dp ; dw3_ph = 1.0_dp
        do i=1,self%wdim
            w1(i,:,:) = self%kbwin%apod_memo_dp( real( win(1,1)+i-1,kind=dp) - loc(1) )
            w2(:,i,:) = self%kbwin%apod_memo_dp( real( win(1,2)+i-1,kind=dp) - loc(2) )
            w3(:,:,i) = self%kbwin%apod_memo_dp( real( win(1,3)+i-1,kind=dp) - loc(3) )
            dapod_tmp(1) = self%kbwin%dapod_memo( real( win(1,1)+i-1,kind=dp) - loc(1) )
            dapod_tmp(2) = self%kbwin%dapod_memo( real( win(1,2)+i-1,kind=dp) - loc(2) )
            dapod_tmp(3) = self%kbwin%dapod_memo( real( win(1,3)+i-1,kind=dp) - loc(3) )
            dw1_t (i,:,:) = - dapod_tmp(1) * dRdangle(1,1)
            dw2_t (:,i,:) = - dapod_tmp(2) * dRdangle(2,1)
            dw3_t (:,:,i) = - dapod_tmp(3) * dRdangle(3,1)
            dw1_p (i,:,:) = - dapod_tmp(1) * dRdangle(1,2)
            dw2_p (:,i,:) = - dapod_tmp(2) * dRdangle(2,2)
            dw3_p (:,:,i) = - dapod_tmp(3) * dRdangle(3,2)
            dw1_ph(i,:,:) = - dapod_tmp(1) * dRdangle(1,3)
            dw2_ph(:,i,:) = - dapod_tmp(2) * dRdangle(2,3)
            dw3_ph(:,:,i) = - dapod_tmp(3) * dRdangle(3,3)
        end do
        wt  = dw1_t  * w2 * w3  + w1 * dw2_t  * w3  + w1 * w2 * dw3_t
        wp  = dw1_p  * w2 * w3  + w1 * dw2_p  * w3  + w1 * w2 * dw3_p
        wph = dw1_ph * w2 * w3  + w1 * dw2_ph * w3  + w1 * w2 * dw3_ph
        w   = w1 * w2 * w3
        N   = sum( w  *  self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) )
        D   = sum( w )
        D2   = D**2
        !! FIXME D2 uninitialized
        res(1) = cmplx( ( sum( wt*self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum(wt) ) / D2 )
        res(2) = cmplx( ( sum( wp*self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum(wp) ) / D2 )
        res(3) = cmplx( ( sum( wph*self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wph ) ) / D2 )
        fcomp  = cmplx( N / D )
    end function fdf_interp_fcomp_memo

    !>  \brief  is a destructor of expanded matrices (imgpolarizer AND expanded projection of)
    subroutine kill_expanded( self )
        class(projector), intent(inout) :: self !< projector instance
        if( allocated(self%cmat_exp) )deallocate(self%cmat_exp)
        self%ldim_exp        = 0
        self%expanded_exists = .false.
    end subroutine kill_expanded

end module simple_projector

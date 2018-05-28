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

logical, parameter :: MEMOIZEKB = .false.
integer, parameter :: NKBPOINTS = 10

type, extends(image) :: projector
    private
    type(kbinterpol)      :: kbwin                   !< window function object
    integer               :: ldim_exp(3,2) = 0       !< expanded FT matrix limits
    complex, allocatable  :: cmat_exp(:,:,:)         !< expanded FT matrix
    real,    allocatable  :: polweights_mat(:,:,:)   !< polar weights matrix for the image to polar transformer
    integer, allocatable  :: polcyc1_mat(:,:,:)      !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)      !< image cyclic adresses for the image to polar transformer
    logical, allocatable  :: is_in_mask(:,:,:)       !< neighbour matrix for the shape mask projector
    integer               :: wdim   = 0              !< dimension of K-B window
    integer               :: iwinsz = 0              !< integer half-window size
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
    procedure          :: dfproject_polar
    procedure          :: fdf_project_polar
    procedure, private :: fdf_project_polar_memo
    procedure          :: interp_fcomp
    procedure, private :: interp_fcomp_memo
    procedure, private :: dinterp_fcomp
    procedure, private :: fdf_interp_fcomp
    procedure, private :: fdf_interp_fcomp_memo
    ! DESTRUCTOR
    procedure :: kill_expanded
end type projector

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor of the expanded Fourier matrix
    subroutine expand_cmat( self, alpha )
        class(projector), intent(inout) :: self
        real,             intent(in)    :: alpha !< oversampling factor
        integer, allocatable :: cyck(:), cycm(:), cych(:)
        integer :: h, k, m, phys(3), logi(3), lims(3,2), ldim(3)
        real    :: winsz
        call self%kill_expanded
        ldim = self%get_ldim()
        if( .not.self%is_ft() ) stop 'volume needs to be FTed before call; expand_cmat; simple_projector'
        if( ldim(3) == 1      ) stop 'only for volumes; expand_cmat; simple_projector'
        self%kbwin         = kbinterpol(KBWINSZ, alpha)
        self%iwinsz        = ceiling(self%kbwin%get_winsz() - 0.5)
        if( MEMOIZEKB ) call self%kbwin%memoize(NKBPOINTS)
        self%wdim          = self%kbwin%get_wdim()
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
            &stop 'expanded fourier matrix does not exist; simple_projector::reset_expanded'
        self%cmat_exp = cmplx(0.,0.)
    end subroutine reset_expanded

    ! FOURIER PROJECTORS

    !> \brief  extracts a Fourier plane from the expanded FT matrix of a volume (self)
    subroutine fproject( self, e, fplane )
        class(projector), intent(inout) :: self
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real        :: loc(3)
        integer     :: h, k, sqarg, sqlp, lims(3,2), phys(3), ldim(3)
        complex(sp) :: comp
        lims   = self%loop_lims(2)
        sqlp   = (maxval(lims(:,2)))**2
        ldim   = fplane%get_ldim()
        call fplane%zero_and_flag_ft()
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(h,k,sqarg,loc,phys,comp) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                loc  = matmul(real([h,k,0]), e%get_mat())
                comp = self%interp_fcomp(loc)
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
        real        :: loc(3)
        integer     :: h, k, sqarg, sqlp, lims(3,2), phys(3), ldim(3)
        complex(sp) :: comp
        lims   = self%loop_lims(2)
        sqlp   = (maxval(lims(:,2)))**2
        ldim   = fplane%get_ldim()
        call fplane%zero_and_flag_ft()
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                loc  = matmul(real([h,k,0]), e%get_mat())
                comp = self%interp_fcomp(loc)
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
        real    :: vec(3), loc(3)
        if( MEMOIZEKB )then
            call self%fproject_polar_memo(iref, e, pftcc, iseven)
            return
        endif
        pdim = pftcc%get_pdim()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = pftcc%get_coord(irot,k)
                vec(3)  = 0.
                loc     = matmul(vec,e%get_mat())
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
        real    :: vec(3), loc(3)
        pdim = pftcc%get_pdim()
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = pftcc%get_coord(irot,k)
                vec(3)  = 0.
                loc     = matmul(vec,e%get_mat())
                call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp_memo(loc), iseven)
            end do
        end do
    end subroutine fproject_polar_memo

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
    function interp_fcomp( self, loc )result( comp )
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
            w(i,:,:) = w(i,:,:) * self%kbwin%apod( real(win(1,1)+i-1)-loc(1) )
            w(:,i,:) = w(:,i,:) * self%kbwin%apod( real(win(1,2)+i-1)-loc(2) )
            w(:,:,i) = w(:,:,i) * self%kbwin%apod( real(win(1,3)+i-1)-loc(3) )
        end do
        ! SUM( kernel x components )
        comp = sum( w * self%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2),win(1,3):win(2,3)) ) / sum( w )
    end function interp_fcomp

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
        real(dp)                                                 :: N, D, D2 !numerator, denominator
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
        D   = D**2
        res(1) = cmplx( ( sum( wt  * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wt  ) ) / D2 )
        res(2) = cmplx( ( sum( wp  * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wp  ) ) / D2 )
        res(3) = cmplx( ( sum( wph * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wph ) ) / D2 )
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
        real(dp)                                                 :: N, D, D2 !numerator, denominator
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
        D   = D**2
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
        real(dp)                                                 :: N, D, D2 !numerator, denominator        
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
        D   = D**2
        res(1) = cmplx( ( sum( wt  * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wt  ) ) / D2 )
        res(2) = cmplx( ( sum( wp  * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wp  ) ) / D2 )
        res(3) = cmplx( ( sum( wph * self%cmat_exp(win(1,1):win(2,1),win(1,2):win(2,2),win(1,3):win(2,3)) ) * D - N * sum( wph ) ) / D2 )
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

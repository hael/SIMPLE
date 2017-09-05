! projection of 3D volumes in the Fourier domain by convolution interpolation
! to generate band-pass limited Cartesian and polar 2D Fourier transforms 
module simple_projector
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs       ! use all in there
use simple_kbinterpol, only: kbinterpol
use simple_image,      only: image
use simple_ori,        only: ori
use simple_oris,       only: oris
use simple_params,     only: params
use simple_syslib,     only: alloc_errchk
use simple_math,       only: sqwin_3d
implicit none

public :: projector
private

type, extends(image) :: projector
    private
    type(kbinterpol)      :: kbwin                               !< window function object
    integer               :: ldim_exp(3,2) = 0                   !< expanded FT matrix limits
    complex, allocatable  :: cmat_exp(:,:,:)                     !< expanded FT matrix
    real,    allocatable  :: polweights_mat(:,:,:)               !< polar weights matrix for the image to polar transformer
    integer, allocatable  :: polcyc1_mat(:,:,:)                  !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)                  !< image cyclic adresses for the image to polar transformer
    logical, allocatable  :: is_in_mask(:,:,:)                   !< neighbour matrix for the shape mask projector
    real                  :: winsz      = KBWINSZ                !< window half-width
    real                  :: alpha      = KBALPHA                !< oversampling ratio
    real                  :: harwin     = real(ceiling(KBWINSZ)) !< rounded window half-width
    real                  :: harwin_exp = 1.0                    !< rounded window half-width in expanded routines
    integer               :: wdim       = 2*ceiling(1.0) + 1     !< harwin_exp is argument to ceiling
    logical               :: expanded_exists=.false.             !< indicates FT matrix existence
  contains
    ! CONSTRUCTORS
    procedure          :: expand_cmat
    ! SETTERS
    procedure          :: reset_expanded
    ! INTERPOLATOR
    procedure          :: extr_gridfcomp
    ! FOURIER PROJECTORS
    procedure          :: fproject
    procedure          :: fproject_polar
    procedure          :: interp_fcomp
    ! DESTRUCTOR
    procedure          :: kill_expanded
end type projector

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor of the expanded Fourier matrix
    subroutine expand_cmat( self )
        use simple_math, only: cyci_1d
        class(projector), intent(inout) :: self
        integer, allocatable :: cyck(:), cycm(:), cych(:)
        integer :: h, k, m, phys(3), logi(3)
        integer :: lims(3,2), ldim(3)
        call self%kill_expanded
        ldim = self%get_ldim()
        if( .not.self%is_ft() ) stop 'volume needs to be FTed before call; expand_cmat; simple_projector'
        if( ldim(3) == 1      ) stop 'only for volumes; expand_cmat; simple_projector'
        self%kbwin         = kbinterpol(KBWINSZ, KBALPHA)
        lims               = self%loop_lims(3)
        self%ldim_exp(:,2) = maxval(abs(lims)) + ceiling(self%harwin_exp)
        self%ldim_exp(:,1) = -self%ldim_exp(:,2)
        allocate( self%cmat_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),&
                                &self%ldim_exp(2,1):self%ldim_exp(2,2),&
                                &self%ldim_exp(3,1):self%ldim_exp(3,2)),&
                                &cych(self%ldim_exp(1,1):self%ldim_exp(1,2)),&
                                &cyck(self%ldim_exp(2,1):self%ldim_exp(2,2)),&
                                &cycm(self%ldim_exp(3,1):self%ldim_exp(3,2)), stat=alloc_stat)
        if(alloc_stat/=0)call alloc_errchk("In: expand_cmat; simple_projector", alloc_stat)
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
    subroutine fproject( self, e, fplane, lp )
        class(projector), intent(inout) :: self
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real, optional,   intent(in)    :: lp
        real    :: loc(3)
        integer :: h, k, sqarg, sqlp, lims(3,2), logi(3), phys(3)
        ! init
        if( present(lp) )then
            lims = self%loop_lims(1,lp)
            sqlp = fplane%get_find(lp)**2
        else
            lims = self%loop_lims(2) 
            sqlp = (maxval(lims(:,2)))**2
        endif
        fplane = cmplx(0.,0.)
        !$omp parallel do collapse(2) schedule(static) default(shared)&
        !$omp private(h,k,sqarg,loc,logi,phys) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = dot_product([h,k],[h,k])
                if(sqarg > sqlp)cycle
                ! address
                logi = [h, k, 0]
                loc  = matmul(real(logi), e%get_mat())
                ! set fourier component
                phys = self%comp_addr_phys(logi)
                call fplane%set_fcomp(logi,phys,self%interp_fcomp(loc))
            end do
        end do
        !$omp end parallel do
    end subroutine fproject

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fproject_polar( self, iref, e, pftcc, serial )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_math,             only: deg2rad
        class(projector),        intent(inout) :: self   !< projector object
        integer,                 intent(in)    :: iref   !< which reference
        class(ori),              intent(inout) :: e      !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc  !< object that holds the polar image
        logical, optional,       intent(in)    :: serial !< thread or serial
        integer :: irot, k, ldim(3), pdim(3), ldim_polft(3)
        real    :: vec(3), loc(3)
        logical :: l_serial = .false.
        if( present(serial) )l_serial = serial
        ldim = self%get_ldim()
        if(ldim(3) == 1)stop 'only for interpolation from 3D images; fproject_polar_1; simple_projector'
        ldim_polft(1:2) = ldim(1:2)
        ldim_polft(3)   = 1
        pdim   = pftcc%get_pdim(.false.)
        if( l_serial )then
            ! this is the serial version of the threaded version just below
            do irot=1,pdim(1)
                do k=pdim(2),pdim(3)
                    vec(:2) = pftcc%get_coord(irot,k)
                    vec(3)  = 0.
                    loc     = matmul(vec,e%get_mat())
                    call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp(loc))
                end do
            end do
        else
            ! threaded version
            !$omp parallel do collapse(2) schedule(static) default(shared)&
            !$omp private(irot,k,vec,loc) proc_bind(close)
            do irot=1,pdim(1)
                do k=pdim(2),pdim(3)
                    vec(:2) = pftcc%get_coord(irot,k)
                    vec(3)  = 0.
                    loc     = matmul(vec,e%get_mat())
                    call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp(loc))
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine fproject_polar

    ! INTERPOLATORS
    
    !> \brief  extracts a Fourier component from a transform (self) by gridding
    !!         not expanded
    function extr_gridfcomp( self, loc ) result( comp_sum )
        use simple_math, only: cyci_1d
        class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        complex, allocatable :: comps(:,:,:)
        real,    allocatable :: w(:,:,:)
        integer, allocatable :: cyc1(:), cyc2(:), cyc3(:)
        complex  :: comp_sum
        integer  :: i, j, m, wdim, incr, lims(3,2), win(3,2), logi(3), phys(3)
        real     :: harwin_here
        ! this is neeed in case expand_cmat hasn't been called
        self%kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        ! init
        harwin_here = self%harwin_exp
        win         = sqwin_3d(loc(1),loc(2),loc(3), harwin_here)
        wdim        = win(1,2) - win(1,1) + 1
        lims        = self%loop_lims(3)
        if( self%is_2d() )then
            ! 2D
            ! init
            allocate( cyc1(wdim), cyc2(wdim), w(wdim, wdim, 1), comps(wdim, wdim, 1))
            comps = cmplx(0.,0.)
            w     = 1.
            do i=1,wdim
                incr = i-1
                ! circular addresses
                cyc1(i)  = cyci_1d(lims(1,:), win(1,1)+incr)
                cyc2(i)  = cyci_1d(lims(2,:), win(2,1)+incr)
                ! interpolation kernel matrix
                w(i,:,1) = w(i,:,1) * self%kbwin%apod( real(win(1,1)+incr)-loc(1) )
                w(:,i,1) = w(:,i,1) * self%kbwin%apod( real(win(2,1)+incr)-loc(2) )
            enddo
            ! fetch fourier components
            do i=1, wdim
                do j=1, wdim
                    if(w(i,j,1) == 0.)cycle
                    logi = [cyc1(i), cyc2(j), 0]
                    phys = self%comp_addr_phys(logi)
                    comps(i,j,1) = self%get_fcomp(logi, phys)
                end do
            end do
            ! SUM( kernel x components )
            comp_sum = sum(w * comps)
            deallocate( w,comps,cyc1,cyc2 )
        else
            ! 3D
            ! init
            allocate( cyc1(wdim), cyc2(wdim), cyc3(wdim), w(wdim, wdim, wdim),comps(wdim, wdim, wdim))
            comps = cmplx(0.,0.)
            w     = 1.
            do i=1, wdim
                incr = i-1
                ! circular addresses
                cyc1(i)  = cyci_1d(lims(1,:), win(1,1)+incr)
                cyc2(i)  = cyci_1d(lims(2,:), win(2,1)+incr)
                cyc3(i)  = cyci_1d(lims(3,:), win(3,1)+incr)
                ! interpolation kernel matrix
                w(i,:,:) = w(i,:,:) * self%kbwin%apod( real(win(1,1)+incr)-loc(1) )
                w(:,i,:) = w(:,i,:) * self%kbwin%apod( real(win(2,1)+incr)-loc(2) )
                w(:,:,i) = w(:,:,i) * self%kbwin%apod( real(win(3,1)+incr)-loc(3) )
            enddo
            ! fetch fourier components
            do i=1, wdim
                do j=1, wdim
                    if(all(w(i,j,:) == 0.))cycle
                    do m=1, wdim
                        if(w(i,j,m) == 0.)cycle
                        logi = [cyc1(i), cyc2(j), cyc3(m)]
                        phys = self%comp_addr_phys(logi)
                        comps(i,j,m) = self%get_fcomp(logi, phys)
                    enddo
                end do
            end do
            ! SUM( kernel x components )
            comp_sum = sum(w * comps)
            deallocate( w,comps,cyc1,cyc2,cyc3 )
        endif
    end function extr_gridfcomp

    !>  \brief is to interpolate from the expanded complex matrix 
    function interp_fcomp( self, loc )result( comp )
        class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        complex :: comp
        real    :: w(1:self%wdim,1:self%wdim,1:self%wdim)
        integer :: i, win(3,2)
        ! interpolation kernel window
        win  = sqwin_3d(loc(1), loc(2), loc(3), self%harwin_exp)
        ! interpolation kernel matrix
        w = 1.
        do i=1,self%wdim
            w(i,:,:) = w(i,:,:) * self%kbwin%apod( real(win(1,1)+i-1)-loc(1) )
            w(:,i,:) = w(:,i,:) * self%kbwin%apod( real(win(2,1)+i-1)-loc(2) )
            w(:,:,i) = w(:,:,i) * self%kbwin%apod( real(win(3,1)+i-1)-loc(3) )
        end do
        ! SUM( kernel x components )
        comp = sum( w * self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2),win(3,1):win(3,2)) )
    end function interp_fcomp
    

    !>  \brief  is a destructor of expanded matrices (imgpolarizer AND expanded projection of)
    subroutine kill_expanded( self )
        class(projector), intent(inout) :: self !< projector instance
        if( allocated(self%cmat_exp) )deallocate(self%cmat_exp)    
        self%ldim_exp        = 0
        self%expanded_exists = .false.
    end subroutine kill_expanded

end module simple_projector

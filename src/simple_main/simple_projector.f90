!>  \brief  SIMPLE projector class
module simple_projector
use simple_defs       ! use all in there
use simple_kbinterpol ! use all in there
use simple_image,     only: image
use simple_ori,       only: ori
use simple_oris,      only: oris
use simple_params,    only: params
use simple_gridding,  only: prep4cgrid
use simple_jiffys,    only: alloc_err, progress
use simple_math,      only: recwin_3d, euclid, sqwin_3d
implicit none

public :: projector
private

logical, parameter :: DEBUG = .true.

type, extends(image) :: projector
    private
    integer               :: ldim_exp(3,2)=0         !< expanded FT matrix limits
    complex, allocatable  :: cmat_exp(:,:,:)         !< expanded FT matrix
    real,    allocatable  :: polweights_mat(:,:,:)   !< polar weights matrix for the image to polar transformer
    integer, allocatable  :: polcyc1_mat(:,:,:)      !< image cyclic adresses for the image to polar transformer
    integer, allocatable  :: polcyc2_mat(:,:,:)      !< image cyclic adresses for the image to polar transformer
    logical, allocatable  :: is_in_mask(:,:,:)       !< neighbour matrix for the shape mask projector
    real                  :: winsz     = 1.5         !< window half-width
    real                  :: alpha     = 2.          !< oversampling ratio
    real                  :: harwin    = 1.          !< rounded window half-width
    real                  :: harwin_exp= 1.          !< rounded window half-width in expanded routines
    logical               :: expanded_exists=.false. !< indicates FT matrix existence
  contains
    ! CONSTRUCTORS
    procedure          :: expand_cmat
    procedure          :: compress_cmat
    ! GETTERS
    procedure          :: get_harwin
    procedure          :: get_harwin_exp
    procedure          :: exp_exists
    ! SETTERS
    procedure          :: reset_expanded
    ! FOURIER PROJECTORS
    procedure          :: fproject
    procedure          :: fproject_expanded
    procedure          :: fproject_polar
    procedure, private :: fproject_polar_expanded
    procedure          :: img2polarft
    procedure          :: extr_gridfcomp
    procedure, private :: interp_fcomp_expanded
    ! IMAGE TO POLAR TRANSFORMER
    procedure          :: init_imgpolarizer
    procedure          :: imgpolarizer
    procedure, private :: kill_imgpolarizer
    procedure          :: kill_expanded
    ! MASK REAL-SPACE PROJECTION
    procedure          :: init_env_rproject
    procedure          :: env_rproject
    procedure          :: kill_env_rproject

end type projector

contains

    ! CONSTRUCTOR

     !>  \brief  is a constructor of the expanded Fourier matrix
    subroutine expand_cmat( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: cyci_1d
        class(projector), intent(inout) :: self
        integer, allocatable :: cyck(:), cycm(:), cych(:)
        integer :: h, k, m, alloc_stat, phys(3), logi(3)
        integer :: lims(3,2), ldim(3)
        ldim = self%get_ldim()
        if( .not.self%is_ft() ) stop 'volume needs to be FTed before call; expand_cmat; simple_image'
        if( ldim(3) == 1      ) stop 'only for volumes; expand_cmat; simple_image'
        self%winsz  = get_kb_winsz()
        self%harwin = real(ceiling(self%winsz))
        self%alpha  = get_kb_alpha()
        lims               = self%loop_lims(3)
        self%ldim_exp(:,2) = maxval(abs(lims)) + ceiling(self%harwin_exp)
        self%ldim_exp(:,1) = -self%ldim_exp(:,2)
        if( allocated(self%cmat_exp) ) deallocate(self%cmat_exp)
        allocate( self%cmat_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),&
                                &self%ldim_exp(2,1):self%ldim_exp(2,2),&
                                &self%ldim_exp(3,1):self%ldim_exp(3,2)),&
                                &cych(self%ldim_exp(1,1):self%ldim_exp(1,2)),&
                                &cyck(self%ldim_exp(2,1):self%ldim_exp(2,2)),&
                                &cycm(self%ldim_exp(3,1):self%ldim_exp(3,2)), stat=alloc_stat)
        call alloc_err("In: expand_cmat; simple_projector", alloc_stat)
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
        !$omp parallel do collapse(3) schedule(auto) default(shared) private(h,k,m,logi,phys)
        do h = self%ldim_exp(1,1),self%ldim_exp(1,2)
            do k = self%ldim_exp(2,1),self%ldim_exp(2,2)
                do m = self%ldim_exp(3,1),self%ldim_exp(3,2)
                    logi = [cych(h),cyck(k),cycm(m)]
                    phys = self%comp_addr_phys(logi)
                    self%cmat_exp(h,k,m) = self%get_fcomp(logi,phys)
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate( cych,cyck,cycm )
        self%expanded_exists = .true.
    end subroutine expand_cmat

    !>  \brief converts the expanded matrix to standard imaginary representation
    subroutine compress_cmat( self )
        class(projector), intent(inout) :: self
        integer :: h, k, m, logi(3), phys(3)
        integer :: lims(3,2)
        if( .not. self%expanded_exists )then
            stop 'expanded complex matrix does not exist; simple_projector :: compress_cmat'
        endif
        lims = self%loop_lims(2) ! excluding redundant Friedel mates
        self = cmplx(0.,0.)
        !$omp parallel do collapse(3) schedule(auto) default(shared) private(h,k,m,logi,phys)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                do m = lims(3,1),lims(3,2)
                    logi = [h,k,m]
                    phys = self%comp_addr_phys(logi)
                    call self%set_fcomp(logi,phys,self%cmat_exp(h,k,m))
                end do
            end do 
        end do
        !$omp end parallel do     
    end subroutine compress_cmat

    ! GETTERS

    !>  \brief  returns the hard window size for projection
    function get_harwin( self ) result( winsz )
        class(projector), intent(inout) :: self
        real :: winsz
        winsz = self%harwin
    end function get_harwin

    !>  \brief  returns the hard window size for projection from expanded matrix
    function get_harwin_exp( self ) result( winsz )
        class(projector), intent(inout) :: self
        real :: winsz
        winsz = self%harwin_exp
    end function get_harwin_exp

    !>  \brief  indicates expanded fourier matix existance
    function exp_exists( self )result( answer )
        class(projector), intent(inout) :: self
        logical :: answer
        answer = self%expanded_exists
    end function exp_exists

    ! SETTERS

    !>  \brief  sets the expanded Fourier matrix to zero
    subroutine reset_expanded( self )
        class(projector), intent(inout) :: self
        if( .not.self%expanded_exists )&
            &stop 'expanded fourier matrix does not exist; simple_image::reset_expanded'
        self%cmat_exp = cmplx(0.,0.)
    end subroutine reset_expanded

    ! FOURIER PROJECTORS

    !> \brief  extracts a Fourier plane from a volume (self)
    subroutine fproject( self, e, fplane, lp )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real, optional,   intent(in)    :: lp
        real    :: vec(3), loc(3)
        integer :: h, k, sqarg, sqlp, lims(3,2), logi(3), phys(3)
        ! init
        fplane = cmplx(0.,0.)
        if( present(lp) )then
            lims = self%loop_lims(1,lp)
        else
            lims = self%loop_lims(2) ! Nyqvist default low-pass limit
        endif
        sqlp = (maxval(lims(:,2)))**2
        !$omp parallel do collapse(2) schedule(auto) default(shared) private(h,k,sqarg,vec,loc,logi,phys)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = h*h+k*k
                if( sqarg <= sqlp )then
                    ! address
                    vec(1) = real(h)
                    vec(2) = real(k)
                    vec(3) = 0.
                    loc = matmul(vec,e%get_mat())
                    ! set fourier component
                    logi = [h,k,0]
                    phys = self%comp_addr_phys(logi)
                    call fplane%set_fcomp(logi,phys,self%extr_gridfcomp(loc))
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine fproject

    !> \brief  extracts a Fourier plane from the expanded FT matrix of a volume (self)
    subroutine fproject_expanded( self, e, fplane, lp )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real, optional,   intent(in)    :: lp
        real    :: vec(3), loc(3)
        integer :: h, k, sqarg, sqlp, wdim, lims(3,2), logi(3), phys(3)
        ! init
        lims = self%loop_lims(2) 
        if( present(lp) )then
            lims = self%loop_lims(1,lp)
            sqlp = fplane%get_find(lp)**2
        else
            sqlp = (maxval(lims(:,2)))**2
        endif
        fplane = cmplx(0.,0.)
        wdim = 2*ceiling(self%harwin_exp) + 1 ! interpolation kernel window size
        !$omp parallel do collapse(2) schedule(auto) default(shared) private(h,k,sqarg,vec,loc,logi,phys)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                sqarg = h*h+k*k
                if( sqarg <= sqlp )then
                    ! address
                    vec(1) = real(h)
                    vec(2) = real(k)
                    vec(3) = 0.
                    loc    = matmul(vec,e%get_mat())
                    ! set fourier component
                    logi = [h,k,0]
                    phys = self%comp_addr_phys(logi)
                    call fplane%set_fcomp(logi,phys,self%interp_fcomp_expanded(wdim, loc))
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine fproject_expanded

    !> \brief  extracts a polar FT from a volume (self)
    subroutine fproject_polar( self, iref, e, pftcc, expanded )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_params, only: params
        class(projector),        intent(inout) :: self  !< projector instance
        integer,                 intent(in)    :: iref  !< logical reference index [1,nrefs]
        class(ori),              intent(inout) :: e     !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc !< polarft_corrcalc object to be filled
        logical,optional,        intent(in)    :: expanded
        integer :: irot, k, ldim_fvol(3), ldim_pft(3), pdim(3)
        real    :: vec(3), loc(3)
        logical :: l_exp
        l_exp = .false.
        if( present(expanded) )then
            l_exp = expanded
            if( l_exp .and. .not. self%exp_exists() )&
                &stop 'expanded fourier components matrix does not exists; fproject_polar_2; simple_projector'
        endif
        ldim_fvol = self%get_ldim()
        ldim_pft  = pftcc%get_ldim()
        if( .not. all(ldim_fvol(1:2) == ldim_pft(1:2)) ) stop 'nonconforming logical dimensions; fproject_polar_2; simple_projector'
        pdim = pftcc%get_pdim(.false.)
        if( l_exp )then
            call self%fproject_polar_expanded(iref, e, pftcc)
        else
            !$omp parallel do collapse(2) schedule(auto) default(shared) private(irot,k)
            do irot=1,pdim(1)
                do k=pdim(2),pdim(3)
                    vec(:2) = pftcc%get_coord(irot,k) 
                    vec(3)  = 0.
                    loc     = matmul(vec,e%get_mat())
                    call pftcc%set_ref_fcomp(iref, irot, k, self%extr_gridfcomp(loc))
                end do
            end do
            !$omp end parallel do
        endif
        call pftcc%memoize_sqsum_ref(iref)
    end subroutine fproject_polar

    !> \brief  extracts a polar FT from a volume's expanded FT (self)
    subroutine fproject_polar_expanded( self, iref, e, pftcc )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_math,             only: deg2rad
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector),        intent(inout) :: self  !< projector object
        integer,                 intent(in)    :: iref  !< which reference
        class(ori),              intent(inout) :: e     !< orientation
        class(polarft_corrcalc), intent(inout) :: pftcc !< object that holds the polar image
        integer :: irot, k, wdim, ldim(3), pdim(3), ldim_polft(3)
        real    :: vec(3), loc(3)
        ldim = self%get_ldim()
        if( ldim(3) == 1 )        stop 'only for interpolation from 3D images; fproject_polar_1; simple_projector'
        if( .not. self%is_ft() )  stop 'volume needs to be FTed before call; fproject_polar_1; simple_projector'
        ldim_polft(1:2) = ldim(1:2)
        ldim_polft(3)   = 1
        pdim = pftcc%get_pdim(.false.)
        wdim = 2*ceiling(self%harwin_exp) + 1 ! interpolation kernel window size
        !$omp parallel do collapse(2) schedule(auto) default(shared) private(irot,k,vec,loc)
        do irot=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2) = pftcc%get_coord(irot,k)
                vec(3)  = 0.
                loc     = matmul(vec,e%get_mat())
                call pftcc%set_ref_fcomp(iref, irot, k, self%interp_fcomp_expanded( wdim, loc ))
            end do
        end do
        !$omp end parallel do
    end subroutine fproject_polar_expanded
    
    !> \brief  transfers a 2D Cartesian image (self) to polar Fourier
    subroutine img2polarft( self, ind, pftcc, isptcl )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use gnufor2
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector),        intent(inout) :: self   !< instance
        integer,                 intent(in)    :: ind    !< logical particle index [fromp,top]
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        logical, optional,       intent(in)    :: isptcl !< to indicate whether particle or reference
        complex, allocatable :: pft(:,:)
        integer :: i, k, ldim_img(3), ldim_pft(3), pdim(3), alloc_stat
        real    :: vec(3)
        logical :: iisptcl
        iisptcl = .true.
        if( present(isptcl) ) iisptcl = isptcl
        ldim_img = self%get_ldim()
        if( ldim_img(3) > 1 )      stop 'only for interpolation from 2D images; img2polarft_2; simple_projector'
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; img2polarft_2; simple_projector'
        ldim_pft = pftcc%get_ldim()
        if( .not. all(ldim_img == ldim_pft) )then
            print *, 'ldim_img: ', ldim_img
            print *, 'ldim_pft: ', ldim_pft
            stop 'logical dimensions do not match; img2polarft_2; simple_projector'
        endif
        if( .not. self%is_ft() ) stop 'image needs to FTed before this operation; simple_projector::img2polarft_2'
        pdim = pftcc%get_pdim(iisptcl)
        allocate( pft(pdim(1),pdim(2):pdim(3)), stat=alloc_stat )
        call alloc_err("In: img2polarft_2; simple_projector", alloc_stat)
        !$omp parallel do collapse(2) schedule(auto) default(shared) private(i,k,vec)
        do i=1,pdim(1)
            do k=pdim(2),pdim(3)
                vec(:2)  = pftcc%get_coord(i,k)
                vec(3)   = 0.
                pft(i,k) = self%extr_gridfcomp(vec)
            end do
        end do
        !$omp end parallel do
        if( iisptcl )then
            call pftcc%set_ptcl_pft(ind, pft)
        else
            call pftcc%set_ref_pft(ind, pft)
        endif
        ! kill the remains
        deallocate(pft)
    end subroutine img2polarft

    ! INTERPOLATORS
    
    !> \brief  extracts a Fourier component from a transform (self) by gridding
    function extr_gridfcomp( self, loc ) result( comp_sum )
        use simple_math, only: cyci_1d
        class(projector), intent(inout) :: self
        real,             intent(in)    :: loc(3)
        real,    allocatable :: w1(:), w2(:), w3(:)
        integer, allocatable :: cyc1(:), cyc2(:), cyc3(:)
        integer :: alloc_stat, i, j, m
        integer :: lims(3,2), win(3,2), logi(3), phys(3)
        complex :: comp, comp_sum, zero
        real    :: harwin_here
        harwin_here = 2.
        zero        = cmplx(0.,0.)
        comp_sum    = zero
        win         = recwin_3d(loc(1),loc(2),loc(3),harwin_here)
        lims        = self%loop_lims(3)
        allocate( w1(win(1,1):win(1,2)), w2(win(2,1):win(2,2)), &
            & cyc1(win(1,1):win(1,2)), cyc2(win(2,1):win(2,2)), stat=alloc_stat)
        call alloc_err("In: extr_gridfcomp; simple_projector, 1", alloc_stat)
        do i=win(1,1),win(1,2)
            w1(i)   = kb_apod(real(i)-loc(1))
            cyc1(i) = cyci_1d( lims(1,:),i )
        end do
        do j=win(2,1),win(2,2)
            w2(j)   = kb_apod(real(j)-loc(2))
            cyc2(j) = cyci_1d( lims(2,:),j )
        end do
        if( self%is_2d() )then
            ! 2D
            do i=win(1,1),win(1,2)
                if( w1(i) == 0. ) cycle
                do j=win(2,1),win(2,2)
                    if( w2(j) == 0. ) cycle
                    logi = [cyc1(i),cyc2(j),0]
                    phys = self%comp_addr_phys(logi)
                    comp = self%get_fcomp(logi,phys)
                    if( comp .eq. zero ) cycle
                    comp_sum = comp_sum+comp*w1(i)*w2(j)
                end do
            end do
            deallocate( w1,w2,cyc1,cyc2 )
        else
            ! 3D
            allocate( w3(win(3,1):win(3,2)), cyc3(win(3,1):win(3,2)), stat=alloc_stat)
            call alloc_err("In: extr_gridfcomp; simple_projector, 2", alloc_stat)
            do m=win(3,1),win(3,2)
                w3(m)   = kb_apod(real(m)-loc(3))
                cyc3(m) = cyci_1d( lims(3,:),m )
            end do
            do i=win(1,1),win(1,2)
                if( w1(i) == 0. ) cycle
                do j=win(2,1),win(2,2)
                    if( w2(j) == 0. ) cycle
                    do m=win(3,1),win(3,2)
                        if( w3(m) == 0. ) cycle
                        logi = [cyc1(i),cyc2(j),cyc3(m)]
                        phys = self%comp_addr_phys(logi)
                        comp = self%get_fcomp(logi,phys)
                        if( comp .eq. zero ) cycle
                        comp_sum = comp_sum+comp*w1(i)*w2(j)*w3(m)
                    end do
                end do
            end do
            deallocate(w1,w2,w3,cyc1,cyc2,cyc3)
        endif
    end function extr_gridfcomp

    !>  \brief is to interpolate from the expanded complex matrix 
    function interp_fcomp_expanded( self, wdim, loc )result( comp )
        class(projector), intent(inout) :: self
        integer,          intent(in)    :: wdim
        real,             intent(in)    :: loc(3)
        complex :: comp
        real    :: w(1:wdim,1:wdim,1:wdim)
        integer :: i, wlen, win(3,2)
        ! interpolation kernel window
        win  = sqwin_3d(loc(1), loc(2), loc(3), self%harwin_exp)
        wlen = wdim**3
        ! interpolation kernel matrix
        w = 1.
        do i=1,wdim
            w(i,:,:) = w(i,:,:) * kb_apod( real(win(1,1)+i-1)-loc(1) )
            w(:,i,:) = w(:,i,:) * kb_apod( real(win(2,1)+i-1)-loc(2) )
            w(:,:,i) = w(:,:,i) * kb_apod( real(win(3,1)+i-1)-loc(3) )
        end do
        ! SUM( kernel x components )
        comp = dot_product( reshape(w,(/wlen/)), reshape(self%cmat_exp( win(1,1):win(1,2),&
            &win(2,1):win(2,2),win(3,1):win(3,2)), (/wlen/)) )
    end function interp_fcomp_expanded
    
    ! IMAGE TO POLAR FT TRANSFORMER

    !> \brief  initialises the image polarizer
    subroutine init_imgpolarizer( self, pftcc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math,             only: sqwin_2d, cyci_1d
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        real, allocatable :: w(:,:)
        real              :: loc(2)
        integer           :: pdim(3), win(2,2), lims(3,2)
        integer           :: i, k, l, wdim, wlen, alloc_stat, cnt
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; init_imgpolarizer; simple_projector'
        call self%kill_imgpolarizer
        wdim = 2*ceiling(self%harwin_exp) + 1
        wlen = wdim**2
        pdim = pftcc%get_pdim(.true.)
        lims = self%loop_lims(3)
        allocate( self%polcyc1_mat(1:pdim(1), pdim(2):pdim(3), 1:wdim),&
                  &self%polcyc2_mat(1:pdim(1), pdim(2):pdim(3), 1:wdim),&
                  &self%polweights_mat(1:pdim(1), pdim(2):pdim(3), 1:wlen),&
                  &w(1:wdim,1:wdim), stat=alloc_stat)
        call alloc_err('in simple_projector :: init_imgpolarizer', alloc_stat)
        !$omp parallel do schedule(auto) default(shared) private(i,k,l,w,loc,cnt,win)
        do i=1,pdim(1)
            do k=pdim(2),pdim(3)
                ! polar coordinates
                loc = pftcc%get_coord(i,k)
                win = sqwin_2d(loc(1), loc(2), self%harwin_exp)
                w   = 1.
                cnt = 0
                do l=1,wdim
                    cnt = cnt + 1
                    ! interpolation weights
                    w(l,:) = w(l,:) * kb_apod( real(win(1,1)+l-1)-loc(1) )
                    w(:,l) = w(:,l) * kb_apod( real(win(2,1)+l-1)-loc(2) )
                    ! cyclic addresses
                    self%polcyc1_mat(i, k, cnt) = cyci_1d(lims(1,:), win(1,1)+l-1)
                    self%polcyc2_mat(i, k, cnt) = cyci_1d(lims(2,:), win(2,1)+l-1)
                end do
                self%polweights_mat(i,k,:) = reshape(w,(/wlen/))
            enddo
        enddo
        !$omp end parallel do
        deallocate(w)
    end subroutine init_imgpolarizer

    !> \brief  creates the polar Fourier transform
    subroutine imgpolarizer( self, pftcc, img_ind, isptcl )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: sqwin_2d, cyci_1d
        use simple_polarft_corrcalc, only: polarft_corrcalc
        class(projector),        intent(inout) :: self   !< projector instance
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        integer,                 intent(in)    :: img_ind
        logical, optional,       intent(in)    :: isptcl
        complex, allocatable :: pft(:,:), comps(:,:)
        integer :: i, k, l, m, alloc_stat, windim, vecdim, addr_l
        integer :: lims(3,2), ldim_img(3), ldim_pft(3), pdim(3), logi(3), phys(3)
        logical :: iisptcl
        if( .not. allocated(self%polweights_mat) )&
        &stop 'the imgpolarizer has not been initialized!; simple_projector :: imgpolarizer'
        iisptcl = .true.
        if( present(isptcl) ) iisptcl = isptcl
        ldim_img = self%get_ldim()
        if( ldim_img(3) > 1 )      stop 'only for interpolation from 2D images; imgpolarizer; simple_projector'
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; imgpolarizer; simple_projector'
        ldim_pft = pftcc%get_ldim()
        if( .not. all(ldim_img == ldim_pft) )then
            print *, 'ldim_img: ', ldim_img
            print *, 'ldim_pft: ', ldim_pft
            stop 'logical dimensions do not match; imgpolarizer; simple_projector'
        endif
        if( .not.self%is_ft() ) stop 'image needs to FTed before this operation; simple_projector :: imgpolarizer'
        pdim   = pftcc%get_pdim(iisptcl)
        windim = 2*ceiling(self%harwin_exp) + 1
        vecdim = windim**2
        allocate( pft(pdim(1),pdim(2):pdim(3)), comps(1:windim,1:windim), stat=alloc_stat )
        call alloc_err("In: imgpolarizer; simple_projector", alloc_stat)
        lims = self%loop_lims(3)
        !$omp parallel do collapse(2) schedule(auto) default(shared) private(i,k,l,m,logi,phys,comps,addr_l)
        do i=1,pdim(1)
            do k=pdim(2),pdim(3)
                do l=1,windim
                    addr_l = self%polcyc1_mat(i,k,l)
                    do m=1,windim
                        logi = [addr_l,self%polcyc2_mat(i,k,m),0]
                        phys = self%comp_addr_phys(logi)
                        comps(l,m) = self%get_fcomp(logi,phys)
                    enddo
                enddo
                pft(i,k) = dot_product(self%polweights_mat(i,k,:), reshape(comps,(/vecdim/)))
            end do
        end do
        !$omp end parallel do
        if( iisptcl )then
            call pftcc%set_ptcl_pft(img_ind, pft)
        else
            call pftcc%set_ref_pft(img_ind, pft)
        endif
        ! kill the remains
        deallocate(pft, comps)
    end subroutine imgpolarizer


    ! REAL-SPACE PROJECTOR

    !>  \brief  
    subroutine init_env_rproject( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self   !< projector instance
        real, allocatable :: rmat(:,:,:)
        real              :: thresh
        integer           :: ldim(3),i, ii, jj, j, k, orig(3)
        call self%kill_env_rproject
        ldim = self%get_ldim()
        if( ldim(3) == 1 )          stop 'only for Volumes; env_rproject; simple_projector'
        if( .not. self%even_dims() )stop 'even dimensions assumed; env_rproject; simple_projector'
        if( self%is_ft() )          stop 'real space only; env_rproject; simple_projector'
        ! init
        call self%norm_bin          ! ensures [0;1] range
        rmat   = self%get_rmat()
        thresh = 0.9999             ! prior soft masking is discarded
        orig = ldim/2+1
        allocate( self%is_in_mask(1-orig(1):ldim(1)-orig(1),&
                                 &1-orig(2):ldim(2)-orig(2),&
                                 &1-orig(3):ldim(3)-orig(3)) )
        self%is_in_mask = .false.
        !$omp parallel do default(shared) private(i,ii,jj,j,k)
        do i=1,ldim(1)-1
            ii = i-orig(1)
            do j=1,ldim(2)-1
                jj = j-orig(2)
                do k=1,ldim(3)-1
                    ! if any of the 8 neighbors is in hard mask, value is set to true
                    self%is_in_mask(ii, jj, k-orig(3)) = any(rmat(i:i+1, j:j+1, k:k+1) >= thresh)
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate(rmat)            
    end subroutine init_env_rproject

    !>  \brief  Envelope projection effector
    subroutine env_rproject(self, e, img, maxrad)
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self   !< projector instance
        class(ori),       intent(inout)    :: e      !< Euler angle
        type(image),      intent(inout) :: img    !< resulting projection image
        real,             intent(in)    :: maxrad !< project inside this radius
        real              :: incr_i(3), incr_j(3), incr_k(3), ray_k(3), corner(3), mmaxrad
        integer           :: orig(3), ldim(3),  lims(2), inds(3), i, j, k, sqmaxrad
        if( .not.allocated(self%is_in_mask) )stop 'the envelope projector has not been initialized'
        ldim = self%get_ldim()
        if( ldim(3) == 1 )          stop 'only for Volumes; env_rproject; simple_projector'
        if( .not. self%even_dims() )stop 'even dimensions assumed; env_rproject; simple_projector'
        if( self%is_ft() )          stop 'real space only; env_rproject; simple_projector'
        ! init
        img      = 0.
        orig     = ldim/2+1
        mmaxrad  = min(maxrad,real(ldim(1))/2.-1.)
        sqmaxrad = nint(mmaxrad**2)
        lims(1)  = orig(1) - ceiling(mmaxrad)
        lims(2)  = orig(2) + ceiling(mmaxrad)
        incr_i   = matmul([1., 0., 0.], e%get_mat())
        incr_j   = matmul([0., 1., 0.], e%get_mat())
        incr_k   = matmul([0., 0., 1.], e%get_mat())
        corner   = matmul(-real([mmaxrad+1,mmaxrad+1,mmaxrad+1]), e%get_mat())
        !$omp parallel do collapse(2) default(shared) private(j,i,k,ray_k,inds)
        do i=lims(1),lims(2) 
            do j=lims(1),lims(2)
                if( (i-orig(1))**2+(j-orig(2))**2 > sqmaxrad )cycle
                ray_k = corner + real(i-lims(1)+1)*incr_i + real(j-lims(1)+1)*incr_j
                do k = lims(1),lims(2)
                    ray_k = ray_k + incr_k
                    inds  = floor(ray_k)
                    if(dot_product(inds,inds) > sqmaxrad) cycle
                    if( self%is_in_mask(inds(1), inds(2), inds(3)) )then
                        call img%set([i,j,1], 1.)
                        exit
                    endif                    
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine env_rproject

    ! DESTRUCTORS

    !>  \brief  is a detructor 
    subroutine kill_imgpolarizer( self )
        class(projector), intent(inout) :: self !< projector instance
        if( allocated(self%polweights_mat) ) deallocate(self%polweights_mat)
        if( allocated(self%polcyc1_mat)    ) deallocate(self%polcyc1_mat)
        if( allocated(self%polcyc2_mat)    ) deallocate(self%polcyc2_mat)
    end subroutine kill_imgpolarizer

    !>  \brief  is a detructor 
    subroutine kill_expanded( self )
        class(projector), intent(inout) :: self !< projector instance
        call self%kill_imgpolarizer
        if( allocated(self%cmat_exp) )deallocate(self%cmat_exp)    
        self%ldim_exp        = 0
        self%expanded_exists = .false.
    end subroutine kill_expanded

    !>  \brief  
    subroutine kill_env_rproject( self )
        class(projector), intent(inout) :: self !< projector instance
        if( allocated(self%is_in_mask) )deallocate(self%is_in_mask)          
    end subroutine kill_env_rproject

end module simple_projector

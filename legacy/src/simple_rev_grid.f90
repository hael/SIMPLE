!>  \brief  SIMPLE reverse gridding class
module simple_rev_grid
use simple_winfuns, only: winfuns
use simple_defs     ! singleton
implicit none

public :: rev_grid, test_rev_grid
private

type :: rev_grid
    private
    type(winfuns)         :: wfuns         !< window functions object
    character(len=STDLEN) :: wfun_str='kb' !< kernel, string descriptor
    real                  :: harwin=2.     !< rounded window half-width
    real                  :: winsz=1.5     !< window half-width
    real                  :: alpha=2.      !< oversampling ratio
    integer               :: xdimpd        !< padded Fourier dim
    integer               :: boxpd         !< padded box size
    real                  :: halfbox       !< half box size (assuming square imgs)
    integer               :: box           !< box size
    integer               :: xdim          !< Fourier size
  contains
    procedure :: prep4rev_grid_real
    procedure :: prep4rev_grid_imag
    procedure :: interpolate_real
    procedure :: interpolate_imag
    procedure :: pre_rev_grid_corr
    procedure :: rt
end type

interface rev_grid
    module procedure constructor
end interface

contains

    !> \brief  is a constructor
    function constructor( box, wfun, winsz, alpha ) result( self )
        integer, intent(in)                         :: box
        character(len=STDLEN), intent(in), optional :: wfun   !< window function, string descriptor
        real, intent(in), optional                  :: winsz  !< window half-width
        real, intent(in), optional                  :: alpha  !< oversampling ratio
        type(rev_grid)                              :: self !< instance
        if( present(wfun) )then
            self%wfun_str = wfun
        endif
        if( present(winsz) )then
            self%winsz = winsz
            self%harwin = real(ceiling(self%winsz))
        endif
        if( present(alpha) )then
            self%alpha = alpha
        endif
        self%wfuns   = winfuns(self%wfun_str,self%winsz,self%alpha)
        self%xdimpd  = round2even(self%alpha*real(box/2))
        self%xdim    = box/2
        self%boxpd   = 2*self%xdimpd
        self%halfbox = real(self%boxpd)/2.
        self%box     = box
    end function
    
    !> \brief  prepares an image for reverse gridding interpolation
    function prep4rev_grid_real( self, img ) result( img4grid )
        use enfors_image, only: image
        class(rev_grid), intent(inout) :: self
        class(image), intent(inout)    :: img
        integer                        :: lims(3,2), h, k, l, alloc_stat
        real, allocatable              :: w(:)
        real                           :: const
        type(image)                    :: img4grid
        ! FT & pad
        call img%fwd_ft
        if( img%is_2d() )then
            call img4grid%new([self%boxpd,self%boxpd,1], img%get_smpd())
        else
            call img4grid%new([self%boxpd,self%boxpd,self%boxpd], img%get_smpd())
        endif
        call img%pad(img4grid)
        ! calculate the instrument function
        lims = img4grid%loop_lims(2)
        allocate( w(-self%xdimpd:self%xdimpd), stat=alloc_stat)
        call alloc_err("In: prep4rev_grid; simple_rev_grid", alloc_stat)
         w = 0.
        const = 1./(self%winsz*real(self%xdimpd))
        do k=-self%xdimpd,self%xdimpd
            w(k) = self%wfuns%eval_instr(real(k)*const)
        end do
        ! divide by the instrument function
        do h=lims(1,1),lims(1,2)
            if( w(h) == 0. ) cycle
            do k=lims(2,1),lims(2,2)
                if( w(k) == 0. ) cycle
                do l=lims(3,1),lims(3,2)
                    if( w(l) == 0. ) cycle
                    call img4grid%div([h,k,l], w(h)*w(k)*w(l))
                end do
            end do
        end do
        ! reverse FT
        call img4grid%bwd_ft
        deallocate(w)
    end function
    
    !> \brief  for pre-reverse gridding correction
!    subroutine pre_rev_grid_corr( self, img )
!        use enfors_image, only: image
!        class(rev_grid), intent(inout) :: self !< instance
!        class(image), intent(inout)    :: img
!        integer :: ldim(3), i, j, k, alloc_stat
!        real    :: ci, cj, ck, w, arg
!        real, allocatable :: w1(:), w2(:), w3(:)
!        if( img%is_ft() )then
!            stop 'need real image for pre-gridding correction; pre_rev_grid_corr; simple_projector'
!        endif
!        ldim = img%get_ldim()
!        allocate( w1(ldim(1)), w2(ldim(2)), w3(ldim(3)), stat=alloc_stat )
!        call alloc_err("In: pre_rev_grid_corr; simple_projector", alloc_stat)
!        ci = -real(ldim(1))/2.
!        do i=1,ldim(1)
!            w1(i) = self%wfuns%eval_instr(ci/real(ldim(1)))
!            ci = ci+1.
!        end do
!        cj = -real(ldim(2))/2.
!        do j=1,ldim(2)
!            w2(j) = self%wfuns%eval_instr(cj/real(ldim(2)))
!            cj = cj+1.
!        end do
!        ck = -real(ldim(3))/2.
!        do k=1,ldim(3)
!            w3(k) = self%wfuns%eval_instr(ck/real(ldim(3)))
!            ck = ck+1.
!        end do
!        do i=1,ldim(1)
!            if( w1(i) == 0. ) cycle
!            do j=1,ldim(2)
!                if( w2(j) == 0. ) cycle
!                do k=1,ldim(3)
!                    if( w3(k) == 0. ) cycle
!                    w = w1(i)*w2(j)*w3(k)
!                    call img%div([i,j,k],w )
!                end do
!            end do
!        end do
!        deallocate(w1,w2,w3)
!    end subroutine
    
    !> \brief  prepares an image for reverse gridding interpolation
    function prep4rev_grid_imag( self, img ) result( img4grid )
        use enfors_image, only: image
        class(rev_grid), intent(inout) :: self
        class(image), intent(inout)    :: img
        integer                        :: lims(3,2), h, k, l, alloc_stat, ldim(3)
        real, allocatable              :: w(:)
        real                           :: const
        type(image)                    :: tmp, img4grid
        tmp = img
        call self%pre_rev_grid_corr(tmp)
        if( img%is_2d() )then
            call img4grid%new([self%boxpd,self%boxpd,1], img%get_smpd())
        else
            call img4grid%new([self%boxpd,self%boxpd,self%boxpd], img%get_smpd())
        endif
        call tmp%pad(img4grid)
        call img4grid%fwd_ft
        call tmp%kill
    end function
    
    !> \brief  interpolates an arbitrarily located pixel value
    function interpolate_real( self, img, x, y ) result( pix )
        use enfors_image, only: image
        class(rev_grid), intent(inout) :: self
        class(image), intent(inout)    :: img
        real, intent(in)               :: x, y
        integer                        :: win(2,2), alloc_stat, inds(2), i, j, lims(2,2)
        real, allocatable              :: w1(:), w2(:)
        real                           :: loc(2), w, pix
        if( img%is_3d() )then
            stop 'Only 2d imags 4 now; interpolate; simple_rev_grid'
        endif
        lims(1,1) = 1
        lims(1,2) = self%boxpd
        lims(2,1) = 1
        lims(2,2) = self%boxpd
        loc(1) = x+self%halfbox
        loc(2) = y+self%halfbox
        win = recwin_2d(loc(1), loc(2), self%harwin)
        allocate( w1(win(1,1):win(1,2)), w2(win(2,1):win(2,2)), stat=alloc_stat )
        call alloc_err("In: interpolate; simple_rev_grid", alloc_stat)
        do i=win(1,1),win(1,2)
            w1(i) = self%wfuns%eval_apod(real(i)-loc(1))
        end do
        do j=win(2,1),win(2,2)
            w2(j) = self%wfuns%eval_apod(real(j)-loc(2))
        end do
        pix = 0.
        do i=win(1,1),win(1,2)
            if( w1(i) == 0. ) cycle
            do j=win(2,1),win(2,2)
                if( w2(j) == 0. ) cycle
                inds = cyci_2d(lims,i,j)
                w = w1(i)*w2(j)
                pix = pix+w*img%get(inds(1),inds(2))
            end do
        end do    
        deallocate(w1,w2)
    end function
    
    !>  \brief  rotation of image by reverse gridding interpolation (from Spider/Sparx)
    function rt( self, img_in, ang ) result( img_rot )
        use enfors_image, only: image
        class(rev_grid), intent(inout) :: self
        class(image), intent(inout)    :: img_in
        real, intent(in)               :: ang
        type(image)                    :: img_rot 
        integer                        :: ldim(3), i, j, h, k, lims(3,2)
        real                           :: mat(2,2), u(2), loc(2)
        ldim = img_in%get_ldim()
        call img_rot%new(ldim, img_in%get_smpd())
        mat = rotmat2d( ang )
        if( img_in%is_ft() )then
            call img_rot%set_ft(.true.)
            lims = img_in%loop_lims(2)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    u(1) = real(h)
                    u(2) = real(k)
                    loc = matmul(u,mat)
                    call img_rot%set_fcomp([h,k,0], self%interpolate_imag(img_in,loc(1),loc(2)))
                end do
            end do
            call img_rot%bwd_ft
        else
            u(1) = -real(ldim(1))/2.
            do i=1,ldim(1)
                u(2) = -real(ldim(2))/2.
                do j=1,ldim(2)
                    loc = matmul(u,mat)
                    call img_rot%set(self%interpolate_real(img_in,loc(1),loc(2)), i, j)
                    u(2) = u(2)+1
                end do
                u(1) = u(1)+1
            end do
        endif
    end function
    
    subroutine test_rev_grid
    
    ! 2 dos: have to check so that rotations are defined as in rtsq
    
        use enfors_image, only: image
        type(image)    :: img, img2, img4grid, img_rot
        type(rev_grid) :: rg
        call img%new([100,100,1], 1.)
        call img2%new([100,100,1], 1.)
        call tester('kb')
        
        contains
        
            subroutine tester( wfun )
                character(len=*), intent(in) :: wfun
                real, allocatable :: corrs(:), res(:)
                character(len=STDLEN) :: wfun_str
                integer :: j
                wfun_str = wfun
                rg = rev_grid(100,wfun_str)
                call img%square(20)
!                call img%ran
!                call img%vis
                img4grid = rg%prep4rev_grid_imag(img)
                img_rot = rg%rt(img4grid, 30.)
                call img_rot%bwd_ft
                call img_rot%clip(img)
                img4grid = rg%prep4rev_grid_imag(img_rot)
                img_rot = rg%rt(img4grid, -30.)
                call img_rot%bwd_ft
                call img_rot%clip(img2)
                call img2%vis
                call img%fsc(img2, res, corrs )
                write(*,*) '>>> TESTING WFUN:', wfun
                do j=1,size(res) 
                    write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
                end do
                write(*,*) 'corr:', img%corr(img2)
            end subroutine
            
    end subroutine

end module simple_rev_grid

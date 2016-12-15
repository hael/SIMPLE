module simple_shalgn
use simple_image,  only: image
use simple_oris,   only: oris
use simple_defs    ! singleton
use simple_params, only: params
implicit none
 
public :: shalgn, test_shalgn
private

type shalgn
    type(image)              :: avg, img      !> average and image
    type(image), allocatable :: refs(:,:)     !> reference images
    class(oris), pointer     :: o_ptr=>null() !> pointer to object that holds orientations
    character(len=STDLEN)    :: stk           !> input stack
    character(len=STDLEN)    :: mode          !> alignment mode = movie|particle          
    real                     :: smpd=0., hp=0., lp=0., msk=0.
    integer                  :: ldim(3)=0, nimgs=0, sh=0, trsstep=1, edge=0
    logical                  :: exists=.false.
  contains
    procedure :: new
    procedure :: align
    procedure :: kill
end type

interface shalgn
    module procedure constructor
end interface

integer, parameter :: MITS=20

contains

    !>  \brief  constructs the shalgn obj
    function constructor( p, o, shalgnmode ) result( self )
        class(params), intent(in)       :: p
        class(oris), intent(in), target :: o
        character(len=*), intent(in)    :: shalgnmode
        type(shalgn)                    :: self
        call self%new(p, o, shalgnmode)
    end function

    !>  \brief  constructs the shalgn obj
    subroutine new( self, p, o, shalgnmode )
        use simple_jiffys, only: alloc_err
        class(shalgn), intent(inout)    :: self
        class(params), intent(in)       :: p
        class(oris), intent(in), target :: o
        character(len=*), intent(in)    :: shalgnmode
        integer                         :: alloc_stat
        ! set constants:
        self%mode    = shalgnmode  ! mode of action
        self%smpd    = p%smpd      ! sampling distance (A)
        self%hp      = p%hp        ! high-pass limit (A)
        self%lp      = p%lp        ! low-pass limit (A)
        self%edge    = p%edge      ! size of edge 4 rectanguar mask with cos falloff
        self%msk     = p%msk       ! radius of spherical mask
        self%stk     = p%stk       ! stack filename
        self%ldim    = p%ldim      ! logical dimension of image
        self%sh      = nint(p%trs) ! integer shift range
        self%trsstep = p%trsstep   ! step size of move
        self%o_ptr   => o          ! pointer to oris
        self%nimgs   = self%o_ptr%get_noris()
        ! allocate what is needed:
        if( self%mode .eq. 'particle' )then
            allocate( self%refs(-self%sh:self%sh,-self%sh:self%sh), stat=alloc_stat )
        else if( self%mode .eq. 'movie' )then
            allocate( self%refs(-1:1,-1:1), stat=alloc_stat )
        else
            stop 'Unsupported mode; new; simple_shalgn'
        endif
        call alloc_err('new; simple_shalgn', alloc_stat)
        call self%avg%new(self%ldim,p%smpd)
        call self%img%new(self%ldim,p%smpd)
        ! indicate existence
        self%exists = .true.
    end subroutine
    
    !>  \brief  the master shiftalignment subroutine
    subroutine align( self )
        use simple_jiffys, only: progress
        class(shalgn), intent(inout) :: self
        integer                      :: j, i
        real                         :: corr, corr_prev, x, y, xnew, ynew
        logical                      :: did_change, refine
        real                         :: movieshifts(-1:1,-1:1,2)
        refine = .false.
        call self%o_ptr%zero_shifts
        call self%o_ptr%set_all('corr',-1.)
        write(*,'(A)') '>>> REFINEMENT'
        do j=1,MITS
            did_change = .false.
            call make_avg
            if( self%mode .eq. 'particle' ) call make_refs(self%sh)
            corr_prev = self%o_ptr%get_avg('corr')
            write(*,'(A)') '>>> ALIGNMENT'
            do i=1,self%nimgs
                call progress(i, self%nimgs)
                x = self%o_ptr%get(i, 'x')
                y = self%o_ptr%get(i, 'y')
                call self%img%read(self%stk, i)
                if( self%mode .eq. 'particle' )then
                    call self%img%mask(self%msk, 'soft') ! apply a soft-edged mask before ft
                else
                    ! Could apply rect mask with cosine fall-off to micrograph
                endif
                call self%img%fwd_ft
                if( self%mode .eq. 'particle' )then
                    call sh_srch_particle(xnew, ynew, corr)
                else
                    call sh_srch_movie(i, xnew, ynew, corr)
                endif
                if( abs(xnew-x) > TINY .or. abs(ynew-y) > TINY )then
                    call self%o_ptr%set(i, 'x',    xnew)
                    call self%o_ptr%set(i, 'y',    ynew)
                    call self%o_ptr%set(i, 'corr', corr)
                    did_change = .true.
                endif
            end do
            corr = self%o_ptr%get_avg('corr')
            write(*,"(1X,A,1X,I3,1X,A,1X,F7.4)") 'Iteration:', j, 'Correlation:', corr
            if( (corr > corr_prev .and. did_change) .or. j < 5 )then
                cycle 
            else
                exit
            endif
        end do

        contains

            subroutine make_avg
                integer :: i, ish, jsh
                write(*,'(A)') '>>> MAKING NEW AVERAGE'
                self%avg = cmplx(0.,0.)
                do i=1,self%nimgs
                    call progress(i, self%nimgs)
                    call self%img%read(self%stk, i)
                    x = self%o_ptr%get(i, 'x')
                    y = self%o_ptr%get(i, 'y')
                    call self%img%fwd_ft
                    call self%img%shift(-x, -y)
                    call self%avg%add(self%img)
                end do
                call self%avg%bwd_ft
                if( self%mode .eq. 'particle' )then
                    ! calculate the rotational average
                    call self%avg%roavg(5.,self%img)
                    self%avg = self%img
                    call self%avg%mask(self%msk, 'soft')
                else
                    ! Could apply rect mask with cosine fall-off to micrograph
                endif
                ! leave average in FTed state
                call self%avg%fwd_ft
            end subroutine
            
            subroutine make_refs( sh, imgind )
                integer, intent(in)           :: sh
                integer, intent(in), optional :: imgind
                integer :: ish, jsh
                real    :: x, y
                if( present(imgind) )then ! movie mode
                    if( imgind == 0 )then
                        x = self%o_ptr%get(imgind, 'x')
                        y = self%o_ptr%get(imgind, 'y')
                    else
                        x = self%o_ptr%get(imgind, 'x')
                        y = self%o_ptr%get(imgind, 'y')
                        x = x+self%o_ptr%get(imgind-1, 'x')
                        y = y+self%o_ptr%get(imgind-1, 'y')
                    endif
                    ! make sure that sh == 1
                    if( sh /= 1 ) stop 'only shift steps of one allowed in movie mode; make_refs; align; simple_shalgn'
                else ! particle mode
                    x = 0.
                    y = 0.
                endif
                do ish=-sh,sh
                    do jsh=-sh,sh
                        self%refs(ish,jsh) = self%avg
                        if( present(imgind) )then ! movie mode
                            movieshifts(ish,jsh,1) = x+real(ish)*real(self%trsstep)
                            movieshifts(ish,jsh,2) = y+real(jsh)*real(self%trsstep)
                            call self%refs(ish,jsh)%shift(movieshifts(ish,jsh,1), movieshifts(ish,jsh,2) , lp_dyn=self%lp)
                        else ! particle mode
                            call self%refs(ish,jsh)%shift(real(ish), real(jsh), lp_dyn=self%lp)
                        endif
                    end do
                end do
            end subroutine
            
            subroutine sh_srch_particle( x, y, corr )
                real, intent(out) :: x, y, corr
                integer :: ish, jsh, xint, yint, xrange(2), yrange(2)
                real :: r
                xint = nint(self%o_ptr%get(i, 'x'))
                yint = nint(self%o_ptr%get(i, 'y'))
                x = real(xint)
                y = real(yint)
                xrange(1) = max(-self%sh,xint-self%trsstep)
                xrange(2) = min( self%sh,xint+self%trsstep)
                yrange(1) = max(-self%sh,yint-self%trsstep)
                yrange(2) = min( self%sh,yint+self%trsstep)
                corr = self%refs(xint,yint)%corr(self%img,self%lp,self%hp)            
                do ish=xrange(1),xrange(2),self%trsstep
                    do jsh=yrange(1),yrange(2),self%trsstep
                        r = self%refs(ish,jsh)%corr(self%img,self%lp,self%hp)
                        if( r > corr )then
                            corr = r
                            x = real(ish)
                            y = real(jsh)
                        endif
                    end do
                end do
            end subroutine
            
            subroutine sh_srch_movie( i, x, y, corr )
                integer, intent(in) :: i
                real, intent(out)   :: x, y, corr
                integer :: ish, jsh
                real    :: r
                ! make references
                call make_refs(1, i)
                ! get previous correlation
                corr = self%o_ptr%get(i, 'corr')
                ! calculate new correlations
                do ish=-1,1
                    do jsh=-1,1
                        r = self%refs(ish,jsh)%corr(self%img,self%lp,self%hp)
                        if( r > corr )then
                            corr = r
                            x = movieshifts(ish,jsh,1)
                            y = movieshifts(ish,jsh,2)
                        endif
                    end do
                end do
            end subroutine

    end subroutine
    
    ! UNIT TEST
    
    subroutine test_shalgn
        use simple_image,  only: image
        use simple_oris,   only: oris
        use simple_params, only: params
        use simple_rnd,    only: ran3
        use simple_procimgfile, only: shift_imgfile
        type(image)  :: img
        type(shalgn) :: shal
        type(oris)   :: o
        type(params) :: p
        integer      :: i
        integer, parameter :: nimgs=500
        character(len=STDLEN), parameter :: fname='test_shalgn_gaus.spi'
        character(len=STDLEN), parameter :: fname_aligned='test_shalgn_gaus_aligned.spi'
        real         :: x, y
        write(*,'(a)') '**info(simple_shalgn_unit_test): testing everything'
        o = oris(nimgs)
        call img%new([100,100,1], 2.)
        do i=1,nimgs
            call img%gauimg(10)
            x = real(nint(ran3()*14.-7.))
            y = real(nint(ran3()*14.-7.))
            call o%set(i, 'x', x)
            call o%set(i, 'y', y)
            call img%fwd_ft
            call img%shift(x, y)
            call img%bwd_ft
            call img%write(fname, i)
        end do
        call o%write('test_shalgn_facit.txt')
        p%stk  = fname
        p%ldim = [100,100,1] 
        p%trs  = 8.
        p%msk  = 45.
        p%smpd = 2.
        p%lp   = 10.
        p%trsstep=2
        call shal%new(p, o, 'particle')
        call shal%align
        call shal%kill
        call o%write('test_shalgn_algndoc.txt')
        call shift_imgfile(fname, fname_aligned, o, 2.)
        write(*,'(a)') 'SIMPLE_SHALGN_UNIT_TEST COMPLETED WITHOUT TERMINAL BUGS ;-)'
        write(*,'(a)') 'PLEASE, INSPECT THE RESULTS'
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(shalgn), intent(inout) :: self
        integer :: ish, jsh
        if( self%exists )then
            do ish=-self%sh,self%sh
                do jsh=-self%sh,self%sh
                    call self%refs(ish,jsh)%kill
                end do
            end do
            call self%avg%kill
            call self%img%kill
        endif
    end subroutine
    
end module simple_shalgn

!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple optimisation module: Fourier shift search
module simple_ft_shsrch
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_image,       only: image
implicit none

public :: ft_shsrch_init, ft_shsrch_reset_ptrs, ft_shsrch_minimize, test_ft_shsrch
private

type(opt_factory)             :: ofac              !< optimizer factory
type(opt_spec)                :: ospec             !< optimizer specification object
class(optimizer), pointer     :: nlopt=>null()     !< pointer to nonlinear optimizer
class(image), pointer         :: reference=>null() !< reference pft
class(image), pointer         :: particle=>null()  !< particle pft
type(image)                   :: tmpimg            !< temporary image for corr calc
real                          :: lp, hp            !< low- and high-pass limits
character(len=:), allocatable :: opt_str           !< optimiser string descriptor
integer                       :: nrestarts=3       !< number of optimisation restarts
real, parameter               :: TOL=1e-4          !< tolerance parameter
integer, parameter            :: MAXITS=30         !< maximum number of iterations

contains
    !> Initialise Fourier shift search
    subroutine ft_shsrch_init( ref, ptcl, lims, lp_in, hp_in, opt, nrestarts_in )
        class(image), target,       intent(in) :: ref, ptcl
        real,                       intent(in) :: lims(2,2), lp_in, hp_in
        character(len=*), optional, intent(in) :: opt
        integer, optional,          intent(in) :: nrestarts_in
        ! set pointers
        reference => ref 
        particle  => ptcl        
        ! set resolution limits
        lp = lp_in
        hp = hp_in
        ! set opt_str & nrestarts
        if( allocated(opt_str) ) deallocate(opt_str)
        if( present(opt) )then
            allocate(opt_str, source=opt)
        else
            allocate(opt_str, source='simplex')
        endif
        nrestarts = 1
        if( present(nrestarts_in) ) nrestarts = nrestarts_in        
        ! make optimizer spec
        call ospec%specify(opt_str, 2, ftol=TOL, gtol=TOL,&
        limits=lims, nrestarts=nrestarts ) ! maxits=MAXITS
        ! set optimizer cost function
        call ospec%set_costfun(ft_shsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
    end subroutine ft_shsrch_init
    
    subroutine ft_shsrch_reset_ptrs( ref, ptcl )
        class(image), intent(in), target :: ref, ptcl
        ! re-set pointers
        reference => ref 
        particle  => ptcl
    end subroutine ft_shsrch_reset_ptrs
    
    function ft_shsrch_cost( vec, D ) result( cost )
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real :: cost
        call tmpimg%copy(particle)
        cost = -reference%corr_shifted(tmpimg, -vec, lp, hp)
    end function ft_shsrch_cost
    
    function ft_shsrch_minimize( prev_corr, prev_shift ) result( cxy )
        real, optional, intent(in) :: prev_corr, prev_shift(2)
        real :: cxy(3), maxshift, maxlim, lims(2,2)
        if( present(prev_shift) )then
            maxshift = maxval(abs(prev_shift))
            maxlim   = maxval(abs(ospec%limits))
            if( maxshift > maxlim )then
                ! re-specify the limits
                lims(:,1) = -maxshift
                lims(:,2) = maxshift
                call ospec%specify(opt_str, 2, ftol=TOL, gtol=TOL,&
                limits=lims, nrestarts=nrestarts ) ! maxits=MAXITS
                ! set optimizer cost function
                call ospec%set_costfun(ft_shsrch_cost)
                ! generate optimizer object with the factory
                call ofac%new(ospec, nlopt)
            endif
            ospec%x = prev_shift
        else
            ospec%x = 0.
        endif
        call nlopt%minimize(ospec, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = ospec%x ! shift
        if( present(prev_corr) )then
            if( abs(cxy(1)-prev_corr) <= TOL )then
                cxy(1)  = prev_corr
                if( present(prev_shift) ) cxy(2:) = prev_shift
            endif
        endif
        call tmpimg%kill
    end function ft_shsrch_minimize
    
    subroutine test_ft_shsrch
        use simple_rnd, only: ran3
        type(image)     :: img_ref, img_ptcl
        real, parameter :: TRS=5.0
        real            :: cxy(3), x, y, lims(2,2)
        integer         :: i
        lims(:,1) = -TRS
        lims(:,2) = TRS
        call img_ref%new([100,100,1], 2.)
        call img_ptcl%new([100,100,1], 2.)
        call img_ref%square(20)
        call img_ref%fwd_ft
        call ft_shsrch_init(img_ref, img_ptcl, lims, 6., 100., 'simplex', 3)
        do i=1,100
            x = ran3()*2*TRS-TRS
            y = ran3()*2*TRS-TRS
            call img_ref%shift(x, y, imgout=img_ptcl)
            cxy = ft_shsrch_minimize()
            if( cxy(1) < 0.995 )then
                stop 'shift alignment failed; test_ft_shsrch; simple_ft_shsrch'
            endif
            call img_ptcl%shift(-cxy(2), -cxy(3))
        end do
        write(*,'(a)') 'SIMPLE_FT_SHSRCH_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ft_shsrch
    
end module simple_ft_shsrch

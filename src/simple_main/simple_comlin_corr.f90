!==Class simple_comlin_corr
!
! simple_comlin_corr implements a decorator pattern for selecting common line correlators. 
! The strategy is to have a non-instantatiable class (singleton) in which pointers are assigned to 
! data required for calculating common line correlations. Note that in several instances the ponited 
! to variables are modified by the method calls (see comments in the code). The code is distributed 
! with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification 
! is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2009-06-12.
! 
!==Changes are documented below
!
!* deugged and incorporated in the _SIMPLE_ library, HE 2009-06-25
!* refurnished according to the new OO standard, HE 2012-06-15
!
module simple_comlin_corr
use simple_build, only: build
use simple_ori,   only: ori
implicit none

type(ori), private               :: o_old
class(build), pointer, private   :: bp=>null()
integer, pointer, private        :: pbootarr(:)=>null(), pstatearr(:)=>null()
integer, pointer, private        :: pstate=>null(), pptcl=>null()
integer                          :: varwhich
real, pointer, private           :: plp_dyn=>null()
logical, private                 :: assoc=.false., decor=.false.

interface pcorr_comlin
    module procedure pcorr_comlin_1
    module procedure pcorr_comlin_2
    module procedure pcorr_comlin_3
    module procedure pcorr_comlin_4
end interface

contains

    !>  \brief  initializes the comlin_corr singleton
    subroutine comlin_corr_init( b, ptcl, lp_dyn )
        use simple_build, only: build
        class(build), intent(in), target :: b
        integer, intent(in), target      :: ptcl
        real, intent(in), target         :: lp_dyn
        bp => b
        pptcl => ptcl
        plp_dyn => lp_dyn
        assoc = .true.
    end subroutine

    !>  \brief  decorates additional functionality to the comlin_corr singleton
    subroutine comlin_corr_decor( bootarr_in, statearr_in, state_in )
        integer, intent(in), optional, target :: bootarr_in(:), statearr_in(:), state_in
        if( assoc )then
            if( present(bootarr_in) )  pbootarr => bootarr_in
            if( present(statearr_in) )then
                if( present(state_in) )then
                    pstatearr => statearr_in
                    pstate    => state_in
                else
                    write(*,'(a)') 'Optional argument state_in required for operation!'
                    write(*,'(a)') 'decor_comlin_corr; simple_comlin_corr'
                endif
            endif
            decor = .true.
         else
            write(*,'(a)') 'Associate comlin_corr before decorating additional functionality!'
            write(*,'(a)') 'decor_comlin_corr; simple_comlin_corr'
            stop
        endif
    end subroutine

    !>  \brief  is for nullifying the pointers associated by decor_comlin_corr
    subroutine comlin_corr_undecor
        pstate    => null()
        pbootarr  => null()
        pstatearr => null()
        decor     = .false.
    end subroutine

    !>  \brief  is for nullifying the bootarr pointer
    subroutine comlin_corr_undecor_boot
        pbootarr => null()
    end subroutine
    
    !>  \brief  is for calculating the per-ptcl continuous common line correlation coefficient
    !!          modifies the clins data structure (make sure to store the one_ptcl contribution)
    function pcorr_comlin_1() result( corr )
        real :: corr
        corr = -1.
        if( assoc )then
            if( associated(pbootarr) )then
                if( associated(pstatearr) )then
                    corr = bp%clins%pcorr( pptcl, plp_dyn, pbootarr, pstatearr, pstate )
                else
                    corr = bp%clins%pcorr( pptcl, plp_dyn, pbootarr )
                endif
            else if( associated(pstatearr) )then
                corr = bp%clins%pcorr( pptcl, plp_dyn, pstatearr, pstate )
            else
                corr = bp%clins%pcorr( pptcl, plp_dyn )
            endif
        else
            write(*,'(a)') 'Associate comlin_corr pointers before calculating correlation!'
            write(*,'(a)') 'pcorr_comlin; simple_comlin_corr'
            stop
        endif
    end function
    
    !>  \brief  is for evaluating new orientations in a
    function pcorr_comlin_2( e1, e2, e3, x, y ) result( corr )
        real, intent(in) :: e1, e2, e3
        real, intent(in), optional :: x, y
        real :: x_here, y_here, corr
        x_here = 0.
        y_here = 0.
        if( present(x) ) x_here = x
        if( present(y) ) y_here = y
        ! update orientations structure
        o_old = bp%a%get_ori(pptcl)
        call bp%a%set_euler(pptcl, [e1,e2,e3])
        call bp%a%set(pptcl, 'x', x_here)
        call bp%a%set(pptcl, 'y', y_here)
        ! calculate correlation
        corr = pcorr_comlin()
        ! put old orientation back
        call bp%a%set_ori(pptcl,o_old)
    end function
    
    !>  \brief  is for evaluating new orientations in a
    function pcorr_comlin_3( var ) result( corr )
        real, intent(in) :: var
        real :: corr
        ! update orientations structure
        o_old = bp%a%get_ori(pptcl)
        select case( varwhich )
            case(1)
                call bp%a%e1set( pptcl, var )
            case(2)
                call bp%a%e2set( pptcl, var )
            case(3)
                call bp%a%e3set( pptcl, var )
            case(4)
                call bp%a%set( pptcl, 'x', var )
            case(5)
                call bp%a%set( pptcl, 'y', var )
            case DEFAULT
                write(*,*) 'Unknown case: ', varwhich
                stop 'simple_comlin_corr :: pcorr_comlin_3'
        end select
        ! calculate correlation
        corr = pcorr_comlin()
        ! put old orientation back
        call bp%a%set_ori(pptcl,o_old)
    end function
    
    !>  \brief  is for calculating the per-ptcl continuous pairwise common line correlation coefficient
    function pcorr_comlin_4( iptcl, jptcl ) result( corr )
        integer, intent(in) :: iptcl, jptcl
        real :: corr
        corr = -1.
        if( assoc )then
            corr = bp%clins%pcorr(iptcl, jptcl, plp_dyn)
        else
            write(*,'(a)') 'Associate comlin_corr pointers before calculating correlation!'
            write(*,'(a)') 'pcorr_comlin; simple_comlin_corr'
            stop
        endif
    end function
    
    !>  \brief  is for monitoring the score in sgd
    function pscore_comlin_sgd( vec, i, N, D ) result( score )
        integer, intent(in) :: i, N, D
        real, intent(in)    :: vec(N,D)
        real                :: score
        pptcl = i ! pointed to ptcl value is modified
        score = pcorr_comlin_2(vec(i,1), vec(i,2), vec(i,3), vec(i,4), vec(i,5))
    end function
    
    !>  \brief  is for calcuating the gradient of the continuous negative common line correlation coefficient 
    !!          given the ptcl _i_, the variable selector _j_, and the solution vector.
    !!          Used in comlin_srch for the stochastic gradient descent-based alignment
    function pgrad_comlin_sgd( vec, i, j, N, D ) result( grad )
        use simple_math, only: numderiv
        integer, intent(in) :: N, D, i, j
        real, intent(in)    :: vec(N,D)
        real, parameter     :: eh=5., xh=1.
        real                :: grad, err
        pptcl    = i ! pointed to ptcl value is modified
        varwhich = j ! to indicate which variable is being differentiated
        if( varwhich <= 3 )then
            call numderiv( pcorr_comlin_3, vec(i,j), eh, err, grad )
        else
            call numderiv( pcorr_comlin_3, vec(i,j), xh, err, grad )
        endif
    end function
    
    !>  \brief  is for calcuating the continuous negative common line correlation coefficient 
    !!          given the ptcl _i_, the integer alignment vector _vec_ and its dimensions (_N_,_L_).
    !!          Used in comlin_srch for the annealing
    function pcost_comlin_sa( vec, i, N, L ) result( cost )
        integer, intent(in) :: i, N, L, vec(N,L)
        real                :: cost
        pptcl = i ! pointed to ptcl value is modified
        cost = -pcorr_comlin(bp%e%e1get(vec(i,1)), bp%e%e2get(vec(i,1)), 0.)
    end function
    
    !>  \brief  is for calcuating the continuous negative common line correlation coefficient 
    !!          given the alignment vector _vec_ and its size _D_. The function is used by
    !!          simple_rfree_search for reference-free multiple conformations alignment.
    !!          modifies the clins data structure (make sure to store the one_ptcl contribution)
    function pcost_comlin_de( vec, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real                :: cost
        cost = 1.
        if( decor )then
            pstate = int(vec(6)) ! pointed to state variable is modified
            cost = -pcorr_comlin(vec(1), vec(2), vec(3), vec(4), vec(5))
        else
            write(*,'(a)') 'Decorate statearr functionality before calculating correlation!'
            write(*,'(a)') 'pcorr_comlin_de; simple_comlin_corr'
            stop
        endif
    end function
    
end module simple_comlin_corr

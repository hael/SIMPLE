!==Class simple_comlin_srch
!
! The simple_comlin_srch singleton provides the reference free alignment/heterogeneity 
! analysis functionality in the _SIMPLE_ library. The code is distributed with the hope 
! that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification 
! is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-06-11.
module simple_comlin_srch
use simple_comlin_corr ! singleton
use simple_build,      only: build
use simple_sa_opt,     only: sa_opt
use simple_ran_tabu,   only: ran_tabu
use simple_jiffys,     only: progress, alloc_err
use simple_rnd,        only: ran_iarr
implicit none

public :: comlin_srch_init, comlin_srch_corr, comlin_srch_corr_distr, comlin_srch_set_states,&
comlin_srch_anneal, comlin_srch_state, comlin_srch_symaxis
private

class(build), pointer :: bp=>null()      !< pointer to builder
type(sa_opt)          :: metro           !< simulated annealing object 
type(ran_tabu)        :: rt              !< for making non-equal random numbers
integer, allocatable  :: isol(:,:)       !< labeling solution for the annealing
integer, allocatable  :: bootarr(:)      !< random sample array
integer, allocatable  :: statearr(:)     !< state assignment array
integer, pointer      :: ptcl=>null()    !< ptcl index
real, pointer         :: dynlim=>null()  !< dynamic lowpass limit
integer               :: astep           !< in-plane angle step size for symmetry search
integer               :: nptcls=0        !< nr of particles
integer               :: nspace=0        !< nr of projection directions
integer               :: N_boot=100      !< random sample size
integer               :: state=0         !< state assignment variable
integer               :: nbest=0         !< de pop size
real                  :: ares=0.         !< angular resolution
real                  :: trs=0.          !< origin shift half interval
real                  :: eullims(3,2)=0. !< Euler limits
real                  :: olims(6,2)=0.   !< optimization limits
real                  :: hp=0.           !< fixed high-pass limit
real                  :: lp=0.           !< fixed low-pass limit

real, parameter    :: cost_err = 0.0001, T_initial = 1000000., T_update = 0.9, T_lowlim = 0.5
real, parameter    :: trs_lrate = 0.05, eul_lrate = 0.1
real, parameter    :: lrates(5) = [eul_lrate,eul_lrate,eul_lrate,trs_lrate,trs_lrate]
integer, parameter :: T_level_max_rearr = 500
logical, parameter :: limvars(5) = [.false.,.false.,.false.,.true.,.true.]

contains

    !>  \brief  is a constructor
    subroutine comlin_srch_init( b, p )
        use simple_build,  only: build
        use simple_params, only: params
        class(build), intent(in), target  :: b
        class(params), intent(in), target :: p
        integer                           :: alloc_stat
        logical                           :: cyhere(6)
        ! allocate:
        allocate( statearr(p%nptcls), isol(nptcls,1), stat=alloc_stat )
        call alloc_err('init_comlin_srch', alloc_stat) 
        statearr = 1
        ! set constants & pointers:
        astep  = p%astep
        ares   = p%ares
        nptcls = p%nptcls
        nspace = p%nspace
        nbest  = p%nbest
        hp     = p%hp
        lp     = p%lp
        trs    = p%trs
        bp     => b
        ptcl   => p%ptcl
        dynlim => p%lp_dyn
        ! make comlin_corr functionality:
        call comlin_corr_init( b, ptcl, dynlim )
        ! make random number functionality
        rt = ran_tabu( p%nptcls )
        ! make annealing object
        call metro%new([nspace], [.true.], nptcls, 1, pcost_comlin_sa, update)
        call metro%set_statefun(statefun)
        ! make constants for differential evolution:
        cyhere     = .false.
        eullims    = p%eullims
        olims(1,1) = p%eullims(1,1)
        olims(1,2) = p%eullims(1,2)
        olims(2,1) = p%eullims(2,1)
        olims(2,2) = p%eullims(2,2)
        olims(3,1) = 0.
        olims(3,2) = 0.0001
        olims(4,1) = -p%trs
        olims(4,2) = p%trs
        olims(5,1) = -p%trs
        olims(5,2) = p%trs
        if( p%nstates > 1 )then
            olims(6,1) = 1.
            olims(6,2) = real(p%nstates)+0.9999
        else
            olims(6,1) = 1.
            olims(6,2) = 1.9999
        endif
    end subroutine
    
    ! GETTER/SETTERS
    
    !>  \brief  is for calculating the joint common line correlation,
    !!          parameterization is with respect to state
    function comlin_srch_corr( s ) result( corr )
        integer, intent(in), optional :: s
        real :: corr
        integer :: ss, i
        ss = 1
        if( present(s) )then
            ss = s
            ! decorate statearr functionality to the comlin_corr object
            call comlin_corr_decor( statearr_in=statearr, state_in=state )
            state = ss
        endif
        ! undecorate any bootstrap decoration
        call comlin_corr_undecor_boot
        corr = 0.
        do i=1,nptcls
            ptcl = i
            corr = corr+pcorr_comlin()
        end do
        corr = corr/real(nptcls)
    end function

    !>  \brief  is for calculating all individual common line correlations
    !!           (for statistical analysis of the distribution)
    function comlin_srch_corr_distr() result( corrs )
        real, allocatable :: corrs(:)
        integer           :: i, alloc_stat
        allocate( corrs(nptcls), stat=alloc_stat )
        call alloc_err('simple_comlin_srch; comlin_srch_corr_distr', alloc_stat) 
        ! decorate statearr functionality to the comlin_corr object
        call comlin_corr_decor(statearr_in=statearr, state_in=state)
        ! undecorate any bootstrap decoration
        call comlin_corr_undecor_boot
        do i=1,nptcls
            ptcl = i
            corrs(i) = pcorr_comlin()
        end do
    end function
    
    !>  \brief  is for setting the state array
    subroutine comlin_srch_set_states
        integer :: i
        do i=1,nptcls
            statearr(i) = nint(bp%a%get(i, 'state'))
        end do
    end subroutine
    
    ! OPTIMIZERS
    
    !>  \brief  is for finding the symmetry axis give an aligned set of images
    subroutine comlin_srch_symaxis( orientation_best, doprint )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(ori), intent(inout) :: orientation_best
        logical, intent(in)       :: doprint
        integer    :: i, j, k
        type(ori)  :: orientation
        type(oris) :: a_copy
        real       :: corr, corr_best
        integer    :: ntot, cnt, lims(3,2)
        if( doprint ) write(*,'(A)') '>>> GLOBAL SYMMETRY AXIS SEARCH'
        a_copy    = bp%a
        corr_best = -1.
        ntot      = 24624
        cnt       = 0
        lims(1,1) = 0
        lims(1,2) = 359
        lims(2,1) = 0
        lims(2,2) = 180
        lims(3,1) = 0
        lims(3,2) = 359
        call orientation%new
        do i=lims(1,1),lims(1,2),10
            call orientation%e1set(real(i))
            do j=lims(2,1),lims(2,2),10
                call orientation%e2set(real(j))
                do k=lims(3,1),lims(3,2),10
                    cnt = cnt+1
                    if( doprint ) call progress(cnt, ntot)
                    call orientation%e3set(real(k))
                    bp%a = a_copy
                    call bp%a%rot(orientation)
                    call bp%se%apply2all(bp%a)
                    corr = comlin_srch_corr()
                    call orientation%set('corr',corr)
                    if( corr > corr_best )then
                        corr_best = corr
                        orientation_best = orientation
                    endif
                end do
            end do
        end do
        if( doprint )then
            write(*,'(A)') '>>> FOUND FIRST SYMMETRY AXIS ORIENTATION'
            call orientation_best%print
            write(*,'(A)') '>>> REFINED SYMMETRY AXIS SEARCH'
        endif
        lims(1,1) = max(nint(orientation_best%e1get()-10.),lims(1,1))
        lims(1,2) = min(nint(orientation_best%e1get()+10.),lims(1,2))
        lims(2,1) = max(nint(orientation_best%e2get()-10.),lims(2,1))
        lims(2,2) = min(nint(orientation_best%e2get()+10.),lims(2,2))
        lims(3,1) = max(nint(orientation_best%e3get()-10.),lims(3,1))
        lims(3,2) = min(nint(orientation_best%e3get()+10.),lims(3,2))
        ntot      = (lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)*(lims(3,2)-lims(3,1)+1)
        cnt       = 0 
        do i=lims(1,1),lims(1,2)
            call orientation%e1set(real(i))
            do j=lims(2,1),lims(2,2)
                call orientation%e2set(real(j))
                do k=lims(3,1),lims(3,2)
                    cnt = cnt+1
                    if( doprint ) call progress(cnt, ntot)
                    call orientation%e3set(real(k))
                    bp%a = a_copy
                    call bp%a%rot(orientation)
                    call bp%se%apply2all(bp%a)
                    corr = comlin_srch_corr()
                    call orientation%set('corr',corr)
                    if( corr > corr_best )then
                        corr_best = corr
                        orientation_best = orientation
                    endif
                end do
            end do
        end do
        if( doprint )then
            write(*,'(A)') '>>> FOUND REFINED SYMMETRY AXIS ORIENTATION'
            call orientation_best%print
        endif
        bp%a = a_copy
        call a_copy%kill
    end subroutine
    
    !>  \brief  is the simulated annealing based algorithm "Reference-free 3D Alignment in a Discrete angular space"
    !!          Elmlund H. et al. "A new cryo-EM single-particle _ab_ _initio_ reconstruction method visualizes 
    !!          secondary structure elements in an ATP-fueled AAA+ motor" J. Mol. Biol. 2008 Jan 25;375(4):934-47.
    subroutine comlin_srch_anneal
        integer :: i
        real    :: cost
        write(*,'(A)') '>>> REFERENCE-FREE ALIGNMENT IN A DISCRETE ANGULAR SPACE'
        ! anneal, the cost function takes care of modifying the ptcl var
        call metro%minimize(cost_err, T_initial, T_update, T_level_max_rearr, T_lowlim, isol, cost )
        do i=1,nptcls
            call bp%a%set_euler(i, bp%e%get_euler(isol(i,1)))
            call bp%a%set(i, 'x',       0.)
            call bp%a%set(i, 'y',       0.)
            call bp%a%set(i, 'corr', -cost)
            call bp%a%set(i, 'state',   1.)
        end do          
    end subroutine
    
    !>  \brief  is for assinging conformational states to aligned ptcls by GRASP
    subroutine comlin_srch_state( slim )
        integer, intent(in) :: slim
        integer :: it, imax(1), statearr_prev(nptcls), statearr_best(nptcls), i, j, s
        real :: euls(3), x, y, corrs(slim), corr_this, corr_prev, corr_best
        if( slim > 1 )then
            ! init best corr
            corr_best = -1.
            ! decorate statearr functionality to the comlin_corr object
            call comlin_corr_decor( statearr_in=statearr, state_in=state )
            ! undecorate any bootstrap decoration
            call comlin_corr_undecor_boot
            do j=1,5
                write(*,'(A,1X,I2,A)') '>>> STATE ASSIGNMENT, ROUND:', j, ' OF 5'
                ! assign random states
                call ran_iarr( statearr, slim )
                ! init
                corr_this = -1.
                it = 0
                do ! local search
                    it = it+1
                    corr_prev = corr_this
                    statearr_prev = statearr
                    corr_this = 0.
                    do i=1,nptcls
                        ptcl = i ! modifying the ptcl var in comlin_corr
                        do s=1,slim
                            state = s ! modifying the state var in comlin_corr
                            euls = bp%a%get_euler(ptcl)
                            x = bp%a%get(ptcl,'x')
                            y = bp%a%get(ptcl,'y')
                            corrs(state) = pcorr_comlin( euls(1), euls(2), euls(3), x, y )
                        end do
                        imax = maxloc(corrs)
                        corr_this = corr_this+corrs(imax(1))
                        statearr(i) = imax(1)    
                    end do
                    corr_this = corr_this/real(nptcls)       
                    if( corr_this > corr_prev .or. it <= 30 ) then
                        write(*,"(1X,A,1X,I7,1X,A,1X,F7.4)") 'Iteration:', it, 'Correlation:', corr_this
                        cycle
                    else
                        exit
                    endif
                end do
                if( corr_this > corr_best )then
                    corr_best = corr_this
                    statearr_best = statearr_prev
                endif
            end do
            do i=1,nptcls
                call bp%a%set(i, 'state', real(statearr_best(i)))
                statearr(i) = statearr_best(i)
            end do
         endif
    end subroutine
    
    ! PRIVATE STUFF

    !>  \brief  is for updating bp%a when a configuration is accepted in RAD
    subroutine update( vec, i, L ) 
        integer, intent(in) :: i, L, vec(L)
        call bp%a%set_ori(i, bp%e%get_ori(vec(1)))
    end subroutine
    
    !>  \brief  is for updating bp%a when a configuration is accepted in the stochastic gradient descent
    subroutine update2( vec, i, D )
        integer, intent(in) :: D, i
        real, intent(in)    :: vec(D)
        call bp%a%set_euler(i, vec(:3))
        call bp%a%set(i, 'x', vec(4))
        call bp%a%set(i, 'y', vec(5))
    end subroutine
    
    !>  \brief  is for temperature dependent update of the lowpass limit used for calculating the common line
    !!          correlation. Also updates the bootstrap sample during annealing. The strategy is to have a small 
    !!          sample size and low resolution cost function at high temperature and increase the sample size 
    !!          during the annealing simultaneously with increased resolution (to speed things up and improve 
    !!          the accuracy of the low temperature alignment)
    subroutine statefun( T, T_init )
        real, intent(in) :: T, T_init
        real             :: fract, resfact
        integer          :: alloc_stat
        fract = log(T)/log(T_init)
        if( fract > 0. )then
            resfact = fract*(hp/2.-lp)
        else
            resfact = 0.
        endif
        dynlim = lp+resfact
        ! bootstrap sample update
        if( T > 1. ) then
            N_boot = min(nptcls,min(100,int((1./log(T))*100.)))
        else
            N_boot = min(nptcls,100)
        endif
        if( allocated(bootarr) ) deallocate(bootarr)
        allocate( bootarr(N_boot), stat=alloc_stat )
        call alloc_err( 'statefun; simple_comlin_srch', alloc_stat)
        call rt%reset
        call rt%ne_ran_iarr(bootarr) ! bootstrap sample update
        ! decorate bootstrap estimation functionality to the comlin_corr object:
        call comlin_corr_decor( bootarr_in=bootarr )
    end subroutine

end module simple_comlin_srch

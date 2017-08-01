module simple_oasis_srch
use simple_ori,         only: ori
use simple_oris,        only: oris
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_build,       only: build
use simple_params,      only: params
use simple_image,       only: image
use simple_defs         ! singleton
use simple_cmdline      ! singleton
implicit none

public ::
private

class(params), pointer    :: pp=>null()          !< pointer 2 params
class(build), pointer     :: bp=>null()          !< pointer 2 build
logical, parameter        :: debug=.true.        !< debug flag
type(opt_factory)         :: ofac                !< optimizer factory
type(opt_spec)            :: ospec               !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()       !< pointer to nonlinear optimizer
integer                   :: nrefs_eval=0        !< number of references evaluated
integer                   :: state_glob=1        !< previous state parameter
type(ori)                 :: o_best              !< best orientation identified by the search
type(oris)                :: oriset              !< set of probabilistic orientations
real                      :: dfx_prev = 0., dfy_prev = 0., angast_prev = 0.

contains
    
    !>  \brief  is a constructor
    subroutine init_oasis_srch( b, p, how2srch )
        use simple_jiffys, only: alloc_err
        class(build),  intent(in), target :: b
        class(params), intent(in), target :: p
        character(len=*), intent(in) :: how2srch
        integer, intent(in)          :: nstates
        integer :: s, alloc_stat
        ! initialize
        bp        => b
        pp        => p
        nrefs_eval = 0
        state_glob = 1
        
        ! make optimiser spec
        if( how2srch .eq. 'soft' )then
            call ospec%specify(p%opt, p%ndim, ftol=1e-4, gtol=1e-4, npeaks=IMPORTANT,&
            limits=p%optlims(:p%ndim,:), cyclic=p%cyclic(:p%ndim), nrestarts=1, verbose=.true.)
        else if( how2srch .eq. 'hard' )then
            call ospec%specify(p%opt, p%ndim, ftol=1e-4, gtol=1e-4,&
            limits=p%optlims(:p%ndim,:), cyclic=p%cyclic(:p%ndim), nrestarts=3, verbose=.true.)
        else
            stop 'unknown flag value how2srch; simple_oasis_srch::new'
        endif
        ! set optimiser cost function
        call ospec%set_costfun(cost_volimg_oasis_srch)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        if( debug ) write(*,*) '*** created new oasis_srch object'
    end subroutine
    
    !>  \brief  is to get the number of probabilistic orientations
    function noris_oasis_srch() result( noris )
        integer :: noris
        noris = oriset%get_noris()
    end function
    
    !>  \brief  is to get one probabilistic orientation
    function get_ori_oasis_srch( iori ) result( o )
        integer, intent(in) :: iori
        type(ori) :: o
        if( iori < 1 .or. iori > oriset%get_noris() )then
            write(*,*) 'iori is out of range'
            write(*,*) 'iori:', iori, 'noris:', oriset%get_noris()
            stop 'In: simple_oasis_srch::get_ori'
        endif
        o = oriset%get_ori(iori)
    end function
    
    !>  \brief  is for generating the abs(CTF) image needed for correlation matching
    subroutine prep_ctfabsimg_oasis_srch( iptcl, o )
        use simple_math, only: euclid
        integer, intent(in)       :: iptcl
        class(ori), intent(inout) :: o
        real                      :: x, y, dfx, dfy, angast, dist
        integer                   :: s, j
        if( pp%ctf .ne. 'no' )then
            ! parse ori
            x = o%get('x')                            
            y = o%get('y')
            s = nint(o%get('state'))
            if( pp%ctfmode .eq. 'astig' )then ! astigmatic CTF
                dfx = bp%a%get(iptcl,'dfx')
                dfy = bp%a%get(iptcl,'dfy')
                angast = bp%a%get(iptcl,'angast')
            else if( pp%ctfmode .eq. 'noastig' )then
                dfx = bp%a%get(iptcl,'dfx')
                dfy = dfx
                angast = 0.
            else
                stop 'unknown ctf mode; simple_oasis_srch::prep_ctfabsimg'
            endif
            dist = euclid([dfx,dfy,angast],[dfx_prev,dfy_prev,angast_prev])  
            if( dist < 0.001 )then
                ! CTF parameters are the same as for the previous particle & no update is needed
            else
                ! CTF parameters have changed and the abs(CTF) image needs to be updated
                call bp%tfun%ctf2img(bp%img_ctf, dfx, 'abs', dfy, angast)
            endif
            dfx_prev    = dfx
            dfy_prev    = dfy
            angast_prev = angast
        endif
        if( debug ) write(*,*) '*** finished oasis_srch::prep_ctfabsimg'
    end subroutine
    
    !>  \brief  calculates the vol vs. img correlation for the input orientation/state
    function corr_volimg_oasis_srch( orientation, state ) result( corr )
        class(ori), intent(inout)     :: orientation
        integer, intent(in), optional :: state
        real                          :: x, y, corr
        ! set state variable
        state_glob = 1
        if( present(state) )then
            state_glob = state
        else
            if( pp%nstates > 1 ) state_glob = nint(orientation%get('state'))
        endif
        ! project volume (OpenMP para)
        call bp%proj%fproject(bp%refvols(state_glob), orientation, bp%img_ref, lp_dyn=pp%lp)
        x = orientation%get('x')
        y = orientation%get('y')
        if( pp%ctf .ne. 'no' )then
            ! multiply the reference section with the abs(CTF) (OpenMP para)
            call bp%img_ref%mul(bp%img_ctf, pp%lp)
        endif
        ! shift the reference section (OpenMP para)
        call bp%img%shift(x, y, imgout=bp%img_sh, lp_dyn=pp%lp)
        ! correlate image with reference section (OpenMP para)
        ! correlation function NOT normalized
!         corr = bp%img_ref%corr_memoized(bp%img_sh, lp_dyn=pp%lp, hp_dyn=pp%hp)
        corr = bp%img_ref%corr(bp%img_sh, lp_dyn=pp%lp, hp_dyn=pp%hp)
    end function
    
    !>  \brief  calculates the vol vs. img cost for input vector
    !!          (remember to set state beforehand via p%ori_glob)
    function corr_volimg_oasis_srch( vec, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real                :: cost
        type(ori)           :: o
        call o%new
        call o%set_euler(vec(1:3))
        call o%set('x',vec(4))
        call o%set('x',vec(5))
        cost = -corr_volimg_oasis_srch(glob_obj%pp%ori_glob)
    end function
    
    !>  \brief  is the OASIS search routine
    subroutine srch_hard_oasis_srch( self, trial )
        use simple_jiffys, only: int2str
        use simple_math,   only: rad2deg
        class(oasis_srch), intent(inout) :: self
        type(ori), intent(inout)         :: trial
        real    :: x(5), x_old, y_old, lowest_cost, corr
        real    :: sumd, avgd, sdevd, mind, maxd
        integer :: ipeaks
        logical :: found_better
        ! extract trial solution & communicate it to OASIS via ospec%x
        x(1:3)  = trial%get_euler()
        x(4)    = trial%get('x')
        x(5)    = trial%get('y')
        x_old   = x(4) 
        y_old   = x(5)
        ospec%x = x
        if( debug ) call trial%print
        ! optimize
        call nlopt%minimize(ospec, lowest_cost)
        if( debug ) write(*,*) '*** oasis found optimum: ', ospec%x, lowest_cost
        ! gather stats
        nrefs_eval = ospec%nevals
        if( debug ) write(*,*) '*** nr of references evaluated: ', nrefs_eval
        call oriset%new(1)
        call oriset%set_ori(1,trial) ! transfer of parameters
        call oriset%set_euler(1,  ospec%x(1:3))
        call oriset%set(1, 'x',   ospec%x(4))
        call oriset%set(1, 'y',   ospec%x(5))
        call oriset%set(1,'corr', -lowest_cost)
        call oriset%set(1,'w',    1.)
        call oriset%set(1,'nmin', 1.)
        if( debug ) write(*,*) '***********************'
        if( debug ) call oriset%print(1)
        if( debug ) write(*,*) '*** finished oasis_srch::srch_hard'
    end subroutine
    
    !>  \brief  is the OASIS search routine
!     subroutine srch( self, trial )
!         use simple_jiffys, only: int2str
!         use simple_math,   only: rad2deg
!         class(oasis_srch), intent(inout) :: self
!         type(ori), intent(inout)         :: trial
!         real    :: w, wsum, x(5), x_old, y_old, lowest_cost, corr
!         real    :: sumd, avgd, sdevd, mind, maxd
!         integer :: ipeaks
!         logical :: found_better
!         ! extract trial solution & communicate it to OASIS via ospec%x
!         x(1:3) = trial%get_euler()
!         x(4)   = trial%get('x')
!         x(5)   = trial%get('y')
!         x_old  = x(4)
!         y_old  = x(5)
!         ospec%x = x
!         ! optimize to find IMPORTANT number of peaks
!         call oas%minimize(ospec, lowest_cost)
!         if( debug ) write(*,*) '*** oasis found optimum: ', ospec%x, lowest_cost
!         ! gather stats
!         nrefs_eval = ospec%nevals
!         if( debug ) write(*,*) '*** nr of references evaluated: ', nrefs_eval
!         if( debug ) write(*,*) '*** nr of minimas sought: ',       ospec%npeaks
!         if( debug ) write(*,*) '*** nr of minimas found: ',        ospec%peakcnt
        ! update oriset (set orientations, correlations & weights)
!         if( ospec%peakcnt == 0 )then
!
!             ospec%peakcnt = 1
!
!             call oriset%new(1)
!             call oriset%set_ori(1,trial)
!             call oriset%set(1,'corr',-lowest_cost)
!             call oriset%set(1,'w',1.)
!             call oriset%set(1,'nmin',0.)
!         else
!             found_better = .true.
!             call oriset%new(ospec%peakcnt)
!             wsum = 0.
!             do ipeaks=1,ospec%peakcnt
!                 ! transfer of parameters
!                 call oriset%set_ori(ipeaks, trial)
!                 ! set peak information
!                 corr = -ospec%peaks(ipeaks, pp%ndim+1)
!                 call oriset%set(ipeaks, 'corr', corr)
!                 call oriset%set_euler(ipeaks, ospec%peaks(ipeaks, 1:3))
!                 ! shifts must be obtained by vector addition
!                 call oriset%set(ipeaks, 'x', x_old-ospec%peaks(ipeaks, 4)) ! update x-shift, revshsgn to fit shift convention
!                 call oriset%set(ipeaks, 'y', y_old-ospec%peaks(ipeaks, 5)) ! update y-shift, revshsgn to fit shift convention
!                 w = 0.
!                 if( corr > 0. ) w = exp(corr)
!                 call oriset%set(ipeaks, 'ow', w)
!                 wsum = wsum+w
!             end do
!             ! normalise weights
!             if( wsum > 0. )then
!                 do ipeaks=1,ospec%peakcnt
!                     w = oriset%get(ipeaks, 'ow')
!                     call oriset%set(ipeaks, 'ow', w/wsum)
!                 end do
!             endif
!         endif
        ! set ouput statistics
!         o_best = oriset%get_ori(ospec%peakcnt)
!         call oriset%set(ospec%peakcnt, 'dist', rad2deg(trial.euldist.o_best))
!         call oriset%set(ospec%peakcnt, 'nmin', real(ospec%peakcnt))
!         if( found_better )then
!             call oriset%set(ospec%peakcnt, 'nbetter', 1.)
!         else
!             call oriset%set(ospec%peakcnt, 'nbetter', 0.)
!         endif
!         if( ospec%peakcnt >= 2 )then
!             call oriset%diststat(sumd, avgd, sdevd, mind, maxd)
!             call oriset%set(ospec%peakcnt, 'sdev', sdevd)
!         else
!             call oriset%set(ospec%peakcnt, 'sdev', 0.)
!         endif
!         if( debug ) write(*,*) '*** finished oasis_srch::prep_ctfabsimg'
!     end subroutine

end module

! OASIS IMPORTANCE SAMPLING ROUTINES

!>  \brief  initialization step for the OASIS minimization
! subroutine oasis_init( self, ospec, sdevs, probs )
!     use simple_rnd, only: ran3
!     class(pori), intent(inout)  :: self
!     class(opt_spec), intent(in) :: ospec
!     real, intent(in), optional  :: sdevs(:), probs(:)
!     real, allocatable :: x(:)
!     logical :: arezero(ngaup)
!     integer :: i
!     ! extract the previous solution
!     x = extract_prev()
!     ! test if best point is set
!     arezero = .false.
!     do i=1,ngaup
!         if( x(i) == 0. ) arezero(i) = .true.
!     end do
!     ! generate initial vector
!     if( all(arezero) )then
!         do i=1,ngaup
!             ! initialize each variable by randomized bounds
!             x(i) = ospec%limits(i,1)+ran3()*(ospec%limits(i,2)-ospec%limits(i,1))
!         end do
!         ! update the parameters in ori object
!         call update_ori_part(x)
!     endif
!     if( sum(sdevs) > 0. )then
!         ! sdevs are already read
!     else
!         ! initialize by setting all the sdevs to the half the interval
!         do i=1,ngaup
!             sdevs(i) = (ospec%limits(i,2)-ospec%limits(i,1))/2.
!             if( i > 3 )then
!                 if( sdevs(i) < 1. ) sdevs(i) = 1.
!             endif
!         end do
!     endif
!     if( present(sdevs) )then
!         do i=1,size(sdevs)
!             sdevs(i) = sdevs(i)
!         end do
!     endif
!     if( allocated(probs) )then
!         if( present(probs) )then
!             if( size(probs) == size(probs) )then
!                 probs = probs
!             else
!                 stop 'size mismatch between probs arrays; oasis_init; simple_pori'
!             endif
!         else
!             if( sum(probs) > 0. )then
!                 ! probs are already read
!             else
!                 probs = 1./real(size(probs))
!             endif
!         endif
!     endif
!     if( report ) write(*,*) 'oasis initialized'
! end subroutine
!
! !>  \brief  importance sampling step 4 the OASIS minimization
! subroutine oasis_mc( self, ospec, refine )
!     use simple_rnd,      only: gasdev, multinomal
!     use simple_opt_subs, only: linmin, check_and_correct_vec
!     class(pori), intent(inout)   :: self
!     class(opt_spec), intent(in)  :: ospec
!     character(len=*), intent(in) :: refine
!     real, allocatable :: x(:)
!     real    :: y, y2, x_old(ngaup)
!     logical :: corrected, found_better
!     integer :: i, j, state, state_best
!     ! extract the previous solution
!     x = extract_prev()
!     ! initialize best-so-far cost
!     yb = ospec%costfun(x, ngaup)
!     if( report ) write(*,*) 'cost of previous solution:', yb
!     ! add the 'old' solution to the list of local minima
!     call locopts%new
!     if( allocated(probs) )then
!         call add_minimum(x, yb, state=nint(optr%get('state')))
!     else
!         call add_minimum(x, yb)
!     endif
!     ! associate costfun in line minimizer spec
!     call spec_linmin%set_costfun(ospec%costfun)
!     ! initialize found-better-indicator
!     found_better =.false.
!     do i=1,spec_linmin%maxits                                ! [1,nsample]
!         x_old = x                                                 ! store old parameter vector
!         do j=1,ngaup                                         ! change the components using Gaussian random sampling
!             x(j) = gasdev(x_old(j), sdevs(j))                ! no limits, since limits give infinite loop at times
!         end do
!         state = 1
!         if( allocated(probs) )then                           ! need to sample multinomal state distribution
!             state = multinomal(probs)                        ! this is the new state
!         endif
!         call ori_glob%set('state', real(state))              ! update ori_glob
!         y = ospec%costfun(x, ngaup)                          ! score the new solution vector
!         if( y <= yb )then                                    ! update the model if a better solution is found
!             found_better = .true.                                 ! indicate better solution found
!             if( report ) write(*,*) 'mc found better:', y
!             yb = y                                           ! update the best cost
!             if( allocated(probs) )then
!                 statecnts(state) = statecnts(state)+1   ! update the state counters
!                 state_best = state                                ! set best state
!             endif
!             if( refine .ne. 'het' .and. refine .ne. 'stoch' )then
!                 spec_linmin%x  = x                           ! set solution in linmin spec
!                 spec_linmin%xi = x-x_old                     ! estimate the search direction
!                 call linmin(spec_linmin,y2)                  ! line minimize
!                 if( y2 < yb )then                            ! better solution found
!                     yb = y2                                  ! update the best cost
!                     x = spec_linmin%x                        ! update solution
!                     if( report )then
!                         write(*,*) 'linmin found better:', y2
!                     endif
!                 endif
!             endif
!             call update_gaup(x)                              ! update Gaussian part of the model
!             call update_probs                                ! update the multinomal part of the model
!             ! correct if cyclic
!             call check_and_correct_vec(spec_linmin, x, corrected)
!             if( corrected )then ! shift the Euler means if x was corrected
!                 do j=1,3
!                     call vars(j)%reset_mean(x(j))
!                 end do
!             endif
!             ! add local minimum to the list of local minima
!             if( allocated(probs) )then
!                 call add_minimum(x, yb, state=state)
!             else
!                 call add_minimum(x, yb)
!             endif
!         else
!             x = x_old ! put back the old point
!         endif
!     end do
!     ! update the parameters in ori object
!     if( found_better )then
!         if( allocated(probs) )then
!             call update_ori_part(x, state_best)
!         else
!             call update_ori_part(x)
!         endif
!     endif
! end subroutine
!
! ! NLP-REFINEMENT PROCEDURE
!
! !>  \brief  for refining a local minimum further using NLP
! subroutine refine_minimum( self, opt, ospec, which )
!     use simple_optimizer, only: optimizer
!     use simple_ori,       only: ori
!     use simple_opt_subs,  only: check_and_correct_vec
!     class(pori), intent(inout)      :: self
!     class(optimizer), intent(inout) :: opt
!     class(opt_spec), intent(inout)  :: ospec
!     integer, intent(in), optional   :: which
!     integer :: pos, j, s
!     real    :: yb, y
!     logical :: corrected
!     if( present(which) )then
!         if( which >= 1 .and. which <= locopts%size() )then
!             pos = which
!         else
!             write(*,*) 'which=', which
!             write(*,*) 'locopts%size=',  locopts%size()
!             stop 'which out of range; refine_minimum; simple_pori'
!         endif
!     else
!         pos = locopts%size()
!     endif
!     if( allocated(probs) )then
!         call get_minimum(ospec%x, yb, state=s, which=pos)
!         call ori_glob%set('state', real(s))
!     else
!         call get_minimum(ospec%x, yb, which=pos)
!     endif
!     ! minimize
!     call opt%minimize(ospec, y)
!     if( y < yb )then ! update the local minimum
!         if( report ) write(*,*) 'NLP found better:', y
!         ! first update the Gaussian model
!         call update_gaup(ospec%x)
!         ! then, correct if cyclic
!         call check_and_correct_vec(ospec, ospec%x, corrected)
!         ! replace the local minimum
!         call locopts%set(pos, rarr=[ospec%x,yb,0.])
!         if( y < yb )then ! if better than previous best
!             if( report ) write(*,*) 'NLP found better best:', y
!             ! update best cost
!             yb = y
!             if( corrected )then ! shift the Euler means if x was corrected
!                 do j=1,3
!                     call vars(j)%reset_mean(ospec%x(j))
!                 end do
!             endif
!             call update_ori_part(ospec%x)
!         endif
!     endif
!     if( report ) write(*,*) 'refined minimum'
! end subroutine

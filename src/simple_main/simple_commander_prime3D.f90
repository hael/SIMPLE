!==Class simple_commander_prime3D
!
! This class contains the set of concrete prime3D commanders of the SIMPLE library. This class provides the glue between the reciver 
! (main reciever is simple_exec program) and the abstract action, which is simply execute (defined by the base class: simple_commander_base). 
! Later we can use the composite pattern to create MacroCommanders (or workflows)
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2016
!
module simple_commander_prime3D
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_commander_base, only: commander_base
use simple_qsys_funs       ! use all in there
use simple_filehandling    ! use all in there
use simple_jiffys          ! use all in there
implicit none

public :: resrange_commander
public :: npeaks_commander
public :: nspace_commander
public :: shellweight3D_commander
public :: prime3D_init_commander
public :: multiptcl_init_commander
public :: het_init_commander
public :: prime3D_commander
public :: cont3D_commander
public :: check3D_conv_commander
private

type, extends(commander_base) :: resrange_commander
  contains
    procedure :: execute      => exec_resrange
end type resrange_commander
type, extends(commander_base) :: npeaks_commander
  contains
    procedure :: execute      => exec_npeaks
end type npeaks_commander
type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander
type, extends(commander_base) :: shellweight3D_commander
 contains
   procedure :: execute      => exec_shellweight3D
end type shellweight3D_commander
type, extends(commander_base) :: prime3D_init_commander 
  contains
    procedure :: execute      => exec_prime3D_init
end type prime3D_init_commander
type, extends(commander_base) :: multiptcl_init_commander
  contains
    procedure :: execute      => exec_multiptcl_init
end type multiptcl_init_commander
type, extends(commander_base) :: het_init_commander
  contains
    procedure :: execute      => exec_het_init
end type het_init_commander
type, extends(commander_base) :: prime3D_commander
  contains
    procedure :: execute      => exec_prime3D
end type prime3D_commander
type, extends(commander_base) :: cont3D_commander
  contains
    procedure :: execute      => exec_cont3D
end type cont3D_commander
type, extends(commander_base) :: check3D_conv_commander
  contains
    procedure :: execute      => exec_check3D_conv
end type check3D_conv_commander

contains

    subroutine exec_resrange( self, cline )
        use simple_hadamard3D_matcher, only: prime3D_find_resrange
        class(resrange_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        if( cline%defined('box') .or. cline%defined('moldiam') )then
            p = params(cline)                     ! parameters generated
            call b%build_general_tbox(p, cline)   ! general objects built
            call b%build_hadamard_prime3D_tbox(p) ! prime objects built
            call prime3D_find_resrange( b, p, p%lp, p%lpstop )
            p%lpstart = p%lp
            call cline%set('lpstart',p%lpstart)   ! for reporting
            call cline%set('lpstop',p%lpstop)     ! for reporting
            write(*,'(A,1X,F5.1)') '>>> LP START:', p%lpstart
            write(*,'(A,2X,F5.1)') '>>> LP STOP:',  p%lpstop
            write(*,'(A,2X,F5.1)') '>>> HP:',       p%hp
        else
            stop 'need either box size or moldiam to estimate resrange; simple_resrange'
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RESRANGE NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_resrange
    
    subroutine exec_npeaks( self, cline )
        use simple_oris, only: oris
        class(npeaks_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(build)  :: b
        type(params) :: p
        p = params(cline) ! parameters generated
        call b%build_general_tbox(p, cline)
        p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        write(*,'(A,1X,I4)') '>>> NPEAKS:', p%npeaks
        ! end gracefully
        call simple_end('**** SIMPLE_NPEAKS NORMAL STOP ****')
    end subroutine exec_npeaks
    
    subroutine exec_nspace(self,cline)
        use simple_math, only: resang
        use simple_oris, only: oris
        class(nspace_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(oris)   :: o
        type(params) :: p
        real         :: ares
        integer      :: i
        p = params(cline) ! parameters generated
        do i=500,5000,500
            o = oris(i)
            call o%spiral
            ares = o%find_angres()
            write(*,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, p%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_shellweight3D(self,cline)
        use simple_cartft_corrcalc, only: cartft_corrcalc
        use simple_cont3D_matcher   ! use all in there
        class(shellweight3D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(params)          :: p
        type(build)           :: b
        type(cartft_corrcalc) :: cftcc
        p = params(cline)                    ! constants & derived constants produced
        ! make sure boxmatch .eq. box
        p%boxmatch = p%box
        call b%build_general_tbox(p, cline)  ! general objects built
        call cont3D_shellweight(b, p, cline)
        call simple_end('**** SIMPLE_SHELLWEIGHT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_shellweight3D
    
    subroutine exec_prime3D_init( self, cline )
        use simple_hadamard3D_matcher, only: gen_random_model, prime3D_find_resrange
        class(prime3D_init_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)       :: p
        type(build)        :: b
        integer, parameter :: MAXIMGS=1000
        p = params(cline) ! parameters generated
        if( p%l_xfel )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('mskfile') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D_init'
            endif
        endif
        if( p%ctf .ne. 'no')then
            if( .not. cline%defined('deftab') )&
            &stop 'need texfile with defocus/astigmatism values for ctf .ne. no mode exec'
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        ! determine resolution range
        if( cline%defined('lp') ) call prime3D_find_resrange( b, p, p%lp, p%lpstop )
        ! determine the number of peaks
        if( .not. cline%defined('npeaks') ) p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
        ! generate the random model
        if( cline%defined('nran') )then
            call gen_random_model( b, p, p%nran )
        else
            if( p%nptcls > MAXIMGS )then
                 call gen_random_model( b, p, MAXIMGS )
            else
                call gen_random_model( b, p )
            endif
        endif
        ! end gracefully
        call qsys_job_finished( p, cline%get_carg('prg') )
        call simple_end('**** SIMPLE_PRIME3D_INIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prime3D_init

    subroutine exec_het_init( self, cline )
        use simple_prime_srch, only: prime_srch
        use simple_rec_master, only: exec_rec_master
        use simple_oris,       only: oris
        use simple_stat               ! use all in there
        use simple_hadamard3D_matcher ! use all in there
        use simple_math               ! use all in there
        use simple_hadamard_common    ! use all in there
        class(het_init_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(params)            :: p
        type(build)             :: b
        type(prime_srch)        :: srch_common !< functionalities common to primesrch2D/3D
        type(oris)              :: best_individual, winner, loser
        type(oris), allocatable :: individuals(:), individuals_opt(:)
        real,       allocatable :: probs(:,:), fitness(:)
        real                    :: best_fitness
        integer                 :: ind, it, i, alloc_stat, noris, loc(1), winner_ind
        character(len=STDLEN), parameter :: fitfun    = 'avg'
        real,                  parameter :: prob_incr = 0.05
        integer,               parameter :: tour_pop  = 3 ! Selection tournament size
        p = params(cline) ! parameters generated
        if( p%nstates < 2 ) stop 'Nonsensical to have nstates < 2; simple_multiptcl_init'
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        call b%build_rec_tbox(p)
        ! Init Indidviduals
        noris = b%a%get_noris()
        allocate( individuals(tour_pop), individuals_opt(tour_pop) ) 
       ! Init state probability vector
        allocate( probs( b%a%get_noris(),p%nstates), stat=alloc_stat)
        allocate( fitness(tour_pop) )
        probs = 1./real(p%nstates)
        if( cline%defined('oritab2') )then
            ! TODO
        endif
        ! init best candidate
        call set_bp_range( b, p, cline )
        best_individual = b%a
        call best_individual%set_all2single( 'corr',-1.)
        !call best_individual%set_all( 'w',    1.)
        ! init for correlations
        srch_common = prime_srch(p, p%nspace*p%nstates, round2even(twopi*real(p%ring2)))
        ! Main loop
        best_fitness = -1.
        do it=p%startit,p%maxits
            write(*,'(A,I6)') '>>> ITERATION ', it
            ! Generation, optimization & fitness
            do i=1,tour_pop
                write(*,'(A,I6)')  '>>> INDIVIDUAL ',i
                call gen_individual( individuals(i) )
                ! ptftcc
                if( it==p%startit .and. i==1 )then
                    call preppftcc4align( b, p, cline )
                    call primesrch3D%kill
                else
                    call prep_refs_pftcc4align( b, p, cline )
                endif
                ! local_search
                individuals_opt(i) = individuals(i)
                call local_search( individuals_opt(i) )
                ! Fitness
                fitness(i) = calc_fitness( individuals_opt(i), fitfun )
            enddo
            ! Tournament selection
            write(*,'(A)')  '>>> PROBABILITY DISTRIBUTION UPDATE'
            loc = maxloc( fitness )
            winner_ind = loc(1)
            winner = individuals( winner_ind )
            print *,'winner fitness', fitness( winner_ind )
            do i=1,tour_pop
                !if( )
                print *,'looser fitness', fitness( i )
                loser = individuals( ind )
                call update_probs( winner, loser )
            enddo
            ! survival of the fitest
            if( fitness(ind) > best_fitness )then
                write(*,'(A)')  '>>> FOUND NEW BEST INDIVIDUAL'            
                best_individual = individuals( ind )
                best_fitness    = fitness( tour_pop )
                print *, calc_fitness( individuals_opt( ind ), fitfun )
                call best_individual%write( 'best_individual_'//int2str_pad(it,3)//'.txt' )
                call individuals_opt(ind)%write( 'best_individual_opt'//int2str_pad(it,3)//'.txt' )
                write(*,'(A,F9.6)')  '>>> CURERENT BEST FITNESS: ', best_fitness           
            endif
            call individuals(ind)%write( 'winner_'//int2str_pad(it,3)//'.txt' )
            call individuals_opt(ind)%write( 'winner_opt'//int2str_pad(it,3)//'.txt' )
        enddo
        ! output
        call best_individual%write( p%outfile )
        call rec_individual( best_individual )
        ! end gracefully
        call simple_end('**** SIMPLE_HET_INIT NORMAL STOP ****', print_simple=.false.)
        contains

            subroutine rec_individual( individual )
                class(oris), intent(inout) :: individual
                type(oris)  :: os
                os  = b%a
                b%a = individual
                call exec_rec_master(b, p, cline )
                b%a = os
            end subroutine rec_individual

            function calc_fitness( individual, method )result( val )
                class(oris),           intent(inout) :: individual
                character(len=STDLEN), intent(in)    :: method
                real, allocatable :: vals(:)
                real    :: val, corrs(p%nstates)
                integer :: iptcl, s, pops( p%nstates )
                val = 0.
                select case( method )
                    case( 'avg' )
                        val = individual%get_avg('corr')
                    case( 'median' )
                        vals = individual%get_all( 'corr' )
                        val  = median( vals )
                        deallocate( vals )
                    case DEFAULT
                        stop 'unknown fitness function'
                end select
            end function

            subroutine gen_individual( individual )
                use simple_rnd,      only: multinomal
                use simple_oris,     only: oris
                class(oris), intent(inout) :: individual
                integer              :: i, iptcl, s
                ! draws individuals from distribution
                individual = best_individual
                write(*,'(A)')  '>>> DRAWING FROM DISTRIBUTION'
                do iptcl = 1,noris
                    s = multinomal( probs( iptcl,:) )
                    call individual%set( iptcl, 'state', real(s) )
                    call progress(iptcl,noris)
                enddo
                ! reconstruction
                call rec_individual( individual )
            end subroutine

            subroutine local_search( individual )
                use simple_ori, only: ori
                use simple_rnd, only: irnd_uni
                class(oris), intent(inout) :: individual
                type(ori) :: o
                real      :: corrs(p%nstates)
                integer   :: iptcl,ref,roind, s, state, proj, loc(1)
                write(*,'(A)')'>>> CORRELATIONS & LOCAL SEARCH'
                do iptcl= 1,noris
                    o     = individual%get_ori( iptcl )
                    state = nint( o%get( 'state' ) )
                    if( state==0 )cycle
                    roind = srch_common%roind( 360.-o%e3get() )
                    proj  = b%e%find_closest_proj( o )
                    corrs = 0.
                    do s = 1,p%nstates   
                        ref = (s-1)*p%nspace + proj                        ! state projdir
                        call preprefs4align(b, p, iptcl, pftcc, ref=ref)
                        corrs(s) = pftcc%corr( ref, iptcl, roind)   ! correlation
                    enddo
                    loc = maxloc(corrs)
                    s   = loc(1)
                    call individual%set(iptcl, 'state', real(s) )                ! updates solution
                    call individual%set(iptcl, 'corr', corrs(s) )                ! updates solution
                    call progress( iptcl, noris )
                enddo
            end subroutine

            subroutine update_probs( won, lost )
                use simple_ori, only: ori
                type(oris), intent(inout) :: won, lost
                integer :: s, iptcl, state, noris
                real :: incr, decr
                incr = prob_incr
                decr = prob_incr / real(p%nstates-1)
                do iptcl = 1,noris
                    state = nint( won%get(iptcl,'state') )
                    if( state /= nint( lost%get(iptcl,'state') ) )then
                        if( state == nint( best_individual%get(iptcl,'state') ) )then
                            ! veer towards best candidate
                            do s=1,p%nstates
                                if( s==state)then
                                    probs( iptcl, s ) = min( 1., probs( iptcl, s )+incr )
                                else
                                    probs( iptcl, s ) = max( 0., probs( iptcl, s )-decr ) 
                                endif 
                            enddo                                   
                        else
                            ! veer away
                            do s=1,p%nstates
                                if( s==state)then
                                    probs( iptcl, s ) = max( 0., probs( iptcl, s )-incr )
                                else
                                    probs( iptcl, s ) = min( 1., probs( iptcl, s )+decr ) 
                                endif 
                            enddo  
                        endif
                        probs( iptcl,: ) = probs( iptcl,: )/sum(probs( iptcl,: )) ! for numerical stability
                    endif
                enddo
            end subroutine

    end subroutine exec_het_init

    ! subroutine exec_het_init( self, cline )
    !     use simple_pftcc_shsrch      ! singleton
    !     use simple_math              ! singleton
    !     use simple_hadamard_common    ! singleton
    !     use simple_prime_srch
    !     use simple_hadamard3D_matcher ! singleton
    !     use simple_rec_master, only: exec_rec_master
    !     use simple_oris,       only: oris
    !     class(het_init_commander), intent(inout) :: self
    !     class(cmdline),            intent(inout) :: cline
    !     type(params)       :: p
    !     type(build)        :: b
    !     type(prime_srch)   :: srch_common             !< functionalities common to primesrch2D/3D
    !     type(oris)         :: prev_a
    !     integer            :: it, iptcl, alloc_stat
    !     logical            :: converged
    !     character(len=STDLEN), parameter :: method = 'shc'
    !     p = params(cline) ! parameters generated
    !     if( p%nstates < 2 ) stop 'Nonsensical to have nstates < 2; simple_multiptcl_init'
    !     call b%build_general_tbox(p, cline)   ! general objects built
    !     call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
    !     if( cline%defined('lp') )then
    !         call b%build_rec_tbox(p)
    !     else
    !         call b%build_eo_rec_tbox(p)
    !     endif
    !     ! init best candidate
    !     call set_bp_range( b, p, cline )
    !     if( cline%defined('oritab2') )then
    !         prev_a = oris( b%a%get_noris() )
    !         call prev_a%read( p%oritab2 )
    !         do iptcl = 1,b%a%get_noris()
    !             call b%a%set( iptcl,'state', prev_a%get( iptcl,'state' ) )
    !         enddo
    !         call prev_a%kill
    !     else
    !         call b%a%rnd_states( p%nstates )
    !     endif
    !     call b%a%write( 'candidate_00.txt')
    !     ! init for correlations
    !     srch_common = prime_srch(p, p%nspace*p%nstates, round2even(twopi*real(p%ring2)))
    !     ! Main loop
    !     do it=p%startit,p%maxits
    !         prev_a = b%a
    !         call update_candidate
    !         if( it==p%startit )then
    !             call preppftcc4align( b, p, cline )
    !             call primesrch3D%kill
    !         else
    !             call prep_refs_pftcc4align( b, p, cline )
    !         endif
    !         call search( method )
    !         ! output
    !         call b%a%write( 'candidate_'//int2str_pad(it,2)//'.txt')
    !         write(*,'(A,I4)')  '>>> ITERATION ',it
    !         converged = b%conv%check_conv3D()
    !         !if( converged )exit
    !     enddo
    !     ! output
    !     call b%a%write( p%outfile )
    !     call update_candidate
    !     ! end gracefully
    !     call simple_end('**** SIMPLE_HET_INIT NORMAL STOP ****')

    !     contains

    !         subroutine update_candidate
    !             use simple_oris, only: oris
    !             use simple_rnd,      only: multinomal
    !             use simple_ran_tabu, only: ran_tabu
    !             type(oris)           :: os, backup
    !             type(ran_tabu)       :: rt
    !             integer, allocatable :: state_inds(:), proj_inds(:), vec(:), proj_discard_inds(:)
    !             logical, allocatable :: vail(:)
    !             real                 :: proj_probs( p%nspace ), probs( p%nstates )
    !             integer              :: i, iptcl, s, proj, proj_pop, minpop, n_discard, ind, proj_ind
    !             integer              :: proj_pops( p%nspace ), pop( p%nstates ), loc(1), state_pop
    !             logical :: discarded
    !             backup = b%a
    !             ! determines population of ptcls to keep per state
    !             pop = 0
    !             do s = 1,p%nstates
    !                 pop( s ) = backup%get_statepop( s )
    !             enddo
    !             if( any(pop<1) )stop 'Empty state in gen_candidate'
    !             loc    = minloc( pop )
    !             minpop = pop( loc(1) )
    !             ! discards ptcls per state
    !             print *,'init pop:', pop
    !             do s=1,p%nstates
    !                 if( s==loc(1) )cycle
    !                 if( pop(s)==0 )cycle
    !                 ! fetch state
    !                 ! calculates projdir distribution
    !                 proj_pops  = 0
    !                 proj_probs = 0.
    !                 state_inds = backup%get_state( s )
    !                 os = oris( pop(s) )
    !                 do i = 1,pop(s)
    !                     call os%set_ori( i, backup%get_ori(state_inds(i)) )
    !                 enddo
    !                 do proj = 1,p%nspace
    !                     proj_pops( proj ) =  os%get_clspop( proj )
    !                 enddo
    !                 proj_probs = real( proj_pops ) / real( pop(s) )
    !                 ! draws number of ptcls per projdir to ignore in reconstruction
    !                 n_discard = pop(s) - minpop
    !                 allocate( proj_discard_inds(n_discard) )
    !                 do i = 1,n_discard
    !                     proj_discard_inds(i) = multinomal( proj_probs )
    !                 enddo
    !                 ! discard ptcls
    !                 do proj = 1,p%nspace
    !                     proj_pop = proj_pops( proj )
    !                     if( proj_pop==0 )cycle
    !                     n_discard = min( proj_pop, count( proj_discard_inds==proj ) )
    !                     if( n_discard==0 )cycle
    !                     proj_inds = os%get_cls_pinds( proj )
    !                     allocate( vec(proj_pop) )
    !                     vec = 1
    !                     vec(1:n_discard) = 0
    !                     rt = ran_tabu( proj_pop )
    !                     call rt%shuffle( vec )
    !                     call rt%kill
    !                     do i = 1,proj_pop
    !                         if( vec(i)==0 )then
    !                             iptcl = state_inds( proj_inds(i) )
    !                             call b%a%set( iptcl, 'state', 0. ) ! only b%a updated
    !                         endif
    !                     enddo                        
    !                     deallocate( vec, proj_inds )
    !                 enddo
    !                 deallocate( proj_discard_inds, state_inds )
    !             enddo
    !             pop = 0
    !             do s = 1,p%nstates
    !                 pop( s ) = b%a%get_statepop( s )
    !             enddo
    !             print *,'final pop', pop
    !             ! reconstruction
    !             call exec_rec_master(b, p, cline )
    !             b%a = backup
    !             call backup%kill
    !         end subroutine

    !         subroutine search( method )
    !             use simple_ori, only: ori
    !             use simple_rnd, only: irnd_uni, ran3
    !             character(len=*) :: method
    !             type(ori) :: o
    !             real      :: prev_corr, corr, frac, corrs(p%nstates), x, y
    !             real      :: target_corr
    !             integer   :: iptcl,ref,roind, s, state, proj, loc(1), pops(p%nstates), cnt, nran
    !             logical   :: avail( p%nstates )
    !             write(*,'(A)')'>>> CORRELATIONS CALCULATION'
    !             !call b%a%set_all( 'corr',-1.)
    !             corrs = 0.
    !             pops  = 0
    !             nran  =0
    !             do iptcl = 1, b%a%get_noris()
    !                 call preprefs4align(b, p, iptcl, pftcc)
    !                 o     = b%a%get_ori( iptcl )
    !                 state = nint( o%get( 'state' ) )
    !                 roind = srch_common%roind( 360.-o%e3get() )
    !                 proj  = b%e%find_closest_proj( o )
    !                 ref   = (state-1)*p%nspace + proj               ! previous reference
    !                 prev_corr   = pftcc%corr( ref, iptcl, roind )   ! previous correlation
    !                 select case( method )
    !                     case( 'shc' )
    !                         target_corr = prev_corr     
    !                     case( 'prev_shc' )
    !                         target_corr = b%a%get( iptcl, 'corr' )
    !                         if( target_corr<=0. )target_corr = prev_corr 
    !                     case DEFAULT
    !                         stop 'Unknown method for target correlation'
    !                 end select
    !                 target_corr = prev_corr     
    !                 ! for stats
    !                 corrs( state ) = corrs( state ) + prev_corr
    !                 pops( state ) = pops( state ) + 1
    !                 ! SHC
    !                 corr  = -1.
    !                 avail = .true.
    !                 avail( state ) = .false.
    !                 cnt = 0
    !                 do while( any(avail) )
    !                     s = irnd_uni( p%nstates )
    !                     if( avail(s) )then
    !                         avail( s ) = .false.
    !                         cnt = cnt + 1
    !                         ref  = (s-1)*p%nspace + proj             ! state projdir
    !                         corr = pftcc%corr( ref, iptcl, roind )   ! correlation
    !                         if( corr >= target_corr )exit            ! first improvement
    !                     endif
    !                 enddo
    !                 if( corr >= target_corr )then                    ! first improvement
    !                     ! new solution
    !                     frac = 100. * real( cnt ) / real( p%nstates )
    !                 else
    !                     ! keep previous
    !                     frac = 100.
    !                     s    = state
    !                     corr = prev_corr
    !                 endif
    !                 if( ran3()*100. > frac )then
    !                     ! randomization
    !                     s    = irnd_uni( p%nstates )
    !                     ref  = (s-1)*p%nspace + proj             ! state projdir
    !                     corr = pftcc%corr( ref, iptcl, roind )   ! correlation
    !                     nran = nran+1
    !                 endif          
    !                 if( s==state )then
    !                     call b%a%set(iptcl, 'mi_state', 1. )                ! updates solution
    !                 else
    !                     call b%a%set(iptcl, 'mi_state', 0. )                ! updates solution
    !                 endif
    !                 call b%a%set(iptcl, 'frac',  frac )                ! updates solution
    !                 call b%a%set(iptcl, 'state', real(s) )             ! updates solution
    !                 call b%a%set(iptcl, 'corr',  corr )            ! updates solution
    !                 call progress( iptcl, b%a%get_noris() )
    !             enddo
    !             print *, 'prev corrs:', corrs/real(pops)
    !             print *, 'nran:', nran
    !         end subroutine
    ! end subroutine exec_het_init

    subroutine exec_multiptcl_init( self, cline )
        use simple_rec_master, only: exec_rec_master
        class(multiptcl_init_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        p = params(cline) ! constants & derived constants produced
        call b%build_general_tbox(p, cline)
        if( cline%defined('state2split') )then
            if( cline%defined('oritab') )then
                p%nstates = b%a%get_nstates()
                call b%a%split_state(p%state2split)
                p%nstates = p%nstates+1
            else
                stop 'Need oritab to be defined when state2split is defined on command line; simple_multiptcl_init'
            endif
        else
            if( p%nstates < 2 ) stop 'Nonsensical to have nstates < 2; simple_multiptcl_init'
            call b%a%rnd_states(p%nstates)
        endif
        if( p%norec .eq. 'no' )then
            if( cline%defined('lp') )then
                call b%build_rec_tbox(p)
            else
                call b%build_eo_rec_tbox(p)
            endif
            call exec_rec_master(b, p, cline, 'startvol')
        endif
        if( p%zero .eq. 'yes' ) call b%a%set_all2single('corr', 0.)
        call b%a%write('multiptcl_startdoc.txt')
        ! end gracefully
        call simple_end('**** SIMPLE_MULTIPTCL_INIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_multiptcl_init
    
    subroutine exec_prime3D( self, cline )
        use simple_hadamard3D_matcher, only: prime3D_exec, prime3D_find_resrange
        use simple_strings,            only: str_has_substr
        class(prime3D_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        integer           :: i, startit, s, n_states, cnt
        logical           :: update_res=.false., converged=.false.
        real              :: lpstart, lpstop
        p = params(cline) ! parameters generated
        if( p%l_xfel )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('automsk') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
            endif
        endif
        if( str_has_substr(p%refine,'neigh') .or. str_has_substr(p%refine,'qcont') &
            &.or. p%refine.eq.'shift' .or. p%refine.eq.'anneal' )then
            if( .not. cline%defined('oritab') )then
                stop 'need oritab input for execution of prime3D with refine mode'
            endif
        endif
        call b%build_general_tbox(p, cline)   ! general objects built
        if( .not. cline%defined('eo') ) p%eo = 'no' ! default
        if( p%eo .eq. 'yes' ) p%dynlp = 'no'    
        if( cline%defined('lp') .or. cline%defined('find')&
        .or. p%eo .eq. 'yes' .or. p%dynlp .eq. 'yes' )then
            ! alles ok!
        else
           stop 'need a starting low-pass limit (set lp or find)!'
        endif
        call b%build_hadamard_prime3D_tbox(p) ! prime3D objects built
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            if( cline%defined('find') )then
                p%lp = b%img%get_lp(p%find)
            endif
            call prime3D_exec(b, p, cline, 0, update_res, converged) ! partition or not, depending on 'part'
        else
            if( p%dynlp .eq. 'yes' )then
                call prime3D_find_resrange( b, p, lpstart, lpstop ) ! determine resolution range
                if( cline%defined('lpstart') )then
                    p%lp = p%lpstart
                else
                    p%lp = lpstart
                endif
                if( cline%defined('lpstop') )then
                    ! defined aleady
                else
                    p%lpstop = lpstop
                endif
            endif
            p%find = int((real(p%box-1)*p%smpd)/p%lp)
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            do i=startit,p%maxits
                call prime3D_exec(b, p, cline, i, update_res, converged)
                if( update_res )then
                    ! dynamic low-pass
                    p%find = p%find+p%fstep
                    p%lp = max(p%lpstop,b%img%get_lp(p%find))
                endif
                if( converged )exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME3D NORMAL STOP ****')
    end subroutine exec_prime3D

    subroutine exec_cont3D( self, cline )
        use simple_cont3D_matcher, only: cont3D_exec
        class(cont3D_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(params) :: p
        type(build)  :: b
        integer       :: i, startit
        logical       :: converged=.false.
        p = params(cline) ! parameters generated
        if( p%xfel .eq. 'yes' )then
            if( cline%defined('msk') .or. cline%defined('mw') .or.&
            cline%defined('nvox') .or. cline%defined('automsk') )then
                stop 'no mask allowed when processing XFEL patterns; simple_prime3D'
            endif
        endif
        call b%build_general_tbox(p, cline)
        call b%build_cont3D_tbox(p)
        if( cline%defined('part') )then
            if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
            call cont3D_exec(b, p, cline, 0, converged) ! partition or not, depending on 'part'
        else
            startit = 1
            if( cline%defined('startit') ) startit = p%startit
            do i=startit,p%maxits
                call cont3D_exec(b, p, cline, i, converged)
                if(converged) exit
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONT3D NORMAL STOP ****')
    end subroutine exec_cont3D
    
    subroutine exec_check3D_conv( self, cline )
        use simple_math,    only: rad2deg, get_lplim
        class(check3D_conv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(params)      :: p
        type(build)       :: b
        real, allocatable :: maplp(:)
        integer           :: istate, loc(1)
        logical           :: here, limset, converged, update_res
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        ! nstates consistency check
        ! Now incompatible with empty states
        !if( cline%defined('nstates') )then
        !    if( p%nstates /= b%a%get_nstates() ) stop 'Inconsistent number of states between command-line and oritab'
        !endif
        limset = .false.
        if( p%eo .eq. 'yes' )then
            allocate( maplp(p%nstates) )
            maplp = 0.
            do istate=1,p%nstates
                if( b%a%get_statepop( istate ) == 0 )cycle ! empty state
                p%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
                inquire(file=p%fsc, exist=here)
                if( here )then
                    b%fsc(istate,:) = file2rarr(p%fsc)
                    maplp(istate)   = max(b%img%get_lp(get_lplim(b%fsc(istate,:))),2.*p%smpd)
                else
                    write(*,*) 'Tried to open the fsc file: ', trim(p%fsc)
                    stop 'but it does not exist!'
                endif
            enddo
            loc     = maxloc( maplp )
            p%state = loc(1)            ! state with worst low-pass
            p%lp    = maplp( p%state )  ! worst lp
            p%fsc   =  'fsc_state'//int2str_pad(p%state,2)//'.bin'
            deallocate(maplp)
            limset = .true.
        endif
        ! Let find override the command line input lp (if given)
        if( .not. limset .and. cline%defined('find') )then
            p%lp = b%img%get_lp(p%find)
            limset = .true.
        endif
        ! Method for setting lp with lowest priority is lp on the command line
        if( cline%defined('lp') ) limset = .true.
        ! If we arrived here and the limit wasn't set: fall over
        if( limset )then
            ! we are happy
        else
            ! we fall over
            stop 'No method available to set low-pass limit! ABORTING...'
        endif
        ! calculate angular threshold
        p%athres = rad2deg(atan(p%lp/(p%moldiam/2.)))
        ! check convergence
        if( cline%defined('update_res') )then
            update_res = .false.
            if( cline%get_carg('update_res').eq.'yes' )update_res = .true.
            converged = b%conv%check_conv3D( update_res )
        else
            converged = b%conv%check_conv3D()
        endif
        ! reports convergence, shift activation, resolution update and
        ! fraction of search space scanned to the distr commander
        if( p%doshift )then
            call cline%set('trs', p%trs)        ! activates shift serach
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        if( update_res )then
            call cline%set('update_res', 'yes') ! fourier index to be updated in distr commander
        else
            call cline%set('update_res', 'no')
        endif
        call cline%set('frac', b%conv%get('frac'))
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK3D_CONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check3D_conv

end module simple_commander_prime3D

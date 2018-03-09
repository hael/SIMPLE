module simple_prime3D_srch_tester
include 'simple_lib.f08'
use simple_strategy2D3D_common    ! singleton
use simple_strategy3D_matcher ! singleton
use simple_cmdline            ! singleton
use simple_strategy3D_srch,       only: strategy3D_srch
use simple_oris,               only: oris
use simple_ori,                only: ori
use simple_build,              only: build
use simple_params,             only: params
use simple_image,              only: image
implicit none

public :: exec_prime3D_srch_test

private
#include "simple_local_flags.inc"

! module global constants
! integer,           parameter :: NPROJS      = 15
! integer,           parameter :: NSTATES     = 3
! character(len=32), parameter :: ptclsname   = 'prime3D_srch_test_ptcls.mrc'
! character(len=32), parameter :: orisname    = 'prime3D_srch_test_oris.txt'
! character(len=32), parameter :: outfilename = 'prime3D_srch_test_algndoc.txt'
! real,              parameter :: LPLIM       = 10.

! ! module global objects
! type(build)               :: b
! type(params)              :: p
! type(prime3D_srch)        :: primesrch3D
! type(oris)                :: o_ptcls
! type(image), allocatable  :: imgs_ptcls(:)


contains

    subroutine exec_prime3D_srch_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
    !     type(cmdline) :: cline_local
    !     ! refine = no; nstates=1;
    !     write(*,*)'>>> REFINE=NO; NSTATES=1'
    !     cline_local = cline
    !     call cline_local%set('refine', 'no')
    !     call cline_local%set('nstates',1.)
    !     call setup_testenv( cline_local, be_verbose )
    !     call test_prep4srch
    !     call test_prepcorr4srch
    !     call test_prep_reforis( cline_local )
    !     call shutdown_testenv
    !     ! refine = shc; nstates=1;
    !     write(*,*)'>>> REFINE=SHC; NSTATES=1'
    !     cline_local = cline
    !     call cline_local%set('refine', 'shc')
    !     call cline_local%set('nstates',1.)
    !     call setup_testenv( cline_local, be_verbose )
    !     call test_prep4srch
    !     call test_prepcorr4srch
    !     call test_prep_reforis( cline_local )
    !     call shutdown_testenv
    !     ! ! refine = no; nstates=4;
    !     write(*,*)'>>> REFINE=NO; NSTATES=',NSTATES
    !     cline_local = cline
    !     call cline_local%set('refine', 'no')
    !     call cline_local%set('nstates',real(NSTATES))
    !     call setup_testenv( cline_local, be_verbose )
    !     call test_prep4srch
    !     call test_prepcorr4srch
    !     call test_prep_reforis( cline_local )
    !     call shutdown_testenv
    !     ! ! refine = neigh; nstates=1;
    !     ! write(*,*)'>>> REFINE=NEIGH; NSTATES=1'
    !     ! cline_local = cline
    !     ! call cline_local%set('refine', 'neigh')
    !     ! call cline_local%set('nstates',1.)
    !     ! call setup_testenv( cline_local, be_verbose )
    !     ! call test_prep_reforis( cline_local )
    !     ! call shutdown_testenv
    !     ! ! refine = shcneigh; nstates=1;
    !     ! write(*,*)'>>> REFINE=SHCNEIGH; NSTATES=1'
    !     ! cline_local = cline
    !     ! call cline_local%set('refine', 'shcneigh')
    !     ! call cline_local%set('nstates',1.)
    !     ! call setup_testenv( cline_local, be_verbose )
    !     ! call test_prep_reforis( cline_local )
    !     ! call shutdown_testenv
    end subroutine exec_prime3D_srch_test

    ! subroutine setup_testenv( cline, be_verbose )
    !     use simple_strings,        only: str_has_substr
    !     use simple_projector_hlev, only: projvol
    !     class(cmdline),    intent(inout) :: cline
    !     logical, optional, intent(in)    :: be_verbose
    !     !type(ori)          :: o
    !     integer            :: i   !, state ,noris

    !     if( present(be_verbose) ) verbose = be_verbose
    !     verbose=.true.
    !     ! it is assumed that vol1, smpd, msk are part of the inputted command line
    !     ! setting the remainder of the command line up in here
    !     if( str_has_substr( cline%get_carg('refine'),'neigh'))then
    !         call cline%set('nnn', 5.)
    !     endif
    !     call cline%set('nspace',  real(NPROJS))
    !     call cline%set('nptcls',  real(NPROJS))
    !     call cline%set('lp',      LPLIM       )
    !     call cline%set('ctf',     'no'        )
    !     call cline%set('outfile', outfilename )
    !     if( str_has_substr(cline%get_carg('refine'),'shc') )then
    !         call cline%set('npeaks', 1.)
    !     else
    !         call cline%set('npeaks', real(MAXNPEAKS))
    !     endif
    !     if( nint(cline%get_rarg('nstates'))==NSTATES )then
    !         call cline%set('vol2',cline%get_carg('vol1'))
    !         call cline%set('vol3',cline%get_carg('vol1'))
    !         call cline%set('vol4',cline%get_carg('vol1'))
    !     endif
    !     VerbosePrint('cline keys set')
    !     ! generate orientations
    !     call o_ptcls%new(NPROJS)
    !     call o_ptcls%spiral
    !     !call o_ptcls%rnd_inpls
    !     call o_ptcls%set_all2single('lp',LPLIM)
    !     call o_ptcls%write( orisname )
    !     VerbosePrint 'orientations generated'
    !     ! create parameters and build
    !     p = params(cline)                     ! parameters generated
    !     ! instantiates builder
    !     call b%build_general_tbox(p, cline)   ! general objects built
    !     ! set resolution range
    !     p%kfromto(1) = 2
    !     p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
    !     ! simulate images
    !     call b%vol%read(p%vols(1))
    !     call b%vol%fwd_ft
    !     call b%vol%expand_cmat(p%alpha)
    !     imgs_ptcls = projvol(b%vol, o_ptcls, p)
    !     call b%vol%kill_expanded
    !     do i=1,NPROJS
    !         call imgs_ptcls(i)%write(ptclsname, i)
    !     end do
    !     call cline%set('stk',  ptclsname)
    !     call cline%set('oritab', orisname )
    !     p = params(cline)                     ! parameters generated
    !     call b%build_general_tbox(p, cline)   ! general objects built
    !     call b%build_hadamard_prime3D_tbox(p) ! prime objects built
    !     ! CALCULATE ANGULAR THRESHOLD
    !     p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
    !     ! DETERMINE THE NUMBER OF PEAKS
    !     if( .not. cline%defined('npeaks') )then
    !         select case(p%refine)
    !             case('no', 'neigh', 'adasym')
    !                 p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
    !                 if( str_has_substr(p%refine,'adasym') ) p%npeaks = p%npeaks * p%nsym
    !             case DEFAULT
    !                 p%npeaks = 1
    !         end select
    !     endif
    !     ! GENERATE PROJECTIONS (POLAR FTs)
    !     call b%a%calc_hard_weights(p%frac)
    !     b%a = o_ptcls
    !     call pftcc%new(NPROJS, p)
    !     call b%img_match%init_polarizer( pftcc, p%alpha )
    !     do i=1,NPROJS
    !         call imgs_ptcls(i)%clip(b%img_match)
    !         call b%img_match%mask(p%msk,'soft')
    !         call b%img_match%fwd_ft
    !         call b%img_match%polarize(pftcc, i, isptcl=.true.)
    !         call b%img_match%polarize(pftcc, i, isptcl=.false.)
    !     enddo
    !     ! call preppftcc4align( b, p, cline )
    !     ! The pftcc & primesrch3D objects are now globally available in the module
    !     ! because of the use simple_hadamard3D_matcher statement in the top
    !     ! now instantiatable, so create it
    !     ! call primesrch3D%new(pftcc, b%a, b%e, p)
    !     VerbosePrint 'end setup_testenv'
    ! end subroutine setup_testenv

    ! subroutine test_prep4srch
    !     type(ori) :: o, o_saved
    !     real      :: e3, shvec(2), x, y
    !     integer   :: nrefs, nrots, ind,i, state, proj
    !     nrefs = primesrch3D%get_ntotrefs()
    !     nrots = primesrch3D%get_nrots()
    !     ! refine=no,shc; states=1
    !     ! do i=1,NPROJS
    !     !     call primesrch3D%new(i, pftcc, b%a, b%e, p)
    !     !     o = b%a%get_ori(i)
    !     !     o_saved = o
    !     !     call o%rnd_inpl( p%trs )
    !     !     call o%set('proj',real(i))
    !     !     state = nint(o%get('state'))
    !     !     e3    = o%e3get()
    !     !     ind   = pftcc%get_roind( 360.-e3 )
    !     !     x     = o%get('x')
    !     !     y     = o%get('y')
    !     !     proj  = b%e%find_closest_proj( o, 1 )
    !     !     call b%a%set_ori(i, o)
    !     !     call primesrch3D%prep4srch(p%lp )
    !     !     if(state.ne.primesrch3D%get_prevstate())stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 1'
    !     !     shvec = primesrch3D%get_prevshvec()
    !     !     if( x.ne.shvec(1) )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 2'
    !     !     if( y.ne.shvec(2) )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 3'
    !     !     if( ind.ne.primesrch3D%get_prevroind() )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 4'
    !     !     if( proj.ne.primesrch3d%get_prevref() )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 5'
    !     !     call b%a%set_ori(i,o_saved)
    !     ! enddo
    !     ! other cases
    !     VerbosePrint 'end setup_prep4srch'
    ! end subroutine test_prep4srch

    ! subroutine test_prepcorr4srch
    !     use simple_rnd,        only: ran3
    !     type(ori) :: o
    !     real      :: prev_corr, corr
    !     integer   :: nrefs, nrots,i,j !, ind,state, proj
    !     nrefs = primesrch3D%get_ntotrefs()
    !     nrots = primesrch3D%get_nrots()
    !     ! if( p%refine.eq.'no' .and. p%ctf.eq.'no' )then
    !     !     ! refine=no; states=1
    !     !     do j=1,10
    !     !         do i=1,p%nptcls
    !     !             call primesrch3D%new(i, pftcc, b%a, b%e, p)
    !     !             o = b%a%get_ori(i)
    !     !             prev_corr = ran3()
    !     !             call b%a%set(i,'corr',prev_corr)
    !     !             call primesrch3D%prep4srch(p%lp )
    !     !             corr = primesrch3D%get_prevcorr()
    !     !             if( p%nstates==1 )then
    !     !                 if( abs(2.* corr - prev_corr) < 0.0001 )then
    !     !                     print *, 'corr = ', corr
    !     !                     stop 'Failed in test_prepcorr4srch::simple_prime3D_srch_tester 1'
    !     !                 endif
    !     !             else

    !     !                 if( corr < 0.99 )then
    !     !                     print *, 'corr = ', corr
    !     !                     stop 'Failed in test_prepcorr4srch::simple_prime3D_srch_tester 2'
    !     !                 endif
    !     !             endif
    !     !         enddo
    !     !     enddo
    !     ! else
    !     !     do i=1,p%nptcls
    !     !         call primesrch3D%new(i, pftcc, b%a, b%e, p)
    !     !         call primesrch3D%prep4srch(p%lp )
    !     !         corr = primesrch3D%get_prevcorr()
    !     !         if( corr < 0.99 )then
    !     !             print *, 'corr = ', corr
    !     !             stop 'Failed in test_prepcorr4srch::simple_prime3D_srch_tester 3'
    !     !         endif
    !     !     enddo
    !     ! endif
    !     VerbosePrint 'end setup_prepcorr4srch'
    ! end subroutine test_prepcorr4srch

    ! subroutine test_prep_reforis( cline)
    !     class(cmdline),    intent(inout) :: cline
    !     type(oris) :: test_os
    !     type(ori)  :: o, oref, orefs
    !     integer    :: iptcl,i, ref,s, ind
    !     !real,allocatable :: vals(:), test_vals(:)
    !     !integer,allocatable :: ivals(:), test_ivals(:)
    !     integer,allocatable :: srch_order(:),nnmat(:,:)
    !     select case(p%refine)
    !         case('no','shc')
    !             ! if( p%nstates==1 )then
    !             !
    !             ! else
    !             !     test_os = oris( p%nspace*p%nstates )
    !             !     do iptcl=1,p%nptcls
    !             !         call primesrch3D%new(iptcl, pftcc, b%a, b%e, p)
    !             !         call primesrch3D%prep4srch(p%lp)
    !             !         test_os = primesrch3D%get_o_refs( p%nspace*p%nstates )
    !             !         do ref=1,p%nspace
    !             !             oref = test_os%get_ori(ref)
    !             !             do s=1,NSTATES-1
    !             !                 orefs = test_os%get_ori(s*NPROJS+ref)
    !             !                 if( (orefs.euldist.oref)>.001)stop 'Failed test_prep_reforis 21'
    !             !             enddo
    !             !         enddo
    !             !         call test_os%nearest_neighbors( NSTATES, nnmat )
    !             !         do ref=1,p%nspace
    !             !             oref = test_os%get_ori(ref)
    !             !             do s=1,NSTATES
    !             !                 orefs = test_os%get_ori(nnmat(ref,s))
    !             !                 if( (orefs.euldist.oref)>.001 .or. nint(oref%get('class')).ne.nint(orefs%get('class')) &
    !             !                     & .or. nint(orefs%get('state'))>NSTATES )stop 'Failed test_prep_reforis 22'
    !             !             enddo
    !             !         enddo
    !             !         deallocate(nnmat)
    !             !         call test_os%kill
    !             !     enddo
    !             ! endif
    !             ! srch_order = primesrch3D%get_srch_order()
    !             ! if( minval(srch_order)<1 )stop 'Failed test_prep_reforis 23'
    !             ! if( maxval(srch_order)>primesrch3D%get_ntotrefs() )stop 'Failed test_prep_reforis 24'
    !             ! do i=1,primesrch3D%get_ntotrefs()
    !             !     if( count(srch_order==i) /= 1 )stop 'Failed test_prep_reforis 25'
    !             ! enddo
    !             ! deallocate( srch_order )
    !         case('neigh','shcneigh')
    !             ! if( p%nstates==1 )then
    !             !     do iptcl=1,p%nptcls
    !             !         call primesrch3D%new(iptcl, pftcc, b%a, b%e, p)
    !             !         o = b%a%get_ori( iptcl )
    !             !         call primesrch3D%prep4srch(p%lp, b%nnmat)
    !             !         srch_order = primesrch3D%get_srch_order()
    !             !         if( minval(srch_order)<1 )stop 'Failed test_prep_reforis 32'
    !             !         if( maxval(srch_order)>primesrch3D%get_ntotrefs() )stop 'Failed test_prep_reforis 33'
    !             !         do ind=1,p%nnn
    !             !             ref = srch_order(ind)
    !             !             if( count(srch_order==ref) /= 1 )stop 'Failed test_prep_reforis 34'
    !             !         enddo
    !             !         deallocate( srch_order )
    !             !     enddo
    !             ! endif
    !         case DEFAULT
    !             stop 'not implemented yet'
    !     end select
    !     VerbosePrint 'end prep_reforis'
    ! end subroutine test_prep_reforis

    ! subroutine shutdown_testenv
    !     integer :: i
    !     call b%kill_general_tbox
    !     call b%kill_hadamard_prime3D_tbox
    !     do i=1,NPROJS
    !         call imgs_ptcls(i)%kill
    !     end do
    !     deallocate(imgs_ptcls)
    ! end subroutine shutdown_testenv

end module simple_prime3D_srch_tester

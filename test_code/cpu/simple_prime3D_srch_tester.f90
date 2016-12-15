module simple_prime3D_srch_tester
use simple_hadamard_common    ! singleton
use simple_hadamard3D_matcher ! singleton
use simple_defs               ! singleton
use simple_cuda_defs          ! singleton
use simple_jiffys,            ! singleton
use simple_cmdline            ! singleton
use simple_math               ! singleton
use simple_timing             ! singleton
use simple_oris,               only: oris
use simple_ori,                only: ori
use simple_build,              only: build
use simple_params,             only: params
use simple_image,              only: image
implicit none

public :: exec_prime3D_srch_test

private

! module global constants
integer,           parameter :: NPROJS      = 15
integer,           parameter :: NSTATES     = 4
integer,           parameter :: MAXNPEAKS   = 10
character(len=32), parameter :: ptclsname   = 'prime3D_srch_test_ptcls.mrc'
character(len=32), parameter :: orisname    = 'prime3D_srch_test_oris.txt'
character(len=32), parameter :: outfilename = 'prime3D_srch_test_algndoc.txt'
real,              parameter :: LPLIM       = 10.

! module global objects
type(build)               :: b
type(params)              :: p
type(oris)                :: o_ptcls
type(image), allocatable  :: imgs_ptcls(:)
logical                   :: verbose=.false.

contains

    subroutine exec_prime3D_srch_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        type(cmdline) :: cline_local
        ! refine = no; nstates=1; 
        write(*,*)'>>> REFINE=NO; NSTATES=1'
        cline_local = cline
        call cline_local%set('refine', 'no')
        call cline_local%set('nstates',1.)
        call setup_testenv( cline_local, be_verbose )
        call test_ctfparms
        call test_stochastic_weights
        call test_sort_shifted_npeaks
        call test_prep4srch
        call test_prepcorr4srch
        call test_prep_reforis
        call shutdown_testenv
        ! refine = shc; nstates=1;
        write(*,*)'>>> REFINE=SHC; NSTATES=1'
        cline_local = cline
        call cline_local%set('refine', 'shc')
        call cline_local%set('nstates',1.)
        call setup_testenv( cline_local, be_verbose )
        call test_ctfparms
        call test_prep4srch
        call test_prepcorr4srch
        call test_prep_reforis
        call shutdown_testenv        
        ! refine = no; nstates=4; 
        write(*,*)'>>> REFINE=NO; NSTATES=',NSTATES
        cline_local = cline
        call cline_local%set('refine', 'no')
        call cline_local%set('nstates',real(NSTATES))
        call setup_testenv( cline_local, be_verbose )
        call test_ctfparms
        call test_stochastic_weights
        call test_sort_shifted_npeaks
        call test_prep4srch
        call test_prepcorr4srch
        call test_prep_reforis
        call shutdown_testenv
        ! refine = neigh; nstates=1;
        write(*,*)'>>> REFINE=NEIGH; NSTATES=1'
        cline_local = cline
        call cline_local%set('refine', 'neigh')
        call cline_local%set('nstates',1.)
        call setup_testenv( cline_local, be_verbose )
        call test_prep_reforis
        call shutdown_testenv
        ! refine = shcneigh; nstates=1;
        write(*,*)'>>> REFINE=SHCNEIGH; NSTATES=1'
        cline_local = cline
        call cline_local%set('refine', 'shcneigh')
        call cline_local%set('nstates',1.)
        call setup_testenv( cline_local, be_verbose )
        call test_prep_reforis
        call shutdown_testenv
        ! refine = qcontneigh; nstates=1;
        write(*,*)'>>> REFINE=QCONTNEIGH; NSTATES=1'
        cline_local = cline
        call cline_local%set('refine', 'qcontneigh')
        call cline_local%set('nstates',1.)
        call setup_testenv( cline_local, be_verbose )
        call test_prep_reforis
        call shutdown_testenv
        ! refine = qcontneigh; nstates=4;
        write(*,*)'>>> REFINE=QCONTNEIGH; NSTATES=4'
        cline_local = cline
        call cline_local%set('refine', 'qcontneigh')
        call cline_local%set('nstates',real(NSTATES))
        call setup_testenv( cline_local, be_verbose )
        call test_prep_reforis
        call shutdown_testenv
    end subroutine exec_prime3D_srch_test

    subroutine setup_testenv( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        type(ori)          :: o
        integer            :: i, state,noris
        if( present(be_verbose) ) verbose = be_verbose
        verbose=.true.
        ! it is assumed that vol1, smpd, msk are part of the inputted command line
        ! setting the remainder of the command line up in here
        if( str_has_substr( cline%get_carg('refine'),'neigh'))then
            call cline%set('nnn', 5.)
        endif
        call cline%set('nspace',  real(NPROJS))
        call cline%set('nptcls',  real(NPROJS))
        call cline%set('lp',      LPLIM       )
        call cline%set('ctf',     'no'        )
        call cline%set('outfile', outfilename )
        if( str_has_substr(cline%get_carg('refine'),'shc') )then
            call cline%set('npeaks', 1.)
        else
            call cline%set('npeaks', real(MAXNPEAKS))
        endif
        if( nint(cline%get_rarg('nstates'))==NSTATES )then
            call cline%set('vol2',cline%get_carg('vol1'))
            call cline%set('vol3',cline%get_carg('vol1'))
            call cline%set('vol4',cline%get_carg('vol1'))
        endif
        if( verbose )print *,'cline keys set'
        ! generate orientations
        call o_ptcls%new(NPROJS)
        call o_ptcls%spiral
        call o_ptcls%rnd_inpls
        call o_ptcls%set_all2single('lp',LPLIM)
        call o_ptcls%write( orisname )
        if( verbose )print *,'orientations generated'
        ! create parameters and build
        p = params(cline)                     ! parameters generated
        p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
        ! instantiates builder
        call b%build_general_tbox(p, cline)   ! general objects built
        ! set resolution range
        p%kfromto(1) = 2
        p%kfromto(2) = b%img%get_find(p%lp)
        ! simulate images
        call b%vol%read(p%vols(1))
        imgs_ptcls = b%proj%projvol(b%vol, o_ptcls, p)
        do i=1,NPROJS
            call imgs_ptcls(i)%write(ptclsname, i)
        end do
        call cline%set('stk',  ptclsname)
        call cline%set('oritab', orisname )
        p = params(cline)                     ! parameters generated
        p%boxmatch = p%box                    !!!!!!!!!!!!!!!!!! 4 NOW
        call b%build_general_tbox(p, cline)   ! general objects built
        call b%build_hadamard_prime3D_tbox(p) ! prime objects built
        ! CALCULATE ANGULAR THRESHOLD
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            select case(p%refine)
                case('no', 'neigh','qcont', 'qcontneigh')
                    p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
                case DEFAULT
                    p%npeaks = 1
            end select
        endif
        ! GENERATE PROJECTIONS (POLAR FTs)
        call b%a%calc_hard_ptcl_weights(p%frac)
        call preppftcc4align( b, p, cline )
        ! The pftcc & primesrch3D objects are now globally available in the module
        ! because of the use simple_hadamard3D_matcher statement in the top
        if( verbose )print *,'end setup_testenv'        
    end subroutine setup_testenv

    subroutine test_prep4srch
        use simple_prime_srch, only: prime_srch
        type(prime_srch) :: srch_common
        type(ori) :: o
        real      :: class,e3, shvec(2), x, y
        integer   :: nrefs, nrots, ind,i, state, proj
        nrefs = primesrch3D%get_ntotrefs()
        nrots = primesrch3D%get_nrots()
        ! refine=no,shc; states=1
        srch_common = prime_srch( p, nrefs, nrots)
        do i=1,NPROJS
            o = b%a%get_ori(i)
            call o%rnd_inpl( p%trs )
            call o%set('class',real(i))
            state = nint(o%get('state'))
            e3    = o%e3get()
            ind   = srch_common%roind( 360.-e3 )
            x     = o%get('x')
            y     = o%get('y')
            proj  = b%e%find_closest_proj( o, 1 )
            call primesrch3D%prep4srch( o )
            if(state.ne.primesrch3D%get_prevstate())stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 1'
            shvec = primesrch3D%get_prevshvec()
            if( x.ne.shvec(1) )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 2'
            if( y.ne.shvec(2) )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 3'
            if( ind.ne.primesrch3D%get_prevroind() )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 4'
            if( proj.ne.primesrch3d%get_prevref() )stop 'Failed simple_prime3D_srch_tester:: test_prep4srch 5'
        enddo

        ! other cases
        if( verbose )print *,'end setup_prep4srch'
    end subroutine test_prep4srch

    subroutine test_prepcorr4srch
        use simple_prime_srch, only: prime_srch
        use simple_rnd,        only: ran3
        type(prime_srch) :: srch_common
        type(ori) :: o
        real      :: prev_corr
        integer   :: nrefs, nrots, ind,i, state, proj, j
        nrefs = primesrch3D%get_ntotrefs()
        nrots = primesrch3D%get_nrots()
        srch_common = prime_srch( p, nrefs, nrots)
        if( p%refine.eq.'no' .and. p%ctf.eq.'no' )then
            ! refine=no; states=1
            do j=1,10
                do i=1,p%nptcls
                    o = b%a%get_ori(i)
                    prev_corr = ran3()
                    call o%set('corr',prev_corr)
                    call primesrch3D%prep4srch( o )
                    call primesrch3D%prep_corr4srch(  pftcc, i, p%lp, o )
                    if( p%nstates==1 )then
                        if( abs(2.*primesrch3D%get_prevcorr()-prev_corr-1.)>0.005 ) &
                            & stop 'Failed in test_prepcorr4srch::simple_prime3D_srch_tester 1'
                    else
                        if( primesrch3D%get_prevcorr()<0.995 ) &
                            & stop 'Failed in test_prepcorr4srch::simple_prime3D_srch_tester 2'
                    endif
                enddo
            enddo
        else
            do i=1,p%nptcls
                o = b%a%get_ori(i)
                prev_corr = ran3()
                call o%set('corr',prev_corr)
                call primesrch3D%prep4srch( o )
                call primesrch3D%prep_corr4srch(  pftcc, i, p%lp, o )
                if( primesrch3D%get_prevcorr()<0.995 ) &
                    & stop 'Failed in test_prepcorr4srch::simple_prime3D_srch_tester 3'
            enddo
        endif
        call srch_common%kill
        if( verbose )print *,'end setup_prepcorr4srch'
    end subroutine test_prepcorr4srch

    subroutine test_stochastic_weights
        type(oris) :: test_os
        real       :: wcorr, corrs(p%npeaks)
        integer :: i
        wcorr   = 0.
        test_os = oris( p%npeaks )
        call test_os%set_all2single('corr',0.66666)
        call test_os%set(p%npeaks,'corr',-.00054)
        call primesrch3D%stochastic_weights( wcorr, test_os )
        corrs = test_os%get_all('corr')
        do i=1,p%npeaks
            print *,i,test_os%get(i,'corr'),test_os%get(i,'ow')
        enddo
        print *,wcorr
        call test_os%kill
    end subroutine test_stochastic_weights

    subroutine test_sort_shifted_npeaks
        type(oris) :: test_os
        type(ori)  :: o
        integer :: i,j,IND
        do j=1,5
            IND = MAXNPEAKS-j
            test_os = oris( p%npeaks )
            do i=1,p%npeaks
                call o%rnd_ori
                call o%set('class',real(i))
                call o%set('corr',real(i)/(4.*real(p%npeaks)))
                call test_os%set_ori(i,o)
            enddo
            call test_os%set(IND,'corr',.5+real(j)/10.)
            call primesrch3D%sort_shifted_npeaks( test_os )
            if( nint(test_os%get(p%npeaks,'class')).ne.IND)print *,'Failed in simple_prime3D_srch_tester::test_sort_shifted_npeaks'
            call test_os%kill
        enddo
    end subroutine test_sort_shifted_npeaks

    subroutine test_ctfparms
        use simple_rnd, only: ran3
        real      :: ctf_vals( p%npeaks,3 )
        type(ori) :: o
        real      :: dfx,dfy,angast
        integer   :: i, j
        p%ctf         = 'flip'
        p%tfplan%flag = p%ctf
        p%tfplan%mode = 'astig'
        do j=1,3
            call primesrch3D%new( b%a, b%e, p )
            do i=1,p%npeaks
                call o%rnd_ori
                dfx    = (ran3()-.5)*2.
                dfy    = (ran3()-.5)*2.
                angast = (ran3()-.5)*179.99
                call o%set('dfx'   , dfx)
                call o%set('dfy'   , dfy)
                call o%set('angast', angast)
                call primesrch3D%prep_ctfparms( o )
                if( primesrch3D%get_angast() .ne. angast)stop 'Failed in simple_prime3D_srch_tester::test_ctfparms 1'
                if( primesrch3D%get_dfx() .ne. dfx)stop 'Failed in simple_prime3D_srch_tester::test_ctfparms 2'
                if( primesrch3D%get_dfy() .ne. dfy)stop 'Failed in simple_prime3D_srch_tester::test_ctfparms 3'
            enddo
        enddo
        p%ctf         = 'no'
        p%tfplan%flag = p%ctf
        p%tfplan%mode = 'no'
        call primesrch3D%set_ctf('no')
    end subroutine test_ctfparms

    subroutine test_prep_reforis
        type(oris) :: test_os
        type(ori)  :: o, oref, orefs
        integer    :: iptcl,i, ref,s, ind
        real,allocatable :: vals(:), test_vals(:)
        integer,allocatable :: ivals(:), test_ivals(:), srch_order(:),nnmat(:,:)
        real :: euldist
        select case(p%refine)
            case('no','shc')
                if( p%nstates==1 )then
                    test_os = oris( p%nspace )
                    do iptcl=1,p%nptcls
                        o = b%a%get_ori( iptcl )
                        call preprefs4align(b, p, iptcl, pftcc) 
                        call primesrch3D%prep4srch( o )
                        test_os = primesrch3D%get_o_refs( p%nspace )
                        ! e1
                        vals = b%e%get_all('e1')
                        test_vals = test_os%get_all('e1')
                        if( any(abs(vals-test_vals)>TINY) )print *,'failed test_prep_reforis 1'
                        ! e2
                        vals = b%e%get_all('e2')
                        test_vals = test_os%get_all('e2')
                        if( any(abs(vals-test_vals)>TINY) )print *,'failed test_prep_reforis 2'
                        ! e3
                        vals = b%e%get_all('e3')
                        test_vals = test_os%get_all('e3')
                        if( any(abs(vals-test_vals)>TINY) )print *,'failed test_prep_reforis 3'
                    enddo
                else
                    test_os = oris( p%nspace*p%nstates )
                    do iptcl=1,p%nptcls
                        o = b%a%get_ori( iptcl )
                        call preprefs4align(b, p, iptcl, pftcc) 
                        call primesrch3D%prep4srch( o )
                        test_os = primesrch3D%get_o_refs( p%nspace*p%nstates )
                        do ref=1,p%nspace
                            oref = test_os%get_ori(ref)
                            do s=1,NSTATES-1
                                orefs = test_os%get_ori(s*NPROJS+ref)
                                if( (orefs.euldist.oref)>.001)stop 'Failed test_prep_reforis 21'
                            enddo
                        enddo
                        nnmat = test_os%nearest_neighbors( NSTATES )
                        do ref=1,p%nspace
                            oref = test_os%get_ori(ref)
                            do s=1,NSTATES
                                orefs = test_os%get_ori(nnmat(ref,s))
                                if( (orefs.euldist.oref)>.001 .or. nint(oref%get('class')).ne.nint(orefs%get('class')) &
                                    & .or. nint(orefs%get('state'))>NSTATES )stop 'Failed test_prep_reforis 22'
                            enddo
                        enddo
                        deallocate(nnmat)
                        call test_os%kill
                    enddo
                endif
                srch_order = primesrch3D%get_srch_order()
                if( minval(srch_order)<1 )stop 'Failed test_prep_reforis 23'
                if( maxval(srch_order)>primesrch3D%get_ntotrefs() )stop 'Failed test_prep_reforis 24'
                do i=1,primesrch3D%get_ntotrefs()
                    if( count(srch_order==i) /= 1 )stop 'Failed test_prep_reforis 25'
                enddo
                deallocate( srch_order )
            case('neigh','shcneigh')
                if( p%nstates==1 )then
                    do iptcl=1,p%nptcls
                        o = b%a%get_ori( iptcl )
                        call preprefs4align(b, p, iptcl, pftcc) 
                        call primesrch3D%prep4srch( o, b%nnmat )
                        srch_order = primesrch3D%get_srch_order()
                        if( minval(srch_order)<1 )stop 'Failed test_prep_reforis 32'
                        if( maxval(srch_order)>primesrch3D%get_ntotrefs() )stop 'Failed test_prep_reforis 33'
                        do ind=1,p%nnn
                            ref = srch_order(ind)
                            if( count(srch_order==ref) /= 1 )stop 'Failed test_prep_reforis 34'
                        enddo
                        deallocate( srch_order )
                    enddo
                endif
            case('qcontneigh')
                if( p%nstates==1 )then
                    do iptcl=1,p%nptcls
                        o = b%a%get_ori( iptcl )
                        call preprefs4align(b, p, iptcl, pftcc) 
                        call primesrch3D%prep4srch( o )
                        test_os = primesrch3D%get_o_refs( p%nnn )
                        do ref=1,p%nnn
                            oref = test_os%get_ori(ref)
                            if( deg2rad(o.euldist.oref) > 2.*p%athres )stop 'Failed test_prep_reforis 41'
                        enddo
                        call test_os%kill
                    enddo
                else
                    do iptcl=1,p%nptcls
                        o = b%a%get_ori( iptcl )
                        call preprefs4align(b, p, iptcl, pftcc) 
                        call primesrch3D%prep4srch( o )
                        test_os = primesrch3D%get_o_refs( p%nnn*p%nstates )
                        do ref=1,p%nnn*p%nstates
                            oref = test_os%get_ori(ref)
                            if( deg2rad(o.euldist.oref) > 2.*p%athres )stop 'Failed test_prep_reforis 51'
                        enddo
                        call test_os%kill
                    enddo
                endif      
            case DEFAULT
                stop'not implemented yet'
        end select
        if( verbose )print *,'end prep_reforis'
    end subroutine test_prep_reforis

    subroutine shutdown_testenv
        integer :: i
        call b%kill_general_tbox
        call b%kill_hadamard_prime3D_tbox
        do i=1,NPROJS
            call imgs_ptcls(i)%kill
        end do
        deallocate(imgs_ptcls)
    end subroutine shutdown_testenv

end module simple_prime3D_srch_tester



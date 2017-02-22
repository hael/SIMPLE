module simple_pftcc_contsrch
use simple_params,           only: params
use simple_build,            only: build
use simple_cmdline,          only: cmdline
use simple_opt_factory,      only: opt_factory
use simple_opt_spec,         only: opt_spec
use simple_optimizer,        only: optimizer
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_projector,        only: projector
use simple_ctf,              only: ctf
use simple_image,            only: image
use simple_ori,              only: ori
use simple_sym,              only: sym
use simple_defs
implicit none

public :: pftcc_contsrch_init, pftcc_contsrch_set_state, pftcc_contsrch_minimize
private

type(opt_factory)                :: ofac              !< optimizer factory
type(opt_spec)                   :: ospec             !< optimizer specification object
class(projector),        pointer :: proj_ptr=>null()  !< pointer to projector
class(params),           pointer :: p_ptr=>null()     !< pointer to params
class(build),            pointer :: b_ptr=>null()     !< pointer to builder
class(optimizer),        pointer :: nlopt=>null()     !< pointer to nonlinear optimizer
class(image),            pointer :: pimg=>null()      !< pointer to image
type(polarft_corrcalc)           :: pftcc
type(ctf)                        :: tfun              !<
type(ori)                        :: o_glob            !< global orientation
type(image),         allocatable :: refvols(:)
real,                allocatable :: ctfmat(:,:)
real                             :: dfx=0., dfy=0., angast=0.
integer                          :: state=1           !< state to evaluate
integer                          :: nstates=0        !< s
contains

    subroutine pftcc_contsrch_init( b, p, cline, img, opt_str, nrestarts )
        use simple_hadamard_common, only: preprefvol
        use simple_math,            only: round2even
        class(build),            target, intent(in) :: b
        class(params),           target, intent(in) :: p
        class(cmdline),                  intent(inout) :: cline
        class(image),            target, intent(in) :: img
        character(len=*),                intent(in) :: opt_str
        integer,                         intent(in) :: nrestarts
        real    :: lims_here(5,2)
        integer :: state
        ! set pointers
        b_ptr => b
        p_ptr => p
        proj_ptr  => b%proj
        pimg      => img
        nstates  = p_ptr%nstates
        ! init pftcc
        if( p%l_xfel )then
            call pftcc%new(1, [1,1], [p_ptr%boxmatch,p_ptr%boxmatch,1],&
            p_ptr%kfromto, p_ptr%ring2, p_ptr%ctf, isxfel='yes')
        else
            call pftcc%new(1, [1,1], [p_ptr%boxmatch,p_ptr%boxmatch,1],&
            p_ptr%kfromto, p_ptr%ring2, p_ptr%ctf)
        endif
        ! make optimizer spec
        lims_here = p%optlims(:5,:)
        call ospec%specify(opt_str, 5, ftol=1e-4, gtol=1e-4, limits=lims_here, nrestarts=nrestarts)
        ! set optimizer cost function
        call ospec%set_costfun(pftcc_contsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        ! volumes prep
        allocate( refvols(nstates) )
        do state=1,nstates
            call preprefvol( b_ptr, p_ptr, cline, state )
            refvols(state) = b%vol
        end do
    end subroutine pftcc_contsrch_init
    
    subroutine pftcc_contsrch_set_state( state_in )
        integer, intent(in) :: state_in
        state = state_in
    end subroutine pftcc_contsrch_set_state
    
    function pftcc_contsrch_get_nevals() result( nevals )
        integer :: nevals
        nevals = ospec%nevals
    end function pftcc_contsrch_get_nevals
    
    function pftcc_contsrch_cost( vec, D ) result( cost )
        use simple_ori, only: ori
        use simple_math, only: rad2deg
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        type(ori) :: o
        real      :: cost
        integer   :: i
        ! enforce the barrier constraint for the shifts
        do i=4,5
            if( vec(i) < ospec%limits(i,1) .or. vec(i) > ospec%limits(i,2) )then
                cost = 1.
                return
            endif
        end do
        ! calculate cost
        o = o_glob
        call o%set_euler(vec(1:3))
        ! shift vector
        call o%set('state', real(state))
        call proj_ptr%fproject_polar( 1, refvols(state), o, p_ptr, pftcc )
        if( p_ptr%ctf .ne. 'no' )then
            call pftcc%apply_ctf(tfun, dfx, dfy=dfy, angast=angast, ref=1,&
                &rot=1, ctfmat=ctfmat)
        endif
        cost = -pftcc%corr( 1, 1, 1, vec(4:5) )
    end function pftcc_contsrch_cost
    
    subroutine pftcc_contsrch_minimize( o )
        use simple_math, only: rad2deg, enforce_cyclic_limit
        class(ori), intent(inout) :: o
        real      :: corr, cost, dist
        real      :: prev_shvec(2)
        ! extract pft from ptcl
        call b_ptr%proj%img2polarft(1, pimg, pftcc)
        ! init CTF
        if( p_ptr%ctf .ne. 'no' )then
            ! CTF parms
            dfx = o%get('dfx')
            if( p_ptr%tfplan%mode.eq.'astig' )then ! astigmatic CTF
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else if( p_ptr%tfplan%mode.eq.'noastig' )then
                dfy    = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; simple_pftcc_contsrch%minimize'
            endif
            ! init CTF object
            tfun   = ctf( p_ptr%smpd, o%get('kv'), o%get('cs'), o%get('fraca'))
            ! fills and stores per ptcl polar CTF mat
            ctfmat = pftcc%create_polar_ctfmat(tfun, dfx, dfy, angast, pftcc%get_refsz())
        endif
        ! initial shift vector
        prev_shvec = o%get_shift()
        ! copy the input orientation
        o_glob = o
        ! initialise optimiser
        ospec%x = 0.
        ospec%x(1:3) = o%get_euler()
        ! search
        call nlopt%minimize(ospec, cost)
        ! report
        corr = -cost
        call o%set('corr', corr)
        ! shifts must be obtained by vector addition after rotation
        call o%set_euler(ospec%x(1:3))
        call o%set_shift( prev_shvec + ospec%x(4:5) )
        ! distance
        dist = 0.5*rad2deg(o_glob.euldist.o)+0.5*rad2deg(o_glob.inpldist.o)
        call o%set('dist',dist)
        if( p_ptr%ctf .ne. 'no' )then
            ! clean CTF exit
            dfx      = 0.
            dfy      = 0.
            angast   = 0.
            deallocate( ctfmat)
        endif
        ! clean exit
        state    = 1
    end subroutine pftcc_contsrch_minimize

end module simple_pftcc_contsrch

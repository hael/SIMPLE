module simple_pftcc_contsrch
use simple_params,           only: params
use simple_build,            only: build
use simple_cmdline,          only: cmdline
use simple_opt_factory,      only: opt_factory
use simple_opt_spec,         only: opt_spec
use simple_optimizer,        only: optimizer
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_projector,        only: projector
use simple_ori,              only: ori
use simple_defs
implicit none

public :: pftcc_contsrch_init, pftcc_contsrch_set_state, pftcc_contsrch_minimize
private

type(opt_factory)                :: ofac            !< optimizer factory
type(opt_spec)                   :: ospec           !< optimizer specification object
class(params),           pointer :: p_ptr =>null()  !< pointer to params
class(optimizer),        pointer :: nlopt =>null()  !< pointer to nonlinear optimizer
class(projector),        pointer :: pimg  =>null()  !< pointer to projector
type(polarft_corrcalc),  pointer :: ppftcc=>null()  !< polar FT calculator
type(ori)                        :: o_glob          !< global orientation
type(projector),     allocatable :: refvols(:)      !< reference volumes
integer                          :: state = 1       !< state to evaluate
logical,               parameter :: debug = .false.

contains

    subroutine pftcc_contsrch_init(b, p, cline, pftcc, img, opt_str, nrestarts)
        use simple_hadamard_common, only: preprefvol
        use simple_jiffys,          only: alloc_err, progress
        class(build),            target, intent(inout) :: b
        class(params),           target, intent(in)    :: p
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        class(cmdline),                  intent(inout) :: cline
        class(projector),        target, intent(in)    :: img
        character(len=*),                intent(in)    :: opt_str
        integer,                         intent(in)    :: nrestarts
        real    :: lims_here(5,2)
        integer :: s, alloc_stat
        ! set pointers
        p_ptr  => p
        ppftcc => pftcc
        pimg   => img
        ! init
        state  = 1
        ! make optimizer spec
        lims_here = p%optlims(:5,:)
        if(p%pgrp .ne. 'c1')lims_here(1:3,:) = b%se%srchrange() ! symmetry
        call ospec%specify(opt_str, 5, ftol=1e-4, gtol=1e-4, limits=lims_here, nrestarts=nrestarts)
        ! set optimizer cost function
        call ospec%set_costfun(pftcc_contsrch_cost)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        ! volumes prep
        if(allocated(refvols))deallocate(refvols)
        allocate( refvols(p_ptr%nstates), stat=alloc_stat )
        call alloc_err("In: pftcc_contsrch_init, simple_pftcc_contsrch", alloc_stat)
        write(*,'(A)') '>>> PREPARING VOLUMES'
        do s = 1, p_ptr%nstates
            call preprefvol( b, p_ptr, cline, s, doexpand=.false. )
            refvols(s) = b%vol
            call refvols(s)%expand_cmat
            call progress(s,p_ptr%nstates)
        end do
        call b%vol%kill_expanded
        if(debug)write(*,*)'pftcc_contsrch_init done'
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
            if(vec(i) < ospec%limits(i,1) .or. vec(i) > ospec%limits(i,2))then
                cost = 1.
                return
            endif
        end do
        ! calculate cost
        o = o_glob
        call o%set_euler(vec(1:3))
        call refvols(state)%fproject_polar(1, o, ppftcc, expanded=.true.)
        if(p_ptr%ctf .ne. 'no')call ppftcc%apply_ctf_single(1, 1)
        cost = -ppftcc%corr(1, 1, 1, vec(4:5))
    end function pftcc_contsrch_cost
    
    subroutine pftcc_contsrch_minimize( o )
        use simple_math, only: rad2deg
        use simple_oris, only: oris
        class(ori), intent(inout) :: o
        type(oris) :: a
        real :: corr, cost, dist, dist_inpl, prev_corr, frac
        real :: prev_shvec(2), dfx, dfy, angast
        ! extract pft from ptcl
        call pimg%img2polarft(1, ppftcc, isptcl=.true.)
        ! init CTF
        if(p_ptr%ctf .ne. 'no')then
            ! CTF parms
            dfx = o%get('dfx')
            if( p_ptr%tfplan%mode.eq.'astig' )then
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else if( p_ptr%tfplan%mode.eq.'noastig' )then
                dfy    = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; simple_pftcc_contsrch%minimize'
            endif
            call a%new(1)
            call a%set_ori(1, o)
            call ppftcc%create_polar_ctfmats(p_ptr%smpd, a)
            call a%kill
        endif
        ! initial shift vector
        prev_shvec = o%get_shift()
        ! copy the input orientation
        call o%set('state', real(state)) ! from pftcc_srch_set_state
        o_glob = o
        ! previous correlation
        call refvols(state)%fproject_polar(1, o, ppftcc, expanded=.true.)
        if(p_ptr%ctf .ne. 'no')call ppftcc%apply_ctf_single(1, 1)
        prev_corr = ppftcc%corr(1, 1, 1, [0.,0.])
        ! initialise optimiser
        ospec%x      = 0.
        ospec%x(1:3) = o%get_euler()
        ! search
        call nlopt%minimize(ospec, cost)
        corr = -cost
        ! report
        if(corr < prev_corr)then
            ! no improvement
            corr      = prev_corr
            dist      = 0.
            dist_inpl = 0.
            frac      = 100.
        else
            ! improvement
            call o%set_euler(ospec%x(1:3))
            ! shifts must be obtained by vector addition
            call o%set_shift(prev_shvec + ospec%x(4:5))
            ! distance
            dist_inpl = rad2deg(o_glob.inpldist.o)
            dist      = rad2deg(o_glob.euldist.o)
            frac      = 100.*(180.-(.5*dist+.5*dist_inpl)**2.)/180.
        endif
        ! sets new values
        call o%set('corr',      corr)
        call o%set('ow',        1.)
        call o%set('dist_inpl', dist_inpl)
        call o%set('dist',      dist)
        call o%set('mi_class',  1.)
        call o%set('mi_inpl',   1.)
        call o%set('mi_state',  1.)
        call o%set('mi_joint',  1.)
        call o%set('frac',      frac)
        call o%set('sdev',      0.)
        ! clean exit
        state = 1
    end subroutine pftcc_contsrch_minimize

end module simple_pftcc_contsrch

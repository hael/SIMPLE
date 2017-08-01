module simple_onflyft_srch
use simple_build,       only: build
use simple_params,      only: params
use simple_latent_alig, only: latent_alig
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_ori,         only: ori
use simple_defs         ! singleton
implicit none

public :: onflyft_srch_init, onflyft_srch_get_ori, onflyft_srch_align, onflyft_srch_kill,&
onflyft_srch_get_ori_best, onflyft_srch_corr_volimg
private

class(params), pointer    :: pp=>null()        !< pointer to parameters class
class(build), pointer     :: bp=>null()        !< pointer to build
type(opt_factory)         :: ofac              !< optimizer factory
type(opt_spec)            :: ospec             !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()     !< pointer to nonlinear optimizer
logical                   :: ctfastig=.false.  !< astigmatic CTF or not
logical                   :: phaseflip=.false. !< phaseflipped images or not
logical                   :: ctfrefine=.false. !< refine CTF or not
logical, parameter        :: report=.false.    !< reporting mode

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine onflyft_srch_init( b, p )
        use simple_jiffys, only: alloc_err
        class(build), intent(in), target     :: b
        class(params), intent(inout), target :: p
        ! kill pre-existing       
        call onflyft_srch_kill
        ! set constants
        pp => p
        bp => b
        ! make sure that the reference volumes are FTed
        do s=1,pp%nstates
            call bp%refvols(s)%fwd_ft
        end do
        ! set ctf modes 
        ctfastig  = .false.
        ctfrefine = .false.
        if( pp%ctf .eq. 'refine' )then
            ctfrefine = .true.
            if( bp%a%isthere('dfy') )then
                ctfastig = .true.
            endif
        endif
        phaseflip = .false.
        if( pp%phaseflip .eq. 'yes' ) phaseflip = .true.
        ! make optimizer spec
        if( p%diversify .eq. 'yes' )then
            call ospec%specify('simplex', pp%ndim, ftol=1e-4, gtol=1e-4,&
            limits=pp%optlims(:pp%ndim,:), cyclic=pp%cyclic(:pp%ndim), nrestarts=1)
        else
            call ospec%specify('powell', pp%ndim, ftol=1e-4, gtol=1e-4,&
            limits=pp%optlims(:pp%ndim,:), cyclic=pp%cyclic(:pp%ndim), nrestarts=1)
        endif
        ! set optimizer cost function
        call ospec%set_costfun(onflyft_srch_cost_volimg)
        ! set optimizer gradient function
        call ospec%set_gcostfun(onflyft_srch_gcost_volimg)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        if(report) write(*,*) 'did onflyft_srch_align_init'
    end subroutine

    ! ALIGNER ROUTINE
    
    !>  \brief  is the master aligner routine
    subroutine onflyft_srch_align( i, trial )
        use simple_ran_tabu, only: ran_tabu
        use simple_jiffys,   only: alloc_err
        integer, intent(in)                         :: i
        type(pori), intent(inout), target, optional :: trial
        class(ori), pointer                         :: optr=>null()
        type(ran_tabu)                              :: rt
        integer :: i, nbetter, alloc_stat, irnd
        ! make sure b%img FTed (assumes that b%img has been read, shifted, filtered, and masked)
        call bp%img%fwd_ft
        ! associate pointer
        if( present(trial) )then
            optr => trial
        else
            optr = bp%a%get_ori_ptr(i)
        endif
        ! check that the type is correct
        select type(optr)
            type is(pori)
                ! go ahead
            type is(ori)
                stop 'unsupported type: ori; onflyft_srch_align'
            class DEFAULT
                stop 'unsupported type; onflyft_srch_align'
        end select
        ! OASIS MC
        call optr%oasis_mc(ctfrefine, onflyft_srch_cost_volimg, pp%refine)
        ! NLP refinement
        nbetter = optr%get_nbetter()
        if( pp%diversify .eq. 'yes' )then
            ! refine randomly select ones (with simplex randomized bounds)
            if( nbetter <= 3 )then
                do i=1,3
                    call optr%refine_minimum(nlopt, ospec, i)
                end do
            else
                call alloc_err("In: onflyft_srch_align", alloc_stat)
                rt = ran_tabu(nbetter)
                do i=1,3
                    irnd = rt%irnd()
                    call rt%insert(irnd)
                    call optr%refine_minimum(nlopt, ospec, irnd)
                end do
                call rt%kill
            endif
        else
            ! refine the best one (with powell)
            call optr%refine_minimum(nlopt, ospec, i)
        endif
        ! calculate orientation spread
        call optr%calc_orispread
        ! calculate orientation weights
        call optr%calc_weights
        ! decimate solution poulation with binary random tournament
        call optr%binary_tournament
        ! we're done, nullify the pointer
        optr => null
        if(report) write(*,*) 'did onflyft_srch_align'
    end subroutine
    
    !>  \brief  calculates the vol vs. img correlation for input orientation
    function onflyft_srch_corr_volimg( orientation, state ) result( corr )
        use simple_math, only: deg2rad
        class(ori), intent(inout)     :: orientation
        integer, intent(in), optional :: state
        real                          :: x, y, corr
        integer                       :: sstate
        sstate = 1 
        if( present(state) )then
            sstate = state
        else
            if( pp%nstates > 1 ) sstate = nint(orientation%get('state'))
        endif
        call bp%proj%fproject(bp%refvols(sstate), orientation, bp%img_ref, lp_dyn=pp%lp) ! project
        x = orientation%get('x')                                                         ! get x shift
        y = orientation%get('y')                                                         ! get y shift
        if( pp%ctf .eq. 'yes' .and. pp%ndim <= 5 )then ! multiply reference with CTF**2
            call bp%img_ref%mul(bp%img_ctfsq, pp%lp)
        endif
        if( pp%ndim > 5 )then ! CTF refinement
            bp%img_sh = bp%img
            ! multiply ptcl & ref with CTF
            if( ctfastig )then
                if( phaseflip )then
                    call bp%tfun%apply(bp%img_sh, orientation%get('dfx'), 'abs',&
                    orientation%get('dfy'), deg2rad(orientation%get('angast')))
                else
                    call bp%tfun%apply(bp%img_sh, orientation%get('dfx'), 'ctf',&
                    orientation%get('dfy'), deg2rad(orientation%get('angast')))
                endif
                call bp%tfun%ctf2img(bp%img_ctfsq, orientation%get('dfx'), 'square',&
                orientation%get('dfy'), deg2rad(orientation%get('angast')))
                call bp%img_ref%mul(bp%img_ctfsq, pp%lp)
            else
                if( phaseflip )then
                    call bp%tfun%apply(bp%img_sh, orientation%get('dfx'), 'abs')
                else
                    call bp%tfun%apply(bp%img_sh, orientation%get('dfx'), 'ctf')
                endif
                call bp%tfun%ctf2img(bp%img_ctfsq, orientation%get('dfx'), 'square')
                call bp%img_ref%mul(bp%img_ctfsq, pp%lp)
            endif
            call bp%img_sh%shift(x, y)          ! shift ptcl, no rev sign here (stupid convention, but what the heck)
            if( allocated(bp%wiener) ) call bp%img%apply_filter(bp%wiener(sstate,:))
            call bp%img_sh%bwd_ft               ! reverse FT ptcl
            call bp%img_sh%mask(pp%msk, 'soft') ! mask ptcl
            call bp%img_sh%fwd_ft               ! fwd FT ptcl
        else
            call bp%img%shift(x, y, imgout=bp%img_sh)
        endif
        corr = bp%img_ref%corr(bp%img_sh,pp%lp,pp%hp) ! correlate
    end function
    
    !>  \brief  calculates the vol vs. img cost for input vector
    !!          (remember to set state beforehand)
    function onflyft_srch_cost_volimg( vec, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real                :: cost
        call pp%ori_glob%set_euler([vec(1),vec(2),vec(3)])
        call pp%ori_glob%set('x',vec(4))
        call pp%ori_glob%set('y',vec(5))
        if( D > 5 )then
            call pp%ori_glob%set('dfx',vec(6))
            if( D > 6 ) call pp%ori_glob%set('dfy',vec(7))
        endif
        cost = -onflyft_srch_corr_volimg(pp%ori_glob)
    end function
    
    !>  \brief  calculates the vol vs. img cost gradient for input vector
    function onflyft_srch_gcost_volimg( vec, D ) result( grad )
        use simple_math, only: numderiv
        integer, intent(in) :: D
        real, intent(inout) :: vec(D)
        real                :: grad(D),err(D),h(D)
        integer             :: which_dim, i
        h(1:3) = 5.
        h(4:5) = 1.
        do i=1,D
            which_dim = i
            call numderiv(func_onedim, vec(i), h(i), err(i), grad(i))
        end do

        contains
        
            !>  \brief  one-dimensional version of the multidimensional func
            function func_onedim( t ) result( r )
                real, intent(in) :: t
                real :: r, s
                ! store the old parameter value
                s = vec(which_dim)
                ! replace
                vec(which_dim) = t
                ! evaluate func
                r = onflyft_srch_cost_volimg( vec, D ) 
                ! put old val back
                vec(which_dim) = s
            end function

    end function
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine onflyft_srch_kill
        pp   => null()
        bp   => null()
        call ofac%kill
        call ospec%kill
        nlopt => null()
    end subroutine

end module simple_onflyft_srch
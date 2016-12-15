module simple_onflyft_srch
use simple_build,       only: build
use simple_params,      only: params
use simple_latent_alig, only: latent_alig
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_ori,         only: ori
use simple_defs         ! singleton
use simple_cmdline      ! singleton
implicit none

public :: onflyft_srch_init, onflyft_srch_align,  onflyft_srch_corr_volimg, onflyft_srch_kill
private

class(params), pointer    :: pp=>null()        !< pointer to parameters class
class(build), pointer     :: bp=>null()        !< pointer to build
type(opt_factory)         :: ofac              !< optimizer factory
type(opt_spec)            :: ospec             !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()     !< pointer to nonlinear optimizer
logical, parameter        :: report=.false.    !< reporting mode

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine onflyft_srch_init( b, p )
        use simple_jiffys, only: alloc_err
        class(build), intent(in), target     :: b
        class(params), intent(inout), target :: p
        integer :: s
        ! kill pre-existing
        call onflyft_srch_kill
        ! set constants
        pp => p
        bp => b
        ! make sure that the reference volumes are FTed
        do s=1,pp%nstates
            call bp%refvols(s)%fwd_ft
        end do
        ! make optimizer spec
        call ospec%specify('simplex', pp%ndim, ftol=1e-4, gtol=1e-4,&
        limits=pp%optlims(:pp%ndim,:), cyclic=pp%cyclic(:pp%ndim), nrestarts=1)
        ! set optimizer cost function
        call ospec%set_costfun(onflyft_srch_cost_volimg)
        ! set optimizer gradient function
        call ospec%set_gcostfun(onflyft_srch_gcost_volimg)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
    end subroutine

    ! ALIGNER ROUTINE
    
    !>  \brief  is the master aligner routine
    subroutine onflyft_srch_align( i, po, trial )
        use simple_ran_tabu, only: ran_tabu
        use simple_pori,     only: pori
        integer, intent(in)                :: i
        type(pori), intent(inout)          :: po
        type(ori), intent(inout), optional :: trial
        type(ori)                          :: o, o_old
        type(ran_tabu)                     :: rt
        character(len=STDLEN)              :: dig
        integer                            :: j, nbetter, irnd, fromto(2)
        logical                            :: here
        real                               :: corr, corr_old, sdevs(pp%ndim)
        integer, parameter                 :: NREFINE=2
        real                               :: angsdev
        ! make sure b%img FTed (assumes that b%img has been read, shifted, filtered, and masked)
        call bp%img%fwd_ft
        ! define pp%porifile name
        if( defined_cmd_arg('part') )then
            write(dig,*) pp%part
            pp%porifile='poris_part'//trim(adjustl(dig))//'.bin'
        endif
        ! check for pp%porifile
        inquire(file=pp%porifile, exist=here)
!        if( here )then
!        else
!            if(report) write(*,*) '*** onflyft_srch ***: trying to create porifile'
!            call create_porifile
!        endif
        ! set o
        if( present(trial) )then
            o = trial
            o_old = bp%a%get_ori(i)
            corr_old = onflyft_srch_corr_volimg(o_old, nint(o_old%get('state'))) 
        else
            o = bp%a%get_ori(i)
        endif
        ! create new pori object
!         call po%new(bp%porispec,o,pp%ori_glob)
        if( report ) write(*,*) '*** onflyft_srch ***: created new pori object'
        if( here )then
            ! do nothing
        else
            fromto(1) = 1
            fromto(2) = pp%nptcls
!            if( defined_cmd_arg('part') )then
!                fromto(1) = pp%fromp
!                fromto(2) = pp%top
!            else
!                fromto(1) = 1
!                fromto(2) = pp%nptcls
!            endif
            angsdev = max(5.,pp%tres)
            do j=fromto(1),fromto(2)
                o = bp%a%get_ori(j)
                if( report ) write(*,*) '*** onflyft_srch ***: got ori'
!                 call po%new(bp%porispec,o,pp%ori_glob)
                if( report ) write(*,*) '*** onflyft_srch ***: created new pori'
                if( pp%refine .eq. 'het' )then
!                     call po%oasis_init(ospec, sdevs=[angsdev,angsdev,angsdev])
                else
!                     call po%oasis_init(ospec) ! BUGS OUT WHEN EO=YES BUT NOT WHEN EO=NO & LP=X
                endif
                if( report ) write(*,*) '*** onflyft_srch ***: initialized oasis'
                call po%write(pp%porifile, j)
                if( report ) write(*,*) '*** onflyft_srch ***: wrote pori 2 file'
            end do
            if(report) write(*,*) '*** onflyft_srch ***: created porifile'
        endif
        ! re-create pori object from file
!         call po%read(pp%porifile, i, o, pp%ori_glob)
        if( report ) write(*,*) '*** onflyft_srch ***: re-created pori object from file'
        ! OASIS MC
!         call po%oasis_mc(ospec, pp%refine)
        if( report ) write(*,*) '*** onflyft_srch ***: did oasis mc'
        ! NLP refinement
        nbetter = po%get_nbetter()
        call o%set('nmin', real(nbetter))
        if( report ) write(*,*) '*** onflyft_srch ***: NBETTER(onflyft_srch)=', nbetter
        if( pp%refine .ne. 'het' .and. pp%refine .ne. 'stoch' )then
            ! refine randomly select ones (with simplex randomized bounds)
            if( nbetter <= NREFINE )then
                do j=1,nbetter
!                     call po%refine_minimum(nlopt, ospec, j)
                end do
            else
                rt = ran_tabu(nbetter)
                do j=1,NREFINE
                    irnd = rt%irnd()
                    call rt%insert(irnd)
!                     call po%refine_minimum(nlopt, ospec, irnd)
                end do
                call rt%kill
            endif
        endif
        ! calculate orientation weights
        call po%calc_weights
        ! decimate solution poulation with binary random tournament
!!!!!!!!! call po%binary_tournament
        ! check if alignment improved
        if( present(trial) )then
            ! check if optimization attempt was successful
            corr = o%get('corr')
            if( corr > corr_old )then
                call bp%a%set_ori(i,o)
            endif
        else
            sdevs = po%get_sdevs()
            call o%set('sdev', sum(sdevs(1:3))/3.)
            call bp%a%set_ori(i,o)
        endif
        ! update the porifile
        call po%write(pp%porifile, i)
        ! the orientation parameters are written to file in matcher
        if(report) write(*,*) '*** onflyft_srch ***: did onflyft_srch_align'
        
!        contains
!            
!            subroutine create_porifile
!                integer :: i, fromto(2)
!                real    :: angsdev
!                angsdev = max(5.,pp%tres)
!                if( defined_cmd_arg('part') )then
!                    fromto(1) = pp%fromp
!                    fromto(2) = pp%top
!                else
!                    fromto(1) = 1
!                    fromto(2) = pp%nptcls
!                endif
!                do i=fromto(1),fromto(2)
!                    o = bp%a%get_ori(i)
!                    if( report ) write(*,*) '*** onflyft_srch ***: got ori'
!                    call po%new(bp%porispec,o,pp%ori_glob)
!                    if( report ) write(*,*) '*** onflyft_srch ***: created new pori'
!                    if( pp%refine .eq. 'het' )then
!                        call po%oasis_init(ospec, sdevs=[angsdev,angsdev,angsdev])
!                    else
!                        call po%oasis_init(ospec) ! BUGS OUT WHEN EO=YES BUT NOT WHEN EO=NO & LP=X
!                    endif
!                    if( report ) write(*,*) '*** onflyft_srch ***: initialized oasis'
!                    call po%write(pp%porifile, i)
!                    if( report ) write(*,*) '*** onflyft_srch ***: wrote pori 2 file'
!                end do
!                if(report) write(*,*) '*** onflyft_srch ***: created porifile'
!            end subroutine
        
    end subroutine
    
    !>  \brief  calculates the vol vs. img correlation for input orientation
    function onflyft_srch_corr_volimg( orientation, state ) result( corr )
        use simple_math, only: deg2rad
        class(ori), intent(inout)     :: orientation
        integer, intent(in), optional :: state
        real                          :: x, y, corr
        integer                       :: sstate
        if( report ) write(*,*) '*** onflyft_srch_corr_volimg ***: inside'
        sstate = 1 
        if( present(state) )then
            sstate = state
        else
            if( pp%nstates > 1 ) sstate = nint(orientation%get('state'))
        endif
        if( report )then
            write(*,*) '*** onflyft_srch_corr_volimg *** ldim refvol:', bp%refvols(sstate)%get_ldim()
                write(*,*) '*** onflyft_srch_corr_volimg *** ldim img_ref:', bp%img_ref%get_ldim()
            write(*,*) '*** onflyft_srch_corr_volimg *** lp_dyn:', pp%lp
        endif
        call bp%proj%fproject(bp%refvols(sstate), orientation, bp%img_ref, lp_dyn=pp%lp) ! project
        if( report ) write(*,*) '*** onflyft_srch ***: projected'
        x = orientation%get('x')                                                         ! get x shift
        y = orientation%get('y')                                                         ! get y shift
        if( pp%ndim > 5 )then ! CTF refinement
            if( pp%ctf .eq. 'flip' ) stop 'Cannot use phase-flipped images for CTF refinement;  onflyft_srch_corr_volimg'     
            bp%img_sh = bp%img
            ! multiply ptcl with CTF
            if( pp%ctfmode .eq. 'astig' )then
                call bp%tfun%apply(bp%img_sh, orientation%get('dfx'), 'ctf',&
                orientation%get('dfy'), orientation%get('angast'))
            else if( pp%ctfmode .eq. 'noastig' )then
                call bp%tfun%apply(bp%img_sh, orientation%get('dfx'), 'ctf')
            endif
            call bp%img_sh%shift(x, y)          ! shift ptcl, no rev sign here (stupid convention, but what the heck)
            call bp%img_sh%bwd_ft               ! reverse FT ptcl
            ! clip image if needed
            if( pp%boxmatch < pp%box ) call bp%img_sh%clip_inplace([pp%boxmatch,pp%boxmatch,1]) ! SQUARE DIMS ASSUMED
            call bp%img_sh%mask(pp%msk, 'soft') ! mask ptcl
            call bp%img_sh%fwd_ft               ! fwd FT ptcl
        else
            call bp%img%shift(x, y, imgout=bp%img_sh)
            if( report ) write(*,*) '*** onflyft_srch ***: shifted'
        endif
        corr = bp%img_ref%corr(bp%img_sh,pp%lp,pp%hp) ! correlate
        if( report ) write(*,*) '*** onflyft_srch ***: correlated'
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

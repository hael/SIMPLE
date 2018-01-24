module simple_ctffit
#include "simple_lib.f08"
use simple_image,       only: image
use simple_ctf,         only: ctf
use simple_opt_spec,    only: opt_spec
use simple_de_opt,      only: de_opt
use simple_simplex_opt, only: simplex_opt
implicit none

public :: ctffit_init, ctffit_x_validated_fit, ctffit_kill
private

class(image), pointer :: ppspec_all              ! all micrograph powerspec
class(image), pointer :: ppspec_lower            ! lower half of micrograph powerspec
class(image), pointer :: ppspec_upper            ! upper half of micrograph powerspec
class(image), pointer :: ppspec_ref              ! powerspec fitting reference
type(image)           :: pspec_ctf               ! CTF powerspec
type(image)           :: imgmsk                  ! mask image
type(ctf)             :: tfun                    ! transfer function object
type(opt_spec)        :: ospec_de                ! optimiser specification differential evolution (DE)
type(opt_spec)        :: ospec_simplex           ! optimiser specification Nelder-Mead (N-M)
type(de_opt)          :: diffevol                ! DE search object
type(simplex_opt)     :: simplexsrch             ! N-M search object
logical, allocatable  :: cc_msk(:,:,:)           ! corr mask
logical               :: l_phaseplate = .false.  ! Volta phase-plate flag
integer               :: ndim         = 3        ! # optimisation dims
real                  :: df_min       = 0.5      ! close 2 focus limit
real                  :: df_max       = 5.0      ! far from focus limit
real                  :: sxx          = 0.       ! memoized corr term
real                  :: limits(4,2)  = 0.       ! barrier limits

contains

    subroutine ctffit_init( pspec_all, pspec_lower, pspec_upper, smpd, kV, Cs, amp_contr, dfrange, resrange, phaseplate )
        class(image), target, intent(inout) :: pspec_all   !< all micrograph powerspec
        class(image), target, intent(inout) :: pspec_lower !< lower half of micrograph powerspec
        class(image), target, intent(inout) :: pspec_upper !< upper half of micrograph powerspec
        real,                 intent(in)    :: smpd        !< sampling distance
        real,                 intent(in)    :: kV          !< acceleration voltage
        real,                 intent(in)    :: Cs          !< constant
        real,                 intent(in)    :: amp_contr   !< amplitude contrast
        real,                 intent(in)    :: dfrange(2)  !< defocus range, [30.0,5.0] default
        real,                 intent(in)    :: resrange(2) !< resolution range, [30.0,5.0] default
        character(len=*),     intent(in)    :: phaseplate  !< Volta phase-plate images (yes|no)
        real    :: hp, lp
        integer :: ldim(3)
        ! set constants
        ppspec_all   => pspec_all
        ppspec_lower => pspec_lower
        ppspec_upper => pspec_upper
        ldim = pspec_all%get_ldim()
        if( dfrange(1) < dfrange(2) )then
            df_min = dfrange(1)
            df_max = dfrange(2)
        else
            stop 'invalid defocus range; simple_ctffit :: new'
        endif
        if( resrange(1) > resrange(2) )then
            hp = resrange(1)
            lp = resrange(2)
        else
            stop 'invalid resolution range; simple_ctffit :: new'
        endif
        select case(trim(phaseplate))
            case('yes')
                l_phaseplate = .true.
                ndim = 4
            case DEFAULT
                ndim = 3
        end select
        ! construct CTF object
        tfun = ctf(smpd, kV, Cs, amp_contr)
        ! prepare powerspectra
        call ppspec_all%dampen_central_cross
        call ppspec_all%subtr_backgr(hp)
        call ppspec_lower%dampen_central_cross
        call ppspec_lower%subtr_backgr(hp)
        call ppspec_upper%dampen_central_cross
        call ppspec_upper%subtr_backgr(hp)
        call pspec_ctf%new(ldim, smpd)
        ! generate correlation mask
        call imgmsk%new(ldim, smpd)
        call imgmsk%resmsk(hp, lp)
        cc_msk = imgmsk%bin2logical()
        ! contruct optimisers
        call seed_rnd
        limits        = 0.
        limits(1:2,1) = df_min
        limits(1:2,2) = df_max
        limits(3,1)   = 0.
        limits(3,2)   = twopi  ! miminise in radians so that the df:s are roughly on the same scale
        if( l_phaseplate )then
            limits(4,1) = 0.
            limits(4,2) = 3.15
        endif
        call ospec_de%specify('de', ndim, limits=limits(1:ndim,:), maxits=400)
        call ospec_simplex%specify('simplex', ndim, limits=limits(1:ndim,:), maxits=80, nrestarts=3)
        if( l_phaseplate )then
            call ospec_de%set_costfun(ctffit_cost_phaseplate)
            call ospec_simplex%set_costfun(ctffit_cost_phaseplate_barrier)
        else
            call ospec_de%set_costfun(ctffit_cost)
            call ospec_simplex%set_costfun(ctffit_cost_barrier)
        endif
        call diffevol%new(ospec_de)
        call simplexsrch%new(ospec_de)
    end subroutine ctffit_init

    subroutine ctffit_x_validated_fit( dfx, dfy, angast, phshift, dferr, cc, diagfname )
        real,             intent(out) :: dfx, dfy, angast, phshift, dferr, cc
        character(len=*), intent(in)  :: diagfname
        type(image) :: pspec_half_n_half
        real        :: dfavg, dfavg_lower, dfavg_upper, dfx_lower, dfx_upper, dfy_lower
        real        :: dfy_upper, angast_lower, angast_upper, phshift_lower, phshift_upper
        real        :: cc_lower, cc_upper
        ! determine parameters based on lower half of micrograph
        call init_srch( ppspec_lower )
        call srch(   dfx_lower, dfy_lower, angast_lower, phshift_lower, cc_lower )
        call refine( dfx_lower, dfy_lower, angast_lower, phshift_lower, cc_lower )
        ! determine parameters based on upper half of micrograph
        call init_srch( ppspec_upper )
        call srch(   dfx_upper, dfy_upper, angast_upper, phshift_upper, cc_upper )
        call refine( dfx_upper, dfy_upper, angast_upper, phshift_upper, cc_upper )
        ! check which solution fits the global powerspec best
        call init_srch( ppspec_all )
        if( l_phaseplate )then
            call tfun%ctf2pspecimg(pspec_ctf, dfx_lower, dfy_lower, angast_lower, add_phshift=phshift_lower)
            cc_lower = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
            call tfun%ctf2pspecimg(pspec_ctf, dfx_upper, dfy_upper, angast_upper, add_phshift=phshift_lower)
            cc_upper = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
        else
            call tfun%ctf2pspecimg(pspec_ctf, dfx_lower, dfy_lower, angast_lower)
            cc_lower = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
            call tfun%ctf2pspecimg(pspec_ctf, dfx_upper, dfy_upper, angast_upper)
            cc_upper = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
        endif
        if( cc_lower >= cc_upper )then
            dfx     = dfx_lower
            dfy     = dfy_lower
            angast  = angast_lower
            phshift = phshift_lower
        else
            dfx     = dfx_upper
            dfy     = dfy_upper
            angast  = angast_upper
            phshift = phshift_upper
        endif
        ! refine the best solution on the global powerspec
        call refine( dfx, dfy, angast, phshift, cc )
        ! calculate error
        dfavg       = (dfx + dfy / 2)
        dfavg_lower = (dfx_lower + dfy_lower / 2)
        dfavg_upper = (dfx_upper + dfy_upper / 2)
        dferr       = maxval([abs(dfavg - dfavg_lower), abs(dfavg - dfavg_upper), abs(dfavg_lower - dfavg_upper)])
        ! make a half-n-half diagnostic
        if( l_phaseplate )then
            call tfun%ctf2pspecimg(pspec_ctf, dfx, dfy, angast, add_phshift=phshift)
        else
            call tfun%ctf2pspecimg(pspec_ctf, dfx, dfy, angast)
        endif
        call pspec_ctf%norm
        call ppspec_ref%norm
        call pspec_ctf%mul(imgmsk)
        pspec_half_n_half = ppspec_ref%before_after(pspec_ctf, cc_msk)
        call pspec_half_n_half%write(trim(diagfname), 1)
        call pspec_half_n_half%kill
    end subroutine ctffit_x_validated_fit

    subroutine init_srch( which_pspec )
        class(image), target :: which_pspec
        ppspec_ref => which_pspec
        ! memoize reference corr components
        call ppspec_ref%prenorm4real_corr(sxx, cc_msk)
    end subroutine init_srch

    subroutine srch( dfx, dfy, angast, phshift, cc )
        real, intent(out) :: dfx, dfy, angast, phshift, cc
        real              :: cost, df, cost_lowest, dfstep
        class(*), pointer :: fun_self => null()
        dfstep = (df_max - df_min) / 100.
        if( l_phaseplate )then
            ! do a first grid search assuming:
            ! no astigmatism
            ! pi half phase shift
            df = df_min
            cost_lowest = ctffit_cost_phaseplate(fun_self, [df,df,0.,PIO2], ndim)
            do while( df <= df_max )
                cost = ctffit_cost_phaseplate(fun_self, [df,df,0.,PIO2], ndim)
                if( cost < cost_lowest )then
                    cost_lowest = cost
                    ospec_de%x = [df,df,0.,PIO2]
                endif
                df = df + dfstep
            end do
        else
            ! do a first grid search assuming:
            ! no astigmatism
            df = df_min
            cost_lowest = ctffit_cost(fun_self, [df,df,0.], ndim)
            do while( df <= df_max )
                cost = ctffit_cost(fun_self, [df,df,0.], ndim)
                if( cost < cost_lowest )then
                    cost_lowest = cost
                    ospec_de%x  = [df,df,0.]
                endif
                df = df + dfstep
            end do
        endif
        ! refinement by DE (Differential Evolution)
        call diffevol%minimize(ospec_de, fun_self, cost)
        dfx     = ospec_de%x(1)
        dfy     = ospec_de%x(2)
        angast  = rad2deg(ospec_de%x(3))
        cc      = -cost
        phshift = 0.
        if( l_phaseplate ) phshift = ospec_de%x(4)
    end subroutine srch

    subroutine refine( dfx, dfy, angast, phshift, cc )
        real, intent(inout) :: dfx, dfy, angast, phshift, cc
        real                :: cost
        class(*), pointer   :: fun_self => null()
        if( l_phaseplate )then
            ospec_simplex%x = [dfx,dfy,angast,phshift]
        else
            ospec_simplex%x = [dfx,dfy,angast]
        endif
        ! refinement with unconstrained Nelder-Mead: critical to accuracy
        call simplexsrch%minimize(ospec_simplex, fun_self, cost)
        ! report solution
        dfx     = ospec_simplex%x(1)
        dfy     = ospec_simplex%x(2)
        angast  = rad2deg(ospec_simplex%x(3))
        cc      = -cost
        phshift = 0.
        if( l_phaseplate ) phshift = ospec_de%x(4)
    end subroutine refine

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    function ctffit_cost( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: cost
        call tfun%ctf2pspecimg(pspec_ctf, vec(1), vec(2), rad2deg(vec(3)))
        cost = -ppspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
    end function ctffit_cost

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    ! with barrier constraint on the astigmatism angle [0,twopi]
    function ctffit_cost_barrier( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: cost
        if( vec(3) < limits(3,1) .or. vec(3) > limits(3,2) )then
            cost = 1.0
            return
        endif
        call tfun%ctf2pspecimg(pspec_ctf, vec(1), vec(2), rad2deg(vec(3)))
        cost = -ppspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
    end function ctffit_cost_barrier

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    function ctffit_cost_phaseplate( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        ! vec(4) is additional phase shift (in radians)
        call tfun%ctf2pspecimg(pspec_ctf, vec(1), vec(2), rad2deg(vec(3)), add_phshift=vec(4))
        cost = -ppspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
    end function ctffit_cost_phaseplate

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    ! with barrier constraint on the astigmatism angle [0,twopi]
    ! with barrier constraint on the phase change      [0,3.15]
    function ctffit_cost_phaseplate_barrier( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        ! vec(4) is additional phase shift (in radians)
        if( vec(3) < limits(3,1) .or. vec(3) > limits(3,2) .or. vec(4) < limits(4,1) .or. vec(4) > limits(4,2) )then
            cost = 1.0
            return
        endif
        call tfun%ctf2pspecimg(pspec_ctf, vec(1), vec(2), rad2deg(vec(3)), add_phshift=vec(4))
        cost = -ppspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
    end function ctffit_cost_phaseplate_barrier

    ! with barrier constraint on the astigmatism angle [0,twopi]

    subroutine ctffit_kill
        ppspec_ref   => null()
        ppspec_all   => null()
        ppspec_lower => null()
        ppspec_upper => null()
        call pspec_ctf%kill
        call ospec_de%kill
        call ospec_simplex%kill
        call diffevol%kill
        call simplexsrch%kill
        call imgmsk%kill
        if( allocated(cc_msk) ) deallocate(cc_msk)
    end subroutine ctffit_kill

end module simple_ctffit

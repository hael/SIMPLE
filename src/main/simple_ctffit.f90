module simple_ctffit
#include "simple_lib.f08"
use simple_image,    only: image
use simple_ctf,      only: ctf
use simple_opt_spec, only: opt_spec
use simple_de_opt,   only: de_opt
implicit none

public :: ctffit_init, ctffit_x_validated_fit, ctffit_kill
private

class(image), pointer :: ppspec_all              ! all micrograph powerspec
class(image), pointer :: ppspec_lower            ! lower half of micrograph powerspec
class(image), pointer :: ppspec_upper            ! upper half of micrograph powerspec
class(image), pointer :: ppspec_ref              ! powerspec fitting reference
class(image), pointer :: ppspec_ref_roavg        ! rotationally averaged powerspec fitting reference
type(image)           :: pspec_ctf               ! CTF powerspec
type(image)           :: pspec_ctf_roavg         ! rotationally averaged CTF powerspec
type(image)           :: pspec_all_roavg         ! rotationally averaged all micrograph powerspec
type(image)           :: pspec_lower_roavg       ! rotationally averaged lower half of micrograph powerspec
type(image)           :: pspec_upper_roavg       ! rotationally averaged upper half of micrograph powerspec
type(image)           :: imgmsk                  ! mask image
type(ctf)             :: tfun                    ! transfer function object
type(ctf)             :: tfun_roavg              ! rotationally averaged transfer function object
type(opt_spec)        :: ospec_de                ! optimiser specification differential evolution (DE)
type(de_opt)          :: diffevol                ! DE search object
logical, allocatable  :: cc_msk(:,:,:)           ! corr mask
logical               :: l_phaseplate = .false.  ! Volta phase-plate flag
integer               :: ndim         = 3        ! # optimisation dims
real                  :: df_min       = 0.5      ! close 2 focus limit
real                  :: df_max       = 5.0      ! far from focus limit
real                  :: sxx          = 0.       ! memoized corr term
real                  :: hp           = 0.       ! high-pass limit
real                  :: lp           = 0.       ! low-pass limit
real                  :: sxx_roavg    = 0.       ! memoized corr term, rotationally averaged power spec
real                  :: limits(4,2)  = 0.       ! barrier limits

integer, parameter :: IARES = 10, NSTEPS = 100

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
                l_phaseplate = .false.
                ndim = 3
        end select
        ! construct CTF objects
        tfun       = ctf(smpd, kV, Cs, amp_contr)
        tfun_roavg = ctf(smpd, kV, Cs, amp_contr)
        ! prepare powerspectra
        call ppspec_all%dampen_central_cross
        call ppspec_all%subtr_backgr(hp)
        call ppspec_lower%dampen_central_cross
        call ppspec_lower%subtr_backgr(hp)
        call ppspec_upper%dampen_central_cross
        call ppspec_upper%subtr_backgr(hp)
        ! prepare rotationally averaged power spectra
        call ppspec_all%roavg(IARES,   pspec_all_roavg)
        call ppspec_lower%roavg(IARES, pspec_lower_roavg)
        call ppspec_upper%roavg(IARES, pspec_upper_roavg)
        ! prepare CTF power spectra
        call pspec_ctf%new(ldim, smpd)
        call pspec_ctf_roavg%new(ldim, smpd)
        ! generate correlation mask
        call imgmsk%new(ldim, smpd)
        call imgmsk%resmsk(hp, lp)
        cc_msk = imgmsk%bin2logical()
        ! contruct DE optimiser
        call seed_rnd
        limits        = 0.
        limits(1:2,1) = df_min
        limits(1:2,2) = df_max
        limits(3,1)   = -180.
        limits(3,2)   = 180.
        if( l_phaseplate )then
            limits(4,1) = 0.
            limits(4,2) = 3.15
        endif
        call ospec_de%specify('de', ndim, limits=limits(1:ndim,:), maxits=400, ftol=1e-4)
        if( l_phaseplate )then
            call ospec_de%set_costfun(ctffit_cost_phaseplate)
        else
            call ospec_de%set_costfun(ctffit_cost)
        endif
        call diffevol%new(ospec_de)
    end subroutine ctffit_init

    subroutine ctffit_x_validated_fit( dfx, dfy, angast, phshift, dferr, cc, ctfscore, diagfname )
        real,             intent(out) :: dfx, dfy, angast, phshift, dferr, cc, ctfscore
        character(len=*), intent(in)  :: diagfname
        real, allocatable :: corrs(:)
        type(image)       :: pspec_half_n_half
        real              :: dfavg, dfavg_lower, dfavg_upper, dfx_lower, dfx_upper, dfy_lower
        real              :: dfy_upper, angast_lower, angast_upper, phshift_lower, phshift_upper
        real              :: cc_lower, cc_upper, df_avg
        integer           :: filtsz, k, hpfind, lpfind
        ! determine parameters based on lower half of micrograph
        call init_srch( ppspec_lower, pspec_lower_roavg )
        call grid_srch(   dfx_lower, dfy_lower, angast_lower, phshift_lower, cc_lower )
        call refine( dfx_lower, dfy_lower, angast_lower, phshift_lower, cc_lower )
        ! determine parameters based on upper half of micrograph
        call init_srch( ppspec_upper, pspec_upper_roavg )
        call grid_srch(   dfx_upper, dfy_upper, angast_upper, phshift_upper, cc_upper )
        call refine( dfx_upper, dfy_upper, angast_upper, phshift_upper, cc_upper )
        ! check which solution fits the global powerspec best
        if( l_phaseplate )then
            call ctf2pspecimg(tfun, pspec_ctf, dfx_lower, dfy_lower, angast_lower, add_phshift=phshift_lower)
            cc_lower = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
            call ctf2pspecimg(tfun, pspec_ctf, dfx_upper, dfy_upper, angast_upper, add_phshift=phshift_lower)
            cc_upper = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
        else
            call ctf2pspecimg(tfun, pspec_ctf, dfx_lower, dfy_lower, angast_lower)
            cc_lower = ppspec_all%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
            call ctf2pspecimg(tfun, pspec_ctf, dfx_upper, dfy_upper, angast_upper)
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
        call init_srch( ppspec_all, pspec_all_roavg )
        call refine( dfx, dfy, angast, phshift, cc )
        ! calculate error
        dfavg       = (dfx + dfy / 2)
        dfavg_lower = (dfx_lower + dfy_lower / 2)
        dfavg_upper = (dfx_upper + dfy_upper / 2)
        dferr       = maxval([abs(dfavg - dfavg_lower), abs(dfavg - dfavg_upper), abs(dfavg_lower - dfavg_upper)])
        ! make a half-n-half diagnostic
        999 if( l_phaseplate )then
            call ctf2pspecimg(tfun, pspec_ctf, dfx, dfy, angast, add_phshift=phshift)
        else
            call ctf2pspecimg(tfun, pspec_ctf, dfx, dfy, angast)
        endif
        call pspec_ctf%norm
        call ppspec_ref%norm
        call pspec_ctf%mul(imgmsk)
        pspec_half_n_half = ppspec_ref%before_after(pspec_ctf, cc_msk)
        call pspec_half_n_half%write(trim(diagfname), 1)
        call pspec_half_n_half%kill
        ! calculate CTF score diagnostic
        df_avg = dfx + dfy / 2.0
        call ctf2pspecimg(tfun, pspec_ctf_roavg, df_avg, df_avg, 0.)
        hpfind = pspec_all_roavg%get_find(hp)
        filtsz = pspec_lower_roavg%get_filtsz()
        call pspec_all_roavg%mask(real(filtsz), 'soft', real(hpfind))
        call pspec_all_roavg%norm
        call pspec_ctf_roavg%mask(real(filtsz), 'soft', real(hpfind))
        allocate(corrs(filtsz))
        call pspec_all_roavg%frc_pspec(pspec_ctf_roavg, corrs)
        call pspec_all_roavg%write('pspec_all_roavg.mrc', 1)
        call pspec_ctf_roavg%write('pspec_ctf_roavg.mrc', 1)
        corrs(1:hpfind) = 1.0
        ctfscore        = real(count(corrs > 0.)) / real(filtsz)
    end subroutine ctffit_x_validated_fit

    subroutine init_srch( which_pspec, which_pspec_roavg )
        class(image), target :: which_pspec, which_pspec_roavg
        ppspec_ref       => which_pspec
        ppspec_ref_roavg => which_pspec_roavg
        ! memoize reference corr components
        call ppspec_ref%prenorm4real_corr(sxx, cc_msk)
        call ppspec_ref_roavg%prenorm4real_corr(sxx_roavg, cc_msk)
    end subroutine init_srch

    ! CONSIDER TASK PARALLELISM HERE
    subroutine grid_srch( dfx, dfy, angast, phshift, cc )
        real, intent(out) :: dfx, dfy, angast, phshift, cc
        real              :: cost, df, cost_lowest, dfstep, df_best
        class(*), pointer :: fun_self => null()
        dfstep = (df_max - df_min) / real(NSTEPS)
        if( l_phaseplate )then
            ! do a first grid search assuming:
            ! no astigmatism
            ! pi half phase shift
            df      = df_min
            df_best = df
            cost_lowest = ctffit_cost_phaseplate_roavg(fun_self, [df,df,0.,PIO2], ndim)
            do while( df <= df_max )
                cost = ctffit_cost_phaseplate_roavg(fun_self, [df,df,0.,PIO2], ndim)
                if( cost < cost_lowest )then
                    cost_lowest = cost
                    df_best     = df
                endif
                df = df + dfstep
            end do
        else
            ! do a first grid search assuming:
            ! no astigmatism
            df = df_min
            df_best = df
            cost_lowest = ctffit_cost_roavg(fun_self, [df,df,0.], ndim)
            do while( df <= df_max )
                cost = ctffit_cost_roavg(fun_self, [df,df,0.], ndim)
                if( cost < cost_lowest )then
                    cost_lowest = cost
                    df_best     = df
                endif
                df = df + dfstep
            end do
        endif
        dfx     = df_best
        dfy     = df_best
        angast  = 0.
        phshift = 0.
        if( l_phaseplate ) phshift = PIO2
        cc      = -cost_lowest
    end subroutine grid_srch

    subroutine refine( dfx, dfy, angast, phshift, cc )
        real, intent(inout) :: dfx, dfy, angast, phshift, cc
        real                :: cost
        class(*), pointer   :: fun_self => null()
        if( l_phaseplate )then
            ospec_de%x = [dfx,dfy,angast,phshift]
        else
            ospec_de%x = [dfx,dfy,angast]
        endif
        call diffevol%minimize(ospec_de, fun_self, cost)
        dfx    = ospec_de%x(1)
        dfy    = ospec_de%x(2)
        angast = ospec_de%x(3)
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
        real ::  cost
        call ctf2pspecimg(tfun, pspec_ctf, vec(1), vec(2), vec(3))
        cost = -ppspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
    end function ctffit_cost

    ! this cost function considers only the rotationally averaged model
    function ctffit_cost_roavg( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        call ctf2pspecimg(tfun, pspec_ctf_roavg, vec(1), vec(1), 0.)
        cost = -ppspec_ref_roavg%real_corr_prenorm(pspec_ctf_roavg, sxx_roavg, cc_msk)
    end function ctffit_cost_roavg

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    function ctffit_cost_phaseplate( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        ! vec(4) is additional phase shift (in radians)
        call ctf2pspecimg(tfun, pspec_ctf, vec(1), vec(2), vec(3), add_phshift=vec(4))
        cost = -ppspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
    end function ctffit_cost_phaseplate

    ! this cost function considers only the rotationally averaged model
    function ctffit_cost_phaseplate_roavg( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        ! vec(4) is additional phase shift (in radians)
        call ctf2pspecimg(tfun, pspec_ctf_roavg, vec(1), vec(1), 0., add_phshift=vec(4))
        cost = -ppspec_ref_roavg%real_corr_prenorm(pspec_ctf_roavg, sxx_roavg, cc_msk)
    end function ctffit_cost_phaseplate_roavg

    !>  \brief  is for making a CTF power-spec image
    subroutine ctf2pspecimg( tfun, img, dfx, dfy, angast, add_phshift )
        use simple_image, only: image
        class(ctf),     intent(inout) :: tfun        !< CTF object
        class(image),   intent(inout) :: img         !< image (output)
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real, optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer :: lims(3,2),h,mh,k,mk,phys(3),ldim(3),inds(3)
        real    :: ang, tval, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3), res, wght
        ! initialize
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call tfun%init(dfx, dfy, angast)
        img      = 0.
        lims     = img%loop_lims(3)
        mh       = maxval(lims(1,:))
        mk       = maxval(lims(2,:))
        inds     = 1
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        !$omp parallel do collapse(2) default(shared) private(h,hinv,k,kinv,inds,spaFreqSq,ang,tval,phys) &
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                inds(1)   = min(max(1,h+mh+1),ldim(1))
                inds(2)   = min(max(1,k+mk+1),ldim(2))
                inds(3)   = 1
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = tfun%eval(spaFreqSq, ang, aadd_phshift)
                tval      = min(1.,max(tval * tval,0.001))
                tval      = sqrt(tval)
                call img%set(inds, tval)
            end do
        end do
        !$omp end parallel do
    end subroutine ctf2pspecimg

    subroutine ctffit_kill
        ppspec_all       => null()
        ppspec_lower     => null()
        ppspec_upper     => null()
        ppspec_ref       => null()
        ppspec_ref_roavg => null()
        call pspec_ctf%kill
        call pspec_ctf_roavg%kill
        call pspec_all_roavg%kill
        call pspec_lower_roavg%kill
        call pspec_upper_roavg%kill
        call imgmsk%kill
        call ospec_de%kill
        call diffevol%kill
        call imgmsk%kill
        if( allocated(cc_msk) ) deallocate(cc_msk)
    end subroutine ctffit_kill

end module simple_ctffit

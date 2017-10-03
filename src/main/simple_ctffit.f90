module simple_ctffit
#include "simple_lib.f08"
use simple_image,    only: image
use simple_ctf,      only: ctf
use simple_opt_spec, only: opt_spec
use simple_de_opt,   only: de_opt
implicit none

public :: ctffit_init, ctffit_srch
private

type(image)          :: pspec_ref
type(image)          :: pspec_ctf
type(ctf)            :: tfun
type(opt_spec)       :: ospec
type(de_opt)         :: diffevol
logical, allocatable :: cc_msk(:,:,:)
integer              :: ldim(3)  = [0,0,0]
real                 :: df_min   = 0.5
real                 :: df_max   = 5.0
real                 :: hp       = 30.0
real                 :: lp       = 5.0
real                 :: sxx      = 0. 

contains

  	subroutine ctffit_init( pspec, smpd, kV, Cs, amp_contr, dfrange, resrange )
        class(image),   intent(in) :: pspec       !< powerspectrum
        real,           intent(in) :: smpd        !< sampling distance
        real,           intent(in) :: kV          !< acceleration voltage
        real,           intent(in) :: Cs          !< constant
        real,           intent(in) :: amp_contr   !< amplitude contrast
        real, optional, intent(in) :: dfrange(2)  !< defocus range, [30.0,5.0] default
        real, optional, intent(in) :: resrange(2) !< resolution range, [30.0,5.0] default
        type(image) :: tmpimg
        real        :: limits(3,2)
        ! set constants
        if( present(dfrange) )then
            if( dfrange(1) < dfrange(2) )then
        		    df_min = dfrange(1)
        		    df_max = dfrange(2)
            else
                stop 'invalid defocuis range; simple_ctffit :: new'
            endif
        endif
        if( present(resrange) )then
          	if( resrange(1) > resrange(2) )then
          		  hp = resrange(1)
          		  lp = resrange(2)
          	else
          		  stop 'invalid resolution range; simple_ctffit :: new'
          	endif
        endif
        ! construct CTF object
        tfun = ctf(smpd, kV, Cs, amp_contr)
        ! prepare powerspectra
        pspec_ref = pspec
        ldim = pspec_ref%get_ldim()
        call pspec_ref%dampen_central_cross
        call pspec_ref%subtr_backgr(hp)
        call pspec_ctf%new(ldim, smpd)
        ! generate correlation mask
        call tmpimg%new(ldim, smpd)
        call tmpimg%resmsk(hp, lp)
        cc_msk = tmpimg%bin2logical()
        call tmpimg%kill
        ! memoize reference corr components
        call pspec_ref%prenorm4real_corr(sxx, cc_msk)
        ! contruct optimiser
        limits(1:2,1) = df_min
        limits(1:2,2) = df_max
        limits(3,1)   = 0.
        limits(3,2)   = twopi ! miminise in radians so that the df:s are roughly on the same scale
        call ospec%specify('de', 3, limits=limits)
        call ospec%set_costfun(ctffit_cost)
        call diffevol%new(ospec)
  	end subroutine ctffit_init

  	subroutine ctffit_srch( dfx, dfy, angast, cc )
    		real, intent(out) :: dfx, dfy, angast, cc
    		real :: cost
    		ospec%x = 0. ! automatic initialisation within the DE
        ! optimisation by DE (Differential Evolution)
    		call diffevol%minimize(ospec, cost)
    		dfx    = ospec%x(1)
    		dfy    = ospec%x(2)
    		angast = rad2deg(ospec%x(3))
    		cc     = -cost
  	end subroutine ctffit_srch

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
  	function ctffit_cost( vec, D ) result( cost )
    		integer, intent(in) :: D
    		real,    intent(in) :: vec(D)
    		real :: cost
    		call tfun%ctf2pspecimg(pspec_ctf, vec(1), vec(2), rad2deg(vec(3)))
    		cost = -pspec_ref%real_corr_prenorm(pspec_ctf, sxx, cc_msk)
  	end function ctffit_cost

end module simple_ctffit

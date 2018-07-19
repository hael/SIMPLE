! Cartesian volume-volume alignment based on band-pass limited cross-correlation
module simple_vol_srch
include 'simple_lib.f08'
use simple_image,       only: image
use simple_opt_spec,    only: opt_spec
use simple_opt_simplex, only: opt_simplex
use simple_projector,   only: projector
implicit none

public :: vol_srch_init, vol_srch_minimize
private

type(opt_spec)    :: ospec         !< optimizer specification object
type(opt_simplex) :: nlopt         !< optimizer object
integer           :: nrestarts = 3 !< simplex restarts (randomized bounds)
real              :: hp, lp, trs   !< srch ctrl params
type(image)       :: vref          !< reference volume
type(projector)   :: vtarg         !< target volume (subjected to rotation and shift)

contains

    subroutine vol_srch_init( vol_ref, vol_target, hp_in, lp_in, trs_in, angres, e_start, nrestarts_in )
        use simple_ori, only: ori
        class(image),      intent(inout) :: vol_ref, vol_target
        real,              intent(in)    :: hp_in, lp_in, trs_in, angres
        class(ori),        intent(in)    :: e_start
        integer, optional, intent(in)    :: nrestarts_in
        integer :: ldim(3), ldim_pd(3), boxpd
        real    :: smpd, lims(6,2), eul(3)
        ! set corr/srch ctrl params
        hp  = hp_in
        lp  = lp_in
        trs = trs_in
        ! set nrestarts
        nrestarts = 3
        if( present(nrestarts_in) ) nrestarts = nrestarts_in
        ! prepare vref and vtarg for corr calc
        ldim    = vol_ref%get_ldim()
        smpd    = vol_ref%get_smpd()
        boxpd   = 2 * round2even(KBALPHA * real(ldim(1) / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        call vref%new(ldim_pd, smpd)
        call vtarg%new(ldim_pd, smpd)
        call vol_ref%pad(vref)
        call vol_target%pad(vtarg)
        call vref%fft
        call vtarg%fft
        call vtarg%expand_cmat(KBALPHA)
        ! make optimizer spec
        eul         =  e_start%get_euler()
        lims(1:3,1) =  eul - angres
        lims(1:3,2) =  eul + angres
        lims(4:6,1) = -trs
        lims(4:6,2) =  trs
        call ospec%specify('simplex', 6, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=nrestarts, maxits=100)
        call ospec%set_costfun(vol_srch_costfun)
        ! generate the simplex optimizer object
        call nlopt%new(ospec)
    end subroutine vol_srch_init

    function vol_srch_costfun( fun_self, vec, D ) result( cost )
        use simple_ori, only: euler2m
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        integer :: sh,h,k,l,lims(3,2),logi(3),phys(3),sqlp,sqhp,sqarg
        real    :: loc(3),rmat(3,3),corr,cost,sumasq,sumbsq
        complex :: comp_targ, comp_ref
        lims   = vref%loop_lims(1,lp)
        rmat   = euler2m(vec(1:3))
        sqlp   = (maxval(lims(:,2)))**2
        sqhp   = max(2,vref%get_find(hp))**2
        corr   = 0.
        sumasq = 0.
        sumbsq = 0.
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sqarg,logi,sh,phys,loc,comp_targ,comp_ref)&
        !$omp reduction(+:corr,sumasq,sumbsq) schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sqarg = h*h + k*k + l*l
                    if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                        logi      = [h,k,l]
                        sh        = nint(hyp(real(h),real(k),real(l)))
                        phys      = vref%comp_addr_phys(logi)
                        loc       = matmul(real(logi), rmat)
                        comp_targ = vtarg%interp_fcomp(loc) * vtarg%oshift(loc, vec(4:6))
                        comp_ref  = vref%get_fcomp(logi, phys)
                        ! corr is real part of the complex mult btw ref and targ*
                        corr      = corr + real(comp_ref * conjg(comp_targ))
                        sumasq    = sumasq + csq(comp_targ)
                        sumbsq    = sumbsq + csq(comp_ref)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( sumasq > 0. .and. sumbsq > 0. )then
            corr = corr / sqrt(sumasq * sumbsq)
        else
            corr = 0.
        endif
        cost = -corr
    end function vol_srch_costfun

    function vol_srch_minimize( ) result( ceulxyz )
        real              :: cost_init, cost, ceulxyz(7)
        class(*), pointer :: fun_self => null()
        cost_init = vol_srch_costfun(fun_self, ospec%x, ospec%ndim)

        print *, 'cost_init: ', cost_init

        call nlopt%minimize(ospec, fun_self, cost)

        print *, 'cost: ', cost

        if( cost < cost_init )then
            ceulxyz(1)  = -cost    ! correlation
            ceulxyz(2:) =  ospec%x ! Eulers and shift
        else
            ceulxyz(1)  = -1.      ! to indicate that better solution wasn't found
            ceulxyz(2:) =  0.
        endif
    end function vol_srch_minimize

end module simple_vol_srch

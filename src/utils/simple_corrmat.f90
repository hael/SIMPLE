!@descr: for calculation of correlation matrices
module simple_corrmat
use simple_pftc_srch_api
use simple_defs
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_builder,          only: builder
implicit none

public :: calc_cartesian_corrmat, calc_inpl_invariant_cc_nomirr
private

interface calc_cartesian_corrmat
    module procedure calc_cartesian_corrmat_1
    module procedure calc_cartesian_corrmat_2
end interface calc_cartesian_corrmat

contains
    
    subroutine calc_cartesian_corrmat_1( imgs, corrmat, msk, lp )
        type(image),       intent(inout) :: imgs(:)
        real, allocatable, intent(out)   :: corrmat(:,:)
        real, optional,    intent(in)    :: msk, lp
        integer, allocatable :: pairs(:,:)
        real    :: rmsk
        integer :: iptcl, jptcl, ipair, cnt, nptcls, ntot
        nptcls = size(imgs)
        rmsk   = 0.9 * real(imgs(1)%get_box()) / 2.0
        if( present(msk) ) rmsk = msk
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nptcls,nptcls), source=1.0)
        ntot = (nptcls*(nptcls-1))/2
        allocate(pairs(ntot,2))
        cnt = 0
        do iptcl = 1, nptcls-1
            do jptcl = iptcl+1, nptcls
                cnt = cnt + 1
                pairs(cnt,1) = iptcl
                pairs(cnt,2) = jptcl
            end do
        end do
        if( present(lp) )then
            ! Fourier correlation
            ! apply a soft-edged mask & fwd FT
            call imgs(1)%memoize_mask_coords
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl=1,nptcls
                call imgs(iptcl)%mask2D_soft(rmsk)
                call imgs(iptcl)%fft()
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(ipair) schedule(static) proc_bind(close)
            do ipair = 1,ntot
                corrmat(pairs(ipair,1),pairs(ipair,2)) = imgs(pairs(ipair,1))%corr(imgs(pairs(ipair,2)), lp_dyn=lp)
            end do
            !$omp end parallel do
        else
            ! Real-space correlation
            !$omp parallel do default(shared) private(ipair) schedule(static) proc_bind(close)
            do ipair=1,ntot
                corrmat(pairs(ipair,1),pairs(ipair,2)) = imgs(pairs(ipair,1))%real_corr(imgs(pairs(ipair,2)))
            end do
            !$omp end parallel do
        endif
        ! symmetrize matrix
        do ipair = 1,ntot
            corrmat(pairs(ipair,2),pairs(ipair,1)) = corrmat(pairs(ipair,1),pairs(ipair,2))
        end do
        ! cleanup
        deallocate(pairs)
    end subroutine calc_cartesian_corrmat_1

    subroutine calc_cartesian_corrmat_2( imgs_sel, imgs_orig, corrmat, msk, lp )
        type(image),       intent(inout) :: imgs_sel(:), imgs_orig(:)
        real, allocatable, intent(out)   :: corrmat(:,:)
        real, optional,    intent(in)    :: msk, lp
        real    :: rmsk
        integer :: iptcl, isel, norig, nsel
        logical :: err
        ! set const
        norig = size(imgs_orig)
        nsel  = size(imgs_sel)
        rmsk  = 0.9 * real(imgs_sel(1)%get_box()) / 2.0
        if( present(msk) ) rmsk = msk
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nsel,norig))
        ! Correlation marix calculation
        if( present(lp) )then
            ! prep imgs for Fourier corrcalc
            ! apply a soft-edged mask & fwd FT
            call imgs_sel(1)%memoize_mask_coords
            do iptcl=1,nsel
                call imgs_sel(iptcl)%mask2D_soft(rmsk)
                call imgs_sel(iptcl)%fft()
            end do
            do iptcl=1,norig
                call imgs_orig(iptcl)%mask2D_soft(rmsk)
                call imgs_orig(iptcl)%fft()
            end do
            ! Fourier correlation
            do isel=1,nsel
                do iptcl=1,norig
                    corrmat(isel,iptcl) = imgs_sel(isel)%corr(imgs_orig(iptcl),lp_dyn=lp)
                end do
            end do
        else
            ! Real space correlations
            ! prenormalization to zero mean and unit variance
            !$omp parallel do default(shared) private(isel,err) schedule(static) proc_bind(close)
            do isel=1,nsel
                call imgs_sel(isel)%prenorm4real_corr(err)
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(iptcl,err) schedule(static) proc_bind(close)
            do iptcl=1,norig
                call imgs_orig(iptcl)%prenorm4real_corr(err)
            end do
            !$omp end parallel do
            ! compute correlations
            !$omp parallel do collapse(2) default(shared) private(isel,iptcl) schedule(static) proc_bind(close)
            do isel=1,nsel
                do iptcl=1,norig
                    corrmat(isel,iptcl) = imgs_sel(isel)%real_corr_prenorm(imgs_orig(iptcl))
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine calc_cartesian_corrmat_2

    function calc_inpl_invariant_cc_nomirr( params, hp, lp, trs, imgs ) result( ccmat )
        class(parameters), intent(in)    :: params
        real,              intent(in)    :: hp, lp, trs
        class(image),      intent(inout) :: imgs(:)
        integer,      parameter     :: MAXITS_SH = 60
        real,         allocatable   :: inpl_corrs(:)
        type(pftc_shsrch_grad)      :: grad_shsrch_obj(nthr_glob)
        type(builder)               :: build
        real,         allocatable   :: ccmat(:,:)
        complex,      allocatable   :: pft(:,:)
        integer :: ldim(3), box, kfromto(2), ithr, i, j, loc, nrots, irot, irot_seed, nimgs
        real    :: smpd, lims(2,2), lims_init(2,2), joint_lims(3,2), cxy(3)
        real(dp) :: rotind_frac
        nimgs      = size(imgs)
        ldim       = imgs(1)%get_ldim()
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! initialize
        call build%pftc%new(params, nimgs, [1,nimgs], kfromto)
        call imgs(1)%memoize4polarize(build%pftc%get_pdim_srch())
        ! in-plane search object objects for parallel execution
        lims(:,1)      = -trs
        lims(:,2)      =  trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        nrots = build%pftc%get_nrots()
        joint_lims(1:2,:) = lims
        joint_lims(3,:) = [1.-2., real(nrots)+2.]
        do ithr = 1, nthr_glob
            if( trim(params%inpl_cont) == 'yes' )then
                call grad_shsrch_obj(ithr)%new_joint(build, joint_lims, MAXITS_SH)
            else
                call grad_shsrch_obj(ithr)%new_legacy(build, lims, lims_init=lims_init, &
                    &maxits=MAXITS_SH)
            endif
        end do
        !$omp parallel do default(shared)  private(i,ithr,pft) proc_bind(close) schedule(static)
        do i = 1, nimgs
            pft = build%pftc%allocate_pft()
            call imgs(i)%fft()
            call imgs(i)%polarize(pft)
            call build%pftc%set_ref_pft(i, pft, iseven=.true.)
            call build%pftc%cp_even_ref2ptcl(i, i)
            call imgs(i)%ifft
            deallocate(pft)
        end do
        !$omp end parallel do
        call build%pftc%memoize_refs
        call build%pftc%memoize_ptcls
        ! register imgs
        allocate(inpl_corrs(nrots), ccmat(nimgs,nimgs))
        ccmat = 1. ! takes care of diagonal elements
        !$omp parallel do private(i,j,ithr,inpl_corrs,loc,irot,irot_seed,cxy,rotind_frac)&
        !$omp default(shared) schedule(dynamic) proc_bind(close)
        do i = 1, nimgs - 1
            do j = i + 1, nimgs
                ithr = omp_get_thread_num() + 1
                call build%pftc%gen_objfun_vals(j, i, [0.,0.], inpl_corrs)
                loc  = maxloc(inpl_corrs,dim=1)
                irot = loc
                call grad_shsrch_obj(ithr)%set_indices(j, i)
                if( trim(params%inpl_cont) == 'yes' )then
                    ! seeded local polish of the pre-selected cell
                    irot_seed = irot
                    cxy = grad_shsrch_obj(ithr)%minimize_joint(irot, [0.,0.], &
                        &sh_rot=.true., rotind_frac=rotind_frac, irot_in=irot_seed)
                else
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                endif
                if( irot == 0 )then ! no improved solution found, put back the old one
                    cxy(1) = inpl_corrs(loc)
                    cxy(2) = 0.
                    cxy(3) = 0.
                    irot   = loc
                endif
                ccmat(j,i) = cxy(1)
                ccmat(i,j) = ccmat(j,i)
            end do
        end do
        !$omp end parallel do
        ! destruct
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call build%pftc%kill
    end function calc_inpl_invariant_cc_nomirr

end module simple_corrmat

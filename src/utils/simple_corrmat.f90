!@descr: for calculation of correlation matrices
module simple_corrmat
use simple_pftc_srch_api
use simple_defs
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
implicit none

public :: calc_cartesian_corrmat, calc_inpl_invariant_cc_nomirr
private

interface calc_cartesian_corrmat
    module procedure calc_cartesian_corrmat_1
    module procedure calc_cartesian_corrmat_2
end interface calc_cartesian_corrmat

type(image)          :: mskimg
integer, allocatable :: pairs(:,:)
integer              :: nptcls, ntot, npix, norig, nsel
    
contains
    
    subroutine calc_cartesian_corrmat_1( imgs, corrmat, msk, lp )
        type(image),       intent(inout) :: imgs(:)
        real, allocatable, intent(out)   :: corrmat(:,:)
        real, optional,    intent(in)    :: msk, lp
        integer :: iptcl, jptcl, ipair, cnt
        nptcls = size(imgs)
        ! prep imgs for corrcalc
        call imgs(1)%memoize_mask_coords
        do iptcl=1,nptcls
            if( present(lp) )then
                ! if( .not. present(msk) ) THROW_HARD('need mask radius (msk) 4 Fourier corr calc!')
                ! apply a soft-edged mask
                call imgs(iptcl)%mask2D_soft(msk)
                ! Fourier transform
                call imgs(iptcl)%fft()
            endif
        end do
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nptcls,nptcls))
        corrmat = 1.
        ntot = (nptcls*(nptcls-1))/2
        if( present(lp) )then ! Fourier correlation
            cnt = 0
            do iptcl=1,nptcls-1
                do jptcl=iptcl+1,nptcls
                    cnt = cnt+1
                    call progress(cnt, ntot)
                    corrmat(iptcl,jptcl) = imgs(iptcl)%corr(imgs(jptcl), lp_dyn=lp)
                    corrmat(jptcl,iptcl) = corrmat(iptcl,jptcl)
                end do
            end do
        else ! Real-space correlation
            ! first make the pairs to balance the parallel section
            allocate(pairs(ntot,2))
            cnt = 0
            do iptcl=1,nptcls-1
                do jptcl=iptcl+1,nptcls
                    cnt = cnt+1
                    pairs(cnt,1) = iptcl
                    pairs(cnt,2) = jptcl
                end do
            end do
            !$omp parallel do default(shared) private(ipair) schedule(static) proc_bind(close)
            do ipair=1,ntot
                corrmat(pairs(ipair,1),pairs(ipair,2))=&
                imgs(pairs(ipair,1))%real_corr(imgs(pairs(ipair,2)))
                corrmat(pairs(ipair,2),pairs(ipair,1)) = corrmat(pairs(ipair,1),pairs(ipair,2))
            end do
            !$omp end parallel do
            deallocate(pairs)
            call mskimg%kill
        endif
    end subroutine calc_cartesian_corrmat_1

    subroutine calc_cartesian_corrmat_2( imgs_sel, imgs_orig, corrmat, msk, lp )
        type(image),       intent(inout) :: imgs_sel(:), imgs_orig(:)
        real, allocatable, intent(out)   :: corrmat(:,:)
        real, optional,    intent(in)    :: msk, lp
        integer :: iptcl, isel
        logical :: doftcalc, domsk
        ! set const
        norig    = size(imgs_orig)
        nsel     = size(imgs_sel)
        ! set operation modes
        domsk    = present(msk)
        doftcalc = present(lp)
        ! prep sel imgs for corrcalc
        call imgs_sel(1)%memoize_mask_coords
        do iptcl=1,nsel
            if( doftcalc )then
                ! if( .not. domsk ) THROW_HARD('need mask radius (msk) 4 Fourier corr calc!')
                ! apply a soft-edged mask
                call imgs_sel(iptcl)%mask2D_soft(msk)
                ! Fourier transform
                call imgs_sel(iptcl)%fft()
            endif
        end do
        ! prep orig imgs for corrcalc
        do iptcl=1,norig
            if( doftcalc )then
                ! apply a soft-edged mask
                call imgs_orig(iptcl)%mask2D_soft(msk)
                ! Fourier transform
                call imgs_orig(iptcl)%fft()
            endif
        end do
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nsel,norig))
        if( doftcalc )then ! Fourier correlation
            do isel=1,nsel
                call progress(isel,nsel)
                do iptcl=1,norig
                    corrmat(isel,iptcl) = imgs_sel(isel)%corr(imgs_orig(iptcl),lp_dyn=lp)
                end do
            end do
        else ! Real-space correlation
            if( domsk)then
                call mskimg%disc(imgs_sel(1)%get_ldim(), imgs_sel(1)%get_smpd(), msk, npix)
            endif
            !$omp parallel do collapse(2) default(shared) private(isel,iptcl) schedule(static) proc_bind(close)
            do isel=1,nsel
                do iptcl=1,norig
                    corrmat(isel,iptcl) = imgs_sel(isel)%real_corr(imgs_orig(iptcl))
                end do
            end do
            !$omp end parallel do
            call mskimg%kill
        endif
    end subroutine calc_cartesian_corrmat_2

    function calc_inpl_invariant_cc_nomirr( hp, lp, trs, imgs ) result( ccmat )
        real,         intent(in)    :: hp, lp, trs
        class(image), intent(inout) :: imgs(:)
        integer,      parameter     :: MAXITS_SH = 60
        real,         allocatable   :: inpl_corrs(:)
        type(pftc_shsrch_grad)      :: grad_shsrch_obj(nthr_glob)
        type(polarft_calc)          :: pftc
        real,         allocatable   :: ccmat(:,:)
        complex,      allocatable   :: pft(:,:)
        integer :: ldim(3), box, kfromto(2), ithr, i, j, k, loc, nrots, irot, nimgs
        real    :: smpd, lims(2,2), lims_init(2,2), cxy(3)
        nimgs      = size(imgs)
        ldim       = imgs(1)%get_ldim()
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! initialize
        call pftc%new(nimgs, [1,nimgs], kfromto)
        call imgs(1)%memoize4polarize(pftc%get_pdim(), KBALPHA)
        ! in-plane search object objects for parallel execution
        lims(:,1)      = -trs
        lims(:,2)      =  trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier='yes',&
            &maxits=MAXITS_SH, opt_angle=.true.)
        end do
        pft = pftc%allocate_pft()
        !$omp parallel do default(shared)  private(i,pft) proc_bind(close) schedule(static)
        do i = 1, nimgs
            call imgs(i)%fft()
            call imgs(i)%polarize(pft)
            call pftc%set_ref_pft(i, pft, iseven=.true.)
            call pftc%cp_even_ref2ptcl(i, i)
            call imgs(i)%ifft
        end do
        !$omp end parallel do
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        ! register imgs
        nrots = pftc%get_nrots()
        allocate(inpl_corrs(nrots), ccmat(nimgs,nimgs))
        ccmat = 1. ! takes care of diagonal elements
        !$omp parallel do private(i,j,ithr,inpl_corrs,loc,irot,cxy)&
        !$omp default(shared) schedule(dynamic) proc_bind(close)
        do i = 1, nimgs - 1
            do j = i + 1, nimgs
                ithr = omp_get_thread_num() + 1
                call pftc%gen_objfun_vals(j, i, [0.,0.], inpl_corrs)
                loc  = maxloc(inpl_corrs,dim=1)
                irot = loc
                call grad_shsrch_obj(ithr)%set_indices(j, i)
                cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
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
        call pftc%kill
        if( allocated(pft) ) deallocate(pft)
    end function calc_inpl_invariant_cc_nomirr

end module simple_corrmat

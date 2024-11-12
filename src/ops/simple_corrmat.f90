! for calculation of correlation matrices
module simple_corrmat
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
implicit none

public :: calc_cartesian_corrmat, calc_inplane_invariant_corrmat
private
#include "simple_local_flags.inc"

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
        do iptcl=1,nptcls
            if( present(lp) )then
                if( .not. present(msk) ) THROW_HARD('need mask radius (msk) 4 Fourier corr calc!')
                ! apply a soft-edged mask
                call imgs(iptcl)%mask(msk, 'soft')
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
        do iptcl=1,nsel
            if( doftcalc )then
                if( .not. domsk ) THROW_HARD('need mask radius (msk) 4 Fourier corr calc!')
                ! apply a soft-edged mask
                call imgs_sel(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs_sel(iptcl)%fft()
            endif
        end do
        ! prep orig imgs for corrcalc
        do iptcl=1,norig
            if( doftcalc )then
                ! apply a soft-edged mask
                call imgs_orig(iptcl)%mask(msk, 'soft')
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

    subroutine calc_inplane_invariant_corrmat( imgs, hp, lp, corrmat)
        use simple_polarizer,         only: polarizer
        use simple_polarft_corrcalc,  only: polarft_corrcalc
        use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
        class(image),         intent(inout)  :: imgs(:)
        real,                 intent(in)     :: hp, lp
        real, allocatable,    intent(out)    :: corrmat(:,:)
        type(pftcc_shsrch_grad), allocatable :: grad_shsrch_obj(:)
        type(polarizer)                      :: polartransform
        type(polarft_corrcalc)               :: pftcc
        real, allocatable :: inpl_corrs(:)
        integer :: n, i, j, ithr, nrots, loc(1), irot
        real    :: lims(2,2), lims_init(2,2), cxy(3)
        n = size(imgs)
        ! resolution limits
        params_glob%kfromto(1) = max(2, calc_fourier_index(hp, params_glob%box, params_glob%smpd))
        params_glob%kfromto(2) =        calc_fourier_index(lp, params_glob%box, params_glob%smpd)
        ! initialize pftcc, polarizer
        call pftcc%new(n, [1,n], params_glob%kfromto)
        call polartransform%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call polartransform%init_polarizer(pftcc, params_glob%alpha)
        ! in-plane search object
        lims(:,1)      = -params_glob%trs
        lims(:,2)      =  params_glob%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        allocate(grad_shsrch_obj(nthr_glob))
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
            &maxits=params_glob%maxits_sh, opt_angle=.true.)
        end do
        
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, n
            call imgs(i)%fft()
            call polartransform%polarize(pftcc, imgs(i), i, isptcl=.false., iseven=.true.)
            call pftcc%cp_even_ref2ptcl(i, i)
            call imgs(i)%ifft()
        end do
        !$omp end parallel do
        call pftcc%memoize_refs
        call pftcc%memoize_ptcls       
        ! generate nearest neighbor structure
        nrots = pftcc%get_nrots()
        allocate(inpl_corrs(nrots), corrmat(n,n), source=0.)
        corrmat = 1.
        !$omp parallel do default(shared) private(i,j,ithr,inpl_corrs,loc,irot,cxy) schedule(dynamic) proc_bind(close)
        do i = 1, n - 1
            do j = i + 1, n
                ithr = omp_get_thread_num() + 1
                call pftcc%gencorrs(i, j, inpl_corrs)
                loc  = maxloc(inpl_corrs)
                irot = loc(1) 
                call grad_shsrch_obj(ithr)%set_indices(i, j)
                cxy = grad_shsrch_obj(ithr)%minimize(irot=irot)
                if( irot > 0 )then
                    corrmat(i,j)= cxy(1)
                else
                    corrmat(i,j) = inpl_corrs(loc(1))
                endif
                corrmat(j,i) = corrmat(i,j)
            enddo
        enddo
        !$omp end parallel do
        ! destruct
        call polartransform%kill
        call pftcc%kill
    end subroutine calc_inplane_invariant_corrmat

end module simple_corrmat

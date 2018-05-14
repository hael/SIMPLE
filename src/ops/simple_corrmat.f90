! for calculation of cartesian correlation matrices

module simple_corrmat
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_error, only: allocchk
use simple_jiffys, only: progress
use simple_image,  only: image

implicit none

public :: calc_cartesian_corrmat!, calc_roinv_corrmat
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
        do iptcl=1,nptcls
            if( present(lp) )then
                if( .not. present(msk) ) stop 'need mask radius (msk) 4 Fourier corr calc!'
                ! apply a soft-edged mask
                call imgs(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs(iptcl)%fft()
            endif
        end do
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nptcls,nptcls), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: calc_cartesian_corrmat_1; simple_corrmat, 1',alloc_stat)
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
            allocate(pairs(ntot,2), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: calc_cartesian_corrmat_1; simple_corrmat, 2',alloc_stat)
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
                if( .not. domsk ) stop 'need mask radius (msk) 4 Fourier corr calc!'
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
        allocate(corrmat(nsel,norig), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: calc_cartesian_corrmat_2; simple_corrmat, 1',alloc_stat)
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

    ! subroutine calc_roinv_corrmat( pftcc, corrmat )
    !     ! assumes particles/references in pftcc identical sets
    !     use simple_polarft_corrcalc, only: polarft_corrcalc
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     real, allocatable,       intent(out)   :: corrmat(:,:)
    !     real    :: corrs(pftcc%get_nrots())
    !     integer :: iref, iptcl, nptcls, nrefs, loc(1)
    !     nptcls = pftcc%get_nptcls()
    !     nrefs  = pftcc%get_nrefs()
    !     if( nptcls /= nrefs ) stop 'nptcls == nrefs in pftcc required; simple_corrmat :: calc_roinv_corrmat'
    !     if( allocated(corrmat) ) deallocate(corrmat)
    !     allocate(corrmat(nptcls,nptcls), stat=alloc_stat)
    !     if(alloc_stat/=0)call allocchk('In: calc_roinv_corrmat; simple_corrmat')
    !     corrmat = 1.
    !     !$omp parallel do default(shared) schedule(guided) private(iref,iptcl,corrs,loc) proc_bind(close)
    !     do iref=1,nptcls - 1
    !         do iptcl=iref + 1,nptcls
    !             ! rotational corr
    !             call pftcc%gencorrs(iref, iptcl, corrs)
    !             loc  = maxloc(corrs)
    !             corrmat(iref,iptcl) = corrs(loc(1))
    !             ! symmetrize
    !             corrmat(iptcl,iref) = corrmat(iref,iptcl)
    !         end do
    !     end do
    !     !$omp end parallel do
    ! end subroutine calc_roinv_corrmat

end module simple_corrmat

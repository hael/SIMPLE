module simple_corrmat
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,  only: image
use simple_jiffys, only: alloc_err, progress
implicit none

public :: calc_cartesian_corrmat
private

interface calc_cartesian_corrmat
    module procedure calc_cartesian_corrmat_1
    module procedure calc_cartesian_corrmat_2
end interface calc_cartesian_corrmat

type(image)          :: mskimg
integer, allocatable :: pairs(:,:)
integer              :: nptcls, ntot, npix, norig, nsel
logical              :: debug=.false.
    
contains
    
    subroutine calc_cartesian_corrmat_1( imgs, corrmat, msk, lp )
        type(image),       intent(inout) :: imgs(:)
        real, allocatable, intent(out)   :: corrmat(:,:)
        real, optional,    intent(in)    :: msk, lp 
        integer :: iptcl, jptcl, ipair, alloc_stat, cnt
        nptcls = size(imgs)
        ! prep imgs for corrcalc
        do iptcl=1,nptcls
            if( present(lp) )then
                if( .not. present(msk) ) stop 'need mask radius (msk) 4 Fourier corr calc!'
                ! apply a soft-edged mask
                call imgs(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs(iptcl)%fwd_ft
            endif
        end do
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nptcls,nptcls), stat=alloc_stat)
        call alloc_err('In: calc_cartesian_corrmat_1; simple_corrmat, 1', alloc_stat)
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
            call alloc_err('In: calc_cartesian_corrmat_1; simple_corrmat, 2', alloc_stat)
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
        integer :: iptcl, isel, alloc_stat, cnt
        norig = size(imgs_orig)
        nsel  = size(imgs_sel)
        ! prep sel imgs for corrcalc
        do iptcl=1,nsel
            if( present(lp) )then
                if( .not. present(msk) ) stop 'need mask radius (msk) 4 Fourier corr calc!'
                ! apply a soft-edged mask
                call imgs_sel(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs_sel(iptcl)%fwd_ft
            endif
        end do
        ! prep orig imgs for corrcalc
        do iptcl=1,norig
            if( present(lp) )then
                ! apply a soft-edged mask
                call imgs_orig(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs_orig(iptcl)%fwd_ft
            endif
        end do
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nsel,norig), stat=alloc_stat)
        call alloc_err('In: calc_cartesian_corrmat_2; simple_corrmat, 1', alloc_stat)
        if( present(lp) )then ! Fourier correlation
            do isel=1,nsel
                call progress(isel,nsel)
                do iptcl=1,norig
                    corrmat(isel,iptcl) = imgs_sel(isel)%corr(imgs_orig(iptcl),lp_dyn=lp)
                end do
            end do
        else ! Real-space correlation
            if( present(msk) ) call mskimg%disc(imgs_sel(1)%get_ldim(), imgs_sel(1)%get_smpd(), msk, npix)
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

end module simple_corrmat

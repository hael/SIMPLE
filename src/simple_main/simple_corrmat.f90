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
    
    subroutine calc_cartesian_corrmat_2( imgs_sel, imgs_orig, corrmat, msk, lp, scale )
        use simple_magic_boxes, only: autoscale
        type(image),       intent(inout) :: imgs_sel(:), imgs_orig(:)
        real, allocatable, intent(out)   :: corrmat(:,:)
        real, optional,    intent(in)    :: msk, lp, scale
        integer :: iptcl, isel, alloc_stat, cnt, ldim(3), box_sc
        real    :: smpd, smpd_target, smpd_sc, sscale
        logical :: doscale, doftcalc, domsk
        ! set const
        norig    = size(imgs_orig)
        nsel     = size(imgs_sel)
        ldim     = imgs_sel(1)%get_ldim()
        smpd     = imgs_sel(1)%get_smpd()
        ! set operation modes
        domsk    = present(msk)
        doftcalc = present(lp)
        doscale  = present(scale) 
        if( doscale )then
            smpd_target = smpd/scale
            call autoscale(ldim(1), smpd, smpd_target, box_sc, smpd_sc, sscale)
        endif
        ! prep sel imgs for corrcalc
        do iptcl=1,nsel
            if( doftcalc )then
                if( .not. domsk ) stop 'need mask radius (msk) 4 Fourier corr calc!'
                ! apply a soft-edged mask
                call imgs_sel(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs_sel(iptcl)%fwd_ft
                ! scale
                call imgs_sel(iptcl)%clip_inplace([box_sc,box_sc,1])
            else
                ! scale
                call imgs_sel(iptcl)%fwd_ft
                call imgs_sel(iptcl)%clip_inplace([box_sc,box_sc,1])
                call imgs_sel(iptcl)%bwd_ft
            endif
        end do
        ! prep orig imgs for corrcalc
        do iptcl=1,norig
            if( doftcalc )then
                ! apply a soft-edged mask
                call imgs_orig(iptcl)%mask(msk, 'soft')
                ! Fourier transform
                call imgs_orig(iptcl)%fwd_ft
                ! scale
                call imgs_orig(iptcl)%clip_inplace([box_sc,box_sc,1])
            else
                ! scale
                call imgs_orig(iptcl)%fwd_ft
                call imgs_orig(iptcl)%clip_inplace([box_sc,box_sc,1])
                call imgs_orig(iptcl)%bwd_ft
            endif
        end do
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nsel,norig), stat=alloc_stat)
        call alloc_err('In: calc_cartesian_corrmat_2; simple_corrmat, 1', alloc_stat)
        if( doftcalc )then ! Fourier correlation
            do isel=1,nsel
                call progress(isel,nsel)
                do iptcl=1,norig
                    corrmat(isel,iptcl) = imgs_sel(isel)%corr(imgs_orig(iptcl),lp_dyn=lp)
                end do
            end do
        else ! Real-space correlation
            if( domsk)then
                if( doscale )then
                    call mskimg%disc(imgs_sel(1)%get_ldim(), imgs_sel(1)%get_smpd(), sscale*msk, npix)
                else
                    call mskimg%disc(imgs_sel(1)%get_ldim(), imgs_sel(1)%get_smpd(), msk, npix)
                endif
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

end module simple_corrmat

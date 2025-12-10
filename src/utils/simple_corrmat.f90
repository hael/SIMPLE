! for calculation of correlation matrices
module simple_corrmat
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: calc_cartesian_corrmat, calc_inpl_invariant_fm
private

interface calc_cartesian_corrmat
    module procedure calc_cartesian_corrmat_1
    module procedure calc_cartesian_corrmat_2
end interface calc_cartesian_corrmat

interface calc_inpl_invariant_fm
    module procedure calc_inpl_invariant_fm_1
    module procedure calc_inpl_invariant_fm_2
end interface calc_inpl_invariant_fm

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
                ! if( .not. present(msk) ) THROW_HARD('need mask radius (msk) 4 Fourier corr calc!')
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
                ! if( .not. domsk ) THROW_HARD('need mask radius (msk) 4 Fourier corr calc!')
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

    subroutine calc_inpl_invariant_fm_1( imgs, hp, lp, trs, corrmat, l_srch_mirr )
        use simple_pftc_shsrch_fm
        use simple_polarizer,        only: polarizer
        use simple_polarft_calc, only: polarft_calc
        class(image),          intent(inout) :: imgs(:)
        real,                  intent(in)    :: hp, lp, trs
        real,    allocatable,  intent(inout) :: corrmat(:,:)
        type(image),           allocatable   :: ccimgs(:,:)
        type(pftc_shsrch_fm), allocatable   :: fm_correlators(:)
        logical, optional,     intent(in)    :: l_srch_mirr
        type(polarizer)        :: polartransform
        type(polarft_calc) :: pftc
        real, parameter :: TRS_STEPSZ = 1.0
        integer :: n, i, j, ithr, ldim(3), box, kfromto(2)
        real    :: offset(2), offsetm(2), ang, angm, smpd, cc, ccm
        logical :: ll_srch_mirr
        ll_srch_mirr = .true.
        if( present(l_srch_mirr) ) ll_srch_mirr = l_srch_mirr
        n          = size(imgs)
        ldim       = imgs(1)%get_ldim()
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! initialize pftc, polarizer
        call pftc%new(n, [1,n], kfromto)
        call polartransform%new([box,box,1], smpd)
        call polartransform%init_polarizer(pftc, KBALPHA)
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(n,n), source=-1.)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, n
            call imgs(i)%fft()
            call polartransform%polarize(pftc, imgs(i), i, isptcl=.false., iseven=.true.)
            call pftc%cp_even_ref2ptcl(i, i)
            call pftc%mirror_ref_pft(i) ! controlled by even/odd flag (see below)
        end do
        !$omp end parallel do
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        ! correlation matrix calculation
        allocate(fm_correlators(nthr_glob),ccimgs(nthr_glob,2))
        do i = 1,nthr_glob
            call fm_correlators(i)%new(trs,TRS_STEPSZ,opt_angle=.false.)
            call ccimgs(i,1)%new(ldim, smpd, wthreads=.false.)
            call ccimgs(i,2)%new(ldim, smpd, wthreads=.false.)
        enddo
        if( ll_srch_mirr )then
            !$omp parallel do default(shared) private(i,j,cc,ccm,ithr,offset,offsetm,ang,angm)&
            !$omp schedule(dynamic) proc_bind(close)
            do i = 1, n - 1
                ithr = omp_get_thread_num()+1
                corrmat(i,i) = 1.
                do j = i + 1, n
                    ! reference to particle
                    call pftc%set_eo(i,.true.)
                    call fm_correlators(ithr)%calc_phasecorr(j, i, imgs(j), imgs(i),&
                        &ccimgs(ithr,1), ccimgs(ithr,2), cc, rotang=ang, shift=offset)
                    ! mirrored reference to particle
                    call pftc%set_eo(i,.false.) ! switch mirror
                    call fm_correlators(ithr)%calc_phasecorr(j, i, imgs(j), imgs(i),&
                    &ccimgs(ithr,1), ccimgs(ithr,2), ccm, mirror=.true., rotang=angm, shift=offsetm)
                    ! higher correlation wins
                    if( ccm > cc )then
                        cc     = ccm
                        ang    = angm
                        offset = offsetm
                    endif
                    corrmat(i,j) = cc
                    corrmat(j,i) = cc
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,j,cc,ithr,offset,ang)&
            !$omp schedule(dynamic) proc_bind(close)
            do i = 1, n - 1
                ithr = omp_get_thread_num()+1
                corrmat(i,i) = 1.
                do j = i + 1, n
                    ! reference to particle
                    call pftc%set_eo(i,.true.)
                    call fm_correlators(ithr)%calc_phasecorr(j, i, imgs(j), imgs(i),&
                        &ccimgs(ithr,1), ccimgs(ithr,2), cc, rotang=ang, shift=offset)
                    corrmat(i,j) = cc
                    corrmat(j,i) = cc
                enddo
            enddo
            !$omp end parallel do
        endif
        corrmat(n,n) = 1. ! leftover
        ! BWD FT
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, n
            call imgs(i)%ifft()
        end do
        !$omp end parallel do
        ! tidy
        call pftc%kill
        call polartransform%kill_polarizer
        call polartransform%kill
        do i = 1,nthr_glob
            call fm_correlators(i)%kill
            call ccimgs(i,1)%kill
            call ccimgs(i,2)%kill
        enddo
    end subroutine calc_inpl_invariant_fm_1

    subroutine calc_inpl_invariant_fm_2( refimgs, imgs, hp, lp, trs, corrmat )
        use simple_pftc_shsrch_fm
        use simple_polarizer,        only: polarizer
        use simple_polarft_calc, only: polarft_calc
        class(image),          intent(inout) :: refimgs(:), imgs(:)
        real,                  intent(in)    :: hp, lp, trs
        real,    allocatable,  intent(inout) :: corrmat(:,:)
        type(image),           allocatable   :: ccimgs(:,:)
        type(pftc_shsrch_fm), allocatable   :: fm_correlators(:)
        type(polarizer)        :: polartransform
        type(polarft_calc) :: pftc
        real, parameter :: TRS_STEPSZ = 1.0
        integer :: n, i, ithr, ldim(3), box, kfromto(2), iref, nrefs
        real    :: offset(2), offsetm(2), ang, angm, smpd, cc, ccm
        n          = size(imgs)
        nrefs      = size(refimgs)
        ldim       = imgs(1)%get_ldim()
        box        = ldim(1)
        smpd       = imgs(1)%get_smpd()
        kfromto(1) = max(2, calc_fourier_index(hp, box, smpd))
        kfromto(2) =        calc_fourier_index(lp, box, smpd)
        ! initialize pftc, polarizer
        call pftc%new(nrefs, [1,n], kfromto)
        call polartransform%new([box,box,1], smpd)
        call polartransform%init_polarizer(pftc, KBALPHA)
        if( allocated(corrmat) ) deallocate(corrmat)
        allocate(corrmat(nrefs,n), source=-1.)
        !$omp parallel default(shared) private(i,iref) proc_bind(close)
        !$omp do schedule(static) 
        do iref = 1, nrefs
            call refimgs(iref)%fft()
            call polartransform%polarize(pftc, refimgs(iref), iref, isptcl=.false., iseven=.true.)
        end do
        !$omp end do
        !$omp do schedule(static) 
        do i = 1, n
            call imgs(i)%fft()
            call polartransform%polarize(pftc, imgs(i), i, isptcl=.true., iseven=.true.)
        end do
        !$omp end do
        !$omp end parallel
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        ! correlation matrix calculation
        allocate(fm_correlators(nthr_glob),ccimgs(nthr_glob,2))
        do i = 1,nthr_glob
            call fm_correlators(i)%new(trs,TRS_STEPSZ,opt_angle=.false.)
            call ccimgs(i,1)%new(ldim, smpd, wthreads=.false.)
            call ccimgs(i,2)%new(ldim, smpd, wthreads=.false.)
        enddo
        do iref = 1, nrefs
            !$omp parallel do default(shared) private(i,cc,ccm,ithr,offset,offsetm,ang,angm)&
            !$omp schedule(static) proc_bind(close)
            do i = 1, n
                ithr = omp_get_thread_num() + 1
                ! reference to particle
                call pftc%set_eo(i,.true.)
                call fm_correlators(ithr)%calc_phasecorr(iref, i, refimgs(iref), imgs(i),&
                    &ccimgs(ithr,1), ccimgs(ithr,2), cc, rotang=ang, shift=offset)
                ! mirrored reference to particle
                call pftc%mirror_ref_pft(iref) ! switch mirror
                call fm_correlators(ithr)%calc_phasecorr(iref, i, refimgs(iref), imgs(i),&
                &ccimgs(ithr,1), ccimgs(ithr,2), ccm, mirror=.true., rotang=angm, shift=offsetm)
                ! higher correlation wins
                if( ccm > cc )then
                    cc     = ccm
                    ang    = angm
                    offset = offsetm
                endif
                corrmat(iref,i) = cc
                call pftc%mirror_ref_pft(iref) ! switch mirror
            enddo
            !$omp end parallel do
        enddo
        ! BWD FT
        !$omp parallel default(shared) private(i,iref) proc_bind(close)
        !$omp do schedule(static) 
        do iref = 1, nrefs
            call refimgs(iref)%ifft()
        end do
        !$omp end do
        !$omp do schedule(static) 
        do i = 1, n
            call imgs(i)%ifft()
        end do
        !$omp end do
        !$omp end parallel
        ! tidy
        call pftc%kill
        call polartransform%kill_polarizer
        call polartransform%kill
        do i = 1,nthr_glob
            call fm_correlators(i)%kill
            call ccimgs(i,1)%kill
            call ccimgs(i,2)%kill
        enddo
    end subroutine calc_inpl_invariant_fm_2

end module simple_corrmat
    
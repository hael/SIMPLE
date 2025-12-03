! optimization(search)-based filter (uniform/nonuniform)
module simple_opt_filter
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image, image_ptr
use simple_parameters, only: params_glob
use simple_butterworth
implicit none
#include "simple_local_flags.inc"

public :: nonuni_filt3D, estimate_lplim, estimate_lplims2D
private

interface estimate_lplim
    module procedure estimate_lplim_1, estimate_lplim_2
end interface estimate_lplim


contains

    ! 3D optimization(search)-based nonuniform filter, paralellized version
    subroutine nonuni_filt3D(odd, even, mskimg, lpstop)
        class(image),           intent(inout) :: odd, even
        class(image), optional, intent(inout) :: mskimg
        real,         optional, intent(in)    :: lpstop
        type(image)          :: odd_copy_rmat, odd_copy_cmat, even_copy_rmat, even_copy_cmat, weights_img,&
                               &diff_img_opt_odd, diff_img_opt_even, diff_img_odd, diff_img_even, odd_filt, even_filt
        integer              :: k,l,m, box, ldim(3), find_start, find_stop, iter_no
        integer              :: filtsz, cutoff_find, lb(3), ub(3), smooth_ext
        real                 :: rad, find_stepsz, val, smpd
        type(image_ptr)      :: pdiff_odd, pdiff_even, pdiff_opt_odd, pdiff_opt_even, pweights
        integer, parameter   :: CHUNKSZ = 20
        real,    pointer     :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_odd_filt(:,:,:), rmat_even_filt(:,:,:)
        real,    allocatable :: fsc(:), cur_fil(:)
        type(string) :: fsc_fname
        ldim       = odd%get_ldim()
        filtsz     = odd%get_filtsz()
        smooth_ext = params_glob%smooth_ext
        box        = ldim(1)
        fsc_fname  = params_glob%fsc
        smpd       = even%get_smpd()
        find_start = calc_fourier_index(params_glob%lpstart_nonuni, box, even%get_smpd())
        find_stop  = find_start
        if( present(lpstop) )then
            find_stop = calc_fourier_index(lpstop, box, smpd)
        else
            ! calculate Fourier index low-pass limit for search based on FSC
            if( .not.file_exists(fsc_fname) ) then
                THROW_HARD('FSC file: '//fsc_fname%to_char()//' not found')
            else
                ! retrieve FSC and calculate low-pass limit
                fsc       = file2rarr(fsc_fname)
                find_stop = min(get_find_at_crit(fsc, 0.1),calc_fourier_index(params_glob%lpstop, box, smpd)) ! little overshoot, filter function anyway applied in polarft_calc
            endif
        endif
        find_stepsz = real(find_stop - find_start)/(params_glob%nsearch - 1)
        if( find_start >= find_stop ) THROW_HARD('nonuni_filt3D: starting Fourier index is larger than ending Fourier index!')
        call       weights_img%new(ldim, smpd)
        call      diff_img_odd%new(ldim, smpd)
        call     diff_img_even%new(ldim, smpd)
        call  diff_img_opt_odd%new(ldim, smpd)
        call diff_img_opt_even%new(ldim, smpd)
        call       weights_img%get_mat_ptrs(pweights)
        call      diff_img_odd%get_mat_ptrs(pdiff_odd)
        call     diff_img_even%get_mat_ptrs(pdiff_even)
        call  diff_img_opt_odd%get_mat_ptrs(pdiff_opt_odd)
        call diff_img_opt_even%get_mat_ptrs(pdiff_opt_even)
        call  odd_copy_rmat%copy(odd)
        call even_copy_rmat%copy(even)
        call  odd_copy_cmat%copy(odd)
        call even_copy_cmat%copy(even)
        call       odd_filt%copy(odd)
        call      even_filt%copy(even)
        call even_copy_cmat%fft
        call odd_copy_cmat%fft
        allocate(cur_fil(box), source=0.)
        if( present(mskimg) )then
            call bounds_from_mask3D(mskimg%bin2logical(), lb, ub)
        else
            lb = (/ 1, 1, 1/)
            ub = (/ box, box, box /)
        endif
        do k = 1, 3
            if( lb(k) < smooth_ext + 1 )   lb(k) = smooth_ext+1
            if( ub(k) > box - smooth_ext ) ub(k) = box - smooth_ext
        enddo
        call weights_img%zero_and_unflag_ft()
        do k = -smooth_ext, smooth_ext
            do l = -smooth_ext, smooth_ext
                do m = -smooth_ext, smooth_ext
                    rad = hyp(k,l,m)
                    val = -rad/(smooth_ext + 1) + 1.
                    if( val > 0 ) call weights_img%set_rmat_at(box/2+k+1, box/2+l+1, box/2+m+1, val)
                enddo
            enddo
        enddo
        call weights_img%fft()
        ! searching for the best fourier index from here and generating the optimized filter
        pdiff_opt_odd%rmat  = huge(val)
        pdiff_opt_even%rmat = huge(val)
        call  odd%get_rmat_ptr(rmat_odd)
        call even%get_rmat_ptr(rmat_even)
        rmat_odd  = 0.
        rmat_even = 0.
        call  odd_filt%get_rmat_ptr( rmat_odd_filt)
        call even_filt%get_rmat_ptr(rmat_even_filt)
        do iter_no = 1, params_glob%nsearch
            cutoff_find = nint(find_start + (iter_no - 1)*find_stepsz)
            ! filtering odd/even
            call  odd_filt%copy_fast( odd_copy_cmat)
            call even_filt%copy_fast(even_copy_cmat)
            call butterworth_filter( odd_filt, cutoff_find, cur_fil)
            call butterworth_filter(even_filt, cutoff_find, cur_fil)
            call even_filt%ifft
            call odd_filt%ifft
            call  odd_filt%sqeuclid_matrix(even_copy_rmat, diff_img_odd)
            call even_filt%sqeuclid_matrix( odd_copy_rmat, diff_img_even)
            ! do the non-uniform, i.e. optimizing at each voxel
            call diff_img_even%fft
            call diff_img_odd%fft
            !$omp parallel workshare
            pdiff_odd% cmat = pdiff_odd %cmat * pweights%cmat
            pdiff_even%cmat = pdiff_even%cmat * pweights%cmat
            !$omp end parallel workshare
            call diff_img_even%ifft
            call diff_img_odd%ifft
            !$omp parallel do collapse(3) default(shared) private(k,l,m) schedule(dynamic,CHUNKSZ) proc_bind(close)
            do m = lb(3),ub(3)
                do l = lb(2),ub(2)
                    do k = lb(1),ub(1)
                        if( pdiff_odd%rmat(k,l,m) < pdiff_opt_odd%rmat(k,l,m) )then
                            rmat_odd(          k,l,m) = rmat_odd_filt( k,l,m)
                            pdiff_opt_odd%rmat(k,l,m) = pdiff_odd%rmat(k,l,m)
                        endif
                        if( pdiff_even%rmat(k,l,m) < pdiff_opt_even%rmat(k,l,m) )then
                            rmat_even(          k,l,m) = rmat_even_filt( k,l,m)
                            pdiff_opt_even%rmat(k,l,m) = pdiff_even%rmat(k,l,m)
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
        enddo
        deallocate(cur_fil)
        call odd_copy_rmat%kill
        call odd_copy_cmat%kill
        call even_copy_rmat%kill
        call even_copy_cmat%kill
        call odd_filt%kill
        call even_filt%kill
        call weights_img%kill
        call diff_img_odd%kill
        call diff_img_even%kill
        call diff_img_opt_odd%kill
        call diff_img_opt_even%kill
    end subroutine nonuni_filt3D

    subroutine estimate_lplim_1( odd, even, mskimg, kfromto, best_ind, odd_filt_out )
        class(image),           intent(inout) :: odd, even, mskimg
        integer,                intent(in)    :: kfromto(2)
        integer,                intent(out)   :: best_ind
        class(image), optional, intent(out)   :: odd_filt_out
        real,    allocatable :: cur_fil(:)
        logical, allocatable :: lmask(:,:,:)
        type(image) :: odd_copy_cmat, odd_filt
        integer     :: box, ldim(3), cutoff_find
        real        :: smpd, cur_cost, best_cost
        if( kfromto(1) == kfromto(2) )then  
            best_ind = kfromto(1)
            return
        endif
        if( odd%is_ft() .or. even%is_ft() ) THROW_HARD('Input even & odd has to be in real-space representation')
        if( mskimg%is_ft() ) THROW_HARD('Input mskimg has to be in real-space representation')
        ldim  = odd%get_ldim()
        box   = ldim(1)
        smpd  = odd%get_smpd()
        lmask = mskimg%bin2logical()
        call odd_copy_cmat%copy(odd)
        call odd_copy_cmat%fft
        call odd_filt%new(ldim, smpd)
        call odd_filt%set_ft(.true.)
        allocate(cur_fil(box), source=0.)
        best_ind  = kfromto(1)
        best_cost = huge(best_cost)
        if( ldim(3) > 1 )then ! 3D, so parallelize
            call odd_filt%set_wthreads(.true.)
            call even%set_wthreads(.true.)
        endif
        do cutoff_find = kfromto(1), kfromto(2)
            ! apply BW kernel to odd
            call odd_filt%copy_fast(odd_copy_cmat)
            call butterworth_filter(odd_filt, cutoff_find, cur_fil)
            call odd_filt%ifft
            ! calculate squared Euclidean distance
            cur_cost = odd_filt%sqeuclid(even, lmask)
            if( cur_cost <= best_cost )then
                best_cost = cur_cost
                best_ind  = cutoff_find
            endif
        enddo
        if( present(odd_filt_out) )then
            call odd_filt_out%new(ldim, smpd)
            call odd_filt_out%set_ft(.true.)
            call odd_filt_out%copy_fast(odd_copy_cmat)
            call butterworth_filter(odd_filt_out, best_ind, cur_fil)
            call odd_filt_out%ifft
        endif
        call odd_copy_cmat%kill
        call odd_filt%kill
    end subroutine estimate_lplim_1

    subroutine estimate_lplim_2( odd, even, mskimg, lprange, lpopt, odd_filt_out )
        class(image),           intent(inout) :: odd, even, mskimg
        real,                   intent(in)    :: lprange(2)
        real,                   intent(out)   :: lpopt
        class(image), optional, intent(out)   :: odd_filt_out
        integer :: ldim(3), box, kfromto(2), best_ind
        real    :: smpd
        ldim       = odd%get_ldim()
        box        = ldim(1)
        smpd       = odd%get_smpd()
        kfromto(1) = calc_fourier_index(lprange(1), box, smpd)
        kfromto(2) = calc_fourier_index(lprange(2), box, smpd)
        call estimate_lplim_1( odd, even, mskimg, kfromto, best_ind, odd_filt_out )
        lpopt = calc_lowpass_lim(best_ind, box, smpd)
    end subroutine estimate_lplim_2

    subroutine estimate_lplims2D( odd, even, mskrad_px, lprange, lpsopt, odd_filt_out )
        class(image),                       intent(inout) :: odd(:), even(:)
        real,                               intent(in)    :: mskrad_px, lprange(2)
        real,                  allocatable, intent(inout) :: lpsopt(:)
        type(image), optional, allocatable, intent(inout) :: odd_filt_out(:)
        type(image), allocatable :: masks(:)
        real    :: smpd
        integer :: iptcl, box, ldim(3), nptcls, npix
        write(logfhandle,'(A)') '>>> 2D UNIFORM FILTERING FOR LP ESTIMATION'
        ! init
        ldim    = odd(1)%get_ldim()
        ldim(3) = 1 ! because we operate on stacks
        box     = ldim(1)
        smpd    = odd(1)%get_smpd()
        nptcls  = size(odd)
        if( allocated(lpsopt) ) deallocate(lpsopt)
        allocate(masks(nptcls), lpsopt(nptcls))
        lpsopt = 0.
        call masks(1)%disc(ldim, smpd, mskrad_px, npix)
        do iptcl = 2, nptcls
            call masks(iptcl)%copy(masks(1))
        end do
        if( present(odd_filt_out) )then
            do iptcl = 1, size(odd_filt_out)
                call odd_filt_out(iptcl)%kill
            end do
            deallocate(odd_filt_out)
            allocate(odd_filt_out(nptcls))
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = 1, nptcls
                call estimate_lplim(odd(iptcl), even(iptcl), masks(iptcl), lprange, lpsopt(iptcl), odd_filt_out(iptcl))
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = 1, nptcls
                call estimate_lplim(odd(iptcl), even(iptcl), masks(iptcl), lprange, lpsopt(iptcl))
            enddo
            !$omp end parallel do
        endif
        do iptcl = 1, nptcls
            call masks(iptcl)%kill
        end do
        deallocate(masks)
    end subroutine estimate_lplims2D

end module simple_opt_filter

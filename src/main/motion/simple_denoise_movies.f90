!@descr: diffusion-map denoising helpers for real-space image stacks
module simple_denoise_movies
use simple_core_module_api
use simple_diff_map_graphs,  only: diffmap_graph, build_euclidean_knn_graph
use simple_diff_map_denoise, only: graph_nystrom_residual_preimage
use simple_diffusion_maps,   only: embed_graph, auto_ndiff_from_eigengap
use simple_image,            only: image
implicit none
private
#include "simple_local_flags.inc"

public :: diffmap_denoise_image_stack
public :: DEFAULT_DIFFMAP_STACK_KNN

integer, parameter :: DEFAULT_DIFFMAP_STACK_KNN = 6

contains

    subroutine diffmap_denoise_image_stack(frames, den_imgs, k_nn)
        type(image),              intent(inout) :: frames(:)
        type(image), allocatable, intent(out)   :: den_imgs(:)
        integer,     optional,    intent(in)    :: k_nn
        type(diffmap_graph) :: graph
        type(image)         :: avg_img
        real,   allocatable :: center(:), features(:,:), coords(:,:), eigvals(:), frame_vec(:)
        integer :: ldim(3), nframes_here, npix, iframe, rank_scan, rank_keep, k_nn_eff
        real    :: rms, smpd, feature_avg
        call validate_realspace_2d_stack(frames, ldim, smpd)
        nframes_here = size(frames)
        k_nn_eff = DEFAULT_DIFFMAP_STACK_KNN
        if( present(k_nn) ) k_nn_eff = k_nn
        npix = product(ldim)
        allocate(center(npix), features(npix,nframes_here))
        center = 0.0
        do iframe = 1,nframes_here
            frame_vec          = frames(iframe)%serialize()
            features(:,iframe) = frame_vec
            center             = center + frame_vec
            deallocate(frame_vec)
        end do
        center = center / real(nframes_here)
        !$omp parallel do default(shared) private(iframe,rms,feature_avg) schedule(static) proc_bind(close)
        do iframe = 1,nframes_here
            features(:,iframe) = features(:,iframe) - center    ! residuals
            feature_avg = sum(features(:,iframe)) / real(npix)
            features(:,iframe) = features(:,iframe) - feature_avg
            rms = sqrt(sum(features(:,iframe)**2) / real(npix))
            if( rms > TINY ) features(:,iframe) = features(:,iframe) / rms
        end do
        !$omp end parallel do
        call build_euclidean_knn_graph(features, min(max(2, k_nn_eff), nframes_here - 1), graph)
        if( graph%n /= nframes_here ) THROW_HARD('frame-stack graph size mismatch')
        rank_scan = frame_diffmap_rank_scan(nframes_here)
        call embed_graph(graph, rank_scan, coords, eigvals)
        rank_keep = select_frame_diffmap_rank(eigvals, size(coords, 1))
        call avg_img%new(ldim, smpd, wthreads=.false.)
        call avg_img%unserialize(center)
        call graph_nystrom_residual_preimage(frames, avg_img, graph, den_imgs, rank_keep)
        call avg_img%kill()
        call graph%kill()
        deallocate(center, features, coords, eigvals)
    end subroutine diffmap_denoise_image_stack

    subroutine validate_realspace_2d_stack(frames, ldim, smpd)
        type(image), intent(in) :: frames(:)
        integer, intent(out) :: ldim(3)
        real,    intent(out) :: smpd
        integer :: iframe, ldim_here(3), nframes_here
        nframes_here = size(frames)
        if( nframes_here < 3 ) THROW_HARD('at least 3 images are required for stack graph denoising')
        ldim = frames(1)%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('stack denoising expects 2D images')
        smpd = frames(1)%get_smpd()
        do iframe = 1,nframes_here
            ldim_here = frames(iframe)%get_ldim()
            if( any(ldim_here /= ldim) ) THROW_HARD('stack denoising image dimensions do not match')
            if( frames(iframe)%is_ft() ) THROW_HARD('stack denoising expects real-space images')
            if( abs(frames(iframe)%get_smpd() - smpd) > 1.e-6 * max(1.0, abs(smpd)) )then
                THROW_HARD('stack denoising image sampling distances do not match')
            endif
        end do
    end subroutine validate_realspace_2d_stack

    integer function frame_diffmap_rank_scan(nframes) result(rank_scan)
        integer, intent(in) :: nframes
        rank_scan = min(24, max(1, nframes - 2))
        if( nframes > 3 ) rank_scan = max(2, rank_scan)
    end function frame_diffmap_rank_scan

    integer function select_frame_diffmap_rank(eigvals, ncoords) result(rank_keep)
        real,    intent(in) :: eigvals(:)
        integer, intent(in) :: ncoords
        integer :: n
        n = min(size(eigvals), max(1, ncoords))
        if( n < 1 )then
            rank_keep = 1
            return
        endif
        rank_keep = auto_ndiff_from_eigengap(eigvals(:n))
        rank_keep = min(max(frame_diffmap_min_rank(n), rank_keep), n)
    end function select_frame_diffmap_rank

    integer function frame_diffmap_min_rank(max_neigs) result(rank_min)
        integer, intent(in) :: max_neigs
        rank_min = 1
        if( max_neigs >= 2 ) rank_min = 2
    end function frame_diffmap_min_rank

end module simple_denoise_movies

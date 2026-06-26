!@descr: diffusion-map denoising helpers for motion-correction movie stacks
module simple_denoise_movies
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, build_euclidean_knn_graph
use simple_diff_map_denoise, only: graph_nystrom_residual_preimage
use simple_diffusion_maps,  only: embed_graph, auto_ndiff_from_eigengap
use simple_image,           only: image
implicit none
private
#include "simple_local_flags.inc"

public :: diffmap_denoise_movie_tiles
public :: diffmap_denoise_movie_stack
public :: DEFAULT_DIFFMAP_MOVIE_KNN

integer, parameter :: DEFAULT_DIFFMAP_MOVIE_KNN = 6

contains

    subroutine diffmap_denoise_movie_tiles(frames, ntiles, overlap_pct, den_imgs, k_nn)
        type(image), intent(in) :: frames(:)
        integer,      intent(in) :: ntiles
        real,         intent(in) :: overlap_pct
        type(image), allocatable, intent(out) :: den_imgs(:)
        integer, optional, intent(in) :: k_nn
        type(image), allocatable :: tile_frames(:), den_tile(:)
        real, allocatable :: weights(:,:)
        integer, allocatable :: xstarts(:), ystarts(:)
        integer :: ldim(3)
        integer :: nx, ny, nframes_here, ix0, iy0, ix, iy, k_nn_eff, tile_x, tile_y
        real :: smpd
        call validate_movie_frame_stack(frames, ldim, smpd)
        nx = ldim(1)
        ny = ldim(2)
        nframes_here = size(frames)
        k_nn_eff = DEFAULT_DIFFMAP_MOVIE_KNN
        if( present(k_nn) ) k_nn_eff = k_nn
        call init_image_stack(den_imgs, nframes_here, ldim, smpd)
        allocate(weights(nx,ny), source=0.0)
        call make_ntile_starts(nx, ntiles, overlap_pct, tile_x, xstarts)
        call make_ntile_starts(ny, ntiles, overlap_pct, tile_y, ystarts)
        call init_image_stack(tile_frames, nframes_here, [tile_x,tile_y,1], smpd)
        do iy = 1,size(ystarts)
            iy0 = ystarts(iy)
            do ix = 1,size(xstarts)
                ix0 = xstarts(ix)
                call extract_tile_stack(frames, ix0, iy0, tile_x, tile_y, tile_frames)
                call diffmap_denoise_movie_stack(tile_frames, den_tile, k_nn_eff)
                call accumulate_tile_stack(den_tile, ix0, iy0, tile_x, tile_y, den_imgs)
                weights(ix0:ix0+tile_x-1,iy0:iy0+tile_y-1) = weights(ix0:ix0+tile_x-1,iy0:iy0+tile_y-1) + 1.0
                call kill_image_stack(den_tile)
            end do
        end do
        if( any(weights <= 0.0) ) THROW_HARD('tile stitching left uncovered pixels')
        call normalize_by_weights(den_imgs, weights)
        call kill_image_stack(tile_frames)
        deallocate(weights, xstarts, ystarts)
    end subroutine diffmap_denoise_movie_tiles

    subroutine diffmap_denoise_movie_stack(frames, den_imgs, k_nn)
        type(image), intent(inout) :: frames(:)
        type(image), allocatable, intent(out) :: den_imgs(:)
        integer, optional, intent(in) :: k_nn
        type(diffmap_graph) :: graph
        type(image) :: avg_img
        real, allocatable :: center(:), residuals(:,:), features(:,:), coords(:,:), eigvals(:), frame_vec(:)
        integer :: ldim(3), nframes_here, npix, iframe, rank_scan, rank_keep, k_nn_eff
        real :: rms
        real :: smpd
        call validate_movie_frame_stack(frames, ldim, smpd)
        nframes_here = size(frames)
        k_nn_eff = DEFAULT_DIFFMAP_MOVIE_KNN
        if( present(k_nn) ) k_nn_eff = k_nn
        npix = product(ldim)
        allocate(center(npix), residuals(npix,nframes_here), features(npix,nframes_here))
        center = 0.0
        do iframe = 1,nframes_here
            frame_vec = frames(iframe)%serialize()
            residuals(:,iframe) = frame_vec
            center = center + frame_vec
            deallocate(frame_vec)
        end do
        center = center / real(nframes_here)
        !$omp parallel do default(shared) private(iframe,rms) schedule(static) proc_bind(close)
        do iframe = 1,nframes_here
            residuals(:,iframe) = residuals(:,iframe) - center
            features(:,iframe)  = residuals(:,iframe)
            features(:,iframe)  = features(:,iframe) - sum(features(:,iframe)) / real(npix)
            rms = sqrt(sum(features(:,iframe)**2) / real(npix))
            if( rms > TINY ) features(:,iframe) = features(:,iframe) / rms
        end do
        !$omp end parallel do
        call build_euclidean_knn_graph(features, min(max(2, k_nn_eff), nframes_here - 1), 'none', graph)
        if( graph%n /= nframes_here ) THROW_HARD('frame-stack graph size mismatch')
        rank_scan = frame_diffmap_rank_scan(nframes_here)
        call embed_graph(graph, rank_scan, coords, eigvals)
        rank_keep = select_frame_diffmap_rank(eigvals, size(coords, 1))
        call avg_img%new(ldim, smpd, wthreads=.false.)
        call avg_img%unserialize(center)
        call graph_nystrom_residual_preimage(frames, avg_img, graph, den_imgs, rank_keep)
        call avg_img%kill()
        call graph%kill()
        deallocate(center, residuals, features, coords, eigvals)
    end subroutine diffmap_denoise_movie_stack

    subroutine validate_movie_frame_stack(frames, ldim, smpd)
        type(image), intent(in) :: frames(:)
        integer, intent(out) :: ldim(3)
        real,    intent(out) :: smpd
        integer :: iframe, ldim_here(3), nframes_here
        nframes_here = size(frames)
        if( nframes_here < 3 ) THROW_HARD('at least 3 frames are required for frame-stack graph denoising')
        ldim = frames(1)%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('movie denoising expects 2D image frames')
        smpd = frames(1)%get_smpd()
        do iframe = 1,nframes_here
            ldim_here = frames(iframe)%get_ldim()
            if( any(ldim_here /= ldim) ) THROW_HARD('movie denoising frame dimensions do not match')
            if( frames(iframe)%is_ft() ) THROW_HARD('movie denoising expects real-space frames')
            if( abs(frames(iframe)%get_smpd() - smpd) > 1.e-6 * max(1.0, abs(smpd)) )then
                THROW_HARD('movie denoising frame sampling distances do not match')
            endif
        end do
    end subroutine validate_movie_frame_stack

    subroutine init_image_stack(frames, nframes_here, ldim, smpd)
        type(image), allocatable, intent(out) :: frames(:)
        integer, intent(in) :: nframes_here, ldim(3)
        real,    intent(in) :: smpd
        integer :: iframe
        allocate(frames(nframes_here))
        do iframe = 1,nframes_here
            call frames(iframe)%new(ldim, smpd, wthreads=.false.)
            call frames(iframe)%zero()
        end do
    end subroutine init_image_stack

    subroutine kill_image_stack(frames)
        type(image), allocatable, intent(inout) :: frames(:)
        integer :: iframe
        if( .not. allocated(frames) ) return
        do iframe = 1,size(frames)
            if( frames(iframe)%exists() ) call frames(iframe)%kill()
        end do
        deallocate(frames)
    end subroutine kill_image_stack

    subroutine extract_tile_stack(frames, ix0, iy0, tile_x, tile_y, tile_frames)
        type(image), intent(in) :: frames(:)
        integer, intent(in) :: ix0, iy0, tile_x, tile_y
        type(image), intent(inout) :: tile_frames(:)
        integer :: iframe, k, l, ip, jp
        !$omp parallel do default(shared) private(iframe,k,l,ip,jp) schedule(static) proc_bind(close)
        do iframe = 1,size(frames)
            do l = iy0,iy0 + tile_y - 1
                jp = l - iy0 + 1
                do k = ix0,ix0 + tile_x - 1
                    ip = k - ix0 + 1
                    call tile_frames(iframe)%set([ip,jp,1], frames(iframe)%get([k,l,1]))
                end do
            end do
            call tile_frames(iframe)%set_ft(.false.)
        end do
        !$omp end parallel do
    end subroutine extract_tile_stack

    subroutine accumulate_tile_stack(tile_frames, ix0, iy0, tile_x, tile_y, out_frames)
        type(image), intent(in) :: tile_frames(:)
        integer, intent(in) :: ix0, iy0, tile_x, tile_y
        type(image), intent(inout) :: out_frames(:)
        real :: val
        integer :: iframe, k, l, ip, jp
        !$omp parallel do default(shared) private(iframe,k,l,ip,jp,val) schedule(static) proc_bind(close)
        do iframe = 1,size(out_frames)
            do l = iy0,iy0 + tile_y - 1
                jp = l - iy0 + 1
                do k = ix0,ix0 + tile_x - 1
                    ip = k - ix0 + 1
                    val = out_frames(iframe)%get([k,l,1]) + tile_frames(iframe)%get([ip,jp,1])
                    call out_frames(iframe)%set([k,l,1], val)
                end do
            end do
            call out_frames(iframe)%set_ft(.false.)
        end do
        !$omp end parallel do
    end subroutine accumulate_tile_stack

    subroutine make_ntile_starts(n, ntiles, overlap_pct, tile, starts)
        integer, intent(in) :: n, ntiles
        real,    intent(in) :: overlap_pct
        integer, intent(out) :: tile
        integer, allocatable, intent(out) :: starts(:)
        real :: overlap_frac, step_frac, denom, start_step
        integer :: max_start, i
        if( n < 4 ) THROW_HARD('frame dimension too small for frame-tile graph denoising')
        if( ntiles < 1 ) THROW_HARD('ntiles must be positive for frame-tile graph denoising')
        if( overlap_pct < 0.0 .or. overlap_pct >= 100.0 ) THROW_HARD('overlap percentage must be in [0,100)')
        if( ntiles == 1 )then
            tile = n
            allocate(starts(1))
            starts = 1
            return
        endif
        overlap_frac = overlap_pct / 100.0
        step_frac    = 1.0 - overlap_frac
        denom        = 1.0 + real(ntiles - 1) * step_frac
        tile         = ceiling(real(n) / denom)
        if( tile < 4 ) THROW_HARD('ntiles/overlap produce tiles smaller than four pixels')
        if( tile > n ) THROW_HARD('ntiles/overlap produce tiles larger than the frame dimension')
        max_start = n - tile + 1
        if( max_start < ntiles ) THROW_HARD('ntiles/overlap require duplicate tile starts')
        allocate(starts(ntiles))
        starts(1)      = 1
        starts(ntiles) = max_start
        start_step = real(max_start - 1) / real(ntiles - 1)
        do i = 2,ntiles - 1
            starts(i) = 1 + nint(real(i - 1) * start_step)
        end do
        do i = 2,ntiles
            if( starts(i) <= starts(i - 1) ) starts(i) = starts(i - 1) + 1
        end do
        do i = ntiles - 1,1,-1
            if( starts(i) >= starts(i + 1) ) starts(i) = starts(i + 1) - 1
        end do
    end subroutine make_ntile_starts

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

    subroutine normalize_by_weights(frames, weights)
        type(image), intent(inout) :: frames(:)
        real, intent(in) :: weights(:,:)
        real :: val
        integer :: iframe, k, l
        !$omp parallel do default(shared) private(iframe,k,l,val) schedule(static) proc_bind(close)
        do iframe = 1,size(frames)
            do l = 1,size(weights,2)
                do k = 1,size(weights,1)
                    val = frames(iframe)%get([k,l,1]) / weights(k,l)
                    call frames(iframe)%set([k,l,1], val)
                end do
            end do
            call frames(iframe)%set_ft(.false.)
        end do
        !$omp end parallel do
    end subroutine normalize_by_weights

end module simple_denoise_movies

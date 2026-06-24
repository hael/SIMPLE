program simple_test_frame_tile_diffmap_denoise
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, build_euclidean_knn_graph
use simple_diffusion_maps,  only: embed_graph, auto_ndiff_from_eigengap
use simple_image,           only: image
implicit none
#include "simple_local_flags.inc"

integer, parameter :: DEFAULT_TILE = 112, DEFAULT_STRIDE = 56
integer, parameter :: DEFAULT_KNN = 6
type(string) :: movie_dir
type(string) :: den_movie
type(string), allocatable :: movie_files(:)
real, allocatable :: raw(:,:,:), den(:,:,:)
integer :: tile_arg, stride_arg, knn_arg, imovie
real :: raw_temporal, den_temporal, avg_delta

call parse_args(movie_dir, tile_arg, stride_arg, knn_arg)
if( movie_dir .eq. '' )then
    call print_usage()
    stop
endif
if( .not. dir_exists(movie_dir) )then
    write(logfhandle,'(A,A)') 'Input movie stack folder not found: ', movie_dir%to_char()
    call print_usage()
    THROW_HARD('missing input movie stack folder')
endif

write(logfhandle,'(A,A)') 'frame-tile diffmap input folder: ', movie_dir%to_char()
call list_movie_stacks(movie_dir, movie_files)
write(logfhandle,'(A,I8)') 'frame-tile diffmap movie count:   ', size(movie_files)
write(logfhandle,'(A,I8,A,I8,A,I8)') &
    'frame-tile diffmap params: tile=', tile_arg, ' stride=', stride_arg, &
    ' k_nn=', knn_arg

do imovie = 1,size(movie_files)
    call read_mrc_stack(movie_files(imovie), raw)
    call denoise_frame_tiles(raw, tile_arg, stride_arg, knn_arg, den)
    den_movie = den_movie_name(movie_files(imovie))
    call write_mrc_stack(den_movie, movie_files(imovie), den)

    raw_temporal = temporal_residual_var(raw)
    den_temporal = temporal_residual_var(den)
    avg_delta    = frame_average_rms_delta(raw, den)

    write(logfhandle,'(A,I8,A,A)') 'frame-tile diffmap movie ', imovie, ': ', movie_files(imovie)%to_char()
    write(logfhandle,'(A,A)') 'frame-tile diffmap output: ', den_movie%to_char()
    write(logfhandle,'(A,I8)') 'frame-tile diffmap frames: ', size(raw,3)
    write(logfhandle,'(A,ES12.4,A,ES12.4,A,F8.2)') &
        'frame-tile temporal residual variance: raw=', raw_temporal, ' den=', den_temporal, &
        ' reduction_pct=', 100.0 * (raw_temporal - den_temporal) / max(raw_temporal, TINY)
    write(logfhandle,'(A,ES12.4)') 'frame average RMS delta raw-vs-den=', avg_delta
    if( den_temporal > raw_temporal * (1.0 + 1.0e-5) )then
        write(logfhandle,'(A)') 'frame-tile warning: denoised temporal residual variance increased'
    endif
    deallocate(raw, den)
end do

contains

    subroutine parse_args(folder, tile, stride, k_nn)
        type(string), intent(out) :: folder
        integer,      intent(out) :: tile, stride, k_nn
        character(len=STDLEN) :: arg, key, val
        integer :: i, eqpos
        folder = ''
        tile   = DEFAULT_TILE
        stride = DEFAULT_STRIDE
        k_nn   = DEFAULT_KNN
        do i = 1,command_argument_count()
            call get_command_argument(i, arg)
            eqpos = index(arg, '=')
            if( eqpos < 2 ) cycle
            key = lowercase(trim(arg(:eqpos-1)))
            val = trim(arg(eqpos+1:))
            select case(trim(key))
                case('dir', 'folder', 'movie_dir')
                    folder = trim(val)
                case('tile')
                    read(val,*) tile
                case('stride')
                    read(val,*) stride
                case('k_nn')
                    read(val,*) k_nn
                case DEFAULT
                    write(logfhandle,'(A,A)') 'Ignoring unknown argument: ', trim(key)
            end select
        end do
    end subroutine parse_args

    subroutine print_usage()
        write(logfhandle,'(A)') 'Usage: simple_test_frame_tile_diffmap_denoise ' // &
            'dir=/path/to/mrc_movie_stacks [tile=112] [stride=56] [k_nn=6]'
        write(logfhandle,'(A)') 'Aliases: folder= and movie_dir= are accepted instead of dir=.'
        write(logfhandle,'(A)') 'Writes each output stack next to the input movie as original_name_den.mrc.'
    end subroutine print_usage

    subroutine list_movie_stacks(folder, files)
        type(string), intent(in) :: folder
        type(string), allocatable, intent(out) :: files(:)
        logical, allocatable :: keep(:)
        type(string) :: bname
        character(len=STDLEN) :: bname_lower
        integer :: ifile
        call simple_list_files_regexp(folder, '\.[mM][rR][cC]$', files)
        if( .not. allocated(files) ) allocate(files(0))
        allocate(keep(size(files)), source=.true.)
        do ifile = 1,size(files)
            bname = basename(files(ifile))
            bname_lower = lowercase(bname%to_char())
            if( len_trim(bname_lower) >= len('_den.mrc') )then
                keep(ifile) = bname_lower(len_trim(bname_lower)-len('_den.mrc')+1:len_trim(bname_lower)) /= '_den.mrc'
            endif
            call bname%kill
        end do
        files = pack(files, keep)
        deallocate(keep)
        if( size(files) < 1 )then
            write(logfhandle,'(A,A)') 'No input .mrc movie stacks found in ', folder%to_char()
            THROW_HARD('no input .mrc movie stacks')
        endif
        call lex_sort(files)
    end subroutine list_movie_stacks

    subroutine read_mrc_stack(stkname, stack)
        type(string),      intent(in)  :: stkname
        real, allocatable, intent(out) :: stack(:,:,:)
        type(image) :: img
        integer :: ldim_here(3), nframes_here, iframe
        real, allocatable :: rmat(:,:,:)
        call find_ldim_nptcls(stkname, ldim_here, nframes_here)
        if( nframes_here < 3 )then
            write(logfhandle,'(A,A,A,I8)') 'Movie stack needs at least 3 frames: ', stkname%to_char(), '; found ', nframes_here
            THROW_HARD('not enough frames in input movie stack')
        endif
        ldim_here(3) = 1
        allocate(stack(ldim_here(1),ldim_here(2),nframes_here))
        call img%new(ldim_here, find_img_smpd(stkname))
        do iframe = 1,nframes_here
            call img%read(stkname, iframe)
            rmat = img%get_rmat()
            stack(:,:,iframe) = rmat(:,:,1)
            deallocate(rmat)
        end do
        call img%kill()
    end subroutine read_mrc_stack

    subroutine denoise_frame_tiles(frames, tile, stride, k_nn, out_frames)
        real,    intent(in)  :: frames(:,:,:)
        integer, intent(in)  :: tile, stride, k_nn
        real, allocatable, intent(out) :: out_frames(:,:,:)
        real, allocatable :: weights(:,:)
        real, allocatable :: den_tile(:,:,:)
        integer, allocatable :: xstarts(:), ystarts(:)
        integer :: nx, ny, nframes_here, ix0, iy0, ix, iy
        nx = size(frames, 1)
        ny = size(frames, 2)
        nframes_here = size(frames, 3)
        if( nframes_here < 3 ) THROW_HARD('at least 3 frames are required for frame-tile graph denoising')
        if( tile < 4 .or. stride < 1 ) THROW_HARD('invalid tile/stride for frame-tile denoising')
        if( tile > nx .or. tile > ny ) THROW_HARD('tile exceeds frame dimensions')
        allocate(out_frames(nx,ny,nframes_here), source=0.0)
        allocate(weights(nx,ny), source=0.0)
        call make_tile_starts(nx, tile, stride, xstarts)
        call make_tile_starts(ny, tile, stride, ystarts)
        do iy = 1,size(ystarts)
            iy0 = ystarts(iy)
            do ix = 1,size(xstarts)
                ix0 = xstarts(ix)
                call denoise_one_tile(frames(ix0:ix0+tile-1,iy0:iy0+tile-1,:), k_nn, den_tile)
                out_frames(ix0:ix0+tile-1,iy0:iy0+tile-1,:) = &
                    out_frames(ix0:ix0+tile-1,iy0:iy0+tile-1,:) + den_tile
                weights(ix0:ix0+tile-1,iy0:iy0+tile-1) = weights(ix0:ix0+tile-1,iy0:iy0+tile-1) + 1.0
                deallocate(den_tile)
            end do
        end do
        if( any(weights <= 0.0) ) THROW_HARD('tile stitching left uncovered pixels')
        call normalize_by_weights(out_frames, weights)
        deallocate(weights, xstarts, ystarts)
    end subroutine denoise_frame_tiles

    subroutine make_tile_starts(n, tile, stride, starts)
        integer, intent(in) :: n, tile, stride
        integer, allocatable, intent(out) :: starts(:)
        integer :: max_start, nregular, nstarts, i
        max_start = n - tile + 1
        nregular  = (max_start - 1) / stride + 1
        nstarts   = nregular
        if( 1 + (nregular - 1) * stride /= max_start ) nstarts = nstarts + 1
        allocate(starts(nstarts))
        do i = 1,nregular
            starts(i) = 1 + (i - 1) * stride
        end do
        if( nstarts > nregular ) starts(nstarts) = max_start
    end subroutine make_tile_starts

    subroutine denoise_one_tile(tile_stack, k_nn, den_tile)
        real,    intent(in)  :: tile_stack(:,:,:)
        integer, intent(in)  :: k_nn
        real, allocatable, intent(out) :: den_tile(:,:,:)
        type(diffmap_graph), target :: graph
        real, allocatable :: center(:), residuals(:,:), features(:,:), coords(:,:), eigvals(:)
        integer :: nx, ny, nframes_here, npix, iframe, rank_scan, rank_keep
        real :: rms
        nx = size(tile_stack, 1)
        ny = size(tile_stack, 2)
        nframes_here = size(tile_stack, 3)
        if( nframes_here < 3 ) THROW_HARD('at least 3 frames are required for frame-tile graph denoising')
        npix = nx * ny
        allocate(center(npix), residuals(npix,nframes_here), features(npix,nframes_here))
        center = 0.0
        do iframe = 1,nframes_here
            center = center + reshape(tile_stack(:,:,iframe), [npix])
        end do
        center = center / real(nframes_here)
        do iframe = 1,nframes_here
            residuals(:,iframe) = reshape(tile_stack(:,:,iframe), [npix]) - center
            features(:,iframe)  = residuals(:,iframe)
            features(:,iframe)  = features(:,iframe) - sum(features(:,iframe)) / real(npix)
            rms = sqrt(sum(features(:,iframe)**2) / real(npix))
            if( rms > TINY ) features(:,iframe) = features(:,iframe) / rms
        end do
        call build_euclidean_knn_graph(features, min(max(2, k_nn), nframes_here - 1), 'none', graph)
        if( graph%n /= nframes_here ) THROW_HARD('frame-tile graph size mismatch')
        rank_scan = frame_diffmap_rank_scan(nframes_here)
        call embed_graph(graph, rank_scan, coords, eigvals)
        rank_keep = select_frame_diffmap_rank(eigvals, size(coords, 1))
        call diffmap_nystrom_residual_preimage(residuals, center, graph, coords, rank_keep, nx, ny, den_tile)
        call graph%kill()
        deallocate(center, residuals, features, coords, eigvals)
    end subroutine denoise_one_tile

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

    subroutine diffmap_nystrom_residual_preimage(residuals, center, graph, coords, rank_keep, nx, ny, den_tile)
        real, intent(in) :: residuals(:,:), center(:)
        real, intent(in) :: coords(:,:)
        type(diffmap_graph), intent(in) :: graph
        integer, intent(in) :: rank_keep, nx, ny
        real, allocatable, intent(out) :: den_tile(:,:,:)
        real, allocatable :: den_residuals(:,:)
        integer, allocatable :: landmark_inds(:)
        integer :: i, j, p, n, nlands, rank_used
        real :: sigma2, w, wsum
        n = size(residuals, 2)
        if( graph%n /= n .or. size(coords, 2) /= n ) THROW_HARD('diffusion preimage frame count mismatch')
        rank_used = min(max(1, rank_keep), size(coords, 1))
        allocate(den_residuals(size(residuals, 1), n), source=0.0)
        nlands = max(2, min(max(1, n / 3), 20))
        allocate(landmark_inds(nlands))
        do i = 1,nlands
            landmark_inds(i) = 1 + ((i - 1) * (n - 1)) / max(1, nlands - 1)
        end do
        sigma2 = diffusion_coord_bandwidth2(coords, rank_used, graph)
        do i = 1,n
            wsum = 0.0
            do j = 1,nlands
                if( landmark_inds(j) == i ) cycle
                w = diffusion_coord_weight(coords, rank_used, i, landmark_inds(j), sigma2)
                if( w <= real(DTINY) ) cycle
                den_residuals(:,i) = den_residuals(:,i) + w * residuals(:,landmark_inds(j))
                wsum = wsum + w
            end do
            do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                j = graph%colind(p)
                if( j < 1 .or. j > n ) cycle
                if( j == i .or. graph%w(p) <= real(DTINY) ) cycle
                w = diffusion_coord_weight(coords, rank_used, i, j, sigma2)
                if( w <= real(DTINY) ) cycle
                den_residuals(:,i) = den_residuals(:,i) + w * residuals(:,j)
                wsum = wsum + w
            end do
            if( wsum > real(DTINY) )then
                den_residuals(:,i) = den_residuals(:,i) / wsum
            else
                den_residuals(:,i) = residuals(:,i)
            endif
        end do
        call residuals_to_tile_stack(den_residuals, center, nx, ny, den_tile)
        deallocate(den_residuals, landmark_inds)
    end subroutine diffmap_nystrom_residual_preimage

    real function diffusion_coord_bandwidth2(coords, rank_used, graph) result(sigma2)
        real, intent(in) :: coords(:,:)
        integer, intent(in) :: rank_used
        type(diffmap_graph), intent(in) :: graph
        integer :: i, j, p, n, cnt
        real :: d2
        n = size(coords, 2)
        sigma2 = 0.0
        cnt = 0
        do i = 1,n
            do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                j = graph%colind(p)
                if( j < 1 .or. j > n .or. j == i ) cycle
                d2 = diffusion_coord_dist2(coords, rank_used, i, j)
                if( d2 <= real(DTINY) ) cycle
                sigma2 = sigma2 + d2
                cnt = cnt + 1
            end do
        end do
        if( cnt > 0 )then
            sigma2 = sigma2 / real(cnt)
        else
            sigma2 = 1.0
        endif
        sigma2 = max(sigma2, 1.0e-6)
    end function diffusion_coord_bandwidth2

    real function diffusion_coord_weight(coords, rank_used, i, j, sigma2) result(w)
        real, intent(in) :: coords(:,:), sigma2
        integer, intent(in) :: rank_used, i, j
        real :: d2
        d2 = diffusion_coord_dist2(coords, rank_used, i, j)
        w = exp(-0.5 * d2 / max(sigma2, 1.0e-6))
    end function diffusion_coord_weight

    real function diffusion_coord_dist2(coords, rank_used, i, j) result(d2)
        real, intent(in) :: coords(:,:)
        integer, intent(in) :: rank_used, i, j
        integer :: k
        d2 = 0.0
        do k = 1,rank_used
            d2 = d2 + (coords(k,i) - coords(k,j))**2
        end do
    end function diffusion_coord_dist2

    subroutine residuals_to_tile_stack(residuals, center, nx, ny, den_tile)
        real, intent(in) :: residuals(:,:), center(:)
        integer, intent(in) :: nx, ny
        real, allocatable, intent(out) :: den_tile(:,:,:)
        integer :: iframe
        allocate(den_tile(nx,ny,size(residuals, 2)))
        do iframe = 1,size(residuals, 2)
            den_tile(:,:,iframe) = reshape(center + residuals(:,iframe), [nx,ny])
        end do
    end subroutine residuals_to_tile_stack

    subroutine normalize_by_weights(stack, weights)
        real, intent(inout) :: stack(:,:,:)
        real, intent(in)    :: weights(:,:)
        integer :: iframe
        do iframe = 1,size(stack, 3)
            stack(:,:,iframe) = stack(:,:,iframe) / weights
        end do
    end subroutine normalize_by_weights

    function den_movie_name(stkname) result(outname)
        type(string), intent(in) :: stkname
        type(string) :: outname
        type(string) :: stk_dir, den_base
        stk_dir = get_fpath(stkname)
        den_base = append2basename(stkname, '_den')
        outname = stk_dir%to_char()//den_base%to_char()
        call stk_dir%kill
        call den_base%kill
    end function den_movie_name

    subroutine write_mrc_stack(outstk, instk, stack)
        type(string), intent(in) :: outstk, instk
        real,         intent(in) :: stack(:,:,:)
        type(image) :: img
        integer :: iframe, ldim_here(3)
        ldim_here = [size(stack, 1), size(stack, 2), 1]
        call img%new(ldim_here, find_img_smpd(instk))
        do iframe = 1,size(stack, 3)
            call img%set_rmat(stack(:,:,iframe:iframe), .false.)
            call img%write(outstk, iframe, del_if_exists=(iframe == 1))
        end do
        call img%kill()
    end subroutine write_mrc_stack

    real function temporal_residual_var(stack) result(var)
        real, intent(in) :: stack(:,:,:)
        real, allocatable :: avg(:,:)
        integer :: iframe
        allocate(avg(size(stack,1),size(stack,2)), source=0.0)
        do iframe = 1,size(stack, 3)
            avg = avg + stack(:,:,iframe)
        end do
        avg = avg / real(size(stack, 3))
        var = 0.0
        do iframe = 1,size(stack, 3)
            var = var + sum((stack(:,:,iframe) - avg)**2)
        end do
        var = var / real(size(stack))
        deallocate(avg)
    end function temporal_residual_var

    real function frame_average_rms_delta(a, b) result(rms)
        real, intent(in) :: a(:,:,:), b(:,:,:)
        real, allocatable :: avga(:,:), avgb(:,:)
        integer :: iframe
        allocate(avga(size(a,1),size(a,2)), avgb(size(a,1),size(a,2)), source=0.0)
        do iframe = 1,size(a, 3)
            avga = avga + a(:,:,iframe)
            avgb = avgb + b(:,:,iframe)
        end do
        avga = avga / real(size(a, 3))
        avgb = avgb / real(size(b, 3))
        rms = sqrt(sum((avga - avgb)**2) / real(size(avga)))
        deallocate(avga, avgb)
    end function frame_average_rms_delta

end program simple_test_frame_tile_diffmap_denoise

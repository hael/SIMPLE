program simple_test_extract_dead_pixels
use simple_core_module_api
use simple_motion_correct_utils, only: extract_outliers
use simple_image, only: image
implicit none
#include "simple_local_flags.inc"

character(len=STDLEN)         :: movie_file
character(len=STDLEN)         :: gainref_file
character(len=:), allocatable :: smpd_char
real                          :: smpd
integer                       :: slen, argc, iarg, n_movies
integer                       :: n_low_common, n_high_common
integer, allocatable          :: low_prod(:,:), high_prod(:,:), low_cur(:,:), high_cur(:,:)
integer                       :: n_common_total, n_common_prev, n_same_streak, n_processed
integer, parameter            :: N_CONVERGED = 5
logical, allocatable          :: low_ref(:,:), high_ref(:,:), low_here(:,:), high_here(:,:)
logical, allocatable          :: low_common_mask(:,:), high_common_mask(:,:)
type(image)                   :: gain_ref, gain_flip_x, gain_flip_y, gain_flip_xy
real, allocatable             :: gain_ref_rmat(:,:,:), gain_flip_x_rmat(:,:,:), gain_flip_y_rmat(:,:,:), gain_flip_xy_rmat(:,:,:)
real                          :: gain_sum_low_ref, gain_sum_high_ref, gain_sum_low_x, gain_sum_high_x
real                          :: gain_sum_low_y, gain_sum_high_y, gain_sum_low_xy, gain_sum_high_xy
real                          :: best_high_sum, best_low_sum
real                          :: gain_smpd
integer                       :: gain_ldim(3), ifoo
integer                       :: best_idx

argc = command_argument_count()
if( argc < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_extract_dead_pixels smpd gainref movie1 [movie2 ...]'
    write(logfhandle,'(a)') 'smpd   : sampling distance in Angstrom/pixel'
    write(logfhandle,'(a)') 'gainref: gain reference image'
    write(logfhandle,'(a)') 'movieN : input movie stacks (non-EER currently supported)'
    stop
endif

call get_command_argument(1, length=slen)
allocate(character(slen) :: smpd_char)
call get_command_argument(1, smpd_char)
read(smpd_char, *) smpd
call get_command_argument(2, gainref_file)

n_movies = argc - 2
n_processed  = 0
n_same_streak = 0
n_common_prev = -1

do iarg=3, argc
    call get_command_argument(iarg, movie_file)
    call extract_outliers(string(trim(movie_file)), smpd, low_here, high_here)
    n_processed = n_processed + 1
    write(logfhandle,'(a,1x,a)') '>>> Processed movie:', trim(movie_file)
    write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> LOW/HIGH outliers:', count(low_here), '/', count(high_here)

    if( .not. allocated(low_ref) )then
        allocate(low_ref(size(low_here,1), size(low_here,2)), source=low_here)
        allocate(high_ref(size(high_here,1), size(high_here,2)), source=high_here)
        allocate(low_prod(size(low_here,1), size(low_here,2)))
        allocate(high_prod(size(high_here,1), size(high_here,2)))
        low_prod  = merge(1, 0, low_here)
        high_prod = merge(1, 0, high_here)
    else
        if( size(low_here,1) /= size(low_ref,1) .or. size(low_here,2) /= size(low_ref,2) )then
            THROW_HARD('Movie dimensions differ from reference movie: '//trim(movie_file))
        endif
        allocate(low_cur(size(low_here,1), size(low_here,2)))
        allocate(high_cur(size(high_here,1), size(high_here,2)))
        low_cur  = merge(1, 0, low_here)
        high_cur = merge(1, 0, high_here)
        low_prod  = low_prod  * low_cur
        high_prod = high_prod * high_cur
        deallocate(low_cur, high_cur)
    endif

    if( allocated(low_here) )  deallocate(low_here)
    if( allocated(high_here) ) deallocate(high_here)
    n_low_common  = count(low_prod  /= 0)
    n_high_common = count(high_prod /= 0)
    n_common_total = n_low_common + n_high_common
    if( n_common_total == n_common_prev )then
        n_same_streak = n_same_streak + 1
    else
        n_same_streak = 0
        n_common_prev = n_common_total
    endif
    write(logfhandle,'(a,1x,i0)') '>>> COMMON LOW-OUTLIER PIXELS (non-zero product):', n_low_common
    write(logfhandle,'(a,1x,i0)') '>>> COMMON HIGH-OUTLIER PIXELS (non-zero product):', n_high_common
    write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> COMMON OUTLIER COUNT STREAK:', n_same_streak, 'target=', N_CONVERGED
    write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> MOVIES PROCESSED:', n_processed, 'of', n_movies
    if( n_same_streak >= N_CONVERGED )then
        write(logfhandle,'(a,1x,i0,1x,a)') '>>> CONVERGED after', n_processed, 'movies (no change for 5 micrographs)'
        exit
    endif
    
end do

n_low_common  = count(low_prod  /= 0)
n_high_common = count(high_prod /= 0)
write(logfhandle,'(a,1x,i0)') '>>> FINAL COMMON LOW-OUTLIER PIXELS (non-zero product):', n_low_common
write(logfhandle,'(a,1x,i0)') '>>> FINAL COMMON HIGH-OUTLIER PIXELS (non-zero product):', n_high_common
write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> FINAL MOVIES PROCESSED:', n_processed, 'of', n_movies

allocate(low_common_mask(size(low_prod,1), size(low_prod,2)), source=(low_prod /= 0))
allocate(high_common_mask(size(high_prod,1), size(high_prod,2)), source=(high_prod /= 0))

if( .not. file_exists(string(trim(gainref_file))) )then
    THROW_HARD('Gain reference file does not exist: '//trim(gainref_file))
endif
call find_ldim_nptcls(string(trim(gainref_file)), gain_ldim, ifoo)
gain_smpd = find_img_smpd(string(trim(gainref_file)))
gain_ldim(3) = 1

call gain_ref%new(gain_ldim, gain_smpd, wthreads=.false.)
call gain_ref%read(string(trim(gainref_file)))
gain_ref_rmat = gain_ref%get_rmat()

if( size(gain_ref_rmat,1) /= size(low_common_mask,1) .or. size(gain_ref_rmat,2) /= size(low_common_mask,2) )then
    THROW_HARD('Gain reference dimensions differ from outlier map dimensions')
endif

call gain_flip_x%new(gain_ldim, gain_smpd, wthreads=.false.)
call gain_flip_x%read(string(trim(gainref_file)))
call gain_flip_x%flip('X')
gain_flip_x_rmat = gain_flip_x%get_rmat()

call gain_flip_y%new(gain_ldim, gain_smpd, wthreads=.false.)
call gain_flip_y%read(string(trim(gainref_file)))
call gain_flip_y%flip('Y')
gain_flip_y_rmat = gain_flip_y%get_rmat()

call gain_flip_xy%new(gain_ldim, gain_smpd, wthreads=.false.)
call gain_flip_xy%read(string(trim(gainref_file)))
call gain_flip_xy%flip('XY')
gain_flip_xy_rmat = gain_flip_xy%get_rmat()

gain_sum_low_ref  = sum(gain_ref_rmat(:,:,1),     mask=low_common_mask)
gain_sum_high_ref = sum(gain_ref_rmat(:,:,1),     mask=high_common_mask)
gain_sum_low_x    = sum(gain_flip_x_rmat(:,:,1),  mask=low_common_mask)
gain_sum_high_x   = sum(gain_flip_x_rmat(:,:,1),  mask=high_common_mask)
gain_sum_low_y    = sum(gain_flip_y_rmat(:,:,1),  mask=low_common_mask)
gain_sum_high_y   = sum(gain_flip_y_rmat(:,:,1),  mask=high_common_mask)
gain_sum_low_xy   = sum(gain_flip_xy_rmat(:,:,1), mask=low_common_mask)
gain_sum_high_xy  = sum(gain_flip_xy_rmat(:,:,1), mask=high_common_mask)

write(logfhandle,'(a)') '>>> GAIN SUMS ON FINAL COMMON OUTLIER MAPS'
write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> GAIN[UNCHANGED] low/high:', gain_sum_low_ref, '/', gain_sum_high_ref
write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> GAIN[FLIP_X]    low/high:', gain_sum_low_x,   '/', gain_sum_high_x
write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> GAIN[FLIP_Y]    low/high:', gain_sum_low_y,   '/', gain_sum_high_y
write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> GAIN[FLIP_XY]   low/high:', gain_sum_low_xy,  '/', gain_sum_high_xy

best_idx      = 1
best_high_sum = gain_sum_high_ref
best_low_sum  = gain_sum_low_ref

if( gain_sum_high_x < best_high_sum .or. (gain_sum_high_x == best_high_sum .and. gain_sum_low_x > best_low_sum) )then
    best_idx      = 2
    best_high_sum = gain_sum_high_x
    best_low_sum  = gain_sum_low_x
endif
if( gain_sum_high_y < best_high_sum .or. (gain_sum_high_y == best_high_sum .and. gain_sum_low_y > best_low_sum) )then
    best_idx      = 3
    best_high_sum = gain_sum_high_y
    best_low_sum  = gain_sum_low_y
endif
if( gain_sum_high_xy < best_high_sum .or. (gain_sum_high_xy == best_high_sum .and. gain_sum_low_xy > best_low_sum) )then
    best_idx      = 4
    best_high_sum = gain_sum_high_xy
    best_low_sum  = gain_sum_low_xy
endif

select case(best_idx)
case(1)
    call gain_ref%write(string('best_gainref_highsum_lowlow.mrc'))
    write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> WROTE BEST GAIN: unchanged high/low=', best_high_sum, '/', best_low_sum
case(2)
    call gain_flip_x%write(string('best_gainref_highsum_lowlow.mrc'))
    write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> WROTE BEST GAIN: flip_x high/low=', best_high_sum, '/', best_low_sum
case(3)
    call gain_flip_y%write(string('best_gainref_highsum_lowlow.mrc'))
    write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> WROTE BEST GAIN: flip_y high/low=', best_high_sum, '/', best_low_sum
case(4)
    call gain_flip_xy%write(string('best_gainref_highsum_lowlow.mrc'))
    write(logfhandle,'(a,1x,f12.5,1x,a,1x,f12.5)') '>>> WROTE BEST GAIN: flip_xy high/low=', best_high_sum, '/', best_low_sum
end select

if( allocated(low_ref) )  deallocate(low_ref)
if( allocated(high_ref) ) deallocate(high_ref)
if( allocated(low_prod) ) deallocate(low_prod)
if( allocated(high_prod) ) deallocate(high_prod)
if( allocated(low_common_mask) ) deallocate(low_common_mask)
if( allocated(high_common_mask) ) deallocate(high_common_mask)
if( allocated(gain_ref_rmat) ) deallocate(gain_ref_rmat)
if( allocated(gain_flip_x_rmat) ) deallocate(gain_flip_x_rmat)
if( allocated(gain_flip_y_rmat) ) deallocate(gain_flip_y_rmat)
if( allocated(gain_flip_xy_rmat) ) deallocate(gain_flip_xy_rmat)
call gain_ref%kill()
call gain_flip_x%kill()
call gain_flip_y%kill()
call gain_flip_xy%kill()
if( allocated(smpd_char) ) deallocate(smpd_char)

end program simple_test_extract_dead_pixels

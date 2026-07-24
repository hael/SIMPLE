program simple_test_search_gain_flips
use simple_core_module_api
use simple_image,                          only: image
use simple_motion_gain_analysis_helpers,   only: read_movies_and_sum_frames
use simple_motion_gain_analysis,           only: gain_flip_analyzer, GAIN_FIRST_ANALYSIS_AT
implicit none
#include "simple_local_flags.inc"

character(len=STDLEN) :: gainref_file, smpd_char, movie_file
real                  :: smpd_avg
integer               :: argc, iarg, n_movies, n_movies_cap
integer               :: ibeg, iend, n_movies_batch, n_movies_done
integer               :: total_frames_batch, total_frames_all
integer, parameter    :: MAX_MOVIES = 200
integer, parameter    :: BATCH_MOVIES = 10
logical               :: ran_analysis
type(string), allocatable :: movie_fnames(:)
type(image)                 :: sum_img
type(gain_flip_analyzer)    :: analyzer

argc = command_argument_count()
if( argc < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_search_gain_flips gainref smpd movie1 [movie2 ...]'
    write(logfhandle,'(a)') 'gainref: gain reference image'
    write(logfhandle,'(a)') 'smpd   : sampling distance in Angstrom/pixel'
    write(logfhandle,'(a)') 'movieN : input movie stacks (non-EER currently supported)'
    stop
endif

call get_command_argument(1, gainref_file)
call get_command_argument(2, smpd_char)
read(smpd_char, *) smpd_avg

n_movies     = argc - 2
n_movies_cap = min(n_movies, MAX_MOVIES)

write(logfhandle,'(a,1x,f10.5)') '>>> Input smpd:', smpd_avg
write(logfhandle,'(a,1x,i0)') '>>> Total input movies:', n_movies
write(logfhandle,'(a,1x,i0)') '>>> Max movies to process:', n_movies_cap

allocate(movie_fnames(n_movies_cap))
do iarg=1,n_movies_cap
    call get_command_argument(iarg + 2, movie_file)
    movie_fnames(iarg) = string(trim(movie_file))
    write(logfhandle,'(a,1x,a,1x,a,1x,f10.5)') '>>> smpd_movie for', trim(movie_file), '=', smpd_avg
end do

call analyzer%new(string(trim(gainref_file)), smpd_avg)
n_movies_done    = 0
total_frames_all = 0
do ibeg=1,n_movies_cap,BATCH_MOVIES
    iend = min(ibeg + BATCH_MOVIES - 1, n_movies_cap)

    call read_movies_and_sum_frames(movie_fnames(ibeg:iend), smpd_avg, sum_img, n_movies_batch, total_frames_batch)
    n_movies_done    = n_movies_done + n_movies_batch
    total_frames_all = total_frames_all + total_frames_batch

    write(logfhandle,'(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0)') '>>> Batch processed movies', ibeg, 'to', iend, 'frames=', total_frames_batch
    call analyzer%analyze_if_due(sum_img, total_frames_batch, n_movies_batch, ran_analysis)
    call sum_img%kill()
    if( ran_analysis .and. analyzer%get_converged() ) exit
end do

if( total_frames_all < 1 )then
    THROW_HARD('No frames were read from input movies')
endif

if( n_movies_done < GAIN_FIRST_ANALYSIS_AT )then
    write(logfhandle,'(a,1x,i0,1x,a,1x,i0,1x,a)') '>>> Stopped after', n_movies_done, &
    'movies; first analysis requires', GAIN_FIRST_ANALYSIS_AT, 'movies'
endif
if( .not. analyzer%get_converged() .and. n_movies_done >= n_movies_cap )then
    if( n_movies_cap == MAX_MOVIES )then
        write(logfhandle,'(a,1x,i0)') '>>> Reached hard stop at movies:', MAX_MOVIES
    else
        write(logfhandle,'(a,1x,i0)') '>>> Reached end of provided movies:', n_movies_done
    endif
endif

call analyzer%kill()
deallocate(movie_fnames)

end program simple_test_search_gain_flips

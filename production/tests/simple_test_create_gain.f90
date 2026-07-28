program simple_test_create_gain
use simple_core_module_api
use simple_image,                          only: image
use simple_motion_gain_helpers,            only: read_movies_and_sum_frames, normalized_inverse_average_intensity
implicit none
#include "simple_local_flags.inc"

character(len=STDLEN) :: gainref_file, smpd_char, movie_file
real                  :: smpd_avg, avg_value
integer               :: argc, iarg, n_movies, n_movies_cap
integer               :: n_movies_read, total_frames_all
integer, parameter    :: MAX_MOVIES = 200
type(string), allocatable :: movie_fnames(:)
type(string)                :: gainref_out, gainref_ext, gainref_base, flip_out
type(image)                 :: sum_img
type(image)                 :: gain_img
type(image)                 :: gain_flip

argc = command_argument_count()
if( argc < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_create_gain gainref smpd movie1 [movie2 ...]'
    write(logfhandle,'(a)') 'gainref: output gain reference image'
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

call read_movies_and_sum_frames(movie_fnames, smpd_avg, sum_img, n_movies_read, total_frames_all)
write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> Processed movies:', n_movies_read, 'frames=', total_frames_all

if( total_frames_all < 1 )then
    THROW_HARD('No frames were read from input movies')
endif

call normalized_inverse_average_intensity(sum_img, total_frames_all, gain_img, avg_value)
write(logfhandle,'(a,1x,g12.5)') '>>> Average intensity:', avg_value
gainref_out = string(trim(gainref_file))
call gain_img%write(gainref_out, del_if_exists=.true.)
gainref_ext  = fname2ext(gainref_out)
gainref_base = get_fbody(gainref_out, gainref_ext, separator=.true.)

call gain_flip%new(gain_img%get_ldim(), gain_img%get_smpd(), wthreads=.false.)
call gain_flip%copy(gain_img)
call gain_flip%flip('X')
flip_out = gainref_base%to_char()//'_flipX.mrc'
call gain_flip%write(flip_out, del_if_exists=.true.)

call gain_flip%copy(gain_img)
call gain_flip%flip('Y')
flip_out = gainref_base%to_char()//'_flipY.mrc'
call gain_flip%write(flip_out, del_if_exists=.true.)

call gain_flip%copy(gain_img)
call gain_flip%flip('XY')
flip_out = gainref_base%to_char()//'_flipXY.mrc'
call gain_flip%write(flip_out, del_if_exists=.true.)

write(logfhandle,'(a,1x,g12.5)') '>>> Global avg_value:', avg_value
write(logfhandle,'(a,1x,a)') '>>> Wrote gain image:', trim(gainref_file)
write(logfhandle,'(a)') '>>> Wrote gain flips: _flipX, _flipY, _flipXY'

call gain_flip%kill()
call gain_img%kill()
call sum_img%kill()

if( n_movies_cap == MAX_MOVIES .and. n_movies_read >= n_movies_cap )then
    write(logfhandle,'(a,1x,i0)') '>>> Reached hard stop at movies:', MAX_MOVIES
else
    write(logfhandle,'(a,1x,i0)') '>>> Reached end of provided movies:', n_movies_read
endif
deallocate(movie_fnames)

end program simple_test_create_gain

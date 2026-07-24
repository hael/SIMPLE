program simple_test_extract_dead_pixels
use simple_core_module_api
use simple_image, only: image
use simple_defs, only: PI
implicit none
#include "simple_local_flags.inc"

character(len=STDLEN)         :: gainref_file
character(len=STDLEN)         :: smpd_char
character(len=STDLEN)         :: movie_file
type(image)                   :: avg_img, frame
type(image)                   :: gain_ref, gain_flip_x, gain_flip_y, gain_flip_xy
type(image)                   :: tmp_avg, tmp_gain
real                          :: smpd_avg, smpd_movie, corr_ref, corr_x, corr_y, corr_xy
integer                       :: argc, iarg, n_movies, nframes, iframe, ifoo
integer                       :: ldim_ref(3), ldim_cur(3)
integer                       :: total_frames
integer                       :: best_idx
integer                       :: prev_best_idx
logical                       :: have_ref_shape
real                          :: best_corr
real                          :: prev_best_corr, prev_delta_r
real                          :: second_corr, delta_r, z_best, z_second, z_gap, fisher_den, r_clip
real                          :: corrs(4)
real                          :: z_ref, z_x, z_y, z_xy, zsep_ref, zsep_x, zsep_y, zsep_xy
real, parameter               :: R_EPS = 1.e-6
real, parameter               :: LP_RES_A = 200.
real                          :: corr_len_pix, corr_area_pix
integer                       :: iori, npix_eff, npix_tot
integer                       :: n_movies_cap, n_movies_processed, analysis_count, conv_streak
integer, parameter            :: FIRST_ANALYSIS_AT = 10
integer, parameter            :: ANALYSIS_STEP = 5
integer, parameter            :: MAX_MOVIES = 200
integer, parameter            :: CONV_STREAK_REQ = 2
real,    parameter            :: CONV_EPS_CORR = 1.e-2
real,    parameter            :: CONV_EPS_DELTA = 1.e-2
logical                       :: have_prev_analysis, converged

argc = command_argument_count()
if( argc < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_extract_dead_pixels gainref smpd movie1 [movie2 ...]'
    write(logfhandle,'(a)') 'gainref: gain reference image'
    write(logfhandle,'(a)') 'smpd   : sampling distance in Angstrom/pixel'
    write(logfhandle,'(a)') 'movieN : input movie stacks (non-EER currently supported)'
    stop
endif

call get_command_argument(1, gainref_file)
call get_command_argument(2, smpd_char)
read(smpd_char, *) smpd_avg
n_movies = argc - 2
n_movies_cap = min(n_movies, MAX_MOVIES)
total_frames = 0
have_ref_shape = .false.
have_prev_analysis = .false.
converged = .false.
n_movies_processed = 0
analysis_count = 0
conv_streak = 0
write(logfhandle,'(a,1x,f10.5)') '>>> Input smpd:', smpd_avg
write(logfhandle,'(a,1x,i0)') '>>> Total input movies:', n_movies
write(logfhandle,'(a,1x,i0)') '>>> Max movies to process:', n_movies_cap

if( .not. file_exists(string(trim(gainref_file))) )then
    THROW_HARD('Gain reference file does not exist: '//trim(gainref_file))
endif

do iarg=3, argc
    if( n_movies_processed >= n_movies_cap ) exit
    call get_command_argument(iarg, movie_file)
    if( fname2format(string(trim(movie_file))) == 'K' )then
        THROW_HARD('EER movies are not supported by this test: '//trim(movie_file))
    endif
    call find_ldim_nptcls(string(trim(movie_file)), ldim_cur, nframes)
    smpd_movie = smpd_avg
    write(logfhandle,'(a,1x,a,1x,a,1x,f10.5)') '>>> smpd_movie for', trim(movie_file), '=', smpd_movie
    
    ldim_cur(3) = 1
    if( .not. have_ref_shape )then
        ldim_ref = ldim_cur
        call avg_img%new(ldim_ref, smpd_avg, wthreads=.false.)
        call avg_img%zero()
        have_ref_shape = .true.

        call find_ldim_nptcls(string(trim(gainref_file)), ldim_cur, ifoo)
        ldim_cur(3) = 1
        if( ldim_cur(1) /= ldim_ref(1) .or. ldim_cur(2) /= ldim_ref(2) )then
            THROW_HARD('Gain reference dimensions differ from movie average dimensions')
        endif

        call gain_ref%new(ldim_ref, smpd_avg, wthreads=.false.)
        call gain_ref%read(string(trim(gainref_file)))
        call gain_ref%bp(0., LP_RES_A)
        call gain_flip_x%new(ldim_ref, smpd_avg, wthreads=.false.)
        call gain_flip_x%copy(gain_ref)
        call gain_flip_x%flip('X')
        call gain_flip_y%new(ldim_ref, smpd_avg, wthreads=.false.)
        call gain_flip_y%copy(gain_ref)
        call gain_flip_y%flip('Y')
        call gain_flip_xy%new(ldim_ref, smpd_avg, wthreads=.false.)
        call gain_flip_xy%copy(gain_ref)
        call gain_flip_xy%flip('XY')

        call tmp_avg%new(ldim_ref, smpd_avg, wthreads=.false.)
        call tmp_gain%new(ldim_ref, smpd_avg, wthreads=.false.)
    else
        if( ldim_cur(1) /= ldim_ref(1) .or. ldim_cur(2) /= ldim_ref(2) )then
            THROW_HARD('Movie dimensions differ from reference: '//trim(movie_file))
        endif
    endif

    call frame%new(ldim_ref, smpd_avg, wthreads=.false.)
    do iframe=1,nframes
        call frame%read(string(trim(movie_file)), iframe)
        call avg_img%add_workshare(frame)
        total_frames = total_frames + 1
    enddo
    call frame%kill()
    n_movies_processed = n_movies_processed + 1
    write(logfhandle,'(a,1x,a,1x,a,1x,i0)') '>>> Added movie', trim(movie_file), 'frames=', nframes

    if( n_movies_processed >= FIRST_ANALYSIS_AT .and. mod(n_movies_processed - FIRST_ANALYSIS_AT, ANALYSIS_STEP) == 0 )then
        analysis_count = analysis_count + 1
        write(logfhandle,'(a)') '>>> ------------------------------------------------------------'
        write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> ANALYSIS', analysis_count, 'at movie count', n_movies_processed

        call tmp_avg%copy(avg_img)
        call tmp_avg%div(real(total_frames))
        call tmp_avg%bp(0., LP_RES_A)
        write(logfhandle,'(a,1x,i0)') '>>> TOTAL FRAMES AVERAGED:', total_frames
        write(logfhandle,'(a,1x,f8.2,a)') '>>> Applied low-pass filter to average:', LP_RES_A, 'A'
        call tmp_avg%write(string('average_all_movies.mrc'))

        call tmp_gain%copy(gain_ref)
        corr_ref = tmp_avg%corr(tmp_gain)
        call tmp_gain%copy(gain_flip_x)
        corr_x = tmp_avg%corr(tmp_gain)
        call tmp_gain%copy(gain_flip_y)
        corr_y = tmp_avg%corr(tmp_gain)
        call tmp_gain%copy(gain_flip_xy)
        corr_xy = tmp_avg%corr(tmp_gain)

        write(logfhandle,'(a)') '>>> CORRELATION TABLE: average_all_movies vs gain variants'
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_unchanged :', corr_ref
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_flip_x    :', corr_x
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_flip_y    :', corr_y
        write(logfhandle,'(a,1x,f12.6)') '>>> gain_flip_xy   :', corr_xy

        corrs = [corr_ref, corr_x, corr_y, corr_xy]
        best_idx = 1
        best_corr = corr_ref
        if( corr_x < best_corr )then
            best_idx = 2
            best_corr = corr_x
        endif
        if( corr_y < best_corr )then
            best_idx = 3
            best_corr = corr_y
        endif
        if( corr_xy < best_corr )then
            best_idx = 4
            best_corr = corr_xy
        endif

        second_corr = huge(second_corr)
        do iori=1,4
            if( iori == best_idx ) cycle
            if( corrs(iori) < second_corr ) second_corr = corrs(iori)
        enddo
        delta_r = second_corr - best_corr

        npix_tot = ldim_ref(1) * ldim_ref(2)
        corr_len_pix = max(1.0, 0.5 * LP_RES_A / max(smpd_avg, 1.e-6))
        corr_area_pix = PI * corr_len_pix * corr_len_pix
        npix_eff = max(4, int(real(npix_tot) / corr_area_pix))
        r_clip = max(-1. + R_EPS, min(1. - R_EPS, best_corr))
        z_best = 0.5 * log((1. + r_clip) / (1. - r_clip))
        r_clip = max(-1. + R_EPS, min(1. - R_EPS, second_corr))
        z_second = 0.5 * log((1. + r_clip) / (1. - r_clip))

        r_clip = max(-1. + R_EPS, min(1. - R_EPS, corr_ref))
        z_ref = 0.5 * log((1. + r_clip) / (1. - r_clip))
        r_clip = max(-1. + R_EPS, min(1. - R_EPS, corr_x))
        z_x = 0.5 * log((1. + r_clip) / (1. - r_clip))
        r_clip = max(-1. + R_EPS, min(1. - R_EPS, corr_y))
        z_y = 0.5 * log((1. + r_clip) / (1. - r_clip))
        r_clip = max(-1. + R_EPS, min(1. - R_EPS, corr_xy))
        z_xy = 0.5 * log((1. + r_clip) / (1. - r_clip))

        if( npix_eff > 3 )then
            fisher_den = sqrt(2. / real(npix_eff - 3))
            z_gap = (z_second - z_best) / fisher_den
            zsep_ref = (z_ref - z_best) / fisher_den
            zsep_x   = (z_x   - z_best) / fisher_den
            zsep_y   = (z_y   - z_best) / fisher_den
            zsep_xy  = (z_xy  - z_best) / fisher_den
        else
            fisher_den = 0.
            z_gap = 0.
            zsep_ref = 0.
            zsep_x   = 0.
            zsep_y   = 0.
            zsep_xy  = 0.
        endif

        write(logfhandle,'(a,1x,f12.6)') '>>> CONFIDENCE delta_r (2nd-best - best):', delta_r
        write(logfhandle,'(a,1x,f12.6)') '>>> CONFIDENCE fisher_z_gap:', z_second - z_best
        write(logfhandle,'(a,1x,f10.3,1x,a,1x,f10.3)') '>>> CONFIDENCE corr_len_pix/area_pix:', corr_len_pix, '/', corr_area_pix
        write(logfhandle,'(a,1x,f12.6,1x,a,1x,i0,1x,a,1x,i0)') '>>> CONFIDENCE z_score:', z_gap, 'npix_eff=', npix_eff, 'npix_tot=', npix_tot
        write(logfhandle,'(a)') '>>> SCORES FOR ALL FLIPS (corr, fisher_z, zsep_vs_best)'
        write(logfhandle,'(a,1x,f12.6,1x,f12.6,1x,f12.6)') '>>> gain_unchanged :', corr_ref, z_ref, zsep_ref
        write(logfhandle,'(a,1x,f12.6,1x,f12.6,1x,f12.6)') '>>> gain_flip_x    :', corr_x,   z_x,   zsep_x
        write(logfhandle,'(a,1x,f12.6,1x,f12.6,1x,f12.6)') '>>> gain_flip_y    :', corr_y,   z_y,   zsep_y
        write(logfhandle,'(a,1x,f12.6,1x,f12.6,1x,f12.6)') '>>> gain_flip_xy   :', corr_xy,  z_xy,  zsep_xy

        select case(best_idx)
        case(1)
            call gain_ref%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): unchanged, corr=', best_corr
        case(2)
            call gain_flip_x%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): flip_x, corr=', best_corr
        case(3)
            call gain_flip_y%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): flip_y, corr=', best_corr
        case(4)
            call gain_flip_xy%write(string('best_gainref_by_corr.mrc'))
            write(logfhandle,'(a,1x,f12.6)') '>>> BEST ORIENTATION (most negative corr): flip_xy, corr=', best_corr
        end select

        if( have_prev_analysis )then
            if( best_idx == prev_best_idx .and. abs(best_corr - prev_best_corr) <= CONV_EPS_CORR .and. abs(delta_r - prev_delta_r) <= CONV_EPS_DELTA )then
                conv_streak = conv_streak + 1
            else
                conv_streak = 0
            endif
        endif
        prev_best_idx = best_idx
        prev_best_corr = best_corr
        prev_delta_r = delta_r
        have_prev_analysis = .true.
        write(logfhandle,'(a,1x,i0,1x,a,1x,i0)') '>>> CONVERGENCE streak:', conv_streak, 'required=', CONV_STREAK_REQ

        if( conv_streak >= CONV_STREAK_REQ )then
            converged = .true.
            write(logfhandle,'(a,1x,i0)') '>>> CONVERGED at movie count:', n_movies_processed
            exit
        endif
    endif
end do

if( total_frames < 1 )then
    THROW_HARD('No frames were read from input movies')
endif
if( n_movies_processed < FIRST_ANALYSIS_AT )then
    write(logfhandle,'(a,1x,i0,1x,a)') '>>> Stopped after', n_movies_processed, 'movies; first analysis requires 10 movies'
endif
if( .not. converged .and. n_movies_processed >= n_movies_cap )then
    if( n_movies_cap == MAX_MOVIES )then
        write(logfhandle,'(a,1x,i0)') '>>> Reached hard stop at movies:', MAX_MOVIES
    else
        write(logfhandle,'(a,1x,i0)') '>>> Reached end of provided movies:', n_movies_processed
    endif
endif

if( have_ref_shape )then
    call avg_img%kill()
    call gain_ref%kill()
    call gain_flip_x%kill()
    call gain_flip_y%kill()
    call gain_flip_xy%kill()
    call tmp_avg%kill()
    call tmp_gain%kill()
endif

end program simple_test_extract_dead_pixels

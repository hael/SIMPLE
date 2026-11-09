!@descr: batched FLEX preparation of full-field time-series tracking images
module simple_tseries_prep_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,            only: logfhandle
use simple_defs_fname,      only: MRC_EXT
use simple_error,           only: simple_exception
use simple_fileio,          only: fopen, fclose, fileiochk
use simple_image,           only: image
use simple_imgarr_utils,    only: dealloc_imgarr
use simple_parameters,      only: parameters
use simple_sp_project,      only: sp_project
use simple_string,          only: string
use simple_string_utils,    only: int2str, int2str_pad
use simple_syslib,          only: simple_mkdir, simple_abspath, file_exists
use simple_tseries_denoise,  only: flex_denoise_imgarr
implicit none
private
public :: tseries_prep4tracking_exec
#include "simple_local_flags.inc"

contains

    !> Owns one project read and mic-segment write; original frame inputs are read-only.
    subroutine tseries_prep4tracking_exec( params )
        class(parameters), intent(in) :: params
        real,    parameter :: SMPD_TARGET = 1.0
        integer, parameter :: GRPFREQ = 10, NFRAMESGRP = 5, BATCHSZ = 100
        character(len=*), parameter :: FILT_DIR = 'filtered/'
        type(string), allocatable :: framenames(:)
        type(image),  allocatable :: frames(:), frames_sc(:), denoised_frames_sc(:)
        type(string)     :: fname
        type(image)      :: avgimg, global_avgimg
        type(sp_project) :: spproj
        real    :: smpd, smpd_sc, scale, frame_smpd
        integer :: ldim(3), ldim_sc(3), i, j, n, s, iframe, nframes, last, numlen, nbuf
        integer :: common_divisor, other_divisor, remainder, scaled_divisor
        call spproj%read(params%projfile)
        nframes = spproj%get_nframes()
        if( nframes == 0 ) THROW_HARD('no movie frames to process; tseries_prep4tracking')
        if( nframes /= spproj%os_mic%get_noris() ) THROW_HARD('missing frame rows; tseries_prep4tracking')
        smpd = spproj%os_mic%get(1, 'smpd')
        if( .not. ieee_is_finite(smpd) .or. smpd <= 0. ) THROW_HARD('invalid frame sampling; tseries_prep4tracking')
        ldim = [spproj%os_mic%get_int(1,'xdim'), spproj%os_mic%get_int(1,'ydim'), 1]
        if( any(ldim(1:2) < 2) .or. any(mod(ldim(1:2),2) /= 0) ) THROW_HARD('even 2D frame dimensions required; tseries_prep4tracking')
        allocate(framenames(nframes))
        do i = 1, nframes
            framenames(i) = spproj%os_mic%get_str(i, 'frame')
            if( .not. file_exists(framenames(i)) ) THROW_HARD('missing frame: '//framenames(i)%to_char())
            if( spproj%os_mic%get_int(i,'xdim') /= ldim(1) .or. spproj%os_mic%get_int(i,'ydim') /= ldim(2) )then
                THROW_HARD('inconsistent frame dimensions; tseries_prep4tracking')
            endif
            frame_smpd = spproj%os_mic%get(i, 'smpd')
            if( .not. ieee_is_finite(frame_smpd) ) THROW_HARD('nonfinite frame sampling; tseries_prep4tracking')
            if( abs(frame_smpd - smpd) > 1.e-4*smpd ) THROW_HARD('inconsistent frame sampling; tseries_prep4tracking')
        enddo
        call simple_mkdir(FILT_DIR)
        numlen = len(int2str(nframes))
        scale = min(1.0, smpd / SMPD_TARGET)
        common_divisor = ldim(1)
        other_divisor = ldim(2)
        do while( other_divisor /= 0 )
            remainder = modulo(common_divisor, other_divisor)
            common_divisor = other_divisor
            other_divisor = remainder
        enddo
        ! A common even multiplier preserves aspect ratio and one isotropic sampling value.
        scaled_divisor = 2*max(1,min(common_divisor/2,nint(scale*real(common_divisor)/2.)))
        ldim_sc(1:2) = (ldim(1:2)/common_divisor)*scaled_divisor
        ldim_sc(3) = 1
        scale   = real(ldim_sc(1)) / real(ldim(1))
        smpd_sc = smpd / scale
        write(logfhandle,'(A,F8.3)') '>>> SCALING FACTOR   : ', scale
        write(logfhandle,'(A,2I6)')  '>>> SCALED DIMENSIONS: ', ldim_sc(1:2)
        nbuf = min(BATCHSZ, nframes)
        allocate(frames(nbuf), frames_sc(nbuf))
        do i = 1, nbuf
            call frames(i)%new(ldim, smpd, wthreads=.false.)
            call frames_sc(i)%new(ldim_sc, smpd_sc, wthreads=.false.)
        enddo
        call avgimg%new(ldim_sc, smpd_sc)
        call global_avgimg%new(ldim, smpd)
        do iframe = 1, nframes, BATCHSZ
            last = min(iframe + BATCHSZ - 1, nframes)
            n = last - iframe + 1
            do i = 1, n
                call frames(i)%read(framenames(iframe+i-1))
                call global_avgimg%add(frames(i), w=1./real(nframes))
                call frames(i)%fft
                call frames(i)%clip(frames_sc(i))
                call frames_sc(i)%ifft
                call frames_sc(i)%subtract_background(45.0, mode='default')
            enddo
            call flex_denoise_imgarr(frames_sc(1:n), denoised_frames_sc, params%neigs, params%halfwin, params%eo_pool == 'yes')
            do i = 1, n, GRPFREQ
                call avgimg%zero
                do j = i, min(i + NFRAMESGRP - 1, n)
                    call avgimg%add(denoised_frames_sc(j))
                enddo
                call avgimg%NLmean2D(patch_size=7, search_radius=11, gaussian_patch=.true.)
                call avgimg%bp(45.0, 2.*smpd_sc)
                call avgimg%norm
                s = iframe + i - 1
                fname = FILT_DIR//'frame_'//int2str_pad(s,numlen)//MRC_EXT
                call avgimg%write(fname)
                ! Retain these return-handoff fields for Nanotation's exported image paths.
                ! Contract: doc/policies/tseries_tracking_external_gui_contract.md
                call spproj%os_mic%set(s, 'track_fname', simple_abspath(fname))
                call spproj%os_mic%set(s, 'fromt', s)
                call spproj%os_mic%set(s, 'tot', min(s + GRPFREQ - 1, nframes))
            enddo
            call dealloc_imgarr(denoised_frames_sc)
        enddo
        call write_outliers
        call spproj%write_segment_inside('mic', params%projfile)
        call spproj%kill
        call dealloc_imgarr(frames)
        call dealloc_imgarr(frames_sc)
        call avgimg%kill
        call global_avgimg%kill
        call framenames(:)%kill
        deallocate(framenames)
    contains

        !> Reports global-average pixels beyond six standard deviations.
        subroutine write_outliers
            real, parameter :: OUTLIER_SIGMA = 6.0
            real    :: ave, sdev, maxv, minv, v, dead_t, hot_t
            integer :: ix, iy, ndead, nhot, funit, io_stat, count_out
            call global_avgimg%write(string('global_avgimg.mrc'))
            call global_avgimg%stats(ave, sdev, maxv, minv)
            dead_t = ave - OUTLIER_SIGMA*sdev
            hot_t  = ave + OUTLIER_SIGMA*sdev
            ndead = 0
            nhot = 0
            count_out = 0
            call fopen(funit, string('outliers.txt'), 'replace', 'write', iostat=io_stat, form='formatted')
            call fileiochk('cannot write outliers.txt; tseries_prep4tracking', io_stat)
            do iy = 1, ldim(2)
                do ix = 1, ldim(1)
                    v = global_avgimg%get([ix,iy,1])
                    if( v < dead_t )then
                        ndead = ndead + 1
                    elseif( v > hot_t )then
                        nhot = nhot + 1
                    else
                        cycle
                    endif
                    count_out = count_out + 1
                    write(funit,'(3I6)') count_out, ix, iy
                enddo
            enddo
            call fclose(funit)
            write(logfhandle,'(A,I6)') '>>> # DEAD PIXELS: ', ndead
            write(logfhandle,'(A,I6)') '>>> # HOT PIXELS : ', nhot
        end subroutine write_outliers

    end subroutine tseries_prep4tracking_exec

end module simple_tseries_prep_strategy

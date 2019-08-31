! time series tracker intended for movies of nanoparticles spinning in solution

module simple_tseries_tracker
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_tvfilter
implicit none

public :: init_tracker, track_particle, tracker_get_nnn
public :: tracker_get_nnfname, tracker_get_stkname, kill_tracker
private
#include "simple_local_flags.inc"

real,    parameter :: EPS     = 0.1
real,    parameter :: SPECW   = 0.0001
logical, parameter :: DOPRINT = .true.
integer, parameter :: CENRATE = 15, NNN = 8

integer,                   allocatable :: particle_locations(:,:)
character(len=LONGSTRLEN), pointer     :: framenames(:)
character(len=:),          allocatable :: dir, fbody
type(image)               :: frame_img      ! individual frame image
type(image)               :: frame_avg      ! average over time window
type(image)               :: frame_avg_filt ! filered version of time window average actually used for tracking
type(image)               :: reference, tmp_img, ptcl_target
type(image)               :: neigh_imgs(NNN), backgr_img, pspec
type(tvfilter)            :: tv
character(len=LONGSTRLEN) :: neighstknames(NNN), stkname
integer                   :: ldim(3), nframes, nx, ny, neigh_cnt
real                      :: sxx
logical                   :: l_neg

contains

    subroutine init_tracker( boxcoord, fnames, dir_in, fbody_in )
        integer,                           intent(in) :: boxcoord(2)
        character(len=LONGSTRLEN), target, intent(in) :: fnames(:)
        character(len=*),                  intent(in)  :: dir_in, fbody_in
        character(len=:), allocatable :: fname
        integer :: n, i
        ! set constants
        dir    = trim(dir_in)
        fbody  = trim(fbody_in)
        l_neg  = .false.
        if( trim(params_glob%neg) .eq. 'yes' ) l_neg = .true.
        select case(trim(params_glob%filter))
            case('no','nlmean')
                ! all good
            case('tv')
                call tv%new()
            case DEFAULT
                THROW_HARD('Unsupported filter in init_tracker!')
        end select
        ! names & dimensions
        nframes = size(fnames)
        framenames => fnames(:)
        call find_ldim_nptcls(framenames(1),ldim,n)
        if( n == 1 .and. ldim(3) == 1 )then
            ! all ok
        else
            write(logfhandle,*) 'ldim(3): ', ldim(3)
            write(logfhandle,*) 'nframes: ', n
            THROW_HARD('init_tracker; assumes one frame per file')
        endif
        nx = ldim(1) - params_glob%box
        ny = ldim(2) - params_glob%box
        ! construct
        allocate(particle_locations(nframes,2), source=0)
        call frame_img%new(ldim, params_glob%smpd)
        call frame_avg%new(ldim, params_glob%smpd)
        call frame_avg_filt%new(ldim, params_glob%smpd)
        call tmp_img%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call reference%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call ptcl_target%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call backgr_img%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call pspec%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        do i=1,NNN
            call neigh_imgs(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd)
            fname = trim(dir)//'/'//trim(fbody)//'_background_nn'//int2str(i)//'.mrcs'
            neighstknames(i) = trim(fname)
        end do
        particle_locations(:,1) = boxcoord(1)
        particle_locations(:,2) = boxcoord(2)
    end subroutine init_tracker

    subroutine track_particle
        character(len=:), allocatable :: fname
        integer :: funit, io_stat, iframe, xind, yind, i
        integer :: pos(2), pos_refined(2)
        logical :: outside
        ! extract first reference
        call update_frame_images(1)
        pos = particle_locations(1,:)
        call update_reference(1, pos)
        call write_background_images_and_update_pspec(1, pos)
        ! init pspec
        call pspec%zero_and_unflag_ft
        ! init neigh counter
        neigh_cnt = 0
        ! track
        write(logfhandle,'(a)') ">>> TRACKING PARTICLE"
        do iframe=2,nframes
            call progress(iframe,nframes)
            ! update frame images & refine position
            call update_frame_images(iframe)
            call refine_position( pos, pos_refined )
            ! update position & reference
            pos = pos_refined
            ! set position and propagate fwd
            particle_locations(iframe:,1) = pos(1)
            particle_locations(iframe:,2) = pos(2)
            call update_reference(iframe, pos)
            call write_background_images_and_update_pspec(iframe, pos)
        end do
        ! write box file and tracked particles
        fname = trim(dir)//'/'//trim(fbody)//'.box'
        call fopen(funit, status='REPLACE', action='WRITE', file=fname,iostat=io_stat)
        call fileiochk("tseries tracker ; write_tracked_series ", io_stat)
        stkname = trim(dir)//'/'//trim(fbody)//'.mrc'
        do iframe=1,nframes
            xind = particle_locations(iframe,1)
            yind = particle_locations(iframe,2)
            write(funit,'(I7,I7,I7,I7,I7)') xind, yind, params_glob%box, params_glob%box, -3
            call frame_img%read(framenames(iframe),1)
            call frame_img%window_slim([xind,yind,1], params_glob%box, ptcl_target, outside)
            if( l_neg ) call ptcl_target%neg()
            call ptcl_target%write(stkname, iframe)
        end do
        call fclose(funit, errmsg="tseries tracker ; write_tracked_series end")
        ! average and write power spectrum for CTF estimation
        call pspec%div(SPECW * real(neigh_cnt))
        fname = trim(dir)//'/'//trim(fbody)//'_pspec4ctf_estimation.mrc'
        call pspec%write(fname, 1)
    end subroutine track_particle

    subroutine update_reference( iframe, pos )
        integer, intent(in) :: iframe, pos(2)
        real    :: xyz(3)
        logical :: outside
        call frame_avg_filt%window_slim(pos, params_glob%box, tmp_img, outside)
        call tmp_img%prenorm4real_corr(sxx)
        if( iframe == 1 )then
            reference = tmp_img
        else
            call reference%mul(1.0 - EPS)
            call reference%add(tmp_img, EPS)
        endif
        if( mod(iframe,CENRATE) == 0 )then
            ! center the reference
            if( l_neg ) call reference%neg()
            xyz = reference%calc_shiftcen(params_glob%cenlp)
            if( l_neg ) call reference%neg()
        endif
    end subroutine update_reference

    subroutine write_background_images_and_update_pspec( iframe, pos )
        integer, intent(in) :: iframe, pos(2)
        integer :: neigh(NNN,2), i
        logical :: outside
        call identify_neighbours
        do i=1,NNN
            call frame_img%window_slim(neigh(i,:), params_glob%box, backgr_img, outside)
            if( outside ) call backgr_img%zero
            if( l_neg )   call backgr_img%neg()
            call backgr_img%write(neighstknames(i), iframe)
            if( .not. outside )then
                neigh_cnt = neigh_cnt + 1
                if( l_neg ) call backgr_img%neg()
                call backgr_img%norm()
                call backgr_img%zero_edgeavg
                call backgr_img%fft()
                call backgr_img%ft2img('sqrt', tmp_img)
                call pspec%add(tmp_img, SPECW)
                call backgr_img%zero_and_unflag_ft
                call tmp_img%zero_and_unflag_ft
            endif
        end do

        contains

            subroutine identify_neighbours
                ! neigh east
                neigh(1,1) = pos(1) + params_glob%box
                neigh(1,2) = pos(2)
                ! neigh west
                neigh(2,1) = pos(1) - params_glob%box
                neigh(2,2) = pos(2)
                ! neigh north
                neigh(3,1) = pos(1)
                neigh(3,2) = pos(2) + params_glob%box
                ! neigh south
                neigh(4,1) = pos(1)
                neigh(4,2) = pos(2) - params_glob%box
                ! neigh north/east
                neigh(5,1) = pos(1) + params_glob%box
                neigh(5,2) = pos(2) + params_glob%box
                ! neigh north/west
                neigh(6,1) = pos(1) - params_glob%box
                neigh(6,2) = pos(2) + params_glob%box
                ! neigh south/east
                neigh(7,1) = pos(1) + params_glob%box
                neigh(7,2) = pos(2) - params_glob%box
                ! neigh south/west
                neigh(8,1) = pos(1) - params_glob%box
                neigh(8,2) = pos(2) - params_glob%box
            end subroutine identify_neighbours

    end subroutine write_background_images_and_update_pspec

    subroutine update_frame_images( iframe )
        integer, intent(in) :: iframe
        integer :: fromto(2), i
        real    :: w
        fromto(1) = iframe - params_glob%nframesgrp/2
        fromto(2) = iframe + params_glob%nframesgrp/2 - 1
        do while(fromto(1) < 1)
            fromto = fromto + 1
        end do
        do while(fromto(2) > nframes)
            fromto = fromto - 1
        end do
        ! create frame_avg
        call frame_avg%zero
        w = 1. / real(params_glob%nframesgrp)
        do i=fromto(1),fromto(2)
            call frame_img%read(framenames(i),1)
            call frame_avg%add(frame_img, w)
        end do
        call frame_avg%div(real(params_glob%nframesgrp))
        ! frame_img should correspond to the current frame
        call frame_img%read(framenames(iframe),1)
        ! make filtered average
        call frame_avg_filt%copy(frame_avg)
        select case(trim(params_glob%filter))
            case('tv')
                call frame_avg_filt%fft()
                call tv%apply_filter(frame_avg_filt, 5.)
                ! always low-pass filter
                call frame_avg_filt%bp(0., params_glob%lp)
            case('nlmean')
                call frame_avg_filt%nlmean
                call frame_avg_filt%fft()
                ! always low-pass filter
                call frame_avg_filt%bp(0., params_glob%lp)
            case DEFAULT
                ! defaults to just low-pass filtering
                call frame_avg_filt%fft()
                call frame_avg_filt%bp(0., params_glob%lp)
        end select
        call frame_avg_filt%ifft()
    end subroutine update_frame_images

    subroutine refine_position( pos, pos_refined )
        integer, intent(in)  :: pos(2)
        integer, intent(out) :: pos_refined(2)
        integer :: xind, yind, xrange(2), yrange(2)
        real    :: corr, target_corr
        logical :: outside
        ! set srch range
        xrange(1) = max(0,  pos(1) - params_glob%offset)
        xrange(2) = min(nx, pos(1) + params_glob%offset)
        yrange(1) = max(0,  pos(2) - params_glob%offset)
        yrange(2) = min(ny, pos(2) + params_glob%offset)
        ! extract image, correlate, find peak
        corr = -1
        do xind=xrange(1),xrange(2)
            do yind=yrange(1),yrange(2)
                call frame_avg_filt%window_slim([xind,yind,1], params_glob%box, ptcl_target, outside)
                target_corr = reference%real_corr_prenorm(ptcl_target, sxx)
                if( target_corr > corr )then
                    pos_refined = [xind,yind]
                    corr = target_corr
                endif
            end do
        end do
    end subroutine refine_position

    integer function tracker_get_nnn()
        tracker_get_nnn = NNN
    end function tracker_get_nnn

    character(len=LONGSTRLEN) function tracker_get_nnfname(i)
        integer, intent(in) :: i
        if( i<1 .or. i>NNN )THROW_HARD('neighbour index out of range')
        tracker_get_nnfname = trim(neighstknames(i))
    end function tracker_get_nnfname

    character(len=LONGSTRLEN) function tracker_get_stkname()
        tracker_get_stkname = trim(stkname)
    end function tracker_get_stkname

    subroutine kill_tracker
        integer :: i
        deallocate(particle_locations)
        framenames => null()
        do i=1,NNN
            call neigh_imgs(i)%kill
        end do
        call frame_img%kill
        call frame_avg%kill
        call frame_avg_filt%kill
        call reference%kill
        call tmp_img%kill
        call ptcl_target%kill
        call backgr_img%kill
        call pspec%kill
        call tv%kill
    end subroutine kill_tracker

end module simple_tseries_tracker

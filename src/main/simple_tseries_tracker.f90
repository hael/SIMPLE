! time series tracker intended for movies of nanoparticles spinning in solution

module simple_tseries_tracker
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: init_tracker, track_particle, write_tracked_series, kill_tracker
private
#include "simple_local_flags.inc"

integer,                   allocatable :: particle_locations(:,:)
character(len=LONGSTRLEN), allocatable :: framenames(:)
real,                      parameter   :: EPS=0.1
logical,                   parameter   :: DOPRINT=.true.
integer,                   parameter   :: CENRATE=15, NNN=8
type(image)           :: frame_img, reference, tmp_img, ptcl_target
type(image)           :: neigh_imgs_mean(NNN), diff_img
integer               :: ldim(3), nframes, box, nx, ny, offset
real                  :: smpd, sxx, lp, cenlp, sumw
character(len=3)      :: neg
logical               :: l_neg

contains

    subroutine init_tracker(  boxcoord )
        integer, intent(in) :: boxcoord(2)
        integer :: n, i
        ! set constants
        box    = params_glob%box
        offset = params_glob%offset
        smpd   = params_glob%smpd
        lp     = params_glob%lp
        cenlp  = params_glob%cenlp
        neg    = params_glob%neg
        l_neg  = .false.
        if( params_glob%neg .eq. 'yes' ) l_neg = .true.
        call read_filetable(params_glob%filetab, framenames)
        nframes = size(framenames)
        call find_ldim_nptcls(framenames(1),ldim,n)
        if( n == 1 .and. ldim(3) == 1 )then
            ! all ok
        else
            write(logfhandle,*) 'ldim(3): ', ldim(3)
            write(logfhandle,*) 'nframes: ', n
            THROW_HARD('init_tracker; assumes one frame per file')
        endif
        nx = ldim(1) - box
        ny = ldim(2) - box
        ! construct
        allocate(particle_locations(nframes,2), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: simple_tseries_tracker :: init_tracker",alloc_stat)
        particle_locations = 0
        call frame_img%new(ldim, smpd)
        call tmp_img%new([box,box,1], smpd)
        call reference%new([box,box,1], smpd)
        call ptcl_target%new([box,box,1], smpd)
        call diff_img%new([box,box,1], smpd)
        do i=1,NNN
            call neigh_imgs_mean(i)%new([box,box,1], smpd)
        end do
        particle_locations(:,1) = boxcoord(1)
        particle_locations(:,2) = boxcoord(2)
    end subroutine init_tracker

    subroutine track_particle
        integer :: pos(2), pos_refined(2), iframe
        ! extract first reference
        call update_frame(1)
        pos = particle_locations(1,:)
        call update_reference(1, pos)
        sumw = 1.0
        call update_background_images(1, pos)
        ! track
        write(logfhandle,'(a)') ">>> TRACKING PARTICLE"
        do iframe=2,nframes
            sumw = sumw + 1.0
            call progress(iframe,nframes)
            ! update frame & refine position
            call update_frame(iframe)
            call refine_position( pos, pos_refined )
            ! update position & reference
            pos = pos_refined
            ! set position and propagate fwd
            particle_locations(iframe:,1) = pos(1)
            particle_locations(iframe:,2) = pos(2)
            call update_reference(iframe, pos)
            call update_background_images(iframe, pos)
        end do
    end subroutine track_particle

    subroutine write_tracked_series( fbody )
        character(len=*), intent(in) :: fbody
        integer :: funit, io_stat, iframe, xind, yind, i
        logical :: outside
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(fbody)//'.box',iostat=io_stat)
        call fileiochk("tseries tracker ; write_tracked_series ", io_stat)
        do iframe=1,nframes
            xind = particle_locations(iframe,1)
            yind = particle_locations(iframe,2)
            write(funit,'(I7,I7,I7,I7,I7)') xind, yind, box, box, -3
            call frame_img%read(framenames(iframe),1)
            call frame_img%window_slim([xind,yind,1], box, reference, outside)
            if( l_neg ) call reference%neg()
            call reference%write(trim(fbody)//'.mrc', iframe)
        end do
        call fclose(funit, errmsg="tseries tracker ; write_tracked_series end")
        do i=1,NNN
            if( l_neg ) call neigh_imgs_mean(i)%neg()
            call neigh_imgs_mean(i)%write(trim(fbody)//'_background_nn'//int2str(i)//'.mrc')
        end do
    end subroutine write_tracked_series

    subroutine update_reference( iframe, pos )
        integer, intent(in) :: iframe, pos(2)
        real    :: xyz(3)
        logical :: outside
        call frame_img%window_slim(pos, box, tmp_img, outside)
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
            xyz = reference%calc_shiftcen(cenlp)
            if( l_neg ) call reference%neg()
        endif
    end subroutine update_reference

    subroutine update_background_images( iframe, pos )
        integer, intent(in) :: iframe, pos(2)
        integer :: neigh(NNN,2), i
        logical :: outside
        call identify_neighbours
        do i=1,NNN
            call frame_img%window_slim(neigh(i,:), box, diff_img, outside)
            if( .not. outside )then
                call diff_img%subtr(neigh_imgs_mean(i))
                call diff_img%div(sumw)
                call neigh_imgs_mean(i)%add(diff_img)
            endif
        end do

        contains

            subroutine identify_neighbours
                ! neigh east
                neigh(1,1) = pos(1) + box
                neigh(1,2) = pos(2)
                ! neigh west
                neigh(2,1) = pos(1) - box
                neigh(2,2) = pos(2)
                ! neigh north
                neigh(3,1) = pos(1)
                neigh(3,2) = pos(2) + box
                ! neigh south
                neigh(4,1) = pos(1)
                neigh(4,2) = pos(2) - box
                ! neigh north/east
                neigh(5,1) = pos(1) + box
                neigh(5,2) = pos(2) + box
                ! neigh north/west
                neigh(6,1) = pos(1) - box
                neigh(6,2) = pos(2) + box
                ! neigh south/east
                neigh(7,1) = pos(1) + box
                neigh(7,2) = pos(2) - box
                ! neigh south/west
                neigh(8,1) = pos(1) - box
                neigh(8,2) = pos(2) - box
            end subroutine identify_neighbours

    end subroutine update_background_images

    subroutine update_frame( iframe )
        integer, intent(in) :: iframe
        call frame_img%read(framenames(iframe),1)
        call frame_img%fft()
        call frame_img%bp(0., lp)
        call frame_img%ifft()
    end subroutine update_frame

    subroutine refine_position( pos, pos_refined )
        integer, intent(in)  :: pos(2)
        integer, intent(out) :: pos_refined(2)
        integer     :: xind, yind, xrange(2), yrange(2)
        real        :: corr, target_corr
        logical     :: outside
        ! set srch range
        xrange(1) = max(0,  pos(1) - offset)
        xrange(2) = min(nx, pos(1) + offset)
        yrange(1) = max(0,  pos(2) - offset)
        yrange(2) = min(ny, pos(2) + offset)
        ! extract image, correlate, find peak
        corr = -1
        do xind=xrange(1),xrange(2)
            do yind=yrange(1),yrange(2)
                call frame_img%window_slim([xind,yind,1], box, ptcl_target, outside)
                target_corr = reference%real_corr_prenorm(ptcl_target, sxx)
                if( target_corr > corr )then
                    pos_refined = [xind,yind]
                    corr = target_corr
                endif
            end do
        end do
    end subroutine refine_position

    subroutine kill_tracker
        deallocate(particle_locations, framenames, stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("simple_tseries_tracker::kill_tracker dealloc",alloc_stat)
        call frame_img%kill
        call reference%kill
        call tmp_img%kill
        call ptcl_target%kill
    end subroutine kill_tracker

end module simple_tseries_tracker

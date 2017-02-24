module simple_tseries_tracker
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_image,        only: image
use simple_jiffys,       only: alloc_err, find_ldim_nptcls
use simple_strings,      only: int2str
use simple_filehandling, only: read_filetable
implicit none

public :: init_tracker, track_particle, write_tracked_series, kill_tracker
private

integer,               allocatable :: particle_locations(:,:)
character(len=STDLEN), allocatable :: framenames(:)
logical, parameter     :: DEBUG=.true., DOPRINT=.true.
type(image)            :: frame_img, reference
integer                :: ldim(3), nframes, box, nx, ny, offset
real                   :: smpd, sxx, lp

contains

    subroutine init_tracker( filetabname, boxcoord, box_in, offset_in, smpd_in, lp_in  )
        character(len=*), intent(in) :: filetabname
        integer,          intent(in) :: boxcoord(2), box_in, offset_in
        real,             intent(in) :: smpd_in, lp_in
        integer :: alloc_stat, n
        ! set constants
        box    = box_in
        offset = offset_in
        smpd   = smpd_in
        lp     = lp_in
        call read_filetable(filetabname, framenames)
        nframes = size(framenames)  
        call find_ldim_nptcls(framenames(1),ldim,n)
        if( n == 1 .and. ldim(3) == 1 )then
            ! all ok
        else
            write(*,*) 'ldim(3): ', ldim(3)
            write(*,*) 'nframes: ', n
            stop 'simple_tseries_tracker :: init_tracker; assumes one frame per file' 
        endif
        nx = ldim(1) - box
        ny = ldim(2) - box
        ! construct
        allocate(particle_locations(nframes,2), stat=alloc_stat)
        call alloc_err("In: simple_tseries_tracker :: init_tracker", alloc_stat)
        particle_locations = 0
        call frame_img%new(ldim, smpd)
        call reference%new([box,box,1], smpd)
        particle_locations(1,1) = boxcoord(1)
        particle_locations(1,2) = boxcoord(2)
    end subroutine init_tracker

    subroutine track_particle
        use simple_jiffys, only: progress
        integer :: pos(2), pos_refined(2), iframe
        ! extract first reference
        call update_frame(1)
        pos = particle_locations(1,:)
        call update_reference(pos)
        ! track
        write(*,'(a)') ">>> TRACKING PARTICLE"
        do iframe=2,nframes
            call progress(iframe,nframes)
            ! update frame & refine position
            call update_frame(iframe)
            call refine_position( pos, pos_refined )
            ! update position & reference
            pos = pos_refined
            particle_locations(iframe,:) = pos
            call update_reference(pos)
        end do
    end subroutine track_particle

    subroutine write_tracked_series( fbody, neg )
        use simple_filehandling, only: get_fileunit
        character(len=*), intent(in) :: fbody
        character(len=*), intent(in) :: neg
        integer :: funit, iframe, xind, yind
        funit = get_fileunit()
        open(unit=funit, status='REPLACE', action='WRITE', file=trim(fbody)//'.box')
        do iframe=1,nframes
            xind = particle_locations(iframe,1)
            yind = particle_locations(iframe,2)
            write(funit,'(I7,I7,I7,I7,I7)') xind, yind, box, box, -3
            call frame_img%read(framenames(iframe),1)
            call frame_img%window_slim([xind,yind,1], box, reference)
            if( neg .eq. 'yes' ) call reference%neg
            call reference%norm
            call reference%write(trim(fbody)//'.mrc', iframe)
        end do
        close(funit)
    end subroutine write_tracked_series

    subroutine update_reference( pos )
        integer, intent(in) :: pos(2)
        call frame_img%window_slim(pos, box, reference)
        call reference%prenorm4real_corr(sxx)
    end subroutine update_reference

    subroutine update_frame( iframe )
        integer, intent(in) :: iframe
        call frame_img%read(framenames(iframe),1)
        call frame_img%fwd_ft
        call frame_img%bp(0., lp)
        call frame_img%bwd_ft
    end subroutine update_frame

    subroutine refine_position( pos, pos_refined )
        integer, intent(in)  :: pos(2)
        integer, intent(out) :: pos_refined(2)
        type(image), allocatable :: target_imgs(:,:)
        real,        allocatable :: target_corrs(:,:)
        integer :: xind, yind, xrange(2), yrange(2)
        real    :: corr
        ! set srch range
        xrange(1) = max(0,  pos(1) - offset)
        xrange(2) = min(nx, pos(1) + offset)
        yrange(1) = max(0,  pos(2) - offset)
        yrange(2) = min(ny, pos(2) + offset)
        ! allocate
        allocate(target_imgs(xrange(1):xrange(2),yrange(1):yrange(2)),&
                 target_corrs(xrange(1):xrange(2),yrange(1):yrange(2)))
        ! extract image matrix
        do xind=xrange(1),xrange(2)
            do yind=yrange(1),yrange(2)
                call target_imgs(xind,yind)%new([box,box,1], smpd)
                call frame_img%window_slim([xind,yind,1], box, target_imgs(xind,yind))
            end do
        end do
        ! correlate
        do xind=xrange(1),xrange(2)
            do yind=yrange(1),yrange(2)
                target_corrs(xind,yind) =&
                &reference%real_corr_prenorm(target_imgs(xind,yind), sxx)
            end do
        end do
        ! find peak
        pos_refined = pos
        corr        = target_corrs(pos(1),pos(2))
        do xind=xrange(1),xrange(2)
            do yind=yrange(1),yrange(2)
                call target_imgs(xind,yind)%kill
                if( target_corrs(xind,yind) > corr )then
                    pos_refined = [xind,yind]
                    corr = target_corrs(xind,yind)
                endif
            end do
        end do
        deallocate(target_imgs, target_corrs)
    end subroutine refine_position

    subroutine kill_tracker
        deallocate(particle_locations, framenames)
        call frame_img%kill
        call reference%kill
    end subroutine kill_tracker

end module simple_tseries_tracker
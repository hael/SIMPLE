! time series tracker intended for movies of nanoparticles spinning in solution
module simple_tseries_tracker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_tvfilter
implicit none

public :: init_tracker, track_particle, tracker_get_nnn, detect_outliers
public :: tracker_get_nnfname, tracker_get_stkname, kill_tracker
private
#include "simple_local_flags.inc"

real,    parameter :: SPECW    = 0.0001
real,    parameter :: TVLAMBDA = 10.
real,    parameter :: BFACTOR  = 5.
integer, parameter :: NNN      = 8

type(image),               allocatable :: ptcls(:), ptcls_saved(:)
type(tvfilter),            allocatable :: tv(:)
real,                      allocatable :: particle_locations(:,:)
character(len=LONGSTRLEN), pointer     :: framenames(:)
character(len=:),          allocatable :: dir, fbody
type(image)               :: frame_img      ! individual frame image
type(image)               :: frame_avg      ! average over time window
type(image)               :: reference, ptcl_target
type(image)               :: neigh_imgs(NNN), backgr_imgs(NNN), tmp_imgs(NNN), pspec
character(len=LONGSTRLEN) :: neighstknames(NNN), stkname
integer                   :: ldim(3), nframes, neigh_cnt
logical                   :: l_neg

integer :: cnt4debug = 0

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
                allocate(tv(params_glob%nframesgrp))
                do i=1,params_glob%nframesgrp
                    call tv(i)%new()
                enddo
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
        ! construct
        allocate(particle_locations(nframes,2), source=0.)
        call frame_img%new(ldim, params_glob%smpd)
        call frame_avg%new(ldim, params_glob%smpd)
        call reference%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call ptcl_target%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        call pspec%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        do i=1,NNN
            call neigh_imgs(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd,wthreads=.false.)
            call backgr_imgs(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd,wthreads=.false.)
            call tmp_imgs(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd,wthreads=.false.)
            fname = trim(dir)//'/'//trim(fbody)//'_background_nn'//int2str(i)//'.mrcs'
            neighstknames(i) = trim(fname)
        end do
        allocate(ptcls(params_glob%nframesgrp),ptcls_saved(params_glob%nframesgrp))
        do i=1,params_glob%nframesgrp
            call ptcls(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd, wthreads=.false.)
            call ptcls_saved(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd, wthreads=.false.)
        end do
        particle_locations(:,1) = real(boxcoord(1))
        particle_locations(:,2) = real(boxcoord(2))
    end subroutine init_tracker

    subroutine detect_outliers
        logical, allocatable :: outliers(:,:)
        integer              :: deadhot(2), iframe, cnt, ncured
        write(logfhandle,'(a)') ">>> DETECTING OUTLIERS"
        call frame_avg%zero_and_unflag_ft
        cnt = 0
        do iframe = 1,nframes,10
            cnt = cnt + 1
            call frame_img%read(framenames(iframe),1)
            call frame_avg%add(frame_img,w=0.01)
        enddo
        call frame_avg%div(0.01*real(cnt))
        call frame_avg%cure_outliers(ncured, 6., deadhot, outliers)
        write(logfhandle,'(a,1x,i7)') '>>> # DEAD PIXELS:', deadhot(1)
        write(logfhandle,'(a,1x,i7)') '>>> # HOT  PIXELS:', deadhot(2)
    end subroutine detect_outliers

    subroutine track_particle( fname_forctf )
        use simple_motion_align_nano, only: motion_align_nano
        character(len=:), allocatable, intent(inout) :: fname_forctf
        type(motion_align_nano)       :: aligner
        character(len=:), allocatable :: fname
        real,             allocatable :: ini_shifts(:,:), opt_shifts(:,:)
        integer :: first_pos(3), funit, io_stat, iframe, xind,yind, i, last_frame, cnt, nrange
        real    :: xyz(3), pos(2)
        logical :: outside
        ! init neigh counter
        neigh_cnt = 0
        call frame_avg%zero_and_unflag_ft
        call reference%zero_and_flag_ft
        write(logfhandle,'(a)') ">>> TRACKING PARTICLES"
        cnt4debug = 0
        do iframe = 1,nframes,params_glob%nframesgrp
            cnt4debug = cnt4debug + 1
            ! read frames & extract particles
            last_frame = min(iframe+params_glob%nframesgrp-1,nframes)
            nrange     = last_frame-iframe+1
            allocate(ini_shifts(nrange,2),source=0.)
            first_pos = nint([particle_locations(iframe,1),particle_locations(iframe,2),1.])
            cnt = 0
            do i=iframe,last_frame
                cnt = cnt + 1
                call ptcls(cnt)%zero_and_unflag_ft
                call frame_img%read(framenames(i),1)
                call frame_img%window_slim(first_pos, params_glob%box, ptcls(cnt), outside)
                if( cnt==1 ) call frame_avg%add(frame_img,w=0.01)
            enddo
            ! prep images: norm, mask, filter, FFT
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1,nrange
                ini_shifts(i,:) = particle_locations(iframe,:) - real(first_pos(1:2))
                call ptcls(i)%norm
                select case(trim(params_glob%filter))
                    case('nlmean')
                        call ptcls(i)%nlmean
                    case('tv')
                        call tv(i)%apply_filter(ptcls(i), TVLAMBDA)
                    case DEFAULT
                        ! nothing to do
                end select
                call ptcls_saved(i)%copy(ptcls(i))
                call ptcls_saved(i)%fft
                call ptcls(i)%mask(real(params_glob%box/2)-3.,'soft')
                call ptcls(i)%fft
            enddo
            !$omp end parallel do
            ! align
            call aligner%new(ptcls)
            call aligner%set_reslims(params_glob%hp, params_glob%lp)
            call aligner%set_bfactor(BFACTOR)
            call aligner%set_trs(real(params_glob%offset))
            if( iframe == 1 )then
                call aligner%align
            else
                ! providing previous reference to mitigate drift
                call aligner%align(reference, ini_shifts=ini_shifts)
            endif
            call aligner%get_opt_shifts(opt_shifts)
            ! generate reference
            call reference%zero_and_flag_ft
            do i = 1,nrange
                call ptcls_saved(i)%shift([-opt_shifts(i,1),-opt_shifts(i,2),0.])
                call reference%add(ptcls_saved(i),w=1./real(nrange))
            enddo
            call reference%ifft
            call reference%mask(real(params_glob%box/2)-3.,'soft')
            ! call reference%write('reference.mrc',cnt4debug)
            if( l_neg ) call reference%neg
            xyz = center_reference()
            if( l_neg ) call reference%neg
            select case(trim(params_glob%filter))
                case('nlmean')
                    call reference%nlmean
                    call reference%fft
                case('tv')
                    call reference%fft
                    call tv(1)%apply_filter(reference, TVLAMBDA)
                case DEFAULT
                    call reference%fft
            end select
            call reference%shift(xyz)
            ! updates shifts
            pos = particle_locations(iframe,:)
            cnt = 0
            do i=iframe,last_frame
                cnt = cnt + 1
                particle_locations(i,1) = pos(1) - opt_shifts(cnt,1)
                particle_locations(i,2) = pos(2) - opt_shifts(cnt,2)
            enddo
            particle_locations(last_frame+1:,1) = particle_locations(last_frame,1)
            particle_locations(last_frame+1:,2) = particle_locations(last_frame,2)
            ! cleanup
            call aligner%kill
            deallocate(ini_shifts)
        enddo
        write(logfhandle,'(a)') ">>> WRITING PARTICLES, NEIGHBOURS & SPECTRUM"
        ! Second pass write box file, tracked particles, neighbours & update spectrum
        call pspec%zero_and_unflag_ft
        call ptcl_target%zero_and_unflag_ft
        fname = trim(dir)//'/'//trim(fbody)//'.box'
        call fopen(funit, status='REPLACE', action='WRITE', file=fname,iostat=io_stat)
        call fileiochk("tseries tracker ; write_tracked_series ", io_stat)
        stkname = trim(dir)//'/'//trim(fbody)//'.mrc'
        do iframe=1,nframes
            xind = nint(particle_locations(iframe,1))
            yind = nint(particle_locations(iframe,2))
            ! read
            call frame_img%read(framenames(iframe),1)
            ! box
            write(funit,'(I7,I7,I7,I7,I7)') xind, yind, params_glob%box, params_glob%box, -3
            ! particle
            call frame_img%window_slim([xind,yind,1], params_glob%box, ptcl_target, outside)
            call ptcl_target%norm
            if( l_neg ) call ptcl_target%neg()
            call ptcl_target%write(stkname, iframe)
            ! neighbors & spectrum
            call write_background_images_and_update_pspec(iframe, [xind,yind])
        end do
        call fclose(funit, errmsg="tseries tracker ; write_tracked_series end")
        ! average and write power spectrum for CTF estimation
        call pspec%div(SPECW * real(neigh_cnt))
        fname_forctf = trim(dir)//'/'//trim(fbody)//'_pspec4ctf_estimation.mrc'
        call pspec%dampen_pspec_central_cross
        call pspec%write(fname_forctf, 1)
        ! trajectory
        call write_trajectory
    end subroutine track_particle

    subroutine write_trajectory
        integer :: xind,yind,i,j,iframe,sz
        call frame_avg%norm
        do iframe = nframes,1,-1
            xind = nint((particle_locations(iframe,1)+1.) + real(params_glob%box)/2.+1.)
            yind = nint((particle_locations(iframe,2)+1.) + real(params_glob%box)/2.+1.)
            sz   = 1
            if( iframe==1 ) sz=3
            do j = yind-sz,yind+sz
                do i = xind-sz,xind+sz
                    call frame_avg%set([i,j,1],-3.)
                enddo
            enddo
            if( iframe==1 )then
                do j = yind-1,yind+1
                    do i = xind-1,xind+1
                        call frame_avg%set([i,j,1],3.)
                    enddo
                enddo
            endif
        enddo
        call frame_avg%write(trim(dir)//'/'//trim(fbody)//'_trajectory.mrc')
    end subroutine write_trajectory

    subroutine write_background_images_and_update_pspec( iframe, pos )
        integer, intent(in) :: iframe, pos(2)
        integer :: neigh(NNN,2), i
        logical :: outside(NNN)
        call identify_neighbours
        !$omp parallel default(shared) private(i) proc_bind(close)
        !$omp do schedule(static)
        do i=1,NNN
            call frame_img%window_slim(neigh(i,:), params_glob%box, backgr_imgs(i), outside(i))
            if( outside(i) ) call backgr_imgs(i)%zero
            if( l_neg )      call backgr_imgs(i)%neg()
        enddo
        !$omp end do nowait
        !$omp single
        do i=1,NNN
            call backgr_imgs(i)%write(neighstknames(i), iframe)
        enddo
        neigh_cnt = neigh_cnt + count(.not.outside)
        !$omp end single
        !$omp do schedule(static)
        do i=1,NNN
            if( .not. outside(i) )then
                if( l_neg ) call backgr_imgs(i)%neg()
                call backgr_imgs(i)%norm()
                call backgr_imgs(i)%zero_edgeavg
                call backgr_imgs(i)%fft()
                call backgr_imgs(i)%ft2img('sqrt', tmp_imgs(i))
                call backgr_imgs(i)%zero_and_unflag_ft
            endif
        end do
        !$omp end do
        !$omp end parallel
        do i=1,NNN
            if( .not. outside(i) )then
                call pspec%add(tmp_imgs(i), SPECW)
            endif
            call tmp_imgs(i)%zero_and_unflag_ft
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

    ! PRIVATE FUNCTIONS

    function center_reference( )result( shift )
        use simple_binimage, only: binimage
        use simple_segmentation, only: otsu_robust_fast
        type(binimage)          :: img, tmp, tmpcc
        real,    pointer     :: rmat(:,:,:), rmat_cc(:,:,:)
        integer, allocatable :: sz(:)
        real                 :: shift(3), thresh(3)
        integer              :: loc, ldim
        ldim = params_glob%box
        img  = reference
        ! low-pass
        call img%bp(0., params_glob%cenlp)
        call img%ifft()
        ! mask
        call img%mask(0.45*real(params_glob%box), 'hard')
        ! median filtering
        call img%real_space_filter(5,'median')
        ! thresholding
        call img%get_rmat_ptr(rmat)
        where(rmat(1:ldim,1:ldim,1) < TINY) rmat(1:ldim,1:ldim,1) = 0.
        tmp = img
        ! binarize
        call otsu_robust_fast(tmp,is2d=.true.,noneg=.true.,thresh=thresh)
        ! median filtering again
        call tmp%real_space_filter(3,'median')
        ! identify biggest connected component
        call tmp%find_ccs(tmpcc)
        sz  = tmpcc%size_ccs()
        loc = maxloc(sz,dim=1)
        ! set to zero all the other connected components
        call tmpcc%get_rmat_ptr(rmat_cc)
        where(abs(rmat_cc(1:ldim,1:ldim,1)-loc) < TINY)
            rmat_cc(1:ldim,1:ldim,1) = 1.
        elsewhere
            rmat_cc(1:ldim,1:ldim,1) = 0.
        endwhere
        ! keep only thresholded values within largest CC
        call img%mul(tmpcc)
        call img%masscen(shift)
        ! cleanup
        call tmp%kill
        call tmpcc%kill
        call img%kill
    end function center_reference

    ! GETTERS

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
        framenames => null()
        do i=1,NNN
            call neigh_imgs(i)%kill
            call tmp_imgs(i)%kill
            call backgr_imgs(i)%kill
        enddo
        do i=1,params_glob%nframesgrp
            call ptcls(i)%kill
            call ptcls_saved(i)%kill
            if(trim(params_glob%filter).eq.'tv') call tv(i)%kill
        end do
        if(trim(params_glob%filter).eq.'tv') deallocate(tv)
        deallocate(particle_locations,ptcls,ptcls_saved)
        call frame_img%kill
        call frame_avg%kill
        call reference%kill
        call ptcl_target%kill
        call pspec%kill
    end subroutine kill_tracker

end module simple_tseries_tracker

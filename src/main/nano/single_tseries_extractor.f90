!@descr: Extractions routine for nanoparticles time series intended for frames and given per frame coordinates
module single_tseries_extractor
use simple_core_module_api
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_sp_project, only: sp_project
implicit none

public :: init_trajectory_extractor, extract_trajectory, kill_trajectory_extractor
private
#include "simple_local_flags.inc"

integer, parameter :: NNN = 8

class(parameters), pointer     :: p_ptr => null()
class(sp_project), pointer     :: p_spproj => null()
type(string),      allocatable :: frame_names(:)
integer,           allocatable :: particle_locations(:,:)
type(string)                   :: dir, fbody
type(image)                    :: frame_img, ptcl_target, pspec, pspec_nn
type(image)                    :: backgr_imgs(NNN), tmp_imgs(NNN)
type(string)                   :: stkname
integer                        :: ldim(3), nframes, box
logical                        :: l_neg

contains

    subroutine init_trajectory_extractor( params, spproj, dir_in, fbody_in )
        class(parameters), target, intent(in) :: params
        class(sp_project), target, intent(in) :: spproj
        class(string),             intent(in) :: dir_in, fbody_in
        integer      :: n, i
        call kill_trajectory_extractor
        ! set pointers
        p_ptr    => params
        p_spproj => spproj 
        ! set constants
        dir   = dir_in
        fbody = fbody_in
        l_neg = trim(p_ptr%neg) .eq. 'yes'
        ! names & dimensions
        if( spproj%os_mic%isthere('intg') )then
            nframes = spproj%get_nintgs()
            allocate(frame_names(nframes))
            do i = 1, nframes
                frame_names(i) = spproj%os_mic%get_str(i,'intg')
            enddo
        else
            nframes = spproj%get_nframes()
            allocate(frame_names(nframes))
            do i = 1, nframes
                frame_names(i) = spproj%os_mic%get_str(i,'frame')
            enddo
        endif
        call find_ldim_nptcls(frame_names(1), ldim, n)
        if( n /= 1 .or. ldim(3) /= 1 )then
            write(logfhandle,*) 'ldim(3): ', ldim(3)
            write(logfhandle,*) 'nframes: ', n
            THROW_HARD('init_tracker; assumes one frame per file')
        endif
        ! init images
        box = p_ptr%box_extract
        call frame_img%new(ldim, p_ptr%smpd)
        call ptcl_target%new([box,box,1], p_ptr%smpd)
        call pspec%new([box,box,1], p_ptr%smpd)
        call pspec_nn%new([box,box,1], p_ptr%smpd)
        do i=1,NNN
            call backgr_imgs(i)%new([box,box,1], p_ptr%smpd,wthreads=.false.)
            call tmp_imgs(i)%new([box,box,1], p_ptr%smpd,wthreads=.false.)
        end do
    end subroutine init_trajectory_extractor

    subroutine extract_trajectory()
        type(oris)            :: track_os
        integer, allocatable  :: track_totrack(:), track2frame(:)
        integer :: ldim_sc(2), iframe, xind,yind, i,  cnt, nrange, noutside, nl
        integer :: first_frame, last_frame, ntrack, first, last
        real    :: x, y, xp, yp, xscale, yscale, delta, T
        write(logfhandle,'(A)') ">>> READING TRAJECTORY"
        nl = nlines(p_ptr%infile)
        call track_os%new(nl, is_ptcl=.false.)
        call track_os%read(p_ptr%infile)
        ntrack = track_os%get_noris()
        ! dimensions
        ldim_sc(1) = track_os%get_int(1,'xdim')
        ldim_sc(2) = track_os%get_int(1,'ydim')
        xscale     = real(ldim_sc(1))/real(ldim(1))
        yscale     = real(ldim_sc(2))/real(ldim(2))
        call build_mapping(track_os, track2frame, track_totrack)
        first_frame = minval(track2frame)
        last_frame  = maxval(track_totrack)
        write(logfhandle,'(A)') ">>> GENERATING PARTICLES COORDINATES"
        allocate(particle_locations(2,first_frame:last_frame))
        T = 0.5*real(box)
        do i = 1,ntrack
            first = track2frame(i)
            if( first < first_frame .or. first > last_frame )then
                THROW_HARD("frame index in coordinates file is out of range")
            endif
            last = track_totrack(i)
            ! 0-based center coordinates -> EMAN 0-based upper left corner coordinates
            x = track_os%get(i,'x') / xscale - T
            y = track_os%get(i,'y') / yscale - T
            if( i == ntrack ) then
                particle_locations(1, first:last) = nint(x)
                particle_locations(2, first:last) = nint(y)
            else
                nrange = last - first + 1
                particle_locations(1, first) = nint(x)
                particle_locations(2, first) = nint(y)
                if( nrange > 1 )then
                    xp = track_os%get(i+1,'x') / xscale - T
                    yp = track_os%get(i+1,'y') / yscale - T
                    do iframe = first, last
                        delta = real(iframe-first)/real(nrange-1)
                        particle_locations(1, iframe) = nint(x + (xp-x) * delta)
                        particle_locations(2, iframe) = nint(y + (yp-y) * delta)
                    enddo
                    particle_locations(1, last) = nint(xp)
                    particle_locations(2, last) = nint(yp)
                endif
            endif
            call progress_gfortran( i, ntrack)
        enddo
        ! outputs to be provided
        ! write(funit,'(I7,I7,I7,I7,I7)') xind, yind, p_ptr%box, p_ptr%box, -3
        ! write(funit2,'(2F12.6)') particle_locations(1:2,iframe) * p_ptr%smpd ! in angstroms
        write(logfhandle,'(A)') ">>> GENERATING PARTICLES AND BACKGROUND POWER SPECTRA"
        call pspec%zero_and_unflag_ft
        call ptcl_target%zero_and_unflag_ft
        stkname = dir%to_char()//'/'//fbody%to_char()//MRC_EXT
        nrange  = last_frame - first_frame + 1
        cnt     = 0
        do iframe = first_frame, last_frame
            cnt  = cnt + 1
            xind = particle_locations(1,iframe)
            yind = particle_locations(2,iframe)
            if( xind<0 .and. xind>-10 ) xind = 0
            if( yind<0 .and. yind>-10 ) yind = 0
            if( xind>ldim(1)-1 .and. xind<ldim(1)+10 ) xind = ldim(1)-1
            if( yind>ldim(2)-1 .and. yind<ldim(2)+10 ) yind = ldim(2)-1
            ! read frame, extract, normalize, flip signs, write particle
            call frame_img%read(frame_names(iframe),1)
            call frame_img%window([xind,yind,1], box, ptcl_target, noutside)
            call ptcl_target%norm
            if( l_neg ) call ptcl_target%neg()
            call ptcl_target%write(stkname, cnt)
            ! neighbors & spectrum
            call update_background_pspec([xind,yind])
            call pspec%add(pspec_nn, w=1./real(nrange))
            call pspec_nn%write(string(dir%to_char()//'/'//fbody%to_char()//'_background_pspec.mrc'),cnt)
            call progress_gfortran( iframe, nrange)
        end do
        ! average and write power spectrum for CTF estimation
        call pspec%dampen_pspec_central_cross
        call pspec%write(string(dir%to_char()//'/'//'pspec4ctf_estimation.mrc'))
        ! cleanup
        call track_os%kill
        deallocate(track2frame, track_totrack)
    end subroutine extract_trajectory

    subroutine build_mapping( track_os, track2frame, track_totrack )
        type(oris),           intent(in)    :: track_os
        integer, allocatable, intent(inout) :: track2frame(:), track_totrack(:)
        type(string)              :: absfname, fname, trackfname
        character(len=LONGSTRLEN) :: str, extstr, numstr
        integer                   :: i, j, iframe, fromt, tot, n
        n = track_os%get_noris()
        allocate(track2frame(n), track_totrack(n))
        do i = 1, n
            j = track_os%get_int(i,'index')
            absfname  = track_os%get_str(i,'filename')
            fname     = basename(absfname)
            extstr    = fname%to_char()
            call split_str(extstr, MRC_EXT, numstr)
            call split_str(numstr, '_', str)
            iframe = str2int(numstr)
            if( p_spproj%os_mic%isthere(iframe,'track_fname') )then
                fromt = p_spproj%os_mic%get_int(iframe,'fromt')
                if( iframe /= fromt ) THROW_HARD('frame index in coordinates file is out of range')
                trackfname = p_spproj%os_mic%get_str(iframe,'track_fname')
                if( trackfname /= absfname ) THROW_HARD('frame names are different')
            else
                THROW_HARD('Fatal frame indexing error')
            endif
            tot = p_spproj%os_mic%get_int(iframe,'tot')
            if( tot < fromt .or. tot > nframes ) THROW_HARD('TOT frame index in coordinates file is out of range')
            track2frame(i)   = fromt
            track_totrack(i) = tot
        enddo
    end subroutine build_mapping

    subroutine update_background_pspec( pos )
        integer, intent(in) :: pos(2)
        integer :: neigh(NNN,2), i, neigh_cnt
        logical :: outside(NNN)
        call identify_neighbours
        neigh_cnt = 0
        call pspec_nn%zero_and_unflag_ft
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static) reduction(+:neigh_cnt)
        do i=1,NNN
            call backgr_imgs(i)%zero_and_unflag_ft
            call tmp_imgs(i)%zero_and_unflag_ft
            call frame_img%window_slim(neigh(i,:), box, backgr_imgs(i), outside(i))
            if( .not.outside(i) )then
                neigh_cnt = neigh_cnt + 1
                if( l_neg ) call backgr_imgs(i)%neg
                call backgr_imgs(i)%norm
                call backgr_imgs(i)%zero_edgeavg
                call backgr_imgs(i)%fft
                call backgr_imgs(i)%ft2img('log', tmp_imgs(i))
            endif
        enddo
        !$omp end parallel do
        do i=1,NNN
            if( .not. outside(i) ) call pspec_nn%add(tmp_imgs(i), 1./real(neigh_cnt))
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

    end subroutine update_background_pspec

    subroutine kill_trajectory_extractor
        integer :: i
        do i=1,NNN
            call tmp_imgs(i)%kill
            call backgr_imgs(i)%kill
        enddo
        if( allocated(particle_locations) ) deallocate(particle_locations)
        call frame_img%kill
        call ptcl_target%kill
        call pspec%kill
        call pspec_nn%kill
        nullify(p_ptr, p_spproj)
    end subroutine kill_trajectory_extractor

end module single_tseries_extractor

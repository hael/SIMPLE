!@descr: optics-group assignment of micrographs: threshold clustering of beam-image shifts within tilt groups
!==============================================================================
! MODULE: simple_optics_groups
!
! PURPOSE:
!   Assigns every micrograph of a project to an optics group, and rebuilds the
!   optics segment, from the micrographs' beam-image shifts (shiftx, shifty):
!     1. tilt groups: the distinct 'tiltgrp' values when beam tilt is used,
!        otherwise one group;
!     2. within each tilt group, h_clust (simple_starproject_utils) merges
!        shifts closer than tilt_thres;
!     3. each resulting cluster is an optics group: os_mic gets 'ogid', and
!        os_optics gets one row per group with its population and centroid.
!
!   Until now this ran inside starproject_stream%stream_export_optics, as a
!   side effect of writing optics.star, followed by a hidden write of the
!   whole project. Here it only changes os_mic and os_optics in memory.
!
! HOME:
!   In src/main/stream/shared for now; it belongs beside the project's optics
!   routines (simple_sp_project_optics). The batch assign_optics_groups
!   (starproject%assign_optics) clusters the same quantity, read from the EPU
!   XML files, with the same h_clust, and could call this once those values
!   are on os_mic.
!
! TESTS:
!   simple_optics_groups_tester
!==============================================================================
module simple_optics_groups
use simple_defs,              only: logfhandle
use simple_string_utils,      only: int2str
use simple_math,              only: elim_dup
use simple_ori,               only: ori
use simple_sp_project,        only: sp_project
use simple_starproject_utils, only: h_clust
implicit none

public :: assign_optics_groups
private

contains

    !> Assigns optics groups to every micrograph of @p spproj, whatever its state, and
    !! rebuilds os_optics. Group ids are @p group_offset + 1, + 2, ...
    subroutine assign_optics_groups( spproj, tilt_thres, l_beamtilt, group_offset )
        class(sp_project), intent(inout) :: spproj
        real,              intent(in)    :: tilt_thres
        logical,           intent(in)    :: l_beamtilt
        integer,           intent(in)    :: group_offset
        type(ori)            :: template
        real,    allocatable :: tiltgrps(:), tilts_uniq(:), shiftxs(:), shiftys(:), states(:)
        real,    allocatable :: shifts(:,:), centroids(:,:), group_cx(:), group_cy(:)
        integer, allocatable :: tiltind(:), allinds(:), members(:), labels(:), populations(:), group_pops(:)
        integer :: nmics, ntilt, itilt, nmembers, ngroups, igroup, imic, iref
        nmics = spproj%os_mic%get_noris()
        if( nmics == 0 ) return
        allinds = [(imic, imic=1,nmics)]
        ! 1. tilt groups
        allocate(tiltind(nmics), source=1)
        ntilt = 1
        if( l_beamtilt )then
            tiltgrps = spproj%os_mic%get_all('tiltgrp')
            call elim_dup(tiltgrps, tilts_uniq)
            ntilt = size(tilts_uniq)
            do imic = 1,nmics
                tiltind(imic) = findloc(tilts_uniq, tiltgrps(imic), 1)
            enddo
        endif
        write(logfhandle,'(A,I8)') '>>> # TILT GROUPS ASSIGNED : ', ntilt
        ! 2. shift clusters within each tilt group
        write(logfhandle,'(A,F8.2)') '>>> CLUSTERING TILT GROUPS USING SHIFTS AND THRESHOLD : ', tilt_thres
        shiftxs = spproj%os_mic%get_all('shiftx')
        shiftys = spproj%os_mic%get_all('shifty')
        allocate(group_cx(0), group_cy(0), group_pops(0))
        ngroups = 0
        do itilt = 1,ntilt
            write(logfhandle,'(A,I8)') '      CLUSTERING TILT GROUP : ', itilt
            members  = pack(allinds, tiltind == itilt)
            nmembers = size(members)
            if( nmembers == 0 ) cycle
            allocate(shifts(nmembers,2), labels(nmembers))
            shifts(:,1) = shiftxs(members)
            shifts(:,2) = shiftys(members)
            call h_clust(shifts, tilt_thres, labels, centroids, populations)
            do imic = 1,nmembers
                call spproj%os_mic%set(members(imic), 'ogid', real(group_offset + ngroups + labels(imic)))
            enddo
            group_cx   = [group_cx,   centroids(:,1)]
            group_cy   = [group_cy,   centroids(:,2)]
            group_pops = [group_pops, populations]
            ngroups    = ngroups + size(populations)
            write(logfhandle,'(A,I8)') '        # SHIFT GROUPS ASSIGNED : ', size(populations)
            deallocate(shifts, labels)
        enddo
        ! 3. the optics segment; CTF constants from the first accepted micrograph, else the first
        states = spproj%os_mic%get_all('state')
        iref   = findloc(states > 0., .true., 1)
        if( iref == 0 ) iref = 1
        call template%new(.false.)
        call template%set('smpd',   spproj%os_mic%get(iref, 'smpd'))
        call template%set('cs',     spproj%os_mic%get(iref, 'cs'))
        call template%set('kv',     spproj%os_mic%get(iref, 'kv'))
        call template%set('fraca',  spproj%os_mic%get(iref, 'fraca'))
        call template%set('state',  1.0)
        call template%set('pop',    0.0)
        call template%set('ogid',   0.0)
        call template%set('opcx',   0.0)
        call template%set('opcy',   0.0)
        call template%set('ogname', 'opticsgroup')
        call spproj%os_optics%new(ngroups, is_ptcl=.false.)
        do igroup = 1,ngroups
            call spproj%os_optics%append(igroup, template)
            call spproj%os_optics%set(igroup, 'ogid',   real(group_offset + igroup))
            call spproj%os_optics%set(igroup, 'pop',    real(group_pops(igroup)))
            call spproj%os_optics%set(igroup, 'opcx',   group_cx(igroup))
            call spproj%os_optics%set(igroup, 'opcy',   group_cy(igroup))
            call spproj%os_optics%set(igroup, 'ogname', 'opticsgroup'//int2str(group_offset + igroup))
        enddo
        call template%kill()
    end subroutine assign_optics_groups

end module simple_optics_groups

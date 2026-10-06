!@descr: optics-group assignment of micrographs: single-linkage grouping of beam-image shifts within tilt groups
!==============================================================================
! MODULE: simple_optics_groups
!
! PURPOSE:
!   Assigns every micrograph of a project to an optics group, and rebuilds the
!   optics segment, from the micrographs' beam-image shifts (shiftx, shifty):
!     1. tilt groups: the distinct 'tiltgrp' values when beam tilt is used,
!        otherwise one group;
!     2. within each tilt group, single linkage: two micrographs are linked
!        when their shifts lie within tilt_thres, and each connected set is
!        an optics group (link_shifts: a grid of cells tilt_thres/sqrt(2)
!        wide and a union-find, no distance matrix, memory linear in the
!        micrographs; equal shifts, as without dir_meta, cost one pass);
!     3. os_mic gets 'ogid', and os_optics one row per group, in id order,
!        with its population and centroid (the mean shift of its members).
!   Every call groups all micrographs again (stream optics assignment does so
!   on every pass), so the groups do not depend on the order micrographs
!   arrived in. Called with last_ogid, the ids are kept across calls: a group
!   takes the id most of its micrographs had (each id once, the group sharing
!   most first) and a new group the next id never given, so ids can have
!   gaps once groups merge (kept_ids).
!
! HOME:
!   In src/main/stream/shared for now; it belongs beside the project's optics
!   routines (simple_sp_project_optics). The batch assign_optics_groups
!   (starproject%assign_optics) clusters the same quantity, read from the EPU
!   XML files, with h_clust (simple_starproject_utils), which keeps a distance
!   matrix; it could call this once those values are on os_mic.
!
! TESTS:
!   simple_optics_groups_tester
!==============================================================================
module simple_optics_groups
use simple_defs,          only: logfhandle
use simple_string_utils,  only: int2str
use simple_math,          only: elim_dup
use simple_srch_sort_loc, only: hpsort
use simple_ori,           only: ori
use simple_sp_project,    only: sp_project
implicit none

public :: assign_optics_groups
private

contains

    !> Assigns optics groups to every micrograph of @p spproj, whatever its state, and rebuilds
    !! os_optics. Without @p last_ogid, group ids are @p group_offset + 1, + 2, ... in the order of
    !! each group's first micrograph. With @p last_ogid, the highest id given so far (updated), the
    !! ids the micrographs had ('ogid' > 0) are kept as kept_ids says, and @p group_offset is not
    !! used.
    subroutine assign_optics_groups( spproj, tilt_thres, l_beamtilt, group_offset, last_ogid )
        class(sp_project), intent(inout) :: spproj
        real,              intent(in)    :: tilt_thres
        logical,           intent(in)    :: l_beamtilt
        integer,           intent(in)    :: group_offset
        integer, optional, intent(inout) :: last_ogid
        type(ori)            :: template
        real,    allocatable :: tiltgrps(:), tilts_uniq(:), shiftxs(:), shiftys(:), states(:)
        real,    allocatable :: group_cx(:), group_cy(:)
        integer, allocatable :: tiltind(:), labels(:), prev_ogids(:), group_ids(:), group_pops(:), rows(:), allinds(:)
        integer :: nmics, ntilt, itilt, ngroups, ngroups_tilt, igroup, imic, iref, irow
        nmics = spproj%os_mic%get_noris()
        if( nmics == 0 ) return
        allinds = [(imic, imic=1,nmics)]
        ! the ids the micrographs had, before they are given new ones
        allocate(prev_ogids(nmics), source=0)
        if( present(last_ogid) )then
            do imic = 1,nmics
                if( spproj%os_mic%isthere(imic, 'ogid') ) prev_ogids(imic) = max(0, spproj%os_mic%get_int(imic, 'ogid'))
            enddo
        endif
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
        ! 2. single-linkage groups within each tilt group
        write(logfhandle,'(A,F8.2)') '>>> GROUPING TILT GROUPS BY SHIFTS WITHIN : ', tilt_thres
        shiftxs = spproj%os_mic%get_all('shiftx')
        shiftys = spproj%os_mic%get_all('shifty')
        allocate(labels(nmics), source=0)
        ngroups = 0
        do itilt = 1,ntilt
            call link_shifts(pack(allinds, tiltind == itilt), shiftxs, shiftys, tilt_thres, labels, ngroups, ngroups_tilt)
            write(logfhandle,'(A,I8,A,I8)') '      TILT GROUP ', itilt, ' # SHIFT GROUPS ASSIGNED : ', ngroups_tilt
        enddo
        ! populations and centroids
        allocate(group_pops(ngroups), source=0)
        allocate(group_cx(ngroups), group_cy(ngroups), source=0.)
        do imic = 1,nmics
            igroup             = labels(imic)
            group_pops(igroup) = group_pops(igroup) + 1
            group_cx(igroup)   = group_cx(igroup) + shiftxs(imic)
            group_cy(igroup)   = group_cy(igroup) + shiftys(imic)
        enddo
        group_cx = group_cx / real(max(1, group_pops))
        group_cy = group_cy / real(max(1, group_pops))
        ! 3. the ids
        if( present(last_ogid) )then
            group_ids = kept_ids(labels, prev_ogids, ngroups, last_ogid)
        else
            group_ids = [(group_offset + igroup, igroup=1,ngroups)]
        endif
        do imic = 1,nmics
            call spproj%os_mic%set(imic, 'ogid', real(group_ids(labels(imic))))
        enddo
        ! 4. the optics segment, a row per group in id order; CTF constants from the first accepted
        ! micrograph, else the first
        rows = [(igroup, igroup=1,ngroups)]
        call hpsort(rows, id_lt)
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
        do irow = 1,ngroups
            igroup = rows(irow)
            call spproj%os_optics%append(irow, template)
            call spproj%os_optics%set(irow, 'ogid',   real(group_ids(igroup)))
            call spproj%os_optics%set(irow, 'pop',    real(group_pops(igroup)))
            call spproj%os_optics%set(irow, 'opcx',   group_cx(igroup))
            call spproj%os_optics%set(irow, 'opcy',   group_cy(igroup))
            call spproj%os_optics%set(irow, 'ogname', 'opticsgroup'//int2str(group_ids(igroup)))
        enddo
        call template%kill()

    contains

        logical function id_lt( g1, g2 )
            integer, intent(in) :: g1, g2
            id_lt = group_ids(g1) < group_ids(g2)
        end function id_lt

    end subroutine assign_optics_groups

    !> Single linkage of the shifts @p xs, @p ys of @p members (micrograph indices): two are linked
    !! when their shifts lie within @p thres, and each connected set gets the next label after
    !! @p ngroups (updated), in the order of the sets' first members; @p nnew is the number of sets.
    !! The shifts are binned into cells thres/sqrt(2) wide, so the members of a cell are linked
    !! without comparison and a shift is compared only with those of the cells up to two away; two
    !! cells already in one set are not compared, and their first link ends their comparison. A
    !! threshold that is not positive links equal shifts only.
    subroutine link_shifts( members, xs, ys, thres, labels, ngroups, nnew )
        integer, intent(in)    :: members(:)
        real,    intent(in)    :: xs(:), ys(:), thres
        integer, intent(inout) :: labels(:), ngroups
        integer, intent(out)   :: nnew
        integer(kind=8), allocatable :: cx(:), cy(:)
        integer,         allocatable :: parent(:), order(:), cell_first(:), cell_last(:), set_label(:)
        real    :: h
        integer :: n, i, k, nb, ncells, dx, dy, iroot
        logical :: l_cells_linked
        n    = size(members)
        nnew = 0
        if( n == 0 ) return
        ! with a positive threshold, any two shifts of a cell lie within it
        l_cells_linked = thres > 0.
        h = 1.
        if( l_cells_linked ) h = thres / sqrt(2.)
        allocate(cx(n), cy(n), parent(n), order(n))
        do i = 1,n
            cx(i)     = floor(xs(members(i)) / h, kind=8)
            cy(i)     = floor(ys(members(i)) / h, kind=8)
            parent(i) = i
            order(i)  = i
        enddo
        call hpsort(order, cell_lt)
        ! the cells: runs of the same cell in that order
        allocate(cell_first(n), cell_last(n))
        ncells = 0
        do i = 1,n
            if( i > 1 )then
                if( cx(order(i)) == cx(order(i-1)) .and. cy(order(i)) == cy(order(i-1)) )then
                    cell_last(ncells) = i
                    cycle
                endif
            endif
            ncells             = ncells + 1
            cell_first(ncells) = i
            cell_last(ncells)  = i
        enddo
        ! within a cell, then each pair of cells once: the neighbours after a cell in that order
        do k = 1,ncells
            if( l_cells_linked )then
                do i = cell_first(k)+1,cell_last(k)
                    call unite(order(cell_first(k)), order(i))
                enddo
            else
                call link_cells(k, k)
            endif
        enddo
        do k = 1,ncells
            do dx = 0,2
                do dy = -2,2
                    if( dx == 0 .and. dy <= 0 ) cycle
                    nb = find_cell(cx(order(cell_first(k))) + dx, cy(order(cell_first(k))) + dy)
                    if( nb > 0 ) call link_cells(k, nb)
                enddo
            enddo
        enddo
        ! the labels, in the order of the sets' first members
        allocate(set_label(n), source=0)
        do i = 1,n
            iroot = root(i)
            if( set_label(iroot) == 0 )then
                nnew             = nnew + 1
                set_label(iroot) = ngroups + nnew
            endif
            labels(members(i)) = set_label(iroot)
        enddo
        ngroups = ngroups + nnew

    contains

        logical function cell_lt( p1, p2 )
            integer, intent(in) :: p1, p2
            if( cx(p1) /= cx(p2) )then
                cell_lt = cx(p1) < cx(p2)
            else
                cell_lt = cy(p1) < cy(p2)
            endif
        end function cell_lt

        ! the cell (qx, qy), 0 when no shift falls in it (binary search over the cells in order)
        integer function find_cell( qx, qy )
            integer(kind=8), intent(in) :: qx, qy
            integer         :: lo, hi, mid
            integer(kind=8) :: mx, my
            find_cell = 0
            lo = 1
            hi = ncells
            do while( lo <= hi )
                mid = (lo + hi) / 2
                mx  = cx(order(cell_first(mid)))
                my  = cy(order(cell_first(mid)))
                if( mx == qx .and. my == qy )then
                    find_cell = mid
                    return
                endif
                if( mx < qx .or. (mx == qx .and. my < qy) )then
                    lo = mid + 1
                else
                    hi = mid - 1
                endif
            enddo
        end function find_cell

        ! the members of cells @p k1 and @p k2 within the threshold are linked
        subroutine link_cells( k1, k2 )
            integer, intent(in) :: k1, k2
            integer :: i1, i2, a, b
            if( l_cells_linked .and. root(order(cell_first(k1))) == root(order(cell_first(k2))) ) return
            do i1 = cell_first(k1),cell_last(k1)
                a = order(i1)
                do i2 = cell_first(k2),cell_last(k2)
                    b = order(i2)
                    if( k1 == k2 .and. i2 <= i1 ) cycle
                    if( (xs(members(a)) - xs(members(b)))**2 + (ys(members(a)) - ys(members(b)))**2 > thres**2 ) cycle
                    call unite(a, b)
                    ! two cells whose members are linked among themselves: one link joins them
                    if( l_cells_linked ) return
                enddo
            enddo
        end subroutine link_cells

        integer function root( i )
            integer, intent(in) :: i
            root = i
            do while( parent(root) /= root )
                parent(root) = parent(parent(root))
                root         = parent(root)
            enddo
        end function root

        subroutine unite( a, b )
            integer, intent(in) :: a, b
            integer :: ra, rb
            ra = root(a)
            rb = root(b)
            if( ra == rb ) return
            if( ra < rb )then
                parent(rb) = ra
            else
                parent(ra) = rb
            endif
        end subroutine unite

    end subroutine link_shifts

    !> The ids of the @p ngroups groups of @p labels, kept from the ids the micrographs had
    !! (@p prev_ogids, 0: none): the pairs (group, previous id) are taken in order of shared
    !! micrographs, most first (then the lower group, then the lower id), each group and each
    !! previous id once. A group left without one takes the next id after @p last_ogid (updated,
    !! and never below a previous id), so an id is never given to two groups.
    function kept_ids( labels, prev_ogids, ngroups, last_ogid ) result( ids )
        integer, intent(in)    :: labels(:), prev_ogids(:), ngroups
        integer, intent(inout) :: last_ogid
        integer, allocatable :: ids(:), mics(:), pair_group(:), pair_id(:), pair_count(:), pairs(:)
        logical, allocatable :: id_taken(:)
        integer :: nmics, npairs, i, ip, igroup
        allocate(ids(ngroups), source=0)
        nmics = size(labels)
        if( nmics > 0 ) last_ogid = max(last_ogid, maxval(prev_ogids))
        ! the pairs present and their counts: the micrographs with a previous id, in (group, id) order
        mics = pack([(i, i=1,nmics)], prev_ogids > 0)
        call hpsort(mics, pair_lt)
        allocate(pair_group(size(mics)), pair_id(size(mics)), pair_count(size(mics)))
        npairs = 0
        do i = 1,size(mics)
            if( npairs > 0 )then
                if( labels(mics(i)) == pair_group(npairs) .and. prev_ogids(mics(i)) == pair_id(npairs) )then
                    pair_count(npairs) = pair_count(npairs) + 1
                    cycle
                endif
            endif
            npairs             = npairs + 1
            pair_group(npairs) = labels(mics(i))
            pair_id(npairs)    = prev_ogids(mics(i))
            pair_count(npairs) = 1
        enddo
        pairs = [(ip, ip=1,npairs)]
        call hpsort(pairs, pair_first)
        allocate(id_taken(max(1, last_ogid)), source=.false.)
        do i = 1,npairs
            ip = pairs(i)
            if( ids(pair_group(ip)) /= 0 ) cycle
            if( id_taken(pair_id(ip)) ) cycle
            ids(pair_group(ip))     = pair_id(ip)
            id_taken(pair_id(ip)) = .true.
        enddo
        do igroup = 1,ngroups
            if( ids(igroup) /= 0 ) cycle
            last_ogid   = last_ogid + 1
            ids(igroup) = last_ogid
        enddo

    contains

        logical function pair_lt( m1, m2 )
            integer, intent(in) :: m1, m2
            if( labels(m1) /= labels(m2) )then
                pair_lt = labels(m1) < labels(m2)
            else
                pair_lt = prev_ogids(m1) < prev_ogids(m2)
            endif
        end function pair_lt

        ! the order pairs are taken in: most shared micrographs, then the lower group, then the lower id
        logical function pair_first( p1, p2 )
            integer, intent(in) :: p1, p2
            if( pair_count(p1) /= pair_count(p2) )then
                pair_first = pair_count(p1) > pair_count(p2)
            else if( pair_group(p1) /= pair_group(p2) )then
                pair_first = pair_group(p1) < pair_group(p2)
            else
                pair_first = pair_id(p1) < pair_id(p2)
            endif
        end function pair_first

    end function kept_ids

end module simple_optics_groups

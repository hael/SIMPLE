!@descr: routines for distributed SIMPLE execution
module simple_map_reduce
use simple_defs
use simple_fileio
use simple_jiffys
use simple_srch_sort_loc
use simple_string
use simple_string_utils
implicit none
#include "simple_local_flags.inc"

contains

    !>  \brief  for generating balanced partitions of nobjs objects
    function split_nobjs_even( nobjs, nparts, szmax ) result( parts )
        integer,           intent(in)  :: nobjs, nparts
        integer, optional, intent(out) :: szmax
        integer, allocatable :: parts(:,:)
        integer :: nobjs_per_part, leftover, istop, istart, ipart, sszmax
        allocate(parts(nparts,2))
        nobjs_per_part = nobjs/nparts
        leftover = nobjs-nobjs_per_part*nparts
        istop  = 0
        istart = 0
        sszmax = 0
        do ipart=1,nparts
            if( ipart == nparts )then
                istart = istop+1
                istop  = nobjs
            else
                if( leftover == 0 )then
                    istart = istop+1
                    istop  = istart+nobjs_per_part-1
                else
                    istop    = ipart*(nobjs_per_part+1)
                    istart   = istop-(nobjs_per_part+1)+1
                    leftover = leftover-1
                endif
            endif
            parts(ipart,1) = istart
            parts(ipart,2) = istop
            sszmax         = max(sszmax,istop - istart + 1)
        end do
        if( present(szmax) ) szmax = sszmax
    end function split_nobjs_even

    !>  \brief  contiguous partitions of size(l_active) objects that balance the active ones: part k
    !!          ends on the last active object of the k-th share of split_nobjs_even(count(l_active),
    !!          nparts), so inactive objects join the part of the next active object (trailing ones
    !!          the last part) and every part keeps at least one object. With every object active the
    !!          result is split_nobjs_even; with none active, or fewer objects than parts, it falls back
    !!          to it. szmax returns the largest number of active objects in a part.
    function split_nobjs_active( l_active, nparts, szmax ) result( parts )
        logical,           intent(in)  :: l_active(:)
        integer,           intent(in)  :: nparts
        integer, optional, intent(out) :: szmax
        integer, allocatable :: parts(:,:), shares(:,:), active_inds(:)
        integer              :: nobjs, nactive, ipart, i, last
        nobjs   = size(l_active)
        nactive = count(l_active)
        if( nactive == 0 .or. nobjs < nparts )then
            parts = split_nobjs_even(nobjs, nparts)
        else
            shares      = split_nobjs_even(nactive, nparts)
            active_inds = pack([(i, i=1,nobjs)], l_active)
            allocate(parts(nparts,2))
            last = 0
            do ipart = 1, nparts
                parts(ipart,1) = last + 1
                if( ipart == nparts )then
                    last = nobjs
                else
                    ! shares(ipart,2) is the number of active objects in parts 1..ipart
                    last = active_inds(shares(ipart,2))
                    last = min(max(last, parts(ipart,1)), nobjs - (nparts - ipart))
                endif
                parts(ipart,2) = last
            end do
        endif
        if( present(szmax) )then
            szmax = 0
            do ipart = 1, nparts
                szmax = max(szmax, count(l_active(parts(ipart,1):parts(ipart,2))))
            end do
        endif
    end function split_nobjs_active

    !>  \brief  for generating balanced partitions for pairwise calculations on nobjs ojects
    subroutine split_pairs_in_parts( nobjs, nparts )
        integer, intent(in)  :: nobjs  !< number objects to analyse in pairs
        integer, intent(in)  :: nparts !< number of partitions (nodes) for parallel execution
        integer              :: npairs, cnt, funit, i, j
        integer              :: ipart, io_stat, numlen
        integer, allocatable :: pairs(:,:), parts(:,:)
        type(string) :: fname
        ! generate all pairs
        npairs = (nobjs*(nobjs-1))/2
        allocate( pairs(npairs,2))
        cnt = 0
        do i=1,nobjs-1
            do j=i+1,nobjs
                cnt = cnt+1
                pairs(cnt,1) = i
                pairs(cnt,2) = j
            end do
        end do
        ! generate balanced partitions of pairs
        parts  = split_nobjs_even(npairs, nparts)     ! realloc lhs
        numlen = len(int2str(nparts))
        ! write the partitions
        do ipart=1,nparts
            call progress(ipart,nparts)
            fname = 'pairs_part' // int2str_pad(ipart,numlen) // '.bin'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, access='STREAM',iostat=io_stat)
            call fileiochk('mapreduce ;split_pairs_in_parts '//fname%to_char(), io_stat)
            write(unit=funit,pos=1,iostat=io_stat) pairs(parts(ipart,1):parts(ipart,2),:)
            ! Check if the write was successful
            if( io_stat .ne. 0 )&
                call fileiochk('mapreduce ;split_pairs_in_parts writing to '//fname%to_char(), io_stat)
            call fclose(funit)
            call fname%kill
        end do
        deallocate(pairs, parts)
    end subroutine split_pairs_in_parts

end module simple_map_reduce

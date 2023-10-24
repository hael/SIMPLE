module simple_srch_sort_loc
use simple_defs
use simple_error, only: simple_exception
implicit none

interface find
    module procedure find_1, find_2
end interface

interface hpsort
    module procedure hpsort_1, hpsort_2, hpsort_3, hpsort_4, hpsort_5
end interface

interface locate
   module procedure locate_1, locate_2
end interface

interface reverse
    module procedure reverse_iarr, reverse_rarr
end interface

contains

    !>   for finding closest element in an ordered list
    subroutine find_1( arr, n, x, j, dist1 )
        integer, intent(in)  :: n           !< size of list
        real,    intent(in)  :: arr(n), x   !< list and search value
        real,    intent(out) :: dist1
        integer, intent(out) :: j
        real                 :: dist2
        j = max(1,locate_1( arr, n, x ))
        dist1 = arr(j)-x
        if( j < n ) then
            dist2 = arr(j+1)-x
            if( abs(dist1) >= abs(dist2) )then
                j = j+1
                dist1 = dist2
            endif
        endif
    end subroutine find_1

    !>   for finding closest element in an ordered list
    subroutine find_2( arr, n, x, j, dist1 )
        integer, intent(in)  :: n         !< size of list
        integer, intent(in)  :: arr(n), x !< list and search value
        integer, intent(out) :: dist1
        integer, intent(out) :: j
        integer              :: dist2
        j = max(1,locate_2( arr, n, x ))
        dist1 = arr(j)-x
        if( j < n ) then
            dist2 = arr(j+1)-x
            if( abs(dist1) >= abs(dist2) )then
                j = j+1
                dist1 = dist2
            endif
        endif
    end subroutine find_2

    ! heapsort from numerical recipes (largest last)
    subroutine hpsort_1( rarr, iarr )
        real,    intent(inout) :: rarr(:)
        integer, intent(inout) :: iarr(:)
        integer :: i, ir, j, l, ia, n
        real    :: ra
        n = size(rarr)
        if( n /= size(iarr) )&
        &call simple_exception('nonconforming array sizes; hpsort_1', __FILENAME__ , __LINE__)
        if( n < 2 ) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
                ia = iarr(l)
            else
                ra = rarr(ir)
                ia = iarr(ir)
                rarr(ir) = rarr(1)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    iarr(1) = ia
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j) < rarr(j+1)) j = j+1
                endif
                if(ra < rarr(j))then
                    rarr(i) = rarr(j)
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
                iarr(i) = ia
            end do
        end do
    end subroutine hpsort_1

    ! heapsort from numerical recepies (largest last)
    subroutine hpsort_2( iarr )
        integer, intent(inout) :: iarr(:)
        integer :: i, ir, j, l, ra, n
        n = size(iarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = iarr(l)
            else
                ra = iarr(ir)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(iarr(j) < iarr(j+1)) j = j+1
                endif
                if(ra < iarr(j))then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ra
            end do
        end do
    end subroutine hpsort_2

    ! heapsort from numerical recepies (largest last)
    subroutine hpsort_3( iarr, p1_lt_p2 )
        integer, intent(inout) :: iarr(:)
        interface
            function p1_lt_p2( p1, p2 ) result( val )
                integer, intent(in) :: p1, p2
                logical :: val
            end function p1_lt_p2
        end interface
        integer :: i, ir, j, l, ra, n
        n = size(iarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = iarr(l)
            else
                ra = iarr(ir)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(p1_lt_p2(iarr(j),iarr(j+1))) j = j+1
                endif
                if(p1_lt_p2(ra,iarr(j)))then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ra
            end do
        end do
    end subroutine hpsort_3

    ! heapsort from numerical recepies (largest last)
    subroutine hpsort_4( rarr )
        real, intent(inout) :: rarr(:)
        integer :: i, ir, j, l, n
        real    :: ra
        n = size(rarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
            else
                ra = rarr(ir)
                rarr(ir) = rarr(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j) < rarr(j+1)) j = j+1
                endif
                if(ra < rarr(j))then
                    rarr(i) = rarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
            end do
        end do
    end subroutine hpsort_4

    ! heapsort from numerical recepies (largest last)
    subroutine hpsort_5( rarr, rarr2 )
        real, intent(inout) :: rarr(:), rarr2(:)
        integer :: i, ir, j, l, n
        real    :: ra, ra2
        n = size(rarr)
        if( n /= size(rarr2) )&
        &call simple_exception('nonconforming array sizes; hpsort_5', __FILENAME__ , __LINE__)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
                ra2 = rarr2(l)
            else
                ra = rarr(ir)
                ra2 = rarr2(ir)
                rarr(ir) = rarr(1)
                rarr2(ir) = rarr2(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    rarr2(1) = ra2
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j) < rarr(j+1)) j = j+1
                endif
                if(ra < rarr(j))then
                    rarr(i) = rarr(j)
                    rarr2(i) = rarr2(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
                rarr2(i) = ra2
            end do
        end do
    end subroutine hpsort_5
    
    !>   given an array arr(1:n), and given value x, locate returns a value j such that x is
    !!         between arr(j) and arr(j+1). arr(1:n) must be monotonic, either increasing or decreasing.
    !!         j=0 or j=n is returned to indicate that x is out of range, from numerical recepies
    pure function locate_1( arr, n, x ) result( j )
        integer, intent(in) :: n
        real, intent(in)    :: arr(n), x
        integer             :: jl, jm, ju, j
        jl = 0                ! initialize lower
        ju = n+1              ! and upper limits
        do while( ju-jl > 1 ) ! if we are not yet done
            jm = (ju+jl)/2    ! compute a midpoint
            if((arr(n) >= arr(1)).eqv.(x >= arr(jm))) then
                jl = jm       ! and replace either the lower limit
            else
                ju = jm       ! or the upper limit, as appropriate
            endif
        end do
        if( abs(x-arr(1)) < TINY )then
            j = 1
        else if( abs(x-arr(n)) < TINY )then
            j = n-1
        else
            j = jl
        endif
    end function locate_1

    !>   given an array arr(1:n), and given value x, locate returns a value j such that x is
    ! between arr(j) and arr(j+1). arr(1:n) must be monotonic, either increasing or decreasing.
    ! j=0 or j=n is returned to indicate that x is out of range, from numerical recepies
    pure function locate_2( arr, n, x ) result( j )
        integer, intent(in) :: n             !< size of list
        integer, intent(in) :: arr(n), x     !< list and search value
        integer             :: jl, jm, ju, j
        jl = 0                ! initialize lower
        ju = n+1              ! and upper limits
        do while( ju-jl > 1 ) ! if we are not yet done
            jm = (ju+jl)/2    ! compute a midpoint
            if((arr(n) >= arr(1)).eqv.(x >= arr(jm))) then
                jl = jm       ! and replace either the lower limit
            else
                ju = jm       ! or the upper limit, as appropriate
            endif
        end do
        if( x == arr(1) )then
            j = 1
        else if( x == arr(n) )then
            j = n-1
        else
            j = jl
        endif
    end function locate_2

    function maxnloc( rarr, n ) result( loc )
        real,    intent(in) :: rarr(:)
        integer, intent(in) :: n
        real    :: arr(n), val
        integer :: loc(n), i, j, sz
        logical :: val_gt(n)
        sz = size(rarr)
        if( sz < n )&
        &call simple_exception('cannot identify more maxima than elements in the array; maxnloc', __FILENAME__ , __LINE__)
        loc = (/(i,i=1,n)/)
        arr = rarr(:n)
        call hpsort(arr, loc)
        call reverse(loc)
        if( sz == n ) return
        do i=n+1,sz
            val    = rarr(i)
            val_gt = val > rarr(loc)
            if( any(val_gt) )then
                do j=1,n
                    if( val_gt(j) )then
                        if( j == 1 )then
                            loc = [i,loc(1:n-1)]
                            exit
                        else if( j == n )then
                            loc(n) = i
                            exit
                        else
                            ! 1:j-1 untouched
                            ! insert i @ location j
                            ! j+1:n-1 forms the second part of the index array
                            loc = [loc(1:j-1),i,loc(j:n-1)]
                            exit
                        endif
                    endif
                end do
            endif
        enddo
    end function maxnloc

    !> sort and return the 3 lowest
    function min3( rarr ) result(min_3)
        real, intent(in) :: rarr(:)
        real :: min_3(3)
        integer :: j, n
        real    :: ra
        n = size(rarr)
        if( n < 4)then
            min_3(:n) = rarr
            return
        end if
        min_3 = rarr(:3)
        call hpsort(min_3)
        do j=4,n
            if(rarr(j) < min_3(3))then
                ra = rarr(j)
                if(ra < min_3(2))then
                    if(ra < min_3(1))then
                        min_3 = (/ ra,  min_3(1), min_3(2) /)
                    else
                        min_3 = (/ min_3(1), ra, min_3(2) /)
                    end if
                else
                    min_3 = (/ min_3(1), min_3(2), ra /)
                end if
            end if
        end do
    end function min3

    function minnloc( rarr, n ) result( loc )
        real,    intent(in) :: rarr(:)
        integer, intent(in) :: n
        real    :: arr(n), val
        integer :: loc(n), i, j, sz
        logical :: val_lt(n)
        sz = size(rarr)
        if( sz < n )&
        &call simple_exception('cannot identify more maxima than elements in the array; minnloc', __FILENAME__ , __LINE__)
        loc = (/(i,i=1,n)/)
        arr = rarr(:n)
        call hpsort(arr, loc)
        if( sz == n ) return
        do i=n+1,sz
            val    = rarr(i)
            val_lt = val < rarr(loc)
            if( any(val_lt) )then
                do j=1,n
                    if( val_lt(j) )then
                        if( j == 1 )then
                            loc = [i,loc(1:n-1)]
                            exit
                        else if( j == n )then
                            loc(n) = i
                            exit
                        else
                            ! 1:j-1 untouched
                            ! insert i @ location j
                            ! j+1:n-1 forms the second part of the index array
                            loc = [loc(1:j-1),i,loc(j:n-1)]
                            exit
                        endif
                    endif
                end do
            endif
        enddo
    end function minnloc

    function peakfinder( vals ) result( peakpos )
        real,    intent(in)  :: vals(:)
        logical, allocatable :: peakpos(:)
        integer :: n, i
        n = size(vals)
        allocate(peakpos(n), source=.false.)
        if( vals(1) > vals(2) )  peakpos(1) = .true.
        do i=2,n-1
            if( vals(i) >= vals(i-1) .and. vals(i) >= vals(i+1) ) peakpos(i) = .true.
        end do
        if( vals(n) > vals(n-1) ) peakpos(n) = .true.
    end function peakfinder

    !>   reverses an integer array
    subroutine reverse_iarr( iarr )
        integer, intent(inout) :: iarr(:) !< array for modification
        integer                :: i, j, iswap, sz, en
        sz = size(iarr,1)
        if( sz < 2 )then
            return
        endif
        if( mod(sz,2) == 0 )then
            en = sz/2+1
        else
            en = (sz+1)/2+1
        endif
        j = 0
        do i = sz,en,-1
            j = j+1
            iswap   = iarr(j)
            iarr(j) = iarr(i)
            iarr(i) = iswap
        end do
    end subroutine reverse_iarr

    !>   reverses a real array
    subroutine reverse_rarr( rarr )
        real, intent(inout) :: rarr(:) !< array for modification
        integer             :: i, j, sz, en
        real                :: rswap
        sz = size(rarr,1)
        if( sz < 2 )then
            return
        endif
        if( mod(sz,2) == 0 )then
            en = sz/2+1
        else
            en = (sz+1)/2+1
        endif
        j = 0
        do i = sz,en,-1
            j = j+1
            rswap   = rarr(j)
            rarr(j) = rarr(i)
            rarr(i) = rswap
        end do
    end subroutine reverse_rarr

    !>   for selecting kth largest, array is modified
    real function selec(k,n,arr)
        integer, intent(in)    :: k,n
        real,    intent(inout) :: arr(:)
        integer :: i,ir,j,l,mid
        real    :: a,temp
        l = 1
        ir = n
    22  if (ir-l.le.1) then
            if (ir-1.eq.1) then
                if (arr(ir).lt.arr(l)) then
                    temp = arr(l)
                    arr(l) = arr(ir)
                    arr(ir) = temp
                endif
            endif
            selec = arr(k)
            return
        else
            mid = (l+ir)/2
            temp = arr(mid)
            arr(mid) = arr(l+1)
            arr(l+1) = temp
            if (arr(l).gt.arr(ir)) then
                temp = arr(l)
                arr(l) = arr(ir)
                arr(ir) = temp
            endif
            if (arr(l+1).gt.arr(ir)) then
                temp = arr(l+1)
                arr(l+1) = arr(ir)
                arr(ir) = temp
            endif
            if (arr(l).gt.arr(l+1)) then
                temp = arr(l)
                arr(l) = arr(l+1)
                arr(l+1) = temp
            endif
            i = l+1
            j = ir
            a = arr(l+1)
    23      continue
            i = i+1
            if (arr(i).lt.a) goto 23
    24      continue
            j = j-1
            if (arr(j).gt.a) goto 24
            if (j.lt.i)      goto 25
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            goto 23
    25      arr(l+1) = arr(j)
            arr(j) = a
            if (j.ge.k) ir = j-1
            if (j.le.k) l = i
        endif
        goto 22
    end function selec
    
    ! recursively sort the columns of a 2D table, w.r.t the sum of the best nrows/ncols
    ! entries of each column
    ! for example:
    !    original table:
    !                   r1          r2          r3
    !           ------------------------------------
    !             c1  |  2           0           6
    !             c2  |  4           0           7
    !             c3  |  3           1           0
    !             c4  |  3           0           2
    !             c5  |  0           4           7
    !             c6  |  6           4           2
    !   sorted table (smallest to largest):
    !                   r2          r1          r3
    !           ------------------------------------
    !             c1  |  0           2           6
    !             c4  |  0           3           2
    !             c3  |  1           3           0
    !             c6  |  4           6           2
    !             c5  |  4           0           7
    !             c2  |  0           4           7
    !
    ! based on this uniformly sorted table, one can cluster:
    !            (c2, c5) -> r3, (c6, c3) -> r1, (c1, c4) -> r2
    subroutine cluster_sort( nrows, ncols, tab, to_ii, to_ir, cur_id, cur_ir )
        integer, intent(in)    :: nrows, ncols
        real,    intent(in)    :: tab(nrows, ncols)
        integer, intent(in)    :: to_ii
        integer, intent(in)    :: to_ir
        integer, intent(inout) :: cur_id(nrows)
        integer, intent(inout) :: cur_ir(ncols)
        integer :: ir, ip, tmp_id(to_ii), tmp_i, orig_id(to_ii), ir_best, num
        real    :: best_sum, sum_prob
        num = int(nrows/ncols)
        if( to_ii <= num ) return
        best_sum = 0.
        orig_id  = cur_id(1:to_ii)
        ir_best  = to_ir
        do ir = 1,to_ir
            tmp_id = (/(ip, ip=1,to_ii)/)
            call hpsort_ind(tmp_id, tab(orig_id,cur_ir(ir)))
            sum_prob = sum(tab(orig_id(tmp_id(to_ii-num+1:to_ii)),cur_ir(ir)))
            if( sum_prob > best_sum )then
                best_sum        = sum_prob
                cur_id(1:to_ii) = orig_id(tmp_id)
                ir_best         = ir
            endif
        enddo
        ! swapping the last with the current best ir
        tmp_i           = cur_ir(to_ir)
        cur_ir(to_ir)   = cur_ir(ir_best)
        cur_ir(ir_best) = tmp_i
        call cluster_sort( nrows, ncols, tab, to_ii - num, to_ir - 1, cur_id, cur_ir )
    end subroutine cluster_sort

    subroutine uniform_sort( nrows, ncols, tab, to_ii, cur_id, cur_ir, mask_ir )
        integer,            intent(in)    :: nrows, ncols
        real,               intent(in)    :: tab(nrows, ncols)
        integer,            intent(in)    :: to_ii
        integer,            intent(inout) :: cur_id(nrows)
        integer,            intent(inout) :: cur_ir(nrows)
        logical,            intent(inout) :: mask_ir(ncols)
        integer :: ir, ip, tmp_i, max_ind_ir, max_ind_ip
        real    :: max_ir(ncols), max_ip(ncols)
        if( to_ii < 1 ) return
        if( .not.(any(mask_ir)) ) mask_ir = .true.
        do ir = 1, ncols
            if( mask_ir(ir) )then
                max_ir(ir) = 0.
                do ip = 1, to_ii
                    if( tab(cur_id(ip), ir) > max_ir(ir) )then
                        max_ir(ir) = tab(cur_id(ip), ir)
                        max_ip(ir) = ip
                    endif
                enddo
            endif
        enddo
        max_ind_ir = maxloc(max_ir, dim=1, mask=mask_ir)
        max_ind_ip = max_ip(max_ind_ir)
        ! swapping max_ind_ip
        tmp_i               = cur_id(to_ii)
        cur_id(to_ii)       = cur_id(max_ind_ip)
        cur_id(max_ind_ip)  = tmp_i
        cur_ir(to_ii)       = max_ind_ir
        mask_ir(max_ind_ir) = .false.
        call uniform_sort(nrows, ncols, tab, to_ii-1, cur_id, cur_ir, mask_ir)
    end subroutine uniform_sort

    ! sorting rarr, but only keep the sorted indeces
    subroutine hpsort_ind( iarr, rarr )
        integer, intent(inout) :: iarr(:)
        real,    intent(in)    :: rarr(:)
        integer :: i, ir, j, l, ia, n
        real    :: ra
        n = size(rarr)
        if( n /= size(iarr) )&
        &call simple_exception('nonconforming array sizes; hpsort_6', __FILENAME__ , __LINE__)
        if( n < 2 ) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ia = iarr(l)
                ra = rarr(ia)
            else
                ia = iarr(ir)
                ra = rarr(ia)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ia
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(iarr(j)) < rarr(iarr(j+1))) j = j+1
                endif
                if(ra < rarr(iarr(j)))then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ia
            end do
        end do
    end subroutine hpsort_ind
end module simple_srch_sort_loc
    
module simple_strategy_tree_helpers
use simple_rnd, only: ran3
implicit none

real,    parameter :: INVALID_CORR        = -huge(1.0)
real,    parameter :: INVALID_CORR_THRESH = INVALID_CORR / 2.0
integer, parameter :: MAX_NTREES          = 2500
integer, parameter :: MAX_NPEAKS          = 64
integer, parameter :: MAX_TREE_REFS       = 1024

contains    
    
    integer function choose_next_child_prob( left_idx, right_idx, corr_left, corr_right ) result(inode_next)
        integer, intent(in) :: left_idx, right_idx
        real,    intent(in) :: corr_left, corr_right
        real :: cmax
        real :: p_left
        real :: p_right
        inode_next = 0
        if( left_idx == 0 .and. right_idx == 0 ) return
        if( left_idx == 0 )then
            inode_next = right_idx
            return
        endif
        if( right_idx == 0 )then
            inode_next = left_idx
            return
        endif
        if( is_invalid_corr(corr_left) .and. is_invalid_corr(corr_right) )then
            if( sample_two(1.0, 1.0) == 1 )then
                inode_next = left_idx
            else
                inode_next = right_idx
            endif
            return
        endif
        cmax = max(corr_left, corr_right)
        if( is_invalid_corr(corr_left) )then
            p_left = 0.0
        else
            p_left = exp(corr_left - cmax)
        endif
        if( is_invalid_corr(corr_right) )then
            p_right = 0.0
        else
            p_right = exp(corr_right - cmax)
        endif
        if( sample_two(p_left, p_right) == 1 )then
            inode_next = left_idx
        else
            inode_next = right_idx
        endif
    end function choose_next_child_prob

    integer function choose_next_child_greedy( left_idx, right_idx, corr_left, corr_right ) result(inode_next)
        integer, intent(in) :: left_idx, right_idx
        real,    intent(in) :: corr_left, corr_right
        inode_next = 0
        if( left_idx == 0 .and. right_idx == 0 ) return
        if( left_idx == 0 )then
            inode_next = right_idx
            return
        endif
        if( right_idx == 0 )then
            inode_next = left_idx
            return
        endif
        if( is_invalid_corr(corr_left) .and. is_invalid_corr(corr_right) )then
            if( sample_two(1.0, 1.0) == 1 )then
                inode_next = left_idx
            else
                inode_next = right_idx
            endif
            return
        endif
        if( is_invalid_corr(corr_left) )then
            inode_next = right_idx
            return
        endif
        if( is_invalid_corr(corr_right) )then
            inode_next = left_idx
            return
        endif
        inode_next = merge(left_idx, right_idx, corr_left >= corr_right)
    end function choose_next_child_greedy

    logical pure function is_invalid_corr( corr ) result(invalid)
        real, intent(in) :: corr
        invalid = corr <= INVALID_CORR_THRESH
    end function is_invalid_corr

    integer function sample_two( p1, p2 ) result(which)
        real, intent(in) :: p1, p2
        real :: psum
        real :: r
        psum = p1 + p2
        if( psum <= 0.0 )then
            which = merge(1, 2, ran3() < 0.5)
            return
        endif
        r = ran3()
        which = merge(1, 2, r < p1 / psum)
    end function sample_two

end module simple_strategy_tree_helpers
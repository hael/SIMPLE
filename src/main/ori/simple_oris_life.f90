submodule (simple_oris) simple_oris_life
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module function constructor( n, is_ptcl ) result( self )
        integer, intent(in) :: n
        logical, intent(in) :: is_ptcl
        type(oris) :: self
        call self%new(n, is_ptcl)
    end function constructor

    module subroutine new_1( self, n, is_ptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        logical,     intent(in)    :: is_ptcl
        integer :: i
        call self%kill
        self%n = n
        allocate( self%o(self%n) )
        do i=1,n
            call self%o(i)%new(is_ptcl)
        end do
    end subroutine new_1

    module subroutine new_2( self, o_arr )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_arr(:)
        integer :: i
        call self%kill
        self%n = size(o_arr)
        allocate( self%o(self%n) )
        do i=1,self%n
            call self%o(i)%copy(o_arr(i))
        end do
    end subroutine new_2

    module subroutine reallocate( self, new_n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: new_n
        type(oris) :: tmp
        integer    :: old_n, i
        logical    :: is_ptcl
        if( self%n == 0 )     THROW_HARD('cannot reallocate non-existing oris; reallocate')
        if( new_n <= self%n ) THROW_HARD('reallocation to smaller size not supported; reallocate')
        is_ptcl = self%o(1)%is_particle()
        old_n = self%n
        tmp   = self
        call self%new(new_n, is_ptcl)
        do i=1,old_n
            self%o(i) = tmp%o(i)
        end do
        call tmp%kill
    end subroutine reallocate

    module function extract_subset_1( self, from, to ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: from, to
        type(oris) :: self_sub
        integer    :: n, cnt, i
        n  = to - from + 1
        call self_sub%new(n, self%is_particle())
        cnt = 0
        do i = from, to
            cnt = cnt + 1
            call self_sub%o(cnt)%copy(self%o(i))
        end do
    end function extract_subset_1

    module function extract_subset_2( self, inds ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: inds(:)
        type(oris) :: self_sub
        integer    :: n, i
        n  = size(inds)
        call self_sub%new(n, self%is_particle())
        do i = 1,n
            call self_sub%o(i)%copy(self%o(inds(i)))
        end do
    end function extract_subset_2

    module subroutine kill_chash( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( allocated(self%o) )then
            do i=1,self%n
                call self%o(i)%kill_chash
            end do
        endif
    end subroutine kill_chash

    module subroutine kill( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( allocated(self%o) )then
            do i=1,self%n
                call self%o(i)%kill
            end do
            deallocate( self%o )
        endif
        self%n = 0
    end subroutine kill

end submodule simple_oris_life
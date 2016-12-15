program permutations
implicit none
integer, parameter :: N = 20
integer, dimension (1 : N) :: permutation
 
  call permute (1)
 
contains

    recursive subroutine permute (pos)
        integer, intent (in) :: pos
        integer              :: value
        if( pos > N )then
            write(*,*) permutation
        else
        do value = 1, N
            if( .not. any (permutation(:pos-1) .eq. value) )then
                permutation(pos) = value
                call permute(pos+1)
            end if
        end do
    end if
 
  end subroutine
 
end program permutations
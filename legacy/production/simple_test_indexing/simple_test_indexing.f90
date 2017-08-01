program test_indexing
use simple_jiffys
use simple_timing
integer, parameter :: n_i=10
integer, parameter :: n_j=35
integer, parameter :: n_k=40
integer :: i, j, k, l, cnt
cnt = 0
do j=1,n_j
    do i=1,n_i
        cnt = cnt+1
        if( cnt /= gf2(j, i) ) stop 'ERROR, gf2'
    end do
end do
cnt = 0
do k=1,n_k
    do j=1,n_j
        do i=1,n_i
            cnt = cnt+1
            if( cnt /= gf3(k, j, i) ) stop 'ERROR, gf3'
        end do
    end do
end do
cnt = 0
do l=1,100
    do k=1,n_k
        do j=1,n_j
            do i=1,n_i
                cnt = cnt+1
                if( cnt /= gf4(l, k, j, i) ) stop 'ERROR, gf4'
            end do
        end do
    end do
end do

contains

    function gf1( i ) result( cnt )
        integer, intent(in) :: i
        integer :: cnt
        cnt = i
    end function
    
    function gf2( j, i ) result( cnt )
        integer, intent(in) :: j, i
        integer :: cnt
        cnt = (j-1)*n_i+gf1(i)
    end function
    
    function gf3( k, j, i ) result( cnt )
        integer, intent(in) :: k, j, i
        integer :: cnt
        cnt = (k-1)*n_i*n_j+gf2( j, i )
    end function
    
    function gf4( l, k, j, i ) result( cnt )
        integer, intent(in) :: l, k, j, i
        integer :: cnt
        cnt = (l-1)*n_i*n_j*n_k+gf3( k, j, i )
    end function
    
end program

program simple_test_shelliter


integer              :: hrange(2), krange(2), h, k, hl, hr, kl, kr
logical, allocatable :: finds_in_shell(:,:), binmap(:,:)

hrange(1) = -10
hrange(2) =  10
krange(1) = -10
krange(2) =  10

allocate( finds_in_shell(hrange(1):hrange(2),krange(1):krange(2)),&
          binmap(hrange(1):hrange(2),krange(1):krange(2)))
binmap = .false.
binmap(0,0) = .true.
do while( .not. all(binmap) )
    call identify_shell
    print *, '*************************'
    call print_shell
    where( finds_in_shell .eqv. .true. ) binmap = .true.
end do

contains

    subroutine identify_shell
        finds_in_shell = .false.
        do h=hrange(1),hrange(2)  
            hl = max(hrange(1),h-1)
            hr = min(hrange(2),h+1)
            do k=krange(1),krange(2)
                if( .not. binmap(h,k) )then
                    kl = max(krange(1),k-1)
                    kr = min(krange(2),k+1)
                    if( any(binmap(hl:hr,kl:kr)) ) finds_in_shell(h,k) = .true.
                end if
            end do
        end do
    end subroutine identify_shell

    subroutine print_shell
        do h=hrange(1),hrange(2) 
            do k=krange(1),krange(2)
                if( finds_in_shell(h,k) ) print *, h, k
            end do
        end do
    end subroutine print_shell

end program simple_test_shelliter

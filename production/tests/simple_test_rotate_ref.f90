program simple_test_rotate_ref
implicit none
integer, parameter :: NP = 100, NK = 300, NR = 200
complex :: ref(NP, NK), ref_rot(NP, NK), fast_ref_rot(NP, NK)
real    :: a(NP, NK), b(NP, NK), start, finish, norm_all, fast_all
integer :: i
call random_number(a)
call random_number(b)
ref      = cmplx(a,b)
norm_all = 0.
fast_all = 0.
do i = 1,NR
    call cpu_time(start)
    call rotate_ref(ref, i, ref_rot)
    call cpu_time(finish)
    norm_all = norm_all + (finish - start)
    call cpu_time(start)
    call fast_rotate_ref(ref, i, fast_ref_rot)
    call cpu_time(finish)
    fast_all = fast_all + (finish - start)
    if( .not.(all(ref_rot .eq. fast_ref_rot)) )then
        print *, 'FAILED'
        stop
    endif
enddo
print *, 'current timing = ', norm_all
print *, '-----------'
print *, 'improved timing = ', fast_all
print *, '-----------'
print *, 'PASSED'

contains

    subroutine rotate_ref( ref_in, irot, ref_rot_out )
        complex, intent(in)  :: ref_in(NP, NK)
        integer, intent(in)  :: irot
        complex, intent(out) :: ref_rot_out(NP, NK)
        integer :: rot, jrot
        do jrot = 1,NP
            rot = jrot - (irot - 1) ! reverse rotation
            if( rot < 1 ) rot = rot + NR
            if( rot > NP )then
                ref_rot_out(jrot,:) = conjg(ref_in(rot-NP,:))
            else
                ref_rot_out(jrot,:) = ref_in(rot,:)
            endif
        enddo
    end subroutine rotate_ref

    subroutine rotate_ref_2( ref_in, irot, ref_rot_out )
        complex, intent(in)  :: ref_in(NP, NK)
        integer, intent(in)  :: irot
        complex, intent(out) :: ref_rot_out(NP, NK)
        integer :: rot, jrot
        do jrot = 1,irot-1
            rot = jrot - (irot - 1) + NR ! reverse rotation
            if( rot > NP )then
                ref_rot_out(jrot,:) = conjg(ref_in(rot-NP,:))
            else
                ref_rot_out(jrot,:) = ref_in(rot,:)
            endif
        enddo
        do jrot = irot,NP
            rot = jrot - (irot - 1) ! reverse rotation
            if( rot > NP )then
                ref_rot_out(jrot,:) = conjg(ref_in(rot-NP,:))
            else
                ref_rot_out(jrot,:) = ref_in(rot,:)
            endif
        enddo
    end subroutine rotate_ref_2

    subroutine rotate_ref_t( ref_in, irot, ref_rot_out )
        complex, intent(in)  :: ref_in(NK, NP)
        integer, intent(in)  :: irot
        complex, intent(out) :: ref_rot_out(NK, NP)
        integer :: rot, jrot
        do jrot = 1,NP
            rot = jrot - (irot - 1) ! reverse rotation
            if( rot < 1 ) rot = rot + NR
            if( rot > NP )then
                ref_rot_out(:,jrot) = conjg(ref_in(:,rot-NP))
            else
                ref_rot_out(:,jrot) = ref_in(:,rot)
            endif
        enddo
    end subroutine rotate_ref_t

    subroutine fast_rotate_ref( ref_in, irot, ref_rot_out )
        complex, intent(in)  :: ref_in(NP, NK)
        integer, intent(in)  :: irot
        complex, intent(out) :: ref_rot_out(NP, NK)
        integer :: mid
        if( irot == 1 )then
            ref_rot_out = ref_in
        elseif( irot >= 2 .and. irot <= NP )then
            mid = NP - irot + 1
            ref_rot_out(   1:irot-1,:) = conjg(ref_in(mid+1:NP, :))
            ref_rot_out(irot:NP,    :) =       ref_in(    1:mid,:)
        elseif( irot == NP + 1 )then
            ref_rot_out = conjg(ref_in)
        else
            mid = NR - irot + 1
            ref_rot_out(irot-NP:NP,       :)  = conjg(ref_in(    1:mid,:))
            ref_rot_out(      1:irot-NP-1,:) =        ref_in(mid+1:NP, :)
        endif
    end subroutine fast_rotate_ref

end program simple_test_rotate_ref
    

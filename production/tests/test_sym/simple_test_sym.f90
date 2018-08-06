program simple_test_ctf
include 'simple_lib.f08'
use simple_sym,   only: sym
use simple_oris,  only: oris
use simple_ori,   only: ori
use simple_timer
implicit none
type(sym)  :: symop
type(ori)  :: o
type(oris) :: os, os_c1
character(len=3), allocatable :: sym_subgrps(:)
integer :: i, j
! symmetry subgroups
write(*,'(A)')'>>> SYMMETRY SUBGROUPS FOR C1'
symop = sym('c1')
call write_subgrp(symop)
write(*,'(A)')'>>> SYMMETRY SUBGROPUS FOR D14'
symop = sym('d14')
call write_subgrp(symop)
write(*,'(A)')'>>> SYMMETRY SUBGROPUS FOR C17'
symop = sym('c17')
call write_subgrp(symop)
write(*,'(A)')'>>> SYMMETRY SUBGROPUS FOR C32'
symop = sym('c32')
call write_subgrp(symop)
write(*,'(A)')'>>> SYMMETRY SUBGROPUS FOR C34'
symop = sym('c34')
call write_subgrp(symop)
write(*,'(A)')'>>> SYMMETRY SUBGROPUS FOR D8'
symop = sym('d8')
call write_subgrp(symop)
write(*,'(A)')'>>> SYMMETRY SUBGROPUS FOR O'
symop = sym('O')
call write_subgrp(symop)

symop = sym('D2')
call write_subgrp(symop)
o = symop%get_symori(1)
call o%print_ori
o = symop%get_symori(2)
call o%print_ori
o = symop%get_symori(3)
call o%print_ori
o = symop%get_symori(4)
call o%print_ori

call symspiral('c3')
call symspiral('t')
call symspiral('o')
call symspiral('i')

contains

    subroutine write_subgrp( se )
        class(sym), intent(inout) :: se
        type(sym) :: symop, tmp
        integer   :: i
        do i=1, se%get_nsubgrp()
            tmp = se%get_subgrp(i)
            write(*,'(I3,1X,A3)') i, tmp%get_pgrp()
        enddo
    end subroutine

    subroutine symspiral( pgrp )
        character(len=*) :: pgrp
        type(oris) :: spiral, tmp
        type(sym)  :: symop
        integer :: i,j
        symop = sym(pgrp, incl_mirror=.false.)
        call tmp%new(200)
        call spiral%new(200)
        call spiral%spiral(symop%get_nsym(), symop%srchrange())
        call spiral%write(pgrp//'_oris.txt')
        call spiral%write2bild(pgrp//'.bild')
        call spiral%mirror2d
        call spiral%write(pgrp//'_mirr_oris.txt')
        call spiral%write2bild(pgrp//'_mirr.bild')
        do i=1,symop%get_nsym()
            do j=1,spiral%get_noris()
                o = symop%apply(spiral%get_ori(j),i)
                call tmp%set_ori(j,o)
            enddo
            call tmp%write2bild(pgrp//'_mirr_'//int2str(i)//'.bild')
        enddo
        call tmp%kill
        call symop%kill
        call spiral%kill
    end subroutine symspiral

end program simple_test_ctf

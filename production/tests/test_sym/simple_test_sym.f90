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

call symspiral('c1')
call symspiral('d7')
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
        type(oris) :: spiral
        type(sym)  :: symop
        integer :: i,j
        symop = sym(pgrp)
        call spiral%new(1000)
        call symop%build_refspiral(spiral)
        call spiral%write2bild(pgrp//'.bild')
        call symop%kill
        call spiral%kill
    end subroutine symspiral

end program simple_test_ctf

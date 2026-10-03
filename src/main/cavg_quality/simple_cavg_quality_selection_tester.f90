!@descr: unit tests for writing the selected/rejected class-average stacks (simple_cavg_quality_selection)
module simple_cavg_quality_selection_tester
use simple_test_utils
use simple_string,                 only: string
use simple_string_utils,           only: int2str
use simple_fileio,                 only: file_exists, del_file
use simple_syslib,                 only: get_process_id
use simple_imghead,                only: find_ldim_nptcls
use simple_image,                  only: image
use simple_cavg_quality_selection, only: write_cavg_selection_stacks
implicit none
private
public :: run_all_cavg_quality_selection_tests

integer, parameter :: BOX = 8, NCLS = 4

contains

    subroutine run_all_cavg_quality_selection_tests()
        write(*,'(A)') '**** running all class-average selection tests ****'
        call test_selection_stacks_split_and_order()
        call test_selection_stacks_replace_existing()
    end subroutine run_all_cavg_quality_selection_tests

    !> classes 1 and 3 selected: the selected stack holds them in order, the rejected stack 2 and 4
    subroutine test_selection_stacks_split_and_order()
        type(image)  :: imgs(NCLS)
        type(string) :: fsel, frej
        write(*,'(A)') 'test_selection_stacks_split_and_order'
        call make_imgs(imgs)
        call stack_names(fsel, frej)
        call write_cavg_selection_stacks(imgs, [1, 0, 1, 0], fsel, frej)
        call check_stack(fsel, [1., 3.], 'selected')
        call check_stack(frej, [2., 4.], 'rejected')
        call cleanup(imgs, fsel, frej)
    end subroutine test_selection_stacks_split_and_order

    !> a second call replaces the stacks instead of appending to them
    subroutine test_selection_stacks_replace_existing()
        type(image)  :: imgs(NCLS)
        type(string) :: fsel, frej
        write(*,'(A)') 'test_selection_stacks_replace_existing'
        call make_imgs(imgs)
        call stack_names(fsel, frej)
        call write_cavg_selection_stacks(imgs, [1, 1, 1, 0], fsel, frej)
        call write_cavg_selection_stacks(imgs, [0, 1, 0, 0], fsel, frej)
        call check_stack(fsel, [2.],         'selected after rewrite')
        call check_stack(frej, [1., 3., 4.], 'rejected after rewrite')
        call cleanup(imgs, fsel, frej)
    end subroutine test_selection_stacks_replace_existing

    ! class i is a constant image of value i
    subroutine make_imgs( imgs )
        type(image), intent(inout) :: imgs(NCLS)
        integer :: icls
        do icls = 1,NCLS
            call imgs(icls)%new([BOX,BOX,1], 1.0, wthreads=.false.)
            imgs(icls) = real(icls)
        enddo
    end subroutine make_imgs

    subroutine stack_names( fsel, frej )
        type(string), intent(inout) :: fsel, frej
        fsel = 'cavg_selection_test_'//int2str(get_process_id())//'_selected.mrc'
        frej = 'cavg_selection_test_'//int2str(get_process_id())//'_rejected.mrc'
    end subroutine stack_names

    ! the stack holds one image per expected value, in order
    subroutine check_stack( fname, values, label )
        type(string),     intent(in) :: fname
        real,             intent(in) :: values(:)
        character(len=*), intent(in) :: label
        type(image) :: img
        integer     :: ldim(3), n, i
        call assert_true(file_exists(fname), label//' stack is written')
        if( .not. file_exists(fname) ) return
        call find_ldim_nptcls(fname, ldim, n)
        call assert_int(size(values), n, label//' stack size')
        if( n /= size(values) ) return
        call img%new([BOX,BOX,1], 1.0, wthreads=.false.)
        do i = 1,n
            call img%read(fname, i)
            call assert_real(values(i), img%get_rmat_at(1,1,1), 1.e-6, label//' image '//int2str(i))
        enddo
        call img%kill
    end subroutine check_stack

    subroutine cleanup( imgs, fsel, frej )
        type(image),  intent(inout) :: imgs(NCLS)
        type(string), intent(in)    :: fsel, frej
        integer :: icls
        do icls = 1,NCLS
            call imgs(icls)%kill
        enddo
        if( file_exists(fsel) ) call del_file(fsel)
        if( file_exists(frej) ) call del_file(frej)
    end subroutine cleanup

end module simple_cavg_quality_selection_tester

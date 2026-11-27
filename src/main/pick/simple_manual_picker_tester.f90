!@descr: manual-picking tests for selection, box resizing, lower-origin viewport geometry and Python-compatible box output
module simple_manual_picker_tester
use simple_defs,          only: dp, logfhandle
use simple_manual_picker, only: manual_picker
use simple_test_utils,    only: assert_true, assert_int, enter_fixture, leave_fixture, tests_failed
use simple_string,        only: string
use simple_image,         only: image
use simple_stack_display, only: stack_display
implicit none
private
public :: run_all_manual_picker_tests

contains

    subroutine run_all_manual_picker_tests()
        write(logfhandle,'(a)') '**** running all manual picker tests ****'
        call test_selection_and_resize
        call test_viewport
        call test_box_output
        call test_display_percentiles
    end subroutine run_all_manual_picker_tests

    subroutine test_selection_and_resize()
        class(manual_picker), allocatable :: p
        real(dp), allocatable :: xy(:,:)
        write(logfhandle,'(a)') 'test_selection_and_resize'
        allocate(p)
        call p%new([101,79], 11)
        call p%click([0._dp,0._dp])
        xy = p%get_centres()
        call assert_true(all(xy(:,1) == [5.5_dp,5.5_dp]), 'a lower-left pick keeps the odd box inside the image')
        call p%click([40._dp,30._dp])
        call p%click([101._dp,79._dp])
        xy = p%get_centres()
        call assert_true(all(xy(:,3) == [95.5_dp,73.5_dp]), 'a top-right pick is clamped to the non-square image')
        call p%click([40._dp,30._dp])
        call assert_int(3, p%count_boxes(), 'clicking an existing box selects it without duplication')
        call assert_int(2, p%get_selected(), 'the middle box is selected')
        call p%move_selected([95._dp,73._dp])
        call p%click([95._dp,73._dp])
        call assert_int(3, p%get_selected(), 'the last-added overlapping box is selected first')
        call p%delete_selected()
        call assert_int(2, p%count_boxes(), 'deletion removes only the selected box')
        call assert_int(0, p%get_selected(), 'deletion clears selection')
        call p%set_box(40)
        xy = p%get_centres()
        call assert_true(all(xy(:,1) == [20._dp,20._dp]), 'slider enlargement moves a lower-edge box inward')
        call assert_true(all(xy(:,2) == [81._dp,59._dp]), 'slider enlargement moves an upper-edge box inward')
        call p%set_box(5)
        call assert_true(all(p%get_centres() == xy), 'slider shrink preserves centres and insertion order')
        call assert_int(5, p%get_box(), 'odd box sizes survive slider changes')
        call p%set_box(999)
        call assert_int(79, p%get_box(), 'box size cannot exceed the shorter image dimension')
        call p%click([-1._dp,20._dp])
        call assert_int(2, p%count_boxes(), 'clicking the backdrop does not add a box')
        call p%kill
        call p%kill
        call assert_int(0, p%count_boxes(), 'kill is idempotent')
    end subroutine test_selection_and_resize

    subroutine test_viewport()
        class(manual_picker), allocatable :: p
        real(dp) :: screen(2), xy(2), anchor(2)
        real(dp), parameter :: tol = 1.e-10_dp
        write(logfhandle,'(a)') 'test_viewport'
        allocate(p)
        call p%new([100,60],10)
        call p%fit(400,300)
        screen = p%image_to_screen([0._dp,60._dp])
        call assert_true(all(abs(screen-[0._dp,30._dp]) < tol), 'fit centres a 400x240 image in a 400x300 viewport')
        screen = p%image_to_screen([25._dp,10._dp])
        call assert_true(all(abs(screen-[100._dp,230._dp]) < tol), 'lower-origin y is flipped exactly once for GTK')
        xy = p%screen_to_image(screen)
        call assert_true(all(abs(xy-[25._dp,10._dp]) < tol), 'screen picking recovers image coordinates')
        anchor = p%screen_to_image([123._dp,187._dp])
        call p%zoom([123._dp,187._dp],1.2_dp)
        xy = p%screen_to_image([123._dp,187._dp])
        call assert_true(all(abs(xy-anchor) < tol), 'wheel zoom keeps the image coordinate under the pointer fixed')
        screen = p%image_to_screen([25._dp,10._dp])
        call p%pan([17._dp,-23._dp])
        call assert_true(all(abs(p%image_to_screen([25._dp,10._dp])-screen-[17._dp,-23._dp]) < tol), &
            &'right drag pans by screen pixels independently of zoom')
        call p%kill
    end subroutine test_viewport

    subroutine test_box_output()
        class(manual_picker), allocatable :: p
        type(string) :: cwd, root
        integer :: nfail0, status, unit, values(4)
        character(len=1024) :: message
        write(logfhandle,'(a)') 'test_box_output'
        nfail0 = tests_failed
        call enter_fixture('manual_picker', cwd, root)
        allocate(p)
        call p%new([101,79],11)
        call p%click([16._dp,17._dp])
        call p%click([80.5_dp,60.5_dp])
        call p%save('positions.box',status,message)
        call assert_int(0,status,'box save succeeds')
        if( status == 0 )then
            open(newunit=unit,file='positions.box',status='old',action='read')
            read(unit,*,iostat=status) values
            call assert_int(0,status,'first box row has four integer fields')
            call assert_true(all(values == [10,12,11,11]), 'Python ties-to-even rounding gives lower-left (10,12), not centre coordinates')
            read(unit,*,iostat=status) values
            call assert_true(all(values == [75,55,11,11]), 'box output keeps insertion order and the MRC y direction')
            read(unit,*,iostat=status) values
            call assert_true(status < 0,'one box produces exactly one row')
            close(unit)
        endif
        call p%save('missing/positions.box',status,message)
        call assert_true(status /= 0 .and. len_trim(message) > 0,'an unwritable destination returns an actionable error')
        call assert_int(2,p%count_boxes(),'failed save keeps every pick in memory')
        call p%kill
        call p%new([32,31],5)
        call p%save('positions.box',status,message)
        open(newunit=unit,file='positions.box',status='old',action='read')
        read(unit,*,iostat=status) values
        call assert_true(status < 0,'saving zero picks replaces stale rows with an empty file, as the Python picker does')
        close(unit)
        call p%kill
        call leave_fixture(cwd,root,nfail0)
    end subroutine test_box_output

    subroutine test_display_percentiles()
        type(image) :: imgs(1)
        class(stack_display), allocatable :: disp
        real :: ramp(10,10,1)
        real(dp) :: bounds(2)
        integer :: i
        write(logfhandle,'(a)') 'test_display_percentiles'
        allocate(disp)
        ramp = reshape([(real(i),i=0,99)],[10,10,1])
        call imgs(1)%new([10,10,1],1.)
        call imgs(1)%set_rmat(ramp,.false.)
        call disp%new(imgs,[.true.],auto_percentiles=[1._dp,99._dp])
        bounds = disp%get_window()
        call assert_true(all(abs(bounds-[0.99_dp,98.01_dp]) < 1.e-10_dp), &
            &'picker contrast starts at NumPy linear 1st/99th percentiles of a 0..99 ramp')
        call disp%kill
        call imgs(1)%kill
    end subroutine test_display_percentiles

end module simple_manual_picker_tester

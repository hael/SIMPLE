!@descr: unit tests for simple_epu_xml: the beam shift of EPU metadata files, written in the layouts EPU uses
module simple_epu_xml_tester
use simple_epu_xml,    only: read_epu_beamshift
use simple_test_utils, only: assert_true, assert_false, assert_real
implicit none
private
public :: run_all_epu_xml_tests

character(len=*), parameter :: FNAME = 'epu_xml_tester.xml'
real,             parameter :: TOL   = 1.e-6

contains

    subroutine run_all_epu_xml_tests()
        write(*,'(A)') '**** running all EPU XML tests ****'
        call test_single_line_document()
        call test_indented_document()
        call test_longer_name_first()
        call test_missing_beamshift()
        call test_missing_component()
        call test_missing_file()
    end subroutine run_all_epu_xml_tests

    !> EPU's usual layout: the whole document on one line, ImageShift (also a:_x and a:_y) before
    !! BeamShift and BeamTilt after it
    subroutine test_single_line_document()
        real    :: x, y
        logical :: l_ok
        write(*,'(A)') 'test_single_line_document'
        call write_doc([character(len=700) :: &
            &'<?xml version="1.0" encoding="utf-8"?><MicroscopeImage xmlns="http://schemas.datacontract.org/2004/07/'//&
            &'Fei.SharedObjects"><microscopeData><optics><ImageShift xmlns:a="http://schemas.datacontract.org/2004/07/'//&
            &'Fei.Types"><a:_x>0.75</a:_x><a:_y>-0.5</a:_y></ImageShift><BeamShift xmlns:a="http://schemas.datacontract'//&
            &'.org/2004/07/Fei.Types"><a:_x>-0.012345</a:_x><a:_y>0.006789</a:_y></BeamShift><BeamTilt xmlns:a="x">'//&
            &'<a:_x>9</a:_x><a:_y>9</a:_y></BeamTilt></optics></microscopeData></MicroscopeImage>'])
        call read_epu_beamshift(FNAME, x, y, l_ok)
        call assert_true(l_ok, 'single line: BeamShift found')
        call assert_real(-0.012345, x, TOL, 'single line: x is BeamShift''s, not ImageShift''s')
        call assert_real( 0.006789, y, TOL, 'single line: y')
        call delete_doc
    end subroutine test_single_line_document

    !> an indented document with line breaks around the values and exponent notation
    subroutine test_indented_document()
        real    :: x, y
        logical :: l_ok
        write(*,'(A)') 'test_indented_document'
        call write_doc([character(len=60) :: '<MicroscopeImage>', '  <BeamShift xmlns:a="urn:fei">', &
            &'    <a:_x>', '      1.25E-03', '    </a:_x>', '    <a:_y>-2.5e-4</a:_y>', '  </BeamShift>', &
            &'</MicroscopeImage>'])
        call read_epu_beamshift(FNAME, x, y, l_ok)
        call assert_true(l_ok, 'indented: BeamShift found')
        call assert_real( 1.25e-3, x, TOL, 'indented: x across line breaks')
        call assert_real(-2.5e-4,  y, TOL, 'indented: y in exponent notation')
        call delete_doc
    end subroutine test_indented_document

    !> an element whose name begins with BeamShift is not BeamShift
    subroutine test_longer_name_first()
        real    :: x, y
        logical :: l_ok
        write(*,'(A)') 'test_longer_name_first'
        call write_doc([character(len=200) :: '<doc><BeamShiftCalibration><a:_x>5</a:_x><a:_y>5</a:_y>'//&
            &'</BeamShiftCalibration><BeamShift><a:_x>0.1</a:_x><a:_y>0.2</a:_y></BeamShift></doc>'])
        call read_epu_beamshift(FNAME, x, y, l_ok)
        call assert_true(l_ok, 'longer name: BeamShift found')
        call assert_real(0.1, x, TOL, 'longer name: x')
        call assert_real(0.2, y, TOL, 'longer name: y')
        call delete_doc
    end subroutine test_longer_name_first

    subroutine test_missing_beamshift()
        real    :: x, y
        logical :: l_ok
        write(*,'(A)') 'test_missing_beamshift'
        call write_doc([character(len=80) :: '<doc><ImageShift><a:_x>1</a:_x><a:_y>2</a:_y></ImageShift></doc>'])
        call read_epu_beamshift(FNAME, x, y, l_ok)
        call assert_false(l_ok, 'no BeamShift: not found')
        call assert_true(x == 0. .and. y == 0., 'no BeamShift: shifts zero')
        call delete_doc
    end subroutine test_missing_beamshift

    subroutine test_missing_component()
        real    :: x, y
        logical :: l_ok
        write(*,'(A)') 'test_missing_component'
        call write_doc([character(len=60) :: '<doc><BeamShift><a:_x>0.3</a:_x></BeamShift></doc>'])
        call read_epu_beamshift(FNAME, x, y, l_ok)
        call assert_false(l_ok, 'no a:_y: not found')
        call assert_true(x == 0. .and. y == 0., 'no a:_y: shifts zero')
        call delete_doc
    end subroutine test_missing_component

    subroutine test_missing_file()
        real    :: x, y
        logical :: l_ok
        write(*,'(A)') 'test_missing_file'
        call read_epu_beamshift('epu_xml_tester_absent.xml', x, y, l_ok)
        call assert_false(l_ok, 'missing file: not found')
    end subroutine test_missing_file

    subroutine write_doc( lines )
        character(len=*), intent(in) :: lines(:)
        integer :: funit, i
        open(newunit=funit, file=FNAME, status='replace', action='write')
        do i = 1, size(lines)
            write(funit,'(A)') trim(lines(i))
        end do
        close(funit)
    end subroutine write_doc

    subroutine delete_doc
        integer :: funit
        open(newunit=funit, file=FNAME, status='old')
        close(funit, status='delete')
    end subroutine delete_doc

end module simple_epu_xml_tester

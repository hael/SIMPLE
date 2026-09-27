!@descr: library tests of the abinitio3D_addon report's docking branch (simple_abinitio3D_addon_report): a moved frame is docked and its rotation recovered
! An asymmetric map of four blobs against its copy rotated by 90 degrees about z
! ((x,y) -> (-y,x) about the box centre): the in-frame correlation falls below
! the docking floor, the report docks the copy onto the map and recovers the
! rotation. The docking search takes seconds, so this is a library sub-suite
! (lib_reconstruction); the rest of the report is unit-tested in
! simple_abinitio3D_addon_report_tester.
module simple_abinitio3D_addon_dock_tester
use simple_core_module_api
use simple_image,                   only: image
use simple_abinitio3D_addon_report, only: abinitio3D_addon_report
use simple_test_utils
implicit none

public :: run_all_abinitio3D_addon_dock_tests
private
#include "simple_local_flags.inc"

real,             parameter :: SMPD    = 2.
real,             parameter :: MSKDIAM = 100.
integer,          parameter :: BOX     = 64
integer,          parameter :: SEED    = 20260927
character(len=*), parameter :: F_MAP   = 'tmp_addon_dock_map.mrc'
character(len=*), parameter :: F_ROT   = 'tmp_addon_dock_rot.mrc'

contains

    subroutine run_all_abinitio3D_addon_dock_tests()
        write(logfhandle,'(A)') '**** running all abinitio3D add-on docking tests ****'
        call write_fixtures()
        call test_rotated_frame()
        call del_file(F_MAP)
        call del_file(F_ROT)
    end subroutine run_all_abinitio3D_addon_dock_tests

    !> the map and its copy rotated by 90 degrees about z
    subroutine write_fixtures()
        real, parameter :: CENTRES(3,4) = reshape([-10.,4.,2., 9.,-6.,0., 2.,11.,-8., -3.,-9.,7.], [3,4])
        real, allocatable :: rmat(:,:,:), rot(:,:,:), noise(:,:,:)
        real    :: d2
        integer :: i, j, k, b, c
        call set_fixed_seed(SEED)
        allocate(rmat(BOX,BOX,BOX), rot(BOX,BOX,BOX), source=0.)
        allocate(noise(BOX,BOX,BOX))
        call random_number(noise)
        c = BOX/2 + 1
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, 4
                        d2 = (real(i-c) - CENTRES(1,b))**2 + (real(j-c) - CENTRES(2,b))**2 + (real(k-c) - CENTRES(3,b))**2
                        rmat(i,j,k) = rmat(i,j,k) + real(b) * exp(-0.5 * d2 / 9.)
                    enddo
                enddo
            enddo
        enddo
        rmat = rmat + 0.02 * (noise - 0.5)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    ! rot(x,y) = map(y,-x): index 2c-i leaves the box only at i=1
                    if( 2*c - i <= BOX ) rot(i,j,k) = rmat(j, 2*c - i, k)
                enddo
            enddo
        enddo
        call write_map(F_MAP, rmat)
        call write_map(F_ROT, rot)
    end subroutine write_fixtures

    subroutine write_map( fname, rmat )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: rmat(:,:,:)
        type(image) :: vol
        call vol%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(fname), del_if_exists=.true.)
        call vol%kill
    end subroutine write_map

    subroutine test_rotated_frame()
        type(abinitio3D_addon_report) :: report
        write(logfhandle,'(A)') '-- test_rotated_frame'
        call report%new(1, 8, SMPD, MSKDIAM)
        call report%compare_maps(1, string(F_MAP), string(F_ROT))
        call assert_true(report%get_corr(1) < 0.9, 'a rotated map does not correlate in the base frame')
        call assert_real(90., report%get_dock_angle(1), 10., 'docking recovers the 90 degree rotation')
        call report%kill
    end subroutine test_rotated_frame

end module simple_abinitio3D_addon_dock_tester

!@descr: for all masks tests
module simple_commanders_test_masks
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_msk_routines
  contains
    procedure :: execute      => exec_test_msk_routines
end type commander_test_msk_routines

type, extends(commander_base) :: commander_test_nano_mask
  contains
    procedure :: execute      => exec_test_nano_mask
end type commander_test_nano_mask

type, extends(commander_base) :: commander_test_score_volume_shape
    contains
        procedure :: execute      => exec_volume_shape_descriptors
end type commander_test_score_volume_shape

contains

subroutine exec_volume_shape_descriptors( self, cline )
    use simple_image,     only: image
    use simple_image_bin, only: image_bin
    use simple_imghead,   only: find_img_smpd
    class(commander_test_score_volume_shape), intent(inout) :: self
    class(cmdline),                             intent(inout) :: cline
    type(parameters) :: p
    type(image)      :: vol
    type(image_bin)  :: mskvol
    integer          :: ldim(3), nimgs
    real             :: header_smpd
    if( command_argument_count() < 2 )then
        write(logfhandle,'(a)') 'Usage: simple_test_exec test=score_volume_shape vol1=volume.mrc [smpd=xx] [lp=20.0]'
        stop
    endif
    call cline%parse_oldschool
    if( .not. cline%defined('lp') ) call cline%set('lp', 15.0)
    if( .not. cline%defined('smpd') )then
        header_smpd = find_img_smpd(cline%get_carg('vol1'))
        if( header_smpd <= 0.0 ) THROW_HARD('Could not determine smpd from volume header')
        call cline%set('smpd', header_smpd)
        write(logfhandle,'(A,F6.3,A)') '>>> PIXEL SIZE: ', header_smpd, ' A'
    endif
    call cline%checkvar('vol1', 1)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%vols(1), ldim, nimgs)
    call vol%new(ldim, p%smpd)
    call vol%read(p%vols(1))
    write(logfhandle,'(A,F9.1,A)') '>>> MSKDIAM: ', p%mskdiam, ' A'
    call mskvol%vol_shape_descr(vol, p%lp, p%msk)
    call mskvol%kill_bimg
    call vol%kill
    call simple_end('**** SIMPLE_TEST_VOLUME_SHAPE_DESCRIPTORS NORMAL STOP ****')
end subroutine exec_volume_shape_descriptors

!> the real-space mask routines inside OpenMP loops give the same result as serially.
!! The mask coordinates are memoised in module variables keyed on the box; the
!! routines THROW inside a parallel region if the memo does not match, so the
!! threaded path (memoise once outside, mask many inside) is a contract of its own.
!! The single-thread semantics are pinned by the `masks` sub-suite of unit_image;
!! this case needs threads (library tier, lib_masks) and asserts parallel == serial.
subroutine exec_test_msk_routines( self, cline )
    !$ use omp_lib
    use simple_test_utils,   only: begin_test_suite, end_test_suite, assert_true, assert_real, report_summary
    use simple_image,        only: image, unmemoize_mask_coords
    class(commander_test_msk_routines), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter :: NIMGS = 8, BOX2 = 128, BOX3 = 48
    real,    parameter :: SMPD = 1.0
    type(image), allocatable :: stk(:), ref(:)
    integer :: i, nthr
    logical :: test_failed
    call begin_test_suite('mask routines, parallel equals serial')
    nthr = 1
    !$ nthr = omp_get_max_threads()
    write(logfhandle,'(A,I3,A)') 'running with ', nthr, ' OpenMP threads (OMP_NUM_THREADS)'
    call assert_true(nthr > 1, 'the threaded path is exercised with more than one thread')
    allocate(stk(NIMGS), ref(NIMGS))
    call check_2d('mask2D_soft',    1)
    call check_2d('mask2D_softavg', 2)
    call check_2d('mask2D_hard',    3)
    call check_3d('mask3D_soft',    1)
    call check_3d('mask3D_softavg', 2)
    call check_3d('mask3D_hard',    3)
    do i = 1,NIMGS
        call stk(i)%kill
        call ref(i)%kill
    end do
    deallocate(stk, ref)
    call unmemoize_mask_coords
    call end_test_suite
    call report_summary(failed=test_failed)
    if( test_failed ) error stop 1
    call simple_end('**** SIMPLE_TEST_MSK_ROUTINES NORMAL STOP ****')

    contains

        subroutine fill( imgs, ldim )
            type(image), intent(inout) :: imgs(:)
            integer,     intent(in)    :: ldim(3)
            integer :: j
            do j = 1,size(imgs)
                call imgs(j)%new(ldim, SMPD)
                call imgs(j)%gauran(0.0, 1.0)
            end do
        end subroutine fill

        subroutine mask_2d( img, which, mskrad )
            type(image),      intent(inout) :: img
            integer,          intent(in)    :: which
            real,             intent(in)    :: mskrad
            select case(which)
                case(1); call img%mask2D_soft(mskrad)
                case(2); call img%mask2D_softavg(mskrad)
                case(3); call img%mask2D_hard(mskrad)
            end select
        end subroutine mask_2d

        subroutine mask_3d( img, which, mskrad )
            type(image),      intent(inout) :: img
            integer,          intent(in)    :: which
            real,             intent(in)    :: mskrad
            select case(which)
                case(1); call img%mask3D_soft(mskrad)
                case(2); call img%mask3D_softavg(mskrad)
                case(3); call img%mask3D_hard(mskrad)
            end select
        end subroutine mask_3d

        subroutine check_2d( name, which )
            character(len=*), intent(in) :: name
            integer,          intent(in) :: which
            real    :: mskrad, maxdiff
            integer :: j
            mskrad = real(BOX2)/3.0
            call unmemoize_mask_coords
            call fill(stk, [BOX2,BOX2,1])
            do j = 1,NIMGS
                call ref(j)%copy(stk(j))
                call mask_2d(ref(j), which, mskrad)   ! serial reference
            end do
            call unmemoize_mask_coords
            call stk(1)%memoize_mask_coords          ! once, outside the parallel region
            !$omp parallel do default(shared) private(j) proc_bind(close) schedule(static)
            do j = 1,NIMGS
                call mask_2d(stk(j), which, mskrad)
            end do
            !$omp end parallel do
            maxdiff = 0.0
            do j = 1,NIMGS
                maxdiff = max(maxdiff, maxval(abs(stk(j)%get_rmat() - ref(j)%get_rmat())))
            end do
            call assert_real(0.0, maxdiff, 1.0e-6, name//': parallel result equals serial')
        end subroutine check_2d

        subroutine check_3d( name, which )
            character(len=*), intent(in) :: name
            integer,          intent(in) :: which
            real    :: mskrad, maxdiff
            integer :: j
            mskrad = real(BOX3)/3.0
            call unmemoize_mask_coords
            call fill(stk, [BOX3,BOX3,BOX3])
            do j = 1,NIMGS
                call ref(j)%copy(stk(j))
                call mask_3d(ref(j), which, mskrad)
            end do
            call unmemoize_mask_coords
            call stk(1)%memoize_mask_coords
            !$omp parallel do default(shared) private(j) proc_bind(close) schedule(static)
            do j = 1,NIMGS
                call mask_3d(stk(j), which, mskrad)
            end do
            !$omp end parallel do
            maxdiff = 0.0
            do j = 1,NIMGS
                maxdiff = max(maxdiff, maxval(abs(stk(j)%get_rmat() - ref(j)%get_rmat())))
            end do
            call assert_real(0.0, maxdiff, 1.0e-6, name//': parallel result equals serial')
        end subroutine check_3d

end subroutine exec_test_msk_routines

subroutine exec_test_nano_mask( self, cline )
    use simple_image,     only: image
    use simple_image_msk, only: automask2D
    use simple_parameters, only: parameters
    class(commander_test_nano_mask),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    ! constants
    character(len=*), parameter :: DEFAULT_STK='selected.spi'
    real,             parameter :: DEFAULT_SMPD=0.358, DEFAULT_MSKDIAM=100.
    integer,          parameter :: NGROW=3, WINSZ=1, EDGE=12
    ! variables
    type(parameters)            :: params
    type(image),    allocatable :: imgs(:)
    real,           allocatable :: diams(:), shifts(:,:)
    integer                     ::  n, i, ldim(3)
    ! setup parameters
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)') 'Usage: simple_test_exec test=nano_mask stk=<images> smpd=<pixel size> mskdiam=<diameter>'
        write(logfhandle,'(a)') 'No input keywords provided; running the default test with selected.spi'
        call cline%set('stk',      DEFAULT_STK)
        call cline%set('smpd',     DEFAULT_SMPD)
        call cline%set('mskdiam',  DEFAULT_MSKDIAM)
    else
        call cline%parse_oldschool
    endif
    call cline%set('amsklp',  20.)
    call cline%set('automsk', 'no')
    call cline%set('part',    1.)
    call cline%checkvar('stk',     1)
    call cline%checkvar('smpd',    2)
    call cline%checkvar('mskdiam', 3)
    call cline%check
    call params%new(cline)
    ! read images
    call find_ldim_nptcls(params%stk, ldim, n)
    allocate(imgs(n))
    do i = 1, n
        call imgs(i)%new(ldim, params%smpd)
        call imgs(i)%read(params%stk, i)
    end do
    ! mask
    call automask2D(params, imgs, NGROW, WINSZ, EDGE, diams, shifts)
    call simple_end('**** SIMPLE_TEST_NANO_MASK_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_nano_mask

end module simple_commanders_test_masks

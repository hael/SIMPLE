!@descr: for all masks tests
module simple_commanders_test_masks
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

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

!@descr: input/output tests run by hand on user data (mrc2jpeg, mrc_validate)
! The hermetic I/O tests are the stack I/O (unit_core), binoris and STAR sub-suites (unit_project).
module simple_commanders_test_io
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_mrc2jpeg
  contains
    procedure :: execute      => exec_test_mrc2jpeg
end type commander_test_mrc2jpeg

type, extends(commander_base) :: commander_test_mrc_validate
  contains
    procedure :: execute      => exec_test_mrc_validate
end type commander_test_mrc_validate

contains

subroutine exec_test_mrc2jpeg( self, cline )
    use simple_jpg
    use simple_image,      only: image
    use simple_cmdline,    only: cmdline
    use simple_parameters, only: parameters
    class(commander_test_mrc2jpeg),    intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    type(string), allocatable :: micname(:)
    type(image)               :: microg 
    type(string)              :: outputfile, fbody
    integer                   :: i, j, nfiles, ldim(3), ifoo, ldim_refs(3)
    type(parameters)          :: p
    call cline%parse_oldschool
    call cline%checkvar('filetab', 1)
    call cline%checkvar('smpd',    2)
    call cline%check
    call p%new(cline)
    call read_filetable(p%filetab, micname)
    nfiles=size(micname)
    do i=1,nfiles
        fbody = get_fbody(basename(micname(i)),'mrc')
        call find_ldim_nptcls(micname(i), ldim, ifoo)
        ldim_refs = [ldim(1), ldim(2), 1]
        if(ldim(3)==1)then
            call microg%new([ldim(1), ldim(2), 1], p%smpd)
            call microg%read(micname(i))
            outputfile = fbody//'.jpeg'
            write(logfhandle,'(a)') '>>> WRITING '//outputfile%to_char()
            call microg%write_jpg(outputfile)
            call microg%kill()
        else
            call microg%new([ldim(1), ldim(2), 1], p%smpd)
            do j=1,ldim(3)
                call microg%read(micname(i),j)
                outputfile = fbody//int2str_pad(j,3)//'.jpeg'
                write(logfhandle,'(a)') '>>> WRITING '//outputfile%to_char()
                call microg%write_jpg(outputfile)
            enddo
            call microg%kill()
        endif
    enddo
    call simple_end('**** SIMPLE_TEST_MRC2JPEG_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mrc2jpeg

subroutine exec_test_mrc_validate( self, cline )
    use simple_image, only: image
    class(commander_test_mrc_validate), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(image)        :: vol
    type(string)       :: vol_file
    real               :: smpd
    integer            :: ldim(3), ifoo
    if( .not. cline%defined('vol')  ) THROW_HARD('The vol keyword is required; expected vol=volume.mrc')
    if( .not. cline%defined('smpd') ) THROW_HARD('The smpd keyword is required; expected smpd=<Angstrom per voxel>')
    vol_file = cline%get_carg('vol')
    smpd     = cline%get_rarg('smpd')
    call find_ldim_nptcls(vol_file, ldim, ifoo)
    write(logfhandle,'(a,3(1x,i0))') 'Input volume dimensions:', ldim
    call vol%new(ldim, smpd)
    call vol%read(vol_file)
    call vol%write(string('vol_simple.mrc'))
    call vol%kill
    call simple_end('**** SIMPLE_TEST_MRC_VALIDATE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mrc_validate

end module simple_commanders_test_io

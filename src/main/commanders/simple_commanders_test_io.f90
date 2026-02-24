!@descr: for all input/output tests
module simple_commanders_test_io
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_imgfile
  contains
    procedure :: execute      => exec_test_imgfile
end type commander_test_imgfile

type, extends(commander_base) :: commander_test_inside_write
  contains
    procedure :: execute      => exec_test_inside_write
end type commander_test_inside_write

type, extends(commander_base) :: commander_test_io
  contains
    procedure :: execute      => exec_test_io
end type commander_test_io

type, extends(commander_base) :: commander_test_io_parallel
  contains
    procedure :: execute      => exec_test_io_parallel
end type commander_test_io_parallel

type, extends(commander_base) :: commander_test_mrc2jpeg
  contains
    procedure :: execute      => exec_test_mrc2jpeg
end type commander_test_mrc2jpeg

type, extends(commander_base) :: commander_test_mrc_validation
  contains
    procedure :: execute      => exec_test_mrc_validation
end type commander_test_mrc_validation

type, extends(commander_base) :: commander_test_stack_io
  contains
    procedure :: execute      => exec_test_stack_io
end type commander_test_stack_io

type, extends(commander_base) :: commander_test_star_export
  contains
    procedure :: execute      => exec_test_star_export
end type commander_test_star_export

type, extends(commander_base) :: commander_test_starfile_test
  contains
    procedure :: execute      => exec_test_starfile_test
end type commander_test_starfile_test

contains

subroutine exec_test_imgfile( self, cline )
    use simple_image,   only: image
    use simple_imgfile, only: imgfile
    class(commander_test_imgfile),    intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    integer     :: ldim(3), i, j
    real        :: smpd=2., corr
    type(image) :: img
    type(image) :: imgs(20)
    logical     :: ft=.false.
    ! SELF-CONSISTENCY TESTS
    ! create a square
    ldim = [120,120,1]
    call img%new(ldim, smpd)
    call img%square(20)
    call img%write(string('squares_mrc.mrc'),1)
    ! write stacks of 5 squares
    do i=1,5
        if( ft ) call img%fft()
        call img%write(string('squares_spider.spi'),i)
        call img%write(string('squares_mrc.mrc'),i)
    end do
    ! create a cube
    ldim = [120,120,120]
    call img%new(ldim, smpd)
    call img%square(20)
    ! write volume files
    do i=1,5
        if( ft ) call img%fft()
        call img%write(string('cube_spider.spi'))
        call img%write(string('cube_mrc.mrc'))
    end do
    ! convert the cubes from SPIDER to MRC & vice versa
    do i=1,5
        call img%read(string('cube_spider.spi'))
        if( ft ) call img%ifft()
        call img%write(string('cube_spider_converted.mrc'))
        call img%read(string('cube_mrc.mrc'))
        if( ft ) call img%ifft()
        call img%write(string('cube_mrc_converted.spi'))
    end do
    ! test SPIDER vs. MRC & converted vs. nonconverted
    do i=1,4
        call imgs(i)%new(ldim, smpd)
        call imgs(i)%read(string('cube_spider.spi'))
        call imgs(i)%read(string('cube_spider_converted.mrc'))
        call imgs(i)%read(string('cube_mrc.mrc'))
        call imgs(i)%read(string('cube_mrc_converted.spi'))
        if( ft ) call imgs(i)%ifft()
    end do
    do i=1,3
        do j=i+1,4
            corr = imgs(i)%corr(imgs(j))
            if( corr < 0.99999 )then
                THROW_HARD('SPIDER vs. MRC & converted vs. nonconverted test failed')
            endif
        end do
    end do
    call simple_end('**** SIMPLE_TEST_IMGFILE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_imgfile

subroutine exec_test_inside_write( self, cline )
    class(commander_test_inside_write), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(sp_project) :: spproj
    integer          :: i
    call spproj%os_mic%new(10, is_ptcl=.false.)
    do i=1,10
        call spproj%os_mic%set(i, 'intg', 'intg')
    end do
    call spproj%os_ptcl2D%new(100, is_ptcl=.true.)
    call spproj%os_ptcl2D%rnd_oris
    call spproj%write(string('original_proj.simple'))
    call spproj%write(string('updated_proj.simple'))
    call spproj%os_stk%new(10, is_ptcl=.false.)
    do i=1,10
        call spproj%os_stk%set(i, 'stk', 'stk')
    end do
    call spproj%write_segment_inside('stk', string('updated_proj.simple'))
    call simple_end('**** SIMPLE_TEST_INSIDE_WRITE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_inside_write

subroutine exec_test_io( self, cline )
    use simple_image,    only: image
    use simple_stack_io, only: stack_io
    class(commander_test_io),    intent(inout) :: self
    class(cmdline),              intent(inout) :: cline
    type(image),      allocatable :: imgs(:)
    type(string)                  :: benchfname
    type(image)                   :: vol
    type(stack_io)                :: stkio_r, stkio_w
    integer,          parameter   :: NVOLRWS   = 10
    integer,          parameter   :: NSTKRWS   = 10
    integer,          parameter   :: NPTCLS    = 1024 * 10
    integer,          parameter   :: BOX       = 512
    integer,          parameter   :: ONE_M     = 1024**2
    integer(kind=8),  parameter   :: NSTKBYTES = NSTKRWS * NPTCLS * BOX * BOX * 4
    integer(kind=8),  parameter   :: NVOLBYTES = NVOLRWS * BOX    * BOX * BOX * 4
    real,             parameter   :: SMPD      = 1.0
    integer(timer_int_kind)       ::  t_stk_w,  t_stk_r,  t_vol_w,  t_vol_r
    real(timer_int_kind)          :: rt_stk_w, rt_stk_r, rt_vol_w, rt_vol_r, rt_tot
    integer :: iptcl, i, fnr
    real    :: mb_per_s_stk_w, mb_per_s_stk_r, mb_per_s_vol_w, mb_per_s_vol_r, mb_per_s_w, mb_per_s_r

    print *, 'simulating a stack of '//int2str(NPTCLS)//' particles'
    allocate(imgs(NPTCLS))
    do iptcl = 1, NPTCLS
        call imgs(iptcl)%new([BOX,BOX,1], SMPD)
        call imgs(iptcl)%ran()
    end do

    print *, 'simulating a volume'
    call vol%new([BOX,BOX,BOX], SMPD)
    call vol%ran()
    call vol%write(string('random_vol.mrc'))

    print *, 'writing the stack '//int2str(NSTKRWS)//' times'
    rt_tot  = 0.
    t_stk_w = tic()
    do i = 1, NSTKRWS
        call stkio_w%open(string('stack_of_random_imgs.mrcs'), SMPD, 'write', box=BOX)
        do iptcl = 1, NPTCLS
            call stkio_w%write(iptcl, imgs(iptcl))
        end do
        call stkio_w%close
    end do
    rt_stk_w = toc(t_stk_w)
    rt_tot   = rt_tot + rt_stk_w

    print *, 'reading the stack '//int2str(NSTKRWS)//' times'
    t_stk_r = tic()
    do i = 1, NSTKRWS
        call stkio_r%open(string('stack_of_random_imgs.mrcs'), SMPD, 'read')
        do iptcl = 1, NPTCLS
            call stkio_r%read(iptcl, imgs(iptcl))
        end do
        call stkio_r%close
    end do
    rt_stk_r = toc(t_stk_r)
    rt_tot   = rt_tot + rt_stk_r

    print *, 'writing the volume '//int2str(NVOLRWS)//' times'
    t_vol_w = tic()
    do i = 1, NVOLRWS
        call vol%write(string('random_vol.mrc'))
    end do
    rt_vol_w = toc(t_vol_w)
    rt_tot   = rt_tot + rt_vol_w

    print *, 'reading the volume '//int2str(NVOLRWS)//' times'
    t_vol_r = tic()
    do i = 1, NVOLRWS
        call vol%read(string('random_vol.mrc'))
    end do
    rt_vol_r = toc(t_vol_r)
    rt_tot   = rt_tot + rt_vol_r

    ! calc MB / s
    mb_per_s_stk_w = real(real(NSTKBYTES,dp)             / real(ONE_M,dp) / rt_stk_w)
    mb_per_s_stk_r = real(real(NSTKBYTES,dp)             / real(ONE_M,dp) / rt_stk_r)
    mb_per_s_vol_w = real(real(NVOLBYTES,dp)             / real(ONE_M,dp) / rt_vol_w)
    mb_per_s_vol_r = real(real(NVOLBYTES,dp)             / real(ONE_M,dp) / rt_vol_r)
    mb_per_s_w     = real(real(NSTKBYTES + NVOLBYTES,dp) / real(ONE_M,dp) / (rt_stk_w + rt_vol_w))
    mb_per_s_r     = real(real(NSTKBYTES + NVOLBYTES,dp) / real(ONE_M,dp) / (rt_stk_r + rt_vol_r))

    ! write benchmark stats
    benchfname = 'SIMPLE_IO_BENCH.txt'
    call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
    write(fnr,'(a)') '*** TIMINGS (s) ***'
    write(fnr,'(a,1x,f9.2)') 'stack  writes : ', rt_stk_w
    write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', rt_stk_r
    write(fnr,'(a,1x,f9.2)') 'volume writes : ', rt_vol_w
    write(fnr,'(a,1x,f9.2)') 'volume reads  : ', rt_vol_r
    write(fnr,'(a,1x,f9.2)') 'total  time   : ', rt_tot
    write(fnr,'(a)') ''
    write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
    write(fnr,'(a,1x,f9.2)') 'stack  writes : ', (rt_stk_w / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', (rt_stk_r / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') 'volume writes : ', (rt_vol_w / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') 'volume reads  : ', (rt_vol_r / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') '% accounted for : ', ((rt_stk_w+rt_stk_r+rt_vol_w+rt_vol_r)/rt_tot) * 100.
    write(fnr,'(a)') ''
    write(fnr,'(a)') '*** READ/WRITE SPEEDS (MB/s) ***'
    write(fnr,'(a,1x,f9.2)') 'stack  writes : ', mb_per_s_stk_w
    write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', mb_per_s_stk_r
    write(fnr,'(a,1x,f9.2)') 'volume writes : ', mb_per_s_vol_w
    write(fnr,'(a,1x,f9.2)') 'volume reads  : ', mb_per_s_vol_r
    write(fnr,'(a,1x,f9.2)') 'total  writes : ', mb_per_s_w
    write(fnr,'(a,1x,f9.2)') 'total  reads  : ', mb_per_s_r
    call fclose(fnr)
    call simple_end('**** SIMPLE_TEST_IO_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_io

subroutine exec_test_io_parallel( self, cline )
    use simple_image,    only: image
    use simple_stack_io, only: stack_io
    class(commander_test_io_parallel),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(image),      allocatable :: imgs(:)
    type(string)                  :: benchfname
    type(image)                   :: vol
    type(stack_io)                :: stkio_r, stkio_w
    integer,          parameter   :: NVOLRWS   = 10
    integer,          parameter   :: NSTKRWS   = 10
    integer,          parameter   :: NPTCLS    = 1024 * 10
    integer,          parameter   :: BOX       = 512
    integer,          parameter   :: ONE_M     = 1024**2
    integer(kind=8),  parameter   :: NSTKBYTES = NSTKRWS * NPTCLS * BOX * BOX * 4
    integer(kind=8),  parameter   :: NVOLBYTES = NVOLRWS * BOX    * BOX * BOX * 4
    real,             parameter   :: SMPD      = 1.0
    integer(timer_int_kind)       ::  t_stk_w,  t_stk_r,  t_vol_w,  t_vol_r
    real(timer_int_kind)          :: rt_stk_w, rt_stk_r, rt_vol_w, rt_vol_r, rt_tot
    integer :: iptcl, i, fnr
    real    :: mb_per_s_stk_w, mb_per_s_stk_r, mb_per_s_vol_w, mb_per_s_vol_r, mb_per_s_w, mb_per_s_r
    print *, 'simulating a stack of '//int2str(NPTCLS)//' particles'
    allocate(imgs(NPTCLS))
    do iptcl = 1, NPTCLS
        call imgs(iptcl)%new([BOX,BOX,1], SMPD)
        call imgs(iptcl)%ran()
    end do
    print *, 'simulating a volume'
    call vol%new([BOX,BOX,BOX], SMPD)
    call vol%ran()
    call vol%write(string('random_vol.mrc'))
    print *, 'writing the stack '//int2str(NSTKRWS)//' times'
    rt_tot  = 0.
    t_stk_w = tic()
    do i = 1, NSTKRWS
        call stkio_w%open(string('stack_of_random_imgs.mrcs'), SMPD, 'write', box=BOX)
        do iptcl = 1, NPTCLS
            call stkio_w%write(iptcl, imgs(iptcl))
        end do
        call stkio_w%close
    end do
    rt_stk_w = toc(t_stk_w)
    rt_tot   = rt_tot + rt_stk_w
    print *, 'reading the stack '//int2str(NSTKRWS)//' times'
    t_stk_r = tic()
    do i = 1, NSTKRWS
        call stkio_r%open(string('stack_of_random_imgs.mrcs'), SMPD, 'read')
        do iptcl = 1, NPTCLS
            call stkio_r%read(iptcl, imgs(iptcl))
        end do
        call stkio_r%close
    end do
    rt_stk_r = toc(t_stk_r)
    rt_tot   = rt_tot + rt_stk_r
    print *, 'writing the volume '//int2str(NVOLRWS)//' times'
    t_vol_w = tic()
    do i = 1, NVOLRWS
        call vol%write(string('random_vol.mrc'))
    end do
    rt_vol_w = toc(t_vol_w)
    rt_tot   = rt_tot + rt_vol_w
    print *, 'reading the volume '//int2str(NVOLRWS)//' times'
    t_vol_r = tic()
    do i = 1, NVOLRWS
        call vol%read(string('random_vol.mrc'))
    end do
    rt_vol_r = toc(t_vol_r)
    rt_tot   = rt_tot + rt_vol_r
    ! calc MB / s
    mb_per_s_stk_w = real(real(NSTKBYTES,dp)             / real(ONE_M,dp) / rt_stk_w)
    mb_per_s_stk_r = real(real(NSTKBYTES,dp)             / real(ONE_M,dp) / rt_stk_r)
    mb_per_s_vol_w = real(real(NVOLBYTES,dp)             / real(ONE_M,dp) / rt_vol_w)
    mb_per_s_vol_r = real(real(NVOLBYTES,dp)             / real(ONE_M,dp) / rt_vol_r)
    mb_per_s_w     = real(real(NSTKBYTES + NVOLBYTES,dp) / real(ONE_M,dp) / (rt_stk_w + rt_vol_w))
    mb_per_s_r     = real(real(NSTKBYTES + NVOLBYTES,dp) / real(ONE_M,dp) / (rt_stk_r + rt_vol_r))
    ! write benchmark stats
    benchfname = 'SIMPLE_IO_BENCH.txt'
    call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
    write(fnr,'(a)') '*** TIMINGS (s) ***'
    write(fnr,'(a,1x,f9.2)') 'stack  writes : ', rt_stk_w
    write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', rt_stk_r
    write(fnr,'(a,1x,f9.2)') 'volume writes : ', rt_vol_w
    write(fnr,'(a,1x,f9.2)') 'volume reads  : ', rt_vol_r
    write(fnr,'(a,1x,f9.2)') 'total  time   : ', rt_tot
    write(fnr,'(a)') ''
    write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
    write(fnr,'(a,1x,f9.2)') 'stack  writes : ', (rt_stk_w / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', (rt_stk_r / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') 'volume writes : ', (rt_vol_w / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') 'volume reads  : ', (rt_vol_r / rt_tot) * 100.
    write(fnr,'(a,1x,f9.2)') '% accounted for : ', ((rt_stk_w+rt_stk_r+rt_vol_w+rt_vol_r)/rt_tot) * 100.
    write(fnr,'(a)') ''
    write(fnr,'(a)') '*** READ/WRITE SPEEDS (MB/s) ***'
    write(fnr,'(a,1x,f9.2)') 'stack  writes : ', mb_per_s_stk_w
    write(fnr,'(a,1x,f9.2)') 'stack  reads  : ', mb_per_s_stk_r
    write(fnr,'(a,1x,f9.2)') 'volume writes : ', mb_per_s_vol_w
    write(fnr,'(a,1x,f9.2)') 'volume reads  : ', mb_per_s_vol_r
    write(fnr,'(a,1x,f9.2)') 'total  writes : ', mb_per_s_w
    write(fnr,'(a,1x,f9.2)') 'total  reads  : ', mb_per_s_r
    call fclose(fnr)
    call simple_end('**** SIMPLE_TEST_IO_PARALLEL_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_io_parallel

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

subroutine exec_test_mrc_validation( self, cline )
    use simple_atoms, only: atoms
    use simple_image, only: image
    class(commander_test_mrc_validation), intent(inout) :: self
    class(cmdline),                       intent(inout) :: cline
    character(len=STDLEN)         :: vol_file
    character(len=:), allocatable :: smpd_char
    type(image)                   :: vol 
    real                          :: smpd
    integer                       :: ldim(3), ifoo, slen
    if( command_argument_count() /= 2 )then
        write(logfhandle,'(a)') 'ERROR! Usage: simple_test_mrc_validate vol.mrc smpd'
        write(logfhandle,'(a)') 'vol.mrc : volume' 
        write(logfhandle,'(a)') 'smpd    : SMPD value in Angstrom per voxel ' 
    else
        call get_command_argument(1, vol_file)
        call get_command_argument(2, length=slen)
        allocate(character(slen) :: smpd_char)
        call get_command_argument(2, smpd_char)
        read(smpd_char, *) smpd
    endif
    call find_ldim_nptcls(string(trim(vol_file)),  ldim, ifoo)
    print *, trim(vol_file), ldim
    call vol%new(ldim, smpd)
    call vol%read(string(trim(vol_file)))
    call vol%write(string('vol_simple.mrc'))
    call vol%kill
    call simple_end('**** SIMPLE_TEST_MRC_VALIDATION_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mrc_validation

subroutine exec_test_stack_io( self, cline )
    use simple_stack_io, only: stack_io
    use simple_image,    only: image
    class(commander_test_stack_io),    intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    type(stack_io)              :: stkio_r, stkio_w
    type(image)                 :: img
    character(len=*), parameter :: stkname = 'cavgs_iter030_ranked.mrc'
    real,             parameter :: smpd    = 1.3
    integer                     :: nptcls, iptcl, ldim(3)
    call stkio_r%open(string(stkname), smpd, 'read', bufsz=100)
    call stkio_w%open(string('outstk_written.mrc'), smpd, 'write', box=256, is_ft=.false.)
    nptcls = stkio_r%get_nptcls()
    ldim   = stkio_r%get_ldim()
    call img%new(ldim, smpd)
    ! read
    do iptcl = 1, nptcls
        call stkio_r%read(iptcl, img)
        call img%write(string('outstk_read.mrc'), iptcl)
    end do
    call stkio_r%close
    ! write
    do iptcl = 1, nptcls
        call img%read(string(stkname), iptcl)
        call stkio_w%write(iptcl, img)
    end do
    call stkio_w%close
    ! readwrite
    call stkio_r%open(string(stkname), smpd, 'read', bufsz=100)
    call stkio_w%open(string('outstk_read_written.mrc'), smpd, 'write', box=256, is_ft=.false.)
    nptcls = stkio_r%get_nptcls()
    ldim   = stkio_r%get_ldim()
    call img%new(ldim, smpd)
    do iptcl = 1, nptcls
        call stkio_r%read(iptcl, img)
        call stkio_w%write(iptcl, img)
    end do
    call stkio_r%close
    call stkio_w%close
    call simple_end('**** SIMPLE_TEST_STACK_IO_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_stack_io

subroutine exec_test_star_export( self, cline )
    use simple_sp_project
    use simple_starfile
    use simple_parameters
    class(commander_test_star_export),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    character(len=*), parameter :: projfile   = 'test.simple'
    character(len=*), parameter :: opticsfile = 'optics.simple'
    type(sp_project)            :: spproj
    integer(timer_int_kind)     :: ms0
    real(timer_int_kind)        :: ms_complete
    if(.not. file_exists(string(projfile))) THROW_HARD(projfile // " does not exist")
    ms0 = tic()
    call spproj%read(string(projfile))
    ms_complete = toc(ms0)
    print *,'read project file in : ', ms_complete; call flush(6)
    ms0 = tic()
    call spproj%write_mics_star()
    ms_complete = toc(ms0)
    print *,'write_mics_star file in : ', ms_complete; call flush(6)
    ms0 = tic()
    call spproj%write_ptcl2D_star()
    ms_complete = toc(ms0)
    print *,'write_ptcl2D_table file in : ', ms_complete; call flush(6)
    call simple_end('**** SIMPLE_TEST_STAR_EXPORT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_star_export

subroutine exec_test_starfile_test( self, cline )
    use simple_starproject_tester
    use simple_starfile_wrappers
    class(commander_test_starfile_test), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    type(starfile_table_type)     :: sfile
    logical                       :: aresult
    character(len=:), allocatable :: retrieved_string
    real(kind=8)                  :: amt, aml
    integer(C_long)               :: object_id, num_objects
    ! step 1: write star-file
    ! alloc and open output file
    call starfile_table__new(sfile)
    call starfile_table__open_ofile(sfile, "outputfile.star")
    call starfile_table__clear(sfile)
    ! first segment ("name1")
    call starfile_table__setName(sfile, "name1")
    ! first segment is list
    call starfile_table__setislist(sfile, .true.)
    call starfile_table__addObject(sfile)
    ! add a comment for good measure
    call starfile_table__setComment(sfile, "this_is_a_comment")
    ! add 3 fields, 1 string and 2 doubles
    call starfile_table__setValue_string(sfile, EMDL_MICROGRAPH_NAME, "this_is_a_string")
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL, 0.12345678901234567890d0)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_EARLY, 42._8)
    ! write first segment
    call starfile_table__write_ofile(sfile)
    call starfile_table__clear(sfile)
    ! create second segment ("name2")
    call starfile_table__setName(sfile, "name2")
    call starfile_table__addObject(sfile)
    ! this one is not a list
    call starfile_table__setIsList(sfile, .false.)
    ! 4 double values
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 0.12345678901234567890d0)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 456._8)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 789._8)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 101112345678._8)
    ! write second segment
    call starfile_table__write_ofile(sfile)
    call starfile_table__clear(sfile)
    ! create third segment ("name3")
    call starfile_table__setName(sfile, "name3")
    ! this one is a list again
    call starfile_table__setIsList(sfile, .true.)
    call starfile_table__addObject(sfile)
    call starfile_table__setValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, 523._8)
    ! write and deallocate
    call starfile_table__write_ofile(sfile)
    call starfile_table__close_ofile(sfile)
    call starfile_table__delete(sfile)
    ! reallocate and read first segment
    call starfile_table__new(sfile)
    call starfile_table__read(sfile, string("outputfile.star"), "name1")
    ! see if we can retrieve string correctly
    aresult = starfile_table__getValue_string(sfile, EMDL_MICROGRAPH_NAME, retrieved_string)
    write (*,*) 'result=', aresult, ' retrieved_string=', retrieved_string
    ! this one should fail
    aresult = starfile_table__getValue_string(sfile, EMDL_MLMODEL_REF_IMAGE, retrieved_string)
    if (aresult) then
       write (*,*) 'result=', aresult, ' retrieved_string=', retrieved_string
    else
       write (*,*) 'result=', aresult
    end if
    ! now read the other two fields; the first one should go through, the other should fail
    aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL, amt)
    write (*,*) 'aresult = ', aresult, ' ; amt = ', amt
    aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, aml)
    write (*,*) 'aresult = ', aresult, ' ; aml = ', aml
    write (*,*) '-----------------'
    ! now read in second segment
    call starfile_table__read(sfile, string("outputfile.star"), "name2")
    ! iterate through list
    object_id = starfile_table__firstobject(sfile)
    num_objects = starfile_table__numberofobjects(sfile)
    do while( (object_id < num_objects) .and. (object_id >= 0) )
          aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, aml)
          write (*,*) 'aresult = ', aresult, ' ; aml = ', aml
       object_id = starfile_table__nextobject(sfile)
    end do
    ! now read in third segment
    call starfile_table__read(sfile, string("outputfile.star"), "name3")
    write (*,*) '-----------------'
    aresult = starfile_table__getValue_double(sfile, EMDL_MICROGRAPH_ACCUM_MOTION_LATE, aml)
    write (*,*) 'aresult = ', aresult, ' ; aml = ', aml
    ! deallocate
    call starfile_table__delete(sfile)
    ! chatgpt generated tests
    call run_all_starproject_tests
    call simple_end('**** SIMPLE_TEST_STARFILE_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_starfile_test

end module simple_commanders_test_io

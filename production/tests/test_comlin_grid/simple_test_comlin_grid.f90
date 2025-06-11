program simple_test_comlin_grid
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_projector,  only: projector
use simple_oris
use simple_ori
implicit none
integer,          parameter   :: ORI_IND = 15, NPLANES = 100
character(len=:), allocatable :: cmd
real,             allocatable :: pspec(:)
type(parameters) :: p
type(cmdline)    :: cline
type(image)      :: vol, noise, ptcl, ptcl_pad, roavg
type(projector)  :: vol_pad
type(oris)       :: spiral
type(ori)        :: o1
integer          :: rc, ifoo, iind
real             :: ave, sdev, maxv, minv, masscen(3), sh(2)
logical          :: mrc_exists
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_comlin_grid smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
    write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Converting .pdb to .mrc...'
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('vol1'   , '1JYX.mrc')
    call cline%set('mskdiam', 180.)
    call cline%set('lp'   ,   3.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('vol1',    3)
call cline%checkvar('mskdiam', 4)
call cline%checkvar('lp',      5)
call cline%check
call p%new(cline)
print *, 'box = ', p%box
print *, 'msk = ', p%msk
call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
call vol%new(p%ldim, p%smpd)
call vol%read(p%vols(1))
call vol%stats('foreground', ave, sdev, maxv, minv)
call spiral%new(NPLANES, is_ptcl=.false.)
call spiral%spiral
call ptcl%new(    [p%box,   p%box,   1],       p%smpd)
call noise%new(   [p%box,   p%box ,  1],       p%smpd)
call vol_pad%new( [p%boxpd, p%boxpd, p%boxpd], p%smpd)
call ptcl_pad%new([p%boxpd, p%boxpd, 1],       p%smpd)
call vol%pad(vol_pad)
call vol_pad%fft
call vol_pad%expand_cmat(p%alpha)
call spiral%get_ori(ORI_IND, o1)
call vol_pad%fproject(o1,ptcl_pad)
call ptcl_pad%ifft
! add noise in a small center region of the vol
call noise%gauran(0., 0.01 * sdev)
call ptcl%zero_and_unflag_ft
call ptcl_pad%clip(ptcl)
call ptcl%add(noise)

end program simple_test_comlin_grid

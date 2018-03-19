!==Program simple_select
!
! <select/begin> is a program for selecting files in table based on image correlation matching. <select/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_select
include 'simple_lib.f08'
use simple_params,  only: params
use simple_build,   only: build
use simple_imgfile, only: imgfile
use simple_image,   only: image
use simple_cmdline, only: cmdline
use simple_corrmat  ! singleton
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
type(image), allocatable           :: imgs_sel(:), imgs_all(:)
type(image)                        :: stk3_img
character(len=STDLEN), allocatable :: imgnames(:)
integer                            :: iimg, isel, nall, nsel, loc(1), ios, funit, ldim(3), ifoo, lfoo(3)
integer, allocatable               :: selected(:)
real, allocatable                  :: correlations(:,:)
logical                            :: debug=.false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_SELECT stk=<all_imgs.ext> stk2=<selected_imgs.ext>'
    write(*,'(a)',advance='no') ' [stk3<stk2selectfrom.ext>] [filetab=<table2selectfrom.txt>]'
    write(*,'(a)',advance='no') ' [outfile=<selected_lines.txt>] [outstk=outstk.ext] [dir=<move files 2 here{selected}>]'
    write(*,'(a)') ' [nthr=<nr of OpenMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',     1)
call cline%checkvar('stk2',    2)
if( cline%defined('stk3') .or. cline%defined('filetab')  )then
    ! all good
else
    stop 'Need either stk3 or filetab on the command line!'
endif
if( .not. cline%defined('outfile') )then
    call cline%set('outfile', 'selected_lines.txt')
endif
call cline%check
p = params(cline)            ! parameters generated
call b%build_general_tbox(p, cline) ! general objects built
! find number of selected cavgs
call find_ldim_nptcls(p%stk2, lfoo, nsel)
! find number of original cavgs
call find_ldim_nptcls(p%stk, lfoo, nall)
! read images
allocate(imgs_sel(nsel), imgs_all(nall))
do isel=1,nsel
    call imgs_sel(isel)%new([p%box,p%box,1], p%smpd)
    call imgs_sel(isel)%read(p%stk2, isel)
end do
do iimg=1,nall
    call imgs_all(iimg)%new([p%box,p%box,1], p%smpd)
    call imgs_all(iimg)%read(p%stk, iimg)
end do
write(*,'(a)') '>>> CALCULATING CORRELATIONS'
call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
! find selected
allocate(selected(nsel))
do isel=1,nsel
    loc = maxloc(correlations(isel,:))
    selected(isel) = loc(1)
    if( debug ) print *, 'selected: ', loc(1), ' with corr: ', correlations(isel,loc(1))
end do
if( cline%defined('filetab') )then
    ! read filetable
    call read_filetable(p%filetab, imgnames)
    if( size(imgnames) /= nall ) stop 'nr of entries in filetab and stk not consistent'

    call fopen(funit,p%outfile, status="replace", action="write", iostat=ios, access="sequential")
    if( ios /= 0 )then
        write(*,*) "Error opening file name", trim(adjustl(p%outfile))
        stop
    endif
    call simple_mkdir(trim(adjustl(p%dir)))
    ! write outoput & move files
    do isel=1,nsel
        write(funit,'(a)') trim(adjustl(imgnames(selected(isel))))
        call simple_rename(trim(adjustl(imgnames(selected(isel)))), trim(adjustl(p%dir)))
    end do
    close(funit)
endif
if( cline%defined('stk3') )then
    call find_ldim_nptcls(p%stk3, ldim, ifoo)
    ldim(3) = 1
    call stk3_img%new(ldim,1.)
    do isel=1,nsel
        call stk3_img%read(p%stk3, selected(isel))
        call stk3_img%write(p%outstk, isel)
    end do
endif
call simple_end('**** SIMPLE_SELECT NORMAL STOP ****')
end program simple_select

!==Program simple_join_ptcls
!
! <join_ptcls/begin> is a program for joining multiple image stacks produced during 
! parallel simple\_extr\_ptcls execution into one <join_ptcls/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_join_ptcls
use simple_cmdline  ! singleton
use simple_defs     ! singleton
use simple_jiffys   ! singleton
use simple_build,   only: build
use simple_params,  only: params
use simple_oris,    only: oris
implicit none
type(params)          :: p
type(build)           :: b
integer               :: frame, cnt, nimgs, i, j, prog_cnt, prog_n, numlen, lfoo(3)
character(len=STDLEN) :: framestack, sumstack, framestack_glob, sumstack_glob
type(oris)            :: o_merge
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_JOIN_PTCLS stk=<sumstack_part1.ext> nframes=<nr movie frames>'
    write(*,'(a)') ' npart=<nr partitions>'
    stop
endif
call parse_cmdline
call cmdcheckvar('stk',     1)
call cmdcheckvar('nframes', 2)
call cmdcheckvar('npart',   3)
call cmdcheck
p = params() ! parameters generated
call b%build_general_tbox(p,do3d=.false.)
write(*,'(a)') '>>> JOINING THE FRAMESTACKS'
prog_cnt = 0
prog_n   = p%nframes*p%npart
numlen   = len(int2str(p%nframes))
do frame=1,p%nframes
    framestack_glob = 'framestack'//int2str_pad(frame,numlen)//p%ext
    cnt = 0
    do i=1,p%npart
        prog_cnt = prog_cnt+1
        call progress(prog_cnt,prog_n)
        framestack = 'framestack'//int2str_pad(frame,numlen)//'_part'//int2str_pad(i,p%numlen)//p%ext
        ! get number of images from stack
        call find_ldim_nptcls(framestack, lfoo, nimgs)
        do j=1,nimgs
            cnt = cnt+1
            call b%img%read(trim(adjustl(framestack)), j)
            call b%img%write(trim(adjustl(framestack_glob)), cnt)
        end do
        call del_binfile(framestack)
    end do
end do
write(*,'(a)') '>>> JOINING THE SUMSTACKS & MERGING THE DOCS'
sumstack_glob = 'sumstack'//p%ext
cnt = 0
! delete possible previous merged parameter file
call del_txtfile('extr_ptcls_params_merged.txt')
do i=1,p%npart
    call progress(i,p%npart)
    sumstack = 'sumstack'//'_part'//int2str_pad(i,p%numlen)//p%ext
    ! get number of images from stack
    call find_ldim_nptcls(sumstack, lfoo, nimgs)
    do j=1,nimgs
        cnt = cnt+1
        call b%img%read(trim(adjustl(sumstack)), j)
        call b%img%write(trim(adjustl(sumstack_glob)), cnt)
    end do
    call del_binfile(sumstack)
    call o_merge%merge_files('extr_ptcls_params_merged.txt',&
    'extr_ptcls_params'//'_part'//int2str_pad(i,p%numlen)//'.txt')
    call del_txtfile('extr_ptcls_params'//'_part'//int2str_pad(i,p%numlen)//'.txt')
end do
call o_merge%write('extr_ptcls_params_merged.txt')
call simple_end('**** SIMPLE_JOIN_PTCLS NORMAL STOP ****')
end program simple_join_ptcls

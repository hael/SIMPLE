program simple_test_ori
include 'simple_lib.f08'
use simple_ori,  only: ori
use simple_oris, only: oris
use simple_binoris
use simple_binoris_io,   only: binwrite_oritab
use simple_sp_project

implicit none
type(sp_project) :: spproj, spproj2
type(oris)                :: os, os2
type(ori)                 :: o, o2
type(binoris)             :: bos
character(len=LONGSTRLEN) :: str
integer :: i
! call o%new(is_ptcl=.true.)
! str = o%pparms2str()
! write(*,'(a)') trim(str) ! works as intended
! str = o%ori2str()
! write(*,'(a)') trim(str) ! works as intended
! call o%set('mic', 'micrograph.mrc')
! call o%set('ctf', 'yes')
! str = o%ori2str()
! write(*,'(a)') trim(str)
! call o2%str2ori(str, is_ptcl=.true.)
! write(*,'(a)') '***************************'
! call o2%print_ori ! works as intended

! CLEVER compiler
! call spproj%os_ptcl2D%new(10, is_ptcl=.true.)
! do i=1,10
!     call spproj%os_ptcl2D%set(i, 'ctf', 'no')
!     call spproj%os_ptcl2D%set(i, 'corr', 0.6458)
!     call spproj%os_ptcl2D%set(i, 'x', 5.0)
! end do
! call binwrite_oritab('test.simple', spproj, spproj%os_ptcl2D, [1,10], isegment=3)

! call bos%open('/Processing/bgal/ori_refac2/2_cleanup2D/algndoc_2.simple')
! call bos%read_header
! do i = 1,13
!     print *,bos%header(i)%fromto,bos%header(i)%n_bytes_per_record
! enddo
! call bos%close



call spproj2%read_segment('ptcl2D', '/Processing/bgal/ori_refac2/2_cleanup2D/algndoc_1.simple')
call spproj2%os_ptcl2D%print_(1)
call spproj2%os_ptcl2D%print_(10)



end program simple_test_ori

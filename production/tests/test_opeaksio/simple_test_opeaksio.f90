program simple_test_opeaksio
! use simple_oris, only: oris
! use simple_binoris_io
! use simple_strings
! use simple_defs_fname
! implicit none
!
!
! type(oris) :: o_peaks(5)
! integer :: i
!
! do i=1,5
!     call o_peaks(i)%new(3)
!     call o_peaks(i)%set_all2single('iptcl', real(i))
!     call o_peaks(i)%set_euler(1, [real(i),real(i),real(i)])
!     call o_peaks(i)%set_euler(2, [real(i),real(i),real(i)])
!     call o_peaks(i)%set_euler(3, [real(i),real(i),real(i)])
!     call o_peaks(i)%set(1, 'ow', real(i))
!     call o_peaks(i)%set(2, 'ow', real(i)-0.2)
!     call o_peaks(i)%set(3, 'ow', real(i)-0.4)
!     call o_peaks(i)%write('peakset'//int2str(i)//'.txt')
! end do
!
! call binwrite_o_peaks('o_peaks_file.simple', [1,5], o_peaks, PTCL3D_SEG)
!
! ! do i=1,5
! !     call o_peaks(i)%kill
! ! end do
!
! ! call binread_o_peaks('o_peaks_file.simple', [1,5], o_peaks, PTCL3D_SEG)
!
! ! do i=1,5
! !     call o_peaks(i)%write('peakset_after'//int2str(i)//'.txt')
! ! end do
!
! ! wipe all odds, they should now be the same
! call o_peaks(1)%kill
! call o_peaks(3)%kill
! call o_peaks(5)%kill
! ! replace evens, they should be replaced in peakset_after files
! do i=2,4,2
!     call o_peaks(i)%set_euler(1, [6.,6.,6.])
!     call o_peaks(i)%set_euler(2, [6.,6.,6.])
!     call o_peaks(i)%set_euler(3, [6.,6.,6.])
!     call o_peaks(i)%set(1, 'ow', 6.)
!     call o_peaks(i)%set(2, 'ow', 6.-0.2)
!     call o_peaks(i)%set(3, 'ow', 6.-0.4)
! end do
!
! call binwrite_o_peaks('o_peaks_file.simple', [1,5], o_peaks, PTCL3D_SEG)
! call binread_o_peaks('o_peaks_file.simple', [1,5], o_peaks, PTCL3D_SEG)
!
! do i=1,5
!     call o_peaks(i)%write('peakset_after'//int2str(i)//'.txt')
! end do



end program simple_test_opeaksio

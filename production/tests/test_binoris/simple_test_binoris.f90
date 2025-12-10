program simple_test_binoris
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_commanders_project, only: commander_new_project
use simple_sp_project
implicit none
type(commander_new_project) :: xnew_project
type(oris)                  :: a1, os_peak1, os_peak2
type(cmdline)               :: cline
type(binoris)               :: bos
logical                     :: mask(5), test_passed
type(string)                :: fname, projname
#include "simple_local_flags.inc"
test_passed = .true.
fname       = 'spproject.simple'
projname    = 'spproject'
mask        = [.true.,.false.,.true.,.false.,.true.]
call a1%new(5, is_ptcl=.false.)
call a1%rnd_oris
call a1%write(string('oris1_facit.txt'))
call os_peak1%new(3, is_ptcl=.false.)
call os_peak1%rnd_oris()
call os_peak1%rnd_states(2)
call os_peak1%rnd_corrs()
call os_peak1%set_all2single('ow', 0.3)
call os_peak1%set_all2single('proj', 3.0)
call os_peak1%write(string('os_peak1_facit.txt'))
call os_peak2%new(5, is_ptcl=.false.)
call os_peak2%rnd_oris()
call os_peak2%rnd_states(2)
call os_peak2%rnd_corrs()
call os_peak2%set_all2single('ow', 0.5)
call os_peak2%set_all2single('proj', 5.0)
call os_peak2%write(string('os_peak2_facit.txt'))
call cline%set('projname', projname)
call cline%set('mkdir',        'no')
call cline%check()
call xnew_project%execute_safe(cline)
print *,'>>> OPEN BINORIS'
call bos%open(fname)
print *,'>>> WRITE HEADER '
call bos%write_header()
 print *,'>>> WRITE SEGMENT '
 call bos%write_segment(1, a1)
 print *,'>>> ADD SEGMENT '
 call bos%add_segment(1, a1)
 print *,'>>> UPDATE BYTE RANGES '
 call bos%update_byte_ranges()
! print *,'>>> READ FIRST SEGMENT RECORD '
! call bos%read_first_segment_record(PROJINFO_SEG, o)
! print *,'>>> READ SEGMENT '
! call bos%read_segment(PROJINFO_SEG, o)
! print *,'>>> READ RECORD '
! call bos%read_record()
! print *,'>>> GET SEGMENTS INFO '
! call bos%get_segments_info()
! print *,'>>> GET N SEGMENTS '
! call bos%get_n_segments()
! print *,'>>> GET FROM TO'
! call box%get_fromto()
! print *,'>>> GET N RECORDS '
! call bos%get_n_records()
! print *,'>>> GET N BYTES PER RECORD '
! call bos%get_n_bytes_per_record()
! print *,'>>> GET N BYTES TOTAL '
! call bos%get_n_bytes_tot()
! print *,'>>> GET FIRST DATA BYTE '
! call bos%get_first_data_byte()
print *,'>>> IS OPENED '
if(.not. bos%is_opened()) test_passed = .false.
print *,'>>> CLOSE BINORIS'
call bos%close()
if( test_passed )then
    print *, '>>> TEST PASSED'
else
    THROW_HARD('>>> TEST FAILED')
endif
! do i=1,5
!     call psrch3D(i)%set_o_peaks(os_peak1)
! end do
! call binwrite_oritab('oris_fill_in.bin', a1, [1,5], psrch3D)
! do i=1,5
!     if( mask(i) ) call psrch3D(i)%set_o_peaks(os_peak2)
! end do
! call binwrite_oritab('oris_mixed.bin', a1, [1,5], psrch3D, mask, 'oris_fill_in.bin')
! call a2%new(5)
! call binread_oritab('oris_mixed.bin', a2, [1,5], psrch3D)
! call a2%write('oris_mixed.txt')
! do i=1,5
!     os_peak1 = psrch3D(i)%get_o_peaks()
!     call os_peak1%write('os_peak'//int2str(i)//'.txt')
! end do
! call a%new(10)
! call a%rnd_oris
! call a%write('oris1_facit.txt')
! call binwrite_oritab('oris1_facit.bin', a, [1,10])
! call a2%new(10)
! call binread_oritab('oris1_facit.bin', a2, [1,10])
! call a2%write('oris1.txt')
! call a%new(10)
! call a%rnd_oris
! call a%rnd_ctf(300., 2.7, 0.1, 2.0, 0.5)
! call a%rnd_states(2)
! call a%write('oris2_facit.txt')
! call binwrite_oritab('oris2_facit.bin', a, [1,10])
! call a2%new(10)
! call binread_ctfparams_and_state('oris2_facit.bin', a2, [1,10])
! call a2%write('deftab.txt')
! call os_peak%new(3)
! call os_peak%rnd_oris()
! call os_peak%rnd_states(2)
! call os_peak%rnd_corrs()
! call os_peak%set_all2single('ow', 0.5)
! call os_peak%set_all2single('proj', 5.0)
! call os_peak%write('os_peak_facit.txt')
! call a%write('oris_full_facit.txt')
! do ipart=1,10
!     call psrch3D(ipart)%set_o_peaks(os_peak)
! end do
! call binwrite_oritab('oris_full.bin', a, [1,10], psrch3D )
! call a2%new(10)
! call binread_oritab('oris_full.bin', a2, [1,10], psrch3D)
! call os_peak%kill
! do ipart=1,10
!     os_peak = psrch3D(ipart)%get_o_peaks()
!     call os_peak%write('os_peak'//int2str(ipart)//'.txt')
! end do
! call a2%write('oris_full.txt')
! ! test the fill-in functionality
! mask = .false.
! mask(1:5) = .true.
! call binwrite_oritab('oris_full_fill_in.bin', a, [1,10], psrch3D, mask, 'oris_full.bin')
! call a2%new(10)
! call binread_oritab('oris_full_fill_in.bin', a2, [1,10], psrch3D)
! call a2%write('oris_full_fill_in.txt')
! parts = split_nobjs_even(10,5)
! do ipart=1,5
!     call bos%new(a, parts(ipart,:))
!     call bos%open('bos'//int2str(ipart)//'.bin', del_if_exists=.true.)
!     call bos%write_header
!     do iptcl=parts(ipart,1),parts(ipart,2)
!         call bos%write_record(iptcl, a)
!     end do
!     call bos%kill
! end do
! call a%kill
! call a%new(10)
! do ipart=1,5
!     call bos%open('bos'//int2str(ipart)//'.bin')
!     do iptcl=parts(ipart,1),parts(ipart,2)
!         call bos%read_record(iptcl, a)
!     end do
!     call bos%kill
! end do
! call a%write('oris1.txt')
! call a%new(10)
! call a%rnd_oris
! call a%write('oris2_facit.txt')
! call os_peak%new(3)
! call os_peak%rnd_oris()
! call os_peak%rnd_states(2)
! call os_peak%rnd_corrs()
! call os_peak%set_all2single('ow', 0.5)
! call os_peak%set_all2single('proj', 5.0)
! call os_peak%write('os_peak_facit.txt')
! do ipart=1,5
!     call bos%new(a, parts(ipart,:), os_peak)
!     call bos%open('bos'//int2str(ipart)//'.bin', del_if_exists=.true.)
!     call bos%write_header
!     do iptcl=parts(ipart,1),parts(ipart,2)
!         call bos%write_record(iptcl, a, os_peak)
!     end do
!     call bos%kill
! end do
! call a%kill
! call a%new(10)
! do ipart=1,5
!     call bos%open('bos'//int2str(ipart)//'.bin')
!     do iptcl=parts(ipart,1),parts(ipart,2)
!         call bos%read_record(iptcl, a, os_peak)
!         call os_peak%write('os_peak'//int2str(ipart)//'.txt')
!     end do
!     call bos%kill
! end do
! call a%write('oris2.txt')
end program simple_test_binoris

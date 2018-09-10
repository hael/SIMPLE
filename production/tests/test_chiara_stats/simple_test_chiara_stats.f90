program simple_test_chiara_stats
include 'simple_lib.f08'
use simple_oris
use gnufor2
use simple_picker_chiara
use simple_micops
use simple_image
use simple_stackops
use simple_math
use simple_test_chiara_try_mod
use simple_edge_detector, only : automatic_thresh_sobel
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline

type(image)        :: imgwin_bin, imgwin
integer            :: n_image
real               :: ave, sdev, maxv, minv, med !stats
real               :: aveB, sdevB, maxvB, minvB, medB !stats
logical            :: errout, erroutB
character (len = *), parameter :: fmt_1 = "(a)" !I/O
character(len = 100) :: iom
integer :: status
type(oris)        :: os
integer           :: nl, i
real, allocatable :: fore_avg(:), fore_avg_bin(:), back_avg(:), back_avg_bin(:)
real, allocatable :: fore_sdev(:), fore_sdev_bin(:), back_sdev(:), back_sdev_bin(:)
real, allocatable :: fore_minv(:), back_minv(:), fore_maxv(:), back_maxv(:)
real, allocatable :: indices(:)
call imgwin_bin%new([70,70,1], 1.)
call imgwin%new([70,70,1], 1.)
open(unit = 17, access = 'sequential', action = 'readwrite',file = "StatisticsWindowsFore.txt", form = 'formatted', iomsg = iom, iostat = status, position = 'append', status = 'replace')
open(unit = 18, file = "StatisticsWindowsBack.txt")
open(unit = 19, file = "StatisticsBinaryFore.txt")
open(unit = 20, file = "StatisticsBinaryBack.txt")
do n_image = 1, 74
    call imgwin_bin%read('/home/lenovoc30/Desktop/MassCenter/try1/centered_particles_BIN.mrc', n_image)
    call imgwin%read('/home/lenovoc30/Desktop/MassCenter/try1/centered_particles.mrc', n_image)
    !FOREGROUND
    call imgwin%stats('foreground', ave, sdev, maxv, minv)
    call imgwin_bin%stats('foreground', aveB, sdevB, maxvB, minvB)
    write(unit = 17, fmt = "(a,i0,4(a,f0.0))") &
        &'image=', n_image,' ave=', ave, ' sdev=', sdev,' maxv=', maxv, ' minv=', minv
    write(unit = 19, fmt = "(a,i0,4(a,f0.0))") &
        &'Bimage=', n_image, ' aveB=', aveB, ' sdevB=', sdevB, ' maxvB=', maxvB, ' minvB=', minvB
    !BACKGROUND STATS
    call imgwin%stats('background', ave, sdev, maxv, minv)
    call imgwin_bin%stats('background', aveB, sdevB, maxvB, minvB)
    write(unit = 18, fmt = "(a,i0,4(a,f0.0))") &
        &'image=', n_image, ' back-ave=', ave, ' back-sdev=', sdev,' back-maxv=', maxv, ' back-minv=', minv
    write(unit = 20, fmt = "(a,i0,4(a,f0.0))") &
        &'Bimage=', n_image,' back-aveB=', aveB, ' back-sdevB=', sdevB, ' back-maxvB=', maxvB, ' back-minvB=', minvB
enddo
close(17, status = "keep")
close(18, status = "keep")
close(19, status = "keep")
close(20, status = "keep")
nl = nlines('/home/lenovoc30/Desktop/MassCenter/try1/StatisticsWindowsFore.txt')
call os%new(nl)
call os%read('StatisticsWindowsFore.txt')
allocate(fore_avg(nl), fore_sdev(nl),fore_minv(nl), fore_maxv(nl))
allocate(indices(nl))
fore_avg     = os%get_all('ave')
! fore_avg_bin = os%get_all('aveB')
! back_avg     = os%get_all('back-ave')
! back_avg_bin = os%get_all('back-aveB')
!
! fore_sdev     = os%get_all('sdev')
! fore_sdev_bin = os%get_all('sdevB')
! back_sdev     = os%get_all('back-sdev')
! back_sdev_bin = os%get_all('back-sdevB')
!
! fore_minv     = os%get_all('minv')
! fore_maxv     = os%get_all('maxv')
! back_minv     = os%get_all('back-minv')
! back_maxv     = os%get_all('back-maxv')

allocate(fore_avg(nl), fore_avg_bin(nl), back_avg(nl), back_avg_bin(nl))
allocate(fore_sdev(nl), fore_sdev_bin(nl), back_sdev(nl), back_sdev_bin(nl))
allocate(fore_minv(nl), back_minv(nl), fore_maxv(nl), back_maxv(nl))
allocate(indices(nl))

do i = 1, nl
    indices(i) = real(i)
enddo
call plot(indices,fore_avg)
! call plot(indices,fore_avg_bin, color1 = 'black')
! call plot(indices,back_avg, color1 = 'red')
print *,  fore_avg
end program simple_test_chiara_stats

module simple_magic_boxes
use simple_math, only: find
implicit none

public :: find_magic_box, print_magic_box_range, find_magic_boxes4scale, autoscale
private

integer, parameter :: NSZS=97
integer :: boxsizes(NSZS) = [32, 36, 40, 48, 52, 56, 64, 66, 70, 72, 80, 84, 88, 100, 104, 108, 112, 120, 128, 130, 132,&
140, 144, 150, 160, 162, 168, 176, 180, 182, 192, 200, 208, 216, 220, 224, 240, 256, 264, 288, 300, 308, 320, 324, 336,&
338, 352, 364, 384, 400, 420, 432, 448, 450, 462, 480, 486, 500, 504, 512, 520, 528, 546, 560, 576, 588, 600, 640, 648,&
650, 660, 672, 686, 700, 702, 704, 720, 726, 728, 750, 768, 770, 784, 800, 810, 840, 882, 896, 910, 924, 936, 972, 980,&
1008, 1014, 1020, 1024]

contains

    function find_magic_box( trial_box ) result( best_box )
        integer, intent(in) :: trial_box
        integer :: best_box, dist, ind
        call find(boxsizes, NSZS, trial_box, ind, dist)
        best_box = boxsizes(ind)
    end function find_magic_box

    subroutine print_magic_box_range( smpd, diam )
        real, intent(in) :: smpd, diam
        integer :: ind_start, ind_stop, dist, i
        call find(boxsizes, NSZS, nint((1.5*diam)/smpd), ind_start, dist)
        call find(boxsizes, NSZS, nint((2.0*diam)/smpd), ind_stop, dist)
        do i=ind_start,ind_stop
            write(*,'(i5)') boxsizes(i)
        end do
        if( ind_stop == NSZS ) write(*,'(a)') 'WARNING! box size may underestimated (max value in list is 1024)'
    end subroutine print_magic_box_range

    function find_magic_boxes4scale( orig_box, scales ) result( boxes )
        integer, intent(in) :: orig_box
        real,    intent(in) :: scales(2)
        integer :: boxes(2), ind1, ind2, dist
        call find(boxsizes, NSZS, nint(real(orig_box)*scales(1)), ind1, dist)
        call find(boxsizes, NSZS, nint(real(orig_box)*scales(2)), ind2, dist)
        boxes(1) = boxsizes(ind1)
        boxes(2) = boxsizes(ind2)
    end function find_magic_boxes4scale

    subroutine autoscale( box_in, smpd_in, box_new, smpd_new, scale )
        integer, intent(in)  :: box_in
        real,    intent(in)  :: smpd_in
        integer, intent(out) :: box_new
        real,    intent(out) :: smpd_new, scale
        real,    parameter   :: SMPD_TARGET = 3.7
        if( smpd_in < SMPD_TARGET )then
            ! ok
        else
            scale    = 1.0
            box_new  = box_in
            smpd_new = smpd_in
            write(*,*) 'Inputted smpd < SMPD_TARGET, no scaling done; simple_magic_boxes :: autoscale'
            return
        endif
        scale    = smpd_in/SMPD_TARGET
        box_new  = find_magic_box(nint(scale*real(box_in)))
        scale    = real(box_new)/real(box_in)
        smpd_new = smpd_in/scale
    end subroutine autoscale

end module simple_magic_boxes

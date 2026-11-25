!@descr: the class-average CNN's input: a class average and its halves resampled to the network's pixel size and box, then scaled by the average's statistics, as scripts/sscnn/cnn_data.py does
! The resampling is numpy's Fourier crop or pad with the transforms' origin at the centre pixel, as one exact operator
! per axis, out = Re(R img R^T); the centre stays the centre pixel at odd and even sizes alike. The image class's
! pad and clip are not used: they need even boxes and weight the spectrum against aliasing.
module simple_cavg_cnn_input
use simple_defs,  only: sp, dp, DTWOPI
use simple_error, only: simple_exception
implicit none

public :: cavg_cnn_input, run_all_cavg_cnn_input_tests
private
#include "simple_local_flags.inc"

type :: cavg_cnn_input
    private
    integer               :: n          = 0, m = 0, box = 0
    logical               :: l_resample = .false.
    real(dp), allocatable :: rr(:,:), ri(:,:)       !< (m, n) the operator R's real and imaginary parts
    real(dp), allocatable :: rrt(:,:), rit(:,:)     !< (n, m) their transposes
    logical               :: exists     = .false.
  contains
    procedure          :: new
    procedure          :: prepare
    procedure          :: kill
    procedure, private :: resample
end type cavg_cnn_input

interface
    module subroutine run_all_cavg_cnn_input_tests()
    end subroutine run_all_cavg_cnn_input_tests
end interface

contains

    !> images of n pixels at smpd to the network's box at smpd_out: m = nint(n smpd / smpd_out) pixels after the
    !! Fourier crop or pad, then a crop or pad about the centre pixel to box
    subroutine new( self, n, smpd, box, smpd_out )
        class(cavg_cnn_input), intent(inout) :: self
        integer,               intent(in)    :: n, box
        real,                  intent(in)    :: smpd, smpd_out
        real(dp) :: theta, cs, sn
        integer  :: l, j, q, qlo, qhi
        call self%kill
        if( n < 1 .or. box < 1 .or. smpd <= 0. .or. smpd_out <= 0. ) THROW_HARD('invalid sizes; cavg_cnn_input%new')
        self%n   = n
        self%box = box
        self%m   = nint(real(n, dp) * real(smpd, dp) / real(smpd_out, dp))
        if( self%m < 1 ) THROW_HARD('the images shrink to nothing at the network pixel size; cavg_cnn_input%new')
        self%l_resample = self%m /= n
        if( self%l_resample )then
            ! the frequencies both grids hold: those of the smaller one, centred on index size/2 (numpy's fftshift)
            if( self%m <= n )then
                qlo = -(self%m / 2)
                qhi = self%m - 1 - self%m / 2
            else
                qlo = -(n / 2)
                qhi = n - 1 - n / 2
            endif
            allocate(self%rr(self%m, n), self%ri(self%m, n))
            do j = 0,n-1
                do l = 0,self%m-1
                    theta = DTWOPI * (real(l - self%m/2, dp) / real(self%m, dp) - real(j - n/2, dp) / real(n, dp))
                    cs = 0.d0
                    sn = 0.d0
                    do q = qlo,qhi
                        cs = cs + cos(real(q, dp) * theta)
                        sn = sn + sin(real(q, dp) * theta)
                    enddo
                    self%rr(l+1,j+1) = cs / real(self%m, dp)
                    self%ri(l+1,j+1) = sn / real(self%m, dp)
                enddo
            enddo
            self%rrt = transpose(self%rr)
            self%rit = transpose(self%ri)
        endif
        self%exists = .true.
    end subroutine new

    !> raw(n, n, channels) (channel 1 the class average, then its halves) to the network's input x(box, box,
    !! channels): resampled, then all channels minus the average's mean inside the inscribed circle and over its
    !! standard deviation there; a constant average is empty and its input zero
    subroutine prepare( self, raw, x, empty )
        class(cavg_cnn_input), intent(in)  :: self
        real(sp),              intent(in)  :: raw(:,:,:)
        real(sp),              intent(out) :: x(:,:,:)
        logical,               intent(out) :: empty
        real(dp), allocatable :: t(:,:,:)
        real(dp)              :: mean, sdev, r2, c
        integer               :: ch, i, j, nin
        if( .not. self%exists ) THROW_HARD('not constructed; cavg_cnn_input%prepare')
        if( size(raw,1) /= self%n .or. size(raw,2) /= self%n ) THROW_HARD('raw images of another size; cavg_cnn_input%prepare')
        if( size(x,1) /= self%box .or. size(x,2) /= self%box .or. size(x,3) /= size(raw,3) )&
            &THROW_HARD('input of another shape; cavg_cnn_input%prepare')
        mean  = sum(real(raw(:,:,1), dp)) / real(self%n**2, dp)
        empty = sum((real(raw(:,:,1), dp) - mean)**2) == 0.d0
        if( empty )then
            x = 0.
            return
        endif
        allocate(t(self%box, self%box, size(raw,3)))
        do ch = 1,size(raw,3)
            call self%resample(real(raw(:,:,ch), dp), t(:,:,ch))
        enddo
        ! the statistics of the average inside the circle of radius box/2 about the centre pixel
        c    = real(self%box / 2, dp)
        r2   = (real(self%box, dp) / 2.d0)**2
        mean = 0.d0
        nin  = 0
        do j = 1,self%box
            do i = 1,self%box
                if( (real(i-1, dp) - c)**2 + (real(j-1, dp) - c)**2 < r2 )then
                    mean = mean + t(i,j,1)
                    nin  = nin + 1
                endif
            enddo
        enddo
        mean = mean / real(nin, dp)
        sdev = 0.d0
        do j = 1,self%box
            do i = 1,self%box
                if( (real(i-1, dp) - c)**2 + (real(j-1, dp) - c)**2 < r2 ) sdev = sdev + (t(i,j,1) - mean)**2
            enddo
        enddo
        sdev = sqrt(sdev / real(nin, dp))
        if( sdev <= 0.d0 ) sdev = 1.d0
        do ch = 1,size(raw,3)
            x(:,:,ch) = real((t(:,:,ch) - mean) / sdev, sp)
        enddo
        deallocate(t)
    end subroutine prepare

    subroutine kill( self )
        class(cavg_cnn_input), intent(inout) :: self
        if( allocated(self%rr) ) deallocate(self%rr, self%ri, self%rrt, self%rit)
        self%n   = 0
        self%m   = 0
        self%box = 0
        self%l_resample = .false.
        self%exists     = .false.
    end subroutine kill

    !> img(n, n) to out(box, box): the operator, then a crop or pad about the centre pixel (index size/2 from 0)
    subroutine resample( self, img, out )
        class(cavg_cnn_input), intent(in)  :: self
        real(dp),              intent(in)  :: img(:,:)
        real(dp),              intent(out) :: out(:,:)
        real(dp), allocatable :: t(:,:)
        integer               :: mm, a
        if( self%l_resample )then
            t = matmul(matmul(self%rr, img), self%rrt) - matmul(matmul(self%ri, img), self%rit)
        else
            t = img
        endif
        mm = size(t,1)
        if( mm >= self%box )then
            a   = mm / 2 - self%box / 2
            out = t(a+1:a+self%box, a+1:a+self%box)
        else
            a   = self%box / 2 - mm / 2
            out = 0.d0
            out(a+1:a+mm, a+1:a+mm) = t
        endif
    end subroutine resample

end module simple_cavg_cnn_input

module simple_ft_expanded_tester
include 'simple_lib.f08'
use simple_ft_expanded,  only: ft_expanded
use simple_ftexp_shsrch, only: ftexp_shsrch
use simple_image,        only: image
implicit none

public :: exec_ft_expanded_test
private
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: LDIM(3)=[240,240,1], SQRAD=60, NTST=50, NNOISY=20
real,    parameter :: SMPD=1.77, TRS=10., HP=100.0, LP=8., SNR=0.2

! module global variables
type(ft_expanded)        :: ftexp_img
type(image)              :: img, img_shifted
type(image), allocatable :: noisy_imgs(:)
integer                  :: x, y

contains

    subroutine exec_ft_expanded_test
        type(image)       :: img, img2
        type(ft_expanded) :: ftexp
        call setup_testenv
        call test_shifted_correlator
        call profile_corrs
    end subroutine exec_ft_expanded_test

    subroutine setup_testenv
        integer :: i
        call img%new(LDIM, SMPD)
        call img%square(SQRAD)
        call ftexp_img%new(img, HP, LP, .true.)
        call img_shifted%new(LDIM, SMPD)
        ! seed the random number generator
        call seed_rnd
        ! generate noisy images
        if( allocated(noisy_imgs) )then
            do i=1,NNOISY
                call noisy_imgs(i)%kill
            end do
            deallocate(noisy_imgs)
        else
            allocate(noisy_imgs(NNOISY))
            do i=1,NNOISY
                noisy_imgs(i) = img
                call noisy_imgs(i)%add_gauran(SNR)
            end do
        endif
    end subroutine setup_testenv

    subroutine test_shifted_correlator
        real    :: dist, corr, corravg, distavg
        integer :: itst
        real    :: shvec(3)
        write(logfhandle,*) 'testing ft_expanded :: shifted correlator'
        corravg = 0.
        distavg = 0.
        do itst=1,NTST
            x = nint(ran3()*2*TRS-TRS)
            y = nint(ran3()*2*TRS-TRS)
            img_shifted = img
            call img_shifted%shift([real(x),real(y),0.])
            call find_shift(dist, corr)
            corravg = corravg + corr
            distavg = distavg + dist
        end do
        corravg = corravg/real(NTST)
        distavg = distavg/real(NTST)
        if( corravg > 0.999 .and. distavg < 0.0001 )then
            ! test passed
        else
            THROW_HARD('****ft_expanded_tester FAILURE :: test_shifted_correlator')
        endif

        contains

            subroutine find_shift( dist, corr_best )
                real, intent(out)  :: dist, corr_best
                type(ft_expanded)  :: ftexp_trial
                type(ftexp_shsrch) :: ftexp_shsrch_trial
                real    :: shvec(3), corr
                integer :: ysh_best, xsh_best, xsh, ysh
                call ftexp_trial%new(img_shifted, HP, LP, .true.)
                call ftexp_shsrch_trial%new(ftexp_img, ftexp_trial, TRS)
                corr_best = -huge(corr)
                do xsh=nint(-TRS),nint(TRS)
                    do ysh=nint(-TRS),nint(TRS)
                        shvec(1) = real(xsh)
                        shvec(2) = real(ysh)
                        shvec(3) = 0.
                        corr = real(ftexp_shsrch_trial%corr_shifted_8(dble(shvec)))
                        if( corr > corr_best )then
                            corr_best = corr
                            xsh_best  = xsh
                            ysh_best  = ysh
                        endif
                    end do
                end do
                call ftexp_img%corr_normalize(ftexp_trial, corr_best)
                dist = euclid(real([xsh_best,ysh_best]), real([-x,-y]))
                call ftexp_trial%kill
                call ftexp_shsrch_trial%kill
            end subroutine find_shift

    end subroutine test_shifted_correlator

    subroutine profile_corrs
#ifdef INTEL
        use ifport
#endif

        integer, parameter   :: NTSTS=1000, NTHR=8
        integer              :: itst
        type(image)          :: img_ref, img_ptcl
        type(ft_expanded)    :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch)   :: ftexp_shsrch1
        real, allocatable    :: shvecs(:,:)
        real(4)    :: corr, actual, delta, shvec(3), tarray(2)
        call img_ref%new([4096,4096,1],SMPD)
        call img_ref%ran
        call img_ref%fft()
        call ftexp_ref%new(img_ref, HP, LP, .true.)
        call img_ptcl%new([4096,4096,1],SMPD)
        call img_ptcl%ran
        call img_ptcl%fft()
        call ftexp_ptcl%new(img_ptcl, HP, LP, .true.)
        call ftexp_shsrch1%new(ftexp_ref, ftexp_ptcl, TRS)
        allocate(shvecs(NTSTS,3))
        do itst=1,NTSTS
            shvecs(itst,1) = ran3()*2*TRS-TRS
            shvecs(itst,2) = ran3()*2*TRS-TRS
            shvecs(itst,3) = 0.
        end do
        actual = etime( tarray )
        write(logfhandle,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        delta = dtime( tarray )
        write(logfhandle,'(A,F9.2)') 'Relative cpu-time:', delta

        write(logfhandle,'(a)') '>>> PROFILING STANDARD CORRELATOR'
        do itst=1,NTST
            corr = img_ref%corr_shifted(img_ptcl, shvecs(itst,:), lp_dyn=LP)
        end do
        actual = etime( tarray )
        write(logfhandle,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        delta = dtime( tarray )
        write(logfhandle,'(A,F9.2)') 'Relative cpu-time:', delta
        write(logfhandle,'(a)') '>>> PROFILING FTEXP CORRELATOR'
        do itst=1,NTST
            corr = real(ftexp_shsrch1%corr_shifted_8(dble(shvecs(itst,:))))
        end do
        actual = etime( tarray )
        write(logfhandle,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        delta = dtime( tarray )
        write(logfhandle,'(A,F9.2)') 'Relative cpu-time:', delta
    end subroutine profile_corrs

end module simple_ft_expanded_tester

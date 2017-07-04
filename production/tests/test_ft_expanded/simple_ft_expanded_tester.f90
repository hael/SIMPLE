module simple_ft_expanded_tester
use simple_ft_expanded, only: ft_expanded
use simple_image,       only: image
implicit none

public :: exec_ft_expanded_test
private

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
        use simple_rnd, only: seed_rnd
        integer :: i
        call img%new(LDIM, SMPD)
        call img%square(SQRAD)
        call ftexp_img%new(img, HP, LP)
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
        use simple_rnd
        real    :: dist, corr, corravg, distavg
        integer :: itst
        real    :: shvec(3)
        write(*,*) 'testing ft_expanded :: shifted correlator'
        corravg = 0.
        distavg = 0.
        do itst=1,NTST
            x = nint(ran3()*2*TRS-TRS)
            y = nint(ran3()*2*TRS-TRS)
            call img%shift(real(x), real(y), imgout=img_shifted)
            call find_shift(dist, corr)
            corravg = corravg + corr
            distavg = distavg + dist
        end do
        corravg = corravg/real(NTST)
        distavg = distavg/real(NTST)
        if( corravg > 0.999 .and. distavg < 0.0001 )then
            ! test passed
        else
            stop '****ft_expanded_tester FAILURE :: test_shifted_correlator'
        endif

        contains

            subroutine find_shift( dist, corr_best )
                use simple_math, only: euclid
                real, intent(out) :: dist, corr_best
                type(ft_expanded) :: ftexp_trial
                real    :: shvec(3), corr
                integer :: ysh_best, xsh_best, xsh, ysh
                call ftexp_trial%new(img_shifted, HP, LP)
                corr_best = -1.
                do xsh=nint(-TRS),nint(TRS)
                    do ysh=nint(-TRS),nint(TRS)
                        shvec(1) = real(xsh)
                        shvec(2) = real(ysh)
                        shvec(3) = 0.
                        corr = ftexp_img%corr_shifted(ftexp_trial, shvec)
                        if( corr > corr_best )then
                            corr_best = corr
                            xsh_best  = xsh
                            ysh_best  = ysh
                        endif
                    end do
                end do
                dist = euclid(real([xsh_best,ysh_best]), real([-x,-y]))
            end subroutine find_shift

    end subroutine test_shifted_correlator

    subroutine profile_corrs
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_rnd
        use simple_syscalls
        integer, parameter   :: NTSTS=1000, NTHR=8
        integer              :: itst
        type(image)          :: img_ref, img_ptcl
        type(ft_expanded)    :: ftexp_ref, ftexp_ptcl
        real, allocatable    :: shvecs(:,:)
        real    :: corr, actual, delta, shvec(3)
!$      call omp_set_num_threads(NTHR)
        call img_ref%new([4096,4096,1],SMPD)
        call img_ref%ran
        call img_ref%fwd_ft
        call ftexp_ref%new(img_ref, HP, LP)
        call img_ptcl%new([4096,4096,1],SMPD)
        call img_ptcl%ran
        call img_ptcl%fwd_ft
        call ftexp_ptcl%new(img_ptcl, HP, LP)
        allocate(shvecs(NTSTS,3))
        do itst=1,NTSTS
            shvecs(itst,1) = ran3()*2*TRS-TRS
            shvecs(itst,2) = ran3()*2*TRS-TRS
            shvecs(itst,3) = 0.
        end do
        actual = getabscpu(.true.)
        delta  = getdiffcpu(.true.) 
        write(*,'(a)') '>>> PROFILING STANDARD CORRELATOR' 
        do itst=1,NTST
            corr = img_ref%corr_shifted(img_ptcl, shvecs(itst,:), lp_dyn=LP)
        end do
        actual = getabscpu(.true.)
        delta  = getdiffcpu(.true.)
        write(*,'(a)') '>>> PROFILING FTEXP CORRELATOR' 
        !$omp parallel do schedule(auto) default(shared) private(itst)
        do itst=1,NTST
            corr = ftexp_ref%corr_shifted(ftexp_ptcl, shvecs(itst,:))
        end do
        !$omp end parallel do
        actual = getabscpu(.true.)
        delta  = getdiffcpu(.true.)
    end subroutine profile_corrs

end module simple_ft_expanded_tester

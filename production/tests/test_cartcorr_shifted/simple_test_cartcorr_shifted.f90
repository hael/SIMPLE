program simple_test_cartcorr_shifted
include 'simple_lib.f08'
use simple_cartft_corrcalc,   only: cartft_corrcalc
use simple_builder,           only: builder
use simple_image,             only: image
use simple_parameters,        only: parameters
use simple_cmdline,           only: cmdline
use simple_cftcc_shsrch_grad, only: cftcc_shsrch_grad
implicit none
type(cmdline)           :: cline
type(builder)           :: b
type(parameters)        :: p
type(cartft_corrcalc)   :: cftcc
type(cftcc_shsrch_grad) :: cftcc_shsrch
type(image)             :: noise_img
real, parameter         :: SHMAG=5.0
integer                 :: i, noise_n, noise_i
real                    :: grad(2), cxy(3), lims(2,2)
integer(timer_int_kind) :: t_tot
real,    parameter      :: NOISE_MIN = .3, NOISE_MAX = .7, NOISE_DEL = 0.1
real                    :: ave, sdev, maxv, minv, noise_lvl, correct_sh(2)
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_cartcorr_shifted stk=<particles.ext>'
    write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] mskdiam=zz [verbose=<yes|no{no}>]'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('stk',     1)
call cline%checkvar('smpd',    2)
call cline%checkvar('mskdiam', 3)
call cline%check
call p%new(cline)
p%kfromto(1) = 1
p%kfromto(2) = p%box/2 - 5
call b%build_general_tbox(p, cline)
call cftcc%new_dev([1,p%nptcls])
do i = 1,p%nptcls
    call b%img%read(p%stk, i)
    call b%img%fft
    call cftcc%set_ptcl(i, b%img)
end do
! t_tot = tic()
! do i = 1,1
!     call b%img%read(p%stk, i)
!     call b%img%fft
!     call b%img%shift2Dserial([-0.5,-0.5]) 
!     call cftcc%set_ref(b%img)
!     call srch_shifts(0.5, i)
! end do
! print *, 'corr_shifted timing = ', toc(t_tot)
! t_tot = tic()
! do i = 1,1
!     call b%img%read(p%stk, i)
!     call b%img%fft
!     call b%img%shift2Dserial([-0.5,-0.5]) 
!     call cftcc%set_ref(b%img)
!     call srch_shifts_ad(0.5, i)
! end do
! print *, 'corr_shifted_ad timing = ', toc(t_tot)
! lbfgsb shift search
lims(1,1) = -6.
lims(1,2) =  6.
lims(2,1) = -6.
lims(2,2) =  6.
call cftcc_shsrch%new(lims)
noise_n  = int((NOISE_MAX - NOISE_MIN)/NOISE_DEL + 1.)
call noise_img%new(p%ldim, p%smpd)
do noise_i = 1, noise_n
    noise_lvl = NOISE_MIN + (noise_i - 1)*NOISE_DEL
    print *, 'noise = ', noise_lvl
    do i = 1, p%nptcls
        correct_sh(1) = 2.*(ran3()-0.5) * 5
        correct_sh(2) = 2.*(ran3()-0.5) * 5
        call b%img%read(p%stk, i)
        call b%img%mask(p%msk, 'soft')
        call b%img%stats('foreground', ave, sdev, maxv, minv)
        call noise_img%gauran(0., noise_lvl * sdev)
        call noise_img%mask(1.5 * p%msk, 'soft')
        call b%img%add(noise_img)
        call b%img%write("stk_img_noisy_"//int2str(noise_i)//".mrc", i)
        call b%img%fft
        call b%img%shift2Dserial(correct_sh) 
        call cftcc%set_ref(b%img)
        call cftcc_shsrch%set_pind(i)
        cxy = cftcc_shsrch%minimize()
        print *, 'iptcl = ', i, ': minimized shift = ', cxy(2:), '; corr = ', cxy(1), '; correct shift = ', correct_sh
    enddo
enddo

contains

    subroutine srch_shifts_ad( shstep, iptcl )
        real,    intent(in) :: shstep
        integer, intent(in) :: iptcl
        real, allocatable   :: srch_space(:,:) 
        integer :: cnt, i
        real    :: x, y, corr
        cnt = 0
        x   = -SHMAG
        do while( x <= SHMAG )
            y = -SHMAG
            do while( y <= SHMAG )
                cnt = cnt + 1
                y = y + shstep
            end do
            x = x + shstep
        end do
        allocate(srch_space(cnt,2), source=0.)
        cnt = 0
        x   = -SHMAG; 
        do while( x <= SHMAG )
            y = -SHMAG
            do while( y <= SHMAG )
                cnt = cnt + 1
                srch_space(cnt,1) = x
                srch_space(cnt,2) = y
                y = y + shstep
            end do
            x = x + shstep
        end do
        do i = 1,cnt
            corr = cftcc%corr_shifted_ad(iptcl, srch_space(i,:), grad)
            print *, srch_space(i,1), srch_space(i,2), corr, grad(1), grad(2)
        end do
    end subroutine srch_shifts_ad

    subroutine srch_shifts( shstep, iptcl )
        real,    intent(in) :: shstep
        integer, intent(in) :: iptcl
        real, allocatable   :: srch_space(:,:) 
        integer :: cnt, i
        real    :: x, y, corr
        cnt = 0
        x   = -SHMAG
        do while( x <= SHMAG )
            y = -SHMAG
            do while( y <= SHMAG )
                cnt = cnt + 1
                y = y + shstep
            end do
            x = x + shstep
        end do
        allocate(srch_space(cnt,2), source=0.)
        cnt = 0
        x   = -SHMAG; 
        do while( x <= SHMAG )
            y = -SHMAG
            do while( y <= SHMAG )
                cnt = cnt + 1
                srch_space(cnt,1) = x
                srch_space(cnt,2) = y
                y = y + shstep
            end do
            x = x + shstep
        end do
        do i = 1,cnt
            corr = cftcc%corr_shifted(iptcl, srch_space(i,:), grad)
            print *, srch_space(i,1), srch_space(i,2), corr, grad(1), grad(2)
        end do
    end subroutine srch_shifts

end program simple_test_cartcorr_shifted

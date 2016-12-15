program simple_test_2Dwiener
use simple_jiffys        ! singleton
use simple_filterer      ! singleton
use simple_estimate_ssnr ! singleton
use simple_defs          ! singleton
use simple_math          ! singleton 
use simple_stat          ! singleton
use simple_cmdline,      only: cmdline
use simple_params,       only: params
use simple_build,        only: build
use simple_image,        only: image
implicit none
type(params)              :: p
type(build)               :: b
type(cmdline)             :: cline
type(image)               :: img_rec
real, allocatable         :: ssnr(:), sqrtssnr(:), ssnr_avg(:), frc(:) 
integer                   :: iptcl, alloc_stat, k, nk
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_TEST_2DWIENER stk=<stack.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' msk=<mask radius(in pixels)> deftab=<input defocus doc>'
    write(*,'(a)',advance='no') ' [ctf=<yes|no|flip|mul{yes}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)',advance='no') ' [fraca=<frac amp contrast{0.07}>] [cs=<spherical aberration'
    write(*,'(a)') ' constant(in mm){2.7}>] [nthr=<nr of OpenMP threads{1}>] [mw=<molecular weight(in kD)>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',    1)
call cline%checkvar('smpd',   2)
call cline%checkvar('msk',    3)
call cline%checkvar('deftab', 4)
if( .not. cline%defined('ctf') ) call cline%set('ctf', 'yes' )
call cline%check
p = params(cline, .false.)
call b%build_general_tbox(p, cline)
allocate(b%imgs(p%nptcls), stat=alloc_stat)
call alloc_err("In: simple_test_2Dwiener, 2", alloc_stat)
do iptcl=1,p%nptcls
    call b%img%read(p%stk, iptcl)
    b%imgs(iptcl) = b%img
end do
call wiener_restore2D(b%imgs, b%a, b%tfun, p%ctfmode, p%ctf, img_rec, p%msk)
call img_rec%write('restored_img.spi')
! ssnr = estimate_ssnr(b%imgs, b%a, p%msk, b%tfun, p%ctfmode, p%ctf)
! frc  = ssnr2fsc(ssnr)
! do k=1,size(ssnr)
!     print *, b%img%get_lp(k), frc(k)
! end do
! allocate(sqrtssnr(size(ssnr)), source=ssnr)
! where(ssnr > 0.)
!     sqrtssnr = sqrt(ssnr)
! else where
!     sqrtssnr = 0.
! end where
!
! do iptcl=1,p%nptcls
!     b%img_tmp = img_rec
!     call prep_ref(b%img_tmp, sqrtssnr)
!     call b%img_tmp%write('references.spi', iptcl)
!     call b%img%read(p%stk, iptcl)
!     call prep_ptcl(b%img, ssnr)
!     call b%img%write('ptcls.spi', iptcl)
! end do

contains
    
    subroutine prep_ptcl( pimg, ssnr )
        class(image), intent(inout) :: pimg
        real, intent(in)            :: ssnr(:)
        type(image)                 :: filter
        real                        :: x, y
        call filter%new([p%box,p%box,1], p%smpd)
        ! prepare the filter
        call b%tfun%ctf2img(filter, b%a%get(iptcl,'dfx'), 'square', b%a%get(iptcl,'dfy'), b%a%get(iptcl,'angast'))
        call filter%apply_filter(ssnr)
        call filter%add(1.0)
        call filter%square_root
        ! prepare the particle image
        call pimg%fwd_ft
        call pimg%shellnorm
        call pimg%apply_filter(filter)
        x = b%a%get(iptcl,'x')
        y = b%a%get(iptcl,'y')
        call b%img%shift(-x, -y)
        call pimg%bwd_ft
        call b%img%rtsq(-b%a%e3get(iptcl), 0., 0.)
        call pimg%mask(p%msk, 'soft')
        call filter%kill
    end subroutine

    subroutine prep_ref( rimg, sqrtssnr )
        class(image), intent(inout) :: rimg
        real, intent(in)            :: sqrtssnr(:)
        ! prepare the reference image
        call rimg%fwd_ft
        call rimg%shellnorm
        call rimg%apply_filter(sqrtssnr)
        call rimg%bwd_ft
        call b%tfun%apply(rimg, b%a%get(iptcl,'dfx'), 'ctf', b%a%get(iptcl,'dfy'), b%a%get(iptcl,'angast'))
        call rimg%bwd_ft
        call rimg%mask(p%msk, 'soft')
    end subroutine

end program simple_test_2Dwiener

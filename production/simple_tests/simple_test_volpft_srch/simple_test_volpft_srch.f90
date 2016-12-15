program simple_test_volpft_srch
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_build,           only: build
use simple_params,          only: params
use simple_ori,             only: ori
use simple_image,           only: image
use simple_rnd,             only: ran3
use simple_math,            only: euclid
use simple_cmdline,         only: cmdline
use simple_volpft_srch      ! singleton
use simple_defs             ! singleton
use simple_gridding         ! singleton
use simple_jiffys           ! singleton
implicit none
type(params)          :: p
type(build)           :: b
type(cmdline)         :: cline
type(volpft_corrcalc) :: vpftcc
type(ori)             :: o_best, e, ranori
type(image)           :: vol_ref
real                  :: corr_best, shvec(3), x, y, z, dist, sumdist, sherr, xf, yf, zf
integer, parameter    :: NTESTS=10, NPEAKS=3
real, parameter       :: TRS=5.
integer               :: itest
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_VOLOPS vol1=<reference.ext> vol2=<target.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' msk=<mask radius(in pixels)> snr=<signal to noise ratio> [hp=<high-pass limit{100}>]'
    write(*,'(a)') ' [lp=<low-pass limit{20}>] [nthr=<nr of openMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('vol1', 1)
call cline%checkvar('vol2', 2)
call cline%checkvar('smpd', 3)
call cline%checkvar('msk',  4)
call cline%checkvar('snr',  5)
call cline%check
p = params(cline,checkdistr=.false.) ! constants & derived constants produced, mode=2
call b%build_general_tbox(p,cline)   ! general objects built
! deal with reference 
call vol_ref%new([p%box,p%box,p%box], p%smpd)
call vol_ref%read(p%vols(1))
call vol_ref%add_gauran(p%snr)
call vol_ref%mask(p%msk,'soft')
call vol_ref%fwd_ft
! deal with target
call b%vol%new([p%box,p%box,p%box], p%smpd)
call b%vol%read(p%vols(2))
call b%vol%add_gauran(p%snr)
call b%vol%mask(p%msk,'soft')
call b%vol%fwd_ft
! create correlator object
call vpftcc%new(vol_ref,b%vol,p%hp,p%lp)
! initialise srch object
call volpft_srch_init(vpftcc, 5.)
sumdist = 0.
sherr   = 0.
print *, 'running ', NTESTS, ' tests'
do itest=1,NTESTS
    call progress(itest,NTESTS)
    call ranori%rnd_ori
    x = ran3()*2*TRS-TRS
    y = ran3()*2*TRS-TRS
    z = ran3()*2*TRS-TRS
    call ranori%set('x',x)
    call ranori%set('y',y)
    call ranori%set('z',z)
    call vpftcc%extract_ref(ranori)
    call vpftcc%shift_orig_ref([x,y,z])
    call volpft_6dimsrch(NPEAKS, corr_best, o_best)
    dist    = o_best.euldist.ranori
    sumdist = sumdist+dist
    xf = o_best%get('x')
    yf = o_best%get('y')
    zf = o_best%get('z')
    sherr = sherr+euclid([x,y,z],[-xf,-yf,-zf])
end do
dist  = sumdist/real(NTESTS)
sherr = sherr/real(NTESTS*3)
print *, 'angular error: ', dist
print *, 'shift error: ', sherr
end program simple_test_volpft_srch

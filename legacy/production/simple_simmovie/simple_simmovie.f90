!==Program simple_simmovie
!
! <simmovie/begin> is a program for simulating a DDD movie. Input is a set of projection images to place. Movie frames
! are then generated related by randomly shifting the base image and applying three different noise sources. Shot and
! detector noise is simulated as well as fixed pattern noise due to dead or hot pixels. This is necessary to simulate
! since any movie alignment procedure must overcome the correlation peak bias in at (0,0) due to this fixed pattern
! noise. <simmovie/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_simmovie
use simple_jiffys,      ! singleton
use simple_cmdline,     only: cmdline
use simple_build,       only: build
use simple_params,      only: params
use simple_ori,         only: ori
use simple_math,        only: deg2rad, gen_ptcl_pos
use simple_image,       only: image
use simple_rnd,         only: ran3
use simple_procimgfile, only: stats_imgfile
implicit none
type(build)          :: b
type(params)         :: p
type(cmdline)        :: cline
type(image)          :: base_image, shifted_base_image
real                 :: snr_pink, snr_detector, ave, sdev, var, med, fracarea, x, y, sherr, dfx, dfy, deferr, angast
integer              :: i, ptclarea, mgrapharea, fixed_frame, alloc_stat
integer, allocatable :: ptcl_positions(:,:)
real, allocatable    :: shifts(:,:)
logical              :: here, debug=.false.
if( command_argument_count() < 1 )then
    write(*,'(a)', advance='no') 'SIMPLE_SIMMOVIE stk=<projs2place.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' msk=<mask radius(in pixels)> xdim=<x dimension(in pixles)> ydim=<y'
    write(*,'(a)', advance='no') ' dimension(in pixles)> snr=<per frame signal2noise ratio> [nframes=<number of'
    write(*,'(a)', advance='no') ' frames{30}>] [trs=<shift'
    write(*,'(a)', advance='no') ' error(in pixles){3.}>] [kv=<acceleration voltage(in kV){300.}>] [fraca=<frac'
    write(*,'(a)', advance='no') ' amp contrast{0.07}>] [cs=<spherical aberrationconstant(in mm){2.7}>]'
    write(*,'(a)', advance='no') ' [defocus=<defocus(in microns){3.}>] [bfac=<bfactor(in A**2){200}>]'
    write(*,'(a)') ' [vis=<yes|no{no}>]'
    stop
endif
call cline%parse
call cline%checkvar('stk',  1)
call cline%checkvar('msk',  2)
call cline%checkvar('xdim', 3)
call cline%checkvar('ydim', 4)
call cline%checkvar('snr',  5)
if( .not. cline%defined('trs') )         call cline%set('trs', 3.)
if( .not. cline%defined('ctf') )         call cline%set('ctf', 'yes' )
if( .not. cline%defined('bfac') )        call cline%set('bfac', 200.)
if( .not. cline%defined('nframes') )     call cline%set('nframes', 30.)
call cline%set('wfun', 'kb')
call cline%set('winsz', 1.5)
call cline%set('alpha', 2.)
call cline%set('eo', 'no')
call cline%check
call cline%set('prg', 'simmovie')
p = params(cline) ! parameters generated
if( p%box == 0 ) stop 'box=0, something is fishy!'
call b%build_general_tbox(p, cline) ! build generated
! set fixed frame
fixed_frame = nint(real(p%nframes)/2.)
! remake the alignment doc
call b%a%new(1)
! check the fractional area occupied by particles & generate particle positions
ptclarea   = p%box*p%box*p%nptcls
mgrapharea = p%xdim*p%ydim
fracarea   = real(ptclarea)/real(mgrapharea)
write(*,'(a,1x,f7.3)') 'Fraction of area occupied by ptcls:', fracarea
if( fracarea > 0.55 )then
    write(*,'(A)') 'It is not recommended that more than 55% of the micrograph area is occupied with particles!'
    write(*,'(A)') 'Please, reduce the number of projection images to place!'
    stop
endif
write(*,'(a)') '>>> GENERATING PARTICLE POSITIONS'
ptcl_positions = gen_ptcl_pos(p%nptcls, p%xdim, p%ydim, p%box)
! calculate stack background stats
inquire(file=p%stk, exist=here)
if( .not. here )then
    write(*,*) 'file does not exis:', p%stk
    stop
endif
call stats_imgfile(p%stk, 'background', ave, sdev, var, med, p%msk)
write(*,'(a)') '>>> STACK BACKGROUND STATISTICS'
write(*,'(a,1x,f7.4)') 'ave: ', ave
write(*,'(a,1x,f7.4)') 'sdev:', sdev
write(*,'(a,1x,f7.4)') 'med: ', med
! make a base image by inserting the projections at ptcl_positions and setting the background to med
call base_image%new([p%xdim,p%ydim,1], p%smpd, backgr=med)
do i=1,p%nptcls
    call b%img%read(p%stk, i)
    call b%img%insert(ptcl_positions(i,:), base_image)
end do
if( p%vis .eq. 'yes' ) call base_image%vis
if( debug ) write(*,'(a)') 'inserted projections'
! calculate snr:s
snr_pink     = p%snr/0.2
snr_detector = p%snr/0.8
if( debug ) write(*,'(a)') 'calculated SNR:s'
! set CTF parameters
deferr = ran3()*0.2
if( ran3() < 0.5 )then
    dfx = p%defocus-deferr
else
    dfx = p%defocus+deferr
endif
deferr = ran3()*0.2
if( ran3() < 0.5 )then
    dfy = p%defocus-deferr
else
    dfy = p%defocus+deferr
endif
angast = ran3()*360.
if( debug ) write(*,'(a)') 'did set CTF parameters'
if( debug ) write(*,'(a)') 'initialized shifts'
! generate shifts
allocate( shifts(p%nframes,2), stat=alloc_stat )
call alloc_err('In: simple_simmovie; shifts', alloc_stat)
x = 0.
y = 0.
do i=1,p%nframes
    ! generate random shifts
    sherr = ran3()*p%trs
    if( ran3() > 0.5 )then
        x = x+sherr
    else
        x = x-sherr
    endif
    sherr = ran3()*p%trs
    if( ran3() > 0.5 )then
        y = y+sherr
    else
        y = y-sherr
    endif
    shifts(i,1) = x
    shifts(i,2) = y
end do
! put the central frame in the series at x,y = (0,0) & fill-up b%a
do i=1,p%nframes
    shifts(i,:) = shifts(i,:)-shifts(fixed_frame,:)
    call b%a%set(1, 'x'//int2str(i), shifts(i,1))
    call b%a%set(1, 'y'//int2str(i), shifts(i,2))
end do
! make and open a stack for the movie frames
if( debug ) write(*,'(a)') 'made stack for output movie frames'
write(*,'(a)') '>>> GENERATING MOVIE FRAMES'
do i=1,p%nframes
    call progress(i,p%nframes)
    ! shift base image
    call base_image%fwd_ft
    if( i > 1 ) then
        call base_image%shift(shifts(i,1),shifts(i,2),imgout=shifted_base_image)
    else
        shifted_base_image = base_image
    endif
    call shifted_base_image%bwd_ft
    ! add pink noise
    call shifted_base_image%add_gauran(snr_pink)
    ! multiply with CTF
    call shifted_base_image%fwd_ft
    if( p%neg .eq. 'yes' )then
        call b%tfun%apply(shifted_base_image, dfx, 'neg', dfy, angast, p%bfac)
    else
        call b%tfun%apply(shifted_base_image, dfx, 'ctf', dfy, angast, p%bfac)
    endif
    call shifted_base_image%bwd_ft
    ! add the detector noise
    call shifted_base_image%add_gauran(snr_detector)
    if( p%vis .eq. 'yes' ) call shifted_base_image%vis
    call shifted_base_image%write('simmovie'//p%ext, i)
    ! set orientation parameters in object
    call b%a%set(1, 'dfx',    dfx)
    call b%a%set(1, 'dfy',    dfy)
    call b%a%set(1, 'angast', angast)
end do
if( debug ) write(*,'(a)') 'generated movie'
! generate the optimal average
base_image = 0.
write(*,'(a)') '>>> GENERATING OPTIMAL AVERAGE'
do i=1,p%nframes
    call progress(i,p%nframes)
    call shifted_base_image%read('simmovie'//p%ext, i)
    x = b%a%get(1, 'x'//int2str(i))
    y = b%a%get(1, 'y'//int2str(i))
    call shifted_base_image%shift(-x,-y)
    call base_image%add(shifted_base_image)
end do
if( debug ) write(*,'(a,1x,f7.4)') 'constant 4 division:', real(p%nframes)
call base_image%div(real(p%nframes))
if( debug ) write(*,'(a)') 'generated optimal average'
call base_image%write('optimal_movie_average'//p%ext, 1)
if( p%vis .eq. 'yes' ) call base_image%vis
! output orientations
call b%a%write('simmovie_params.txt')
call simple_end('**** SIMPLE_SIMMOVIE NORMAL STOP ****')
end program simple_simmovie

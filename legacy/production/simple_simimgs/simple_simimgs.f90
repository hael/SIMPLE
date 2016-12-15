!==Program simple_simimgs
!
! <simimgs/begin> is a program for simulating cryo-EM images. It is not a very sophisticated simulator, but it is nevertheless 
! useful for testing purposes. It does not do any multi-slice simulation and it cannot be used for simulating molecules 
! containing heavy atoms. It does not even accept a PDB file as an input. Input is 
! a cryo-EM map, which we usually generate from a PDB file using EMAN's program \texttt{pdb2mrc}. \prgname{simple\_simimgs}
! then projects the volume using Fourier interpolation, applies 20\% of the total noise to the images (pink noise), Fourier 
! transforms them, and multiplies them with astigmatic CTF and B-factor. The images are inverse FTed before the remaining 80\% 
! of the noise (white noise) is added. <simimgs/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_simimgs
use simple_jiffys,   ! singleton
use simple_cmdline,  only: cmdline
use simple_build,    only: build
use simple_params,   only: params
use simple_ori,      only: ori
use simple_rnd,      only: ran3
use simple_math,     only: deg2rad
use simple_gridding, only: prep4cgrid
use simple_timing
implicit none 
type(build)   :: b
type(ori)     :: orientation
type(params)  :: p
type(cmdline) :: cline
real          :: snr_pink, snr_detector, bfac, bfacerr, dfx, dfy, angast
integer       :: i, cnt, ntot
logical       :: debug=.false.
!Start of the execution commands
call timestamp()
call start_Alltimers_cpu()
!parameter
if( command_argument_count() < 1 )then
    write(*,'(a)', advance='no') 'SIMPLE_SIMIMGS vol1=<invol.ext> smpd=<sampling distance(in A)> msk=<mask radius(in pixels)>'
    write(*,'(a)', advance='no') ' nptcls=<number of particles> snr=<signal2noise ratio> [sherr=<shift error(in pixels){2}>]'
    write(*,'(a)', advance='no') ' [ctf=<yes|no|flip|mul{yes}>] [kv=<acceleration voltage(in kV){300.}>] [fraca=<frac amp'
    write(*,'(a)', advance='no') ' contrast{0.07}>] [cs=<spherical aberration constant(in mm){2.7}>] [defocus=<defocus'
    write(*,'(a)', advance='no') '(in microns){3.0}>] [deferr=<defocus error(in microns){1.0}>] [astigerr=<astigmatism'
    write(*,'(a)', advance='no') ' error(in microns){0.1}>] [bfac=<bfactor(in A**2){0}>] [bfacerr=<bfactor error'
    write(*,'(a)', advance='no') '(in A**2){0}>] [oritab=<input alignment doc>] [outfile=<output alignment doc{simoris.txt}>]' 
    write(*,'(a)', advance='no') ' [outstk=<ouput stack{simimgs.ext}>] [nthr=<nr of OpenMP threads{1}>] [single=<yes|no{no}>]'
    write(*,'(a)') ' [ndiscrete=<nr of discrete orientations>]'
    stop
endif
call cline%parse
call cline%checkvar('vol1',   1)
call cline%checkvar('smpd',   2)
call cline%checkvar('msk',    3)
call cline%checkvar('nptcls', 4)
call cline%checkvar('snr',    5)
call cline%set('nspace', cline%get_rarg('nptcls'))
if( .not. cline%defined('sherr') .and. .not. cline%defined('oritab') ) call cline%set('sherr', 2.)
if( .not. cline%defined('ctf') )      call cline%set('ctf', 'yes' )
if( .not. cline%defined('astigerr') ) call cline%set('astigerr', 0.1)
if( .not. cline%defined('bfacerr') )  call cline%set('bfacerr', 0.)
call cline%set('wfun', 'kb')
call cline%set('winsz', 1.5)
call cline%set('alpha', 2.)
call cline%set('eo', 'no')
call cline%check
call cline%set('prg', 'simimgs')
p = params(cline, .false.) ! parameters generated
call b%build_general_tbox(p, cline)
if( .not. cline%defined('outstk') ) p%outstk = 'simimgs'//p%ext
if( cline%defined('part') )then
    if( .not. cline%defined('outfile') ) stop 'need unique output file for parallel jobs'
else
    if( .not. cline%defined('outfile') ) p%outfile = 'simoris.txt'
endif
if( p%box == 0 ) stop 'box=0, something is fishy! Perhaps forgotten to input volume or stack?'
! generate orientation/CTF parameters
if( cline%defined('ndiscrete') )then
    if( p%ndiscrete > 0 )then
        call b%a%rnd_oris_discrete(p%ndiscrete, p%nsym, p%eullims)
    endif
    call b%a%rnd_inpls(p%trs)
else if( .not. cline%defined('oritab') .and. p%single .eq. 'no' )then
    call b%a%rnd_oris(p%sherr)
endif
if( debug )then
    write(*,*) 'CTF parameters used'
    write(*,*) 'kv = ', p%kv
    write(*,*) 'cs = ', p%cs
    write(*,*) 'fraca = ', p%fraca
endif
if( .not. b%a%isthere('dfx') )then
    if( p%ctf  .ne. 'no' ) call b%a%rnd_ctf(p%defocus, p%deferr, p%astigerr)
endif
if( debug ) write(*,'(A)') '>>> DONE GENERATING ORIENTATION/CTF PARAMETERS'
call b%a%write(p%outfile)
! calculate snr:s
snr_pink = p%snr/0.2
snr_detector = p%snr/0.8
if( debug ) write(*,'(A)') '>>> DONE CALCULATING SNR:S'
! prepare for image generation
call b%vol%read(p%vols(1))
if( debug ) write(*,'(A)') '>>> DID READ VOL'
call prep4cgrid(b%vol, b%vol_pad, p%msk, wfuns=b%proj%get_wfuns())
if( debug ) write(*,'(A)') '>>> DONE PREPARING FOR IMAGE GENERATION'
write(*,'(A)') '>>> GENERATING IMAGES'
cnt = 0
ntot = p%top-p%fromp+1
do i=p%fromp,p%top
    cnt = cnt+1
    call progress(cnt,ntot)
    ! zero images
    b%img_pad = cmplx(0.,0.)
    b%img = 0.
    ! extract ori
    orientation = b%a%get_ori(i)
    ! project vol
    call b%proj%fproject(b%vol_pad, orientation, b%img_pad)
    ! shift
    call b%img_pad%shift(orientation%get('x'),orientation%get('y'))
    ! back FT
    call b%img_pad%bwd_ft
    ! add pink noise
    if( p%snr < 3. ) call b%img_pad%add_gauran(snr_pink)
    call b%img_pad%fwd_ft
    ! apply ctf/bfactor
    bfacerr = ran3()*p%bfacerr
    if( ran3() < 0.5 )then
        bfac = p%bfac-bfacerr
    else
        bfac = p%bfac+bfacerr
    endif
    if( orientation%isthere('dfx') .and. orientation%isthere('dfy') )then
        dfx = orientation%get('dfx')
        dfy = orientation%get('dfy')
        angast = orientation%get('angast')
        if( cline%defined('bfac') )then
            call b%tfun%apply(b%img_pad,dfx, 'ctf', dfy, angast, bfac=bfac)
        else
            call b%tfun%apply(b%img_pad,dfx, 'ctf', dfy, angast)
        endif
    else if( orientation%isthere('dfx') )then
        dfx = orientation%get('dfx')
        dfy = dfx
        angast = 0.
        if( cline%defined('bfac') )then
            call b%tfun%apply(b%img_pad, orientation%get('dfx'), 'ctf', bfac=bfac)
        else
            call b%tfun%apply(b%img_pad, orientation%get('dfx'), 'ctf')
        endif
    else
        if( cline%defined('bfac') ) call b%img_pad%apply_bfac(bfac)
    endif
    ! add detector noise
    call b%img_pad%bwd_ft
    if( p%snr < 3. ) call b%img_pad%add_gauran(snr_detector)
    if( p%ctf .eq. 'flip' )then
        ! simulate phase-flipped images
        call b%img_pad%fwd_ft
        call b%tfun%apply(b%img_pad, dfx, 'flip', dfy, angast)
        call b%img_pad%bwd_ft
    else if( p%ctf .eq. 'mul' )then
        ! simulate CTF multiplied images
        call b%img_pad%fwd_ft
        call b%tfun%apply(b%img_pad, dfx, 'ctf', dfy, angast)
        call b%img_pad%bwd_ft
    endif
    ! clip
    call b%img_pad%clip(b%img)
    ! write to stack
    if( cline%defined('part') )then
        call b%img%write('simimgs_part'//int2str_pad(p%part,p%numlen)//p%ext, cnt)
    else
        call b%img%write(p%outstk, i)
    endif
end do
call simple_end('**** SIMPLE_SIMIMGS NORMAL STOP ****')
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
!shutting down the timers
call stop_Alltimers_cpu()

end program simple_simimgs

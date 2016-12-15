program simple_postproc_vol
use simple_masker,       only: automask
use simple_build,        only: build
use simple_params,       only: params
use simple_cmdline,      only: cmdline
use simple_estimate_ssnr ! singleton
use simple_defs          ! singleton
use simple_jiffys        ! singleton
implicit none
type(params)      :: p
type(build)       :: b
type(cmdline)     :: cline
real, allocatable :: fsc(:), spec_count3D(:), tmparr(:), pssnr_ctfsq3D(:), pssnr3D(:), sqrtssnr(:), optlp(:)
integer           :: k, state=1
if( command_argument_count() < 4 )then
    write(*,'(a)',advance='no') 'SIMPLE_POSTPROC_VOL vol1=<invol.ext> smpd=<sampling distance(in A)>'
    write(*,'(a)',advance='no') ' fsc=<fsc_state01.bin> msk=<mask radius(in pixels)> [ctfsqspec=<ctfsqspec_state01.bin>]'
    write(*,'(a)',advance='no') ' [mw=<molecular weight(in kD)>] [bfac=<bfactor(in A**2){200.}>] [automsk=<yes|no{no}>]'
    write(*,'(a)',advance='no') ' [amsklp=<low-pass limit(in A){20}>] [edge=<edge size for softening molecular'
    write(*,'(a)',advance='no') ' envelope(in pixels){14}>] [binwidth=<number of layers before softening molecular'
    write(*,'(a)') ' envelope(in pixels){1}>]'
    stop
endif
call cline%parse
call cline%checkvar('vol1', 1)
call cline%checkvar('smpd', 2)
call cline%checkvar('fsc',  3)
call cline%check
call cline%set('prg', 'postproc_vol')
p = params(cline, checkdistr=.false.) ! constants & derived constants produced, mode=2
call b%build_general_tbox(p, cline)   ! general objects built
if( file_exists(p%fsc) )then
    fsc = file2rarr(p%fsc)
    ! allocate CTF**2-dependent PSSNR term
    allocate(pssnr_ctfsq3D(size(fsc)))
    if( file_exists(p%ctfsqspec) )then
        ! get the count spectrum
        spec_count3D = b%vol%spectrum('count')
        ! get the ctfsq spectrum
        tmparr = file2rarr(p%ctfsqspec)
        ! calculate the CTF**2-dependent component of the PSSNR
        where( tmparr > 1e-6 )
            pssnr_ctfsq3D = spec_count3D/tmparr
        else where
            pssnr_ctfsq3D = 0.
        end where
    else
        pssnr_ctfsq3D = 1.
    endif
    allocate(pssnr3D(size(fsc)), sqrtssnr(size(fsc)))
    pssnr3D = estimate_pssnr3D(p%avr, fsc)*pssnr_ctfsq3D    
    optlp = ssnr2optlp(pssnr3D)
else
    write(*,*) 'FSC file: ', trim(p%fsc), ' not in cwd'
    stop
endif
call b%vol%read(p%vols(state))
call b%vol%fwd_ft
if( cline%defined('bfac') )then
    call b%vol%apply_bfac(p%bfac)
endif
call b%vol%apply_filter(optlp)
call b%vol%bwd_ft
p%vols_msk(state) = add2fbody(trim(p%vols(state)), p%ext, 'msk')
if( p%automsk .eq. 'yes' )then
    p%masks(state) = 'automask_state'//int2str_pad(state,2)//p%ext
    call automask(b, p, cline, b%vol, b%mskvol, p%vols_msk(state), p%masks(state))
else
    call b%vol%mask(p%msk, 'soft')
endif
p%outvol = add2fbody(trim(p%vols(state)), p%ext, 'pproc')
call b%vol%write(p%outvol)
call simple_end('**** SIMPLE_POSTPROC_VOL NORMAL STOP ****')
end program simple_postproc_vol

module simple_classrefine_tester
use simple_build,       only: build
use simple_params,      only: params
use simple_cmdline,     only: cmdline
use simple_simulator,   only: simimg
use simple_ori,         only: ori
use simple_classrefine, only: classrefine
use simple_rnd          ! singleton
use simple_defs         ! singleton
implicit none

public :: exec_classrefine_test
private

! module global constants
integer, parameter :: NIMGS=50
real, parameter    :: TRS=5.0, KV=300., CS=2.7, FRACA=0.07, DEFOCUS=2.5, DFERR=1.0, ASTIGERR=0.2, BFAC=50.
real, parameter    :: ROERR=10.0, SNR=0.2, SNR_PINK=SNR/0.2, SNR_DETECTOR=SNR/0.8, LPLIM=20.
character(len=3),      parameter :: CTFFLAG='yes'
character(len=STDLEN), parameter :: CTFMODE='astig'

! module global variables
type(params)      :: p
type(build)       :: b
type(classrefine) :: crefine
logical           :: verbose=.false.

contains

    subroutine exec_classrefine_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        call setup_testenv( cline, be_verbose )
        write(*,*) '****classrefine_test, init'
        call crefine%refine_master(LPLIM)
        b%img = crefine%get_avg()
        call b%img%write('classrefine_tester_avg.mrc', 1)
        write(*,*) '****classrefine_test, completed'
    end subroutine exec_classrefine_test

    subroutine setup_testenv( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        integer   :: i
        real      :: e3, x, y
        type(ori) :: o
        verbose = .false.
        call cline%set('ctf', CTFFLAG)
        call cline%set('trs', TRS)
        ! create parameters and build
        p = params(cline)                   ! parameters generated
        call b%build_general_tbox(p, cline) ! general objects built
        p%tfplan%mode = CTFMODE
        ! read image
        call b%img%read(p%stk, 1)
        ! simulate orientations
        call b%a%new(NIMGS)
        call b%a%rnd_ctf(KV, CS, FRACA, DEFOCUS, DFERR, ASTIGERR)
        do i=1,NIMGS
            e3 = ran3()*2.*ROERR-ROERR
            x  = ran3()*2.0*TRS-TRS
            y  = ran3()*2.0*TRS-TRS
            call b%a%set(i, 'x', x)
            call b%a%set(i, 'y', y)
            call b%a%e3set(i, e3)
            call b%a%set(i, 'class', 1.0)
        end do
        ! simulate images
        do i=1,NIMGS
            b%img_copy = b%img
            o = b%a%get_ori(i)
            call simimg(b%img_copy, o, b%tfun, CTFFLAG, SNR, SNR_PINK, SNR_DETECTOR, BFAC)
            call b%img_copy%write('classrefine_tester_simimgs.mrc', i)
        end do
        ! zero the shifts (to provide a more realistic situation)
        do i=1,NIMGS
            call b%a%set(i, 'x', 0.0)
            call b%a%set(i, 'y', 0.0)
        end do
        call b%a%write('classrefine_tester_simoris.txt')
        call crefine%new(b, p, cline, 'classrefine_tester_simimgs.mrc', 1, p%msk, p%ctf)
    end subroutine setup_testenv

end module simple_classrefine_tester

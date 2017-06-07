module simple_wiener2D_tester
use simple_jiffys         ! singleton
use simple_filterer       ! singleton
use simple_params,        only: params
use simple_build,         only: build
use simple_image,         only: image
use simple_oris,          only: oris
use simple_cmdline,       only: cmdline
use simple_commander_sim, only: simimgs_commander
implicit none

public :: exec_wiener2D_test
private

! module global constants
integer, parameter           :: NPROJS=100
real, parameter              :: TRS=5.0, DEF=2.5, DFERR=0.3
real, parameter              :: SNR=0.2, BFAC=50., LP=8.
character(len=32), parameter :: ptclsname = 'wiener2D_test_ptcls.mrc'
character(len=32), parameter :: orisname  = 'wiener2D_test_oris.txt'

! module global variables
type(build)              :: b
type(params)             :: p
type(oris)               :: oset
type(image)              :: img_rec
type(image), allocatable :: imgs(:), img_ref(:)
type(simimgs_commander)  :: xsimimgs
type(cmdline)            :: cline_here
integer                  :: ldim(3)
logical                  :: verbose=.false.

contains

    subroutine exec_wiener2D_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        call setup_testenv( cline, be_verbose )
        write(*,*) '****wiener2D_test, init'
        call test_wiener_restore2D
        call shutdown_testenv
        write(*,*) '****wiener2D_test, completed'
    end subroutine exec_wiener2D_test

    subroutine setup_testenv( cline, be_verbose )
        use simple_projector_hlev, only: projvol_expanded
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        type(oris) :: o_single
        integer    :: iproj
        verbose = .false.
        if( present(be_verbose) ) verbose = be_verbose
        cline_here = cline
        ! it is assumed that vol1, smpd, msk are part of the inputted command line
        ! setting the rest of the command line up in here
        call cline_here%set('nptcls',    real(NPROJS))
        call cline_here%set('snr',       SNR         )
        call cline_here%set('ndiscrete', 1.0         )
        call cline_here%set('outstk',    ptclsname   )
        call cline_here%set('outfile',   orisname    )
        call cline_here%set('nspace',    real(NPROJS))
        call cline_here%set('sherr',     2.0         )
        call cline_here%set('ctf',       'yes'       )
        call cline_here%set('dferr',     DFERR       ) 
        call cline_here%set('bfac',      BFAC        )
        call cline_here%set('bfacerr',   0.          )
        call cline_here%set('wfun',      'kb'        )
        call cline_here%set('winsz',     1.5         )
        call cline_here%set('alpha',     2.0         )
        call cline_here%set('eo',        'no'        )
        call cline_here%set('trs',       'TRS'       )
        ! simulate images
        call xsimimgs%execute(cline_here)
        ! make params
        call cline_here%set('stk', ptclsname)
        p = params(cline_here)
        call b%build_general_tbox(p, cline_here) ! general objects built
        ! read back in images and orientations
        allocate(imgs(NPROJS))
        ldim = [p%box,p%box,1]
        do iproj=1,NPROJS
            call imgs(iproj)%new(ldim, p%smpd)
            call imgs(iproj)%read(ptclsname, iproj)
        end do
        call oset%new(NPROJS)
        call oset%read(orisname)
        ! create reconstructed image
        call img_rec%new([p%box,p%box,1], p%smpd)
        ! create reference
        call o_single%new(1) 
        call o_single%set_euler(1, [0.,0.,0.])
        call b%vol%read(p%vols(1))
        img_ref = projvol_expanded(b%vol, o_single,  p)
        call o_single%kill
    end subroutine setup_testenv

    subroutine test_wiener_restore2D
        real :: corr
        type(ctfplan) :: tfplan
        if( verbose ) write(*,*) 'testing simple_filterer :: wiener_restore2D'
        tfplan%mode = 'astig'
        tfplan%flag = 'yes'
        call wiener_restore2D(imgs, oset, tfplan, img_rec, p%msk)
        ! check if test passed
        call img_ref(1)%mask(p%msk, 'soft')
        call img_rec%mask(p%msk, 'soft')
        call img_rec%vis
        corr = img_ref(1)%corr(img_rec, LP)
        if( corr > 0.945 )then
            ! the test passed
        else
            write(*,*) 'ref/rec corr: ', corr
            stop '****wiener2D_tester FAILURE simple_filterer :: wiener_restore2D'
        endif 
    end subroutine test_wiener_restore2D

    subroutine shutdown_testenv
        integer :: iproj
        call b%kill_general_tbox
        call oset%kill
        call img_ref(1)%kill
        call img_rec%kill
        do iproj=1,NPROJS
            call imgs(iproj)%kill
        end do
        deallocate(imgs,img_ref)
    end subroutine shutdown_testenv

end module simple_wiener2D_tester
!@descr: sharpening integration tests with closed-form identical and opposite half-map FSCs
module simple_commanders_volops_tester
use simple_core_module_api
use simple_atoms,             only: atoms
use simple_cmdline,           only: cmdline
use simple_image,             only: image
use simple_exec_volume,       only: exec_volume_commander
use simple_test_utils,        only: assert_true, assert_real, assert_int
implicit none
private
public :: run_all_sharpvol_tests

contains

    subroutine run_all_sharpvol_tests()
        integer, parameter :: BOX=64
        real, parameter :: SMPD=2.
        type(atoms) :: model
        type(image) :: vol, result
        class(cmdline), allocatable :: cline
        logical :: did_execute, silent
        real, allocatable :: fsc(:)
        real, pointer :: rmat(:,:,:)
        real :: coords(3,4), positive_peak
        integer :: i
        character(len=*), parameter :: STEM='tmp_sharpvol'
        write(*,'(A)') '**** running sharpvol half-map FSC tests ****'
        coords = reshape([40.,50.,60., 76.,44.,58., 56.,82.,66., 60.,60.,88.], [3,4])
        call model%new(4, dummy=.true.)
        do i=1,4
            call model%set_name(i, ' C  ')
            call model%set_coord(i, coords(:,i))
            call model%set_resnum(i, i)
            call model%set_occupancy(i, 1.)
        enddo
        call model%guess_element()
        call model%writepdb(string(STEM//'.pdb'))
        call model%pdb2mrc(pdbfile=string(STEM//'.pdb'), volfile=string(STEM//'.mrc'), &
            &smpd=SMPD, vol_dim=[BOX,BOX,BOX])
        call model%kill()
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%read(string(STEM//'.mrc'))
        call vol%write(string(STEM//'_even.mrc'))
        call vol%write(string(STEM//'_odd.mrc'))
        allocate(cline)
        call cline%set('prg', 'sharpvol')
        call cline%set('vol1', STEM//'.mrc')
        call cline%set('pdbfile', STEM//'.pdb')
        call cline%set('smpd', SMPD)
        call cline%set('mkdir', 'no')
        call cline%set('mirr', 'no')
        call cline%set('nthr', 1)
        call cline%set('bfac', 0.)
        ! No FSC or low-pass input: the half-maps alone supply the filter.
        did_execute = .false.
        call exec_volume_commander('sharpvol', cline, silent, did_execute)
        call assert_true(did_execute, 'volume router executes sharpvol')
        fsc = file2rarr(string(STEM//'_fsc.bin'))
        call assert_int(BOX/2, size(fsc), 'FSC has one value per Fourier shell')
        call assert_real(0., maxval(abs(fsc-1.)), 1.e-5, 'identical halves have FSC one')
        call result%new([BOX,BOX,BOX], SMPD)
        call result%read(string(STEM//'_pproc.mrc'))
        call result%get_rmat_ptr(rmat)
        positive_peak = maxval(abs(rmat))
        call assert_true(positive_peak > 1.e-6, 'correlated halves retain model-masked density')
        ! Reusing the output path must replace the old curve, not read it as input.
        call cline%set('fsc', STEM//'_fsc.bin')
        call vol%mul(-1.)
        call vol%write(string(STEM//'_odd.mrc'))
        did_execute = .false.
        call exec_volume_commander('sharpvol', cline, silent, did_execute)
        call assert_true(did_execute, 'volume router executes sharpvol')
        fsc = file2rarr(string(STEM//'_fsc.bin'))
        call assert_real(0., maxval(abs(fsc+1.)), 1.e-5, 'opposite halves replace the FSC with minus one')
        call result%read(string(STEM//'_pproc.mrc'))
        call result%get_rmat_ptr(rmat)
        call assert_real(0., maxval(abs(rmat)), positive_peak*1.e-6, 'nonpositive FSC suppresses the sharpened map')
        call result%read(string(STEM//'_lp.mrc'))
        call result%get_rmat_ptr(rmat)
        call assert_real(0., maxval(abs(rmat)), positive_peak*1.e-6, 'nonpositive FSC also suppresses the low-pass map')
        call result%kill()
        call vol%kill()
        call cline%kill()
        deallocate(cline)
        call del_file(STEM//'.pdb')
        call del_file(STEM//'.mrc')
        call del_file(STEM//'_even.mrc')
        call del_file(STEM//'_odd.mrc')
        call del_file(STEM//'_pdb.mrc')
        call del_file(STEM//'_pdb_docked.mrc')
        call del_file(STEM//'_fsc.bin')
        call del_file(STEM//'_lp.mrc')
        call del_file(STEM//'_pproc.mrc')
    end subroutine run_all_sharpvol_tests

end module simple_commanders_volops_tester

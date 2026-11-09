!@descr: tracking-preparation project outputs, batch boundaries, raw means and rectangular scaling
module simple_tseries_prep4tracking_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,                only: CWD_GLOB
use simple_string,              only: string
use simple_cmdline,             only: cmdline
use simple_image,               only: image
use simple_sp_project,          only: sp_project
use single_commanders_tseries,  only: commander_tseries_prep4tracking
use simple_fileio,              only: file_exists, simple_getcwd
use simple_imghead,             only: find_ldim_nptcls, find_img_smpd
use simple_syslib,              only: simple_chdir
use simple_string_utils,        only: int2str_pad
use simple_test_utils,          only: assert_true, assert_false, assert_int, assert_real, &
    &enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_tseries_prep4tracking_tests

! Header arithmetic involves only a multiplication/division in single precision.
real, parameter :: HEADER_TOL = 16.*epsilon(1.)

contains

    subroutine run_all_tseries_prep4tracking_tests()
        use simple_ui, only: make_ui
        write(*,'(A)') '**** running all tseries prep4tracking tests ****'
        call make_ui
        call test_project_outputs(102, 24, 16, 1., 'prep4tracking_batch', .true.)
        call test_project_outputs(2, 64, 48, 0.5, 'prep4tracking_scale', .false.)
    end subroutine run_all_tseries_prep4tracking_tests

    !> Independently planted signal, defects and frame offsets determine the raw mean and tracking table.
    subroutine test_project_outputs( nframes, nx, ny, smpd, tag, defects )
        integer,          intent(in) :: nframes, nx, ny
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: tag
        logical,          intent(in) :: defects
        class(cmdline), allocatable :: cl
        class(sp_project), allocatable :: source, result
        type(commander_tseries_prep4tracking) :: prep
        type(image) :: img
        type(string) :: cwd, root, run_dir, source_file, fname, tracked, cwd_glob_saved
        real, allocatable :: pattern(:,:,:), actual(:,:,:), expected(:,:,:)
        real :: mean_tol
        integer :: nfail0, t, x, y, ntracked, unit, ios, idx, ox, oy, nout, ldim(3), nimgs, outdim(3)
        logical :: saw_hot, saw_dead
        write(*,'(A)') 'test_project_outputs: '//tag
        nfail0 = tests_failed
        cwd_glob_saved = CWD_GLOB
        call enter_fixture(tag, cwd, root)
        CWD_GLOB = root%to_char()
        allocate(cl, source, result, pattern(nx,ny,1), expected(nx,ny,1))
        source_file = root//'/tracking.simple'
        call source%os_mic%new(nframes, is_ptcl=.false.)
        do y = 1, ny
            do x = 1, nx
                pattern(x,y,1) = exp(-real((x-nx/2)**2+(y-ny/2)**2)/(8./smpd**2))
            enddo
        enddo
        if( defects )then
            pattern(1,1,1) = 1000.
            pattern(nx,ny,1) = -1000.
        endif
        call img%new([nx,ny,1], smpd, wthreads=.false.)
        do t = 1, nframes
            fname = root//'/frame_'//int2str_pad(t,3)//'.mrc'
            call img%set_rmat(pattern + 0.01*real(t), .false.)
            call img%write(fname)
            call source%os_mic%set(t, 'frame', fname)
            call source%os_mic%set(t, 'xdim', nx)
            call source%os_mic%set(t, 'ydim', ny)
            call source%os_mic%set(t, 'smpd', smpd)
            call source%os_mic%set(t, 'state', 1)
            call source%os_mic%set(t, 'kv', 300.)
        enddo
        call source%update_projinfo(source_file)
        call source%write(source_file)
        call source%kill
        call cl%set('prg', 'tseries_prep4tracking')
        call cl%set('projfile', source_file)
        call cl%set('nthr', 1)
        call cl%set('neigs', 4)
        call prep%execute(cl)
        call simple_getcwd(run_dir)
        call assert_true(run_dir /= root, 'tracking preparation creates a separate execution directory')
        call result%read(string('tracking.simple'))
        call assert_int(nframes, result%get_nframes(), 'tracking preparation retains the full frame table')
        outdim = [nint(real(nx)*min(1.,smpd)), nint(real(ny)*min(1.,smpd)), 1]
        ntracked = 0
        do t = 1, nframes
            if( mod(t-1,10) == 0 )then
                call assert_true(result%os_mic%isthere(t,'track_fname'), 'every group start receives a tracking record')
                if( result%os_mic%isthere(t,'track_fname') )then
                    ntracked = ntracked + 1
                    tracked = result%os_mic%get_str(t, 'track_fname')
                    call assert_true(file_exists(tracked), 'every tracking record names a written image')
                    if( file_exists(tracked) )then
                        call find_ldim_nptcls(tracked, ldim, nimgs)
                        call assert_true(all(ldim == outdim), 'tracking image header preserves the rectangular aspect ratio')
                        call assert_real(max(1.,smpd), find_img_smpd(tracked), HEADER_TOL, 'tracking header records scaled sampling')
                        call assert_real(real(nx)*smpd, real(ldim(1))*find_img_smpd(tracked), &
                            &HEADER_TOL*real(nx), 'tracking image preserves horizontal field of view')
                        call assert_real(real(ny)*smpd, real(ldim(2))*find_img_smpd(tracked), &
                            &HEADER_TOL*real(ny), 'tracking image preserves vertical field of view')
                        call img%new(outdim, max(1.,smpd), wthreads=.false.)
                        call img%read(tracked)
                        actual = img%get_rmat()
                        call assert_true(all(ieee_is_finite(actual)), 'tracking image contains finite pixels')
                        call assert_true(maxval(actual)-minval(actual) > 0.01, 'tracking image retains planted contrast')
                    endif
                    call assert_int(t, result%os_mic%get_int(t,'fromt'), 'tracking group begins at its original frame')
                    call assert_int(min(t+9,nframes), result%os_mic%get_int(t,'tot'), 'tracking group ends within available frames')
                endif
            else
                call assert_false(result%os_mic%isthere(t,'track_fname'), 'non-start frames receive no spurious tracking record')
            endif
        enddo
        call assert_int((nframes+9)/10, ntracked, 'both full batches and short tails are published')
        call assert_true(file_exists(string('global_avgimg.mrc')), 'the native raw global mean is written')
        if( file_exists(string('global_avgimg.mrc')) )then
            call img%new([nx,ny,1], smpd, wthreads=.false.)
            call img%read(string('global_avgimg.mrc'))
            actual = img%get_rmat()
            expected = pattern + 0.01*real(nframes+1)/2.
            ! Standard first-order SP summation bound, including the final division.
            mean_tol = 2.*real(nframes)*epsilon(1.)
            call assert_true(all(abs(actual-expected) <= mean_tol*max(1.,abs(expected))), &
                &'the raw global image is the independent arithmetic mean of every input frame')
        endif
        if( defects )then
            saw_hot = .false.
            saw_dead = .false.
            nout = 0
            open(newunit=unit, file='outliers.txt', status='old', action='read', iostat=ios)
            call assert_int(0, ios, 'the defect coordinate table is written')
            if( ios == 0 )then
                do
                    read(unit,*,iostat=ios) idx, ox, oy
                    if( ios /= 0 ) exit
                    nout = nout + 1
                    if( ox == 1 .and. oy == 1 ) saw_hot = .true.
                    if( ox == nx .and. oy == ny ) saw_dead = .true.
                enddo
                close(unit)
                call assert_int(2, nout, 'only the planted hot and dead pixels are reported')
                call assert_true(saw_hot .and. saw_dead, 'defect coordinates match the planted one-based pixels')
            endif
        endif
        call source%read(source_file)
        do t = 1, nframes
            call assert_false(source%os_mic%isthere(t,'track_fname'), 'the source project remains unchanged')
        enddo
        call img%new([nx,ny,1], smpd, wthreads=.false.)
        call img%read(root//'/frame_'//int2str_pad(nframes,3)//'.mrc')
        actual = img%get_rmat()
        call assert_true(all(actual == pattern+0.01*real(nframes)), 'the source frame pixels remain unchanged')
        call source%kill
        call result%kill
        call cl%kill
        call img%kill
        call simple_chdir(root)
        CWD_GLOB = root%to_char()
        call leave_fixture(cwd, root, nfail0)
        CWD_GLOB = cwd_glob_saved%to_char()
    end subroutine test_project_outputs

end module simple_tseries_prep4tracking_tester

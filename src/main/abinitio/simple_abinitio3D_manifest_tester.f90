!@descr: unit tests of the abinitio3D run manifest (simple_abinitio3D_manifest)
! Round trip of every record kind (the manifest read back is written again line
! for line; a path-valued input keeps its '/'), the refusals of another schema
! version, a truncated file, a missing or wrong checksum, an unknown field or
! input key, records after the end marker or the checksum and an incomplete run
! with a valid checksum; an input the format cannot hold leaves the manifest
! unpublished instead of stopping the run; the registered manifest of a
! project (a bare name resolves against the project file's own directory, an
! absolute one is kept, the registered run identifier must match); the replay
! of the base run's settings onto a command line; and the frozen-project
! validation (an eligible add-on output is a frozen input, so add-ons chain)
! with its negative cases: an ineligible or foreign manifest, a
! missing state map, a changed particle layout, stack table, optics/CTF
! parameters or final map, and another particle count. In-memory projects, one
! 8-pixel map and one small sigma2 stand-in file whose name has blanks.
module simple_abinitio3D_manifest_tester
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,                only: STDLEN
use simple_string,              only: string
use simple_fileio,              only: del_file, file_exists, fopen, fclose
use simple_image,               only: image
use simple_cmdline,             only: cmdline
use simple_sp_project,          only: sp_project
use simple_sigma2_state_file,   only: sigma2_state_digest_begin, sigma2_state_digest_text
use simple_abinitio3D_manifest, only: abinitio3D_manifest, abinitio3D_stage_record
use simple_test_utils
implicit none
private
public :: run_all_abinitio3D_manifest_tests

character(len=*), parameter :: MAN_FNAME   = 'tmp_abinitio3D_manifest_tester.txt'
character(len=*), parameter :: MAN_FNAME2  = 'tmp_abinitio3D_manifest_tester_2.txt'
character(len=*), parameter :: PROJ_FNAME  = 'tmp_abinitio3D_manifest_tester.simple'
character(len=*), parameter :: VOL_FNAME   = 'tmp_abinitio3D_manifest_tester_vol.mrc'
character(len=*), parameter :: SIGMA_FNAME = 'tmp abinitio3D manifest tester sigma2.bin'
character(len=*), parameter :: RUN_ID      = 'abinitio3D_20260926T120000.000_1'
integer,          parameter :: NPTCLS      = 12
integer,          parameter :: NLADDER     = 8

contains

    subroutine run_all_abinitio3D_manifest_tests()
        write(*,'(A)') '**** running all abinitio3D manifest tests ****'
        call test_round_trip()
        call test_file_refusals()
        call test_registration()
        call test_replay()
        call test_frozen_validation()
        call del_file(VOL_FNAME)
        call del_file(SIGMA_FNAME)
    end subroutine run_all_abinitio3D_manifest_tests

    ! ---- fixtures -------------------------------------------------------------

    !> NLADDER stages with distinct values in every field
    function ladder() result( stages )
        type(abinitio3D_stage_record) :: stages(NLADDER)
        integer :: i
        do i = 1, NLADDER
            stages(i)%lp_planned     = 20. - real(i)
            stages(i)%lp_emitted     = 20. - real(i) - 0.25
            stages(i)%lpstop_emitted = merge(0., 20. - real(i), i >= 6)
            stages(i)%box_crop       = 64 + 8*i
            stages(i)%smpd_crop      = 1.3 * 128. / real(64 + 8*i)
            stages(i)%scale          = real(64 + 8*i) / 128.
            stages(i)%trslim         = 0.1 * real(i)
            stages(i)%frc_crit       = 0.5
            stages(i)%l_autoscale    = i < NLADDER
            stages(i)%l_lpset        = mod(i,2) == 0
        enddo
    end function ladder

    !> a completed eligible abinitio3D run of nstates over spproj, built the way
    !! exec_abinitio3D builds it; its artifacts are those spproj registers
    subroutine make_manifest( man, spproj, nstates )
        type(abinitio3D_manifest), intent(inout) :: man
        type(sp_project),          intent(inout) :: spproj
        integer,                   intent(in)    :: nstates
        type(cmdline) :: cl
        call man%new(RUN_ID, 'abinitio3D', .true., spproj, 'raw')
        call man%set_solution(nstates, 'c3', 128, 1.3, 180., 'independent', 6)
        call man%set_sampling(10000, 8000, 1., .true.)
        call man%set_stage_line(.false., 0., .true., 6.)
        call man%set_ladder(1, 5, ladder())
        call cl%set('pgrp',        'c3')
        call cl%set('rec_backend', 'pcg')
        call cl%set('lpstop',      6.)
        call cl%set('vol1',        '/abs/refs/startvol_state01.mrc')
        call cl%set('nthr',        8)    ! not a manifest input key
        call man%record_inputs(cl)
        call man%record_artifacts(spproj)
        call cl%kill
    end subroutine make_manifest

    !> a particle project of NPTCLS rows in two stacks with CTF parameters, a
    !! registered state-1 map and a registered sigma2 state file
    subroutine make_project( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: i, istk
        call spproj%kill
        call spproj%projinfo%new(1, is_ptcl=.false.)
        call spproj%projinfo%set(1, 'projname', 'tester_project')
        call spproj%os_stk%new(2, is_ptcl=.false.)
        do istk = 1, 2
            call spproj%os_stk%set(istk, 'stk',   'stack_'//char(48+istk)//'.mrcs')
            call spproj%os_stk%set(istk, 'fromp', (istk-1)*NPTCLS/2 + 1)
            call spproj%os_stk%set(istk, 'top',   istk*NPTCLS/2)
            call spproj%os_stk%set(istk, 'box',   64)
            call spproj%os_stk%set(istk, 'smpd',  1.3)
            call spproj%os_stk%set(istk, 'ctf',   'yes')
            call spproj%os_stk%set(istk, 'kv',    300.)
            call spproj%os_stk%set(istk, 'cs',    2.7)
            call spproj%os_stk%set(istk, 'fraca', 0.1)
        enddo
        call spproj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        call spproj%os_ptcl2D%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            istk = merge(1, 2, i <= NPTCLS/2)
            call spproj%os_ptcl3D%set(i, 'stkind', istk)
            call spproj%os_ptcl3D%set(i, 'indstk', i - (istk-1)*NPTCLS/2)
            call spproj%os_ptcl3D%set(i, 'dfx',    1.0 + 0.1*real(i))
            call spproj%os_ptcl3D%set(i, 'dfy',    1.1 + 0.1*real(i))
            call spproj%os_ptcl3D%set(i, 'angast', 10.)
            call spproj%os_ptcl3D%set_state(i, 1)
        enddo
        spproj%os_ptcl2D = spproj%os_ptcl3D
        call write_volume(VOL_FNAME, 1.)
        call spproj%add_vol2os_out(string(VOL_FNAME), 1.3, 1, 'vol')
        call write_lines(SIGMA_FNAME, [character(len=1024) :: 'sigma2 stand-in'])
        call spproj%projinfo%set(1, 'sigma2_state', SIGMA_FNAME)
    end subroutine make_project

    subroutine write_volume( fname, val )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: val
        type(image) :: vol
        real :: rmat(8,8,8)
        rmat = val
        call vol%new([8,8,8], 1.3, wthreads=.false.)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(fname), del_if_exists=.true.)
        call vol%kill
    end subroutine write_volume

    !> a written manifest as lines, for tampering
    subroutine read_lines( fname, lines, n )
        character(len=*),    intent(in)  :: fname
        character(len=1024), intent(out) :: lines(:)
        integer,             intent(out) :: n
        integer :: funit, io_stat
        n = 0
        call fopen(funit, file=string(fname), status='OLD', action='READ', iostat=io_stat)
        if( io_stat /= 0 ) return
        do while( n < size(lines) )
            read(funit,'(A)',iostat=io_stat) lines(n+1)
            if( io_stat /= 0 ) exit
            n = n + 1
        enddo
        call fclose(funit)
    end subroutine read_lines

    subroutine write_lines( fname, lines )
        character(len=*),    intent(in) :: fname
        character(len=1024), intent(in) :: lines(:)
        integer :: funit, io_stat, i
        call fopen(funit, file=string(fname), status='REPLACE', action='WRITE', iostat=io_stat)
        do i = 1, size(lines)
            write(funit,'(A)') trim(lines(i))
        enddo
        call fclose(funit)
    end subroutine write_lines

    !> records followed by the checksum line the writer would give them
    subroutine write_with_checksum( fname, lines )
        character(len=*),    intent(in) :: fname
        character(len=1024), intent(in) :: lines(:)
        character(len=1024) :: str
        integer(int64)     :: checksum
        integer :: i
        checksum = sigma2_state_digest_begin()
        do i = 1, size(lines)
            call sigma2_state_digest_text(checksum, trim(lines(i)))
        enddo
        write(str,'(I0)') checksum
        str = 'checksum '//trim(str)
        call write_lines(fname, [lines, [str]])
    end subroutine write_with_checksum

    !> the text file holds the line
    logical function file_contains( fname, line ) result( l_found )
        character(len=*), intent(in) :: fname, line
        character(len=1024) :: lines(256)
        integer :: n, i
        call read_lines(fname, lines, n)
        l_found = .false.
        do i = 1, n
            if( trim(lines(i)) == trim(line) ) l_found = .true.
        enddo
    end function file_contains

    !> the two text files hold the same lines
    logical function same_text( fname1, fname2 ) result( l_same )
        character(len=*), intent(in) :: fname1, fname2
        character(len=1024) :: lines1(256), lines2(256)
        integer :: n1, n2
        call read_lines(fname1, lines1, n1)
        call read_lines(fname2, lines2, n2)
        l_same = n1 == n2 .and. n1 > 0
        if( l_same ) l_same = all(lines1(1:n1) == lines2(1:n2))
    end function same_text

    ! ---- tests ----------------------------------------------------------------

    subroutine test_round_trip()
        type(abinitio3D_manifest)     :: man, back
        type(abinitio3D_stage_record) :: a, b, stages(NLADDER)
        type(sp_project)      :: spproj
        type(string)          :: path
        character(len=STDLEN) :: msg
        logical :: found, l_same_stages
        integer :: status, i
        write(*,'(A)') 'test_round_trip'
        call make_project(spproj)
        call make_manifest(man, spproj, 2)
        call man%write(string(MAN_FNAME), status, msg)
        call assert_int(0, status, 'a manifest is written: '//trim(msg))
        call assert_false(file_exists(MAN_FNAME//'.tmp'), 'the atomic write leaves no temporary file')
        call back%read(string(MAN_FNAME), status, msg)
        call assert_int(0, status, 'a written manifest reads back: '//trim(msg))
        call assert_string_eq(MAN_FNAME, back%get_fname(), 'the manifest knows the file it was read from')
        ! identity, solution, sampling, stage line, ladder, inputs and
        ! artifacts: the manifest read back is written again record for record
        call back%write(string(MAN_FNAME2), status, msg)
        call assert_true(same_text(MAN_FNAME, MAN_FNAME2), 'every manifest record round trips exactly')
        call assert_char(RUN_ID, trim(back%get_run_id()), 'run identifier round trip')
        call assert_int(2,   back%get_nstates(), 'nstates round trip')
        call assert_int(128, back%get_box(),     'box round trip')
        call assert_real(1.3, back%get_smpd(), 0., 'smpd round trip (exact)')
        call assert_int(5, back%get_last_stage(), 'last stage round trip')
        call assert_int(NLADDER, back%get_nstages(), 'ladder length round trip')
        stages = ladder()
        l_same_stages = .true.
        do i = 1, NLADDER
            a = stages(i)
            b = back%get_stage(i)
            l_same_stages = l_same_stages .and. a%lp_planned == b%lp_planned .and. a%lp_emitted == b%lp_emitted .and. &
                &a%lpstop_emitted == b%lpstop_emitted .and. a%box_crop == b%box_crop .and. a%smpd_crop == b%smpd_crop &
                &.and. a%scale == b%scale .and. a%trslim == b%trslim .and. a%frc_crit == b%frc_crit .and. &
                &(a%l_autoscale .eqv. b%l_autoscale) .and. (a%l_lpset .eqv. b%l_lpset)
        enddo
        call assert_true(l_same_stages, 'every ladder field round trips exactly')
        call assert_char('raw', trim(back%get_ptcl_src()), 'particle source round trip')
        call assert_true(file_contains(MAN_FNAME2, 'input vol1 /abs/refs/startvol_state01.mrc'), &
            &'a path-valued input round trips with its slashes')
        ! artifacts: the registered map, and the sigma2 state, whose path has blanks
        call assert_true(back%matches_artifact('vol', 1, string(VOL_FNAME)), 'the state map is recorded with its digest')
        call assert_false(back%matches_artifact('vol', 2, string(VOL_FNAME)), 'no map is recorded for an unregistered state')
        call back%get_artifact('sigma2_state', 0, path, found)
        call assert_true(found, 'the registered sigma2 state is recorded')
        if( found ) call assert_true(back%matches_artifact('sigma2_state', 0, path), &
            &'an artifact path with blanks round trips to the file it names')
        call del_file(MAN_FNAME)
        call del_file(MAN_FNAME2)
        call man%kill
        call back%kill
        call spproj%kill
        call path%kill
    end subroutine test_round_trip

    subroutine test_file_refusals()
        type(abinitio3D_manifest)  :: man, back
        type(sp_project)           :: spproj
        type(cmdline)              :: cl
        character(len=1024)        :: lines(256), tampered(256)
        character(len=STDLEN)      :: msg
        integer :: status, n, i, iend
        write(*,'(A)') 'test_file_refusals'
        call make_project(spproj)
        call make_manifest(man, spproj, 2)
        call man%write(string(MAN_FNAME), status, msg)
        call read_lines(MAN_FNAME, lines, n)
        iend = 0
        do i = 1, n
            if( trim(lines(i)) == 'end' ) iend = i
        enddo
        call assert_int(n-1, iend, 'the end marker precedes the checksum line')
        ! the tester's checksum is the writer's
        call write_with_checksum(MAN_FNAME, lines(1:iend))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_int(0, status, 'records re-checksummed by the tester read back: '//trim(msg))
        ! another schema version
        tampered(1:n) = lines(1:n)
        tampered(1)   = 'abinitio3D_manifest 2'
        call write_lines(MAN_FNAME, tampered(1:n))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'another schema version is refused')
        ! not a manifest
        tampered(1) = 'abinitio3D_addon_frozen_context 1'
        call write_lines(MAN_FNAME, tampered(1:n))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'a file of another schema is refused')
        ! truncations
        call write_lines(MAN_FNAME, lines(1:n-1))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'a manifest without its checksum is refused')
        call write_lines(MAN_FNAME, lines(1:n/2))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'a manifest truncated in the middle is refused')
        ! checksum: one changed value
        tampered(1:n) = lines(1:n)
        do i = 1, n
            if( lines(i)(1:6) == 'nrows ' ) tampered(i) = 'nrows 13'
        enddo
        call write_lines(MAN_FNAME, tampered(1:n))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'a changed record fails the checksum')
        call assert_true(index(msg, 'checksum') > 0, 'the refusal names the checksum')
        ! unknown field
        tampered(1:n+1) = [lines(1:iend-1), [character(len=1024) :: 'frozen_rows 1 2 3'], lines(iend:n)]
        call write_lines(MAN_FNAME, tampered(1:n+1))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0 .and. index(msg, 'unknown') > 0, 'an unknown field is refused')
        ! an input key the manifest does not record, with a valid checksum
        tampered(1:iend) = [lines(1:iend-2), [character(len=1024) :: 'input nthr 8'], lines(iend-1:iend-1)]
        tampered(iend+1) = lines(iend)
        call write_with_checksum(MAN_FNAME, tampered(1:iend+1))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0 .and. index(msg, 'input key') > 0, 'an unknown input key is refused')
        ! records after the end marker
        tampered(1:n+1) = [lines(1:iend), [character(len=1024) :: 'nrows 12'], lines(iend+1:n)]
        call write_lines(MAN_FNAME, tampered(1:n+1))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'records after the end marker are refused')
        ! records after the checksum line
        call write_lines(MAN_FNAME, [lines(1:n), [character(len=1024) :: 'nrows 12']])
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0 .and. index(msg, 'after its checksum') > 0, 'records after the checksum are refused')
        ! an incomplete run, even with a valid checksum
        tampered(1:iend) = lines(1:iend)
        do i = 1, iend
            if( lines(i)(1:7) == 'status ' ) tampered(i) = 'status incomplete'
        enddo
        call write_with_checksum(MAN_FNAME, tampered(1:iend))
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0 .and. index(msg, 'completed run') > 0, &
            &'a manifest of an incomplete run is refused')
        ! a missing file
        call del_file(MAN_FNAME)
        call back%read(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0, 'a missing manifest is refused')
        ! an input value the format cannot hold: no manifest, no stop
        call make_manifest(man, spproj, 2)
        call cl%set('vol1', 'start vol.mrc')
        call man%record_inputs(cl)
        call man%write(string(MAN_FNAME), status, msg)
        call assert_true(status /= 0 .and. .not. file_exists(MAN_FNAME), &
            &'an input with blanks leaves the manifest unpublished with a status')
        call cl%kill
        call man%kill
        call back%kill
        call spproj%kill
    end subroutine test_file_refusals

    subroutine test_registration()
        type(abinitio3D_manifest) :: man, other, back
        type(sp_project)      :: spproj
        character(len=STDLEN) :: msg
        integer :: status
        write(*,'(A)') 'test_registration'
        call make_project(spproj)
        call make_manifest(man, spproj, 1)
        call man%write(string(MAN_FNAME), status, msg)
        call back%read_registered(spproj, string(PROJ_FNAME), status, msg)
        call assert_true(status /= 0, 'a project without a registration has no manifest')
        ! a bare name, next to the project file
        call man%register(spproj, MAN_FNAME)
        call back%read_registered(spproj, string(PROJ_FNAME), status, msg)
        call assert_int(0, status, 'a registered manifest is read: '//trim(msg))
        call assert_char(RUN_ID, trim(back%get_run_id()), 'the registered manifest is the run''s')
        call back%read_registered(spproj, string('/some/run/proj.simple'), status, msg)
        call assert_true(status /= 0 .and. index(msg, '/some/run/'//MAN_FNAME) > 0, &
            &'a bare name resolves against the project file''s own directory')
        ! an absolute registration
        call man%register(spproj, '/elsewhere/m.txt')
        call back%read_registered(spproj, string(PROJ_FNAME), status, msg)
        call assert_true(status /= 0 .and. index(msg, '/elsewhere/m.txt') > 0, 'an absolute registration is kept')
        ! the same file registered for another run
        call other%new('another_run', 'abinitio3D', .true., spproj, 'raw')
        call other%register(spproj, MAN_FNAME)
        call back%read_registered(spproj, string(PROJ_FNAME), status, msg)
        call assert_true(status /= 0, 'a manifest of another registered run identifier is refused')
        call del_file(MAN_FNAME)
        call man%kill
        call other%kill
        call back%kill
        call spproj%kill
    end subroutine test_registration

    subroutine test_replay()
        type(abinitio3D_manifest) :: man
        type(sp_project) :: spproj
        type(cmdline)    :: cl
        write(*,'(A)') 'test_replay'
        call make_project(spproj)
        call make_manifest(man, spproj, 2)
        call cl%set('lp', 12.)   ! a stale limit the replay drops
        call man%replay(cl)
        call assert_string_eq('pcg', cl%get_carg('rec_backend'), 'a replayed input keeps its value')
        call assert_string_eq('c3',  cl%get_carg('pgrp'),        'the point group of the solution')
        call assert_string_eq('c3',  cl%get_carg('pgrp_start'),  'the search starts in the point group of the solution')
        call assert_false(cl%defined('lp'), 'lp is replayed only when it was on the stage line')
        call assert_real(6., cl%get_rarg('lpstop'), 0., 'lpstop from the stage line shape')
        call assert_real(180., cl%get_rarg('mskdiam'), 0., 'the mask diameter of the solution')
        call assert_int(2, cl%get_iarg('nstates'), 'the state layout of the solution')
        call assert_string_eq('independent', cl%get_carg('multivol_mode'), 'a multi-state solution replays independently')
        call assert_int(5,     cl%get_iarg('nstages'), 'the ladder ends at the base run''s last stage')
        call assert_int(10000, cl%get_iarg('nsample'), 'the effective nsample')
        call assert_false(cl%defined('nthr'), 'a key outside the manifest inputs is not replayed')
        call assert_string_eq('raw', cl%get_carg('ptcl_src'), 'the particle source of the solution is replayed')
        call man%kill
        call spproj%kill
        call cl%kill
    end subroutine test_replay

    subroutine test_frozen_validation()
        type(sp_project)          :: spproj, altered
        type(abinitio3D_manifest) :: man, foreign
        character(len=STDLEN)     :: msg
        integer :: status
        write(*,'(A)') 'test_frozen_validation'
        call make_project(spproj)
        call make_manifest(man, spproj, 1)
        call man%validate_frozen(spproj, status, msg)
        call assert_int(0, status, 'the manifest validates against its own project: '//trim(msg))
        call make_foreign('abinitio3D', .false., 1)
        call foreign%validate_frozen(spproj, status, msg)
        call assert_true(status /= 0, 'an ineligible manifest is refused')
        call make_foreign('abinitio3D_addon', .true., 1)
        call foreign%validate_frozen(spproj, status, msg)
        call assert_int(0, status, 'an eligible add-on output is a frozen input (add-ons chain): '//trim(msg))
        call make_foreign('refine3D_auto', .true., 1)
        call foreign%validate_frozen(spproj, status, msg)
        call assert_true(status /= 0, 'a manifest of another program is refused')
        call make_foreign('abinitio3D', .true., 2)
        call foreign%validate_frozen(spproj, status, msg)
        call assert_true(status /= 0, 'a manifest without a map for every state is refused')
        altered = spproj
        call altered%os_ptcl3D%set(3, 'indstk', 4)
        call man%validate_frozen(altered, status, msg)
        call assert_true(status /= 0, 'a changed particle layout is refused')
        altered = spproj
        call altered%os_stk%set(2, 'box', 72)   ! outside the particle layout digest
        call man%validate_frozen(altered, status, msg)
        call assert_true(status /= 0 .and. index(msg, 'stack table') > 0, 'a changed stack table is refused')
        altered = spproj
        call altered%os_ptcl3D%set(5, 'dfx', 3.3)
        call man%validate_frozen(altered, status, msg)
        call assert_true(status /= 0, 'changed CTF parameters are refused')
        altered = spproj
        call altered%os_ptcl3D%reallocate(NPTCLS+1)
        call man%validate_frozen(altered, status, msg)
        call assert_true(status /= 0, 'another particle count is refused')
        call write_volume(VOL_FNAME, 2.)
        call man%validate_frozen(spproj, status, msg)
        call assert_true(status /= 0, 'a registered map that is not the base run''s map is refused')
        call spproj%kill
        call altered%kill
        call man%kill
        call foreign%kill

    contains

        subroutine make_foreign( program_name, l_eligible, nstates )
            character(len=*), intent(in) :: program_name
            logical,          intent(in) :: l_eligible
            integer,          intent(in) :: nstates
            call foreign%new(RUN_ID, program_name, l_eligible, spproj, 'raw')
            call foreign%set_solution(nstates, 'c3', 128, 1.3, 180., 'independent', 6)
            call foreign%record_artifacts(spproj)
        end subroutine make_foreign

    end subroutine test_frozen_validation

end module simple_abinitio3D_manifest_tester

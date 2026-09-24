!@descr: unit test routines for the project container as a whole (sp_project)
! The phase-shift policy across the record front doors (add_single_movie/add_stk/add_stktab,
! get_micparams/get_ctfparams, text and binary serialisation), the write/read round trip of the mic and
! particle segments, what the project file must not persist, partial alignment documents merged back
! (merge_algndocs) and the GUI's JSON view of a segment (print_segment_json: windows, sorting, histogram and
! plot blocks, absent optionals).
module simple_sp_project_tester
use simple_core_module_api   ! string, oris, ctfparams, del_file, filepath, syslib, PI, ran3, logfhandle
use simple_test_utils        ! assertions etc.
use simple_sp_project,       only: sp_project
use simple_binoris_io,       only: binwrite_oritab
use simple_image,            only: image
use json_kinds,              only: RK, IK
use json_module,             only: json_core, json_value
implicit none
private
public :: run_all_sp_project_tests

integer,          parameter :: NMICS  = 7
integer,          parameter :: NPTCLS = 300
real,             parameter :: TOL    = 1.0e-6
character(len=*), parameter :: PROJ_FILE  = 'sp_project_tester.simple'
character(len=*), parameter :: JSON_FILE  = 'sp_project_tester.json'
character(len=*), parameter :: IMG_FILE   = 'sp_project_tester_record.mrc'
character(len=*), parameter :: STK_FILE   = 'sp_project_tester_stack.mrc'
character(len=*), parameter :: MIC_TXT    = 'sp_project_tester_mic.txt'
character(len=*), parameter :: PHASE_PROJ = 'sp_project_tester_phase.simple'
character(len=*), parameter :: MERGE_PROJ = 'sp_project_tester_merge.simple'

contains

    subroutine run_all_sp_project_tests()
        write(*,'(A)') '**** running all project record tests ****'
        call set_fixed_seed(20260924)   ! the round-trip records are random draws; seed_rnd read /dev/urandom
        call test_phase_records()
        call test_phase_serialisation()
        call test_write_read_roundtrip()
        call test_read_does_not_mutate_projinfo()
        call test_partial_docs_merge()
        call test_print_segment_json()
        call cleanup()
    end subroutine run_all_sp_project_tests

    subroutine cleanup()
        call del_file(string(PROJ_FILE))
        call del_file(string(JSON_FILE))
        call del_file(string(IMG_FILE))
        call del_file(string(STK_FILE))
        call del_file(string(MIC_TXT))
        call del_file(string(PHASE_PROJ))
        call del_file(string(MERGE_PROJ))
    end subroutine cleanup

    !---------------- fixtures ----------------

    function zero_phase_ctf() result( ctfvars )
        type(ctfparams) :: ctfvars
        ctfvars%smpd     = 1.
        ctfvars%kv       = 300.
        ctfvars%cs       = 2.7
        ctfvars%fraca    = 0.1
        ctfvars%dfx      = 1.5
        ctfvars%dfy      = 1.6
        ctfvars%angast   = 5.
        ctfvars%ctfflag  = CTFFLAG_YES
        ctfvars%phshift  = 0.
    end function zero_phase_ctf

    ! writes an 8x8 image to fname (nimgs times) so the add_* front doors find a stack on disk
    subroutine write_dummy_stack( fname, nimgs )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: nimgs
        type(image) :: img
        integer     :: i
        call img%new([8,8,1], 1.)
        do i = 1,nimgs
            call img%write(string(fname), i)
        end do
        call img%kill
    end subroutine write_dummy_stack

    ! a project with NMICS micrographs (CTF + file names) and NPTCLS particles in ptcl2D and ptcl3D
    subroutine make_project( proj )
        type(sp_project), intent(inout) :: proj
        type(string) :: template
        real         :: eullims(3,2)
        integer      :: i
        eullims(1,:) = [0.,360.]
        eullims(2,:) = [0.,180.]
        eullims(3,:) = [0.,360.]
        call proj%os_mic%new(NMICS, is_ptcl=.false.)
        do i = 1,NMICS
            template = 'FoilHole_'//int2str_pad(i,16)
            call proj%os_mic%set(i, 'movie',    template%to_char()//'_fractions.tiff')
            call proj%os_mic%set(i, 'intg',     template%to_char()//'_intg.mrc')
            call proj%os_mic%set(i, 'forctf',   template%to_char()//'_forctf.mrc')
            call proj%os_mic%set(i, 'thumb',    template%to_char()//'_thumb.jpg')
            call proj%os_mic%set(i, 'smpd',     1.34)
            call proj%os_mic%set(i, 'kv',       300.)
            call proj%os_mic%set(i, 'cs',       2.7 )
            call proj%os_mic%set(i, 'fraca',    0.1 )
            call proj%os_mic%set(i, 'dfx',      1.+ran3())
            call proj%os_mic%set(i, 'dfy',      1.+ran3())
            call proj%os_mic%set(i, 'angast',   360.*ran3())
            call proj%os_mic%set(i, 'phshift',  0.)
        enddo
        call proj%os_ptcl2D%new(NPTCLS, .true.)
        do i = 1,NPTCLS
            call proj%os_ptcl2D%set_dfx(i, 1.+ran3())
            call proj%os_ptcl2D%set_dfy(i, 1.+ran3())
            call proj%os_ptcl2D%set(i, 'corr',    ran3())
            call proj%os_ptcl2D%set_class(i, nint(ran3()*100.))
            call proj%os_ptcl2D%set_state(i, nint(ran3()))
        enddo
        proj%os_ptcl3D = proj%os_ptcl2D
        call proj%os_ptcl2D%rnd_oris( trs=5.0, eullims=eullims)
        call proj%os_ptcl2D%delete_entry('z')
        call proj%os_ptcl3D%rnd_oris( trs=5.0, eullims=eullims)
        call proj%os_ptcl3D%delete_entry('class')
        call proj%os_ptcl3D%delete_entry('z')
        do i = 1,NPTCLS
            call proj%os_ptcl3D%set(i, 'proj', nint(1000.*ran3()))
        enddo
        call proj%update_projinfo(string(PROJ_FILE))
    end subroutine make_project

    subroutine assert_segments_equal( os1, os2, n, label )
        class(oris),      intent(inout) :: os1, os2
        integer,          intent(in)    :: n
        character(len=*), intent(in)    :: label
        type(string) :: str1, str2
        integer      :: i, nbad
        call assert_int(n, os2%get_noris(), label//': record count')
        if( os2%get_noris() /= n ) return
        nbad = 0
        do i = 1,n
            str1 = os1%ori2str(i)
            str2 = os2%ori2str(i)
            if( str1 /= str2 ) nbad = nbad + 1
            call str1%kill
            call str2%kill
        enddo
        call assert_int(0, nbad, label//': every record is identical (ori2str)')
    end subroutine assert_segments_equal

    subroutine assert_phase_at( os, i, label, expected )
        class(oris),      intent(in) :: os
        integer,          intent(in) :: i
        character(len=*), intent(in) :: label
        real,             intent(in) :: expected
        call assert_true(os%isthere(i, 'phshift'), label//': phshift present')
        if( os%isthere(i, 'phshift') ) call assert_real(expected, os%get(i, 'phshift'), TOL, label//': phshift value')
    end subroutine assert_phase_at

    !---------------- phase shift policy ----------------

    ! New CTF-bearing records always carry an explicit numerical phase. Zero is the conventional CTF
    ! value, not an absent or provenance-dependent value. A phase beyond pi must survive every mapping
    ! boundary unchanged: folding it into [0,pi) would negate the transfer function, which is what makes
    ! particles fitted either side of the pi wrap cancel each other in 2D class averages. Only a full turn
    ! is the identity, and it must reduce rather than accumulate.
    subroutine test_phase_records()
        type(sp_project) :: proj
        type(ctfparams)  :: ctfvars, mapped
        type(oris)       :: ptcl_oris
        type(string)     :: stknames(1)
        write(*,'(A)') 'test_phase_records'
        call write_dummy_stack(IMG_FILE, 1)
        ! zero phase through the movie and stack front doors
        ctfvars = zero_phase_ctf()
        call proj%add_single_movie(string(IMG_FILE), ctfvars)
        call proj%add_stk(string(IMG_FILE), ctfvars)
        call assert_phase_at(proj%os_mic,    1, 'zero phase os_mic',    0.)
        call assert_phase_at(proj%os_stk,    1, 'zero phase os_stk',    0.)
        call assert_phase_at(proj%os_ptcl2D, 1, 'zero phase os_ptcl2D', 0.)
        call assert_phase_at(proj%os_ptcl3D, 1, 'zero phase os_ptcl3D', 0.)
        call assert_int(1, proj%os_mic%get_noris(),    'add_single_movie adds one micrograph')
        call assert_int(1, proj%os_stk%get_noris(),    'add_stk adds one stack')
        call assert_int(1, proj%os_ptcl2D%get_noris(), 'add_stk adds the stack particles to ptcl2D')
        call assert_int(1, proj%os_ptcl3D%get_noris(), 'add_stk adds the stack particles to ptcl3D')
        call proj%kill
        ! per-particle phases remain authoritative when a stack contains mixed phases; prepending an
        ! existing stack also exercises local input versus global project particle indexing in add_stktab
        call write_dummy_stack(STK_FILE, 2)
        call ptcl_oris%new(2, is_ptcl=.true.)
        call ptcl_oris%set_dfx(1, ctfvars%dfx)
        call ptcl_oris%set_dfy(1, ctfvars%dfy)
        call ptcl_oris%set(1, 'phshift', PI/4.)
        call ptcl_oris%set_dfx(2, ctfvars%dfx)
        call ptcl_oris%set_dfy(2, ctfvars%dfy)
        call ptcl_oris%set(2, 'phshift', PIO2)
        stknames(1) = STK_FILE
        call proj%add_stk(string(IMG_FILE), ctfvars)
        call proj%add_stktab(stknames, ctfvars, ptcl_oris)
        call assert_int(2, proj%os_stk%get_noris(),    'add_stktab appends a second stack')
        call assert_int(3, proj%os_ptcl2D%get_noris(), 'add_stktab appends its two particles after the first stack')
        call assert_phase_at(proj%os_stk,    2, 'heterogeneous os_stk',           0.)
        call assert_phase_at(proj%os_ptcl2D, 2, 'heterogeneous os_ptcl2D first',  PI/4.)
        call assert_phase_at(proj%os_ptcl2D, 3, 'heterogeneous os_ptcl2D second', PIO2)
        call ptcl_oris%kill
        call proj%kill
        ! a phase beyond pi survives every mapping boundary
        mapped = ctfvars
        mapped%phshift = PI + PI/4.
        call proj%add_single_movie(string(IMG_FILE), mapped)
        call proj%add_stk(string(IMG_FILE), mapped)
        call assert_phase_at(proj%os_mic,    1, 'mapped os_mic',    PI + PI/4.)
        call assert_phase_at(proj%os_stk,    1, 'mapped os_stk',    PI + PI/4.)
        call assert_phase_at(proj%os_ptcl2D, 1, 'mapped os_ptcl2D', PI + PI/4.)
        call assert_phase_at(proj%os_ptcl3D, 1, 'mapped os_ptcl3D', PI + PI/4.)
        mapped = proj%get_micparams(1)
        call assert_real(PI + PI/4., mapped%phshift, TOL, 'os_mic -> ctfparams keeps the phase beyond pi')
        call assert_real(300., mapped%kv,    TOL, 'os_mic -> ctfparams carries the voltage')
        call assert_real(2.7,  mapped%cs,    TOL, 'os_mic -> ctfparams carries Cs')
        call assert_real(0.1,  mapped%fraca, TOL, 'os_mic -> ctfparams carries the amplitude contrast')
        call assert_real(0.,   mapped%dfx,   TOL, 'a movie record has no defocus yet (add_single_movie stores optics only)')
        call assert_real(0.,   mapped%angast, TOL, 'a movie record has no astigmatism angle yet')
        mapped = proj%get_ctfparams('stk', 1)
        call assert_real(PI + PI/4., mapped%phshift, TOL, 'os_stk -> ctfparams keeps the phase beyond pi')
        mapped = proj%get_ctfparams('ptcl2D', 1)
        call assert_real(PI + PI/4., mapped%phshift, TOL, 'os_ptcl2D -> ctfparams keeps the phase beyond pi')
        call assert_real(1.6, mapped%dfy, TOL, 'os_ptcl2D -> ctfparams carries dfy')
        mapped = proj%get_ctfparams('ptcl3D', 1)
        call assert_real(PI + PI/4., mapped%phshift, TOL, 'os_ptcl3D -> ctfparams keeps the phase beyond pi')
        call assert_true(mapped%ctfflag == CTFFLAG_YES, 'os_ptcl3D -> ctfparams carries the CTF flag')
        call proj%kill
        ! only a full turn is the identity, and it reduces rather than accumulates
        mapped%phshift = TWOPI + PI + PI/4.
        mapped%dfx     = ctfvars%dfx
        mapped%dfy     = ctfvars%dfy
        call proj%add_single_movie(string(IMG_FILE), mapped)
        call proj%add_stk(string(IMG_FILE), mapped)
        call assert_phase_at(proj%os_ptcl2D, 1, 'wrapped os_ptcl2D', PI + PI/4.)
        call assert_phase_at(proj%os_mic,    1, 'wrapped os_mic',    PI + PI/4.)
        call proj%kill
    end subroutine test_phase_records

    ! serialisation is the final schema boundary: even partially assembled project records materialise
    ! the identity phase rather than persisting an absent field, on the text and the binary route
    subroutine test_phase_serialisation()
        type(sp_project) :: proj, proj_read
        write(*,'(A)') 'test_phase_serialisation'
        call proj%os_mic%new(1, is_ptcl=.false.)
        call proj%os_stk%new(1, is_ptcl=.false.)
        call proj%os_ptcl2D%new(1, is_ptcl=.true.)
        call proj%os_ptcl3D%new(1, is_ptcl=.true.)
        call proj%os_mic%set(1, 'smpd', 1.)
        call proj%os_mic%write(string(MIC_TXT))
        call proj_read%os_mic%new(1, is_ptcl=.false.)
        call proj_read%read_segment('mic', string(MIC_TXT))
        call assert_phase_at(proj_read%os_mic, 1, 'text read os_mic', 0.)
        call proj_read%kill
        call proj%update_projinfo(string(PHASE_PROJ))
        call proj%write(string(PHASE_PROJ))
        call proj_read%read(string(PHASE_PROJ))
        call assert_phase_at(proj_read%os_mic,    1, 'serialised os_mic',    0.)
        call assert_phase_at(proj_read%os_stk,    1, 'serialised os_stk',    0.)
        call assert_phase_at(proj_read%os_ptcl2D, 1, 'serialised os_ptcl2D', 0.)
        call assert_phase_at(proj_read%os_ptcl3D, 1, 'serialised os_ptcl3D', 0.)
        call proj%kill
        call proj_read%kill
    end subroutine test_phase_serialisation

    !---------------- write / read ----------------

    subroutine test_write_read_roundtrip()
        type(sp_project) :: proj, proj_read
        write(*,'(A)') 'test_write_read_roundtrip'
        call make_project(proj)
        call del_file(string(PROJ_FILE))
        call proj%write(string(PROJ_FILE))
        call assert_true(file_exists(string(PROJ_FILE)), 'write creates the project file')
        call proj_read%read(string(PROJ_FILE))
        call assert_segments_equal(proj%os_mic,    proj_read%os_mic,    NMICS,  'os_mic')
        call assert_segments_equal(proj%os_ptcl2D, proj_read%os_ptcl2D, NPTCLS, 'os_ptcl2D')
        call assert_segments_equal(proj%os_ptcl3D, proj_read%os_ptcl3D, NPTCLS, 'os_ptcl3D')
        call assert_int(0, proj_read%os_stk%get_noris(), 'an empty segment stays empty')
        call assert_false(proj_read%compenv%isthere('simple_path'), 'the project does not persist compenv simple_path')
        call assert_true(proj_read%projinfo%isthere(1, 'projname'), 'projinfo carries the project name')
        call assert_string_eq('FoilHole_0000000000000003_intg.mrc', proj_read%os_mic%get_str(3, 'intg'), 'os_mic string values survive')
        call assert_real(proj%os_mic%get(5, 'dfx'), proj_read%os_mic%get(5, 'dfx'), TOL, 'os_mic real values survive')
        call proj%kill
        call proj_read%kill
    end subroutine test_write_read_roundtrip

    ! reading a project through an absolute path from another working directory must not rewrite
    ! the on-disk projinfo (cwd bookkeeping is the caller's business)
    subroutine test_read_does_not_mutate_projinfo()
        type(sp_project) :: proj, probe
        type(string)     :: cwd_orig, probe_dir, proj_abs, projinfo_before, projinfo_after
        integer          :: status
        write(*,'(A)') 'test_read_does_not_mutate_projinfo'
        call make_project(proj)
        call del_file(string(PROJ_FILE))
        call proj%write(string(PROJ_FILE))
        call probe%read_segment('projinfo', string(PROJ_FILE))
        projinfo_before = probe%projinfo%ori2str(1)
        call probe%kill
        call simple_getcwd(cwd_orig)
        proj_abs  = filepath(cwd_orig, string(PROJ_FILE))
        probe_dir = filepath(cwd_orig, string('sp_project_read_probe_'//int2str(get_process_id())))
        call simple_mkdir(probe_dir)
        call simple_chdir(probe_dir, status)
        call assert_int(0, status, 'chdir into the probe directory')
        if( status /= 0 ) return
        call probe%read(proj_abs)
        call assert_int(NMICS, probe%os_mic%get_noris(), 'the project reads through its absolute path')
        call probe%kill
        call simple_chdir(cwd_orig, status)
        call assert_int(0, status, 'chdir back out of the probe directory')
        if( status /= 0 ) return
        call simple_rmdir(probe_dir, status)
        call probe%read_segment('projinfo', string(PROJ_FILE))
        projinfo_after = probe%projinfo%ori2str(1)
        call assert_true(projinfo_before == projinfo_after, 'sp_project%read leaves the on-disk projinfo unchanged')
        call probe%kill
        call proj%kill
    end subroutine test_read_does_not_mutate_projinfo

    !---------------- partial alignment documents ----------------

    ! three partial ptcl2D documents written by binwrite_oritab (as the distributed workflow does) merged
    ! back into an empty project reproduce the source segment record for record
    subroutine test_partial_docs_merge()
        type(sp_project) :: proj, merged
        write(*,'(A)') 'test_partial_docs_merge'
        call make_project(proj)
        call merged%os_ptcl2D%new(NPTCLS, .true.)
        call merged%update_projinfo(string(MERGE_PROJ))
        call binwrite_oritab(string('doc_1.simple'), proj, proj%os_ptcl2D, [  1, 100],    isegment=PTCL2D_SEG)
        call binwrite_oritab(string('doc_2.simple'), proj, proj%os_ptcl2D, [101, 200],    isegment=PTCL2D_SEG)
        call binwrite_oritab(string('doc_3.simple'), proj, proj%os_ptcl2D, [201, NPTCLS], isegment=PTCL2D_SEG)
        call assert_true(file_exists(string('doc_2.simple')), 'partial documents are written')
        call merged%merge_algndocs(NPTCLS, 3, 'ptcl2D', 'doc_', 1)
        call del_file(string('doc_1.simple'))
        call del_file(string('doc_2.simple'))
        call del_file(string('doc_3.simple'))
        call assert_segments_equal(proj%os_ptcl2D, merged%os_ptcl2D, NPTCLS, 'merged os_ptcl2D')
        call assert_int(0, merged%os_ptcl3D%get_noris(), 'merging ptcl2D documents leaves ptcl3D alone')
        call proj%kill
        call merged%kill
    end subroutine test_partial_docs_merge

    !---------------- JSON view ----------------

    ! print_segment_json writes the GUI's view of a segment to logfhandle; divert it to a file and parse
    ! it back: the data window, the indices before/after it, sorting in both directions, and the
    ! histogram and plot blocks
    subroutine test_print_segment_json()
        type(sp_project)          :: proj
        type(json_core)           :: json
        type(json_value), pointer :: root, data, child, block
        real(RK)                  :: dval, prev
        real(RK),    allocatable  :: dvec(:)
        integer(IK), allocatable  :: ivec(:)
        real                      :: dfx_sorted(NMICS)
        integer                   :: order(NMICS), i, unit, logfhandle_saved, n
        logical                   :: found, l_sorted
        write(*,'(A)') 'test_print_segment_json'
        call make_project(proj)
        dfx_sorted = proj%os_mic%get_all('dfx')
        order      = [(i, i=1,NMICS)]
        call hpsort(dfx_sorted, order)
        ! plain window [2,5], no sorting
        call dump_json(proj, 'mic', fromto=[2,5], sort_key='n', sort_asc='yes', hist='no', plot_key='')
        call json%initialize()
        call json%parse(file=JSON_FILE, p=root)
        call assert_false(json%failed(), 'the JSON output parses')
        call json%get(root, 'data', data, found)
        call assert_true(found, 'a data array is present')
        if( found )then
            call assert_int(4, json%count(data), 'the window [2,5] holds four records')
            do i = 1,json%count(data)
                call json%get_child(data, i, child)
                call json%get(child, 'dfx', dval, found)
                call assert_true(found, 'record '//int2str(i)//' carries dfx')
                if( found ) call assert_real(proj%os_mic%get(i+1, 'dfx'), real(dval), 1.e-5, 'record '//int2str(i)//' is micrograph '//int2str(i+1))
            end do
            call json%get_child(data, 1, child)
            call json%get(child, 'intg', block, found)
            call assert_true(found, 'string keys are emitted')
        endif
        call json%get(root, 'indices_pre', ivec, found)
        call assert_true(found, 'indices_pre is present for a window starting after 1')
        if( found )then
            call assert_int(1, size(ivec), 'one index precedes the window')
            call assert_int(1, int(ivec(1)), 'index 1 precedes the window')
        endif
        call json%get(root, 'indices_post', ivec, found)
        call assert_true(found, 'indices_post is present for a window ending before the last record')
        if( found )then
            call assert_int(2, size(ivec), 'two indices follow the window')
            call assert_int(6, int(ivec(1)), 'index 6 follows the window')
            call assert_int(7, int(ivec(2)), 'index 7 follows the window')
        endif
        call json%get(root, 'histogram', block, found)
        call assert_false(found, 'no histogram block unless asked for')
        call json%destroy(root)
        ! sorted ascending on dfx with histogram and plot blocks, whole segment
        call dump_json(proj, 'mic', sort_key='dfx', sort_asc='yes', hist='yes', plot_key='dfy')
        call json%parse(file=JSON_FILE, p=root)
        call json%get(root, 'data', data, found)
        call assert_int(NMICS, json%count(data), 'no window: every record is emitted')
        l_sorted = .true.
        prev     = -huge(prev)
        do i = 1,json%count(data)
            call json%get_child(data, i, child)
            call json%get(child, 'dfx', dval, found)
            if( dval < prev ) l_sorted = .false.
            prev = dval
        end do
        call assert_true(l_sorted, 'records come out in ascending dfx order')
        call json%get_child(data, 1, child)
        call json%get(child, 'dfx', dval, found)
        call assert_real(dfx_sorted(1), real(dval), 1.e-5, 'the first record has the smallest dfx')
        call json%get(root, 'histogram', block, found)
        call assert_true(found, 'a histogram block is present')
        if( found )then
            call json%get(block, 'data', dvec, found)
            call assert_int(20, size(dvec), 'the histogram has 20 bins')
            call assert_real(real(NMICS), real(sum(dvec)), 1.e-5, 'the histogram counts every record')
            call json%get(block, 'labels', dvec, found)
            call assert_int(20, size(dvec), 'the histogram has 20 bin labels')
        endif
        call json%get(root, 'plot', block, found)
        call assert_true(found, 'a plot block is present')
        if( found )then
            call json%get(block, 'data', data, found)
            call assert_int(NMICS, json%count(data), 'the plot has one point per record')
            call json%get_child(data, 1, child)
            call json%get(child, 'x', dval, found)
            call assert_real(proj%os_mic%get(1, 'dfx'), real(dval), 1.e-5, 'plot x is the sort key of the record')
            call json%get(child, 'y', dval, found)
            call assert_real(proj%os_mic%get(1, 'dfy'), real(dval), 1.e-5, 'plot y is the plot key of the record')
        endif
        call json%destroy(root)
        ! descending on dfx with a window: positions 2..5 of the descending order; indices_pre and
        ! indices_post are the records above and below the window as displayed, in display order
        ! (the GUI's reading; they were the ascending head and tail, i.e. swapped, before 2026-09-25)
        call dump_json(proj, 'mic', fromto=[2,5], sort_key='dfx', sort_asc='no', hist='no', plot_key='')
        call json%parse(file=JSON_FILE, p=root)
        call assert_false(json%failed(), 'the descending-window JSON parses')
        call json%get(root, 'data', data, found)
        call assert_int(4, json%count(data), 'descending window [2,5] holds four records')
        l_sorted = .true.
        prev     = huge(prev)
        do i = 1,json%count(data)
            call json%get_child(data, i, child)
            call json%get(child, 'dfx', dval, found)
            if( dval > prev ) l_sorted = .false.
            prev = dval
        end do
        call assert_true(l_sorted, 'records come out in descending dfx order')
        call json%get_child(data, 1, child)
        call json%get(child, 'dfx', dval, found)
        call assert_real(dfx_sorted(NMICS-1), real(dval), 1.e-5, 'the window starts at the second largest dfx')
        n = 0
        call json%get(root, 'indices_pre', ivec, found)
        call assert_true(found, 'descending window: indices_pre present')
        if( found )then
            n = n + size(ivec)
            call assert_int(1, size(ivec), 'descending window: one record above the window')
            if( size(ivec) == 1 ) call assert_int(order(NMICS), int(ivec(1)), 'descending window: the largest dfx is above the window')
        endif
        call json%get(root, 'indices_post', ivec, found)
        call assert_true(found, 'descending window: indices_post present')
        if( found )then
            n = n + size(ivec)
            call assert_int(2, size(ivec), 'descending window: two records below the window')
            if( size(ivec) == 2 )then
                call assert_int(order(2), int(ivec(1)), 'descending window: the second smallest dfx comes first below the window')
                call assert_int(order(1), int(ivec(2)), 'descending window: the smallest dfx is last')
            endif
        endif
        call assert_int(NMICS - 4, n, 'the indices outside the window account for the other records')
        call json%destroy(root)
        ! histogram and plot asked for without a sort key: neither block, no dereference of the absent key
        ! (calculate_histogram/calculate_plot read sort_key unguarded before 2026-09-23)
        call dump_json(proj, 'mic', hist='yes', plot_key='dfy')
        call json%parse(file=JSON_FILE, p=root)
        call assert_false(json%failed(), 'JSON without a sort key parses')
        call json%get(root, 'data', data, found)
        call assert_int(NMICS, json%count(data), 'no sort key: every record is emitted in file order')
        call json%get(root, 'histogram', block, found)
        call assert_false(found, 'no sort key: no histogram block')
        call json%get(root, 'plot', block, found)
        call assert_false(found, 'no sort key: no plot block')
        call json%destroy(root)
        ! a particle segment without a window
        call dump_json(proj, 'ptcl3D', sort_key='n', sort_asc='yes', hist='no', plot_key='')
        call json%parse(file=JSON_FILE, p=root)
        call json%get(root, 'data', data, found)
        call assert_int(NPTCLS, json%count(data), 'ptcl3D: one entry per particle')
        call json%get_child(data, 7, child)
        call json%get(child, 'proj', dval, found)
        call assert_real(proj%os_ptcl3D%get(7, 'proj'), real(dval), 1.e-5, 'ptcl3D: particle values are emitted')
        call json%destroy(root)
        call proj%kill

        contains

            ! absent optionals are passed on as absent
            subroutine dump_json( proj, oritype, fromto, sort_key, sort_asc, hist, plot_key )
                type(sp_project),           intent(inout) :: proj
                character(len=*),           intent(in)    :: oritype
                integer,          optional, intent(in)    :: fromto(2)
                character(len=*), optional, intent(in)    :: sort_key, sort_asc, hist, plot_key
                logfhandle_saved = logfhandle
                open(newunit=unit, file=JSON_FILE, status='replace', action='write')
                logfhandle = unit
                call proj%print_segment_json(oritype, string(PROJ_FILE), fromto=fromto, sort_key=sort_key,&
                    &sort_asc=sort_asc, hist=hist, plot_key=plot_key)
                close(unit)
                logfhandle = logfhandle_saved
            end subroutine dump_json

    end subroutine test_print_segment_json

end module simple_sp_project_tester

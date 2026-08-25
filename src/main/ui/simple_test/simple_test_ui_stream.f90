!@descr: module defining the user interfaces for stream test programs in the simple_test_exec suite
module simple_test_ui_stream
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('stream', 'Stream', 130)
type(ui_program), target :: abinitio2D_stream
type(ui_program), target :: assign_optics
type(ui_program), target :: gen_pickrefs
type(ui_program), target :: master
type(ui_program), target :: pick_extract
type(ui_program), target :: preproc
type(ui_program), target :: sieve_cavgs

contains

    subroutine construct_test_stream_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_abinitio2D_stream(tsttab)
        call new_assign_optics(tsttab)
        call new_gen_pickrefs(tsttab)
        call new_master(tsttab)
        call new_pick_extract(tsttab)
        call new_preproc(tsttab)
        call new_sieve_cavgs(tsttab)
    end subroutine construct_test_stream_programs

    subroutine new_abinitio2D_stream( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call abinitio2D_stream%new('abinitio2D_stream', 'abinitio2D_stream', &
            &'is a test program for abinitio2D_stream', 'simple_test_exec', .false.)
        call add_ui_program('abinitio2D_stream', abinitio2D_stream, tsttab, UI_CATEGORY)
    end subroutine new_abinitio2D_stream

    subroutine new_assign_optics( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call assign_optics%new('assign_optics', 'assign_optics', &
            &'is a test program for assign_optics', 'simple_test_exec', .false.)
        call add_ui_program('assign_optics', assign_optics, tsttab, UI_CATEGORY)
    end subroutine new_assign_optics

    subroutine new_gen_pickrefs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call gen_pickrefs%new('gen_pickrefs', 'gen_pickrefs', &
            &'is a test program for gen_pickrefs', 'simple_test_exec', .false.)
        call add_ui_program('gen_pickrefs', gen_pickrefs, tsttab, UI_CATEGORY)
    end subroutine new_gen_pickrefs

    subroutine new_master( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call master%new('master', 'master', 'is a test program for master', 'simple_test_exec', .false.)
        call add_ui_program('master', master, tsttab, UI_CATEGORY)
    end subroutine new_master

    subroutine new_pick_extract( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call pick_extract%new('pick_extract', 'pick_extract', &
            &'is a test program for pick_extract', 'simple_test_exec', .false.)
        call add_ui_program('pick_extract', pick_extract, tsttab, UI_CATEGORY)
    end subroutine new_pick_extract

    subroutine new_preproc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call preproc%new('preproc', 'preproc', 'is a test program for preproc', 'simple_test_exec', .false.)
        call add_ui_program('preproc', preproc, tsttab, UI_CATEGORY)
    end subroutine new_preproc

    subroutine new_sieve_cavgs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call sieve_cavgs%new('sieve_cavgs', 'sieve_cavgs', &
            &'is a test program for sieve_cavgs', 'simple_test_exec', .false.)
        call add_ui_program('sieve_cavgs', sieve_cavgs, tsttab, UI_CATEGORY)
    end subroutine new_sieve_cavgs

end module simple_test_ui_stream

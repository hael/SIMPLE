!@descr: module defining the user interfaces for stats  programs in the simple_test_exec suite
module simple_test_ui_stats
use simple_ui_modules
implicit none

type(ui_program), target :: class_sample_test
type(ui_program), target :: clustering
type(ui_program), target :: ctf_test
type(ui_program), target :: eo_diff
type(ui_program), target :: extr_frac
type(ui_program), target :: multinomal_test
type(ui_program), target :: pca_all
type(ui_program), target :: pca_imgvar
type(ui_program), target :: sp_project

contains

    subroutine construct_stats_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_class_sample_test(prgtab)
        call new_clustering(prgtab)
        call new_ctf_test(prgtab)
        call new_eo_diff(prgtab)
        call new_extr_frac(prgtab)
        call new_multinomal_test(prgtab)
        call new_pca_all(prgtab)
        call new_pca_imgvar(prgtab)
        call new_sp_project(prgtab)
    end subroutine construct_stats_programs

    subroutine print_stats_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('STADISTICS:', C_UNDERLINED)
        write(logfhandle,'(A)') class_sample_test%name%to_char()
        write(logfhandle,'(A)') clustering%name%to_char()
        write(logfhandle,'(A)') ctf_test%name%to_char()
        write(logfhandle,'(A)') eo_diff%name%to_char()
        write(logfhandle,'(A)') extr_frac%name%to_char()
        write(logfhandle,'(A)') multinomal_test%name%to_char()
        write(logfhandle,'(A)') pca_all%name%to_char()
        write(logfhandle,'(A)') pca_imgvar%name%to_char()
        write(logfhandle,'(A)') sp_project%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_stats_programs

    subroutine new_class_sample_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call class_sample_test%new(&
        &'class_sample_test',&                 ! name
        &'class_sample_test ',&                ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call class_sample_test%add_input(UI_IO, )
        ! parameter input/output
        !call class_sample_test%add_input(UI_IMG, )
        ! alternative inputs
        !call class_sample_test%add_input(UI_PARM, )
        ! search controls
        !call class_sample_test%add_input(UI_SRCH, )
        ! filter controls
        !call class_sample_test%add_input(UI_FILT, )
        ! mask controls
        !call class_sample_test%add_input(UI_MASK, )
        ! computer controls
        !call class_sample_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('class_sample_test', class_sample_test, prgtab)
    end subroutine new_class_sample_test

    subroutine new_clustering( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call clustering%new(&
        &'clustering',&                        ! name
        &'clustering ',&                       ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call clustering%add_input(UI_IO, )
        ! parameter input/output
        !call clustering%add_input(UI_IMG, )
        ! alternative inputs
        !call clustering%add_input(UI_PARM, )
        ! search controls
        !call clustering%add_input(UI_SRCH, )
        ! filter controls
        !call clustering%add_input(UI_FILT, )
        ! mask controls
        !call clustering%add_input(UI_MASK, )
        ! computer controls
        !call clustering%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('clustering', clustering, prgtab)
    end subroutine new_clustering

    subroutine new_ctf_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ctf_test%new(&
        &'ctf_test',&                          ! name
        &'ctf_test ',&                         ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call ctf_test%add_input(UI_IO, )
        ! parameter input/output
        !call ctf_test%add_input(UI_IMG, )
        ! alternative inputs
        !call ctf_test%add_input(UI_PARM, )
        ! search controls
        !call ctf_test%add_input(UI_SRCH, )
        ! filter controls
        !call ctf_test%add_input(UI_FILT, )
        ! mask controls
        !call ctf_test%add_input(UI_MASK, )
        ! computer controls
        !call ctf_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('ctf_test', ctf_test, prgtab)
    end subroutine new_ctf_test

    subroutine new_eo_diff( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call eo_diff%new(&
        &'eo_diff',&                           ! name
        &'eo_diff ',&                          ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call eo_diff%add_input(UI_IO, )
        ! parameter input/output
        !call eo_diff%add_input(UI_IMG, )
        ! alternative inputs
        !call eo_diff%add_input(UI_PARM, )
        ! search controls
        !call eo_diff%add_input(UI_SRCH, )
        ! filter controls
        !call eo_diff%add_input(UI_FILT, )
        ! mask controls
        !call eo_diff%add_input(UI_MASK, )
        ! computer controls
        !call eo_diff%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('eo_diff', eo_diff, prgtab)
    end subroutine new_eo_diff

    subroutine new_extr_frac( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call extr_frac%new(&
        &'extr_frac',&                         ! name
        &'extr_frac ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call extr_frac%add_input(UI_IO, )
        ! parameter input/output
        !call extr_frac%add_input(UI_IMG, )
        ! alternative inputs
        !call extr_frac%add_input(UI_PARM, )
        ! search controls
        !call extr_frac%add_input(UI_SRCH, )
        ! filter controls
        !call extr_frac%add_input(UI_FILT, )
        ! mask controls
        !call extr_frac%add_input(UI_MASK, )
        ! computer controls
        !call extr_frac%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('extr_frac', extr_frac, prgtab)
    end subroutine new_extr_frac

    subroutine new_multinomal_test( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call multinomal_test%new(&
        &'multinomal_test',&                   ! name
        &'multinomal_test ',&                  ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call multinomal_test%add_input(UI_IO, )
        ! parameter input/output
        !call multinomal_test%add_input(UI_IMG, )
        ! alternative inputs
        !call multinomal_test%add_input(UI_PARM, )
        ! search controls
        !call multinomal_test%add_input(UI_SRCH, )
        ! filter controls
        !call multinomal_test%add_input(UI_FILT, )
        ! mask controls
        !call multinomal_test%add_input(UI_MASK, )
        ! computer controls
        !call multinomal_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('multinomal_test', multinomal_test, prgtab)
    end subroutine new_multinomal_test

    subroutine new_pca_all( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pca_all%new(&
        &'pca_all',&                           ! name
        &'pca_all ',&                          ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call pca_all%add_input(UI_IO, )
        ! parameter input/output
        !call pca_all%add_input(UI_IMG, )
        ! alternative inputs
        !call pca_all%add_input(UI_PARM, )
        ! search controls
        !call pca_all%add_input(UI_SRCH, )
        ! filter controls
        !call pca_all%add_input(UI_FILT, )
        ! mask controls
        !call pca_all%add_input(UI_MASK, )
        ! computer controls
        !call pca_all%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('pca_all', pca_all, prgtab)
    end subroutine new_pca_all

    subroutine new_pca_imgvar( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pca_imgvar%new(&
        &'pca_imgvar',&                        ! name
        &'pca_imgvar ',&                       ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call pca_imgvar%add_input(UI_IO, )
        ! parameter input/output
        !call pca_imgvar%add_input(UI_IMG, )
        ! alternative inputs
        !call pca_imgvar%add_input(UI_PARM, )
        ! search controls
        !call pca_imgvar%add_input(UI_SRCH, )
        ! filter controls
        !call pca_imgvar%add_input(UI_FILT, )
        ! mask controls
        !call pca_imgvar%add_input(UI_MASK, )
        ! computer controls
        !call pca_imgvar%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('pca_imgvar', pca_imgvar, prgtab)
    end subroutine new_pca_imgvar

    subroutine new_sp_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call sp_project%new(&
        &'sp_project',&                        ! name
        &'sp_project ',&                       ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call sp_project%add_input(UI_IO, )
        ! parameter input/output
        !call sp_project%add_input(UI_IMG, )
        ! alternative inputs
        !call sp_project%add_input(UI_PARM, )
        ! search controls
        !call sp_project%add_input(UI_SRCH, )
        ! filter controls
        !call sp_project%add_input(UI_FILT, )
        ! mask controls
        !call sp_project%add_input(UI_MASK, )
        ! computer controls
        !call sp_project%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('sp_project', sp_project, prgtab)
    end subroutine new_sp_project

end module simple_test_ui_stats

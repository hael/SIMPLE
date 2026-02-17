!@ descr: module defining the user interfaces for stats programs in the simple_test_exec suite
module simple_test_ui_stats
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_clustering
type(ui_program), target :: simple_test_pca_all
type(ui_program), target :: simple_test_pca_imgvar
type(ui_program), target :: simple_test_class_sample
type(ui_program), target :: simple_test_multinomal
type(ui_program), target :: simple_test_extr_frac
type(ui_program), target :: simple_test_eo_diff
type(ui_program), target :: simple_test_ctf
type(ui_program), target :: simple_test_sp_project

contains

    subroutine construct_stats_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_clustering(prgtab)
        call new_simple_test_pca_all(prgtab)
        call new_simple_test_pca_imgvar(prgtab)
        call new_simple_test_class_sample(prgtab)
        call new_simple_test_multinomal(prgtab)
        call new_simple_test_extr_frac(prgtab)
        call new_simple_test_eo_diff(prgtab)
        call new_simple_test_ctf(prgtab)
        call new_simple_test_sp_project(prgtab)
    end subroutine construct_stats_programs

    subroutine print_stats_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('STADISTICS:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_clustering%name%to_char()
        write(logfhandle,'(A)') simple_test_pca_all%name%to_char()
        write(logfhandle,'(A)') simple_test_pca_imgvar%name%to_char()
        write(logfhandle,'(A)') simple_test_class_sample%name%to_char()
        write(logfhandle,'(A)') simple_test_multinomal%name%to_char()
        write(logfhandle,'(A)') simple_test_extr_frac%name%to_char()
        write(logfhandle,'(A)') simple_test_eo_diff%name%to_char()
        write(logfhandle,'(A)') simple_test_ctf%name%to_char()
        write(logfhandle,'(A)') simple_test_sp_project%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_stats_programs

    subroutine new_simple_test_clustering( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_clustering', simple_test_clustering, prgtab)
    end subroutine new_simple_test_clustering

    subroutine new_simple_test_pca_all( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_pca_all', simple_test_pca_all, prgtab)
    end subroutine new_simple_test_pca_all

    subroutine new_simple_test_pca_imgvar( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_pca_imgvar', simple_test_pca_imgvar, prgtab)
    end subroutine new_simple_test_pca_imgvar

    subroutine new_simple_test_class_sample( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_class_sample', simple_test_class_sample, prgtab)
    end subroutine new_simple_test_class_sample

    subroutine new_simple_test_multinomal( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_multinomal', simple_test_multinomal, prgtab)
    end subroutine new_simple_test_multinomal

    subroutine new_simple_test_extr_frac( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_extr_frac', simple_test_extr_frac, prgtab)
    end subroutine new_simple_test_extr_frac

    subroutine new_simple_test_eo_diff( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_eo_diff', simple_test_eo_diff, prgtab)
    end subroutine new_simple_test_eo_diff

    subroutine new_simple_test_ctf( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_ctf', simple_test_ctf, prgtab)
    end subroutine new_simple_test_ctf

    subroutine new_simple_test_sp_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_sp_project', simple_test_sp_project, prgtab)
    end subroutine new_simple_test_sp_project

end module simple_test_ui_stats

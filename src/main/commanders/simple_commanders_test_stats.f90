!@descr: for all stats tests
module simple_commanders_test_stats
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_class_sample_test
  contains
    procedure :: execute      => exec_test_class_sample_test
end type commander_test_class_sample_test

type, extends(commander_base) :: commander_test_clustering
  contains
    procedure :: execute      => exec_test_clustering
end type commander_test_clustering

type, extends(commander_base) :: commander_test_ctf_test
  contains
    procedure :: execute      => exec_test_ctf_test
end type commander_test_ctf_test

type, extends(commander_base) :: commander_test_eo_diff
  contains
    procedure :: execute      => exec_test_eo_diff
end type commander_test_eo_diff

type, extends(commander_base) :: commander_test_extr_frac
  contains
    procedure :: execute      => exec_test_extr_frac
end type commander_test_extr_frac

type, extends(commander_base) :: commander_test_multinomal_test
  contains
    procedure :: execute      => exec_test_multinomal_test
end type commander_test_multinomal_test

type, extends(commander_base) :: commander_test_pca_all
  contains
    procedure :: execute      => exec_test_pca_all
end type commander_test_pca_all

type, extends(commander_base) :: commander_test_pca_imgvar
  contains
    procedure :: execute      => exec_test_pca_imgvar
end type commander_test_pca_imgvar

type, extends(commander_base) :: commander_test_sp_project
  contains
    procedure :: execute      => exec_test_sp_project
end type commander_test_sp_project

contains

subroutine exec_test_class_sample_test( self, cline )
    class(commander_test_class_sample_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_CLASS_SAMPLE_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_class_sample_test

subroutine exec_test_clustering( self, cline )
    class(commander_test_clustering),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_CLUSTERING_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_clustering

subroutine exec_test_ctf_test( self, cline )
    class(commander_test_ctf_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_CTF_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ctf_test

subroutine exec_test_eo_diff( self, cline )
    class(commander_test_eo_diff),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_EO_DIFF_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_eo_diff

subroutine exec_test_extr_frac( self, cline )
    class(commander_test_extr_frac),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_EXTR_FRAC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_extr_frac

subroutine exec_test_multinomal_test( self, cline )
    class(commander_test_multinomal_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MULTINOMAL_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_multinomal_test

subroutine exec_test_pca_all( self, cline )
    class(commander_test_pca_all),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_PCA_ALL_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_pca_all

subroutine exec_test_pca_imgvar( self, cline )
    class(commander_test_pca_imgvar),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_PCA_IMGVAR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_pca_imgvar

subroutine exec_test_sp_project( self, cline )
    class(commander_test_sp_project),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SP_PROJECT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_sp_project

end module simple_commanders_test_stats

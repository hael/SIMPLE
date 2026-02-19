!@descr: execution of test stats processing commanders
module simple_test_exec_stats
use simple_cmdline,               only: cmdline
use simple_commanders_test_stats, only: commander_test_class_sample_test, commander_test_clustering, &
                                        commander_test_ctf_test, commander_test_eo_diff, &
                                        commander_test_extr_frac, commander_test_multinomal_test, &
                                        commander_test_pca_all, commander_test_pca_imgvar, &
                                        commander_test_sp_project
implicit none

public :: exec_test_stats_commander
private

type(commander_test_class_sample_test) :: xclass_sample_test
type(commander_test_clustering)        :: xclustering
type(commander_test_ctf_test)          :: xctf_test
type(commander_test_eo_diff)           :: xeo_diff
type(commander_test_extr_frac)         :: xextr_frac
type(commander_test_multinomal_test)   :: xmultinomal_test
type(commander_test_pca_all)           :: xpca_all
type(commander_test_pca_imgvar)        :: xpca_imgvar
type(commander_test_sp_project)        :: xsp_project

contains

    subroutine exec_test_stats_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'class_sample_test' )
                call xclass_sample_test%execute(cline)
            case( 'clustering' )
                call xclustering%execute(cline)
            case( 'ctf_test' )
                call xctf_test%execute(cline)
            case( 'eo_diff' )
                call xeo_diff%execute(cline)
            case( 'extr_frac' )
                call xextr_frac%execute(cline)
            case( 'multinomal_test' )
                call xmultinomal_test%execute(cline)
            case( 'pca_all' )
                call xpca_all%execute(cline)
            case( 'pca_imgvar' )
                call xpca_imgvar%execute(cline)
            case( 'sp_project' )
                call xsp_project%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_stats_commander

end module simple_test_exec_stats   

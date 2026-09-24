!@descr: unit tests for the promoted pairwise-neighbour feature of class-average quality (simple_cavg_quality_relations)
! Four classes with symmetric correlation and distance matrices; with k = 2 the feature of the
! first two classes is the mean pairwise distance to their two best-correlated neighbours
! (0.15 and 0.20 by hand), and the feature schema carries a name.
module simple_cavg_quality_relations_tester
use simple_cavg_quality_relations, only: calculate_promoted_feature
use simple_cavg_quality_types,     only: CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1
use simple_test_utils
implicit none
private
public :: run_all_cavg_quality_relations_tests

contains

    subroutine run_all_cavg_quality_relations_tests()
        write(*,'(A)') '**** running all class-average quality relation tests ****'
        call test_promoted_feature()
    end subroutine run_all_cavg_quality_relations_tests

    subroutine test_promoted_feature()
        real, parameter :: TOL = 1.0e-6
        real    :: cc(4,4), distance(4,4), raw(4)
        integer :: class_inds(4)
        write(*,'(A)') 'test_promoted_feature'
        cc = 0.0
        distance = 0.0
        raw = 0.0
        class_inds = [1, 2, 3, 4]
        cc(1,2) = 0.9
        cc(1,3) = 0.7
        cc(1,4) = 0.2
        cc(2,3) = 0.8
        cc(2,4) = 0.1
        cc(3,4) = 0.4
        cc = cc + transpose(cc)
        distance(1,2) = 0.10
        distance(1,3) = 0.20
        distance(1,4) = 0.90
        distance(2,3) = 0.30
        distance(2,4) = 0.80
        distance(3,4) = 0.70
        distance = distance + transpose(distance)
        call calculate_promoted_feature(cc, distance, class_inds, 2, raw)
        call assert_real(0.15, raw(1), TOL, 'relational feature: CC-anchor top-k mean')
        call assert_real(0.20, raw(2), TOL, 'relational feature: per-class neighbour ordering')
        call assert_true(len_trim(CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1) > 0, 'relational feature schema is named')
    end subroutine test_promoted_feature

end module simple_cavg_quality_relations_tester

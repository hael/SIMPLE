program simple_test_cavg_quality_relations
use simple_cavg_quality_relations, only: test_cavg_quality_relations
use simple_test_utils,             only: report_summary, tests_failed
implicit none
call test_cavg_quality_relations
call report_summary
if( tests_failed > 0 ) error stop 1
end program simple_test_cavg_quality_relations

!@descr: smoke test that invokes stream p03 initial analysis commander
program simple_test_stream_initial_analysis
use simple_core_module_api
use simple_cmdline, only: cmdline
use simple_stream_p03_initial_analysis, only: stream_p03_initial_analysis
implicit none

type(stream_p03_initial_analysis) :: commander
type(cmdline)                     :: cline

! Create the minimum watched folder layout expected by the stream stage.
call simple_mkdir('test_initial_analysis')
call simple_chdir('test_initial_analysis')

write(*,*) ">>> CREATING TEST FOLDERS"
call cline%set('prg',        'gen_pickrefs')
call cline%set('dir_target', '../preprocessing')
call cline%set('mkdir',      'yes')

call commander%execute(cline)

end program simple_test_stream_initial_analysis

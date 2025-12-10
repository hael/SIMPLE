program simple_test_ori
include 'simple_lib.f08'
use simple_ori
use simple_sym
use json_kinds
use json_module
use simple_chash
use simple_sp_project
use simple_cmdline
use simple_commanders_project
implicit none
#include "simple_local_flags.inc"
type(json_value), pointer   :: json_ori
type(commander_new_project) :: xnew_project
type(ori)                   :: o_truth, o, o1, compenv_o
type(cmdline)               :: cline
type(sp_project)            :: spproj
type(sym)                   :: pgrpsyms
type(chash)                 :: qdescr
real                        :: vec(3), angle, rotmat(3, 3), euls(3), shvec(2)
type(string)                :: key, projname, projfile, str_projname, str_projrec, str_ori
logical                     :: test_passed
type(string), allocatable   :: keys(:)
test_passed = .true.
projname    = 'str_proj'
call pgrpsyms%new('c1')
call o_truth%new(.true.)
call pgrpsyms%rnd_euler(o_truth)
print *, '>>> MATRIX ', o_truth%get_mat()
call o_truth%get_axis_angle( vec, angle )
print *, '>>> VEC ', vec
print *, '>>> ANGLE ',angle/pi*180.
key         = ""
rotmat(:,:) = 1.
euls        = [1.,2.,3.]
shvec(:)    = 1.
call o%new(.true.)
call o1%new(.true.)
print *, '>>> ORI FROM MATRIX'
call o%ori_from_rotmat(rotmat,.true.)
print *, '>>> REJECT ORI'
call o%reject()
print *, '>>> COPY ORI'
call o%copy(o1)
print *, '>>> MATRIX ',o%get_mat() 
call o%get_axis_angle( vec, angle )
print *, '>>> VECTOR ',vec
print *, '>>> ANGLE ', angle/pi*180.
if(angle /= 0.) THROW_HARD('>>> TEST FAILED angle')
print *, '>>> APPEND IN ORI'
call o%append_ori(o1)
print *, '>>> DELETE ENTRY projname'
call o%delete_entry('projname')
call o%set('projname', 'apoferritin')
str_projname = o%get_str('projname')
print *, '>>> SET PROJNAME ', str_projname%to_char()
call o%set('projrec', 'yes')
str_projrec = o%get_str('projrec')
print *, '>>> SET PROJREC ', str_projrec%to_char()
call o%set('kv', 300)
print *, '>>> KEY kv ',o%get('kv')
call o%set('cs', 2.7d0)
print *, '>>> KEY cs ',o%get('cs')
print *, '>>> KEY cs is there ',o%isthere('cs')
if(.not. o%isthere('cs')) THROW_HARD('TEST FAILED; isthere') 
print *, '>>> KEY projname is character ',o%ischar('projname')
if(.not. o%ischar('projname')) THROW_HARD('TEST FAILED; projname')
str_ori = o%ori2str()
print *,'>>> ORI TO STRING ', str_ori%to_char()
print *,'>>> ORI TO JSON'
call o%ori2json(json_ori, .true.)
print *,'>>> STRLEN_TRIM',o%ori_strlen_trim()
call cline%set('projname', projname)
call cline%set('mkdir',        'no')
call cline%check()
call xnew_project%execute_safe(cline)
projfile=projname%to_char()//'.simple'
call spproj%read_segment('compenv', projfile)
print *,'>>> ORI2CHASH '
call spproj%compenv%get_ori(1, compenv_o)
qdescr = compenv_o%ori2chash()
print *,'>>> CHASH2ORI'
call o%chash2ori(qdescr)
print *, '>>> GET KEYS '
keys=o%get_keys()
print *, '>>> GET CTFVARS ',o%get_ctfvars()
print *, '>>> PRINT ORI'
call o%print_ori()
call o_truth%kill()
call o%kill()
call o1%kill()
if( test_passed )then
    print *, '>>> TEST PASSED'
else
    THROW_HARD('>>> TEST FAILED')
endif
end program simple_test_ori

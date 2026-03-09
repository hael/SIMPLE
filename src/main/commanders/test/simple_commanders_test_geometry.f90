!@descr: for all geometry tests
module simple_commanders_test_geometry
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_angres
  contains
    procedure :: execute      => exec_test_angres
end type commander_test_angres

type, extends(commander_base) :: commander_test_ori_test
  contains
    procedure :: execute      => exec_test_ori_test
end type commander_test_ori_test

type, extends(commander_base) :: commander_test_oris_test
  contains
    procedure :: execute      => exec_test_oris_test
end type commander_test_oris_test

type, extends(commander_base) :: commander_test_sym_test
  contains
    procedure :: execute      => exec_test_sym_test
end type commander_test_sym_test

type, extends(commander_base) :: commander_test_uniform_euler
  contains
    procedure :: execute      => exec_test_uniform_euler
end type commander_test_uniform_euler

type, extends(commander_base) :: commander_test_uniform_rot
  contains
    procedure :: execute      => exec_test_uniform_rot
end type commander_test_uniform_rot

contains

subroutine exec_test_angres( self, cline )
    class(commander_test_angres), intent(inout) :: self
    class(cmdline),               intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ANGRES_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_angres

subroutine exec_test_ori_test( self, cline )
    use json_kinds
    use json_module
    use simple_cmdline,    only: cmdline
    use simple_sp_project, only: sp_project
    use simple_commanders_project_core
    class(commander_test_ori_test), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    type(json_value), pointer   :: json_ori
    type(commander_new_project) :: xnew_project
    type(ori)                   :: o_truth, o, o1, compenv_o
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
    call xnew_project%execute(cline)
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
    call simple_end('**** SIMPLE_TEST_ORI_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ori_test

subroutine exec_test_oris_test( self, cline )
    class(commander_test_oris_test), intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    type(oris) :: o, o1, o2, o_subset, o1_subset
    logical    :: test_passed
    integer    :: inds(10)=[1,2,3,4,5,6,7,8,9,10]
    test_passed=.true.
    print *,'>>> CREATE ORIS '
    call o%new(100, is_ptcl=.false.)
    o1 = oris(100, is_ptcl=.false.)
    print *,'>>> REALLOCATE '
    call o1%reallocate(200)
    if( o1%get_noris() .ne. 200 ) test_passed = .false.
    if( .not. test_passed ) THROW_HARD('reallocate oris failed!')
    print *,'>>> EXTRACT SUBSET'
    o_subset=o%extract_subset(1, 10)
    o1_subset=o1%extract_subset(inds)
    print *,'>>> ORIS ELEMENT EXISTS ', o%exists(1)
    if(.not. o%exists(1)) test_passed=.false.
    call o%rnd_oris(5.)
    print *,'>>> ORIS WRITE'
    call o%write(string('test_oris_rndoris.txt'))
    print *,'>>> ORIS READ'
    call o2%read(string('test_oris_rndoris.txt'))
    print *,'>>> ORIS WRITE 2'
    call o2%write(string('test_oris_rndoris_copy.txt'))
    call o%kill()
    call o1%kill()
    call o_subset%kill()
    call o1_subset%kill()
    if( test_passed )then
       print *, '>>> TEST PASSED'
    else
       THROW_HARD('>>> TEST FAILED')
    endif
    call simple_end('**** SIMPLE_TEST_ORIS_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_oris_test

subroutine exec_test_sym_test( self, cline )
    use simple_core_module_api
    use simple_sym, only: sym_tester
    class(commander_test_sym_test), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    type(sym)                     :: se
    integer                       :: nsubgrps
    call se%new('d6')
    nsubgrps = se%get_nsubgrp()
    ! print *, '# subgroups: ', nsubgrps
    ! do i=1,nsubgrps
    !     print *, se%get_subgrp_descr(i)
    ! end do
    call sym_tester('c1')
    call sym_tester('c4')
    call sym_tester('c5')
    call sym_tester('d2')
    call sym_tester('d7')
    call sym_tester('t')
    call sym_tester('o')
    call sym_tester('i')
    call simple_end('**** SIMPLE_TEST_SYM_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_sym_test

subroutine exec_test_uniform_euler( self, cline )
    class(commander_test_uniform_euler), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    integer, parameter :: N_SAMPLES = 2000
    type(sym)   :: pgrpsyms
    type(ori)   :: o
    integer     :: i
    call pgrpsyms%new('c1')
    call o%new(.true.)    
    call pgrpsyms%rnd_euler(o)
    do i = 1, N_SAMPLES
        call pgrpsyms%rnd_euler(o)
        print "(f20.15, f20.15, f20.15)", o%get_normal()*100
    enddo
    call pgrpsyms%kill
    call o%kill
    call simple_end('**** SIMPLE_TEST_UNIFORM_EULER_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_uniform_euler

subroutine exec_test_uniform_rot( self, cline )
    class(commander_test_uniform_rot),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    integer, parameter  :: N_SAMPLES = 4, N_ARR = 2000
    real,    parameter  :: SIGMA = 1.
    real    :: unif(N_SAMPLES), u1, u2, mat(3,3), q0, q1, q2, q3
    integer :: i, j
    call srand(time())
    do i = 1, N_ARR
        do j = 1, N_SAMPLES
            call rgauss(SIGMA, u1, u2)
            unif(j) = u1
        enddo
        unif = unif/sqrt(sum(unif**2))
        ! print "(f20.15, f20.15, f20.15, f20.15)", unif
        q0 = unif(1)
        q1 = unif(2)
        q2 = unif(3)
        q3 = unif(4)
        ! First row of the rotation matrix
        mat(1,1) = 2 * (q0 * q0 + q1 * q1) - 1;
        mat(1,2) = 2 * (q1 * q2 - q0 * q3);
        mat(1,3) = 2 * (q1 * q3 + q0 * q2);    
        ! Second row of the rotation matrix
        mat(2,1) = 2 * (q1 * q2 + q0 * q3);
        mat(2,2) = 2 * (q0 * q0 + q2 * q2) - 1;
        mat(2,3) = 2 * (q2 * q3 - q0 * q1);
        ! Third row of the rotation matrix
        mat(3,1) = 2 * (q1 * q3 - q0 * q2);
        mat(3,2) = 2 * (q2 * q3 + q0 * q1);
        mat(3,3) = 2 * (q0 * q0 + q3 * q3) - 1;
        print "(f20.15, f20.15, f20.15)", matmul([1., 0., 0.], mat);
        ! print *, mat(1,:)
        ! print *, mat(2,:)
        ! print *, mat(3,:)
    enddo

    call simple_end('**** SIMPLE_TEST_UNIFORM_ROT_WORKFLOW NORMAL STOP ****')
 
    contains

        ! Box—Muller method
        subroutine rgauss( sig, y1, y2 )
            real, intent(in)    :: sig
            real, intent(inout) :: y1, y2
            real :: x1, x2, w
            w = 0. ; x1 = 0. ; x2 = 0.
            do while( ( w .ge. 1.0 ) .or. ( w .eq. 0. ) )
                x1 = 2.0 * rand(0) - 1.0
                x2 = 2.0 * rand(0) - 1.0
                w  = x1 * x1 + x2 * x2
            end do
            w = sig*sqrt( (-2.0 * log( w ) ) / w )
            y1 = x1 * w
            y2 = x2 * w
        end subroutine rgauss

end subroutine exec_test_uniform_rot

end module simple_commanders_test_geometry

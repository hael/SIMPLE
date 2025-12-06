module simple_ori_tester
use simple_ori                 ! the class under test
use simple_string, only: string
use simple_test_utils          ! assertions etc.
use simple_defs                ! dp, STDLEN, etc.
implicit none
private
public :: run_all_ori_tests

real, parameter :: DEG_TOL = 1.0e-3      ! loose angle tolerance in degrees
real, parameter :: EPS     = 1.0e-5      ! small numeric tolerance

contains

    subroutine run_all_ori_tests()
        write(*,'(A)') '**** running all ori tests ****'
        call test_ctor_and_flags()
        call test_euler_basic_set_get()
        call test_euler_matrix_roundtrip()
        call test_normal_and_mat_consistency()
        call test_m2euler_fast_vs_m2euler()
        call test_2Dshift_set_get_and_rounding()
        call test_state_class_proj_eo_sampled_updatecnt()
        call test_get_dfx_dfy_and_setters()
        call test_isthere_and_ischar()
        call test_copy_ptcl_to_ptcl()
        call test_copy_ptcl_to_nonptcl()
        call test_copy_nonptcl_to_ptcl()
        call test_delete_2Dclustering()
        call test_delete_3Dalignment()
        call test_transfer_2Dparams()
        call test_transfer_3Dparams()
        call test_mirror3d_involution()
        call test_mirror2d_involution()
        call test_transp_involution()
        call test_geodesic_metrics()
        call test_euler_compose_vs_compeuler()
        call test_ori2str_and_str2ori_roundtrip()
        ! call report_summary()
    end subroutine run_all_ori_tests

    !---------------- basic lifecycle / flags ----------------

    subroutine test_ctor_and_flags()
        type(ori) :: o1, o2
        logical :: ex, ispt
        write(*,'(A)') 'test_ctor_and_flags'
        call o1%new_ori(.true.)
        ex   = o1%exists()
        ispt = o1%is_particle()
        call assert_true(ex,                    'new_ori(.true.) sets existence')
        call assert_true(ispt,                  'new_ori(.true.) sets is_particle=.true.')
        call o2%new_ori(.false.)
        call assert_true(o2%exists(),           'new_ori(.false.) sets existence')
        call assert_true(.not. o2%is_particle(),'new_ori(.false.) sets is_particle=.false.')
        call o1%kill()
        call assert_true(.not. o1%exists(),     'kill() resets existence=.false.')
    end subroutine test_ctor_and_flags

    !---------------- Euler basic set/get ----------------

    subroutine test_euler_basic_set_get()
        type(ori) :: o
        real :: e(3)
        write(*,'(A)') 'test_euler_basic_set_get'
        call o%new_ori(.false.)
        call o%set_euler([1.0, 2.0, 3.0])
        e = o%get_euler()
        call assert_real(1.0, e(1), DEG_TOL, 'set_euler/get_euler e1')
        call assert_real(2.0, e(2), DEG_TOL, 'set_euler/get_euler e2')
        call assert_real(3.0, e(3), DEG_TOL, 'set_euler/get_euler e3')
        call o%e1set(99.0)
        call o%e2set(98.0)
        call o%e3set(97.0)
        e = o%get_euler()
        call assert_real(99.0, e(1), DEG_TOL, 'e1set then get_euler')
        call assert_real(98.0, e(2), DEG_TOL, 'e2set then get_euler')
        call assert_real(97.0, e(3), DEG_TOL, 'e3set then get_euler')
    end subroutine test_euler_basic_set_get

    !---------------- Euler <-> rotation matrix roundtrip ----------------

    subroutine test_euler_matrix_roundtrip()
        type(ori) :: o
        real :: e_in(3), e_out(3), R(3,3)
        integer :: i
        write(*,'(A)') 'test_euler_matrix_roundtrip'
        e_in = [20.0, 40.0, 60.0]
        call o%new_ori(.false.)
        call o%set_euler(e_in)
        R     = euler2m(e_in)
        e_out = m2euler(R)
        do i = 1,3
            call assert_real(e_in(i), e_out(i), 1.0e-2, 'euler2m/m2euler roundtrip component')
        end do
    end subroutine test_euler_matrix_roundtrip

    !---------------- normal vs matrix consistency ----------------

    subroutine test_normal_and_mat_consistency()
        type(ori) :: o
        real :: e(3), n1(3), n2(3), R(3,3), z(3)
        write(*,'(A)') 'test_normal_and_mat_consistency'
        call o%new_ori(.false.)
        e = [10.0, 50.0, 30.0]
        call o%set_euler(e)
        n1 = o%get_normal()
        R  = o%get_mat()
        z  = [0.0, 0.0, 1.0]
        n2 = matmul(z, R)
        call assert_real(n1(1), n2(1), EPS, 'get_normal vs get_mat(1)')
        call assert_real(n1(2), n2(2), EPS, 'get_normal vs get_mat(2)')
        call assert_real(n1(3), n2(3), EPS, 'get_normal vs get_mat(3)')
    end subroutine test_normal_and_mat_consistency

    !---------------- m2euler_fast vs m2euler ----------------

    subroutine test_m2euler_fast_vs_m2euler()
        type(ori) :: o
        real :: e0(3), R(3,3), e1(3), e2(3)
        write(*,'(A)') 'test_m2euler_fast_vs_m2euler'
        e0 = [35.0, 70.0, 110.0]
        call o%new_ori(.false.)
        call o%set_euler(e0)
        R  = o%get_mat()
        e1 = m2euler(R)
        e2 = m2euler_fast(R)
        call assert_real(e1(1), e2(1), 1.0e-2, 'm2euler_fast vs m2euler e1')
        call assert_real(e1(2), e2(2), 1.0e-2, 'm2euler_fast vs m2euler e2')
        call assert_real(e1(3), e2(3), 1.0e-2, 'm2euler_fast vs m2euler e3')
    end subroutine test_m2euler_fast_vs_m2euler

    !---------------- 2D shift set/get/rounding ----------------

    subroutine test_2Dshift_set_get_and_rounding()
        type(ori) :: o
        real :: v(2), off(2)
        write(*,'(A)') 'test_2Dshift_set_get_and_rounding'
        call o%new_ori(.true.)
        call o%set_shift([1.2, -3.4])
        v = o%get_2Dshift()
        call assert_real(1.2,  v(1), 1.0e-6, 'set_shift/get_2Dshift x')
        call assert_real(-3.4, v(2), 1.0e-6, 'set_shift/get_2Dshift y')
        call o%round_shifts()
        v = o%get_2Dshift()
        call assert_int(1,  nint(v(1)), 'round_shifts x')
        call assert_int(-3, nint(v(2)), 'round_shifts y')
        ! calc_offset2D must be consistent with psi
        call o%set_euler([0.0, 0.0, 90.0])  ! 90° in-plane
        call o%set_shift([1.0, 0.0])        ! shift along x
        call o%calc_offset2D(off)
        call assert_real(0.0, off(1), 1.0e-5, 'calc_offset2D offset(1)')
        call assert_real(1.0, off(2), 1.0e-5, 'calc_offset2D offset(2)')
    end subroutine test_2Dshift_set_get_and_rounding

    !---------------- basic per-particle scalar getters ----------------

    subroutine test_state_class_proj_eo_sampled_updatecnt()
        type(ori) :: o
        write(*,'(A)') 'test_state_class_proj_eo_sampled_updatecnt'
        call o%new_ori(.true.)
        call o%set_state(  2)
        call o%set_class(  5)
        call o%set('proj', 7)
        call o%set('eo',   1)
        call o%set('sampled', 3)
        call o%set('updatecnt', 9)
        call assert_int(2, o%get_state(),     'get_state')
        call assert_int(5, o%get_class(),     'get_class')
        call assert_int(7, o%get_proj(),      'get_proj')
        call assert_int(1, o%get_eo(),        'get_eo')
        call assert_int(3, o%get_sampled(),   'get_sampled')
        call assert_int(9, o%get_updatecnt(), 'get_updatecnt')
    end subroutine test_state_class_proj_eo_sampled_updatecnt

    !---------------- dfx/dfy setters and getters ----------------

    subroutine test_get_dfx_dfy_and_setters()
        type(ori) :: o
        write(*,'(A)') 'test_get_dfx_dfy_and_setters'
        call o%new_ori(.true.)
        call o%set_dfx(1.23)
        call o%set_dfy(4.56)
        call assert_real(1.23, o%get_dfx(), 1.0e-6, 'set_dfx/get_dfx')
        call assert_real(4.56, o%get_dfy(), 1.0e-6, 'set_dfy/get_dfy')
    end subroutine test_get_dfx_dfy_and_setters

    !---------------- isthere / ischar ----------------

    subroutine test_isthere_and_ischar()
        type(ori) :: o
        type(string) :: s
        write(*,'(A)') 'test_isthere_and_ischar'
        call o%new_ori(.true.)
        ! pparms-driven key
        call o%set_state(1)
        call assert_true(o%isthere('state'),     'isthere sees mapped pparms key')
        ! numeric dynamic key
        call o%set('corr', 0.5)
        call assert_true(o%isthere('corr'),      'isthere sees numeric hash key')
        call assert_true(.not. o%ischar('corr'), 'corr is numeric, not char')
        ! char key
        call o%set('tag', 'ABC')
        call assert_true(o%isthere('tag'),       'isthere sees char-key')
        call assert_true(o%ischar('tag'),        'ischar true for chtab key')
        s = o%get_str('tag')
        call assert_char('ABC', s%to_char(),     'get(char) returns correct string')
    end subroutine test_isthere_and_ischar

    !---------------- copy: ptcl -> ptcl ----------------

    subroutine test_copy_ptcl_to_ptcl()
        type(ori)   :: a, b
        type(string):: tag
        write(*,'(A)') 'test_copy_ptcl_to_ptcl'
        call a%new_ori(.true.)
        call a%set_state(3)
        call a%set_class(7)
        call a%set_euler([11.0, 22.0, 33.0])
        call a%set_shift([1.0, 2.0])
        call a%set('corr', 0.77)
        call a%set('tag',  'XYZ')
        call b%new_ori(.true.)
        b = a   ! uses copy via assignment interface
        call assert_true( b%is_particle(),     'copy: target still particle')
        call assert_int(3, b%get_state(),      'copy: state copied')
        call assert_int(7, b%get_class(),      'copy: class copied')
        call assert_real(0.77, b%get('corr'), 1.0e-6, 'copy: corr copied')
        call assert_true(b%ischar('tag'),      'copy: tag char present')
        tag = b%get_str('tag')
        call assert_char('XYZ', tag%to_char(), 'copy: tag value')
    end subroutine test_copy_ptcl_to_ptcl

    !---------------- copy: ptcl -> non-ptcl ----------------

    subroutine test_copy_ptcl_to_nonptcl()
        type(ori) :: a, b
        real :: e(3), sh(2)
        write(*,'(A)') 'test_copy_ptcl_to_nonptcl'
        call a%new_ori(.true.)
        call a%set_state(2)
        call a%set_class(5)
        call a%set_euler([10.0, 20.0, 30.0])
        call a%set_shift([3.0, 4.0])
        call b%new_ori(.false.)
        b = a
        call assert_true(.not. b%is_particle(), 'ptcl->nonptcl: target non-ptcl')
        call assert_int(2, b%get_state(), 'ptcl->nonptcl: state mapped to hash')
        call assert_int(5, b%get_class(), 'ptcl->nonptcl: class mapped to hash')
        e  = b%get_euler()
        sh = b%get_2Dshift()
        call assert_real(10.0, e(1), 1.0e-3, 'ptcl->nonptcl e1')
        call assert_real(20.0, e(2), 1.0e-3, 'ptcl->nonptcl e2')
        call assert_real(30.0, e(3), 1.0e-3, 'ptcl->nonptcl e3')
        call assert_real(3.0,  sh(1), 1.0e-6, 'ptcl->nonptcl shift x')
        call assert_real(4.0,  sh(2), 1.0e-6, 'ptcl->nonptcl shift y')
    end subroutine test_copy_ptcl_to_nonptcl

    !---------------- copy: non-ptcl -> ptcl ----------------

    subroutine test_copy_nonptcl_to_ptcl()
        type(ori) :: a, b
        real :: e(3), sh(2)
        write(*,'(A)') 'test_copy_nonptcl_to_ptcl'

        call a%new_ori(.false.)
        call a%set('state', 4)
        call a%set('class', 8)
        call a%set('e1', 15.0)
        call a%set('e2', 25.0)
        call a%set('e3', 35.0)
        call a%set('x',  -1.0)
        call a%set('y',   2.5)

        call b%new_ori(.true.)
        b = a

        call assert_true(b%is_particle(), 'nonptcl->ptcl: target particle')
        call assert_int(4, b%get_state(), 'nonptcl->ptcl state')
        call assert_int(8, b%get_class(), 'nonptcl->ptcl class')

        e  = b%get_euler()
        sh = b%get_2Dshift()
        call assert_real(15.0, e(1), 1.0e-3, 'nonptcl->ptcl e1')
        call assert_real(25.0, e(2), 1.0e-3, 'nonptcl->ptcl e2')
        call assert_real(35.0, e(3), 1.0e-3, 'nonptcl->ptcl e3')
        call assert_real(-1.0, sh(1), 1.0e-6, 'nonptcl->ptcl shift x')
        call assert_real( 2.5, sh(2), 1.0e-6, 'nonptcl->ptcl shift y')
    end subroutine test_copy_nonptcl_to_ptcl

    !---------------- delete_2Dclustering ----------------

    subroutine test_delete_2Dclustering()
        type(ori) :: o
        real :: sh(2)
        write(*,'(A)') 'test_delete_2Dclustering'

        call o%new_ori(.true.)
        call o%set_class(10)
        call o%set_euler([0.0, 0.0, 45.0])
        call o%set_shift([5.0, -7.0])
        call o%set('corr', 0.9)
        call o%set('frac', 0.5)

        call o%delete_2Dclustering(keepshifts=.false., keepcls=.false.)
        sh = o%get_2Dshift()

        call assert_int(0, o%get_class(), 'delete_2Dclustering clears class')
        call assert_real(0.0, o%e3get(), 1.0e-6, 'delete_2Dclustering clears e3')
        call assert_real(0.0, sh(1), 1.0e-6, 'delete_2Dclustering clears x')
        call assert_real(0.0, sh(2), 1.0e-6, 'delete_2Dclustering clears y')
        call assert_real(0.0, o%get('corr'), 1.0e-6, 'delete_2Dclustering clears corr')
        call assert_real(0.0, o%get('frac'), 1.0e-6, 'delete_2Dclustering clears frac')
    end subroutine test_delete_2Dclustering

    !---------------- delete_3Dalignment ----------------

    subroutine test_delete_3Dalignment()
        type(ori) :: o
        real :: sh(2)
        write(*,'(A)') 'test_delete_3Dalignment'

        call o%new_ori(.true.)
        call o%set('proj', 2)
        call o%set_euler([10.0, 20.0, 30.0])
        call o%set_shift([3.0, -4.0])
        call o%set('corr', 0.7)
        call o%set('frac', 0.4)

        call o%delete_3Dalignment(keepshifts=.false.)
        sh = o%get_2Dshift()

        call assert_int(0, o%get_proj(),  'delete_3Dalignment clears proj')
        call assert_real(0.0, o%e1get(), 1.0e-6, 'delete_3Dalignment clears e1')
        call assert_real(0.0, o%e2get(), 1.0e-6, 'delete_3Dalignment clears e2')
        call assert_real(0.0, o%e3get(), 1.0e-6, 'delete_3Dalignment clears e3')
        call assert_real(0.0, sh(1), 1.0e-6, 'delete_3Dalignment clears x')
        call assert_real(0.0, sh(2), 1.0e-6, 'delete_3Dalignment clears y')
        call assert_real(0.0, o%get('corr'), 1.0e-6, 'delete_3Dalignment clears corr')
        call assert_real(0.0, o%get('frac'), 1.0e-6, 'delete_3Dalignment clears frac')
    end subroutine test_delete_3Dalignment

    !---------------- transfer_2Dparams ----------------

    subroutine test_transfer_2Dparams()
        type(ori) :: src, dst
        real :: sh(2)
        write(*,'(A)') 'test_transfer_2Dparams'

        call src%new_ori(.true.)
        call dst%new_ori(.true.)

        call src%set_class(5)
        call src%set('corr', 0.8)
        call src%set('frac', 0.3)
        call src%set('sampled', 1)
        call src%set('updatecnt', 4)
        call src%set('w', 2.0)
        call src%set('eo', 1)
        call src%set_euler([10.0, 20.0, 30.0])
        call src%set_shift([1.0, 2.0])

        call dst%transfer_2Dparams(src)

        sh = dst%get_2Dshift()

        call assert_int(5, dst%get_class(), 'transfer_2Dparams class')
        call assert_real(0.8, dst%get('corr'), 1.0e-6, 'transfer_2Dparams corr')
        call assert_real(0.3, dst%get('frac'), 1.0e-6, 'transfer_2Dparams frac')
        call assert_int(1, dst%get_sampled(), 'transfer_2Dparams sampled')
        call assert_int(4, dst%get_updatecnt(), 'transfer_2Dparams updatecnt')
        call assert_real(2.0, dst%get('w'), 1.0e-6, 'transfer_2Dparams w')
        call assert_int(1, dst%get_eo(), 'transfer_2Dparams eo')
        call assert_real(10.0, dst%e1get(), 1.0e-3, 'transfer_2Dparams e1')
        call assert_real(20.0, dst%e2get(), 1.0e-3, 'transfer_2Dparams e2')
        call assert_real(30.0, dst%e3get(), 1.0e-3, 'transfer_2Dparams e3')
        call assert_real(1.0, sh(1), 1.0e-6, 'transfer_2Dparams x')
        call assert_real(2.0, sh(2), 1.0e-6, 'transfer_2Dparams y')
    end subroutine test_transfer_2Dparams

    !---------------- transfer_3Dparams ----------------

    subroutine test_transfer_3Dparams()
        type(ori) :: src, dst
        real :: sh(2)
        write(*,'(A)') 'test_transfer_3Dparams'
        call src%new_ori(.true.)
        call dst%new_ori(.true.)
        call src%set('proj', 3)
        call src%set('corr', 0.9)
        call src%set('frac', 0.2)
        call src%set('sampled', 1)
        call src%set('updatecnt', 5)
        call src%set('w', 3.0)
        call src%set('eo', 2)
        call src%set_euler([15.0, 25.0, 35.0])
        call src%set_shift([2.0, -1.0])
        call dst%transfer_3Dparams(src)
        sh = dst%get_2Dshift()
        call assert_int(3, dst%get_proj(), 'transfer_3Dparams proj')
        call assert_real(0.9, dst%get('corr'), 1.0e-6, 'transfer_3Dparams corr')
        call assert_real(0.2, dst%get('frac'), 1.0e-6, 'transfer_3Dparams frac')
        call assert_int(1, dst%get_sampled(), 'transfer_3Dparams sampled')
        call assert_int(5, dst%get_updatecnt(), 'transfer_3Dparams updatecnt')
        call assert_real(3.0, dst%get('w'), 1.0e-6, 'transfer_3Dparams w')
        call assert_int(2, dst%get_eo(), 'transfer_3Dparams eo')
        call assert_real(15.0, dst%e1get(), 1.0e-3, 'transfer_3Dparams e1')
        call assert_real(25.0, dst%e2get(), 1.0e-3, 'transfer_3Dparams e2')
        call assert_real(35.0, dst%e3get(), 1.0e-3, 'transfer_3Dparams e3')
        call assert_real( 2.0, sh(1), 1.0e-6, 'transfer_3Dparams x')
        call assert_real(-1.0, sh(2), 1.0e-6, 'transfer_3Dparams y')
    end subroutine test_transfer_3Dparams

    !---------------- mirror3d is an involution ----------------

    subroutine test_mirror3d_involution()
        type(ori) :: o
        real :: e0(3), e2(3)
        write(*,'(A)') 'test_mirror3d_involution'
        e0 = [10.0, 40.0, 70.0]
        call o%new_ori(.false.)
        call o%set_euler(e0)
        call o%mirror3d()
        call o%mirror3d()
        e2 = o%get_euler()
        call assert_real(e0(1), e2(1), 1.0e-2, 'mirror3d twice returns original e1')
        call assert_real(e0(2), e2(2), 1.0e-2, 'mirror3d twice returns original e2')
        call assert_real(e0(3), e2(3), 1.0e-2, 'mirror3d twice returns original e3')
    end subroutine test_mirror3d_involution

    !---------------- mirror2d is an involution ----------------

    subroutine test_mirror2d_involution()
        type(ori) :: o
        real :: e0(3), e2(3)
        write(*,'(A)') 'test_mirror2d_involution'
        e0 = [50.0, 60.0, 120.0]
        call o%new_ori(.false.)
        call o%set_euler(e0)
        call o%mirror2d()
        call o%mirror2d()
        e2 = o%get_euler()
        call assert_real(e0(1), e2(1), 1.0e-2, 'mirror2d twice returns original e1')
        call assert_real(e0(2), e2(2), 1.0e-2, 'mirror2d twice returns original e2')
        call assert_real(e0(3), e2(3), 1.0e-2, 'mirror2d twice returns original e3')
    end subroutine test_mirror2d_involution

    !---------------- transpose is an involution ----------------

    subroutine test_transp_involution()
        type(ori) :: o
        real :: e0(3), e2(3)
        write(*,'(A)') 'test_transp_involution'
        e0 = [33.0, 77.0, 123.0]
        call o%new_ori(.false.)
        call o%set_euler(e0)
        call o%transp()
        call o%transp()
        e2 = o%get_euler()
        call assert_real(e0(1), e2(1), 1.0e-2, 'transp twice returns original e1')
        call assert_real(e0(2), e2(2), 1.0e-2, 'transp twice returns original e2')
        call assert_real(e0(3), e2(3), 1.0e-2, 'transp twice returns original e3')
    end subroutine test_transp_involution

    !---------------- geodesic metrics ----------------

    subroutine test_geodesic_metrics()
        type(ori) :: a, b
        real :: d_frob, d_trace, e1(3), e2(3)
        write(*,'(A)') 'test_geodesic_metrics'
        call a%new_ori(.false.)
        call b%new_ori(.false.)
        e1 = [0.0, 0.0, 0.0]
        e2 = [0.0, 30.0, 0.0]
        call a%set_euler(e1)
        call b%set_euler(e1)
        d_frob  = a .geod. b
        d_trace = a%geodesic_dist_trace(b)
        call assert_real(0.0, d_frob,  1.0e-6, 'geodesic_frobdev zero for equal')
        call assert_real(0.0, d_trace, 1.0e-6, 'geodesic_dist_trace zero for equal')
        call b%set_euler(e2)
        d_frob  = a .geod. b
        d_trace = a%geodesic_dist_trace(b)
        call assert_true(d_frob  > 0.0, 'geodesic_frobdev > 0 for different')
        call assert_true(d_trace > 0.0, 'geodesic_dist_trace > 0 for different')
    end subroutine test_geodesic_metrics

    !---------------- euler_compose vs compeuler ----------------

    subroutine test_euler_compose_vs_compeuler()
        type(ori) :: o1, o2, o_out
        real :: e1(3), e2(3), e_comp(3), e_out(3)
        integer :: i
        write(*,'(A)') 'test_euler_compose_vs_compeuler'
        e1 = [10.0, 20.0, 30.0]
        e2 = [5.0,  15.0, 25.0]
        call o1%new_ori(.false.)
        call o2%new_ori(.false.)
        call o_out%new_ori(.false.)
        call o1%set_euler(e1)
        call o2%set_euler(e2)
        call euler_compose(e1, e2, e_comp)
        call o1%compose(o2, o_out)
        e_out = o_out%get_euler()
        do i = 1,3
            call assert_real(e_comp(i), e_out(i), 1.0e-2, 'compose vs euler_compose component')
        end do
    end subroutine test_euler_compose_vs_compeuler

    !---------------- ori2str / str2ori roundtrip ----------------

    subroutine test_ori2str_and_str2ori_roundtrip()
        type(ori)    :: o1, o2
        type(string) :: s
        real :: e1(3), e2(3), sh1(2), sh2(2)
        write(*,'(A)') 'test_ori2str_and_str2ori_roundtrip'
        call o1%new_ori(.true.)
        call o1%set_state(2)
        call o1%set_class(5)
        call o1%set_euler([12.0, 34.0, 56.0])
        call o1%set_shift([-1.5, 2.5])
        call o1%set('corr', 0.77)
        call o1%set('tag',  'ABC')
        s = o1%ori2str()
        call o2%str2ori(s%to_char(), is_ptcl=.true.)
        call assert_true(o2%is_particle(), 'str2ori preserves is_ptcl')
        call assert_int(2, o2%get_state(), 'str2ori: state')
        call assert_int(5, o2%get_class(), 'str2ori: class')
        e1  = o1%get_euler()
        e2  = o2%get_euler()
        sh1 = o1%get_2Dshift()
        sh2 = o2%get_2Dshift()
        call assert_real(e1(1), e2(1), 1.0e-3, 'ori2str/str2ori e1')
        call assert_real(e1(2), e2(2), 1.0e-3, 'ori2str/str2ori e2')
        call assert_real(e1(3), e2(3), 1.0e-3, 'ori2str/str2ori e3')
        call assert_real(sh1(1), sh2(1), 1.0e-6, 'ori2str/str2ori shift x')
        call assert_real(sh1(2), sh2(2), 1.0e-6, 'ori2str/str2ori shift y')
        call assert_real(o1%get('corr'), o2%get('corr'), 1.0e-6, 'ori2str/str2ori corr')
        s = o2%get_str('tag')
        call assert_char('ABC', s%to_char(), 'ori2str/str2ori tag')
    end subroutine test_ori2str_and_str2ori_roundtrip

end module simple_ori_tester

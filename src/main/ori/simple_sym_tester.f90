!@descr: unit test routines for the sym class (point-group symmetries)
module simple_sym_tester
use simple_test_utils ! assertions etc.
use simple_defs       ! PI etc.
use simple_ori,       only: ori
use simple_oris,      only: oris
use simple_sym,       only: sym, is_valid_pointgroup
implicit none
private
public :: run_all_sym_tests

real,    parameter :: DEG_TOL   = 1.0e-2 ! angle tolerance in degrees
real,    parameter :: MAT_TOL   = 1.0e-3 ! rotation-matrix element tolerance
integer, parameter :: NSPIRAL   = 400    ! projection directions per reference spiral (must be even)
integer, parameter :: NTRIALS   = 20     ! random orientations per point group
character(len=3), parameter :: PGRPS(9)  = ['c1 ','c2 ','c4 ','c5 ','d2 ','d7 ','t  ','o  ','i  ']
integer,          parameter :: NSYMS(9)  = [ 1,    2,    4,    5,    4,    14,   12,   24,   60  ]

contains

    subroutine run_all_sym_tests()
        write(*,'(A)') '**** running all sym tests ****'
        call test_is_valid_pointgroup()
        call test_construction_and_order()
        call test_eullims()
        call test_subgroups()
        call test_symops_form_a_group()
        call test_apply()
        call test_rnd_euler_within_limits()
        call test_rot_to_asym()
        call test_symrandomize()
        call test_build_refspiral()
        ! call report_summary()
    end subroutine run_all_sym_tests

    !---------------- point-group string validation ----------------

    subroutine test_is_valid_pointgroup()
        write(*,'(A)') 'test_is_valid_pointgroup'
        call assert_true(is_valid_pointgroup('c1'),        'c1 is a valid point group')
        call assert_true(is_valid_pointgroup('c12'),       'c12 is a valid point group')
        call assert_true(is_valid_pointgroup('C3'),        'upper-case C3 is a valid point group')
        call assert_true(is_valid_pointgroup('d2'),        'd2 is a valid point group')
        call assert_true(is_valid_pointgroup('D4'),        'upper-case D4 is a valid point group')
        call assert_true(is_valid_pointgroup('t'),         't is a valid point group')
        call assert_true(is_valid_pointgroup('o'),         'o is a valid point group')
        call assert_true(is_valid_pointgroup('i'),         'i is a valid point group')
        call assert_true(.not. is_valid_pointgroup('c0'),  'c0 is rejected')
        call assert_true(.not. is_valid_pointgroup('d1'),  'd1 is rejected')
        call assert_true(.not. is_valid_pointgroup('c'),   'bare c is rejected')
        call assert_true(.not. is_valid_pointgroup('tt'),  'tt is rejected')
        call assert_true(.not. is_valid_pointgroup('x2'),  'unknown letter is rejected')
    end subroutine test_is_valid_pointgroup

    !---------------- construction, order and classification ----------------

    subroutine test_construction_and_order()
        type(sym) :: se
        integer   :: i
        write(*,'(A)') 'test_construction_and_order'
        do i = 1,size(PGRPS)
            call se%new(trim(PGRPS(i)))
            call assert_int(NSYMS(i), se%get_nsym(),            'get_nsym '//trim(PGRPS(i)))
            call assert_char(trim(PGRPS(i)), trim(se%get_pgrp()),'get_pgrp '//trim(PGRPS(i)))
            call assert_true(se%is_circular()    .eqv. (PGRPS(i)(1:1) == 'c'), 'is_circular '//trim(PGRPS(i)))
            call assert_true(se%is_dihedral()    .eqv. (PGRPS(i)(1:1) == 'd'), 'is_dihedral '//trim(PGRPS(i)))
            call assert_true(se%is_platonic()    .eqv. (index('toi', PGRPS(i)(1:1)) > 0), 'is_platonic '//trim(PGRPS(i)))
            call assert_true(se%is_icosahedral() .eqv. (PGRPS(i)(1:1) == 'i'), 'is_icosahedral '//trim(PGRPS(i)))
            call se%kill
        end do
        ! upper-case input is normalised
        call se%new('C4')
        call assert_char('c4', trim(se%get_pgrp()), 'C4 is normalised to c4')
        call assert_int(4, se%get_nsym(),          'C4 has four operations')
        call se%kill
        ! structure constructor
        se = sym('d3')
        call assert_int(6, se%get_nsym(), 'sym() constructor builds d3')
        call se%kill
    end subroutine test_construction_and_order

    !---------------- Euler-angle limits of the asymmetric unit ----------------

    subroutine test_eullims()
        write(*,'(A)') 'test_eullims'
        call check_lims('c1', 360.0,       180.0)
        call check_lims('c4',  90.0,       180.0)
        call check_lims('c5',  72.0,       180.0)
        call check_lims('d2', 180.0,        90.0)
        call check_lims('d7', 360.0/7.0,    90.0)
        call check_lims('t',  120.0, 70.52877936550931)
        call check_lims('o',   90.0, 54.735610317245346)
        call check_lims('i',   72.0, 37.37736814064969)

        contains

            subroutine check_lims( pgrp, phi_max, theta_max )
                character(len=*), intent(in) :: pgrp
                real,             intent(in) :: phi_max, theta_max
                type(sym) :: se
                real      :: lims(3,2)
                call se%new(pgrp)
                lims = se%get_eullims()
                call assert_real(0.0,       lims(1,1), DEG_TOL, 'eullims phi lower bound '//pgrp)
                call assert_real(phi_max,   lims(1,2), DEG_TOL, 'eullims phi upper bound '//pgrp)
                call assert_real(0.0,       lims(2,1), DEG_TOL, 'eullims theta lower bound '//pgrp)
                call assert_real(theta_max, lims(2,2), DEG_TOL, 'eullims theta upper bound '//pgrp)
                call assert_real(0.0,       lims(3,1), DEG_TOL, 'eullims psi lower bound '//pgrp)
                call assert_real(360.0,     lims(3,2), DEG_TOL, 'eullims psi upper bound '//pgrp)
                call se%kill
            end subroutine check_lims

    end subroutine test_eullims

    !---------------- subgroup tables ----------------

    subroutine test_subgroups()
        type(sym) :: se, sub
        integer   :: i, j, n
        write(*,'(A)') 'test_subgroups'
        do i = 1,size(PGRPS)
            call se%new(trim(PGRPS(i)))
            n = se%get_nsubgrp()
            call assert_true(n >= 1, 'at least one subgroup '//trim(PGRPS(i)))
            sub = se%get_subgrp(1)
            call assert_char(trim(PGRPS(i)), trim(sub%get_pgrp()), 'first subgroup is the group itself '//trim(PGRPS(i)))
            call sub%kill
            if( NSYMS(i) > 1 )then
                call assert_true(se%has_subgrp('c1'), 'c1 is a subgroup '//trim(PGRPS(i)))
                sub = se%get_subgrp(n)
                call assert_char('c1', trim(sub%get_pgrp()), 'last subgroup is c1 '//trim(PGRPS(i)))
                call sub%kill
            endif
            ! Lagrange: the order of every subgroup divides the order of the group
            do j = 1,n
                sub = se%get_subgrp(j)
                call assert_int(0, mod(NSYMS(i), sub%get_nsym()), 'subgroup order divides group order '//trim(PGRPS(i)))
                call sub%kill
            end do
            call se%kill
        end do
        ! spot checks of the tables
        call se%new('d6')
        call assert_true(se%has_subgrp('d3'), 'd6 has subgroup d3')
        call assert_true(se%has_subgrp('c3'), 'd6 has subgroup c3')
        call assert_true(se%has_subgrp('c2'), 'd6 has subgroup c2')
        call assert_true(.not. se%has_subgrp('c4'), 'd6 has no subgroup c4')
        call se%kill
        call se%new('o')
        call assert_int(9, se%get_nsubgrp(),  'o has nine subgroups')
        call assert_true(se%has_subgrp('t'),  'o has subgroup t')
        call assert_true(se%has_subgrp('d4'), 'o has subgroup d4')
        call se%kill
        call se%new('i')
        call assert_true(se%has_subgrp('d5'),       'i has subgroup d5')
        call assert_true(.not. se%has_subgrp('c4'), 'i has no subgroup c4')
        call se%kill
    end subroutine test_subgroups

    !---------------- the symmetry operators form a group ----------------

    subroutine test_symops_form_a_group()
        write(*,'(A)') 'test_symops_form_a_group'
        call check_group('c4')
        call check_group('d2')
        call check_group('d7')
        call check_group('t')
        call check_group('o')
        call check_group('i')

        contains

            subroutine check_group( pgrp )
                character(len=*), intent(in) :: pgrp
                type(sym) :: se
                real, allocatable :: ops(:,:,:)
                real    :: identity(3,3), prod(3,3)
                integer :: n, i, j
                logical :: orthonormal, distinct, closed
                identity = 0.0
                do i = 1,3
                    identity(i,i) = 1.0
                end do
                call se%new(pgrp)
                n = se%get_nsym()
                allocate(ops(3,3,n))
                do i = 1,n
                    call se%get_sym_rmat(i, ops(:,:,i))
                end do
                ! operator 1 is the identity (rot_to_asym relies on it)
                call assert_true(maxval(abs(ops(:,:,1) - identity)) < MAT_TOL, 'symop 1 is the identity '//pgrp)
                ! every operator is a proper rotation
                orthonormal = .true.
                do i = 1,n
                    prod = matmul(ops(:,:,i), transpose(ops(:,:,i)))
                    if( maxval(abs(prod - identity)) > MAT_TOL ) orthonormal = .false.
                    if( abs(det3(ops(:,:,i)) - 1.0) > MAT_TOL )  orthonormal = .false.
                end do
                call assert_true(orthonormal, 'symops are proper rotations '//pgrp)
                ! all operators are distinct
                distinct = .true.
                do i = 1,n-1
                    do j = i+1,n
                        if( maxval(abs(ops(:,:,i) - ops(:,:,j))) < MAT_TOL ) distinct = .false.
                    end do
                end do
                call assert_true(distinct, 'symops are distinct '//pgrp)
                ! closure under composition
                closed = .true.
                do i = 1,n
                    do j = 1,n
                        prod = matmul(ops(:,:,i), ops(:,:,j))
                        if( .not. matches_some_op(prod, ops) ) closed = .false.
                    end do
                end do
                call assert_true(closed, 'symops are closed under composition '//pgrp)
                deallocate(ops)
                call se%kill
            end subroutine check_group

    end subroutine test_symops_form_a_group

    !---------------- apply (ori and Euler-triplet forms agree) ----------------

    subroutine test_apply()
        type(sym) :: se
        type(ori) :: o, o_sym
        real      :: e_in(3), e_sym(3), e_ori(3)
        integer   :: k, i
        logical   :: consistent
        write(*,'(A)') 'test_apply'
        call se%new('c4')
        call o%new_ori(.true.)
        e_in = [20.0, 50.0, 70.0]
        call o%set_euler(e_in)
        call o%set_shift([1.0, -2.0])
        call o%set('corr', 0.5)
        ! the identity operator leaves the orientation unchanged
        call se%apply(o, 1, o_sym)
        e_ori = o_sym%get_euler()
        do i = 1,3
            call assert_real(e_in(i), e_ori(i), DEG_TOL, 'apply with symop 1 is the identity')
        end do
        ! the two forms of apply agree for every operator
        consistent = .true.
        do k = 1,se%get_nsym()
            call se%apply(o, k, o_sym)
            call se%apply(e_in, k, e_sym)
            e_ori = o_sym%get_euler()
            do i = 1,3
                if( abs(e_ori(i) - e_sym(i)) > DEG_TOL ) consistent = .false.
            end do
        end do
        call assert_true(consistent, 'apply(ori) and apply(euls) agree for all symops')
        ! non-Euler parameters are transferred
        call se%apply(o, 2, o_sym)
        call assert_true(o_sym%is_particle(),            'apply preserves is_ptcl')
        call assert_real(0.5, o_sym%get('corr'), 1.0e-6, 'apply transfers corr')
        call assert_real(1.0, o_sym%get('x'),    1.0e-6, 'apply transfers x')
        ! for c4 the second operator is a 90 degree rotation about z: the projection direction moves in phi only
        e_sym = o_sym%get_euler()
        call assert_real(e_in(2), e_sym(2), DEG_TOL, 'c4 symop 2 preserves theta')
        call o%kill
        call o_sym%kill
        call se%kill
    end subroutine test_apply

    !---------------- rnd_euler samples inside the asymmetric unit ----------------

    subroutine test_rnd_euler_within_limits()
        write(*,'(A)') 'test_rnd_euler_within_limits'
        call check_rnd('c4')
        call check_rnd('d2')
        call check_rnd('t')
        call check_rnd('o')
        call check_rnd('i')

        contains

            subroutine check_rnd( pgrp )
                character(len=*), intent(in) :: pgrp
                type(sym) :: se
                type(ori) :: o
                real      :: lims(3,2), e(3)
                integer   :: i
                logical   :: ok
                call se%new(pgrp)
                lims = se%get_eullims()
                call o%new_ori(.false.)
                ok = .true.
                do i = 1,NTRIALS
                    call se%rnd_euler(o)
                    e = o%get_euler()
                    if( .not. within_lims(e, lims) ) ok = .false.
                end do
                call assert_true(ok, 'rnd_euler stays within eullims '//pgrp)
                call o%kill
                call se%kill
            end subroutine check_rnd

    end subroutine test_rnd_euler_within_limits

    !---------------- rot_to_asym maps into the asymmetric unit by a group operation ----------------

    subroutine test_rot_to_asym()
        write(*,'(A)') 'test_rot_to_asym'
        call check_rot('c4')
        call check_rot('d7')
        call check_rot('t')
        call check_rot('o')
        call check_rot('i')

        contains

            subroutine check_rot( pgrp )
                character(len=*), intent(in) :: pgrp
                type(sym) :: se
                type(ori) :: o, orig
                real      :: lims(3,2)
                integer   :: i
                logical   :: in_asu, by_symop
                call se%new(pgrp)
                lims = se%get_eullims()
                call o%new_ori(.false.)
                in_asu   = .true.
                by_symop = .true.
                do i = 1,NTRIALS
                    call o%rnd_euler
                    orig = o
                    call se%rot_to_asym(o)
                    if( .not. within_lims(o%get_euler(), lims) ) in_asu = .false.
                    if( .not. related_by_symop(se, orig, o) )    by_symop = .false.
                end do
                call assert_true(in_asu,   'rot_to_asym lands inside eullims '//pgrp)
                call assert_true(by_symop, 'rot_to_asym applies a group operation '//pgrp)
                call o%kill
                call orig%kill
                call se%kill
            end subroutine check_rot

    end subroutine test_rot_to_asym

    !---------------- symrandomize replaces each orientation by a symmetry mate ----------------

    subroutine test_symrandomize()
        type(sym)  :: se
        type(oris) :: os, os_orig
        type(ori)  :: o, o_orig
        integer    :: i
        logical    :: by_symop
        write(*,'(A)') 'test_symrandomize'
        call se%new('d2')
        call os%new(NTRIALS, is_ptcl=.false.)
        call o%new_ori(.false.)
        do i = 1,NTRIALS
            call se%rnd_euler(o)
            call os%set_ori(i, o)
        end do
        os_orig = os
        call se%symrandomize(os)
        call assert_int(NTRIALS, os%get_noris(), 'symrandomize keeps the number of orientations')
        by_symop = .true.
        do i = 1,NTRIALS
            call os%get_ori(i, o)
            call os_orig%get_ori(i, o_orig)
            if( .not. related_by_symop(se, o_orig, o) ) by_symop = .false.
        end do
        call assert_true(by_symop, 'symrandomize applies group operations only')
        call o%kill
        call o_orig%kill
        call os%kill
        call os_orig%kill
        call se%kill
    end subroutine test_symrandomize

    !---------------- reference spiral: size, redundancy, north pole, mirror pairs ----------------

    subroutine test_build_refspiral()
        integer :: i
        write(*,'(A)') 'test_build_refspiral'
        do i = 1,size(PGRPS)
            call check_spiral(trim(PGRPS(i)))
        end do

        contains

            subroutine check_spiral( pgrp )
                character(len=*), intent(in) :: pgrp
                type(sym)  :: se
                type(oris) :: os
                type(ori)  :: o, north_pole
                real       :: lims(3,2), normals(3,NSPIRAL)
                integer    :: i, j, nredundant, npartner, nhalf
                logical    :: found, mirr_ok, asu_ok
                call se%new(pgrp)
                lims  = se%get_eullims()
                nhalf = NSPIRAL/2
                call os%new(NSPIRAL, is_ptcl=.false.)
                call se%build_refspiral(os)
                call assert_int(NSPIRAL, os%get_noris(), 'build_refspiral keeps the requested size '//pgrp)
                ! no two projection directions coincide (closer than ~0.08 degrees), with one
                ! designed exception: for d/o/i the mirror of the north pole is symmetry-equivalent
                ! to the pole itself, so build_refspiral jitters the pole by <= 0.5 degrees and the
                ! pair (pole, mirror mate) can be arbitrarily close. Such a pair is allowed only
                ! between mirror partners (i, i+nhalf), and at most once.
                do i = 1,NSPIRAL
                    call os%get_ori(i, o)
                    normals(:,i) = o%get_normal()
                end do
                nredundant = 0
                npartner   = 0
                do i = 1,NSPIRAL-1
                    do j = i+1,NSPIRAL
                        if( dot_product(normals(:,i), normals(:,j)) > 1.0 - 1.0e-6 )then
                            if( j == i + nhalf )then
                                npartner = npartner + 1
                            else
                                nredundant = nredundant + 1
                            endif
                            write(*,'(A,2I6,3F9.3,A,3F9.3)') '   near-coincident directions ', i, j, &
                                &os%get_euler(i), ' / ', os%get_euler(j)
                        endif
                    end do
                end do
                call assert_int(0, nredundant,  'build_refspiral has no redundant directions '//pgrp)
                call assert_true(npartner <= 1, 'build_refspiral: at most the jittered pole coincides with its mirror mate '//pgrp)
                ! the north pole is always present
                call north_pole%new_ori(.false.)
                call north_pole%set_euler([0.0, 0.0, 0.0])
                found = .false.
                do i = 1,NSPIRAL
                    call os%get_ori(i, o)
                    if( (o.euldist.north_pole) < 0.01 )then
                        found = .true.
                        exit
                    endif
                end do
                call assert_true(found, 'build_refspiral contains the north pole '//pgrp)
                ! the first half is the un-mirrored asymmetric unit
                asu_ok = .true.
                do i = 1,nhalf
                    if( .not. within_lims(os%get_euler(i), lims) ) asu_ok = .false.
                end do
                call assert_true(asu_ok, 'build_refspiral first half lies inside eullims '//pgrp)
                ! mirror pairs are cross-referenced
                mirr_ok = .true.
                do i = 1,nhalf
                    if( nint(os%get(i,       'mirr')) /= i + nhalf ) mirr_ok = .false.
                    if( nint(os%get(i+nhalf, 'mirr')) /= i )         mirr_ok = .false.
                end do
                call assert_true(mirr_ok, 'build_refspiral flags mirror pairs '//pgrp)
                call o%kill
                call north_pole%kill
                call os%kill
                call se%kill
            end subroutine check_spiral

    end subroutine test_build_refspiral

    !---------------- helpers ----------------

    pure real function det3( m )
        real, intent(in) :: m(3,3)
        det3 = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
            &- m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
            &+ m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))
    end function det3

    pure logical function matches_some_op( R, ops )
        real, intent(in) :: R(3,3), ops(:,:,:)
        integer :: k
        matches_some_op = .false.
        do k = 1,size(ops,3)
            if( maxval(abs(R - ops(:,:,k))) < MAT_TOL )then
                matches_some_op = .true.
                return
            endif
        end do
    end function matches_some_op

    ! phi and theta inside the half-open box of the asymmetric unit (psi is unrestricted)
    pure logical function within_lims( e, lims )
        real, intent(in) :: e(3), lims(3,2)
        within_lims = e(1) >= lims(1,1) - DEG_TOL .and. e(1) < lims(1,2) + DEG_TOL .and. &
                     &e(2) >= lims(2,1) - DEG_TOL .and. e(2) < lims(2,2) + DEG_TOL
    end function within_lims

    ! true when some symmetry operation maps orig onto o (full rotation, not just the direction)
    logical function related_by_symop( se, orig, o )
        class(sym), intent(in) :: se
        class(ori), intent(in) :: orig, o
        type(ori) :: o_k
        integer   :: k
        related_by_symop = .false.
        do k = 1,se%get_nsym()
            call se%apply(orig, k, o_k)
            if( (o_k .geod. o) < MAT_TOL )then
                related_by_symop = .true.
                exit
            endif
        end do
        call o_k%kill
    end function related_by_symop

end module simple_sym_tester

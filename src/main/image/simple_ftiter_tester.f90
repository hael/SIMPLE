!@descr: unit tests for the Fourier index iterator (simple_ftiter): loop limits, logical/physical addressing, resolution conversions
! Replaces the in-module test_ftiter/test_addr (private access, a duplicated block, a stop at the first
! failure). Loop limits and address maps for even, odd and non-square 2D and 3D boxes; the low-pass
! limits, symmetric in the third dimension of a volume since 2026-09-25 (they ran from l = 0, so a
! volume correlation or low-pass product saw half of the half-space); every dimension takes the
! Fourier scale of the first, as the whole Fourier layer indexes isotropically.
module simple_ftiter_tester
use simple_test_utils
use simple_defs
use simple_string_utils, only: int2str
use simple_ftiter,       only: ftiter
implicit none
private
public :: run_all_ftiter_tests

contains

    subroutine run_all_ftiter_tests()
        write(*,'(A)') '**** running all ftiter tests ****'
        call test_logical_limits()
        call test_address_maps()
        call test_resolution()
        call test_lowpass_limits()
        call test_reset()
    end subroutine run_all_ftiter_tests

    function box_tag( ldim ) result( tag )
        integer, intent(in) :: ldim(3)
        character(len=:), allocatable :: tag
        tag = int2str(ldim(1))//'x'//int2str(ldim(2))//'x'//int2str(ldim(3))//': '
    end function box_tag

    ! h from 0 to ldim(1)/2 (Friedel half); k and l from -ldim/2 to ldim/2-1 (even) or -(ldim-1)/2 to
    ! (ldim-1)/2 (odd); a 2D box has l = 0; mode 3 adds the negative h, down to -ldim(1)/2
    subroutine test_logical_limits()
        type(ftiter) :: fit
        write(*,'(A)') 'test_logical_limits'
        call fit%new([100,100,1], 2.)
        call check(fit%loop_lims(2), [0,-50,0, 50,49,0],   '100x100x1: mode 2')
        call check(fit%loop_lims(3), [-50,-50,0, 50,49,0], '100x100x1: mode 3')
        call fit%new([100,100,100], 2.)
        call check(fit%loop_lims(2), [0,-50,-50, 50,49,49], '100x100x100: mode 2')
        call fit%new([101,101,1], 2.)
        call check(fit%loop_lims(2), [0,-50,0, 50,50,0],   '101x101x1: mode 2')
        call check(fit%loop_lims(3), [-50,-50,0, 50,50,0], '101x101x1: mode 3')
        call fit%new([120,90,80], 1.)
        call check(fit%loop_lims(2), [0,-45,-40, 60,44,39], '120x90x80: mode 2')

      contains

        subroutine check( lims, expected, label )
            integer,          intent(in) :: lims(3,2), expected(6)
            character(len=*), intent(in) :: label
            call assert_true(all(lims == reshape(expected, [3,2])), label//' loop limits')
        end subroutine check

    end subroutine test_logical_limits

    ! comp_addr_phys maps the logical Friedel half one to one onto the physical array
    ! (ldim(1)/2+1, ldim(2), ldim(3)); comp_addr_logi inverts it; a negative h lands on the storage of
    ! its Friedel mate; the three forms of comp_addr_phys agree
    subroutine test_address_maps()
        integer, parameter :: NBOX = 6
        integer, parameter :: BOXES(3,NBOX) = reshape([100,100,1, 101,101,1, 120,90,1, 40,40,40, 21,21,21, 24,18,16], [3,NBOX])
        type(ftiter)         :: fit
        integer, allocatable :: hits(:,:,:)
        integer :: ib, h, k, l, lims(3,2), phys(3), phys2(2), logi(3), ldim(3)
        logical :: inverse_ok, forms_ok, friedel_ok, in_bounds
        write(*,'(A)') 'test_address_maps'
        do ib = 1,NBOX
            ldim = BOXES(:,ib)
            call fit%new(ldim, 1.)
            lims = fit%loop_lims(2)
            allocate(hits(ldim(1)/2+1, ldim(2), ldim(3)), source=0)
            inverse_ok = .true.
            forms_ok   = .true.
            friedel_ok = .true.
            in_bounds  = .true.
            do h = lims(1,1),lims(1,2)
                do k = lims(2,1),lims(2,2)
                    do l = lims(3,1),lims(3,2)
                        phys = fit%comp_addr_phys(h,k,l)
                        if( any(phys < 1) .or. any(phys > [ldim(1)/2+1, ldim(2), ldim(3)]) )then
                            in_bounds = .false.
                            cycle
                        endif
                        hits(phys(1),phys(2),phys(3)) = hits(phys(1),phys(2),phys(3)) + 1
                        logi = fit%comp_addr_logi(phys(1),phys(2),phys(3))
                        if( any(logi /= [h,k,l]) ) inverse_ok = .false.
                        if( any(fit%comp_addr_phys([h,k,l]) /= phys) ) forms_ok = .false.
                        if( ldim(3) == 1 )then
                            phys2 = fit%comp_addr_phys(h,k)
                            if( any(phys2 /= phys(1:2)) ) forms_ok = .false.
                        endif
                        if( h > 0 .and. -k >= lims(2,1) .and. -k <= lims(2,2) .and. -l >= lims(3,1) .and. -l <= lims(3,2) )then
                            if( any(fit%comp_addr_phys(-h,-k,-l) /= phys) ) friedel_ok = .false.
                        endif
                    end do
                end do
            end do
            call assert_true(in_bounds,         box_tag(ldim)//'every physical address lies in the half-complex array')
            call assert_true(all(hits == 1),    box_tag(ldim)//'the logical half maps one to one onto the physical array')
            call assert_true(inverse_ok,        box_tag(ldim)//'comp_addr_logi inverts comp_addr_phys')
            call assert_true(forms_ok,          box_tag(ldim)//'the array, (h,k,l) and (h,k) forms agree')
            call assert_true(friedel_ok,        box_tag(ldim)//'(-h,-k,-l) addresses the storage of (h,k,l)')
            deallocate(hits)
        end do
    end subroutine test_address_maps

    ! the Fourier scale is the box in Angstroms, ldim(1)*smpd for even and odd boxes: index = box/res
    subroutine test_resolution()
        type(ftiter) :: fit
        write(*,'(A)') 'test_resolution'
        fit = ftiter([100,100,1], 2.)
        call assert_int(50,     fit%get_lfny(1),         '100 px at 2 A: Nyquist index 50')
        call assert_int(25,     fit%get_find(1, 8.),     '100 px at 2 A: 8 A is index 25')
        call assert_real(8.,    fit%get_lp(1, 25),    1.e-6, '100 px at 2 A: index 25 is 8 A')
        call assert_real(0.125, fit%get_spat_freq(1, 25), 1.e-7, '100 px at 2 A: index 25 is 1/8 per A')
        call assert_int(0,      fit%get_lhp(1),          'no high-pass index by default')
        call assert_true(all(fit%get_ldim() == [100,100,1]), 'get_ldim')
        call assert_real(2.,    fit%get_smpd(), 0.,      'get_smpd')
        fit = ftiter([101,101,1], 2.)
        call assert_int(50,     fit%get_lfny(1),         '101 px at 2 A: Nyquist index 50')
        call assert_real(8.08,  fit%get_lp(1, 25),    1.e-5, '101 px at 2 A: index 25 is 202/25 A')
        ! the second dimension of a non-square box uses the scale of the first
        fit = ftiter([120,90,1], 1.)
        call assert_int(15,     fit%get_find(2, 8.),     '120x90 at 1 A: 8 A is index 15 in y too (scale of x)')
        call assert_int(60,     fit%get_lfny(2),         '120x90 at 1 A: Nyquist index of x in y too')
    end subroutine test_resolution

    ! loop_lims(1, lp): h from 0, k and (for a volume) l symmetric, to the index of lp, clamped to
    ! [3, Nyquist]
    subroutine test_lowpass_limits()
        type(ftiter) :: fit
        write(*,'(A)') 'test_lowpass_limits'
        call fit%new([100,100,1], 2.)
        call assert_true(all(fit%loop_lims(1, 8.)   == reshape([0,-25,0, 25,25,0], [3,2])), '2D, 8 A: index 25')
        call assert_true(all(fit%loop_lims(1, 1.)   == reshape([0,-50,0, 50,50,0], [3,2])), '2D, 1 A: clamped to Nyquist')
        call assert_true(all(fit%loop_lims(1, 150.) == reshape([0,-3,0,  3,3,0],   [3,2])), '2D, 150 A: clamped to index 3')
        call fit%new([100,100,100], 2.)
        call assert_true(all(fit%loop_lims(1, 8.)   == reshape([0,-25,-25, 25,25,25], [3,2])), &
            &'3D, 8 A: the third dimension is symmetric like the second')
    end subroutine test_lowpass_limits

    subroutine test_reset()
        type(ftiter) :: fit
        write(*,'(A)') 'test_reset'
        call fit%new([64,64,1], 1.5)
        call fit%reset
        call assert_true(all(fit%get_ldim() == [1,1,1]), 'reset: unit dimensions')
        call assert_real(0., fit%get_smpd(), 0., 'reset: no sampling distance')
        call assert_int(0, fit%get_lfny(1), 'reset: no Nyquist index')
    end subroutine test_reset

end module simple_ftiter_tester

!@descr: tests of the density-domain order mask and its core signal guard
submodule (simple_nu_filter) simple_nu_filter_ordermask_tester
use simple_test_utils, only: assert_true
implicit none
contains

    !> The selection is serial (no parallel region), so it has no thread-count test. Fixture: a 32^3 box with
    !! a core of 4096 dense voxels (margin 0.9-1.1, label 6 A) and 4096 dense voxels around it (margin 0-0.05,
    !! label 20 A), unit half-map products; the expected values are those of the Phase 0 NumPy dry run.
    module subroutine run_all_nu_order_mask_tests()
        integer, parameter :: n = 32, ncore = 4096, nperi = 4096
        real                 :: margin(n,n,n), locres(n,n,n), even(n,n,n), odd(n,n,n), factor(n,n,n), ref(n,n,n)
        logical              :: dense(n,n,n), core(n,n,n), peri(n,n,n)
        type(nu_order_stats) :: st, base
        real(dp)             :: measured
        integer              :: i, j, k, idx
        dense  = .false.
        core   = .false.
        peri   = .false.
        margin = 0.
        locres = 20.
        even   = 1.
        odd    = 1.
        idx    = 0
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    idx = idx + 1
                    if( idx <= ncore )then
                        dense(i,j,k)  = .true.
                        core(i,j,k)   = .true.
                        locres(i,j,k) = 6.
                        margin(i,j,k) = 0.9 + 0.2 * real(idx-1) / real(ncore-1)
                    else if( idx <= ncore + nperi )then
                        dense(i,j,k)  = .true.
                        peri(i,j,k)   = .true.
                        margin(i,j,k) = 0.05 * real(idx-ncore-1) / real(nperi-1)
                    endif
                enddo
            enddo
        enddo
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(.not.st%l_off, 'membrane order mask is active')
        call assert_true(abs(st%fraction-0.50) < 1.e-6, 'the first qualifying fraction is 0.50')
        call assert_true(st%kept_core >= 0.95, 'core keeps at least 95 percent of cross-half signal')
        call assert_true(sum(factor, mask=peri) / real(nperi) <= 0.6, 'membrane fixture mean attenuation meets frozen bound')
        call assert_true(st%l_peri_signal .and. .not.st%l_fine_signal, &
            &'regional diagnostics distinguish a positive denominator from an empty fine region')
        where( peri ) locres = 19.
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st, coarsest_lp=18.)
        call assert_true(.not.st%l_fine_signal, 'fine-region diagnostic uses the actual coarsest bank cutoff')
        where( peri ) locres = 20.
        measured = sum(real(factor,dp)**2, mask=core) / real(ncore,dp)
        call assert_true(abs(measured-real(st%kept_core,dp)) < 1.e-4_dp*abs(measured), &
            &'reported core retention equals measured masked signal')
        ref  = factor
        base = st
        where( .not.dense ) margin = 0.37
        where( .not.dense ) even   = 0.
        where( .not.dense ) odd    = 0.
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(st%fraction == base%fraction .and. all(factor == ref .or. .not.dense), &
            &'constrained and unconstrained twins agree on dense voxels')
        dense = core
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(.not.st%l_off .and. st%kept_core >= 0.95, 'soluble fixture keeps its ordered core')
        call assert_true(.not.st%l_peri_signal, 'an empty peripheral region has no retention ratio')
        dense  = core .or. peri
        margin = 1.
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(st%l_off .and. all(factor == 1.), 'no fraction satisfying the core guard disables attenuation')
        locres = 20.
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(st%l_off .and. all(factor == 1.), 'too few core voxels disables attenuation')
        where( core ) locres = 6.
        where( core ) odd    = -1.
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(st%l_off .and. all(factor == 1.), 'nonpositive core signal disables attenuation')
        where( core ) odd = 1.
        margin = 0.
        call select_nu_order_mask(margin, locres, dense, even, odd, 6., factor, st)
        call assert_true(st%l_off .and. all(factor == 1.), 'no positive threshold disables attenuation')
    end subroutine run_all_nu_order_mask_tests

end submodule simple_nu_filter_ordermask_tester

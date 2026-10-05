!@descr: unit test routines for the accumulator-domain trailing-reconstruction blend
! Covers the image primitives scale_mats and sum_reduce_mats and the weighting contracts of
! blend_trailing_accumulators (simple_commanders_rec_distr), on tiny synthetic accumulators with known answers.
module simple_accum_blend_tester
use simple_test_utils, only: assert_real, assert_true
use simple_image,      only: image
use simple_oris,       only: population_blend_weights
implicit none
private
public :: run_all_accum_blend_tests

contains

    !> Deterministic test of the accumulator-domain trailing-reconstruction
    !! recurrence on tiny synthetic arrays. Exercises the image-level primitives
    !! the blend is built from (scale_mats, sum_reduce_mats) and verifies the
    !! weighting contracts of simple_commanders_rec_distr::blend_trailing_accumulators:
    !!   - a full-mass chain + fractional partials scaled by u/f restores with
    !!     current-map coefficient exactly u (ufrac_trec contract)
    !!   - a chain start (no chain yet) seeds the chain with the current partials
    !!     normalized by 1/f, so it carries full sampling mass and the next
    !!     iteration's effective update equals the realized fraction f
    !!   - the chain-start iteration itself restores the current sample alone:
    !!     the previous map has no weight in it
    !!   - one particle cohort held for k iterations (refine3D_states cohort schedule)
    !!     reaches the cumulative current-map coefficient 1 - (1-u)^k
    subroutine run_all_accum_blend_tests()
        write(*,'(A)') '**** running all trailing-reconstruction blend tests ****'
        call test_trail_rec_blend()
        call test_trail_rec_population()
        call test_trail_rec_cohort()
    end subroutine run_all_accum_blend_tests

    subroutine test_trail_rec_blend()
        real, parameter :: D_FULL = 2.0  ! full-dataset per-voxel sampling density
        real, parameter :: V_CUR  = 3.0  ! current-map voxel value
        real, parameter :: V_PREV = 1.0  ! previous-map voxel value
        real, parameter :: TOL    = 1.e-4
        ! (realized f, applied u) pairs from the review validation matrix;
        ! the expected restored current-map coefficient is u in every case
        real, parameter :: FU_PAIRS(2,5) = reshape([0.1, 0.5,  &
                                                    1.0, 0.5,  &
                                                    0.5, 0.5,  &
                                                    0.1, 0.0,  &
                                                    0.1, 0.1], [2,5])
        integer :: ipair
        real    :: f, u, restored, w_eff
        write(*,'(A)') 'test_trail_rec_blend'
        ! primitives
        call test_scale_mats_primitive()
        ! chain-mode recurrence: full-mass chain, u/f-scaled partials
        do ipair = 1, size(FU_PAIRS, 2)
            f = FU_PAIRS(1,ipair)
            u = FU_PAIRS(2,ipair)
            call run_recurrence(f, u, restored)
            w_eff = (restored - V_PREV) / (V_CUR - V_PREV)
            call assert_real(u, w_eff, TOL, 'effective current weight equals applied fraction u')
        enddo
        ! chain start: seed = (1/f) * partials must carry full mass, so the next
        ! iteration with no override restores with effective weight f
        f = 0.1
        call run_chain_start_then_update(f, restored)
        w_eff = (restored - V_PREV) / (V_CUR - V_PREV)
        call assert_real(f, w_eff, TOL, 'after a chain start the effective update weight equals realized fraction f')
        ! the chain-start iteration ships the current sample's map
        do ipair = 1, size(FU_PAIRS, 2)
            if( FU_PAIRS(1,ipair) < 0.99 ) call run_chain_start_iteration(FU_PAIRS(1,ipair))
        enddo

    contains

        subroutine make_accum( img, rho, map_value, density )
            type(image),       intent(inout) :: img
            real, allocatable, intent(inout) :: rho(:,:,:)
            real,              intent(in)    :: map_value, density
            integer :: shp(3)
            call img%new([8,8,8], 1.0)
            call img%set_ft(.true.)
            call img%set_cmat(cmplx(map_value * density, 0.))
            shp = img%get_array_shape()
            if( allocated(rho) ) deallocate(rho)
            allocate(rho(shp(1),shp(2),shp(3)), source=density)
        end subroutine make_accum

        real function restored_at_origin( img, rho ) result( val )
            type(image), intent(in) :: img
            real,        intent(in) :: rho(:,:,:)
            val = real(img%get_cmat_at(1,1,1)) / rho(1,1,1)
        end function restored_at_origin

        subroutine test_scale_mats_primitive()
            type(image)       :: img
            real, allocatable :: rho(:,:,:)
            call make_accum(img, rho, V_CUR, D_FULL)
            call img%scale_mats(rho, 0.25)
            call assert_real(0.25 * V_CUR * D_FULL, real(img%get_cmat_at(1,1,1)), TOL, 'scale_mats scales cmat')
            call assert_real(0.25 * D_FULL,         rho(1,1,1),                   TOL, 'scale_mats scales rho')
            call img%kill
            deallocate(rho)
        end subroutine test_scale_mats_primitive

        !> chain mode: cur partials carry mass f*D; scale by u/f, decay full-mass
        !! chain by (1-u), sum, restore; also assert the blended mass stays D
        subroutine run_recurrence( f_in, u_in, restored_val )
            real, intent(in)  :: f_in, u_in
            real, intent(out) :: restored_val
            type(image)       :: cur, chain
            real, allocatable :: rho_cur(:,:,:), rho_chain(:,:,:)
            call make_accum(cur,   rho_cur,   V_CUR,  f_in * D_FULL)
            call make_accum(chain, rho_chain, V_PREV, D_FULL)
            call cur%scale_mats(rho_cur, u_in / f_in)
            call chain%scale_mats(rho_chain, 1.0 - u_in)
            call cur%sum_reduce_mats(chain, rho_cur, rho_chain)
            call assert_real(D_FULL, rho_cur(1,1,1), TOL, 'blended chain keeps full sampling mass')
            restored_val = restored_at_origin(cur, rho_cur)
            call cur%kill
            call chain%kill
            deallocate(rho_cur, rho_chain)
        end subroutine run_recurrence

        !> chain-start seeding then one no-override update at realized fraction f
        subroutine run_chain_start_then_update( f_in, restored_val )
            real, intent(in)  :: f_in
            real, intent(out) :: restored_val
            type(image)       :: cur, seed
            real, allocatable :: rho_cur(:,:,:), rho_seed(:,:,:)
            ! iteration 1: fractional partials of the previous map, normalized to full mass
            call make_accum(seed, rho_seed, V_PREV, f_in * D_FULL)
            call seed%scale_mats(rho_seed, 1.0 / f_in)
            call assert_real(D_FULL, rho_seed(1,1,1), TOL, 'chain-start seed carries full sampling mass')
            ! iteration 2: current partials at realized f, chain decayed by (1-f)
            call make_accum(cur, rho_cur, V_CUR, f_in * D_FULL)
            call seed%scale_mats(rho_seed, 1.0 - f_in)
            call cur%sum_reduce_mats(seed, rho_cur, rho_seed)
            restored_val = restored_at_origin(cur, rho_cur)
            call cur%kill
            call seed%kill
            deallocate(rho_cur, rho_seed)
        end subroutine run_chain_start_then_update

        !> The chain-start iteration as blend_trailing_accumulators runs it: the
        !! current partials (map V_CUR, mass f*D) are scaled by 1/f and written as
        !! the chain, then scaled back by f and restored. The written chain has
        !! full mass and the current map; the restored map and its density are
        !! those of the current partials alone, so the previous map (V_PREV) has
        !! no weight in this iteration whatever f is.
        subroutine run_chain_start_iteration( f_in )
            real, intent(in)  :: f_in
            type(image)       :: cur, ref
            real, allocatable :: rho_cur(:,:,:), rho_ref(:,:,:)
            call make_accum(cur, rho_cur, V_CUR, f_in * D_FULL)
            call make_accum(ref, rho_ref, V_CUR, f_in * D_FULL)
            call cur%scale_mats(rho_cur, 1.0 / f_in)
            call assert_real(D_FULL, rho_cur(1,1,1), TOL, 'chain start writes a full-mass chain')
            call assert_real(V_CUR, restored_at_origin(cur, rho_cur), TOL, 'chain start writes the current map as the chain')
            call cur%scale_mats(rho_cur, f_in)
            call assert_real(rho_ref(1,1,1), rho_cur(1,1,1), TOL, &
                &'chain-start iteration restores with the current sample''s own density')
            call assert_real(real(ref%get_cmat_at(1,1,1)), real(cur%get_cmat_at(1,1,1)), TOL, &
                &'chain-start iteration restores the current sample''s own sums')
            call assert_real(1.0, (restored_at_origin(cur, rho_cur) - V_PREV) / (V_CUR - V_PREV), TOL, &
                &'chain-start iteration ships the current sample''s map (previous map weight 0)')
            call cur%kill
            call ref%kill
            deallocate(rho_cur, rho_ref)
        end subroutine run_chain_start_iteration

    end subroutine test_trail_rec_blend

    !> Population rule on the accumulator primitives: one unit of sampling density per
    !! particle, a chain representing M particles of the previous map and current partials
    !! of n particles of the current map. For each population change the blended density
    !! equals N, the population now represented, and the restored current-map coefficient
    !! is the applied fraction u; without a population change the weights are the former
    !! ones (1 - u on the chain). The former rule is kept as a control: it loses density
    !! when first-time particles join. Expected values follow from s = u/f, w = (1-u)*N/M.
    subroutine test_trail_rec_population()
        real, parameter :: V_CUR  = 3.0
        real, parameter :: V_PREV = 1.0
        real, parameter :: TOL    = 1.e-4
        ! N, n, M: first-time rows, deactivation, re-activation, no change, n = N
        integer, parameter :: CASES(3,5) = reshape([130, 40, 100,  &
                                                     80, 20, 100,  &
                                                    120, 10, 100,  &
                                                    100, 25, 100,  &
                                                     60, 60,  90], [3,5])
        real    :: s, w, mnew, dens, restored, f, u
        integer :: icase
        write(*,'(A)') 'test_trail_rec_population'
        do icase = 1, size(CASES, 2)
            f = real(CASES(2,icase)) / real(CASES(1,icase))
            ! default u = f, then a ufrac_trec-like override
            u = f
            call population_blend_weights(CASES(1,icase), CASES(2,icase), real(CASES(3,icase)), s, w, mnew)
            call blend(s, w, CASES(2,icase), CASES(3,icase), dens, restored)
            call assert_real(real(CASES(1,icase)), dens, TOL, 'population rule: blended density = N')
            call assert_real(real(CASES(1,icase)), mnew, TOL, 'population rule: recorded M = N')
            call assert_real(u, (restored - V_PREV) / (V_CUR - V_PREV), TOL, 'population rule: current-map coefficient = f')
            if( CASES(1,icase) == CASES(3,icase) ) call assert_real(1.0 - f, w, TOL, 'no population change: chain weight 1 - f')
            if( CASES(2,icase) < CASES(1,icase) )then
                u = 0.5 * f
                call population_blend_weights(CASES(1,icase), CASES(2,icase), real(CASES(3,icase)), s, w, mnew, ufrac=u)
                call blend(s, w, CASES(2,icase), CASES(3,icase), dens, restored)
                call assert_real(real(CASES(1,icase)), dens, TOL, 'population rule with ufrac: blended density = N')
                call assert_real(u, (restored - V_PREV) / (V_CUR - V_PREV), TOL, 'population rule: current-map coefficient = u')
            endif
        enddo
        ! control: the former weights (s = 1, chain 1 - f) on the first case (N 130, n 40, M 100:
        ! 30 first-time particles)
        call blend(1.0, 1.0 - 40./130., 40, 100, dens, restored)
        call assert_real(40. + (1. - 40./130.) * 100., dens, TOL, 'former rule: density n + (1 - f) M')
        call assert_true(dens < 130. - 1., 'former rule loses density when first-time particles join')

    contains

        subroutine blend( s_in, w_in, n_in, m_in, dens_out, restored_out )
            real,    intent(in)  :: s_in, w_in
            integer, intent(in)  :: n_in, m_in
            real,    intent(out) :: dens_out, restored_out
            type(image)       :: cur, chain
            real, allocatable :: rho_cur(:,:,:), rho_chain(:,:,:)
            call make_unit_accum(cur,   rho_cur,   V_CUR,  real(n_in))
            call make_unit_accum(chain, rho_chain, V_PREV, real(m_in))
            call cur%scale_mats(rho_cur, s_in)
            call chain%scale_mats(rho_chain, w_in)
            call cur%sum_reduce_mats(chain, rho_cur, rho_chain)
            dens_out     = rho_cur(1,1,1)
            restored_out = real(cur%get_cmat_at(1,1,1)) / rho_cur(1,1,1)
            call cur%kill
            call chain%kill
            deallocate(rho_cur, rho_chain)
        end subroutine blend

        subroutine make_unit_accum( img, rho, map_value, density )
            type(image),       intent(inout) :: img
            real, allocatable, intent(inout) :: rho(:,:,:)
            real,              intent(in)    :: map_value, density
            integer :: shp(3)
            call img%new([8,8,8], 1.0)
            call img%set_ft(.true.)
            call img%set_cmat(cmplx(map_value * density, 0.))
            shp = img%get_array_shape()
            if( allocated(rho) ) deallocate(rho)
            allocate(rho(shp(1),shp(2),shp(3)), source=density)
        end subroutine make_unit_accum

    end subroutine test_trail_rec_population

    !> One cohort (refine3D_states cohort_sampling) reconstructed for k iterations of a frequency block:
    !! every iteration its partials (cohort map, mass f*D) are scaled by u/f and the full-mass chain by
    !! 1-u. The even and odd chains and their density must follow the closed form: map coefficient of the
    !! cohort 1 - (1-u)^k, earlier maps fading as (1-u)^k, density D throughout.
    subroutine test_trail_rec_cohort()
        real,    parameter :: D_FULL = 2.0
        real,    parameter :: V_PREV(2) = [1.0, 2.0]  ! even, odd previous maps
        real,    parameter :: V_COH(2)  = [3.0, 5.0]  ! even, odd cohort maps
        real,    parameter :: TOL    = 1.e-4
        integer, parameter :: NITS   = 4
        real,    parameter :: F_COH  = 0.25           ! realized fraction of the cohort; u = f
        type(image)       :: cur, chain
        real, allocatable :: rho_cur(:,:,:), rho_chain(:,:,:)
        real    :: coef
        integer :: ieo, k
        write(*,'(A)') 'test_trail_rec_cohort'
        do ieo = 1, 2
            call make_unit_density_accum(chain, rho_chain, V_PREV(ieo), D_FULL)
            do k = 1, NITS
                call make_unit_density_accum(cur, rho_cur, V_COH(ieo), F_COH * D_FULL)
                call cur%scale_mats(rho_cur, 1.0)                ! u/f = 1 for u = f
                call chain%scale_mats(rho_chain, 1.0 - F_COH)
                call cur%sum_reduce_mats(chain, rho_cur, rho_chain)
                ! the blend is the next chain
                call chain%copy(cur)
                rho_chain = rho_cur
                call cur%kill
                coef = 1.0 - (1.0 - F_COH)**k
                call assert_real(D_FULL, rho_chain(1,1,1), TOL, 'cohort chain keeps full sampling mass')
                call assert_real(V_PREV(ieo) + coef * (V_COH(ieo) - V_PREV(ieo)), &
                    &real(chain%get_cmat_at(1,1,1)) / rho_chain(1,1,1), TOL, &
                    &'cohort held k iterations: map coefficient 1 - (1-u)^k')
            end do
            call chain%kill
            deallocate(rho_cur, rho_chain)
        end do

    contains

        subroutine make_unit_density_accum( img, rho, map_value, density )
            type(image),       intent(inout) :: img
            real, allocatable, intent(inout) :: rho(:,:,:)
            real,              intent(in)    :: map_value, density
            integer :: shp(3)
            call img%new([8,8,8], 1.0)
            call img%set_ft(.true.)
            call img%set_cmat(cmplx(map_value * density, 0.))
            shp = img%get_array_shape()
            if( allocated(rho) ) deallocate(rho)
            allocate(rho(shp(1),shp(2),shp(3)), source=density)
        end subroutine make_unit_density_accum

    end subroutine test_trail_rec_cohort

end module simple_accum_blend_tester

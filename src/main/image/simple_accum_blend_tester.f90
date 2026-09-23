!@descr: unit test routines for the accumulator-domain trailing-reconstruction blend
! Covers the image primitives scale_mats and sum_reduce_mats and the weighting contracts of
! blend_trailing_accumulators (simple_commanders_rec_distr), on tiny synthetic accumulators with known answers.
module simple_accum_blend_tester
use simple_test_utils, only: assert_real
use simple_image,      only: image
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
    !!   - a bootstrap seed normalized by 1/f carries full sampling mass, so the
    !!     next iteration's effective update equals the realized fraction f
    subroutine run_all_accum_blend_tests()
        write(*,'(A)') '**** running all trailing-reconstruction blend tests ****'
        call test_trail_rec_blend()
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
        ! bootstrap: seed = (1/f) * partials must carry full mass, so the next
        ! iteration with no override restores with effective weight f
        f = 0.1
        call run_bootstrap_then_update(f, restored)
        w_eff = (restored - V_PREV) / (V_CUR - V_PREV)
        call assert_real(f, w_eff, TOL, 'post-bootstrap effective update weight equals realized fraction f')

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

        !> bootstrap seeding then one no-override update at realized fraction f
        subroutine run_bootstrap_then_update( f_in, restored_val )
            real, intent(in)  :: f_in
            real, intent(out) :: restored_val
            type(image)       :: cur, seed
            real, allocatable :: rho_cur(:,:,:), rho_seed(:,:,:)
            ! iteration 1: fractional partials of the previous map, normalized to full mass
            call make_accum(seed, rho_seed, V_PREV, f_in * D_FULL)
            call seed%scale_mats(rho_seed, 1.0 / f_in)
            call assert_real(D_FULL, rho_seed(1,1,1), TOL, 'bootstrap seed carries full sampling mass')
            ! iteration 2: current partials at realized f, chain decayed by (1-f)
            call make_accum(cur, rho_cur, V_CUR, f_in * D_FULL)
            call seed%scale_mats(rho_seed, 1.0 - f_in)
            call cur%sum_reduce_mats(seed, rho_cur, rho_seed)
            restored_val = restored_at_origin(cur, rho_cur)
            call cur%kill
            call seed%kill
            deallocate(rho_cur, rho_seed)
        end subroutine run_bootstrap_then_update

    end subroutine test_trail_rec_blend

end module simple_accum_blend_tester

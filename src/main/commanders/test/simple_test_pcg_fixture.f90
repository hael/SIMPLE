!@descr: PCG test fixture shared by testers: the relative differences of two open raw accumulators (B and D)
module simple_test_pcg_fixture
use simple_core_module_api,   only: dp
use simple_error,             only: simple_exception
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: compare_raw
#include "simple_local_flags.inc"

contains

    !> relative L2 differences of two open raw accumulators over the rows the fold reads (h >= 0 and the
    !! wrapped Nyquist mate at the lower h bound), the rows a persisted accumulator carries
    subroutine compare_raw( op, other, b_relerr, d_relerr )
        class(reconstructor_pcg), intent(in)  :: op, other
        real,                    intent(out) :: b_relerr, d_relerr
        complex, allocatable :: b(:,:,:), b_other(:,:,:)
        real,    allocatable :: d(:,:,:), d_other(:,:,:)
        real(dp) :: bnum, bden, dnum, dden
        integer  :: hlo
        call op%get_raw_accum(b, d)
        call other%get_raw_accum(b_other, d_other)
        if( any(lbound(b) /= lbound(b_other)) .or. any(ubound(b) /= ubound(b_other)) )then
            THROW_HARD('raw PCG comparison geometry mismatch')
        endif
        hlo  = lbound(b,1)
        bnum = sum(abs(cmplx(b(hlo,:,:) - b_other(hlo,:,:), kind=dp))**2) + &
              &sum(abs(cmplx(b(0:,:,:)  - b_other(0:,:,:),  kind=dp))**2)
        bden = sum(abs(cmplx(b_other(hlo,:,:), kind=dp))**2) + sum(abs(cmplx(b_other(0:,:,:), kind=dp))**2)
        dnum = sum(real(d(hlo,:,:) - d_other(hlo,:,:), dp)**2) + sum(real(d(0:,:,:) - d_other(0:,:,:), dp)**2)
        dden = sum(real(d_other(hlo,:,:), dp)**2) + sum(real(d_other(0:,:,:), dp)**2)
        b_relerr = real(sqrt(bnum / max(1.0_dp, bden)))
        d_relerr = real(sqrt(dnum / max(1.0_dp, dden)))
    end subroutine compare_raw

end module simple_test_pcg_fixture

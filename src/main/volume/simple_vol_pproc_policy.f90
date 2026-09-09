!@descr: per-state mask artifact compatibility check shared by volume assembly, postprocess and the abinitio final rec
module simple_vol_pproc_policy
use simple_core_module_api
implicit none

private

public :: state_mask_is_compatible

contains

    subroutine state_mask_is_compatible( mskfile_state, box, smpd, exists, compatible )
        class(string), intent(in)  :: mskfile_state
        integer,       intent(in)  :: box
        real,          intent(in)  :: smpd
        logical,       intent(out) :: exists
        logical,       intent(out) :: compatible
        real    :: smpd_mask
        integer :: ldim_mask(3), nptcls_mask
        exists      = file_exists(mskfile_state)
        compatible  = .false.
        if( .not. exists ) return
        call find_ldim_nptcls(mskfile_state, ldim_mask, nptcls_mask)
        smpd_mask  = find_img_smpd(mskfile_state)
        compatible = all(ldim_mask == box) .and.  abs(smpd_mask - smpd) <= 1.e-6
    end subroutine state_mask_is_compatible

end module simple_vol_pproc_policy

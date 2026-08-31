!@descr: workflow-neutral frequency-stage planning for refine3D wrappers
module simple_refine3D_stage_plan
use simple_error,     only: simple_exception
use simple_estimate_ssnr, only: lpstages_setlims
use simple_math_ft,   only: calc_fourier_index, calc_lowpass_lim
use simple_type_defs, only: lp_crop_inf
implicit none
#include "simple_local_flags.inc"

public :: refine3D_stage_plan_entry, plan_refine3D_frequency_stages
private

type :: refine3D_stage_plan_entry
    integer :: stage      = 0
    integer :: first_iter = 0
    integer :: last_iter  = 0
    integer :: nits       = 0
    type(lp_crop_inf) :: lpinfo
end type refine3D_stage_plan_entry

contains

    subroutine plan_refine3D_frequency_stages( box, smpd, lpstart, lpstop, total_nits, block_nits, &
        &first_iter, master_lpinfo, stages )
        integer,                           intent(in)  :: box, total_nits, block_nits, first_iter
        real,                              intent(in)  :: smpd, lpstart, lpstop
        type(lp_crop_inf),                 intent(in)  :: master_lpinfo
        type(refine3D_stage_plan_entry), allocatable, intent(out) :: stages(:)
        integer :: find_start, find_stop, istage, nstages, nits_remaining
        real    :: rfind, rfind_incr
        type(lp_crop_inf) :: stage_lpinfo(1)
        if( box < 8 ) THROW_HARD('frequency-stage planning requires box >= 8')
        if( smpd <= 0. ) THROW_HARD('frequency-stage planning requires positive sampling')
        if( total_nits < 1 ) THROW_HARD('frequency-stage planning requires at least one iteration')
        if( block_nits < 1 ) THROW_HARD('frequency-stage planning requires a positive block length')
        if( first_iter < 1 ) THROW_HARD('frequency-stage planning requires first_iter >= 1')
        find_start = max(5,       calc_fourier_index(lpstart, box, smpd))
        find_stop  = min(box/2-2, calc_fourier_index(lpstop,  box, smpd))
        if( find_stop < find_start )then
            THROW_HARD('frequency-stage planning requires lpstop to be at least as high-resolution as lpstart')
        endif
        nstages    = ceiling(real(total_nits) / real(block_nits))
        if( nstages > 2 )then
            rfind_incr = real(find_stop - find_start) / real(nstages - 2)
        else
            rfind_incr = real(find_stop - find_start)
        endif
        allocate(stages(nstages))
        rfind = real(find_start) - rfind_incr
        nits_remaining = total_nits
        do istage = 1,nstages
            rfind = rfind + rfind_incr
            stages(istage)%stage      = istage
            stages(istage)%first_iter = first_iter + (istage - 1) * block_nits
            stages(istage)%nits       = min(block_nits, nits_remaining)
            stages(istage)%last_iter  = stages(istage)%first_iter + stages(istage)%nits - 1
            stages(istage)%lpinfo     = master_lpinfo
            stages(istage)%lpinfo%lp  = max(calc_lowpass_lim(nint(rfind), box, smpd), lpstop)
            call lpstages_setlims(box, 1, smpd, stages(istage)%lpinfo%lp, &
                &stages(istage)%lpinfo%lp, stage_lpinfo)
            stages(istage)%lpinfo%box_crop    = stage_lpinfo(1)%box_crop
            stages(istage)%lpinfo%smpd_crop   = stage_lpinfo(1)%smpd_crop
            stages(istage)%lpinfo%scale       = stage_lpinfo(1)%scale
            stages(istage)%lpinfo%trslim      = stage_lpinfo(1)%trslim
            stages(istage)%lpinfo%l_autoscale = stage_lpinfo(1)%l_autoscale
            stages(istage)%lpinfo%l_lpset = .true.
            nits_remaining = nits_remaining - stages(istage)%nits
        enddo
    end subroutine plan_refine3D_frequency_stages

end module simple_refine3D_stage_plan

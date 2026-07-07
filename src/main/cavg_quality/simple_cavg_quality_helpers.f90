!@descr: helper utilities for class-average quality workflows
module simple_cavg_quality_helpers
use simple_string,             only: string
use simple_cavg_quality_types, only: CAVG_REJECT_REASON_NONE, CAVG_REJECT_REASON_POP_NONPOS, CAVG_REJECT_REASON_POP_LOWFRAC, &
                                     CAVG_REJECT_REASON_BAD_PIXELS, CAVG_REJECT_REASON_NO_COMPONENT, CAVG_REJECT_REASON_MASK_GEOMETRY, &
                                     CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW
implicit none
private

public :: cavg_rejection_reason_string

contains

    function cavg_rejection_reason_string( reason_code ) result( reason_text )
        integer, intent(in) :: reason_code
        type(string)        :: reason_text
        select case( reason_code )
            case( CAVG_REJECT_REASON_NONE )
                reason_text = 'none'
            case( CAVG_REJECT_REASON_POP_NONPOS )
                reason_text = 'pop_nonpos'
            case( CAVG_REJECT_REASON_POP_LOWFRAC )
                reason_text = 'pop_lowfrac'
            case( CAVG_REJECT_REASON_BAD_PIXELS )
                reason_text = 'bad_pixels'
            case( CAVG_REJECT_REASON_NO_COMPONENT )
                reason_text = 'no_component'
            case( CAVG_REJECT_REASON_MASK_GEOMETRY )
                reason_text = 'mask_geometry'
            case( CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW )
                reason_text = 'bp_center_edge_low'
            case default
                reason_text = 'unknown_reason_code'
        end select
    end function cavg_rejection_reason_string

end module simple_cavg_quality_helpers
!@descr: helper utilities for class-average quality workflows
module simple_cavg_quality_helpers
use simple_string,             only: string
use simple_cavg_quality_types, only: CAVG_REJECT_REASON_NONE, CAVG_REJECT_REASON_POP, &
                                     CAVG_REJECT_REASON_BAD_PIXELS, CAVG_REJECT_REASON_NO_COMPONENT, CAVG_REJECT_REASON_MASK_GEOMETRY, &
                                     CAVG_REJECT_REASON_MODEL
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
            case( CAVG_REJECT_REASON_POP )
                reason_text = 'low population'
            case( CAVG_REJECT_REASON_BAD_PIXELS )
                reason_text = 'bad pixels'
            case( CAVG_REJECT_REASON_NO_COMPONENT )
                reason_text = 'no mask component'
            case( CAVG_REJECT_REASON_MASK_GEOMETRY )
                reason_text = 'mask geometry issues'
            case( CAVG_REJECT_REASON_MODEL )
                reason_text = 'rejected by learned model'
            case default
                reason_text = 'unknown_reason_code'
        end select
    end function cavg_rejection_reason_string

end module simple_cavg_quality_helpers
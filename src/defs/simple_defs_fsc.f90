!@descr: file prefixes for phase-randomized FSC diagnostics
module simple_defs_fsc
implicit none
private

public :: FSC_UNMASKED_FBODY, FSC_MASKED_FBODY, FSC_RANDOMIZED_FBODY

character(len=*), parameter :: FSC_UNMASKED_FBODY   = 'fscu_state'
character(len=*), parameter :: FSC_MASKED_FBODY     = 'fsct_state'
character(len=*), parameter :: FSC_RANDOMIZED_FBODY = 'fscn_state'

end module simple_defs_fsc

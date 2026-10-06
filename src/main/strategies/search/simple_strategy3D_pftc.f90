!@descr: abstract parent of the polar (PFTC) 3D search strategies, owning their polar search object
! The representation-neutral strategy3D holds the particle spec and the four deferred methods
! only; the polar strategies (greedy, greedy_inpl, greedy_sub, shc, eval, prob) extend this
! type, which adds the strategy3D_srch object s they drive.
! A Cartesian strategy extends strategy3D directly and carries no polar state (plan section 6.3).
module simple_strategy3D_pftc
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_srch
implicit none

public :: strategy3D_pftc
private

type, abstract, extends(strategy3D) :: strategy3D_pftc
    type(strategy3D_srch) :: s
end type strategy3D_pftc

end module simple_strategy3D_pftc

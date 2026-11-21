!@descr: file-name constants of the order reference mask (automsk=order); outside simple_defs_fname, which the
! core module API re-exports, so that adding a name here recompiles only its users
module simple_defs_automask
implicit none
private
public :: ORDER_ENVMASK_FBODY, ORDER_ENVMASK_SUFFIX
character(len=*), parameter :: ORDER_ENVMASK_FBODY  = 'order_envmask3D_state' !< refinement artifact, per state
character(len=*), parameter :: ORDER_ENVMASK_SUFFIX = '_order_envmask'        !< nu_filt3D output beside its average
end module simple_defs_automask

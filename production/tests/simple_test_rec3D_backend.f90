program simple_test_rec3D_backend
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_rec3D_strategy, only: rec3D_strategy, rec3D_inmem_strategy, rec3D_pcg_inmem_strategy, &
    &create_rec3D_strategy, rec3D_backend_id, rec3D_backend_is_wired, &
    &REC3D_BACKEND_INVALID, REC3D_BACKEND_GRIDDING, REC3D_BACKEND_PCG
implicit none

type(parameters) :: params
type(cmdline)    :: cline
class(rec3D_strategy), allocatable :: strategy

if( trim(params%rec_backend) /= 'gridding' )then
    error stop 'rec_backend default must preserve gridding'
endif
if( params%maxits_pcg /= 2 )then
    error stop 'maxits_pcg default must be 2'
endif
if( rec3D_backend_id('gridding') /= REC3D_BACKEND_GRIDDING )then
    error stop 'gridding backend resolution failed'
endif
if( .not. rec3D_backend_is_wired(REC3D_BACKEND_GRIDDING) )then
    error stop 'gridding backend must remain wired'
endif
if( rec3D_backend_id('pcg') /= REC3D_BACKEND_PCG )then
    error stop 'pcg backend resolution failed'
endif
if( .not. rec3D_backend_is_wired(REC3D_BACKEND_PCG) )then
    error stop 'pcg backend must resolve to the shared-memory production strategy'
endif
if( rec3D_backend_id('invalid') /= REC3D_BACKEND_INVALID )then
    error stop 'invalid backend resolution failed'
endif

strategy = create_rec3D_strategy(cline)
select type(strategy)
    type is(rec3D_inmem_strategy)
    class default
        error stop 'default backend must construct the gridding in-memory strategy'
end select
deallocate(strategy)

call cline%set('rec_backend', 'pcg')
strategy = create_rec3D_strategy(cline)
select type(strategy)
    type is(rec3D_pcg_inmem_strategy)
    class default
        error stop 'pcg backend must construct the PCG in-memory strategy'
end select
deallocate(strategy)
call cline%kill

write(*,'(A)') 'REC3D BACKEND SELECTOR TEST PASSED'

end program simple_test_rec3D_backend

program simple_test_clineparse
use simple_cmdline, only: cmdline
use simple_chash,   only: chash
use simple_defs     ! singleton
implicit none
type(cmdline) :: cline
type(chash)   :: hash

call cline%set('prg'       , 'prime2D'   )
call cline%set('stk'       , 'ptcls.mrc' )
call cline%set('smpd'      , 1.77        )
call cline%set('msk'       , 120.        )
call cline%set('ncls'      , 400.        )
call cline%set('split_mode', 'even'      )
call cline%set('nthr'      , 16.         )
call cline%set('nparts'    , 20.         )

call cline%gen_job_descr(hash)
call hash%print_key_val_pairs

call cline%gen_job_descr(hash, prg='prime2D_init')
call hash%print_key_val_pairs


end program simple_test_clineparse
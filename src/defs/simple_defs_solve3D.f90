!@descr: file names of the solve3D workflows (solve3D_cavgs, solve3D_chunks)
module simple_defs_solve3D
implicit none

! solve3D_cavgs
character(len=*), parameter :: CAVGS_SIGMA2_STATE_FNAME = 'solve3D_cavgs_sigma2_state.bin'

! solve3D_chunks
character(len=*), parameter :: CHUNKS_SOURCE_PROJ = 'source.simple'
character(len=*), parameter :: CHUNKS_BASE_PROJ   = 'base.simple'

end module simple_defs_solve3D

program simple_test_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_picker_utils, only: picker_utils
use simple_image,        only: image
implicit none

character(len=*), parameter :: micname = '/home/elmlundho/cache/relion_tut/20170629_00021_frameImage_intg.mrc'
character(len=*), parameter :: boxrefs = '/home/elmlundho/cache/relion_tut/boxrefs.mrc'
real,             parameter :: SMPD    = 0.885, MOLDIAM = 180.
type(image)                 :: micimg
type(image),    allocatable :: refs(:)
type(picker_utils)          :: putils
integer                     :: ldim(3), ifoo, nthr, nptcls, ldim_refs(3), nrefs, iref
character(len=LONGSTRLEN)   :: boxname_out, cwd

!$ nthr = omp_get_max_threads()
!$ call omp_set_num_threads(nthr)
nthr_glob = nthr
call simple_getcwd(cwd)
allocate(CWD_GLOB, source=trim(cwd))

call find_ldim_nptcls(micname, ldim, ifoo)
call find_ldim_nptcls(boxrefs, ldim_refs, nrefs)
ldim_refs(3) = 1
allocate(refs(nrefs))
do iref = 1,nrefs
    call refs(iref)%new(ldim_refs, SMPD)
    call refs(iref)%read(boxrefs, iref)
    ! call refs(iref)%write('foo.mrc', iref)
end do
call micimg%new(ldim, SMPD)
call micimg%read(micname)
call putils%new(micimg, SMPD, MOLDIAM, 'black', 'test')
call putils%set_refs(refs, 200.)

! call putils%exec_gaupicker(micimg, SMPD, MOLDIAM, 'black', micname, boxname_out, nptcls)

! print *, 'boxname_out ', trim(boxname_out)
! print *, 'nptcls      ', nptcls

end program simple_test_picker_utils

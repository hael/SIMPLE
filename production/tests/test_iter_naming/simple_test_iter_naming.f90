program test_iter_naming
include 'simple_lib.f08'
implicit none
character(len=*), parameter :: STARTVOL = 'recvol_state01_SIM.mrc'
character(len=*), parameter :: RECVOL   = 'recvol_state01.mrc'
character(len=*), parameter :: ATOMS    = 'recvol_state01_ATMS.pdb'
character(len=STDLEN)       :: vol_iter, atms_iter, sim_iter, fbody
integer :: i
do i = 1, 5
    fbody = get_fbody(basename(RECVOL), 'mrc')
    vol_iter  = trim(fbody)//'_iter'//int2str_pad(i,3)//'.mrc'
    atms_iter = trim(fbody)//'_iter'//int2str_pad(i,3)//'_ATMS.pdb'
    sim_iter  = trim(fbody)//'_iter'//int2str_pad(i,3)//'_SIM.mrc'
    print *, trim(vol_iter)
    print *, trim(atms_iter)
    print *, trim(sim_iter)
end do
end program test_iter_naming

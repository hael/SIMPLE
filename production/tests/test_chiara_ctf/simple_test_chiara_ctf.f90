program simple_test_chiara_ctf
use simple_ctf,   only: ctf
use simple_image, only: image
implicit none

real,    parameter :: smpd=1.2, kV=300, Cs=2.7, amp_contr=0.1, dfx=1.54, dfy =1.72, angast=30.
integer, parameter :: box=256

type(image) :: ctf_image, img4viz
type(ctf)   :: tfun

call ctf_image%new([box,box,1], smpd)
call img4viz%new([box,box,1], smpd)
tfun = ctf(smpd, kV, Cs, amp_contr)
call tfun%ctf2img(ctf_image, dfx, dfy, angast)
call ctf_image%ft2img('real', img4viz)
call img4viz%vis

end program simple_test_chiara_ctf

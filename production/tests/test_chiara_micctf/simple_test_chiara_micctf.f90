program simple_test_chiara_micctf
  include 'simple_lib.f08'
  use simple_oris
  use gnufor2
  use simple_image
  use simple_ctf
  use simple_stackops
  use simple_math
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline
type(image) :: img, img_out, mic
type(ctf)   :: tfun
type(ctfparams) :: ctfparms
real, allocatable :: rmat(:,:,:)
integer :: i, winsz, m
logical ::yes_no
real :: matrix(7,7,1), thresh
integer :: vec(10), ldim(3), nptcls
real :: vec_r(10), smpd
real :: xhist_i(5), yhist_i(5), xhist_r1(5), yhist_r1(5), xhist(7), yhist(7)

ctfparms%smpd   = 1.32
ctfparms%cs     = 2.7
ctfparms%kv     = 300
ctfparms%fraca  = 0.1
ctfparms%dfx    = 2.62365627
ctfparms%dfy    = 2.60851598
ctfparms%angast = -29.8392296
call find_ldim_nptcls('/home/lenovoc30/Desktop/ANTERGOS/forctf/0001_forctf.mrc', ldim, nptcls)
call mic%new(ldim, ctfparms%smpd)
call mic%read('/home/lenovoc30/Desktop/ANTERGOS/forctf/0001_forctf.mrc')
call mic%clip_inplace([3600,3600,1])
call mic%write('mic_clip.mrc')
call mic%bp(100., 0.)
call mic%write('mic_hp.mrc')
call mic%fft
tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)
call tfun%apply_serial(mic,'ctf',ctfparms)
call mic%ifft
call mic%write('mic_mul.mrc')
! watch the results in ~/Desktop/ANTERGOS/forctf/test
end program simple_test_chiara_micctf

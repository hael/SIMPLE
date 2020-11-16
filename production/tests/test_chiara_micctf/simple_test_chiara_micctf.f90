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
  use simple_tvfilter,    only: tvfilter

type(image)     :: mic
type(ctf)       :: tfun
type(ctfparams) :: ctfparms
type(tvfilter)  :: tvf
integer         :: ldim(3), nptcls
real :: lambda

ctfparms%smpd   = 1.32
ctfparms%cs     = 2.7
ctfparms%kv     = 300
ctfparms%fraca  = 0.1
ctfparms%dfx    = 2.62365627
ctfparms%dfy    = 2.60851598
ctfparms%angast = -29.8392296
call find_ldim_nptcls('/home/chiara/Desktop/Chiara/ANTERGOS/forctf/0001_forctf.mrc', ldim, nptcls)
call mic%new(ldim, ctfparms%smpd)
call mic%read('/home/chiara/Desktop/Chiara/ANTERGOS/forctf/0001_forctf.mrc')
! I want to compare WienerLikeRestoration and TVfiltering fo the running time.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WIENER-LIKE RESTORATION
! call mic%fft
! tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)
! call tfun%wienerlike_restoration(mic, ctfparms)
! !call mic%bp(500., 0.) !hp
! call mic%ifft
! call mic%write('mic_wienerlike_restored.mrc')
! watch the results in /home/chiara/Desktop/Chiara/ANTERGOS/forctf/test or /test_2/
! other results in /home/chiara/Desktop/Chiara/WienerLikeRestoration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TV FILTERING
lambda = 100.
call tvf%new()
call tvf%apply_filter(mic, lambda)
call mic%write('mic_tv_filtered.mrc')
end program simple_test_chiara_micctf

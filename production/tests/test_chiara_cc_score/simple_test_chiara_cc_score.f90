program try_cc_score
  include 'simple_lib.f08'
  use simple_image
  use simple_oris
  use gnufor2
  use simple_ctf,  only : ctf
  implicit none
  type(oris)        :: os
  integer           :: nl, i
  real, allocatable :: cc90(:), indices(:), ctf_estimatecc(:)
  nl = nlines('/home/lenovoc30/Desktop/ANTERGOS/forctf/test/1_ctf_estimate/info.txt')
  write(logfhandle,*) "nl = ", nl
  call os%new(nl, is_ptcl=.false.)
  call os%read('info.txt')
  allocate(cc90(nl), ctf_estimatecc(nl), indices(nl), source = 0.)
  cc90           = os%get_all('cc90')
  ctf_estimatecc = os%get_all('ctf_estimatecc')
  do i = 1, nl
    indices(i) = real(i)
  enddo
  call plot(indices,cc90)
  call plot(indices,ctf_estimatecc, color1 = 'black')
end program try_cc_score

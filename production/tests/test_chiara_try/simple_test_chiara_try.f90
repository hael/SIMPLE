program simple_test_chiara_try
    !$ use omp_lib
    !$ use omp_lib_kinds
  include 'simple_lib.f08'
  use simple_powerspec_analysis
  use gnufor2
  use simple_ctf
  use simple_micops
  use simple_image
  use simple_stackops
  use simple_math
  use simple_segmentation
  use simple_parameters, only: parameters
  use simple_cmdline,    only: cmdline
  type(image)       :: img, img_cc
  real, allocatable :: rmat(:,:,:), rmat_t(:,:,:)
  integer :: i, ldim(3), nptcls, box, nthr
  type(ctf)       :: tfun
  type(ctfparams) :: ctfparms
  logical :: yes_no
  real :: smpd

    call img%new([512,512,1],1.)
    call img%ellipse([256,256],[20.,20.], 'yes')
    call img%ellipse([352,312],[10.,5.], 'yes')
    call img%ellipse([160,200],[10.,5.], 'yes')
    call img%write('ToyImage2.mrc')
    call img%find_connected_comps(img_cc)
    call img_cc%write('ToyImage2CC.mrc')
    rmat = img_cc%get_rmat()
    do i = 1, int(maxval(rmat))
      yes_no = is_symmetric(img_cc, i)
      write(logfhandle,*) 'CC ', i, ' is symmetric ', yes_no
    enddo
end program simple_test_chiara_try

! ctfparms%smpd   = 1.32
! ctfparms%cs     = 2.7
! ctfparms%kv     = 300
! ctfparms%fraca  = 0.1
! ctfparms%dfx    = 2.62365627
! ctfparms%dfy    = 2.60851598
! ctfparms%angast = -29.8392296
! call find_ldim_nptcls('/home/chiara/Desktop/Chiara/ANTERGOS/forctf/0001_forctf.mrc', ldim, nptcls)
! call mic%new(ldim, ctfparms%smpd)
! tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)


! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])

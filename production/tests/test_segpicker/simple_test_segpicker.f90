program simple_test_segpicker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters
use simple_pickseg
! use simple_image,        only: image
! use simple_tvfilter,     only: tvfilter
! use simple_segmentation, only: otsu_robust_fast, otsu_img
! use simple_binimage,     only: binimage
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: mic       = '/home/elmlundho/cache/NanoX/NP_Minyoung_dry/intgs/FoilHole_6388490_Data_6383078_31_20240207_130359_EER_intg.mrc'
type(parameters), target    :: params
type(pickseg) :: picker
! character(len=*), parameter :: pcontrast = 'black'
! character(len=*), parameter :: extracted = 'extracted.mrc'
! real,             parameter :: SHRINK    = 4.
! real,             parameter :: lp        = 10.
! real,             parameter :: lambda    = 3.
! real,             parameter :: nsig      = 1.5
! logical,          parameter :: L_WRITE   = .false.
! integer, allocatable :: sz(:)
! real,    allocatable :: diams(:)
! integer              :: ldim(3), nframes, ldim_shrink(3), i, sz_min, sz_max, nptcls, boxcoord(2), box
! real                 :: smpd, smpd_shrink
! logical              :: outside
! real                 :: otsu_t, px(3)
! type(image)          :: img_mic, img_win
! type(binimage)       :: img_mic_shrink, img_cc
! type(tvfilter)       :: tvf
! type(stats_struct)   :: sz_stats, diam_stats

params_glob => params
params_glob%pcontrast = 'black'
params_glob%lp        = 10.
params_glob%nsig      = 1.5
call picker%pick( mic )

! call find_ldim_nptcls(mic, ldim, nframes, smpd)
! print *, 'ldim(1)    ', ldim(1)
! print *, 'ldim(2)    ', ldim(2)
! print *, 'ldim(3)    ', ldim(3)
! print *, 'nframes    ', nframes
! print *, 'smpd_found ', smpd

! ! shrink micrograph
! ldim_shrink(1) = round2even(real(ldim(1))/SHRINK)
! ldim_shrink(2) = round2even(real(ldim(2))/SHRINK)
! ldim_shrink(3) = 1
! smpd_shrink    = smpd * SHRINK
! call img_mic%new(ldim, smpd)
! call img_mic%read(mic)
! call img_mic%mul(real(product(ldim))) ! to prevent numerical underflow when performing FFT
! call img_mic%fft
! call img_mic_shrink%new_bimg(ldim_shrink, smpd_shrink)
! call img_mic_shrink%set_ft(.true.)
! call img_mic%clip(img_mic_shrink)
! select case(trim(pcontrast))
!     case('black')
!         ! flip contrast (assuming black particle contrast on input)
!         call img_mic_shrink%mul(-1.)
!     case('white')
!         ! nothing to do
!     case DEFAULT
!         THROW_HARD('uknown pcontrast parameter, use (black|white)')
! end select

! ! low-pass filter micrograph
! call img_mic_shrink%bp(0., lp)
! call img_mic_shrink%ifft
! call img_mic%ifft
! if( L_WRITE ) call img_mic_shrink%write('mic_shrink_lp.mrc')

! ! TV denoising
! call tvf%new()
! call tvf%apply_filter(img_mic_shrink, lambda)
! call tvf%kill
! if( L_WRITE ) call img_mic_shrink%write('mic_shrink_lp_tv.mrc')
! call otsu_img(img_mic_shrink, otsu_t)
! call img_mic_shrink%set_imat
! if( L_WRITE ) call img_mic_shrink%write_bimg('mic_shrink_lp_tv_bin.mrc')
! call img_mic_shrink%erode
! call img_mic_shrink%erode
! if( L_WRITE ) call img_mic_shrink%write_bimg('mic_shrink_lp_tv_bin_erode.mrc')

! ! identify connected components
! call img_mic_shrink%find_ccs(img_cc)
! if( L_WRITE ) call img_cc%write_bimg('mic_shrink_lp_tv_bin_erode_cc.mrc')
! call img_cc%get_nccs(nptcls)
! print *, 'nptcls before elimination: ', nptcls
! ! eliminate connected components that are too large or too small
! sz = img_cc%size_ccs()
! call calc_stats(real(sz), sz_stats)
! print *, 'avg size: ', sz_stats%avg
! print *, 'med size: ', sz_stats%med
! print *, 'sde size: ', sz_stats%sdev
! print *, 'min size: ', sz_stats%minv
! print *, 'max size: ', sz_stats%maxv
! sz_min = nint(sz_stats%avg - sz_stats%sdev)
! sz_max = nint(sz_stats%avg + sz_stats%sdev)
! call img_cc%elim_ccs([sz_min,sz_max])
! call img_cc%get_nccs(nptcls)
! ! print *, 'nptcls after  elimination: ', nptcls
! sz = img_cc%size_ccs()
! call calc_stats(real(sz), sz_stats)
! print *, 'avg size: ', sz_stats%avg
! print *, 'med size: ', sz_stats%med
! print *, 'sde size: ', sz_stats%sdev
! print *, 'min size: ', sz_stats%minv
! print *, 'max size: ', sz_stats%maxv
! allocate(diams(nptcls), source=0.)
! call calc_stats(diams, diam_stats)
! do i = 1, nptcls
!     call img_cc%diameter_cc(i, diams(i))
! end do
! call calc_stats(diams, diam_stats)
! print *, 'avg diam: ', diam_stats%avg
! print *, 'med diam: ', diam_stats%med
! print *, 'sde diam: ', diam_stats%sdev
! print *, 'min diam: ', diam_stats%minv
! print *, 'max diam: ', diam_stats%maxv
! box = find_magic_box(2 * nint(diam_stats%med/smpd))
! call img_win%new([box,box,1], smpd)
! do i = 1, nptcls
!     px       = center_mass_cc(i)
!     boxcoord = nint((real(SHRINK)*px(1:2))-real(box)/2.)
!     call img_mic%window_slim(boxcoord, box, img_win, outside)
!     call img_win%write(extracted, i)
! end do

! contains

!     function center_mass_cc( i_cc ) result( px )
!         integer,           intent(in)    :: i_cc
!         real :: px(3)
!         integer, allocatable :: pos(:,:)
!         integer, allocatable :: imat_cc(:,:,:)
!         imat_cc = int(img_cc%get_rmat())
!         where(imat_cc .ne. i_cc) imat_cc = 0
!         call get_pixel_pos(imat_cc,pos)
!         px(1) = sum(pos(1,:))/real(size(pos,dim = 2))
!         px(2) = sum(pos(2,:))/real(size(pos,dim = 2))
!         px(3) = 1.
!         if(allocated(imat_cc)) deallocate(imat_cc)
!     end function center_mass_cc

end program simple_test_segpicker
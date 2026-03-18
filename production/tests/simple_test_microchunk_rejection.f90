
program main
  use simple_defs,             only: COSMSKHALFWIDTH
  use simple_image,            only: image
  use simple_string,           only: string
  use simple_fileio,           only: basename, swap_suffix, del_file
  use simple_imghead,          only: get_mrcfile_info
  use simple_stack_io,         only: stack_io
  use simple_imgarr_utils,     only: dealloc_imgarr
  use simple_stream_microchunk_utils, only: calc_rejection_score, reject_outliers, reject_auto

  implicit none

  type :: microchunk_cavg
    type(string) :: path
    integer      :: msk
  end type microchunk_cavg
  
  type(image), allocatable :: cavg_imgs(:)
  real,        allocatable :: rmat(:,:,:), stack(:,:,:), scores(:)
  logical,     allocatable :: l_rejected(:)
  type(image)              :: img
  type(string)             :: stkname
  type(stack_io)           :: stkio_r
 ! type(class_scores_t)     :: results
  type(microchunk_cavg)    :: microchunk_cavgs(7)
  integer                  :: icls, ichunk, istk, ldim(3)
  real                     :: smpd, threshold, snr, mskrad
  threshold = 0.25
  microchunk_cavgs(1)%path = 'microchunk_rejection/SLC/SLC_cavgs_1.mrc'
  microchunk_cavgs(1)%msk  = 156
  microchunk_cavgs(2)%path = 'microchunk_rejection/SLC/SLC_cavgs_2.mrc'
  microchunk_cavgs(2)%msk  = 156
  microchunk_cavgs(3)%path = 'microchunk_rejection/SLC/SLC_cavgs_3.mrc'
  microchunk_cavgs(3)%msk  = 156
  microchunk_cavgs(4)%path = 'microchunk_rejection/SLC/SLC_cavgs_4.mrc'
  microchunk_cavgs(4)%msk  = 156
  microchunk_cavgs(5)%path = 'microchunk_rejection/SLC/SLC_cavgs_5.mrc'
  microchunk_cavgs(5)%msk  = 156
  microchunk_cavgs(6)%path = 'microchunk_rejection/Betagal/BGAL_cavgs_1.mrc'
  microchunk_cavgs(6)%msk  = 200
  microchunk_cavgs(7)%path = 'microchunk_rejection/Betagal/BGAL_cavgs_3.mrc'
  microchunk_cavgs(7)%msk  = 200
  do ichunk=1, size(microchunk_cavgs)
    call get_mrcfile_info(microchunk_cavgs(ichunk)%path, ldim, 'M', smpd, .false.)
    allocate(cavg_imgs(ldim(3)))
    allocate(l_rejected(ldim(3)))
    l_rejected = .false.
    call stkio_r%open(microchunk_cavgs(ichunk)%path, smpd, 'read', bufsz=1)
    do icls=1, ldim(3)
      call cavg_imgs(icls)%new([ldim(1), ldim(2), 1], smpd, wthreads=.false.)
      call stkio_r%read(icls, cavg_imgs(icls))
    end do
    call stkio_r%close()
    ! outlier rejection
    mskrad = min(real(ldim(2)/2) - COSMSKHALFWIDTH - 1., 0.5 * real(microchunk_cavgs(ichunk)%msk)/smpd)
    call reject_outliers(cavg_imgs, mskrad, l_rejected)
    istk = 0
    stkname = swap_suffix(basename(microchunk_cavgs(ichunk)%path), '_rejected_outliers.mrc', '.mrc')
    call del_file(stkname)
    call stkio_r%open(microchunk_cavgs(ichunk)%path, smpd, 'read', bufsz=1)
    do icls=1, ldim(3)
      if( l_rejected(icls) ) then 
        istk = istk + 1
        call img%new([ldim(1), ldim(2), 1], smpd, wthreads=.false.)
        call stkio_r%read(icls, img)
        call img%write(stkname, istk)
        call img%kill()
      endif
    end do
    call stkio_r%close()
    call reject_auto(cavg_imgs, l_rejected)
    istk = 0
    stkname = swap_suffix(basename(microchunk_cavgs(ichunk)%path), '_rejected_auto.mrc', '.mrc')
    call del_file(stkname)
    call stkio_r%open(microchunk_cavgs(ichunk)%path, smpd, 'read', bufsz=1)
    do icls=1, ldim(3)
      if( l_rejected(icls) ) then 
        istk = istk + 1
        call img%new([ldim(1), ldim(2), 1], smpd, wthreads=.false.)
        call stkio_r%read(icls, img)
        call img%write(stkname, istk)
        call img%kill()
      endif
    end do
    istk = 0
    stkname = swap_suffix(basename(microchunk_cavgs(ichunk)%path), '_selected.mrc', '.mrc')
    call del_file(stkname)
    call stkio_r%open(microchunk_cavgs(ichunk)%path, smpd, 'read', bufsz=1)
    do icls=1, ldim(3)
      if( .not.l_rejected(icls) ) then 
        istk = istk + 1
        call img%new([ldim(1), ldim(2), 1], smpd, wthreads=.false.)
        call stkio_r%read(icls, img)
        call img%write(stkname, istk)
        call img%kill()
      endif
    end do
    call stkio_r%close()
    call dealloc_imgarr(cavg_imgs)
    deallocate(l_rejected)
  end do

end program main

module simple_nano_utils
  use simple_image,    only : image
  use simple_binimage, only : binimage
  ! use simple_syslib,   only : get_process_id
  use simple_rnd
  use simple_math
  use simple_defs ! singleton
  ! use simple_error!,  only: , simple_exception
  implicit none

  public :: remove_graphene_peaks
  private

contains

  ! This subroutine generates a mask that identifies the graphene
  ! peaks in the spectrum.
  subroutine generate_peakmask(img_spec, msk)
    type(image),          intent(inout) :: img_spec
    integer, allocatable, intent(inout) :: msk(:,:,:)
    integer, parameter :: WIDTH   =2 ! width of the shell for peak identification
    integer, parameter :: N_LAYERS=2 ! nb of layers to grow
    real,    parameter :: UP_LIM  =2.14, L_LIM=1.23 ! resolution where to fine peaks (2 bands)
    real,    pointer   :: rmat(:,:,:)
    real, allocatable  :: x(:),rmat_aux(:,:,:)
    real               :: thresh, smpd
    integer            :: ldim(3),i,j,h,k,sh,nframe,cnt,UFlim,LFlim,lims(2)
    type(binimage)     :: img_mask, img_mask_roavg
    nthr_glob=8
    if(allocated(msk)) deallocate(msk)
    ldim = img_spec%get_ldim()
    smpd = img_spec%get_smpd()
    call img_spec%get_rmat_ptr(rmat)
    allocate(rmat_aux(1:ldim(1),1:ldim(2),1), source= 0.)
    rmat_aux(1:ldim(1),1:ldim(2),1)=rmat(1:ldim(1),1:ldim(2),1)
    UFlim   = calc_fourier_index(max(2.*smpd,UP_LIM),ldim(1),smpd)
    lims(1) = UFlim-WIDTH
    lims(2) = UFlim+WIDTH
    cnt = 0
    do i = 1, ldim(1)
      do j = 1, ldim(2)
        h   = -int(ldim(1)/2) + i - 1
        k   = -int(ldim(2)/2) + j - 1
        sh  =  nint(hyp(real(h),real(k)))
        if(sh < lims(2) .and. sh > lims(1)) then
          !do nothing
        else
          rmat_aux(i,j,1) = 0. ! set to zero everything outside the considered shell
        endif
      enddo
    enddo
    x = pack(rmat_aux, abs(rmat_aux) > TINY)
    call otsu(x, thresh)
    deallocate(x)
    allocate(msk(ldim(1),ldim(2),1), source =0)
    where(rmat(1:ldim(1),1:ldim(2),1)>=thresh .and. abs(rmat_aux(:,:,1)) > TINY)  msk(:,:,1) = 1
    ! Do the same thing for the other band
    ! reset
    rmat_aux(:,:,1) = rmat(1:ldim(1),1:ldim(2),1)
    LFlim   = calc_fourier_index(max(2.*smpd,L_LIM),ldim(1),smpd)
    lims(1) = LFlim-WIDTH
    lims(2) = LFlim+WIDTH
    do i = 1, ldim(1)
      do j = 1, ldim(2)
        h   = -int(ldim(1)/2) + i - 1
        k   = -int(ldim(2)/2) + j - 1
        sh  =  nint(hyp(real(h),real(k)))
        if(sh < lims(2) .and. sh > lims(1)) then
          ! do nothing
        else
          rmat_aux(i,j,1) = 0.
        endif
      enddo
    enddo
    x = pack(rmat_aux, abs(rmat_aux) > 0.) ! pack where it's not 0
    call otsu(x, thresh)
    where(rmat(1:ldim(1),1:ldim(2),1)>=thresh  .and. abs(rmat_aux(:,:,1)) > TINY) msk(:,:,1) = 1
    ! Debug
    call img_mask%new_bimg(ldim, smpd)
    call img_mask%set_imat(msk)
    call img_mask%update_img_rmat()
    call img_mask_roavg%new_bimg(ldim,smpd)
    call img_mask%roavg(60,img_mask_roavg)
    call img_mask_roavg%get_rmat_ptr(rmat)
    where(rmat<TINY)
      rmat(:,:,:) = 0.
    elsewhere
      rmat(:,:,:) = 1.
    endwhere
    call img_mask_roavg%set_imat()
    call img_mask_roavg%grow_bins(N_LAYERS)
    call img_mask_roavg%get_imat(msk)
    call img_mask_roavg%kill_bimg
    call img_mask%kill_bimg
  end subroutine generate_peakmask


  ! This subroutine operates on raw_img on the basis of the
  ! provided logical mask. For each pxl in the msk, it subsitutes
  ! its value with the avg value in the neighbours of the pxl that
  ! do not belong to the msk. The neighbours are in a shpere of radius RAD.
  subroutine filter_peaks(raw_img,msk)
    type(image), intent(inout) :: raw_img
    integer,     intent(in)    :: msk(:,:,:)
    complex, pointer   :: cmat(:,:,:)
    real,    parameter :: RAD = 4.
    real        :: rabs, ang
    integer     :: nptcls,i,j,ii,jj,ldim(3),cnt,hp,h,k,kp,phys(3)
    ldim  = raw_img%get_ldim()
    call raw_img%norm()
    call raw_img%fft()
    call raw_img%get_cmat_ptr(cmat)
    do i = ldim(1)/2+1, ldim(1)
      do j = 1, ldim(2)
        if(msk(i,j,1)>0) then
          cnt = 0 ! number of pxls in the sphere that are not flagged by msk
          rabs = 0.
          do ii = 1, ldim(1)
            h  = ii-(ldim(1)/2+1)
            do jj = 1, ldim(2)
              if(ii/=i .and. jj/=j) then
                if(euclid(real([i,j]), real([ii,jj])) <= RAD) then
                  if(msk(ii,jj,1)<1) then
                    cnt = cnt + 1
                    k = jj-(ldim(2)/2+1)
                    phys = raw_img%comp_addr_phys(h,k,0)
                    rabs = rabs + sqrt(csq(cmat(phys(1),phys(2),1))) ! magnitude
                  endif
                endif
              endif
            enddo
          enddo
          h = i-(ldim(1)/2+1)
          k = j-(ldim(2)/2+1)
          phys = raw_img%comp_addr_phys(h,k,0)
          if(cnt == 0) then
            ! random number
            rabs = ran3()
          else
            rabs = rabs/real(cnt)
          endif
          ang = TWOPI*ran3()
          cmat(phys(1),phys(2),1) = rabs * cmplx(cos(ang), sin(ang))
        endif
      enddo
    enddo
    call raw_img%ifft()
  end subroutine filter_peaks

  subroutine remove_graphene_peaks(raw_img, spectrum, outfile)
    type(image), intent(inout) :: raw_img
    type(image), intent(inout) :: spectrum
    character(len=*), optional, intent(in) :: outfile
    integer, allocatable       :: msk(:,:,:)
    call generate_peakmask(spectrum, msk)
    call filter_peaks(raw_img,msk)
    if(present(outfile)) call raw_img%write(outfile)
  end subroutine remove_graphene_peaks
end module simple_nano_utils

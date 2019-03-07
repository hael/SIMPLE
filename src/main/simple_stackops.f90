! stack image processing routines for SPIDER/MRC files
module simple_stackops
include 'simple_lib.f08'
use simple_image,   only: image
use simple_oris,    only: oris
implicit none

public :: make_pattern_stack, binarize_stack, prepare_stack
public :: acf_stack, make_avg_stack, stats_imgfile, frameavg_stack
private
#include "simple_local_flags.inc"

contains

    !>  \brief  is for raising exception
    subroutine raise_exception( n, ldim, routine )
        integer, intent(in) :: n        !< num of images
        integer, intent(in) :: ldim(3)  !< logical dimensions
        character(len=*)    :: routine  !< Error message caller
        if( n < 1 .or. any(ldim == 0) )then
            write(logfhandle,*) routine
            write(logfhandle,*) 'The input stack is corrupt!'
            write(logfhandle,*) 'Number of images: ', n
            write(logfhandle,*) 'Logical dimensions: ', ldim
            THROW_HARD('procimgfile exception')
        endif
    end subroutine raise_exception

    !>  \brief  is for making a stack of vectors for PCA analysis
    !! \param fnameStack,fnamePatterns filenames for stack and pattern
    subroutine make_pattern_stack( fnameStack, fnamePatterns, l_mask, D, recsz, avg )
        character(len=*),  intent(in)  :: fnameStack, fnamePatterns
        logical,           intent(in)  :: l_mask(:,:,:) !< true for pixels to extract
        integer,           intent(out) :: D, recsz      !< record size
        real, allocatable, intent(out) :: avg(:)        !< frame stack average
        type(image)        :: img
        real, allocatable  :: pcavec(:)
        integer            :: n, fnum, ier, i, ldim(3), ldim_mask(3)
        logical            :: err
        call find_ldim_nptcls(fnameStack, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'make_pattern_stack' )
        ldim_mask = [size(l_mask, dim=1),size(l_mask, dim=2),size(l_mask, dim=3)]
        if( .not. all(ldim_mask .eq. ldim) )then
            write(logfhandle,*) 'ERROR! nonconforming matrix sizes'
            write(logfhandle,*) 'ldim of image: ', ldim
            write(logfhandle,*) 'ldim of mask : ', ldim_mask
            THROW_HARD('make_pattern_stack')
        endif
        ! build and initialise objects
        call img%new(ldim,1.)
        D = count(l_mask)
        allocate(pcavec(D), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('make_pattern_stack; simple_procimgfile, 1', alloc_stat)
        pcavec = 0.
        inquire(iolength=recsz) pcavec
        deallocate(pcavec)
        if( allocated(avg) ) deallocate(avg)
        allocate(avg(D), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('make_pattern_stack; simple_procimgfile, 2', alloc_stat)
        avg = 0.
        ! extract patterns and write to file
        call fopen(fnum, status='replace', action='readwrite', file=fnamePatterns,&
             access='direct', form='unformatted', recl=recsz, iostat=ier)
        call fileiochk('make_pattern_stack; simple_procimgfile', ier)
        write(logfhandle,'(a)') '>>> MAKING PATTERN STACK'
        do i=1,n
            call progress(i,n)
            call img%read(fnameStack, i)
            pcavec = img%serialize(l_mask)
            avg = avg + pcavec
            write(fnum,rec=i) pcavec
            deallocate(pcavec)
        end do
        avg = avg/real(n)
        allocate(pcavec(D), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('make_pattern_stack; simple_procimgfile, 3', alloc_stat)
        do i=1,n
            read(fnum,rec=i) pcavec
            pcavec = pcavec-avg
            write(fnum,rec=i) pcavec
        end do
        deallocate(pcavec)
        call fclose(fnum,errmsg='make_pattern_stack; simple_procimgfile')
        call img%kill
    end subroutine make_pattern_stack

    !>   make_avg_imgfile is for making an average of the entire stack
    !! \param fname stack filename
    !! \param avgname output filename
    subroutine make_avg_stack( fname, avgname, smpd )
        character(len=*), intent(in) :: fname, avgname
        real,             intent(in) :: smpd              !< sampling distance, resolution
        type(image)                  :: avg, img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'make_avg_imgfile' )
        call avg%new(ldim, smpd)
        call img%new(ldim, smpd)
        avg = 0.
        write(logfhandle,'(a)') '>>> MAKING GLOBAL STACK AVERAGE'
        do i=1,n
            call progress(i,n)
            call img%read(fname, i)
            call avg%add(img)
        end do
        call avg%div(real(n))
        call avg%write(avgname,1)
    end subroutine make_avg_stack

    !>  \brief  is for making frame averages of dose-fractionated image series
    !! \param fname2process  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param navg number of averages
    subroutine frameavg_stack( fname2process, fname, navg, smpd )
        character(len=*), intent(in) :: fname2process, fname
        integer,          intent(in) :: navg
        real,             intent(in) :: smpd
        type(image) :: img, avg
        integer     :: i, n, cnt, cnt2, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'frameavg_imgfile' )
        call img%new(ldim,smpd)
        call avg%new(ldim,smpd)
        cnt = 0
        cnt2 = 0
        write(logfhandle,'(a)') '>>> AVERAGING FRAMES'
        do i=1,n
            call progress(i,n)
            cnt = cnt+1
            call img%read(fname2process, i)
            if( cnt <= navg )then
                call avg%add(img)
            endif
            if( cnt == navg )then
                cnt2 = cnt2+1
                call avg%write(fname, cnt2)
                cnt = 0
                avg = 0.
            endif
        end do
        call img%kill
        call avg%kill
    end subroutine frameavg_stack

    !>  \brief  is for calculating the acf of a stack
    !! \param fname2acf  output filename
    !! \param fname  input filename
    subroutine acf_stack( fname2acf, fname )
        character(len=*), intent(in) :: fname2acf, fname
        type(image)   :: img
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2acf, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'acf_imgfile' )
        call img%new(ldim,1.)
        write(logfhandle,'(a)') '>>> CALCULATING ACF:S OF THE IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2acf, i)
            call img%acf
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine acf_stack

    !>  \brief  is for calculating image statistics
    subroutine stats_imgfile( fname, os, msk )
        character(len=*), intent(in)  :: fname
        class(oris),      intent(out) :: os
        real, optional,   intent(in)  :: msk
        real              :: ave, sdev, med, minv, maxv, spec
        real, allocatable :: spectrum(:)
        type(image)       :: img
        integer           :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'stats_imgfile' )
        call img%new(ldim,1.)
        call os%new(n)
        write(logfhandle,'(a)') '>>> CALCULATING STACK STATISTICS'
        do i=1,n
            call progress(i,n)
            call img%read(fname,i)
            call img%stats('foreground', ave, sdev, maxv, minv, med=med, msk=msk)
            call img%fft()
            call img%spectrum('power',spectrum)
            spec = sum(spectrum)/real(size(spectrum))
            call os%set(i, 'ave',   ave)
            call os%set(i, 'sdev',  sdev)
            call os%set(i, 'maxv',  maxv)
            call os%set(i, 'minv',  minv)
            call os%set(i, 'med',   med)
            call os%set(i, 'spec',  spec)
            call os%set(i, 'dynrange', maxv - minv)
        end do
        call img%kill
    end subroutine stats_imgfile

    ! This subroutine calculates the logaritm of all the micrographs in a stack.
    ! It is meant to be applied to power spectra images, in order to decrease
    ! the influence of the central pixel.
    subroutine calc_log_stack( fname2process, fname, smpd )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd
        real, allocatable :: rmat(:,:,:)
        type(image)       :: img
        integer           :: n, ldim(3), i
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'calc_log_imgfile' )
        call img%new(ldim,smpd)
        do i=1,n
            call img%read(fname2process, i)
            call img%scale_pixels([1.,real(ldim(1))*2.+1.]) !not to have problems in calculating the log
            rmat = img%get_rmat()
            rmat = log(rmat)
            call img%set_rmat(rmat)
            call img%write(fname, i)
            rmat = 0.
        end do
        call img%kill
        deallocate(rmat)
    end subroutine calc_log_stack

    ! This subroutine calculates the logaritm of all the micrographs in a stack.
    ! It is meant to be applied to power spectra images, in order to decrease
    ! the influence of the central pixel.
    ! Edge detaction is more challengin in close to focus powerspectra
    ! images. There are more appropriate settings for this case, in
    ! order to set them, flag the variable 'close_to_focus'
    subroutine apply_tvf_stack( fname2process, fname, smpd, lambda, close_to_focus)
        use simple_tvfilter
        character(len=*),  intent(in) :: fname2process, fname
        real,              intent(in) :: smpd
        real,              intent(in) :: lambda ! parameter for tv filtering
        logical, optional, intent(in) :: close_to_focus(:)
        type(image)       :: img
        integer           :: n, ldim(3), i
        type(tvfilter)    :: tvf
        real              :: llambda !to override var lambda in case of close_to_focus imaging
        logical, allocatable :: cclose_to_focus(:)
        call find_ldim_nptcls(fname2process, ldim, n)
        if(present(close_to_focus) .and. size(close_to_focus) .ne. n) THROW_HARD('Wrong dim in close_to_focus input; apply_tvf_stack')
        allocate (cclose_to_focus(n), source = .false.)
        if(present(close_to_focus)) cclose_to_focus = close_to_focus
        ldim(3) = 1
        call tvf%new()
        call raise_exception( n, ldim, 'apply_tvf_imgfile' )
        call img%new(ldim,smpd)
        llambda = lambda
        do i = 1, n
            if(cclose_to_focus(i)) llambda = 10.
            call img%read(fname2process, i)
            !call img%scale_pixels([1.,real(ldim(1))*2.+1.]) !not to have problems in calculating the log
            call tvf%apply_filter(img, llambda)
            call img%write(fname, i)
            llambda = lambda !restore input
        enddo
        call img%kill
    end subroutine apply_tvf_stack

    ! This subroutine takes in input a stack of images and performs
    ! a series of operations on it in order to prepare it to be
    ! binarised. It is written for power spectra stacks.
    ! Procedure:
    !                -) logarithm of the image (after pixel scaling)
    !                -) background subtraction
    !                -) median filter
    subroutine prepare_stack( fname2process, fname, smpd, lp, close_to_focus)
        use simple_procimgfile,   only : subtr_backgr_imgfile, real_filter_imgfile
        character(len=*),  intent(in) :: fname2process, fname
        real,              intent(in) :: smpd, lp
        logical, optional, intent(in) :: close_to_focus(:)
        integer     :: winsz !for median filtering
        integer     :: n, ldim(3), i
        type(image) :: img
        real        :: lambda !for tv filtering
        logical, allocatable :: cclose_to_focus(:)
        call find_ldim_nptcls(fname2process, ldim, n)
        if(present(close_to_focus) .and. size(close_to_focus) .ne. n) THROW_HARD('Wrong dim in close_to_focus input; prepare_stack')
        allocate (cclose_to_focus(n), source = .false.)
        if(present(close_to_focus)) cclose_to_focus = close_to_focus
        ldim(3) = 1
        call img%new(ldim,smpd)
        lambda  = 1 !10. for close to focus images
        winsz   = 7 !just for close to focus images
        call raise_exception( n, ldim, 'calc_log_imgfile' )
        call calc_log_stack( fname2process, fname, smpd ) !logaritm of the images in the stack
        call apply_tvf_stack(fname,fname,smpd,lambda,cclose_to_focus)
        call subtr_backgr_imgfile( fname, fname, smpd, lp )
        write(logfhandle,'(a)') '>>> REAL-SPACE FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname, i)
            if(cclose_to_focus(i)) then !otherwise do not median filter
                call img%real_space_filter(winsz, 'median')
            endif
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine prepare_stack

    !This subroutine performs Canny automatic edge detection on a stack
    !of images.
    !If discard_borders is set to true, the external frame of the binary
    !image is set to 0, to take care of border effects.
    !default: discard_borders=.false.
    subroutine binarize_stack(fname2process,fname,smpd,discard_borders)
      use simple_segmentation, only : canny
      character(len=*),  intent(in) :: fname2process, fname
      real,              intent(in) :: smpd
      logical, optional, intent(in) :: discard_borders
      real, allocatable :: rmat_bin(:,:,:)
      logical     :: ddiscard_borders
      integer     :: n, ldim(3), i
      type(image) :: img
      call find_ldim_nptcls(fname2process, ldim, n)
      ldim(3) = 1
      ddiscard_borders = .false.
      if(present(discard_borders)) ddiscard_borders = discard_borders
      call img%new(ldim,smpd)
      do i = 1, n
          call img%read(fname2process, i)
          call canny(img)
          if(ddiscard_borders) then
              rmat_bin = img%get_rmat()
              rmat_bin(1:int(ldim(1)/10), :, 1) = 0.       !bottom horizontal border
              rmat_bin(ldim(1)-int(ldim(1)/10):ldim(1),:, 1) = 0. !top horizontal border
              rmat_bin(:,1:int(ldim(2)/10), 1) = 0.        !bottom vertical border
              rmat_bin(:,ldim(2)-int(ldim(2)/10):ldim(2), 1) = 0.  !top vertical border
              call img%set_rmat(rmat_bin)
          endif
          call img%write(fname,i)
      enddo
      call img%kill
    end subroutine binarize_stack
end module simple_stackops

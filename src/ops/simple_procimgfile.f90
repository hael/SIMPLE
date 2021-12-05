! stack image processing routines for SPIDER/MRC files
module simple_procimgfile
include 'simple_lib.f08'
use simple_image,   only: image
use simple_oris,    only: oris
use simple_tvfilter
implicit none

!! Basic Operations
public :: copy_imgfile, diff_imgfiles, pad_imgfile, resize_imgfile, clip_imgfile, mirror_imgfile
public :: random_selection_from_imgfile, resize_and_clip_imgfile
public :: subtr_backgr_imgfile, random_cls_from_imgfile, selection_from_tseries_imgfile
public :: roavg_imgfile
!! Normalisation
public :: norm_bin_imgfile, norm_imgfile, norm_ext_imgfile, noise_norm_imgfile
public :: shellnorm_imgfile, matchfilt_imgfile
!! Morphological
public :: neg_imgfile, bin_imgfile
public :: mask_imgfile, taper_edges_imgfile
!! Filters
public :: ft2img_imgfile, masscen_imgfile, cure_imgfile, apply_bfac_imgfile
public :: shift_imgfile, bp_imgfile, shrot_imgfile, add_noise_imgfile, nlmean_imgfile, corrfilt_imgfile
public :: corrfilt_tseries_imgfile, tvfilt_tseries_imgfile, real_filter_imgfile, sharpen_imgfile, phase_rand_imgfile
public :: apply_ctf_imgfile, tvfilter_imgfile
private


interface subtr_backgr_imgfile
    module procedure subtr_backgr_imgfile_1
    module procedure subtr_backgr_imgfile_2
end interface

interface corrfilt_imgfile
    module procedure corrfilt_imgfile_1
    module procedure corrfilt_imgfile_2
end interface
#include "simple_local_flags.inc"

contains

    !>  \brief  is for raising exception
    subroutine raise_exception_imgfile( n, ldim, routine )
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
    end subroutine raise_exception_imgfile

    !>  \brief  is for copying image file
    !! \param fname2copy,fname Filenames to copy
    subroutine copy_imgfile( fname2copy, fname, smpd, fromto )
        character(len=*), intent(in) :: fname2copy, fname
        real,             intent(in) :: smpd       !< sampling distance
        integer,          intent(in) :: fromto(2)  !< range
        type(image) :: img
        integer     :: n, i, cnt, ldim(3)
        call find_ldim_nptcls(fname2copy, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile(n, ldim, 'copy_imgfile')
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> COPYING IMAGES'
        cnt = 0
        do i=fromto(1),fromto(2)
            cnt = cnt+1
            call progress(cnt, fromto(2)-fromto(1)+1)
            call img%read(fname2copy, i)
            call img%write(fname, cnt)
        end do
        call img%kill
    end subroutine copy_imgfile

    subroutine diff_imgfiles( fname, fname2subtr, fname_out, smpd, fromto )
        character(len=*),  intent(in) :: fname, fname2subtr, fname_out
        real,              intent(in) :: smpd       !< sampling distance
        integer, optional, intent(in) :: fromto(2)  !< ranges
        type(image) :: img2, diffimg
        integer     :: n, i, cnt, ldim(3), ldim2(3),n2
        call find_ldim_nptcls(fname, ldim, n)
        call find_ldim_nptcls(fname2subtr, ldim2, n2)
        ldim(3) = 1; ldim2(3)=1
        call raise_exception_imgfile(n, ldim, 'diff_imgfiles')
        if ( n /= n2 .or. ldim(1) /=ldim2(1)  .or. ldim(2) /= ldim2(2)) &
            THROW_HARD ('procimgfile exception; diff_imgfiles mismatch')
        call diffimg%new(ldim,smpd)
        call img2%new(ldim,smpd)
        if( n >= 1 )then
            write(logfhandle,'(a)') '>>> SUBTRACTING IMAGES'
            if( present(fromto) )then
                cnt = 0
                do i=fromto(1),fromto(2)
                    cnt = cnt+1
                    call progress(cnt, fromto(2)-fromto(1)+1)
                    call diffimg%read(fname, i)
                    call img2%read(fname2subtr, i)
                    diffimg = diffimg - img2
                    call diffimg%write(fname_out, cnt)
                end do
            else
                do i=1,n
                   call progress(i, n)
                    call diffimg%read(fname, i)
                    call img2%read(fname2subtr, i)
                    diffimg = diffimg - img2
                    call diffimg%write(fname_out, cnt)
                end do
            endif
        endif
        call diffimg%kill
        call img2%kill
    end subroutine diff_imgfiles

    !>  \brief pad_imgfile is for padding
    !! \param fname2pad,fname filename strings
    !! \param ldim_pad logical dimension of padding
    !! \param smpd sampling distance
    !!
    subroutine pad_imgfile( fname2pad, fname, ldim_pad, smpd )
        character(len=*), intent(in) :: fname2pad, fname
        integer,          intent(in) :: ldim_pad(3)
        real,             intent(in) :: smpd
        type(image)                  :: img, img_pad
        integer                      :: n, i, ldim(3)
        real                         :: ave, sdev, maxv, minv, med
        call find_ldim_nptcls(fname2pad, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'pad_imgfile' )
        if( ldim_pad(1) >= ldim(1) .and. ldim(2) >= ldim(2)&
             .and. ldim(3) >= ldim(3) )then
            call img%new(ldim,smpd)
            call img_pad%new(ldim_pad,smpd)
            write(logfhandle,'(a)') '>>> PADDING IMAGES'
            do i=1,n
                call progress(i,n)
                call img%read(fname2pad, i)
                if( img%is_ft() )then
                    call img%pad(img_pad) ! FT state preserved
                else
                    ! get background statistics
                    call img%stats('background', ave, sdev, maxv, minv, med=med)
                    call img%pad(img_pad, backgr=med) ! FT state preserved
                endif
                call img_pad%write(fname, i)
            end do
            call img%kill
            call img_pad%kill
        end if
    end subroutine pad_imgfile

    !>  \brief resize_imgfile is for resizing
    !! \param fname2resize,fname string filenames
    !! \param smpd sampling distance
    !! \param ldim_new new logical dimension
    !! \param smpd_new new sampling distance
    !! \param fromptop index range for copying
    !!
    subroutine resize_imgfile( fname2resize, fname, smpd, ldim_new, smpd_new, fromptop )
        character(len=*),  intent(in)  :: fname2resize, fname
        real,              intent(in)  :: smpd
        integer,           intent(in)  :: ldim_new(3)
        real,              intent(out) :: smpd_new
        integer, optional, intent(in)  :: fromptop(2)
        type(image) :: img, img_resized, blank_img
        integer     :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2resize, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'resize_imgfile' )
        call img%new(ldim,smpd,wthreads=.false.)
        call img_resized%new(ldim_new,smpd,wthreads=.false.) ! this sampling distance will be overwritten
        write(logfhandle,'(a)') '>>> RESIZING IMAGES'
        if( present(fromptop) )then
            prange = fromptop
        else
            prange(1) = 1
            prange(2) = n
        endif
        sz  = prange(2)-prange(1)+1
        if( all(prange == 0) )then
            call simple_touch(fname)
        else
            cnt = 0
            do i=prange(1),prange(2)
                cnt = cnt+1
                call progress(cnt,sz)
                call img%read(fname2resize, i)
                call img%fft()
                if( ldim_new(1) <= ldim(1) .and. ldim_new(2) <= ldim(2)&
                     .and. ldim_new(3) <= ldim(3) )then
                    call img%clip(img_resized)
                else
                    call img%pad(img_resized)
                endif
                call img_resized%ifft()
                call img_resized%write(fname, cnt)
            end do
            smpd_new = img_resized%get_smpd()
        endif
        call img%kill
        call img_resized%kill
        call blank_img%kill
    end subroutine resize_imgfile

    subroutine clip_imgfile( fname2clip, fname, ldim_clip, smpd )
        character(len=*), intent(in) :: fname2clip, fname
        integer,          intent(in) :: ldim_clip(3)
        real,             intent(in) :: smpd
        type(image) :: img, img_clip
        integer     :: n, i, ldim(3)
        call find_ldim_nptcls(fname2clip, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'clip_imgfile' )
        ! do the work
        if( ldim_clip(1) <= ldim(1) .and. ldim_clip(2) <= ldim(2)&
             .and. ldim_clip(3) <= ldim(3) )then
            call img%new(ldim,smpd)
            call img_clip%new(ldim_clip,smpd)
            write(logfhandle,'(a)') '>>> CLIPPING IMAGES'
            do i=1,n
                call progress(i,n)
                call img%read(fname2clip, i)
                call img%clip(img_clip) ! FT state preserved
                call img_clip%write(fname, i)
            end do
            call img%kill
            call img_clip%kill
        end if
    end subroutine clip_imgfile

    subroutine roavg_imgfile( fname2roavg, fname, angstep, smpd )
        character(len=*), intent(in) :: fname2roavg, fname
        integer,          intent(in) :: angstep
        real,             intent(in) :: smpd
        type(image) :: img, img_roavg
        integer     :: n, i, ldim(3)
        call find_ldim_nptcls(fname2roavg, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'roavg_imgfile' )
        call img%new(ldim,smpd)
        call img_roavg%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> ROTATIONALLY AVERAGING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2roavg, i)
            call img%roavg(angstep, img_roavg)
            call img_roavg%write(fname, i)
        end do
        call img%kill
        call img_roavg%kill
    end subroutine roavg_imgfile

    !>  \brief mirror imgfile
    subroutine mirror_imgfile( fname2mirr, fname, mirr_flag, smpd )
        character(len=*), intent(in) :: fname2mirr, fname
        character(len=1), intent(in) :: mirr_flag
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: n, i, ldim(3)
        call find_ldim_nptcls(fname2mirr, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'clip_imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> MIRRORING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2mirr, i)
            call img%mirror(mirr_flag)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine mirror_imgfile

    !>  \brief  resize_and_clip_imgfile is for resizing and clipping
    !! \param fname2resize  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param ldim_new new logical dimension
    !! \param ldim_clip clipped logical dimension
    !! \param smpd_new new sampling distance
    !! \param fromptop  index range for copying
    !!
    subroutine resize_and_clip_imgfile( fname2resize, fname, smpd, ldim_new, ldim_clip, smpd_new, fromptop )
        character(len=*),  intent(in)  :: fname2resize, fname
        real,              intent(in)  :: smpd
        integer,           intent(in)  :: ldim_new(3)
        integer,           intent(in)  :: ldim_clip(3)
        real,              intent(out) :: smpd_new
        integer, optional, intent(in)  :: fromptop(2)
        type(image) :: img, img_resized, img_clip
        integer     :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2resize, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'resize_imgfile' )
        call img%new(ldim,smpd)
        call img_resized%new(ldim_new,smpd) ! this sampling distance will be overwritten
        write(logfhandle,'(a)') '>>> RESIZING IMAGES'
        if( present(fromptop) )then
            prange = fromptop
        else
            prange(1) = 1
            prange(2) = n
        endif
        sz  = prange(2)-prange(1)+1
        cnt = 0
        do i=prange(1),prange(2)
            cnt = cnt+1
            call progress(cnt,sz)
            call img%read(fname2resize, i)
            call img%fft()
            call img%clip(img_resized)
            call img_resized%ifft()
            call img_clip%new(ldim_clip,img_resized%get_smpd())
            if( ldim_clip(1) <= ldim_new(1) .and. ldim_clip(2) <= ldim_new(2)&
                 .and. ldim_clip(3) <= ldim_new(3) )then
                call img_resized%clip(img_clip)
            else
                call img_resized%pad(img_clip)
            endif
            call img_clip%write(fname, cnt)
        end do
        smpd_new = img_resized%get_smpd()
        call img%kill
        call img_resized%kill
        call img_clip%kill
    end subroutine resize_and_clip_imgfile

    subroutine subtr_backgr_imgfile_1( fname2subtr, fname, smpd, lp )
        character(len=*), intent(in) :: fname2subtr, fname
        real,             intent(in) :: smpd, lp
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2subtr, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'subtr_backgr_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SUBTRACTING BACKGROUND FROM IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2subtr, i)
            call img%subtr_backgr(lp)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine subtr_backgr_imgfile_1

    subroutine subtr_backgr_imgfile_2( fname2subtr, fname, smpd, lp, mask )
        character(len=*),  intent(in) :: fname2subtr, fname
        real,              intent(in) :: smpd, lp
        logical,           intent(in) :: mask(:) !which imgs of the stack to consider
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2subtr, ldim, n)
        ldim(3) = 1
        if(size(mask) .ne. n) THROW_HARD('Wrong dim in input mask; subtr_backgr_imgfile')
        call raise_exception_imgfile( n, ldim, 'subtr_backgr_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SUBTRACTING BACKGROUND FROM IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2subtr, i)
            if(mask(i)) then
                call img%subtr_backgr(lp)
            endif
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine subtr_backgr_imgfile_2

    subroutine norm_bin_imgfile( fname2norm, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> BIN NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%norm_bin
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine norm_bin_imgfile

    subroutine norm_imgfile( fname2norm, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%norm()
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine norm_imgfile

    !>  \brief  is for normalization [mean,stdev] ->  [new mean, new stdev]
    !! \param fname2norm  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param avg,sdev new average and standard deviation
    subroutine norm_ext_imgfile( fname2norm, avg, sdev, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in)  :: avg, sdev, smpd
        type(image)   :: img
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%norm_ext(avg, sdev)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine norm_ext_imgfile

    !>  \brief  is for noise normalization
    !! \param fname2norm  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param msk masking threshold
    subroutine noise_norm_imgfile( fname2norm, msk, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: msk, smpd
        type(image)          :: img
        integer              :: i, n, ldim(3)
        logical, allocatable :: lmsk(:,:,:)
        real                 :: sdev_noise
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'noise_norm_imgfile' )
        call img%disc(ldim, smpd, msk, lmsk)
        write(logfhandle,'(a)') '>>> NOISE NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%noise_norm(lmsk, sdev_noise)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine noise_norm_imgfile

    !>  \brief  is for noise normalization
    !! \param fname2norm  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    subroutine shellnorm_imgfile( fname2norm, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'noise_norm_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SHELL NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%shellnorm()
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine shellnorm_imgfile

    !>  \brief  is for noise normalization
    !! \param fname2norm  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    subroutine matchfilt_imgfile( fname2filt, fname, frcs_fname, smpd )
        use simple_class_frcs,    only: class_frcs
        use simple_estimate_ssnr, only: fsc2optlp_sub
        character(len=*), intent(in) :: fname2filt, fname, frcs_fname
        real,             intent(in) :: smpd
        type(class_frcs)  :: frcs
        type(image)       :: img
        real, allocatable :: frc(:), filter(:)
        integer           :: i, n, ldim(3)
        call frcs%read(frcs_fname)
        call find_ldim_nptcls(fname2filt, ldim, n)
        ldim(3) = 1
        if( frcs%get_nprojs().ne.n )then
            write(logfhandle,*) '# imgs: ',n
            write(logfhandle,*) '# class_frcs: ',frcs%get_nprojs()
            THROW_HARD('inconsistent dimensions; matchfilt_imgfile')
        endif
        call raise_exception_imgfile( n, ldim, 'matchfilt_imgfile' )
        call img%new(ldim,smpd)
        if( frcs%get_filtsz().ne.img%get_filtsz() )then
            write(logfhandle,*) 'img filtsz: ',img%get_filtsz()
            write(logfhandle,*) 'frcs filtsz: ',frcs%get_filtsz()
            THROW_HARD('Inconsistent filter dimensions; matchfilt_imgfile')
        endif
        write(logfhandle,'(a)') '>>> SHELL NORMALIZING AND FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filt, i)
            frc = frcs%get_frc(i, ldim(1))
            allocate(filter(size(frc)), source=1.)
            call fsc2optlp_sub(frcs%get_filtsz(), frc, filter)
            where( filter < TINY )filter = 0.
            call img%fft()
            call img%shellnorm_and_apply_filter(frc)
            call img%ifft()
            call img%write(fname, i)
            deallocate(frc,filter)
        end do
        call img%kill
    end subroutine matchfilt_imgfile

    !>  \brief  is for contrast inversion
    !! \param fname2neg  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    subroutine neg_imgfile( fname2neg, fname, smpd )
        character(len=*), intent(in) :: fname2neg, fname
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2neg, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'neg_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> INVERTING IMAGE CONTRAST'
        do i=1,n
            call progress(i,n)
            call img%read(fname2neg, i)
            call img%neg()
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine neg_imgfile

    !>  \brief  is for centering based on center of mass
    !! \param fname2masscen  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param  msk threshold
    subroutine masscen_imgfile( fname2masscen, fname, smpd, lp, msk )
        character(len=*), intent(in) :: fname2masscen, fname
        real,             intent(in) :: smpd, lp, msk
        type(image) :: img
        integer     :: i, n, ldim(3)
        real        :: xyz(3)
        call find_ldim_nptcls(fname2masscen, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'masscen_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> CENTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2masscen, i)
            xyz = img%calc_shiftcen(lp, msk)
            print *, xyz(:2)
            call img%shift(xyz)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine masscen_imgfile

    !>  \brief  is for curing
    !! \param fname2cure  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    subroutine cure_imgfile( fname2cure, fname, smpd )
        character(len=*), intent(in) :: fname2cure, fname
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: i, n, n_nans, io_stat,filnum, ldim(3)
        real        :: maxv, minv, ave, sdev
        call find_ldim_nptcls(fname2cure, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'cure_imgfile' )
        call fopen(filnum,'cure_stats.txt',status='replace',iostat=io_stat)
        call fileiochk("cure_imgfile error", io_stat)
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> CURING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2cure, i)
            call img%cure(maxv, minv, ave, sdev, n_nans)
            write(filnum,'(A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,I9)') 'MAX:', maxv,&
                 'MIN:', minv, 'AVE:', ave, 'SDEV:', sdev, 'NANS:', n_nans
            call img%write(fname, i)
        end do
        call fclose(filnum,errmsg="cure_imgfile close error")
        call img%kill
    end subroutine cure_imgfile

    !>  \brief  is for adding noise
    !! \param fname2process  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param snr signal to noise ratio
    subroutine add_noise_imgfile( fname2process, fname, snr, smpd )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: snr, smpd
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'add_noise_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> ADDING NOISE TO THE IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            call img%add_gauran(snr)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine add_noise_imgfile

    !>  \brief  is for Fourier transformation visualization
    !! \param fname2process  output filename
    !! \param fname  input filename
    !! \param which file type
    subroutine ft2img_imgfile( fname2process, fname, which )
        character(len=*), intent(in) :: fname2process, fname
        character(len=*), intent(in) :: which
        type(image)   :: img, img2
        integer       :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'ft2img_imgfile' )
        call img%new(ldim,1.)
        if( ldim(3) == 1 )then
            do i=1,n
                call img%read(fname2process, i)
                call img%ft2img(which, img2)
                call img2%write(fname, i)
            end do
        else
            call img%read(fname2process)
            call img%ft2img(which, img2)
            call img2%write(fname)
        endif
        call img%kill
        call img2%kill
    end subroutine ft2img_imgfile

    !>  \brief  is for origin shifting a stack according to info in o
    !! \param fname2shift  output filename
    !! \param fname  input filename
    !! \param o orientation object used for shifting
    !! \param mul multiplier
    !! \param round optional flag for rounding
    !! \param smpd  sampling distance
    subroutine shift_imgfile( fname2shift, fname, o, smpd, mul, round )
        use simple_oris, only: oris
        character(len=*),  intent(in)    :: fname2shift, fname
        class(oris),       intent(inout) :: o
        real,              intent(in)    :: smpd
        real,    optional, intent(in)    :: mul
        logical, optional, intent(in)    :: round
        type(image)   :: img
        integer       :: n, i, ldim(3)
        real          :: x, y, xhere, yhere
        logical       :: rround
        call find_ldim_nptcls(fname2shift, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'shift_imgfile' )
        rround = .false.
        if( present(round) ) rround = round
        if( n /= o%get_noris() ) THROW_HARD('inconsistent nr entries; shift_imgfile')
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SHIFTING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2shift, i)
            x = o%get(i, 'x')
            y = o%get(i, 'y')
            if( present(mul) )then
                xhere = x*mul
                yhere = y*mul
            else
                xhere = x
                yhere = y
            endif
            if( rround )then
                xhere = real(nint(xhere))
                yhere = real(nint(yhere))
            endif
            call img%shift([-xhere,-yhere,0.])
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine shift_imgfile

    !>  \brief  is for band-pass filtering
    subroutine bp_imgfile( fname2filter, fname, smpd, hp, lp, width )
        character(len=*), intent(in) :: fname2filter, fname
        real,             intent(in) :: smpd, hp, lp
        real, optional,   intent(in) :: width
        type(image) :: img
        integer     :: n, i, ldim(3)
        real        :: wwidth
        wwidth = 10.0
        if( present(width) ) wwidth = width
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'bp_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> BAND-PASS FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filter, i)
            call img%bp(hp,lp,width=width)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine bp_imgfile

    !>  \brief  is for applying B-factor
    subroutine apply_bfac_imgfile( fname2filter, fname, smpd, bfac )
        character(len=*), intent(in) :: fname2filter, fname
        real,             intent(in) :: smpd, bfac
        type(image) :: img
        integer     :: n, i, ldim(3)
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'apply_bfac_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> APPLYING B-FACTOR TO THE IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filter, i)
            call img%apply_bfac(bfac)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine apply_bfac_imgfile

    !>  \brief real_filter_imgfile is for real-space filtering
    subroutine real_filter_imgfile( fname2filter, fname, smpd, which_filter, winsz )
        character(len=*), intent(in) :: fname2filter, fname
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: which_filter
        integer,          intent(in) :: winsz
        type(image)     :: img
        integer         :: n, i, ldim(3)
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'real_filter_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> REAL-SPACE FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filter, i)
            call img%real_space_filter(winsz, which_filter)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine real_filter_imgfile

    subroutine sharpen_imgfile( fname2filter, fname, smpd, bfac )
        character(len=*), intent(in) :: fname2filter, fname
        real,             intent(in) :: smpd, bfac
        type(image)     :: img
        integer         :: n, i, ldim(3)
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'sharpen_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SHARPENING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filter, i)
            call img%fft
            call img%apply_bfac(bfac)
            call img%ifft
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine sharpen_imgfile

    !>  \brief  is for phase randomization
    subroutine phase_rand_imgfile( fname2process, fname, smpd, lp )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, lp
        type(image) :: img
        integer :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'phase_rand_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> PHASE RANDOMIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            call img%phase_rand(lp)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine phase_rand_imgfile

    !>  \brief  is for apply the tv filter to an image file
    subroutine tvfilter_imgfile( fname2process, fname, smpd, lambda )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, lambda
        type(tvfilter) :: tv
        type(image)    :: img
        integer        :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'tvfilter_imgfile' )
        call img%new(ldim,smpd)
        call tv%new()
        write(logfhandle,'(a)') '>>> APPLYING TV FILTER TO IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            call tv%apply_filter(img,lambda)
            call img%write(fname, i)
        end do
        call img%kill
        call tv%kill
    end subroutine tvfilter_imgfile

    !>  \brief  is for apply the non-local mean filter to an image file
    subroutine nlmean_imgfile( fname2process, fname, smpd, noise_sdev )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd
        real, optional,   intent(in) :: noise_sdev
        type(image)    :: img
        integer        :: n, i, ldim(3)
        logical        :: noise_sdev_present
        noise_sdev_present = present(noise_sdev)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'nlmean_imgfile' )
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> APPLYING NLMEAN FILTER TO IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            if( noise_sdev_present )then
                call img%nlmean(sdev_noise=noise_sdev)
            else
                call img%nlmean
            endif
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine nlmean_imgfile

    ! convolution-based filtering globally over the entire time-series to make maximum
    ! use of the redundancu in the data. Discovered by mistake and we are surprised it works
    ! since SIMPLE is not supposed to support non-square transforms
    subroutine corrfilt_tseries_imgfile( fname2process, fname, smpd, sigma, lp )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, sigma, lp
        type(image)    :: imgproc, gau, imgresult
        integer        :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        call raise_exception_imgfile( n, ldim, 'corrfilt_tseries_imgfile' )
        call gau%new(ldim, smpd)
        call imgproc%new(ldim,smpd)
        call imgresult%new(ldim,smpd)
        call imgresult%set_ft(.true.)
        call gau%gauimg3D(sigma, sigma, sigma, cutoff=5.)
        call gau%fft
        write(logfhandle,'(a)') '>>> APPLYING CORR FILTER TO TIME-SERIES'
        call imgproc%read(fname2process)
        call imgproc%fft
        call imgproc%phase_corr(gau,imgresult,lp)
        call imgresult%write(fname)
        call imgproc%kill
        call gau%kill
        call imgresult%kill
    end subroutine corrfilt_tseries_imgfile

    subroutine tvfilt_tseries_imgfile( fname2process, fname, smpd, lambda )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: lambda, smpd
        type(tvfilter) :: tv
        type(image)    :: img
        integer        :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        call raise_exception_imgfile( n, ldim, 'tvfilt_tseries_imgfile' )
        call tv%new()
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> APPLYING TV FILTER TO TIME-SERIES'
        call img%read(fname2process)
        call tv%apply_filter_3d(img,lambda)
        call img%write(fname)
        call img%kill
        call tv%kill
    end subroutine tvfilt_tseries_imgfile

    !>  \brief  is for apply the correlation filter to an image file
    !    fname2process: images to be processed
    !    fname        : output images
    subroutine corrfilt_imgfile_1( fname2process, fname, smpd, sigma, lp, isvol )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, sigma, lp
        logical,          intent(in) :: isvol
        type(image)    :: imgproc, gau, imgresult
        integer        :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        call raise_exception_imgfile( n, ldim, 'corrfilt_imgfile_1' )
        if( .not. isvol ) ldim(3) = 1
        call gau%new(ldim, smpd)
        call imgproc%new(ldim,smpd)
        call imgresult%new(ldim,smpd)
        call imgresult%set_ft(.true.)
        if(ldim(3) == 1) then !2D case
            call gau%gauimg2D(sigma, sigma, cutoff=5.)
            call gau%fft
            write(logfhandle,'(a)') '>>> APPLYING CORR FILTER TO IMAGES'
            do i=1,n
                call progress(i,n)
                call imgproc%read(fname2process, i)
                call imgproc%fft
                call imgproc%phase_corr(gau,imgresult,lp)
                call imgresult%write(fname, i)
                call imgresult%zero_and_flag_ft
            end do
        else !3D case
            call gau%gauimg3D(sigma, sigma, sigma, cutoff=5.)
            call gau%fft
            write(logfhandle,'(a)') '>>> APPLYING CORR FILTER TO VOLUME'
            call imgproc%read(fname2process)
            call imgproc%fft
            call imgproc%phase_corr(gau,imgresult,lp)
            call imgresult%write(fname)
        endif
        call imgproc%kill
        call gau%kill
        call imgresult%kill
    end subroutine corrfilt_imgfile_1

    subroutine corrfilt_imgfile_2( fname2process, fname, element, smpd, lp )
        use simple_atoms, only : atoms
        character(len=*), intent(in) :: fname2process, fname
        character(len=2), intent(in) :: element
        real,             intent(in) :: smpd, lp
        type(image)    :: imgproc, img_atom, imgresult
        type(atoms)    :: atom
        integer        :: n, i, ldim(3)
        real           :: cutoff
        call find_ldim_nptcls(fname2process, ldim, n)
        call raise_exception_imgfile( n, ldim, 'corrfilt_imgfile_2' )
        if(ldim(3) == 1) then !2D case
            THROW_HARD('Convolution with atom is possible just in 3D; corrfilt_imgfile_2')
        else !3D case
            write(logfhandle,'(a)') '>>> APPLYING CORR FILTER TO VOLUME'
            call img_atom%new(ldim, smpd)
            call imgproc%new(ldim,smpd)
            call imgresult%new(ldim,smpd)
            call imgresult%set_ft(.true.)
            cutoff = 12.*smpd
            call imgproc%read(fname2process)
            call imgproc%fft
            call atom%new(1)
            call atom%set_element(1,element)
            call atom%set_coord(1,smpd*(real(ldim)/2.)) !DO NOT NEED THE +1
            call atom%convolve(img_atom, cutoff)
            call img_atom%fft
            call imgproc%phase_corr(img_atom,imgresult,lp)
            call imgresult%write(fname)
        endif
        call imgproc%kill
        call img_atom%kill
        call imgresult%kill
    end subroutine corrfilt_imgfile_2

    !>  \brief  is for applying CTF
    subroutine apply_ctf_imgfile( fname2process, fname, o, smpd, mode, bfac )
        use simple_oris,  only: oris
        use simple_ctf,   only: ctf
        character(len=*), intent(in)    :: fname2process, fname
        class(oris),      intent(inout) :: o
        real,             intent(in)    :: smpd
        character(len=*), intent(in)    :: mode
        real, optional,   intent(in)    :: bfac
        type(image)   :: img
        integer       :: i, n, ldim(3)
        type(ctf)     :: tfun
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3)    = 1
        call raise_exception_imgfile( n, ldim, 'apply_ctf_imgfile' )
        ! do the work
        if( n /= o%get_noris() ) THROW_HARD('inconsistent nr entries; apply_ctf_imgfile')
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> APPLYING CTF TO IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            tfun = ctf(smpd, o%get(i,'kv'), o%get(i,'cs'), o%get(i,'fraca'))
            if( o%isthere('dfy') )then ! astigmatic CTF
                call tfun%apply(img, o%get(i,'dfx'), mode, dfy=o%get(i,'dfy'), angast=o%get(i,'angast'), bfac=bfac)
            else ! non-astigmatic CTF
                call tfun%apply(img, o%get(i,'dfx'), mode, bfac=bfac)
            endif
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine apply_ctf_imgfile

    !>  \brief  is for origin shifting and rotating a stack
    !!          according to info in o
    !! \param fname2process  output filename
    !! \param fname  input filename
    !! \param o orientation object used for shifting and rotating
    !! \param mul multiplier
    subroutine shrot_imgfile( fname2process, fname, o, smpd, mul )
        use simple_oris,  only: oris
        character(len=*), intent(in)    :: fname2process, fname
        class(oris),      intent(inout) :: o
        real,             intent(in)    :: smpd
        real, optional,   intent(in)    :: mul
        type(image)   :: img, img_rot
        integer       :: n, i, ldim(3)
        real          :: x, y
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'shrot_imgfile' )
        ! do the work
        if( n /= o%get_noris() )then
            write(logfhandle,*) 'nr of entries in stack: ', n
            write(logfhandle,*) 'nr of entries in oris object: ', o%get_noris()
            THROW_HARD('inconsistent nr entries; shrot_imgfile')
        endif
        call img%new(ldim,smpd)
        call img_rot%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SHIFTING AND ROTATING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            call img%fft()
            x = o%get(i, 'x')
            y = o%get(i, 'y')
            if( present(mul) )then
                call img%shift([-x*mul,-y*mul,0.])
            else
                call img%shift([-x,-y,0.])
            endif
            call img%ifft()
            call img%rtsq(-o%e3get(i), 0., 0., img_rot)
            call img_rot%write(fname, i)
        end do
        call img%kill
        call img_rot%kill
    end subroutine shrot_imgfile

    !>  \brief  is for applying a soft circular mask to all images in stack
    !! \param fname2mask  output filename
    !! \param fname  input filename
    !! \param mskrad masking radius in pixels
    !! \param smpd sampling distance
    !! \param inner,width masking parameters
    !! \param which  file type
    subroutine mask_imgfile( fname2mask, fname, mskrad, smpd, inner, width, which )
        character(len=*),           intent(in) :: fname2mask, fname
        real,                       intent(in) :: mskrad
        real,                       intent(in) :: smpd
        real,             optional, intent(in) :: inner, width
        character(len=*), optional, intent(in) :: which
        type(image)           :: img
        integer               :: n, i, ldim(3)
        character(len=STDLEN) :: msktype
        msktype = 'soft'
        if( present(which) )msktype=which
        call find_ldim_nptcls(fname2mask, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'mask_imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> MASKING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2mask, i)
            call img%norm()
            call img%mask(mskrad, which, inner=inner, width=width)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine mask_imgfile

    !>  \brief  is for tapering edges of all images in stack
    !! \param fname2mask  output filename
    !! \param fname  input filename
    !! \param smpd sampling distance
    subroutine taper_edges_imgfile( fname2mask, fname, smpd)
        character(len=*),           intent(in) :: fname2mask, fname
        real,                       intent(in) :: smpd
        type(image)           :: img
        integer               :: n, i, ldim(3)
        call find_ldim_nptcls(fname2mask, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'taper_edges_imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> TAPERING EDGES OF IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2mask, i)
            call img%taper_edges
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine taper_edges_imgfile

    !>  \brief  is for binarizing all images in the stack
    !! \param fname2process  output filename
    !! \param fname  input filename
    !! \param thres binary threshold
    !! \param smpd sampling distance
    subroutine bin_imgfile( fname2process, fname, smpd, thres )
        use simple_segmentation, only: otsu_robust_fast
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd
        real, optional,   intent(in) :: thres
        type(image) :: img
        integer     :: n, i, ldim(3)
        real        :: thresh(3)
        logical     :: didft
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'bin_imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> BINARIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            didft = .false.
            if( img%is_ft() )then
                call img%ifft()
                didft = .true.
            endif
            if( present(thres) )then
                call img%binarize(thres)
            else
                call otsu_robust_fast(img, is2d=.true., noneg=.false., thresh=thresh)
            endif
            if( didft ) call img%fft()
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine bin_imgfile

    !>  random_selection_from_imgfile
    subroutine random_selection_from_imgfile( spproj, fname, box, nran )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: fname
        integer,           intent(in)    :: nran, box
        type(ran_tabu)                   :: rt
        type(image)                      :: img, img_scaled
        character(len=:), allocatable    :: stkname
        logical,          allocatable    :: mask(:)
        real    :: smpd
        integer :: nptcls, box_ori, ldim(3), i, ii, ldim_scaled(3), ind
        logical :: doscale
        ! dimensions
        nptcls      = spproj%get_nptcls()
        smpd        = spproj%os_stk%get(1,'smpd')
        box_ori     = nint(spproj%os_stk%get(1,'box'))
        ldim        = [box_ori,box_ori,1]
        ldim_scaled = [box,box,1]
        doscale     = box /= box_ori
        if( doscale )call img_scaled%new(ldim_scaled,smpd) ! this sampling distance will be overwritten
        call raise_exception_imgfile( nptcls, ldim, 'random_selection_from_imgfile' )
        call img%new(ldim,smpd)
        ! state mask
        if( spproj%os_ptcl2D%isthere('state') )then
            mask = spproj%os_ptcl2D%get_all('state') > 0.5
        else
            allocate(mask(nptcls), source=.true.)
        endif
        if( count(mask) < nran )THROW_HARD('Insufficient images; random_selection_from_imgfile')
        rt = ran_tabu(nptcls)
        do i=1,nptcls
            if( .not. mask(i) ) call rt%insert(i)
        end do
        ! copy
        write(logfhandle,'(a)') '>>> RANDOMLY SELECTING IMAGES'
        do i = 1,nran
            call progress(i, nran)
            ii = rt%irnd()
            call rt%insert(ii)
            call spproj%get_stkname_and_ind('ptcl2D', ii, stkname, ind)
            call img%read(stkname, ind)
            call img%norm()
            if( doscale )then
                call img%fft()
                if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2)&
                     .and. ldim_scaled(3) <= ldim(3) )then
                    call img%clip(img_scaled)
                else
                    call img%pad(img_scaled)
                endif
                call img_scaled%ifft()
                call img_scaled%write(fname, i)
            else
                call img%write(fname, i)
            endif
        end do
        call rt%kill
        call img%kill
        call img_scaled%kill
    end subroutine random_selection_from_imgfile

    !>  selection from imgfile of time series
    subroutine selection_from_tseries_imgfile( spproj, fname, box, nsel )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: fname
        integer,           intent(in)    :: nsel, box
        type(image)                      :: img, img_scaled
        character(len=:), allocatable    :: stkname
        real,             allocatable    :: states(:)
        integer,          allocatable    :: inds(:), inds_packed(:)
        real    :: smpd
        integer :: nptcls, box_ori, ldim(3), i, n
        integer :: ldim_scaled(3), ind, cnt, stepsz
        logical :: doscale
        ! dimensions
        nptcls      = spproj%get_nptcls()
        smpd        = spproj%os_stk%get(1,'smpd')
        box_ori     = nint(spproj%os_stk%get(1,'box'))
        ldim        = [box_ori,box_ori,1]
        ldim_scaled = [box,box,1]
        doscale     = box /= box_ori
        if( doscale )call img_scaled%new(ldim_scaled,smpd) ! this sampling distance will be overwritten
        call raise_exception_imgfile( nptcls, ldim, 'selection_from_tseries_imgfile' )
        call img%new(ldim,smpd)
        ! copy
        write(logfhandle,'(a)') '>>> SELECTING IMAGES FROM TIME-SERIES'
        ! extract the indices of the nonzero states
        allocate(inds(nptcls), source=(/(i,i=1,nptcls)/))
        if( spproj%os_ptcl2D%get_noris() == nptcls )then
            if( spproj%os_ptcl2D%isthere('state') )then
                states = spproj%os_ptcl2D%get_all('state')
            endif
        endif
        if( .not. allocated(states) ) allocate(states(nptcls), source=1.0)
        inds_packed = pack(inds, mask=states > 0.5)
        n      = size(inds_packed)
        stepsz = n / nsel
        cnt    = 0
        ! selection
        do i = 1,n,stepsz
            cnt = cnt + 1
            if( cnt > nsel ) exit
            call spproj%get_stkname_and_ind('ptcl2D', inds_packed(i), stkname, ind)
            call img%read(stkname, ind)
            call img%norm()
            if( doscale )then
                call img%fft()
                if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2)&
                     .and. ldim_scaled(3) <= ldim(3) )then
                    call img%clip(img_scaled)
                else
                    call img%pad(img_scaled)
                endif
                call img_scaled%ifft()
                call img_scaled%write(fname, cnt)
            else
                call img%write(fname, cnt)
            endif
        end do
        deallocate(states, inds, inds_packed)
        call img%kill
        call img_scaled%kill
    end subroutine selection_from_tseries_imgfile

    !>  random_selection_from_imgfile
    !! \param fname2selfrom  output filename
    !! \param fname  input filename
    !! \param ncls number of classes
    subroutine random_cls_from_imgfile( spproj, fname, ncls )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: fname
        integer,           intent(in)    :: ncls
        type(image)                      :: img, img_avg
        character(len=:), allocatable    :: stkname
        real               :: smpd
        integer            :: nptcls, ldim(3), i, j, ii, ind, box
        integer, parameter :: navg = 100
        ! dimensions
        nptcls = spproj%get_nptcls()
        smpd   = spproj%os_stk%get(1,'smpd')
        box    = nint(spproj%os_stk%get(1,'box'))
        ldim   = [box,box,1]
        call raise_exception_imgfile( nptcls, ldim, 'random_selection_from_imgfile' )
        call img%new(ldim,smpd)
        call img_avg%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> RANDOMLY GENERATING CLUSTER CENTERS'
        do i = 1,ncls
            call progress(i, ncls)
            img_avg = 0.
            do j = 1,navg
                ii = irnd_uni(nptcls)
                call spproj%get_stkname_and_ind('ptcl2D', ii, stkname, ind)
                call img%read(stkname, ind)
                call img%norm()
                call img_avg%add(img)
            enddo
            call img_avg%norm()
            call img_avg%write(fname, i)
        enddo
        call img%kill
        call img_avg%kill
    end subroutine random_cls_from_imgfile

end module simple_procimgfile

! stack image processing routines for SPIDER/MRC files
module simple_procimgfile
include 'simple_lib.f08'
use simple_image,   only: image
use simple_oris,    only: oris
implicit none

private ! :: raise_exception_imgfile

public :: copy_imgfile, diff_imgfiles, make_pattern_stack, pad_imgfile, resize_imgfile, clip_imgfile
public :: resize_and_clip_imgfile, resize_imgfile_double, norm_bin_imgfile, norm_imgfile, norm_ext_imgfile
public :: noise_norm_imgfile, shellnorm_imgfile, stats_imgfile, neg_imgfile, masscen_imgfile, cure_imgfile
public :: acf_imgfile, add_noise_imgfile, frameavg_imgfile, ft2img_imgfile,shift_imgfile,bp_imgfile
public :: real_filter_imgfile,phase_rand_imgfile, apply_ctf_imgfile,shrot_imgfile,mask_imgfile, taper_edges_imgfile
public :: bin_imgfile, make_avg_imgfile, random_selection_from_imgfile

#include "simple_local_flags.inc"
contains

    !>  \brief  is for raising exception
    subroutine raise_exception_imgfile( n, ldim, routine )
        integer, intent(in) :: n        !< num of images
        integer, intent(in) :: ldim(3)  !< logical dimensions
        character(len=*)    :: routine  !< Error message caller
        if( n < 1 .or. any(ldim == 0) )then
            write(*,*) routine
            write(*,*) 'The input stack is corrupt!'
            write(*,*) 'Number of images: ', n
            write(*,*) 'Logical dimensions: ', ldim
            call simple_stop ('procimgfile exception')
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
        write(*,'(a)') '>>> COPYING IMAGES'
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
            call simple_stop ('procimgfile exception ; diff_imgfiles mismatch')
        call diffimg%new(ldim,smpd)
        call img2%new(ldim,smpd)
        if( n >= 1 )then
            write(*,'(a)') '>>> SUBTRACTING IMAGES'
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

    !>  \brief  is for making a stack of normalized vectors for PCA analysis
    !! \param fnameStack,fnamePatterns filenames for stacka and pattern
    subroutine make_pattern_stack( fnameStack, fnamePatterns, mskrad, D, recsz, avg, otab, hfun )
        character(len=*),            intent(in)    :: fnameStack, fnamePatterns
        real,                        intent(in)    :: mskrad       !< mask radius
        integer,                     intent(out)   :: D, recsz     !< record size
        real, allocatable, optional, intent(out)   :: avg(:)       !< frame stack average
        class(oris),       optional, intent(inout) :: otab         !< oris table
        character(len=*),  optional, intent(in)    :: hfun         !< which normalise pca vec
        type(image)        :: img
        real, allocatable  :: pcavec(:)
        real               :: x, y
        integer            :: n, fnum, ier, i, ldim(3)
        logical            :: err

        call find_ldim_nptcls(fnameStack, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'make_pattern_stack' )
        ! build and initialise objects
        call img%new(ldim,1.)
        D = img%get_npix(mskrad)
        allocate(pcavec(D), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('make_pattern_stack; simple_procimgfile, 1', alloc_stat)
        pcavec = 0.
        inquire(iolength=recsz) pcavec
        deallocate(pcavec)
        if( present(avg) )then
            allocate(avg(D), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('make_pattern_stack; simple_procimgfile, 2', alloc_stat)
            avg = 0.
        endif
        ! extract patterns and write to file
        call fopen(fnum, status='replace', action='readwrite', file=fnamePatterns,&
             access='direct', form='unformatted', recl=recsz, iostat=ier)
        call fileiochk('make_pattern_stack; simple_procimgfile', ier)
        write(*,'(a)') '>>> MAKING PATTERN STACK'
        do i=1,n
            call progress(i,n)
            call img%read(fnameStack, i)
            if( present(otab) )then
                ! shift image
                call img%fft()
                x = otab%get(i, 'x')
                y = otab%get(i, 'y')
                call img%shift([-x,-y,0.])
                ! rotate image
                call img%ifft()
                call img%rtsq(-otab%e3get(i), 0., 0.)
            else
                if( img%is_ft() ) call img%ifft()
            endif
            call img%serialize(pcavec, mskrad)
            err = .false.
            if( present(hfun) )then
                call normalize_sigm(pcavec)
            else
                call normalize(pcavec, err)
            endif
            if( err ) write(*,'(a,i7)') 'WARNING: variance zero! image nr: ', i
            if( present(avg) ) avg = avg+pcavec
            write(fnum,rec=i) pcavec
            if( debug .or. global_debug )then
                call check4nans(pcavec)
                call img%serialize(pcavec, mskrad)
                call img%write('unserialized.spi', i)
            endif
            deallocate(pcavec)
        end do
        if( present(avg) )then
            avg = avg/real(n)
            allocate(pcavec(D), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('make_pattern_stack; simple_procimgfile, 3', alloc_stat)
            do i=1,n
                read(fnum,rec=i) pcavec
                pcavec = pcavec-avg
                write(fnum,rec=i) pcavec
            end do
            deallocate(pcavec)
        endif
        call fclose(fnum,errmsg='make_pattern_stack; simple_procimgfile')
        call img%kill
    end subroutine make_pattern_stack

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
            write(*,'(a)') '>>> PADDING IMAGES'
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
        type(image) :: img, img_resized
        integer     :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2resize, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'resize_imgfile' )
        call img%new(ldim,smpd)
        call img_resized%new(ldim_new,smpd) ! this sampling distance will be overwritten
        write(*,'(a)') '>>> RESIZING IMAGES'
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
        call img%kill
        call img_resized%kill
    end subroutine resize_imgfile

    !>  \brief clip_imgfile is for clipping
    !! \param fname2clip output filename
    !! \param fname input filename
    !! \param ldim_clip clipped logical dimension
    !! \param smpd sampling distance
    !!
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
            write(*,'(a)') '>>> CLIPPING IMAGES'
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
        write(*,'(a)') '>>> RESIZING IMAGES'
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

    !>  \brief  is for resizing and clipping
     !! \param fname2resize  output filename
    !! \param fname1,fname2  input filenames
    !! \param smpd  sampling distance
    !! \param ldims_new new logical dimension
    !! \param smpds_new new sampling distance
    !! \param fromptop  index range for copying
    !!
    subroutine resize_imgfile_double( fname2resize, fname1, fname2, smpd, ldims_new, smpds_new, fromptop )
        character(len=*),  intent(in)  :: fname2resize, fname1, fname2
        real,              intent(in)  :: smpd
        integer,           intent(in)  :: ldims_new(2,3)
        real,              intent(out) :: smpds_new(2)
        integer, optional, intent(in)  :: fromptop(2)
        type(image) :: img, img_resized1, img_resized2
        integer     :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2resize, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'resize_imgfile' )
        call img%new(ldim,smpd)
        call img_resized1%new(ldims_new(1,:),smpd) ! this sampling distance will be overwritten
        call img_resized2%new(ldims_new(2,:),smpd) ! this sampling distance will be overwritten
        write(*,'(a)') '>>> RESIZING IMAGES'
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
            call img%clip(img_resized1)
            call img%clip(img_resized2)
            call img_resized1%ifft()
            call img_resized2%ifft()
            call img_resized1%write(fname1, cnt)
            call img_resized2%write(fname2, cnt)
        end do
        smpds_new(1) = img_resized1%get_smpd()
        smpds_new(2) = img_resized2%get_smpd()
        call img%kill
        call img_resized1%kill
        call img_resized2%kill
    end subroutine resize_imgfile_double

    subroutine norm_bin_imgfile( fname2norm, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: smpd
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> BIN NORMALIZING IMAGES'
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
        write(*,'(a)') '>>> NORMALIZING IMAGES'
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
        write(*,'(a)') '>>> NORMALIZING IMAGES'
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
        type(image)      :: img
        integer          :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'noise_norm_imgfile' )
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> NOISE NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%noise_norm(msk)
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
        write(*,'(a)') '>>> SHELL NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%shellnorm()
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine shellnorm_imgfile

    !>  \brief  is for calculating image statistics
    subroutine stats_imgfile( fname, os, msk )
        use simple_oris, only: oris
        character(len=*), intent(in)  :: fname
        class(oris),      intent(out) :: os
        real, optional,   intent(in)  :: msk
        real              :: ave, sdev, med, minv, maxv, spec
        real, allocatable :: spectrum(:)
        type(image)       :: img
        integer           :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'stats_imgfile' )
        call img%new(ldim,1.)
        call os%new_clean(n)
        write(*,'(a)') '>>> CALCULATING STACK STATISTICS'
        do i=1,n
            call progress(i,n)
            call img%read(fname,i)
            call img%stats('foreground', ave, sdev, maxv, minv, med=med, msk=msk)
            call img%fft()
            !$ allocate(spectrum(1))
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
        write(*,'(a)') '>>> INVERTING IMAGE CONTRAST'
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
        write(*,'(a)') '>>> CENTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2masscen, i)
            xyz = img%center(lp, msk)
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
        write(*,'(a)') '>>> CURING IMAGES'
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

    !>  \brief  is for calculating the acf of a stack
    !! \param fname2acf  output filename
    !! \param fname  input filename
    subroutine acf_imgfile( fname2acf, fname )
        character(len=*), intent(in) :: fname2acf, fname
        type(image)   :: img
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2acf, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'acf_imgfile' )
        call img%new(ldim,1.)
        write(*,'(a)') '>>> CALCULATING ACF:S OF THE IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2acf, i)
            call img%acf
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine acf_imgfile

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
        write(*,'(a)') '>>> ADDING NOISE TO THE IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            call img%add_gauran(snr)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine add_noise_imgfile

    !>  \brief  is for making frame averages of dose-fractionated image series
    !! \param fname2process  output filename
    !! \param fname  input filename
    !! \param smpd  sampling distance
    !! \param navg number of averages
    subroutine frameavg_imgfile( fname2process, fname, navg, smpd )
        character(len=*), intent(in) :: fname2process, fname
        integer,          intent(in) :: navg
        real,             intent(in) :: smpd
        type(image) :: img, avg
        integer     :: i, n, cnt, cnt2, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'frameavg_imgfile' )
        call img%new(ldim,smpd)
        call avg%new(ldim,smpd)
        cnt = 0
        cnt2 = 0
        write(*,'(a)') '>>> AVERAGING FRAMES'
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
    end subroutine frameavg_imgfile

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
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; shift; simple_procimgfile'
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> SHIFTING IMAGES'
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
    !> bp_imgfile
    !! \param fname2filter output filename
    !! \param fname input filename
    !! \param smpd sampling distance
    !! \param hp high-pass cutoff freq
    !! \param lp low-pass cutoff freq
    !!
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
        write(*,'(a)') '>>> BAND-PASS FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filter, i)
            call img%bp(hp,lp,width=width)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine bp_imgfile

    !>  \brief real_filter_imgfile is for real-space filtering
    !! \param fname2filter input filename
    !! \param fname output filename
    !! \param smpd sampling distance
    !! \param which_filter filter tag
    !! \param winsz window size
    !!
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
        write(*,'(a)') '>>> REAL-SPACE FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2filter, i)
            call img%real_space_filter(winsz, which_filter)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine real_filter_imgfile

    !>  \brief  is for phase randomization
    !! \param fname2process output filename
    !! \param fname input filename
    !! \param smpd sampling distance
    !! \param lp low-pass cutoff freq
    !!
    subroutine phase_rand_imgfile( fname2process, fname, smpd, lp )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, lp
        type(image) :: img
        integer :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'phase_rand_imgfile' )
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> PHASE RANDOMIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            call img%phase_rand(lp)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine phase_rand_imgfile

    !>  \brief  is for applying CTF
    !! \param fname2process output filename
    !! \param fname input filename
    !! \param smpd sampling distance
    !! \param mode,bfac ctf mode and bfac args
    !! \param o orientation object used for ctf
    !!
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
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; apply_ctf; simple_procimgfile'
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> APPLYING CTF TO IMAGES'
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
            write(*,*) 'nr of entries in stack: ', n
            write(*,*) 'nr of entries in oris object: ', o%get_noris()
            stop 'inconsistent nr entries; shrot; simple_procimgfile'
        endif
        call img%new(ldim,smpd)
        call img_rot%new(ldim,smpd)
        write(*,'(a)') '>>> SHIFTING AND ROTATING IMAGES'
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
    !! \param mskrad masking radius
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
        write(*,'(a)') '>>> MASKING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2mask, i)
            call img%norm()
            call img%mask(mskrad, which, inner, width)
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
        write(*,'(a)') '>>> TAPERING EDGES OF IMAGES'
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
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd
        real, optional,   intent(in) :: thres
        type(image) :: img
        integer     :: n, i, ldim(3)
        logical     :: didft
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'bin_imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> BINARIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            didft = .false.
            if( img%is_ft() )then
                call img%ifft()
                didft = .true.
            endif
            if( present(thres) )then
                call img%bin(thres)
            else
                call img%bin_kmeans
            endif
            if( didft ) call img%fft()
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine bin_imgfile

    !>   make_avg_imgfile is for making an average of the entire stack
    !! \param fname stack filename
    !! \param avgname output filename
    subroutine make_avg_imgfile( fname, avgname, smpd )
        character(len=*), intent(in) :: fname, avgname
        real,             intent(in) :: smpd              !< sampling distance, resolution
        type(image)                  :: avg, img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'make_avg_imgfile' )
        call avg%new(ldim, smpd)
        call img%new(ldim, smpd)
        avg = 0.
        write(*,'(a)') '>>> MAKING GLOBAL STACK AVERAGE'
        do i=1,n
            call progress(i,n)
            call img%read(fname, i)
            call avg%add(img)
        end do
        call avg%div(real(n))
        call avg%write(avgname,1)
    end subroutine make_avg_imgfile

    !>  random_selection_from_imgfile
    !! \param fname2selfrom  output filename
    !! \param fname  input filename
    !! \param nran number of random samples
    !! \param smpd sampling distance
    subroutine random_selection_from_imgfile( spproj, fname, box, nran)
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
        if( count(mask) < nran )stop 'Insufficient images; simple_procimgfile%random_selection_from_imgfile'
        rt = ran_tabu(nptcls)
        do i=1,nptcls
            if( .not. mask(i) ) call rt%insert(i)
        end do
        ! copy
        write(*,'(a)') '>>> RANDOMLY SELECTING IMAGES'
        do i = 1,nran
            call progress(i, nran)
            ii = rt%irnd()
            call rt%insert(ii)
            call spproj%get_stkname_and_ind('ptcl2D', ii, stkname, ind)
            call img%read(stkname, ind)
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

end module simple_procimgfile

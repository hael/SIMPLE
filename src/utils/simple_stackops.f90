! stack image processing routines for SPIDER/MRC files
module simple_stackops
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: make_pcavec_stack, prep_imgfile4movie
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
    !! \param fnameStack,fnamePCAvecs filenames for stack and pattern
    subroutine make_pcavec_stack( fnameStack, fnamePCAvecs, D, recsz, avg, l_mask )
        class(string),     intent(in)  :: fnameStack, fnamePCAvecs
        integer,           intent(out) :: D, recsz      !< record size
        real, allocatable, intent(out) :: avg(:)        !< frame stack average
        logical, optional, intent(in)  :: l_mask(:,:,:) !< true for pixels to extract
        type(image)          :: img
        real,    allocatable :: pcavec(:)
        logical, allocatable :: ll_mask(:,:,:)
        integer              :: n, fnum, ier, i, ldim(3), ldim_mask(3)
        call find_ldim_nptcls(fnameStack, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'make_pcavec_stack' )
        if( present(l_mask) )then
            ldim_mask = [size(l_mask, dim=1),size(l_mask, dim=2),size(l_mask, dim=3)]
            if( .not. all(ldim_mask .eq. ldim) )then
                write(logfhandle,*) 'ERROR! nonconforming matrix sizes'
                write(logfhandle,*) 'ldim of image: ', ldim
                write(logfhandle,*) 'ldim of mask : ', ldim_mask
                THROW_HARD('make_pcavec_stack')
            endif
            allocate(ll_mask(ldim_mask(1),ldim_mask(2),ldim_mask(3)), source=l_mask)
        else
            ldim_mask = ldim
            allocate(ll_mask(ldim_mask(1),ldim_mask(2),ldim_mask(3)), source=.true.)
        endif
        ! build and initialise objects
        call img%new(ldim,1.)
        D = count(ll_mask)
        allocate(pcavec(D), source=0.)
        inquire(iolength=recsz) pcavec
        deallocate(pcavec)
        if( allocated(avg) ) deallocate(avg)
        allocate(avg(D), source=0.)
        ! extract patterns and write to file
        call fopen(fnum, status='replace', action='readwrite', file=fnamePCAvecs,&
             access='direct', form='unformatted', recl=recsz, iostat=ier)
        call fileiochk('make_pcavec_stack; simple_procimgstk', ier)
        write(logfhandle,'(a)') '>>> MAKING PCAVEC STACK'
        do i=1,n
            call progress(i,n)
            call img%read(fnameStack, i)
            pcavec = img%serialize(ll_mask)
            avg = avg + pcavec
            write(fnum,rec=i) pcavec
            deallocate(pcavec)
        end do
        avg = avg/real(n)
        allocate(pcavec(D))
        do i=1,n
            read(fnum,rec=i) pcavec
            pcavec = pcavec - avg
            write(fnum,rec=i) pcavec
        end do
        deallocate(pcavec)
        call fclose(fnum)
        call img%kill
    end subroutine make_pcavec_stack

    !>   make_avg_imgfile is for making an average of the entire stack
    !! \param fname stack filename
    !! \param avgname output filename
    subroutine make_avg_stack( fname, avgname, smpd )
        class(string), intent(in) :: fname, avgname
        real,          intent(in) :: smpd              !< sampling distance, resolution
        type(image) :: avg, img
        integer     :: i, n, ldim(3)
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
        class(string), intent(in) :: fname2process, fname
        integer,       intent(in) :: navg
        real,          intent(in) :: smpd
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
        class(string), intent(in) :: fname2acf, fname
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
        class(string),  intent(in)  :: fname
        class(oris),    intent(out) :: os
        real, optional, intent(in)  :: msk
        real              :: ave, sdev, med, minv, maxv, spec
        real, allocatable :: spectrum(:)
        type(image)       :: img
        integer           :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception( n, ldim, 'stats_imgfile' )
        call img%new(ldim,1.)
        call os%new(n, is_ptcl=.false.)
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

    !>  \brief  is for calculating image statistics
    subroutine prep_imgfile4movie( fname, smpd )
        class(string), intent(in) :: fname
        real,          intent(in) :: smpd
        type(string) :: fname_out, str_cat, str_mp4
        type(image)  :: img
        integer      :: i, ldim(3), numlen, funit, iostat, nptcls
        call find_ldim_nptcls(fname, ldim, nptcls)
        ldim(3) = 1
        call img%new(ldim, smpd)
        ! jpeg output
        numlen = len(int2str(nptcls))
        do i=1,nptcls
            call progress(i, nptcls)
            fname_out = fname_new_ext(basename(fname), string(int2str_pad(i,numlen)//'.jpg'))
            call img%read(fname, i)
            call img%write_jpg(fname_out, quality=100, norm=.false.)
        end do
        call img%kill
        ! script ouput
        str_cat = fname_new_ext(basename(fname),string('*.jpg'))
        str_mp4 = fname_new_ext(basename(fname),string('mp4'))
        call fopen(funit, status='replace', action='write', file=string('makemovie.script'),iostat=iostat)
        write(funit,*)'# Increase the argument to the first occurence of argument -r:v for a faster movie'
        write(funit,*)'cat '//str_cat%to_char()//' | ffmpeg -y -f image2pipe -vcodec mjpeg -r:v 5 '&
            &//'-i - -vcodec:v libx264 -preset veryslow -r:v 25 '//str_mp4%to_char()
        call fclose(funit)
        iostat = simple_chmod('makemovie.script','+x')
        call fname_out%kill
        call str_cat%kill
        call str_mp4%kill
        write(logfhandle,'(A)')'>>> Execute the generated script: ./makemovie.script'
    end subroutine prep_imgfile4movie

end module simple_stackops

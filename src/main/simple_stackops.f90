! stack image processing routines for SPIDER/MRC files
module simple_stackops
include 'simple_lib.f08'
use simple_image,   only: image
use simple_oris,    only: oris
implicit none
private
#include "simple_local_flags.inc"

public :: make_pattern_stack, acf_stack, make_avg_stack, stats_imgfile, frameavg_stack
contains

    !>  \brief  is for raising exception
    subroutine raise_exception( n, ldim, routine )
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
    end subroutine raise_exception

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
        call raise_exception( n, ldim, 'make_pattern_stack' )
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
        write(*,'(a)') '>>> MAKING GLOBAL STACK AVERAGE'
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
        write(*,'(a)') '>>> CALCULATING ACF:S OF THE IMAGES'
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
        write(*,'(a)') '>>> CALCULATING STACK STATISTICS'
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

end module simple_stackops

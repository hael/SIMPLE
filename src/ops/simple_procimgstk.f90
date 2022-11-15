! stack image processing routines for SPIDER/MRC files
module simple_procimgstk
include 'simple_lib.f08'
use simple_image,    only: image
use simple_stack_io, only: stack_io
use simple_tvfilter
implicit none
#include "simple_local_flags.inc"

contains

    subroutine raise_exception_imgfile( n, ldim, routine )
        integer, intent(in) :: n
        integer, intent(in) :: ldim(3)
        character(len=*)    :: routine
        if( n < 1 .or. any(ldim == 0) )then
            write(logfhandle,*) routine
            write(logfhandle,*) 'The input stack is corrupt!'
            write(logfhandle,*) 'Number of images: ', n
            write(logfhandle,*) 'Logical dimensions: ', ldim
            THROW_HARD('procimgfile exception')
        endif
    end subroutine raise_exception_imgfile

    subroutine copy_imgfile( fname2copy, fname, smpd, fromto )
        character(len=*), intent(in) :: fname2copy, fname
        real,             intent(in) :: smpd
        integer,          intent(in) :: fromto(2)
        type(stack_io) :: stkio_r, stkio_w
        type(image)    :: img
        integer        :: n, i, cnt, ldim(3)
        call find_ldim_nptcls(fname2copy, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile(n, ldim, 'copy_imgfile')
        call stkio_r%open(fname2copy, smpd, 'read')
        call stkio_w%open(fname, smpd, 'write', box=ldim(1), is_ft=.false.)
        call img%new(ldim,smpd)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> COPYING IMAGES'
        cnt = 0
        do i=fromto(1),fromto(2)
            cnt = cnt+1
            call progress(cnt, fromto(2)-fromto(1)+1)
            call stkio_r%read(i, img)
            call stkio_w%write(cnt, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine copy_imgfile

    subroutine pad_imgfile( fname2pad, fname, ldim_pad, smpd )
        character(len=*), intent(in) :: fname2pad, fname
        integer,          intent(in) :: ldim_pad(3)
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w, stkio_w2
        type(image)                  :: img, img_pad
        integer                      :: n, i, ldim(3)
        real                         :: ave, sdev, maxv, minv, med
        call find_ldim_nptcls(fname2pad, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'pad_imgfile' )
        if( ldim_pad(1) >= ldim(1) .and. ldim(2) >= ldim(2) )then
            call stkio_r%open(fname2pad, smpd, 'read')
            call stkio_w%open(fname,     smpd, 'write', box=ldim_pad(1), is_ft=.false.)
            call stkio_w2%open(fname,    smpd, 'write', box=ldim_pad(1), is_ft=.true.)
            call img%new(ldim,smpd)
            call img_pad%new(ldim_pad,smpd)
            write(logfhandle,'(a)') '>>> PADDING IMAGES'
            do i=1,n
                call progress(i,n)
                call stkio_r%read(i, img)
                if( img%is_ft() )then
                    call img%pad(img_pad) ! FT state preserved
                    call stkio_w2%write(i, img_pad)
                else
                    ! get background statistics
                    call img%stats('background', ave, sdev, maxv, minv, med=med)
                    call img%pad(img_pad, backgr=med) ! FT state preserved
                    call stkio_w%write(i, img_pad)
                endif
            end do
            call stkio_r%close
            call stkio_w%close
            call stkio_w2%close
            call img%kill
            call img_pad%kill
        end if
    end subroutine pad_imgfile

    subroutine scale_imgfile( fname2scale, fname, smpd, ldim_new, smpd_new, fromptop )
        character(len=*),  intent(in)  :: fname2scale, fname
        real,              intent(in)  :: smpd
        integer,           intent(in)  :: ldim_new(3)
        real,              intent(out) :: smpd_new
        integer, optional, intent(in)  :: fromptop(2)
        type(stack_io) :: stkio_r, stkio_w
        type(image) :: img, img_scaled
        integer     :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2scale, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'scale_imgfile' )
        call stkio_r%open(fname2scale, smpd, 'read')
        call stkio_w%open(fname, smpd, 'write', box=ldim_new(1))
        call img%new(ldim,smpd,wthreads=.false.)
        call img_scaled%new(ldim_new,smpd,wthreads=.false.) ! this sampling distance will be overwritten
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> SCALING IMAGES'
        if( present(fromptop) )then
            prange = fromptop
        else
            prange(1) = 1
            prange(2) = n
        endif
        sz = prange(2) - prange(1) + 1
        if( sz == 0 )then
            THROW_WARN('no images to operate on given input range')
        else
            cnt = 0
            do i=prange(1),prange(2)
                cnt = cnt+1
                call progress(cnt,sz)
                call stkio_r%read(i, img)
                call img%fft()
                if( ldim_new(1) <= ldim(1) .and. ldim_new(2) <= ldim(2)&
                     .and. ldim_new(3) <= ldim(3) )then
                    call img%clip(img_scaled)
                else
                    call img%pad(img_scaled)
                endif
                call img_scaled%ifft()
                call stkio_w%write(cnt, img_scaled)
            end do
            smpd_new = img_scaled%get_smpd()
        endif
        call stkio_r%close
        call stkio_w%close
        call img%kill
        call img_scaled%kill
    end subroutine scale_imgfile

    subroutine clip_imgfile( fname2clip, fname, ldim_clip, smpd )
        character(len=*), intent(in) :: fname2clip, fname
        integer,          intent(in) :: ldim_clip(3)
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img, img_clip
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2clip, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'clip_imgfile' )
        if( ldim_clip(1) <= ldim(1) .and. ldim_clip(2) <= ldim(2) )then
            call stkio_r%open(fname2clip, smpd, 'read')
            call stkio_w%open(fname,      smpd, 'write', box=ldim_clip(1))
            call img%new(ldim,smpd)
            call img_clip%new(ldim_clip,smpd)
            write(logfhandle,'(a)') '>>> CLIPPING IMAGES'
            do i=1,n
                call progress(i,n)
                call stkio_r%read(i, img)
                call img%clip(img_clip) ! FT state preserved
                call stkio_w%write(i, img_clip)
            end do
            call stkio_r%close
            call stkio_w%close
            call img%kill
            call img_clip%kill
        end if
    end subroutine clip_imgfile

    subroutine roavg_imgfile( fname2roavg, fname, angstep, smpd )
        character(len=*), intent(in) :: fname2roavg, fname
        integer,          intent(in) :: angstep
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img, img_roavg
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2roavg, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'roavg_imgfile' )
        call stkio_r%open(fname2roavg, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        call img_roavg%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> ROTATIONALLY AVERAGING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%roavg(angstep, img_roavg)
            call stkio_w%write(i, img_roavg)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
        call img_roavg%kill
    end subroutine roavg_imgfile

    subroutine mirror_imgfile( fname2mirr, fname, mirr_flag, smpd )
        character(len=*), intent(in) :: fname2mirr, fname
        character(len=1), intent(in) :: mirr_flag
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2mirr, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'clip_imgfile' )
        call stkio_r%open(fname2mirr, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> MIRRORING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%mirror(mirr_flag)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine mirror_imgfile

    subroutine scale_and_clip_imgfile( fname2scale, fname, smpd, ldim_new, ldim_clip, smpd_new, fromptop )
        character(len=*),  intent(in)  :: fname2scale, fname
        real,              intent(in)  :: smpd
        integer,           intent(in)  :: ldim_new(3)
        integer,           intent(in)  :: ldim_clip(3)
        real,              intent(out) :: smpd_new
        integer, optional, intent(in)  :: fromptop(2)
        type(stack_io)                 :: stkio_r, stkio_w
        type(image)                    :: img, img_scaled, img_clip
        integer                        :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2scale, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'scale_and_clip_imgfile' )
        call stkio_r%open(fname2scale, smpd, 'read')
        call stkio_w%open(fname, smpd, 'write', box=ldim_clip(1))
        call img%new(ldim,smpd)
        call img_scaled%new(ldim_new,smpd) ! this sampling distance will be overwritten
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> SCALING IMAGES'
        if( present(fromptop) )then
            prange = fromptop
        else
            prange(1) = 1
            prange(2) = n
        endif
        sz = prange(2) - prange(1) + 1
        if( sz == 0 )then
            THROW_WARN('no images to operate on given input range')
        else
            cnt = 0
            do i=prange(1),prange(2)
                cnt = cnt+1
                call progress(cnt,sz)
                call stkio_r%read(i, img)
                call img%fft()
                call img%clip(img_scaled)
                call img_scaled%ifft()
                call img_clip%new(ldim_clip,img_scaled%get_smpd())
                if( ldim_clip(1) <= ldim_new(1) .and. ldim_clip(2) <= ldim_new(2) )then
                    call img_scaled%clip(img_clip)
                else
                    call img_scaled%pad(img_clip)
                endif
                call stkio_w%write(cnt, img_clip)
            end do
            smpd_new = img_scaled%get_smpd()
        endif
        call stkio_r%close
        call stkio_w%close
        call img%kill
        call img_scaled%kill
        call img_clip%kill
    end subroutine scale_and_clip_imgfile

    subroutine norm_imgfile( fname2norm, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        call stkio_r%open(fname2norm, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%norm()
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine norm_imgfile

    subroutine noise_norm_imgfile( fname2norm, msk, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: msk, smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: i, n, ldim(3)
        logical, allocatable         :: lmsk(:,:,:)
        real                         :: sdev_noise
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'noise_norm_imgfile' )
        call stkio_r%open(fname2norm, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%disc(ldim, smpd, msk, lmsk)
        write(logfhandle,'(a)') '>>> NOISE NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%noise_norm(lmsk, sdev_noise)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine noise_norm_imgfile

    subroutine shellnorm_imgfile( fname2norm, fname, smpd )
        character(len=*), intent(in) :: fname2norm, fname
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'shellnorm_imgfile' )
        call stkio_r%open(fname2norm, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> SHELL NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%shellnorm()
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine shellnorm_imgfile

    subroutine matchfilt_imgfile( fname2filt, fname, frcs_fname, smpd )
        use simple_class_frcs, only: class_frcs
        character(len=*), intent(in) :: fname2filt, fname, frcs_fname
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(class_frcs)             :: frcs
        type(image)                  :: img
        real, allocatable            :: frc(:), filter(:)
        integer                      :: i, n, ldim(3)
        call frcs%read(frcs_fname)
        call find_ldim_nptcls(fname2filt, ldim, n)
        ldim(3) = 1
        if( frcs%get_nprojs().ne.n )then
            write(logfhandle,*) '# imgs: ',n
            write(logfhandle,*) '# class_frcs: ',frcs%get_nprojs()
            THROW_HARD('inconsistent dimensions; matchfilt_imgfile')
        endif
        call raise_exception_imgfile( n, ldim, 'matchfilt_imgfile' )
        call stkio_r%open(fname2filt, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        if( frcs%get_filtsz().ne.img%get_filtsz() )then
            write(logfhandle,*) 'img filtsz: ',img%get_filtsz()
            write(logfhandle,*) 'frcs filtsz: ',frcs%get_filtsz()
            THROW_HARD('Inconsistent filter dimensions; matchfilt_imgfile')
        endif
        write(logfhandle,'(a)') '>>> SHELL NORMALIZING AND FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            frc = frcs%get_frc(i, ldim(1))
            allocate(filter(size(frc)), source=1.)
            call fsc2optlp_sub(frcs%get_filtsz(), frc, filter)
            where( filter < TINY )filter = 0.
            call img%fft()
            call img%apply_filter(filter)
            call img%shellnorm
            call img%ifft()
            call stkio_w%write(i, img)
            deallocate(frc,filter)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine matchfilt_imgfile

    subroutine neg_imgfile( fname2neg, fname, smpd )
        character(len=*), intent(in) :: fname2neg, fname
        real,             intent(in) :: smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname2neg, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'neg_imgfile' )
        call stkio_r%open(fname2neg, smpd, 'read')
        call stkio_w%open(fname,     smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> INVERTING IMAGE CONTRAST'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%neg()
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine neg_imgfile

    subroutine masscen_imgfile( fname2masscen, fname, smpd, lp, msk )
        character(len=*), intent(in) :: fname2masscen, fname
        real,             intent(in) :: smpd, lp, msk
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: i, n, ldim(3)
        real                         :: xyz(3)
        call find_ldim_nptcls(fname2masscen, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'masscen_imgfile' )
        call stkio_r%open(fname2masscen, smpd, 'read')
        call stkio_w%open(fname,         smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> CENTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            xyz = img%calc_shiftcen(lp, msk)
            print *, xyz(:2)
            call img%shift(xyz)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine masscen_imgfile

    subroutine add_noise_imgfile( fname2process, fname, snr, smpd )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: snr, smpd
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'add_noise_imgfile' )
        call stkio_r%open(fname2process, smpd, 'read')
        call stkio_w%open(fname, smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> ADDING NOISE TO THE IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%add_gauran(snr)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine add_noise_imgfile

    subroutine bp_imgfile( fname2filter, fname, smpd, hp, lp, width )
        character(len=*), intent(in) :: fname2filter, fname
        real,             intent(in) :: smpd, hp, lp
        real, optional,   intent(in) :: width
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: n, i, ldim(3)
        real                         :: wwidth
        wwidth = 10.0
        if( present(width) ) wwidth = width
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'bp_imgfile' )
        call stkio_r%open(fname2filter, smpd, 'read')
        call stkio_w%open(fname,        smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> BAND-PASS FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%bp(hp,lp,width=width)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine bp_imgfile

    subroutine real_filter_imgfile( fname2filter, fname, smpd, which_filter, winsz )
        character(len=*), intent(in) :: fname2filter, fname
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: which_filter
        integer,          intent(in) :: winsz
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'real_filter_imgfile' )
        call stkio_r%open(fname2filter, smpd, 'read')
        call stkio_w%open(fname, smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> REAL-SPACE FILTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%real_space_filter(winsz, which_filter)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine real_filter_imgfile

    subroutine phase_rand_imgfile( fname2process, fname, smpd, lp )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, lp
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'phase_rand_imgfile' )
        call stkio_r%open(fname2process, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> PHASE RANDOMIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%phase_rand(lp)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine phase_rand_imgfile

    subroutine tvfilter_imgfile( fname2process, fname, smpd, lambda )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd, lambda
        type(stack_io)               :: stkio_r, stkio_w
        type(tvfilter)               :: tv
        type(image)                  :: img
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'tvfilter_imgfile' )
        call stkio_r%open(fname2process, smpd, 'read')
        call stkio_w%open(fname,         smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        call tv%new()
        write(logfhandle,'(a)') '>>> APPLYING TV FILTER TO IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call tv%apply_filter(img,lambda)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
        call tv%kill
    end subroutine tvfilter_imgfile

    subroutine nlmean_imgfile( fname2process, fname, smpd, noise_sdev )
        character(len=*), intent(in) :: fname2process, fname
        real,             intent(in) :: smpd
        real, optional,   intent(in) :: noise_sdev
        type(stack_io)               :: stkio_r, stkio_w
        type(image)                  :: img
        integer                      :: n, i, ldim(3)
        logical                      :: noise_sdev_present
        noise_sdev_present = present(noise_sdev)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'nlmean_imgfile' )
        call stkio_r%open(fname2process, smpd, 'read')
        call stkio_w%open(fname,         smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> APPLYING NLMEAN FILTER TO IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            if( noise_sdev_present )then
                call img%nlmean(sdev_noise=noise_sdev)
            else
                call img%nlmean
            endif
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine nlmean_imgfile

    subroutine apply_ctf_imgfile( fname2process, fname, o, smpd, mode, bfac )
        use simple_ctf,   only: ctf
        character(len=*), intent(in)    :: fname2process, fname
        class(oris),      intent(inout) :: o
        real,             intent(in)    :: smpd
        character(len=*), intent(in)    :: mode
        real, optional,   intent(in)    :: bfac
        type(stack_io)                  :: stkio_r, stkio_w
        type(image)                     :: img
        integer                         :: i, n, ldim(3)
        type(ctf)                       :: tfun
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3)    = 1
        call raise_exception_imgfile( n, ldim, 'apply_ctf_imgfile' )
        if( n /= o%get_noris() ) THROW_HARD('inconsistent nr entries; apply_ctf_imgfile')
        call stkio_r%open(fname2process, smpd, 'read')
        call stkio_w%open(fname,         smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> APPLYING CTF TO IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            tfun = ctf(smpd, o%get(i,'kv'), o%get(i,'cs'), o%get(i,'fraca'))
            if( o%isthere('dfy') )then ! astigmatic CTF
                call tfun%apply(img, o%get_dfx(i), mode, dfy=o%get_dfy(i), angast=o%get(i,'angast'), bfac=bfac)
            else ! non-astigmatic CTF
                call tfun%apply(img, o%get_dfx(i), mode, bfac=bfac)
            endif
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine apply_ctf_imgfile

    subroutine mask_imgfile( fname2mask, fname, mskrad, smpd, inner, width, which )
        character(len=*),           intent(in) :: fname2mask, fname
        real,                       intent(in) :: mskrad
        real,                       intent(in) :: smpd
        real,             optional, intent(in) :: inner, width
        character(len=*), optional, intent(in) :: which
        type(stack_io)                         :: stkio_r, stkio_w
        type(image)                            :: img
        integer                                :: n, i, ldim(3)
        character(len=STDLEN)                  :: msktype
        msktype = 'soft'
        if( present(which) )msktype=which
        call find_ldim_nptcls(fname2mask, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'mask_imgfile' )
        call stkio_r%open(fname2mask, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> MASKING IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%norm()
            call img%mask(mskrad, which, inner=inner, width=width)
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine mask_imgfile

    subroutine taper_edges_imgfile( fname2mask, fname, smpd)
        character(len=*),           intent(in) :: fname2mask, fname
        real,                       intent(in) :: smpd
        type(stack_io)                         :: stkio_r, stkio_w
        type(image)                            :: img
        integer                                :: n, i, ldim(3)
        call find_ldim_nptcls(fname2mask, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'taper_edges_imgfile' )
        call stkio_r%open(fname2mask, smpd, 'read')
        call stkio_w%open(fname,      smpd, 'write', box=ldim(1))
        call img%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> TAPERING EDGES OF IMAGES'
        do i=1,n
            call progress(i,n)
            call stkio_r%read(i, img)
            call img%taper_edges
            call stkio_w%write(i, img)
        end do
        call stkio_r%close
        call stkio_w%close
        call img%kill
    end subroutine taper_edges_imgfile

   subroutine shift_imgfile( fname2shift, fname, o, smpd, mul )
       character(len=*),  intent(in)    :: fname2shift, fname
       class(oris),       intent(inout) :: o
       real,              intent(in)    :: smpd
       real,    optional, intent(in)    :: mul
       type(stack_io) :: stkio_r, stkio_w
       type(image)    :: img
       integer        :: n, i, ldim(3)
       real           :: x, y, xhere, yhere
       call find_ldim_nptcls(fname2shift, ldim, n)
       ldim(3) = 1
       call raise_exception_imgfile( n, ldim, 'shift_imgfile' )
       if( n /= o%get_noris() ) THROW_HARD('inconsistent nr entries; shift_imgfile')
       call stkio_r%open(fname2shift, smpd, 'read')
       call stkio_w%open(fname,       smpd, 'write', box=ldim(1))
       call img%new(ldim,smpd)
       write(logfhandle,'(a)') '>>> SHIFTING IMAGES'
       do i=1,n
           call progress(i,n)
           call stkio_r%read(i, img)
           x = o%get(i, 'x')
           y = o%get(i, 'y')
           if( present(mul) )then
               xhere = x * mul
               yhere = y * mul
           else
               xhere = x
               yhere = y
           endif
           call img%shift([-xhere,-yhere,0.])
           call stkio_w%write(i, img)
       end do
       call stkio_r%close
       call stkio_w%close
       call img%kill
   end subroutine shift_imgfile

    subroutine random_selection_from_imgfile( spproj, fname, box, nran )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: fname
        integer,           intent(in)    :: nran, box
        type(stack_io)                   :: stkio_w
        type(ran_tabu)                   :: rt
        type(image)                      :: img, img_scaled
        character(len=:), allocatable    :: stkname
        logical,          allocatable    :: mask(:)
        real                             :: smpd
        integer                          :: nptcls, box_ori, ldim(3), i, ii, ldim_scaled(3), ind
        logical                          :: doscale
        ! dimensions
        nptcls      = spproj%get_nptcls()
        smpd        = spproj%os_stk%get(1,'smpd')
        box_ori     = nint(spproj%os_stk%get(1,'box'))
        ldim        = [box_ori,box_ori,1]
        ldim_scaled = [box,box,1]
        doscale     = box /= box_ori
        call raise_exception_imgfile( nptcls, ldim, 'random_selection_from_imgfile' )
        if( doscale )then
            call img_scaled%new(ldim_scaled,smpd) ! this sampling distance will  be overwritten
            call stkio_w%open(fname, smpd, 'write', box=ldim_scaled(1))
        else
            call stkio_w%open(fname, smpd, 'write', box=ldim(1))
        endif
        call img%new(ldim,smpd)
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
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> RANDOMLY SELECTING IMAGES'
        do i = 1,nran
            call progress(i, nran)
            ii = rt%irnd()
            call rt%insert(ii)
            call spproj%get_stkname_and_ind('ptcl2D', ii, stkname, ind)
            call img%read(stkname, ind)
            call img%norm()
            if( doscale )then
                call img%fft()
                if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2) )then
                    call img%clip(img_scaled)
                else
                    call img%pad(img_scaled)
                endif
                call img_scaled%ifft()
                call stkio_w%write(i, img_scaled)
            else
                call stkio_w%write(i, img)
            endif
        end do
        call stkio_w%close
        call rt%kill
        call img%kill
        call img_scaled%kill
    end subroutine random_selection_from_imgfile

    subroutine selection_from_tseries_imgfile( spproj, fname, box, nsel )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: fname
        integer,           intent(in)    :: nsel, box
        type(stack_io)                   :: stkio_r, stkio_w
        type(image)                      :: img, img_scaled
        character(len=:), allocatable    :: stkname
        real,             allocatable    :: states(:)
        integer,          allocatable    :: inds(:), inds_packed(:)
        real                             :: smpd
        integer                          :: nptcls, box_ori, ldim(3), i, n
        integer                          :: ldim_scaled(3), ind, cnt, stepsz
        logical                          :: doscale
        nptcls      = spproj%get_nptcls()
        smpd        = spproj%os_stk%get(1,'smpd')
        box_ori     = nint(spproj%os_stk%get(1,'box'))
        ldim        = [box_ori,box_ori,1]
        ldim_scaled = [box,box,1]
        doscale     = box /= box_ori
        call raise_exception_imgfile( nptcls, ldim, 'selection_from_tseries_imgfile' )
        if( doscale )then
            call img_scaled%new(ldim_scaled,smpd) ! this sampling distance will  be overwritten
            call stkio_w%open(fname, smpd, 'write', box=ldim_scaled(1))
        else
            call stkio_w%open(fname, smpd, 'write', box=ldim(1))
        endif
        call img%new(ldim,smpd)
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
        do i = 1,n,stepsz
            cnt = cnt + 1
            if( cnt > nsel ) exit
            call spproj%get_stkname_and_ind('ptcl2D', inds_packed(i), stkname, ind)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, ldim) )then
                call stkio_r%close
                call stkio_r%open(stkname, smpd, 'read')
            endif
            call stkio_r%read(ind, img)
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
                call stkio_w%write(cnt, img_scaled)
            else
                call stkio_w%write(cnt, img)
            endif
        end do
        deallocate(states, inds, inds_packed)
        call stkio_r%close
        call stkio_w%close
        call img%kill
        call img_scaled%kill
    end subroutine selection_from_tseries_imgfile

    subroutine random_cls_from_imgfile( spproj, fname, ncls )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: fname
        integer,           intent(in)    :: ncls
        type(stack_io)                   :: stkio_r, stkio_w
        type(image)                      :: img, img_avg
        character(len=:), allocatable    :: stkname
        real                             :: smpd
        integer                          :: nptcls, ldim(3), i, j, ii, ind, box
        integer, parameter               :: navg = 100
        ! dimensions
        nptcls = spproj%get_nptcls()
        smpd   = spproj%os_stk%get(1,'smpd')
        box    = nint(spproj%os_stk%get(1,'box'))
        ldim   = [box,box,1]
        call raise_exception_imgfile( nptcls, ldim, 'random_selection_from_imgfile' )
        call stkio_w%open(fname, smpd, 'write', box=box)
        call img%new(ldim,smpd)
        call img_avg%new(ldim,smpd)
        write(logfhandle,'(a)') '>>> RANDOMLY GENERATING CLUSTER CENTERS'
        do i = 1,ncls
            call progress(i, ncls)
            img_avg = 0.
            do j = 1,navg
                ii = irnd_uni(nptcls)
                call spproj%get_stkname_and_ind('ptcl2D', ii, stkname, ind)
                if( .not. stkio_r%stk_is_open() )then
                    call stkio_r%open(stkname, smpd, 'read')
                else if( .not. stkio_r%same_stk(stkname, ldim) )then
                    call stkio_r%close
                    call stkio_r%open(stkname, smpd, 'read')
                endif
                call stkio_r%read(ind, img)
                call img%norm()
                call img_avg%add(img)
            enddo
            call img_avg%norm()
            call stkio_w%write(i, img_avg)
        enddo
        call stkio_r%close
        call stkio_w%close
        call img%kill
        call img_avg%kill
    end subroutine random_cls_from_imgfile

end module simple_procimgstk

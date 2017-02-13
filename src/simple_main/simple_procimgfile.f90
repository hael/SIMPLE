!>  \brief  SIMPLE stack image processing routines for SPIDER/MRC files
module simple_procimgfile
use simple_defs
use simple_image,  only: image
use simple_jiffys, only: find_ldim_nptcls, progress
implicit none

private :: raise_exception_imgfile
public

contains
    
    !>  \brief  is for raising exception
    subroutine raise_exception_imgfile( n, ldim, routine )
        integer, intent(in) :: n, ldim(3)
        character(len=*)    :: routine
        if( n < 1 .or. any(ldim == 0))then
            write(*,*) routine
            write(*,*) 'The input stack is corrupt!'
            write(*,*) 'Number of images: ', n
            write(*,*) 'Logical dimensions: ', ldim
            stop
        endif
    end subroutine raise_exception_imgfile

    !> \brief  for checking if file is completely written to disk
    logical function file_written2disk( fname )
        use simple_image, only: image
        character(len=*), intent(in) :: fname
        integer     :: ldim(3), nframes, iframe
        type(image) :: testimg
        logical     :: read_failure
        ! get_number of frames and dim from stack
        ! this assumes that the header IS written to disk
        ! it should be because it is small
        call find_ldim_nptcls(fname, ldim, nframes, endconv=endconv)
        ldim(3) = 1
        call testimg%new(ldim, 1.0)
        do iframe=1,nframes
            call testimg%read(fname, iframe, read_failure=read_failure)
            if( read_failure ) exit
        end do
        call testimg%kill
        file_written2disk = .not. read_failure
    end function file_written2disk

    !>  \brief  is for copying
    subroutine copy_imgfile( fname2copy, fname, fromto , smpd_in)
        character(len=*), intent(in)  :: fname2copy, fname !< filenames
        integer, intent(in), optional :: fromto(2)         !< range
        real, intent(in), optional    :: smpd_in
        type(image)                   :: img
        integer                       :: n, i, cnt, ldim(3)
        !local variable
        real :: smpd
        call find_ldim_nptcls(fname2copy, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile(n, ldim, 'copy_imgfile')
        if ( .not. present(smpd_in) ) then
           smpd = 1.0 
        else
           smpd = smpd_in
        end if
        ! do the work
        call img%new(ldim,smpd)
        if( n >= 1 )then
            write(*,'(a)') '>>> COPYING IMAGES'
            if( present(fromto) )then
                cnt = 0
                do i=fromto(1),fromto(2)
                    cnt = cnt+1
                    call progress(cnt, fromto(2)-fromto(1)+1)
                    call img%read(fname2copy, i)
                    call img%write(fname, cnt)
                end do
            else                
                do i=1,n
                    call progress(i, n)
                    call img%read(fname2copy, i)
                    call img%write(fname, cnt)
                end do
            endif
        endif
        call img%kill
    end subroutine copy_imgfile
    
    !>  \brief  is for making a stack of normalized vectors for PCA analysis    
    subroutine make_pattern_stack( fnameStack, fnamePatterns, mskrad, D, recsz, avg, otab, hfun )
        use simple_filehandling, only: get_fileunit, fopen_err
        use simple_stat,         only: normalize, normalize_sigm
        use simple_oris,         only: oris
        use simple_math,         only: check4nans
        use simple_jiffys,       only: alloc_err
        character(len=*),            intent(in)    :: fnameStack, fnamePatterns !< filenames
        real,                        intent(in)    :: mskrad                    !< mask radius
        integer,                     intent(out)   :: D, recsz                  !< D=number of pixels and recsz=record size 
        real, allocatable, optional, intent(out)   :: avg(:)                    !< average
        class(oris),       optional, intent(inout) :: otab                      !< orientation table
        character(len=*),  optional, intent(in)    :: hfun                      !< hidden unit function descriptor
        type(image)        :: img
        real, allocatable  :: pcavec(:)
        real               :: x, y
        integer            :: n, fnum, ier, i, alloc_stat, ldim(3)
        logical            :: err
        logical, parameter :: debug=.false.        
        call find_ldim_nptcls(fnameStack, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'make_pattern_stack' )
        ! build and initialise objects
        call img%new(ldim,1.)
        D = img%get_npix(mskrad)
        allocate(pcavec(D), stat=alloc_stat)
        call alloc_err('make_pattern_stack; simple_procimgfile, 1', alloc_stat)
        pcavec = 0.
        inquire(iolength=recsz) pcavec
        deallocate(pcavec)
        if( present(avg) )then
            allocate(avg(D), stat=alloc_stat)
            call alloc_err('make_pattern_stack; simple_procimgfile, 2', alloc_stat)
            avg = 0.
        endif
        ! extract patterns and write to file
        fnum = get_fileunit()
        open(unit=fnum, status='replace', action='readwrite', file=fnamePatterns,&
        access='direct', form='unformatted', recl=recsz, iostat=ier)
        call fopen_err('make_pattern_stack; simple_procimgfile', ier)
        write(*,'(a)') '>>> MAKING PATTERN STACK'
        do i=1,n
            call progress(i,n)
            call img%read(fnameStack, i)
            if( present(otab) )then
                ! shift image
                call img%fwd_ft
                x = otab%get(i, 'x')
                y = otab%get(i, 'y')
                call img%shift(-x, -y)
                ! rotate image
                call img%bwd_ft
                call img%rtsq(-otab%e3get(i), 0., 0.)
            else
                if( img%is_ft() ) call img%bwd_ft
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
            if( debug )then
                call check4nans(pcavec)
                call img%serialize(pcavec, mskrad)
                call img%write('unserialized.spi', i)
            endif
            deallocate(pcavec)
        end do
        if( present(avg) )then
            avg = avg/real(n)
            allocate(pcavec(D), stat=alloc_stat)
            call alloc_err('make_pattern_stack; simple_procimgfile, 3', alloc_stat)
            do i=1,n
                read(fnum,rec=i) pcavec
                pcavec = pcavec-avg
                write(fnum,rec=i) pcavec
            end do
            deallocate(pcavec)
        endif
        close(unit=fnum)
        call img%kill
    end subroutine make_pattern_stack
    
    !>  \brief  is for padding
    subroutine pad_imgfile( fname2pad, fname, ldim_pad )
        character(len=*), intent(in) :: fname2pad, fname !< filenames
        integer, intent(in)          :: ldim_pad(3)      !< desired dimensions
        type(image)                  :: img, img_pad
        integer                      :: n, i, ldim(3)
        real                         :: ave, sdev, var, med
        call find_ldim_nptcls(fname2pad, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'pad_imgfile' )
        ! do the work
        if( ldim_pad(1) >= ldim(1) .and. ldim(2) >= ldim(2)&
        .and. ldim(3) >= ldim(3) )then
            call img%new(ldim,1.)
            call img_pad%new(ldim_pad,1.)
            write(*,'(a)') '>>> PADDING IMAGES'
            do i=1,n
                call progress(i,n)
                call img%read(fname2pad, i)
                if( img%is_ft() )then
                    call img%pad(img_pad) ! FT state preserved
                else
                    call img%stats('background', ave, sdev, var, med=med) ! get background statistics
                    call img%pad(img_pad, backgr=med) ! FT state preserved
                endif
                call img_pad%write(fname, i)
            end do
            call img%kill
            call img_pad%kill
        end if
    end subroutine pad_imgfile
    
    !>  \brief  is for resizing
    subroutine resize_imgfile( fname2resize, fname, smpd, ldim_new, fromptop )
        character(len=*),  intent(in) :: fname2resize, fname !< filenames
        real,              intent(in) :: smpd                !< original sampling distance
        integer,           intent(in) :: ldim_new(3)         !< desired dimensions
        integer, optional, intent(in) :: fromptop(2)         !< particle range
        type(image)                   :: img, img_resized
        integer                       :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2resize, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'resize_imgfile' )
        ! do the work
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
            call img%fwd_ft
            if( ldim_new(1) <= ldim(1) .and. ldim_new(2) <= ldim(2)&
            .and. ldim_new(3) <= ldim(3) )then
                call img%clip(img_resized)
            else
                call img%pad(img_resized)
            endif
            call img_resized%bwd_ft
            call img_resized%write(fname, cnt)
        end do
        call img%kill
        call img_resized%kill
    end subroutine resize_imgfile
    
    !>  \brief  is for clipping
    subroutine clip_imgfile( fname2clip, fname, ldim_clip )
        character(len=*), intent(in) :: fname2clip, fname !< filenames
        integer, intent(in)          :: ldim_clip(3)      !< desired dimensions
        type(image)                  :: img, img_clip
        integer                      :: n, i, ldim(3)
        call find_ldim_nptcls(fname2clip, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'clip_imgfile' )
        ! do the work
        if( ldim_clip(1) <= ldim(1) .and. ldim_clip(2) <= ldim(2)&
        .and. ldim_clip(3) <= ldim(3) )then
            call img%new(ldim,1.)
            call img_clip%new(ldim_clip,1.)
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
    
    !>  \brief  is for resizing and clipping
    subroutine resize_and_clip_imgfile( fname2resize, fname, smpd, ldim_new, ldim_clip, fromptop )
        character(len=*),  intent(in) :: fname2resize, fname !< filenames
        real,              intent(in) :: smpd                !< original sampling distance
        integer,           intent(in) :: ldim_new(3)         !< scaled dimensions
        integer,           intent(in) :: ldim_clip(3)        !< clipped (final) dimensions
        integer, optional, intent(in) :: fromptop(2)         !< particle range
        type(image)                   :: img, img_resized, img_clip
        integer                       :: n, i, ldim(3), prange(2), cnt, sz
        call find_ldim_nptcls(fname2resize, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'resize_imgfile' )
        ! do the work
        call img%new(ldim,1.)
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
            call img%fwd_ft
            call img%clip(img_resized)
            call img_resized%bwd_ft
            call img_clip%new(ldim_clip,img_resized%get_smpd())
            if( ldim_clip(1) <= ldim_new(1) .and. ldim_clip(2) <= ldim_new(2)&
            .and. ldim_clip(3) <= ldim_new(3) )then
                call img_resized%clip(img_clip)
            else
                call img_resized%pad(img_clip)
            endif
            call img_clip%write(fname, cnt)
        end do
        call img%kill
        call img_resized%kill
        call img_clip%kill
    end subroutine resize_and_clip_imgfile
    
    !>  \brief  is for normalization
    subroutine norm_imgfile( fname2norm, fname, hfun )
        character(len=*), intent(in) :: fname2norm, fname !< filenames
        character(len=*), intent(in), optional :: hfun    !< hidden unit function for initialization
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%norm(hfun)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine norm_imgfile
    
    !>  \brief  is for normalization
    subroutine norm_ext_imgfile( fname2norm, avg, sdev, fname )
        character(len=*), intent(in) :: fname2norm, fname !< filenames
        real, intent(in)             :: avg, sdev         !< average, standard deviation
        type(image)   :: img
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'norm_imgfile' )
        ! do the work
        call img%new(ldim,1.)
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
    subroutine noise_norm_imgfile( fname2norm, msk, fname )
        character(len=*), intent(in) :: fname2norm, fname !< filenames
        real, intent(in) :: msk
        type(image)      :: img
        integer          :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'noise_norm_imgfile' )
        ! do the work
        call img%new(ldim,1.)
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
    subroutine shellnorm_imgfile( fname2norm, fname )
        character(len=*), intent(in) :: fname2norm, fname !< filenames
        type(image)      :: img
        integer          :: i, n, ldim(3)
        call find_ldim_nptcls(fname2norm, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'noise_norm_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> SHELL NORMALIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2norm, i)
            call img%shellnorm
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine shellnorm_imgfile

    !>  \brief  is for calculating statistics
    subroutine stats_imgfile( fname, which, ave, sdev, var, med, msk )
        character(len=*), intent(in) :: fname, which        !< filename, selector
        real, intent(out)            :: ave, sdev, var, med !< ouputs
        real, intent(in), optional   :: msk                 !< optional mask parameter
        real        :: ave_loc, sdev_loc, var_loc, med_loc
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'stats_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> CALCULATING STACK STATISTICS'
        ave  = 0.
        sdev = 0.
        var  = 0.
        med  = 0.
        do i=1,n
            call progress(i,n)
            call img%read(fname,i)
            call img%stats(which, ave_loc, sdev_loc, var_loc, msk=msk, med=med_loc) 
            ave  = ave  + ave_loc
            sdev = sdev + sdev_loc
            var  = var  + var_loc
            med  = med  + med_loc
        end do
        ave  = ave/real(n)
        sdev = sdev/real(n)
        var  = var/real(n)
        med  = med/real(n)
        call img%kill
    end subroutine stats_imgfile
    
    !>  \brief  is for contrast inversion
    subroutine neg_imgfile( fname2neg, fname )
        character(len=*), intent(in) :: fname2neg, fname !< filenames
        type(image)   :: img
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2neg, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'neg_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> INVERTING IMAGE CONTRAST'
        do i=1,n
            call progress(i,n)
            call img%read(fname2neg, i)
            call img%neg
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine neg_imgfile
    
    !>  \brief  is for centering based on center of mass
    subroutine masscen_imgfile( fname2masscen, fname, smpd, lp, neg, msk, tres )
        character(len=*), intent(in) :: fname2masscen, fname
        real,             intent(in) :: smpd, lp
        character(len=*), intent(in) :: neg
        real, optional,   intent(in) :: msk, tres
        type(image) :: img
        integer     :: i, n, ldim(3)
        real        :: xyz(3)
        call find_ldim_nptcls(fname2masscen, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'masscen_imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        write(*,'(a)') '>>> CENTERING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2masscen, i)
            xyz = img%center(lp, neg, msk, tres)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine masscen_imgfile
    
    !>  \brief  is for curing
    subroutine cure_imgfile( fname2cure, fname )
        use simple_image,        only: image
        use simple_filehandling, only: get_fileunit
        character(len=*), intent(in) :: fname2cure, fname
        type(image) :: img
        integer     :: i, n, n_nans, filnum, ldim(3)
        real        :: maxv, minv, ave, sdev
        call find_ldim_nptcls(fname2cure, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'cure_imgfile' )
        ! do the work
        filnum = get_fileunit()
        open(unit=filnum, status='replace', file='cure_stats.txt')
        call img%new(ldim,1.)
        write(*,'(a)') '>>> CURING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2cure, i)
            call img%cure(maxv, minv, ave, sdev, n_nans)
            write(filnum,'(A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,f12.2,1X,A,1X,I9)') 'MAX:', maxv,&
            'MIN:', minv, 'AVE:', ave, 'SDEV:', sdev, 'NANS:', n_nans
            call img%write(fname, i)
        end do
        close(filnum)
        call img%kill
    end subroutine cure_imgfile
    
     !>  \brief  is for calculating the acf of a stack
    subroutine acf_imgfile( fname2acf, fname )
        character(len=*), intent(in) :: fname2acf, fname
        type(image)   :: img
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2acf, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'acf_imgfile' )
        ! do the work
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
    subroutine add_noise_imgfile( fname2process, fname, snr )
        character(len=*), intent(in) :: fname2process, fname
        real, intent(in) :: snr
        type(image)      :: img
        integer          :: i, n, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'add_noise_imgfile' )
        ! do the work
        call img%new(ldim,1.)
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
    subroutine frameavg_imgfile( fname2process, fname, navg )
        character(len=*), intent(in) :: fname2process, fname
        integer, intent(in) :: navg
        type(image)         :: img, avg
        integer :: i, n, cnt, cnt2, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'frameavg_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        call avg%new(ldim,1.)
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
    
    !>  \brief  is for swapping the sign of the imaginary part of the FT
    subroutine signswap_aimag_imgfile( fname2process, fname )
        character(len=*), intent(in) :: fname2process, fname
        type(image)   :: img
        integer       :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'signswap_aimag_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        do i=1,n
            call img%read(fname2process, i)
            call img%signswap_aimag
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine signswap_aimag_imgfile
    
    !>  \brief  is for swapping the sign of the real part of the FT
    subroutine signswap_real_imgfile( fname2process, fname )
        character(len=*), intent(in) :: fname2process, fname
        type(image)   :: img
        integer       :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'signswap_real_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        do i=1,n
            call img%read(fname2process, i)
            call img%signswap_real
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine signswap_real_imgfile
    
    !>  \brief  is for Fourier transformation visualization
    subroutine ft2img_imgfile( fname2process, fname, which )
        character(len=*), intent(in) :: fname2process, fname
        character(len=*), intent(in) :: which
        type(image)   :: img, img2
        integer       :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'ft2img_imgfile' )
        ! do the work
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
    subroutine shift_imgfile( fname2shift, fname, o, smpd, mul, round )
        use simple_oris,  only: oris
        character(len=*), intent(in)  :: fname2shift, fname
        class(oris), intent(inout)    :: o
        real, intent(in)              :: smpd
        real, intent(in), optional    :: mul
        logical, intent(in), optional :: round
        type(image)   :: img
        integer       :: n, i, ldim(3)
        real          :: x, y, xhere, yhere
        logical       :: rround
        call find_ldim_nptcls(fname2shift, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'shift_imgfile' )
        ! do the work
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
            call img%shift(-xhere, -yhere)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine shift_imgfile
    
    !>  \brief  is for band-pass filtering
    subroutine bp_imgfile( fname2filter, fname, smpd, hp, lp )
        character(len=*), intent(in) :: fname2filter, fname
        real, intent(in) :: smpd, hp, lp
        type(image)      :: img
        integer          :: n, i, ldim(3)
        real, parameter  :: width=10.
        call find_ldim_nptcls(fname2filter, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'bp_imgfile' )
        ! do the work
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
    
    !>  \brief  is for phase randomization
    subroutine phase_rand_imgfile( fname2process, fname, smpd, lp )
        character(len=*), intent(in) :: fname2process, fname
        real, intent(in) :: smpd, lp
        type(image)      :: img
        integer :: n, i, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'phase_rand_imgfile' )
        ! do the work
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
    subroutine apply_ctf_imgfile( fname2process, fname, o, smpd, mode, bfac )
        use simple_math,  only: deg2rad
        use simple_oris,  only: oris
        use simple_ctf,   only: ctf
        character(len=*), intent(in)    :: fname2process, fname !< filenames
        class(oris),      intent(inout) :: o                    !< orientations object
        real,             intent(in)    :: smpd                 !< sampling distance, a/pix
        character(len=*), intent(in)    :: mode                 !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: bfac                 !< bfactor
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
    
    !>  \brief  is for visualization of the CTF power spectrum
    subroutine ctf2imgfile( fname2process, fname, o, smpd, tfun, bfac )
        use simple_math,  only: deg2rad
        use simple_oris,  only: oris
        use simple_ctf,   only: ctf
        character(len=*), intent(in) :: fname2process, fname !< filenames
        class(oris), intent(inout)   :: o                    !< orientations object
        real, intent(in)             :: smpd                 !< sampling distance, a/pix
        class(ctf), intent(inout)    :: tfun                 !< transfer function 
        real, intent(in), optional   :: bfac                 !< bfactor
        type(image)   :: img, img2
        integer       :: i, n, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'ctf2imgfile' )
        ! do the work
        call img%new(ldim,smpd)
        call img2%new(ldim,smpd)
        write(*,'(a)') '>>> MAKING CTF POWSPEC IMAGES'
        do i=1,n
            call progress(i,n)
            if( o%isthere('dfy') )then ! astigmatic CTF   
                call tfun%ctf2img(img, o%get(i,'dfx'), 'square', dfy=o%get(i,'dfy'), angast=o%get(i,'angast'), bfac=bfac)
            else
                call tfun%ctf2img(img, o%get(i,'dfx'), 'square', bfac=bfac)
            endif
            call img%ft2img('real', img2)
            call img2%write(fname, i)
        end do
        call img%kill
        call img2%kill
    end subroutine ctf2imgfile
    
    !>  \brief  is for origin shifting and rotating a stack
    !!          according to info in o 
    subroutine shrot_imgfile( fname2process, fname, o, smpd, mul )
        use simple_oris,  only: oris
        character(len=*), intent(in) :: fname2process, fname !< filenames
        class(oris), intent(inout)   :: o                    !< orientations
        real, intent(in)             :: smpd                 !< sampling distance
        real, intent(in), optional   :: mul                  !< optional multiplicaton factor
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
            call img%fwd_ft
            x = o%get(i, 'x')
            y = o%get(i, 'y')
            if( present(mul) )then
                call img%shift(-x*mul, -y*mul)
            else
                call img%shift(-x, -y)
            endif
            call img%bwd_ft
            call img%rtsq(-o%e3get(i), 0., 0., img_rot)
            call img_rot%write(fname, i)
        end do
        call img%kill
        call img_rot%kill
    end subroutine shrot_imgfile
    
    !>  \brief  is for applying a soft circular mask to all images in stack
    subroutine mask_imgfile( fname2mask, fname, mskrad, inner, width, which )
        character(len=*), intent(in)           :: fname2mask, fname !< filenames
        real,             intent(in)           :: mskrad            !< mask radius
        real,             intent(in), optional :: inner, width      !< inner & outer mask radii (optional)
        character(len=*), intent(in), optional :: which
        type(image)           :: img
        integer               :: n, i, ldim(3)
        character(len=STDLEN) :: msktype
        msktype = 'soft'
        if( present(which) )msktype=which
        call find_ldim_nptcls(fname2mask, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'mask_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> MASKING IMAGES' 
        do i=1,n
            call progress(i,n)
            call img%read(fname2mask, i)
            call img%mask(mskrad, which, inner, width) ! FT state preserved 
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine mask_imgfile
    
    !>  \brief  is for binarizing all images in the stack
    subroutine bin_imgfile( fname2process, fname, tres )
        character(len=*), intent(in) :: fname2process, fname !< filenames
        real, intent(in), optional   :: tres                 !< optional treshold
        type(image) :: img
        integer     :: n, i, ldim(3)
        logical     :: didft
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'bin_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> BINARIZING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2process, i)
            didft = .false.
            if( img%is_ft() )then
                call img%bwd_ft
                didft = .true.
            endif
            if( present(tres) )then
                call img%bin(tres)
            else
                call img%bin
            endif
            if( didft ) call img%fwd_ft
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine bin_imgfile
    
    !>  \brief  is for making an average of the entire stack
    subroutine make_avg_imgfile( fname, avgname )
        character(len=*), intent(in) :: fname, avgname !< filenames
        type(image)                  :: avg, img
        integer                      :: i, n, ldim(3)
        call find_ldim_nptcls(fname, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'make_avg_imgfile' )
        ! do the work
        call avg%new(ldim, 1.)
        call img%new(ldim, 1.)
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
    
    !>  \brief  is for making class averages according to the 
    !!          classification info in o
    subroutine make_cavgs_imgfile( fname2process, fname, o, minp, list )
        use simple_oris,  only: oris
        use simple_sll,   only: sll
        character(len=*), intent(in)     :: fname2process, fname !< filenames
        class(oris), intent(inout)       :: o                    !< orientations
        integer, intent(in), optional    :: minp                 !< minimum number of particles          
        type(sll), intent(out), optional :: list                 !< optional singly linked list
        type(image)                      :: avg
        logical                          :: empty
        integer :: i, ncls, ncls_here, pop, ldim(3), n
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile(ncls, ldim, 'make_cavgs_imgfile' )
        ! extract the number of classes from the oris structure
        ncls = o%get_ncls()
        ! build
        if( present(list) ) call list%new
        ncls_here = 0
        write(*,'(a)') '>>> MAKING CLASS AVERAGES'
        do i=1,ncls
            call progress(i,ncls)
            pop = o%get_clspop(i)
            if( present(minp) )then
                if( pop >= minp )then
                    if( present(list) ) call list%add(iarr=[i])
                    ncls_here = ncls_here+1
                    call make_cavg_imgfile(fname2process, fname, o, i, avg, empty)
                    call avg%write(fname, ncls_here)
                endif   
            else
                if( pop > 0 )then
                    if( present(list) ) call list%add(iarr=[i])
                endif
                call make_cavg_imgfile(fname2process, fname, o, i, avg, empty)
                call avg%write(fname, i)
            endif     
        end do
    end subroutine make_cavgs_imgfile
    
    !>  \brief  is for making one class average according to the 
    !!          classification info in o and the inputted class
    subroutine make_cavg_imgfile( fname2process, fname, o, class, avg, empty, smpd )  
        use simple_oris,  only: oris
        character(len=*), intent(in) :: fname2process, fname !< filenames
        class(oris), intent(inout)   :: o                    !< orientations
        integer, intent(in)          :: class                !< class index
        type(image), intent(inout)   :: avg                  !< output average
        logical, intent(out)         :: empty                !< empty indicator
        real, intent(in), optional   :: smpd                 !< optional sampling distance
        type(image)                  :: img
        real                         :: ssmpd
        integer :: j, ncls, pop, n, mycls, ldim(3)
        call find_ldim_nptcls(fname2process, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'make_cavg_imgfile' )
        ! do the work
        ssmpd = 1.
        if( present(smpd) ) ssmpd = smpd
        ncls = o%get_ncls()
        pop = o%get_clspop(class) 
        empty = .false.
        call avg%new(ldim, ssmpd)
        if( pop < 2 )then
            empty = .true.
        else
            call img%new(ldim, ssmpd)
            avg = 0.         
            do j=1,n
                mycls = nint(o%get(j, 'class'))
                if( mycls == class )then
                    call img%read(fname2process, j)
                    call avg%add(img)
                endif
            end do
            call avg%div(real(pop))
            call img%kill
        endif
    end subroutine make_cavg_imgfile

    subroutine mirror_imgfile( fname2mirr, fname, mirr )
        character(len=*), intent(in) :: fname, fname2mirr, mirr
        type(image) :: img
        integer     :: i, n, ldim(3)
        call find_ldim_nptcls(fname2mirr, ldim, n)
        ldim(3) = 1
        call raise_exception_imgfile( n, ldim, 'mirror_imgfile' )
        ! do the work
        call img%new(ldim,1.)
        write(*,'(a)') '>>> MIRRORING IMAGES'
        do i=1,n
            call progress(i,n)
            call img%read(fname2mirr, i)
            call img%mirror(mirr)
            call img%write(fname, i)
        end do
        call img%kill
    end subroutine mirror_imgfile

end module simple_procimgfile

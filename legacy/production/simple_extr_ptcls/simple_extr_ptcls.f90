!==Program simple_extr_ptcls
!
! <extr_ptcls/begin> is a program that extracts particle images from DDD movies or integrated movies. Boxfiles are assumed to be in EMAN format
! but we provide a conversion script (\texttt{relion2emanbox.pl}) for *.star files containing particle coordinates obtained with Relion.
! The program creates one stack per movie frame as well as a stack of corrected framesums. In addition to single-particle image stacks, the program 
! produces a parameter file \texttt{extr\_ptcls\_params.txt} that can be used in conjunction with other SIMPLE programs. We obtain CTF parameters 
! with CTFFIND4 using the script (\texttt{exec_ctffind.pl}) but if you have already obtained CTF parameters from CTFFIND4, please see section 
! "CTF parameters convention" in the manual.<extr_ptcls/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2011-08-16.
!
program simple_extr_ptcls
use simple_defs       ! singleton
use simple_jiffys     ! singleton
use simple_cmdline,   only: cmdline
use simple_build,     only: build
use simple_params,    only: params
use simple_nrtxtfile, only: nrtxtfile
use simple_imgfile,   only: imgfile
use simple_image,     only: image
use simple_math,      only: euclid, hpsort, median
use simple_oris,      only: oris
use simple_stat,      only: moment
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
integer                            :: nmovies, nboxfiles, funit_movies, funit_box, nframes, frame, pind
integer                            :: i, j, k, alloc_stat, ldim(3), box_current, movie, ndatlines, nptcls
integer                            :: npix, npix_backgr, npix_tot, cnt, niter, nnn, fromto(2), orig_box
integer                            :: movie_ind, numlen, ntot, lfoo(3), ifoo
type(nrtxtfile)                    :: boxfile
character(len=STDLEN)              :: mode, framestack, sumstack
character(len=STDLEN), allocatable :: movienames(:), boxfilenames(:)
real, allocatable                  :: boxdata(:,:), nndists(:), noise_pixels(:), pixels(:)
integer, allocatable               :: nninds(:), pinds(:)
real                               :: x, y, dfx, dfy, angast, ctfres, med, ave,sdev, var
type(image)                        :: img_frame, mskimg
type(oris)                         :: outoris, local_oris
logical                            :: err
logical, parameter                 :: debug = .false.
if( command_argument_count() < 2 )then
    write(*,'(a)',advance='no') 'SIMPLE_EXTR_PTCLS filetab=<movies.txt> boxtab=<boxfiles.txt> smpd=<sampling'
    write(*,'(a)',advance='no') ' distance(in A)> msk=<rough estimate particle radius(pixels)> [oritab=<align'
    write(*,'(a)',advance='no') 'ment/CTF doc>] [neg=<yes|no{yes}>] [mapnn=<map nearest neighs(yes|no){no}>]'
    write(*,'(a)') ' [noise_norm=<yes|no{no}>] [outside=<yes|no{yes}>]'
    stop
endif
call cline%parse
call cline%checkvar('filetab', 1)
call cline%checkvar('boxtab',  2)
call cline%checkvar('smpd',    3)
call cline%checkvar('msk',     4)
if( .not. cline%defined('neg') )then
    call cline%set('neg', 'yes')
endif
call cline%check
p = params(cline, checkdistr=.false.) ! constants & derived constants produced
! check file inout existence
if( .not. file_exists(p%filetab) ) stop 'inputted filetab does not exist in cwd'
if( .not. file_exists(p%boxtab) )  stop 'inputted boxtab does not exist in cwd'
nmovies = nlines(p%filetab)
if( debug ) write(*,*) 'nmovies: ', nmovies
nboxfiles = nlines(p%boxtab)
if( debug ) write(*,*) 'nboxfiles: ', nboxfiles
if( nmovies /= nboxfiles ) stop 'number of entries in inputted files do not match!'
funit_movies = get_fileunit()
open(unit=funit_movies, status='old', file=p%filetab)
funit_box = get_fileunit()
open(unit=funit_box, status='old', file=p%boxtab)
allocate( movienames(nmovies), boxfilenames(nmovies), stat=alloc_stat )
call alloc_err('In: simple_extr_ptcls; boxdata etc., 1', alloc_stat)

! remove output file
call del_txtfile('extr_ptcls_params.txt')

! read the filenames
do movie=1,nmovies
    read(funit_movies,'(a256)') movienames(movie)
    read(funit_box,'(a256)') boxfilenames(movie)
end do
close(funit_movies)
close(funit_box)
if( debug ) write(*,*) 'read the filenames'

! fill up local oris
local_oris = oris(nmovies)
if( cline%defined('oritab') ) call local_oris%read(p%oritab)

! determine loop range
fromto(1) = 1
fromto(2) = nmovies
ntot = fromto(2)-fromto(1)+1
if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)

! count the number of particles & find ldim
nptcls = 0
do movie=fromto(1),fromto(2)
    if( file_exists(boxfilenames(movie)) )then 
        call boxfile%new(boxfilenames(movie), 1)
        ndatlines = boxfile%get_ndatalines()
        nptcls = nptcls+ndatlines
        call boxfile%kill
    else
        write(*,*) 'WARNING! The inputted boxfile (below) does not exist'
        write(*,*) trim(boxfilenames(movie))
    endif
end do
call find_ldim_nptcls(movienames(1), ldim, ifoo)
if( debug ) write(*,*) 'number of particles: ', nptcls

! create frame
call img_frame%new([ldim(1),ldim(2),1], p%smpd)

! initialize
call outoris%new(nptcls)
pind = 0
if( debug ) write(*,*) 'made outoris'

! loop over exposures (movies)
niter = 0
do movie=fromto(1),fromto(2)
    
    ! show progress
    if( niter > 1 )then
        call progress(niter,ntot)
    endif
    
    ! get movie index (parsing the number string from the filename)
    call fname2ind(trim(adjustl(movienames(movie))), movie_ind)

    ! process boxfile
    ndatlines = 0
    if( file_exists(boxfilenames(movie)) )then
        call boxfile%new(boxfilenames(movie), 1)
        ndatlines = boxfile%get_ndatalines()        
    endif
    if( ndatlines == 0 ) cycle
    
    ! update iteration counter (has to be after the cycle statements or suffer bug!!!)
    niter = niter+1
    
    ! read box data
    allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()),&
    nninds(ndatlines), nndists(ndatlines), pinds(ndatlines), stat=alloc_stat )
    call alloc_err('In: simple_extr_ptcls; boxdata etc., 2', alloc_stat)
    do j=1,ndatlines
        call boxfile%readNextDataLine(boxdata(j,:))
        orig_box = nint(boxdata(j,3))
        if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
            stop 'Only square windows are currently allowed!'
        endif
    end do
    
    ! create particle index list and set movie index
    if( .not. cline%defined('box') ) p%box = nint(boxdata(1,3)) ! use the box size from the box file
    do j=1,ndatlines
        if( box_inside(ldim, nint(boxdata(j,1:2)), p%box) )then
            pind     = pind+1
            pinds(j) = pind
            call outoris%set(pinds(j), 'movie', real(movie_ind))
        else
            pinds(j) = 0
        endif
    end do
    
    ! create nearest neighbor structure
    if( p%mapnn .eq. 'yes' )then
        if( ndatlines > 1 )then
            do i=1,ndatlines
                do j=1,ndatlines
                    if( i == j )then
                        nndists(j) = huge(x)
                    else
                        nndists(j) = euclid(boxdata(i,1:2), boxdata(j,1:2))
                    endif
                end do
                nninds = pinds
                call hpsort(ndatlines, nndists, nninds)
                nnn = min(20,ndatlines-1)
                do j=1,nnn
                    if( pinds(i) > 0 .and. nndists(j) > 0 )then
                        call outoris%set(pinds(i), 'nnn',             real(nnn))        ! number of nearest neighbors
                        call outoris%set(pinds(i), 'nn'//int2str(j),  real(nninds(j)))  ! nearest neighbor particle index
                        call outoris%set(pinds(i), 'dnn'//int2str(j), real(nndists(j))) ! nearest neighbor distance
                    endif
                end do
            end do
        endif
    endif
    
    ! check box parsing
    if( .not. cline%defined('box') )then
        if( niter == 1 )then
            p%box = nint(boxdata(1,3))   ! use the box size from the box file
        else
            box_current = nint(boxdata(1,3))
            if( box_current /= p%box )then
                write(*,*) 'box_current: ', box_current, 'box in params: ', p%box
                stop 'inconsistent box sizes in box files'
            endif
        endif
        if( p%box == 0 )then
            write(*,*) 'ERROR, box cannot be zero!'
            stop
        endif
    endif
    if( debug ) write(*,*) 'did check box parsing'
    
    ! get number of frames from stack
    call find_ldim_nptcls(movienames(movie), lfoo, nframes )
    numlen  = len(int2str(nframes))
    if( debug ) write(*,*) 'number of frames: ', nframes
    
    ! build general objects
    if( niter == 1 )then
        call b%build_general_tbox(p,cline,do3d=.false.)
        call mskimg%disc([p%box,p%box,1], p%smpd, p%msk, npix)
        call mskimg%bin_inv
        npix_backgr = p%box*p%box-npix
    endif
    npix_tot = npix_backgr*ndatlines
    if( npix_tot == 0 )then
        write(*,*) "ERROR, total nr of pixels cannot be zero"
        write(*,*) "npix_backgr: ", npix_backgr
        write(*,*) "ndatlines: ", ndatlines
        stop
    endif
    allocate( noise_pixels(npix_tot), stat=alloc_stat )
    call alloc_err('In: simple_extr_ptcls; boxdata, 3', alloc_stat)

    ! extract ctf info
    if( b%a%isthere('dfx') )then
        dfx = b%a%get(movie,'dfx')
        ctfres = b%a%get(movie,'ctfres')
        angast= 0.
        if( b%a%isthere('dfy') )then ! astigmatic CTF
            if( .not. b%a%isthere('angast') ) stop 'need angle of astigmatism for CTF correction'
            dfy = b%a%get(movie,'dfy')
            angast = b%a%get(movie,'angast')
        endif
        ! set CTF info in outoris
        do i=1,ndatlines
            if( pinds(i) > 0 )then
                call outoris%set(pinds(i), 'dfx', dfx)
                call outoris%set(pinds(i), 'ctfres', ctfres)
                if( b%a%isthere('dfy') )then
                    call outoris%set(pinds(i), 'angast', angast)
                    call outoris%set(pinds(i), 'dfy', dfy)
                endif
            endif
        end do
        if( debug ) write(*,*) 'did set CTF parameters dfx/dfy/angast/ctfres: ', dfx, dfy, angast, ctfres
    endif
    
    ! loop over frames
    do frame=1,nframes
        
        ! read frame
        call img_frame%read(movienames(movie),frame)
        
        if( nframes > 1 )then
            ! shift frame according to global shift (drift correction)
            if( b%a%isthere(movie, 'x'//int2str(frame)) .and. b%a%isthere(movie, 'y'//int2str(frame))  )then
                call img_frame%fwd_ft
                x = b%a%get(movie, 'x'//int2str(frame))
                y = b%a%get(movie, 'y'//int2str(frame))
                call img_frame%shift(-x,-y)
                call img_frame%bwd_ft
            else
                write(*,*) 'no shift parameters available for alignment of frames'
                stop 'use simple_unblur_movies if you want to integrate frames'
            endif
        endif
        
        ! extract the particle images & normalize
        if( nframes > 1 )then
            framestack = 'framestack'//int2str_pad(frame,numlen)//p%ext
        else
            framestack = 'sumstack'//p%ext
        endif
        cnt = 0
        do j=1,ndatlines ! loop over boxes
            if( pinds(j) > 0 )then
                ! modify coordinates if change in box (shift by half)
                if( debug ) print *, 'original coordinate: ', boxdata(j,1:2) 
                if( orig_box /= p%box ) boxdata(j,1:2) = boxdata(j,1:2)-real(p%box-orig_box)/2.
                if( debug ) print *, 'shifted coordinate: ', boxdata(j,1:2) 
                ! extract the window    
                call img_frame%window(nint(boxdata(j,1:2)), p%box, b%img)
                if( p%neg .eq. 'yes' ) call b%img%neg
                pixels = b%img%extr_pixels(mskimg)
                do k=1,size(pixels)
                    cnt = cnt+1
                    noise_pixels(cnt) = pixels(k)
                end do
                deallocate(pixels)
                call b%img%write(trim(adjustl(framestack)), pinds(j))
            endif
        end do
        med = median(noise_pixels(:cnt))
        call moment(noise_pixels(:cnt), ave, sdev, var, err)
        ! normalize the framestack according to the median and
        ! standard devation of the noise over the entire frame
        do j=1,ndatlines ! loop over boxes
            if( pinds(j) > 0 )then
                call b%img%read(trim(adjustl(framestack)), pinds(j))
                if( p%noise_norm .eq. 'yes' )then
                    call b%img%norm_ext(med, sdev)
                else
                    call b%img%norm
                endif
                call b%img%write(trim(adjustl(framestack)), pinds(j))
            endif
        end do
    end do
    
    if( nframes > 1 )then
        ! create the sum and normalize it
        sumstack = 'sumstack'//p%ext
        cnt = 0
        do j=1,ndatlines ! loop over boxes
            if( pinds(j) > 0 )then
                b%img_copy = 0.
                do frame=1,nframes ! loop over frames
                    framestack = 'framestack'//int2str_pad(frame,numlen)//p%ext
                    call b%img%read(trim(adjustl(framestack)), pinds(j))
                    call b%img_copy%add(b%img)
                end do
                call b%img_copy%write(trim(adjustl(sumstack)), pinds(j))
                pixels = b%img_copy%extr_pixels(mskimg)
                do k=1,size(pixels)
                    cnt = cnt+1
                    noise_pixels(cnt) = pixels(k)
                end do
                deallocate(pixels)
            endif
        end do
        med = median(noise_pixels(:cnt))
        call moment(noise_pixels(:cnt), ave, sdev, var, err)
        ! normalize the sumstack according to the median and
        ! standard devation of the noise over the entire framesum
        do j=1,ndatlines ! loop over boxes
            if( pinds(j) > 0 )then
                call b%img%read(trim(adjustl(sumstack)), pinds(j))
                if( p%noise_norm .eq. 'yes' )then
                    call b%img%norm_ext(med, sdev)
                else
                    call b%img%norm
                endif
                call b%img%write(trim(adjustl(sumstack)), pinds(j))
            endif
        end do
    endif
    
    ! write output
    do j=1,ndatlines ! loop over boxes
        if( pinds(j) > 0 )then
            call outoris%write(pinds(j), 'extr_ptcls_params.txt')
        endif
    end do
    
    ! destruct
    call boxfile%kill
    deallocate(boxdata, nninds, nndists, noise_pixels, pinds)
end do
call simple_end('**** SIMPLE_EXTR_PTCLS NORMAL STOP ****')

contains
    
    function box_inside( ldim, coord, box ) result( inside )
        integer, intent(in) :: ldim(3), coord(2), box
        integer             :: fromc(2), toc(2)
        logical             :: inside
        if( p%outside .eq. 'yes' )then
            inside = .true.
            return
        endif
        fromc = coord+1     ! compensate for the c-range that starts at 0
        toc = fromc+(box-1) ! the lower left corner is 1,1
        inside = .true.
        if( toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
    end function
    
end program simple_extr_ptcls

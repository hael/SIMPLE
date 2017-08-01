!==Program simple_volops
!
! <volops/begin> provides standard single-particle image processing routines that are applied to MRC or SPIDER volumes.
! <comment/begin> If you input two volumes and the sampling distance, the FSC is calculated between the volumes. The 
! FSC plot is written to STDOUT together with resolution estimates at $FSC=0.5$ and $FSC=0.143$. The volumes subjected 
! to FSC calculation should be masked with a soft-edged (not hard-edged) mask and they should not have been subjected 
! to any "auto" or threshold masking. If \texttt{phrand} and \texttt{lp} are given, the Fourier phases of the input 
! volume \texttt{vol1} are randomized. \texttt{msk} is used for spherical masking with a soft (cosine edge) fall-off. 
! \texttt{lp} and \texttt{hp} are the low-pass and high-pass limits used for filtering. To add noise to a volume, give
! the desired signal-to-noise ratio via \texttt{snr}. Give \texttt{center=yes} and \texttt{lp} to center the input volume 
! according to center of mass. The 3D origin shift vector is found by low-pass filtering the volume to \texttt{lp}, 
! binarizing the density, identifying the center of mass, and calculating the vector needed to place the center of mass 
! in the center of the box. \texttt{soften=yes} applies first a real-space low-pass filter then softened by a cosine edge of
! pixel width \texttt{edge}. \texttt{mskfile} is used for masking a volume using an externally generated mask. 
! \texttt{countvox=yes} counts the number of foreground voxels (the binarization method is k-means). \texttt{newbox} and 
! \texttt{scale} are used for resizing the volume. \texttt{msktype} controls the mask type (\texttt{hard} or \texttt{soft}). 
! \texttt{inner} controls the radius of the inner mask with fall-off \texttt{width}. \texttt{cube} is used to generate a binary
! cube (4 testing purposes). \texttt{e1,e2,e3} is the Euler triplet used to rotate the input volume using Kaiser-Bessel 
! interpolation in Fourier space. \texttt{corner} is used for filling in the corners of the box with binary cubes (4 testing
! purposes). \texttt{neg} inverts the contrast of the input volume by multiplication with $-1$ in Fourier space. 
! \texttt{voltab} and \texttt{voltab2} are used to give text files with the names of volume files that are correlated and the
! nearest neighbour structure of the comparison is written to STDOUT.<comment/end> <volops/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_volops
use simple_defs         ! singleton
use simple_jiffys,      ! singleton
use simple_cmdline,     only: cmdline
use simple_build,       only: build
use simple_params,      only: params
use simple_image,       only: image
use simple_ori,         only: ori
use simple_math,        only: get_resolution, hpsort, round2odd, nvoxfind
use simple_sll,         only: sll
use simple_timing
implicit none
type(params)                       :: p
type(build)                        :: b
type(cmdline)                      :: cline
type(ori)                          :: o
type(image)                        :: vol2, mskimg
real, allocatable                  :: res(:), corrs(:), spec(:)
integer, allocatable               :: order(:)
character(len=STDLEN), allocatable :: signs(:), files(:)
real                               :: fsc05, fsc0143, ave, sdev, var, med, xyz(3)
integer                            :: j, npix, alloc_stat, i, ntot, nvgrp1, nvgrp2, k
logical                            :: here
!Start of the execution commands
call timestamp()
call start_Alltimers_cpu()
!parameter
if( command_argument_count() < 1 )then
   write(*,'(a)', advance='no') 'SIMPLE_VOLOPS [vol1=<invol.ext>] [vol2=<invol2.ext>]'
   write(*,'(a)', advance='no') ' [smpd=<sampling distance(in A)>] [outvol=<outvol.ext>]'
   write(*,'(a)', advance='no') ' [nthr=<nr of openMP threads{1}>] [phrand=<yes|no{no}>]'
   write(*,'(a)', advance='no') ' [msk=<mask radius(in pixels)>] [lp=<low-pass limit{20}>]'
   write(*,'(a)', advance='no') ' [hp=<high-pass limit{100}>] [snr=<signal-to-noise ratio>]'
   write(*,'(a)', advance='no') ' [center=<yes|no{no}>] [soften=<yes|no{no}>] [guinier=<yes|no{no}>]'
   write(*,'(a)', advance='no') ' [bfac=<bfactor(in A**2){200.}>] [edge=<edge size for softening'
   write(*,'(a)', advance='no') ' molecular envelope(in pixels){12}>] [mskfile=<mask.ext>]'
   write(*,'(a)', advance='no') ' [countvox=<yes|no{no}>] [newbox=<scaled box>] [scale=<scale factor{1}>]'
   write(*,'(a)', advance='no') ' [msktype=<hard|soft{soft}>] [inner=<inner mask radius(in pixels)>]'
   write(*,'(a)', advance='no') ' [width=<pixels falloff inner mask{10}>] [cube=<side (in pixels){0}>]'
   write(*,'(a)', advance='no') ' [e1=<1st Euler{0}>] [e2=<2nd Euler{0}>] [e3=<3d Euler{0}>]'
   write(*,'(a)', advance='no') ' [corner=<corner size{0}>] [neg=<yes|no{no}>]'
   write(*,'(a)', advance='no') ' [voltab=<file table>] [voltab2=<file table>] [bin=<yes|no{no}>]'
   write(*,'(a)', advance='no') ' [nvox=<nr of voxels{0}>] [xsh=<x shift(pixels){0}>]'
   write(*,'(a)') ' [ysh=<y shift(pixels){0}>] [zsh=<z shift(pixels){0}>] [norm=<yes|no{no}>]'
   stop  
endif
call cline%parse
call cline%set('prg', 'volops')
p = params(cline,checkpara=.false.)         ! constants & derived constants produced, mode=2
call b%build_general_tbox(p, cline)         ! general objects built
call b%vol%new([p%box,p%box,p%box], p%smpd) ! reallocate vol (boxmatch issue)
if( cline%defined('voltab') )then
    if( cline%defined('voltab2') )then
        if( .not. cline%defined('box')  ) stop 'need box size for corr calc!'
        if( .not. cline%defined('smpd') ) stop 'need sampling distance (smpd) for corr calc!'
        if( .not. cline%defined('lp')   ) stop 'need low-pass limit for corr calc!'
        if( .not. cline%defined('hp')   ) stop 'need high-pass limit for corr calc! (set it to half the particle diam)'
        write(*,*) '>>> ASSUMING VOLUMES TO BE COMPARED HAVE BEEN MASKED WITH A SOFT, NOT HARD OR THRESHOLD, MASK'
        call vol2%new([p%box,p%box,p%box],p%smpd)
        inquire(FILE=p%voltab, EXIST=here)
        if( .not. here )then
            write(*,*) 'File does not exist! ', trim(adjustl(p%voltab))
            stop 
        endif
        inquire(FILE=p%voltab2, EXIST=here)
        if( .not. here )then
            write(*,*) 'File does not exist! ', trim(adjustl(p%voltab2))
            stop 
        endif
        nvgrp1 = nlines(p%voltab)
        nvgrp2 = nlines(p%voltab2)
        ntot   = nvgrp1+nvgrp2
        signs  = make_signs_two_groups(nvgrp1, nvgrp2)
        files  = merge_txtfiles(p%voltab, p%voltab2)
        allocate( corrs(ntot), order(ntot), stat=alloc_stat )
        call alloc_err('In: simple_volops', alloc_stat)
        write(*,'(a)') '>>> PRINTING NEAREST NEIGHBOR STRUCTURE'
        do i=1,ntot
            order=(/(j,j=1,ntot)/)
            call b%vol%read(files(i))
            do j=1,ntot
                if( i == j )then
                    corrs(j) = -1.
                else
                    call vol2%read(files(j))
                    corrs(j) = b%vol%corr(vol2,lp_dyn=p%lp,hp_dyn=p%hp)
                endif
            end do
            call hpsort(ntot, corrs, order)
            write(*,'(1X,A,A,A,3X)', advance="no") '((',trim(adjustl(signs(i))),'))'
            do j=ntot,2,-1
                write(*,'(A,1X,F7.4,3X)', advance="no") trim(adjustl(signs(order(j)))), corrs(j) 
            end do
            write(*,*) ''
        end do
    else
        stop 'need second table (voltab2) for volume comparison!'
    endif
    goto 999
endif
inquire(FILE=p%vols(1), EXIST=here)
if( here ) call b%vol%read(p%vols(1))
if( here )then
    if( p%norm  .eq. 'yes' )then
        call b%vol%shellnorm
        spec = b%vol%spectrum('power')
        do k=1,size(spec)
            print *, k, spec(k)
        end do
    else if( p%center .eq. 'yes' )then
        if( .not. cline%defined('smpd') ) stop 'smpd needed 4 centering'
        if( .not. cline%defined('lp') )   stop 'low-pass limit needed 4 centering'
        xyz = b%vol%center(p%lp,p%msk)
        write(*,*) 'Identified center (xyz): ', xyz(1), xyz(2), xyz(3)
        if( cline%defined('clip') ) goto 998
    else if( cline%defined('e1') .or. (cline%defined('e2') .or. cline%defined('e3')) )then
        call o%new
        call o%set_euler([p%e1,p%e2,p%e3])
        vol2 = b%proj%rotvol(b%vol, o, p)
        b%vol = vol2
    else if( p%phrand .eq. 'yes' )then
        if( .not. cline%defined('smpd') ) stop 'smpd needed 4 phase randomization'
        if( .not. cline%defined('lp') )   stop 'low-pass limit needed 4 phase randomization'
        call b%vol%phase_rand(p%lp)
    else if( p%vols(2) .ne. '' )then
        if( .not. cline%defined('smpd') ) stop 'need smpd 4 FSC calculation!'
        call vol2%copy(b%vol)
        call vol2%read(p%vols(2))
        if( cline%defined('msk') )then
            call mskimg%disc(vol2%get_ldim(), vol2%get_smpd(), p%msk, npix )
            write(*,'(A,1X,F7.3)') '>>> CORRELATION:', b%vol%real_corr(vol2, mskimg)
        else if( cline%defined('lp') )then
            write(*,'(A,1X,F7.3)') '>>> CORRELATION:', b%vol%corr(vol2, p%lp)
        else
            call b%vol%fsc(vol2, res, corrs)
            do j=1,size(res) 
                write(*,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
            end do
            call get_resolution( corrs, res, fsc05, fsc0143 )
            write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', fsc0143
            write(*,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', fsc05
        endif
    else if( p%guinier .eq. 'yes' )then
        p%bfac = b%vol%guinier_bfac(p%hp, p%lp)
        write(*,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', p%bfac
    else if( p%mskfile .ne. '' )then
        call mskimg%new([p%box,p%box,p%box], p%smpd)
        call mskimg%read(p%mskfile)
        b%vol = b%vol*mskimg
    else if( p%countvox .eq. 'yes' )then
        call b%vol%bin
        write(*,'(A,1X,I7)') '>>> NR OF FOREGROUND VOXELS:', b%vol%nforeground()
    else if( p%bin.eq. 'yes' )then
        call b%vol%bin
    else if( cline%defined('newbox') .or. cline%defined('scale') )then
        call vol2%new([p%newbox,p%newbox,p%newbox],p%smpd)
        call b%vol%fwd_ft
        call vol2%set_ft(.true.)
        if( p%newbox < p%box )then
            call b%vol%clip(vol2)
        else if( p%newbox > p%box )then
            call b%vol%pad(vol2)
        else
            vol2 = b%vol
        endif
        b%vol = vol2
        call b%vol%bwd_ft
        p%box = p%newbox
        if( cline%defined('clip') ) goto 998
    else if( cline%defined('xsh') )then
        call b%vol%shift(p%xsh,p%ysh,p%zsh)
    else if( cline%defined('corner') )then
        call b%vol%corner(p%corner,vol2)
        b%vol = vol2
    endif
else if( p%cube > 0 )then
    if( p%box < p%cube ) stop 'need box > cube side'
    call b%vol%square(p%cube/2)
else
    stop 'ERROR: vol1 does not exist!'
endif
998 if( cline%defined('clip') )then
    call vol2%new([p%clip,p%clip,p%clip],p%smpd)
    if( p%clip < p%box )then 
        call b%vol%clip(vol2)
    else
        if( cline%defined('msk') )then
            call b%vol%stats( 'background', ave, sdev, var, med, p%msk ) 
        else
            call b%vol%stats( 'background', ave, sdev, var, med ) 
        endif
        call b%vol%pad(vol2, backgr=med)
    endif
    b%vol = vol2
endif
if( cline%defined('neg') )   call b%vol%neg
if( cline%defined('snr') )   call b%vol%add_gauran(p%snr)
if( cline%defined('hp') )then
    if( .not. cline%defined('smpd') ) stop 'smpd needed 4 filtering'
    call b%vol%bp(p%hp,0.)
endif
if( cline%defined('lp') .and. p%center .eq. 'no' .and. p%phrand .eq. 'no'  )then
    if( .not. cline%defined('smpd') ) stop 'smpd needed 4 filtering'
    call b%vol%bp(0.,p%lp)
endif
if( cline%defined('bfac') )then
    if( .not. cline%defined('smpd') ) stop 'smpd needed 4 bfactor application'
    call b%vol%apply_bfac(p%bfac)
endif
if( p%soften .eq. 'yes' )then
    call b%vol%cos_edge(p%edge)
endif
if( cline%defined('msk') )then
    if( cline%defined('inner') )then
        if( cline%defined('width') )then
            call b%vol%mask(p%msk, p%msktype, inner=p%inner, width=p%width)
        else
            call b%vol%mask(p%msk, p%msktype, inner=p%inner)
        endif
    else
        call b%vol%mask(p%msk, p%msktype)
    endif
endif
if( cline%defined('mirr') ) call b%vol%mirror(p%mirr)
if( p%outvol .ne. '' ) call b%vol%write(p%outvol, del_if_exists=.true.)
999 call simple_end('**** SIMPLE_VOLOPS NORMAL STOP ****')

!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
!shutting down the timers
call stop_Alltimers_cpu()

contains

    !>  \brief  for generating unique character strings for two groups
    function make_signs_two_groups(ng1, ng2 ) result( signs )
        integer, intent(in)   :: ng1, ng2
        character(len=STDLEN), allocatable :: signs(:)
        integer :: alloc_stat, i, cnt
        if( allocated(signs) ) deallocate(signs)
        allocate( signs(ng1+ng2), stat=alloc_stat)
        call alloc_err('make_signs_two_groups; simple_jiffys', alloc_stat)
        cnt = 0
        do i=1,ng1
            cnt = cnt+1
            signs(cnt) = '*'//int2str(i)//'*'
        end do
        do i=1,ng2
            cnt = cnt+1
            signs(cnt) = '@'//int2str(i)//'@'
        end do
    end function
    
end program simple_volops

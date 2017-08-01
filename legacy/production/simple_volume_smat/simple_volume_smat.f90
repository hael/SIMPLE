!==Program simple_volume_smat
!
! <volume\_smat/begin> is a program for creating a similarity matrix based on volume2volume
! correlation. <clin\_smat/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!
program simple_volume_smat
use simple_defs     ! singleton
use simple_jiffys,  ! singleton
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_image,   only: image
implicit none
type(params), target               :: p
type(cmdline)                      :: cline
integer                            :: funit, io_stat, cnt, npairs, npix, nvols
integer                            :: ivol, jvol, ldim(3), alloc_stat, ipair, ifoo
real, allocatable                  :: corrmat(:,:), corrs(:)
integer, allocatable               :: pairs(:,:)
type(image)                        :: vol1, vol2, mskimg
character(len=STDLEN), allocatable :: vollist(:)
character(len=:), allocatable      :: fname
logical, parameter                 :: debug=.false.
if( command_argument_count() < 3 )then
    write(*,'(a)', advance='no') 'SIMPLE_VOLUME_SMAT vollist=<list of volumes> smpd=<sampling distance(in A)>'
    write(*,'(a)', advance='no') ' [lp=<low-pass limit(in A)>] [msk=<mask radius(in pixels)>] [hp=<high-pass'
    write(*,'(a)') ' limit(in A)>] [nthr=<nr of OpenMP threads{1}>]'
    stop
endif
call cline%parse
call cline%checkvar('vollist', 1)
call cline%checkvar('smpd',    2)
call cline%check    ! checks args and prints cmdline to cmdline.dat
p = params(cline, .false.) ! constants & derived constants produced
nvols  = nlines(p%vollist)
npairs = (nvols*(nvols-1))/2
! read in list of volumes
allocate(vollist(nvols))
funit = get_fileunit()
open(unit=funit, status='old', file=p%vollist)
do ivol=1,nvols
    read(funit,'(a256)') vollist(ivol)
    if( debug ) write(*,*) 'read volume: ', vollist(ivol)
end do
close(unit=funit)
! find logical dimension & make volumes for matching
call find_ldim_nptcls(vollist(1), ldim, ifoo)
if( debug ) write(*,*) 'found logical dimension: ', ldim
call vol1%new(ldim,p%smpd)
call vol2%new(ldim,p%smpd)
if( debug ) write(*,*) 'allocated volumes'
if( cline%defined('part') )then
    npairs = p%top-p%fromp+1
    if( debug ) print *, 'allocating this number of similarities: ', npairs
    allocate(corrs(p%fromp:p%top), pairs(p%fromp:p%top,2), stat=alloc_stat)
    call alloc_err('In: simple_comlin_smat, 1', alloc_stat)
    ! read the pairs
    funit = get_fileunit()
    allocate(fname, source='pairs_part'//int2str_pad(p%part,p%numlen)//'.bin')
    if( .not. file_exists(fname) )then
        write(*,*) 'file: ', fname, 'does not exist!'
        write(*,*) 'If all pair_part* are not in cwd, please execute simple_split_pairs to generate the required files'
        stop 'I/O error; simple_comlin_smat'
    endif
    open(unit=funit, status='OLD', action='READ', file=fname, access='STREAM')
    if( debug ) print *, 'reading pairs in range: ', p%fromp, p%top
    read(unit=funit,pos=1,iostat=io_stat) pairs(p%fromp:p%top,:)
    ! Check if the read was successful
    if( io_stat .ne. 0 )then
        write(*,'(a,i0,2a)') '**ERROR(simple_volume_smat): I/O error ', io_stat, ' when reading file: ', fname
        stop 'I/O error; simple_comlin_smat'
    endif
    close(funit)
    deallocate(fname)
    ! make real-space mask if needed
    if( .not. cline%defined('lp') .and. cline%defined('msk') )then 
        call mskimg%disc(vol1%get_ldim(), p%smpd, p%msk, npix)
    endif
    ! calculate the similarities
    cnt = 0
    do ipair=p%fromp,p%top
        cnt = cnt+1
        call progress(cnt, npairs)
        ivol = pairs(ipair,1)
        jvol = pairs(ipair,2)
        call vol1%read(vollist(ivol))
        call vol2%read(vollist(jvol))
        if( cline%defined('lp') )then
            if( cline%defined('msk') )then
                ! apply a soft-edged mask
                call vol1%mask(p%msk, 'soft')
                call vol2%mask(p%msk, 'soft')
            endif
            corrs(ipair) = vol1%corr(vol2,lp_dyn=p%lp,hp_dyn=p%hp)
        else
            if( cline%defined('msk') )then
                corrs(ipair) = vol1%real_corr(vol2, mskimg)
            else
                corrs(ipair) = vol1%real_corr(vol2)
            endif
        endif
    end do
    if( debug ) print *, 'did set this number of similarities: ', cnt
    ! write the similarities
    funit = get_fileunit()
    allocate(fname, source='similarities_part'//int2str_pad(p%part,p%numlen)//'.bin')
    open(unit=funit, status='REPLACE', action='WRITE', file=fname, access='STREAM')
    write(unit=funit,pos=1,iostat=io_stat) corrs(p%fromp:p%top)
    ! Check if the write was successful
    if( io_stat .ne. 0 )then
        write(*,'(a,i0,2a)') '**ERROR(simple_volume_smat): I/O error ', io_stat, ' when writing to: ', fname
        stop 'I/O error; simple_comlin_smat'
    endif
    close(funit)
    deallocate(fname, corrs)
else
    ! generate similarity matrix
    allocate(corrmat(nvols,nvols))
    corrmat = 1.
    cnt = 0
    do ivol=1,nvols-1
        do jvol=ivol+1,nvols
            cnt = cnt+1
            call progress(cnt, npairs)
            call vol1%read(vollist(ivol))
            if( debug ) write(*,*) 'read vol1: ', vollist(ivol)
            call vol2%read(vollist(jvol))
            if( debug ) write(*,*) 'read vol1: ', vollist(ivol)
            corrmat(ivol,jvol) = vol1%corr(vol2,lp_dyn=p%lp,hp_dyn=p%hp)
            if( debug ) write(*,*) 'corr ', ivol, jvol, corrmat(ivol,jvol)
            corrmat(jvol,ivol) = corrmat(ivol,jvol)
        end do
    end do
    funit = get_fileunit()
    open(unit=funit, status='REPLACE', action='WRITE', file='vol_smat.bin', access='STREAM')
    write(unit=funit,pos=1,iostat=io_stat) corrmat
    if( io_stat .ne. 0 )then
        write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to clin_smat.bin'
        stop 'I/O error; simple_volume_smat'
    endif
    close(funit)
endif
call simple_end('**** SIMPLE_VOL_SMAT NORMAL STOP ****')
end program simple_volume_smat

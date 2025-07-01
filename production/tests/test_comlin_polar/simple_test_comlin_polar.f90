program simple_test_comlin_polar
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_image,            only: image
use simple_projector,        only: projector
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_polarizer,        only: polarizer
use simple_comlin,           only: polar_comlin_pfts, gen_polar_comlins
implicit none
integer,          parameter   :: NPLANES = 300, ORI_IND1 = 10, ORI_IND2 = 55
character(len=:), allocatable :: cmd
type(image),      allocatable :: pad_fplanes(:)
complex,          allocatable :: cmat(:,:), pfts(:,:,:)
complex,          pointer     :: ref_ptrs_even(:,:,:), ref_ptrs_odd(:,:,:)
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_polarizer
type(cmdline)                 :: cline
type(image)                   :: noise, ptcl, ptcl_pad, fplanes(NPLANES), vol, img, fplane_polar
type(oris)                    :: spiral
type(ori)                     :: o1, o2
type(projector)               :: vol_pad
type(polar_fmap)              :: pcomlines(NPLANES,NPLANES)
integer :: ifoo, rc, i, lims(3,2), ithr, box, errflg, kfromto(2)
real    :: ave, sdev, maxv, minv
logical :: mrc_exists
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_comlin_polar smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
    write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Converting .pdb to .mrc...'
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('vol1'   , '1JYX.mrc')
    call cline%set('mskdiam', 180.)
    call cline%set('lp'   ,   3.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('vol1',    3)
call cline%checkvar('mskdiam', 4)
call cline%checkvar('lp',      5)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
call vol%new(p%ldim, p%smpd)
call vol%read(p%vols(1))
call vol%stats('foreground', ave, sdev, maxv, minv)
call spiral%new(NPLANES, is_ptcl=.false.)
call spiral%spiral
call ptcl%new(    [p%box,   p%box,   1],       p%smpd)
call noise%new(   [p%box,   p%box ,  p%box],   p%smpd)
call img%new(     [p%box,   p%box,   1],       p%smpd)
call vol_pad%new( [p%boxpd, p%boxpd, p%boxpd], p%smpd)
call ptcl_pad%new([p%boxpd, p%boxpd, 1],       p%smpd)
! add noise in a small center region of the vol
call noise%gauran(0., 5. * sdev)
call noise%mask(1.5 * p%msk, 'soft')
call vol%add(noise)
call vol%pad(vol_pad)
call vol_pad%fft
call vol_pad%expand_cmat(p%alpha)
call spiral%get_ori(ORI_IND1, o1)
call vol_pad%fproject(o1,ptcl_pad)
call ptcl_pad%ifft
call ptcl_pad%clip(ptcl)
call ptcl%write('reproj_com_reprojcom.mrc', 1)
call ptcl%fft
lims = vol%loop_lims(3)
allocate(pad_fplanes(p%nthr))
do ithr = 1, p%nthr
    call pad_fplanes(ithr)%new([p%boxpd, p%boxpd, 1], p%smpd)
enddo
!$omp parallel do default(shared) private(i,ithr,o2)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    ithr = omp_get_thread_num() + 1
    ! fplanes are used for optimization below
    call spiral%get_ori(i, o2)
    call fplanes(i)%new([p%box, p%box, 1], p%smpd)
    call vol_pad%fproject(o2,pad_fplanes(ithr))
    call pad_fplanes(ithr)%ifft
    call pad_fplanes(ithr)%clip(fplanes(i))
enddo
!$omp end parallel do
do i = 1, spiral%get_noris()
    call fplanes(i)%write('fplanes.mrc', i)
    call fplanes(i)%fft
enddo
! polar stuffs
kfromto = p%kfromto
call img_polarizer%new([p%box,p%box,1],p%smpd, wthreads=.false.)
call pftcc%new(NPLANES, [1,NPLANES], kfromto)
call img_polarizer%init_polarizer(pftcc, p%alpha)
call fplane_polar%new([p%box,p%box,1],1.0)
do i = 1, spiral%get_noris()
    call img_polarizer%polarize(pftcc, fplanes(i), i, isptcl=.false., iseven=.true.)
    call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
    call fplane_polar%fft
    call fplane_polar%set_cmat(cmat)
    call fplane_polar%shift_phorig()
    call fplane_polar%ifft
    call fplane_polar%write('fplanes_polar.mrc', i)
enddo
! Testing polar line generation
call pftcc%get_refs_ptr(ref_ptrs_even, ref_ptrs_odd)
allocate(pfts(pftcc%get_pftsz(), kfromto(1):kfromto(2), NPLANES), source=cmplx(0.,0.))
call gen_polar_comlins(pftcc, spiral, pcomlines)
call polar_comlin_pfts(pcomlines, ref_ptrs_even, pfts)
i = ORI_IND1
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 1)
call pftcc%set_ref_pft(i, pfts(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 2)
where( sqrt(real(pfts(:,:,i)*conjg(pfts(:,:,i)))) > TINY ) ref_ptrs_even(:,:,i) = (ref_ptrs_even(:,:,i) + pfts(:,:,i))/2.
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 3)
i = ORI_IND2
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 4)
call pftcc%set_ref_pft(i, pfts(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 5)
where( sqrt(real(pfts(:,:,i)*conjg(pfts(:,:,i)))) > TINY ) ref_ptrs_even(:,:,i) = (ref_ptrs_even(:,:,i) + pfts(:,:,i))/2.
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 6)
end program simple_test_comlin_polar
program simple_test_comlin_polar
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_image,            only: image
use simple_projector,        only: projector
use simple_comlin,           only: polar_comlin_map
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_polarizer,        only: polarizer
implicit none
integer,          parameter   :: NPLANES = 300, ORI_IND1 = 10, ORI_IND2 = 1
character(len=:), allocatable :: cmd
type(image),      allocatable :: pad_fplanes(:)
complex,          allocatable :: cmat(:,:), pft(:,:), pft_line(:)
complex,          pointer     :: ref_ptrs_even(:,:,:), ref_ptrs_odd(:,:,:)
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_polarizer
type(cmdline)                 :: cline
type(image)                   :: noise, ptcl, ptcl_pad, fplanes(NPLANES), vol, img, fplane_polar
type(oris)                    :: spiral
type(ori)                     :: o1, o2
type(projector)               :: vol_pad
integer :: ifoo, rc, i, j, lims(3,2), ithr, box, irot, kind, errflg, irot_l, irot_r, kfromto(2), pftsz
real    :: ave, sdev, maxv, minv, hk1(2), hk2(2), rotmats(3,3,NPLANES), loc1_3D(3), loc2_3D(3), denom,&
          &a1, a2, b1, b2, line_xyz(3), irot_real, invmats(3,3,NPLANES), k_real, w1, w2
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
call spiral%spiral(no_ends=.true.)
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
do i = 1, NPLANES
    call spiral%get_ori(i, o2)
    rotmats(:,:,i) = o2%get_mat()
    call matinv(rotmats(:,:,i), invmats(:,:,i), 3, errflg)
enddo
! Testing polar line generation
pftsz   = pftcc%get_pftsz()
irot    = 5
kind    = kfromto(1) + 2
hk1     = pftcc%get_coord(irot,kind)
irot    = 16
kind    = kfromto(1) + 7
hk2     = pftcc%get_coord(irot,kind)
call pftcc%get_refs_ptr(ref_ptrs_even, ref_ptrs_odd)
allocate( pft(pftsz, kfromto(1):kfromto(2)), pft_line(kfromto(1):kfromto(2)), source=cmplx(0.,0.))
!
i       = ORI_IND1
loc1_3D = matmul([hk1(1), hk1(2), 0.], rotmats(:,:,i))
loc2_3D = matmul([hk2(1), hk2(2), 0.], rotmats(:,:,i))
denom   = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
a1      = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
b1      = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
! images generated by weighted polar common lines
pft     = cmplx(0.,0.)
do j = 1, NPLANES
    if( j == i )cycle
    ! getting the 3D common line
    loc1_3D       = matmul([hk1(1), hk1(2), 0.], rotmats(:,:,j))
    loc2_3D       = matmul([hk2(1), hk2(2), 0.], rotmats(:,:,j))
    denom         = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
    a2            = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
    b2            = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
    line_xyz(1:2) = [1., -(a1-a2)/(b1-b2)]
    line_xyz(3)   = a2*line_xyz(1) + b2*line_xyz(2)
    ! reproject the 3D common line to the polar line on the j-th reference
    line_xyz      = matmul(line_xyz, invmats(:,:,j))
    call pftcc%get_polar_coord(line_xyz(1:2), irot_real, k_real)
    if( irot_real < 1. ) irot_real = irot_real + real(pftsz)
    ! compute the interpolated polar common line, between irot_j and irot_j+1
    irot_l = floor(irot_real)
    irot_r = irot_l + 1
    w2     = irot_real - real(irot_l)
    if( irot_l > pftsz ) irot_l = irot_l - pftsz
    if( irot_r > pftsz ) irot_r = irot_r - pftsz
    pft_line = (1.-w2) * ref_ptrs_even(irot_l,:,j) + w2 * ref_ptrs_even(irot_r,:,j)
    ! reproject the 3D common line to the polar line on the i-th reference
    line_xyz(1:2) = [1., -(a1-a2)/(b1-b2)]
    line_xyz(3)   = a1*line_xyz(1) + b1*line_xyz(2)
    line_xyz      = matmul(line_xyz, invmats(:,:,i))
    call pftcc%get_polar_coord(line_xyz(1:2), irot_real, k_real)
    if( irot_real < 1. ) irot_real = irot_real + real(pftsz)
    ! extrapolate the interpolated polar common line to irot_i and irot_i+1 of i-th reference
    irot_l = floor(irot_real)
    irot_r = irot_l + 1
    w1     = irot_real - real(irot_l)
    if( irot_l > pftsz ) irot_l = irot_l - pftsz
    if( irot_r > pftsz ) irot_r = irot_r - pftsz
    pft(irot_l,:) = pft(irot_l,:) + (1.-w1) * pft_line
    pft(irot_r,:) = pft(irot_r,:) +     w1  * pft_line
enddo
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i), iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 1)
call pftcc%set_ref_pft(i, pft, iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 2)
call pftcc%set_ref_pft(i, ref_ptrs_even(:,:,i) + pft, iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 3)
end program simple_test_comlin_polar
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
integer,          parameter   :: NPLANES = 100, ORI_IND1 = 10, ORI_IND2 = 1
character(len=:), allocatable :: cmd
type(polar_fmap), allocatable :: polar_coords(:)
type(image),      allocatable :: pad_fplanes(:)
complex,          allocatable :: cmat(:,:), pft(:,:), cur_pft(:,:), pft_line(:)
type(polar_fmap)              :: polar_map(NPLANES)
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_polarizer
type(cmdline)                 :: cline
type(image)                   :: noise, ptcl, ptcl_pad, fplanes(NPLANES), vol, img, fplane_polar
type(oris)                    :: spiral
type(ori)                     :: o1, o2
type(projector)               :: vol_pad
integer :: ifoo, rc, i, j, f_ind, lims(3,2), ithr, box, pdim(3), ndim, ori_pcoord(2), tar_pcoord(2),&
          &irot, kind, errflg, irot1, irot2, kfromto(2)
real    :: ave, sdev, maxv, minv, hk(2), e1_rotmat(3,3), loc1_3D(3), loc2_3D(3), denom, a1, a2, b1, b2, line_xyz(3),&
          &irot_real, e1_inv(3,3), k_real, w1, w2
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
call img_polarizer%new([p%box,p%box,1],p%smpd, wthreads=.false.)
call pftcc%new(NPLANES, [1,NPLANES], p%kfromto)
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
pdim = pftcc%get_pdim()
! print *, 'pdim = ', pdim
ndim = pdim(1) * 2 * (pdim(3) - pdim(2))
allocate(polar_coords(p%nthr))
do ithr = 1, p%nthr
    allocate(polar_coords(ithr)%tar_find(  ndim),&
            &polar_coords(ithr)%ori_inds(2,ndim),&
            &polar_coords(ithr)%tar_inds(2,ndim))
enddo
!$omp parallel do default(shared) private(i,ithr)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    ithr = omp_get_thread_num() + 1
    ! common line mapping is independent of fplanes below
    call polar_comlin_map(lims, i, spiral, pftcc, polar_map(i), polar_coords(ithr))
    ! print *, 'i = ', i, '; npoints = ', polar_map(i)%n_points
enddo
!$omp end parallel do
! images generated by polar common lines only
i       = ORI_IND1
kfromto = pftcc%get_kfromto()
allocate(pft(pftcc%get_pftsz(), kfromto(1):kfromto(2)), cur_pft(pftcc%get_pftsz(), kfromto(1):kfromto(2)))
call pftcc%get_ref_pft( i, .true., pft )
pft = cmplx(0.,0.)
do j = 1, polar_map(i)%n_points
    f_ind      = polar_map(i)%tar_find(j)
    if( f_ind == i )cycle
    call pftcc%get_ref_pft( f_ind, .true., cur_pft )
    ori_pcoord = polar_map(i)%ori_inds(:,j)
    tar_pcoord = polar_map(i)%tar_inds(:,j)
    pft(ori_pcoord(1), ori_pcoord(2)) = pft(ori_pcoord(1), ori_pcoord(2)) + cur_pft(tar_pcoord(1), tar_pcoord(2))
enddo
call pftcc%set_ref_pft(i, pft, iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 1)
! Testing polar line generation
e1_rotmat = o1%get_mat()
irot      = 5
kind      = p%kfromto(1) + 2
hk        = pftcc%get_coord(irot,kind)
loc1_3D   = matmul([hk(1), hk(2), 0.], e1_rotmat)
irot      = 16
kind      = p%kfromto(1) + 7
hk        = pftcc%get_coord(irot,kind)
loc2_3D   = matmul([hk(1), hk(2), 0.], e1_rotmat)
denom     = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
a1        = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
b1        = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
! images generated by weighted polar common lines
i   = ORI_IND1
call pftcc%get_ref_pft( i, .true., pft )
pft = cmplx(0.,0.)
allocate(pft_line(kfromto(1):kfromto(2)), source=cmplx(0.,0.))
do j = 1, NPLANES
    if( j == i )cycle
    call spiral%get_ori(j, o2)
    e1_rotmat = o2%get_mat()
    irot      = 5
    kind      = p%kfromto(1) + 2
    hk        = pftcc%get_coord(irot,kind)
    loc1_3D   = matmul([hk(1), hk(2), 0.], e1_rotmat)
    irot      = 16
    kind      = p%kfromto(1) + 7
    hk        = pftcc%get_coord(irot,kind)
    loc2_3D   = matmul([hk(1), hk(2), 0.], e1_rotmat)
    denom     = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
    a2        = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
    b2        = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
    !
    line_xyz(1:2) = [1., -(a1-a2)/(b1-b2)]
    line_xyz(3)   = a1*line_xyz(1) + b1*line_xyz(2)
    e1_rotmat     = o1%get_mat()
    call matinv(e1_rotmat, e1_inv, 3, errflg)
    line_xyz      = matmul(line_xyz, e1_inv)
    call pftcc%get_polar_coord(line_xyz(1:2), irot_real, k_real)
    if( irot_real < 1. ) irot_real = irot_real + real(pftcc%get_pftsz())
    irot1 = floor(irot_real)
    w1    = irot_real - real(irot1)
    !
    line_xyz(1:2) = [1., -(a1-a2)/(b1-b2)]
    line_xyz(3)   = a2*line_xyz(1) + b2*line_xyz(2)
    e1_rotmat     = o2%get_mat()
    call matinv(e1_rotmat, e1_inv, 3, errflg)
    line_xyz      = matmul(line_xyz, e1_inv)
    call pftcc%get_polar_coord(line_xyz(1:2), irot_real, k_real)
    if( irot_real < 1. ) irot_real = irot_real + real(pftcc%get_pftsz())
    irot2 = floor(irot_real)
    w2    = irot_real - real(irot2)
    !
    call pftcc%get_ref_pft( j, .true., cur_pft )
    if( irot2 > pftcc%get_pftsz() ) irot2 = irot2 - pftcc%get_pftsz()
    pft_line = (1.-w2) * cur_pft(irot2,:)
    irot2    = irot2 + 1
    if( irot2 > pftcc%get_pftsz() ) irot2 = irot2 - pftcc%get_pftsz()
    pft_line = pft_line + w2 * cur_pft(irot2,:)
    if( irot1 > pftcc%get_pftsz() ) irot1 = irot1 - pftcc%get_pftsz()
    pft(irot1,:) = pft(irot1,:) + (1.-w1) * pft_line
    irot1        = irot1 + 1
    if( irot1 > pftcc%get_pftsz() ) irot1 = irot1 - pftcc%get_pftsz()
    pft(irot1,:) = pft(irot1,:) +     w1  * pft_line
enddo
call pftcc%set_ref_pft(i, pft, iseven=.true.)
call pftcc%polar2cartesian(i, .true., cmat, box, box_in=p%box)
call fplane_polar%fft
call fplane_polar%set_cmat(cmat)
call fplane_polar%shift_phorig()
call fplane_polar%ifft
call fplane_polar%write('polar_comlin.mrc', 2)
end program simple_test_comlin_polar
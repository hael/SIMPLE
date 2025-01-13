program simple_test_common_lines
include 'simple_lib.f08'
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_image,            only: image
use simple_projector,        only: projector
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_polarizer,        only: polarizer
implicit none
integer,          parameter   :: NPLANES = 50, ORI_IND = 15
character(len=:), allocatable :: cmd
complex,          allocatable :: cmat(:,:)
type(fplan_map)               :: coord_map(NPLANES)
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_polarizer
type(cmdline)                 :: cline
type(image)                   :: noise, ptcl, ptcl_pad, fplanes(NPLANES), pad_fplanes(NPLANES), ptcl_polar
type(oris)                    :: spiral
type(ori)                     :: o1, o2
type(projector)               :: vol_pad, vol
integer                       :: ifoo, rc, errflg, i, j, ori_phys(3), target_phys(3), f_ind, box
real                          :: res_fsc05, res_fsc0143, ave, sdev, maxv, minv, med
complex                       :: diff
logical                       :: mrc_exists
real                          :: vec(1,3), A(3,3), vec_A(1,3), A_inv(3,3), inv_vec_A(1,3), total_costs(NPLANES)
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_common_lines smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
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
call vol_pad%new( [p%boxpd, p%boxpd, p%boxpd], p%smpd)
call ptcl_pad%new([p%boxpd, p%boxpd, 1],       p%smpd)
call noise%new(   [p%boxpd, p%boxpd, 1],       p%smpd)
call vol%pad(vol_pad)
print *, 'lims = ', vol%loop_lims(2)
print *, 'ldim = ', vol%get_ldim()
call vol_pad%fft
call vol_pad%expand_cmat(p%alpha)
call spiral%get_ori(ORI_IND, o1)
call vol_pad%fproject(o1,ptcl_pad)
print *, 'lims_pad = ', vol_pad%loop_lims(2)
print *, 'ldim_pad = ', vol_pad%get_ldim()
call ptcl_pad%ifft
! add noise in a small center region of the vol
call noise%gauran(0., 0.05 * sdev)
call noise%mask(p%msk, 'soft')
call ptcl_pad%add(noise)
call ptcl_pad%clip(ptcl)
call ptcl%write('test_images.mrc', 1)
do i = 1, spiral%get_noris()
    call spiral%get_ori(i, o2)
    call fplanes(i)%new(    [p%box,   p%box,   1], p%smpd)
    call pad_fplanes(i)%new([p%boxpd, p%boxpd, 1], p%smpd)
    call vol_pad%fproject(o2,pad_fplanes(i))
    call pad_fplanes(i)%ifft
    call pad_fplanes(i)%clip(fplanes(i))
    call fplanes(i)%fft
    call vol%fproject_map(i, spiral, coord_map(i))
enddo
call ptcl_pad%fft
call ptcl%fft
total_costs = 0.
do i = 1, spiral%get_noris()
    do j = 1, coord_map(i)%n_points
        f_ind          = coord_map(i)%target_find(j)
        if( f_ind == i ) cycle
        ori_phys       = coord_map(i)%ori_phys(:,j)
        target_phys    = coord_map(i)%target_phys(:,j)
        diff           = ptcl%get_cmat_at(ori_phys) - fplanes(f_ind)%get_cmat_at(target_phys)
        total_costs(i) = total_costs(i) + diff * conjg(diff)
    enddo
    total_costs(i) = sqrt(total_costs(i) / real(coord_map(i)%n_points))
enddo
print *, 'original ptcl ind = ', ORI_IND
print *, 'searched ptcl ind = ', minloc(total_costs)
print *, total_costs



! call vol_pad%fproject_map(ORI_IND, spiral, coord_map)
! do i = 1, size(coord_map)
!     ori_phys    = coord_map(i)%ori_phys
!     target_phys = coord_map(i)%target_phys
!     f_ind       = coord_map(i)%target_find
!     call fplane2_pad%set_cmat_at(ori_phys, fplanes(f_ind)%get_cmat_at(target_phys))
! enddo
! polar stuffs
call img_polarizer%new([p%box,p%box,1],p%smpd, wthreads=.false.)
call pftcc%new(NPLANES, [1,NPLANES], p%kfromto)
call img_polarizer%init_polarizer(pftcc, p%alpha)
call ptcl%fft
call img_polarizer%cartesian2polar(pftcc, ptcl, ORI_IND, isptcl=.false., iseven=.true.)
call pftcc%polar2cartesian(ORI_IND, .true., cmat, box, box_in=p%box)
call ptcl_polar%new([box,box,1],1.0)
call ptcl_polar%set_cmat(cmat)
call ptcl_polar%shift_phorig()
call ptcl_polar%ifft
call ptcl_polar%write('test_images.mrc', 2)
call img_polarizer%polarize(pftcc, ptcl, ORI_IND, isptcl=.false., iseven=.true.)
call pftcc%polar2cartesian(ORI_IND, .true., cmat, box, box_in=p%box)
call ptcl_polar%fft
call ptcl_polar%set_cmat(cmat)
call ptcl_polar%shift_phorig()
call ptcl_polar%ifft
call ptcl_polar%write('test_images.mrc', 3)
! testing
A(1,:)   = [1., 2., 3.]
A(2,:)   = [4., 5., 6.]
A(3,:)   = [7., 9., 10.]
vec(1,:) = [1., 2., 0.]
vec_A = matmul(vec, A)
call matinv(A, A_inv, 3, errflg)
inv_vec_A = matmul(vec_A, A_inv)
end program simple_test_common_lines
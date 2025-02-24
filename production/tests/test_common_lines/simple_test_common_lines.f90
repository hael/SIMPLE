program simple_test_common_lines
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_projector,  only: projector
use simple_comlin,     only: comlin_map
implicit none
integer,          parameter   :: NPLANES = 50, ORI_IND = 15
character(len=:), allocatable :: cmd
type(fplan_map)  :: all_coords
type(fplan_map)  :: coord_map(NPLANES)
type(parameters) :: p
type(cmdline)    :: cline
type(image)      :: noise, ptcl, ptcl_pad, fplanes(NPLANES), pad_fplane, vol
type(oris)       :: spiral
type(ori)        :: o1, o2
type(projector)  :: vol_pad
integer          :: ifoo, rc, i, j, ori_phys(3), target_phys(3), f_ind, lims(3,2)
real             :: ave, sdev, maxv, minv, total_costs(NPLANES)
complex          :: diff
logical          :: mrc_exists

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
print *, 'box  = ', p%box
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
call pad_fplane%new([p%boxpd, p%boxpd, 1], p%smpd)
lims = vol%loop_lims(2)
allocate(all_coords%target_find(  (lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
        &all_coords%ori_phys(   3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
        &all_coords%target_phys(3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)))
do i = 1, spiral%get_noris()
    ! common line mapping is independent of fplanes below
    call comlin_map(lims, i, spiral, coord_map(i), all_coords)
    ! fplanes are used for optimization below
    call spiral%get_ori(i, o2)
    call fplanes(i)%new([p%box, p%box, 1], p%smpd)
    call vol_pad%fproject(o2,pad_fplane)
    call pad_fplane%ifft
    call pad_fplane%clip(fplanes(i))
    call fplanes(i)%fft
enddo
call ptcl%fft
total_costs = 0.
do i = 1, spiral%get_noris()
    do j = 1, coord_map(i)%n_points
        f_ind          = coord_map(i)%target_find(j)
        if( f_ind == i ) cycle
        ori_phys       = coord_map(i)%ori_phys(:,j)
        target_phys    = coord_map(i)%target_phys(:,j)
        diff           = ptcl%get_cmat_at(ori_phys) - fplanes(f_ind)%get_cmat_at(target_phys)
        total_costs(i) = total_costs(i) + real(diff * conjg(diff))
    enddo
    total_costs(i) = sqrt(total_costs(i) / real(coord_map(i)%n_points))
enddo
print *, 'truth    ptcl ind = ', ORI_IND
print *, 'searched ptcl ind = ', minloc(total_costs)
print *, total_costs

! call vol_pad%fproject_map(ORI_IND, spiral, coord_map)
! do i = 1, size(coord_map)
!     ori_phys    = coord_map(i)%ori_phys
!     target_phys = coord_map(i)%target_phys
!     f_ind       = coord_map(i)%target_find
!     call fplane2_pad%set_cmat_at(ori_phys, fplanes(f_ind)%get_cmat_at(target_phys))
! enddo

end program simple_test_common_lines
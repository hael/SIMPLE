program simple_test_common_lines
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_projector,  only: projector
use simple_comlin,     only: comlin_map
implicit none
integer,          parameter   :: NPLANES = 50, ORI_IND1 = 15, ORI_IND2 = 43
character(len=:), allocatable :: cmd
type(fplan_map),  allocatable :: all_coords(:)
type(image),      allocatable :: pad_fplanes(:)
type(fplan_map)  :: coord_map(NPLANES)
type(parameters) :: p
type(cmdline)    :: cline
type(image)      :: noise, ptcl, ptcl_pad, fplanes(NPLANES), vol, img
type(oris)       :: spiral
type(ori)        :: o1, o2
type(projector)  :: vol_pad
integer          :: ifoo, rc, i, j, k, ori_phys(3), target_phys(3), f_ind, lims(3,2), ithr, cnts(NPLANES,NPLANES), min_i, min_k
real             :: ave, sdev, maxv, minv, total_costs(NPLANES), pair_costs(NPLANES,NPLANES), minval
complex          :: diff, val
logical          :: mrc_exists, found(NPLANES,NPLANES)

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
call noise%new(   [p%box,   p%box ,  1],       p%smpd)
call img%new(     [p%box,   p%box,   1],       p%smpd)
call vol_pad%new( [p%boxpd, p%boxpd, p%boxpd], p%smpd)
call ptcl_pad%new([p%boxpd, p%boxpd, 1],       p%smpd)
call vol%pad(vol_pad)
call vol_pad%fft
call vol_pad%expand_cmat(p%alpha)
call spiral%get_ori(ORI_IND1, o1)
call vol_pad%fproject(o1,ptcl_pad)
call ptcl_pad%ifft
! add noise in a small center region of the vol
call noise%gauran(0., 0.1 * sdev)
call noise%mask(p%msk, 'soft')
call ptcl_pad%clip(ptcl)
call ptcl%add(noise)
call ptcl%write('reproj_com_reprojcom.mrc', 1)
call ptcl%fft
lims = vol%loop_lims(2)
allocate(pad_fplanes(p%nthr),all_coords(p%nthr))
do ithr = 1, p%nthr
    call pad_fplanes(ithr)%new([p%boxpd, p%boxpd, 1], p%smpd)
    allocate(all_coords(ithr)%target_find(  (lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%ori_phys(   3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%ori_four(   3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%target_phys(3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%target_four(3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)))
enddo
!$omp parallel do default(shared) private(i,ithr,o2)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    ithr = omp_get_thread_num() + 1
    ! common line mapping is independent of fplanes below
    call comlin_map(lims, i, spiral, coord_map(i), all_coords(ithr))
    ! fplanes are used for optimization below
    call spiral%get_ori(i, o2)
    call fplanes(i)%new([p%box, p%box, 1], p%smpd)
    call vol_pad%fproject(o2,pad_fplanes(ithr))
    call pad_fplanes(ithr)%ifft
    call pad_fplanes(ithr)%clip(fplanes(i))
    call fplanes(i)%fft
enddo
!$omp end parallel do
total_costs = 0.
!$omp parallel do default(shared) private(i,j,f_ind,ori_phys,target_phys,diff)&
!$omp proc_bind(close) schedule(static)
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
!$omp end parallel do
print *, 'truth    ptcl ind = ', ORI_IND1
print *, 'searched ptcl ind = ', minloc(total_costs)
! common line detection between two fplanes
pair_costs = 0.
cnts       = 0
found      = .false.
!$omp parallel do default(shared) private(i,j,k,f_ind,ori_phys,target_phys,diff)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    do k = 1, spiral%get_noris()
        if( k == i )then
            pair_costs(i,k) = huge(ave)
            cnts(i,k)       = 1
            cycle
        endif
        do j = 1, coord_map(i)%n_points
            f_ind           = coord_map(i)%target_find(j)
            if( f_ind /= k ) cycle
            ori_phys        = coord_map(i)%ori_phys(:,j)
            target_phys     = coord_map(i)%target_phys(:,j)
            if( .not.(check_kfromto(coord_map(i)%ori_four(1:2,j)) .and. check_kfromto(coord_map(i)%target_four(1:2,j))) ) cycle
            diff            = fplanes(ORI_IND1)%get_cmat_at(ori_phys) - fplanes(ORI_IND2)%get_cmat_at(target_phys)
            pair_costs(i,k) = pair_costs(i,k) + real(diff * conjg(diff))
            cnts(i,k)       = cnts(i,k) + 1
            found(i,k)      = .true.
        enddo
    enddo
enddo
!$omp end parallel do
minval = huge(ave)
do i = 1, spiral%get_noris()
    do k = 1, spiral%get_noris()
        if( k == i .or. .not.found(i,k) .or. cnts(i,k) < p%kfromto(2)-p%kfromto(1) )cycle
        if( pair_costs(i,k)/real(cnts(i,k)) < minval )then
            minval = pair_costs(i,k)/real(cnts(i,k))
            min_i  = i
            min_k  = k
        endif
    enddo
enddo
print *, 'BOX                    = ', p%box
print *, 'kfromto                = ', p%kfromto
print *, 'ORI_IND1 ind           = ', ORI_IND1
print *, 'ORI_IND2 ind           = ', ORI_IND2
print *, 'searched ORI_IND1 ind  = ', min_i
print *, 'searched ORI_IND2 ind  = ', min_k
print *, 'minval                 = ', minval
print *, 'cnt                    = ', cnts(min_i,min_k)
print *, 'min(ORI_IND1,ORI_IND2) = ', pair_costs(ORI_IND1,ORI_IND2)/real(cnts(ORI_IND1,ORI_IND2))
print *, 'cnt                    = ', cnts(ORI_IND1,ORI_IND2)
! images generated by common lines only
i = minloc(total_costs, dim=1)
call img%zero_and_flag_ft
do j = 1, coord_map(i)%n_points
    f_ind       = coord_map(i)%target_find(j)
    ori_phys    = coord_map(i)%ori_phys(:,j)
    target_phys = coord_map(i)%target_phys(:,j)
    if( .not.(check_kfromto(coord_map(i)%ori_four(1:2,j)) .and. check_kfromto(coord_map(i)%target_four(1:2,j))) ) cycle
    val         = img%get_cmat_at(ori_phys)
    call img%set_cmat_at(ori_phys, val + fplanes(f_ind)%get_cmat_at(target_phys))
enddo
call img%ifft
call img%write('reproj_com_reprojcom.mrc', 2)
! noisy img + img from common lines
call img%zero_and_flag_ft
call img%set_cmat(ptcl%get_cmat())
do j = 1, coord_map(i)%n_points
    f_ind       = coord_map(i)%target_find(j)
    ori_phys    = coord_map(i)%ori_phys(:,j)
    target_phys = coord_map(i)%target_phys(:,j)
    if( .not.(check_kfromto(coord_map(i)%ori_four(1:2,j)) .and. check_kfromto(coord_map(i)%target_four(1:2,j))) ) cycle
    val         = img%get_cmat_at(ori_phys)
    call img%set_cmat_at(ori_phys, val + fplanes(f_ind)%get_cmat_at(target_phys))
enddo
call img%ifft
call img%write('reproj_com_reprojcom.mrc', 3)

contains

    function check_kfromto(addr) result(okay)
        integer, intent(in) :: addr(2)
        logical :: okay
        integer :: sh
        sh   = nint(hyp(addr(1),addr(2)))
        okay = .false.
        if( sh >= p%kfromto(1) .and. sh <= p%kfromto(2) ) okay = .true.
    end function check_kfromto

end program simple_test_common_lines
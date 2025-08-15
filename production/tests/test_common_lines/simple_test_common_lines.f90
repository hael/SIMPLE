program simple_test_common_lines
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_image,      only: image
use simple_projector,  only: projector
use simple_comlin,     only: comlin_map, fplan_map
implicit none
integer,          parameter   :: NPLANES = 20, ORI_IND1 = 3, ORI_IND2 = 15, NTHETAS = 18
character(len=:), allocatable :: cmd
type(fplan_map),  allocatable :: all_coords(:)
type(image),      allocatable :: pad_fplanes(:)
real,             allocatable :: thetas(:)
type(fplan_map)  :: coord_map(NPLANES)
type(parameters) :: p
type(cmdline)    :: cline
type(image)      :: noise, ptcl, ptcl_pad, fplanes(NPLANES), vol, img, rot_fplane
type(oris)       :: spiral
type(ori)        :: o1, o2
type(projector)  :: vol_pad
integer          :: ifoo, rc, i, j, k, ori_phys(3), tar_phys(3), f_ind, lims(3,2), ithr, cnts(NPLANES,NPLANES),&
                   &min_i, min_k, itheta, ori_four(2), tar_four(2), errflg, rot_four(2), itheta1, itheta2, jtheta, theta_ori, theta_tar
real             :: ave, sdev, maxv, minv, total_costs(NPLANES), all_costs(NTHETAS,NPLANES,NPLANES), minval, &
                   &shifts(3), rval, theta, pair_costs(NPLANES,NPLANES), found_thetas(NPLANES,NPLANES), mat(2,2), mat_inv(2,2),&
                   &theta1, theta2, the_costs(NTHETAS,NTHETAS)
complex          :: diff, val, val1, val2
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
lims = vol%loop_lims(3)
allocate(pad_fplanes(p%nthr),all_coords(p%nthr))
do ithr = 1, p%nthr
    call pad_fplanes(ithr)%new([p%boxpd, p%boxpd, 1], p%smpd)
    allocate(all_coords(ithr)%tar_find(  (lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%ori_phys(3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%ori_four(2,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%tar_phys(3,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)),&
            &all_coords(ithr)%tar_four(2,(lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)))
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
enddo
!$omp end parallel do
do i = 1, spiral%get_noris()
    call fplanes(i)%write('fplanes.mrc', i)
    call fplanes(i)%fft
enddo
total_costs = 0.
!$omp parallel do default(shared) private(i,j,f_ind,ori_phys,tar_phys,diff)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    do j = 1, coord_map(i)%n_points
        f_ind          = coord_map(i)%tar_find(j)
        if( f_ind == i ) cycle
        ori_phys       = coord_map(i)%ori_phys(:,j)
        tar_phys       = coord_map(i)%tar_phys(:,j)
        diff           = ptcl%get_cmat_at(ori_phys) - fplanes(f_ind)%get_cmat_at(tar_phys)
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
! rotating first fplane
itheta1    = 5
itheta     = 5
thetas     = real((/(i,i=0,360/20-1)/)) * 20.
theta1     = thetas(itheta1)
call rot_fplane%new([p%box, p%box, 1], p%smpd)
call rotate_img(theta1, fplanes(ORI_IND1), rot_fplane)
call fplanes(ORI_IND1)%ifft
call fplanes(ORI_IND1)%write('ori_fplane1.mrc', 1)
call fplanes(ORI_IND1)%fft
call rot_fplane%ifft
call rot_fplane%write('rot_fplane1.mrc', 1)
call rot_fplane%fft
! cost function check
i = ORI_IND1
k = ORI_IND2
itheta = 1
minval = 0.
do j = 1, coord_map(i)%n_points
    f_ind    = coord_map(i)%tar_find(j)
    if( f_ind /= k .or. f_ind == i ) cycle
    ori_four = coord_map(i)%ori_four(:,j)
    tar_four = coord_map(i)%tar_four(:,j)
    if( .not.(check_kfromto(ori_four) .and. check_kfromto(tar_four)) ) cycle
    tar_phys = coord_map(i)%tar_phys(:,j)
    call rotmat2d(thetas(itheta), mat)
    rot_four      = nint(matmul(real(ori_four),mat))
    if( .not.(check_kfromto(rot_four)) ) cycle
    ori_phys(1:2) = fplanes(ORI_IND1)%comp_addr_phys(rot_four(1),rot_four(2))
    ori_phys(3)   = 1
    val1          = fplanes(ORI_IND1)%get_cmat_at(ori_phys)
    val2          = fplanes(ORI_IND2)%get_cmat_at(tar_phys)
    if( ori_four(1) < 0 ) val1 = conjg(val1)
    if( tar_four(1) < 0 ) val2 = conjg(val2)
    diff          = val1   - val2
    minval        = minval + real(diff * conjg(diff))
enddo
print *, 'cost without rotation  = ', minval
! rotate ORI_IND1
call fplanes(ORI_IND1)%copy_fast(rot_fplane)
! cost function check after rotation
i      = ORI_IND1
k      = ORI_IND2
minval = 0.
all_costs(:,i,k) = 0.
do j = 1, coord_map(i)%n_points
    f_ind    = coord_map(i)%tar_find(j)
    if( f_ind /= k .or. f_ind == i ) cycle
    ori_four = coord_map(i)%ori_four(:,j)
    tar_four = coord_map(i)%tar_four(:,j)
    if( .not.(check_kfromto(ori_four) .and. check_kfromto(tar_four)) ) cycle
    tar_phys = coord_map(i)%tar_phys(:,j)
    do itheta = 1, size(thetas)
        call rotmat2d(thetas(itheta), mat)
        rot_four = nint(matmul(real(ori_four),mat))
        if( .not.(check_kfromto(rot_four)) ) cycle
        ori_phys(1:2)         = fplanes(ORI_IND1)%comp_addr_phys(rot_four(1),rot_four(2))
        ori_phys(3)           = 1
        val1                  = fplanes(ORI_IND1)%get_cmat_at(ori_phys)
        val2                  = fplanes(ORI_IND2)%get_cmat_at(tar_phys)
        if( ori_four(1) < 0 ) val1 = conjg(val1)
        if( tar_four(1) < 0 ) val2 = conjg(val2)
        diff                  = val1 - val2
        all_costs(itheta,i,k) = all_costs(itheta,i,k) + real(diff * conjg(diff))
    enddo
enddo
itheta = minloc(all_costs(:,i,k), dim=1)
print *, 'cost with rotation     = ', all_costs(itheta,i,k)
print *, 'truth theta1           = ', thetas(itheta1)
print *, 'found theta1           = ', thetas(itheta)
! both rotated fplanes
itheta2 = 13
theta2  = thetas(itheta2)
call rot_fplane%kill
call rot_fplane%new([p%box, p%box, 1], p%smpd)
call rotate_img(theta2, fplanes(ORI_IND2), rot_fplane)
call fplanes(ORI_IND2)%ifft
call fplanes(ORI_IND2)%write('ori_fplane2.mrc', 1)
call fplanes(ORI_IND2)%fft
call rot_fplane%ifft
call rot_fplane%write('rot_fplane2.mrc', 1)
call rot_fplane%fft
call fplanes(ORI_IND2)%copy_fast(rot_fplane)
! cost function check after rotation
i         = ORI_IND1
k         = ORI_IND2
the_costs = 0.
do j = 1, coord_map(i)%n_points
    f_ind    = coord_map(i)%tar_find(j)
    if( f_ind /= k .or. f_ind == i ) cycle
    ori_four = coord_map(i)%ori_four(:,j)
    tar_four = coord_map(i)%tar_four(:,j)
    if( .not.(check_kfromto(ori_four) .and. check_kfromto(tar_four)) ) cycle
    do itheta = 1, size(thetas)
        call rotmat2d(thetas(itheta), mat)
        rot_four = nint(matmul(real(ori_four),mat))
        if( .not.(check_kfromto(rot_four)) ) cycle
        ori_phys(1:2)         = fplanes(ORI_IND1)%comp_addr_phys(rot_four(1),rot_four(2))
        ori_phys(3)           = 1
        val1                  = fplanes(ORI_IND1)%get_cmat_at(ori_phys)
        if( ori_four(1) < 0 ) val1 = conjg(val1)
        do jtheta = 1, size(thetas)
            call rotmat2d(thetas(jtheta), mat)
            rot_four = nint(matmul(real(tar_four),mat))
            if( .not.(check_kfromto(rot_four)) ) cycle
            tar_phys(1:2) = fplanes(ORI_IND2)%comp_addr_phys(rot_four(1),rot_four(2))
            tar_phys(3)   = 1
            val2          = fplanes(ORI_IND2)%get_cmat_at(tar_phys)
            if( tar_four(1) < 0 ) val2 = conjg(val2)
            diff                     = val1 - val2
            the_costs(itheta,jtheta) = the_costs(itheta,jtheta) + real(diff * conjg(diff))
        enddo
    enddo
enddo
minval    = huge(ave)
theta_ori = 1
theta_tar = 1
do itheta = 1, size(thetas)
    do jtheta = 1, size(thetas)
        if( the_costs(itheta,jtheta) < minval )then
            theta_ori = itheta
            theta_tar = jtheta
            minval    = the_costs(itheta,jtheta)
        endif
    enddo
enddo
print *, 'cost of both rotations = ', minval
print *, 'truth theta1           = ', thetas(itheta1)
print *, 'truth theta2           = ', thetas(itheta2)
print *, 'found theta1           = ', thetas(theta_ori)
print *, 'found theta2           = ', thetas(theta_tar)
! over all projection directions
!$omp parallel do collapse(2) default(shared) private(i,j,k,f_ind,ori_phys,tar_phys,diff,itheta,ori_four,tar_four,mat,rot_four,val1,val2)&
!$omp proc_bind(close) schedule(static)
do i = 1, spiral%get_noris()
    do k = 1, spiral%get_noris()
        if( k == i )then
            pair_costs(i,k) = huge(ave)
            cnts(i,k)       = 1
            cycle
        endif
        all_costs(:,i,k) = 0.
        do j = 1, coord_map(i)%n_points
            f_ind    = coord_map(i)%tar_find(j)
            if( f_ind /= k .or. f_ind == i ) cycle
            ori_four = coord_map(i)%ori_four(:,j)
            tar_four = coord_map(i)%tar_four(:,j)
            if( .not.(check_kfromto(ori_four) .and. check_kfromto(tar_four)) ) cycle
            tar_phys = coord_map(i)%tar_phys(:,j)
            do itheta = 1, size(thetas)
                call rotmat2d(thetas(itheta), mat)
                rot_four = nint(matmul(real(ori_four),mat))
                if( .not.(check_kfromto(rot_four)) ) cycle
                ori_phys(1:2)         = fplanes(ORI_IND1)%comp_addr_phys(rot_four(1),rot_four(2))
                ori_phys(3)           = 1
                val1                  = fplanes(ORI_IND1)%get_cmat_at(ori_phys)
                val2                  = fplanes(ORI_IND2)%get_cmat_at(tar_phys)
                if( ori_four(1) < 0 ) val1 = conjg(val1)
                if( tar_four(1) < 0 ) val2 = conjg(val2)
                diff                  = val1 - val2
                all_costs(itheta,i,k) = all_costs(itheta,i,k) + real(diff * conjg(diff))
            enddo
        enddo
        itheta            = minloc(all_costs(:,i,k), dim=1)
        found_thetas(i,k) = thetas(itheta)
        do j = 1, coord_map(i)%n_points
            f_ind    = coord_map(i)%tar_find(j)
            if( f_ind /= k .or. f_ind == i ) cycle
            ori_four = coord_map(i)%ori_four(:,j)
            tar_four = coord_map(i)%tar_four(:,j)
            if( .not.(check_kfromto(ori_four) .and. check_kfromto(tar_four)) ) cycle
            tar_phys        = coord_map(i)%tar_phys(:,j)
            call rotmat2d(thetas(itheta), mat)
            rot_four        = nint(matmul(real(ori_four),mat))
            if( .not.(check_kfromto(rot_four)) ) cycle
            ori_phys(1:2)   = fplanes(ORI_IND1)%comp_addr_phys(rot_four(1),rot_four(2))
            ori_phys(3)     = 1
            val1            = fplanes(ORI_IND1)%get_cmat_at(ori_phys)
            val2            = fplanes(ORI_IND2)%get_cmat_at(tar_phys)
            if( ori_four(1) < 0 ) val1 = conjg(val1)
            if( tar_four(1) < 0 ) val2 = conjg(val2)
            diff            = val1 - val2
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
        if( k == i .or. .not.found(i,k) .or. .not.found(k,i) .or. cnts(i,k) > p%box )cycle
        rval = min(sqrt(pair_costs(i,k))/real(cnts(i,k)), sqrt(pair_costs(k,i))/real(cnts(k,i)))
        if( rval < minval )then
            minval = rval
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
print *, 'theta                  = ', found_thetas(min_i,min_k)
print *, 'cnt                    = ', cnts(min_i,min_k)
print *, 'min(ORI_IND1,ORI_IND2) = ', sqrt(pair_costs(ORI_IND1,ORI_IND2))/real(cnts(ORI_IND1,ORI_IND2))
print *, 'min(ORI_IND2,ORI_IND1) = ', sqrt(pair_costs(ORI_IND2,ORI_IND1))/real(cnts(ORI_IND2,ORI_IND1))
! images generated by common lines only
i = minloc(total_costs, dim=1)
call img%zero_and_flag_ft
do j = 1, coord_map(i)%n_points
    f_ind    = coord_map(i)%tar_find(j)
    ori_phys = coord_map(i)%ori_phys(:,j)
    tar_phys = coord_map(i)%tar_phys(:,j)
    if( .not.(check_kfromto(coord_map(i)%ori_four(1:2,j)) .and. check_kfromto(coord_map(i)%tar_four(1:2,j))) ) cycle
    val      = img%get_cmat_at(ori_phys)
    call img%set_cmat_at(ori_phys, val + fplanes(f_ind)%get_cmat_at(tar_phys))
enddo
call img%ifft
call img%write('reproj_com_reprojcom.mrc', 2)
! noisy img + img from common lines
call img%zero_and_flag_ft
call img%set_cmat(ptcl%get_cmat())
do j = 1, coord_map(i)%n_points
    f_ind    = coord_map(i)%tar_find(j)
    ori_phys = coord_map(i)%ori_phys(:,j)
    tar_phys = coord_map(i)%tar_phys(:,j)
    if( .not.(check_kfromto(coord_map(i)%ori_four(1:2,j)) .and. check_kfromto(coord_map(i)%tar_four(1:2,j))) ) cycle
    val      = img%get_cmat_at(ori_phys)
    call img%set_cmat_at(ori_phys, val + fplanes(f_ind)%get_cmat_at(tar_phys))
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

    subroutine rotate_img( e3, img_in, img_out )
        real,        intent(in)    :: e3
        type(image), intent(inout) :: img_in, img_out
        complex :: fcompl, fcompll
        real    :: mat(2,2), loc(2), dist(2), kw
        integer :: h, k, logi_lims(3,2), cyc_lims(3,2), cyc_limsR(2,2), win_corner(2), l, ll, m, mm, phys(2)
        logi_lims      = img_in%loop_lims(3)
        cyc_lims       = img_in%loop_lims(3)
        cyc_limsR(:,1) = cyc_lims(1,:)
        cyc_limsR(:,2) = cyc_lims(2,:)
        call rotmat2d(-e3, mat)
        call img_out%zero_and_flag_ft
        do h = logi_lims(1,1),logi_lims(1,2)
            do k = logi_lims(2,1),logi_lims(2,2)
                ! Rotation
                loc        = matmul(real([h,k]),mat)
                win_corner = nint(loc)
                phys       = img_in%comp_addr_phys(win_corner(1), win_corner(2))
                fcompl     = img_in%get_cmat_at(phys(1),phys(2),1)
                phys       = img_in%comp_addr_phys(h,k)
                call img_out%set_cmat_at(phys(1),phys(2),1, fcompl)
            end do
        end do
    end subroutine rotate_img

    subroutine rotate_img_interp( e3, img_in, img_out )
        real,        intent(in)    :: e3
        type(image), intent(inout) :: img_in, img_out
        complex :: fcompl, fcompll
        real    :: mat(2,2), loc(2), dist(2), kw
        integer :: h, k, logi_lims(3,2), cyc_lims(3,2), cyc_limsR(2,2), win_corner(2), l, ll, m, mm, phys(2)
        logi_lims      = img_in%loop_lims(2)
        cyc_lims       = img_in%loop_lims(3)
        cyc_limsR(:,1) = cyc_lims(1,:)
        cyc_limsR(:,2) = cyc_lims(2,:)
        call rotmat2d(-e3, mat)
        call img_out%zero_and_flag_ft
        do h = logi_lims(1,1),logi_lims(1,2)
            do k = logi_lims(2,1),logi_lims(2,2)
                ! Rotation
                loc        = matmul(real([h,k]),mat)
                win_corner = floor(loc) ! bottom left corner
                dist       = loc - real(win_corner)
                ! Bi-linear interpolation
                l       = cyci_1d(cyc_limsR(:,1), win_corner(1))
                ll      = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
                m       = cyci_1d(cyc_limsR(:,2), win_corner(2))
                mm      = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
                ! l, bottom left corner
                phys    = img_in%comp_addr_phys(l,m)
                kw      = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
                fcompl  = kw * img_in%get_cmat_at(phys(1), phys(2),1)
                ! l, bottom right corner
                phys    = img_in%comp_addr_phys(l,mm)
                kw      = (1.-dist(1))*dist(2)
                fcompl  = fcompl + kw * img_in%get_cmat_at(phys(1), phys(2),1)
                if( l < 0 ) fcompl = conjg(fcompl) ! conjugation when required!
                ! ll, upper left corner
                phys    = img_in%comp_addr_phys(ll,m)
                kw      = dist(1)*(1.-dist(2))
                fcompll = kw * img_in%get_cmat_at(phys(1), phys(2),1)
                ! ll, upper right corner
                phys    = img_in%comp_addr_phys(ll,mm)
                kw      = dist(1)*dist(2)
                fcompll = fcompll + kw * img_in%get_cmat_at(phys(1), phys(2),1)
                if( ll < 0 ) fcompll = conjg(fcompll) ! conjugation when required!
                ! update with interpolated values
                phys    = img_in%comp_addr_phys(h,k)
                call img_out%set_cmat_at(phys(1),phys(2),1, fcompl + fcompll)
            end do
        end do
    end subroutine rotate_img_interp

end program simple_test_common_lines
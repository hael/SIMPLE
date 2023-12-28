program simple_test_3D_core_finder
include 'simple_lib.f08'
use simple_image,        only: image
implicit none
#include "simple_local_flags.inc"
real,    parameter  :: smpd=0.358, np_rad_A=50.0, shell_size_A=1.0
real                :: suma, suma2, sumax, sumay
type(image)                 :: img_part1, img_part2, dists_img
real,    allocatable        :: rvec1(:), rvec2(:)
logical, allocatable        :: mask(:,:,:), shell_mask(:,:,:)
real(kind=c_float), pointer :: rmat_img_part1(:,:,:)=>null()
real(kind=c_float), pointer :: rmat_img_part2(:,:,:)=>null()
real(kind=c_float), pointer :: rmat_dists_img(:,:,:)=>null()
real                        :: np_rad_vox=np_rad_A/smpd, shell_size_vox=shell_size_A/smpd
real                        :: rad_core, rad_core_vox, corr
character(len=256)          :: fn_img_part1, fn_img_part2
integer                     :: ldim1(3), ldim2(3), ldim_refs(3), n, nshells, ifoo
real                        :: mean, mean1, sdev1, var1, mean2, sdev2, var2, dist_lbound, dist_ubound
logical                     :: err
if( command_argument_count() /= 2 )then
    write(logfhandle,'(a)')  'Usage: simple_test_3D_core_finder vol_part1.mrc vol_part2.mrc radius_core'
    write(logfhandle,'(a)')  'vol_part1.mrc: 3D volume density map of part1'
    write(logfhandle,'(a)')  'vol_part2.mrc: 3D volume density map of part2', NEW_LINE('a')
    stop
else
    call get_command_argument(1, fn_img_part1)
    call get_command_argument(2, fn_img_part2)
endif
! create and read in volumes
call find_ldim_nptcls(fn_img_part1, ldim1, ifoo)
call find_ldim_nptcls(fn_img_part2, ldim2, ifoo)
if( ldim1(1) /= ldim2(1) .or. ldim1(2) /= ldim2(2) .or. ldim1(3) /= ldim2(3) )then
    THROW_HARD('****Non-conforming dimensions of 3D density maps vol1 and vol2')
endif
ldim_refs = [ldim1(1), ldim1(2), ldim1(3)]
call img_part1%new(ldim_refs, smpd)
call img_part2%new(ldim_refs, smpd)
call dists_img%new(ldim_refs, smpd)
call img_part1%read(fn_img_part1)
call img_part2%read(fn_img_part2)
call dists_img%cendist
call dists_img%get_rmat_ptr( rmat_dists_img )
! allocations and parameter calculations
rad_core_vox  = rad_core/smpd
nshells       = int(np_rad_vox / shell_size_vox)
allocate(mask(ldim_refs(1), ldim_refs(2), ldim_refs(3)),&
  &shell_mask(ldim_refs(1), ldim_refs(2), ldim_refs(3)), source=.true.)
! mask negative values
call img_part1%get_rmat_ptr( rmat_img_part1 )
call img_part2%get_rmat_ptr( rmat_img_part2 )
where( rmat_img_part1(:ldim_refs(1),:ldim_refs(2),:ldim_refs(3)) < 0. ) mask = .false.
where( rmat_img_part2(:ldim_refs(1),:ldim_refs(2),:ldim_refs(3)) < 0. ) mask = .false.
! thresholding
! call img_part1%zero_below(0.)
! call img_part2%zero_below(0.)
! take averages of shells out to the NP radius
do n = 0, nshells
    dist_lbound = real(n)*shell_size_vox
    dist_ubound = dist_lbound + shell_size_vox
    where( (rmat_dists_img(:ldim_refs(1),:ldim_refs(2),:ldim_refs(3)) > dist_lbound) .and.&
          &(rmat_dists_img(:ldim_refs(1),:ldim_refs(2),:ldim_refs(3)) < dist_ubound) .and. mask)
        shell_mask = .true.
    else where
        shell_mask = .false.
    end where
    if( count(shell_mask) < 3 )then
        corr = 0.
        mean = 0.
    else
        rvec1 = img_part1%serialize(shell_mask)
        rvec2 = img_part2%serialize(shell_mask)
        call moment(rvec1, mean1, sdev1, var1, err)
        call moment(rvec2, mean2, sdev2, var2, err)
        mean = mean1 - mean2
        corr = pearsn_serial(rvec1, rvec2)
        ! corr = img_part1%real_corr(img_part2, shell_mask) ! works too
    endif
    print *, n, dist_lbound, dist_ubound, corr, mean
enddo
call img_part1%kill()
call img_part2%kill()
call dists_img%kill()
end program simple_test_3D_core_finder

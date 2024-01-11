program simple_test_picker_comp
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_picker_utils, only: picker_utils
use simple_pickgau, only: pickgau, read_mic_raw
use simple_linalg

implicit none
character(len = 200) :: micname 
real, parameter :: SMPD = 1.72, MOLDIAM = 180., SMPD_SHRINK = 4., MOLDIAM_MAX = 180.
integer, parameter :: OFFSET = 3
integer :: nptcls, nthr, i, iostat, cnt, nlines_pickgau, nlines_picker_utils
integer, allocatable :: coords_pickgau(:,:), coords_picker_utils(:,:)
logical :: match
character(len=LONGSTRLEN) :: boxname_out, cwd
character(len=:), allocatable :: dir_out
type(picker_utils) :: putils
type(pickgau) :: pgau, pgau_refine

!$ nthr = omp_get_max_threads()
!$ call omp_set_num_threads(nthr)
nthr_glob = nthr
dir_out = PATH_HERE
call simple_getcwd(cwd)
allocate(CWD_GLOB, source=trim(cwd))

! Testing to make sure the file opens correctly
micname = trim( CWD_GLOB // trim('/simulate_movie_1_intg.mrc'))
!micname = trim( CWD_GLOB // trim('/SLC22A6_intg_1.mrc'))

open(unit=16,file=micname,iostat=iostat,status='old')
if (iostat .eq. 0) then
    print *, 'File ' // trim(micname) // ' exists'
else
    print *, 'An error occured when trying to open the file'
end if
close(16)

! single moldiam pick, picker_utils
print *, 'EXECUTING PICK FROM SIMPLE_PICKER_UTILS'
call putils%new(micname=micname, pcontrast='black', smpd=SMPD, moldiam=MOLDIAM)
print *, 'CREATED NEW PICKER OBJECT'
call putils%exec_picker(boxname_out=boxname_out, nptcls=nptcls, dir_out=dir_out)
print *, 'EXEC_PICKER IS COMPLETE'
call putils%kill

print *, ' '

! single moldiam pick, pickgau
print *, 'EXECUTING PICK FROM SIMPLE_PICKGAU'
call read_mic_raw(micname=micname,smpd=SMPD)
print *, 'READ MICROGRAPH'
call pgau%new_gaupicker(pcontrast='black', smpd_shrink=SMPD_SHRINK, moldiam=MOLDIAM, moldiam_max=MOLDIAM_MAX, offset=OFFSET)
print *, 'CREATED NEW PICKGAU OBJECT'
call pgau_refine%new_gaupicker(pcontrast='black', smpd_shrink=(SMPD_SHRINK / 2), moldiam=MOLDIAM, moldiam_max=MOLDIAM_MAX, offset=1)
print *, 'CREATED REFINED PICKGAU OBJECT'
call pgau%gaupick(pgau_refine)
print *, 'GAUPICK IS COMPLETE'
call pgau%kill
call pgau_refine%kill
print *, ' '

! Check that box files match
! after_center_filter
nlines_pickgau = count_lines('pickgau_after_center_filter.box')
nlines_picker_utils = count_lines('picker_utils_after_center_filter.box')
allocate(coords_pickgau(nlines_pickgau,2))
allocate(coords_picker_utils(nlines_picker_utils,2))
call extract_coords('pickgau_after_center_filter.box', nlines_pickgau, coords_pickgau)
call extract_coords('picker_utils_after_center_filter.box', nlines_picker_utils, coords_picker_utils)
match = compare_coords(coords_pickgau, coords_picker_utils, nlines_pickgau, nlines_picker_utils)
if (match) then
    print *, 'CENTER_FILTER BOX FILES MATCH'
else
    print *, 'CENTER_FILTER BOX FILES DO NOT MATCH'
end if
deallocate(coords_pickgau)
deallocate(coords_picker_utils)

!reset
nlines_pickgau = 0
nlines_picker_utils = 0

! after_distance_filter
nlines_pickgau = count_lines('pickgau_after_distance_filter.box')
nlines_picker_utils = count_lines('picker_utils_after_distance_filter.box')
allocate(coords_pickgau(nlines_pickgau,2))
allocate(coords_picker_utils(nlines_picker_utils,2))
call extract_coords('pickgau_after_distance_filter.box', nlines_pickgau, coords_pickgau)
call extract_coords('picker_utils_after_distance_filter.box', nlines_picker_utils, coords_picker_utils)
match = compare_coords(coords_pickgau, coords_picker_utils, nlines_pickgau, nlines_picker_utils)
if (match) then
    print *, 'DISTANCE_FILTER BOX FILES MATCH'
else
    print *, 'DISTANCE_FILTER BOX FILES DO NOT MATCH'
end if
deallocate(coords_pickgau)
deallocate(coords_picker_utils)

!reset
nlines_pickgau = 0
nlines_picker_utils = 0

! after_remove_outliers
nlines_pickgau = count_lines('pickgau_after_remove_outliers.box')
nlines_picker_utils = count_lines('picker_utils_after_remove_outliers.box')
allocate(coords_pickgau(nlines_pickgau,2))
allocate(coords_picker_utils(nlines_picker_utils,2))
call extract_coords('pickgau_after_remove_outliers.box', nlines_pickgau, coords_pickgau)
call extract_coords('picker_utils_after_remove_outliers.box', nlines_picker_utils, coords_picker_utils)
match = compare_coords(coords_pickgau, coords_picker_utils, nlines_pickgau, nlines_picker_utils)
if (match) then
    print *, 'REMOVE_OUTLIERS BOX FILES MATCH'
else
    print *, 'REMOVE_OUTLIERS BOX FILES DO NOT MATCH'
end if
deallocate(coords_pickgau)
deallocate(coords_picker_utils)

!reset
nlines_pickgau = 0
nlines_picker_utils = 0

! after_refine_upscaled
nlines_pickgau = count_lines('pickgau_after_refine_upscaled.box')
nlines_picker_utils = count_lines('picker_utils_after_refine_positions.box')
allocate(coords_pickgau(nlines_pickgau,2))
allocate(coords_picker_utils(nlines_picker_utils,2))
call extract_coords('pickgau_after_refine_upscaled.box', nlines_pickgau, coords_pickgau)
call extract_coords('picker_utils_after_refine_positions.box', nlines_picker_utils, coords_picker_utils)
match = compare_coords(coords_pickgau, coords_picker_utils, nlines_pickgau, nlines_picker_utils)
!match = compare_coords(coords_pickgau, coords_picker_utils, nlines_pickgau, nlines_pickgau)
if (match) then
    print *, 'REFINED BOX FILES MATCH'
else
    print *, 'REFINED BOX FILES DO NOT MATCH'
end if
call find_closest(  coords_picker_utils, coords_pickgau, nlines_picker_utils , nlines_pickgau)
deallocate(coords_pickgau)
deallocate(coords_picker_utils)

print *, ' '
print *, 'TEST COMPLETE!'

contains

function count_lines(filename) result (cnt)
    character(len=*),     intent(in) :: filename
    integer :: cnt, iostat
    open(unit=20, file=filename, status='old', iostat=iostat)
    if (iostat == 0) then
        !file has been opened successfully
        cnt = 0
        do 
            read(20, *, iostat=iostat)
            if (iostat /= 0) EXIT
            cnt = cnt + 1
        end do
    else
        print *, 'CANNOT OPEN FILE ' // filename
        RETURN
    end if
    close(20)
end function count_lines

subroutine extract_coords(filename, cnt, coords)
    character(len=*),     intent(in) :: filename
    integer,              intent(in) :: cnt
    integer,              intent(out) :: coords(cnt,2)
    character(len=80) :: line
    integer :: iostat, i
    coords = 0
    open(unit=21, file=filename, status='old', iostat=iostat)
    do i = 1, cnt
        read(21, *, iostat=iostat) coords(i,1), coords(i,2)
        if (iostat /= 0) EXIT
    end do
    close(21)
end subroutine extract_coords

function compare_coords(coords_1, coords_2, length_1, length_2) result(coords_same)
    integer, intent(in) :: length_1, length_2
    integer, intent(in) :: coords_1(length_1,2), coords_2(length_2,2)
    logical :: coords_same 
    logical, allocatable :: mask(:)
    integer :: i, j, val_1(2), val_2(2), cnt
    coords_same = .false.
    cnt = 0
    if (size(coords_1) == size(coords_2)) then
        allocate(mask(size(coords_1, dim=1)))
        do i=1,size(coords_1,dim=1)
            val_1 = coords_1(i,:)
            do j=1,size(coords_2,dim=1)
                val_2 = coords_2(j,:)
                if (val_1(1) == val_2(1) .and. val_1(2) == val_2(2)) then 
                    mask(i) = .true.
                    cnt = cnt + 1
                end if
            end do
        end do
        if (count(mask) == size(coords_1, dim=1)) then
            coords_same = .true.
        else
            coords_same = .false.
            print *, size(coords_1, dim=1) - cnt , " DIFFERENCES BETWEEN FILES"
        end if
        deallocate(mask)
    else
        print *, 'COMPARE_COORDS: FILE SIZES DO NOT MATCH'
        print *, 'SIZE COORDS_1 = ', size(coords_1, dim=1), ' BUT SIZE COORDS_2 = ', size(coords_2, dim=1)
    end if
end function compare_coords

subroutine find_closest(coords_1, coords_2, length_1, length_2, filename)
    integer,                    intent(in) :: length_1, length_2
    integer,                    intent(in) :: coords_1(length_1,2), coords_2(length_2,2)
    character(len=*), optional, intent(in) :: filename
    integer :: i, j, closest_coord(2), min_loc(1)
    real, allocatable :: dists(:)
    character(len=1) :: diff
    real :: min_dist

    if (present(filename)) then
        open(unit=22, file=filename)
    else
        open(unit=22, file='combined_coords.csv')
    end if

    do i = 1, size(coords_1,dim=1)
        allocate(dists(size(coords_2, dim=1)))
        do j = 1, size(coords_2,dim=1)
            dists(j) = euclid(real(coords_1(i,:)),real(coords_2(j,:)))
        end do
        min_loc = minloc(dists)
        min_dist = minval(dists)
        closest_coord = coords_2(min_loc(1),:)
        write(22,'(I7,I7,I7,I7,f6.2)') coords_1(i,:), closest_coord(1), closest_coord(2), min_dist
        deallocate(dists)
    end do

    close(22)
end subroutine find_closest

end program simple_test_picker_comp
program simple_test_picker_comp
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_picker_utils
use simple_pickgau, only: pickgau, read_mic_raw, gaupick_multi
use simple_linalg
use simple_image

implicit none
character(len = 200) :: micname 
real, parameter :: SMPD = 1.72, MOLDIAM = 100., SMPD_SHRINK = 4., MOLDIAM_MAX = 200., SMPD_SHRINK2 = 2.
integer, parameter :: OFFSET = 3
integer :: nthr, iostat, idiam, nmoldiams, loc_max(1)
real :: moldiams(11), mic_stats(11,5), moldiams_optimal(2)
character(len=LONGSTRLEN) :: cwd
character(len=:), allocatable :: dir_out
real, allocatable :: max_smds(:,:)
real, allocatable :: max_ksstats(:,:)
real, allocatable :: max_a_peaks(:,:)
type(image) :: ref_imgs(2)
type(pickgau) :: gaup, gaup_refine


!$ nthr = omp_get_max_threads()
!$ call omp_set_num_threads(nthr)
nthr_glob = nthr
dir_out = PATH_HERE
call simple_getcwd(cwd)
allocate(CWD_GLOB, source=trim(cwd))

! Testing to make sure the file opens correctly
micname = trim( CWD_GLOB // trim('/short_linker_intg.mrc'))

open(unit=16,file=micname,iostat=iostat,status='old')
if (iostat .eq. 0) then
    print *, 'File ' // trim(micname) // ' exists'
else
    print *, 'An error occured when trying to open the file'
end if
close(16)


print *, ' '

! ! single moldiam pick, pickgau
! print *, 'EXECUTING PICK FROM SIMPLE_PICKGAU'
! call read_mic_raw(micname=micname,smpd=SMPD)
! print *, 'READ MICROGRAPH'
! call pgau%new_gaupicker(pcontrast='black', smpd_shrink=SMPD_SHRINK, moldiam=MOLDIAM, moldiam_max=MOLDIAM_MAX, offset=OFFSET)
! print *, 'CREATED NEW PICKGAU OBJECT'
! call pgau_refine%new_gaupicker(pcontrast='black', smpd_shrink=(SMPD_SHRINK / 2), moldiam=MOLDIAM, moldiam_max=MOLDIAM_MAX, offset=1)
! print *, 'CREATED REFINED PICKGAU OBJECT'
! call pgau%gaupick(pgau_refine)
! print *, 'GAUPICK IS COMPLETE'
! call pgau%kill
! call pgau_refine%kill
! print *, ' '


! multiple moldiam pick, pickgau
call read_mic_raw(micname=micname,smpd=SMPD)
! print *, 'EXECUTING GAUPICK_MULTI FROM SIMPLE_PICKGAU'
! moldiams = [100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.]
! call gaupick_multi( pcontrast='black', smpd_shrink=SMPD_SHRINK, moldiams=moldiams, offset=OFFSET, mic_stats=mic_stats )
! open(unit=30,file='mic_stats.csv')
! do idiam = 1,11
!     write(30,'(5(f7.3))') mic_stats(idiam,:)
! end do

! ! find local and absolute maxima
! ! local maxima
! nmoldiams = size(moldiams)
! max_smds = find_local_maxima(mic_stats(:,1),mic_stats(:,2),nmoldiams)
! max_ksstats = find_local_maxima(mic_stats(:,1),mic_stats(:,3),nmoldiams)
! max_a_peaks = find_local_maxima(mic_stats(:,1),mic_stats(:,4),nmoldiams)

! print *, 'local max smds occur at ', max_smds(:,1)
! print *, 'local max ksstats occur at ', max_ksstats(:,1)
! print *, 'local max a_peaks occur at ', max_a_peaks(:,1)

! ! absolute maxima
! loc_max = maxloc(mic_stats(:,2))
! print *, 'absolute max smd occurs at ', mic_stats(loc_max(1),1)

! loc_max = maxloc(mic_stats(:,3))
! print *, 'absolute max ksstat occurs at ', mic_stats(loc_max(1),1)

! loc_max = maxloc(mic_stats(:,4))
! print *, 'absolute max a_peak occurs at ', mic_stats(loc_max(1),1)

moldiams_optimal = [100.,200.]
call gaup%new_gaupicker_multi('black', SMPD_SHRINK, moldiams_optimal, offset=OFFSET, roi=.false.)
call gaup%gaupick

contains

function count_lines( filename ) result ( cnt )
    character(len=*),     intent(in) :: filename
    integer :: cnt, iostat
    open(unit=20, file=filename, status='old', iostat=iostat)
    cnt = 0
    if( iostat == 0 )then
        !file has been opened successfully
        do 
            read(20, *, iostat=iostat)
            if( iostat /= 0 ) exit 
            cnt = cnt + 1
        end do
    else
        print *, 'Cannot open file' // filename
        return 
    endif
    close(20)
end function count_lines

subroutine extract_coords( filename, cnt, coords )
    character(len=*),     intent(in)  :: filename
    integer,              intent(in)  :: cnt
    integer,              intent(out) :: coords(cnt,2)
    integer :: iostat, i
    coords = 0
    open(unit=21, file=filename, status='old', iostat=iostat)
    do i = 1, cnt
        read(21, *, iostat=iostat) coords(i,1), coords(i,2)
        if (iostat /= 0) EXIT
    end do
    close(21)
end subroutine extract_coords

function compare_coords( coords_1, coords_2, length_1, length_2 ) result( coords_same )
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

subroutine find_closest( coords_1, coords_2, length_1, length_2, filename )
    integer,                    intent(in) :: length_1, length_2
    integer,                    intent(in) :: coords_1(length_1,2), coords_2(length_2,2)
    character(len=*), optional, intent(in) :: filename
    integer :: i, j, closest_coord(2), min_loc(1)
    real, allocatable :: dists(:)
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

module simple_starproject_utils
include 'simple_lib.f08'
!$ use omp_lib
use simple_sp_project, only: sp_project
use simple_oris, only: oris
use simple_ori,  only: ori
implicit none

public :: LEN_LINE, LEN_FLAG, stk_map, star_flag, star_data, star_file, tilt_info
public :: enable_rlnflag, enable_splflag, enable_splflags, get_rlnflagindex, center_boxes
public :: split_dataline, h_clust
public :: VERBOSE_OUTPUT
private
#include "simple_local_flags.inc"

integer, parameter :: LEN_LINE = 2048
integer, parameter :: LEN_FLAG = 64

logical            :: VERBOSE_OUTPUT =.false.

type stk_map
    character(LEN=2056) :: stkpath
    integer             :: stkmax
end type stk_map

type star_flag
    character(LEN=LEN_FLAG) :: rlnflag
    character(LEN=LEN_FLAG) :: rlnflag2   = ''
    character(LEN=LEN_FLAG) :: splflag
    character(LEN=LEN_FLAG) :: splflag2   = ''
    integer                 :: ind        = 0
    logical                 :: present    = .false.
    logical                 :: string     = .false.
    logical                 :: int        = .false.
    real                    :: mult       = 0
    logical                 :: imagesplit = .false.
    logical                 :: addstk     = .false.
end type star_flag
    
type star_data
    type(star_flag), allocatable :: flags(:)
    integer                      :: flagscount = 0
    integer                      :: datastart  = 0
    integer                      :: dataend    = 0 
end type star_data

type star_file
    character(len=LONGSTRLEN) :: filename
    character(len=LEN_LINE)   :: rootdir
    type(star_data)           :: optics
    type(star_data)           :: stacks
    type(star_data)           :: micrographs
    type(star_data)           :: particles2D
    type(star_data)           :: particles3D
    type(star_data)           :: clusters2D
    integer, allocatable      :: opticsmap(:)
    integer, allocatable      :: stkmap(:,:) ! (stkid : z)
    integer                   :: stkptclcount
    logical                   :: initialised = .false.
end type star_file

type tilt_info
    character(len=LONGSTRLEN) :: basename
    integer                   :: initialtiltgroupid = 1
    integer                   :: finaltiltgroupid   = 1
    real                      :: tiltx = 0.0
    real                      :: tilty = 0.0
    real                      :: smpd  = 0.0
    real                      :: kv    = 0.0
    real                      :: fraca = 0.0
    real                      :: cs    = 0.0
    real                      :: box   = 0.0
end type tilt_info

contains

    subroutine enable_rlnflag(rlnflag, flags, flagindex)
        type(star_flag), intent(inout) :: flags(:)
        character(LEN=LEN_FLAG) :: rlnflag
        integer                 :: flagindex, i
        do i = 1, size(flags)
            if(index(flags(i)%rlnflag, trim(adjustl(rlnflag))) > 0 .AND. len_trim(flags(i)%splflag) > 0) then
                flags(i)%present = .true.
                flags(i)%ind = flagindex
                if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), char(9), "mapping ", flags(i)%rlnflag, " => ", trim(adjustl(flags(i)%splflag))
                exit
            end if
        end do
    end subroutine enable_rlnflag

    subroutine enable_splflag(splflag, flags)
        type(star_flag),  intent(inout) :: flags(:)
        character(LEN=*), intent(in)    :: splflag
        integer :: i
        do i = 1, size(flags)
            if((index(flags(i)%splflag, trim(adjustl(splflag))) > 0 .OR. index(flags(i)%splflag2, trim(adjustl(splflag))) > 0) .AND. len_trim(flags(i)%rlnflag) > 0) then
                if(flags(i)%present .eqv. .false.) then
                    flags(i)%present = .true.
                    if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), char(9), "mapping ", flags(i)%splflag, " => ", trim(adjustl(flags(i)%rlnflag))
                end if
                exit
            end if
        end do
    end subroutine enable_splflag

    subroutine enable_splflags(sporis, flags)
        class(oris),              intent(inout) :: sporis
        type(star_flag),          intent(inout) :: flags(:)
        type(ori)                               :: testori
        character(len=XLONGSTRLEN), allocatable :: keys(:)
        integer                                 :: iori, ikey, testcount
        ! find 1st non state 0 ori
        testcount = 0
        do iori = 1,sporis%get_noris()
            if(sporis%get_state(iori) > 0) then
                call sporis%get_ori(iori, testori)
                keys = testori%get_keys()
                do ikey=1, size(keys)
                    call enable_splflag(trim(adjustl(keys(ikey))), flags)
                end do
                testcount = testcount + 1
                if(testcount == 10) then ! Do 10 tests to counter xpos or ypos being 0
                    exit
                end if
            end if
        end do
    end subroutine enable_splflags

    subroutine get_rlnflagindex(rlnflag, flags, flagindex)
        type(star_flag),  intent(in)  :: flags(:)
        character(LEN=*), intent(in)  :: rlnflag
        integer,          intent(out) :: flagindex
        integer :: i
        flagindex = 0
        do i = 1, size(flags)
            if(index(flags(i)%rlnflag, trim(adjustl(rlnflag))) > 0) then
                flagindex = flags(i)%ind
                exit
            end if
        end do
    end subroutine get_rlnflagindex

    subroutine split_dataline(line, splitline)
        character(len=LEN_LINE), intent(in)    :: line
        character(len=LEN_LINE), intent(inout) :: splitline(:)
        integer :: iend, istart, flagid
        flagid = 1
        iend   = 1
        istart = 1
        do while (flagid <= size(splitline))
            do while (line(istart:istart + 1) .eq. " ")         
                istart = istart + 1           
            end do          
            iend = index(line(istart + 1:), ' ') + istart
            splitline(flagid) = trim(adjustl(line(istart:iend)))
            istart = iend + 1
            flagid = flagid + 1
        end do
    end subroutine split_dataline
    
    subroutine center_boxes(spproj, sporis)
        class(sp_project),        intent(inout) :: spproj
        class(oris),              intent(inout) :: sporis
        integer                                 :: iori
        real                                    :: stkind, stkbox, xpos, ypos
        do iori=1, sporis%get_noris()
            if(sporis%get_state(iori) > 0) then
                stkind = sporis%get(iori, "stkind")
                stkbox = spproj%os_stk%get(int(stkind), "box")
                xpos = sporis%get(iori, "xpos") + (stkbox / 2)
                ypos = sporis%get(iori, "ypos") + (stkbox / 2)
                call sporis%set(iori, "xpos", xpos)
                call sporis%set(iori, "ypos", ypos)
            end if
        end do
    end subroutine center_boxes
    
    ! distance threshold based yerarchical clustering
    ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
    subroutine h_clust(data_in, thresh,labels, centroids, populations)
        real,                 intent(in)  :: data_in(:,:)   ! input data, point coords
        real,                 intent(in)  :: thresh         ! threshold for class merging
        integer,              intent(out) :: labels(:)      ! labels of the elements in vec
        real,    allocatable, intent(out) :: centroids(:,:) ! centroids of the classes
        integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
        real,    allocatable :: mat(:,:)                    ! pariwise distances matrix
        logical, allocatable :: mask(:), outliers(:)
        integer :: N, i, j, cnt, ncls, filnum, io_stat
        integer :: index(1), loc1(1), loc2(1)
        real    :: d
        if( size(data_in, dim = 2) .ne. 2 )then
            THROW_HARD('Input data should be two dimensional!; h_clust')
        endif
        N = size(data_in, dim = 1) ! number of data points
        ! 1) calc all the couples of distances, using euclid dist
        allocate(mat(N,N), source = 0.)
        do i = 1, N-1
            do j = i+1, N
                mat(i,j) = sqrt((data_in(i,1)-data_in(j,1))**2 + (data_in(i,2)-data_in(j,2))**2) ! pariwise euclidean distance
                mat(j,i) = mat(i,j)
            enddo
        enddo
        ! 2) Generate binary clusters
        allocate(mask(N),     source = .true. )
        allocate(outliers(N), source = .false.)
        ncls = 0
        do i = 1, N
            if( mask(i) )then ! if it's not already been clustered
                mask(i) = .false.
                ! find the index of the couple
                d = minval(mat(i,:), mask)
                index(:) = minloc(mat(i,:), mask)
                ncls = ncls + 1
                ! assign labels
                labels(i) = ncls
                if(d <= thresh) then ! if it's not an outlier (it has a couple)
                    labels(index(1)) = labels(i)
                    mask(index(1)) = .false. ! index(1) has already been clustered
                else
                    outliers(i) = .true.
                endif
            endif
        enddo
        ! 3) Calculate centroids
        allocate(centroids(ncls,2), source = 0.)
        mask = .true. ! reset
        do i = 1, ncls
            ! find the first member of the class
            loc1(:) = minloc(abs(labels-i))
            if(.not. outliers(loc1(1))) then
                mask(loc1(1)) = .false.
                ! find the second member of the class
                loc2(:) = minloc(abs(labels-i), mask)
                mask(loc2(1)) = .false.
                centroids(i,1) = (data_in(loc1(1),1) + data_in(loc2(1),1))/2.
                centroids(i,2) = (data_in(loc1(1),2) + data_in(loc2(1),2))/2.
            else ! the class has just one element
                loc1(:) = minloc(abs(labels-i))
                centroids(i,1) = data_in(loc1(1),1)
                centroids(i,2) = data_in(loc1(1),2)
                mask(loc1(1)) = .false.
            endif
        enddo
        mask  = .true. ! reset
        ! 4) merge classes
        do i = 1, ncls-1
            do j = i+1, ncls
                if(sqrt((centroids(i,1)-centroids(j,1))**2+(centroids(i,2)-centroids(j,2))**2) <= thresh) then ! merge classes
                    ! change label to class j
                    where(labels == j) labels = i
                endif
            enddo
        enddo
        ! 5) Reoder labels
        cnt = 0
        do i = 1, ncls
            if( any(labels== i) )then !there is a class labelled i
                cnt = cnt + 1
                where(labels == i) labels = cnt
            endif
        enddo
        ! 6) recalculate centroids
        deallocate(centroids)
        ncls = maxval(labels) ! the nr of classes is maxval(labels)
        allocate(centroids(ncls,2), source = 0.)
        allocate(populations(ncls), source = 0 )
        mask = .true. ! reset
        do i = 1, ncls
            populations(i) = count(labels == i) ! counts the nr of elements in the class
            ! find the all the cnt member of the class and update the centroids
            do j = 1, populations(i)
                loc1(:) = minloc(abs(labels-i), mask)
                mask(loc1(1)) = .false. ! do not consider this member of the class anymore
                centroids(i,1) = centroids(i,1)+ data_in(loc1(1),1)
                centroids(i,2) = centroids(i,2)+ data_in(loc1(1),2)
            enddo
            centroids(i,:) = centroids(i,:)/real(populations(i))
        enddo
        ! output generation
        do i = 1, ncls
            if( VERBOSE_OUTPUT ) write(logfhandle, *) char(9), char(9), char(9), 'Class ', i, 'Population ', populations(i)
        enddo
    end subroutine h_clust

end module simple_starproject_utils

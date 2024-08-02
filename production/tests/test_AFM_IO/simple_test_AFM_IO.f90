program AFM_File_IO
use iso_c_binding
include 'simple_lib.f08'
use simple_AFM_image
! use simple_AFM_image 
!AFM IO
type :: AFM_image
    type(image), allocatable :: img_array(:)
    character(len = 20), allocatable :: img_names(:)
end type AFM_image 

type    :: bin_header5
    INTEGER(KIND=2) :: version 
    INTEGER(KIND=2) :: checksum 
    INTEGER(KIND=4) :: wfmSize 
    INTEGER(KIND=4) :: formulaSize 
    INTEGER(KIND=4) :: noteSize 
    INTEGER(KIND=4) :: dataEUnitsSize 
    INTEGER(KIND=4) :: dimEUnitsSize(4) 
    INTEGER(KIND=4) :: dimLabelsSize(4)
    INTEGER(KIND=4) :: sIndicesSize 
    INTEGER(KIND=4) :: optionsSize1 
    INTEGER(KIND=4) :: optionsSize2
end type

type    :: wave_header5
    INTEGER(KIND=4)    :: Structure_Padding
    INTEGER(C_INT32_T) :: creationDate
    INTEGER(C_INT32_T) :: modDate
    INTEGER(KIND=4)    :: npnts 
    INTEGER(KIND=2)    :: type 
    INTEGER(KIND=2)    :: dlock 
    CHARACTER          :: whpad(6)
    INTEGER(KIND=2)    :: whVersion 
    CHARACTER          :: bname(32)
    INTEGER(KIND=4)    :: whpad2
    INTEGER(KIND=4)    :: DataFolder_Padding
    INTEGER(KIND=4)    :: nDim(4)
    REAL(KIND=8)       :: sfA(4)
    REAL(KIND=8)       :: sfB(4)
    CHARACTER          :: dataUnits(4)
    CHARACTER          :: dimUnits(4,4)
    INTEGER(KIND=2)    :: fsvalid
    INTEGER(KIND=2)    :: whpad3
    REAL(KIND=8)       :: topFullScale, botFullScale
    INTEGER(KIND = 4)  :: dataEUnits, dimEUnits(4), dimLabels(4), waveNoteH
    INTEGER(KIND=4)    :: whUnused(16)
    INTEGER(KIND=2)    :: aModified
    INTEGER(KIND=2)    :: wModified 
    INTEGER(KIND=2)    :: sModified  
    CHARACTER          :: useBits
    CHARACTER          :: kindBits
    INTEGER(KIND=4)    :: formula_pointer
    INTEGER(KIND=4)    :: depID
    INTEGER(KIND=2)    :: whpad4, srcFldr
    CHARACTER(KIND=4)  :: fileName
    INTEGER(KIND=4)    :: sIndices 
end type

! interface
!     include 'simple_lib.f08'
!     use simple_image
!     type    :: AFM
!         type(image_stack)      :: AFM_stack
!         character(len = 30), pointer   :: names(:) => null()
!     contains 
!         procedure   :: set_names
!     end type AFM
! end interface 


call read_ibw('/Users/atifao/Downloads/IBW/Cob_450016.ibw')
! Can add integer/complex support later. 
! open(newunit = data, file = fn_in, status = 'old', access='stream')
! read(in, pos = 385) first_entry 
! delete file at some point
contains 

    subroutine read_ibw(fn_in) 
        character(len=*),           intent(in)  :: fn_in 
        integer         :: in, check, real_type
        integer                          :: real_type1, data, total_bytes, bytes_read, iter_ind, prop_ind, img_ind, i
        real(kind = 4), allocatable          :: Rank3_Data_4byte(:, :, :, :)
        real(kind = 8), allocatable          :: Rank3_Data_8byte(:, :, :, :)
        character(:), allocatable                            :: channel_info
        character(len = 10)    :: iteration(2), properties(4)
        type(AFM_image)       :: AFM
        type(bin_header5)   :: binheader
        type(wave_header5)  :: waveheader
        iteration = [character(len = 10) :: 'Trace', 'Retrace']
        properties = [character(len = 10) :: 'Height', 'Amplitude', 'Phase', 'ZSensor' ]
        if(index(fn_in, '.ibw') == 0) then 
            print *, 'Error: only .ibw files are supported'
        end if
        open(newunit = check, file = fn_in, status = 'old', access='stream')
        read(check) binheader%version
        if(binheader%version > 5) then
            open(newunit = in, file = fn_in, status = 'old', access='stream', convert='swap')
        else
            open(newunit = in, file = fn_in, status = 'old', access='stream')
        end if
        read(in)binheader, waveheader
        if(binheader%version /= 5) then
            print *, "Error: only version 5 files are supported"
            stop 
        end if 
        allocate(Rank3_Data_4byte(waveheader%nDim(1) ,waveheader%nDim(2), 1, waveheader%nDim(3)))
        ! select case (waveheader1%type)
        ! case (2) 
        !     real_type1 = 4
        ! case (4)
        !     real_type1 = 8
        !     Rank3_Data_8byte = dble(Rank3_Data_4byte) 
        ! case default
        !     print *, "Error: only float data is currently supported"
        !     stop
        ! end select
        open(newunit = data, file = fn_in, status = 'old', access='stream')
        read(data, pos = 385) Rank3_Data_4byte
        inquire(data, pos = bytes_read)
        allocate(character(binheader%dimLabelsSize(3)) :: channel_info)
        allocate(AFM%img_names(size(iteration) * size(properties)))
        read(data, pos = binheader%noteSize + bytes_read) channel_info
        print *, channel_info 
        ! getting image names 
        do prop_ind = 1, size(properties)
            do iter_ind = 1, size(iteration)
                if( index(channel_info, trim(properties(prop_ind)) // trim(iteration(iter_ind)) ) /= 0 .AND. iter_ind == 2 ) then 
                    AFM%img_names(iter_ind*prop_ind) = trim(properties(prop_ind)) // trim(iteration(iter_ind))
                else if( index(channel_info, trim(properties(prop_ind)) // trim(iteration(iter_ind)) ) /= 0 .AND. iter_ind == 1) then
                    AFM%img_names(2*prop_ind - iter_ind) = trim(properties(prop_ind)) // trim(iteration(iter_ind)) 
                end if 
            end do 
        end do     
        print *, AFM%img_names 

        allocate(AFM%img_array(waveheader%nDim(3)))
        do img_ind = 1, waveheader%nDim(3)
            call AFM%img_array(img_ind)%new([waveheader%nDim(1), waveheader%nDim(2), 1], real(waveheader%sfA(1)) * 10**10)
        end do 
        ! print *,  img_array(1)%get_smpd()
        print *, AFM%img_array(1)%get_ldim()
        print *, waveheader%nDim(3)
        print *, size(Rank3_Data_4byte, dim=1), size(Rank3_Data_4byte, dim=2), size(Rank3_Data_4byte, dim=3)
        call AFM%img_array(1)%set_rmat(Rank3_Data_4byte(:, :, 1, :), .false.)
        do img_ind = 1, waveheader%nDim(3)
            call AFM%img_array(img_ind)%set_rmat(Rank3_Data_4byte(:, :, img_ind, :), .false.)
        end do 
        ! print *, img_array(1)%get_rmat()
        ! might need to transpose or rotate
        ! call img_array(3)%vis()
        call AFM%img_array(findloc(index(AFM%img_names, 'HeightTrace'),1, dim = 1))%vis()
    end subroutine read_ibw

    ! subroutine get_property(AFM, property)
    !     type(AFM_image), intent(in)      :: AFM
    !     character(len = *), intent(in)                     :: property
    !     AFM%img_array(findloc(index(AFM%img_names, property), 1, dim = 1))
    ! end subroutine get_property 

end program AFM_File_IO
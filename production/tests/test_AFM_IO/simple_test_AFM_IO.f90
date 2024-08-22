program AFM_File_IO
use iso_c_binding
include 'simple_lib.f08'
use simple_image
use simple_segmentation
use simple_hash
use simple_ced_filter,         only: ced_filter_2D
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
type(AFM_image) ::  Cob16
type(image)     ::  HeightTrace, Height_Avg
type(image)     ::  HeightRetrace, Retrace_Copy
type(image)     ::  Sum
integer         :: i, j, fu
real, allocatable   :: rmat_test(:,:,:)
real            :: h_shift(2), shvec(3), a, iterx, itery 
real,    pointer     :: J11_rmat(:,:,:)=>null()
! call J11%new([1024, 1024, 1], 5.0)
! call J11%get_rmat_ptr(J11_rmat)
! allocate(J11_rmat(1024, 1024, 1))
! print *, J11_rmat
! character(len=*), parameter :: OUT_FILE = '/Users/atifao/Downloads/IBW/16_im.txt' ! Output file.
! character(len=*), parameter :: PLT_FILE = 'plot.plt' ! Gnuplot file.
! real       :: test_mat(3,3), test_sing(3), test_vec(3,3)
! test_mat = reshape((/ 2, 1, 3, 3, 7, 2, 6, 9, -1 /), shape(test_mat))
! call svdcmp_sp(test_mat, test_sing, test_vec)

! print *, test_sing
! add padding to non square images 
call read_ibw('/Users/atifao/Downloads/IBW/Cob_450016.ibw', Cob16)
call zero_padding(Cob16)

call get_AFM(Cob16, 'HeightTrace', HeightTrace)
HeightRetrace = Cob16%img_array(findloc(index(Cob16%img_names, 'HeightRetrace'),1, dim = 1))
call align_avg(Cob16)
call Height_Avg%copy(HeightTrace)
call Retrace_Copy%copy(HeightRetrace)
! print *, HeightTrace%mean()
! call ced_filter_2D(HeightRetrace, 0.7)
! call HeightTrace%vis()
! print *, HeightTrace%mean()

call HeightTrace%fft()
call HeightRetrace%fft()

! call HeightTrace%fcorr_shift(HeightRetrace, 20., h_shift, .true.) 
! have to store height trace since it is destroyed on fcorr_shift output
shvec = [h_shift(1), h_shift(2), 0.]
! print *, shvec 
call HeightRetrace%fft()
call HeightRetrace%shift(shvec)
call HeightRetrace%ifft()
! call HeightTrace%ifft()

! call HeightRetrace%vis()
! call Trace_Copy%vis()
Height_Avg = HeightRetrace * 0.5 + Height_Avg* 0.5
! call Height_Avg%vis()
! call Trace_Copy%vis()
! call softmin_avg(Trace_Copy, HeightRetrace)
call otsu_img(Height_Avg, positive = .false.)
! call Height_Avg%vis()

call otsu_img(Retrace_Copy, positive = .false.)
! call Retrace_Copy%vis()
! do i = 1, Trace_Copy%get_ldim(1)
!     do j = 1, Trace_Copy%ldim(2)
!         -beta**(-1)*(exp(-beta*x) + exp(-beta*y)) + log(2)
!     end do 
! end do 

! shvec(3) = 1
! do iterx = 9.38, 10.09, 0.01
!     do itery = -1.09, -0.61, 0.01
!         shvec(1) = iterx
!         shvec(2) = itery
!         a = HeightTrace%corr_shifted(HeightRetrace, shvec)
!         print *, iterx, itery, a
!     end do 
! end do 

!need to adjust make function sub-pixel accurate. 


! print *, h_shift

! normalize retrace to max value of normalized retrace

! Subtract = HeightTrace + PhaseTrace
! call Subtract%vis()
! rmat_test = HeightTrace%get_rmat()

! open (action = 'write', file = OUT_FILE, unit = fu)
!     write (fu, *) rmat_test
! close (fu)
! call HeightTrace%vis()
! call HeightRetrace%vis()
! Can add integer/complex support later. 
! open(newunit = data, file = fn_in, status = 'old', access='stream')
! read(in, pos = 385) first_entry 
! delete file at some point
!normalize and then align =
contains 
    subroutine read_ibw(fn_in, AFM)
        character(len=*),           intent(in)  :: fn_in 
        integer         :: in, check, real_type
        integer                          :: real_type1, data, total_bytes, bytes_read, iter_ind, prop_ind, img_ind, i
        real(kind = 4), allocatable          :: Rank3_Data_4byte(:, :, :, :)
        character(:), allocatable                            :: channel_info
        character(len = 10)    :: iteration(2), properties(4)
        type(AFM_image), intent(out)       :: AFM
        type(bin_header5)   :: binheader
        type(wave_header5)  :: waveheader
        iteration = [character(len = 10) :: 'Trace', 'Retrace']
        properties = [character(len = 10) :: 'Height', 'Amplitude', 'Phase', 'ZSensor' ]
        if(index(fn_in, '.ibw') == 0) then 
            print *, 'Error: only .ibw files are supported'
            stop 
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
        open(newunit = data, file = fn_in, status = 'old', access='stream')
        read(data, pos = 385) Rank3_Data_4byte
        inquire(data, pos = bytes_read)
        allocate(character(binheader%dimLabelsSize(3)) :: channel_info)
        allocate(AFM%img_names(size(iteration)*size(properties)))
        read(data, pos = binheader%noteSize + bytes_read) channel_info
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
        allocate(AFM%img_array(waveheader%nDim(3)))
        do img_ind = 1, waveheader%nDim(3)
            call AFM%img_array(img_ind)%new([waveheader%nDim(1), waveheader%nDim(2), 1], real(waveheader%sfA(1)) * 10**10)
            call AFM%img_array(img_ind)%set_rmat(Rank3_Data_4byte(:, :, img_ind, :), .false.)
        end do 
        deallocate(Rank3_Data_4byte)
    end subroutine read_ibw

    ! Zero Padding to make it square (maybe set to intensity average)
    subroutine zero_padding(AFM_pad)
        type(AFM_image), intent(inout)  :: AFM_pad
        integer        :: img_ind
        integer :: dim(3)
        do img_ind = 1, size(AFM_pad%img_array)
            dim = AFM_pad%img_array(img_ind)%get_ldim()
            if( dim(1) /= dim(2) ) then 
                call AFM_pad%img_array(img_ind)%pad_inplace([maxval(dim), maxval(dim), 1])
            end if 
        end do 
    end subroutine zero_padding 
  
    subroutine destriping(AFM)
        type(AFM_image), intent(inout)  :: AFM 
    end subroutine destriping

    subroutine matrix_log(rmat, rmat_eigvl, rmat_eigvc)
        real, allocatable   :: rmat(:, :), rmat_eigvl(:, :), rmat_eigvc(:, :)
    end subroutine matrix_log

    ! alignment subroutine for all measurements.
    subroutine align_avg(Align_AFM)
        type(AFM_image), intent(inout)     :: Align_AFM
        real, allocatable                  :: shifts(:, :)
        integer                            :: num_avg, avg_ind
        call Align_AFM%img_array(1)%vis()
        num_avg = size(Align_AFM%img_array) / 2
        allocate(shifts(2, num_avg))
        if(modulo(size(shifts), 2) /= 0) then 
            print *, 'One or more traces are not present'
        end if 
        do avg_ind = 1, size(Align_AFM%img_array)
            call Align_AFM%img_array(avg_ind)%fft()
        end do
        do avg_ind = 2, size(Align_AFM%img_array), 2
            call Align_AFM%img_array(avg_ind - 1)%fcorr_shift(Align_AFM%img_array(avg_ind), 20., shifts(:, avg_ind / 2), .true.)
        end do 
        do avg_ind = 1, num_avg 
            call Align_AFM%img_array(2*avg_ind)%fft()
            call Align_AFM%img_array(2*avg_ind)%shift([shifts(1, avg_ind), shifts(2, avg_ind), 0.])
            call Align_AFM%img_array(2*avg_ind)%ifft()
        end do
        do avg_ind = 2, size(Align_AFM%img_array), 2
           Align_AFM%img_array(avg_ind - 1) =  Align_AFM%img_array(avg_ind) * 0.5 + Align_AFM%img_array(avg_ind - 1) * 0.5
           call Align_AFM%img_array(avg_ind - 1)%vis()
        end do 
    end subroutine align_avg

    !change this to function... 
    subroutine get_AFM(AFM_Hash, key, image_at_key)
        type(AFM_image), intent(in) :: AFM_Hash
        character(*), intent(in)    :: key 
        type(image), intent(out)  :: image_at_key
        image_at_key = AFM_Hash%img_array(findloc(index(AFM_Hash%img_names, key),1, dim = 1))
    end subroutine get_AFM  

    function fget_AFM(AFM_Hash, key) result(image_at_key)
        type(AFM_image), intent(in) :: AFM_Hash
        character(*), intent(in)    :: key 
        type(image)     :: image_at_key
        image_at_key = AFM_Hash%img_array(findloc(index(AFM_Hash%img_names, key),1, dim = 1))
    end function fget_AFM  

    subroutine softmin_avg(trace, retrace, beta)
        type(image), intent(inout) :: trace 
        type(image), intent(in) :: retrace
        real, optional, intent(in) :: beta 
        real                       :: beta_d, int_t, int_r, int_avg
        integer                    :: dim(3), xdim, ydim
        beta_d = 0.0000000002
        if(present(beta)) beta_d = beta 
        dim = trace%get_ldim()
        do xdim = 1, dim(1)
            do ydim = 1, dim(2)
                int_t = trace%get_rmat_at(xdim, ydim, 1)
                int_r = retrace%get_rmat_at(xdim, ydim, 1)
                int_avg = beta_d**(-1.) * log(2.) * exp(abs(int_t - int_r)*(-beta_d)**2.)
                print *, int_t, int_r, int_avg 
                call trace%set_rmat_at(xdim, ydim, 1, int_avg)
            end do 
        end do
    end subroutine softmin_avg
    ! print *,  img_array(1)%get_smpd()
    ! print *, AFM%img_array(1)%get_ldim()
    ! print *, waveheader%nDim(3)
    ! print *, size(Rank3_Data_4byte, dim=1), size(Rank3_Data_4byte, dim=2), size(Rank3_Data_4byte, dim=3)
    ! call AFM%img_array(1)%set_rmat(Rank3_Data_4byte(:, :, 1, :), .false.)
    ! do img_ind = 1, 8
    !     print *, AFM%img_array(1)%real_corr(AFM%img_array(img_ind))
    ! end do     
    ! print *, AFM%img_array(1)%real_corr(AFM%img_array(1))
    ! print *, img_array(1)%get_rmat()
    ! call img_array(3)%vis()
    ! call AFM%img_array(findloc(index(AFM%img_names, 'HeightTrace'),1, dim = 1))%vis()
    ! call AFM%img_array(findloc(index(AFM%img_names, 'HeightRetrace'),1, dim = 1))%vis()
    ! HeightTrace = AFM%img_array(findloc(index(AFM%img_names, 'HeightTrace'),1, dim = 1))
    ! call otsu_img(HeightTrace, positive = .false.)
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
    ! real(kind = 8), allocatable          :: Rank3_Data_8byte(:, :, :, :)
end program AFM_File_IO
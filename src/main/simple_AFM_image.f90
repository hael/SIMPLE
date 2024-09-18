module simple_AFM_image 
use iso_c_binding
include 'simple_lib.f08'
use simple_image
use simple_ced_filter,         only: ced_filter_2D
use simple_pickseg
use simple_parameters
use simple_segmentation
use simple_binimage
use simple_neighs
implicit none 

logical     :: L_DEBUG = .true.
type :: AFM_image
    type(image), allocatable :: img_array(:)
    character(len = 50), allocatable :: img_names(:)
contains
    procedure   :: read_ibw
    procedure   :: pick_valid
    procedure   :: get_AFM
    procedure   :: zero_padding
    procedure   :: align_avg
end type AFM_image 
contains
    subroutine read_ibw(AFM, fn_in)
        class(AFM_image), intent(out)       :: AFM
        character(len=*),      intent(in)  :: fn_in 
        integer         :: in, check, real_type
        integer                          :: real_type1, data, total_bytes, bytes_read, iter_ind, prop_ind, img_ind, i
        real(kind = 4), allocatable          :: Rank3_Data_4byte(:, :, :, :)
        character(:), allocatable                            :: channel_info
        character(len = 10)    :: iteration(2), properties(4)
    
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
            call AFM%img_array(img_ind)%new([waveheader%nDim(1), waveheader%nDim(2), 1], real(waveheader%sfA(1)) * 10.**10.)
            call AFM%img_array(img_ind)%set_rmat(Rank3_Data_4byte(:, :, img_ind, :), .false.)
            call AFM%img_array(img_ind)%norm_minmax()
            call AFM%img_array(img_ind)%bin_inv
        end do 
        deallocate(Rank3_Data_4byte)
        
    end subroutine read_ibw

    subroutine pick_valid(AFM_in)
        class(AFM_image), intent(inout)    :: AFM_in 
        type(image)                        :: HeightTrace, HeightRetrace, AvgHeight
        CHARACTER(len=255)                 :: cwd
        type(pickseg)                      :: avg_p, trace_p, retrace_p
        type(image)                        :: avg_slim, trace_slim, retrace_slim 
        integer                            :: ldim_box(3), box_iter, search_iter, neighbor_iter, box_count
        real                               :: smpd_box, neighbor_corr(8), coord_corr(3,8), max_corr(20)
        integer, allocatable               :: pickpos(:, :), val_center_r(:, :), val_center_t(:, :)
        integer                            :: coord_test(2), center_x(20), center_y(20)
        logical                            :: outs
        integer                            :: neighbor(3, 8), nsiz, center(3)
        real                               :: corr_r, corr_t, val_score_t, val_score_r 
        real, allocatable                  :: corr_final_t(:), corr_final_r(:)
        
        call get_AFM(AFM_in, 'AvgHeight', AvgHeight)
        call get_AFM(AFM_in, 'HeightTrace', HeightTrace)
        call get_AFM(AFM_in, 'HeightRetrace', HeightRetrace)

        call AvgHeight%norm_minmax()
        call HeightTrace%norm_minmax()
        call HeightRetrace%norm_minmax()
        
        call getcwd(cwd)
        call AvgHeight%write(trim(cwd) // 'avg.mrc')
        call avg_p%pick(trim(cwd) // 'avg.mrc')
        call HeightTrace%write(trim(cwd) // 'trace.mrc')
        call trace_p%pick(trim(cwd) // 'trace.mrc')
        call HeightRetrace%write(trim(cwd) // 'retrace.mrc')
        call retrace_p%pick(trim(cwd) // 'retrace.mrc')

        ldim_box = [avg_p%box_raw, avg_p%box_raw, 1]
        smpd_box =  AvgHeight%get_smpd()
        call avg_slim%new(ldim_box, smpd_box )
        call trace_slim%new(ldim_box, smpd_box)
        call retrace_slim%new(ldim_box, smpd_box)

        box_count = 0 
        ! allocate(val_centers(3, avg_p%get_nboxes()))
        allocate(val_center_r(3, avg_p%get_nboxes()))
        allocate(val_center_t(3, avg_p%get_nboxes()))
        val_center_t = 0
        val_center_r = 0
        call avg_p%get_positions(pickpos) 
        allocate(corr_final_r(avg_p%get_nboxes()))
        allocate(corr_final_t(avg_p%get_nboxes()))
        do box_iter = 1, avg_p%get_nboxes() 
            
            coord_test = pickpos(box_iter, :)
            call AvgHeight%window_slim(coord_test, avg_p%box_raw, avg_slim, outs)
            call HeightTrace%window_slim(coord_test, avg_p%box_raw, trace_slim, outs)
            call HeightRetrace%window_slim(coord_test, avg_p%box_raw, retrace_slim, outs)

            corr_r =  avg_slim%real_corr(retrace_slim)
            corr_t = avg_slim%real_corr(trace_slim)
            center = [coord_test(1), coord_test(2), 1]
            center_x(1) = coord_test(1)
            center_y(1) = coord_test(2)
            ! print *, 'box: ', box_iter 
            if( corr_r < 0.1 .or. corr_t < 0.1) then 
                ! print *, 'box is outside the image'
                cycle
            end if 
            box_count = box_count + 1 

            call nn_val(HeightTrace, corr_final_t, val_center_t)                                                 
            call nn_val(HeightRetrace, corr_final_r, val_center_r)
            
        end do 
        
        if(L_DEBUG) then 
            if(sum(corr_final_t) /  box_count > sum(corr_final_r) / box_count) then 
                print *, 'retrace is more noisy than trace'
            else
                print *, 'trace is more noisy than retrace'
            end if 
            print *, sum(corr_final_t) /  box_count, sum(corr_final_r) / box_count
            print *, val_center_t, val_center_r
            do box_iter = 73, 75
                    call AvgHeight%window_slim(pickpos(box_iter, :), avg_p%box_raw, avg_slim, outs)
                    call avg_slim%vis()
                    call HeightTrace%window_slim(val_center_t(:2, box_iter), avg_p%box_raw, trace_slim, outs)
                    call trace_slim%vis()
                    call HeightRetrace%window_slim(val_center_r(:2, box_iter), avg_p%box_raw, retrace_slim, outs)
                    call retrace_slim%vis()
            end do
        end if 
        
        contains 
            subroutine nn_val(im_iter, corr_final, val_centers)
                type(image), intent(in)     :: im_iter
                real, intent(out) :: corr_final(avg_p%get_nboxes())
                integer, intent(out)    :: val_centers(3, avg_p%get_nboxes())
                ! output validated centers
                do search_iter = 2, 20
                    call neigh_8_1(AvgHeight%get_ldim(), center, neighbor, nsiz)
                    do neighbor_iter = 1, nsiz
                        call im_iter%window_slim([neighbor(1, neighbor_iter), neighbor(2, neighbor_iter)], avg_p%box_raw, trace_slim, outs)
                        coord_corr(:, neighbor_iter) = [real(neighbor(1, neighbor_iter)), real(neighbor(2, neighbor_iter)), avg_slim%real_corr(trace_slim)]
                        neighbor_corr(neighbor_iter) = avg_slim%real_corr(trace_slim)
                    end do 
                    center = [ int(coord_corr(1, maxloc(neighbor_corr))), int(coord_corr(2, maxloc(neighbor_corr))), 1 ] 
                    ! print *, center 
                    max_corr(search_iter) = maxval(neighbor_corr)
                    if(search_iter > 2 .and. max_corr(search_iter - 1) < max_corr(search_iter)) then 
                        ! print *, center, max_corr(search_iter)
                        val_centers(:, box_iter) = center 
                        corr_final(box_iter) = max_corr(search_iter)
                        exit
                    end if 
                end do
            end subroutine 
    end subroutine pick_valid  

    subroutine get_AFM(AFM_Hash, key, image_at_key)
        class(AFM_image), intent(in) :: AFM_Hash
        character(*), intent(in)    :: key 
        type(image), intent(out)  :: image_at_key
        image_at_key = AFM_Hash%img_array(findloc(index(AFM_Hash%img_names, key),1, dim = 1))
    end subroutine get_AFM

    subroutine zero_padding(AFM_pad)
        class(AFM_image), intent(inout)  :: AFM_pad
        integer        :: img_ind
        integer :: dim(3)
        do img_ind = 1, size(AFM_pad%img_array)
            dim = AFM_pad%img_array(img_ind)%get_ldim()
            if( dim(1) /= dim(2) ) then 
                call AFM_pad%img_array(img_ind)%pad_inplace([maxval(dim), maxval(dim), 1])
            end if 
        end do 
    end subroutine zero_padding
    
    subroutine align_avg(AFM_in, Align_AFM)
        class(AFM_image), intent(in)     :: AFM_in 
        class(AFM_image), intent(out)    :: Align_AFM
        real, allocatable                  :: shifts(:, :)
        integer                            :: num_avg, avg_ind, prop_ind, new_size, count, num_mic, tr_ind, retr_ind, i 
        character(len = 50)                :: new_name 
        num_mic = size(AFM_in%img_array)
        new_size = int(num_mic*1.5)
        allocate(Align_AFM%img_array(new_size))
        allocate(Align_AFM%img_names(new_size))
        num_avg = num_mic/2
        allocate(shifts(2, num_avg))
        if(modulo(size(shifts), 2) /= 0) then 
            print *, 'One or more traces are not present'
        end if 
        do i = 1, num_mic
            Align_AFM%img_array(i) = AFM_in%img_array(i)
            Align_AFM%img_names(i) = AFM_in%img_names(i)
        end do 
        count = 0
        do i = num_mic + 1, new_size
            Align_AFM%img_array(i) = AFM_in%img_array(i - num_mic + count)
            count = count + 1
        end do 
    
        do avg_ind = 1, size(Align_AFM%img_array)
            call Align_AFM%img_array(avg_ind)%fft()
        end do

        count = 0
        do avg_ind = num_mic + 1, new_size
            count = count + 1
            call Align_AFM%img_array(avg_ind)%fcorr_shift(Align_AFM%img_array(avg_ind - num_mic + count), 20., shifts(:, count), .true.)
        end do 
        
        do avg_ind = 1, num_avg
            call Align_AFM%img_array(2*avg_ind)%fft()
            call Align_AFM%img_array(2*avg_ind)%shift([shifts(1, avg_ind), shifts(2, avg_ind), 0.])
            call Align_AFM%img_array(2*avg_ind)%ifft()
        end do
        
        do avg_ind = 1, new_size 
            if(Align_AFM%img_array(avg_ind)%is_ft()) then
                call Align_AFM%img_array(avg_ind)%ifft()
            end if    
        end do 
        count = 0 
        do avg_ind = num_mic + 1, new_size 
           tr_ind = avg_ind - num_mic + count
           retr_ind = avg_ind - num_mic + 1 + count
           Align_AFM%img_array(avg_ind) =  Align_AFM%img_array(tr_ind) * 0.5 + Align_AFM%img_array(retr_ind) * 0.5
           new_name = AFM_in%img_names(tr_ind)
           Align_AFM%img_names(avg_ind) = 'Avg' // new_name(1:len_trim(new_name)  - 5)
           count = count + 1
        end do 
    end subroutine align_avg

end module simple_AFM_image 

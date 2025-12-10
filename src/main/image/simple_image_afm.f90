module simple_image_afm 
use iso_c_binding
include 'simple_lib.f08'
use simple_aff_prop,            only: aff_prop
use simple_class_frcs,          only: class_frcs
use simple_cmdline,             only: cmdline
use simple_pftc_shsrch_fm,     only: pftc_shsrch_fm
use simple_polarft_calc,    only: polarft_calc
use simple_polarizer,           only: polarizer
use simple_corrmat
use simple_fileio
use simple_ftiter
use simple_gauss2Dfit
use simple_image
use simple_image_bin
use simple_neighs
use simple_parameters
use simple_parameters
use simple_pickseg
use simple_segmentation
use simple_srch_sort_loc
use simple_syslib
implicit none 
#include "simple_local_flags.inc"

logical  :: L_DEBUG = .false.
type :: image_afm
    private
    type(image),  allocatable :: img_array(:)
    type(string), allocatable :: img_names(:)
    type(string)              :: stack_string
contains
    procedure :: align_avg
    procedure :: get_AFM
    procedure :: pick_valid
    procedure :: read_ibw
    procedure :: zero_padding
end type image_afm 

contains

    subroutine read_ibw( AFM, fn_in )
        class(image_afm), intent(out) :: AFM
        class(string),    intent(in)  :: fn_in 
        real(kind=4), allocatable :: Rank3_Data_4byte(:, :, :, :)
        character(len=:), allocatable :: channel_info
        character(len=10) :: iteration(2), properties(4)
        integer :: in, check
        integer :: data, bytes_read, iter_ind, prop_ind, img_ind
        type :: bin_header5
            integer(kind=2) :: version 
            integer(kind=2) :: checksum 
            integer(kind=4) :: wfmSize 
            integer(kind=4) :: formulaSize 
            integer(kind=4) :: noteSize 
            integer(kind=4) :: dataEUnitsSize 
            integer(kind=4) :: dimEUnitsSize(4) 
            integer(kind=4) :: dimLabelsSize(4)
            integer(kind=4) :: sIndicesSize 
            integer(kind=4) :: optionsSize1 
            integer(kind=4) :: optionsSize2
        end type
        type :: wave_header5
            integer(kind=4)    :: structure_padding
            integer(C_INT32_T) :: creationDate
            integer(C_INT32_T) :: modDate
            integer(kind=4)    :: npnts 
            integer(kind=2)    :: type 
            integer(kind=2)    :: dlock 
            character          :: whpad(6)
            integer(kind=2)    :: whVersion 
            character          :: bname(32)
            integer(kind=4)    :: whpad2
            integer(kind=4)    :: DataFolder_Padding
            integer(kind=4)    :: nDim(4)
            real(kind=8)       :: sfA(4)
            real(kind=8)       :: sfB(4)
            character          :: dataUnits(4)
            character          :: dimUnits(4,4)
            integer(kind=2)    :: fsvalid
            integer(kind=2)    :: whpad3
            real(kind=8)       :: topFullScale, botFullScale
            integer(kind=4)    :: dataEUnits, dimEUnits(4), dimLabels(4), waveNoteH
            integer(kind=4)    :: whUnused(16)
            integer(kind=2)    :: aModified
            integer(kind=2)    :: wModified 
            integer(kind=2)    :: sModified  
            character          :: useBits
            character          :: kindBits
            integer(kind=4)    :: formula_pointer
            integer(kind=4)    :: depID
            integer(kind=2)    :: whpad4, srcFldr
            character(kind=4)  :: fileName
            integer(kind=4)    :: sIndices 
        end type
        type(bin_header5)  :: binheader
        type(wave_header5) :: waveheader  
        iteration  = [character(len = 10) :: 'Trace', 'Retrace']
        properties = [character(len = 10) :: 'Height', 'Amplitude', 'Phase', 'ZSensor' ]
        if( fn_in%substr_ind('.ibw') == 0 ) THROW_HARD('Only .ibw files are supported')
        open(newunit = check, file = fn_in%to_char(), status = 'old', access='stream')
        read(check) binheader%version
        if( binheader%version > 5 )then
#if USE_AFM
            open(newunit = in, file = fn_in%to_char(), status = 'old', access='stream', convert='swap')
#else
            open(newunit = in, file = fn_in%to_char(), status = 'old', access='stream')
#endif
        else
            open(newunit = in, file = fn_in%to_char(), status = 'old', access='stream')
        endif
        read(in)binheader, waveheader
        if( binheader%version /= 5 ) THROW_HARD('Only version 5 files are supported')
        allocate(Rank3_Data_4byte(waveheader%nDim(1) ,waveheader%nDim(2), 1, waveheader%nDim(3)))
        open(newunit = data, file = fn_in%to_char(), status = 'old', access='stream')
        read(data, pos = 385) Rank3_Data_4byte
        inquire(data, pos = bytes_read)
        allocate(character(binheader%dimLabelsSize(3)) :: channel_info)
        allocate(AFM%img_names(size(iteration)*size(properties)))
        read(data, pos = binheader%noteSize + bytes_read) channel_info
        ! getting image names 
        do prop_ind = 1, size(properties)
            do iter_ind = 1, size(iteration)
                if(     index(channel_info, trim(properties(prop_ind))//trim(iteration(iter_ind))) /= 0 .and. iter_ind == 2 )then 
                    AFM%img_names(iter_ind*prop_ind)     = trim(properties(prop_ind))//trim(iteration(iter_ind))
                elseif( index(channel_info, trim(properties(prop_ind))//trim(iteration(iter_ind))) /= 0 .and. iter_ind == 1 )then
                    AFM%img_names(2*prop_ind - iter_ind) = trim(properties(prop_ind))//trim(iteration(iter_ind)) 
                endif 
            enddo 
        enddo     
        allocate(AFM%img_array(waveheader%nDim(3)))
        do img_ind = 1, waveheader%nDim(3)
            call AFM%img_array(img_ind)%new([waveheader%nDim(1), waveheader%nDim(2), 1], real(waveheader%sfA(1)) * 10.**10.)
            call AFM%img_array(img_ind)%set_rmat(Rank3_Data_4byte(:, :, img_ind, :), .false.)
            call AFM%img_array(img_ind)%norm_minmax()
        enddo 
        deallocate(Rank3_Data_4byte)
    end subroutine read_ibw

    subroutine pick_valid( AFM_in, outname, avg_p )
        class(image_afm), intent(inout) :: AFM_in 
        class(string),    intent(in)    :: outname
        type(pickseg),    intent(inout) :: avg_p   
        integer, allocatable :: pickpos(:, :), val_center_r(:, :), val_center_t(:, :)
        real,    allocatable :: corr_final_t(:), corr_final_r(:)
        type(image)   :: HeightTrace, HeightRetrace, AvgHeight
        type(pickseg) :: trace_p, retrace_p
        type(image)   :: avg_slim, trace_slim, retrace_slim 
        integer       :: ldim_box(3), box_iter, search_iter, neighbor_iter, box_count
        real          :: smpd_box, neighbor_corr(8), coord_corr(3,8), max_corr(20)
        integer       :: coord_test(2), center_x(20), center_y(20)
        logical       :: outs
        integer       :: neighbor(3, 8), nsiz, center(3)
        real          :: corr_r, corr_t
        call get_AFM(AFM_in, 'AvgHeight',     AvgHeight)
        call get_AFM(AFM_in, 'HeightTrace',   HeightTrace)
        call get_AFM(AFM_in, 'HeightRetrace', HeightRetrace)
        call AvgHeight%norm_minmax()
        call HeightTrace%norm_minmax()
        call HeightRetrace%norm_minmax()
        call AvgHeight%write(outname//'avg.mrc')
        call avg_p%pick(outname//'avg.mrc',.true.)
        call HeightTrace%write(outname//'trace.mrc')
        call trace_p%pick(outname//'trace.mrc',.true.)
        call HeightRetrace%write(outname//'retrace.mrc')
        call retrace_p%pick(outname//'retrace.mrc',.true.)
        ldim_box = [avg_p%get_boxsize(), avg_p%get_boxsize(), 1]
        smpd_box = AvgHeight%get_smpd()
        call avg_slim%new(ldim_box, smpd_box )
        call trace_slim%new(ldim_box, smpd_box)
        call retrace_slim%new(ldim_box, smpd_box)
        box_count = 0 
        allocate(val_center_r(3, avg_p%get_nboxes()))
        allocate(val_center_t(3, avg_p%get_nboxes()))
        val_center_t = 0
        val_center_r = 0
        call avg_p%get_positions(pickpos) 
        allocate(corr_final_r(avg_p%get_nboxes()))
        allocate(corr_final_t(avg_p%get_nboxes()))
        do box_iter = 1, avg_p%get_nboxes() 
            coord_test = pickpos(box_iter, :)
            call AvgHeight%window_slim(coord_test, avg_p%get_boxsize(), avg_slim, outs)
            call HeightTrace%window_slim(coord_test, avg_p%get_boxsize(), trace_slim, outs)
            call HeightRetrace%window_slim(coord_test, avg_p%get_boxsize(), retrace_slim, outs)
            corr_r      =  avg_slim%real_corr(retrace_slim)
            corr_t      = avg_slim%real_corr(trace_slim)
            center      = [coord_test(1), coord_test(2), 1]
            center_x(1) = coord_test(1)
            center_y(1) = coord_test(2)
            ! print *, 'box: ', box_iter 
            if( corr_r < 0.1 .or. corr_t < 0.1 )then 
                ! print *, 'box is outside the image'
                cycle
            endif 
            box_count = box_count + 1 
            call nn_val(HeightTrace, corr_final_t, val_center_t)                                                 
            call nn_val(HeightRetrace, corr_final_r, val_center_r)
            ! check if line from hough transform is located within largest cc. if not, keep the box. 
        enddo 
        if( L_DEBUG )then 
            if( sum(corr_final_t) /  box_count > sum(corr_final_r) / box_count )then 
                print *, 'retrace is more noisy than trace'
            else
                print *, 'trace is more noisy than retrace'
            endif 
            print *, sum(corr_final_t) /  box_count, sum(corr_final_r) / box_count
            print *, val_center_t, val_center_r
            do box_iter = 73, 75
                call AvgHeight%window_slim(pickpos(box_iter, :), avg_p%get_boxsize(), avg_slim, outs)
                call avg_slim%vis()
                call HeightTrace%window_slim(val_center_t(:2, box_iter), avg_p%get_boxsize(), trace_slim, outs)
                call trace_slim%vis()
                call HeightRetrace%window_slim(val_center_r(:2, box_iter), avg_p%get_boxsize(), retrace_slim, outs)
                call retrace_slim%vis()
            enddo
        endif 
    
        contains 

            subroutine nn_val( im_iter, corr_final, val_centers )
                type(image), intent(in)  :: im_iter
                real,        intent(out) :: corr_final(avg_p%get_nboxes())
                integer,     intent(out) :: val_centers(3, avg_p%get_nboxes())
                ! output validated centers
                do search_iter = 2, 20
                    call neigh_8_1(AvgHeight%get_ldim(), center, neighbor, nsiz)
                    do neighbor_iter = 1, nsiz
                        call im_iter%window_slim([neighbor(1, neighbor_iter), neighbor(2, neighbor_iter)], avg_p%get_boxsize(), trace_slim, outs)
                        coord_corr(:, neighbor_iter) = [real(neighbor(1, neighbor_iter)), real(neighbor(2, neighbor_iter)), avg_slim%real_corr(trace_slim)]
                        neighbor_corr(neighbor_iter) = avg_slim%real_corr(trace_slim)
                    enddo 
                    center = [ int(coord_corr(1, maxloc(neighbor_corr))), int(coord_corr(2, maxloc(neighbor_corr))), 1 ] 
                    ! print *, center 
                    max_corr(search_iter) = maxval(neighbor_corr)
                    if( search_iter > 2 .and. max_corr(search_iter - 1) < max_corr(search_iter) )then 
                        ! print *, center, max_corr(search_iter)
                        val_centers(:, box_iter) = center 
                        corr_final(box_iter)     = max_corr(search_iter)
                        exit
                    endif 
                enddo
            end subroutine nn_val

    end subroutine pick_valid  

    subroutine get_AFM( AFM_Hash, key, image_at_key )
        class(image_afm), intent(in)  :: AFM_Hash
        character(*),     intent(in)  :: key 
        type(image),      intent(out) :: image_at_key
        image_at_key = AFM_Hash%img_array(findloc_str(AFM_Hash%img_names, key))
    end subroutine get_AFM

    subroutine zero_padding( AFM_pad )
        class(image_afm), intent(inout) :: AFM_pad
        integer :: img_ind
        integer :: dim(3)
        do img_ind = 1, size(AFM_pad%img_array)
            dim = AFM_pad%img_array(img_ind)%get_ldim()
            if( dim(1) /= dim(2) )then 
                call AFM_pad%img_array(img_ind)%pad_inplace([maxval(dim), maxval(dim), 1])
            endif
        enddo
    end subroutine zero_padding
    
    subroutine align_avg( AFM_in, Align_AFM )
        class(image_afm), intent(in)  :: AFM_in 
        class(image_afm), intent(out) :: Align_AFM
        real, allocatable :: shifts(:, :)
        integer           :: num_avg, avg_ind, new_size, count, num_mic, tr_ind, retr_ind, i 
        type(string)      :: new_name 
        num_mic  = size(AFM_in%img_array)
        new_size = int(num_mic*1.5)
        allocate(Align_AFM%img_array(new_size))
        allocate(Align_AFM%img_names(new_size))
        num_avg = num_mic/2
        allocate(shifts(2, num_avg))
        if(modulo(size(shifts), 2) /= 0) THROW_HARD( 'One or more traces are not present')
        do i = 1, num_mic
            Align_AFM%img_array(i) = AFM_in%img_array(i)
            Align_AFM%img_names(i) = AFM_in%img_names(i)
        enddo 
        count = 0
        do i = num_mic + 1, new_size
            Align_AFM%img_array(i) = AFM_in%img_array(i - num_mic + count)
            count = count + 1
        enddo
        do avg_ind = 1, size(Align_AFM%img_array)
            call Align_AFM%img_array(avg_ind)%fft()
        enddo
        count = 0
        do avg_ind = num_mic + 1, new_size
            count = count + 1
            call Align_AFM%img_array(avg_ind)%fcorr_shift(Align_AFM%img_array(avg_ind - num_mic + count), 20., shifts(:, count), .true.)
        enddo 
        do avg_ind = 1, num_avg
            call Align_AFM%img_array(2*avg_ind)%fft()
            call Align_AFM%img_array(2*avg_ind)%shift([shifts(1, avg_ind), shifts(2, avg_ind), 0.])
            call Align_AFM%img_array(2*avg_ind)%ifft()
        enddo
        do avg_ind = 1, new_size 
            if( Align_AFM%img_array(avg_ind)%is_ft() )then
                call Align_AFM%img_array(avg_ind)%ifft()
            endif    
        enddo 
        count = 0 
        do avg_ind = num_mic + 1, new_size 
            tr_ind                       = avg_ind - num_mic + count
            retr_ind                     = avg_ind - num_mic + 1 + count
            Align_AFM%img_array(avg_ind) = Align_AFM%img_array(tr_ind) * 0.5 + Align_AFM%img_array(retr_ind) * 0.5
            new_name                     = AFM_in%img_names(tr_ind)
            Align_AFM%img_names(avg_ind) = 'Avg'//new_name%to_char([1,new_name%strlen_trim()- 5])
            count = count + 1
        enddo 
    end subroutine align_avg

    ! Hough transform to identify lines in any angle. Output logical array of line positions to discard picks 
    subroutine hough_lines( img_in, theta_range, mask )
        class(image),   intent(inout) :: img_in
        real,           intent(out)   :: mask(:,:,:)
        real, optional, intent(in)    :: theta_range(2)
        type(image) :: img_denoised
        type(image) :: img_edge 
        real,    parameter   :: min_theta = -PI/2., theta_step = PI/180., rad_step = 1
        real,    allocatable :: angles(:), rad(:), curr_rads(:), sins(:), coss(:), emat(:,:,:), rmat(:,:,:), imat(:,:,:)
        integer, allocatable :: accumulator(:,:), line_pos(:,:)
        real    :: curr_rad, theta_range_def(2), smpd
        integer :: dims(3), diagonal, a_grid, r_grid, i, x, y, t, r, curr_rad_r, end_px, pix_cnt
        integer :: draw, line_num, traversed, min_line = 5
        logical :: debug_m = .false.
        if( present(theta_range) )then 
            theta_range_def = theta_range
        else
            theta_range_def = [-PI/2., PI/2.]
        endif
        dims = img_in%get_ldim()
        smpd = img_in%get_smpd()
        ! pre-processing
        call canny(img_in, img_edge)
        allocate(rmat(dims(1), dims(2), dims(3)))
        allocate(emat(dims(1), dims(2), dims(3)))
        allocate(imat(dims(1), dims(2), dims(3)))
        emat     = img_edge%get_rmat()
        ! max radius
        diagonal = ceiling(sqrt(real(dims(1))**2. + real(dims(2))**2.))
        r_grid   = 2*diagonal
        allocate(rad(r_grid))
        allocate(curr_rads(r_grid))
        rad    = [(-diagonal + (i - 1)*rad_step, i = 1, r_grid)]
        a_grid = nint(abs(theta_range_def(2) - theta_range_def(1))*theta_step**(-1)) + 1
        allocate(angles(a_grid))
        angles = [(theta_range_def(1) + (i - 1)*theta_step, i = 1, size(angles))]
        allocate(sins(size(angles)))
        allocate(coss(size(angles)))
        do i = 1, size(angles)
            sins(i) = sin(angles(i))
            coss(i) = cos(angles(i))
        enddo 
        allocate(accumulator(size(rad), size(angles)))
        accumulator = 0
        do x = 1, dims(1)
            do y = 1, dims(2)
                if( emat(x, y, 1) > 0. )then 
                    if( debug_m )then
                        print *, 'coordinate:', x,y
                    endif 
                    do t = 1, size(angles)
                        curr_rad   = x*coss(t) + y*sins(t)
                        curr_rads  = curr_rad
                        curr_rad_r = minloc(abs(rad - curr_rads), 1)
                        if( debug_m )then
                            write(logfhandle,'(2(A,f8.4))') 'radius:', curr_rad, 'angle:', angles(t)
                        endif 
                        accumulator(curr_rad_r,t) = accumulator(curr_rad_r,t) + 1                
                    enddo
                endif 
            enddo 
        enddo  
        call img_denoised%new(dims, smpd)
        allocate(line_pos(dims(1), dims(2)))
        line_pos = 0 
        ! finding local maxima
        do t = 1, size(angles)
            do r = 1, size(rad)
                if( accumulator(r, t) > 10 .and. angles(t) > PI/2. - 0.01 .and. angles(t) < PI/2. + 0.01 )then 
                    if( debug_m )then
                        print *, angles(t), accumulator(r,t), rad(r)*coss(t), rad(r)*sins(t)
                    endif 
                    do draw = 0, accumulator(r,t)
                        line_pos(1, nint(rad(r)*sins(t))) = accumulator(r,t)
                    enddo 
                endif 
            enddo 
        enddo 
        ! line alignment + multiple lines in row
        do line_num = 1, size(line_pos(1, :))
            if( line_pos(1, line_num) > 0 .and. line_pos(1, line_num) < dims(1) )then 
                x         = 1
                traversed = 0
                if( debug_m )then
                    print *, 'y = ', line_num
                endif 
                do while( line_pos(1, line_num) > 0 .and. x < dims(1) - min_line )
                    if( x < dims(1) .and. nint(sum(emat(x:x + min_line, line_num, 1))) == size(emat(x:x + min_line, line_num, 1))) then 
                        do while( emat(x, line_num, 1) > 0 .and. x < dims(1) - min_line ) 
                            traversed = traversed + 1
                            x = x + 1
                        enddo 
                        line_pos(1, line_num)             = line_pos(1, line_num) - traversed 
                        if( debug_m )then
                            print *, x, traversed 
                        endif 
                        line_pos(x - traversed, line_num) = traversed 
                        traversed                         = 0
                    endif 
                    x = x + 1
                enddo
            endif  
        enddo 
        ! connect close lines, final line detection
        pix_cnt = 5
        do x = 2, dims(1)
            do y = 1, dims(2)
                if( line_pos(x,y) + x < dims(1) )then
                    end_px = line_pos(x, y) + x
                    if( sum(line_pos(end_px:end_px + pix_cnt, y)) < 1. .and. end_px + pix_cnt < dims(1) )then 
                        line_pos(x, y) = line_pos(x, y) + sum(line_pos(end_px:end_px + pix_cnt, y))      
                    endif  
                    if( line_pos(x, y) > 50 )then 
                        do draw = 1, line_pos(x,y)
                            call img_denoised%set_rmat_at(x + draw, y, 1, 1.0)
                        enddo
                    endif 
                endif 
            enddo 
        enddo
        emat = img_denoised%get_rmat()
        ! rmat = img_in%get_rmat()
        mask = 0.
        where( emat > 0.5 )
            mask = 1. 
        end where
    end subroutine hough_lines 

    subroutine mask42D( img_in, AFM_pick, bin_cc, hough_mask, pick_vec )
        type(pickseg),  intent(in)    :: AFM_pick
        real,           intent(in)    :: hough_mask(:,:,:)
        type(image),    intent(inout) :: img_in 
        type(image_bin), intent(inout) :: bin_cc 
        type(image),    intent(inout) :: pick_vec(:)
        integer, allocatable :: pos(:, :) 
        real,    allocatable :: im_in_rmat(:,:,:), im_win_rmat(:,:,:), bin_win_rmat(:,:,:), area(:,:)
        type(image)    :: img_win, bin_win, lin_win, lines
        type(image_bin) :: bin_erode 
        real           :: smpd, new_cen(3)
        integer        :: ldim(3), ldim_pick(3), i, windim(3), j, counts, pad = 150
        logical        :: outside
        ldim   = img_in%get_ldim()
        smpd   = img_in%get_smpd()
        windim = [AFM_pick%get_boxsize(), AFM_pick%get_boxsize(), 1]
        call bin_cc%new(AFM_pick%get_ldim(), AFM_pick%get_smpd_shrink())
        call bin_cc%read(string('mic_shrink_lp_tv_bin_erode_cc.mrc'))
        call lines%new(ldim, smpd)
        call lines%set_rmat(hough_mask, .false.)
        ldim_pick = AFM_pick%get_ldim()
        if( ldim_pick(1) /= ldim_pick(2) )then 
            call bin_cc%pad_inplace( [maxval(AFM_pick%get_ldim()), maxval(AFM_pick%get_ldim()), 1])
        endif
        call bin_cc%pad_inplace([ldim(1) + pad, ldim(2) + pad, 1])
        call img_in%pad_inplace([ldim(1) + pad, ldim(2) + pad, 1])
        allocate(im_in_rmat(ldim(1),ldim(2),ldim(3)))
        im_in_rmat = img_in%get_rmat()
        where( nint(bin_cc%get_rmat()) < 1. )
            im_in_rmat = 0.
        end where 
        call img_in%set_rmat(im_in_rmat,.false.)
        allocate(pos(AFM_pick%get_nboxes(), 2))
        call img_win%new(windim,AFM_pick%get_smpd_shrink())
        call bin_win%new(windim,AFM_pick%get_smpd_shrink())
        call lin_win%new(windim,AFM_pick%get_smpd_shrink())
        call bin_erode%new(windim,AFM_pick%get_smpd_shrink())
        allocate(im_win_rmat(windim(1),windim(2),windim(3)))
        allocate(bin_win_rmat(windim(1),windim(2),windim(3)))
        allocate(area(AFM_pick%get_nboxes(),AFM_pick%get_nboxes()))
        ! allocate(pick_vec(AFM_pick%get_nboxes()))
        ! make copy, pad the copy in place and window slim on that 
        ! pad = 2*windim(1)
        area(:,:) = 0.
        counts    = 0
        do i = 1, AFM_pick%get_nboxes()
            call AFM_pick%get_positions(pos, i)
            outside = .false. 
            if( minval(pos(i,:)) < 0. ) outside = .true. 
            ! accounting for positions outside mic. 
            if( pos(i,1) < 0 .and. pos(i,2) < 0 )then
                pos(i,1) = pos(i,1) + pad 
                pos(i,2) = pos(i,2) + pad 
            elseif( pos(i,1) < 0 )then 
                pos(i,1) = pos(i,1) + pad 
            elseif( pos(i,2) < 0 )then  
                pos(i,2) = pos(i,2) + pad 
            endif 
            call bin_cc%window_slim(pos(i,:), AFM_pick%get_boxsize(), bin_win, outside)
            do j = 1, AFM_pick%get_nboxes()
                if( count(nint(bin_win%get_rmat()) == j) > 300 .and. area(i,j) /= 1. )then 
                    area(i,j) = count(nint(bin_win%get_rmat()) == j)
                endif
            enddo 
            if( i > 1 ) area(:i-1,maxloc(area(i,:),1)) = 1.
            area(i+1:,maxloc(area(i,:),1)) = 1.
            bin_win_rmat = 0.
            where(nint(bin_win%get_rmat()) == maxloc(area(i,:),1))
                bin_win_rmat = 1.
            end where 
            call bin_win%set_rmat(bin_win_rmat,.false.)
            call bin_win%masscen(new_cen)
            pos(i,:) = pos(i,:) + nint(new_cen)
            call bin_cc%window_slim(pos(i,:), AFM_pick%get_boxsize(), bin_win, outside)
            bin_win_rmat = bin_win%get_rmat()
            where(nint(bin_win%get_rmat()) /= maxloc(area(i,:),1))
                bin_win_rmat = 0.
            end where 
            call bin_win%set_rmat(bin_win_rmat,.false.)
            call img_in%window_slim(pos(i,:), AFM_pick%get_boxsize(), img_win, outside)
            im_win_rmat = img_win%get_rmat()
            where(bin_win%get_rmat() < 1.)
                im_win_rmat = 0.
            end where
            call img_win%set_rmat(im_win_rmat, .false.)
            call lines%window_slim(pos(i,:), AFM_pick%get_boxsize(), lin_win, outside)
            ! throw some hough lines away. can adjust parameters.
            if( lin_win%mean() > 0.01 .and. lin_win%real_corr(img_win) > 0. )then   
                cycle
            endif 
            ! can change to each indiv. particle radius 
            ! call img_win%mask(msk_rad,'soft') 
            if( img_win%mean() > 0.)then  
                pick_vec(i) = img_win
            endif 
        enddo 
    end subroutine mask42D

    function per_pix_var( img1, img2 ) result( var )
        type(image), intent(in) :: img1, img2
        real, allocatable :: rmat1(:,:,:), rmat2(:,:, :)
        integer :: ldim(3)
        real    :: var 
        ldim = img1%get_ldim()
        allocate(rmat1(ldim(1), ldim(2), ldim(3)), rmat2(ldim(1), ldim(2), ldim(3)), source = 0.)
        rmat1 = img1%get_rmat()
        rmat2 = img2%get_rmat()
        rmat1 = abs(rmat1 - rmat2)
        var   = norm2(rmat1)
    end function per_pix_var

    ! searches for pick from one micrograph in another similar micrograph
    subroutine pick_search( trace_picks, trace, retrace, retrace_centers )
        type(image),   intent(inout) :: trace, retrace
        type(pickseg), intent(in)    :: trace_picks  
        integer,       intent(out)   :: retrace_centers(:,:)
        integer, allocatable :: coords(:,:)
        type(image) :: trace_slim, retrace_slim
        integer     :: search_iter, neighbor_iter, box_iter, center(3), nsiz, neighbor(3,8), dim(3), dim_shrink(3)
        logical     :: outside
        real        :: neighbor_corr(8), coord_corr(3,8), max_corr(20), smpd, smpd_shrink, shrink = 1.
        call trace_picks%get_positions(coords)
        dim  = trace%get_ldim()
        smpd = trace%get_smpd()
        dim_shrink(1) = round2even(real(dim(1))/shrink)
        dim_shrink(2) = round2even(real(dim(2))/shrink)
        dim_shrink(3) = 1
        smpd_shrink = smpd * shrink
        call trace%fft()
        call trace%clip_inplace(dim_shrink)
        call trace%ifft()
        call trace%fft()
        call trace%clip_inplace(dim_shrink)
        call trace%ifft()
        call trace_slim%new([trace_picks%get_boxsize(), trace_picks%get_boxsize(), 1], smpd_shrink)
        call retrace_slim%new([trace_picks%get_boxsize(), trace_picks%get_boxsize(), 1], smpd_shrink)
        do box_iter = 1, trace_picks%get_nboxes() 
            center  = [coords(box_iter, 1), coords(box_iter, 2), 1]
            print *, 'trace center', center
            if( center(1) < 0 ) center(1) = center(1) + dim_shrink(1)/4
            if( center(2) < 0 ) center(2) = center(2) + dim_shrink(2)/4
            call trace%window_slim([center(1), center(2)], trace_picks%get_boxsize(), trace_slim, outside)
            do search_iter = 2, 10
                call neigh_8_1(trace%get_ldim(), center, neighbor, nsiz)
                do neighbor_iter = 1, nsiz
                    call retrace%window_slim([neighbor(1, neighbor_iter), neighbor(2, neighbor_iter)], trace_picks%get_boxsize(), retrace_slim, outside)
                    coord_corr(:, neighbor_iter) = [real(neighbor(1, neighbor_iter)), real(neighbor(2, neighbor_iter)), trace_slim%real_corr(retrace_slim)]
                    neighbor_corr(neighbor_iter) = retrace_slim%real_corr(trace_slim)
                enddo 
                center = [ int(coord_corr(1, maxloc(neighbor_corr))), int(coord_corr(2, maxloc(neighbor_corr))), 1 ] 
                max_corr(search_iter) = maxval(neighbor_corr)
                if( search_iter > 2 .and. max_corr(search_iter - 1) < max_corr(search_iter) )then 
                    retrace_centers(:, box_iter) = [center(1), center(2)]
                    exit
                elseif( search_iter > 8 )then 
                    retrace_centers(:, box_iter) = [coords(box_iter, 1), coords(box_iter, 2), 1]
                    exit 
                endif 
            enddo
            print *, 'center found:', retrace_centers(:,box_iter)
        enddo  
    end subroutine pick_search

    subroutine mask_incl_retrace( trace_pick, trace, retrace, trace_vec, retrace_vec ) 
        type(pickseg), intent(inout) :: trace_pick
        type(image),   intent(inout) :: trace, retrace
        type(image),   intent(inout) :: trace_vec(:), retrace_vec(:) 
        integer, allocatable :: coords(:,:)
        integer        :: ldim(3), windim(3)
        real           :: smpd
        type(image_bin) :: retrace_image_bin
        call trace_pick%get_positions(coords)
        ldim   = trace_pick%get_ldim()
        smpd   = trace_pick%get_smpd_shrink()
        windim = trace_pick%get_ldim()
        call retrace%mul(real(product(ldim)))
        call retrace%fft()
        call retrace_image_bin%new_bimg(ldim, smpd)
        call retrace_image_bin%set_ft(.true.)
        call retrace%clip(retrace_image_bin)
        call retrace_image_bin%bp(0., params_glob%lp)
        call retrace%ifft()
        call retrace_image_bin%ifft()
        call otsu_img(retrace_image_bin)
        call trace%pad_inplace([2*ldim(1), 2*ldim(2), 1])
        call retrace_image_bin%pad_inplace([2*ldim(1), 2*ldim(2), 1])
        call retrace_image_bin%erode()
        call retrace_image_bin%erode()
        call retrace_image_bin%dilate()
        call retrace_image_bin%vis()
        ! pad images for negative box pos.
        ! take trace positions as center of retrace boxes.
        ! mask in the same way
        ! need to apply same procedures to generate mic_shrink_lp_tv_bin_erode_cc.mrc
        ! band pass 
        ! probably get rid of TV denoising 
        ! otsu 
        ! erode, erode, dilate
        ! find ccs labels the connected components on the processed micrograph.
        ! mask other particles. 
    end subroutine mask_incl_retrace

    subroutine gau_fit( im_in, ncls )
        type(image),  intent(in) :: im_in 
        integer,      intent(in) :: ncls ! number of gaussians
        type(image), allocatable :: clus_stk(:), gau2D_stk(:)
        real,        allocatable :: rmat(:,:,:), euc_dist(:,:), rmat_labl(:,:,:)
        integer,     allocatable :: labels(:), centers(:), x(:), y(:), npnts(:), temp_vec(:)
        logical,     allocatable :: mask(:)
        type(image) :: gau2D_sum
        real        :: smpd, rand, cen_gauss(2), cov_gauss(2,2), corr
        integer     :: ldim(3), i, j, ndat, maxits = 100, iter, l, counts, sum 
        call im_in%vis()
        ldim = im_in%get_ldim()
        ldim = [ldim(1),ldim(2),1]
        allocate(rmat(ldim(1),ldim(2),1))
        rmat = im_in%get_rmat()
        ndat = count(rmat(:,:,1) > 0.1)
        allocate(x(ndat),y(ndat), source = 0)
        ndat = 0 
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if( rmat(i,j,1) > 0.1 )then 
                    ndat = ndat + 1
                    x(ndat) = i 
                    y(ndat) = j
                endif 
            enddo 
        enddo 
        ! 2D kmeans 
        allocate(euc_dist(ndat,ncls), source = 0.)
        allocate(labels(ndat), centers(ncls), source = 0)
        allocate(mask(ndat), source = .false.)
        allocate(npnts(ncls))
        ! rand init
        do i = 1, ncls
            call random_number(rand)
            centers(i) = nint(ndat * rand)
        enddo
        do i = 1, ncls - 1
            npnts(i) = floor(real(ndat/ncls))
        enddo 
        npnts(ncls) = npnts(ncls - 1) + mod(ndat,ncls)
        do iter = 1, maxits
            do i = 1, ncls 
                do j = 1, ndat 
                    euc_dist(j,i) = sqrt(real((x(j) - x(centers(i)))**2 + (y(j) - y(centers(i))))**2)
                enddo
            enddo 
            do i = 1, ncls 
                allocate(temp_vec(npnts(i)))
                temp_vec = minnloc(euc_dist(:,i),npnts(i))
                do j = 1, npnts(i)
                    labels(temp_vec(j)) = i
                enddo 
                deallocate(temp_vec)
            enddo 
            ! calc new center
            do i = 1, ncls
                sum = 0
                ! need sum of indices
                do j = 1, ndat 
                    if(labels(j) == i ) sum = sum + j
                enddo 
                centers(i) = nint(real(sum / npnts(i)))
            enddo 
        enddo 
        ! map labels back to image
        counts = 0
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if(rmat(i,j,1) > 0.1) then 
                    rmat(i,j,1) = labels(counts)
                    counts = counts + 1
                else
                    rmat(i,j,1) = 0. 
                endif
            enddo 
        enddo  
        allocate(clus_stk(ncls))
        allocate(gau2D_stk(ncls))
        allocate(rmat_labl(ldim(1), ldim(2), 1), source = 0.)
        call gau2D_sum%new(ldim, smpd)
        do l = 1, ncls
            rmat_labl = rmat
            do i = 1, ldim(1)
                do j = 1, ldim(2)
                    if( nint(rmat_labl(i,j,1)) /= l) rmat_labl(i,j,1) = 0.
                enddo 
            enddo
            call clus_stk(l)%new(ldim, smpd)
            call clus_stk(l)%set_rmat(rmat_labl, .false.)
            call gauss2Dfit(clus_stk(l), cen_gauss, cov_gauss, corr, gau2D_stk(l))
            call otsu_img(gau2D_stk(l))
            gau2D_stk(l) = gau2D_stk(l)*l
            call gau2D_sum%add(gau2D_stk(l))
        enddo 
        call gau2D_sum%vis()
    end subroutine gau_fit

end module simple_image_afm 

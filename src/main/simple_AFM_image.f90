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

type :: AFM_image
    type(image), allocatable :: img_array(:)
    character(len = 50), allocatable :: img_names(:)
contains
    procedure   :: pick_valid1
    procedure   :: get_AFM1
end type AFM_image 
contains
    subroutine pick_valid1(AFM_in)
        class(AFM_image), intent(inout)    :: AFM_in 
        type(image)                        :: HeightTrace, HeightRetrace, AvgHeight
        CHARACTER(len=255)                 :: cwd
        type(pickseg)                      :: avg_p, trace_p, retrace_p
        type(image)                        :: avg_slim, trace_slim, retrace_slim 
        integer                            :: ldim_box(3), box_iter, search_iter, neighbor_iter, box_count
        real                               :: smpd_box, neighbor_corr(8), coord_corr(3,8), max_corr(20)
        integer, allocatable               :: pickpos(:, :), val_centers(:, :)   
        integer                            :: coord_test(2), center_x(20), center_y(20)
        logical                            :: outs = .true. 
        integer                            :: neighbor(3, 8), nsiz, center(3)
        real                               :: corr_r, corr_t, val_score_t, val_score_r 
        ! real, allocatable                  :: corr_final(:)   
        real, allocatable                  :: corr_final_t(:), corr_final_r(:)
        

        call get_AFM1(AFM_in, 'HeightRetrace', HeightRetrace)
        call get_AFM1(AFM_in, 'HeightTrace', HeightTrace)
        call get_AFM1(AFM_in, 'AvgHeight', AvgHeight)

        call HeightRetrace%norm_minmax()
        call HeightTrace%norm_minmax()
        call AvgHeight%norm_minmax()


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
        ! allocate(corr_final(avg_p%get_nboxes()))
        allocate(val_centers(3, avg_p%get_nboxes()))
        call avg_p%get_positions(pickpos) 
        allocate(corr_final_r(avg_p%get_nboxes()))
        allocate(corr_final_t(avg_p%get_nboxes()))
        do box_iter = 1, avg_p%get_nboxes() 
            !  print *, 'box: ', box_iter 
            coord_test = pickpos(box_iter, :)
            call AvgHeight%window_slim(coord_test, avg_p%box_raw, avg_slim, outs)
            call HeightTrace%window_slim(coord_test, avg_p%box_raw, trace_slim, outs)
            call HeightRetrace%window_slim(coord_test, avg_p%box_raw, retrace_slim, outs)

            corr_r =  avg_slim%real_corr(retrace_slim)
            corr_t = avg_slim%real_corr(trace_slim)
            center = [coord_test(1), coord_test(2), 1]
            center_x(1) = coord_test(1)
            center_y(1) = coord_test(2)
            
            if( corr_r < 0.1 .or. corr_t < 0.1) then 
                ! print *, 'box is outside the image'
                cycle
            end if 
            box_count = box_count + 1 

            call nn_val(HeightTrace, corr_final_t)
            
            
            call nn_val(HeightRetrace, corr_final_r)
            
        end do 
        if(sum(corr_final_t) /  box_count > sum(corr_final_r) / box_count) then 
            print *, 'retrace is more noisy than trace'
        else
            print *, 'trace is more noisy than retrace'
        end if 
        contains 
            subroutine nn_val(im_iter, corr_final)
                type(image), intent(in)     :: im_iter
                real, intent(out) :: corr_final(avg_p%get_nboxes())
                ! output validated centers. add debug for some print, vis. 
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
    end subroutine pick_valid1     

    subroutine get_AFM1(AFM_Hash, key, image_at_key)
        class(AFM_image), intent(in) :: AFM_Hash
        character(*), intent(in)    :: key 
        type(image), intent(out)  :: image_at_key
        image_at_key = AFM_Hash%img_array(findloc(index(AFM_Hash%img_names, key),1, dim = 1))
    end subroutine get_AFM1

end module simple_AFM_image 
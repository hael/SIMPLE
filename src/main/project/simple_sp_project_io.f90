submodule(simple_sp_project) simple_sp_project_io
implicit none
#include "simple_local_flags.inc"
contains

    ! Printers

    module subroutine print_info( self, fname )
        class(sp_project),     intent(inout) :: self
        class(string),         intent(in)    :: fname
        type(binoris_seginfo), allocatable   :: hinfo(:)
        integer,               allocatable   :: seginds(:)
        type(string) :: projfile, record
        integer :: i
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//fname%to_char()//' not supported; sp_project :: print_info')
        endif
        projfile = fname
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            call self%bos%get_segments_info(seginds, hinfo)
            if( allocated(hinfo) )then
                do i = 1,size(hinfo)
                    call self%bos%read_record(seginds(i), hinfo(i)%first_data_byte, record)
                    write(logfhandle,'(a)') format_str(format_str('SEGMENT '//int2str_pad(seginds(i),2)//' of type: '//segment2oritype(seginds(i)), C_BOLD), C_UNDERLINED)
                    write(logfhandle,'(a)') format_str(segment2info(seginds(i), hinfo(i)%n_records), C_ITALIC)
                    write(logfhandle,'(a)') format_str('first record:', C_BOLD)//' '//record%to_char()
                end do
            endif
            call self%bos%close
        else
            THROW_HARD('projfile: '//projfile%to_char()//' nonexistent; print_info')
        endif
    end subroutine print_info

    module subroutine print_info_json( self, fname )
        class(sp_project), intent(inout)   :: self
        class(string),     intent(in)      :: fname
        type(string),          allocatable :: keys(:)
        type(binoris_seginfo), allocatable :: hinfo(:)
        integer,               allocatable :: seginds(:)
        type(json_value),      pointer     :: json_root, json_seg, json_real_keys, json_char_keys
        type(oris)      :: vol_oris
        type(json_core) :: json
        type(string)    :: projfile, record
        type(ori)       :: seg_ori
        logical         :: is_ptcl = .false.
        integer         :: i, j, noris
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//fname%to_char()//' not supported; sp_project :: print_info_json')
        endif
        projfile = fname
        call json%initialize(no_whitespace=.true.)
        call json%create_object(json_root,'')
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            call self%bos%get_segments_info(seginds, hinfo)
            if( allocated(hinfo) )then
                do i = 1,size(hinfo)
                    if(hinfo(i)%n_records .gt. 0) then
                        ! initialise json
                        call json%create_object(json_seg, segment2oritype(seginds(i)))
                        call json%add(json_seg, 'n', int(hinfo(i)%n_records))
                        call json%add(json_seg, 'info',segment2info(seginds(i), hinfo(i)%n_records))
                        call json%create_array(json_real_keys, 'numeric_keys')
                        call json%create_array(json_char_keys, 'character_keys')
                        ! read 1st record and convert to ori
                        if(seginds(i) == PTCL2D_SEG .or. seginds(i) == PTCL3D_SEG) is_ptcl = .true.
                        call self%bos%read_record(seginds(i), hinfo(i)%first_data_byte, record)
                        call seg_ori%str2ori(record%to_char(), is_ptcl)
                        ! get keys and test if real or character
                        keys = seg_ori%get_keys()
                        do j = 1, size(keys)
                            if(seg_ori%ischar(keys(j)%to_char())) then
                                call json%add(json_char_keys, '', keys(j)%to_char())
                            else
                                call json%add(json_real_keys, '', keys(j)%to_char())
                            end if
                        end do 
                        ! add to json
                        call json%add(json_seg, json_real_keys)
                        call json%add(json_seg, json_char_keys)
                        call json%add(json_root, json_seg)
                        ! add vols section
                        if(seginds(i) == OUT_SEG) then
                            call self%os_out%new(int(hinfo(i)%n_records), is_ptcl=.false.)
                            call self%bos%read_segment(OUT_SEG, self%os_out)
                            call self%get_all_vols( vol_oris )
                            noris = vol_oris%get_noris()
                            if( noris > 0 ) then
                                call json%create_object(json_seg, 'vols')
                                call json%add(json_seg, 'n', noris)
                                call json%add(json_root, json_seg)
                            end if
                        end if
                        ! clean up
                        is_ptcl = .false.
                        call keys%kill
                    end if
                end do
            endif
            call self%bos%close
        else
            THROW_HARD('projfile: '//projfile%to_char()//' nonexistent; print_info_json')
        endif
        call json%print(json_root, logfhandle)
        write(logfhandle,*) !need a newline else the footer is on same line as json
    end subroutine print_info_json

    module subroutine print_segment( self, oritype, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        integer, optional, intent(in)    :: fromto(2)
        type(string) :: str
        integer      :: ffromto(2), iori, noris
        logical      :: fromto_present
        fromto_present = present(fromto)
        if( fromto_present ) ffromto = fromto
        select case(trim(oritype))
            case('mic')
                noris = self%os_mic%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_mic%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No mic-type oris available to print; sp_project :: print_segment'
                endif
            case('stk')
                noris = self%os_stk%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_stk%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No stk-type oris available to print; sp_project :: print_segment'
                endif
            case('ptcl2D')
                noris = self%os_ptcl2D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_ptcl2D%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No ptcl2D-type oris available to print; sp_project :: print_segment'
                endif
            case('cls2D')
                noris = self%os_cls2D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_cls2D%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No cls2D-type oris available to print; sp_project :: print_segment'
                endif
            case('cls3D')
                noris = self%os_cls3D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_cls3D%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No cls3D-type oris available to print; sp_project :: print_segment'
                endif
            case('ptcl3D')
                noris = self%os_ptcl3D%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_ptcl3D%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No ptcl3D-type oris available to print; sp_project :: print_segment'
                endif
            case('out')
                noris = self%os_out%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_out%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No out-type oris available to print; sp_project :: print_segment'
                endif
            case('optics')
                noris = self%os_optics%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%os_optics%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No optics-type oris available to print; sp_project :: print_segment'
                endif
            case('projinfo')
                noris = self%projinfo%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%projinfo%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No projinfo-type oris available to print; sp_project :: print_segment'
                endif
            case('jobproc')
                noris = self%jobproc%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%jobproc%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No jobproc-type oris available to print; sp_project :: print_segment'
                endif
            case('compenv')
                noris = self%compenv%get_noris()
                if( noris > 0 )then
                    if( .not. fromto_present ) ffromto = [1,noris]
                    do iori=ffromto(1),ffromto(2)
                        str = self%compenv%ori2str(iori)
                        write(logfhandle,'(a)') str%to_char()
                        call str%kill
                    end do
                else
                    write(logfhandle,*) 'No compenv-type oris available to print; sp_project :: print_segment'
                endif
            case DEFAULT
                THROW_HARD('unsupported oritype flag; print_segment')
        end select
    end subroutine print_segment

    module subroutine print_segment_json( self, oritype, projfile, fromto, sort_key, sort_asc, hist, nran, boxes, plot_key )
        class(sp_project),           intent(inout) :: self
        character(len=*),            intent(in)    :: oritype
        class(string),               intent(in)    :: projfile
        character(len=*),  optional, intent(in)    :: sort_key, sort_asc, hist, plot_key
        integer,           optional, intent(in)    :: fromto(2), nran
        logical,           optional, intent(in)    :: boxes
        type(json_core)                            :: json
        type(json_value),  pointer                 :: json_root, json_data, json_hist, json_ori, json_pre, json_post, json_plot
        type(oris)                                 :: vol_oris, fsc_oris
        type(ori)                                  :: tmp_ori
        real,              allocatable             :: projections(:,:)
        integer,           allocatable             :: indices(:), indices_pre(:), indices_post(:)
        real                                       :: rnd
        integer                                    :: ffromto(2), iori, noris, boxsize
        logical                                    :: fromto_present, sort, sort_ascending, copy_oris, l_boxes
        fromto_present = present(fromto)
        if( fromto_present ) ffromto = fromto
        l_boxes = .false.
        if(present(boxes)) l_boxes = boxes
        sort    = .false.
        if( present(sort_key) ) then
            if(sort_key .ne. '' .and. sort_key .ne. 'n') sort = .true.
        end if
        sort_ascending = .true.
        if( present(sort_asc) ) then
            if(sort_asc .eq. "no") sort_ascending = .false.
        end if
        call json%initialize(no_whitespace=.true.)
        call json%create_object(json_root,'')
        call json%create_array(json_data, 'data')
        select case(trim(oritype))
            case('mic')
                noris = self%os_mic%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%os_mic)
                    do iori=1, size(indices)
                        call self%os_mic%ori2json(indices(iori), json_ori, boxes=l_boxes)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(allocated(indices_pre)) then
                        call json%create_array(json_pre, 'indices_pre')
                        do iori=1, size(indices_pre)
                            call json%add(json_pre, '', indices_pre(iori))
                        enddo
                        call json%add(json_root, json_pre)
                        deallocate(indices_pre)
                    end if
                    if(allocated(indices_post)) then
                        call json%create_array(json_post, 'indices_post')
                        do iori=1, size(indices_post)
                            call json%add(json_post, '', indices_post(iori))
                        enddo
                        call json%add(json_root, json_post)
                        deallocate(indices_post)
                    end if
                    if(present(hist))     call calculate_histogram(self%os_mic)
                    if(present(plot_key)) call calculate_plot(self%os_mic)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No mic-type oris available to print; sp_project :: print_segment_json'
                endif
            case('stk')
                noris = self%os_stk%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%os_stk)
                    do iori=1, size(indices)
                        call self%os_stk%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%os_stk)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No stk-type oris available to print; sp_project :: print_segment_json'
                endif
            case('ptcl2D')
                copy_oris=.false.
                noris = self%os_ptcl2D%get_noris()
                if( noris > 0 )then
                    call self%read_segment('stk', projfile)
                    if(present(nran)) then
                        if(nran .gt. 0) then
                            ! randomly choose nran particles
                            call random_seed()
                            allocate(indices(nran))
                            iori = 1
                            do while (iori <= nran)   
                                call random_number(rnd)
                                if(self%os_ptcl2D%get(ceiling(noris * rnd), "state") .gt. 0.0) then
                                    indices(iori) = ceiling(noris * rnd)
                                    iori = iori + 1
                                end if
                            end do
                            call self%set_ptcl2D_thumb(projfile, indices, boxsize)
                            do iori=1, size(indices)
                                call self%os_ptcl2D%get_ori(indices(iori), tmp_ori)
                                call tmp_ori%ori2json(json_ori)
                                call json%add(json_ori, 'n', iori) 
                                call json%add(json_ori, "thumbn",   size(indices))
                                call json%add(json_ori, "thumbdim", JPEG_DIM)
                                call json%add(json_ori, "thumbidx", iori)
                                call json%add(json_ori, "box",      boxsize)
                                call json%add(json_data, json_ori)
                            end do
                            copy_oris = .true.
                        end if
                    end if
                    if(.not. copy_oris) then
                        call calculate_indices(self%os_ptcl2D)
                        do iori=1, size(indices)
                            call self%os_ptcl2D%print(iori)
                            call self%os_ptcl2D%ori2json(indices(iori), json_ori)
                            call json%add(json_data, json_ori)
                        end do
                    end if
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%os_ptcl2D)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No ptcl2D-type oris available to print; sp_project :: print_segment_json'
                endif
            case('cls2D')
                noris = self%os_cls2D%get_noris()
                if( noris > 0 )then
                    if(.not. self%os_cls2D%isthere(1, "thumb")) then
                        ! create thumb
                        call self%read_segment('out', projfile)
                        call self%set_cavgs_thumb(projfile)
                        call self%write_segment_inside('cls2D', fname=projfile)
                    end if
                    call calculate_indices(self%os_cls2D)
                    do iori=1, size(indices)
                        if(self%os_cls2D%get_state(indices(iori)) .gt. 0) then
                            call self%os_cls2D%ori2json(indices(iori), json_ori)
                            call json%add(json_data, json_ori)
                        end if
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist))     call calculate_histogram(self%os_cls2D)
                    if(present(plot_key)) call calculate_plot(self%os_cls2D)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No cls2D-type oris available to print; sp_project :: print_segment_json'
                endif
            case('cls3D')
                noris = self%os_cls3D%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%os_cls3D)
                    do iori=1, size(indices)
                        call self%os_cls3D%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%os_cls3D)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No cls3D-type oris available to print; sp_project :: print_segment_json'
                endif
            case('ptcl3D')
                noris = self%os_ptcl3D%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%os_ptcl3D)
                    do iori=1, size(indices)
                        call self%os_ptcl3D%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%os_ptcl3D)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No ptcl3D-type oris available to print; sp_project :: print_segment_json'
                endif
            case('out')
                noris = self%os_out%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%os_out)
                    do iori=1, size(indices)
                        call self%os_out%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%os_out)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No out-type oris available to print; sp_project :: print_segment_json'
                endif
            case('optics')
                noris = self%os_optics%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%os_optics)
                    do iori=1, size(indices)
                        call self%os_optics%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%os_optics)
                    call self%read_segment('mic', projfile)
                    call calculate_optics_plot()
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No optics-type oris available to print; sp_project :: print_segment_json'
                endif
            case('projinfo')
                noris = self%projinfo%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%projinfo)
                    do iori=1, size(indices)
                        call self%projinfo%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%projinfo)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No projinfo-type oris available to print; sp_project :: print_segment_json'
                endif
            case('jobproc')
                noris = self%jobproc%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%jobproc)
                    do iori=1, size(indices)
                        call self%jobproc%ori2json(indices(iori), json_ori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%jobproc)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No jobproc-type oris available to print; sp_project :: print_segment_json'
                endif
            case('compenv')
                noris = self%compenv%get_noris()
                if( noris > 0 )then
                    call calculate_indices(self%compenv)
                    do iori=1, size(indices)
                        call self%compenv%ori2json(indices(iori), json_ori)
                        call json%add(json_root, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%compenv)
                    deallocate(indices)
                else
                    write(logfhandle,*) 'No compenv-type oris available to print; sp_project :: print_segment_json'
                endif
            case('vol')
                call self%get_all_vols( vol_oris )
                call self%get_all_fscs( fsc_oris )
                noris = vol_oris%get_noris() 
                if( noris > 0 ) then
                    call self%read_segment('ptcl3D', projfile)
                    if(self%os_ptcl3D%get_noris() .gt. 0) call get_projections(noris)
                    do iori=1, noris
                        call vol_oris%ori2json(iori, json_ori)
                        call add_fsc(iori)
                        if(self%os_ptcl3D%get_noris() .gt. 0) call add_oriplot(iori)
                        call json%add(json_data, json_ori)
                    end do
                    call json%add(json_root, json_data)
                    if(present(hist)) call calculate_histogram(self%jobproc)
                    if(allocated(projections))       deallocate(projections)
                else
                    write(logfhandle,*) 'No volumes available to print; sp_project :: print_segment_json'
                endif
                call vol_oris%kill
                call fsc_oris%kill
            case DEFAULT
                THROW_HARD('unsupported oritype flag; print_segment_json')
        end select
        call json%print(json_root, logfhandle)
        write(logfhandle, *)
        call tmp_ori%kill
        contains

            subroutine calculate_histogram( seg_oris )
                use, intrinsic :: iso_c_binding 
                include 'simple_lib.f08'
                type(oris),        intent(in)  :: seg_oris
                type(histogram)                :: histgrm
                type(json_value),  pointer     :: data, labels
                real,              allocatable :: rvec(:)
                integer                        :: n_bins, i
                n_bins = 20
                if(hist .eq. 'yes') then
                    if((.not. sort_key .eq. '') .and. (.not. sort_key .eq. 'n')) then
                        call json%create_object(json_hist,'histogram')
                        call json%create_array(data,     "data")
                        call json%create_array(labels,   "labels")
                        rvec = seg_oris%get_all(sort_key)
                        call histgrm%new(n_bins, rvec)
                        do i=1, n_bins
                            call json%add(data,   '', dble(histgrm%get(i)))
                            call json%add(labels, '', dble(histgrm%get_x(i)))
                        end do
                        call json%add(json_hist, data)
                        call json%add(json_hist, labels)
                        call json%add(json_root, json_hist)
                        call histgrm%kill()
                    end if
                end if
            end subroutine calculate_histogram

            subroutine calculate_plot( seg_oris )
                type(oris),        intent(in)  :: seg_oris
                type(json_value),  pointer     :: data, xy
                real,              allocatable :: rvecx(:), rvecy(:)
                integer                        :: i
                if((.not. sort_key .eq. '') .and. (.not. plot_key .eq. '')) then
                    call json%create_object(json_plot, 'plot')
                    call json%create_array(data,     "data")
                    rvecx = seg_oris%get_all(sort_key)
                    rvecy = seg_oris%get_all(plot_key)
                    if(size(rvecx) .eq. size(rvecy)) then
                        do i = 1, size(rvecx)
                            call json%create_object(xy, 'xy')
                            if(sort_key == 'n') then
                                call json%add(xy, 'x', dble(i))
                            else
                                call json%add(xy, 'x', dble(rvecx(i)))
                            endif
                            call json%add(xy, 'y', dble(rvecy(i)))
                            call json%add(data, xy)
                        enddo
                    endif
                    call json%add(json_plot, data)
                    call json%add(json_root, json_plot)
                endif
            end subroutine calculate_plot

            subroutine calculate_optics_plot()
                type(json_value),  pointer :: optics_plot, datasets, dataset, data, xy
                integer                    :: i, j
                if(self%os_optics%get_noris() .eq. 0)   return
                if(self%os_mic%get_noris()    .eq. 0)   return
                if(.not. self%os_mic%isthere('ogid'))   return
                if(.not. self%os_mic%isthere('shiftx')) return
                if(.not. self%os_mic%isthere('shifty')) return
                call json%create_object(optics_plot, 'assignments')
                call json%add(optics_plot, 'type', 'plot_scatter')
                call json%create_array(datasets, 'datasets')
                do i = 1, self%os_optics%get_noris()
                    call json%create_object(dataset, 'dataset')
                    call json%create_array(data, 'data')
                    call json%add(dataset, 'label', 'optics group ' // int2str(i))
                    do j = 1, self%os_mic%get_noris()
                        if(self%os_mic%get(j, 'ogid') == i) then
                            call json%create_object(xy, 'xy')
                            call json%add(xy, 'x', dble(self%os_mic%get(j, 'shiftx')))
                            call json%add(xy, 'y', dble(self%os_mic%get(j, 'shifty')))
                            call json%add(data, xy)
                        end if
                    end do
                    call json%add(dataset,  data)
                    call json%add(datasets, dataset)
                end do
                call json%add(optics_plot, datasets)
                call json%add(json_root,   optics_plot)
            end subroutine calculate_optics_plot

            subroutine calculate_indices( seg_oris )
                type(oris), intent(in)  :: seg_oris
                integer,    allocatable :: order(:)
                integer                 :: i
                if( .not. fromto_present ) ffromto = [1,noris]
                if( .not. sort_ascending)  ffromto = [noris - ffromto(2) + 1, noris - ffromto(1) + 1]
                allocate(indices(ffromto(2) - ffromto(1) + 1))
                if(ffromto(1) .gt. 1)     allocate(indices_pre(ffromto(1) - 1))
                if(ffromto(2) .lt. noris) allocate(indices_post(noris - fromto(2)))
                if(sort) then
                    order = sort_oris(seg_oris)
                    indices(:) = order(ffromto(1):ffromto(2))
                    if(allocated(indices_pre))  indices_pre(:)  = order(:ffromto(1) - 1)
                    if(allocated(indices_post)) indices_post(:) = order(ffromto(2) + 1:)
                    deallocate(order)
                else
                    do i=1, size(indices)
                        indices(i) = ffromto(1) + i - 1
                    end do
                    if(allocated(indices_pre)) then
                        do i=1, size(indices_pre)
                            indices_pre(i) = i
                        end do
                    end if
                    if(allocated(indices_post)) then
                        do i=1, size(indices_post)
                            indices_post(i) = ffromto(2) + i
                        end do
                    end if
                end if
                if(.not. sort_ascending) call reverse(indices)
            end subroutine calculate_indices

            subroutine add_fsc( iori_l )
                integer,          intent(in)  :: iori_l
                type(json_value), pointer     :: fsc_json, datasets, dataset, data, labels
                type(string) :: fscfile
                real,             allocatable :: fsc(:), res(:)
                real                          :: smpd_l, box_l, fsc05, fsc0143
                integer                       :: ifsc, fsc05_crossed_bin, fsc0143_crossed_bin
                logical                       :: fsc05_crossed, fsc0143_crossed
                if(.not. vol_oris%get_noris() .eq. fsc_oris%get_noris()) return
                if(.not. vol_oris%isthere(iori_l, "smpd")) return
                if(.not. fsc_oris%isthere(iori_l, "fsc")) return
                if(.not. fsc_oris%isthere(iori_l, "box")) return
                call fsc_oris%getter(iori_l, "fsc", fscfile)
                smpd_l = vol_oris%get(iori_l, "smpd")
                box_l  = fsc_oris%get(iori_l, "box")
                if(.not. file_exists(fscfile)) THROW_HARD('fsc file doesnt exist; print_segment_json')
                fsc = file2rarr(fscfile)
                res = get_resarr(int(box_l), smpd_l)
                call json%create_object(fsc_json, 'fsc')
                call json%add(fsc_json, 'type', "plot_bar")
                call json%create_array(datasets, "datasets")
                call json%create_array(data,     "data")
                call json%create_array(labels,   "labels")
                call json%create_object(dataset, "dataset")
                fsc05_crossed   = .false.
                fsc0143_crossed = .false.
                do ifsc=1, size(fsc)
                    if(.not. fsc05_crossed) then
                        if(fsc(ifsc) .gt. 0.5) then
                            fsc05 = res(ifsc)
                        else
                            fsc05_crossed = .true.
                            fsc05_crossed_bin = ifsc
                        end if
                    end if
                    if(.not. fsc0143_crossed) then
                        if(fsc(ifsc) .gt. 0.143) then
                            fsc0143 = res(ifsc)
                        else
                            fsc0143_crossed = .true.
                            fsc0143_crossed_bin = ifsc
                        end if
                    end if
                    call json%add(data,   '', dble(fsc(ifsc)))
                    call json%add(labels, '', dble(res(ifsc)))
                end do
                call json%add(dataset, 'borderColor', "rgba(30, 144, 255, 0.5)")
                call json%add(dataset, 'pointStyle', .false.)
                call json%add(dataset, 'cubicInterpolationMode', 'monotone')
                call json%add(dataset, 'tension', dble(0.4))
                call json%add(dataset, data)
                call json%add(datasets, dataset)
                call json%add(fsc_json, datasets)
                call json%add(fsc_json, labels)
                call json%add(json_ori, fsc_json)
                call json%add(json_ori, 'fsc05',   dble(fsc05))
                call json%add(json_ori, 'fsc0143', dble(fsc0143))
                call json%add(json_ori, 'fsc05bin',   fsc05_crossed_bin)
                call json%add(json_ori, 'fsc0143bin', fsc0143_crossed_bin)
                call fscfile%kill
            end subroutine add_fsc

            subroutine add_oriplot( iori_l )
                integer, intent(in)       :: iori_l
                type(json_value), pointer :: oriplot_json, datasets, dataset, data, xy
                integer, allocatable      :: projection_counts(:,:)
                integer                   :: iptcl, state, proj, iproj, idataset, max(2)
                if(.not. self%os_ptcl3D%isthere("state")) return
                allocate(projection_counts(size(projections, 1), 2))
                projection_counts = 0
                do iptcl=1, self%os_ptcl3D%get_noris()
                    state = self%os_ptcl3D%get_int(iptcl, "state")
                    proj  = self%os_ptcl3D%get_int(iptcl, "proj")
                    if(state .ne. iori_l) cycle
                    projection_counts(proj, 1) = projection_counts(proj, 1) + 1
                end do
                max = maxval(projection_counts, 1)
                do iproj = 1, size(projection_counts, 1)
                    if(projection_counts(iproj, 1) .gt. ceiling(max(1)/10.0)) then
                        projection_counts(iproj, 2) = 1
                    else if (projection_counts(iproj, 1) .gt. ceiling(max(1)/100.0)) then
                        projection_counts(iproj, 2) = 2
                    else if (projection_counts(iproj, 1) .gt. ceiling(max(1)/1000.0)) then
                        projection_counts(iproj, 2) = 3
                    else
                        projection_counts(iproj, 2) = 4
                    end if
                end do
                call json%create_object(oriplot_json, 'orientations')
                call json%add(oriplot_json, 'type', "plot_scatter")
                call json%create_array(datasets, 'datasets')
                do idataset = 1, 4
                    call json%create_object(dataset, 'dataset')
                    call json%create_array(data, 'data')
                    call json%add(dataset, 'label', 'logfold population ' // int2str(idataset))
                    if(idataset .eq. 1) call json%add(dataset, 'backgroundColor', "rgb(255, 99,  71 )")
                    if(idataset .eq. 2) call json%add(dataset, 'backgroundColor', "rgb(255, 215, 0  )")
                    if(idataset .eq. 3) call json%add(dataset, 'backgroundColor', "rgb(60,  179, 113)")
                    if(idataset .eq. 4) call json%add(dataset, 'backgroundColor', "rgb(30,  144, 255)")
                    do iproj = 1, size(projection_counts, 1)
                        if(projection_counts(iproj, 2) .ne. idataset) cycle
                        if(projection_counts(iproj, 1) .eq. 0) cycle
                        call json%create_object(xy, 'xy')
                        call json%add(xy, 'x', dble(projections(iproj, 2)))
                        call json%add(xy, 'y', dble(projections(iproj, 3)))
                        call json%add(data, xy)
                    end do
                    call json%add(dataset,      data)
                    call json%add(datasets,     dataset)
                end do
                call json%add(oriplot_json, datasets)
                call json%add(json_ori,     oriplot_json)
                if(allocated(projection_counts)) deallocate(projection_counts)
            end subroutine add_oriplot

            function sort_oris( seg_oris ) result( arr )
                type(oris), intent(in)  :: seg_oris
                integer,    allocatable :: arr(:)
                real,       allocatable :: sort_vals(:)
                if(.not. seg_oris%isthere(sort_key)) THROW_HARD('invalid sort key; print_segment_json')
                if(seg_oris%ischar(1, sort_key))     THROW_HARD('sort values are characters; print_segment_json')
                sort_vals = seg_oris%get_all(sort_key)
                allocate(arr(size(sort_vals)))
                do iori=1,size(sort_vals)
                    arr(iori) = iori
                end do
                call hpsort(sort_vals, arr)
                deallocate(sort_vals)
            end function sort_oris

            subroutine get_projections( noris_l )
                integer, intent(in) :: noris_l
                real                :: minproj, maxproj, e1, e2
                integer             :: iproj, iptcl, proj, state
                if(.not. self%os_ptcl3D%isthere("state")) return
                if(.not. self%os_ptcl3D%isthere("proj"))  return
                if(.not. self%os_ptcl3D%isthere("e1"))    return
                if(.not. self%os_ptcl3D%isthere("e2"))    return
                call self%os_ptcl3D%minmax("proj", minproj, maxproj)
                allocate(projections(int(maxproj), 3))
                do iproj=1, int(maxproj)
                    projections(iproj, 1) = 0.0 !active
                    projections(iproj, 2) = 0.0 !e1
                    projections(iproj, 3) = 0.0 !e2
                end do
                do iptcl=1, self%os_ptcl3D%get_noris()
                    state = self%os_ptcl3D%get_int(iptcl, "state")
                    proj  = self%os_ptcl3D%get_int(iptcl, "proj")
                    if(state .eq. 0) cycle
                    if(proj  .eq. 0) cycle
                    if(projections(proj, 1) .eq. 0.0) then
                        e1 = self%os_ptcl3D%get(iptcl, "e1")
                        e2 = self%os_ptcl3D%get(iptcl, "e2")
                        projections(proj, 1) = 1.0
                        projections(proj, 2) = e1
                        projections(proj, 3) = e2
                    end if
                end do
            end subroutine get_projections

    end subroutine print_segment_json

    ! Readers

    module subroutine read( self, fname, wthreads )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        logical, optional, intent(in)    :: wthreads
        type(string) :: projfile
        integer :: isegment
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('format of: '//fname%to_char()//' not supported; read')
        endif
        projfile = fname
        if( .not. file_exists(projfile) )then
            THROW_HARD('file: '// projfile%to_char() //' does not exist; read')
        endif
        call self%bos%open(projfile)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment, wthreads=wthreads)
        end do
        call self%bos%close
        call self%update_projinfo(fname)
        call self%write_segment_inside('projinfo', fname)
    end subroutine read

    module subroutine read_non_data_segments( self, fname )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        type(string) :: projfile
        integer :: iseg
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//fname%to_char()//' not supported; read_non_data_segments')
        endif
        projfile = fname
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            do iseg=11,MAXN_OS_SEG
                call self%segreader(iseg)
            end do
            call self%bos%close
        else
            THROW_HARD('projfile: '// projfile%to_char() //' nonexistent; read_non_data_segments')
        endif
    end subroutine read_non_data_segments

    module subroutine read_ctfparams_state_eo( self, fname )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        integer :: isegment
        if( .not. file_exists(fname) )then
            THROW_HARD('inputted file: '//fname%to_char()//' does not exist; read_ctfparams_state_eo')
        endif
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//fname%to_char()//' not supported; read_ctfparams_state_eo')
        endif
        call self%bos%open(fname)
        do isegment=1,self%bos%get_n_segments()
            call self%segreader(isegment, only_ctfparams_state_eo=.true.)
        end do
        call self%bos%close
    end subroutine read_ctfparams_state_eo

    module subroutine read_mic_stk_ptcl2D_segments( self, fname, wthreads )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        logical, optional, intent(in)    :: wthreads
        type(string) :: projfile
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('format of: '//fname%to_char()//' not supported; read')
        endif
        projfile = fname
        if( .not. file_exists(projfile) )then
            THROW_HARD('file: '// projfile%to_char()//' does not exist; read')
        endif
        call self%bos%open(projfile)
        call self%segreader(MIC_SEG, wthreads=wthreads)
        call self%segreader(STK_SEG, wthreads=wthreads)
        call self%segreader(PTCL2D_SEG, wthreads=wthreads)
        call self%bos%close
    end subroutine read_mic_stk_ptcl2D_segments

    module subroutine read_segment( self, oritype, fname, fromto, wthreads )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        class(string),     intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: wthreads
        integer :: isegment
        if( .not. file_exists(fname) )then
            THROW_HARD('inputted file: '//fname%to_char()//' does not exist; read_segment')
        endif
        select case(fname2format(fname))
            case('O')
                ! *.simple project file
                isegment = oritype2segment(oritype)
                call self%bos%open(fname)
                call self%segreader(isegment,fromto=fromto,wthreads=wthreads)
                call self%bos%close
            case('T')
                ! *.txt plain text ori file
                select case(trim(oritype))
                    case('mic')
                        call self%os_mic%read(fname)
                    case('stk')
                        call self%os_stk%read(fname)
                    case('ptcl2D')
                        call self%os_ptcl2D%read(fname, fromto=fromto)
                    case('cls2D')
                        call self%os_cls2D%read(fname)
                    case('cls3D')
                        call self%os_cls3D%read(fname,  fromto=fromto)
                    case('ptcl3D')
                        call self%os_ptcl3D%read(fname, fromto=fromto)
                    case('out')
                        call self%os_out%read(fname)
                    case('optics')
                        call self%os_optics%read(fname)
                    case('projinfo')
                        call self%projinfo%read(fname)
                    case('jobproc')
                        call self%jobproc%read(fname)
                    case('compenv')
                        call self%compenv%read(fname)
                    case DEFAULT
                        THROW_HARD('unsupported oritype flag; read_segment')
                end select
            case DEFAULT
                THROW_HARD('file format of: '//fname%to_char()//' not supported; read_segment')
        end select
    end subroutine read_segment

    module subroutine segreader( self, isegment, fromto, only_ctfparams_state_eo, wthreads )
        class(sp_project),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional,          intent(in)    :: fromto(2)
        logical, optional,          intent(in)    :: only_ctfparams_state_eo, wthreads
        integer :: n
        n = self%bos%get_n_records(isegment)
        select case(isegment)
            case(MIC_SEG)
                call self%os_mic%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_mic,    wthreads=wthreads)
            case(STK_SEG)
                call self%os_stk%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_stk,    only_ctfparams_state_eo=only_ctfparams_state_eo, wthreads=wthreads)
            case(PTCL2D_SEG)
                call self%os_ptcl2D%new(n, is_ptcl=.true.)
                call self%bos%read_segment(isegment, self%os_ptcl2D, only_ctfparams_state_eo=only_ctfparams_state_eo, wthreads=wthreads)
            case(CLS2D_SEG)
                call self%os_cls2D%new(n,  is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%os_cls3D%new(n,  is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_cls3D)
            case(PTCL3D_SEG)
                call self%os_ptcl3D%new(n, is_ptcl=.true.)
                call self%bos%read_segment(isegment, self%os_ptcl3D, only_ctfparams_state_eo=only_ctfparams_state_eo, wthreads=wthreads)
            case(OUT_SEG)
                call self%os_out%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_out)
            case(OPTICS_SEG)
                call self%os_optics%new(n,    is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%os_optics)
            case(PROJINFO_SEG)
                call self%projinfo%new(n,  is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%jobproc%new(n,   is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%compenv%new(n,   is_ptcl=.false.)
                call self%bos%read_segment(isegment, self%compenv)
        end select
    end subroutine segreader

    module subroutine read_segments_info( self, fname, seginds, seginfos )
        class(sp_project),                  intent(inout) :: self
        class(string),                      intent(in)  :: fname
        type(binoris_seginfo), allocatable, intent(out) :: seginfos(:)
        integer,               allocatable, intent(out) :: seginds(:)
        type(string) :: projfile
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//fname%to_char()//' not supported; sp_project :: read_segments_info')
        endif
        projfile = fname
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            call self%bos%get_segments_info(seginds, seginfos)
            call self%bos%close
        else
            THROW_HARD('projfile: '//projfile%to_char()//' nonexistent; print_info')
        endif
    end subroutine read_segments_info

    !>  Convenience funtion for checking # of movies/mics, stacks and particles
    module subroutine read_data_info( self, fname, nmics, nstks, nptcls )
        class(sp_project),   intent(inout) :: self
        class(string),       intent(in)    :: fname
        integer,             intent(out)   :: nmics, nstks, nptcls
        type(binoris_seginfo), allocatable :: seginfos(:)
        integer,               allocatable :: seginds(:)
        integer :: i,n2D, n3D
        call self%read_segments_info(fname, seginds, seginfos)
        nmics  = 0
        nstks  = 0
        n2D    = 0
        n3D    = 0
        do i = 1,size(seginds)
            select case(seginds(i))
            case(MIC_SEG)
                nmics = int(seginfos(i)%n_records)
            case(STK_SEG)
                nstks = int(seginfos(i)%n_records)
            case(PTCL2D_SEG)
                n2D = int(seginfos(i)%n_records)
            case(PTCL3D_SEG)
                n3D = int(seginfos(i)%n_records)
            end select
        enddo
        ! if( n2D /= n3D )then
        !     THROW_WARN('Inconsistent # of particles in the 2D/3D segments; read_data_info: '//trim(fname))
        ! endif
        nptcls = n2D
    end subroutine read_data_info

    ! Writers

    module subroutine write( self, fname, fromto, isegment, tempfile )
        class(sp_project),                    intent(inout) :: self
        class(string),              optional, intent(in)    :: fname
        integer,                    optional, intent(in)    :: fromto(2)
        integer(kind(ENUM_ORISEG)), optional, intent(in)    :: isegment
        logical,                    optional, intent(in)    :: tempfile
        type(string) :: projfile, tmpfile
        integer(kind(ENUM_ORISEG))    :: iseg
        logical :: l_tmp
        l_tmp = .false.
        if( present(tempfile) ) l_tmp = tempfile
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                THROW_HARD('file format of: '//fname%to_char()//' not supported; write')
            endif
            projfile = fname
        else
            call self%projinfo%getter(1, 'projfile', projfile)
        endif
        if( l_tmp )then
            tmpfile = swap_suffix(projfile, string(METADATA_EXT), string('.tmp'))
            call self%bos%open(tmpfile, del_if_exists=.true.)
        else
            call self%bos%open(projfile, del_if_exists=.true.)
        endif
        if( present(isegment) )then
            call self%segwriter(isegment, fromto)
        else
            do iseg=1,MAXN_OS_SEG
                call self%segwriter(iseg, fromto)
            end do
        endif
        ! update header
        call self%bos%write_header
        call self%bos%close
        if( l_tmp ) call simple_rename(tmpfile, projfile)
    end subroutine write

    module subroutine write_segment_inside( self, oritype, fname, fromto )
        class(sp_project),          intent(inout) :: self
        character(len=*),           intent(in)    :: oritype
        class(string),    optional, intent(in)    :: fname
        integer,          optional, intent(in)    :: fromto(2)
        type(string) :: projfile
        integer(kind(ENUM_ORISEG)) :: iseg
        if( present(fname) )then
            if( fname2format(fname) .ne. 'O' )then
                THROW_HARD('file format of: '//fname%to_char()//' not supported; sp_project :: write')
            endif
            projfile = fname
        else
            call self%projinfo%getter(1, 'projfile', projfile)
        endif
        if( file_exists(projfile) )then
            iseg = oritype2segment(oritype)
            call self%bos%open(projfile)
            call self%segwriter_inside(iseg, fromto)
        else
            call self%write(fname, fromto)
        endif
        ! no need to update header (taken care of in binoris object)
        call self%bos%close
    end subroutine write_segment_inside

    module subroutine write_non_data_segments( self, fname )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        type(string) :: projfile
        integer :: iseg
        if( fname2format(fname) .ne. 'O' )then
            THROW_HARD('file format of: '//fname%to_char()//' not supported; sp_project :: write_non_data_segments')
        endif
        projfile = fname
        if( file_exists(projfile) )then
            call self%bos%open(projfile)
            do iseg=11,MAXN_OS_SEG
                call self%segwriter(iseg)
            end do
            ! update header
            call self%bos%write_header
            call self%bos%close
        else
            THROW_HARD('projfile: '//projfile%to_char()//' nonexistent; write_non_data_segments')
        endif
    end subroutine write_non_data_segments

    module subroutine write_segment2txt( self, oritype, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        class(string),     intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        select case(fname2format(fname))
            case('O')
                THROW_HARD('write_segment2txt is not supported for *.simple project files; write_segment2txt')
            case('T')
                ! *.txt plain text ori file
                select case(trim(oritype))
                    case('mic')
                        if( self%os_mic%get_noris() > 0 )then
                            call self%os_mic%write(fname)
                        else
                            THROW_WARN('no mic-type oris available to write; write_segment2txt')
                        endif
                    case('stk')
                        if( self%os_stk%get_noris() > 0 )then
                            call self%os_stk%write(fname)
                        else
                            THROW_WARN('no stk-type oris available to write; write_segment2txt')
                        endif
                    case('ptcl2D')
                        if( self%os_ptcl2D%get_noris() > 0 )then
                            call self%os_ptcl2D%write(fname, fromto)
                        else
                            THROW_WARN('no ptcl2D-type oris available to write; write_segment2txt')
                        endif
                    case('cls2D')
                        if( self%os_cls2D%get_noris() > 0 )then
                            call self%os_cls2D%write(fname)
                        else
                            THROW_WARN('no cls2D-type oris available to write; write_segment2txt')
                        endif
                    case('cls3D')
                        if( self%os_cls3D%get_noris() > 0 )then
                            call self%os_cls3D%write(fname,  fromto)
                        else
                            THROW_WARN('no cls3D-type oris available to write; write_segment2txt')
                        endif
                    case('ptcl3D')
                        if( self%os_ptcl3D%get_noris() > 0 )then
                            call self%os_ptcl3D%write(fname, fromto)
                        else
                            THROW_WARN('no ptcl3D-type oris available to write; write_segment2txt')
                        endif
                    case('out')
                        if( self%os_out%get_noris() > 0 )then
                            call self%os_out%write(fname)
                        else
                            THROW_WARN('no out-type oris available to write; write_segment2txt')
                        endif
                    case('optics')
                        if( self%os_optics%get_noris() > 0 )then
                            call self%os_optics%write(fname)
                        else
                            THROW_WARN('no optics-type oris available to write; write_segment2txt')
                        endif
                    case('projinfo')
                        if( self%projinfo%get_noris() > 0 )then
                            call self%projinfo%write(fname, fromto)
                        else
                            THROW_WARN('no projinfo-type oris available to write; write_segment2txt')
                        endif
                    case('jobproc')
                        if( self%jobproc%get_noris() > 0 )then
                            call self%jobproc%write(fname)
                        else
                            THROW_WARN('no jobproc-type oris available to write; write_segment2txt')
                        endif
                    case('compenv')
                        if( self%compenv%get_noris() > 0 )then
                            call self%compenv%write(fname)
                        else
                            THROW_WARN('no compenv-type oris available to write; write_segment2txt')
                        endif
                    case DEFAULT
                        THROW_HARD('unsupported oritype flag; write_segment2txt')
                end select
            case DEFAULT
                THROW_HARD('file format of: '//fname%to_char()//'not supported; write_segment2txt')
        end select
    end subroutine write_segment2txt

    module subroutine segwriter( self, isegment, fromto )
        class(sp_project), intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional, intent(in)    :: fromto(2)
        select case(isegment)
            case(MIC_SEG)
                call self%bos%write_segment(isegment, self%os_mic, fromto)
            case(STK_SEG)
                call self%bos%write_segment(isegment, self%os_stk)
            case(PTCL2D_SEG)
                call self%bos%write_segment(isegment, self%os_ptcl2D, fromto)
            case(CLS2D_SEG)
                call self%bos%write_segment(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%bos%write_segment(isegment, self%os_cls3D, fromto)
            case(PTCL3D_SEG)
                call self%bos%write_segment(isegment, self%os_ptcl3D, fromto)
            case(OUT_SEG)
                call self%bos%write_segment(isegment, self%os_out)
            case(OPTICS_SEG)
                call self%bos%write_segment(isegment, self%os_optics)
            case(PROJINFO_SEG)
                call self%bos%write_segment(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%bos%write_segment(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%bos%write_segment(isegment, self%compenv)
        end select
    end subroutine segwriter

    module subroutine segwriter_inside( self, isegment, fromto )
        class(sp_project),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional,          intent(in)    :: fromto(2)
        select case(isegment)
            case(MIC_SEG)
                call self%bos%write_segment_inside(isegment, self%os_mic, fromto)
            case(STK_SEG)
                call self%bos%write_segment_inside(isegment, self%os_stk)
            case(PTCL2D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_ptcl2D, fromto)
            case(CLS2D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_cls2D)
            case(CLS3D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_cls3D, fromto)
            case(PTCL3D_SEG)
                call self%bos%write_segment_inside(isegment, self%os_ptcl3D, fromto)
            case(OUT_SEG)
                call self%bos%write_segment_inside(isegment, self%os_out)
            case(OPTICS_SEG)
                call self%bos%write_segment_inside(isegment, self%os_optics)
            case(PROJINFO_SEG)
                call self%bos%write_segment_inside(isegment, self%projinfo)
            case(JOBPROC_SEG)
                call self%bos%write_segment_inside(isegment, self%jobproc)
            case(COMPENV_SEG)
                call self%bos%write_segment_inside(isegment, self%compenv)
        end select
    end subroutine segwriter_inside

    module subroutine write_mics_star( self, fname )
        class(sp_project),       intent(inout) :: self
        class(string), optional, intent(in)    :: fname
        type(starfile) :: star
        type(string)   :: l_fname
        if( self%os_mic%get_noris() == 0 ) return
        if(present(fname)) then 
            l_fname = fname
        else
            l_fname = MICS_STAR_BODY // STAR_EXT
        end if
        call star%init(l_fname, verbose=.true.)
        call star%write_optics_table(self%os_optics)
        call star%write_mics_table(self%os_mic)
        call star%complete()
    end subroutine write_mics_star

    module subroutine write_ptcl2D_star( self, fname )
        class(sp_project),       intent(inout) :: self
        class(string), optional, intent(in)    :: fname
        type(starfile) :: star
        type(string)   :: l_fname
        if( self%os_mic%get_noris() == 0 ) return
        if(present(fname)) then 
            l_fname = fname
        else
            l_fname = PTCL2D_STAR_BODY // STAR_EXT
        end if
        call star%init(l_fname, verbose=.true.)
        call star%write_optics_table(self%os_optics)
        call star%write_ptcl2D_table(self%os_ptcl2D, self%os_stk, mics_oris=self%os_mic)
        call star%complete()
    end subroutine write_ptcl2D_star

    module subroutine write_optics_map( self, fname_prefix )
        class(sp_project),          intent(inout) :: self
        character(len=*),           intent(in)    :: fname_prefix
        type(sp_project) :: spproj_optics
        type(nrtxtfile)  :: map_file
        real             :: mapline(2)
        integer          :: imic
        call map_file%new(string(fname_prefix//TXT_EXT), 2, 2)
        do imic=1, self%os_mic%get_noris()
            if(self%os_mic%isthere(imic, 'importind') .and. self%os_mic%isthere(imic, 'ogid')) then
               mapline(1) = self%os_mic%get_int(imic, 'importind')
               mapline(2) = self%os_mic%get_int(imic, 'ogid')
               call map_file%write(mapline)
            end if
        end do
        call map_file%kill
        call spproj_optics%os_optics%copy(self%os_optics, is_ptcl=.false.)
        call spproj_optics%write(string(fname_prefix//METADATA_EXT))
        call spproj_optics%kill
    end subroutine write_optics_map

    !------ Private Non-type-bound helpers ------

    module function segment2info( iseg, n_records ) result( info )
        integer(kind(ENUM_ORISEG)), intent(in) :: iseg
        integer(kind=8),            intent(in) :: n_records
        character(len=:), allocatable :: info, nrecs_str
        nrecs_str = int2str(int(n_records,kind=4))
        info = nrecs_str//' record(s) of '
        select case(iseg)
            case(MIC_SEG)
                info = info//'movie and micrograph (integrated movie) info, one per movie/micrograph'
            case(STK_SEG)
                info = info//'stack (extracted particles) info, one per stack of particles'
            case(PTCL2D_SEG)
                info = info//'2D information generated by cluster2D, one per particle'
            case(CLS2D_SEG)
                info = info//'data generated by cluster2D, one per 2D cluster'
            case(CLS3D_SEG)
                info = info//'3D information for class averages, one per class'
            case(PTCL3D_SEG)
                info = info//'3D information, one per particle'
            case(OUT_SEG)
                info = info//'critical project outputs: class averages, 3D volumes, FSC/FRC files etc.'
            case(OPTICS_SEG)
                info = info//'optics group information.'
            case(PROJINFO_SEG)
                info = info//'information about the project, project name etc.'
            case(JOBPROC_SEG)
                info = info//'all command-lines executed throughout the project'
            case(COMPENV_SEG)
                info = info//'computing environment specifications, queue system, memory per job etc.'
            case DEFAULT
                write(logfhandle,*) 'iseg: ', iseg
                THROW_HARD('unsupported segment of kind(ENUM_ORISEG); segment2oritype')
        end select
    end function segment2info

    module function segment2oritype( iseg ) result( oritype )
        integer(kind(ENUM_ORISEG)), intent(in) :: iseg
        character(len=:), allocatable :: oritype
        select case(iseg)
            case(MIC_SEG)
                oritype = 'mic'
            case(STK_SEG)
                oritype = 'stk'
            case(PTCL2D_SEG)
                oritype = 'ptcl2D'
            case(CLS2D_SEG)
                oritype = 'cls2D'
            case(CLS3D_SEG)
                oritype = 'cls3D'
            case(PTCL3D_SEG)
                oritype = 'ptcl3D'
            case(OUT_SEG)
                oritype = 'out'
            case(OPTICS_SEG)
                oritype = 'optics'
            case(PROJINFO_SEG)
                oritype = 'projinfo'
            case(JOBPROC_SEG)
                oritype = 'jobproc'
            case(COMPENV_SEG)
                oritype = 'compenv'
            case DEFAULT
                write(logfhandle,*) 'iseg: ', iseg
                THROW_HARD('unsupported segment of kind(ENUM_ORISEG); segment2oritype')
        end select
    end function segment2oritype

end submodule simple_sp_project_io

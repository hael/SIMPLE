module simple_stream_utils
include 'simple_lib.f08'
use simple_cmdline,             only: cmdline
use simple_default_clines,      only: set_automask2D_defaults
use simple_image,               only: image
use simple_parameters,          only: parameters, params_glob
use simple_stack_io,            only: stack_io
use simple_stream_communicator, only: stream_http_communicator 
use simple_gui_utils
use simple_nice
use simple_qsys_funs
implicit none
#include "simple_local_flags.inc"

contains

    !> To deal with dynamic user input diring streaming
    subroutine update_user_params( cline_here, httpcom )
        type(cmdline),                            intent(inout) :: cline_here
        type(stream_http_communicator), optional, intent(inout) :: httpcom
        type(oris) :: os
        character(kind=CK,len=:), allocatable :: interactive, ring
        real       :: tilt_thres, beamtilt, astigthreshold, ctfresthreshold, icefracthreshold
        real(kind=dp) :: icefracthreshold_dp
        real       :: moldiam_refine
        integer    :: moldiam_refine_int, astigthreshold_int, ctfresthreshold_int, moldiam_int, moldiam_ring_int
        logical    :: found
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS) )then
            call os%read(string(USER_PARAMS))
            if( os%isthere(1,'tilt_thres') ) then
                tilt_thres = os%get(1,'tilt_thres')
                if( abs(tilt_thres-params_glob%tilt_thres) > 0.001) then
                     if(tilt_thres < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO LOW: ',tilt_thres
                     else if(tilt_thres > 1) then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO HIGH: ',tilt_thres
                     else
                         params_glob%tilt_thres = tilt_thres
                         params_glob%updated    = 'yes'
                         call cline_here%set('tilt_thres', params_glob%tilt_thres)
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES UPDATED TO: ',tilt_thres
                     endif
                endif
            endif
            if( os%isthere(1,'beamtilt') ) then
                beamtilt = os%get(1,'beamtilt')
                if( beamtilt .eq. 1.0 ) then
                    params_glob%beamtilt = 'yes'
                    params_glob%updated  = 'yes'
                    call cline_here%set('beamtilt', params_glob%beamtilt)
                    write(logfhandle,'(A)')'>>> OPTICS ASSIGNMENT UDPATED TO USE BEAMTILT'
                else if( beamtilt .eq. 0.0 ) then
                    params_glob%beamtilt = 'no'
                    params_glob%updated  = 'yes'
                    call cline_here%set('beamtilt', params_glob%beamtilt)
                    write(logfhandle,'(A)')'>>> OPTICS ASSIGNMENT UDPATED TO IGNORE BEAMTILT'
                else
                    write(logfhandle,'(A,F8.2)')'>>> OPTICS UPDATE INVALID BEAMTILT VALUE: ',beamtilt
                endif
            endif
            if( os%isthere(1,'astigthreshold') ) then
                astigthreshold = os%get(1,'astigthreshold')
                if( abs(astigthreshold-params_glob%astigthreshold) > 0.001) then
                     if(astigthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD TOO LOW: ',astigthreshold
                     else if(astigthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM TOO HIGH: ',astigthreshold
                     else
                         params_glob%astigthreshold = astigthreshold
                         params_glob%updated    = 'yes'
                         call cline_here%set('astigthreshold', params_glob%astigthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD UPDATED TO: ',astigthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'ctfresthreshold') ) then
                ctfresthreshold = os%get(1,'ctfresthreshold')
                if( abs(ctfresthreshold-params_glob%ctfresthreshold) > 0.001) then
                     if(ctfresthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD TOO LOW: ',ctfresthreshold
                     else if(ctfresthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION TOO HIGH: ',ctfresthreshold
                     else
                         params_glob%ctfresthreshold = ctfresthreshold
                         params_glob%updated    = 'yes'
                         call cline_here%set('ctfresthreshold', params_glob%ctfresthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD UPDATED TO: ',ctfresthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'icefracthreshold') ) then
                icefracthreshold = os%get(1,'icefracthreshold')
                if( abs(icefracthreshold-params_glob%icefracthreshold) > 0.001) then
                     if(icefracthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD TOO LOW: ',icefracthreshold
                     else if(icefracthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION TOO HIGH: ',icefracthreshold
                     else
                         params_glob%icefracthreshold = icefracthreshold
                         params_glob%updated    = 'yes'
                         call cline_here%set('icefracthreshold', params_glob%icefracthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD UPDATED TO: ',icefracthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'moldiam_refine') ) then
                moldiam_refine = os%get(1,'moldiam_refine')
                if( abs(moldiam_refine-params_glob%moldiam_refine) > 0.001) then
                     if(moldiam_refine < 10)then
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO LOW: ' , moldiam_refine
                     else if(moldiam_refine > 1000) then
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO HIGH: ', moldiam_refine
                     else
                         params_glob%moldiam_refine = moldiam_refine
                         params_glob%updated        = 'yes'
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER UPDATED TO: ', moldiam_refine
                     endif
                endif
            endif
            call del_file(USER_PARAMS)
        endif
        call os%kill
        ! nice
        if(present(httpcom)) then
            if( httpcom%arg_associated() )then
                ! moldiam_refine
                call httpcom%get_json_arg('moldiam_refine', moldiam_refine_int, found)
                if(found) then
                    if( abs(real(moldiam_refine_int)-params_glob%moldiam_refine) > 0.001) then
                        if(moldiam_refine_int < 10 .and. moldiam_refine_int > 0)then
                            write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO LOW: ' , real(moldiam_refine_int)
                        else if(moldiam_refine_int > 1000) then
                            write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO HIGH: ', real(moldiam_refine_int)
                        else
                            params_glob%moldiam_refine = real(moldiam_refine_int)
                            params_glob%updated        = 'yes'
                            write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER UPDATED TO: ', real(moldiam_refine_int)
                        end if
                    end if
                end if
                ! moldiam
                call httpcom%get_json_arg('moldiam', moldiam_int, found)
                if(found) then
                    if( abs(real(moldiam_int)-params_glob%moldiam) > 0.001) then
                        if(moldiam_int < 20 .and. moldiam_int > 0)then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR DIAMETER TOO LOW: ' , real(moldiam_int)
                        else if(moldiam_int > 1000) then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR DIAMETER TOO HIGH: ', real(moldiam_int)
                        else
                            params_glob%moldiam = real(moldiam_int)
                            params_glob%updated        = 'yes'
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR DIAMETER UPDATED TO: ', real(moldiam_int)
                        end if
                    end if
                end if
                ! moldiam_ring
                call httpcom%get_json_arg('moldiam_ring', moldiam_ring_int, found)
                if(found) then
                    if( abs(real(moldiam_ring_int)-params_glob%moldiam_ring) > 0.001) then
                        if(moldiam_ring_int < 20 .and. moldiam_ring_int > 0)then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR RING DIAMETER TOO LOW: ' , real(moldiam_ring_int)
                        else if(moldiam_ring_int > 1000) then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR RING DIAMETER TOO HIGH: ', real(moldiam_ring_int)
                        else
                            params_glob%moldiam_ring = real(moldiam_ring_int)
                            params_glob%updated      = 'yes'
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR RING DIAMETER UPDATED TO: ', real(moldiam_ring_int)
                        end if
                    end if
                end if
                ! interactive
                call httpcom%get_json_arg('interactive', interactive, found)
                if(found) then
                    params_glob%interactive = interactive
                    params_glob%updated        = 'yes'
                    write(logfhandle,'(A,A)')'>>> INTERACTIVE UPDATED TO: ', interactive
                end if
                ! ring
                call httpcom%get_json_arg('ring', ring, found)
                if(found) then
                    params_glob%ring = ring
                    params_glob%updated        = 'yes'
                    write(logfhandle,'(A,A)')'>>> RING UPDATED TO: ', ring
                end if
                ! astigthreshold
                call httpcom%get_json_arg('astigthreshold', astigthreshold_int, found)
                if(found) then
                    if( abs(real(astigthreshold_int)-params_glob%astigthreshold) > 0.001) then
                        params_glob%astigthreshold = real(astigthreshold_int)
                        params_glob%updated        = 'yes'
                        write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD UPDATED TO: ', real(astigthreshold_int)
                    end if
                end if
                ! ctfresthreshold
                call httpcom%get_json_arg('ctfresthreshold', ctfresthreshold_int, found)
                if(found) then
                    if( abs(real(ctfresthreshold_int)-params_glob%ctfresthreshold) > 0.001) then
                        params_glob%ctfresthreshold = real(ctfresthreshold_int)
                        params_glob%updated        = 'yes'
                        write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD UPDATED TO: ', real(ctfresthreshold_int)
                    end if
                end if
                ! icefracthreshold
                call httpcom%get_json_arg('icefracthreshold', icefracthreshold_dp, found)
                if(found) then
                    if( abs(icefracthreshold_dp-params_glob%icefracthreshold) > 0.001) then
                        params_glob%icefracthreshold = real(icefracthreshold_dp)
                        params_glob%updated        = 'yes'
                        write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD UPDATED TO: ', icefracthreshold_dp
                    end if
                end if
            end if
            call httpcom%destroy_arg
        end if
    end subroutine update_user_params

    subroutine wait_for_folder( httpcom, folder, stopmsg )
        use simple_stream_communicator, only:stream_http_communicator
        class(stream_http_communicator), intent(inout) :: httpcom
        class(string),                   intent(in)    :: folder
        character(len=*),                intent(in)    :: stopmsg
        integer :: i
        if(.not. dir_exists(folder)) then
            write(logfhandle, *) ">>> WAITING FOR ", folder%to_char(), " TO BE GENERATED"
            ! wait up to 24 hours
            do i = 1, 8640
                if(dir_exists(folder)) then
                    write(logfhandle, *) ">>> ", folder%to_char(), " FOUND"
                    exit
                endif
                call sleep(10)
                call httpcom%send_jobstats()
                if( httpcom%exit_status() )then
                    ! termination
                    write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                    call httpcom%term()
                    call simple_end(trim(stopmsg))
                    call EXIT(0)
                endif
            end do
        endif
    end subroutine wait_for_folder

    ! Class rejection routine based on image moments & Total Variation Distance
    subroutine class_rejection( os, mask, adjust )
        class(oris),    intent(in)    :: os
        logical,        intent(inout) :: mask(:)
        real, optional, intent(in)    :: adjust
        real,    allocatable :: vals(:), x(:)
        logical, allocatable :: msk(:)
        real    :: eff_mean_thresh, eff_rel_var_thresh, eff_abs_var_thresh
        real    :: eff_tvd_thresh, eff_min_thresh, eff_max_thresh
        integer :: icls, i, n
        logical :: has_mean, has_var, has_tvd, has_minmax
        n   = os%get_noris()
        if( size(mask) /= n )THROW_HARD('Incompatible sizes! class rejection')
        msk = mask
        if( os%isthere('pop') )then
            do icls=1,n
                msk(icls) = os%get(icls,'pop') > 0.5
            enddo
        endif
        mask = msk
        if( count(msk) <= 5 )then
            deallocate(msk)
            return
        endif
        ! Effective threshold
        eff_mean_thresh    = params_glob%stream_mean_threshold
        eff_rel_var_thresh = params_glob%stream_rel_var_threshold
        eff_abs_var_thresh = params_glob%stream_abs_var_threshold
        eff_tvd_thresh     = max(0.001,min(0.999,params_glob%stream_tvd_theshold))
        eff_min_thresh     = -params_glob%stream_minmax_threshold
        eff_max_thresh     =  params_glob%stream_minmax_threshold
        if( present(adjust) )then
            eff_mean_thresh    = adjust * eff_mean_thresh
            eff_rel_var_thresh = adjust * eff_rel_var_thresh
            eff_abs_var_thresh = adjust * eff_abs_var_thresh
            eff_tvd_thresh     = min(0.999, adjust * eff_tvd_thresh)
            eff_min_thresh     = adjust * eff_min_thresh
            eff_max_thresh     = adjust * eff_max_thresh
        endif
        ! selection
        has_mean   = os%isthere('mean')
        has_var    = os%isthere('var')
        has_tvd    = os%isthere('tvd')
        has_minmax = os%isthere('min') .and. os%isthere('max')
        ! Mean
        if( has_mean )then
            vals = os%get_all('mean')
            x    = pack(vals, mask=msk)
            call robust_scaling(x)
            i = 0
            do icls = 1,n
                if( msk(icls) )then
                    i = i+1
                    if( mask(icls) ) mask(icls) = x(i) > eff_mean_thresh
                endif
            enddo
        endif
        ! Variance
        if( has_var )then
            vals = os%get_all('var')
            x    = pack(vals, mask=msk)
            call robust_scaling(x)
            i = 0
            do icls = 1,n
                if( msk(icls) )then
                    i = i+1
                    if( mask(icls) ) mask(icls) = x(i)       < eff_rel_var_thresh
                    if( mask(icls) ) mask(icls) = vals(icls) < eff_abs_var_thresh
                endif
            enddo
        endif
        ! Total Variation Distance
        if( has_tvd )then
            vals = os%get_all('tvd')
            do icls = 1,n
                if( mask(icls) ) mask(icls) = vals(icls) < eff_tvd_thresh
            enddo
        endif
        ! Min/max
        if( has_minmax )then
            do icls = 1,n
                if( mask(icls) )then
                    if(  (os%get(icls,'min') < eff_min_thresh).and.&
                        &(os%get(icls,'max') > eff_max_thresh) )then
                        mask(icls) = .false.
                    endif
                endif
            enddo
        endif
        deallocate(msk)
        if(allocated(vals) ) deallocate(vals, x)
    end subroutine class_rejection

    function stream_datestr() 
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        character(16) :: stream_datestr
        integer,dimension(8) :: values
        ! using keyword arguments
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        write(stream_datestr, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') values(1), '/', values(2), '/', values(3), '_', values(5), ':', values(6)
    end function stream_datestr

    subroutine process_selected_refs( imgfile, smpd, selection, mskdiam, box_for_pick, box_for_extract, nxtiles, nytiles )
        use simple_image_msk, only: automask2d
        class(string),   intent(in)    :: imgfile
        real,            intent(in)    :: smpd
        integer,         intent(in)    :: selection(:)
        real,            intent(out)   :: mskdiam
        integer,         intent(out)   :: box_for_pick, box_for_extract
        integer,         intent(inout) :: nxtiles, nytiles
        type(parameters)               :: params
        type(cmdline)                  :: cline
        type(stack_io)                 :: stkio_r, stkio_w
        type(image),       allocatable :: cavgs(:)
        class(parameters), pointer     :: params_ptr
        real,              allocatable :: diams(:), shifts(:,:)
        logical,           parameter   :: DEBUG = .false.
        real    :: maxdiam, mskrad_in_pix, moldiam
        integer :: ldim(3), icls, ncls, nsel
        nsel = size(selection)
        if( nsel == 0 ) return
        params_ptr => params_glob ! for safe call to automask2D
        nullify(params_glob)
        write(logfhandle,'(A,I6,A)')'>>> USER SELECTED FROM POOL: ', nsel,' clusters'
        write(logfhandle,'(A,A)')'>>> WRITING SELECTED CLUSTERS TO: ', STREAM_SELECTED_REFS // STK_EXT
        ! set defaults
        call set_automask2D_defaults(cline)
        ! parse parameters
        call params%new(cline)
        call find_ldim_nptcls(imgfile, ldim, ncls)
        ldim(3) = 1
        params%msk  = real(ldim(1)/2) - COSMSKHALFWIDTH ! for automasking
        params%smpd = smpd                              ! for automasking
        if( DEBUG )then
            print *, 'imgfile:              ', imgfile%to_char()
            print *, 'file_exists(imgfile): ', file_exists(imgfile)
            print *, 'ldim:                 ', ldim(1), ldim(2), ldim(3)
            print *, 'ncls:                 ', ncls
            print *, 'params%msk:           ', params%msk
            print *, 'params%smpd:          ', params%smpd
            print *, 'nsel:                 ', nsel
        endif
        allocate( cavgs(nsel) )
        do icls = 1, nsel
            call cavgs(icls)%new([ldim(1),ldim(2),1], smpd)
        end do
        call stkio_r%open(imgfile, smpd, 'read', bufsz=ncls)
        call stkio_r%read_whole
        do icls = 1, nsel
            call stkio_r%get_image(selection(icls), cavgs(icls))
        end do
        call automask2D(cavgs, params%ngrow, nint(params%winsz), params%edge, diams, shifts)       
        box_for_pick    = min(round2even(maxval(diams) / smpd + 2. * COSMSKHALFWIDTH), ldim(1))
        moldiam         = smpd * box_for_pick
        mskdiam         = moldiam * MSK_EXP_FAC
        mskrad_in_pix   = (mskdiam / smpd) /2.
        maxdiam         = moldiam + moldiam * BOX_EXP_FAC
        box_for_extract = find_larger_magic_box(round2even(maxdiam / smpd))
        if( DEBUG )then
            print *, 'box_for_pick:    ', box_for_pick
            print *, 'mskdiam (in A):  ', mskdiam
            print *, 'mskrad_in_pix:   ', mskrad_in_pix
            print *, 'box_for_extract: ', box_for_extract
        endif
        call stkio_w%open(string(STREAM_SELECTED_REFS)//STK_EXT, smpd, 'write', box=box_for_extract, bufsz=nsel)
        do icls=1, nsel
            call stkio_r%get_image(selection(icls), cavgs(icls))
            if( ldim(1) > box_for_extract ) then
                call cavgs(icls)%clip_inplace([box_for_extract,box_for_extract,1])
            else
                call cavgs(icls)%pad_inplace([box_for_extract,box_for_extract,1])
            endif
            call cavgs(icls)%mask(mskrad_in_pix, 'softavg')
            call stkio_w%write(icls, cavgs(icls))
            call cavgs(icls)%kill
        end do
        deallocate(cavgs)
        call stkio_w%close
        call stkio_r%close
        ! write jpeg
        call mrc2jpeg_tiled(string(STREAM_SELECTED_REFS)//STK_EXT, string(STREAM_SELECTED_REFS)//JPG_EXT, n_xtiles=nxtiles, n_ytiles=nytiles)
        ! put back pointer to params_glob
        params_glob => params_ptr
        nullify(params_ptr)
    end subroutine process_selected_refs
    
    function get_latest_optics_map_id(optics_dir) result (lastmap)
        class(string), optional, intent(in)   :: optics_dir
        type(string), allocatable :: map_list(:)
        type(string) :: map_str, map_i_str
        integer      :: imap, prefix_len, testmap, lastmap
        lastmap = 0
        if(optics_dir .ne. "") then
            if(dir_exists(optics_dir)) call simple_list_files(optics_dir%to_char()//'/'// OPTICS_MAP_PREFIX //'*'//TXT_EXT, map_list)
        endif
        if(allocated(map_list)) then
            prefix_len = len(optics_dir%to_char() // '/' // OPTICS_MAP_PREFIX) + 1
            do imap=1, size(map_list)
                map_str   = map_list(imap)%to_char([prefix_len,map_list(imap)%strlen_trim()])
                map_i_str = swap_suffix(map_str, "", TXT_EXT)
                testmap   = map_i_str%to_int()
                if(testmap > lastmap) lastmap = testmap
                call map_str%kill
                call map_i_str%kill
            enddo
            deallocate(map_list)
        endif
    end function get_latest_optics_map_id

end module simple_stream_utils

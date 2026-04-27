!@descr: various stream utilities
module simple_stream_utils
use simple_core_module_api
use json_kinds
use json_module
use simple_image,               only: image
use simple_cmdline,             only: cmdline
use simple_qsys_env,            only: qsys_env
use simple_rec_list,            only: rec_list, project_rec
use simple_gui_utils,           only: mrc2jpeg_tiled
use simple_sp_project,          only: sp_project
use simple_default_clines,      only: set_automask2D_defaults
use simple_parameters,          only: parameters
use simple_stack_io,            only: stack_io
use simple_stream_communicator, only: stream_http_communicator 
implicit none
#include "simple_local_flags.inc"

contains

    subroutine terminate_stream( params, msg )
        class(parameters), intent(in) :: params
        character(len=*), intent(in) :: msg
        if(trim(params%async).eq.'yes') then
            if( file_exists(TERM_STREAM) ) call simple_end('**** '//trim(msg)//' ****', print_simple=.false.)
        endif
    end subroutine terminate_stream

    subroutine create_stream_project( spproj, cline, projname )
        type(sp_project),  intent(inout) :: spproj
        type(cmdline),     intent(inout) :: cline
        type(string),         intent(in) :: projname
        if( .not. cline%defined('projfile') ) then
            call cline%set('projname', projname)
            call cline%set('projfile', string(projname%to_char() // METADATA_EXT))
        endif
        call spproj%update_projinfo(cline)
        call spproj%update_compenv(cline)
        call spproj%write()
    end subroutine create_stream_project

    subroutine init_stream_qenv( params, qenv, envvar )
        type(parameters), intent(inout) :: params
        type(qsys_env),   intent(inout) :: qenv
        type(string),        intent(in) :: envvar
        character(len=STDLEN)           :: chunk_part_env
        integer                         :: envlen
        call get_environment_variable(envvar%to_char(), chunk_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(params, 1, exec_bin=string('simple_exec'), qsys_partition=string(trim(chunk_part_env)))
        else
            call qenv%new(params, 1, exec_bin=string('simple_exec'))
        end if
    end subroutine init_stream_qenv 

    subroutine import_new_projects( project_list, projects, n_mics_imported, n_ptcls_imported )
        type(rec_list),              intent(inout) :: project_list
        type(string),   allocatable, intent(inout) :: projects(:)
        integer,                     intent(inout) :: n_mics_imported, n_ptcls_imported
        type(project_rec)                          :: prec
        type(sp_project)                           :: spproj
        type(string)                               :: projabspath
        integer :: iproj, imic
        if( .not.allocated(projects) ) return
        if( size(projects) == 0      ) return
        do iproj=1, size(projects)
            call spproj%read(projects(iproj))
            ! because pick_extract purges state=0 and nptcls=0 mics,
            ! all mics can be assumed associated with particles
            if( spproj%os_mic%get_noris() == 0) then
                write(logfhandle, *) "ERROR: mic noris 0", projects(iproj)%to_char()
                call spproj%kill()
                cycle
            end if
            if( spproj%os_stk%get_noris() == 0) then
                write(logfhandle, *) "ERROR: stk noris 0", projects(iproj)%to_char()
                call spproj%kill()
                cycle
            end if
            projabspath = simple_abspath(projects(iproj))
            do imic = 1, spproj%os_mic%get_noris()
                prec%id          = project_list%size() + 1
                prec%projname    = projabspath
                prec%micind      = imic
                prec%nptcls      = spproj%os_mic%get_int(imic,'nptcls')
                prec%nptcls_sel  = prec%nptcls
                prec%included    = .false.
                n_mics_imported  = n_mics_imported + 1
                n_ptcls_imported = n_ptcls_imported + prec%nptcls
                call project_list%push_back(prec)
            enddo
            ! cleanup
            call spproj%kill()
        end do
    end subroutine import_new_projects

    !> To deal with dynamic user input diring streaming
    subroutine update_user_params( params, cline_here, httpcom )
        class(parameters),                   intent(inout) :: params
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
                if( abs(tilt_thres-params%tilt_thres) > 0.001) then
                     if(tilt_thres < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO LOW: ',tilt_thres
                     else if(tilt_thres > 1) then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO HIGH: ',tilt_thres
                     else
                        params%tilt_thres = tilt_thres
                        params%updated    = 'yes'
                        call cline_here%set('tilt_thres', params%tilt_thres)
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES UPDATED TO: ',tilt_thres
                     endif
                endif
            endif
            if( os%isthere(1,'beamtilt') ) then
                beamtilt = os%get(1,'beamtilt')
                if( beamtilt .eq. 1.0 ) then
                    params%beamtilt = 'yes'
                    params%updated  = 'yes'
                    call cline_here%set('beamtilt', params%beamtilt)
                    write(logfhandle,'(A)')'>>> OPTICS ASSIGNMENT UDPATED TO USE BEAMTILT'
                else if( beamtilt .eq. 0.0 ) then
                    params%beamtilt = 'no'
                    params%updated  = 'yes'
                    call cline_here%set('beamtilt', params%beamtilt)
                    write(logfhandle,'(A)')'>>> OPTICS ASSIGNMENT UDPATED TO IGNORE BEAMTILT'
                else
                    write(logfhandle,'(A,F8.2)')'>>> OPTICS UPDATE INVALID BEAMTILT VALUE: ',beamtilt
                endif
            endif
            if( os%isthere(1,'astigthreshold') ) then
                astigthreshold = os%get(1,'astigthreshold')
                if( abs(astigthreshold-params%astigthreshold) > 0.001) then
                     if(astigthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD TOO LOW: ',astigthreshold
                     else if(astigthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM TOO HIGH: ',astigthreshold
                     else
                        params%astigthreshold = astigthreshold
                        params%updated    = 'yes'
                        call cline_here%set('astigthreshold', params%astigthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD UPDATED TO: ',astigthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'ctfresthreshold') ) then
                ctfresthreshold = os%get(1,'ctfresthreshold')
                if( abs(ctfresthreshold-params%ctfresthreshold) > 0.001) then
                     if(ctfresthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD TOO LOW: ',ctfresthreshold
                     else if(ctfresthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION TOO HIGH: ',ctfresthreshold
                     else
                        params%ctfresthreshold = ctfresthreshold
                        params%updated    = 'yes'
                        call cline_here%set('ctfresthreshold', params%ctfresthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD UPDATED TO: ',ctfresthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'icefracthreshold') ) then
                icefracthreshold = os%get(1,'icefracthreshold')
                if( abs(icefracthreshold-params%icefracthreshold) > 0.001) then
                     if(icefracthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD TOO LOW: ',icefracthreshold
                     else if(icefracthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION TOO HIGH: ',icefracthreshold
                     else
                        params%icefracthreshold = icefracthreshold
                        params%updated    = 'yes'
                        call cline_here%set('icefracthreshold', params%icefracthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD UPDATED TO: ',icefracthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'moldiam_refine') ) then
                moldiam_refine = os%get(1,'moldiam_refine')
                if( abs(moldiam_refine-params%moldiam_refine) > 0.001) then
                     if(moldiam_refine < 10)then
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO LOW: ' , moldiam_refine
                     else if(moldiam_refine > 1000) then
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO HIGH: ', moldiam_refine
                     else
                        params%moldiam_refine = moldiam_refine
                        params%updated        = 'yes'
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
                    if( abs(real(moldiam_refine_int)-params%moldiam_refine) > 0.001) then
                        if(moldiam_refine_int < 10 .and. moldiam_refine_int > 0)then
                            write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO LOW: ' , real(moldiam_refine_int)
                        else if(moldiam_refine_int > 1000) then
                            write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO HIGH: ', real(moldiam_refine_int)
                        else
                            params%moldiam_refine = real(moldiam_refine_int)
                            params%updated        = 'yes'
                            write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER UPDATED TO: ', real(moldiam_refine_int)
                        end if
                    end if
                end if
                ! moldiam
                call httpcom%get_json_arg('moldiam', moldiam_int, found)
                if(found) then
                    if( abs(real(moldiam_int)-params%moldiam) > 0.001) then
                        if(moldiam_int < 20 .and. moldiam_int > 0)then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR DIAMETER TOO LOW: ' , real(moldiam_int)
                        else if(moldiam_int > 1000) then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR DIAMETER TOO HIGH: ', real(moldiam_int)
                        else
                            params%moldiam = real(moldiam_int)
                            params%updated        = 'yes'
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR DIAMETER UPDATED TO: ', real(moldiam_int)
                        end if
                    end if
                end if
                ! moldiam_ring
                call httpcom%get_json_arg('moldiam_ring', moldiam_ring_int, found)
                if(found) then
                    if( abs(real(moldiam_ring_int)-params%moldiam_ring) > 0.001) then
                        if(moldiam_ring_int < 20 .and. moldiam_ring_int > 0)then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR RING DIAMETER TOO LOW: ' , real(moldiam_ring_int)
                        else if(moldiam_ring_int > 1000) then
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR RING DIAMETER TOO HIGH: ', real(moldiam_ring_int)
                        else
                            params%moldiam_ring = real(moldiam_ring_int)
                            params%updated      = 'yes'
                            write(logfhandle,'(A,F8.2)')'>>> MOLECULAR RING DIAMETER UPDATED TO: ', real(moldiam_ring_int)
                        end if
                    end if
                end if
                ! interactive
                call httpcom%get_json_arg('interactive', interactive, found)
                if(found) then
                    params%interactive = interactive
                    params%updated        = 'yes'
                    write(logfhandle,'(A,A)')'>>> INTERACTIVE UPDATED TO: ', interactive
                end if
                ! ring
                call httpcom%get_json_arg('ring', ring, found)
                if(found) then
                    params%ring = ring
                    params%updated        = 'yes'
                    write(logfhandle,'(A,A)')'>>> RING UPDATED TO: ', ring
                end if
                ! astigthreshold
                call httpcom%get_json_arg('astigthreshold', astigthreshold_int, found)
                if(found) then
                    if( abs(real(astigthreshold_int)-params%astigthreshold) > 0.001) then
                        params%astigthreshold = real(astigthreshold_int)
                        params%updated        = 'yes'
                        write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD UPDATED TO: ', real(astigthreshold_int)
                    end if
                end if
                ! ctfresthreshold
                call httpcom%get_json_arg('ctfresthreshold', ctfresthreshold_int, found)
                if(found) then
                    if( abs(real(ctfresthreshold_int)-params%ctfresthreshold) > 0.001) then
                        params%ctfresthreshold = real(ctfresthreshold_int)
                        params%updated        = 'yes'
                        write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD UPDATED TO: ', real(ctfresthreshold_int)
                    end if
                end if
                ! icefracthreshold
                call httpcom%get_json_arg('icefracthreshold', icefracthreshold_dp, found)
                if(found) then
                    if( abs(icefracthreshold_dp-params%icefracthreshold) > 0.001) then
                        params%icefracthreshold = real(icefracthreshold_dp)
                        params%updated        = 'yes'
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

    subroutine wait_for_folder2( folder )
        class(string),                   intent(in)    :: folder
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
            end do
        endif
    end subroutine wait_for_folder2

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

    subroutine process_selected_refs( params, imgfile, smpd, selection, mskdiam, box_for_pick, box_for_extract, nxtiles, nytiles )
        use simple_image_msk, only: automask2d
        class(parameters), intent(inout) :: params
        class(string),   intent(in)    :: imgfile
        real,            intent(in)    :: smpd
        integer,         intent(in)    :: selection(:)
        real,            intent(out)   :: mskdiam
        integer,         intent(out)   :: box_for_pick, box_for_extract
        integer,         intent(inout) :: nxtiles, nytiles
        type(parameters)               :: local_params
        type(cmdline)                  :: cline
        type(stack_io)                 :: stkio_r, stkio_w
        type(image),       allocatable :: cavgs(:)
        real,              allocatable :: diams(:), shifts(:,:)
        logical,           parameter   :: DEBUG = .false.
        real    :: maxdiam, mskrad_in_pix, moldiam
        integer :: ldim(3), icls, ncls, nsel
        nsel = size(selection)
        if( nsel == 0 ) return
        write(logfhandle,'(A,I6,A)')'>>> USER SELECTED FROM POOL: ', nsel,' clusters'
        write(logfhandle,'(A,A)')'>>> WRITING SELECTED CLUSTERS TO: ', STREAM_SELECTED_REFS // STK_EXT
        ! set defaults
        call set_automask2D_defaults(cline)
        ! parse parameters
        call local_params%new(cline)
        call find_ldim_nptcls(imgfile, ldim, ncls)
        ldim(3) = 1
        local_params%msk  = real(ldim(1)/2) - COSMSKHALFWIDTH ! for automasking
        local_params%smpd = smpd                              ! for automasking
        if( DEBUG )then
            print *, 'imgfile:              ', imgfile%to_char()
            print *, 'file_exists(imgfile): ', file_exists(imgfile)
            print *, 'ldim:                 ', ldim(1), ldim(2), ldim(3)
            print *, 'ncls:                 ', ncls
            print *, 'params%msk:           ', local_params%msk
            print *, 'params%smpd:          ', local_params%smpd
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
        call automask2D(local_params, cavgs, local_params%ngrow, nint(local_params%winsz), local_params%edge, diams, shifts)
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
        ! mask memoization
        call cavgs(1)%memoize_mask_coords(box=box_for_extract)
        do icls=1, nsel
            call stkio_r%get_image(selection(icls), cavgs(icls))
            if( ldim(1) > box_for_extract ) then
                call cavgs(icls)%clip_inplace([box_for_extract,box_for_extract,1])
            else
                call cavgs(icls)%pad_inplace([box_for_extract,box_for_extract,1])
            endif
            call cavgs(icls)%mask2D_softavg(mskrad_in_pix)
            call stkio_w%write(icls, cavgs(icls))
            call cavgs(icls)%kill
        end do
        deallocate(cavgs)
        call stkio_w%close
        call stkio_r%close
        ! write jpeg
        call mrc2jpeg_tiled(string(STREAM_SELECTED_REFS)//STK_EXT, string(STREAM_SELECTED_REFS)//JPG_EXT, n_xtiles=nxtiles, n_ytiles=nytiles)
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

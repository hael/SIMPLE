module simple_nice
use, intrinsic :: iso_c_binding 
include 'simple_lib.f08'
use simple_sp_project, only: sp_project
use simple_socket_comm
use simple_histogram
use json_kinds
use json_module
use unix, only : c_pthread_t, c_pthread_mutex_t 
use unix, only : c_pthread_create, c_pthread_join
use unix, only : c_pthread_mutex_init, c_pthread_mutex_destroy
use unix, only : c_pthread_mutex_lock, c_pthread_mutex_unlock
implicit none

real, parameter, dimension(21)  :: astig_hist_bins = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
real, parameter, dimension(19)  :: ctf_res_bins    = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
real, parameter, dimension(21)  :: ice_score_bins  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
real, parameter, dimension(15)  :: res_2d_bins     = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0]

type, private :: nice_thread_comm_message
    type(c_pthread_mutex_t)               :: lock
    logical                               :: terminate
    integer                               :: msg_len, procid
    integer(kind=c_int16_t)               :: port
    character(len=39)                     :: ip
    character(kind=CK,len=:), allocatable :: msg_str
    type(json_value),         pointer     :: answer_json
end type nice_thread_comm_message

type(nice_thread_comm_message), target, private, save :: nice_thread_comm

type, private :: nice_stat_root
    character(len=:), allocatable :: status
    character(len=:), allocatable :: stage
    integer                       :: procid
    logical                       :: user_input = .false.
end type nice_stat_root

type, private :: nice_stat_thumb_image
    character(len=:),                 allocatable :: path
    character(len=:),                 allocatable :: boxfile
    character(len=:),                 allocatable :: uid
    logical                                       :: grid      = .false.
    integer                                       :: spritedim = 0
    integer                                       :: spriten   = 0
    integer                                       :: id        = 0
    integer                                       :: static_id = 0
    real                                          :: scale     = 0.0
    type(nice_stat_thumb_image_meta), allocatable :: metadata(:)
end type nice_stat_thumb_image

type, private :: nice_stat_thumb_image_meta
    character(len=64)             :: label
    integer                       :: type = 0 ! 1:integer, 2:real
    integer                       :: dp   = 1 
    integer,          allocatable :: data_int(:)
    real,             allocatable :: data_real(:)
end type nice_stat_thumb_image_meta

type, private :: nice_plot_doughnut
    character(len=64), allocatable :: labels(:)
    character(len=64), allocatable :: colours(:)
    integer,           allocatable :: data(:)
end type nice_plot_doughnut

type, private :: nice_plot_bar
    character(len=64), allocatable :: labels(:)
    character(len=64), allocatable :: colours(:)
    integer,           allocatable :: data(:)
end type nice_plot_bar

type, private :: nice_view_micrographs
    logical :: active               = .false.
    integer :: movies_imported      = -1
    integer :: movies_processed     = -1
    integer :: micrographs          = -1
    integer :: micrographs_rejected = -1
    integer :: compute_in_use       = -1
    real    :: avg_ctf_resolution   = 0.0
    real    :: avg_ice_score        = 0.0
    real    :: avg_astigmatism      = 0.0
    real    :: cutoff_ctf_res       = 0.0
    real    :: cutoff_ice_score     = 0.0
    real    :: cutoff_astigmatism   = 0.0
    character(16) :: last_movie_imported = ""
    type(nice_stat_thumb_image) :: thumbnail
    type(nice_stat_thumb_image), allocatable :: thumbnail_carousel(:)
    logical,                     allocatable :: thumbnail_carousel_mask(:)
    type(histogram)             :: ctf_res_histogram
    type(histogram)             :: astig_histogram
    type(histogram)             :: ice_score_histogram
    type(ori)                   :: plot_ori
end type nice_view_micrographs

type, private :: nice_view_optics
    logical :: active               = .false.
    integer :: micrographs          = 0
    integer :: micrographs_rejected = 0
    integer :: opticsgroups         = 0
    character(16) :: last_micrograph_imported = ""
    real,  allocatable :: opc(:,:)
end type nice_view_optics

type, private :: nice_view_pick
    logical :: active               = .false.
    integer :: micrographs_imported = 0
    integer :: micrographs_rejected = 0
    integer :: micrographs_picked   = 0
    integer :: gaussian_diameter    = 0
    integer :: suggested_diameter   = 0
    character(16) :: last_micrograph_imported = ""
    integer, allocatable        :: active_search_diameters(:)
    integer, allocatable        :: refined_search_diameters(:)
    integer, allocatable        :: complete_search_diameters(:)
    type(nice_stat_thumb_image) :: thumbnail
    type(nice_stat_thumb_image) :: pickrefs_thumbnail
    type(nice_stat_thumb_image), allocatable :: thumbnail_carousel(:)
    logical,                     allocatable :: thumbnail_carousel_mask(:)
    real,  allocatable :: coords(:,:)
end type nice_view_pick

type, private :: nice_view_cls2D
    logical :: active                    = .false.
    integer :: particles_extracted       = -1
    integer :: particles_imported        = -1
    integer :: iteration                 = -1
    integer :: number_classes            = -1
    integer :: number_classes_rejected   = -1
    integer :: number_particles_assigned = -1
    integer :: number_particles_rejected = -1
    integer :: snapshot_id               = -1
    integer :: boxsizea                  = -1
    real    :: maximum_resolution        = 0.0
    real    :: cutoff_res                = 0.0
    real    :: mskdiam                   = 0.0
    real    :: rejection_params(3)       = 0.0
    character(16) :: last_particles_imported = ""
    character(16) :: last_iteration          = ""
    character(16) :: snapshot_time           = ""
    character(6)  :: cutoff_type             = ""
    type(nice_stat_thumb_image) :: thumbnail
    type(nice_stat_thumb_image) :: pool_rejected_thumbnail
    type(nice_stat_thumb_image) :: chunk_rejected_thumbnail
    type(histogram)             :: res_histogram
    type(histogram)             :: ndev_histogram
    type(ori)                   :: plot_ori
    type(json_value),  pointer  :: snapshot_json
end type nice_view_cls2D

type, private :: nice_view_ini3D
    logical       :: active        = .false.
    integer       :: stage         = -1
    integer       :: number_states = -1
    real          :: lp            = 0.0
    character(16) :: last_stage_completed = ""
    type(oris)    :: vol_oris
    type(oris)    :: fsc_oris
end type nice_view_ini3D

type, private :: nice_view_vols
    logical   :: active      = .false.
    integer   :: number_vols = -1
    type(ori) :: plot_ori
end type nice_view_vols

type, public :: simple_nice_communicator
    integer,                        public  :: pid, current_checksum
    logical,                        private :: remote_active
    integer,                        private :: port
    character(len=39),              private :: ip
    type(simple_socket),            private :: socket
    type(json_core),                private :: stat_json
    type(json_value),   pointer,    private :: stat_json_root
    type(c_pthread_t),              private :: comm_thread
    logical,                        public  :: exit, stop
    type(nice_stat_root),           public  :: stat_root
    type(nice_view_micrographs),    public  :: view_micrographs
    type(nice_view_optics),         public  :: view_optics
    type(nice_view_pick),           public  :: view_pick
    type(nice_view_cls2D),          public  :: view_cls2D
    type(nice_view_ini3D),          public  :: view_ini3D
    type(nice_view_vols),           public  :: view_vols
    type(json_value),   pointer,    public  :: update_arguments

    contains
        procedure, private :: init_1, init_2
        generic            :: init => init_1, init_2
        procedure, public  :: terminate
        procedure, public  :: cycle
        procedure, private :: start_comm_thread
        procedure, private :: terminate_comm_thread
        procedure, private :: generate_stat_json
        procedure, private :: text_data_object_1 
        procedure, private :: text_data_object_2
        procedure, private :: text_data_object_3
        generic            :: text_data_object => text_data_object_1, text_data_object_2, text_data_object_3
        procedure, private :: doughnut_plot_object
        procedure, private :: bar_plot_object
        generic            :: plot_object => doughnut_plot_object, bar_plot_object
        procedure, private :: fsc_object
        procedure, private :: image_thumbnail_object
        procedure, private :: calculate_checksum
        procedure, public  :: update_micrographs
        procedure, public  :: update_optics
        procedure, public  :: update_pick
        procedure, public  :: update_cls2D
        procedure, public  :: update_ini3D
        procedure, public  :: update_vols
        procedure, public  :: update_from_project
        procedure, public  :: get_real_keys_json
        procedure, public  :: get_stat_json

end type simple_nice_communicator

    contains

    subroutine init_1(this, procid, serveraddr)
        class(simple_nice_communicator), intent(inout) :: this
        character(*),                    intent(in)    :: serveraddr
        integer,                         intent(in)    :: procid
        real,                            allocatable   :: histogram_rvec(:)
        character(len=STDLEN)                          :: ip, port_str
        integer                                        :: io_stat
        if(procid > 0 .and. serveraddr .ne. "") then
            this%remote_active = .true.
            port_str(:) = serveraddr(:)
            call split_str(port_str, ':', ip)
            this%ip   = trim(ip)
            this%port = str2int(trim(port_str), io_stat)
            if(io_stat .gt. 0) then
                this%remote_active = .false.
                write(logfhandle, *) ">>> REMOTE COMMUNICATION TO NICE DISABLED DUE TO MALFORMED ADDRESS STRING ", trim(serveraddr)
            end if
        else
            this%remote_active = .false.
        endif
        this%exit = .false.
        this%stop = .false.
        this%current_checksum = 0
        this%stat_root%status = "running"
        this%stat_root%stage  = "initialising"
        this%stat_root%procid = procid
        !init histograms
        allocate(histogram_rvec(size(ctf_res_bins)))
        histogram_rvec = ctf_res_bins
        call this%view_micrographs%ctf_res_histogram%new(histogram_rvec)
        deallocate(histogram_rvec)
        allocate(histogram_rvec(size(astig_hist_bins)))
        histogram_rvec = astig_hist_bins
        call this%view_micrographs%astig_histogram%new(histogram_rvec)
        deallocate(histogram_rvec)
        allocate(histogram_rvec(size(ice_score_bins)))
        histogram_rvec = ice_score_bins
        call this%view_micrographs%ice_score_histogram%new(histogram_rvec)
        deallocate(histogram_rvec)
        allocate(histogram_rvec(size(res_2d_bins)))
        histogram_rvec = res_2d_bins
        call this%view_cls2D%res_histogram%new(histogram_rvec)
        deallocate(histogram_rvec)
        nullify(this%view_cls2D%snapshot_json)
        ! init update_arguments
        call this%stat_json%create_null(this%update_arguments,'null')
        if(this%remote_active) then
            call this%start_comm_thread()
        end if
    end subroutine init_1

    subroutine init_2(this, procid, serveraddr)
        class(simple_nice_communicator), intent(inout) :: this
        class(string),                   intent(in)    :: serveraddr
        integer,                         intent(in)    :: procid
        call this%init_1(procid, serveraddr%to_char())
    end subroutine init_2

    subroutine terminate(this, stop, failed, export_project)
        class(simple_nice_communicator), intent(inout) :: this
        logical,               optional, intent(in)    :: stop, failed
        type(sp_project),      optional, intent(inout) :: export_project
        if(present(stop)) then
            if(stop) then
                this%stat_root%status = "stopped"
                this%stat_root%stage  = "user commanded stop"
            end if
        else if(present(failed)) then
            if(failed) then
                this%stat_root%status = "failed"
                this%stat_root%stage  = "failed"
            end if   
        else if(present(export_project)) then
            call this%update_from_project(export_project)
            this%stat_root%status = "finished"
            this%stat_root%stage  = "complete" 
        else
            this%stat_root%status = "finished"
            this%stat_root%stage  = "complete"
        end if
        call this%cycle()
        if(this%remote_active) then
            call this%terminate_comm_thread()
        end if
    end subroutine terminate

    subroutine cycle(this)
        class(simple_nice_communicator), intent(inout) :: this
        character(kind=CK,len=:), allocatable          :: json_str
        integer                                        :: rc, checksum, trial_count
        type(json_value),                pointer       :: update_arguments
        logical :: found, exit_argument, trial_exit, found_arguments_1, found_arguments_2, found_arguments_3
        call this%generate_stat_json()
        call this%stat_json%print_to_string(this%stat_json_root, json_str)
        checksum = this%calculate_checksum(json_str)
        if(checksum .ne. this%current_checksum) then
            this%current_checksum = checksum
            if(this%remote_active) then
                trial_exit = .false.
                do trial_count = 1, 10
                    rc = c_pthread_mutex_lock(nice_thread_comm%lock)
                    ! only queue message if msg_str is unallocated
                    if(.not. allocated(nice_thread_comm%msg_str)) then
                        nice_thread_comm%msg_len = len(json_str)
                        nice_thread_comm%msg_str = json_str
                        trial_exit = .true.
                    end if
                    rc = c_pthread_mutex_unlock(nice_thread_comm%lock)
                    if(trial_exit) exit
                    call sleep(1)
                end do
                if(trial_count .eq. 10) then 
                    write(logfhandle, *) "Max nice_thread_comm trials reached. Disabling remote nice communication"
                    this%remote_active = .false.
                    ! write json to file
                end if
            else
                ! write json to file
            endif
        endif
        rc = c_pthread_mutex_lock(nice_thread_comm%lock)
        if(associated(nice_thread_comm%answer_json)) then
            if(this%stat_json%count(nice_thread_comm%answer_json) > 0) then
                call this%stat_json%get(nice_thread_comm%answer_json, 'update', update_arguments, found_arguments_1)
                if(found_arguments_1) then
                    call this%stat_json%clone(update_arguments, this%update_arguments)
                    call this%stat_json%get(this%update_arguments, 'exit', exit_argument, found)
                    if(found) then
                        if(exit_argument) then
                            this%exit = .true.
                        end if
                    end if
                    call this%stat_json%get(this%update_arguments, 'stop', exit_argument, found)
                    if(found) then
                        if(exit_argument) then
                            this%stop = .true.
                        end if
                    end if
                    call this%stat_json%remove_if_present(nice_thread_comm%answer_json, "update")
                end if
                call this%stat_json%get(nice_thread_comm%answer_json, 'update2', update_arguments, found_arguments_2)
                if(found_arguments_2) call this%stat_json%rename(update_arguments, "update")
                call this%stat_json%get(nice_thread_comm%answer_json, 'update3', update_arguments, found_arguments_3)
                if(found_arguments_3) call this%stat_json%rename(update_arguments, "update")
            end if
           ! call this%stat_json%destroy(nice_thread_comm%answer_json)
        end if
        rc = c_pthread_mutex_unlock(nice_thread_comm%lock)
        if(allocated(json_str)) deallocate(json_str)
    end subroutine cycle

    subroutine start_comm_thread(this)
        class(simple_nice_communicator), intent(inout) :: this
        integer                                        :: rc
        nice_thread_comm%terminate = .false.
        nice_thread_comm%procid    = this%stat_root%procid
        nice_thread_comm%ip        = this%ip
        nice_thread_comm%port      = this%port
        rc = c_pthread_mutex_init(nice_thread_comm%lock, c_null_ptr)
        if(rc .ne. 0) then
            write(logfhandle, *) "simple_nice:start_comm_thread failed to init mutex. Remote sync deactivated"
            this%remote_active = .false.
            return
        endif
        rc = c_pthread_create(this%comm_thread, c_null_ptr, c_funloc(remote_heartbeat), c_null_ptr)
        if(rc .ne. 0) then
            write(logfhandle, *) "simple_nice:start_comm_thread failed to create pthread. Remote sync deactivated"
            this%remote_active = .false.
            return
        endif
    end subroutine start_comm_thread

    subroutine terminate_comm_thread(this)
        class(simple_nice_communicator), intent(inout) :: this
        type(c_ptr)                                    :: ptr
        integer                                        :: rc
        rc = c_pthread_mutex_lock(nice_thread_comm%lock)
        nice_thread_comm%terminate = .true.
        rc = c_pthread_mutex_unlock(nice_thread_comm%lock)
        rc = c_pthread_join(this%comm_thread, ptr)
        rc = c_pthread_mutex_destroy(nice_thread_comm%lock)
    end subroutine terminate_comm_thread

    subroutine remote_heartbeat() bind(c)
        integer                                :: rc, procid, heartbeat_count
        logical                                :: terminate, send
        type(simple_socket)                    :: socket
        type(json_core)                        :: json
        type(json_value), pointer              :: ans_json, update_json, update_arguments_1, update_arguments_2, update_arguments_3
        logical                                :: found_update, found_arguments_1, found_arguments_2, found_arguments_3
        character(kind=CK,len=:),  allocatable :: msg_str_loc
        character(kind=CK, len=:), allocatable :: ans_str_loc
        heartbeat_count = 0
        call json%initialize()
        do
            send = .false.
            rc = c_pthread_mutex_lock(nice_thread_comm%lock)
            terminate = nice_thread_comm%terminate
            procid    = nice_thread_comm%procid
            if(nice_thread_comm%msg_len > 0 .and. allocated(nice_thread_comm%msg_str)) then
                msg_str_loc = int2str_pad(procid, 10) // nice_thread_comm%msg_str // char(0)
                nice_thread_comm%msg_len = 0
                deallocate(nice_thread_comm%msg_str)
                send = .true.
                heartbeat_count = 0
            else
                if(heartbeat_count == 10) then
                    msg_str_loc = int2str_pad(procid, 10) // char(0)
                    send = .true.
                    heartbeat_count = 0
                else
                    heartbeat_count = heartbeat_count + 1
                end if
            end if
            rc = c_pthread_mutex_unlock(nice_thread_comm%lock)
            if(send) then
                call socket%open(port=nice_thread_comm%port, ip=trim(nice_thread_comm%ip))
                call socket%set_options
                call socket%send(msg_str_loc)
                call socket%read(ans_str_loc)
                call socket%close
                if(allocated(ans_str_loc)) then
                    rc = c_pthread_mutex_lock(nice_thread_comm%lock)
                   !  call json%destroy(nice_thread_comm%answer_json)
                   ! call json%parse(nice_thread_comm%answer_json, ans_str_loc)
                    call json%parse(ans_json, ans_str_loc)
                    call json%get(ans_json, 'update', update_json, found_update)
                    if(found_update) then
                        if(associated(nice_thread_comm%answer_json)) then
                            call json%get(nice_thread_comm%answer_json, 'update', update_arguments_1, found_arguments_1)
                            call json%get(nice_thread_comm%answer_json, 'update2', update_arguments_2, found_arguments_2)
                            call json%get(nice_thread_comm%answer_json, 'update3', update_arguments_3, found_arguments_3)
                            if(.not. found_arguments_1) then
                                call json%rename(update_json, 'update')
                                call json%add(nice_thread_comm%answer_json, update_json)
                            else if(.not. found_arguments_2) then
                                call json%rename(update_json, 'update2')
                                call json%add(nice_thread_comm%answer_json, update_json)
                            else if(.not. found_arguments_3) then
                                call json%rename(update_json, 'update3')
                                call json%add(nice_thread_comm%answer_json, update_json)
                            else
                                write(logfhandle, *) ">>> MAX NICE UPDATES REACHED"
                            end if
                        else
                            call json%destroy(nice_thread_comm%answer_json)
                            call json%create_object(nice_thread_comm%answer_json, '')
                            call json%add(nice_thread_comm%answer_json, update_json)
                        end if
                    end if
                    rc = c_pthread_mutex_unlock(nice_thread_comm%lock)
                end if
            end if
            if(allocated(msg_str_loc)) deallocate(msg_str_loc)
            if(allocated(ans_str_loc)) deallocate(ans_str_loc)
            if(terminate) exit
            if(.not. send) call sleep(1)
        end do
    end subroutine remote_heartbeat

    subroutine generate_stat_json(this)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value),                pointer       :: views, headlinestats, histograms, plots
        call this%stat_json%initialize(no_whitespace=.true.)
        call this%stat_json%create_object(this%stat_json_root,'')
        call this%stat_json%add(this%stat_json_root, "status",            this%stat_root%status)
        call this%stat_json%add(this%stat_json_root, "stage",             this%stat_root%stage)
        call this%stat_json%add(this%stat_json_root, "user_input",        this%stat_root%user_input)
        if(this%stat_root%user_input) this%stat_root%user_input = .false.
        call this%stat_json%create_object(headlinestats, "headlinestats")
        if(this%view_micrographs%active) call add_headline_micrographs()
        if(this%view_cls2D%active)       call add_headline_cls2D()
        call this%stat_json%create_object(views,         "views")
        if(this%view_micrographs%active) call add_view_micrographs()
        if(this%view_optics%active)      call add_view_optics()
        if(this%view_pick%active)        call add_view_pick()
        if(this%view_cls2D%active)       call add_view_cls2D()
        if(this%view_vols%active)        call add_view_vols()
        if(this%view_ini3D%active)       call add_view_ini3D()
        call this%stat_json%create_object(histograms,     "histograms")
        if(this%view_micrographs%active) call add_histograms_micrographs()
        if(this%view_cls2D%active)       call add_histograms_cls2D()
        call this%stat_json%create_object(plots,     "plots")
        if(this%view_optics%active) call add_plot_optics()
        call this%stat_json%add(this%stat_json_root, headlinestats)
        call this%stat_json%add(this%stat_json_root, views)
        call this%stat_json%add(this%stat_json_root, histograms)
        call this%stat_json%add(this%stat_json_root, plots)

        contains

            subroutine add_headline_micrographs()
                type(json_value), pointer :: movies_headline_section, micrographs_headline_section
                type(json_value), pointer :: movies_imported, last_movie_imported, n_micrographs
                type(json_value), pointer :: avg_ctf_resolution
                type(json_value), pointer :: avg_ice_score
                type(json_value), pointer :: latest_mic_img
                ! movies
                if(this%view_micrographs%movies_imported >= 0) then
                    call this%stat_json%create_object(movies_headline_section, 'movies')
                    call this%text_data_object(movies_imported,     "movies_imported",     this%view_micrographs%movies_imported)
                    call this%stat_json%add(movies_headline_section, movies_imported)
                    if(trim(this%view_micrographs%last_movie_imported) .ne. "") then
                        call this%text_data_object(last_movie_imported, "last_movie_imported", this%view_micrographs%last_movie_imported)
                        call this%stat_json%add(movies_headline_section, last_movie_imported)
                    end if
                    call this%stat_json%add(headlinestats, movies_headline_section)
                end if
                !! CTF ref
                ! micrographs
                if(this%view_micrographs%micrographs >= 0) then
                    call this%stat_json%create_object(micrographs_headline_section, 'micrographs')
                    call this%text_data_object(n_micrographs,        "micrographs",          this%view_micrographs%micrographs)
                    call this%stat_json%add(micrographs_headline_section,          n_micrographs)
                    if(this%view_micrographs%avg_ctf_resolution > 0.) then
                        call this%text_data_object(avg_ctf_resolution, "avg_ctf_resolution", this%view_micrographs%avg_ctf_resolution, dp=1)
                        call this%stat_json%add(micrographs_headline_section, avg_ctf_resolution)
                    end if
                    !! Ice score
                    if(this%view_micrographs%avg_ice_score > 0.) then
                        call this%text_data_object(avg_ice_score, "avg_ice_score", this%view_micrographs%avg_ice_score, dp=2)
                        call this%stat_json%add(micrographs_headline_section, avg_ice_score)
                    end if
                    !! thumbnail
                    if(this%view_micrographs%thumbnail%path .ne. "" .and. this%view_micrographs%thumbnail%id > 0) then
                        call this%image_thumbnail_object(latest_mic_img, this%view_micrographs%thumbnail, name="latest")
                        call this%stat_json%add(micrographs_headline_section, latest_mic_img)
                    end if
                    !! add to sections
                    call this%stat_json%add(headlinestats, micrographs_headline_section)
                end if
            end subroutine add_headline_micrographs

            subroutine add_view_micrographs()
                type(json_value), pointer :: micrographs, thumbnail_carousel
                type(json_value), pointer :: movies_section, movies_doughnut, micrographs_section
                type(json_value), pointer :: movies_imported, last_movie_imported, n_micrographs
                type(json_value), pointer :: micrographs_rejected, avg_ctf_resolution, interactive_plot
                type(json_value), pointer :: avg_ice_score, avg_astigmatism
                type(json_value), pointer :: carousel_img
                type(nice_plot_doughnut)  :: status_plot
                integer                   :: i
                call this%stat_json%create_object(micrographs, 'micrographs')
                ! movies
                if(this%view_micrographs%movies_imported >= 0) then
                    call this%stat_json%create_object(movies_section,          'movies')
                    call this%text_data_object(movies_imported,     "movies_imported",     this%view_micrographs%movies_imported)
                    call this%stat_json%add(movies_section,          movies_imported)
                    if(trim(this%view_micrographs%last_movie_imported) .ne. "") then
                        call this%text_data_object(last_movie_imported, "last_movie_imported", this%view_micrographs%last_movie_imported)
                        call this%stat_json%add(movies_section,          last_movie_imported)
                    end if
                    allocate(status_plot%labels(2))
                    allocate(status_plot%colours(2))
                    allocate(status_plot%data(2))
                    status_plot%labels(1) = "processed"
                    status_plot%labels(2) = "queued"
                    status_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                    status_plot%colours(2) = "rgba(211, 211, 211, 0.5)"
                    status_plot%data(1) = this%view_micrographs%movies_processed
                    status_plot%data(2) = this%view_micrographs%movies_imported - this%view_micrographs%movies_processed 
                    call this%plot_object(movies_doughnut, status_plot)
                    call this%stat_json%add(movies_section,          movies_doughnut)
                    call this%stat_json%add(micrographs,             movies_section)
                end if
                ! micrographs
                if(this%view_micrographs%micrographs >= 0) then
                    call this%stat_json%create_object(micrographs_section,          'micrographs')
                    !!n_micrographs
                    call this%text_data_object(n_micrographs,        "micrographs",          this%view_micrographs%micrographs)
                    call this%stat_json%add(micrographs_section,          n_micrographs)
                    !! rejected micrographs
                    if(this%view_micrographs%micrographs_rejected >= 0) then
                        call this%text_data_object(micrographs_rejected, "micrographs_rejected", this%view_micrographs%micrographs_rejected)
                        call this%stat_json%add(micrographs_section,          micrographs_rejected)
                    end if
                    !! CTF ref
                    if(this%view_micrographs%avg_ctf_resolution > 0.) then
                        call this%text_data_object(avg_ctf_resolution, "avg_ctf_resolution", this%view_micrographs%avg_ctf_resolution, dp=1)
                        call this%stat_json%add(micrographs_section,          avg_ctf_resolution)
                    end if
                    !! Ice score
                    if(this%view_micrographs%avg_ice_score > 0.) then
                        call this%text_data_object(avg_ice_score, "avg_ice_score", this%view_micrographs%avg_ice_score, dp=2)
                        call this%stat_json%add(micrographs_section,          avg_ice_score)
                    end if
                    !! Astigmatism
                    if(this%view_micrographs%avg_astigmatism > 0.) then
                        call this%text_data_object(avg_astigmatism, "avg_astigmatism", this%view_micrographs%avg_astigmatism, dp=2)
                        call this%stat_json%add(micrographs_section, avg_astigmatism)
                    end if
                    !! add to sections
                    call this%stat_json%add(micrographs, micrographs_section)
                end if

                ! thumbnail carousel
                if (allocated(this%view_micrographs%thumbnail_carousel)) then
                    if(size(this%view_micrographs%thumbnail_carousel) > 0) then
                        call this%stat_json%create_object(thumbnail_carousel, 'latest')
                        call this%stat_json%add(thumbnail_carousel, "type", "thumbnail_carousel")
                        this%view_micrographs%thumbnail_carousel = pack(this%view_micrographs%thumbnail_carousel, this%view_micrographs%thumbnail_carousel_mask)
                        if (allocated(this%view_micrographs%thumbnail_carousel_mask)) deallocate(this%view_micrographs%thumbnail_carousel_mask)
                        allocate(this%view_micrographs%thumbnail_carousel_mask(size(this%view_micrographs%thumbnail_carousel)))
                        this%view_micrographs%thumbnail_carousel_mask = .true.
                        do i=1, size(this%view_micrographs%thumbnail_carousel)
                            call this%image_thumbnail_object(carousel_img, this%view_micrographs%thumbnail_carousel(i), name="thumbnail_" // int2str(i))
                            call this%stat_json%add(thumbnail_carousel, carousel_img)
                        end do
                        call this%stat_json%add(micrographs, thumbnail_carousel)
                    end if
                end if
                ! interactive plot
                call create_interactive_plot(this%view_micrographs%plot_ori, interactive_plot)
                call this%stat_json%add(micrographs, interactive_plot)
                ! add to views
                call this%stat_json%add(views, micrographs)
                ! deallocate
                if(allocated(status_plot%labels))  deallocate(status_plot%labels)
                if(allocated(status_plot%colours)) deallocate(status_plot%colours)
                if(allocated(status_plot%data))    deallocate(status_plot%data)
            end subroutine add_view_micrographs

            subroutine add_histograms_micrographs()
                type(json_value), pointer :: ctf_res_hist, ice_score_hist, astig_hist
                type(nice_plot_bar)       :: ctf_res_plot, ice_score_plot, astig_plot
                character(len=:), allocatable              :: val_tmp
                integer                   :: i, n
                !! ctf res
                n = size(ctf_res_bins)
                allocate(ctf_res_plot%labels(n))
                allocate(ctf_res_plot%data(n))
                allocate(ctf_res_plot%colours(1))
                ctf_res_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                do i=1, n
                    ctf_res_plot%labels(i) = int2str(int(ctf_res_bins(i))) // 'Å'
                    ctf_res_plot%data(i)   = this%view_micrographs%ctf_res_histogram%get(i)
                end do
                if(this%view_micrographs%cutoff_ctf_res .gt. 0.0) then
                    call this%plot_object(ctf_res_hist, ctf_res_plot, name="ctf_resolution", value=trim(int2str(int(this%view_micrographs%cutoff_ctf_res))) // 'Å')
                else
                    call this%plot_object(ctf_res_hist, ctf_res_plot, name="ctf_resolution")
                end if
                call this%stat_json%add(histograms, ctf_res_hist)
                !! ice score
                n = size(ice_score_bins)
                allocate(ice_score_plot%labels(n))
                allocate(ice_score_plot%data(n))
                allocate(ice_score_plot%colours(1))
                ice_score_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                do i=1, n
                    write(ice_score_plot%labels(i), "(F4.2)") ice_score_bins(i)
                    ice_score_plot%data(i)   = this%view_micrographs%ice_score_histogram%get(i)
                end do
                if(this%view_micrographs%cutoff_ice_score .gt. 0.0) then
                    allocate(character(len=4) :: val_tmp)
                    write(val_tmp, "(F4.2)") this%view_micrographs%cutoff_ice_score
                    call this%plot_object(ice_score_hist, ice_score_plot, name="ice_score", value=val_tmp)
                    deallocate(val_tmp)
                else
                    call this%plot_object(ice_score_hist, ice_score_plot, name="ice_score")
                end if
                call this%stat_json%add(histograms, ice_score_hist)
                !! astig
                n = size(astig_hist_bins)
                allocate(astig_plot%labels(n))
                allocate(astig_plot%data(n))
                allocate(astig_plot%colours(1))
                astig_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                do i=1, n
                    astig_plot%labels(i) = int2str(int(astig_hist_bins(i))) // '%'
                    astig_plot%data(i)   = this%view_micrographs%astig_histogram%get(i)
                end do
                call this%plot_object(astig_hist, astig_plot, name="astigmatism")
                if(this%view_micrographs%cutoff_astigmatism .gt. 0.0) then
                    call this%plot_object(astig_hist, astig_plot, name="astigmatism", value=trim(int2str(int(this%view_micrographs%cutoff_astigmatism))) // '%')
                else
                    call this%plot_object(astig_hist, astig_plot, name="astigmatism")
                end if
                call this%stat_json%add(histograms, astig_hist)
                ! deallocate
                if(allocated(ctf_res_plot%labels))    deallocate(ctf_res_plot%labels)
                if(allocated(ctf_res_plot%data))      deallocate(ctf_res_plot%data)
                if(allocated(ctf_res_plot%colours))   deallocate(ctf_res_plot%colours)
                if(allocated(ice_score_plot%labels))  deallocate(ice_score_plot%labels)
                if(allocated(ice_score_plot%data))    deallocate(ice_score_plot%data)
                if(allocated(ice_score_plot%colours)) deallocate(ice_score_plot%colours)
                if(allocated(astig_plot%labels))      deallocate(astig_plot%labels)
                if(allocated(astig_plot%data))        deallocate(astig_plot%data)
                if(allocated(astig_plot%colours))     deallocate(astig_plot%colours)
            end subroutine add_histograms_micrographs

            subroutine add_histograms_cls2D()
                type(json_value),     pointer :: res_hist
                type(nice_plot_bar)           :: res_plot
                integer                       :: i, n
                !! res
                n = size(res_2d_bins)
                allocate(res_plot%labels(n))
                allocate(res_plot%data(n))
                allocate(res_plot%colours(1))
                res_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                do i=1, n
                    res_plot%labels(i) = int2str(int(res_2d_bins(i))) // 'Å'
                    res_plot%data(i)   = this%view_cls2D%res_histogram%get(i)
                end do
                if(this%view_cls2D%cutoff_res .gt. 0.0) then
                    call this%plot_object(res_hist, res_plot, name="class_resolution", value=trim(int2str(int(this%view_cls2D%cutoff_res))) // 'Å', extra_1=trim(this%view_cls2D%cutoff_type))
                else
                    call this%plot_object(res_hist, res_plot, name="class_resolution")
                end if
                call this%stat_json%add(histograms, res_hist)
                ! deallocate
                if(allocated(res_plot%labels))    deallocate(res_plot%labels)
                if(allocated(res_plot%data))      deallocate(res_plot%data)
                if(allocated(res_plot%colours))   deallocate(res_plot%colours)
            end subroutine add_histograms_cls2D

            subroutine add_view_optics()
                type(json_value), pointer :: optics
                type(json_value), pointer :: micrographs_section, optics_groups_section
                type(json_value), pointer :: micrographs_imported, last_micrograph_imported, optics_groups
                call this%stat_json%create_object(optics, 'optics')
                ! micrographs
                call this%stat_json%create_object(micrographs_section, 'micrographs')
                call this%text_data_object(micrographs_imported,     "micrographs_imported",     this%view_optics%micrographs)
                call this%stat_json%add(micrographs_section, micrographs_imported)
                if(trim(this%view_optics%last_micrograph_imported) .ne. "") then
                    call this%text_data_object(last_micrograph_imported, "last_micrograph_imported", this%view_optics%last_micrograph_imported)
                    call this%stat_json%add(micrographs_section, last_micrograph_imported)
                end if
                call this%stat_json%add(optics, micrographs_section)
                ! optics groups
                call this%stat_json%create_object(optics_groups_section, 'optics_groups')
                call this%text_data_object(optics_groups,     "optics_groups_assigned",     this%view_optics%opticsgroups)
                call this%stat_json%add(optics_groups_section, optics_groups)
                call this%stat_json%add(optics, optics_groups_section)
                ! add to views
                call this%stat_json%add(views, optics)
            end subroutine add_view_optics

            subroutine create_interactive_plot(seg_ori, interactive_plot)
                type(ori),                 intent(in)    :: seg_ori
                type(json_value), pointer, intent(inout) :: interactive_plot
                type(json_value), pointer :: keys
                call this%stat_json%create_object(interactive_plot, 'interactive_plot')
                call this%stat_json%add(interactive_plot, 'type', 'plot_interactive')
                call this%get_real_keys_json(keys, seg_ori)
                call this%stat_json%add(interactive_plot, keys)
            end subroutine create_interactive_plot

            subroutine add_plot_optics()
                type(json_value), pointer :: optics_plot, datasets, dataset, data, xy
                integer                   :: i, j
                call this%stat_json%create_object(optics_plot, 'assignments')
                call this%stat_json%add(optics_plot, 'type', 'plot_scatter')
                call this%stat_json%create_array(datasets,     'datasets')
                do i = 1, this%view_optics%opticsgroups
                    call this%stat_json%create_object(dataset, 'dataset')
                    call this%stat_json%create_array(data, 'data')
                    call this%stat_json%add(dataset, 'label', 'optics group ' // int2str(i))
                    do j = 1, size(this%view_optics%opc, 1)
                        if(this%view_optics%opc(j, 1) == i) then
                            call this%stat_json%create_object(xy, 'xy')
                            call this%stat_json%add(xy, 'x', dble(this%view_optics%opc(j, 2)))
                            call this%stat_json%add(xy, 'y', dble(this%view_optics%opc(j, 3)))
                            call this%stat_json%add(data, xy)
                        end if
                    end do
                    call this%stat_json%add(dataset,  data)
                    call this%stat_json%add(datasets, dataset)
                end do
                call this%stat_json%add(optics_plot, datasets)
                call this%stat_json%add(plots,       optics_plot)
            end subroutine add_plot_optics

            subroutine add_view_pick()
                type(json_value), pointer :: pick
                type(json_value), pointer :: micrographs_section, search_diameters_section, thumbnail_carousel, carousel_img
                type(json_value), pointer :: micrographs_doughnut, last_micrograph_imported, micrographs_imported
                type(json_value), pointer :: active_search, refined_search, complete_search, pick_section, gaussian_diameter
                type(json_value), pointer :: pickrefs_section, pickrefs_image
                type(nice_plot_doughnut)  :: status_plot
                integer                   :: i
                call this%stat_json%create_object(pick, 'pick')
                ! micrographs
                call this%stat_json%create_object(micrographs_section, 'micrographs')
                call this%text_data_object(micrographs_imported, "micrographs_imported",     this%view_pick%micrographs_imported)
                call this%stat_json%add(micrographs_section, micrographs_imported)
                if(trim(this%view_pick%last_micrograph_imported) .ne. "") then
                    call this%text_data_object(last_micrograph_imported, "last_micrograph_imported", this%view_pick%last_micrograph_imported)
                    call this%stat_json%add(micrographs_section, last_micrograph_imported)
                end if
                allocate(status_plot%labels(3))
                allocate(status_plot%colours(3))
                allocate(status_plot%data(3))
                status_plot%labels(1) = "picked"
                status_plot%labels(2) = "rejected"
                status_plot%labels(3) = "queued"
                status_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                status_plot%colours(2) = "rgba(255, 99, 71, 0.75)"
                status_plot%colours(3) = "rgba(211, 211, 211, 0.5)"
                status_plot%data(1) = this%view_pick%micrographs_picked
                status_plot%data(2) = this%view_pick%micrographs_rejected
                status_plot%data(3) = this%view_pick%micrographs_imported - this%view_pick%micrographs_picked - this%view_pick%micrographs_rejected
                call this%plot_object(micrographs_doughnut, status_plot)
                call this%stat_json%add(micrographs_section, micrographs_doughnut)
                call this%stat_json%add(pick, micrographs_section)
                ! pick
                call this%stat_json%create_object(pick_section, 'picking')
                if(this%view_pick%gaussian_diameter .gt. 0) then
                    call this%text_data_object(gaussian_diameter, "gaussian_diameter", this%view_pick%gaussian_diameter)
                    call this%stat_json%add(pick_section, gaussian_diameter)
                    call this%stat_json%add(pick, pick_section)
                end if
                ! pickrefs
                if(this%view_pick%pickrefs_thumbnail%path .ne. "" .and. this%view_pick%pickrefs_thumbnail%id > 0) then
                    call this%stat_json%create_object(pickrefs_section, 'picking templates')
                    call this%image_thumbnail_object(pickrefs_image, this%view_pick%pickrefs_thumbnail)
                    call this%stat_json%add(pickrefs_section, pickrefs_image)
                    call this%stat_json%add(pickrefs_section, "type", "thumbnail_grid")
                    call this%stat_json%add(pick, pickrefs_section)
                end if
                ! thumbnail carousel
                if (allocated(this%view_pick%thumbnail_carousel)) then
                    if(size(this%view_pick%thumbnail_carousel) > 0) then
                         call this%stat_json%create_object(thumbnail_carousel, 'latest')
                         call this%stat_json%add(thumbnail_carousel, "type", "thumbnail_carousel")
                         this%view_pick%thumbnail_carousel = pack(this%view_pick%thumbnail_carousel, this%view_pick%thumbnail_carousel_mask)
                         if (allocated(this%view_pick%thumbnail_carousel_mask)) deallocate(this%view_pick%thumbnail_carousel_mask)
                         allocate(this%view_pick%thumbnail_carousel_mask(size(this%view_pick%thumbnail_carousel)))
                         this%view_pick%thumbnail_carousel_mask = .true.
                         do i=1, size(this%view_pick%thumbnail_carousel)
                             call this%image_thumbnail_object(carousel_img, this%view_pick%thumbnail_carousel(i), name="thumbnail_" // int2str(i))
                             call this%stat_json%add(thumbnail_carousel, carousel_img)
                         end do
                         call this%stat_json%add(pick, thumbnail_carousel)
                    end if
                 end if
                ! pick diameters
                if(allocated(this%view_pick%active_search_diameters) .or. allocated(this%view_pick%refined_search_diameters) .or. allocated(this%view_pick%complete_search_diameters)) then
                    call this%stat_json%create_object(search_diameters_section, 'search_diameters')
                    call this%stat_json%add(search_diameters_section, "type", "search_slider")
                    if(allocated(this%view_pick%active_search_diameters)) then
                        call this%stat_json%create_array(active_search, 'active')
                        do i=1, size(this%view_pick%active_search_diameters)
                             call this%stat_json%add(active_search, '', this%view_pick%active_search_diameters(i))
                        end do
                        call this%stat_json%add(search_diameters_section, active_search)
                    end if 
                    if(allocated(this%view_pick%refined_search_diameters)) then
                        call this%stat_json%create_array(refined_search, 'refined')
                        do i=1, size(this%view_pick%refined_search_diameters)
                             call this%stat_json%add(refined_search, '', this%view_pick%refined_search_diameters(i))
                        end do
                        call this%stat_json%add(search_diameters_section, refined_search)
                    end if 
                    if(allocated(this%view_pick%complete_search_diameters)) then
                        call this%stat_json%create_array(complete_search, 'complete')
                        do i=1, size(this%view_pick%complete_search_diameters)
                             call this%stat_json%add(complete_search, '', this%view_pick%complete_search_diameters(i))
                        end do
                        call this%stat_json%add(search_diameters_section, complete_search)
                    end if
                    if(this%view_pick%suggested_diameter > 0) then
                        call this%stat_json%add(search_diameters_section, "suggested", this%view_pick%suggested_diameter)
                    end if
                    call this%stat_json%add(pick, search_diameters_section)
                end if
                ! add to sections
                call this%stat_json%add(views, pick)
                ! deallocate
                if(allocated(status_plot%labels))  deallocate(status_plot%labels)
                if(allocated(status_plot%colours)) deallocate(status_plot%colours)
                if(allocated(status_plot%data))    deallocate(status_plot%data)
            end subroutine add_view_pick

            subroutine add_headline_cls2D()
                type(json_value), pointer :: particles_section, cls2D_section
                type(json_value), pointer :: particles_imported, particles_assigned, last_iteration_time
                type(json_value), pointer :: latest_classes_image, number_classes
                ! particles
                if(this%view_cls2D%particles_imported > -1) then
                    call this%stat_json%create_object(particles_section, 'particles')
                    call this%text_data_object(particles_imported, "particles_imported", this%view_cls2D%particles_imported)
                    call this%stat_json%add(particles_section, particles_imported)
                    if(this%view_cls2D%number_particles_assigned > -1) then
                        call this%text_data_object(particles_assigned, "particles_assigned", this%view_cls2D%number_particles_assigned)
                        call this%stat_json%add(particles_section, particles_assigned)
                    end if
                    call this%stat_json%add(headlinestats, particles_section)
                end if
                ! 2Dcls2D_section
                 if(this%view_cls2D%number_classes > -1) then
                    call this%stat_json%create_object(cls2D_section, '2D')
                    if(this%view_cls2D%number_classes > -1) then
                        call this%text_data_object(number_classes, "number_classes", this%view_cls2D%number_classes)
                        call this%stat_json%add(cls2D_section, number_classes)
                    end if
                    if(trim(this%view_cls2D%last_iteration) .ne. "") then
                        call this%text_data_object(last_iteration_time, "last_iteration_time", this%view_cls2D%last_iteration)
                        call this%stat_json%add(cls2D_section, last_iteration_time)
                    end if
                    !! thumbnail
                    if(this%view_cls2D%thumbnail%path .ne. "" .and. this%view_cls2D%thumbnail%id > 0) then
                        call this%image_thumbnail_object(latest_classes_image, this%view_cls2D%thumbnail, name="latest")
                        call this%stat_json%add(cls2D_section, latest_classes_image)
                    end if
                    call this%stat_json%add(headlinestats, cls2D_section)
                end if
            end subroutine add_headline_cls2D

            subroutine add_view_cls2D()
                type(json_value), pointer :: cls2D
                type(json_value), pointer :: particles_section, cls2D_section, grid_section, snapshot_section
                type(json_value), pointer :: particles_extracted, particles_imported, last_particles_imported
                type(json_value), pointer :: iteration, number_classes, number_classes_rejected, assignment_doughnut, last_iteration_time, interactive_plot
                type(json_value), pointer :: latest_classes_image, maximum_resolution, chunk_rejected_classes_image, chunk_rejected_grid_section, snapshot_json
                type(json_value), pointer :: pool_rejected_grid_section, pool_rejected_classes_image, snapshot_last, snapshot_last_id, rejection, lpthres, ndev, mskdiam, boxsizea
                type(nice_plot_doughnut)  :: status_plot
                call this%stat_json%create_object(cls2D, 'cls2D')
                ! particles
                if(this%view_cls2D%particles_imported > -1 .or. this%view_cls2D%number_particles_assigned > -1) then
                    call this%stat_json%create_object(particles_section, 'particles')
                    if(this%view_cls2D%particles_extracted > -1) then
                        call this%text_data_object(particles_extracted, "particles_extracted", this%view_cls2D%particles_extracted)
                        call this%stat_json%add(particles_section, particles_extracted)
                    end if
                     if(this%view_cls2D%particles_imported > -1) then
                        call this%text_data_object(particles_imported, "particles_imported", this%view_cls2D%particles_imported)
                        call this%stat_json%add(particles_section, particles_imported)
                    end if
                    if(trim(this%view_cls2D%last_particles_imported) .ne. "") then
                        call this%text_data_object(last_particles_imported, "last_particles_imported", this%view_cls2D%last_particles_imported)
                        call this%stat_json%add(particles_section, last_particles_imported)
                    end if
                    ! assignment plot
                    allocate(status_plot%labels(3))
                    allocate(status_plot%colours(3))
                    allocate(status_plot%data(3))
                    status_plot%labels(1) = "assigned"
                    status_plot%labels(2) = "rejected"
                    status_plot%labels(3) = "queued"
                    status_plot%colours(1) = "rgba(30, 144, 255, 0.5)"
                    status_plot%colours(2) = "rgba(255, 99, 71, 0.75)"
                    status_plot%colours(3) = "rgba(211, 211, 211, 0.5)"
                    status_plot%data(1) = this%view_cls2D%number_particles_assigned
                    status_plot%data(2) = this%view_cls2D%number_particles_rejected
                    status_plot%data(3) = this%view_cls2D%particles_imported - (this%view_cls2D%number_particles_assigned + this%view_cls2D%number_particles_rejected)
                    call this%plot_object(assignment_doughnut, status_plot)
                    call this%stat_json%add(particles_section, assignment_doughnut)
                    call this%stat_json%add(cls2D, particles_section)
                end if
                ! 2D
                if(this%view_cls2D%number_classes > -1) then
                    call this%stat_json%create_object(cls2D_section, '2D')
                    if(this%view_cls2D%iteration > 0) then
                        call this%text_data_object(iteration, "iteration", this%view_cls2D%iteration)
                        call this%stat_json%add(cls2D_section, iteration)
                    end if
                    if(this%view_cls2D%number_classes > -1) then
                        call this%text_data_object(number_classes, "number_classes", this%view_cls2D%number_classes)
                        call this%stat_json%add(cls2D_section, number_classes)
                    end if
                    if(this%view_cls2D%number_classes_rejected > -1) then
                        call this%text_data_object(number_classes_rejected, "number_classes_rejected", this%view_cls2D%number_classes_rejected)
                        call this%stat_json%add(cls2D_section, number_classes_rejected)
                    end if
                    if(this%view_cls2D%maximum_resolution .gt. 0.0) then
                        call this%text_data_object(maximum_resolution, "maximum_resolution", this%view_cls2D%maximum_resolution, dp=1)
                        call this%stat_json%add(cls2D_section, maximum_resolution)
                    end if
                    if(trim(this%view_cls2D%last_iteration) .ne. "") then
                        call this%text_data_object(last_iteration_time, "last_iteration_time", this%view_cls2D%last_iteration)
                        call this%stat_json%add(cls2D_section, last_iteration_time)
                    end if
                    !! rejection parameters
                    if(maxval(this%view_cls2D%rejection_params) .gt. 0.0) then
                        call this%text_data_object(rejection, "rejection",  this%view_cls2D%rejection_params(1), dp=1)
                        call this%stat_json%add(cls2D_section, rejection)
                        call this%text_data_object(lpthres, "lp_threshold", this%view_cls2D%rejection_params(2), dp=1)
                        call this%stat_json%add(cls2D_section, lpthres)
                        call this%text_data_object(ndev, "ndev", this%view_cls2D%rejection_params(3), dp=1)
                        call this%stat_json%add(cls2D_section, ndev)
                    end if
                    if(this%view_cls2D%mskdiam .gt. 0.0) then
                        call this%text_data_object(mskdiam, "mskdiam", round2even(this%view_cls2D%mskdiam))
                        call this%stat_json%add(cls2D_section, mskdiam)
                    end if
                    if(this%view_cls2D%boxsizea .gt. 0) then
                        call this%text_data_object(boxsizea, "boxsizea", this%view_cls2D%boxsizea)
                        call this%stat_json%add(cls2D_section, boxsizea)
                    end if
                    call this%stat_json%add(cls2D, cls2D_section)
                end if
                ! snapshots
                if(this%view_cls2D%snapshot_id .gt. 0) then
                    call this%stat_json%create_object(snapshot_section, 'snapshot')
                    call this%text_data_object(snapshot_last_id, "snapshot_last_id",   this%view_cls2D%snapshot_id)
                    call this%stat_json%add(snapshot_section, snapshot_last_id)
                    call this%text_data_object(snapshot_last,    "snapshot_last_time", this%view_cls2D%snapshot_time)
                    call this%stat_json%add(snapshot_section, snapshot_last)
                    ! snapshot json
                    if(associated(this%view_cls2D%snapshot_json)) then
                        call this%stat_json%clone(this%view_cls2D%snapshot_json, snapshot_json)
                        call this%stat_json%rename(snapshot_json, "", "snapshot_json")
                        call this%stat_json%add(snapshot_section, snapshot_json)
                    end if
                    call this%stat_json%add(cls2D, snapshot_section)
                end if
                !! thumbnail grid
                if(this%view_cls2D%thumbnail%path .ne. "" .and. this%view_cls2D%thumbnail%id > 0) then
                    call this%stat_json%create_object(grid_section, 'latest')
                    call this%image_thumbnail_object(latest_classes_image, this%view_cls2D%thumbnail)
                    call this%stat_json%add(grid_section, latest_classes_image)
                    call this%stat_json%add(grid_section, "type", "thumbnail_grid")
                    call this%stat_json%add(cls2D, grid_section)
                end if
                !! pool rejected
                if(this%view_cls2D%pool_rejected_thumbnail%path .ne. "" .and. this%view_cls2D%pool_rejected_thumbnail%id > 0) then
                    call this%stat_json%create_object(pool_rejected_grid_section, 'pool_rejected_classes')
                    call this%image_thumbnail_object(pool_rejected_classes_image, this%view_cls2D%pool_rejected_thumbnail)
                    call this%stat_json%add(pool_rejected_grid_section, pool_rejected_classes_image)
                    call this%stat_json%add(pool_rejected_grid_section, "type", "thumbnail_grid")
                    call this%stat_json%add(cls2D, pool_rejected_grid_section)
                end if
                !! chunk rejected
                if(this%view_cls2D%chunk_rejected_thumbnail%path .ne. "" .and. this%view_cls2D%chunk_rejected_thumbnail%id > 0) then
                    call this%stat_json%create_object(chunk_rejected_grid_section, 'chunk_rejected_classes')
                    call this%image_thumbnail_object(chunk_rejected_classes_image, this%view_cls2D%chunk_rejected_thumbnail)
                    call this%stat_json%add(chunk_rejected_grid_section, chunk_rejected_classes_image)
                    call this%stat_json%add(chunk_rejected_grid_section, "type", "thumbnail_grid")
                    call this%stat_json%add(cls2D, chunk_rejected_grid_section)
                end if
                ! interactive plot
                call create_interactive_plot(this%view_cls2D%plot_ori, interactive_plot)
                call this%stat_json%add(cls2D, interactive_plot)
                ! add to sections
                call this%stat_json%add(views, cls2D)
                ! deallocate
                if(allocated(status_plot%labels))  deallocate(status_plot%labels)
                if(allocated(status_plot%colours)) deallocate(status_plot%colours)
                if(allocated(status_plot%data))    deallocate(status_plot%data)
            end subroutine add_view_cls2D

            subroutine add_view_ini3D()
                type(json_value), pointer :: ini3D
                type(json_value), pointer :: vols_section, stage, number_states, last_stage_completed, lp, vols, vol_json, fsc_json
                real                      :: fsc05, fsc0143
                integer                   :: noris, iori
                nullify(fsc_json)
                call this%stat_json%create_object(ini3D, 'ini3D')
                ! vols
                if(this%view_ini3D%stage > -1) then
                    call this%stat_json%create_object(vols_section, '3D')
                    call this%text_data_object(stage, "stage", this%view_ini3D%stage)
                    call this%stat_json%add(vols_section, stage)
                    if(trim(this%view_ini3D%last_stage_completed) .ne. "") then
                        call this%text_data_object(last_stage_completed, "last_stage_completed", this%view_ini3D%last_stage_completed)
                        call this%stat_json%add(vols_section, last_stage_completed)
                    end if
                    if(this%view_ini3D%number_states > -1) then
                        call this%text_data_object(number_states, "number_states", this%view_ini3D%number_states)
                        call this%stat_json%add(vols_section, number_states)
                    end if
                    if(this%view_ini3D%lp > 0.0) then
                        call this%text_data_object(lp, "low_pass_limit", this%view_ini3D%lp, dp=1)
                        call this%stat_json%add(vols_section, lp)
                    end if
                    noris = this%view_ini3D%vol_oris%get_noris() 
                    if( noris > 0 ) then
                        call this%stat_json%create_array(vols, 'vols')
                        do iori=1, noris
                            call this%view_ini3D%vol_oris%ori2json(iori, vol_json)
                            call this%fsc_object(iori, fsc_json,  this%view_ini3D%vol_oris,  this%view_ini3D%fsc_oris, fsc05, fsc0143 )
                            if(fsc05 .gt. 0.0) call this%stat_json%add(vol_json, 'fsc05',    dble(fsc05))
                            if(fsc05 .gt. 0.0) call this%stat_json%add(vol_json, 'fsc0143',  dble(fsc0143))
                            if(associated(fsc_json)) call this%stat_json%add(vol_json, fsc_json)
                            call this%stat_json%add(vols, vol_json)
                        end do
               !        if(present(hist)) call calculate_histogram(self%jobproc)
                        call this%stat_json%add(vols_section, vols)
                    end if
                    call this%stat_json%add(ini3D, vols_section)
                end if
                ! interactive plot
    !            call create_interactive_plot(this%view_vols%plot_ori, interactive_plot)
    !            call this%stat_json%add(vols, interactive_plot)
                ! add to sections
                call this%stat_json%add(views, ini3D)
            end subroutine add_view_ini3D

            subroutine add_view_vols()
                type(json_value), pointer :: vols
                type(json_value), pointer :: vols_section, number_vols, interactive_plot
                call this%stat_json%create_object(vols, 'vols')
                ! vols
                if(this%view_vols%number_vols > -1) then
                    call this%stat_json%create_object(vols_section, 'volumes')
                    if(this%view_vols%number_vols > -1) then
                        call this%text_data_object(number_vols, "number_volumes", this%view_vols%number_vols)
                        call this%stat_json%add(vols_section, number_vols)
                    end if
                    call this%stat_json%add(vols, vols_section)
                end if
                ! interactive plot
                call create_interactive_plot(this%view_vols%plot_ori, interactive_plot)
                call this%stat_json%add(vols, interactive_plot)
                ! add to sections
                call this%stat_json%add(views, vols)
            end subroutine add_view_vols


    end subroutine generate_stat_json

    subroutine fsc_object( this, iori_l, fsc_json, vol_oris, fsc_oris, fsc05, fsc0143)
        class(simple_nice_communicator), intent(inout) :: this
        integer,                         intent(in)    :: iori_l
        type(json_value), pointer,       intent(inout) :: fsc_json
        type(oris),                      intent(in)    :: vol_oris, fsc_oris
        real,                            intent(out)   :: fsc05, fsc0143
        type(json_value), pointer     :: datasets, dataset, data, labels !fsc_json
        type(string)                  :: fscfile
        real,             allocatable :: fsc(:), res(:)
        real                          :: smpd_l, box_l
        integer                       :: ifsc
        logical                       :: fsc05_crossed, fsc0143_crossed
        if(.not. vol_oris%get_noris() .eq. fsc_oris%get_noris()) return
        if(.not. vol_oris%isthere(iori_l, "smpd")) return
        if(.not. fsc_oris%isthere(iori_l, "fsc")) return
        if(.not. fsc_oris%isthere(iori_l, "box")) return
        call fsc_oris%getter(iori_l, "fsc", fscfile)
        smpd_l = vol_oris%get(iori_l, "smpd")
        box_l  = fsc_oris%get(iori_l, "box")
        if(.not. file_exists(fscfile)) return
        fsc = file2rarr(fscfile)
        res = get_resarr(int(box_l), smpd_l)
        call this%stat_json%create_object(fsc_json, 'fsc')
        call this%stat_json%add(fsc_json, 'type', "plot_bar")
        call this%stat_json%create_array(datasets, "datasets")
        call this%stat_json%create_array(data,     "data")
        call this%stat_json%create_array(labels,   "labels")
        call this%stat_json%create_object(dataset, "dataset")
        fsc05_crossed   = .false.
        fsc0143_crossed = .false.
        do ifsc=1, size(fsc)
            if(.not. fsc05_crossed) then
                if(fsc(ifsc) .gt. 0.5) then
                    fsc05 = res(ifsc)
                else
                    fsc05_crossed = .true.
                end if
            end if
            if(.not. fsc0143_crossed) then
                if(fsc(ifsc) .gt. 0.143) then
                    fsc0143 = res(ifsc)
                else
                    fsc0143_crossed = .true.
                end if
            end if
            call this%stat_json%add(data,   '', dble(fsc(ifsc)))
            call this%stat_json%add(labels, '', dble(res(ifsc)))
        end do
        call this%stat_json%add(dataset, 'borderColor', "rgba(30, 144, 255, 0.5)")
        call this%stat_json%add(dataset, 'pointStyle', .false.)
        call this%stat_json%add(dataset, 'cubicInterpolationMode', 'monotone')
        call this%stat_json%add(dataset, 'tension', dble(0.4))
        call this%stat_json%add(dataset, data)
        call this%stat_json%add(datasets, dataset)
        call this%stat_json%add(fsc_json, datasets)
        call this%stat_json%add(fsc_json, labels)
        call fscfile%kill
    end subroutine fsc_object

    subroutine update_micrographs(this, movies_imported, movies_processed, last_movie_imported, micrographs, micrographs_rejected, compute_in_use,&
        avg_ctf_resolution, avg_ice_score, avg_astigmatism, thumbnail, thumbnail_id, thumbnail_static_id, carousel, clear_carousel)
        class(simple_nice_communicator), intent(inout) :: this
        integer,               optional, intent(in)    :: movies_imported, movies_processed, micrographs, micrographs_rejected, compute_in_use, thumbnail_id, thumbnail_static_id
        real,                  optional, intent(in)    :: avg_ctf_resolution, avg_ice_score, avg_astigmatism
        logical,               optional, intent(in)    :: last_movie_imported, carousel, clear_carousel
        character(len=*),      optional, intent(in)    :: thumbnail
        type(nice_stat_thumb_image)                    :: new_thumbnail
        integer :: uid, i
        real    :: rnd
        logical :: new
        this%view_micrographs%active = .true.
        if(present(movies_imported))      this%view_micrographs%movies_imported      = movies_imported
        if(present(movies_processed))     this%view_micrographs%movies_processed     = movies_processed
        if(present(micrographs))          this%view_micrographs%micrographs          = micrographs
        if(present(micrographs_rejected)) this%view_micrographs%micrographs_rejected = micrographs_rejected
        if(present(compute_in_use))       this%view_micrographs%compute_in_use       = compute_in_use
        if(present(avg_ctf_resolution))   this%view_micrographs%avg_ctf_resolution   = avg_ctf_resolution
        if(present(avg_ice_score))        this%view_micrographs%avg_ice_score        = avg_ice_score
        if(present(avg_astigmatism))      this%view_micrographs%avg_astigmatism      = avg_astigmatism
        if(present(last_movie_imported) .and. last_movie_imported)  this%view_micrographs%last_movie_imported  = datestr()
        if(present(thumbnail) .and. present(thumbnail_id) .and. present(thumbnail_static_id)) then
            call seed_rnd()
            call random_number(rnd)
            uid = floor(1000000 * rnd)
            if(present(carousel)) then
                if(carousel) then
                    if(.not. allocated(this%view_micrographs%thumbnail_carousel))      allocate(this%view_micrographs%thumbnail_carousel(0))
                    if(.not. allocated(this%view_micrographs%thumbnail_carousel_mask)) allocate(this%view_micrographs%thumbnail_carousel_mask(0))
                    if(present(clear_carousel)) then
                        if(clear_carousel) this%view_micrographs%thumbnail_carousel_mask = .false.
                    end if
                    new = .true.
                    do i=1, size(this%view_micrographs%thumbnail_carousel)
                        if(this%view_micrographs%thumbnail_carousel(i)%static_id == thumbnail_static_id &
                        .and. this%view_micrographs%thumbnail_carousel(i)%id == thumbnail_id) then
                            new = .false.
                            this%view_micrographs%thumbnail_carousel_mask(i) = .true.
                        end if
                    end do
                    if(new) then
                        new_thumbnail%path      = thumbnail
                        new_thumbnail%id        = thumbnail_id
                        new_thumbnail%static_id = thumbnail_static_id
                        new_thumbnail%uid       = "U" // int2str(uid)
                        this%view_micrographs%thumbnail_carousel      = [this%view_micrographs%thumbnail_carousel,      new_thumbnail]
                        this%view_micrographs%thumbnail_carousel_mask = [this%view_micrographs%thumbnail_carousel_mask, .true.]
                    end if
                else 
                    if(thumbnail_id .ne. this%view_micrographs%thumbnail%id) then
                        this%view_micrographs%thumbnail%path      = thumbnail
                        this%view_micrographs%thumbnail%id        = thumbnail_id
                        this%view_micrographs%thumbnail%static_id = thumbnail_static_id
                        this%view_micrographs%thumbnail%uid       = "U" // int2str(uid)
                    end if
                end if
            end if
        end if
    end subroutine update_micrographs

    subroutine update_optics(this, micrographs, assigned, last_micrograph_imported)
        class(simple_nice_communicator), intent(inout) :: this
        integer,               optional, intent(in)    :: micrographs, assigned
        logical,               optional, intent(in)    :: last_micrograph_imported
        this%view_optics%active = .true.
        if(present(micrographs)) this%view_optics%micrographs  = micrographs
        if(present(assigned))    this%view_optics%opticsgroups = assigned
        if(present(last_micrograph_imported)) then
            if(last_micrograph_imported) this%view_optics%last_micrograph_imported  = datestr()
        end if
    end subroutine update_optics

    subroutine update_pick(this, micrographs_imported, micrographs_rejected, micrographs_picked, last_micrograph_imported, thumbnail, thumbnail_id, &
    thumbnail_static_id, carousel, clear_carousel, boxfile, scale, gaussian_diameter, suggested_diameter, pickrefs_thumbnail, pickrefs_thumbnail_id, &
    pickrefs_thumbnail_n_tiles, pickrefs_thumbnail_static_id, pickrefs_thumbnail_scale)
        class(simple_nice_communicator), intent(inout) :: this
        integer,               optional, intent(in)    :: micrographs_imported, micrographs_rejected, micrographs_picked, thumbnail_id, thumbnail_static_id, gaussian_diameter, suggested_diameter
        integer,               optional, intent(in)    :: pickrefs_thumbnail_id, pickrefs_thumbnail_static_id, pickrefs_thumbnail_n_tiles
        logical,               optional, intent(in)    :: last_micrograph_imported, carousel, clear_carousel
        character(len=*),      optional, intent(in)    :: thumbnail, pickrefs_thumbnail, boxfile
        real,                  optional, intent(in)    :: scale, pickrefs_thumbnail_scale
        type(nice_stat_thumb_image)                    :: new_thumbnail
        integer :: uid, i
        real    :: rnd
        logical :: new
        this%view_pick%active = .true.
        if(present(micrographs_imported)) this%view_pick%micrographs_imported = micrographs_imported
        if(present(micrographs_rejected)) this%view_pick%micrographs_rejected = micrographs_rejected
        if(present(micrographs_picked))   this%view_pick%micrographs_picked   = micrographs_picked
        if(present(gaussian_diameter))    this%view_pick%gaussian_diameter    = gaussian_diameter
        if(present(suggested_diameter))   this%view_pick%suggested_diameter   = suggested_diameter
        if(present(last_micrograph_imported)) then
            if(last_micrograph_imported)  this%view_pick%last_micrograph_imported  = datestr()
        end if
        if(present(thumbnail) .and. present(thumbnail_id) .and. present(thumbnail_static_id)) then
            call seed_rnd()
            call random_number(rnd)
            uid = floor(1000000 * rnd)
            if(present(carousel)) then
                if(carousel) then
                    if(.not. allocated(this%view_pick%thumbnail_carousel))      allocate(this%view_pick%thumbnail_carousel(0))
                    if(.not. allocated(this%view_pick%thumbnail_carousel_mask)) allocate(this%view_pick%thumbnail_carousel_mask(0))
                    if(present(clear_carousel)) then
                        if(clear_carousel) this%view_pick%thumbnail_carousel_mask = .false.
                    end if
                    new = .true.
                    do i=1, size(this%view_pick%thumbnail_carousel)
                        if(this%view_pick%thumbnail_carousel(i)%static_id == thumbnail_static_id &
                        .and. this%view_pick%thumbnail_carousel(i)%id == thumbnail_id) then
                            new = .false.
                            this%view_pick%thumbnail_carousel_mask(i) = .true.
                        end if
                    end do
                    if(new) then
                        new_thumbnail%path      = thumbnail
                        new_thumbnail%id        = thumbnail_id
                        new_thumbnail%static_id = thumbnail_static_id
                        new_thumbnail%uid       = "U" // int2str(uid)
                        if(present(boxfile)) new_thumbnail%boxfile = boxfile
                        if(present(scale))   new_thumbnail%scale   = scale
                        this%view_pick%thumbnail_carousel      = [this%view_pick%thumbnail_carousel,      new_thumbnail]
                        this%view_pick%thumbnail_carousel_mask = [this%view_pick%thumbnail_carousel_mask, .true.]
                    end if
                else
                    if(thumbnail_id .ne. this%view_pick%thumbnail%id) then
                        this%view_pick%thumbnail%path      = thumbnail
                        this%view_pick%thumbnail%id        = thumbnail_id
                        this%view_pick%thumbnail%static_id = thumbnail_static_id
                        this%view_pick%thumbnail%uid       = "U" // int2str(uid)
                        if(present(boxfile)) this%view_pick%thumbnail%boxfile = boxfile
                        if(present(scale))   this%view_pick%thumbnail%scale   = scale
                    end if
                end if
            end if
        end if
        if(present(pickrefs_thumbnail) .and. present(pickrefs_thumbnail_id) .and. present(pickrefs_thumbnail_static_id) .and. present(pickrefs_thumbnail_n_tiles)) then
            if(pickrefs_thumbnail_id .ne. this%view_pick%pickrefs_thumbnail%id) then
                call seed_rnd()
                call random_number(rnd)
                uid = floor(1000000 * rnd)
                this%view_pick%pickrefs_thumbnail%path      = pickrefs_thumbnail
                this%view_pick%pickrefs_thumbnail%id        = pickrefs_thumbnail_id
                this%view_pick%pickrefs_thumbnail%static_id = pickrefs_thumbnail_static_id
                this%view_pick%pickrefs_thumbnail%uid       = "U" // int2str(uid)
                this%view_pick%pickrefs_thumbnail%grid      = .true.
                this%view_pick%pickrefs_thumbnail%spriten   = pickrefs_thumbnail_n_tiles
                this%view_pick%pickrefs_thumbnail%spritedim = JPEG_DIM
                if(present(pickrefs_thumbnail_scale))   this%view_pick%pickrefs_thumbnail%scale = pickrefs_thumbnail_scale
                if(allocated(this%view_pick%pickrefs_thumbnail%metadata)) deallocate(this%view_pick%pickrefs_thumbnail%metadata)
                allocate(this%view_pick%pickrefs_thumbnail%metadata(0))
            end if
        end if
    end subroutine update_pick

    subroutine update_cls2D(this, particles_extracted, particles_imported, last_particles_imported, iteration, number_classes, number_classes_rejected, &
    number_particles_assigned, number_particles_rejected, maximum_resolution, last_iteration, thumbnail, thumbnail_id, thumbnail_static_id, &
    thumbnail_n_tiles, stats_mask, stats_resolution, stats_population, scale, pool_rejected_thumbnail, pool_rejected_thumbnail_id, pool_rejected_thumbnail_static_id, &
    pool_rejected_thumbnail_n_tiles, pool_rejected_scale, chunk_rejected_thumbnail, chunk_rejected_thumbnail_id, chunk_rejected_thumbnail_static_id, &
    chunk_rejected_thumbnail_n_tiles, chunk_rejected_scale, snapshot_id, snapshot_time, rejection_params, snapshot_json, snapshot_json_clear)
        class(simple_nice_communicator), intent(inout) :: this
        integer,               optional, intent(in)    :: particles_extracted, particles_imported, iteration, number_classes, number_classes_rejected, thumbnail_n_tiles
        integer,               optional, intent(in)    :: number_particles_assigned, number_particles_rejected, thumbnail_id, thumbnail_static_id, pool_rejected_thumbnail_id
        integer,               optional, intent(in)    :: pool_rejected_thumbnail_static_id, pool_rejected_thumbnail_n_tiles, chunk_rejected_thumbnail_id
        integer,               optional, intent(in)    :: chunk_rejected_thumbnail_static_id, chunk_rejected_thumbnail_n_tiles, snapshot_id
        character(len=*),      optional, intent(in)    :: thumbnail, last_iteration, pool_rejected_thumbnail, chunk_rejected_thumbnail, snapshot_time
        logical,               optional, intent(in)    :: last_particles_imported, snapshot_json_clear
        real,                  optional, intent(in)    :: maximum_resolution
        real,                  optional, intent(in)    :: stats_mask(:), stats_population(:), stats_resolution(:)
        real,                  optional, intent(in)    :: scale, pool_rejected_scale, chunk_rejected_scale
        real,                  optional, intent(in)    :: rejection_params(3)
        type(json_value), pointer, optional, intent(in) :: snapshot_json
        logical, allocatable                           :: msk(:)
        integer, allocatable                           :: id(:)
        type(nice_stat_thumb_image_meta)               :: meta_res
        integer :: uid, i
        real    :: rnd
        this%view_cls2D%active = .true.
        if(present(particles_extracted))       this%view_cls2D%particles_extracted       = particles_extracted
        if(present(particles_imported))        this%view_cls2D%particles_imported        = particles_imported
        if(present(iteration))                 this%view_cls2D%iteration                 = iteration
        if(present(number_classes))            this%view_cls2D%number_classes            = number_classes
        if(present(number_classes_rejected))   this%view_cls2D%number_classes_rejected   = number_classes_rejected
        if(present(number_particles_assigned)) this%view_cls2D%number_particles_assigned = number_particles_assigned
        if(present(number_particles_rejected)) this%view_cls2D%number_particles_rejected = number_particles_rejected
        if(present(maximum_resolution))        this%view_cls2D%maximum_resolution        = maximum_resolution
        if(present(snapshot_id))               this%view_cls2D%snapshot_id               = snapshot_id
        if(present(snapshot_time)) then
            if(trim(snapshot_time) .ne. "") this%view_cls2D%snapshot_time = trim(snapshot_time)
        end if
        if(present(last_iteration)) then
            if(trim(last_iteration) .ne. "") this%view_cls2D%last_iteration = trim(last_iteration)
        end if
        if(present(last_particles_imported)) then
            if (last_particles_imported) this%view_cls2D%last_particles_imported = datestr()
        end if
        if(present(thumbnail) .and. present(thumbnail_id) .and. present(thumbnail_static_id) .and. present(thumbnail_n_tiles)) then
            if(thumbnail_id .ne. this%view_cls2D%thumbnail%id) then
                call seed_rnd()
                call random_number(rnd)
                uid = floor(1000000 * rnd)
                this%view_cls2D%thumbnail%path      = thumbnail
                this%view_cls2D%thumbnail%id        = thumbnail_id
                this%view_cls2D%thumbnail%static_id = thumbnail_static_id
                this%view_cls2D%thumbnail%uid       = "U" // int2str(uid)
                this%view_cls2D%thumbnail%grid      = .true.
                this%view_cls2D%thumbnail%spriten   = thumbnail_n_tiles
                this%view_cls2D%thumbnail%spritedim = JPEG_DIM
                if(present(scale))   this%view_cls2D%thumbnail%scale = scale
                if(allocated(this%view_cls2D%thumbnail%metadata)) deallocate(this%view_cls2D%thumbnail%metadata)
                allocate(this%view_cls2D%thumbnail%metadata(0))
            end if
        end if
        if(present(pool_rejected_thumbnail) .and. present(pool_rejected_thumbnail_id) .and. present(pool_rejected_thumbnail_static_id) .and. present(pool_rejected_thumbnail_n_tiles)) then
            if(pool_rejected_thumbnail_id .ne. this%view_cls2D%pool_rejected_thumbnail%id) then
                call seed_rnd()
                call random_number(rnd)
                uid = floor(1000000 * rnd)
                this%view_cls2D%pool_rejected_thumbnail%path      = pool_rejected_thumbnail
                this%view_cls2D%pool_rejected_thumbnail%id        = pool_rejected_thumbnail_id
                this%view_cls2D%pool_rejected_thumbnail%static_id = pool_rejected_thumbnail_static_id
                this%view_cls2D%pool_rejected_thumbnail%uid       = "U" // int2str(uid)
                this%view_cls2D%pool_rejected_thumbnail%grid      = .true.
                this%view_cls2D%pool_rejected_thumbnail%spriten   = pool_rejected_thumbnail_n_tiles
                this%view_cls2D%pool_rejected_thumbnail%spritedim = JPEG_DIM
                if(present(pool_rejected_scale))   this%view_cls2D%pool_rejected_thumbnail%scale = pool_rejected_scale
                if(allocated(this%view_cls2D%pool_rejected_thumbnail%metadata)) deallocate(this%view_cls2D%pool_rejected_thumbnail%metadata)
                allocate(this%view_cls2D%pool_rejected_thumbnail%metadata(0))
            end if
        end if
        if(present(chunk_rejected_thumbnail) .and. present(chunk_rejected_thumbnail_id) .and. present(chunk_rejected_thumbnail_static_id) .and. present(chunk_rejected_thumbnail_n_tiles)) then
            if(chunk_rejected_thumbnail_id .ne. this%view_cls2D%chunk_rejected_thumbnail%id) then
                call seed_rnd()
                call random_number(rnd)
                uid = floor(1000000 * rnd)
                this%view_cls2D%chunk_rejected_thumbnail%path      = chunk_rejected_thumbnail
                this%view_cls2D%chunk_rejected_thumbnail%id        = chunk_rejected_thumbnail_id
                this%view_cls2D%chunk_rejected_thumbnail%static_id = chunk_rejected_thumbnail_static_id
                this%view_cls2D%chunk_rejected_thumbnail%uid       = "U" // int2str(uid)
                this%view_cls2D%chunk_rejected_thumbnail%grid      = .true.
                this%view_cls2D%chunk_rejected_thumbnail%spriten   = chunk_rejected_thumbnail_n_tiles
                this%view_cls2D%chunk_rejected_thumbnail%spritedim = JPEG_DIM
                if(present(chunk_rejected_scale))   this%view_cls2D%chunk_rejected_thumbnail%scale = chunk_rejected_scale
                if(allocated(this%view_cls2D%chunk_rejected_thumbnail%metadata)) deallocate(this%view_cls2D%chunk_rejected_thumbnail%metadata)
                allocate(this%view_cls2D%chunk_rejected_thumbnail%metadata(0))
            end if
        end if
        if(present(stats_mask)) then
            allocate(id(size(stats_mask)))
            do i = 1, size(stats_mask)
                id(i) = i
            end do
            msk = stats_mask > 0
            meta_res%label = "id"
            meta_res%type  = 1
            meta_res%data_int = pack(id, msk)
            if(.not. allocated(this%view_cls2D%thumbnail%metadata)) allocate(this%view_cls2D%thumbnail%metadata(0))
            this%view_cls2D%thumbnail%metadata = [meta_res, this%view_cls2D%thumbnail%metadata]
        end if
        if(present(stats_mask) .and. present(stats_resolution)) then
            if(size(stats_mask) .eq. size(stats_resolution)) then
                msk = stats_mask > 0
                meta_res%label = "res"
                meta_res%type  = 2
                meta_res%data_real = pack(stats_resolution, msk)
                if(.not. allocated(this%view_cls2D%thumbnail%metadata)) allocate(this%view_cls2D%thumbnail%metadata(0))
                this%view_cls2D%thumbnail%metadata = [meta_res, this%view_cls2D%thumbnail%metadata]
            end if
        end if
        if(present(stats_mask) .and. present(stats_population)) then
            if(size(stats_mask) .eq. size(stats_population)) then
                msk = stats_mask > 0
                meta_res%label = "pop"
                meta_res%type  = 1
                meta_res%data_int = pack(stats_population, msk)
                if(.not. allocated(this%view_cls2D%thumbnail%metadata)) allocate(this%view_cls2D%thumbnail%metadata(0))
                this%view_cls2D%thumbnail%metadata = [meta_res, this%view_cls2D%thumbnail%metadata]
            end if
        end if
        if(present(rejection_params)) this%view_cls2D%rejection_params = rejection_params
        if(present(snapshot_json)) then
            if(associated(snapshot_json)) then
                if(associated(this%view_cls2D%snapshot_json)) then
                    call this%stat_json%destroy(this%view_cls2D%snapshot_json)
                    nullify(this%view_cls2D%snapshot_json)
                end if
                call this%stat_json%clone(snapshot_json, this%view_cls2D%snapshot_json)
            end if
        end if
        if(present(snapshot_json_clear)) then
            if(snapshot_json_clear) then
                if(associated(this%view_cls2D%snapshot_json)) then
                    call this%stat_json%destroy(this%view_cls2D%snapshot_json)
                    nullify(this%view_cls2D%snapshot_json)
                end if
            end if
        end if
    end subroutine update_cls2D

    subroutine update_ini3D(this, stage, number_states, last_stage_completed, lp, vol_oris, fsc_oris)
        class(simple_nice_communicator),           intent(inout) :: this
        type(oris),                      optional, intent(in)    :: vol_oris, fsc_oris
        integer,                         optional, intent(in)    :: stage, number_states
        logical,                         optional, intent(in)    :: last_stage_completed
        real,                            optional, intent(in)    :: lp
        this%view_ini3D%active = .true.
        if(present(stage))         this%view_ini3D%stage         = stage
        if(present(number_states)) this%view_ini3D%number_states = number_states
        if(present(lp))            this%view_ini3D%lp            = lp
        if(present(last_stage_completed)) then
            if(last_stage_completed) this%view_ini3D%last_stage_completed = datestr()
        end if
        if(present(vol_oris)) then
            call this%view_ini3D%vol_oris%kill()
            this%view_ini3D%vol_oris = vol_oris
        end if
        if(present(fsc_oris)) then
            call this%view_ini3D%fsc_oris%kill()
            this%view_ini3D%fsc_oris = fsc_oris
        end if
    end subroutine update_ini3D

    subroutine update_vols(this, number_vols)
        class(simple_nice_communicator),           intent(inout) :: this
        integer,                         optional, intent(in)    :: number_vols
        this%view_vols%active = .true.
        if(present(number_vols)) this%view_vols%number_vols = number_vols
    end subroutine update_vols

    subroutine text_data_object_1(this, ptr, key, val)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        character(*),                    intent(in)    :: key
        integer,                         intent(in)    :: val
        call this%stat_json%create_object(ptr, key)
        call this%stat_json%add(ptr, "value",  val)
        call this%stat_json%add(ptr, "type",  "integer")
    end subroutine text_data_object_1

    subroutine text_data_object_2(this, ptr, key, val, dp)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        character(*),                    intent(in)    :: key
        real,                            intent(in)    :: val
        integer, optional,               intent(in)    :: dp
        call this%stat_json%create_object(ptr, key)
        call this%stat_json%add(ptr, "value",  dble(val))
        call this%stat_json%add(ptr, "type",  "float")
        if(present(dp)) call this%stat_json%add(ptr, "dp",  dp)
    end subroutine text_data_object_2

    subroutine text_data_object_3(this, ptr, key, val)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        character(*),                    intent(in)    :: key
        character(*),                    intent(in)    :: val
        call this%stat_json%create_object(ptr, key)
        call this%stat_json%add(ptr, "value",  val)
        call this%stat_json%add(ptr, "type",  "string")
    end subroutine text_data_object_3

    subroutine image_thumbnail_object(this, ptr, thumb_image, name)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        type(nice_stat_thumb_image),     intent(in)    :: thumb_image
        character(*),     optional,      intent(in)    :: name
        type(json_value), pointer                      :: metadata, metadata_obj
        integer                                        :: i
        if(present(name)) then
            call this%stat_json%create_object(ptr, name)
        else
            call this%stat_json%create_object(ptr, "thumbnail")
        end if
        call this%stat_json%add(ptr, "id",        thumb_image%id)
        call this%stat_json%add(ptr, "static_id", thumb_image%static_id)
        call this%stat_json%add(ptr, "path",      thumb_image%path)
        call this%stat_json%add(ptr, "uid",       thumb_image%uid)
        if(thumb_image%scale .gt. 0.0) call this%stat_json%add(ptr, "scale", dble(thumb_image%scale))
        if(thumb_image%grid .and. thumb_image%spritedim .gt. 0 .and. thumb_image%spriten .gt. 0) then
            call this%stat_json%add(ptr, "type",      "thumbnail_grid")
            call this%stat_json%add(ptr, "spritedim", thumb_image%spritedim)
            call this%stat_json%add(ptr, "spriten",   thumb_image%spriten)
        else
            call this%stat_json%add(ptr, "type",      "thumbnail")
        end if
        if(allocated(thumb_image%boxfile)) call this%stat_json%add(ptr, "boxfile", thumb_image%boxfile)
        if(allocated(thumb_image%metadata)) then
            call this%stat_json%create_object(metadata, "metadata")
            do i=1, size(thumb_image%metadata)
                call this%stat_json%create_object(metadata_obj, trim(thumb_image%metadata(i)%label))
                if(thumb_image%metadata(i)%type == 1) then
                    ! integer
                    call this%stat_json%add(metadata_obj, "type", "integer")
                    call this%stat_json%add(metadata_obj, "data", thumb_image%metadata(i)%data_int)
                else if(thumb_image%metadata(i)%type == 2) then
                    ! real
                    call this%stat_json%add(metadata_obj, "type", "float")
                    call this%stat_json%add(metadata_obj, "dp",   thumb_image%metadata(i)%dp)
                    call this%stat_json%add(metadata_obj, "data", dble(thumb_image%metadata(i)%data_real))
                end if
                call this%stat_json%add(metadata, metadata_obj)
            end do
            call this%stat_json%add(ptr, metadata)
        end if
    end subroutine image_thumbnail_object

    subroutine doughnut_plot_object(this, ptr, plot, name, value)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        type(nice_plot_doughnut),        intent(in)    :: plot
        character(*),     optional,      intent(in)    :: name, value
        type(json_value),                pointer       :: data, colours, labels, datasets, dataset
        integer                                        :: i
        if(present(name)) then
            call this%stat_json%create_object(ptr, name)
        else
            call this%stat_json%create_object(ptr, "plot")
        end if
        call this%stat_json%add(ptr, 'type', "plot_doughnut")
        call this%stat_json%create_array(datasets, "datasets")
        call this%stat_json%create_array(data,     "data")
        call this%stat_json%create_array(colours,  "backgroundColor")
        call this%stat_json%create_array(labels,   "labels")
        do i=1, size(plot%data)
            call this%stat_json%add(data,    '', plot%data(i))
            call this%stat_json%add(colours, '', trim(plot%colours(i)))
            call this%stat_json%add(labels,  '', trim(plot%labels(i)))
        end do
        call this%stat_json%create_object(dataset, "dataset")
        call this%stat_json%add(dataset, data)
        call this%stat_json%add(dataset, colours)
        call this%stat_json%add(datasets, dataset)
        call this%stat_json%add(ptr, datasets)
        call this%stat_json%add(ptr, labels)
    end subroutine doughnut_plot_object

    subroutine bar_plot_object(this, ptr, plot, name, value, extra_1)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        type(nice_plot_bar),             intent(in)    :: plot
        character(*),     optional,      intent(in)    :: name, value, extra_1
        type(json_value),                pointer       :: data, labels, datasets, dataset
        integer                                        :: i
        if(present(name)) then
            call this%stat_json%create_object(ptr, name)
        else
            call this%stat_json%create_object(ptr, "plot")
        end if
        call this%stat_json%add(ptr, 'type', "plot_bar")
        call this%stat_json%create_array(datasets, "datasets")
        call this%stat_json%create_array(data,     "data")
        call this%stat_json%create_array(labels,   "labels")
        do i=1, size(plot%data)
            call this%stat_json%add(data,    '', plot%data(i))
            call this%stat_json%add(labels,  '', trim(plot%labels(i)))
        end do
        call this%stat_json%create_object(dataset, "dataset")
        call this%stat_json%add(dataset, 'backgroundColor', trim(plot%colours(1)))
        call this%stat_json%add(dataset, data)
        call this%stat_json%add(datasets, dataset)
        call this%stat_json%add(ptr, datasets)
        call this%stat_json%add(ptr, labels)
        if(present(value)) call this%stat_json%add(ptr, "value", value)
        if(present(extra_1)) call this%stat_json%add(ptr, "extra_1", extra_1)
    end subroutine bar_plot_object

    function calculate_checksum(this, str)
        class(simple_nice_communicator), intent(inout)   :: this
        character(kind=CK,len=:),allocatable, intent(in) :: str
        integer       :: calculate_checksum, i
        calculate_checksum = 0
        do i = 1, len(str)
            calculate_checksum = calculate_checksum + ichar(str(i:i))
        end do
    end function calculate_checksum

    function datestr() 
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        character(16) :: datestr
        integer,dimension(8) :: values
        ! using keyword arguments
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        write(datestr, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') values(1), '/', values(2), '/', values(3), '_', values(5), ':', values(6)
    end function datestr

    subroutine update_from_project(this, spproj)
        class(simple_nice_communicator), intent(inout) :: this
        type(sp_project),                intent(inout) :: spproj
        type(oris)                                     :: vol_oris
        logical,            allocatable                :: state_mask(:)
        integer,            allocatable                :: pinds(:)
        if(spproj%os_mic%get_noris() .gt. 0) then
            if(allocated(state_mask)) deallocate(state_mask)
            if(allocated(pinds))      deallocate(pinds)
            call spproj%os_mic%mask_from_state(1, state_mask, pinds)
            call spproj%os_mic%get_ori(1, this%view_micrographs%plot_ori)
            call this%update_micrographs(micrographs=size(state_mask), micrographs_rejected=size(state_mask)-count(state_mask))
        end if
        if(spproj%os_cls2D%get_noris() .gt. 0) then
            if(allocated(state_mask)) deallocate(state_mask)
            if(allocated(pinds))      deallocate(pinds)
            call spproj%os_cls2D%mask_from_state(1, state_mask, pinds)
            call spproj%os_cls2D%get_ori(1, this%view_cls2D%plot_ori)
            call this%update_cls2D(number_classes=size(state_mask), number_classes_rejected=size(state_mask)-count(state_mask))
        end if
        call spproj%get_all_vols(vol_oris)
        if(vol_oris%get_noris() .gt. 0) then
            if(allocated(state_mask)) deallocate(state_mask)
            if(allocated(pinds))      deallocate(pinds)
            call vol_oris%mask_from_state(1, state_mask, pinds)
            call vol_oris%get_ori(1, this%view_vols%plot_ori)
            call this%update_vols(number_vols=size(state_mask))
        end if
        call vol_oris%kill
        if(allocated(state_mask)) deallocate(state_mask)
        if(allocated(pinds))      deallocate(pinds)
    end subroutine update_from_project

    subroutine get_real_keys_json(this, ptr, seg_ori)
        class(simple_nice_communicator), intent(inout) :: this
        type(json_value), pointer,       intent(inout) :: ptr
        type(ori),                       intent(in)    :: seg_ori
        type(string), allocatable :: keys(:)
        integer :: j
        call this%stat_json%create_array(ptr, 'keys')
        keys = seg_ori%get_keys()
        do j = 1, size(keys)
            if(.not. seg_ori%ischar(keys(j)%to_char())) then
                call this%stat_json%add(ptr, '', keys(j)%to_char())
            end if
        end do
        call this%stat_json%add(ptr, '', 'n')
    end subroutine get_real_keys_json

    subroutine get_stat_json(this, json)
        class(simple_nice_communicator),       intent(inout) :: this
        type(json_value),             pointer, intent(inout) :: json
        call this%generate_stat_json()
        call this%stat_json%clone(this%stat_json_root, json)
    end subroutine get_stat_json

end module simple_nice
module simple_commander_abinitio3D_stream
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: params_glob
use simple_sp_project,         only: sp_project
use simple_qsys_env,           only: qsys_env
use simple_oris,               only: oris
use simple_stream_utils
implicit none

public :: init_abinitio3D, request_snapshot, add_projects2jobslist, submit_jobs
public :: analysis, update_user_params3D, restart_cleanup, get_nrecs_submitted
public :: check_processes

private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: USER_PARAMS3D = 'stream3D_user_params.txt'
character(len=STDLEN), parameter :: EXECDIR       = 'run_'
character(len=STDLEN), parameter :: DIR_ALGNVOLS  = 'alignedvolumes/'

type(procrecord),    allocatable :: procrecs(:)             ! records of all jobs
real,                allocatable :: algnvols_corrmat(:,:)   ! Volumes alignment correlation matrix
integer                          :: nvols_glob = 0          ! Total number of volumes generated (runs completed)
integer(kind=dp)                 :: time_last_submission

contains

    subroutine init_abinitio3D
        call cleanup_abinitio3D
        ! For volumes alignment & analysis
        call simple_mkdir(DIR_ALGNVOLS)
        params_glob%gridding  = 'yes'
        params_glob%interpfun = 'kb'
    end subroutine init_abinitio3D

    subroutine cleanup_abinitio3D
        nvols_glob = 0
        time_last_submission = 0
        call kill_procrecords(procrecs)
        if(allocated(algnvols_corrmat))deallocate(algnvols_corrmat)
    end subroutine cleanup_abinitio3D

    ! request snapshot generation to 2D process
    ! nptcls: # of particles that will be used in next run of abinitio3D
    subroutine request_snapshot( nptcls )
        integer, intent(in) :: nptcls
        type(oris) :: os
        call os%new(1, is_ptcl=.false.)
        call os%set(1,'nptcls', nptcls)
        call os%write(trim(params_glob%dir_target)//'/'//trim(SNAPSHOT_REQUEST))
        call os%kill
        write(logfhandle,'(A,A)')'>>> REQUESTED ONE SNAPSHOT: ',cast_time_char(simple_gettime())
    end subroutine request_snapshot

    ! Append each project file to the processing queue, no submission performed
    subroutine add_projects2jobslist( projects )
        use simple_procimgstk, only: scale_imgfile
        use simple_class_frcs
        character(len=LONGSTRLEN), allocatable, intent(in) :: projects(:)
        character(len=LONGSTRLEN), allocatable :: files(:)
        character(len=:),          allocatable :: projname, projfile, direxec, id, string
        character(len=:),          allocatable :: path, frc, cavgs, tmpl, tmp1, tmp2
        type(sp_project)          :: spproj
        type(class_frcs)          :: frcs
        character(len=LONGSTRLEN) :: project_fname
        real    :: smpd, smpd_sc
        integer :: i,j,n, jobid, box
        logical :: found
        if( .not.allocated(projects) )return
        do i = 1,size(projects)
            project_fname = trim(projects(i))
            if( .not.file_exists(project_fname) )then
                THROW_WARN('Could not find project: '//trim(project_fname))
                cycle
            endif
            ! prep record & folder
            jobid    = get_nrecs() + 1
            id       = int2str_pad(jobid, params_glob%numlen)
            direxec  = trim(EXECDIR)//id//'/'       ! relative path
            projname = 'abinitio_'  //id
            projfile = projname//trim(METADATA_EXT)
            call simple_mkdir(direxec)
            direxec  = simple_abspath(direxec)//'/' ! absolute path & trailing /
            ! update project
            call spproj%read(project_fname)
            if( spproj%projinfo%get_noris() /= 1) THROW_HARD('Critical error add_projects2jobslist')
            call spproj%projinfo%set(1, 'projname', projname)
            call spproj%projinfo%set(1, 'projfile', projfile)
            call spproj%projinfo%set(1, 'cwd',      direxec)
            call spproj%os_ptcl3D%delete_2Dclustering
            ! fetch frcs & classes
            path = trim(stemname(project_fname))
            tmpl = trim(basename_safe(project_fname))
            tmpl = trim(get_fbody(tmpl, METADATA_EXT, separator=.false.))
            smpd = spproj%get_smpd()
            box  = spproj%get_box()
            call spproj%get_cavgs_stk(tmp1, n, smpd_sc)
            ! call find_ldim_nptcls(string, ldim_sc, n, smpd=smpd_sc)
            ! classes, pad to original dimensions as images were scaled during alignment
            found  = .false.
            string = trim(tmpl)//'_'//trim(CAVGS_ITER_FBODY)
            call simple_list_files_regexp(path, '\.mrc$|\.mrcs$', files)
            if( allocated(files) )then
                do j = 1,size(files)
                    if( str_has_substr(files(j), string) )then
                        if( str_has_substr(files(j),'_even') .or.&
                            &str_has_substr(files(j),'_odd') )cycle
                        string = trim(files(j))
                        found = .true.
                        exit
                    endif
                enddo
                deallocate(files)
            endif
            if( .not.found )then
                THROW_WARN('Cannot find classes at: '//string//'. Skipping')
                call spproj%kill
                return
            endif
            j     = index(string,CAVGS_ITER_FBODY)
            cavgs = trim(direxec)//string(j:len_trim(string))
            call scale_imgfile(string, cavgs, smpd_sc, [box,box,1], smpd)
            call spproj%add_cavgs2os_out(cavgs, smpd, 'cavg')
            tmp1 = add2fbody(string, params_glob%ext, '_even')
            tmp2 = add2fbody(cavgs,  params_glob%ext, '_even')
            call scale_imgfile(tmp1, tmp2, smpd_sc, [box,box,1], smpd)
            tmp1 = add2fbody(string, params_glob%ext, '_odd')
            tmp2 = add2fbody(cavgs,  params_glob%ext, '_odd')
            call scale_imgfile(tmp1, tmp2, smpd_sc, [box,box,1], smpd)
            ! frcs, pad to original dimensions as images were scaled during alignment
            string = trim(path)//trim(tmpl)//'_'//trim(FRCS_FILE)
            if( .not.file_exists(string) )then
                THROW_WARN('Cannot find '//string//'. Skipping')
                call spproj%kill
                return
            endif
            frc  = direxec//trim(FRCS_FILE)
            call frcs%read(string)
            call frcs%pad(smpd, box)
            call frcs%write(frc)
            call spproj%add_frcs2os_out(frc, 'frc2D')
            ! write project
            call spproj%write(direxec//projfile)
            ! new record
            call append_procrecord(procrecs, id, direxec, projfile)
            ! tidy
            call spproj%kill
            call frcs%kill
            write(logfhandle,'(A,A)')'>>> ADDED ONE JOB: ',procrecs(get_nrecs())%folder
        enddo
    end subroutine add_projects2jobslist

    ! submits all abinitio3D jobs to the queue
    subroutine submit_jobs( cline )
        type(cmdline),     intent(in) :: cline
        integer :: i, nsubmitted
        nsubmitted = get_nrecs_submitted()
        if( nsubmitted >= params_glob%maxnruns ) return ! all resources used
        if( nsubmitted == get_nrecs() )         return ! all jobs submitted
        do i = 1,get_nrecs()
            if( procrecs(i)%submitted ) cycle
            if( get_nrecs_submitted() >= params_glob%maxnruns ) exit
            call submit_one_process(cline, procrecs(i))
        enddo
    end subroutine submit_jobs

    ! submits one abinitio3D job to the queue
    subroutine submit_one_process( cline, record )
        type(cmdline),     intent(in)    :: cline
        class(procrecord), intent(inout) :: record
        type(chash)                :: job_descr
        type(qsys_env)             :: qenv
        type(cmdline)              :: cline_abinitio3D
        character(len=XLONGSTRLEN) :: cwd, cwd_old
        ! updates command line
        cline_abinitio3D = cline
        call cline_abinitio3D%set('projfile', record%projfile)
        call cline_abinitio3D%gen_job_descr(job_descr)
        ! submit
        cwd_old = trim(cwd_glob)
        call chdir(record%folder)
        call simple_getcwd(cwd)
        cwd_glob = trim(cwd)
        call qenv%new(params_glob%nparts, exec_bin='simple_exec')
        call qenv%exec_simple_prg_in_queue_async(cline_abinitio3D, 'execscript', 'execlog')
        call chdir(cwd_old)
        cwd_glob = trim(cwd_old)
        time_last_submission = time8()
        ! update global record
        record%submitted = .true.
        ! tidy
        call cline_abinitio3D%kill
        call qenv%kill
        call job_descr%kill
        write(logfhandle,'(A,A)')'>>> SUBMITTED ONE JOB: ',cast_time_char(simple_gettime())
    end subroutine submit_one_process

    subroutine check_processes( ncomplete )
        integer, intent(inout) :: ncomplete
        integer :: i
        ncomplete = 0
        if( .not.allocated(procrecs) )return
        do i = 1,size(procrecs)
            if( procrecs(i)%complete ) cycle
            if( procrecs(i)%submitted )then
                procrecs(i)%complete = file_exists(trim(procrecs(i)%folder)//&
                    &trim(ABINITIO3D_FINISHED))
                if( procrecs(i)%complete )then
                    procrecs(i)%included = .false.
                    ncomplete = ncomplete + 1
                    write(logfhandle,'(A,A)')'>>> ONE JOB COMPLETED AT: ',trim(procrecs(i)%folder)
                endif
            endif
        enddo
    end subroutine check_processes

    subroutine analysis( master_spproj )
        use simple_projector_hlev, only: rotvol
        use simple_image,          only: image
        class(sp_project), intent(inout) :: master_spproj
        type(procrecord),    allocatable :: prevrecs(:)    ! already analyzed
        type(procrecord),    allocatable :: newrecs(:)     ! to analyze
        character(len=:),    allocatable :: projfile, projections, fsc_fname
        type(sp_project) :: spproj
        type(ori)        :: o
        type(image)      :: vol, rvol, reprojs
        real    :: R(3,3), smpd
        integer :: i, j, n_new, n_prev, box, nptcls, fsc_box ! Joe: box & fsc_box are the same?
        n_new = get_nrecs2postprocess()
        if( n_new == 0 ) return
        n_prev = get_nrecs_completed() - n_new
        if( n_prev > 0 ) prevrecs = pack(procrecs, mask=procrecs(:)%included)
        allocate(newrecs(n_new))
        ! Individual analysis
        j = 0
        do i = 1,size(procrecs)
            if( .not.procrecs(i)%complete ) cycle   ! not finished or submitted
            if(      procrecs(i)%included ) cycle   ! previously done
            j        = j + 1
            projfile = procrecs(i)%folder//procrecs(i)%projfile
            ! parameters
            call spproj%read_segment('ptcl3D', projfile)
            ! unfiltered volume
            call spproj%read_segment('out', projfile)
            ! should use postprocessed volume?
            call spproj%get_vol('vol', 1, procrecs(i)%volume, smpd, box)
            call spproj%get_fsc( 1, fsc_fname, fsc_box )
            nptcls = spproj%os_ptcl3D%get_noris(consider_state=.true.)
            call spproj%kill
            ! Individual inertia tensor alignment
            call vol%new([box,box,box], smpd)
            call vol%read(procrecs(i)%volume)
            procrecs(i)%alnvolume = trim(DIR_ALGNVOLS)//'alnvol_'//trim(procrecs(i)%id)//trim(params_glob%ext)
            if( trim(uppercase(params_glob%pgrp)).eq.'C1' )then
                ! align volume inertia tensor principal axes to reference axes
                call vol%calc_principal_axes_rotmat(params_glob%mskdiam, R)
                call o%ori_from_rotmat(R, .false.)
                rvol = rotvol(vol, o, [0.,0.,0.])
                call rvol%write(procrecs(i)%alnvolume)
                call rvol%generate_orthogonal_reprojs(reprojs)
            else
                call vol%write(procrecs(i)%alnvolume)
                call vol%generate_orthogonal_reprojs(reprojs)
            endif
            call reprojs%write_jpg(trim(DIR_ALGNVOLS)//'alnvol_reprojs_'//trim(procrecs(i)%id)//'.jpg')
            call reprojs%kill
            call rvol%kill
            call o%kill
            ! updates all-to-all correlation matrix
            if( (get_nrecs_included()==0) .and. (j==1) )then
                ! first time, nothing to do
            else
                call vol%mirror('x')
                call vol%write('mirr.mrc')
                call dock_new_volume(i)
                call del_file('mirr.mrc')
            endif
            call vol%kill
            ! will not go through this loop again
            procrecs(i)%included = .true.
            ! book-keeping
            newrecs(j) = procrecs(i)
            ! add to global project as state=nvols_glob
            nvols_glob = nvols_glob+1
            call master_spproj%add_vol2os_out(procrecs(i)%alnvolume, smpd, nvols_glob, 'vol', box, pop=nptcls)
            call master_spproj%add_fsc2os_out(fsc_fname, nvols_glob, fsc_box)
            ! TODO align classes to volumes!
        enddo
        ! write global project
        call master_spproj%write_segment_inside('out',params_glob%projfile)
        ! cleanup
        call kill_procrecords(prevrecs)
        call kill_procrecords(newrecs)
    end subroutine analysis

    ! docks volume to all the previous ones
    subroutine dock_new_volume( indnew )
        use simple_dock_vols, only: dock_vols
        integer,          intent(in)  :: indnew
        type(dock_vols)               :: dvols
        real,             allocatable :: tmp(:,:)
        character(len=:), allocatable :: newvol
        integer,          allocatable :: inds(:)
        real    :: eul(3), eul_mirr(3), shift(3), shift_mirr(3), cc, cc_mirr
        integer :: i, ii, n, ncompl, nincl, l
        if( .not.allocated(procrecs) )return
        nincl  = get_nrecs_included()
        ncompl = get_nrecs_completed()
        if( nincl < 1 )  return
        if( ncompl < 1 ) return
        l      = nincl + 1
        n      = size(procrecs)
        inds   = pack((/(i,i=1,n)/), mask=procrecs(:)%complete)
        newvol = trim(procrecs(indnew)%volume)
        write(logfhandle,'(A,A)') '>>> DOCKING VOLUME: ',newvol
        if( allocated(algnvols_corrmat) )then
            allocate(tmp(l,l),source=0.)
            tmp(:nincl,:nincl)    = algnvols_corrmat(:,:)
            deallocate(algnvols_corrmat)
            call move_alloc(tmp, algnvols_corrmat)
            algnvols_corrmat(l,l) = 1.0
        else
            allocate(algnvols_corrmat(2,2),source=0.)
            algnvols_corrmat(1,1) = 1.0
            algnvols_corrmat(2,2) = 1.0
        endif
        do i = 1,nincl
            ii = inds(i)
            if( ii == indnew ) exit
            call dvols%new(newvol, procrecs(ii)%volume, params_glob%smpd, 150., 10.0, params_glob%mskdiam)
            call dvols%srch
            call dvols%get_dock_info(eul, shift, cc)
            call dvols%kill
            call dvols%new('mirr.mrc', procrecs(ii)%volume, params_glob%smpd, 150., 10.0, params_glob%mskdiam)
            call dvols%srch
            call dvols%get_dock_info(eul_mirr, shift_mirr, cc_mirr)
            call dvols%kill
            if( cc_mirr > cc )then
                cc    = cc_mirr
                shift = shift_mirr
                eul   = eul_mirr
            endif
            algnvols_corrmat(l,i) = max(cc,cc_mirr)
            algnvols_corrmat(i,l) = algnvols_corrmat(l,i)
        end do
    end subroutine dock_new_volume

    ! Book-keeping

    !> Updates current parameters with user input
    !  Parameters can be either global ones used by master process (params_glob updated)
    !  AND/OR can be for abinitio3D, in which case cline must be updated!
    subroutine update_user_params3D( cline, updated )
        use simple_sym, only: is_valid_pointgroup
        class(cmdline), intent(inout) :: cline
        logical,        intent(out)   :: updated
        type(oris)                    :: os
        character(len=:), allocatable :: string
        integer :: nptcls
        updated = .false.
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS3D) )then
            call os%read(USER_PARAMS3D)
            if( os%isthere(1,'nptcls') )then
                nptcls = os%get_int(1,'nptcls')
                if( nptcls > 0 )then
                    params_glob%nptcls = nptcls
                    write(logfhandle,'(A,I8)')'>>> # OF PTCLS / SNAPSHOT UPDATED TO: ',params_glob%nptcls
                    updated = .true.
                else
                    THROW_WARN('Invalid NPTCLS!')
                endif
            endif
            if( os%isthere(1,'pgrp') )then
                call os%getter(1,'pgrp', string)
                if( is_valid_pointgroup(string) )then
                    params_glob%pgrp = trim(string)
                    call cline%set('pgrp', params_glob%pgrp) ! command-line updated
                    write(logfhandle,'(A,A)')'>>> PGRP UPDATED TO: ', trim(params_glob%pgrp)
                    updated = .true.
                else
                    THROW_WARN('Invalid PGRP!')
                endif
            endif
            if( os%isthere(1,'pgrp_start') )then
                call os%getter(1,'pgrp_start', string)
                if( is_valid_pointgroup(string) )then
                    params_glob%pgrp_start = trim(string)
                    call cline%set('pgrp_start', params_glob%pgrp_start) ! command-line updated
                    write(logfhandle,'(A,A)')'>>> PGRP_START UPDATED TO: ', trim(params_glob%pgrp_start)
                    updated = .true.
                else
                    THROW_WARN('Invalid PGRP_START!')
                endif
            endif
            call del_file(USER_PARAMS3D)
        endif
        call os%kill
    end subroutine update_user_params3D

    subroutine restart_cleanup
        character(len=LONGSTRLEN), allocatable :: files(:)
        character(len=STDLEN),     allocatable :: folders(:)
        integer :: i
        call del_file(USER_PARAMS3D)
        call del_file(TERM_STREAM)
        call simple_list_files_regexp('.', '\.mrc$|\.mrcs$|\.txt$|\.star$|\.eps$|\.jpeg$|\.jpg$|\.dat$|\.bin$', files)
        if( allocated(files) )then
            do i = 1,size(files)
                call del_file(files(i))
            enddo
            deallocate(files)
        endif
        folders = simple_list_dirs('.')
        if( allocated(folders) )then
            do i = 1,size(folders)
                if( str_has_substr(folders(i),trim(EXECDIR)) ) call simple_rmdir(folders(i))
            enddo
            deallocate(folders)
        endif
    end subroutine restart_cleanup

    ! Total number of jobs that have been analyzed
    integer function get_nrecs_included()
    get_nrecs_included = 0
        if( allocated(procrecs) ) get_nrecs_included = count(procrecs%included)
    end function get_nrecs_included

    ! number completed jobs
    integer function get_nrecs_completed()
        get_nrecs_completed = 0
        if( allocated(procrecs) ) get_nrecs_completed = count(procrecs%complete)
    end function get_nrecs_completed

    ! total number of jobs
    integer function get_nrecs()
        get_nrecs = 0
        if( allocated(procrecs) ) get_nrecs = size(procrecs)
    end function get_nrecs

    ! number of jobs currently submitted
    integer function get_nrecs_submitted()
        get_nrecs_submitted = 0
        if( allocated(procrecs) )then
            get_nrecs_submitted = count(procrecs%submitted) - count(procrecs%complete)
        endif
    end function get_nrecs_submitted

    ! number of jobs to post-process/analyze
    integer function get_nrecs2postprocess()
        get_nrecs2postprocess = 0
        if( allocated(procrecs) ) get_nrecs2postprocess = count(procrecs%complete .and.(.not.procrecs%included))
    end function get_nrecs2postprocess

end module simple_commander_abinitio3D_stream

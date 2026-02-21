!@descr: project commanders for particle-related things
module simple_commanders_project_ptcl
use simple_commanders_api
use simple_stream_watcher, only: stream_watcher
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_import_particles
  contains
    procedure :: execute      => exec_import_particles
end type commander_import_particles

type, extends(commander_base) :: commander_import_boxes
  contains
    procedure :: execute      => exec_import_boxes
end type commander_import_boxes

type, extends(commander_base) :: commander_zero_project_shifts
  contains
    procedure :: execute      => exec_zero_project_shifts
end type commander_zero_project_shifts

type, extends(commander_base) :: commander_prune_project_distr
  contains
    procedure :: execute      => exec_prune_project_distr
end type commander_prune_project_distr

type, extends(commander_base) :: commander_prune_project
  contains
    procedure :: execute      => exec_prune_project
end type commander_prune_project

type, extends(commander_base) :: commander_scale_project_distr
  contains
    procedure :: execute      => exec_scale_project_distr
end type commander_scale_project_distr

type, extends(commander_base) :: commander_split_stack
  contains
    procedure :: execute      => exec_split_stack
end type commander_split_stack

contains

    subroutine exec_import_particles( self, cline )
        class(commander_import_particles), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(string), allocatable      :: stkfnames(:)
        real,         allocatable      :: line(:)
        type(simple_nice_communicator) :: nice_communicator
        type(string)                   :: phaseplate, ctfstr
        type(parameters)               :: params
        type(sp_project)               :: spproj
        type(oris)                     :: os
        type(nrtxtfile)                :: paramfile
        type(ctfparams)                :: ctfvars
        integer                        :: ldim1(3), i, ndatlines, nrecs, n_ori_inputs, nstks
        logical                        :: inputted_oritab, inputted_plaintexttab, inputted_deftab, inputted_stk_den
        logical                        :: l_stktab_per_stk_parms, is_ptcl
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ctf')   ) call cline%set('ctf',   'yes')
        l_stktab_per_stk_parms = .true.
        call params%new(cline)
        ! PARAMETER INPUT MANAGEMENT
        ! parameter input flags
        inputted_oritab       = cline%defined('oritab')
        inputted_deftab       = cline%defined('deftab')
        inputted_plaintexttab = cline%defined('plaintexttab')
        inputted_stk_den      = cline%defined('stk2')
        n_ori_inputs          = count([inputted_oritab,inputted_deftab,inputted_plaintexttab])
        ! exceptions
        if( n_ori_inputs > 1 )then
            THROW_HARD('multiple parameter sources inputted, please use (oritab|deftab|plaintexttab); exec_import_particles')
        endif
        if( cline%defined('stk') .and. cline%defined('stktab') )then
            THROW_HARD('stk and stktab are both defined on command line, use either or; exec_import_particles')
        endif
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            if( trim(params%ctf) .ne. 'no' )then
                ! there needs to be associated parameters of some form
                if( n_ori_inputs < 1 )then
                    THROW_HARD('stk or stktab input requires associated parameter input when ctf .ne. no (oritab|deftab|plaintexttab)')
                endif
            endif
        else
            THROW_HARD('either stk or stktab needed on command line; exec_import_particles')
        endif
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! set is particle flag for correct field parsing
        is_ptcl = .false.
        if( cline%defined('stk') ) is_ptcl = .true.
        ! oris input
        if( inputted_oritab )then
            ndatlines = binread_nlines(params%oritab, params%spproj_iseg)
            call os%new(ndatlines, is_ptcl=is_ptcl )
            call binread_oritab(params%oritab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_deftab )then
            ndatlines = binread_nlines(params%deftab, params%spproj_iseg)
            call os%new(ndatlines, is_ptcl=is_ptcl )
            call binread_ctfparams_state_eo(params%deftab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_plaintexttab )then
            call paramfile%new(params%plaintexttab, 1)
            ndatlines = paramfile%get_ndatalines()
            nrecs     = paramfile%get_nrecs_per_line()
            if( nrecs < 1 .or. nrecs > 4 .or. nrecs == 2 )then
                THROW_HARD('unsupported nr of rec:s in plaintexttab; exec_import_particles')
            endif
            call os%new(ndatlines, is_ptcl=is_ptcl )
            allocate( line(nrecs) )
            do i=1,ndatlines
                call paramfile%readNextDataLine(line)
                select case(params%dfunit)
                    case( 'A' )
                        line(1) = line(1)/1.0e4
                        if( nrecs > 1 )  line(2) = line(2)/1.0e4
                    case( 'microns' )
                        ! nothing to do
                    case DEFAULT
                        THROW_HARD('unsupported dfunit; exec_import_particles')
                end select
                select case(params%angastunit)
                    case( 'radians' )
                        if( nrecs == 3 ) line(3) = rad2deg(line(3))
                    case( 'degrees' )
                        ! nothing to do
                    case DEFAULT
                        THROW_HARD('unsupported angastunit; exec_import_particles')
                end select
                select case(params%phshiftunit)
                    case( 'radians' )
                        ! nothing to do
                    case( 'degrees' )
                        if( nrecs == 4 ) line(4) = deg2rad(line(4))
                    case DEFAULT
                        THROW_HARD('unsupported phshiftunit; exec_import_particles')
                end select
                call os%set_dfx(i, line(1))
                if( nrecs > 1 )then
                    call os%set_dfy(i,       line(2))
                    call os%set(i, 'angast', line(3))
                endif
                if( nrecs > 3 )then
                    call os%set(i, 'phshift', line(4))
                endif
            end do
        endif
        if( cline%defined('stktab') )then
            ! importing from stktab
            call read_filetable(params%stktab, stkfnames)
            nstks = size(stkfnames)
            if( params%mkdir.eq.'yes' )then
                do i=1,nstks
                    if(stkfnames(i)%to_char([1,1]).ne.'/') stkfnames(i) = PATH_PARENT//stkfnames(i)%to_char()
                    if( .not. file_exists(stkfnames(i)) ) THROW_HARD('modified filetable entry '//stkfnames(i)%to_char()//' does not exist')
                enddo
            endif
            l_stktab_per_stk_parms = (os%get_noris() == nstks)
            if( (n_ori_inputs == 1) .and. l_stktab_per_stk_parms )then
                ! sampling distance
                call os%set_all2single('smpd', params%smpd)
                ! acceleration voltage
                if( cline%defined('kv') )then
                    call os%set_all2single('kv', params%kv)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'kv') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks acceleration volatage (kv)'
                            THROW_HARD('provide kv on command line or update input document; exec_import_particles')
                        endif
                    end do
                endif
                ! spherical aberration
                if( cline%defined('cs') )then
                    call os%set_all2single('cs', params%cs)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'cs') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks spherical aberration constant (cs)'
                            THROW_HARD('provide cs on command line or update input document; exec_import_particles')
                        endif
                    end do
                endif
                ! fraction of amplitude contrast
                if( cline%defined('fraca') )then
                    call os%set_all2single('fraca', params%fraca)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'fraca') )then
                            write(logfhandle,*) 'os entry: ', i, ' lacks fraction of amplitude contrast (fraca)'
                            THROW_HARD('provide fraca on command line or update input document; exec_import_particles')
                        endif
                    end do
                endif
                ! phase-plate
                if( cline%defined('phaseplate') )then
                    call os%set_all2single('phaseplate', trim(params%phaseplate))
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'phaseplate') )then
                            call os%set(i, 'phaseplate', 'no')
                        endif
                    end do
                endif
                call os%getter(1, 'phaseplate', phaseplate)
                if( phaseplate .eq. 'yes' )then
                    if( .not. os%isthere(1,'phshift') )then
                        THROW_HARD('phaseplate .eq. yes requires phshift input, currently lacking; exec_import_particles')
                    endif
                endif
                ! ctf flag
                if( cline%defined('ctf') )then
                    call os%set_all2single('ctf', trim(params%ctf))
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'ctf') )then
                            call os%set(i, 'ctf', 'yes')
                        endif
                    end do
                endif
                call os%getter(1, 'ctf', ctfstr)
                if( ctfstr .ne. 'no' )then
                    if( .not. os%isthere(1,'dfx') )then
                        THROW_HARD('ctf .ne. no requires dfx input, currently lacking; exec_import_particles')
                    endif
                endif
            endif
        endif
        if( cline%defined('stk') .or. (cline%defined('stktab').and..not.l_stktab_per_stk_parms) )then
            ctfvars%smpd = params%smpd
            select case(trim(params%ctf))
                case('yes')
                    ctfvars%ctfflag = CTFFLAG_YES
                case('no')
                    ctfvars%ctfflag = CTFFLAG_NO
                case('flip')
                    ctfvars%ctfflag = CTFFLAG_FLIP
                case DEFAULT
                    write(logfhandle,*)
                    THROW_HARD('unsupported ctf flag: '//trim(params%ctf)//'; exec_import_particles')
            end select
            if( ctfvars%ctfflag .ne. CTFFLAG_NO )then
                if( .not. cline%defined('kv')    ) THROW_HARD('kv (acceleration voltage in kV{300}) input required when importing movies; exec_import_particles')
                if( .not. cline%defined('cs')    ) THROW_HARD('cs (spherical aberration constant in mm{2.7}) input required when importing movies; exec_import_particles')
                if( .not. cline%defined('fraca') ) THROW_HARD('fraca (fraction of amplitude contrast{0.1}) input required when importing movies; exec_import_particles')
                if( cline%defined('phaseplate') )then
                    phaseplate = cline%get_carg('phaseplate')
                else
                    phaseplate ='no'
                endif
                ctfvars%kv           = params%kv
                ctfvars%cs           = params%cs
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = phaseplate .eq. 'yes'
            endif
        endif

        ! PROJECT FILE MANAGEMENT
        call spproj%read(params%projfile)

        ! UPDATE FIELDS
        ! add stack if present
        if( cline%defined('stk') )then
            if( n_ori_inputs == 0 .and. trim(params%ctf) .eq. 'no' )then
                ! get number of particles from stack
                call find_ldim_nptcls(params%stk, ldim1, params%nptcls)
                call os%new(params%nptcls, is_ptcl=is_ptcl )
            endif
            ! state = 1 by default
            call os%set_all2single('state', 1.0)
            call os%set_all2single('w',     1.0)
            call spproj%add_single_stk(params%stk, ctfvars, os)
        endif
        ! add list of stacks (stktab) if present
        if( cline%defined('stktab') )then
            ! state = 1 by default
            call os%set_all2single('state', 1.0)
            call os%set_all2single('w',     1.0)
            if( l_stktab_per_stk_parms )then
                ! per stack parameters
                call spproj%add_stktab(stkfnames, os)
            else
                ! per particle parameters
                call spproj%add_stktab(stkfnames, ctfvars, os)
            endif
        endif
        ! WRITE PROJECT FILE
        call spproj%write ! full write since this is guaranteed to be the first import
        call nice_communicator%terminate()
        call simple_end('**** IMPORT_PARTICLES NORMAL STOP ****')
    end subroutine exec_import_particles

    subroutine exec_import_boxes( self, cline )
        class(commander_import_boxes), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        integer                        :: nos_mic, nboxf, i
        type(string), allocatable      :: boxfnames(:)
        type(string)                   :: boxfname, intg, cwd, boxname
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call simple_getcwd(cwd)
        ! project file management
        if( .not. file_exists(params%projfile) )then
            THROW_HARD('project file: '//params%projfile%to_char()//' does not exist! exec_import_boxes')
        endif
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        if(params%reset_boxfiles .eq. 'yes') then
            call spproj%read(params%projfile)
            do i=1,spproj%os_mic%get_noris()
                intg = spproj%os_mic%get_str(i, "intg")
                boxname = swap_suffix(intg, ".box", ".mrc")
                boxfname = cwd //'/'//basename(boxname)
                call spproj%os_mic%set(i, 'boxfile', boxfname)
                call spproj%os_mic%set(i, 'nptcls',  0)
            end do
            call spproj%os_stk%new(0,    is_ptcl=.false.)
            call spproj%os_ptcl2D%new(0, is_ptcl=.true.)
            call spproj%os_cls2D%new(0,  is_ptcl=.false.)
            call spproj%os_cls3D%new(0,  is_ptcl=.false.)
            call spproj%os_ptcl3D%new(0, is_ptcl=.true.)
            call spproj%write(params%projfile)
            call boxname%kill
        else
            call spproj%read(params%projfile)
            ! get boxfiles into os_mic
            call read_filetable(params%boxtab, boxfnames)
            nboxf   = size(boxfnames)
            nos_mic = spproj%os_mic%get_noris()
            if( nboxf /= nos_mic )then
                write(logfhandle,*) '# boxfiles       : ', nboxf
                write(logfhandle,*) '# os_mic entries : ', nos_mic
                THROW_HARD('# boxfiles .ne. # os_mic entries; exec_import_boxes')
            endif
            do i=1,nos_mic
                if( params%mkdir.eq.'yes' )then
                    if(boxfnames(i)%to_char([1,1]).ne.'/') boxfnames(i) = PATH_PARENT//boxfnames(i)%to_char()
                endif
                boxfname = simple_abspath(boxfnames(i))
                call spproj%set_boxfile(i, boxfname)
            end do
            ! write project file
            call spproj%write_segment_inside('mic') ! all that's needed here
        end if
        call nice_communicator%terminate()
        call simple_end('**** IMPORT_BOXES NORMAL STOP ****')
    end subroutine exec_import_boxes

    subroutine exec_zero_project_shifts( self, cline )
        class(commander_zero_project_shifts), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        call cline%set('mkdir', 'yes')
        call params%new(cline, silent=.true.)
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        call spproj%read(params%projfile)
        call spproj%os_ptcl2D%zero_shifts
        call spproj%os_ptcl3D%zero_shifts
        call spproj%write(params%projfile)
        call spproj%kill
        call nice_communicator%terminate()
    end subroutine exec_zero_project_shifts

    subroutine exec_prune_project_distr( self, cline )
        class(commander_prune_project_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(cmdline)                 :: cline_distr
        type(sp_project)              :: spproj
        type(sp_project), allocatable :: spproj_part(:)
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(chash),      allocatable :: part_params(:)
        integer,          allocatable :: parts(:,:), pinds(:)
        real,             allocatable :: states(:)
        type(string)                  :: fname
        logical,          allocatable :: part_mask(:)
        integer :: imic,nmics,cnt,istk,nstks,ipart,nptcls,nparts,iptcl
        integer :: nstks_orig,nptcls_orig,nmics_orig,i
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! sanity checks
        call spproj%read(params%projfile)
        nstks = spproj%get_n_insegment('stk')
        if( nstks == 0 ) THROW_HARD('No stack to process!')
        nptcls = spproj%get_n_insegment('ptcl2D')
        if( nstks == 0 ) THROW_HARD('No particles to process!')
        ! identify the particle indices with state .ne. 0 unless state is given on command line
        states = spproj%os_ptcl2D%get_all('state')
        allocate(pinds(nptcls), source=(/(i,i=1,nptcls)/))
        if( cline%defined('state') )then
            pinds = pack(pinds, mask=abs(states - real(params%state)) < 0.1)
        else
            pinds = pack(pinds, mask=states > 0.5)
        endif
        call arr2txtfile(pinds, string('selected_indices.txt'))
        deallocate(states, pinds)
        nmics       = spproj%get_nintgs()
        nstks_orig  = nstks
        nptcls_orig = nptcls
        nmics_orig  = nmics
        ! DISTRIBUTED EXECUTION
        ! setup the environment for distributed execution
        cline_distr = cline
        call cline_distr%set('prg',     'prune_project')
        call cline_distr%set('nthr',    1)
        call cline_distr%set('oritype', 'stk')
        nparts = min(params%nparts, nstks)
        allocate(spproj_part(nparts),part_mask(nparts))
        parts     = split_nobjs_even(nstks, nparts)
        part_mask = .true.
        allocate(part_params(nparts))
        do ipart=1,nparts
            call part_params(ipart)%new(2)
            call part_params(ipart)%set('fromp',int2str(parts(ipart,1)))
            call part_params(ipart)%set('top',  int2str(parts(ipart,2)))
        end do
        call qenv%new(params, nparts)
        ! prepare job description
        call cline_distr%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR)
        ! ASSEMBLY
        do ipart = 1,nparts
            fname = ALGN_FBODY//int2str(ipart)//METADATA_EXT
            part_mask(ipart) = file_exists(fname)
        enddo
        ! copy updated micrographs
        nstks = 0
        if( nmics > 0 )then
            cnt = 0
            do ipart = 1,nparts
                if( .not.part_mask(ipart) ) cycle
                fname = ALGN_FBODY//int2str(ipart)//METADATA_EXT
                call spproj_part(ipart)%read_segment('mic',fname)
                cnt = cnt + spproj_part(ipart)%os_mic%get_noris()
            enddo
            nmics = cnt
            if( nmics > 0 )then
                call spproj%os_mic%new(nmics, is_ptcl=.false.)
                cnt = 0
                do ipart = 1,nparts
                    if( spproj_part(ipart)%os_mic%get_noris() == 0 ) cycle
                    do imic = 1,spproj_part(ipart)%os_mic%get_noris()
                        cnt = cnt + 1
                        call spproj%os_mic%transfer_ori(cnt, spproj_part(ipart)%os_mic, imic)
                        if( spproj%os_mic%get_int(cnt,'nptcls') > 0 ) nstks = nstks + 1
                    enddo
                    call spproj_part(ipart)%kill
                enddo
                if( nstks /= nmics ) THROW_HARD('Inconsistent number of stacks and micrographs!')
            endif
            write(logfhandle,'(A,I8,A,I8)')'>>> # OF MICROGRAPHS:',nmics_orig,' -> ', nmics
        endif
        ! copy updated stacks
        cnt = 0
        do ipart = 1,nparts
            if( .not.part_mask(ipart) ) cycle
            fname = ALGN_FBODY//int2str(ipart)//METADATA_EXT
            call spproj_part(ipart)%read_segment('stk',fname)
            cnt = cnt + spproj_part(ipart)%os_stk%get_noris()
        enddo
        nstks = cnt
        write(logfhandle,'(A,I8,A,I8)')'>>> # OF STACKS     :',nstks_orig,' -> ', nstks
        if( nstks > 0 )then
            call spproj%os_stk%new(nstks, is_ptcl=.false.)
            cnt = 0
            do ipart = 1,nparts
                if( spproj_part(ipart)%os_stk%get_noris() == 0 ) cycle
                do istk = 1,spproj_part(ipart)%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj%os_stk%transfer_ori(cnt, spproj_part(ipart)%os_stk, istk)
                enddo
                call spproj_part(ipart)%kill
            enddo
            ! copy updated particles 2D segment
            cnt = 0
            do ipart = 1,nparts
                if( .not.part_mask(ipart) ) cycle
                fname = ALGN_FBODY//int2str(ipart)//METADATA_EXT
                call spproj_part(ipart)%read_segment('ptcl2D',fname)
                cnt = cnt + spproj_part(ipart)%os_ptcl2D%get_noris()
            enddo
            nptcls = cnt
            call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
            cnt = 0
            do ipart = 1,nparts
                if( spproj_part(ipart)%os_ptcl2D%get_noris() == 0 ) cycle
                do iptcl = 1,spproj_part(ipart)%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    call spproj%os_ptcl2D%transfer_ori(cnt, spproj_part(ipart)%os_ptcl2D, iptcl)
                enddo
                call spproj_part(ipart)%kill
            enddo
            ! copy updated particles 3D segment
            call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
            cnt = 0
            do ipart = 1,nparts
                if( .not.part_mask(ipart) ) cycle
                fname = ALGN_FBODY//int2str(ipart)//METADATA_EXT
                call spproj_part(ipart)%read_segment('ptcl3D',fname)
                if( spproj_part(ipart)%os_ptcl3D%get_noris() == 0 ) cycle
                do iptcl = 1,spproj_part(ipart)%os_ptcl3D%get_noris()
                    cnt = cnt + 1
                    call spproj%os_ptcl3D%transfer_ori(cnt, spproj_part(ipart)%os_ptcl3D, iptcl)
                enddo
                call spproj_part(ipart)%kill
            enddo
        endif
        write(logfhandle,'(A,I8,A,I8)')'>>> # OF PARTICLES  :',nptcls_orig,' -> ', nptcls
        ! final write
        call spproj%write(params%projfile)
        ! clean up
        call spproj%kill
        call qsys_cleanup
        do ipart = 1,nparts
            call part_params(ipart)%kill
            call spproj_part(ipart)%kill
            call del_file(ALGN_FBODY//int2str(ipart)//METADATA_EXT)
        enddo
        deallocate(spproj_part,part_params)
        ! end gracefully
        call simple_end('**** SIMPLE_PRUNE_PROJECT_DISTR NORMAL STOP ****')
    end subroutine exec_prune_project_distr

    subroutine exec_prune_project( self, cline )
        class(commander_prune_project), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: img
        type(sp_project)     :: spproj, spproj_out
        type(ori)            :: o_stk
        type(string)         :: newstkname,stkname,ext
        logical, allocatable :: stks_mask(:), ptcls_mask(:)
        integer, allocatable :: stkinds(:), stk2mic_inds(:), mic2stk_inds(:)
        type(string)         :: absstkname, stkdir
        real                 :: smpd
        integer              :: iptcl, istk, stk_cnt, nptcls_tot, ptcl_cnt
        integer              :: box, nstks, nstks_tot, fromp, top, fromp_glob, top_glob, nmics_tot
        integer              :: nstks_part, nptcls_part, stkind, nstks_prev, ptcl_glob
        ! init
        call params%new(cline)
        if( params%dir .eq. '' )then
            stkdir = PATH_HERE
        else
            stkdir = params%dir
        endif
        ! particles
        call spproj%read_segment('ptcl2D', params%projfile)
        nptcls_tot = spproj%os_ptcl2D%get_noris()
        allocate(ptcls_mask(nptcls_tot), stkinds(nptcls_tot))
        nptcls_part = 0
        !$omp parallel do proc_bind(close) default(shared) private(iptcl) reduction(+:nptcls_part)
        do iptcl=1,nptcls_tot
            ptcls_mask(iptcl) = spproj%os_ptcl2D%get_state(iptcl) > 0
            if( ptcls_mask(iptcl) )then
                stkinds(iptcl) = spproj%os_ptcl2D%get_int(iptcl,'stkind')
                if( stkinds(iptcl) >= params%fromp .and. stkinds(iptcl) <= params%top )then
                    nptcls_part = nptcls_part+1
                endif
            else
                stkinds(iptcl) = 0
            endif
        enddo
        !$omp end parallel do
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj_out%os_ptcl2D%new(nptcls_part, is_ptcl=.true.)
        call spproj_out%os_ptcl3D%new(nptcls_part, is_ptcl=.true.)
        ! stacks
        call spproj%read_segment('stk', params%projfile)
        nstks_tot = spproj%get_nstks()
        if( nstks_tot == 0 ) THROW_HARD('No images to operate on!')
        allocate(stks_mask(nstks_tot))
        do istk=1,nstks_tot
            stks_mask(istk) = spproj%os_stk%get_state(istk) > 0
            if( count(stkinds==istk) == 0 ) stks_mask(istk) = .false.
        enddo
        nstks = count(stks_mask)
        nstks_part = count(stks_mask(params%fromp:params%top))
        if( nstks_part == 0 )then
            call qsys_job_finished(string('simple_commanders_project_ptcl :: exec_prune_project'))
            return
        endif
        call spproj_out%os_stk%new(nstks_part, is_ptcl=.false.)
        ! micrographs
        call spproj%read_segment('mic', params%projfile)
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot > 0 )then
            call spproj%get_mic2stk_inds(mic2stk_inds, stk2mic_inds)
            call spproj_out%os_mic%new(nstks_part, is_ptcl=.false.)
        endif
        ! new stacks
        box  = spproj%get_box()
        smpd = spproj%get_smpd()
        write(logfhandle,'(A)')'>>> GENERATING STACK(S)'
        call img%new([box,box,1],smpd)
        call simple_mkdir(stkdir)
        nstks_prev = count(stks_mask(:params%fromp-1))
        stkind     = nstks_prev
        stk_cnt    = 0
        if( params%fromp == 1 )then
            top_glob   = 0
        else
            top       = spproj%os_stk%get_top(params%fromp-1)
            top_glob  = count(ptcls_mask(1:top))
        endif
        ptcl_glob  = 0
        do istk=params%fromp,params%top
            if( .not.stks_mask(istk) ) cycle
            stk_cnt = stk_cnt + 1
            stkind  = stkind  + 1
            call spproj%os_stk%get_ori(istk, o_stk)
            call o_stk%getter('stk',stkname)
            ext        = fname2ext(stkname)
            newstkname = stkdir//get_fbody(basename(stkname),ext)//STK_EXT
            fromp = o_stk%get_fromp()
            top   = o_stk%get_top()
            fromp_glob = top_glob+1
            ptcl_cnt   = 0
            do iptcl=fromp,top
                if( .not.ptcls_mask(iptcl) )cycle
                ptcl_glob = ptcl_glob + 1
                top_glob  = top_glob+1
                ptcl_cnt  = ptcl_cnt+1
                ! copy image
                if(spproj%os_ptcl2D%isthere(iptcl, 'indstk') .and. spproj%os_ptcl2D%get(iptcl, 'indstk') > 0.0) then
                    ! write(logfhandle, *) "STK ", spproj%os_ptcl2D%get(iptcl,'indstk')
                    call img%read(stkname, spproj%os_ptcl2D%get_int(iptcl,'indstk'))
                else
                    ! write(logfhandle, *) "STK2 " // int2str(spproj%os_ptcl2D%get_int(iptcl,'indstk'))
                    call img%read(stkname, iptcl-fromp+1)
                endif
                call img%write(newstkname, ptcl_cnt)
                ! update orientations
                call spproj_out%os_ptcl2D%transfer_ori(ptcl_glob, spproj%os_ptcl2D, iptcl)
                call spproj_out%os_ptcl3D%transfer_ori(ptcl_glob, spproj%os_ptcl3D, iptcl)
                call spproj_out%os_ptcl2D%set(ptcl_glob,'stkind',stkind)
                call spproj_out%os_ptcl3D%set(ptcl_glob,'stkind',stkind)
                call spproj_out%os_ptcl2D%set(ptcl_glob,'indstk',ptcl_cnt)
                call spproj_out%os_ptcl3D%set(ptcl_glob,'indstk',ptcl_cnt)
            enddo
            ! update stack
            absstkname = simple_abspath(newstkname)
            call o_stk%set('stk',   absstkname)
            call o_stk%set('fromp', fromp_glob)
            call o_stk%set('top',   top_glob)
            call o_stk%set('nptcls',ptcl_cnt)
            call spproj_out%os_stk%set_ori(stk_cnt, o_stk)
            ! update micrograph
            if( nmics_tot > 0 ) then
                call spproj_out%os_mic%transfer_ori(stk_cnt, spproj%os_mic, stk2mic_inds(istk))
                call spproj_out%os_mic%set(stk_cnt,'nptcls',ptcl_cnt)
            endif
        enddo
        spproj_out%projinfo = spproj%projinfo
        spproj_out%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris() > 0 ) spproj_out%jobproc = spproj%jobproc
        call spproj%kill
        call spproj_out%write(string(ALGN_FBODY)//int2str(params%part)//METADATA_EXT)
        ! cleanup
        call spproj_out%kill
        call img%kill
        call o_stk%kill
        ! end gracefully
        call qsys_job_finished(string('simple_commanders_project_ptcl :: exec_prune_project'))
    end subroutine exec_prune_project

    subroutine exec_scale_project_distr( self, cline )
        class(commander_scale_project_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(chash),      allocatable :: part_params(:)
        integer,          allocatable :: parts(:,:)
        integer,          parameter   :: MAX_NCUNITS = 64
        type(string)     :: dir_target, projfile_sc
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        type(cmdline)    :: cline_scale
        type(parameters) :: params
        type(builder)    :: build
        real             :: smpd, smpd_target
        integer          :: ipart, nparts, nstks, box, newbox, nthr_orig, nparts_orig, ncunits_orig
        logical          :: gen_sc_project
        ! mkdir=yes: a new *_sc project + stacks are generated
        ! mkdir=no : only stacks are scaled
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        gen_sc_project = cline%get_carg('mkdir').eq.'yes'
        ! make parameters and project
        call params%new(cline)
        params%nptcls = 1 ! to avoid excessive memory allocation
        call build%build_spproj(params, cline)
        call build%spproj%read_segment('stk',params%projfile)
        nstks = build%spproj%os_stk%get_noris()
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! copy command line
        cline_scale  = cline
        ! save overridden parameters
        nparts_orig  = params%nparts
        ncunits_orig = params%ncunits
        nthr_orig    = params%nthr
        if( cline%defined('nparts') )then
            nparts = min(params%nparts * params%nthr, nstks)
        else
            nparts = params%nthr
        endif
        call cline_scale%set('nparts', nparts)
        smpd = build%spproj%get_smpd()
        box  = build%spproj%get_box()
        if( gen_sc_project )then
            ! make new project & scales
            smpd_target = max(smpd, smpd * real(box)/real(params%newbox))
            dir_target = filepath(PATH_PARENT,'stack_parts_sc')
            if( cline%defined('dir_target') ) dir_target = cline%get_carg('dir_target')
            call simple_mkdir(dir_target)
            call build%spproj%scale_projfile(smpd_target, projfile_sc, cline, cline_scale, dir=dir_target)
            newbox = cline_scale%get_iarg('newbox')
            if( newbox == box )then
                write(logfhandle,*)'Inconsistent input dimensions: from ',box,' to ',newbox
                THROW_HARD('inconsistent input dimensions; exec_scale_project_distr')
            endif
        else
            newbox = params%newbox
        endif
        ! needs to be re-set
        call cline_scale%set('smpd', smpd)
        call cline_scale%set('box',  box)
        ! setup the environment for distributed execution
        params%nparts       = nparts
        params_glob%nparts  = nparts
        params%ncunits      = min(MAX_NCUNITS, nparts)
        params_glob%ncunits = min(MAX_NCUNITS, nparts)
        params%nthr         = 1
        params_glob%nthr    = 1
        call qenv%new(params, nparts)
        ! prepares stack-based parts
        parts = split_nobjs_even(nstks, nparts)
        allocate(part_params(nparts))
        do ipart=1,nparts
            call part_params(ipart)%new(2)
            call part_params(ipart)%set('fromp',int2str(parts(ipart,1)))
            call part_params(ipart)%set('top',  int2str(parts(ipart,2)))
        end do
        ! prepare job description
        call cline_scale%gen_job_descr(job_descr)
        call job_descr%set('prg',      'scale')
        call job_descr%set('newbox',   int2str(newbox))
        call job_descr%set('autoscale','no')
        call job_descr%set('nthr',     int2str(1))
        call job_descr%set('nparts',   int2str(nparts))
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, part_params=part_params, array=L_USE_SLURM_ARR, extra_params=params)
        ! delete copy in working directory
        if( gen_sc_project ) call del_file(params%projfile)
        ! clean
        call qsys_cleanup
        ! end gracefully
        params_glob%nparts  = nparts_orig
        params_glob%ncunits = ncunits_orig
        params_glob%nthr    = nthr_orig
        call build%spproj%kill
        call simple_end('**** SIMPLE_SCALE_PROJECT_DISTR NORMAL STOP ****')
    end subroutine exec_scale_project_distr

    subroutine exec_split_stack( self, cline )
        class(commander_split_stack), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call cline%set("oritype", 'stk')
        call params%new(cline, silent=.true.)
        call spproj%read(params%projfile)
        if( cline%defined('dir') )then
            call spproj%split_stk(params%nparts, params%dir)
        else
            call spproj%split_stk(params%nparts)
        endif
        call spproj%write(params%projfile)
        call spproj%kill
        call simple_end('**** SPLIT_STACK NORMAL STOP ****')
    end subroutine exec_split_stack

end module simple_commanders_project_ptcl

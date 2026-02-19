!@descr: for producing class averages
module simple_commanders_mkcavgs
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
use simple_new_classaverager
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_make_cavgs_distr
  contains
    procedure :: execute      => exec_make_cavgs_distr
end type commander_make_cavgs_distr

type, extends(commander_base) :: commander_make_cavgs
  contains
    procedure :: execute      => exec_make_cavgs
end type commander_make_cavgs

type, extends(commander_base) :: commander_cavgassemble
  contains
    procedure :: execute      => exec_cavgassemble
end type commander_cavgassemble

type, extends(commander_base) :: commander_write_classes
  contains
    procedure :: execute      => exec_write_classes
end type commander_write_classes

contains

    subroutine exec_make_cavgs_distr( self, cline )
        class(commander_make_cavgs_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)             :: params
        type(builder)                :: build
        type(cmdline)                :: cline_cavgassemble
        type(qsys_env)               :: qenv
        type(chash)                  :: job_descr
        type(commander_make_cavgs)   :: xmk_cavgs_shmem
        type(commander_cavgassemble) :: xcavgassemble
        integer :: ncls_here, nthr_here
        logical :: l_shmem
        call cline%set('wiener', 'full')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('ml_reg')  ) call cline%set('ml_reg',      'no')
        if( (cline%defined('ncls')).and. cline%defined('nspace') )then
            THROW_HARD('NCLS and NSPACE cannot be both defined!')
        endif
        if( cline%defined('nspace') )then
            if( cline%get_carg('oritype').eq.'ptcl2D' )then
                THROW_HARD('NSPACE & PTCL2D are incompatible!')
            endif
            call cline%set('oritype', 'ptcl3D')
        else
            call cline%set('oritype', 'ptcl2D')
        endif
        if( cline%defined('nparts') )then
            l_shmem = cline%get_iarg('nparts') == 1
        else
            l_shmem = .true.
        endif
        ! deal with # threads for the master process
        if( .not.l_shmem ) call set_master_num_threads(nthr_here, string('CLUSTER2D'))
        ! parse parameters & project
        call build%init_params_and_build_spproj(cline, params)
        if( cline%defined('nspace') )then
            ! handled in exec_make_cavgs
        else
            ncls_here = build%spproj_field%get_n('class')
            if( .not. cline%defined('ncls') )then
                call cline%set('ncls', ncls_here)
                params%ncls = ncls_here
            endif
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( l_shmem  )then
            call xmk_cavgs_shmem%execute_safe(cline)
            if(trim(params_glob%async).eq.'yes') call simple_touch(MAKECAVGS_FINISHED)
            return
        endif
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! prepare command lines from prototype master
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg',  'cavgassemble')
        call cline_cavgassemble%set('nthr',  nthr_here)
        if( trim(params%oritype).eq.'ptcl3D' )then
            call cline_cavgassemble%set('ncls', params%nspace)
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
        ! assemble class averages
        call xcavgassemble%execute_safe(cline_cavgassemble)
        ! end
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
        if(trim(params_glob%async).eq.'yes') call simple_touch(MAKECAVGS_FINISHED)
    end subroutine exec_make_cavgs_distr

    subroutine exec_make_cavgs( self, cline )
        class(commander_make_cavgs), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: ncls_here, icls
        logical :: l_shmem
        if( (cline%defined('ncls')).and. cline%defined('nspace') )then
            THROW_HARD('NCLS and NSPACE cannot be both defined!')
        endif
        if( cline%defined('nspace') )then
            if( cline%get_carg('oritype').eq.'ptcl2D' )then
                THROW_HARD('NSPACE & PTCL2D are incompatible!')
            endif
            call cline%set('oritype', 'ptcl3D')
            call cline%set('ncls',    cline%get_iarg('nspace'))
        else
            call cline%set('oritype', 'ptcl2D')
        endif
        ! set shared-memory flag
        l_shmem = set_shmem_flag( cline )
        if( l_shmem .and. .not. cline%defined('refs') ) THROW_HARD('need input refs (filename) for shared-memory execution')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> GENERATING CLUSTER CENTERS'
        ! deal with the orientations
        if( trim(params%oritype).eq.'ptcl3D' )then
            ! 3D class averages
            call build%eulspace%new(params%nspace, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(build%eulspace)
            call build%spproj%os_ptcl3D%set_projs(build%eulspace)
            call build%spproj%os_ptcl3D%proj2class
        else
            ! 2D
            ncls_here = build%spproj_field%get_n('class')
            if( .not. cline%defined('ncls') ) params%ncls = build%spproj_field%get_n('class')
            if( params%l_remap_cls )then
                call build%spproj_field%remap_cls()
                if( cline%defined('ncls') )then
                    if( params%ncls < ncls_here ) THROW_HARD('inputted ncls < ncls_in_oritab not allowed!')
                    if( params%ncls > ncls_here )then
                        call build%spproj_field%expand_classes(params%ncls)
                    endif
                endif
            else if( params%tseries .eq. 'yes' )then
                if( .not. cline%defined('ncls') )then
                    THROW_HARD('# class averages (ncls) need to be part of command line when tseries=yes')
                endif
                call build%spproj_field%ini_tseries(params%ncls, 'class')
                call build%spproj_field%partition_eo
            else if( params%proj_is_class.eq.'yes' )then
                call build%spproj_field%proj2class
            endif
        endif
        ! setup weights
        call build%spproj_field%calc_hard_weights2D(params%frac, params%ncls)
        ! even/odd partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! write
        if( l_shmem )then
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        else
            if( params%part .eq. 1 ) call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! choice of algorithm
        if( l_shmem )then
            if( L_NEW_CAVGER )then
                call cavger_new_new(build)
                call cavger_new_transf_oridat( build%spproj )
                call cavger_new_read_euclid_sigma2
                call cavger_new_assemble_sums( .false. )
                call cavger_new_restore_cavgs( params%frcs )
                call cavger_new_gen2Dclassdoc( build%spproj )
                call cavger_new_write_all(params%refs, params%refs_even, params%refs_odd)
                call cavger_new_kill
            else
                call cavger_new(build)
                call cavger_transf_oridat(build%spproj)
                call cavger_read_euclid_sigma2
                call cavger_assemble_sums( .false. )
                call cavger_merge_eos_and_norm
                call cavger_calc_and_write_frcs_and_eoavg(params%frcs, params%which_iter)
                call cavger_gen2Dclassdoc(build%spproj)
                call cavger_write(params%refs,      'merged')
                call cavger_write(params%refs_even, 'even'  )
                call cavger_write(params%refs_odd,  'odd'   )
                call cavger_kill
            endif
        else
            ! distributed: write partial sums only
            if( L_NEW_CAVGER )then
                call cavger_new_new(build)
                call cavger_new_transf_oridat(build%spproj)
                call cavger_new_read_euclid_sigma2
                call cavger_new_assemble_sums( .false. )
                call cavger_new_readwrite_partial_sums('write')
                call cavger_new_kill
            else
                call cavger_new(build)
                call cavger_transf_oridat(build%spproj)
                call cavger_read_euclid_sigma2
                call cavger_assemble_sums( .false. )
                call cavger_readwrite_partial_sums('write')
                call cavger_kill
            endif
            call qsys_job_finished(string('simple_commanders_cluster2D :: exec_make_cavgs'))
        endif
        ! book-keeping
        if( l_shmem )then
            select case(trim(params%oritype))
                case('ptcl2D')
                    call build%spproj%write_segment_inside('cls2D', params%projfile)
                    call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
                    call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg')
                    call build%spproj%write_segment_inside('out', params%projfile)
                case('ptcl3D')
                    do icls = 1,params%nspace
                        call build%spproj%os_cls3D%set_euler(icls, build%eulspace%get_euler(icls))
                    enddo
                    if( cline%defined('outfile') )then
                        call build%spproj%os_cls3D%write(params%outfile)
                    else
                        call build%spproj%os_cls3D%write(string('cls3D_oris.txt'))
                    endif
                    call build%spproj%write_segment_inside('cls3D', params%projfile)
                    call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc3D')
                    call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg3D')
                    call build%spproj%write_segment_inside('out', params%projfile)
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE: '//trim(params%oritype))
            end select
        endif
        ! end gracefully
        call build%kill_strategy2D_tbox
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs

    subroutine exec_cavgassemble( self, cline )
        use simple_classaverager
        class(commander_cavgassemble), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(starproject)  :: starproj
        type(polarft_calc) :: pftc
        type(string)       :: fname
        real, allocatable  :: states(:)
        real               :: clw
        integer            :: iterstr_start, iterstr_end, iter, io_stat, icls
        integer            :: pftsz, kfromto(2), ncls
        logical            :: l_stream
        if( .not.cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = cline%get_carg('stream')=='yes'
            call cline%set('stream','no')
        endif
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        if( l_stream ) params_glob%stream = 'yes'
        if( cline%defined('which_iter') )then
            params%refs      = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//params%ext%to_char()
            params%refs_even = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//'_even'//params%ext%to_char()
            params%refs_odd  = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//'_odd'//params%ext%to_char()
        else if( .not. cline%defined('refs') )then
            params%refs      = 'start2Drefs'//params%ext%to_char()
            params%refs_even = 'start2Drefs_even'//params%ext%to_char()
            params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
        endif
        if( params%l_polar )then
            fname = 'cavgs_even_part'//int2str_pad(1,params%numlen)//BIN_EXT
            call polaft_dims_from_file_header(fname, pftsz, kfromto, ncls)
            call fname%kill
            call pftc%new(1, [1,1], kfromto)
            if( trim(params%ref_type)=='comlin_hybrid' )then
                call pftc%polar_cavger_new(.true., nrefs=params%ncls)
                call pftc%polar_cavger_calc_pops(build%spproj)
                call build%pgrpsyms%new('c1')
                params%nsym    = build%pgrpsyms%get_nsym()
                params%eullims = build%pgrpsyms%get_eullims()
                call build%eulspace%new(params%ncls, is_ptcl=.false.)
                call build%pgrpsyms%build_refspiral(build%eulspace)
                clw = min(1.0, max(0.0, 1.0-max(0.0, real(params_glob%extr_iter-4)/real(params_glob%extr_lim-3))))
                call pftc%polar_cavger_assemble_sums_from_parts(reforis=build%eulspace, symop=build%pgrpsyms, clin_anneal=clw)
            else
                call pftc%polar_cavger_new(.false., nrefs=params%ncls)
                call pftc%polar_cavger_calc_pops(build%spproj)
                call pftc%polar_cavger_assemble_sums_from_parts
            endif
            call terminate_stream('SIMPLE_CAVGASSEMBLE HARD STOP 1')
            call pftc%polar_cavger_calc_and_write_frcs_and_eoavg(build%clsfrcs, build%spproj_field%get_update_frac(), params%frcs, cline)
            call pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
            call pftc%polar_cavger_gen2Dclassdoc(build_glob%spproj, build%clsfrcs)
            call pftc%kill
            call pftc%polar_cavger_kill
        else
            if( L_NEW_CAVGER )then
                call cavger_new_new(build)
                call cavger_new_transf_oridat( build%spproj )
                call cavger_new_assemble_sums_from_parts()
                ! classdoc gen needs to be after calc of FRCs
                call cavger_new_gen2Dclassdoc(build%spproj) ! populates the cls2D field in project
                ! write references
                call terminate_stream('SIMPLE_CAVGASSEMBLE HARD STOP')
                call cavger_new_write_all(params%refs, params%refs_even, params%refs_odd)
                call cavger_new_kill
            else
                call cavger_new(build)
                call cavger_transf_oridat( build%spproj )
                call cavger_assemble_sums_from_parts()
                call cavger_calc_and_write_frcs_and_eoavg(params%frcs, params%which_iter)
                ! classdoc gen needs to be after calc of FRCs
                call cavger_gen2Dclassdoc(build%spproj) ! populates the cls2D field in project
                ! write references
                call terminate_stream('SIMPLE_CAVGASSEMBLE HARD STOP')
                call cavger_write(params%refs,      'merged')
                call cavger_write(params%refs_even, 'even'  )
                call cavger_write(params%refs_odd,  'odd'   )
                call cavger_kill
            endif
        endif
        ! get iteration from which_iter else from refs filename and write cavgs starfile
        if( cline%defined('which_iter') ) then
            call starproj%export_cls2D(build%spproj, params%which_iter)
        else if( cline%defined('refs') .and. params%refs%substr_ind(CAVGS_ITER_FBODY) > 0 ) then
            iterstr_start = params%refs%substr_ind(CAVGS_ITER_FBODY) + 10
            iterstr_end   = params%refs%substr_ind(params%ext) - 1
            iter = str2int(params%refs%to_char([iterstr_start,iterstr_end]), io_stat)
            call starproj%export_cls2D(build%spproj, iter)
        end if
        ! updates project
        select case(trim(params%oritype))
            case('ptcl2D')
                ! cls2D and state congruent cls3D
                call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
                states = build%spproj%os_cls2D%get_all('state')
                call build%spproj%os_cls3D%set_all('state',states)
                deallocate(states)
                call build%spproj%add_frcs2os_out( string(FRCS_FILE), 'frc2D')
                if( .not. params%l_polar )then ! no Cartesian class averages in the polar version
                    call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg')
                endif
                ! multiple fields updated, do a full write
                call build%spproj%write(params%projfile)
            case('ptcl3D')
                call build%eulspace%new(params%nspace, is_ptcl=.false.)
                call build%pgrpsyms%build_refspiral(build%eulspace)
                do icls = 1,params%ncls
                    call build%spproj%os_cls3D%set_euler(icls, build%eulspace%get_euler(icls))
                enddo
                call build%spproj%write_segment_inside('cls3D', params%projfile)
                if( cline%defined('outfile') )then
                    call build%spproj%os_cls3D%write(params%outfile)
                else
                    call build%spproj%os_cls3D%write(string('cls3D_oris.txt'))
                endif
                call build%spproj%add_frcs2os_out( string(FRCS_FILE), 'frc3D')
                if( .not. params%l_polar )then ! no Cartesian class averages in the polar version
                    call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg3D')
                endif
                call build%spproj%write_segment_inside('out', params%projfile)
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE: '//trim(params%oritype))
        end select
        ! end gracefully
        call starproj%kill
        call build%spproj%kill
        call build%kill_general_tbox
        call build%kill_strategy2D_tbox
        call simple_end('**** SIMPLE_CAVGASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('CAVGASSEMBLE_FINISHED')
    end subroutine exec_cavgassemble

    subroutine exec_write_classes( self, cline )
        class(commander_write_classes), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(image)      :: img_cavg
        type(string)     :: cavgsstk, classname, stkname
        type(image),        allocatable :: imgs_class(:)
        real,               allocatable :: states(:), inpls(:,:)
        integer,            allocatable :: pops(:), pinds(:)
        real(kind=c_float), allocatable :: rmat_rot(:,:,:)
        integer :: ncls, n, ldim(3), icls, pop_max, ind_in_stk, i, cnt
        real    :: smpd
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! get class average stack
        call spproj%get_cavgs_stk(cavgsstk, ncls, smpd)
        call find_ldim_nptcls(cavgsstk, ldim, n)
        ldim(3) = 1
        if( n /= ncls ) THROW_HARD('Incosistent # classes in project file vs cavgs stack; exec_write_classes')
        ! get state flag array
        states = spproj%os_cls2D%get_all('state')
        ! get ncls from ptcl2D field
        n = spproj%os_ptcl2D%get_n('class')
        if( n /= ncls ) THROW_HARD('Incosistent # classes in ptcl2D field of spproj vs cavgs stack; exec_write_classes')
        ! find out maximum population and allocate image arrays accordingly
        call spproj%os_ptcl2D%get_pops(pops, 'class')
        pop_max = maxval(pops)
        write(logfhandle,'(A,I5)') '>>> MAXIMUM CLASS POPULATION: ', pop_max
        allocate(imgs_class(pop_max), inpls(pop_max,3), rmat_rot(ldim(1),ldim(2),1))
        rmat_rot = 0.
        inpls    = 0.
        do i=1,pop_max
            call imgs_class(i)%new(ldim, smpd, wthreads=.false.)
        end do
        call img_cavg%new(ldim, smpd)
        ! loop over classes
        do icls=1,ncls
            if( states(icls) < 0.5 ) cycle
            ! get particle indices of class
            call spproj%os_ptcl2D%get_pinds(icls, 'class', pinds)
            if( .not. allocated(pinds) ) cycle
            ! read the class average
            call img_cavg%read(cavgsstk, icls)
            ! read the images and get the in-plane parameters
            do i=1,size(pinds)
                ! read
                call spproj%get_stkname_and_ind('ptcl2D', pinds(i), stkname, ind_in_stk)
                call imgs_class(i)%read(stkname, ind_in_stk)
                ! get params
                inpls(i,1)  = spproj%os_ptcl2D%e3get(pinds(i))
                inpls(i,2:) = spproj%os_ptcl2D%get_2Dshift(pinds(i))
            end do
            ! rotate the images (in parallel)
            !$omp parallel do default(shared) private(i,rmat_rot) schedule(static) proc_bind(close)
            do i=1,size(pinds)
                call imgs_class(i)%fft
                call imgs_class(i)%shift2Dserial([-inpls(i,2),-inpls(i,3)])
                call imgs_class(i)%ifft
                call imgs_class(i)%rtsq_serial(inpls(i,1), 0., 0., rmat_rot)
                call imgs_class(i)%set_rmat(rmat_rot,.false.)
            end do
            !$omp end parallel do
            ! make a filename for the class
            classname = 'class'//int2str_pad(icls,5)//STK_EXT
            ! write the class average first, followed by the rotated and shifted particles
            call img_cavg%write(classname, 1)
            cnt = 1
            do i=1,size(pinds)
                cnt = cnt + 1
                call imgs_class(i)%write(classname, cnt)
            end do
        end do
        ! destruct
        call spproj%kill
        call img_cavg%kill
        do i=1,size(imgs_class)
            call imgs_class(i)%kill
        end do
        deallocate(imgs_class, inpls, rmat_rot)
        if( allocated(states) ) deallocate(states)
        if( allocated(pops)   ) deallocate(pops)
        if( allocated(pinds)  ) deallocate(pinds)
        ! end gracefully
        call simple_end('**** SIMPLE_WRITE_CLASSES NORMAL STOP ****')
    end subroutine exec_write_classes

end module simple_commanders_mkcavgs

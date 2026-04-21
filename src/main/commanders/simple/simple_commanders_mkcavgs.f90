!@descr: for producing class averages
module simple_commanders_mkcavgs
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
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
        call run_make_cavgs_workflow(cline, from_distr_cmd=.true.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_make_cavgs( self, cline )
        class(commander_make_cavgs), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call run_make_cavgs_workflow(cline, from_distr_cmd=.false.)
    end subroutine exec_make_cavgs

    ! ------------------------------------------------------------------
    ! Unified runtime-polymorphic workflow
    ! ------------------------------------------------------------------

    subroutine run_make_cavgs_workflow( cline, from_distr_cmd )
        use simple_make_cavgs_strategy, only: make_cavgs_strategy, make_cavgs_hooks, create_make_cavgs_strategy
        use simple_cmdline,             only: cmdline
        use simple_parameters,          only: parameters
        class(cmdline), intent(inout) :: cline
        logical,        intent(in)    :: from_distr_cmd
        class(make_cavgs_strategy), allocatable :: strategy
        type(make_cavgs_hooks) :: hooks
        type(parameters) :: params
        ! Ensure distributed scripts see the correct program name.
        call cline%set('prg',    'make_cavgs')
        call cline%set('oritype','ptcl2D')
        ! Provide master-side hook for cavgassemble (avoids strategy importing commander modules).
        hooks%run_cavgassemble => make_cavgs_exec_cavgassemble
        strategy = create_make_cavgs_strategy(cline, hooks, from_distr_cmd=from_distr_cmd)
        call strategy%apply_defaults(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message(), print_simple=.false.)
        ! exec_make_cavgs_distr behavior: async touch marker
        call strategy%after_end(params, cline)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine run_make_cavgs_workflow

    ! ----------------------------------------------------------------------
    ! Hook implementation: run cavgassemble using the existing commander
    ! ----------------------------------------------------------------------

    subroutine make_cavgs_exec_cavgassemble( cline, nthr )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: nthr
        type(commander_cavgassemble) :: xcavgassemble
        type(cmdline)               :: cline_cavgassemble
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg',  'cavgassemble')
        call cline_cavgassemble%set('nthr', nthr)
        call xcavgassemble%execute(cline_cavgassemble)
        call cline_cavgassemble%kill
    end subroutine make_cavgs_exec_cavgassemble

    subroutine exec_cavgassemble( self, cline )
        class(commander_cavgassemble), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)   :: params
        type(builder)      :: build
        type(starproject)  :: starproj
        type(polarft_calc) :: pftc
        type(string)       :: fname
        real, allocatable  :: states(:)
        integer            :: iterstr_start, iterstr_end, iter, io_stat
        integer            :: pftsz, kfromto(2), ncls
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
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
            call polarft_dims_from_file_header(fname, pftsz, kfromto, ncls)
            call fname%kill
            call pftc%new(params, 1, [1,1], kfromto)
            call pftc%polar_cavger_new(.false., nrefs=params%ncls)
            call pftc%polar_cavger_calc_pops(build%spproj)
            call pftc%polar_cavger_assemble_sums_from_parts
            call pftc%polar_cavger_merge_eos_and_norm2D(build%clsfrcs, params%frcs)
            call pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
            call pftc%polar_cavger_gen2Dclassdoc(build%spproj, build%clsfrcs)
            call pftc%kill
            call pftc%polar_cavger_kill
        else
            call cavger_new(params, build)
            call cavger_assemble_sums_from_parts
            call cavger_gen2Dclassdoc
            call terminate_stream(params, 'SIMPLE_CAVGASSEMBLE HARD STOP')
            call cavger_write_all(params%refs, params%refs_even, params%refs_odd)
            call cavger_kill
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

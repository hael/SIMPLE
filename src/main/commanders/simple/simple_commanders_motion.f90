!@descr: motion correction and movie fractionation commanders
module simple_commanders_motion
include "starfile_enum.inc"
use simple_commanders_api
use simple_gui_communicator, only: gui_communicator
implicit none

public :: commander_refine_motion_model
public :: commander_fractionate_movies_distr
public :: commander_fractionate_movies
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_refine_motion_model
  contains
    procedure :: execute       => exec_refine_motion_model
end type commander_refine_motion_model

type, extends(commander_base) :: commander_fractionate_movies_distr
  contains
    procedure :: execute       => exec_fractionate_movies_distr
end type commander_fractionate_movies_distr

type, extends(commander_base) :: commander_fractionate_movies
  contains
    procedure :: execute       => exec_fractionate_movies
end type commander_fractionate_movies

contains

    subroutine exec_refine_motion_model( self, cline )
        use simple_core_module_api,                    only: simple_end
        use simple_refine_motion_model_strategy,      only: refine_motion_model_strategy, &
                                                               create_refine_motion_model_strategy
        use simple_cmdline,                            only: cmdline
        use simple_parameters,                         only: parameters
        class(commander_refine_motion_model), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        class(refine_motion_model_strategy), allocatable :: strategy
        type(parameters)                                   :: params
        type(gui_communicator)                             :: gui_comm
        call cline%set('prg', 'refine_motion_model')
        strategy = create_refine_motion_model_strategy(cline)
        call strategy%apply_defaults(cline)
        call strategy%initialize(params, cline)
        call gui_comm%new(params)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call gui_comm%add_metadata(params%projfile, oritype='ptcl')
        call gui_comm%kill()
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_refine_motion_model

    subroutine exec_fractionate_movies_distr( self, cline )
        use simple_starproject,          only: starproject
        use simple_motion_correct_utils, only: flip_gain
        class(commander_fractionate_movies_distr), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        type(parameters)  :: params
        type(sp_project)  :: spproj
        type(chash)       :: job_descr
        type(qsys_env)    :: qenv
        type(starproject) :: starproj
        integer           :: nmovies
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        if( .not.cline%defined('mcconvention') ) call cline%set('mcconvention', 'simple')
        if( .not.cline%defined('fromf') )        call cline%set('fromf',        1)
        if( .not.cline%defined('tof') )          call cline%set('tof',          0)
        call params%new(cline)
        call spproj%read_segment(params%oritype, params%projfile)
        ! sanity checks
        if( (params%fromf < 1) ) THROW_HARD('Invalid fractions range!')
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 ) THROW_HARD('No movie to process!')
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! processing gain reference if needed
        call flip_gain(cline, params%gainref, params%flipgain)
        ! setup the environment for distributed execution
        params%nparts = min(nmovies, params%nparts)
        call cline%set('nparts', params%nparts)
        call qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call starproj%export_mics(spproj)
        ! cleanup
        call qsys_cleanup(params)
        call spproj%kill
        call starproj%kill
        call simple_end('**** SIMPLE_FRACTIONATE_MOVIES_DISTR NORMAL STOP ****')
    end subroutine exec_fractionate_movies_distr

    subroutine exec_fractionate_movies( self, cline )
        use simple_micrograph_generator
        use simple_fsc, only: plot_fsc
        class(commander_fractionate_movies), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        logical,            parameter :: L_DEBUG = .false.
        type(string)                  :: mic_fname,forctf_fname, ext, mov_fname
        type(string)                  :: mic_fbody, star_fname, background_fname
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(mic_generator)           :: generator
        type(ori)                     :: o
        type(image)                   :: micrograph_dw, micrograph_nodw, mic, background
        real,             allocatable :: frc(:), res(:)
        type(string)                  :: orig_mic
        integer :: nmovies, imov, cnt, n
        call cline%set('mkdir',   'no')
        call cline%set('oritype', 'mic')
        if( .not.cline%defined('mcconvention') ) call cline%set('mcconvention', 'simple')
        if( .not.cline%defined('fromf') )        call cline%set('fromf',        1)
        if( .not.cline%defined('tof') )          call cline%set('tof',          0)
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity checks
        if( (params%fromf < 1) ) THROW_HARD('Invalid fractions range!')
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 ) THROW_HARD('No movie to process!')
        ! Main loop
        cnt = 0
        do imov = params%fromp,params%top
            call spproj%os_mic%get_ori(imov, o)
            if( .not.o%isthere('movie') ) cycle
            if( .not.o%isthere('intg')  ) cycle
            if( o%get_state() == 0 ) cycle
            cnt = cnt + 1
            orig_mic = o%get_str('intg')
            ! new micrograph
            call generator%new(o, params%mcconvention, [params%fromf, params%tof], params%dw=='yes')
            select case(trim(params%mcconvention))
            case('cs')
                call generator%generate_micrographs(micrograph_dw, micrograph_nodw, background=background)
            case DEFAULT
                call generator%generate_micrographs(micrograph_dw, micrograph_nodw)
            end select
            ! file naming
            mov_fname = generator%get_moviename()
            mic_fbody = basename(mov_fname)
            ext       = fname2ext(mic_fbody)
            mic_fbody = get_fbody(mic_fbody, ext)
            select case(trim(params%mcconvention))
            case('simple')
                mic_fname    = mic_fbody//INTGMOV_SUFFIX//params%ext%to_char()
                forctf_fname = mic_fbody//FORCTF_SUFFIX //params%ext%to_char()
            case('motioncorr', 'relion')
                mic_fname    = mic_fbody//params%ext%to_char()
                forctf_fname = mic_fbody//'_noDW'//params%ext%to_char()
            case('cryosparc','cs')
                mic_fname    = mic_fbody//'_patch_aligned_doseweighted'//params%ext%to_char()
                forctf_fname = mic_fbody//'_patch_aligned'             //params%ext%to_char()
            case DEFAULT
                THROW_HARD('Unsupported convention!')
            end select
            star_fname = mic_fbody//STAR_EXT
            ! write
            if( .not.micrograph_dw%exists() )then
                ! doses not defined
                mic_fname = forctf_fname
                call micrograph_nodw%write(mic_fname)
            else
                call micrograph_dw%write(mic_fname)
                call micrograph_nodw%write(forctf_fname)
            endif
            if( background%exists() )then
                background_fname = mic_fbody//'_background'//params%ext%to_char()
                call background%write(background_fname)
            endif
            call generator%write_star(star_fname%to_char())
            ! parameters update
            call o%set('intg',        simple_abspath(mic_fname))
            call o%set('forctf',      simple_abspath(forctf_fname))
            call o%set('mc_starfile', simple_abspath(star_fname))
            call o%set('imgkind',     'mic')
            call o%set('smpd',        micrograph_nodw%get_smpd())
            call o%delete_entry('thumb')
            call spproj%os_mic%set_ori(imov, o)
            if( L_DEBUG )then
                call mic%copy(micrograph_dw)
                call mic%read(orig_mic)
                call mic%fft
                call micrograph_dw%fft
                n = fdim(mic%get_box())-1
                allocate(frc(n))
                res = mic%get_res()
                call mic%fsc(micrograph_dw, frc)
                call plot_fsc(n, frc, res, o%get('smpd'), mic_fbody%to_char())
                deallocate(frc,res)
                call mic%kill
            endif
            ! tidy
            call micrograph_dw%kill
            call micrograph_nodw%kill
            call background%kill
        enddo
        call generator%kill
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, [params%fromp,params%top], isegment=MIC_SEG)
        call spproj%kill
        call qsys_job_finished(params, string('simple_commanders_motion :: exec_fractionate_movies'))
        call simple_end('**** SIMPLE_FRACTIONATE_MOVIES NORMAL STOP ****')
    end subroutine exec_fractionate_movies

end module simple_commanders_motion

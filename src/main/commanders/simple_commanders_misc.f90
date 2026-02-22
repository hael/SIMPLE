!@descr: stuff that didn't fit elsewhere
module simple_commanders_misc
include "starfile_enum.inc"
use simple_commanders_api
use simple_simple_volinterp, only: rotvol
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_print_fsc
  contains
    procedure :: execute       => exec_print_fsc
end type commander_print_fsc

type, extends(commander_base) :: commander_print_magic_boxes
  contains
    procedure :: execute       => exec_print_magic_boxes
end type commander_print_magic_boxes

type, extends(commander_base) :: commander_print_dose_weights
  contains
    procedure :: execute       => exec_print_dose_weights
end type commander_print_dose_weights

type, extends(commander_base) :: commander_kstest
  contains
    procedure :: execute       => exec_kstest
end type commander_kstest

type, extends(commander_base) :: commander_pearsn
  contains
    procedure :: execute       => exec_pearsn
end type commander_pearsn

type, extends(commander_base) :: commander_mkdir
  contains
    procedure :: execute       => exec_mkdir
end type commander_mkdir

type, extends(commander_base) :: commander_fractionate_movies_distr
  contains
    procedure :: execute       => exec_fractionate_movies_distr
end type commander_fractionate_movies_distr

type, extends(commander_base) :: commander_fractionate_movies
  contains
    procedure :: execute       => exec_fractionate_movies
end type commander_fractionate_movies

contains

    !>  for printing the binary FSC files produced by PRIME3D
    subroutine exec_print_fsc( self, cline )
        use simple_fsc,        only: plot_fsc
        use simple_class_frcs, only: class_frcs
        class(commander_print_fsc), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        real,  allocatable :: res(:), fsc(:)
        type(parameters) :: params
        type(image)      :: img
        type(class_frcs) :: frcs
        type(string)     :: tmpl_fname
        integer          :: k,n
        real             :: res0143, res05
        logical :: l_fsc, l_frcs
        call params%new(cline)
        l_fsc  = cline%defined('fsc')
        l_frcs = cline%defined('frcs')
        if( .not.l_fsc .and. .not.l_frcs ) THROW_HARD('FSC or FRCS must be defined!')
        if( l_fsc )then
            if( .not.cline%defined('smpd') .or. .not.cline%defined('box')) then
                THROW_HARD('SMPD and BOX must be defined!')
            endif
            call img%new([params%box,params%box,1], params%smpd)
            res = img%get_res()
            fsc = file2rarr(params%fsc)
            n = size(fsc)
            do k=1,n
                write(logfhandle,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:', res(k), '>>> FSC:', fsc(k)
            end do
            ! get & print resolution
            call get_resolution(fsc, res, res05, res0143)
            write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
            write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
            ! plot
            tmpl_fname = get_fbody(params%fsc,BIN_EXT,separator=.false.)
            call plot_fsc(n, fsc, res, params%smpd, tmpl_fname%to_char())
            call img%kill
        endif
        if( l_frcs )then
            call frcs%read(params%frcs)
            call frcs%plot_frcs(string('frcs'))
            call frcs%kill
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_FSC NORMAL STOP ****')
    end subroutine exec_print_fsc

    !> for printing magic box sizes (fast FFT)
    subroutine exec_print_magic_boxes( self, cline )
        class(commander_print_magic_boxes), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call print_magic_box_range(params%smpd, params%moldiam )
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_MAGIC_BOXES NORMAL STOP ****')
    end subroutine exec_print_magic_boxes

    subroutine exec_print_dose_weights( self, cline )
        class(commander_print_dose_weights), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)  :: params
        real, allocatable :: weights(:,:), res(:)
        integer           :: iframe, k, filtsz
        call params%new(cline)
        call calc_dose_weights(params%nframes, params%box, params%smpd, params%kV, params%total_dose, weights)
        filtsz = size(weights, dim=2)
        res = get_resarr(params%box, params%smpd)
        write(logfhandle,'(A)') 'RESOLUTION, DOSE_WEIGHTS'
        do k = 1,filtsz
            write(logfhandle, '(F7.1,A)', advance='no') res(k), ', '
            do iframe = 1,params%nframes - 1
                write(logfhandle, '(f3.1,A)', advance='no') weights(iframe,k), ', '
            end do
            write(logfhandle, '(f3.1,1X)') weights(iframe,k)
        end do
        write(logfhandle,*)
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_DOSE_WEIGHTS_NORMAL STOP ****')
    end subroutine exec_print_dose_weights

    subroutine exec_kstest( self, cline )
        class(commander_kstest), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)   :: params
        integer            :: ndat1, ndat2
        real               :: ksstat, prob, ave1, sdev1, var, ave2, sdev2
        real, allocatable  :: dat1(:), dat2(:)
        logical            :: err
        call params%new(cline)
        call read_nrs_dat(params%infile,  dat1, ndat1)
        call read_nrs_dat(params%infile2, dat2, ndat2)
        write(logfhandle,'(a)') '>>> STATISTICS OF THE TWO DISTRIBUTIONS'
        call moment(dat1, ave1, sdev1, var, err)
        call moment(dat2, ave2, sdev2, var, err)
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile : ', ave1, sdev1
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile2: ', ave2, sdev2
        write(logfhandle,'(a)') '>>> KOLMOGOROV-SMIRNOV TEST TO DEDUCE EQUIVALENCE OR NON-EQUIVALENCE BETWEEN TWO DISTRIBUTIONS'
        call kstwo(dat1, ndat1, dat2, ndat2, ksstat, prob)
        write(logfhandle,'(a,1x,f4.2)') 'K-S statistic = ', ksstat
        write(logfhandle,'(a,1x,f4.2)') 'P             = ', prob
        write(logfhandle,'(a)') 'P represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution'
        write(logfhandle,'(a)') 'Small P values show that the cumulative distribution functions of the two data sets differ significantly'
        ! end gracefully
        call simple_end('**** SIMPLE_KSTEST NORMAL STOP ****')
    end subroutine exec_kstest

    subroutine exec_pearsn( self, cline )
        class(commander_pearsn), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)   :: params
        integer            :: ndat1, ndat2
        real               :: corr, ave1, sdev1, var, ave2, sdev2
        real, allocatable  :: dat1(:), dat2(:)
        logical            :: err
        call params%new(cline)
        call read_nrs_dat(params%infile,  dat1, ndat1)
        call read_nrs_dat(params%infile2, dat2, ndat2)
        if( ndat1 /= ndat2 ) THROW_HARD('Input distributions not identical')
        call moment(dat1, ave1, sdev1, var, err)
        call moment(dat2, ave2, sdev2, var, err)
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile : ', ave1, sdev1
        write(logfhandle,'(a,1x,f4.2,1x,f4.2)') 'mean & sdev for infile2: ', ave2, sdev2
        write(logfhandle,'(a)') '>>> PEARSON CORRELATION OF THE TWO DISTRIBUTIONS'
        corr = pearsn(dat1, dat2)
        write(logfhandle,'(a,1x,f4.2)') 'CORR = ', corr
        write(logfhandle,'(a)') 'P represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution'
        write(logfhandle,'(a)') 'Small P values show that the cumulative distribution functions of the two data sets differ significantly'
        ! end gracefully
        call simple_end('**** SIMPLE_PEARSN NORMAL STOP ****')
    end subroutine exec_pearsn

    subroutine exec_mkdir( self, cline )
        class(commander_mkdir), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        call cline%set('mkdir', 'yes')
        call params%new(cline)
    end subroutine exec_mkdir

    subroutine exec_fractionate_movies_distr( self, cline )
        use simple_starproject, only: starproject
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
            call generator%new(o, params%mcconvention, [params%fromf, params%tof])
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
                call micrograph_nodw%write(mic_fname)
                forctf_fname = mic_fname
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
        call qsys_job_finished(params, string('simple_commander_mic :: exec_fractionate_movies'))
        call simple_end('**** SIMPLE_FRACTIONATE_MOVIES NORMAL STOP ****')
    end subroutine exec_fractionate_movies

    ! utility routines

    subroutine read_nrs_dat( filename, arr, ndat )
        class(string),     intent(in)  :: filename
        real, allocatable, intent(out) :: arr(:)
        integer,           intent(out) :: ndat
        integer            :: ndatlines, nrecs, i, j, cnt
        real, allocatable  :: line(:)
        type(nrtxtfile)    :: nrsfile
        call nrsfile%new(filename, 1)
        ndatlines = nrsfile%get_ndatalines()
        nrecs     = nrsfile%get_nrecs_per_line()
        ndat = ndatlines * nrecs
        allocate( line(nrecs), arr(ndat) )
        cnt = 0
        do i=1,ndatlines
            call nrsfile%readNextDataLine(line)
            do j=1,nrecs
                cnt = cnt + 1
                arr(cnt) = line(j)
            end do
        end do
        deallocate (line)
        call nrsfile%kill
    end subroutine read_nrs_dat

end module simple_commanders_misc

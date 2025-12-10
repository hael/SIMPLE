! concrete commander: miscallenaous routines
module simple_commanders_misc
include 'simple_lib.f08'
include "starfile_enum.inc"
use simple_binoris_io
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_simple_volinterp, only: rotvol
use simple_sp_project,       only: sp_project
use simple_image,            only: image
use simple_builder,          only: builder
use simple_parameters,       only: parameters
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

type, extends(commander_base) :: commander_comparemc
  contains
    procedure :: execute       => exec_comparemc
end type commander_comparemc

type, extends(commander_base) :: commander_afm
  contains
    procedure :: execute       => exec_afm
end type commander_afm

contains

    subroutine exec_afm( self, cline )
        use simple_pickseg
        use simple_image_afm
        use simple_corrmat
        use simple_aff_prop
        use simple_srch_sort_loc
        class(commander_afm), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        type(parameters), target     :: params
        type(image),     allocatable :: pick_vec(:), ordered_pick_vec(:)
        type(string)    :: directory, test_file, sim_dir
        integer                 :: i, nptcls, temp_ldim(3), pick_dim(3), cropped_dim(3), test_count
        integer, allocatable    :: medoids(:), labels(:), labels_ind(:)
        real,    allocatable    :: corrmat(:, :), var_mat(:, :), labels_real(:)
        logical, allocatable    :: corrmat_mask(:, :)
        type(aff_prop)          :: clus
        real        :: sim_sum
        directory = '/Users/atifao/Downloads/IBW_orig/'
        test_file = '/Users/atifao/Downloads/MRC_T/17.mrc'
        sim_dir   = '/Users/atifao/Downloads/mrc_for_clus.mrc'
        params_glob => params
        params_glob%pcontrast = 'white'
        params_glob%lp  = 10.
        params_glob%nsig  = 1.5 
        call cline%set('objfun','cc')
        call cline%set('ctf',    'no')
        call cline%set('sh_inv',  'yes')
        call cline%set('objfun', 'cc')
        call cline%set('mkdir', 'no')
        call cline%set('lambda', 0.)
        call cline%set('trs',     50.0)
        call cline%set('box',     150)
        call cline%set('smpd',    5.0)
        params_glob%cc_objfun = 0
        params_glob%maxits_sh = 200
        params_glob%shbarrier = 'yes'
        call params%new(cline)
        call find_ldim_nptcls(sim_dir,temp_ldim, nptcls)
        cropped_dim = [temp_ldim(1), temp_ldim(2) , 1]
        params_glob%ldim = cropped_dim
        params_glob%box = cropped_dim(1)
        test_count = temp_ldim(3)
        ! test_count   = 100
        allocate(pick_vec(test_count))
        allocate(ordered_pick_vec(test_count))
        do i = 1, temp_ldim(3) 
            pick_dim = [temp_ldim(1), temp_ldim(2), 1]
            call pick_vec(i)%new(pick_dim, params%smpd, .false.)
            call pick_vec(i)%read(sim_dir, i)
            call pick_vec(i)%clip_inplace(cropped_dim) 
            if(i > test_count - 1) exit 
        end do
        allocate(var_mat(test_count, test_count), source = 0.)
        ! sklearn uses negative squared euclidean distance 
        corrmat = 2.*(1.-corrmat)
        where( corrmat < 0. ) corrmat = 0.
        where( corrmat > 4. ) corrmat = 4.
        corrmat = -corrmat
        ! initial clustering 
        allocate(corrmat_mask(test_count,test_count), source = .true.)
        allocate(labels_ind(test_count))
        do i = 1, test_count
            labels_ind(i) = i
        end do
        allocate(labels_real(test_count))
        call clus%new(test_count, corrmat, pref = median(pack(corrmat, corrmat_mask) / 1.5))
        call clus%propagate(medoids, labels, sim_sum)
        labels_real = real(labels)
        call hpsort_1(labels_real, labels_ind)
        do i = 1, test_count
            call ordered_pick_vec(i)%copy(pick_vec(labels_ind(i)))
            call ordered_pick_vec(i)%ifft()
            call ordered_pick_vec(i)%write(string('ordered_picks.mrc'), i)
        end do 
    end subroutine exec_afm

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

    !>  dsym_cylinder search intended for symmetry of order D
    subroutine exec_dsym_volinit( dsym_os, cylinder)
        use simple_parameters,   only: params_glob
        use simple_image,        only: image
        use simple_segmentation, only: otsu_robust_fast
        class(oris),   intent(inout) :: dsym_os
        class(image),  intent(inout) :: cylinder
        type(image)          :: read_img, img_msk, dist_img, roavg_img, topview
        type(sym)            :: se
        real,    allocatable :: corrs(:), forsort(:), e2(:), radii(:)
        integer, allocatable :: labels(:)
        logical, allocatable :: l_msk(:,:,:)
        integer  :: halfnoris, cnt1, cnt2, i, l, noris
        real     :: cen1, cen2, sum1, sum2, sumvals, sdev_noise
        real     :: minmax(2), width, height, sh1(3), ang, thresh(3)
        if(params_glob%pgrp(1:1).eq.'d' .or. params_glob%pgrp(1:1).eq.'D')then
            ang = 360. / real(se%get_nsym()/2)
        else if(params_glob%pgrp(1:1).eq.'c' .or. params_glob%pgrp(1:1).eq.'C')then
            ang = 360. / real(se%get_nsym())
        else
            THROW_HARD('Not intended for platonic point-groups')
        endif
        ! init
        noris = dsym_os%get_noris()
        allocate(corrs(noris),  source=-1.)
        allocate(radii(noris),e2(noris),  source=0.)
        allocate(labels(noris), source=0)
        call se%new(params_glob%pgrp)
        call cylinder%new( [params_glob%box, params_glob%box, params_glob%box], params_glob%smpd)
        call read_img%new( [params_glob%box, params_glob%box, 1], params_glob%smpd)
        call topview%new( [params_glob%box, params_glob%box, 1], params_glob%smpd)
        call roavg_img%new([params_glob%box, params_glob%box, 1], params_glob%smpd)
        call dist_img%new( [params_glob%box, params_glob%box, 1], params_glob%smpd)
        call dist_img%cendist
        ! prep mask
        call img_msk%new([params_glob%box,params_glob%box,1],params_glob%smpd)
        img_msk = 1.
        call img_msk%mask(params_glob%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        ! centers, calculates self to rotational averages images & radii
        do i = 1, noris
            call read_img%read(params_glob%stk,i)
            call read_img%norm_noise(l_msk, sdev_noise)
            ! center image
            sh1 = read_img%calc_shiftcen(params_glob%cenlp, params_glob%msk)
            call read_img%shift(-sh1)
            call dsym_os%set_shift(i, -sh1)
            ! rotational image
            call read_img%roavg(nint(ang), roavg_img)
            call roavg_img%write(string('roavg.mrc'),i)
            corrs(i) = read_img%real_corr(roavg_img, l_msk)
            ! radii
            call read_img%bp(0., params_glob%cenlp)
            call read_img%mask(params_glob%msk, 'hard')
            call otsu_robust_fast(read_img, is2d=.true., noneg=.false., thresh=thresh)
            call read_img%write(string('bin.mrc'),i)
            call read_img%mul(dist_img)
            minmax = read_img%minmax()
            radii(i) = min(params_glob%msk, minmax(2))
        enddo
        call dsym_os%set_all('corr', corrs)
        ! indentify views along z and x-/y- axes
        forsort = corrs
        call hpsort(forsort)
        halfnoris = nint(real(noris)/2.)
        cen1      = sum(forsort(1:halfnoris)) / real(halfnoris)
        cen2      = sum(forsort(halfnoris+1:noris)) / real(noris-halfnoris)
        sumvals   = sum(forsort)
        ! do 100 iterations of k-means
        do l = 1, 100
            sum1 = 0.
            cnt1 = 0
            do i = 1, noris
                if( (cen1-forsort(i))**2. < (cen2-forsort(i))**2. )then
                    cnt1 = cnt1 + 1
                    sum1 = sum1 + forsort(i)
                endif
            end do
            cnt2 = noris - cnt1
            sum2 = sumvals - sum1
            cen1 = sum1 / real(cnt1)
            cen2 = sum2 / real(cnt2)
        end do
        ! label views
        if(cen1 > cen2)then
            ! cen1: along z-axis
            labels = 2
            where( (cen1-corrs)**2. < (cen2-corrs)**2. )labels = 1
        else
            labels = 1
            where( (cen1-corrs)**2. < (cen2-corrs)**2. )labels = 2
        endif
        e2 = 90.
        where(labels == 1) e2 = 0.
        call dsym_os%set_all('e2', e2)
        call dsym_os%set_all('class', real(labels))
        ! rotates input oris to asymmetric unit
        call se%rotall_to_asym(dsym_os)
        ! estimate height and cylinder radius (0.9 to account for overestimation)
        width  = 0.9 * (sum(radii, mask=(labels==1)) / real(count(labels==1)))
        height = 0.9 * (2. * sum(radii, mask=(labels==2)) / real(count(labels==2)))
        ! dummy top view
        topview = 0.
        do i=1,noris
            if(labels(i)==1)then
                call read_img%read(params_glob%stk,i)
                call read_img%norm_noise(l_msk, sdev_noise)
                call read_img%rtsq(0., -dsym_os%get(i,'x'), -dsym_os%get(i,'y'))
                call topview%add(read_img)
            endif
        enddo
        call topview%div( real(count(labels==1)) )
        call topview%roavg(nint(ang), roavg_img)
        topview = roavg_img
        call topview%norm()
        call topview%mask(params_glob%msk, 'soft')
        cylinder = 0.
        do i=1,params_glob%box
            if( abs( real(i-1)-real(params_glob%box)/2. ) < height/2.)then
                call cylinder%set_slice( i, topview )
            endif
        enddo
        ! cleanup
        deallocate(corrs, e2, radii, labels)
        call dist_img%kill
        call read_img%kill
        call se%kill
        call roavg_img%kill
        call topview%kill
    end subroutine exec_dsym_volinit

    subroutine exec_mkdir( self, cline )
        class(commander_mkdir), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        call cline%set('mkdir', 'yes')
        call params%new(cline)
    end subroutine exec_mkdir

    subroutine exec_fractionate_movies_distr( self, cline )
        use simple_starproject, only: starproject
        use simple_qsys_env,    only: qsys_env
        use simple_qsys_funs
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
        if( .not.cline%defined('interpfun') )    call cline%set('interpfun',    'linear')
        call params%new(cline)
        call spproj%read_segment(params%oritype, params%projfile)
        ! sanity checks
        if( (params%fromf < 1) ) THROW_HARD('Invalid fractions range!')
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 ) THROW_HARD('No movie to process!')
        call spproj%kill
        select case(trim(params%interpfun))
        case('linear', 'nn')
            ! all good
        case DEFAULT
            THROW_HARD('Invalid interpolation scheme: '//trim(params%interpfun))
        end select
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        params%nparts = min(nmovies, params%nparts)
        call cline%set('nparts', params%nparts)
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call starproj%export_mics(spproj)
        ! cleanup
        call qsys_cleanup
        call spproj%kill
        call starproj%kill
        call simple_end('**** SIMPLE_FRACTIONATE_MOVIES_DISTR NORMAL STOP ****')
    end subroutine exec_fractionate_movies_distr

    subroutine exec_fractionate_movies( self, cline )
        use simple_micrograph_generator
        use simple_qsys_funs
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
        logical :: l_bilinear_interp
        call cline%set('mkdir',   'no')
        call cline%set('oritype', 'mic')
        if( .not.cline%defined('mcconvention') ) call cline%set('mcconvention', 'simple')
        if( .not.cline%defined('fromf') )        call cline%set('fromf',        1)
        if( .not.cline%defined('tof') )          call cline%set('tof',          0)
        if( .not.cline%defined('interpfun') )    call cline%set('interpfun',    'linear')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity checks
        if( (params%fromf < 1) ) THROW_HARD('Invalid fractions range!')
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 ) THROW_HARD('No movie to process!')
        select case(trim(params%interpfun))
        case('linear')
            l_bilinear_interp = .true.
        case('nn')
            l_bilinear_interp = .false.
        case DEFAULT
            THROW_HARD('Invalid interpolation scheme: '//trim(params%interpfun))
        end select
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
            call generator%new(o, params%mcconvention, [params%fromf, params%tof], l_bilinear_interp)
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
        call qsys_job_finished(string('simple_commander_mic :: exec_fractionate_movies'))
        call simple_end('**** SIMPLE_FRACTIONATE_MOVIES NORMAL STOP ****')
    end subroutine exec_fractionate_movies

    subroutine exec_comparemc( self, cline )
        use simple_starfile_wrappers
        use CPlot2D_wrapper_module
        integer, PARAMETER :: POLYDIM = 18
        integer, PARAMETER :: NGRID   = 25 ! has to be odd
        class(commander_comparemc), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        real,         allocatable :: rmsds(:), gofs(:,:), avg(:), std(:), npatch(:), rmsds2(:), avg2(:), std2(:)
        real(dp),     allocatable :: poly1(:,:), poly2(:,:), x(:), y(:), t(:), offsets1(:,:), offsets2(:,:)
        type(string), allocatable :: projects_fnames(:), moviedocs(:)
        integer,      allocatable :: indextab(:,:), order(:), nxpatch(:), nypatch(:)
        logical,      allocatable :: mask(:)
        type(string)              :: movie, moviedoc
        type(CPlot2D_type)        :: plot
        type(CDataSet_type)       :: dataSet
        type(CDataPoint_type)     :: point
        type(parameters)          :: params
        type(sp_project)          :: first_spproj, spproj, spproj2
        real     :: smpd, binning
        integer  :: i, j, k, npoints, ldim1(3), ldim2(3), nmics, nprojects, nx, ny, nframes, funit, io_stat, cnt
        logical  :: found
        if(.not.cline%defined('mkdir')) call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call read_filetable( params%infile, projects_fnames )
        if(.not.allocated(projects_fnames))then
            THROW_HARD('FILETAB invalid format')
        endif
        nprojects = size(projects_fnames,1)
        write(logfhandle,'(A,I6)')'>>> # of project files:', nprojects
        do i = 1,nprojects
            if( projects_fnames(i)%to_char([1,1]).ne.'/' )then
                if( trim(params%mkdir).eq.'yes')then
                    projects_fnames(i) = string('../')//projects_fnames(i)
                endif
                projects_fnames(i) = simple_abspath(projects_fnames(i))
            endif
            if( .not.file_exists(projects_fnames(i))) THROW_HARD('Could not find: '//projects_fnames(i)%to_char())
        enddo
        call first_spproj%read_segment(params%oritype, projects_fnames(1))
        nmics = first_spproj%get_nintgs()
        write(logfhandle,'(A,I6)')'>>> # of micrographs  :', nmics
        allocate(indextab(nprojects,nmics),source=0)
        allocate(moviedocs(nmics),gofs(nprojects,nmics),npatch(nprojects),nxpatch(nprojects),nypatch(nprojects))
        nxpatch(1) = nint(first_spproj%os_mic%get(1,'nxpatch'))
        nypatch(1) = nint(first_spproj%os_mic%get(1,'nypatch'))
        gofs = 0.0
        ! working out order of movies in each project file
        do i = 1,nmics
            indextab(1,i) = i
            moviedocs(i)  = first_spproj%os_mic%get_str(i,'mc_starfile')
            gofs(1,i)     = sqrt(first_spproj%os_mic%get(i,'gofx')**2.0+first_spproj%os_mic%get(i,'gofy')**2.0)
        enddo
        do i = 2,nprojects
            call spproj%read_segment(params%oritype, projects_fnames(i))
            do j = 1,nmics
                found = .false.
                do k = 1,nmics
                    moviedoc = basename(spproj%os_mic%get_str(k,'mc_starfile'))
                    if( moviedoc.eq.basename(moviedocs(j)) )then
                        indextab(i,j) = k
                        gofs(i,j) = sqrt(spproj%os_mic%get(k,'gofx')**2.0+spproj%os_mic%get(k,'gofy')**2.0)
                        found = .true.
                    endif
                    if( found ) exit
                enddo
                if( .not.found )then
                    indextab(:,j) = 0
                    THROW_WARN('Could not find corresponding doc: '//moviedoc%to_char())
                endif
            enddo
            nxpatch(i) = nint(spproj%os_mic%get(indextab(i,1),'nxpatch'))
            nypatch(i) = nint(spproj%os_mic%get(indextab(i,1),'nypatch'))
        enddo
        npatch = sqrt(real(nxpatch*nypatch))
        ! order
        allocate(order(nprojects),avg(nprojects),std(nprojects))
        order = (/(i,i=1,nprojects)/)
        call hpsort(npatch, order)
        ! goodness of fit
        do i = 1,nprojects
            j = order(i)
            mask   = indextab(j,:) > 0
            avg(i) = sum(gofs(j,:), mask=mask) / real(count(mask))
            std(i) = sqrt(sum((gofs(j,:)-avg(i))**2.0)/real(count(mask)))
        enddo
        ! csv output
        call fopen(funit, string('goodnessoffit.txt'), action='READWRITE', status='UNKNOWN', iostat=io_stat)
        call fileiochk('could not write goodnessoffit.txt', io_stat)
        write(funit,'(A)')'sp_project,nx,ny,npatches,gof_avg,gof_std'
        do i = 1,nprojects
            write(funit,'(A,A1,I3,A1,I3,A1,F6.1,A1,F8.3,A1,F8.3)') projects_fnames(order(i))%to_char(),&
                &',',nxpatch(order(i)),',',nypatch(order(i)),',',npatch(i),',',avg(i),',',std(i)
        enddo
        call fclose(funit)
        ! eps output
        call CPlot2D__new(plot, 'Goodness of  fit'//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot, real(50*nprojects,dp))
        call CPlot2D__SetYAxisSize(plot, real(50*maxval(avg+2.0*std),dp))
        call CPlot2D__SetXAxisTitle(plot, 'n = sqrt(nx * ny)'//c_null_char)
        call CPlot2D__SetYAxisTitle(plot, 'shifts RMS patch vs poly (pixels) '//c_null_char)
        call CPlot2D__SetDrawLegend(plot, C_FALSE)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet,C_TRUE)
        call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 0.0_c_double)
        do i = 1,nprojects
            call CDataPoint__new2(real(npatch(i),dp), real(avg(i), dp), point)
            call CDataSet__AddDataPoint(dataSet, point)
            call CDataPoint__delete(point)
        enddo
        call CPlot2D__AddDataSet(plot, dataset)
        call CDataSet__delete(dataset)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet,C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
        do i = 1,nprojects
            call CDataPoint__new2(real(npatch(i),dp), real(avg(i)+std(i), dp), point)
            call CDataSet__AddDataPoint(dataSet, point)
            call CDataPoint__delete(point)
        enddo
        call CPlot2D__AddDataSet(plot, dataset)
        call CDataSet__delete(dataset)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet,C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
        do i = 1,nprojects
            call CDataPoint__new2(real(npatch(i),dp), real(avg(i)-std(i), dp), point)
            call CDataSet__AddDataPoint(dataSet, point)
            call CDataPoint__delete(point)
        enddo
        call CPlot2D__AddDataSet(plot, dataset)
        call CDataSet__delete(dataset)
        call CPlot2D__OutputPostScriptPlot(plot, 'goodnessoffit.eps'//C_NULL_CHAR)
        call CPlot2D__delete(plot)
        ! generate dense 2.5D grid
        call parse_movie_star(moviedocs(1), poly1, ldim1, binning, smpd)
        nframes = ldim1(3)
        nx      = NGRID
        ny      = nint(real(NGRID) * real(ldim1(2)) / real(ldim1(1)))
        if( is_even(nx) ) ny = ny + 1
        npoints = nx*ny*nframes
        allocate(x(nx),y(ny),t(nframes),offsets1(npoints,2),offsets2(npoints,2),source=0.d0)
        write(logfhandle,'(A,3I6)')'>>> GRID DIMENSIONS   :', nx, ny, nframes
        do i = 1,nx
            x(i) = -0.5d0 + real(i,dp)/real(nx+1,dp)
        enddo
        do i = 1,ny
            y(i) = -0.5d0 + real(i,dp)/real(ny+1,dp)
        enddo
        do i = 1,nframes
            t(i) = real(i-1,dp)
        enddo
        ! parsing and calculating
        if(allocated(mask))deallocate(mask)
        if(allocated(avg))deallocate(avg)
        if(allocated(std))deallocate(std)
        allocate(mask(nmics),rmsds(nmics),avg(nprojects),std(nprojects),rmsds2(nmics),avg2(nprojects),std2(nprojects))
        avg  = 0.0
        std  = 0.0
        avg2 = 0.0
        std2 = 0.0
        do i = 1,nprojects-1
            call spproj%read_segment(params%oritype, projects_fnames(order(i)))
            call spproj2%read_segment(params%oritype, projects_fnames(order(i+1)))
            mask   = .false.
            rmsds  = 0.0
            rmsds2 = 0.0
            do k = 1,nmics
                if( indextab(order(i),k) == 0 .or. indextab(order(i+1),k) == 0) cycle
                mask(k) = .true.
                moviedoc = spproj%os_mic%get_str(indextab(order(i),k),'mc_starfile')
                call parse_movie_star(moviedoc, poly1, ldim1, binning, smpd)
                moviedoc = spproj2%os_mic%get_str(indextab(order(i+1),k),'mc_starfile')
                call parse_movie_star(moviedoc, poly2, ldim2, binning, smpd)
                call polynomial2shifts(poly1, x,y,t, offsets1)
                call polynomial2shifts(poly2, x,y,t, offsets2)
                rmsds(k)  = sqrt(sum((offsets2-offsets1)**2.0) / real(npoints,dp))
                rmsds2(k) = sqrt(sum(offsets1**2.0) / real(npoints,dp))
            enddo
            cnt = count(mask)
            if( cnt > 0 )then
                avg(i)  = sum(rmsds,mask=mask) / real(cnt)
                std(i)  = sqrt(sum((rmsds-avg(i))**2,mask=mask)/real(cnt))
                avg2(i) = sum(rmsds2,mask=mask) / real(cnt)
                std2(i) = sqrt(sum((rmsds2-avg2(i))**2,mask=mask)/real(cnt))
            endif
        enddo
        i = nprojects
        call spproj%read_segment(params%oritype, projects_fnames(order(i)))
        mask   = .false.
        rmsds2 = 0.0
        do k = 1,nmics
            if( indextab(order(i),k) == 0) cycle
            mask(k) = .true.
            moviedoc = spproj%os_mic%get_str(indextab(order(i),k),'mc_starfile')
            call parse_movie_star(moviedoc, poly1, ldim1, binning, smpd)
            call polynomial2shifts(poly1, x,y,t, offsets1)
            rmsds2(k) = sqrt(sum(offsets1**2.0) / real(npoints,dp))
        enddo
        cnt = count(mask)
        if( cnt > 0 )then
            avg2(i)  = sum(rmsds2,mask=mask) / real(cnt)
            std2(i)  = sqrt(sum((rmsds2-avg2(i))**2,mask=mask)/real(cnt))
        endif
        ! absolute rmsds
        call fopen(funit, string('absolute_rmsd.txt'), action='READWRITE', status='UNKNOWN', iostat=io_stat)
        call fileiochk('could not write incr_rmsd.txt', io_stat)
        write(funit,'(A)')'sp_project,nx,ny,npatches,avg,std'
        do i = 1,nprojects
            write(funit,'(A,A1,I3,A1,I3,A1,F6.1,A1,F8.3,A1,F8.3)') projects_fnames(order(i))%to_char(),&
                &',',nxpatch(order(i)),',',nypatch(order(i)),',',npatch(i),',',avg2(i),',',std2(i)
        enddo
        call fclose(funit)
        call write_plot('Absolute rmsd', 'n = sqrt(nx * ny)', 'shifts rmsd', nprojects, npatch, avg2, std2, 'absolute_rmsd.eps')
        ! incrementals rmsds
        call fopen(funit, string('incremental_rmsd.txt'), action='READWRITE', status='UNKNOWN', iostat=io_stat)
        call fileiochk('could not write incr_rmsd.txt', io_stat)
        write(funit,'(A)')'index,descr,avg,std'
        do i = 1,nprojects-1
            write(funit,'(I3,A1,A,A1,F8.3,A1,F8.3)')i,',',&
                &int2str(nxpatch(order(i)))//'x'//int2str(nypatch(order(i)))//'_vs_'//int2str(nxpatch(order(i+1)))//'x'//int2str(nypatch(order(i+1))),&
                &',',avg(i),',',std(i)
        enddo
        call fclose(funit)
        order = (/(i,i=1,nprojects)/)
        i = nprojects-1
        call write_plot('Incremental rmsd', 'Order', 'incremtal shifts rmsd', i, real(order(:i)), avg(:i), std(:i), 'incremental_rmsd.eps')
        call simple_end('**** SIMPLE_COMPAREMC NORMAL STOP ****')

        contains

            subroutine write_plot(title, abscissa, ordinate, n, x, y, ystd, fname)
                character(len=*), intent(in) :: title,abscissa,ordinate,fname
                integer,          intent(in) :: n
                real,             intent(in) :: x(n), y(n), ystd(n)
                integer :: i
                call CPlot2D__new(plot, trim(title)//C_NULL_CHAR)
                call CPlot2D__SetXAxisSize(plot, real(50*n,dp))
                call CPlot2D__SetYAxisSize(plot, real(250*maxval(y+2.0*ystd),dp))
                call CPlot2D__SetXAxisTitle(plot, abscissa//c_null_char)
                call CPlot2D__SetYAxisTitle(plot, ordinate//c_null_char)
                call CPlot2D__SetDrawLegend(plot, C_FALSE)
                call CDataSet__new(dataSet)
                call CDataSet__SetDrawMarker(dataSet,C_TRUE)
                call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 0.0_c_double)
                call CDataPoint__new2(real(0.0,dp), real(0.0, dp), point)
                call CDataSet__AddDataPoint(dataSet, point)
                call CDataPoint__delete(point)
                do i = 1,n
                    call CDataPoint__new2(real(x(i),dp), real(y(i), dp), point)
                    call CDataSet__AddDataPoint(dataSet, point)
                    call CDataPoint__delete(point)
                enddo
                call CPlot2D__AddDataSet(plot, dataset)
                call CDataSet__delete(dataset)
                call CDataSet__new(dataSet)
                call CDataSet__SetDrawMarker(dataSet,C_FALSE)
                call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
                do i = 1,n
                    call CDataPoint__new2(real(x(i),dp), real(y(i)+ystd(i), dp), point)
                    call CDataSet__AddDataPoint(dataSet, point)
                    call CDataPoint__delete(point)
                enddo
                call CPlot2D__AddDataSet(plot, dataset)
                call CDataSet__delete(dataset)
                call CDataSet__new(dataSet)
                call CDataSet__SetDrawMarker(dataSet,C_FALSE)
                call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
                do i = 1,n
                    call CDataPoint__new2(real(x(i),dp), real(y(i)-ystd(i), dp), point)
                    call CDataSet__AddDataPoint(dataSet, point)
                    call CDataPoint__delete(point)
                enddo
                call CPlot2D__AddDataSet(plot, dataset)
                call CDataSet__delete(dataset)
                call CPlot2D__OutputPostScriptPlot(plot, trim(fname)//C_NULL_CHAR)
                call CPlot2D__delete(plot)
            end subroutine write_plot

            subroutine polynomial2shifts(cs, xs, ys, times, offsets)
                real(dp), intent(in)  :: cs(POLYDIM,2), xs(nx), ys(ny), times(nframes)
                real(dp), intent(out) :: offsets(npoints,2)
                integer :: i,j,k,l
                l = 0
                do i = 1,nx
                    do j = 1,ny
                        do k = 1,nframes
                            l = l + 1
                            offsets(l,1) = sample_poly(cs(:,1), xs(i), ys(j), times(k))
                            offsets(l,2) = sample_poly(cs(:,2), xs(i), ys(j), times(k))
                        enddo
                    enddo
                enddo
            end subroutine polynomial2shifts

            real(dp) function sample_poly(c, x, y, t)
                real(dp), intent(in) :: c(POLYDIM), x, y, t
                real(dp) :: x2, y2, xy, t2, t3
                x2 = x * x
                y2 = y * y
                xy = x * y
                t2 = t * t
                t3 = t2 * t
                sample_poly =               c( 1) * t      + c( 2) * t2      + c( 3) * t3
                sample_poly = sample_poly + c( 4) * t * x  + c( 5) * t2 * x  + c( 6) * t3 * x
                sample_poly = sample_poly + c( 7) * t * x2 + c( 8) * t2 * x2 + c( 9) * t3 * x2
                sample_poly = sample_poly + c(10) * t * y  + c(11) * t2 * y  + c(12) * t3 * y
                sample_poly = sample_poly + c(13) * t * y2 + c(14) * t2 * y2 + c(15) * t3 * y2
                sample_poly = sample_poly + c(16) * t * xy + c(17) * t2 * xy + c(18) * t3 * xy
            end function sample_poly

            integer function parse_int( table, emdl_id )
                class(starfile_table_type) :: table
                integer, intent(in)        :: emdl_id
                logical :: l_ok
                l_ok = starfile_table__getValue_int(table, emdl_id, parse_int)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end function parse_int

            real(dp) function parse_double( table, emdl_id )
                class(starfile_table_type) :: table
                integer, intent(in)        :: emdl_id
                logical  :: l_ok
                l_ok = starfile_table__getValue_double(table, emdl_id, parse_double)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end function parse_double

            subroutine parse_string( table, emdl_id, str )
                class(starfile_table_type)                 :: table
                integer,                       intent(in)  :: emdl_id
                character(len=:), allocatable, intent(out) :: str
                logical :: l_ok
                l_ok = starfile_table__getValue_string(table, emdl_id, str)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end subroutine parse_string

            subroutine parse_movie_star( fname, polynomial, ldim, binning, smpd )
                class(string),         intent(in)  :: fname
                real(dp), allocatable, intent(out) :: polynomial(:,:)
                integer,               intent(out) :: ldim(3)
                real,                  intent(out) :: binning, smpd
                type(string),    allocatable :: names(:)
                type(starfile_table_type)    :: table
                integer          :: i,n,ind,start_frame,motion_model
                integer(C_long)  :: num_objs, object_id
                character(len=:), allocatable :: buffer
                call starfile_table__new(table)
                call starfile_table__getnames(table, fname, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(table, fname, names(i)%to_char() )
                    select case(names(i)%to_char())
                        case('general')
                            ldim(1) = parse_int(table, EMDL_IMAGE_SIZE_X)
                            ldim(2) = parse_int(table, EMDL_IMAGE_SIZE_Y)
                            ldim(3) = parse_int(table, EMDL_IMAGE_SIZE_Z)
                            call parse_string(table, EMDL_MICROGRAPH_MOVIE_NAME, buffer)
                            movie        = buffer
                            binning      = real(parse_double(table, EMDL_MICROGRAPH_BINNING))
                            smpd         = real(parse_double(table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
                            start_frame  = parse_int(table, EMDL_MICROGRAPH_START_FRAME)
                            motion_model = parse_int(table, EMDL_MICROGRAPH_MOTION_MODEL_VERSION)
                            if( motion_model /= 1 )then
                                THROW_HARD('No polynomial found in: '//fname%to_char())
                            else
                                allocate(polynomial(POLYDIM,2),source=0.d0)
                            endif
                            if(allocated(buffer)) deallocate(buffer)
                        case('local_motion_model')
                            object_id = starfile_table__firstobject(table) ! base 0
                            num_objs  = starfile_table__numberofobjects(table)
                            do while( (object_id < num_objs) .and. (object_id >= 0) )
                                ind = parse_int(table, EMDL_MICROGRAPH_MOTION_COEFFS_IDX) + 1
                                if( ind > POLYDIM )then
                                    polynomial(ind-POLYDIM,2) = parse_double(table, EMDL_MICROGRAPH_MOTION_COEFF)
                                else
                                    polynomial(ind,        1) = parse_double(table, EMDL_MICROGRAPH_MOTION_COEFF)
                                endif
                                object_id = starfile_table__nextobject(table)
                            end do
                        case DEFAULT
                            cycle
                    end select
                enddo
                call starfile_table__delete(table)
            end subroutine parse_movie_star

    end subroutine exec_comparemc

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

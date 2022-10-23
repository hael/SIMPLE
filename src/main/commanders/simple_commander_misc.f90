! concrete commander: miscallenaous routines
module simple_commander_misc
include 'simple_lib.f08'
include "starfile/starfile_enum.inc"
use simple_binoris_io
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_projector_hlev, only: rotvol
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_builder,        only: builder
use simple_parameters,     only: parameters
implicit none

public :: masscen_commander
public :: print_fsc_commander
public :: print_magic_boxes_commander
public :: print_dose_weights_commander
public :: kstest_commander
public :: mkdir_commander
public :: remoc_commander
public :: comparemc_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: masscen_commander
  contains
    procedure :: execute      => exec_masscen
end type masscen_commander

type, extends(commander_base) :: print_fsc_commander
  contains
    procedure :: execute       => exec_print_fsc
end type print_fsc_commander

type, extends(commander_base) :: print_magic_boxes_commander
  contains
    procedure :: execute       => exec_print_magic_boxes
end type print_magic_boxes_commander

type, extends(commander_base) :: print_dose_weights_commander
  contains
    procedure :: execute       => exec_print_dose_weights
end type print_dose_weights_commander

type, extends(commander_base) :: kstest_commander
  contains
    procedure :: execute       => exec_kstest
end type kstest_commander

type, extends(commander_base) :: mkdir_commander
  contains
    procedure :: execute       => exec_mkdir
end type mkdir_commander

type, extends(commander_base) :: remoc_commander
  contains
    procedure :: execute       => exec_remoc
end type remoc_commander

type, extends(commander_base) :: comparemc_commander
  contains
    procedure :: execute       => exec_comparemc
end type comparemc_commander

contains

    !> centers base on centre of mass
     subroutine exec_masscen( self, cline )
        use simple_procimgstk, only: masscen_imgfile
        class(masscen_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        params%cenlp = params%lp
        ! center of mass centering
        call masscen_imgfile(params%stk, params%outstk, params%smpd, params%lp, params%msk)
        ! end gracefully
        call simple_end('**** SIMPLE_MASSCEN NORMAL STOP ****')
    end subroutine exec_masscen

    !>  for printing the binary FSC files produced by PRIME3D
    subroutine exec_print_fsc( self, cline )
        use simple_fsc, only: plot_fsc
        class(print_fsc_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)      :: params
        type(image)           :: img
        character(len=STDLEN) :: tmpl_fname
        real,     allocatable :: res(:), fsc(:)
        integer               :: k,n
        real                  :: res0143, res05
        call params%new(cline)
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
        tmpl_fname = get_fbody(params%fsc,trim(BIN_EXT),separator=.false.)
        call plot_fsc(n, fsc, res, params%smpd, tmpl_fname)
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_FSC NORMAL STOP ****')
    end subroutine exec_print_fsc

    !> for printing magic box sizes (fast FFT)
    subroutine exec_print_magic_boxes( self, cline )
        class(print_magic_boxes_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call print_magic_box_range(params%smpd, params%moldiam )
        ! end gracefully
        call simple_end('**** SIMPLE_PRINT_MAGIC_BOXES NORMAL STOP ****')
    end subroutine exec_print_magic_boxes

    subroutine exec_print_dose_weights( self, cline )
        class(print_dose_weights_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)  :: params
        real, allocatable :: weights(:,:), res(:)
        integer           :: iframe, k, filtsz
        call params%new(cline)
        call calc_dose_weights(params%nframes, params%box, params%smpd, params%kV, params%exp_time, params%dose_rate, weights)
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
        class(kstest_commander), intent(inout) :: self
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

        contains

            subroutine read_nrs_dat( filename, arr, ndat )
                character(len=*),  intent(in)  :: filename
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

    end subroutine exec_kstest

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
        real     :: cen1, cen2, sum1, sum2, sumvals, sdev, sdev_noise
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
            call read_img%noise_norm(l_msk, sdev_noise)
            ! center image
            sh1 = read_img%calc_shiftcen(params_glob%cenlp, params_glob%msk)
            call read_img%shift(-sh1)
            call dsym_os%set_shift(i, -sh1)
            ! rotational image
            call read_img%roavg(nint(ang), roavg_img)
            call roavg_img%write('roavg.mrc',i)
            corrs(i) = read_img%real_corr(roavg_img, l_msk)
            ! radii
            call read_img%bp(0., params_glob%cenlp)
            call read_img%mask(params_glob%msk, 'hard')
            call otsu_robust_fast(read_img, is2d=.true., noneg=.false., thresh=thresh)
            call read_img%write('bin.mrc',i)
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
                call read_img%noise_norm(l_msk, sdev_noise)
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
        class(mkdir_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        call cline%set('mkdir', 'yes')
        call params%new(cline)
    end subroutine exec_mkdir

    subroutine exec_remoc( self, cline )
        use simple_starfile_wrappers
        use simple_projector
        class(remoc_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        class(str4arr),   allocatable :: mics(:),mics_from_ptcls(:), star_movies(:)
        character(len=:), allocatable :: fname, movie, vol_even_fname, vol_odd_fname, mask_fname, mic_fname
        real,             allocatable :: fsc(:), optfilter(:)
        type(parameters) :: params
        type(ctfparams)  :: ctf_glob, ctf_mic
        type(oris)       :: os
        type(projector)  :: vol_even, vol_odd
        type(image)      :: vol, vol_mask
        real    :: mov_smpd
        integer :: ldim(3), ldim_mov(3), nmics, nmovies, filtsz
        logical :: l_scale
        call cline%set('mkdir',       'yes')
        if( .not. cline%defined('trs')           ) call cline%set('trs',           20.)
        if( .not. cline%defined('lpstop')        ) call cline%set('lpstop',         5.)
        if( .not. cline%defined('bfac')          ) call cline%set('bfac',          50.)
        if( .not. cline%defined('wcrit')         ) call cline%set('wcrit',   'softmax')
        if( .not. cline%defined('eer_fraction')  ) call cline%set('eer_fraction',  20.)
        if( .not. cline%defined('eer_upsampling')) call cline%set('eer_upsampling', 1.)
        call cline%set('oritype', 'mic')
        call params%new(cline)
        ! Particles
        call parse_particles( params%star_ptcl, os, mics, ctf_glob, params%box, params%nptcls )
        params%smpd = ctf_glob%smpd
        ldim = params%box
        ! Micrographs
        call parse_mics( params%star_mic, mics, star_movies, ctf_mic, mov_smpd )
        params%smpd = ctf_glob%smpd
        nmics   = size(mics)
        nmovies = size(star_movies)
        if( abs(params%smpd-ctf_mic%smpd) > 0.01 ) THROW_HARD('Inconsistent pixel sizes between micrographs and particles!')
        params%scale = mov_smpd / ctf_glob%smpd
        l_scale = abs(params%scale-1.) > 0.01
        if( l_scale )then
            params%alpha = 1./params%scale
            ldim_mov = round2even(real(ldim)/params%scale)
        else
            params%alpha = 1.
            ldim_mov     = ldim
        endif
        ! Model
        call parse_model(params%star_model, vol_even_fname, vol_odd_fname, mask_fname, fsc)
        mask_fname = get_fpath(params%star_model)//'/'//trim(mask_fname)
        filtsz = size(fsc)
        allocate(optfilter(filtsz),source=0.)
        call fsc2optlp_sub(filtsz, fsc, optfilter)
        ! volumes prep
        call vol%new(ldim,params%smpd)
        call vol_mask%new(ldim,params%smpd)
        call vol%read(vol_even_fname)
        call vol_mask%read(mask_fname)
        call vol%mul(vol_mask)
        if( l_scale )then
            call vol_even%new(ldim_mov,mov_smpd)
            call vol%pad(vol_even)
        else
            vol_even = vol
        endif
        call vol%read(vol_odd_fname)
        call vol%mul(vol_mask)
        call vol_mask%kill
        if( l_scale )then
            call vol_odd%new(ldim_mov,mov_smpd)
            call vol%pad(vol_odd)
        else
            vol_odd = vol
        endif
        call vol%kill
        call vol_even%fft
        call vol_odd%fft
        call vol_even%apply_filter(optfilter)
        call vol_odd%apply_filter(optfilter)
        call vol_even%expand_cmat(1.,norm4proj=.true.)
        call vol_odd%expand_cmat(1.,norm4proj=.true.)
        ! end gracefully
        call simple_end('**** SIMPLE_REMOC NORMAL STOP ****')
        contains

            integer function parse_int( table, emdl_id )
                class(starfile_table_type) :: table
                integer, intent(in)        :: emdl_id
                logical :: l_ok
                l_ok = starfile_table__getValue_int(table, emdl_id, parse_int)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end function parse_int

            real function parse_double( table, emdl_id )
                class(starfile_table_type) :: table
                integer, intent(in)        :: emdl_id
                real(dp) :: v
                logical  :: l_ok
                l_ok = starfile_table__getValue_double(table, emdl_id, v)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
                parse_double = real(v)
            end function parse_double

            subroutine parse_string( table, emdl_id, string )
                class(starfile_table_type)                 :: table
                integer,                       intent(in)  :: emdl_id
                character(len=:), allocatable, intent(out) :: string
                logical :: l_ok
                l_ok = starfile_table__getValue_string(table, emdl_id, string)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end subroutine parse_string

            subroutine parse_particles( fname, os, mics, ctf, box, nptcls_out )
                character(len=*),               intent(in)  :: fname
                class(oris),                    intent(out) :: os
                class(str4arr),    allocatable, intent(out) :: mics(:)
                class(ctfparams),               intent(out) :: ctf
                integer,                        intent(out) :: box, nptcls_out
                type(str4arr),    allocatable :: names(:)
                type(str4arr)    ::fname_here
                type(starfile_table_type)     :: star_table
                type(ctfparams)  :: ctfparms
                real             :: euls(3),pos(2),shift(2)
                integer          :: i,n,ind,eo
                integer(C_long)  :: num_objs, object_id
                fname_here%str = trim(fname)//C_NULL_CHAR
                call starfile_table__new(star_table)
                call starfile_table__getnames(star_table, fname_here%str, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(star_table, fname_here%str, names(i)%str )
                    select case(trim(names(i)%str))
                    case('optics')
                        ! ignoring multiple optics groups for now
                        object_id = starfile_table__firstobject(star_table) ! base 0
                        ! ctf
                        ctf%fraca = parse_double(star_table, EMDL_CTF_Q0)
                        ctf%cs    = parse_double(star_table, EMDL_CTF_CS)
                        ctf%kv    = parse_double(star_table, EMDL_CTF_VOLTAGE)
                        ctf%smpd  = parse_double(star_table, EMDL_IMAGE_PIXEL_SIZE)
                        ! dimensions
                        box = parse_int(star_table, EMDL_IMAGE_SIZE)
                    case('particles')
                        object_id  = starfile_table__firstobject(star_table) ! base 0
                        num_objs   = starfile_table__numberofobjects(star_table)
                        nptcls_out = int(num_objs - object_id)
                        call os%new(nptcls_out, is_ptcl=.true.)
                        allocate(mics(nptcls_out))
                        ind = 0
                        do while( (object_id < num_objs) .and. (object_id >= 0) )
                            ind = ind+1
                            ! Micrograph name
                            call parse_string(star_table, EMDL_MICROGRAPH_NAME, mics(ind)%str)
                            ! picking coordinates
                            pos(1) = parse_double(star_table, EMDL_IMAGE_COORD_X) - box/2 - 1
                            pos(2) = parse_double(star_table, EMDL_IMAGE_COORD_Y) - box/2 - 1
                            ! CTF
                            ctfparms%dfx    = parse_double(star_table, EMDL_CTF_DEFOCUSU) / 10000.0 ! Microns
                            ctfparms%dfy    = parse_double(star_table, EMDL_CTF_DEFOCUSU) / 10000.0 ! microns
                            ctfparms%angast = parse_double(star_table, EMDL_CTF_DEFOCUS_ANGLE)
                            ! e/o
                            eo = parse_int(star_table, EMDL_PARTICLE_RANDOM_SUBSET) - 1
                            ! Alignment parameters
                            euls(1)  = parse_double(star_table, EMDL_ORIENT_ROT)
                            euls(2)  = parse_double(star_table, EMDL_ORIENT_TILT)
                            euls(3)  = parse_double(star_table, EMDL_ORIENT_PSI)
                            shift(1) = parse_double(star_table, EMDL_ORIENT_ORIGIN_X_ANGSTROM) / ctf_glob%smpd ! in pixels
                            shift(2) = parse_double(star_table, EMDL_ORIENT_ORIGIN_X_ANGSTROM) / ctf_glob%smpd ! in pixels
                            ! transfer to oris object
                            call os%set_euler(ind,   euls)
                            call os%set_shift(ind,   shift)
                            call os%set_dfx(ind,     ctfparms%dfx)
                            call os%set_dfy(ind,     ctfparms%dfy)
                            call os%set(ind,'angast',ctfparms%angast)
                            call os%set(ind,'xpos',  pos(1))
                            call os%set(ind,'ypos',  pos(2))
                            call os%set(ind,'eo',    real(eo))
                            ! done
                            object_id = starfile_table__nextobject(star_table)
                        end do
                    case DEFAULT
                        THROW_HARD('Invalid table: '//trim(names(i)%str))
                    end select
                enddo
                call starfile_table__delete(star_table)
            end subroutine parse_particles

            subroutine parse_mics( fname, mics, star_movies, ctf, mov_smpd )
                character(len=*),                intent(in) :: fname
                class(str4arr),    allocatable, intent(out) :: mics(:), star_movies(:)
                class(ctfparams),               intent(out) :: ctf
                real,                           intent(out) :: mov_smpd
                type(str4arr),    allocatable :: names(:)
                type(starfile_table_type)     :: star_table
                integer          :: i,n,nmics,ind
                integer(C_long)  :: num_objs, object_id
                call starfile_table__new(star_table)
                call starfile_table__getnames(star_table, trim(fname)//C_NULL_CHAR, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(star_table, trim(fname)//C_NULL_CHAR, names(i)%str )
                    select case(trim(names(i)%str))
                    case('optics')
                        ! ignoring multiple optics groups for now
                        object_id = starfile_table__firstobject(star_table) ! base 0
                        ! ctf
                        ctf%fraca = parse_double(star_table, EMDL_CTF_Q0)
                        ctf%cs    = parse_double(star_table, EMDL_CTF_CS)
                        ctf%kv    = parse_double(star_table, EMDL_CTF_VOLTAGE)
                        ctf%smpd  = parse_double(star_table, EMDL_MICROGRAPH_PIXEL_SIZE)
                        mov_smpd  = parse_double(star_table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE)
                    case('micrographs')
                        object_id = starfile_table__firstobject(star_table) ! base 0
                        num_objs  = starfile_table__numberofobjects(star_table)
                        nmics     = int(num_objs - object_id)
                        allocate(mics(nmics),star_movies(nmics))
                        ind = 0
                        do while( (object_id < num_objs) .and. (object_id >= 0) )
                            ind = ind+1
                            call parse_string(star_table, EMDL_MICROGRAPH_NAME, mics(ind)%str)
                            call parse_string(star_table, EMDL_MICROGRAPH_METADATA_NAME, star_movies(ind)%str)
                            ! done
                            object_id = starfile_table__nextobject(star_table)
                        end do
                    case DEFAULT
                        THROW_HARD('Invalid table: '//trim(names(i)%str))
                    end select
                enddo
                call starfile_table__delete(star_table)
            end subroutine parse_mics

            subroutine parse_model( fname, vol_even, vol_odd, mask, fsc )
                character(len=*),              intent(in)  :: fname
                character(len=:), allocatable, intent(out) :: vol_even, vol_odd, mask
                real,             allocatable, intent(out) :: fsc(:)
                type(str4arr),    allocatable :: names(:)
                type(starfile_table_type)     :: star_table
                integer          :: i,n,ind,nn
                integer(C_long)  :: num_objs, object_id
                call starfile_table__new(star_table)
                call starfile_table__getnames(star_table, trim(fname)//C_NULL_CHAR, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(star_table, trim(fname)//C_NULL_CHAR, names(i)%str )
                    select case(trim(names(i)%str))
                    case('general')
                        object_id = starfile_table__firstobject(star_table)
                        ! e/o volumes, mask
                        call parse_string(star_table, EMDL_POSTPROCESS_UNFIL_HALFMAP1, vol_even)
                        call parse_string(star_table, EMDL_POSTPROCESS_UNFIL_HALFMAP2, vol_odd)
                        call parse_string(star_table, EMDL_MASK_NAME, mask)
                    case('fsc')
                        object_id = starfile_table__firstobject(star_table)
                        num_objs  = starfile_table__numberofobjects(star_table)
                        nn        = int(num_objs - object_id)
                        if(nn /= params%box/2+1 ) THROW_HARD('Inconsistent image & fsc dimensions!')
                        allocate(fsc(params%box/2),source=0.)
                        do while( (object_id < num_objs) .and. (object_id >= 0) )
                            ind = parse_int(star_table, EMDL_SPECTRAL_IDX)
                            if( ind /= 0 ) fsc(ind) = parse_double(star_table, EMDL_POSTPROCESS_FSC_TRUE)
                            ! done
                            object_id = starfile_table__nextobject(star_table)
                        end do
                    case('guinier')
                        ! Ignored for now
                    case DEFAULT
                        THROW_HARD('Invalid table: '//trim(names(i)%str))
                    end select
                enddo
                call starfile_table__delete(star_table)
            end subroutine parse_model

    end subroutine exec_remoc

    subroutine exec_comparemc( self, cline )
        use simple_starfile_wrappers
        use CPlot2D_wrapper_module
        integer, PARAMETER :: POLYDIM = 18
        integer, PARAMETER :: NGRID   = 25 ! has to be odd
        class(comparemc_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        real,                         allocatable :: rmsds(:), gofs(:,:), avg(:), std(:), npatch(:), rmsds2(:), avg2(:), std2(:)
        real(dp),                     allocatable :: poly1(:,:), poly2(:,:), x(:), y(:), t(:), offsets1(:,:), offsets2(:,:)
        character(len=:),             allocatable :: movie, moviedoc
        character(len=LONGSTRLEN),    allocatable :: projects_fnames(:), moviedocs(:)
        integer,                      allocatable :: indextab(:,:), order(:), nxpatch(:), nypatch(:)
        logical,                      allocatable :: mask(:)
        type(CPlot2D_type)    :: plot2D
        type(CDataSet_type)   :: dataSet
        type(CDataPoint_type) :: point
        type(parameters)      :: params
        type(sp_project)      :: first_spproj, spproj, spproj2
        real(dp) :: rmsd, maxrmsd, minrmsd
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
            if( projects_fnames(i)(1:1).ne.'/' )then
                if( trim(params%mkdir).eq.'yes')then
                    projects_fnames(i) = '../'//trim(projects_fnames(i))
                endif
                projects_fnames(i) = simple_abspath(projects_fnames(i), errmsg='simple_fileio::make_relativepath: '//trim(projects_fnames(i)))
            endif
            if( .not.file_exists(projects_fnames(i))) THROW_HARD('Could not find: '//trim(projects_fnames(i)))
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
            moviedocs(i)  = first_spproj%os_mic%get_static(i,'mc_starfile')
            gofs(1,i)     = sqrt(first_spproj%os_mic%get(i,'gofx')**2.0+first_spproj%os_mic%get(i,'gofy')**2.0)
        enddo
        do i = 2,nprojects
            call spproj%read_segment(params%oritype, projects_fnames(i))
            do j = 1,nmics
                found = .false.
                do k = 1,nmics
                    moviedoc = basename(spproj%os_mic%get_static(k,'mc_starfile'))
                    if( trim(moviedoc).eq.trim(basename(moviedocs(j))) )then
                        indextab(i,j) = k
                        gofs(i,j) = sqrt(spproj%os_mic%get(k,'gofx')**2.0+spproj%os_mic%get(k,'gofy')**2.0)
                        found = .true.
                    endif
                    if( found ) exit
                enddo
                if( .not.found )then
                    indextab(:,j) = 0
                    THROW_WARN('Could not find corresponding doc: '//trim(moviedoc))
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
        call fopen(funit, 'goodnessoffit.txt', action='READWRITE', status='UNKNOWN', iostat=io_stat)
        call fileiochk('could not write goodnessoffit.txt', io_stat)
        write(funit,'(A)')'sp_project,nx,ny,npatches,gof_avg,gof_std'
        do i = 1,nprojects
            write(funit,'(A,A1,I3,A1,I3,A1,F6.1,A1,F8.3,A1,F8.3)')trim(projects_fnames(order(i))),&
                &',',nxpatch(order(i)),',',nypatch(order(i)),',',npatch(i),',',avg(i),',',std(i)
        enddo
        call fclose(funit)
        ! eps output
        call CPlot2D__new(plot2D, 'Goodness of  fit'//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, real(50*nprojects,dp))
        call CPlot2D__SetYAxisSize(plot2D, real(50*maxval(avg+2.0*std),dp))
        call CPlot2D__SetXAxisTitle(plot2D, 'n = sqrt(nx * ny)'//c_null_char)
        call CPlot2D__SetYAxisTitle(plot2D, 'shifts RMS patch vs poly (pixels) '//c_null_char)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet,C_TRUE)
        call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 0.0_c_double)
        do i = 1,nprojects
            call CDataPoint__new2(real(npatch(i),dp), real(avg(i), dp), point)
            call CDataSet__AddDataPoint(dataSet, point)
            call CDataPoint__delete(point)
        enddo
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet,C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
        do i = 1,nprojects
            call CDataPoint__new2(real(npatch(i),dp), real(avg(i)+std(i), dp), point)
            call CDataSet__AddDataPoint(dataSet, point)
            call CDataPoint__delete(point)
        enddo
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet,C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
        do i = 1,nprojects
            call CDataPoint__new2(real(npatch(i),dp), real(avg(i)-std(i), dp), point)
            call CDataSet__AddDataPoint(dataSet, point)
            call CDataPoint__delete(point)
        enddo
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        call CPlot2D__OutputPostScriptPlot(plot2D, 'goodnessoffit.eps'//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
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
                moviedoc = trim(spproj%os_mic%get_static(indextab(order(i),k),'mc_starfile'))
                call parse_movie_star(moviedoc, poly1, ldim1, binning, smpd)
                moviedoc = trim(spproj2%os_mic%get_static(indextab(order(i+1),k),'mc_starfile'))
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
            moviedoc = trim(spproj%os_mic%get_static(indextab(order(i),k),'mc_starfile'))
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
        call fopen(funit, 'absolute_rmsd.txt', action='READWRITE', status='UNKNOWN', iostat=io_stat)
        call fileiochk('could not write incr_rmsd.txt', io_stat)
        write(funit,'(A)')'sp_project,nx,ny,npatches,avg,std'
        do i = 1,nprojects
            write(funit,'(A,A1,I3,A1,I3,A1,F6.1,A1,F8.3,A1,F8.3)')trim(projects_fnames(order(i))),&
                &',',nxpatch(order(i)),',',nypatch(order(i)),',',npatch(i),',',avg2(i),',',std2(i)
        enddo
        call fclose(funit)
        call write_plot('Absolute rmsd', 'n = sqrt(nx * ny)', 'shifts rmsd', nprojects, npatch, avg2, std2, 'absolute_rmsd.eps')
        ! incrementals rmsds
        call fopen(funit, 'incremental_rmsd.txt', action='READWRITE', status='UNKNOWN', iostat=io_stat)
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
                call CPlot2D__new(plot2D, trim(title)//C_NULL_CHAR)
                call CPlot2D__SetXAxisSize(plot2D, real(50*n,dp))
                call CPlot2D__SetYAxisSize(plot2D, real(250*maxval(y+2.0*ystd),dp))
                call CPlot2D__SetXAxisTitle(plot2D, abscissa//c_null_char)
                call CPlot2D__SetYAxisTitle(plot2D, ordinate//c_null_char)
                call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
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
                call CPlot2D__AddDataSet(plot2D, dataset)
                call CDataSet__delete(dataset)
                call CDataSet__new(dataSet)
                call CDataSet__SetDrawMarker(dataSet,C_FALSE)
                call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
                do i = 1,n
                    call CDataPoint__new2(real(x(i),dp), real(y(i)+ystd(i), dp), point)
                    call CDataSet__AddDataPoint(dataSet, point)
                    call CDataPoint__delete(point)
                enddo
                call CPlot2D__AddDataSet(plot2D, dataset)
                call CDataSet__delete(dataset)
                call CDataSet__new(dataSet)
                call CDataSet__SetDrawMarker(dataSet,C_FALSE)
                call CDataSet__SetDatasetColor(dataSet, 0.4_c_double, 0.4_c_double, 0.4_c_double)
                do i = 1,n
                    call CDataPoint__new2(real(x(i),dp), real(y(i)-ystd(i), dp), point)
                    call CDataSet__AddDataPoint(dataSet, point)
                    call CDataPoint__delete(point)
                enddo
                call CPlot2D__AddDataSet(plot2D, dataset)
                call CDataSet__delete(dataset)
                call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname)//C_NULL_CHAR)
                call CPlot2D__delete(plot2D)
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

            subroutine parse_string( table, emdl_id, string )
                class(starfile_table_type)                 :: table
                integer,                       intent(in)  :: emdl_id
                character(len=:), allocatable, intent(out) :: string
                logical :: l_ok
                l_ok = starfile_table__getValue_string(table, emdl_id, string)
                if(.not.l_ok) THROW_HARD('Missing value in table!')
            end subroutine parse_string

            subroutine parse_movie_star( fname, polynomial, ldim, binning, smpd )
                character(len=*),      intent(in)  :: fname
                real(dp), allocatable, intent(out) :: polynomial(:,:)
                integer,               intent(out) :: ldim(3)
                real,                  intent(out) :: binning, smpd
                type(str4arr),    allocatable :: names(:)
                type(starfile_table_type)     :: table
                integer          :: i,n,ind,start_frame,motion_model
                integer(C_long)  :: num_objs, object_id
                call starfile_table__new(table)
                call starfile_table__getnames(table, trim(fname)//C_NULL_CHAR, names)
                n = size(names)
                do i =1,n
                    call starfile_table__read(table, trim(fname)//C_NULL_CHAR, names(i)%str )
                    select case(trim(names(i)%str))
                    case('general')
                        ldim(1) = parse_int(table, EMDL_IMAGE_SIZE_X)
                        ldim(2) = parse_int(table, EMDL_IMAGE_SIZE_Y)
                        ldim(3) = parse_int(table, EMDL_IMAGE_SIZE_Z)
                        call parse_string(table, EMDL_MICROGRAPH_MOVIE_NAME, movie)
                        binning      = real(parse_double(table, EMDL_MICROGRAPH_BINNING))
                        smpd         = real(parse_double(table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
                        start_frame  = parse_int(table, EMDL_MICROGRAPH_START_FRAME)
                        motion_model = parse_int(table, EMDL_MICROGRAPH_MOTION_MODEL_VERSION)
                        if( motion_model /= 1 )then
                            THROW_HARD('No polynomial found in: '//trim(fname))
                        else
                            allocate(polynomial(POLYDIM,2),source=0.d0)
                        endif
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

end module simple_commander_misc

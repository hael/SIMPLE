!@descr: experimental commanders under development
module single_commanders_experimental
use simple_commanders_api
use simple_commanders_volops, only: commander_reproject
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cavgsproc_nano
  contains
    procedure :: execute      => exec_cavgsproc_nano
end type commander_cavgsproc_nano

type, extends(commander_base) :: commander_cavgseoproc_nano
  contains
    procedure :: execute      => exec_cavgseoproc_nano
end type commander_cavgseoproc_nano

type, extends(commander_base) :: commander_ptclsproc_nano
  contains
    procedure :: execute      => exec_ptclsproc_nano
end type commander_ptclsproc_nano

type, extends(commander_base) :: commander_tsegmaps_core_finder
  contains
    procedure :: execute      => exec_tsegmaps_core_finder
end type commander_tsegmaps_core_finder

type, extends(commander_base) :: commander_trajectory_make_projavgs
  contains
    procedure :: execute      => exec_trajectory_make_projavgs
end type commander_trajectory_make_projavgs

contains

    subroutine exec_cavgsproc_nano( self, cline )
        use single_commanders_nano3D, only: commander_refine3D_nano
        class(commander_cavgsproc_nano), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)              :: params
        type(commander_refine3D_nano) :: xrefine3D_nano
        type(commander_reproject)     :: xreproject
        type(cmdline)                 :: cline_refine3D_cavgs, cline_reproject
        type(image),      allocatable :: imgs(:)
        type(image)                   :: img_w
        type(sp_project)              :: spproj
        type(string)                  :: cavgs_stk, command_plot, fname
        real,             allocatable :: rstates(:), rad_cc(:,:), rad_dists(:,:)
        logical,          allocatable :: state_mask(:)
        integer,          allocatable :: pinds(:)
        integer :: ncavgs, i, j, cnt, cnt2, tmax, tmin, tstamp
        real    :: smpd
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_cavgsproc_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        ! retrieve cavgs stack
        call spproj%get_cavgs_stk(cavgs_stk, ncavgs, smpd, fail=.false.)
        if( ncavgs /= 0 )then
            ! update cline_refine3D_cavgs accordingly
            call cline_refine3D_cavgs%set('prg',      'refine3D_nano')
            call cline_refine3D_cavgs%set('vol1',      params%vols(1))
            call cline_refine3D_cavgs%set('pgrp',         params%pgrp)
            call cline_refine3D_cavgs%set('mskdiam',   params%mskdiam)
            call cline_refine3D_cavgs%set('nthr',         params%nthr)
            call cline_refine3D_cavgs%set('mkdir',               'no')
            call cline_refine3D_cavgs%set('maxits',                 1)
            call cline_refine3D_cavgs%set('projfile', params%projfile)
            call cline_refine3D_cavgs%set('oritype',          'cls3D')
            call cline_refine3D_cavgs%set('objfun',              'cc')
            call cline_refine3D_cavgs%set('lp',                   1.0)
            call cline_refine3D_cavgs%set('trs',                  5.0)
            call cline_refine3D_cavgs%set('nspace',             10000)
            call cline_refine3D_cavgs%set('center',              'no')
            ! convention for executing shared-memory workflows from within another workflow with a parameters object declared
            call xrefine3D_nano%execute_safe(cline_refine3D_cavgs)
            ! align cavgs
            call spproj%read_segment('cls3D', params%projfile) ! now the newly generated cls3D field will be read...
            ! ...so write out its content
            call spproj%os_cls3D%write(string('cavgs_oris.txt'))
            ! ...and get the state flags
            if( allocated(rstates) ) deallocate(rstates)
            rstates = spproj%os_cls3D%get_all('state')
            ! prepare for re-projection
            call cline_reproject%set('vol1',      params%vols(1))
            call cline_reproject%set('outstk',     'reprojs.mrc')
            call cline_reproject%set('smpd',                smpd)
            call cline_reproject%set('oritab',  'cavgs_oris.txt')
            call cline_reproject%set('pgrp',         params%pgrp)
            call cline_reproject%set('nthr',         params%nthr)
            call xreproject%execute_safe(cline_reproject)
            ! compute radial cross-correlation between cavgs and reproj
            allocate(rad_cc(ncavgs,params%box/2), rad_dists(ncavgs,params%box/2))
            ! write cavgs & reprojections
            allocate(imgs(2*ncavgs), state_mask(ncavgs))
            cnt = 0
            open(unit=25, file="radial_analysis.csv")
            do i = 1,2*ncavgs,2
                cnt = cnt + 1
                if( rstates(cnt) > 0.5 )then
                    call imgs(i    )%new([params%box,params%box,1], smpd)
                    call imgs(i + 1)%new([params%box,params%box,1], smpd)
                    call imgs(i    )%read(cavgs_stk,     cnt)
                    call imgs(i + 1)%read(string('reprojs.mrc'), cnt)
                    call img_w%new([params%box,params%box,1], smpd)
                    ! cavgs images are weighted using radial cross-correlation
                    call spproj%os_ptcl2D%get_pinds(cnt, 'class', pinds)
                    if( .not. allocated(pinds) ) cycle
                    tmax   = maxval(pinds)
                    tmin   = minval(pinds)
                    tstamp = tmin + (tmax-tmin)/2
                    call imgs(i)%radial_cc(imgs(i+1), img_w, smpd, rad_cc(cnt,:), rad_dists(cnt,:)) 
                    do j = 1, size(rad_dists,dim=2)
                        write(25,'(i6, i6, 2F18.6, i6, i6)') cnt, tstamp, rad_dists(cnt,j), rad_cc(cnt,j), tmax-tmin, size(pinds)
                    enddo
                    ! filter out cavgs
                    call imgs(i    )%norm
                    call imgs(i + 1)%norm
                    state_mask(cnt) = .true.
                else
                    state_mask(cnt) = .false.
                endif
                write(25, '(A)') "          "

            enddo
            command_plot = "gnuplot -e "//'"'//"set pm3d map; set zrange[-.4:1]; splot " //"'"//"radial_analysis.csv"//"'"// &
            " u 2:3:4 ; set term png; set xlabel " //"'"//"Time"//"'"// "; set ylabel " //"'"//"Radius({\305})"// &
            "'"// "; set title " //"'"//"Radial Cross-correlation"//"'"// "; set nokey; set output 'radial_analysis.png'; replot" //'"'
            call execute_command_line(command_plot%to_char())
            close(25)
            cnt   = 0
            cnt2  = 1 ! needed because we have omissions
            fname = 'cavgs_vs_reprojections.mrc'
            do i = 1,2*ncavgs,2
                cnt = cnt + 1
                if( state_mask(cnt) )then
                    call imgs(i    )%write(fname, cnt2    )
                    call imgs(i + 1)%write(fname, cnt2 + 1)
                    call imgs(i    )%kill
                    call imgs(i + 1)%kill
                    cnt2 = cnt2 + 2
                endif
            enddo
            deallocate(imgs)
        endif ! end of class average-based validation
        call exec_cmdline('rm -rf fsc* fft* recvol* RES* reprojs_recvol* reprojs* reproject_oris.txt stderrout')
        ! deallocate
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** CAVGSPROC_NANO NORMAL STOP ****')
    end subroutine exec_cavgsproc_nano

    subroutine exec_cavgseoproc_nano( self, cline )
        class(commander_cavgseoproc_nano), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: cavg_even, cavg_odd, img_w
        type(sp_project)     :: spproj
        type(string)         :: cavgs_stk, cavgs_stk_even, cavgs_stk_odd
        type(string)         :: command_plot
        real,    allocatable :: rstates(:), rad_cc(:,:), rad_dists(:,:)
        logical, allocatable :: state_mask(:)
        integer, allocatable :: pinds(:)
        integer :: ncavgs, i, j, cnt, tmax, tmin, tstamp
        real    :: smpd
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_cavgseoproc_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        ! retrieve even and odd cavgs stack
        call spproj%get_cavgs_stk(cavgs_stk, ncavgs, smpd, fail=.false.)
        cavgs_stk_even = add2fbody(cavgs_stk,'.mrc','_even')
        cavgs_stk_odd  = add2fbody(cavgs_stk,'.mrc','_odd')
        if( allocated(rstates) ) deallocate(rstates)
        rstates = spproj%os_cls3D%get_all('state')
        ! compute radial cross-correlation between the even and odd cavgs
        print *,"params", params%box
        allocate(rad_cc(ncavgs,params%box/2), rad_dists(ncavgs,params%box/2), state_mask(ncavgs))
        cnt = 0
        open(unit=25, file="evenodd_radial_analysis.csv")
        print *, "ncavgs", ncavgs
        do i = 1, ncavgs
            cnt = cnt + 1
            if( rstates(cnt) > 0.5 )then
                print *,">>>PROCESSING CLASS ", i
                call cavg_odd%new([params%box,params%box,1], smpd)
                call cavg_odd%new([params%box,params%box,1], smpd)
                call cavg_even%new([params%box,params%box,1], smpd)
                call img_w%new([params%box,params%box,1], smpd)
                call cavg_odd%read(cavgs_stk_even,     i)
                call cavg_even%read(cavgs_stk_odd,     i)
                call spproj%os_ptcl2D%get_pinds(cnt, 'class', pinds)
                if( .not. allocated(pinds) ) cycle
                tmax   = maxval(pinds)
                tmin   = minval(pinds)
                tstamp = tmin + (tmax-tmin)/2
                call cavg_odd%radial_cc(cavg_even, img_w, smpd, rad_cc(cnt,:), rad_dists(cnt,:)) 
                do j = 1, size(rad_dists,dim=2)
                    write(25,'(i6, i6, 2F18.6, i6, i6)') cnt, tstamp, rad_dists(cnt,j), rad_cc(cnt,j), tmax-tmin, size(pinds)
                enddo
            else
                state_mask(cnt) = .false.
            endif
            write(25, '(A)') "          "
        enddo
        command_plot = "gnuplot -e "//'"'//"set pm3d map; set zrange[-.4:1]; splot " //"'"//"evenodd_radial_analysis.csv"//"'"// &
            " u 2:3:4 ; set term png; set xlabel " //"'"//"Time"//"'"// "; set ylabel " //"'"//"Radius({\305})"// &
            "'"// "; set title " //"'"//"Even-Odd Radial Cross-correlation"//"'"// "; set nokey; set output 'radial_analysis.png'; replot" //'"'
        call execute_command_line(command_plot%to_char())
        close(25)
        ! deallocate
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** CAVGSEOPROC_NANO NORMAL STOP ****')
    end subroutine exec_cavgseoproc_nano

    subroutine exec_ptclsproc_nano( self, cline )
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, discrete_read_imgbatch
        class(commander_ptclsproc_nano), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(parameters)     :: params
        type(builder)        :: build
        type(image)          :: cavg_img
        type(sp_project)     :: spproj
        type(string)         :: command_plot
        type(string)         :: cavgs_stk
        integer              :: cnt, iptcl, icavgs, ncavgs, nptcls
        real                 :: smpd, corr
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_ptclsproc_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        call build%init_params_and_build_general_tbox(cline, params)
        ! retrieve number of cavgs stack
        call spproj%get_cavgs_stk(cavgs_stk, ncavgs, smpd, fail=.false.)
        !sz_cls2D = spproj%os_cls2D%get_noris()
        print *, "Number of cavgs ", ncavgs
        do icavgs = 1, ncavgs
            cnt = 0
            call spproj%os_ptcl2D%get_pinds(icavgs, 'class', pinds)
            if( .not. allocated(pinds) ) cycle
            nptcls = size(pinds)
            call cavg_img%new([params%box,params%box,1], smpd)
            call cavg_img%read(cavgs_stk,     icavgs)
            write(logfhandle,'(2(A,I3))')'>>> PROCESSING CLASS ', icavgs, ' NUMBER OF PARTICLES ', nptcls
            if( .not. allocated(pinds) ) cycle
            open(unit=25, file="ptcls_cc_analysis.csv")
            do iptcl = 1, nptcls
                cnt = cnt + 1
                ! compute cross-correlation between particle image and the class average image
                call prepimgbatch(build, nptcls)
                call read_imgbatch(build, [1, nptcls])
                call discrete_read_imgbatch(build, nptcls, pinds(:), [1,nptcls] )
                corr = build%imgbatch(iptcl)%real_corr(cavg_img)
                write(25,'(2i6, 2F18.6)') icavgs, pinds(cnt), real(cnt)/real(nptcls), corr
            enddo
            write(25, '(A)') "          "
            call cavg_img%kill
        enddo
        command_plot = "gnuplot -e "//'"'//"set view map; set zrange[-.4:1]; splot " //"'"//"ptcls_cc_analysis.csv"//"'"// &
            " u 1:2:4 ; set term png; set xlabel " //"'"//"Class average"//"'"// "; set ylabel " //"'"//"Time"// &
            "'"// "; set title " //"'"//"Time Particles Class Analysis"//"'"// "; set nokey; set output 'ptcl_analysis.png'; replot" //'"'
        call execute_command_line(command_plot%to_char())
        call exec_cmdline('rm -rf fsc* fft* stderrout')
        ! deallocate
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** PTCLSPROC_NANO NORMAL STOP ****')
    end subroutine exec_ptclsproc_nano

    subroutine exec_tsegmaps_core_finder( self, cline )
        use simple_opt_mask,  only: estimate_spher_mask
        class(commander_tsegmaps_core_finder), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline !< command line input
        real, parameter   :: RAD_LB = 5.0 ! 5.0 A radial boundary for averaging
        type(parameters)  :: params
        integer           :: ivol, nvols, ldim(3), ifoo, loc(1), best_msk, mskfromto(2)
        real, allocatable :: msk_in_pix(:)
        type(image_msk)   :: mskvol
        type(image)       :: vol, vol_ref, mskimg, mskimg2, mskimg2inv, vol_sum, vol_sum2, vol_w
        type(string), allocatable :: volnames(:)
        call params%new(cline)
        call read_filetable(params%filetab, volnames)
        nvols = size(volnames)
        if( params%mkdir.eq.'yes' )then
            do ivol = 1,nvols
                if(volnames(ivol)%to_char([1,1]).ne.'/') volnames(ivol) = '../'//volnames(ivol)%to_char()
            enddo
        endif
        call find_ldim_nptcls(volnames(1), ldim, ifoo)
        call vol%new(ldim, params%smpd)
        allocate(msk_in_pix(nvols), source=0.)
        do ivol = 1,nvols
            call vol%read(volnames(ivol))
            call mskvol%estimate_spher_mask_diam(vol, AMSKLP_NANO, msk_in_pix(ivol))
            write(logfhandle,*) ivol, 'mask diameter in A: ', 2. * msk_in_pix(ivol) * params%smpd
            call mskvol%kill
        end do
        loc = minloc(msk_in_pix) ! use the smallest NP to drive the radial averaging
        call vol_ref%new(ldim, params%smpd)
        call vol_ref%read(volnames(loc(1)))
        call vol_sum%copy(vol_ref)
        call vol_w%new(ldim, params%smpd)
        vol_w = 1.
        call mskimg%disc(ldim, params%smpd, msk_in_pix(loc(1)))
        mskfromto(1) = int(RAD_LB/params%smpd)
        mskfromto(2) = int(msk_in_pix(loc(1)))
        do ivol = 1,nvols
            if( ivol == loc(1) ) cycle
            call vol%zero_and_unflag_ft
            call vol%read(volnames(ivol))
            call estimate_spher_mask(vol_ref, vol, mskimg, mskfromto, best_msk)
            if( best_msk > mskfromto(1) )then ! agreement beyond the inputted lower radial limit
                call mskimg2%disc(ldim, params%smpd, real(best_msk))
                call vol%mul(mskimg2)
                call vol_sum%add(vol)
                call vol_w%add(mskimg2)
            else
                best_msk = 0
            endif
        enddo
        call vol_sum%div(vol_w)
        call vol_sum%write(string('radial_average_vol_')// int2str_pad(loc(1),2)//'.mrc')
        ! now do the reverse
        call vol_sum2%new(ldim, params%smpd)
        do ivol = 1,nvols
            vol_ref = vol_sum
            call vol%zero_and_unflag_ft
            call vol%read(volnames(ivol))
            call estimate_spher_mask(vol_ref, vol, mskimg, mskfromto, best_msk)
            call mskimg2%disc(ldim, params%smpd, real(best_msk))
            call mskimg2inv%copy(mskimg2)
            call mskimg2inv%bin_inv
            call vol_ref%mul(mskimg2)
            call vol%mul(mskimg2inv)
            call vol_sum2%zero_and_unflag_ft
            call vol_sum2%add(vol)
            call vol_sum2%add(vol_ref)
            call vol_sum2%write(string('core_inserted_vol_')// int2str_pad(ivol,2)//'.mrc')
        end do
    end subroutine exec_tsegmaps_core_finder

    subroutine exec_trajectory_make_projavgs( self, cline )
        use simple_strategy2D3D_common
        class(commander_trajectory_make_projavgs), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters)               :: params
        type(builder)                  :: build
        type(image),       allocatable :: pavgs(:), rot_imgs(:)
        real,              allocatable :: sumw(:), ref_weights(:,:)
        integer,           allocatable :: pinds(:), batches(:,:)
        logical,           allocatable :: ptcl_mask(:)
        real    :: euls_ref(3), euls(3), x, y, sdev_noise, dist, dist_threshold, w
        real    :: spiral_step
        integer :: nbatches, batchsz_max, batch_start, batch_end, batchsz, cnt
        integer :: iptcl, iref, ibatch, nptcls2update, i
        logical :: fall_over
        call cline%set('tseries', 'yes')
        call cline%set('pgrp',    'c1')
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('nspace') ) call cline%set('nspace',  300)
        if( .not. cline%defined('athres') ) call cline%set('athres',  10.)
        call build%init_params_and_build_strategy3D_tbox(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case DEFAULT
                THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        call build%eulspace%set_all2single('state',1.)
        ! allocations
        allocate(pavgs(params%nspace),sumw(params%nspace))
        sumw = 0.0
        !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
        do iref = 1,params%nspace
            call pavgs(iref)%new([params%box,params%box,1],params%smpd,wthreads=.false.)
            call pavgs(iref)%zero_and_unflag_ft
        enddo
        !$omp end parallel do
        ! particles mask and indices
        allocate(ptcl_mask(params%nptcls),source=.true.)
        !$omp parallel do default(shared) private(iptcl) proc_bind(close) schedule(static)
        do iptcl = 1,params%nptcls
            ptcl_mask(iptcl) = build%spproj_field%get_state(iptcl) > 0
        enddo
        !$omp end parallel do
        nptcls2update = count(ptcl_mask)
        allocate(pinds(nptcls2update))
        i = 0
        do iptcl = 1,params%nptcls
            if( ptcl_mask(iptcl) )then
                i        = i+1
                pinds(i) = iptcl
            endif
        enddo
        ! batch prep
        batchsz_max = min(nptcls2update,params_glob%nthr*BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls2update)/real(batchsz_max))
        batches     = split_nobjs_even(nptcls2update, nbatches)
        batchsz_max = maxval(batches(:,2)-batches(:,1)+1)
        call prepimgbatch(build, batchsz_max)
        allocate(ref_weights(batchsz_max,params%nspace),rot_imgs(batchsz_max))
        !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
        do i = 1,batchsz_max
            call rot_imgs(i)%new([params%box,params%box,1],params%smpd,wthreads=.false.)
        enddo
        !$omp end parallel do
        ! angular threshold
        euls_ref       = 0.
        euls           = [0.,params%athres,0.]
        dist_threshold = geodesic_frobdev(euls_ref,euls) / (2.*sqrt(2.))
        spiral_step    = rad2deg(3.809/sqrt(real(params%nspace)))
        do ibatch=1,nbatches
            call progress_gfortran(ibatch,nbatches)
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            ! Weights
            ref_weights = -1.0
            do i = 1,batchsz
                iptcl   = pinds(batch_start+i-1)
                euls    = build%spproj_field%get_euler(iptcl)
                euls(3) = 0.
                do iref = 1,params%nspace
                    if( build%eulspace%get_state(iref) == 0 ) cycle
                    euls_ref    = build%eulspace%get_euler(iref)
                    euls_ref(3) = 0.
                    dist        = geodesic_frobdev(euls_ref,euls) / (2.*sqrt(2.)) ! => [0;1]
                    if( dist > dist_threshold ) cycle
                    w = 1. / (1. + 999.0*dist)
                    ref_weights(i,iref) = w
                enddo
            enddo
            ! Images prep
            call discrete_read_imgbatch(build, batchsz, pinds(batch_start:batch_end), [1,batchsz])
            !$omp parallel do default(shared) private(i,iptcl,x,y,sdev_noise) proc_bind(close) schedule(static)
            do i = 1, batchsz
                iptcl = pinds(batch_start+i-1)
                call build%imgbatch(i)%norm_noise(build%lmsk, sdev_noise)
                x = build%spproj_field%get(iptcl, 'x')
                y = build%spproj_field%get(iptcl, 'y')
                call build%imgbatch(i)%fft
                call build%imgbatch(i)%shift2Dserial(-[x,y])
                call build%imgbatch(i)%ifft
                call rot_imgs(i)%zero_and_flag_ft
                call build%imgbatch(i)%rtsq(-build%spproj_field%e3get(iptcl),0.,0.,rot_imgs(i))
            enddo
            !$omp end parallel do
            ! Projection direction weighted sum
            !$omp parallel do default(shared) private(i,iref,w) proc_bind(close) schedule(static)
            do iref = 1,params%nspace
                if( build%eulspace%get_state(iref) == 0 ) cycle
                do i = 1,batchsz
                    w = ref_weights(i,iref)
                    if( w < TINY ) cycle
                    sumw(iref) = sumw(iref) + w
                    call pavgs(iref)%add(rot_imgs(i),w=w)
                enddo
            enddo
            !$omp end parallel do
        enddo
        ! Weights normalization
        !$omp parallel do default(shared) private(iref,w) proc_bind(close) schedule(static)
        do iref = 1,params%nspace
            if( build%eulspace%get_state(iref) == 0 ) cycle
            if( sumw(iref) > 0.001 )then
                call pavgs(iref)%div(sumw(iref))
            else
                call pavgs(iref)%zero_and_unflag_ft
            endif
        enddo
        !$omp end parallel do
        ! write
        cnt = 0
        do iref = 1,params%nspace
            if( build%eulspace%get_state(iref) == 0 ) cycle
            cnt = cnt + 1
            call pavgs(iref)%write(params%outstk,cnt)
        enddo
        call build%eulspace%write(string('projdirs.txt'))
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_PROJAVGS NORMAL STOP ****')
    end subroutine exec_trajectory_make_projavgs

  end module single_commanders_experimental

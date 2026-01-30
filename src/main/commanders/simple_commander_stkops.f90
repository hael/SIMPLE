!@descr: operations on image stacks
module simple_commanders_stkops
use simple_commander_module_api
use simple_pftc_srch_api
use simple_strategy2D_utils, only: calc_cluster_cavgs_dmat, calc_match_cavgs_dmat
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_convert
  contains
    procedure :: execute      => exec_convert
end type commander_convert

type, extends(commander_base) :: commander_stack
  contains
    procedure :: execute      => exec_stack
end type commander_stack

type, extends(commander_base) :: commander_stackops
  contains
    procedure :: execute      => exec_stackops
end type commander_stackops

type, extends(commander_base) :: commander_cluster_stack
  contains
    procedure :: execute      => exec_cluster_stack
end type commander_cluster_stack

type, extends(commander_base) :: commander_match_stacks
  contains
    procedure :: execute      => exec_match_stacks
end type commander_match_stacks

contains

    !> for converting between SPIDER and MRC formats
    subroutine exec_convert( self, cline )
        class(commander_convert), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: ldim(3), nptcls, iptcl
        if( cline%defined('stk') )then
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            call find_ldim_nptcls(params%stk,ldim,nptcls)
            ldim(3) = 1
            if( ldim(1) /= ldim(2) )then
                call build%img%new(ldim,params%smpd)
            endif
            do iptcl=1,params%nptcls
                call progress(iptcl, params%nptcls)
                call build%img%read(params%stk, iptcl)
                call build%img%write(params%outstk, iptcl)
            end do
        else if( cline%defined('vol1') )then
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
            call build%vol%read(params%vols(1))
            call build%img%write(params%outvol)
        else
            THROW_HARD('either vol1 or stk argument required to execute simple_convert')
        endif
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CONVERT NORMAL STOP ****')
    end subroutine exec_convert

    !>  for stacking individual images or multiple stacks into one
    subroutine exec_stack( self, cline )
        class(commander_stack), intent(inout)  :: self
        class(cmdline),         intent(inout)  :: cline
        type(string), allocatable :: filenames(:)
        type(parameters) :: params
        integer          :: nfiles, ldim(3), ifile, ifoo, cnt
        integer          :: lfoo(3), nimgs, iimg
        type(stack_io)   :: stkio_r, stkio_w
        type(image)      :: img, tmp
        logical          :: l_clip
        call cline%set('mkdir', 'no')
        call params%new(cline)
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        call find_ldim_nptcls(filenames(1),ldim,ifoo)
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        if( ldim(1) /= ldim(2) )then
            ! prepare img and tmp for reading
            call img%new(ldim, params%smpd)
            ! loop over files
            cnt = 0
            do ifile=1,nfiles
                if( .not. file_exists(filenames(ifile)) )then
                    write(logfhandle,*) 'inputted file does not exist: ', filenames(ifile)%to_char()
                endif
                call find_ldim_nptcls(filenames(ifile),lfoo,nimgs)
                do iimg=1,nimgs
                    cnt = cnt + 1
                    call img%read(filenames(ifile), iimg)
                    call img%write(params%outstk, cnt)
                end do
                call progress(ifile, nfiles)
            end do
        else
            params%box = ldim(1)
            ! prepare img and tmp for reading
            call img%new(ldim, params%smpd)
            l_clip = cline%defined('clip')
            if( l_clip )then
                call tmp%new([params%clip,params%clip,1], params%smpd)
                call stkio_w%open(params%outstk, params%smpd, 'write', box=params%clip)
            else
                call stkio_w%open(params%outstk, params%smpd, 'write', box=ldim(1))
            endif
            ! loop over files
            cnt = 0
            do ifile=1,nfiles
                if( .not. file_exists(filenames(ifile)) )then
                    write(logfhandle,*) 'inputted file does not exist: ', filenames(ifile)%to_char()
                endif
                call find_ldim_nptcls(filenames(ifile),lfoo,nimgs)
                do iimg=1,nimgs
                    cnt = cnt+1
                    if( .not. stkio_r%stk_is_open() )then
                        call stkio_r%open(filenames(ifile), params%smpd, 'read')
                    else if( .not. stkio_r%same_stk(filenames(ifile), ldim) )then
                        call stkio_r%close
                        call stkio_r%open(filenames(ifile), params%smpd, 'read')
                    endif
                    call stkio_r%read(iimg, img)
                    if( l_clip )then
                        call img%clip(tmp)
                        call stkio_w%write(cnt, tmp)
                    else
                        call stkio_w%write(cnt, img)
                    endif
                end do
                call progress(ifile, nfiles)
            end do
            call stkio_r%close
            call stkio_w%close
        endif
        call img%kill
        call tmp%kill
        ! end gracefully
        call simple_end('**** SIMPLE_STACK NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_stack 

    subroutine exec_stackops( self, cline )
        use simple_stackops
        use simple_procimgstk
        class(commander_stackops), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters)          :: params
        type(builder)             :: build
        type(ran_tabu)            :: rt
        type(oris)                :: o_here, os_ran
        type(ori)                 :: o, o2
        integer,      allocatable :: pinds(:)
        complex(sp),  pointer     :: ptre(:,:,:), ptro(:,:,:)
        type(string)              :: fname
        type(polarft_calc)        :: pftc
        type(image)               :: img
        integer :: i, s, cnt, nincl, j, k
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! ------------------------------------------------------------------
        ! main control flow: first matching branch wins
        ! ------------------------------------------------------------------
        if( cline%defined('polar') ) then
            call handle_polar_representation()
        else if( cline%defined('nran') ) then
            call handle_random_selection()
        else if( cline%defined('frac') ) then
            call handle_fishing_frac_only()
        else if( (cline%defined('state') .or. cline%defined('class')) .and. .not. cline%defined('frac') ) then
            call handle_state_class_without_frac()
        else if( params%neg .eq. 'yes' ) then
            ! invert contrast
            call neg_imgfile(params%stk, params%outstk, params%smpd)
        else if( params%acf .eq. 'yes' ) then
            ! auto correlation function
            call acf_stack(params%stk, params%outstk)
        else if( params%stats .eq. 'yes' ) then
            ! produce image statistics
            call stats_imgfile(params%stk, build%spproj_field)
            call build%spproj_field%write(string('image_statistics.txt'))
        else if( params%loc_sdev .eq. 'yes' ) then
            ! produce local standard deviations
            call loc_sdev_imgfile(params%stk, params%outstk, nint(params%winsz), params%smpd)
        else if( params%nquanta > 0 ) then
            ! quantize
            call quantize_imgfile(params%stk, params%outstk, params%nquanta, params%smpd)
        else if( params%nframesgrp > 0 ) then
            ! create frame averages
            call frameavg_stack(params%stk, params%outstk, params%nframesgrp, params%smpd)
        else if( params%vis .eq. 'yes' ) then
            ! visualize
            do i=1,params%nptcls
                call build%img%read(params%stk, i)
                call build%img%vis()
            end do
        else if( params%avg .eq. 'yes' ) then
            ! average
            call make_avg_stack(params%stk, params%outstk, params%smpd)
        else if( params%roavg .eq. 'yes' ) then
            ! rotationally average
            call roavg_imgfile(params%stk, params%outstk, params%angstep, params%smpd)
        else if( cline%defined('snr') ) then
            ! add noise
            call add_noise_imgfile(params%stk, params%outstk, params%snr, params%smpd)
        else if( cline%defined('top') .and. .not. cline%defined('part') ) then
            ! copy
            call copy_imgfile(params%stk, params%outstk, params%smpd, fromto=[params%fromp,params%top])
        else if( params%mirr .ne. 'no' ) then
            ! mirror
            call mirror_imgfile(params%stk, params%outstk, params%mirr, params%smpd)
        else if( params%makemovie .eq. 'yes' ) then
            ! prepare for movie
            call prep_imgfile4movie(params%stk, params%smpd)
        else
            ! default
            write(logfhandle,*) 'Nothing to do!'
        end if
        ! end gracefully
        call finish_cleanup()

    contains

        ! ------------------------------------------------------------------
        ! POLAR REPRESENTATION
        ! ------------------------------------------------------------------
        subroutine handle_polar_representation()
            complex, allocatable :: pft(:,:)
            call pftc%new(params%nptcls, [1,params%nptcls], params%kfromto)
            call  build%img_crop%memoize4polarize(pftc%get_pdim(), params%alpha)
            pft = pftc%allocate_pft()
            do i = 1, params%nptcls
                call build%img_crop%read(params%stk, i)
                call build%img_crop%fft
                call build%img_crop%polarize(pft, mask=build%l_resmsk)
                call pftc%set_ref_pft(i, pft, iseven=.true.)
            end do
            call img%new([pftc%get_pftsz(), params%kfromto(2), 1], 1.0)
            do i = 1, params%nptcls
                call pftc%get_ref_pft(i, .true.,  pft)
                img = 0.
                do j = 1, pftc%get_pftsz()
                    do k = params%kfromto(1), params%kfromto(2)
                        call img%set([j,k,1], real(abs(pft(j,k))))
                    end do
                end do
                call img%write(params%outstk, i)
            end do
            call img%kill
            if( allocated(pft) ) deallocate(pft)
        end subroutine handle_polar_representation

        ! ------------------------------------------------------------------
        ! RANDOM SELECTION
        ! ------------------------------------------------------------------
        subroutine handle_random_selection()
            if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> RANDOMLY SELECTING IMAGES'
            allocate( pinds(params%nran) )
            rt = ran_tabu(params%nptcls)
            call rt%ne_ran_iarr(pinds)
            if( cline%defined('oritab') .or. cline%defined('deftab') ) then
                call del_file(params%outfile)
            end if
            call os_ran%new(params%nran, is_ptcl=.false.)
            do i=1,params%nran
                call progress(i, params%nran)
                call build%img%read(params%stk, pinds(i))
                call build%img%write(params%outstk, i)
                if( cline%defined('oritab') .or. cline%defined('deftab') ) then
                    call build%spproj_field%get_ori(pinds(i), o)
                    call os_ran%set_ori(i, o)
                    call o%kill
                end if
            end do
            if( cline%defined('oritab') .or. cline%defined('deftab') ) then
                call build%spproj_field%write(params%outfile, [1,params%nran])
            end if
        end subroutine handle_random_selection

        ! ------------------------------------------------------------------
        ! FISHING EXPEDITIONS – FRAC ONLY
        ! ------------------------------------------------------------------
        subroutine handle_fishing_frac_only()
            if( params%oritab == '' ) &
                THROW_HARD('need input orientation doc for fishing expedition; simple_stackops')
            ! determine how many particles to include
            if( params%frac < 0.99 ) then
                nincl = nint(real(params%nptcls)*params%frac)
            else
                nincl = params%nptcls
            end if
            ! order the particles
            pinds = build%spproj_field%order()
            ! fish the best ones out
            if( cline%defined('state') ) then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'state'))
                    if( s == params%state ) then
                        cnt = cnt+1
                        call build%img%read(params%stk, pinds(i))
                        call build%img%write(params%outstk, cnt)
                    end if
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt, is_ptcl=.true.)
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'state'))
                    if( s == params%state ) then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(pinds(i), o2)
                        call o_here%set_ori(cnt, o2)
                    end if
                end do
                fname = 'extracted_oris_state'//int2str_pad(params%state,2)//TXT_EXT
            else if( cline%defined('class') ) then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'class'))
                    if( s == params%class ) then
                        cnt = cnt+1
                        call build%img%read(params%stk, pinds(i))
                        call build%img%write(params%outstk, cnt)
                    end if
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt, is_ptcl=.true.)
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'class'))
                    if( s == params%class ) then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(pinds(i), o2)
                        call o_here%set_ori(cnt, o2)
                    end if
                end do
                fname = 'extracted_oris_class'//int2str_pad(params%class,5)//TXT_EXT
            else
                o_here = oris(nincl, is_ptcl=.true.)
                do i=1,nincl
                    call progress(i, nincl)
                    call build%img%read(params%stk, pinds(i))
                    call build%img%write(params%outstk, i)
                    call build%spproj_field%get_ori(pinds(i), o2)
                    call o_here%set_ori(i, o2)
                end do
                fname = 'extracted_oris'//TXT_EXT
            end if
            if( cline%defined('outfile') ) then
                call o_here%write(params%outfile, [1,o_here%get_noris()])
            else
                call o_here%write(fname, [1,o_here%get_noris()])
            end if
        end subroutine handle_fishing_frac_only

        ! ------------------------------------------------------------------
        ! FISHING EXPEDITIONS – STATE/CLASS WITHOUT FRAC
        ! ------------------------------------------------------------------
        subroutine handle_state_class_without_frac()
            if( params%oritab == '' ) &
                THROW_HARD('need input orientation doc for fishing expedition; simple_stackops')
            if( cline%defined('state') ) then
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'state'))
                    if( s == params%state ) then
                        cnt = cnt+1
                        call build%img%read(params%stk, i)
                        call build%img%write(params%outstk, cnt)
                    end if
                end do
                ! make orientation structure for the extracted ones
                o_here = oris(cnt, is_ptcl=.true.)
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'state'))
                    if( s == params%state ) then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(i, o2)
                        call o_here%set_ori(cnt, o2)
                    end if
                end do
                fname = 'extracted_oris_state'//int2str_pad(params%state,2)//TXT_EXT
            else if( cline%defined('class') ) then
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'class'))
                    if( s == params%class ) then
                        cnt = cnt+1
                        call build%img%read(params%stk, i)
                        call build%img%write(params%outstk, cnt)
                    end if
                end do
                ! make orientation structure for the extracted ones
                o_here = oris(cnt, is_ptcl=.true.)
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'class'))
                    if( s == params%class ) then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(i, o2)
                        call o_here%set_ori(cnt, o2)
                    end if
                end do
                fname = 'extracted_oris_class'//int2str_pad(params%class,5)//TXT_EXT
            end if
            if( cline%defined('outfile') ) then
                call o_here%write(params%outfile, [1,o_here%get_noris()])
            else
                call o_here%write(fname, [1,o_here%get_noris()])
            end if
        end subroutine handle_state_class_without_frac

        ! ------------------------------------------------------------------
        ! COMMON CLEANUP
        ! ------------------------------------------------------------------
        subroutine finish_cleanup()
            call o_here%kill
            call os_ran%kill
            call o%kill
            call o2%kill
            call rt%kill
            call build%kill_general_tbox
            call simple_end('**** SIMPLE_STACKOPS NORMAL STOP ****')
        end subroutine finish_cleanup

    end subroutine exec_stackops

    subroutine exec_cluster_stack( self, cline )
        use simple_clustering_utils, only: cluster_dmat
        class(commander_cluster_stack), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)              :: params
        type(image),      allocatable :: stk_imgs(:)
        integer,          allocatable :: labels(:), i_medoids(:)
        real,             allocatable :: mm(:, :), dmat(:, :)
        real    :: smpd, mskrad, oa_min, oa_max
        integer :: i, ldim(3), box, nimgs, nclust
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        if( .not. cline%defined('hp')         ) THROW_HARD('Need high-pass limit (hp) on command line!')
        if( .not. cline%defined('lp')         ) THROW_HARD('Need low-pass limit (lp) on command line!')
        if( .not. cline%defined('mskdiam')    ) THROW_HARD('Need mask diameter in A (mskdiam) on command line!')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',         10.)
        if( .not. cline%defined('clust_crit') ) call cline%set('clust_crit', 'fm')
        ! master parameters
        call params%new(cline)
        ! prep
        stk_imgs = read_stk_into_imgarr(params%stk)
        nimgs    = size(stk_imgs)
        smpd     = stk_imgs(1)%get_smpd()
        ldim     = stk_imgs(1)%get_ldim()
        box      = ldim(1)
        mskrad   = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * params%mskdiam/smpd)
        ! create the stuff needed in the loop
        allocate(mm(nimgs,2), source=0.)
        call stk_imgs(1)%memoize_mask_coords
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nimgs
            ! normalization
            call stk_imgs(i)%norm
            ! mask
            call stk_imgs(i)%mask2D_soft(mskrad, backgr=0.)
            ! stash minmax
            mm(i,:) = stk_imgs(i)%minmax(mskrad)
        end do
        !$omp end parallel do
        ! calculate overall minmax
        oa_min      = minval(mm(:,1))
        oa_max      = maxval(mm(:,2))
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! calculate distance matrix
        dmat = calc_cluster_cavgs_dmat(params, stk_imgs, [oa_min,oa_max], params%clust_crit)
        ! cluster
        if( cline%defined('ncls') )then
            nclust = params%ncls
            call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
        else
            call cluster_dmat( dmat, 'aprop', nclust, i_medoids, labels)
        endif
        call arr2file(real(i_medoids), string(CLUST_MEDIODS_FNAME))
        call dealloc_imgarr(stk_imgs)
        stk_imgs = read_stk_into_imgarr(params%stk)
        call write_imgarr( nimgs, stk_imgs, labels, 'cluster', '.mrcs' )
    end subroutine exec_cluster_stack

    subroutine exec_match_stacks( self, cline )
        class(commander_match_stacks), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(image),       allocatable :: stk_imgs_ref(:), stk_imgs_match(:)
        integer,           allocatable :: labels_match(:)
        real,              allocatable :: mm_ref(:, :), dmat(:, :)
        integer :: nmatch, nrefs, ldim(3), i, nclust, imatch, box
        real    :: smpd, oa_minmax(2), mskrad
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        if( .not. cline%defined('hp')         ) THROW_HARD('Need high-pass limit (hp) on command line!')
        if( .not. cline%defined('lp')         ) THROW_HARD('Need low-pass limit (lp) on command line!')
        if( .not. cline%defined('mskdiam')    ) THROW_HARD('Need mask diameter in A (mskdiam) on command line!')
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',         10.)
        if( .not. cline%defined('clust_crit') ) call cline%set('clust_crit', 'cc')
        ! master parameters
        call params%new(cline)
        ! prep
        stk_imgs_ref   = read_stk_into_imgarr(params%stk)
        stk_imgs_match = read_stk_into_imgarr(params%stk2)
        nrefs          = size(stk_imgs_ref)
        nclust         = nrefs
        nmatch         = size(stk_imgs_match)
        smpd           = stk_imgs_ref(1)%get_smpd()
        ldim           = stk_imgs_ref(1)%get_ldim()
        box            = ldim(1)
        mskrad         = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * params%mskdiam/smpd)
        ! create the stuff needed in the loop
        allocate(mm_ref(nrefs,2), source=0.)
        call stk_imgs_ref(1)%memoize_mask_coords
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, nrefs
            ! normalization
            call stk_imgs_ref(i)%norm
            ! mask
            call stk_imgs_ref(i)%mask2D_soft(mskrad, backgr=0.)
            ! stash minmax
            mm_ref(i,:) = stk_imgs_ref(i)%minmax(mskrad)
        end do
        !$omp end parallel do
        ! calculate overall minmax
        oa_minmax(1) = minval(mm_ref(:,1))
        oa_minmax(2) = maxval(mm_ref(:,2))
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! generate matching distance matrix
        dmat = calc_match_cavgs_dmat(params, stk_imgs_ref, stk_imgs_match, oa_minmax, params%clust_crit)
        ! genrate cluster distance matrix
        allocate(labels_match(nmatch), source=0)
        labels_match = minloc(dmat, dim=1)

        print *, 'size(lables_match)       : ', size(labels_match)
        print *, 'nmatch                   : ', nmatch
        print *, 'size(minloc(dmat, dim=1)): ', size(minloc(dmat, dim=1))

        ! write matching references
        do imatch = 1, nmatch
            call stk_imgs_ref(labels_match(imatch))%write(string('matching_refs.mrcs'), imatch)
        end do
        call dealloc_imgarr(stk_imgs_ref)
        call dealloc_imgarr(stk_imgs_match)
    end subroutine exec_match_stacks

end module simple_commanders_stkops

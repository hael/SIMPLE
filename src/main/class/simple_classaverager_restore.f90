!@descr: Routines to perform the classes restoration and processing
submodule (simple_classaverager) simple_classaverager_restore
use simple_imgarr_utils,     only: alloc_imgarr, dealloc_imgarr
use simple_gridding,         only: prep2D_inv_instrfun4mul
use simple_oris,             only: population_blend_weights
use simple_cavg_sums,        only: cavg_sums, CAVG_SUMS_STATE, CAVG_SUMS_CONTRIB, CAVG_SUMS_OK,&
                                  &cavg_contrib_fname
use simple_euclid_sigma2,    only: sigma2_half
use simple_fplane_grid,      only: fplane_grid
implicit none
#include "simple_local_flags.inc"

! restoration state
type(ptcl_record),   allocatable :: precs(:)                  !< Particle records
type(image),         allocatable :: tmp_pad_imgs(:)           !< Temporary images for on-the-fly classes update
type(cavgs_set)                  :: cavgs                     !< Class averages
type(builder),        pointer    :: b_ptr  => null()          !< active builder instance
integer,             allocatable :: eo_pops(:,:)              !< Even/odd class populations
real,                allocatable :: center_offsets(:,:)       !< Class-centering offsets of the references (2,ncls)
integer                          :: ncls       = 0            !< # classes
integer                          :: ldim(3)        = [0,0,0]  !< logical dimension of image
integer                          :: ldim_crop(3)   = [0,0,0]  !< logical dimension of cropped image
integer                          :: ldim_pd(3)     = [0,0,0]  !< logical dimension of image, padded
integer                          :: ldim_croppd(3) = [0,0,0]  !< logical dimension of cropped image, padded
real                             :: smpd       = 0.           !< sampling distance
real                             :: smpd_crop  = 0.           !< cropped sampling distance
logical                          :: l_cropped_ptcls = .false. !< particles handed to cavger_update_sums are at box_crop

contains

    !>  \brief  Constructor
    module subroutine cavger_new( params, build )
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        p_ptr => params
        b_ptr => build
        call cavger_kill
        ncls          = p_ptr%ncls
        ! smpd
        smpd          = p_ptr%smpd
        smpd_crop     = p_ptr%smpd_crop
        ! set ldims
        ldim          = [p_ptr%box,       p_ptr%box,       1]
        ldim_crop     = [p_ptr%box_crop,  p_ptr%box_crop,  1]
        ldim_croppd   = [p_ptr%box_croppd,p_ptr%box_croppd,1]
        ldim_pd       = [p_ptr%boxpd,     p_ptr%boxpd,     1]
        ! instantiate class averages
        call cavgs%new_set(ldim_crop(1:2), ncls)
        ! populations
        allocate(eo_pops(2,ncls),source=0)
        ! class-centering offsets of the references, applied once to the carried sums
        allocate(center_offsets(2,ncls),source=0.)
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers metadata to the instance for a subset of particles
    module subroutine cavger_transf_oridat( nptcls, pinds, updated_only )
        integer,            intent(in) :: nptcls
        integer,            intent(in) :: pinds(nptcls)
        logical,  optional, intent(in) :: updated_only
        class(oris), pointer :: spproj_field
        integer              :: i, iptcl, stkind
        logical              :: l_updated_only
        l_updated_only = .false.
        if( present(updated_only) ) l_updated_only = updated_only
        ! indices
        precs(1:nptcls)%pind = pinds(:)
        if( nptcls < size(precs) ) precs(nptcls+1:)%pind = 0
        ! fetch data from project
        call b_ptr%spproj%ptr2oritype(p_ptr%oritype, spproj_field)
        !$omp parallel do default(shared) private(i,iptcl,stkind) schedule(static) proc_bind(close)
        do i = 1,nptcls
            iptcl = pinds(i)
            if( iptcl == 0 ) cycle
            if( l_updated_only )then
                if( spproj_field%get_updatecnt(iptcl) == 0 )then
                    precs(i)%pind = 0
                    cycle
                endif
            endif
            if( spproj_field%get_state(iptcl)==0 )then
                precs(i)%pind   = 0
            else
                precs(i)%pind   = iptcl
                precs(i)%eo     = spproj_field%get_eo(iptcl)
                precs(i)%ctfparams = b_ptr%spproj%get_ctfparams(p_ptr%oritype,iptcl)
                precs(i)%tfun   = ctf(p_ptr%smpd_crop, precs(i)%ctfparams%kv, precs(i)%ctfparams%cs, precs(i)%ctfparams%fraca)
                precs(i)%class  = spproj_field%get_class(iptcl)
                precs(i)%e3     = spproj_field%e3get(iptcl)
                precs(i)%shift  = spproj_field%get_2Dshift(iptcl)
                call b_ptr%spproj%map_ptcl_ind2stk_ind(p_ptr%oritype, iptcl, stkind, precs(i)%ind_in_stk)
            endif
        end do
        !$omp end parallel do
        nullify(spproj_field)
    end subroutine cavger_transf_oridat

    !>  \brief  for loading sigma2
    module subroutine cavger_read_euclid_sigma2
        type(string) :: fname
        logical :: found
        if( p_ptr%l_ml_reg )then
            call b_ptr%spproj%get_sigma2_state_path(fname, found)
            if( .not. found ) THROW_HARD('particle project has no canonical sigma2 state path')
            ! restoration whitens by the two global spectra; it reads no particle records
            call b_ptr%esig%new(p_ptr, b_ptr%pftc, fname, p_ptr%box)
            call b_ptr%esig%read_global
            call fname%kill
        end if
    end subroutine cavger_read_euclid_sigma2

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population and average correlation
    module subroutine cavger_gen2Dclassdoc()
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(p_ptr%ncls)
        real(dp) :: corrs(p_ptr%ncls)
        real     :: frc05, frc0143
        integer  :: iptcl, icls, pop, nptcls
        select case(trim(p_ptr%oritype))
            case('ptcl2D')
                ptcl_field => b_ptr%spproj%os_ptcl2D
                cls_field  => b_ptr%spproj%os_cls2D
            case('ptcl3D')
                ptcl_field => b_ptr%spproj%os_ptcl3D
                cls_field  => b_ptr%spproj%os_cls3D
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE: '//trim(p_ptr%oritype))
        end select
        nptcls = ptcl_field%get_noris()
        pops  = 0
        corrs = 0.d0
        !$omp parallel do default(shared) private(iptcl,icls) schedule(static)&
        !$omp proc_bind(close) reduction(+:pops,corrs)
        do iptcl=1,nptcls
            if( ptcl_field%get_state(iptcl) == 0 ) cycle
            icls = ptcl_field%get_class(iptcl)
            if( icls<1 .or. icls>p_ptr%ncls )cycle
            pops(icls)  = pops(icls)  + 1
            corrs(icls) = corrs(icls) + real(ptcl_field%get(iptcl,'corr'),dp)
        enddo
        !$omp end parallel do
        where(pops>1)
            corrs = corrs / real(pops)
        elsewhere
            corrs = -1.
        end where
        call cls_field%new(p_ptr%ncls, is_ptcl=.false.)
        do icls=1,p_ptr%ncls
            pop = pops(icls)
            call b_ptr%clsfrcs%estimate_res(icls, frc05, frc0143)
            call ptcl_field%set_field2single('class', icls, 'res', frc0143)
            call cls_field%set(icls, 'class',     icls)
            call cls_field%set(icls, 'pop',       pop)
            call cls_field%set(icls, 'res',       frc0143)
            call cls_field%set(icls, 'corr',      corrs(icls))
            if( pop > 0 )then
                call cls_field%set(icls, 'state', 1.0) ! needs to be default val if no selection has been done
            else
                call cls_field%set(icls, 'state', 0.0) ! exclusion
            endif
        end do
    end subroutine cavger_gen2Dclassdoc

    ! Calculators

    ! Initialize objects for on-the-fly classes update. The sums always start from zero:
    ! workers and the shared-memory matcher accumulate the current sample only, and the
    ! assembly owner blends the carried sums (cavger_commit_carryover)
    module subroutine cavger_init_online( maxbatchsz, cropped_ptcls )
        integer,           intent(in) :: maxbatchsz
        logical, optional, intent(in) :: cropped_ptcls
        ! Whether cavger_update_sums will be fed box_crop particles (from the
        ! downscaled cache) rather than the full-size originals. Set explicitly by
        ! the caller rather than inferred, so the offline assembly path is unaffected.
        l_cropped_ptcls = .false.
        if( present(cropped_ptcls) ) l_cropped_ptcls = cropped_ptcls
        ! Zero sums
        call cavgs%zero_set(.true.)
        ! Work images. box_croppd and boxpd cover the same physical extent, so their
        ! Fourier grids share a spacing and index hp means the same spatial frequency
        ! in both; stack_accumulate_fplane only ever reads |hp| <= 2*nyq of the class
        ! average, which is exactly the extent of the cropped padded grid.
        if( l_cropped_ptcls )then
            call alloc_imgarr(nthr_glob, ldim_croppd, smpd_crop, tmp_pad_imgs)
        else
            call alloc_imgarr(nthr_glob, ldim_pd, smpd, tmp_pad_imgs)
        endif
        ! particle records
        allocate(precs(maxbatchsz))
        precs(:)%pind = 0
        ! populations
        eo_pops(:,:) = 0
        ! Memoization for cropped padded image, will be overwritten during search
        call memoize_ft_maps(ldim_croppd(1:2), p_ptr%smpd_crop)
    end subroutine cavger_init_online

    ! Deallocate objects  on-the-fly classes update
    module subroutine cavger_dealloc_online()
        if( allocated(tmp_pad_imgs))then
            call forget_ft_maps
            call dealloc_imgarr(tmp_pad_imgs)
            deallocate(precs)
        endif
    end subroutine cavger_dealloc_online

    subroutine cavger_update_sums( nptcls, ptcl_imgs )
        integer,      intent(in)    :: nptcls
        class(image), intent(inout) :: ptcl_imgs(nptcls)
        type(fplane_type) :: fplanes(nthr_glob)
        type(fplane_grid) :: grid
        type(ctfparams) :: ctfparms_here
        integer, allocatable :: cls_first(:), cls_order(:)
        integer :: iptcl, ithr, icls, i, j, ih, box_pd
        real    :: crop_factor, shift_here(2)
        ! Memoization for the padded image tmp_pad_imgs was allocated with. The
        ! (ldim,smpd) pair fixes the physical extent, and boxpd*smpd == box_croppd*
        ! smpd_crop, so a given index h denotes the same spatial frequency either way.
        if( l_cropped_ptcls )then
            call memoize_ft_maps(ldim_croppd(1:2), p_ptr%smpd_crop)
        else
            call memoize_ft_maps(ldim_pd(1:2), p_ptr%smpd)
        endif
        crop_factor = real(p_ptr%box_crop) / real(p_ptr%box)
        ! the plane grid of the padded box, once for the call: the spectra's cached one with ML regularization;
        ! its spectra span the shells of the box it is built for (box_crop/2 for cropped particles)
        box_pd = merge(ldim_croppd(1), ldim_pd(1), l_cropped_ptcls)
        if( p_ptr%l_ml_reg )then
            call b_ptr%esig%prep_fplane_grid(box_pd, lbound(b_ptr%esig%sigma2_noise,1))
        else
            call grid%new(box_pd)
        endif
        ! the batch bucketed by class once (a stable counting sort): class icls holds the batch indices
        ! cls_order(cls_first(icls):cls_first(icls+1)-1), in batch order, so each class's sums
        ! accumulate in the same order as a scan over the whole batch
        allocate(cls_first(ncls+1), source=0)
        allocate(cls_order(nptcls))
        do i = 1, nptcls
            if( precs(i)%pind == 0 ) cycle
            icls = precs(i)%class
            if( icls < 1 .or. icls > ncls ) cycle
            cls_first(icls+1) = cls_first(icls+1) + 1
        enddo
        cls_first(1) = 1
        do icls = 1, ncls
            cls_first(icls+1) = cls_first(icls+1) + cls_first(icls)
        enddo
        ! cls_first(icls+1) - 1 is the last slot of class icls; fill from the front with a running cursor
        cls_order = 0
        do i = 1, nptcls
            if( precs(i)%pind == 0 ) cycle
            icls = precs(i)%class
            if( icls < 1 .or. icls > ncls ) cycle
            cls_order(cls_first(icls)) = i
            cls_first(icls) = cls_first(icls) + 1
        enddo
        ! the cursors now point at the next class's first slot: shift them back by one class
        do icls = ncls, 1, -1
            cls_first(icls+1) = cls_first(icls)
        enddo
        cls_first(1) = 1
        !$omp parallel do default(shared) private(icls,j,i,iptcl,ithr,ih,ctfparms_here,shift_here)&
        !$omp schedule(static,1) proc_bind(close)
        do icls = 1, ncls
            do j = cls_first(icls), cls_first(icls+1) - 1
                i     = cls_order(j)
                ithr  = omp_get_thread_num() + 1
                iptcl = precs(i)%pind
                ! particle: normalize, pad & forward FT
                ! The mask has to match the box the particle actually has. Cached
                ! particles were already noise-normalized at the full box before being
                ! Fourier-cropped, and cropping discards most of the noise power, so
                ! re-normalizing here would scale them up by the crop factor while
                ! leaving ctfsq_plane and sigma2 alone.
                if( l_cropped_ptcls )then
                    call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(b_ptr%lmsk_crop, &
                        &tmp_pad_imgs(ithr), renorm=.false.)
                else
                    call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(b_ptr%lmsk, tmp_pad_imgs(ithr))
                endif
                ! shift, CTF and ML regularization in Fourier plane generation
                ctfparms_here = precs(i)%ctfparams
                shift_here    = precs(i)%shift
                if( l_cropped_ptcls )then
                    ! shconst is PI/(ldim/2), so shifts must be in the pixel units of
                    ! the padded box. boxpd shares the original pixel size and needs no
                    ! conversion; box_croppd is in smpd_crop pixels and does. Likewise
                    ! the CTF kernel reads spatial frequency in cycles/pixel of the
                    ! current grid, so the CTF must be told the cropped pixel size.
                    ctfparms_here%smpd = ctfparms_here%smpd / crop_factor != smpd_crop
                    shift_here         = shift_here * crop_factor
                endif
                if( p_ptr%l_ml_reg ) then
                    ih = sigma2_half(precs(i)%eo)
                    call tmp_pad_imgs(ithr)%gen_fplane4rec(p_ptr%smpd_crop, ctfparms_here, shift_here, &
                        &fplanes(ithr), b_ptr%esig%fplane, half=ih)
                else
                    call tmp_pad_imgs(ithr)%gen_fplane4rec(p_ptr%smpd_crop, ctfparms_here, shift_here, &
                        &fplanes(ithr), grid)
                endif
                ! rotation, interpolation and accumulation
                select case(precs(i)%eo)
                case(0,-1)
                    eo_pops(1,icls) = eo_pops(1,icls) + 1
                    call cavgs%even%accumulate_fplane(precs(i)%e3, fplanes(ithr), icls)
                case(1)
                    eo_pops(2,icls) = eo_pops(2,icls) + 1
                    call cavgs%odd%accumulate_fplane(precs(i)%e3, fplanes(ithr), icls)
                end select
            enddo
        enddo
        !$omp end parallel do
        call grid%kill
        deallocate(cls_first, cls_order)
    end subroutine cavger_update_sums

    !>  \brief  is for generating class averages offline
    module subroutine cavger_assemble_sums()
        use simple_matcher_ptcl_io, only: prepimgbatch, ptcl_read_session, killimgbatch, &
            &ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since
        class(oris), pointer :: spproj_field
        type(string)         :: source_stk
        type(ptcl_read_counts) :: reads0
        type(ptcl_read_session) :: reads
        integer :: pinds(READBUFFSZ)
        integer :: iptcl, i, ind_in_stk, nptcls_eff, batchind, ibatch_start, ibatch_end
        integer :: first_stkind, fromp, istk, nptcls_in_stk, last_stkind, ibatch, nbatches, nptcls_tot
        integer(timer_int_kind) :: t_tot, t
        real(timer_int_kind)    :: rt_init, rt_transf, rt_sums
        if( L_BENCH_GLOB )then
            t_tot     = tic()
            reads0    = ptcl_read_counts_now()
            rt_transf = 0.
            rt_sums   = 0.
        endif
        nptcls_tot = 0
        ! fetch data from project
        call b_ptr%spproj%ptr2oritype(p_ptr%oritype, spproj_field)
        ! Initialize temporary arrays
        call cavger_init_online(READBUFFSZ)
        ! Prep for image reading
        call prepimgbatch(p_ptr, b_ptr, READBUFFSZ)
        call reads%new(p_ptr)
        ! Stack & batch loops
        call b_ptr%spproj%map_ptcl_ind2stk_ind(p_ptr%oritype, p_ptr%fromp, first_stkind, ind_in_stk)
        call b_ptr%spproj%map_ptcl_ind2stk_ind(p_ptr%oritype, p_ptr%top,   last_stkind,  ind_in_stk)
        call b_ptr%spproj%get_stkname_and_ind(p_ptr%oritype, p_ptr%fromp, source_stk, ind_in_stk)
        write(logfhandle,'(A,A,A,I0,A,I0,A,A,A,I0)') '>>> MAKE_CAVGS RAW_SOURCE: oritype=', &
            trim(p_ptr%oritype), ' particles=', p_ptr%fromp, '-', &
            p_ptr%top, ' stack=', trim(source_stk%to_char()), ' indstk(first)=', ind_in_stk
        call flush(logfhandle)
        if( L_BENCH_GLOB ) rt_init = toc(t_tot)
        do istk = first_stkind,last_stkind
            ! Particles range in stack
            fromp         = b_ptr%spproj%os_stk%get_fromp(istk)
            nptcls_in_stk = b_ptr%spproj%os_stk%get_top(istk) - fromp + 1   ! # of particles in stack
            ! batch loop
            nbatches = ceiling(real(nptcls_in_stk)/real(READBUFFSZ))
            do ibatch = 1,nbatches
                ibatch_start = (ibatch - 1) * READBUFFSZ + 1                    ! first index in current batch
                ibatch_end   = min(nptcls_in_stk, ibatch_start + READBUFFSZ - 1)! last  index in current batch
                ! identify valid particles
                pinds(:)   = 0                                              ! Global valid particle indices in this batch
                batchind   = 0                                              ! ptcl index in batch
                do i = ibatch_start,ibatch_end
                    iptcl = fromp + i - 1                                   ! Global particle index
                    if( (iptcl < p_ptr%fromp).or.(iptcl > p_ptr%top) ) cycle
                    if( (spproj_field%get_state(iptcl)==0)  )          cycle
                    batchind        = batchind + 1                          ! index in batch
                    pinds(batchind) = iptcl
                enddo
                nptcls_eff = batchind                                       ! # valid particles in batch
                if( nptcls_eff == 0 ) cycle
                nptcls_tot = nptcls_tot + nptcls_eff
                ! Transfer orientation parameters
                if( L_BENCH_GLOB ) t = tic()
                call cavger_transf_oridat( nptcls_eff, pinds(1:nptcls_eff) )
                if( L_BENCH_GLOB ) rt_transf = rt_transf + toc(t)
                ! Read images
                call reads%read(p_ptr, b_ptr, nptcls_eff, pinds(1:nptcls_eff), [1,nptcls_eff])
                ! Interpolate images and update class sums
                if( L_BENCH_GLOB ) t = tic()
                call cavger_update_sums(nptcls_eff, b_ptr%imgbatch(1:nptcls_eff))
                if( L_BENCH_GLOB ) rt_sums = rt_sums + toc(t)
            enddo   ! batch loop
        enddo       ! stack loop
        ! cleanup
        call source_stk%kill
        nullify(spproj_field)
        call reads%kill
        call cavger_dealloc_online
        call killimgbatch(b_ptr)
        if( L_BENCH_GLOB ) call write_assemble_sums_bench(nptcls_tot, ptcl_read_counts_since(reads0), &
            &rt_init, rt_transf, rt_sums, toc(t_tot))
    end subroutine cavger_assemble_sums

    !> The bench file of a part's offline class sums (CAVGSUMS_BENCH_PARTppp.txt)
    subroutine write_assemble_sums_bench( nptcls, reads, rt_init, rt_transf, rt_sums, rt_tot )
        use simple_matcher_ptcl_io, only: ptcl_read_counts, write_ptcl_read_bench_context, write_mem_bench_context
        integer,                intent(in) :: nptcls
        type(ptcl_read_counts), intent(in) :: reads
        real(timer_int_kind),   intent(in) :: rt_init, rt_transf, rt_sums, rt_tot
        type(string) :: fname
        integer :: fnr
        fname = string('CAVGSUMS_BENCH_PART')//int2str_pad(p_ptr%part, max(1, p_ptr%numlen))//TXT_EXT
        call fopen(fnr, FILE=fname, STATUS='REPLACE', action='WRITE')
        write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a,a)')  'cavgsums program                    : ', p_ptr%prg%to_char()
        write(fnr,'(a,i0)') 'cavgsums process partition          : ', p_ptr%part
        write(fnr,'(a,i0)') 'cavgsums worker threads             : ', p_ptr%nthr
        write(fnr,'(a,i0)') 'cavgsums nclasses                   : ', ncls
        write(fnr,'(a,i0)') 'cavgsums particles                  : ', nptcls
        call write_mem_bench_context(fnr, 'cavgsums')
        call write_ptcl_read_bench_context(fnr, 'cavgsums', reads)
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f0.2)') 'cavgsums initialisation             :', rt_init
        write(fnr,'(a,1x,f0.2)') 'cavgsums orientation transfer       :', rt_transf
        write(fnr,'(a,1x,f0.2)') 'cavgsums particle reading           :', reads%rt_read
        write(fnr,'(a,1x,f0.2)') 'cavgsums interpolation and sums     :', rt_sums
        write(fnr,'(a,1x,f0.2)') 'cavgsums total time                 :', rt_tot
        if( rt_tot > 0. )then
            write(fnr,'(a,1x,f0.2)') 'cavgsums % accounted for            :', &
                &((rt_init + rt_transf + reads%rt_read + rt_sums) / rt_tot) * 100.
        endif
        call fclose(fnr)
        call fname%kill
    end subroutine write_assemble_sums_bench

    !>  \brief  merges the even/odd pairs and normalises the sums, merge low resolution
    !    frequencies, calculates & writes FRCs and optionally applies regularization
    module subroutine cavger_restore_cavgs( frcs_fname )
        class(string), intent(in) :: frcs_fname
        real, allocatable :: frcs(:,:)
        type(cavgs_set)   :: cavgs_bak
        type(stack)       :: even_tmp, odd_tmp
        type(image)       :: gridcorr_img
        integer           :: eo_pop(2), icls, ithr, find, pop, filtsz_crop
        logical           :: l_regularize_avg
        ! temporary objects for frc calculation & regularization
        filtsz_crop = fdim(ldim_crop(1))-1
        l_regularize_avg   = p_ptr%l_ml_reg
        allocate(frcs(filtsz_crop,ncls),source=0.0)
        call cavgs_bak%new_set(ldim_crop, ncls)
        call even_tmp%new_stack(ldim_crop, nthr_glob, alloc_ctfsq=.false.)
        call odd_tmp%new_stack( ldim_crop, nthr_glob, alloc_ctfsq=.false.)
        call memoize_ft_maps(ldim_crop(1:2), smpd_crop)
        gridcorr_img = prep2D_inv_instrfun4mul(ldim_crop, ldim_croppd, smpd_crop)
        ! Main loop
        !$omp parallel do default(shared) private(icls,ithr,eo_pop,pop,find)&
        !$omp schedule(static) proc_bind(close)
        do icls = 1,ncls
            ithr   = omp_get_thread_num() + 1
            eo_pop = eo_pops(:,icls)
            pop    = sum(eo_pop)
            if(pop == 0)then
                call cavgs%even%zero_slice(icls, .false.)
                call cavgs%odd%zero_slice(icls, .false.)
                call cavgs%merged%zero_slice(icls, .false.)
            else
                ! even + odd
                cavgs%merged%slices(icls)%ft = .true.
                cavgs%merged%cmat(:,:,icls)  = cavgs%even%cmat(:,:,icls)  + cavgs%odd%cmat(:,:,icls)
                cavgs%merged%ctfsq(:,:,icls) = cavgs%even%ctfsq(:,:,icls) + cavgs%odd%ctfsq(:,:,icls)
                ! backup current classes
                if( l_regularize_avg ) call cavgs_bak%copy_fast(cavgs, icls, .true.)
                ! CTF2 density correction
                if( eo_pop(1) > 1 ) call cavgs%even%ctf_dens_correct(icls)
                if( eo_pop(2) > 1 ) call cavgs%odd%ctf_dens_correct(icls)
                if( pop       > 1 ) call cavgs%merged%ctf_dens_correct(icls)
                ! iFT
                call cavgs%even%ifft(icls)
                call cavgs%odd%ifft(icls)
                call cavgs%merged%ifft(icls)
                ! FRC calculation
                even_tmp%rmat(:,:,ithr)  = cavgs%even%rmat(:,:,icls)
                odd_tmp%rmat(:,:,ithr)   = cavgs%odd%rmat(:,:,icls)
                even_tmp%slices(ithr)%ft = .false.
                odd_tmp%slices(ithr)%ft  = .false.
                call even_tmp%softmask(ithr)
                call odd_tmp%softmask(ithr)
                call even_tmp%fft(ithr)
                call odd_tmp%fft(ithr)
                call even_tmp%frc(odd_tmp, ithr, frcs(:,icls))
                ! ML-regularization: add inverse of noise power to ctfsq & normalize again
                if( l_regularize_avg )then
                    ! add noise power term to denominator
                    call cavgs_bak%even%add_invnoisepower2rho(icls, filtsz_crop, frcs(:,icls))
                    call cavgs_bak%odd%add_invnoisepower2rho(icls, filtsz_crop, frcs(:,icls))
                    if( eo_pop(1) < 3 ) cavgs_bak%even%ctfsq(:,:,icls) = cavgs_bak%even%ctfsq(:,:,icls) + 1.0
                    if( eo_pop(2) < 3 ) cavgs_bak%odd%ctfsq(:,:,icls)  = cavgs_bak%odd%ctfsq(:,:,icls)  + 1.0
                    cavgs_bak%merged%ctfsq(:,:,icls) = cavgs_bak%even%ctfsq(:,:,icls) + cavgs_bak%odd%ctfsq(:,:,icls)
                    ! re-normalize cavg
                    call cavgs_bak%even%ctf_dens_correct(icls)
                    call cavgs_bak%odd%ctf_dens_correct(icls)
                    call cavgs_bak%merged%ctf_dens_correct(icls)
                    ! transfer back cavgs & iFT
                    call cavgs%copy_fast(cavgs_bak, icls, .true.)
                    call cavgs%even%ifft(icls)
                    call cavgs%odd%ifft(icls)
                    call cavgs%merged%ifft(icls)
                endif
                ! average low-resolution info between eo pairs
                find = b_ptr%clsfrcs%estimate_find_for_eoavg(icls, 1)
                call cavgs%merged%fft(icls)
                call cavgs%even%fft(icls)
                call cavgs%odd%fft(icls)
                call cavgs%even%insert_lowres_serial(cavgs%merged, icls, find)
                call cavgs%odd%insert_lowres_serial(cavgs%merged, icls, find)
                ! gridding correction
                call cavgs%merged%ifft(icls)
                call cavgs%even%ifft(icls)
                call cavgs%odd%ifft(icls)
                call gridcorr_img%mul_rmat(cavgs%even%rmat(:ldim_crop(1),:ldim_crop(2),icls:icls))
                call gridcorr_img%mul_rmat(cavgs%odd%rmat(:ldim_crop(1),:ldim_crop(2),icls:icls))
                call gridcorr_img%mul_rmat(cavgs%merged%rmat(:ldim_crop(1),:ldim_crop(2),icls:icls))
            endif
            ! store FRC
            call b_ptr%clsfrcs%set_frc(icls, frcs(:,icls), 1)
        end do
        !$omp end parallel do
        ! write FRCs
        call b_ptr%clsfrcs%write(frcs_fname)
        ! cleanup
        call gridcorr_img%kill
        call forget_ft_maps
        call even_tmp%kill_stack
        call odd_tmp%kill_stack
        call cavgs_bak%kill_set
        deallocate(frcs)
    end subroutine cavger_restore_cavgs

    ! I/O

    module subroutine cavger_write_eo( fname_e, fname_o )
        class(string), intent(in) :: fname_e, fname_o
        call cavgs%even%write(fname_e, .false.)
        call cavgs%odd%write(fname_o, .false.)
    end subroutine cavger_write_eo

    module subroutine cavger_write_all( fname, fname_e, fname_o )
        class(string), intent(in) :: fname, fname_e, fname_o
        call cavger_write_merged( fname)
        call cavger_write_eo( fname_e, fname_o )
    end subroutine cavger_write_all

    module subroutine cavger_write_merged( fname)
        class(string), intent(in) :: fname
        call cavgs%merged%write(fname, .false.)
    end subroutine cavger_write_merged

    module subroutine cavger_read_all()
        if( .not. file_exists(p_ptr%refs) ) THROW_HARD('references (REFS) does not exist in cwd')
        call read_cavgs(p_ptr%refs, 'merged')
        if( file_exists(p_ptr%refs_even) )then
            call read_cavgs(p_ptr%refs_even, 'even')
        else
            call read_cavgs(p_ptr%refs, 'even')
        endif
        if( file_exists(p_ptr%refs_odd) )then
            call read_cavgs(p_ptr%refs_odd, 'odd')
        else
            call read_cavgs(p_ptr%refs, 'odd')
        endif
    end subroutine cavger_read_all

    !>  \brief  submodule utility for reading class averages (image type)
    subroutine read_cavgs( fname, which )
        use simple_imgarr_utils, only: read_stk_into_imgarr
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        class(image), pointer :: pcavgs(:)
        integer               :: ldim_read(3), icls
        ! read
        select case(trim(which))
            case('even')
                cavgs_even = read_stk_into_imgarr(fname)
                pcavgs => cavgs_even
            case('odd')
                cavgs_odd = read_stk_into_imgarr(fname)
                pcavgs => cavgs_odd
            case('merged')
                cavgs_merged = read_stk_into_imgarr(fname)
                pcavgs => cavgs_merged
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
        ldim_read = pcavgs(1)%get_ldim()
        ! scale
        if( any(ldim_read /= ldim_crop) )then
            if( ldim_read(1) > ldim_crop(1) )then
                ! Cropping is not covered
                THROW_HARD('Incompatible cavgs dimensions! ; cavger_read')
            else if( ldim_read(1) < ldim_crop(1) )then
                ! Fourier padding
                !$omp parallel do proc_bind(close) schedule(static) default(shared) private(icls)
                do icls = 1,ncls
                    call pcavgs(icls)%fft
                    call pcavgs(icls)%pad_inplace(ldim_crop)
                    call pcavgs(icls)%ifft
                end do
                !$omp end parallel do
            endif
        endif
        nullify(pcavgs)
    end subroutine read_cavgs

    !>  \brief  writes this worker's current-iteration class sums, accumulated from zero, with
    !!         its class-centering offsets and populations (distributed execution). Workers never
    !!         read carried state; the assembly owner blends (cavger_assemble_sums_from_parts)
    module subroutine cavger_write_contribution( l_frac )
        logical, intent(in) :: l_frac
        type(cavg_sums) :: contrib
        call contrib%new(CAVG_SUMS_CONTRIB, ncls, ldim_crop(1), smpd_crop, part=p_ptr%part)
        call cavgs2sums(contrib)
        call contrib%set_contrib_meta(center_offsets, eo_pops, l_frac)
        call contrib%write(cavg_contrib_fname(p_ptr%part))
        call contrib%kill
    end subroutine cavger_write_contribution

    !>  \brief  shared-memory assembly owner: blends the carried sums into the current ones
    !!         (when l_frac), publishes the new carried set and leaves the blended sums and
    !!         their populations in place for cavger_restore_cavgs
    module subroutine cavger_commit_carryover( l_frac )
        logical, intent(in) :: l_frac
        type(cavg_sums) :: cur
        integer         :: acc_pops(2,ncls)
        call cur%new(CAVG_SUMS_STATE, ncls, ldim_crop(1), smpd_crop)
        call cavgs2sums(cur)
        acc_pops = eo_pops
        call commit_carryover(cur, l_frac, acc_pops)
        call cur%kill
    end subroutine cavger_commit_carryover

    ! Owner blend (doc/refactoring_notes/completed/class_average_and_reconstruct3d_partials_refactoring.md 4.1):
    ! with l_frac and a usable previous set, cur <- s*cur + w*shift(prev), M <- s*n + w*M (N, n on the merged
    ! project); else cur is the new set with M = acc_pops. Publishes it, copies it into cavgs, sets eo_pops to match.
    subroutine commit_carryover( cur, l_frac, acc_pops )
        type(cavg_sums), intent(inout) :: cur
        logical,         intent(in)    :: l_frac
        integer,         intent(in)    :: acc_pops(2,ncls)
        type(cavg_sums)       :: prev
        real,     allocatable :: mprev(:)
        real(dp), allocatable :: mass(:), mass_cur(:)
        integer,  allocatable :: nrep(:), nsmp(:)
        real    :: s(ncls), w(ncls), mnew(ncls)
        integer :: icls, iptcl, status, eo
        logical :: l_blend
        l_blend = .false.
        if( l_frac )then
            call prev%read(string(CAVG_STATE_FILE), CAVG_SUMS_STATE, status)
            if( status == CAVG_SUMS_OK )then
                l_blend = prev%matches(ncls, ldim_crop(1), smpd_crop)
            endif
            if( .not. l_blend ) THROW_WARN('no usable previous class sums; restoring from the current sample only')
        endif
        ! sampling mass of the current sample, per particle, for the carried-mass log
        call cur%class_mass(mass_cur)
        if( l_blend )then
            call prev%get_mrep(mprev)
            call b_ptr%spproj_field%get_group_update_counts('class', ncls, nrep, nsmp)
            do icls = 1, ncls
                if( arg(center_offsets(:,icls)) > CENTHRESH ) call prev%shift_class(icls, center_offsets(:,icls))
            enddo
            call population_blend_weights(nrep, nsmp, mprev, s, w, mnew)
            call cur%blend(prev, s, w)
            ! the restored sums represent the active, updated particles of each class
            eo_pops = 0
            do iptcl = 1, b_ptr%spproj_field%get_noris()
                if( b_ptr%spproj_field%get_state(iptcl) == 0 )     cycle
                if( b_ptr%spproj_field%get_updatecnt(iptcl) <= 0 ) cycle
                icls = b_ptr%spproj_field%get_class(iptcl)
                if( icls < 1 .or. icls > ncls ) cycle
                eo = merge(2, 1, b_ptr%spproj_field%get_eo(iptcl) == 1)
                eo_pops(eo,icls) = eo_pops(eo,icls) + 1
            enddo
        else
            mnew    = real(sum(acc_pops, dim=1))
            eo_pops = acc_pops
        endif
        call cur%set_mrep(mnew)
        call cur%write(string(CAVG_STATE_FILE))
        call sums2cavgs(cur)
        ! one summary line per iteration: counts, weights, and the carried mass per represented
        ! particle relative to the current sample's mass per particle (1 without drift)
        call cur%class_mass(mass)
        if( l_blend )then
            write(logfhandle,'(A,4I9,4F8.4,2ES12.4,F8.4)') '>>> CAVG CARRY-OVER N n MPREV MNEW / W AVG MIN MAX / S / '//&
                &'MASS PER M, PER n / RATIO:', sum(nrep), sum(nsmp), nint(sum(mprev)), nint(sum(mnew)),&
                &sum(w, mask=nrep>0) / real(max(1,count(nrep>0))), minval(w, mask=nrep>0), maxval(w, mask=nrep>0),&
                &sum(s, mask=nsmp>0) / real(max(1,count(nsmp>0))), real(sum(mass)) / max(1., sum(mnew)),&
                &real(sum(mass_cur)) / real(max(1, sum(nsmp))),&
                &(real(sum(mass)) / max(1., sum(mnew))) / max(TINY, real(sum(mass_cur)) / real(max(1, sum(nsmp))))
        else
            write(logfhandle,'(A,I9,ES12.4)') '>>> CAVG CARRY-OVER NONE (CURRENT SAMPLE ONLY) M_NEW / MASS PER M: ', &
                &nint(sum(mnew)), real(sum(mass)) / max(1., sum(mnew))
        endif
        call prev%kill
    end subroutine commit_carryover

    ! copy the module's even/odd accumulators into a cavg_sums container
    subroutine cavgs2sums( sums )
        type(cavg_sums), intent(inout) :: sums
        call sums%set_sums(cavgs%even%cmat, cavgs%odd%cmat, cavgs%even%ctfsq, cavgs%odd%ctfsq)
    end subroutine cavgs2sums

    ! copy a cavg_sums container into the module's even/odd accumulators (Fourier space)
    subroutine sums2cavgs( sums )
        type(cavg_sums), intent(in) :: sums
        call sums%get_sums(cavgs%even%cmat, cavgs%odd%cmat, cavgs%even%ctfsq, cavgs%odd%ctfsq)
        cavgs%even%slices(:)%ft = .true.
        cavgs%odd%slices(:)%ft  = .true.
    end subroutine sums2cavgs

    !>  \brief  Fourier-pads the carried class sums in the current directory once to a larger crop
    !!         box of the same physical extent (the streaming pool's crop-box upsample). Without a
    !!         usable carried set there is nothing to pad: the next iteration is a full update.
    module subroutine cavger_pad_carried_sums( box_crop, smpd_crop )
        integer, intent(in) :: box_crop
        real,    intent(in) :: smpd_crop
        type(cavg_sums) :: carried
        integer         :: status
        call carried%read(string(CAVG_STATE_FILE), CAVG_SUMS_STATE, status)
        if( status == CAVG_SUMS_OK )then
            call carried%pad_to(box_crop, smpd_crop)
            call carried%write(string(CAVG_STATE_FILE))
        endif
        call carried%kill
    end subroutine cavger_pad_carried_sums

    !>  \brief  records the class-centering offset applied to the reference of class icls, so that
    !!         the assembly owner shifts the carried sums once, before blending
    module subroutine cavger_set_center_offset( offset, icls )
        real,    intent(in) :: offset(2)
        integer, intent(in) :: icls
        center_offsets(:,icls) = offset
    end subroutine cavger_set_center_offset

    !>  \brief  distributed assembly owner: sums the workers' current contributions in ascending
    !!         part order, checks that they agree on the class-centering offsets and on the
    !!         carry-over mode, blends and publishes the carried sums, restores the class averages
    !!         and deletes the contributions
    module subroutine cavger_assemble_sums_from_parts
        integer(timer_int_kind) ::  t_init,  t_io,  t_merge_eos_and_norm,  t_tot
        real(timer_int_kind)    :: rt_init, rt_io, rt_merge_eos_and_norm, rt_tot
        type(cavg_sums)      :: total, part_sums
        type(string)         :: benchfname
        real,    allocatable :: offsets(:,:)
        integer, allocatable :: part_pops(:,:)
        real    :: offsets_ref(2,ncls)
        integer :: ipart, fnr, status, acc_pops(2,ncls)
        logical :: l_frac
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call total%new(CAVG_SUMS_STATE, ncls, ldim_crop(1), smpd_crop)
        acc_pops = 0
        l_frac   = .false.
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            t_io    = tic()
        endif
        do ipart = 1, p_ptr%nparts
            call part_sums%read(cavg_contrib_fname(ipart), CAVG_SUMS_CONTRIB, status)
            if( status /= CAVG_SUMS_OK ) THROW_HARD('missing or unreadable class-sum contribution of part '//int2str(ipart))
            if( .not. part_sums%matches(ncls, ldim_crop(1), smpd_crop) )then
                THROW_HARD('class-sum contribution of part '//int2str(ipart)//' does not match the run geometry')
            endif
            call part_sums%get_offsets(offsets)
            if( ipart == 1 )then
                l_frac      = part_sums%get_l_frac()
                offsets_ref = offsets
            else
                if( part_sums%get_l_frac() .neqv. l_frac ) THROW_HARD('workers disagree on the class carry-over mode')
                if( any(abs(offsets - offsets_ref) > 1.e-3) )then
                    THROW_WARN('workers disagree on class-centering offsets; the offsets of part 1 are applied')
                endif
            endif
            call part_sums%get_eo_pops(part_pops)
            acc_pops = acc_pops + part_pops
            call total%accumulate(part_sums)
        enddo
        call part_sums%kill
        center_offsets = offsets_ref
        call commit_carryover(total, l_frac, acc_pops)
        call total%kill
        do ipart = 1, p_ptr%nparts
            call del_file(cavg_contrib_fname(ipart))
        enddo
        if( L_BENCH_GLOB ) rt_io = toc(t_io)
        ! Restoration of e/o/merged classes
        if( L_BENCH_GLOB ) t_merge_eos_and_norm = tic()
        call cavger_restore_cavgs(p_ptr%frcs)
        ! Benchmark
        if( L_BENCH_GLOB )then
            rt_merge_eos_and_norm = toc(t_merge_eos_and_norm)
            rt_tot                = toc(t_tot)
            benchfname = 'CAVGASSEMBLE_BENCH.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'initialisation       :', rt_init
            write(fnr,'(a,1x,f0.2)') 'I/O, sum and blend   :', rt_io
            write(fnr,'(a,1x,f0.2)') 'merge eo-pairs & norm:', rt_merge_eos_and_norm
            write(fnr,'(a,1x,f0.2)') 'total time           :', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f0.2)') 'initialisation        :', (rt_init/rt_tot)               * 100.
            write(fnr,'(a,1x,f0.2)') 'I/O, sum and blend    :', (rt_io/rt_tot)                 * 100.
            write(fnr,'(a,1x,f0.2)') 'merge eo-pairs & norm :', (rt_merge_eos_and_norm/rt_tot) * 100.
            write(fnr,'(a,1x,f0.2)') '% accounted for       :',&
            &((rt_init+rt_io+rt_merge_eos_and_norm)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine cavger_assemble_sums_from_parts

    ! DESTRUCTOR

    !>  \brief  is a destructor
    module subroutine cavger_kill()
        call dealloc_cavgs
        if( allocated(eo_pops)        ) deallocate(eo_pops)
        if( allocated(center_offsets) ) deallocate(center_offsets)
    end subroutine cavger_kill

    !>  \brief submodule private destructor utility
    subroutine dealloc_cavgs
        call dealloc_imgarr(cavgs_even)
        call dealloc_imgarr(cavgs_odd)
        call dealloc_imgarr(cavgs_merged)
        call cavgs%kill_set
        ncls    = 0
        ldim    = 0; ldim_crop   = 0
        ldim_pd = 0; ldim_croppd = 0
        smpd    = 0.;smpd_crop   = 0.
    end subroutine dealloc_cavgs

    ! PUBLIC UTILITIES

    module subroutine transform_ptcls( params, build, spproj, oritype, icls, timgs, pinds, phflip, cavg, imgs_ori, pinds_in, &
        &keep_ft, gridcorr)
        use simple_sp_project,          only: sp_project
        use simple_matcher_ptcl_io,     only: ptcl_read_session, prepimgbatch, killimgbatch
        use simple_memoize_ft_maps
        class(parameters), target,          intent(in)    :: params
        class(builder),                     target, intent(inout) :: build
        class(sp_project),                  intent(inout) :: spproj
        character(len=*),                   intent(in)    :: oritype
        integer,                            intent(in)    :: icls
        type(image),           allocatable, intent(inout) :: timgs(:)
        integer,               allocatable, intent(inout) :: pinds(:)
        logical,     optional,              intent(in)    :: phflip
        type(image), optional,              intent(inout) :: cavg
        type(image), optional, allocatable, intent(inout) :: imgs_ori(:)
        integer,     optional,              intent(in)    :: pinds_in(:)
        logical,     optional,              intent(in)    :: keep_ft   !< return the class-frame Fourier planes (no ifft, no gridding correction)
        type(image), optional,              intent(inout) :: gridcorr  !< the gridding-correction image, for the caller
        class(oris), pointer :: pos
        type(kbinterpol)     :: kbwin
        type(image)          :: img(nthr_glob), gridcorr_img
        type(ptcl_read_session) :: reads
        type(ctfparams)      :: ctfparms
        type(ctf)            :: tfun
        type(string)         :: str
        real,     allocatable :: kbw(:,:)
        integer,  allocatable :: phys_addrh_ori(:,:), phys_addrk_ori(:,:)
        complex :: fcomp, fcompl
        real    :: mat(2,2), shift(2), loc(2), e3, pf2
        integer :: ldim_pd(3), ldim(3),logi_lims(3,2),cyc_lims(3,2),cyc_limsR(2,2),win(2,2)
        integer :: i,iptcl, l,m, pop, h,k, hh,kk,hp,kp, ithr, iwinsz, wdim, physh,physk
        logical :: l_phflip, l_ori_imgs, l_conjg, l_keep_ft
        p_ptr => params
        b_ptr => build
        ! parse inputs
        l_phflip = .false.
        if( present(phflip) ) l_phflip = phflip
        l_keep_ft = .false.
        if( present(keep_ft) ) l_keep_ft = keep_ft
        l_ori_imgs = present(imgs_ori)
        if(present(cavg)) call cavg%kill
        call dealloc_imgarr(timgs)
        if( l_ori_imgs ) call dealloc_imgarr(imgs_ori)
        ! particles indices
        select case(trim(oritype))
            case('ptcl2D')
                str = 'class'
            case('ptcl3D')
                str = 'proj'
            case DEFAULT
                THROW_HARD('ORITYPE not supported!')
        end select
        call spproj%ptr2oritype( oritype, pos )
        if( allocated(pinds) ) deallocate(pinds)
        if( present(pinds_in) )then
            if( size(pinds_in) < 1 ) return
            if( any(pinds_in < 1) .or. any(pinds_in > pos%get_noris()) )then
                THROW_HARD('particle index outside orientation field; transform_ptcls')
            endif
            allocate(pinds, source=pinds_in)
        else
            call pos%get_pinds(icls, str%to_char(), pinds)
        endif
        if( .not.(allocated(pinds)) ) return
        pop = size(pinds)
        if( pop == 0 ) return
        ! Phase flipping sanity check
        if( l_phflip )then
            select case( spproj%get_ctfflag_type(oritype, pinds(1)) )
            case(CTFFLAG_NO)
                THROW_HARD('NO CTF INFORMATION COULD BE FOUND')
            case(CTFFLAG_FLIP)
                THROW_WARN('Images have already been phase-flipped, phase flipping is deactivated')
                l_phflip = .false.
            case(CTFFLAG_YES)
                ! all good
            case DEFAULT
                THROW_HARD('UNSUPPORTED CTF FLAG')
            end select
        endif
        ! interpolation variables
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        allocate(kbw(wdim,wdim),source=0.)
        ! Dimensions and limits
        ldim    = [p_ptr%box,   p_ptr%box,   1]
        ldim_pd = [p_ptr%boxpd, p_ptr%boxpd, 1]
        ! memoization for original size
        call memoize_ft_maps(ldim, p_ptr%smpd)
        phys_addrh_ori = ft_map_phys_addrh
        phys_addrk_ori = ft_map_phys_addrk
        logi_lims      = ft_map_lims
        ! memoization for padded size
        call memoize_ft_maps(ldim_pd, p_ptr%smpd)
        cyc_lims       = ft_map_lims_nr
        cyc_limsR(:,1) = cyc_lims(1,:)
        cyc_limsR(:,2) = cyc_lims(2,:)
        ! Oversampling correction factor
        pf2 = real(OSMPL_PAD_FAC**2)
        ! Gridding correction object
        gridcorr_img = prep2D_inv_instrfun4mul(ldim, ldim_pd, p_ptr%smpd)
        ! transformed output images
        call alloc_imgarr(pop, ldim, p_ptr%smpd, timgs)
        ! temporary objects
        call prepimgbatch(p_ptr, b_ptr, pop)
        !$omp parallel do private(ithr) default(shared) schedule(static) proc_bind(close)
        do ithr = 1, nthr_glob
            call img(ithr)%new(ldim_pd, p_ptr%smpd, wthreads=.false.)
        end do
        !$omp end parallel do
        if( l_ori_imgs ) call alloc_imgarr(pop, ldim, p_ptr%smpd, imgs_ori)
        ! Read all images
        call reads%new(p_ptr)
        call reads%read(p_ptr, b_ptr, pop, pinds(:), [1,pop])
        call reads%kill
        ! Transformation and rotation loop
        !$omp parallel do private(i,ithr,iptcl,shift,e3,ctfparms,tfun,mat,&
        !$omp& h,k,hh,kk,hp,kp,loc,win,l,m,physh,physk,kbw,fcomp,fcompl,l_conjg) &
        !$omp default(shared) schedule(static) proc_bind(close)
        do i = 1,pop
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(i)
            shift = pos%get_2Dshift(iptcl)
            e3    = pos%e3get(iptcl)
            call img(ithr)%zero_and_flag_ft
            call timgs(i)%zero_and_flag_ft
            ! normalisation, padding & forward FT
            call b_ptr%imgbatch(i)%norm_noise_taper_edge_pad_fft(b_ptr%lmsk, img(ithr))
            if( l_ori_imgs )then
                call img(ithr)%ifft
                call img(ithr)%clip(imgs_ori(i))
                call img(ithr)%fft
            endif
            ! optional phase-flipping
            if( l_phflip )then
                ctfparms = spproj%get_ctfparams(oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call img(ithr)%apply_ctf(tfun, 'flip', ctfparms)
            endif
            ! shift
            call img(ithr)%shift2Dserial(-shift)
            ! rotation matrix
            call rotmat2d(-e3, mat)
            do h = logi_lims(1,1),logi_lims(1,2)
                ! padded h-coordinate
                hp = h * OSMPL_PAD_FAC
                do k = logi_lims(2,1),logi_lims(2,2)
                    ! padded k-coordinate
                    kp = k * OSMPL_PAD_FAC
                    ! rotation on the padded lattice
                    loc = matmul(real([hp,kp]),mat)
                    ! interpolation window
                    win(1,:) = nint(loc)
                    win(2,:) = win(1,:) + iwinsz
                    win(1,:) = win(1,:) - iwinsz
                    ! interpolation kernel
                    call kbwin%apod_mat_2d_fast(loc, iwinsz, wdim, kbw)
                    ! interpolation from padded images
                    fcomp = CMPLX_ZERO
                    do l = 1,wdim
                        hh      = win(1,1)+l-1
                        l_conjg = hh < 0
                        hh      = cyci_1d(cyc_limsR(:,1), hh)
                        fcompl  = CMPLX_ZERO
                        do m = 1,wdim
                            kk     = win(1,2)+m-1
                            kk     = cyci_1d(cyc_limsR(:,2), kk)
                            physh  = ft_map_phys_addrh(hh,kk)
                            physk  = ft_map_phys_addrk(hh,kk)
                            fcompl = fcompl + kbw(l,m) * img(ithr)%get_cmat_at(physh,physk,1)
                        enddo
                        fcomp = fcomp + merge(conjg(fcompl), fcompl, l_conjg)
                    end do
                    ! oversampling scaling correction
                    fcomp = pf2 * fcomp
                    ! sets Fourier component
                    physh = phys_addrh_ori(h,k)
                    physk = phys_addrk_ori(h,k)
                    call timgs(i)%set_cmat_at(physh, physk, 1, fcomp)
                enddo
            enddo
            ! backwards FT & gridding correction
            if( .not. l_keep_ft )then
                call timgs(i)%ifft
                call timgs(i)%mul(gridcorr_img)
            endif
        enddo
        !$omp end parallel do
        if( present(cavg) )then
            call cavg%copy(timgs(1))
            do i =2,pop,1
                call cavg%add(timgs(i))
            enddo
            call cavg%div(real(pop))
            if( l_keep_ft )then
                call cavg%ifft
                call cavg%mul(gridcorr_img)
            endif
        endif
        if( present(gridcorr) ) call gridcorr%copy(gridcorr_img)
        ! cleanup
        call killimgbatch(build)
        call forget_ft_maps
        do ithr = 1, nthr_glob
            call img(ithr)%kill
        end do
        call gridcorr_img%kill
        call str%kill
        nullify(pos)
    end subroutine transform_ptcls

end submodule simple_classaverager_restore

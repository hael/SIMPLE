!@descr: Cartesian online/offline 3D reconstruction module
module simple_matcher_3Drec
use simple_core_module_api
use simple_timer
use simple_builder,         only: builder
use simple_cmdline,         only: cmdline
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_kbinterpol,      only: kbinterpol
use simple_math,            only: ceil_div, floor_div
use simple_math_ft,         only: fplane_get_cmplx, fplane_get_ctfsq
use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
use simple_parameters,      only: parameters
use simple_refine3D_fnames, only: refine3D_partial_rec_fbody, refine3D_state_vol_fname
implicit none

public :: init_rec, prep_imgs4rec, update_rec, write_partial_recs, finalize_rec_objs, calc_3Drec, calc_projdir3Drec
private
#include "simple_local_flags.inc"

! Experimental Cartesian reconstruction mode. When enabled, particle Fourier
! planes are inserted with weighted nearest-cell gridding instead of the
! standard full KB splat in reconstructor%insert_plane_oversamp.
logical, parameter :: RECON_USE_GRID_PLANE_NN = .false.

contains

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, allocatable :: pops(:)
        integer :: istate
        call build%spproj_field%get_pops(pops, 'state', maxn=params%nstates)
        do istate = 1, params%nstates
            if( pops(istate) > 0)then
                call build%eorecvols(istate)%new(params, build%spproj)
                call build%eorecvols(istate)%reset_all
            endif
        end do
        deallocate(pops)
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols( params, build )
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer :: istate
        do istate = 1, params%nstates
            call build%eorecvols(istate)%kill
        end do
    end subroutine killrecvols

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( build, fpl, se, o )
        class(builder),     intent(inout) :: build
        class(fplane_type), intent(in)    :: fpl
        class(sym),         intent(inout) :: se
        class(ori),         intent(inout) :: o
        real    :: pw
        integer :: s, eo
        ! state flag
        s = o%get_state()
        if( s == 0 ) return
        ! eo flag
        eo = o%get_eo()
        ! particle-weight
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        if( pw > TINY )then
            if( RECON_USE_GRID_PLANE_NN )then
                call build%eorecvols(s)%grid_plane_nn(se, o, fpl, eo, pw)
            else
                call build%eorecvols(s)%grid_plane(se, o, fpl, eo, pw)
            endif
        endif
    end subroutine grid_ptcl

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( params, build, cline, nptcls, pinds )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable   :: fpls(:)
        integer :: batchlims(2), ibatch, batchsz
        logical :: DEBUG = .false.
        integer(timer_int_kind) :: t, t0
        real(timer_int_kind)    :: t_init, t_read, t_prep, t_grid, t_tot
        if( DEBUG ) t0 = tic()
        ! Initialize objects for recontruction
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        ! Prep batch image objects
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        if( DEBUG ) t_init = toc(t)
        ! gridding batch loop
        if( DEBUG ) then
            t_read = 0.d0
            t_prep = 0.d0
            t_grid = 0.d0
        endif
        do ibatch = 1,nptcls,MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            ! read images
            if( DEBUG ) t = tic()
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            if( DEBUG ) t_read = t_read + toc(t)
            ! preprocess images into padded objects
            if( DEBUG ) t = tic()
            call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz),&
                                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            if( DEBUG ) t_prep = t_prep + toc(t)
            ! insert padded slices into lattice
            if( DEBUG ) t = tic()
            call update_rec(params, build, batchsz, pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            if( DEBUG ) t_grid = t_grid + toc(t)
        end do
        ! Write partial reconstructions and clean up reconstruction objects
        call write_partial_recs(params, build, cline, fpls)
        call finalize_rec_objs(params, build)
        if( DEBUG .and. (params%part==1) )then
            t_tot = toc(t0)
            print *,'Init          : ', t_init
            print *,'Read          : ', t_read
            print *,'Prep          : ', t_prep
            print *,'Grid          : ', t_grid
            print *,'Total rec time: ', t_tot
        endif
    end subroutine calc_3Drec

    !> Volumetric 3D reconstruction from projection-direction pre-averages.
    subroutine calc_projdir3Drec( params, build, cline, nptcls, pinds )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable   :: fpls(:), projdirs(:,:)
        integer,           allocatable   :: eopops(:,:,:), states(:), state_pinds(:)
        type(kbinterpol) :: kbwin
        type(ori) :: orientation
        real      :: w, e3
        integer   :: batchlims(2), batchsz, ibatch, iptcl, iproj, eo, peo, state, pop
        integer   :: i, j, s, iwinsz, wdim
        logical   :: DEBUG = .false.
        integer(timer_int_kind) :: t, t0
        real(timer_int_kind)    :: t_init, t_read, t_prep, t_sum, t_grid, t_tot
        if( DEBUG ) t0 = tic()
        if( nptcls < 1 ) return
        call build%spproj_field%set_projs(build%eulspace)
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz = ceiling(KBWINSZ - 0.5)
        wdim   = kbwin%get_wdim()
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        allocate(projdirs(params%nspace,2), eopops(params%nspace,2,params%nstates), states(nptcls))
        eopops = 0
        states = 0
        !$omp parallel do default(shared) private(i,iptcl,iproj,eo,state,orientation) &
        !$omp schedule(static) proc_bind(close) reduction(+:eopops)
        do i = 1,nptcls
            iptcl = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            if( orientation%isstatezero() ) cycle
            iproj = build%spproj_field%get_int(iptcl, 'proj')
            state = build%spproj_field%get_state(iptcl)
            if( iproj < 1 .or. iproj > params%nspace ) cycle
            if( state < 1 .or. state > params%nstates ) cycle
            eo = build%spproj_field%get_eo(iptcl) + 1
            if( eo < 1 .or. eo > 2 ) cycle
            states(i) = state
            eopops(iproj,eo,state) = eopops(iproj,eo,state) + 1
        end do
        !$omp end parallel do
        if( DEBUG )then
            t_init = toc(t)
            t_read = 0.d0
            t_prep = 0.d0
            t_sum  = 0.d0
            t_grid = 0.d0
        endif
        do s = 1,params%nstates
            if( sum(eopops(:,:,s)) == 0 ) cycle
            if( allocated(state_pinds) ) deallocate(state_pinds)
            state_pinds = pack(pinds, mask=(states == s))
            call zero_projdirs(projdirs)
            do ibatch = 1,size(state_pinds),MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(state_pinds), ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                if( DEBUG ) t = tic()
                call discrete_read_imgbatch(params, build, size(state_pinds), state_pinds, batchlims)
                if( DEBUG ) t_read = t_read + toc(t)
                if( DEBUG ) t = tic()
                call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz), &
                    &state_pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
                if( DEBUG ) t_prep = t_prep + toc(t)
                if( DEBUG ) t = tic()
                !$omp parallel do default(shared) private(j,iproj,eo,i,iptcl,peo,w,e3,pop) &
                !$omp schedule(dynamic) proc_bind(close)
                do j = 1,2*params%nspace
                    if( j <= params%nspace )then
                        iproj = j
                        eo    = 1
                    else
                        iproj = j - params%nspace
                        eo    = 2
                    endif
                    pop = eopops(iproj,eo,s)
                    if( pop == 0 ) cycle
                    do i = batchlims(1),batchlims(2)
                        iptcl = state_pinds(i)
                        if( build%spproj_field%get_int(iptcl, 'proj') /= iproj ) cycle
                        peo = build%spproj_field%get_eo(iptcl) + 1
                        if( peo /= eo ) cycle
                        w = build%spproj_field%get(iptcl, 'w')
                        if( w < TINY ) cycle
                        e3 = build%spproj_field%e3get(iptcl)
                        call accum_projdir_rotated(projdirs(iproj,eo), fpls(i-batchlims(1)+1), e3, w)
                        pop = pop - 1
                        if( pop == 0 ) exit
                    enddo
                enddo
                !$omp end parallel do
                if( DEBUG ) t_sum = t_sum + toc(t)
            enddo
            ! Projection-direction sums are only an acceleration trick: keep
            ! numerator/rho statistics un-restored and let volassemble apply
            ! the ML prior in the same volume-domain step as polar=no.
            if( DEBUG ) t = tic()
            do iproj = 1,params%nspace
                call build%eulspace%get_ori(iproj, orientation)
                call orientation%set_state(s)
                call orientation%set('w', 1.)
                if( eopops(iproj,1,s) > 0 )then
                    call orientation%set('eo', 0)
                    call grid_ptcl(build, projdirs(iproj,1), build%pgrpsyms, orientation)
                endif
                if( eopops(iproj,2,s) > 0 )then
                    call orientation%set('eo', 1)
                    call grid_ptcl(build, projdirs(iproj,2), build%pgrpsyms, orientation)
                endif
            end do
            if( DEBUG ) t_grid = t_grid + toc(t)
        enddo
        call kill_projdirs(projdirs)
        if( allocated(state_pinds) ) deallocate(state_pinds)
        deallocate(projdirs, eopops, states)
        call write_partial_recs(params, build, cline, fpls)
        call finalize_rec_objs(params, build)
        call orientation%kill
        if( DEBUG .and. (params%part==1) )then
            t_tot = toc(t0)
            print *,'Projdir init : ', t_init
            print *,'Projdir read : ', t_read
            print *,'Projdir prep : ', t_prep
            print *,'Projdir sum  : ', t_sum
            print *,'Projdir grid : ', t_grid
            print *,'Projdir total: ', t_tot
        endif

    contains

        subroutine accum_projdir_rotated( projdir, fplane, e3, weight )
            type(fplane_type), intent(inout) :: projdir
            type(fplane_type), intent(in)    :: fplane
            real,              intent(in)    :: e3
            real,              intent(in)    :: weight
            complex :: fcomp
            real    :: loc(2), rmat(2,2), tvalsq
            integer :: fpllims(3,2), h, k, hp, kp, hmax_here, sh
            integer :: nyq_crop, pf
            if( .not. allocated(fplane%cmplx_plane) ) return
            if( .not. allocated(projdir%cmplx_plane) ) call init_projdir(projdir, fplane)
            pf = OSMPL_PAD_FAC
            fpllims = fplane%frlims
            fpllims(1,1) = ceil_div (fplane%frlims(1,1), pf)
            fpllims(1,2) = floor_div(fplane%frlims(1,2), pf)
            fpllims(2,1) = ceil_div (fplane%frlims(2,1), pf)
            fpllims(2,2) = floor_div(fplane%frlims(2,2), pf)
            nyq_crop = fplane%nyq / pf
            call rotmat2D(e3, rmat)
            rmat = real(pf) * rmat
            do k = fpllims(2,1),0
                if( k == 0 )then
                    hmax_here = min(fpllims(1,2), -1)
                else
                    hmax_here = fpllims(1,2)
                endif
                do h = fpllims(1,1),hmax_here
                    sh = nint(hyp(real(h), real(k)))
                    if( sh > nyq_crop ) cycle
                    hp     = h * pf
                    kp     = k * pf
                    fcomp  = fplane_get_cmplx(fplane, hp, kp)
                    tvalsq = fplane_get_ctfsq(fplane, hp, kp)
                    if( abs(real(fcomp)) + abs(aimag(fcomp)) <= TINY .and. tvalsq <= TINY ) cycle
                    ! Keep padded-fplane normalization. Unlike 2D class averages,
                    ! 3D insertion applies the OSMPL_PAD_FAC**2 amplitude factor.
                    loc = matmul(real([h,k]), rmat)
                    call splat_projdir_sample(projdir, loc, fcomp, tvalsq, weight)
                enddo
            enddo
            fcomp  = fplane_get_cmplx(fplane, 0, 0)
            tvalsq = fplane_get_ctfsq(fplane, 0, 0)
            if( abs(real(fcomp)) + abs(aimag(fcomp)) > TINY .or. tvalsq > TINY )then
                call add_projdir_sample(projdir, 0, 0, fcomp, tvalsq, weight)
            endif
        end subroutine accum_projdir_rotated

        subroutine splat_projdir_sample( projdir, loc, fcomp, tvalsq, weight )
            type(fplane_type), intent(inout) :: projdir
            real,              intent(in)    :: loc(2), tvalsq, weight
            complex,           intent(in)    :: fcomp
            real    :: kbw(wdim,wdim), w
            integer :: win(2,2), h, k, hh, kk, i, j
            complex :: val
            call kbwin%apod_mat_2d_fast(loc, iwinsz, wdim, kbw)
            win(1,:) = nint(loc)
            win(2,:) = win(1,:) + iwinsz
            win(1,:) = win(1,:) - iwinsz
            do j = 1,wdim
                kk = win(1,2) + j - 1
                do i = 1,wdim
                    hh = win(1,1) + i - 1
                    w = weight * kbw(i,j)
                    if( kk < 0 )then
                        call add_projdir_sample(projdir, hh, kk, fcomp, tvalsq, w)
                    else if( kk > 0 )then
                        h   = -hh
                        k   = -kk
                        val = conjg(fcomp)
                        call add_projdir_sample(projdir, h, k, val, tvalsq, w)
                    else
                        call add_projdir_sample(projdir, hh, 0, fcomp, tvalsq, w)
                        if( hh /= 0 )then
                            call add_projdir_sample(projdir, -hh, 0, conjg(fcomp), tvalsq, w)
                        endif
                    endif
                enddo
            enddo
        end subroutine splat_projdir_sample

        subroutine add_projdir_sample( projdir, h, k, val, tvalsq, w )
            type(fplane_type), intent(inout) :: projdir
            integer,           intent(in)    :: h, k
            complex,           intent(in)    :: val
            real,              intent(in)    :: tvalsq, w
            if( h < lbound(projdir%cmplx_plane,1) .or. h > ubound(projdir%cmplx_plane,1) ) return
            if( k < lbound(projdir%cmplx_plane,2) .or. k > ubound(projdir%cmplx_plane,2) ) return
            projdir%cmplx_plane(h,k) = projdir%cmplx_plane(h,k) + w * val
            projdir%ctfsq_plane(h,k) = projdir%ctfsq_plane(h,k) + w * tvalsq
        end subroutine add_projdir_sample

        subroutine init_projdir( projdir, fplane )
            type(fplane_type), intent(inout) :: projdir
            type(fplane_type), intent(in)    :: fplane
            integer :: hmin, hmax, kmin, kmax
            hmin = lbound(fplane%cmplx_plane,1)
            hmax = ubound(fplane%cmplx_plane,1)
            kmin = lbound(fplane%cmplx_plane,2)
            kmax = ubound(fplane%cmplx_plane,2)
            projdir%shconst = fplane%shconst
            projdir%frlims  = fplane%frlims
            projdir%nyq     = fplane%nyq
            allocate(projdir%cmplx_plane(hmin:hmax,kmin:kmax), source=cmplx(0.,0.))
            allocate(projdir%ctfsq_plane(hmin:hmax,kmin:kmax), source=0.)
        end subroutine init_projdir

        subroutine zero_projdirs( planes )
            type(fplane_type), intent(inout) :: planes(:,:)
            integer :: i, j
            do j = 1,size(planes,2)
                do i = 1,size(planes,1)
                    if( allocated(planes(i,j)%cmplx_plane) ) planes(i,j)%cmplx_plane = cmplx(0.,0.)
                    if( allocated(planes(i,j)%ctfsq_plane) ) planes(i,j)%ctfsq_plane = 0.
                enddo
            enddo
        end subroutine zero_projdirs

        subroutine kill_projdirs( planes )
            type(fplane_type), intent(inout) :: planes(:,:)
            integer :: i, j
            do j = 1,size(planes,2)
                do i = 1,size(planes,1)
                    if( allocated(planes(i,j)%cmplx_plane) ) deallocate(planes(i,j)%cmplx_plane)
                    if( allocated(planes(i,j)%ctfsq_plane) ) deallocate(planes(i,j)%ctfsq_plane)
                enddo
            enddo
        end subroutine kill_projdirs

    end subroutine calc_projdir3Drec

    !>  Initiates objects required for online volumetric 3d reconstruction
    !>  Does not read images
    subroutine init_rec( params, build, maxbatchsz, fplanes, init_volumes )
        use simple_imgarr_utils, only: alloc_imgarr
        class(parameters),              intent(in)    :: params
        class(builder),                 intent(inout) :: build
        integer,                        intent(in)    :: maxbatchsz
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        logical, optional,              intent(in)    :: init_volumes
        logical :: l_init_volumes
        l_init_volumes = .true.
        if( present(init_volumes) ) l_init_volumes = init_volumes
        ! sanity check for ml_reg
        if( params%l_ml_reg )then
            if( .not. allocated(build%esig%sigma2_noise) )then
                THROW_HARD('build%esig%sigma2_noise is not allocated while ml_reg is enabled; calc_3Drec')
            endif
        endif
        ! init volumes
        if( l_init_volumes ) call preprecvols(params, build)
        ! allocate convenience CTF & memory aligned objects
        if( allocated(fplanes) )  deallocate(fplanes)
        allocate(fplanes(maxbatchsz))
        ! heap of padded images
        call alloc_imgarr(nthr_glob, [params%boxpd, params%boxpd, 1], params%smpd, build%img_pad_heap)
    end subroutine init_rec

    !> Preprocess particle images for online volumetric 3d reconstruction
    subroutine prep_imgs4rec( params, build, nptcls, ptcl_imgs, pinds, fplanes )
        use simple_image, only: image
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ctfparams) :: ctfparms(nthr_glob)
        real      :: shift(2)
        integer   :: iptcl, i, ithr, kfromto(2)
        ! logical/physical address mapping for padded Fourier planes
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        ! gridding batch loop
        kfromto = build%esig%get_kfromto()
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            if( params%l_ml_reg )then
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr),&
                &shift, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl))
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr),&
                &shift, fplanes(i))
            endif
        end do
        !$omp end parallel do
    end subroutine prep_imgs4rec

    !> volumetric 3d reconstruction
    subroutine update_rec( params, build, nptcls, pinds, fplanes )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ori) :: orientation
        integer   :: iptcl, i
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        ! gridding
        do i = 1,nptcls
            iptcl  = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            if( orientation%isstatezero() ) cycle
            call grid_ptcl(build, fplanes(i), build%pgrpsyms, orientation)
        end do
        call orientation%kill
    end subroutine update_rec

    !> volumetric 3d reconstruction
    subroutine write_partial_recs( params, build, cline, fplanes )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters),              intent(inout) :: params
        class(builder),                 intent(inout) :: build
        class(cmdline),                 intent(inout) :: cline
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        integer :: i, s, numlen_part
        ! deallocate convenience objects
        do i = 1,size(fplanes)
            if( allocated(fplanes(i)%cmplx_plane) ) deallocate(fplanes(i)%cmplx_plane)
            if( allocated(fplanes(i)%ctfsq_plane) ) deallocate(fplanes(i)%ctfsq_plane)
        end do
        deallocate(fplanes)
        call dealloc_imgarr(build%img_pad_heap)
        ! write partial reconstructions for downstream volassemble
        numlen_part = max(1, params%numlen)
        do s=1,params%nstates
            if( build%spproj_field%get_pop(s, 'state') == 0 )then
                build%fsc(s,:) = 0.
                cycle
            endif
            call build%eorecvols(s)%compress_exp
            call build%eorecvols(s)%write_eos(refine3D_partial_rec_fbody(s, params%part, numlen_part))
            if( .not. cline%defined('force_volassemble') )then
                params%vols(s) = refine3D_state_vol_fname(s)
                call cline%set('vol'//int2str(s), params%vols(s))
            endif
        end do
    end subroutine write_partial_recs

    subroutine finalize_rec_objs( params, build )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        call dealloc_imgarr(build%img_pad_heap)
        call forget_ft_maps
        call killrecvols(params, build)
    end subroutine finalize_rec_objs

end module simple_matcher_3Drec

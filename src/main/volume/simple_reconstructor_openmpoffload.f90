!@descr: provides one routine for gpu-accelerated reconstruction
module simple_reconstructor_openmpoffload
use simple_core_module_api
use simple_builder,          only: builder
use simple_parameters,       only: parameters
use simple_reconstructor_eo, only: reconstructor_eo
use simple_matcher_ptcl_io,  only: prepimgbatch, discrete_read_imgbatch
use simple_matcher_3Drec,    only: prep_imgs4rec, update_rec, init_rec, write_partial_recs, finalize_rec_objs
use simple_cmdline,          only: cmdline
use simple_math,             only: ceil_div, floor_div
use simple_kbinterpol,       only: apod_kb15_a2
implicit none

public :: calc_3Drec_gpu
private
#include "simple_local_flags.inc"

logical              :: DEBUG = .true.
real(timer_int_kind) :: t_init, t_read, t_prep, t_grid, t_h2d, t_launch, t_overlap, t_wait, t_finalize, t_tot

contains

    !> volumetric 3D reconstruction
    subroutine calc_3Drec_gpu( params, build, cline, nptcls, pinds )
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(cmdline),          intent(inout) :: cline
        integer,                 intent(in)    :: nptcls
        integer,                 intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable :: fpls(:)
        real,              allocatable, target :: symmats(:,:,:), rotmats(:,:,:)
        integer   :: vollims(3,2)
        integer   :: cdim(3), clb(3), jsym, nsym, h_edge, nyq
        integer(timer_int_kind) :: t, t0
#ifndef USE_OPENMP_OFFLOAD
        THROW_HARD('calc_3Drec_gpu is part of the GPU path. Use calc_3Drec instead')
#else
        if( DEBUG ) t0 = tic()
        ! Initialize objects for recontruction
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        ! Prep batch image objects
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        ! 3D limits
        vollims = build%eorecvols(1)%even%loop_lims(2)
        h_edge  = vollims(1,1)
        ! Setup rotation matrices
        allocate(symmats(3,3,build%pgrpsyms%get_nsym()))
        nsym = build%pgrpsyms%get_nsym()
        do jsym = 1, nsym
            call build%pgrpsyms%get_sym_rmat(jsym, symmats(:,:,jsym))
        end do
        allocate(rotmats(2,3,MAXIMGBATCHSZ))
        ! Arrays
        clb  = lbound(build%eorecvols(1)%even%cmat_exp)
        cdim = ubound(build%eorecvols(1)%even%cmat_exp) - clb + 1
        nyq  = build%eorecvols(1)%even%get_lfny(1)
        if( DEBUG ) then
            t_init     = toc(t)
            t_read     = 0.d0
            t_prep     = 0.d0
            t_grid     = 0.d0
            t_h2d      = 0.d0
            t_launch   = 0.d0
            t_overlap  = 0.d0
            t_wait     = 0.d0
            t_finalize = 0.d0
        endif
        ! Gridding interpolation of all particles
        call insert_all_slices(params, build, nptcls, pinds, fpls,&
            & build%eorecvols(1)%even%cmat_exp, build%eorecvols(1)%odd%cmat_exp,&
            & build%eorecvols(1)%even%rho_exp, build%eorecvols(1)%odd%rho_exp,&
            & cdim, clb, h_edge, nyq, symmats, rotmats, nsym)
        if( DEBUG ) t = tic()
        deallocate(symmats, rotmats)
        ! Write partial reconstructions
        call write_partial_recs(params, build, cline, fpls)
        ! Clean up
        call finalize_rec_objs(params, build)
        ! Timings
        if( DEBUG .and. (params%part==1) )then
            t_finalize = toc(t)
            t_tot      = toc(t0)
            print *,'Init          : ', t_init
            print *,'Read          : ', t_read
            print *,'Prep          : ', t_prep
            print *,'Grid          : ', t_grid
            print *,'H2D/map       : ', t_h2d
            print *,'GPU launch    : ', t_launch
            print *,'CPU overlap   : ', t_overlap
            print *,'GPU wait      : ', t_wait
            print *,'Finalize      : ', t_finalize
            print *,'Total rec time: ', t_tot
        endif
#endif
    end subroutine calc_3Drec_gpu

#ifdef USE_OPENMP_OFFLOAD

    subroutine insert_all_slices(params, build, nptcls, pinds, fpls,&
            &cmatexp_e, cmatexp_o, rhoexp_e, rhoexp_o, cdim, clb,&
            &h_edge, nyq, symmats, rotmats, nsym)
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls, cdim(3), clb(3), h_edge, nyq, nsym
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fpls(MAXIMGBATCHSZ)
        complex,   target, intent(inout) :: cmatexp_e(cdim(1),cdim(2),cdim(3))
        complex,   target, intent(inout) :: cmatexp_o(cdim(1),cdim(2),cdim(3))
        real,      target, intent(inout) :: rhoexp_e(cdim(1),cdim(2),cdim(3))
        real,      target, intent(inout) :: rhoexp_o(cdim(1),cdim(2),cdim(3))
        real,      target, intent(in)    :: symmats(3,3,nsym)
        real,      target, intent(inout) :: rotmats(2,3,MAXIMGBATCHSZ)
        complex,    allocatable :: fplanes(:,:,:)
        real,       allocatable :: ctfsqplanes(:,:,:)
        logical,    allocatable :: even(:)
        integer(timer_int_kind) :: t, t_stage
        real(timer_int_kind)    :: dt_h2d, dt_launch, dt_wait
        integer :: fpllims(3,2), fpllims_pd(3,2), cdim2D(2), clb2D(2), batchlims(2)
        integer :: ibatch, batchsz, sz, nbatch
        ! prep first batch
        nbatch = ceiling(real(nptcls)/real(MAXIMGBATCHSZ))
        call prep_batch(1, batchsz, batchlims)
        if( DEBUG ) t = tic()
        fpllims_pd   = fpls(1)%frlims
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), OSMPL_PAD_FAC)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), OSMPL_PAD_FAC)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), OSMPL_PAD_FAC)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), OSMPL_PAD_FAC)
        clb2D  = lbound(fpls(1)%cmplx_plane)
        cdim2D = ubound(fpls(1)%cmplx_plane) - clb2D + 1
        allocate(fplanes(cdim2D(1),cdim2D(2),MAXIMGBATCHSZ),&
            &ctfsqplanes(cdim2D(1),cdim2D(2),MAXIMGBATCHSZ),&
            &even(MAXIMGBATCHSZ))
        if( DEBUG ) t_prep = t_prep + toc(t)

        call update_ptcls_arrays(batchsz, batchlims(1), build%spproj_field, fpls,&
        & nptcls, pinds, even, rotmats, fplanes, ctfsqplanes)

        ! Device storage
        if( DEBUG ) t = tic()
        !$omp target enter data map(alloc: ctfsqplanes(1:cdim2D(1),1:cdim2D(2),1:MAXIMGBATCHSZ),&
        !$omp& fplanes(1:cdim2D(1),1:cdim2D(2),1:MAXIMGBATCHSZ),&
        !$omp& rotmats(1:2,1:3,1:MAXIMGBATCHSZ), even(1:MAXIMGBATCHSZ))
        !$omp target enter data map(to: symmats(1:3,1:3,1:nsym),&
        !$omp& rhoexp_e(1:cdim(1),1:cdim(2),1:cdim(3)), rhoexp_o(1:cdim(1),1:cdim(2),1:cdim(3)),&
        !$omp& cmatexp_e(1:cdim(1),1:cdim(2),1:cdim(3)), cmatexp_o(1:cdim(1),1:cdim(2),1:cdim(3)))
        if( DEBUG ) t_h2d = t_h2d + toc(t)
        do ibatch = 1,nbatch
            ! prep batch meta-data & matrices
            ! call update_ptcls_arrays(batchsz, batchlims(1), build%spproj_field, fpls,&
            ! & nptcls, pinds, even, rotmats, fplanes, ctfsqplanes)
            if( DEBUG )then
                dt_h2d    = 0.d0
                dt_launch = 0.d0
                dt_wait   = 0.d0
                t_stage   = tic()
            endif
            sz = batchsz    ! For use on device
            !$omp target update to(fplanes(1:cdim2D(1),1:cdim2D(2),1:sz),even(1:sz),&
            !$omp& ctfsqplanes(1:cdim2D(1),1:cdim2D(2),1:sz), rotmats(1:2,1:3,1:sz))
            if( DEBUG )then
                dt_h2d = toc(t_stage)
                t_h2d  = t_h2d + dt_h2d
                t_stage = tic()
            endif
            !$omp target data use_device_addr(even, rotmats, symmats, fplanes,&
            !$omp& ctfsqplanes, cmatexp_e, cmatexp_o, rhoexp_e, rhoexp_o)
            call insert_slices_batch(sz, cdim, cdim2D, clb2D, clb, fpllims, nsym,&
                &h_edge, nyq, even, rotmats, symmats, fplanes, ctfsqplanes,&
                &cmatexp_e, cmatexp_o, rhoexp_e, rhoexp_o)
            !$omp end target data
            if( DEBUG )then
                dt_launch = toc(t_stage)
                t_launch  = t_launch + dt_launch
                t_stage   = tic()
            endif
            if( ibatch < nbatch )then
                ! process images for next batch while device is busy
                call prep_batch(ibatch+1, batchsz, batchlims)
                call update_ptcls_arrays(batchsz, batchlims(1), build%spproj_field, fpls,&
                    &nptcls, pinds, even, rotmats, fplanes, ctfsqplanes)
            endif
            if( DEBUG )then
                t_overlap = t_overlap + toc(t_stage)
                t_stage   = tic()
            endif
            !$omp taskwait
            if( DEBUG )then
                dt_wait = toc(t_stage)
                t_wait  = t_wait + dt_wait
                t_grid  = t_grid + dt_h2d + dt_launch + dt_wait
            endif
        enddo
        if( DEBUG ) t = tic()
        !$omp target exit data map(delete:fplanes(1:cdim2D(1),1:cdim2D(2),1:MAXIMGBATCHSZ),&
        !$omp& ctfsqplanes(1:cdim2D(1),1:cdim2D(2),1:MAXIMGBATCHSZ),&
        !$omp& rotmats(1:2,1:3,1:MAXIMGBATCHSZ), even(1:MAXIMGBATCHSZ))
        !$omp target exit data map(from: rhoexp_e(1:cdim(1),1:cdim(2),1:cdim(3)),&
        !$omp& rhoexp_o(1:cdim(1),1:cdim(2),1:cdim(3)), cmatexp_e(1:cdim(1),1:cdim(2),1:cdim(3)),&
        !$omp& cmatexp_o(1:cdim(1),1:cdim(2),1:cdim(3))) map(release: symmats(1:3,1:3,1:nsym))
        deallocate(fplanes, ctfsqplanes, even)
        if( DEBUG ) t_finalize = t_finalize + toc(t)
        contains

            subroutine prep_batch( batchid, sz, lims )
                integer, intent(in) :: batchid
                integer, intent(out):: sz, lims(2)
                integer(timer_int_kind) :: t_local
                lims(1) = (batchid-1)*MAXIMGBATCHSZ + 1
                lims(2) = min(nptcls, batchid*MAXIMGBATCHSZ)
                sz      = lims(2) - lims(1) + 1
                ! read images
                if( DEBUG ) t_local = tic()
                call discrete_read_imgbatch(params, build, nptcls, pinds, lims)
                if( DEBUG ) t_read = t_read + toc(t_local)
                ! preprocess images into padded objects
                if( DEBUG ) t_local = tic()
                call prep_imgs4rec(params, build, sz, build%imgbatch(:sz), pinds(lims(1):lims(2)), fpls(:sz))
                if( DEBUG ) t_prep = t_prep + toc(t_local)
            end subroutine prep_batch

    end subroutine insert_all_slices

    subroutine insert_slices_batch(sz, cdim, cdim2D, clb2D, clb, fpllims, nsym,&
        &h_edge, nyq, even, rotmats, symmats, fplanes, ctfsqplanes,&
        &cmatexp_e, cmatexp_o, rhoexp_e, rhoexp_o)
        integer,         intent(in)    :: sz, cdim(3), cdim2D(2), clb2D(2), clb(3), fpllims(3,2)
        integer,         intent(in)    :: nsym, h_edge, nyq
        logical,         intent(in)    :: even(MAXIMGBATCHSZ)
        real,            intent(in)    :: rotmats(2,3,MAXIMGBATCHSZ)
        real,            intent(in)    :: symmats(3,3,nsym)
        complex,         intent(in)    :: fplanes(cdim2D(1),cdim2D(2),MAXIMGBATCHSZ)
        real,            intent(in)    :: ctfsqplanes(cdim2D(1),cdim2D(2),MAXIMGBATCHSZ)
        complex, target, intent(inout) :: cmatexp_e(cdim(1),cdim(2),cdim(3))
        complex, target, intent(inout) :: cmatexp_o(cdim(1),cdim(2),cdim(3))
        real,    target, intent(inout) :: rhoexp_e(cdim(1),cdim(2),cdim(3))
        real,    target, intent(inout) :: rhoexp_o(cdim(1),cdim(2),cdim(3))
        integer, parameter :: WDIM = 3
        real,      pointer :: rho(:,:,:)
        complex,   pointer :: cmat(:,:,:)
        complex :: comp
        real    :: loc(3), base(3), wx(WDIM), wy(WDIM), wz(WDIM), w_ctfsq
        real    :: pf2, r11, r12, r13, r21, r22, r23, sumx, sumy, sumz, ctfsq, wyz
        integer :: win(3,2), nyqsq, iwinsz, h, i, k, l, m, isym, iy, iz, ky, mz
        integer :: h_sq, k_max_h, k_lo, k_hi, hp, kp, hpb, kpb
        iwinsz = ceiling(KBWINSZ - 0.5)
        nyqsq  = nyq * (nyq + 1)
        pf2    = real(OSMPL_PAD_FAC**2)
        !$omp target teams distribute parallel do nowait collapse(2) proc_bind(close)&
        !$omp& map(to: sz, cdim(1:3), cdim2D(1:2), clb2D(1:2), clb(1:3),&
        !$omp& fpllims(1:3,1:2), nsym, h_edge, nyqsq, pf2, iwinsz)&
        !$omp& has_device_addr(even, rotmats, symmats, fplanes,&
        !$omp& ctfsqplanes, cmatexp_e, cmatexp_o, rhoexp_e, rhoexp_o)&
        !$omp& default(shared) private(h,i,k,l,m,isym,r11,r12,r13,r21,r22,r23,ctfsq,&
        !$omp& win,sumx,sumy,sumz,ky,mz,wx,wy,wz,wyz,base,hp,kp,hpb,kpb,comp,&
        !$omp& loc,w_ctfsq,h_sq,k_max_h,k_lo,k_hi,cmat,rho)
        do i = 1, sz
            do isym = 1, nsym
                if( even(i) )then
                    cmat => cmatexp_e
                    rho  => rhoexp_e
                else
                    cmat => cmatexp_o
                    rho  => rhoexp_o
                endif
                r11 = rotmats(1,1,i)*symmats(1,1,isym) + rotmats(1,2,i)*symmats(2,1,isym) + rotmats(1,3,i)*symmats(3,1,isym)
                r12 = rotmats(1,1,i)*symmats(1,2,isym) + rotmats(1,2,i)*symmats(2,2,isym) + rotmats(1,3,i)*symmats(3,2,isym)
                r13 = rotmats(1,1,i)*symmats(1,3,isym) + rotmats(1,2,i)*symmats(2,3,isym) + rotmats(1,3,i)*symmats(3,3,isym)
                r21 = rotmats(2,1,i)*symmats(1,1,isym) + rotmats(2,2,i)*symmats(2,1,isym) + rotmats(2,3,i)*symmats(3,1,isym)
                r22 = rotmats(2,1,i)*symmats(1,2,isym) + rotmats(2,2,i)*symmats(2,2,isym) + rotmats(2,3,i)*symmats(3,2,isym)
                r23 = rotmats(2,1,i)*symmats(1,3,isym) + rotmats(2,2,i)*symmats(2,3,isym) + rotmats(2,3,i)*symmats(3,3,isym)
                do m = 0, WDIM-1
                    do h = fpllims(1,1)+m, fpllims(1,2), WDIM
                        h_sq = h*h
                        if( h_sq > nyqsq ) cycle
                        k_max_h = int(sqrt(real(nyqsq - h_sq)))
                        k_lo    = max(fpllims(2,1), -k_max_h)
                        k_hi    = min(fpllims(2,2),  k_max_h)
                        loc(1)  = real(h) * r11 + real(k_lo-1) * r21
                        loc(2)  = real(h) * r12 + real(k_lo-1) * r22
                        loc(3)  = real(h) * r13 + real(k_lo-1) * r23
                        hp      = h * OSMPL_PAD_FAC
                        do k = k_lo, k_hi
                            loc(1) = loc(1) + r21
                            loc(2) = loc(2) + r22
                            loc(3) = loc(3) + r23
                            win(:,1) = nint(loc)
                            win(:,2) = win(:,1) + iwinsz
                            win(:,1) = win(:,1) - iwinsz
                            if( win(1,2) < h_edge ) cycle
                            kp = k * OSMPL_PAD_FAC
                            if( kp <= 0 )then
                                hpb   = hp - clb2D(1) + 1
                                kpb   = kp - clb2D(2) + 1
                                comp  = fplanes(hpb,kpb,i)
                                ctfsq = ctfsqplanes(hpb,kpb,i)
                            else
                                hpb   = -hp - clb2D(1) + 1
                                kpb   = -kp - clb2D(2) + 1
                                comp  = conjg(fplanes(hpb,kpb,i))
                                ctfsq = ctfsqplanes(hpb,kpb,i)
                            endif
                            if( abs(real(comp)) + abs(aimag(comp)) <= TINY .and. ctfsq <= TINY ) cycle
                            comp = pf2 * comp
                            base = real(win(:,1)) - loc
                            sumx = 0.0; sumy = 0.0; sumz = 0.0
                            do l = 1, WDIM
                                wx(l) = apod_kb15_a2(base(1) + real(l-1))
                                wy(l) = apod_kb15_a2(base(2) + real(l-1))
                                wz(l) = apod_kb15_a2(base(3) + real(l-1))
                                sumx = sumx + wx(l)
                                sumy = sumy + wy(l)
                                sumz = sumz + wz(l)
                            end do
                            wx = wx / sumx; wy = wy / sumy; wz = wz / sumz
                            win(:,1) = win(:,1) - clb
                            do iz = 1, WDIM
                                mz = win(3,1) + iz
                                do iy = 1, WDIM
                                    ky      = win(2,1) + iy
                                    wyz     = wy(iy) * wz(iz)
                                    w_ctfsq = wyz * ctfsq
                                    cmat(win(1,1)+1, ky, mz) = cmat(win(1,1)+1, ky, mz) + wx(1) * wyz * comp
                                    cmat(win(1,1)+2, ky, mz) = cmat(win(1,1)+2, ky, mz) + wx(2) * wyz * comp
                                    cmat(win(1,1)+3, ky, mz) = cmat(win(1,1)+3, ky, mz) + wx(3) * wyz * comp
                                    rho( win(1,1)+1, ky, mz) = rho( win(1,1)+1, ky, mz) + wx(1) * w_ctfsq
                                    rho( win(1,1)+2, ky, mz) = rho( win(1,1)+2, ky, mz) + wx(2) * w_ctfsq
                                    rho( win(1,1)+3, ky, mz) = rho( win(1,1)+3, ky, mz) + wx(3) * w_ctfsq
                                end do
                            end do
                        end do
                    end do
                end do
            enddo
        end do
        !$omp end target teams distribute parallel do
    end subroutine insert_slices_batch

    subroutine update_ptcls_arrays( bsz, bstart, os, fpls, nptcls, pinds, &
        & even, rotmats, fplanes, ctfsqplanes )
        integer,            intent(in)  :: bsz, bstart, nptcls
        class(oris),        intent(in)  :: os
        class(fplane_type), intent(in)  :: fpls(:)
        integer,            intent(in)  :: pinds(nptcls)
        logical,            intent(out) :: even(MAXIMGBATCHSZ)
        real,               intent(out) :: rotmats(2,3,MAXIMGBATCHSZ)
        complex,            intent(out) :: fplanes(:,:,:)
        real,               intent(out) :: ctfsqplanes(:,:,:)
        integer(timer_int_kind) :: t
        real    :: rmat(3,3)
        integer :: i,j,iptcl
        if( DEBUG ) t = tic()
        !$omp parallel do proc_bind(close) private(i,j,iptcl,rmat) default(shared) schedule(static)
        do i = 1, bsz
            j       = bstart + i - 1
            iptcl   = pinds(j)
            even(i) = os%get_eo(iptcl) == 0
            rmat               = os%get_mat(iptcl)
            rotmats(:,:,i)     = rmat(1:2,:)
            fplanes(:,:,i)     = fpls(i)%cmplx_plane(:,:)
            ctfsqplanes(:,:,i) = fpls(i)%ctfsq_plane(:,:)
        end do
        !$omp end parallel do
        if( DEBUG ) t_prep = t_prep + toc(t)
    end subroutine update_ptcls_arrays

#endif
end module simple_reconstructor_openmpoffload

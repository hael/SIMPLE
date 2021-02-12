module simple_motion_refine_nano
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,                         only: image
use simple_polarizer,                     only: polarizer
use simple_oris,                          only: oris
use simple_polarft_corrcalc
use simple_pftcc_shsrch_grad,             only: pftcc_shsrch_grad
use simple_parameters,                    only: params_glob
use simple_tseries_graphene_subtr
implicit none
! private
#include "simple_local_flags.inc"

public :: motion_refine_nano
!

logical, parameter :: DEBUG   = .true.
real,    parameter :: HALFROT = 12.
integer, parameter :: MAXNFRAMESPERPTCL = 50

type image_ptr
    real,    public, pointer :: rmat(:,:,:)
    complex, public, pointer :: cmat(:,:,:)
end type image_ptr

type :: motion_refine_nano
    class(oris),  pointer                :: pptcl2D
    class(oris),  pointer                :: pcls2D
    class(oris),  pointer                :: pframes
    type(image),             allocatable :: particles(:), padded_imgs(:), rot_imgs(:)
    type(pftcc_shsrch_grad), allocatable :: grad_shsrch_objs(:)
    type(polarizer),         allocatable :: ptcl_polarizers(:)
    type(image),             allocatable :: background_pspecs(:)
    character(LONGSTRLEN),   allocatable :: frames(:)
    character(len=:),        allocatable :: refs
    type(polarizer)                      :: reference
    type(image)                          :: frame
    type(polarft_corrcalc)               :: pftcc
    real,                    allocatable :: shifts_sub(:,:)       ! sub-pixel after re-extraction: ori = discrete + sub
    real,                    allocatable :: shifts(:,:), e3(:)
    logical,                 allocatable :: lresmsk(:)
    real                                 :: msk
    integer                              :: ldim_frame(3), ldim_ptcl(3)
    integer                              :: nframes, nptcls,ncls,nrots
    logical                              :: exists = .false.

    contains
        procedure          :: new
        procedure          :: correct
        procedure, private :: rotate_and_add
        procedure, private :: kernel_correction
        procedure          :: kill
end type motion_refine_nano

contains

    subroutine new( self, os_ptcl2D, os_cls2D, os_frames, refs )
        class(motion_refine_nano), intent(inout) :: self
        class(oris), target,       intent(inout) :: os_ptcl2D, os_cls2D, os_frames
        character(len=*),          intent(in)    :: refs
        integer :: i
        call self%kill
        self%pptcl2D => os_ptcl2D
        self%pcls2D  => os_cls2D
        self%pframes => os_frames
        self%refs    = trim(refs)
        self%nptcls  = self%pptcl2D%get_noris()
        self%nframes = self%pframes%get_noris()
        self%ncls    = self%pcls2D%get_noris()
        allocate(self%frames(self%nframes))
        do i = 1,self%nframes
            self%frames(i) = trim(self%pframes%get_static(i,'frame'))
        enddo
        self%ldim_frame = [nint(self%pframes%get(1,'xdim')), nint(self%pframes%get(1,'ydim')), 1]
        self%ldim_ptcl  = [params_glob%box, params_glob%box, 1]
        call self%frame%new(self%ldim_frame, params_glob%smpd)
        call self%reference%new(self%ldim_ptcl, params_glob%smpd)
        self%msk = real(params_glob%box/2.)-COSMSKHALFWIDTH-5.
        ! resolution limits
        params_glob%kfromto(1) = calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd)
        params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
        params_glob%kstop      = params_glob%kfromto(2)
        params_glob%boxpd      = 2*params_glob%box
        self%lresmsk = calc_graphene_mask(params_glob%box, params_glob%smpd)
        ! convenience objects
        allocate(self%ptcl_polarizers(nthr_glob),self%grad_shsrch_objs(nthr_glob),self%rot_imgs(nthr_glob))
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1,nthr_glob
            call self%rot_imgs(i)%new([params_glob%boxpd,params_glob%boxpd,1],params_glob%smpd,wthreads=.false.)
            call self%ptcl_polarizers(i)%new(self%ldim_ptcl,params_glob%smpd,wthreads=.false.)
        enddo
        !$omp end parallel do
    end subroutine new

    subroutine correct( self, new_stk_fname )
        class(motion_refine_nano), intent(inout) :: self
        character(len=*),          intent(in)    :: new_stk_fname
        type(image)              :: tmpimg, ave, neighbour_imgs(9), tmpimgs(9)
        real,    allocatable :: corrs(:), corrmat(:,:), weights(:),ranked_weights(:)
        logical, allocatable :: lmsk(:,:,:), ptcl_mask(:), ptcl_cls_mask(:)
        integer, allocatable :: pinds(:), ptcl_pos(:,:)
        real    :: shift(2), lims(2,2), lims_init(2,2), ang, sdev_noise, angstep, cxy(3), graphene_pos(2),threshold,wsum
        integer :: fromto(2),i, j, k, icls, iptcl,pop,noutside, irot, hn,npop,ithr, ix,iy, cnt
        integer :: iring2, ikfromto(2), ikstop, inthr_glob
        logical :: valid_neighbour(9)
        inthr_glob = nthr_glob
        iring2     = params_glob%ring2
        ikfromto   = params_glob%kfromto
        ikstop     = params_glob%kstop
        ! shift search limits
        lims(:,1)       = -params_glob%trs
        lims(:,2)       =  params_glob%trs
        lims_init(:,1)  = -2.
        lims_init(:,2)  =  2.
        ! convenienece images init
        call tmpimg%disc(self%ldim_ptcl, params_glob%smpd, self%msk, lmsk)
        call tmpimg%kill
        call ave%new(self%ldim_ptcl, params_glob%smpd)
        allocate(ptcl_mask(self%nptcls),ptcl_cls_mask(self%nptcls),source=.false.)
        hn = (params_glob%nframesgrp-1)/2
        do i = 1,9
            call neighbour_imgs(i)%new(self%ldim_ptcl,params_glob%smpd,wthreads=.false.)
            call tmpimgs(i)%new(self%ldim_ptcl,params_glob%smpd,wthreads=.false.)
        enddo
        ! Main class loop
        do icls = 1,self%ncls
            call progress(icls,self%ncls)
            if( self%pcls2D%get_state(icls) == 0 ) cycle
            pop = nint(self%pcls2D%get(icls,'pop'))
            if( pop == 0 ) cycle
            ptcl_mask     = .false.
            ptcl_cls_mask = .false.
            ! which frames to extract, indexing and further init
            do iptcl = 1,self%nptcls
                if( self%pptcl2D%get_state(iptcl) == 0 ) cycle
                if( nint(self%pptcl2D%get(iptcl,'class')) /= icls ) cycle
                ptcl_cls_mask(iptcl) = .true.
                fromto = [iptcl-hn, iptcl+hn]
                if(fromto(1) < 1          ) fromto = [1, params_glob%nframesgrp]
                if(fromto(2) > self%nptcls) fromto = [self%nptcls-params_glob%nframesgrp+1,self%nptcls]
                do j = fromto(1),fromto(2)
                    ptcl_mask(j) = .true.
                enddo
            enddo
            npop = count(ptcl_mask)
            allocate(pinds(npop),self%particles(npop),self%shifts(npop,2),self%background_pspecs(npop),&
                &self%shifts_sub(npop,2),self%e3(npop), ptcl_pos(npop,2))
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1,npop
                call self%particles(i)%new(self%ldim_ptcl, params_glob%smpd,wthreads=.false.)
                call self%background_pspecs(i)%new(self%ldim_ptcl, params_glob%smpd)
                call self%background_pspecs(i)%zero_and_unflag_ft
            enddo
            !$omp end parallel do
            ! Extract
            i = 0
            do iptcl = 1,self%nptcls
                if( .not.ptcl_mask(iptcl) ) cycle
                i        = i+1
                pinds(i) = iptcl
                self%e3(i)    = 0.
                ptcl_pos(i,:) = nint([self%pptcl2D%get(iptcl,'xpos'), self%pptcl2D%get(iptcl,'ypos')])
                noutside = 0
                call self%frame%read(self%frames(iptcl))
                call self%frame%window([ptcl_pos(i,1),ptcl_pos(i,2),1], params_glob%box, self%particles(i), noutside)
            enddo
            ! prep pftcc
            call self%pftcc%new(1,[1,npop])
            angstep    = abs(self%pftcc%get_rot(2)-self%pftcc%get_rot(1))
            self%nrots = self%pftcc%get_nrots()
            allocate(corrs(self%nrots))
            !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
            do ithr = 1,nthr_glob
                call self%ptcl_polarizers(ithr)%init_polarizer(self%pftcc, KBALPHA)
            enddo
            !$omp end parallel do
            ! prep reference
            call self%reference%zero_and_unflag_ft
            call self%reference%read(self%refs,icls)
            call self%reference%mask(params_glob%msk,'soft')
            ! call self%reference%mask(params_glob%msk,'soft',backgr=0.)
            call self%reference%fft
            call self%reference%init_polarizer(self%pftcc, KBALPHA)
            call self%reference%polarize(self%pftcc,1, .false., .true., mask=self%lresmsk)
            call self%pftcc%memoize_ffts
            !$omp parallel do default(shared) private(i,iptcl,ithr,cxy,irot,sdev_noise,shift) schedule(static) proc_bind(close)
            do i = 1,npop
                iptcl    = pinds(i)
                ithr     = omp_get_thread_num() + 1
                ! prep frame particle
                call self%particles(i)%neg
                call self%particles(i)%subtr_backgr_ramp(lmsk)
                call self%ptcl_polarizers(ithr)%copy(self%particles(i))
                call self%ptcl_polarizers(ithr)%noise_norm(lmsk, sdev_noise)
                call self%ptcl_polarizers(ithr)%mask(self%msk, 'soft', backgr=0.)
                call self%ptcl_polarizers(ithr)%fft
                call self%ptcl_polarizers(ithr)%polarize(self%pftcc,i,.true.,.true.,mask=self%lresmsk)
                ! exhaustive shift/in-plane rotation search
                call self%pftcc%prep_matchfilt(i, 1, 1)
                ithr = omp_get_thread_num() + 1
                call self%grad_shsrch_objs(ithr)%new(lims, lims_init, coarse_init=.true.,maxits=60)
                call self%grad_shsrch_objs(ithr)%set_indices(1, i)
                irot = 0
                cxy  = self%grad_shsrch_objs(ithr)%minimize_exhaustive(irot, halfrot=HALFROT)
                self%e3(i) = self%pftcc%get_rot(irot)
                self%shifts(i,:) = cxy(2:3)
                ! update parameters & prep for re-extraction
                shift                = real(ptcl_pos(i,:)) - self%shifts(i,:)
                ptcl_pos(i,:)        = nint(shift)
                self%shifts_sub(i,:) = real(ptcl_pos(i,:)) - shift
                call self%particles(i)%zero_and_unflag_ft
            enddo
            !$omp end parallel do
            call self%pftcc%kill
            ! Re-extraction, graphene backround subtraction prep
            do i = 1,npop
                iptcl         = pinds(i)
                noutside      = 0
                call self%frame%read(self%frames(iptcl))
                ! frame particles
                call self%frame%window([ptcl_pos(i,1),ptcl_pos(i,2),1], params_glob%box, self%particles(i), noutside)
                ! summing background neighbours power-spectra on the fly
                !$omp parallel do default(shared) collapse(2) private(cnt,j,k,ix,iy,noutside) schedule(static) proc_bind(close)
                do j = -1,1
                    do k = -1,1
                        cnt = 3*(j+1)+k+2
                        if( j==0 .and. k==0)then
                            valid_neighbour(cnt) = .false.
                            cycle
                        endif
                        ix  = ptcl_pos(i,1) + j*params_glob%box
                        iy  = ptcl_pos(i,2) + k*params_glob%box
                        noutside = 0
                        call neighbour_imgs(cnt)%zero_and_unflag_ft
                        call self%frame%window([ix,iy,1], params_glob%box, neighbour_imgs(cnt), noutside)
                        if( noutside > 0 )then
                            valid_neighbour(cnt) = .false.
                        else
                            valid_neighbour(cnt) = .true.
                            call neighbour_imgs(cnt)%neg
                            call neighbour_imgs(cnt)%norm
                            call neighbour_imgs(cnt)%zero_edgeavg
                            call neighbour_imgs(cnt)%fft
                            call neighbour_imgs(cnt)%ft2img('sqrt',tmpimgs(cnt))
                        endif
                    enddo
                enddo
                !$omp end parallel do
                do j = 1,9
                    if( valid_neighbour(j) ) call self%background_pspecs(i)%add(tmpimgs(j))
                enddo
                call self%background_pspecs(i)%div(real(count(valid_neighbour)))
            enddo
            ! preprocess newly extracted particles
            allocate(self%padded_imgs(npop))
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1,npop
                call self%particles(i)%neg
                call self%particles(i)%subtr_backgr_ramp(lmsk)
            enddo
            !$omp end parallel do
            ! graphene background subtraction
            call init_graphene_subtr(params_glob%box, params_glob%smpd)
            do i = 1,npop
                iptcl = pinds(i)
                call ave%zero_and_unflag_ft
                fromto = [i-hn, i+hn]
                if(fromto(1) < 1   ) fromto = [1, params_glob%nframesgrp]
                if(fromto(2) > npop) fromto = [npop-params_glob%nframesgrp+1,npop]
                do j = fromto(1),fromto(2)
                    call ave%add(self%background_pspecs(j))
                enddo
                call ave%div(real(params_glob%nframesgrp))
                call calc_peaks( ave, graphene_pos(1), graphene_pos(2) )
                call remove_lattices(self%particles(i), graphene_pos(1), graphene_pos(2))
            enddo
            call kill_graphene_subtr
            ! restore global parameters
            nthr_glob = inthr_glob
            !$ call omp_set_num_threads(nthr_glob)
            params_glob%ring2   = iring2
            params_glob%kfromto = ikfromto
            params_glob%kstop   = ikstop
            ! rotate, shift, sum
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1,npop
                call self%padded_imgs(i)%new([params_glob%boxpd,params_glob%boxpd,1],params_glob%smpd,wthreads=.false.)
                call self%particles(i)%pad_fft(self%padded_imgs(i))
                call self%padded_imgs(i)%shift2Dserial(-self%shifts_sub(i,:))
            enddo
            !$omp end parallel do
            ! WEIGHING FORK
            if( trim(params_glob%wcrit).ne.'no' )then
                ! user inputted weighing scheme
                ! re-fill pftcc for weights determination
                allocate(corrmat(npop,npop), source=-1.)
                call self%pftcc%new(npop,[1,npop])
                !$omp parallel do default(shared) private(i,ithr) schedule(static) proc_bind(close)
                do i = 1,npop
                    ithr = omp_get_thread_num() + 1
                    call self%particles(i)%fft
                    call self%particles(i)%shift2Dserial(-self%shifts_sub(i,:))
                    call self%particles(i)%ifft
                    call self%ptcl_polarizers(ithr)%copy(self%particles(i))
                    call self%ptcl_polarizers(ithr)%noise_norm(lmsk, sdev_noise)
                    call self%ptcl_polarizers(ithr)%mask(params_glob%msk, 'soft', backgr=0.)
                    call self%ptcl_polarizers(ithr)%fft
                    call self%ptcl_polarizers(ithr)%polarize(self%pftcc,i,.true., .true.,mask=self%lresmsk)
                    call self%ptcl_polarizers(ithr)%polarize(self%pftcc,i,.false.,.true.,mask=self%lresmsk)
                enddo
                !$omp end parallel do
                call self%pftcc%memoize_ffts
                ! correlations for weights
                !$omp parallel do default(shared) private(i,j,corrs) schedule(static) proc_bind(close)
                do i = 1,npop
                    corrmat(i,i) = 1.
                    do j = i+1,npop
                        call self%pftcc%prep_matchfilt(i, j, 1)
                        call self%pftcc%gencorrs(i,j,corrs)
                        corrmat(i,j) = maxval(corrs)
                        corrmat(j,i) = corrmat(i,j)
                    enddo
                enddo
                !$omp end parallel do
                call self%pftcc%kill
                !$omp parallel do default(shared) private(i,iptcl,ithr,j,ang,weights,threshold,ranked_weights,wsum)&
                !$omp schedule(static) proc_bind(close)
                do i = 1,npop
                    iptcl = pinds(i)
                    if( .not.ptcl_cls_mask(iptcl) )cycle
                    ithr   = omp_get_thread_num() + 1
                    ! weights
                    weights = corrs2weights(corrmat(:,i), params_glob%wcrit_enum)
                    if( npop > MAXNFRAMESPERPTCL )then
                        ranked_weights = weights
                        call hpsort(ranked_weights)
                        call reverse(ranked_weights) ! largest first
                        threshold = ranked_weights(MAXNFRAMESPERPTCL + 1)
                        where(weights <= threshold) weights = 0.
                        wsum = sum(weights)
                        if( wsum > TINY )then
                            weights = weights / wsum
                        else
                            weights = 0.
                        endif
                    endif
                    ! rotation and summation
                    call self%rot_imgs(ithr)%zero_and_flag_ft
                    do j = 1,npop
                        if( j == i )then
                            call self%rot_imgs(ithr)%add(self%padded_imgs(i), weights(i))
                        else
                            if( weights(j) < 0.000001 ) cycle
                            ! rotation relative to frame #i
                            ang = self%e3(j) - self%e3(i)
                            if( ang > 180. ) ang = ang - 360.
                            if( abs(ang) > 0.001 )then
                                call self%rotate_and_add(self%padded_imgs(j),ang,weights(j))
                            else
                                call self%rot_imgs(ithr)%add(self%padded_imgs(j),weights(j))
                            endif
                        endif
                    enddo
                    ! prep new particle
                    call self%rot_imgs(ithr)%ifft
                    call self%kernel_correction(self%rot_imgs(ithr))
                    call self%particles(i)%zero_and_unflag_ft
                    call self%rot_imgs(ithr)%clip(self%particles(i))
                    call self%particles(i)%subtr_backgr_ramp(lmsk)
                    call self%particles(i)%noise_norm(lmsk,sdev_noise)
                enddo
                !$omp end parallel do
                deallocate(corrmat)
            else
                ! Defaults to original uniform weights average with nframesgrp frames
                !$omp parallel do default(shared) private(i,iptcl,ithr,fromto,j,ang) schedule(static) proc_bind(close)
                do i = 1,npop
                    iptcl = pinds(i)
                    if( .not.ptcl_cls_mask(iptcl) )cycle
                    ithr   = omp_get_thread_num() + 1
                    fromto = [i-hn, i+hn]
                    if(fromto(1) < 1   ) fromto = [1, params_glob%nframesgrp]
                    if(fromto(2) > npop) fromto = [npop-params_glob%nframesgrp+1,npop]
                    call self%rot_imgs(ithr)%zero_and_flag_ft
                    do j = fromto(1),fromto(2)
                        if( j == i )then
                            call self%rot_imgs(ithr)%add(self%padded_imgs(i))
                        else
                            ! rotation relative to frame #i
                            ang = self%e3(j) - self%e3(i)
                            if( ang > 180. ) ang = ang - 360.
                            if( abs(ang) > 0.001 )then
                                call self%rotate_and_add(self%padded_imgs(j),ang)
                            else
                                call self%rot_imgs(ithr)%add(self%padded_imgs(j))
                            endif
                        endif
                    enddo
                    ! prep new particle
                    call self%rot_imgs(ithr)%ifft
                    call self%kernel_correction(self%rot_imgs(ithr))
                    call self%particles(i)%zero_and_unflag_ft
                    call self%rot_imgs(ithr)%clip(self%particles(i))
                    call self%particles(i)%subtr_backgr_ramp(lmsk)
                    call self%particles(i)%noise_norm(lmsk,sdev_noise)
                enddo
                !$omp end parallel do
            endif
            ! write new particles
            do i = 1,npop
                iptcl = pinds(i)
                if( ptcl_cls_mask(iptcl) )then
                    call self%particles(i)%write(new_stk_fname,iptcl)
                    call self%pptcl2D%e3set(iptcl,360.-self%e3(i))
                    call self%pptcl2D%set(iptcl,'x',0.)
                    call self%pptcl2D%set(iptcl,'y',0.)
                    call self%pptcl2D%set(iptcl,'xpos',real(ptcl_pos(i,1)))
                    call self%pptcl2D%set(iptcl,'ypos',real(ptcl_pos(i,2)))
                endif
            enddo
            ! memory management
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1,npop
                call self%particles(i)%kill
                call self%padded_imgs(i)%kill
                call self%background_pspecs(i)%kill
            enddo
            !$omp end parallel do
            deallocate(self%particles,self%padded_imgs,ptcl_pos)
            deallocate(pinds,self%shifts,self%shifts_sub,corrs,self%e3,self%background_pspecs)
        enddo
        ! zero state==0
        call tmpimg%new(self%ldim_ptcl, params_glob%smpd)
        call tmpimg%zero_and_unflag_ft
        do iptcl = 1,self%nptcls
            if( self%pptcl2D%get_state(iptcl)==0 ) call tmpimg%write(new_stk_fname,iptcl)
        enddo
        call tmpimg%read(new_stk_fname,self%nptcls)
        call tmpimg%write(new_stk_fname,self%nptcls)
        ! cleanup
        call tmpimg%kill
        call ave%kill
        do i = 1,9
            call neighbour_imgs(i)%kill
            call tmpimgs(i)%kill
        enddo
        do ithr = 1,nthr_glob
            call self%rot_imgs(ithr)%kill
            call self%ptcl_polarizers(ithr)%kill_polarizer
            call self%ptcl_polarizers(ithr)%kill
            call self%grad_shsrch_objs(ithr)%kill
        enddo
    end subroutine correct

    subroutine rotate_and_add( self, img, ang, w )
        use simple_kbinterpol
        class(motion_refine_nano), intent(inout) :: self
        class(image),              intent(inout) :: img
        real,                      intent(in)    :: ang
        real,            optional, intent(in)    :: w
        type(image_ptr)          :: pimg
        ! type(kbinterpol)         :: kbwin
        ! real, allocatable        :: w(:,:)
        ! integer, allocatable     :: cyc1(:), cyc2(:)
        complex :: fcomp
        real    :: loc(2), mat(2,2), dist(2), w_here
        integer :: lims(3,2), cyc_lims(3,2), cyc_limsR(2,3), phys_cmat(3), win_corner(2)
        integer :: h,k,l,ll,mm,m, ithr!, iwinsz,wdim,incr
        ithr = omp_get_thread_num() + 1
        w_here = 1.
        if( present(w) ) w_here = w
        call self%rot_imgs(ithr)%get_cmat_ptr(pimg%cmat)
        lims       = self%rot_imgs(ithr)%loop_lims(2)
        cyc_lims   = self%rot_imgs(ithr)%loop_lims(3)
        cyc_limsR(:,1) = cyc_lims(1,:)  ! limits in fortran layered format
        cyc_limsR(:,2) = cyc_lims(2,:)  ! to avoid copy on cyci_1d call
        ! kbwin  = kbinterpol(KBWINSZ, params_glob%alpha)
        ! wdim   = kbwin%get_wdim()
        ! iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        call rotmat2d(ang, mat)
        ! allocate(cyc1(wdim), cyc2(wdim), w(wdim, wdim))
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                loc = matmul(real([h,k]),mat)
                ! bi-linear interpolation
                win_corner = floor(loc) ! bottom left corner
                dist  = loc - real(win_corner)
                l     = cyci_1d(cyc_limsR(:,1), win_corner(1))
                ll    = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
                m     = cyci_1d(cyc_limsR(:,2), win_corner(2))
                mm    = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
                fcomp =         (1.-dist(1))*(1.-dist(2)) * img%get_fcomp2D(l, m)  ! bottom left corner
                fcomp = fcomp + (1.-dist(1))*dist(2)      * img%get_fcomp2D(l, mm) ! bottom right corner
                fcomp = fcomp + dist(1)*(1.-dist(2))      * img%get_fcomp2D(ll,m)  ! upper left corner
                fcomp = fcomp + dist(1)*dist(2)           * img%get_fcomp2D(ll,mm) ! upper right corner
                ! ! kb interpolation
                ! win_corner = nint(loc) - iwinsz
                ! ! weights kernel
                ! w = 1.
                ! do l=1,wdim
                !     incr = l - 1
                !     ! circular addresses
                !     cyc1(l) = cyci_1d(cyc_limsR(:,1), win_corner(1) + incr)
                !     cyc2(l) = cyci_1d(cyc_limsR(:,2), win_corner(2) + incr)
                !     ! interpolation kernel matrix
                !     w(l,:) = w(l,:) * kbwin%apod( real(win_corner(1) + incr) - loc(1) )
                !     w(:,l) = w(:,l) * kbwin%apod( real(win_corner(2) + incr) - loc(2) )
                ! enddo
                ! w = w / sum(w)
                ! ! interpolation
                ! fcomp = cmplx(0.,0.)
                ! do l=1,wdim
                !     do m=1,wdim
                !         if( w(l,m) < TINY ) cycle
                !         fcomp = fcomp + img%get_fcomp2D(cyc1(l),cyc2(m)) * w(l,m)
                !     end do
                ! end do
                ! addition
                phys_cmat = img%comp_addr_phys([h,k,0])
                pimg%cmat(phys_cmat(1),phys_cmat(2),1) = pimg%cmat(phys_cmat(1),phys_cmat(2),1) + w_here*fcomp
            end do
        end do
    end subroutine rotate_and_add

    subroutine kernel_correction( self, img )
        class(motion_refine_nano), intent(in)    :: self
        class(image),              intent(inout) :: img
        real    :: dist(2), pid, pad_sc, sinc, center(2)
        integer :: i,j
        center = real(params_glob%boxpd/2 + 1)
        pad_sc = 1. / real(params_glob%boxpd)
        do i = 1,params_glob%boxpd
            dist(1) = pad_sc*(real(i)-center(1))
            do j = 1,params_glob%boxpd
                dist(2) = pad_sc*(real(j)-center(2))
                pid     = PI*sqrt(sum(dist**2.))
                if( pid > TINY )then
                    sinc = sin(pid) / pid
                    call img%div([i,j,1], sinc*sinc)
                endif
            enddo
        enddo
    end subroutine kernel_correction

    subroutine kill(self)
        class(motion_refine_nano), intent(inout) :: self
        nullify(self%pptcl2D)
        nullify(self%pcls2D)
        nullify(self%pframes)
        call self%frame%kill
        call self%reference%kill_polarizer
        call self%reference%kill
        call self%pftcc%kill
        if(allocated(self%lresmsk))deallocate(self%lresmsk)
        self%exists = .false.
    end subroutine kill

end module simple_motion_refine_nano

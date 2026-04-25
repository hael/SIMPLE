!@descr: submodule for controlling various state-related things in the polarops module
submodule (simple_polarft_calc) simple_polarft_ops_state
use simple_memoize_ft_maps, only: memoize_ft_maps
implicit none
#include "simple_local_flags.inc"

type :: polar_rec_group_plan
    integer :: ngroups = 0
    logical,  allocatable :: ptcl_active(:), ptcl_even(:)
    real(dp), allocatable :: ptcl_w(:)
    integer,  allocatable :: ptcl_proj(:), ptcl_inpl(:), ptcl_group(:), members(:)
    integer,  allocatable :: proj(:), inpl(:), counts(:), offsets(:), next(:)
    real,     allocatable :: ptcl_R(:,:,:), ptcl_euls(:,:), R(:,:,:), euls(:,:)
end type polar_rec_group_plan

type :: polar_rec_ref_plan
    integer :: nrefs = 0
    real, allocatable :: tR(:,:,:), normal(:,:)
end type polar_rec_ref_plan

type :: polar_rec_sym_group_plan
    real, allocatable :: R(:,:,:), normal(:,:)
end type polar_rec_sym_group_plan

type :: polar_rec_pair_plan
    integer :: ntasks = 0
    integer, allocatable :: counts(:), offsets(:), group(:)
    integer, allocatable :: irot0(:), nangle_scan(:), max_irot_delta(:)
    real,    allocatable :: R2d(:,:,:), pol2cart(:,:), phi(:), sin_theta(:)
end type polar_rec_pair_plan

contains

    !> Module initialization
    module subroutine polar_cavger_new( self, l_comlin, nrefs )
        class(polarft_calc), intent(inout) :: self
        logical,             intent(in)    :: l_comlin
        integer,   optional, intent(in)    :: nrefs
        call self%polar_cavger_kill
        self%l_comlin = l_comlin
        if( present(nrefs) )then
            self%ncls = nrefs
        else
            self%ncls = self%nrefs
        endif
        allocate(self%prev_eo_pops(2,self%ncls), self%eo_pops(2,self%ncls), source=0)
        allocate(self%pfts_even(self%pftsz, self%kfromto(1):self%interpklim, self%ncls),&
                &self%pfts_odd( self%pftsz, self%kfromto(1):self%interpklim, self%ncls),&
                &self%ctf2_even(self%pftsz, self%kfromto(1):self%interpklim, self%ncls),&
                &self%ctf2_odd( self%pftsz, self%kfromto(1):self%interpklim, self%ncls),&
                &self%pfts_merg(self%pftsz, self%kfromto(1):self%interpklim, self%ncls))
        call self%polar_cavger_zero_pft_refs
        self%pfts_merg = DCMPLX_ZERO
    end subroutine polar_cavger_new

    module subroutine polar_cavger_zero_pft_refs( self )
        class(polarft_calc), intent(inout) :: self
        !$omp parallel workshare
        self%pfts_even = DCMPLX_ZERO
        self%pfts_odd  = DCMPLX_ZERO
        self%ctf2_even = 0.d0
        self%ctf2_odd  = 0.d0
        !$omp end parallel workshare
    end subroutine polar_cavger_zero_pft_refs

    module subroutine polar_cavger_set_ref_pft( self, icls, which )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: icls
        character(len=*),    intent(in)    :: which
        select case(trim(which))
            case('merged')
                self%pfts_refs_even(:,:,icls) = cmplx(self%pfts_merg(:,self%kfromto(1):self%interpklim,icls),kind=sp)
            case('even')
                self%pfts_refs_even(:,:,icls) = cmplx(self%pfts_even(:,self%kfromto(1):self%interpklim,icls),kind=sp)
            case('odd')
                self%pfts_refs_odd(:,:,icls)  = cmplx(self%pfts_odd(:,self%kfromto(1):self%interpklim,icls),kind=sp)
        end select
    end subroutine polar_cavger_set_ref_pft

    module subroutine polar_cavger_calc_pops( self, spproj )
        class(polarft_calc),       intent(inout) :: self
        class(sp_project), target, intent(in)    :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer :: i, icls, iptcl, eo
        integer :: eo_pops(2,self%ncls)
        logical :: l_3D
        l_3D = .false.
        select case(trim(self%p_ptr%oritype))
        case('ptcl2D')
            ptcl_field => spproj%os_ptcl2D
            cls_field  => spproj%os_cls2D
        case('ptcl3D')
            ptcl_field => spproj%os_ptcl3D
            cls_field  => spproj%os_cls3D
            l_3D       = .true.
        case DEFAULT
            THROW_HARD('Unsupported ORITYPE: '//trim(self%p_ptr%oritype))
        end select
        eo_pops = 0
        !$omp parallel do schedule(guided) proc_bind(close) default(shared)&
        !$omp private(iptcl,eo,icls)&
        !$omp reduction(+:eo_pops)
        do iptcl = 1,ptcl_field%get_noris()
            if( ptcl_field%get_state(iptcl) == 0  ) cycle
            if( ptcl_field%get(iptcl,'w') < SMALL ) cycle
            eo = ptcl_field%get_eo(iptcl)+1
            if( l_3D )then
                icls = ptcl_field%get_proj(iptcl)
            else
                icls = ptcl_field%get_class(iptcl)
            endif
            if( icls < 1 .or. icls > self%ncls ) cycle
            eo_pops(eo,icls) = eo_pops(eo,icls) + 1
        enddo
        !$omp end parallel do
        self%eo_pops      = eo_pops
        self%prev_eo_pops = 0
        if( cls_field%get_noris() == self%ncls )then
            do i = 1,self%ncls
                if( .not.cls_field%isthere(i,'prev_pop_even') ) cycle
                self%prev_eo_pops(1,i) = cls_field%get_int(i,'prev_pop_even')
                self%prev_eo_pops(2,i) = cls_field%get_int(i,'prev_pop_odd')
            enddo
        endif
    end subroutine polar_cavger_calc_pops

    !>  \brief  Updates Fourier components and normalization matrices with new particles
    module subroutine polar_cavger_update_sums( self, nptcls, pinds, spproj, incr_shifts, is3D )
        class(polarft_calc),         intent(inout) :: self
        integer,                     intent(in)    :: nptcls
        integer,                     intent(in)    :: pinds(nptcls)
        class(sp_project),           intent(inout) :: spproj
        real,              optional, intent(in)    :: incr_shifts(2,nptcls)
        logical,           optional, intent(in)    :: is3d
        class(oris), pointer :: spproj_field
        complex(sp) :: rptcl(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: rctf2(self%pftsz,self%kfromto(1):self%interpklim), w
        real(sp)    :: rctf(self%pftsz,self%kfromto(1):self%interpklim)
        real        :: sigma2(self%kfromto(1):self%interpklim), incr_shift(2)
        integer     :: eopops(2,self%ncls), i, icls, iptcl, irot, k
        logical     :: l_ctf, l_even, l_3D, l_shift
        l_3D = .false.
        if( present(is3D) ) l_3D = is3D
        l_shift = .false.
        if( present(incr_shifts) ) l_shift = .true.
        ! retrieve particle info & pointers
        call spproj%ptr2oritype(self%p_ptr%oritype, spproj_field)
        l_ctf = self%is_with_ctf()
        ! update classes
        eopops = 0
        !$omp parallel do default(shared) proc_bind(close) schedule(static) reduction(+:eopops)&
        !$omp private(i,iptcl,w,l_even,icls,irot,incr_shift,rptcl,rctf,rctf2,k,sigma2)
        do i = 1,nptcls
            ! particles parameters
            iptcl = pinds(i)
            if( spproj_field%get_state(iptcl) == 0  ) cycle
            w = real(spproj_field%get(iptcl,'w'),dp)
            if( w < DSMALL ) cycle
            l_even = spproj_field%get_eo(iptcl)==0
            if( l_3D )then
                icls = spproj_field%get_proj(iptcl)
            else
                icls = spproj_field%get_class(iptcl)
            endif
            irot = self%get_roind_fast(spproj_field%e3get(iptcl))
            if( l_shift )then
                incr_shift = incr_shifts(:,i)
                ! weighted restoration
                if( any(abs(incr_shift) > 1.e-6) ) call self%shift_ptcl(iptcl, -incr_shift)
            endif
            ! Particle rotation
            call self%rotate_ptcl(i, irot, rptcl)
            ! Particle weight
            rptcl = real(w) * rptcl
            ! Particle ML regularization
            if( self%p_ptr%l_ml_reg )then
                sigma2 = self%sigma2_noise(self%kfromto(1):self%interpklim,iptcl)
                do k = self%kfromto(1),self%interpklim
                    rptcl(:,k) = rptcl(:,k) / sigma2(k)
                enddo
            endif
            ! Array updates
            if( l_ctf )then
                ! weighted CTF2
                call self%rotate_ctf(i, irot, rctf)
                rctf2 = w * real(rctf,kind=dp)**2
                rptcl = rptcl * rctf    ! PhFlip(X).|CTF|
                ! CTF2 ML regularization
                if( self%p_ptr%l_ml_reg )then
                    do k = self%kfromto(1),self%interpklim
                        rctf2(:,k) = rctf2(:,k) / real(sigma2(k),dp)
                    enddo
                endif
                if( l_even )then
                    !$omp critical
                    self%pfts_even(:,:,icls) = self%pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                    self%ctf2_even(:,:,icls) = self%ctf2_even(:,:,icls) + rctf2
                    !$omp end critical
                else
                    !$omp critical
                    self%pfts_odd(:,:,icls)  = self%pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                    self%ctf2_odd(:,:,icls)  = self%ctf2_odd(:,:,icls)  + rctf2
                    !$omp end critical
                endif
            else
                if( self%p_ptr%l_ml_reg )then
                    ! CTF2=1 & ML regularization
                    do k = self%kfromto(1),self%interpklim
                        rctf2(:,k) = w / real(sigma2(k),dp)
                    enddo
                    if( l_even )then
                        !$omp critical
                        self%pfts_even(:,:,icls) = self%pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                        self%ctf2_even(:,:,icls) = self%ctf2_even(:,:,icls) + rctf2
                        !$omp end critical
                    else
                        !$omp critical
                        self%pfts_odd(:,:,icls)  = self%pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                        self%ctf2_odd(:,:,icls)  = self%ctf2_odd(:,:,icls)  + rctf2
                        !$omp end critical
                    endif
                else
                    if( l_even )then
                        !$omp critical
                        self%pfts_even(:,:,icls) = self%pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                        self%ctf2_even(:,:,icls) = self%ctf2_even(:,:,icls) + w
                        !$omp end critical
                    else
                        !$omp critical
                        self%pfts_odd(:,:,icls)  = self%pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                        self%ctf2_odd(:,:,icls)  = self%ctf2_odd(:,:,icls)  + w
                        !$omp end critical
                    endif
                endif
            endif
            ! total population
            if( l_even )then
                eopops(1,icls) = eopops(1,icls) + 1
            else
                eopops(2,icls) = eopops(2,icls) + 1
            endif
        enddo
        !$omp end parallel do
        self%eo_pops = self%eo_pops + eopops
        ! cleanup
        nullify(spproj_field)
    end subroutine polar_cavger_update_sums

    module subroutine polar_cavger_insert_ptcls_oversamp( self, eulspace, ptcl_field, symop, nptcls, pinds, fpls )
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(polarft_calc),        intent(inout) :: self
        class(oris),                intent(in)    :: eulspace
        class(oris), pointer,       intent(inout) :: ptcl_field
        class(sym),                 intent(in)    :: symop
        integer,                    intent(in)    :: nptcls, pinds(nptcls)
        class(fplane_type), target, intent(inout) :: fpls(nptcls)
        real,     parameter   :: zvec(3) = [0.,0.,1.]                   ! normal vector
        real,     parameter   :: DT      = KBWINSZ                      ! distance threshold
        real(dp), parameter   :: PF2     = real(OSMPL_PAD_FAC**2,dp)    ! Oversampling factor
        type(kbinterpol)      :: kb
        complex(dp) :: fcomp
        integer     :: cl_irots(self%pftsz*(self%interpklim-self%kfromto(1)+1))
        integer     :: cl_shs(  self%pftsz*(self%interpklim-self%kfromto(1)+1))
        integer     :: cl_hhs(  self%pftsz*(self%interpklim-self%kfromto(1)+1))
        integer     :: cl_kks(  self%pftsz*(self%interpklim-self%kfromto(1)+1))
        real(dp)    :: cl_wts(  self%pftsz*(self%interpklim-self%kfromto(1)+1))
        real(dp)    :: pw, ctfsq, w, dkb01 ,dkb02, dkb03, base_w
        type(polar_rec_group_plan)     :: groups
        type(polar_rec_ref_plan)       :: refs
        type(polar_rec_sym_group_plan) :: symgroups
        type(polar_rec_pair_plan)      :: pairs
        real        :: Rproj(3,3), Rsym(3,3), R(3,3), R2d(2,2), pol2cart(2), hk(2), proj_euls(3)
        real        :: normal_proj(3), normal_ptcl(3), rhk(2), euls(3), phi, theta
        real        :: dang, sin_dang, dCL, dz, psi, sin_theta, de3, proj_h, proj_k, max_sin_dang, max_dang
        integer     :: flims(2,3), nsym, iproj, i, iptcl, hh, kk, sh, irot, isym, ncl, icl
        integer     :: irot0, iang, max_irot_delta, nangle_scan
        integer     :: srcproj, inpl_ind, igrp, jgrp, imember, member_start, member_end, itask
        logical     :: l_even, l_self
        ! Interpolation parameters
        kb    = kbinterpol(KBWINSZ, KBALPHA)
        flims = transpose(fpls(1)%frlims)
        dkb01 = real(kb%apod(0.),dp)
        dkb02 = dkb01*dkb01
        dkb03 = dkb01*dkb02
        ! Looping over the un-mirrored asymmetric unit
        refs%nrefs = eulspace%get_noris()/2
        ! geometric precomputations
        allocate(groups%ptcl_active(nptcls), groups%ptcl_even(nptcls), groups%ptcl_w(nptcls))
        allocate(groups%ptcl_proj(nptcls), groups%ptcl_inpl(nptcls), groups%ptcl_group(nptcls))
        allocate(groups%members(nptcls), groups%proj(nptcls), groups%inpl(nptcls), groups%counts(nptcls))
        allocate(groups%offsets(nptcls), groups%next(nptcls), groups%ptcl_R(3,3,nptcls))
        allocate(groups%ptcl_euls(3,nptcls), groups%R(3,3,nptcls), groups%euls(3,nptcls))
        groups%ptcl_active = .false.
        groups%ptcl_even   = .false.
        groups%ptcl_w      = 0.d0
        groups%ptcl_proj   = 0
        groups%ptcl_inpl   = 0
        groups%ptcl_group  = 0
        !$omp parallel do default(shared) private(i,iptcl) schedule(static) proc_bind(close)
        do i = 1,nptcls
            iptcl = pinds(i)
            if( ptcl_field%get_state(iptcl) == 0 ) cycle
            groups%ptcl_w(i) = real(ptcl_field%get(iptcl, 'w'), dp)
            if( groups%ptcl_w(i) < 1.d-6 ) cycle
            groups%ptcl_active(i) = .true.
            groups%ptcl_even(i)   = ptcl_field%get_eo(iptcl) == 0
            groups%ptcl_proj(i)   = ptcl_field%get_proj(iptcl)
            groups%ptcl_euls(:,i) = ptcl_field%get_euler(iptcl)
            groups%ptcl_R(:,:,i)  = euler2m(groups%ptcl_euls(:,i))
            groups%ptcl_inpl(i)   = self%get_roind_fast(groups%ptcl_euls(3,i))
        enddo
        !$omp end parallel do
        ! group particles by source projection and in-plane bin to reuse geometry
        groups%ngroups = 0
        groups%proj    = 0
        groups%inpl    = 0
        groups%counts  = 0
        groups%offsets = 0
        groups%members = 0
        do i = 1,nptcls
            if( .not. groups%ptcl_active(i) ) cycle
            srcproj = groups%ptcl_proj(i)
            if( srcproj < 1 .or. srcproj > self%ncls ) srcproj = -i
            inpl_ind = groups%ptcl_inpl(i)
            igrp = 0
            do jgrp = 1,groups%ngroups
                if( groups%proj(jgrp) /= srcproj ) cycle
                if( groups%inpl(jgrp) /= inpl_ind ) cycle
                de3 = modulo(abs(groups%ptcl_euls(3,i) - groups%euls(3,jgrp)), 360.0)
                if( de3 > 180.0 ) de3 = 360.0 - de3
                if( de3 > 1.e-3 ) cycle
                igrp = jgrp
                exit
            enddo
            if( igrp == 0 )then
                groups%ngroups = groups%ngroups + 1
                igrp = groups%ngroups
                groups%proj(igrp)   = srcproj
                groups%inpl(igrp)   = inpl_ind
                groups%euls(:,igrp) = groups%ptcl_euls(:,i)
                groups%R(:,:,igrp)  = groups%ptcl_R(:,:,i)
            endif
            groups%ptcl_group(i) = igrp
            groups%counts(igrp) = groups%counts(igrp) + 1
        enddo
        member_start = 1
        do igrp = 1,groups%ngroups
            groups%offsets(igrp) = member_start
            groups%next(igrp)    = member_start
            member_start         = member_start + groups%counts(igrp)
        enddo
        do i = 1,nptcls
            if( .not. groups%ptcl_active(i) ) cycle
            igrp = groups%ptcl_group(i)
            imember = groups%next(igrp)
            groups%members(imember) = i
            groups%next(igrp) = imember + 1
        enddo
        allocate(refs%tR(3,3,refs%nrefs), refs%normal(3,refs%nrefs))
        do iproj = 1,refs%nrefs
            proj_euls          = eulspace%get_euler(iproj)
            Rproj              = euler2m(proj_euls)
            refs%tR(:,:,iproj) = transpose(Rproj)
            refs%normal(:,iproj) = matmul(zvec, Rproj)
        enddo
        allocate(symgroups%R(3,3,groups%ngroups), symgroups%normal(3,groups%ngroups))
        ! Main loop
        nsym = symop%get_nsym()
        do isym = 1,nsym
            call symop%get_sym_rmat(isym, Rsym)
            do igrp = 1,groups%ngroups
                if( isym == 1 )then
                    symgroups%R(:,:,igrp) = groups%R(:,:,igrp)
                else
                    symgroups%R(:,:,igrp) = matmul(groups%R(:,:,igrp), Rsym)
                endif
                symgroups%normal(:,igrp) = matmul(zvec, symgroups%R(:,:,igrp))
            enddo
            if( allocated(pairs%counts) ) deallocate(pairs%counts, pairs%offsets)
            if( allocated(pairs%group) )then
                deallocate(pairs%group, pairs%irot0, pairs%nangle_scan, pairs%max_irot_delta)
                deallocate(pairs%R2d, pairs%pol2cart, pairs%phi, pairs%sin_theta)
            endif
            allocate(pairs%counts(refs%nrefs), pairs%offsets(refs%nrefs+1), source=0)
            !$omp parallel do default(shared) private(iproj,igrp,normal_proj,normal_ptcl,l_self)&
            !$omp& proc_bind(close) schedule(static)
            do iproj = 1,refs%nrefs
                normal_proj = refs%normal(:,iproj)
                do igrp = 1,groups%ngroups
                    normal_ptcl = symgroups%normal(:,igrp)
                    l_self = myacos(abs(dot_product(normal_proj, normal_ptcl))) < 1.e-4
                    if( .not. l_self ) pairs%counts(iproj) = pairs%counts(iproj) + 1
                enddo
            enddo
            !$omp end parallel do
            pairs%offsets(1) = 1
            do iproj = 1,refs%nrefs
                pairs%offsets(iproj+1) = pairs%offsets(iproj) + pairs%counts(iproj)
            enddo
            pairs%ntasks = pairs%offsets(refs%nrefs+1) - 1
            allocate(pairs%group(pairs%ntasks),       pairs%irot0(pairs%ntasks))
            allocate(pairs%nangle_scan(pairs%ntasks), pairs%max_irot_delta(pairs%ntasks))
            allocate(pairs%R2d(2,2,pairs%ntasks),     pairs%pol2cart(2,pairs%ntasks))
            allocate(pairs%phi(pairs%ntasks),         pairs%sin_theta(pairs%ntasks))
            !$omp parallel do default(shared) private(iproj,igrp,itask,normal_proj,normal_ptcl,l_self)&
            !$omp& private(R,euls,psi,phi,theta,sin_theta,pol2cart,R2d,irot0)&
            !$omp& private(max_sin_dang,max_dang,max_irot_delta,nangle_scan)&
            !$omp& proc_bind(close) schedule(static)
            do iproj = 1,refs%nrefs
                normal_proj = refs%normal(:,iproj)
                itask = pairs%offsets(iproj)
                do igrp = 1,groups%ngroups
                    normal_ptcl = symgroups%normal(:,igrp)
                    l_self = myacos(abs(dot_product(normal_proj, normal_ptcl))) < 1.e-4
                    if( l_self ) cycle
                    R = matmul(symgroups%R(:,:,igrp), refs%tR(:,:,iproj))
                    euls = m2euler(R)
                    psi = 180.0 - euls(3)
                    phi = euls(1)
                    if( phi > 180.0 )then
                        phi = phi - 180.0
                        psi = psi + 180.0
                    else if( phi < 0.0 )then
                        phi = phi + 180.0
                        psi = psi - 180.0
                    endif
                    theta     = euls(2)
                    sin_theta = sin(deg2rad(theta))
                    pol2cart  = [sin(deg2rad(phi)), -cos(deg2rad(phi))]
                    call rotmat2d(180.0+psi-phi, R2D)
                    irot0 = modulo(nint(phi / self%dang), self%pftsz) + 1
                    max_irot_delta = 0
                    if( abs(sin_theta) < 1.e-6 )then
                        nangle_scan = self%pftsz
                    else
                        max_sin_dang = DT / (abs(sin_theta) * real(self%kfromto(1)*OSMPL_PAD_FAC))
                        if( max_sin_dang >= 1.0 )then
                            nangle_scan = self%pftsz
                        else
                            max_dang = rad2deg(asin(max_sin_dang))
                            max_irot_delta = ceiling(max_dang / self%dang) + 1
                            if( 2*max_irot_delta + 1 >= self%pftsz )then
                                nangle_scan = self%pftsz
                            else
                                nangle_scan = 2*max_irot_delta + 1
                            endif
                        endif
                    endif
                    pairs%group(itask)          = igrp
                    pairs%phi(itask)            = phi
                    pairs%sin_theta(itask)      = sin_theta
                    pairs%pol2cart(:,itask)     = pol2cart
                    pairs%R2d(:,:,itask)        = R2D
                    pairs%irot0(itask)          = irot0
                    pairs%nangle_scan(itask)    = nangle_scan
                    pairs%max_irot_delta(itask) = max_irot_delta
                    itask = itask + 1
                enddo
            enddo
            !$omp end parallel do
            !$omp parallel do default(shared) private(iproj,itask)&
            !$omp& private(sh,irot,R2d,hk,rhk,w,ctfsq,fcomp,hh,kk,sin_theta,dang,dCL)&
            !$omp& private(sin_dang,dz,l_even,phi,pol2cart)&
            !$omp& private(igrp,imember,member_start,member_end)&
            !$omp& private(ncl,icl,cl_irots,cl_shs,cl_hhs,cl_kks,cl_wts,proj_h,proj_k,base_w)&
            !$omp& private(irot0,iang,max_irot_delta,nangle_scan,i,pw)&
            !$omp& proc_bind(close) schedule(static)
            do iproj = 1,refs%nrefs
                do itask = pairs%offsets(iproj), pairs%offsets(iproj+1) - 1
                    igrp           = pairs%group(itask)
                    member_start   = groups%offsets(igrp)
                    member_end     = member_start + groups%counts(igrp) - 1
                    phi            = pairs%phi(itask)
                    sin_theta      = pairs%sin_theta(itask)
                    pol2cart       = pairs%pol2cart(:,itask)
                    R2D            = pairs%R2d(:,:,itask)
                    irot0          = pairs%irot0(itask)
                    nangle_scan    = pairs%nangle_scan(itask)
                    max_irot_delta = pairs%max_irot_delta(itask)
                    ncl = 0
                    do iang = 1, nangle_scan
                        if( nangle_scan == self%pftsz )then
                            irot = iang
                        else
                            irot = modulo(irot0 + iang - max_irot_delta - 2, self%pftsz) + 1
                        endif
                        ! angle betwen CL an current reprojection line
                        dang     = modulo(abs(self%angtab(irot) - phi), 180.0)
                        sin_dang = sin(deg2rad(dang))
                        ! dz increases with sh, so if the first shell is too far, all shells are too far
                        if( abs(sin_theta * sin_dang) * real(self%kfromto(1)*OSMPL_PAD_FAC) > DT ) cycle
                        do sh = self%kfromto(1), self%interpklim
                            proj_h = real(sh*OSMPL_PAD_FAC) * pol2cart(1)
                            proj_k = real(sh*OSMPL_PAD_FAC) * pol2cart(2)
                            ! rejects polar points to far from CL
                            hk(1) = real(OSMPL_PAD_FAC) * self%polar(1,sh,irot)
                            if( abs(proj_h - hk(1)) > DT ) exit
                            hk(2) = real(OSMPL_PAD_FAC) * self%polar(2,sh,irot)
                            if( abs(proj_k - hk(2)) > DT ) exit
                            ! distance from point to common line sin_dang * sqrt(sum(hk**2)) simplifies to:
                            dCL = sin_dang * real(sh*OSMPL_PAD_FAC)
                            ! relative altitude of the point to slice
                            dz  = sin_theta * dCL
                            ! rejects out-of-plane point
                            if( abs(dz) > DT ) exit
                            ! 2D mapping
                            rhk = matmul(hk, R2D)
                            ! NN interpolation with KB weight
                            hh     = nint(rhk(1))
                            kk     = nint(rhk(2))
                            base_w = real(kb%apod(rhk(1)-real(hh)),dp) * real(kb%apod(rhk(2)-real(kk)),dp)
                            base_w = base_w * real(kb%apod(dz),dp)
                            if( base_w < DTINY ) cycle
                            ncl = ncl + 1
                            cl_irots(ncl) = irot
                            cl_shs(ncl)   = sh
                            cl_hhs(ncl)   = cyci_1d(flims(:,1), hh)
                            cl_kks(ncl)   = cyci_1d(flims(:,2), kk)
                            cl_wts(ncl)   = base_w / dkb03
                        enddo
                    enddo
                    if( ncl == 0 ) cycle
                    do imember = member_start, member_end
                        i      = groups%members(imember)
                        pw     = groups%ptcl_w(i)
                        l_even = groups%ptcl_even(i)
                        ! Loop over the PFT and interpolate relevant components
                        do icl = 1,ncl
                            irot  = cl_irots(icl)
                            sh    = cl_shs(icl)
                            hh    = cl_hhs(icl)
                            kk    = cl_kks(icl)
                            w     = pw * cl_wts(icl)
                            fcomp = PF2 * w * cmplx(fplane_get_cmplx(fpls(i), hh,kk), kind=dp)
                            ctfsq =       w * real(fplane_get_ctfsq(fpls(i),  hh,kk), dp)
                            if( l_even )then
                                self%pfts_even(irot,sh,iproj) = self%pfts_even(irot,sh,iproj) + fcomp
                                self%ctf2_even(irot,sh,iproj) = self%ctf2_even(irot,sh,iproj) + ctfsq
                            else
                                self%pfts_odd(irot,sh,iproj)  = self%pfts_odd(irot,sh,iproj)  + fcomp
                                self%ctf2_odd(irot,sh,iproj)  = self%ctf2_odd(irot,sh,iproj)  + ctfsq
                            endif
                        enddo
                    enddo           ! particle
                enddo               ! common-line task
            enddo                   ! slice
            !$omp end parallel do
        enddo                       ! symmetry
    end subroutine polar_cavger_insert_ptcls_oversamp

    module subroutine polar_cavger_kill( self )
        class(polarft_calc), intent(inout) :: self
        if( allocated(self%pfts_even)    ) deallocate(self%pfts_even)
        if( allocated(self%pfts_odd)     ) deallocate(self%pfts_odd)
        if( allocated(self%pfts_merg)    ) deallocate(self%pfts_merg)
        if( allocated(self%ctf2_even)    ) deallocate(self%ctf2_even)
        if( allocated(self%ctf2_odd)     ) deallocate(self%ctf2_odd)
        if( allocated(self%eo_pops)      ) deallocate(self%eo_pops)
        if( allocated(self%prev_eo_pops) ) deallocate(self%prev_eo_pops)
        self%ncls     = 0
        self%l_comlin = .false.
    end subroutine polar_cavger_kill

    ! Estimates the center of the volume based on the distribution of
    ! the individual particles in-plane offsets and map the shifts to both
    ! the particles and the references stored in the pftc
    module subroutine center_3Dpolar_refs( self, algndoc, algnrefs )
        class(polarft_calc), intent(inout) :: self
        class(oris),         intent(inout) :: algndoc
        class(oris),         intent(in)    :: algnrefs
        real    :: R(3,3),offset3D(3), offset2D(3)
        integer :: iref
        ! estimate 3D offset from particle alignement parameters
        call algndoc%calc_avg_offset3D(offset3D, state=1)
        ! report 3D offset to particles 2D offsets
        call algndoc%map3dshift22d(-offset3D, state=1)
        ! report 3D offset to alignment references
        !$omp parallel do proc_bind(close) default(shared) private(iref,R,offset2D)
        do iref = 1,self%ncls
            ! Projection direction rotation matrix
            R = euler2m([algnrefs%e1get(iref), algnrefs%e2get(iref), 0.0])
            ! 3D Shift rotated with respect to projection direction
            offset2D = matmul(R, offset3D)
            ! Apply offset to e/o references
            call self%shift_ref(iref, offset2D(1:2))
        enddo
        !$omp end parallel do
    end subroutine center_3Dpolar_refs

end submodule simple_polarft_ops_state

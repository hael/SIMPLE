!@descr: submodule for controlling various state-related things in the polarops module
submodule (simple_polarft_calc) simple_polarft_ops_state
use simple_memoize_ft_maps, only: memoize_ft_maps
use simple_fgrid_obsfield,  only: fgrid_obs_field_eo
implicit none
#include "simple_local_flags.inc"

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

    ! Direct particle-to-polar insertion routine
    module subroutine polar_cavger_insert_ptcls_direct( self, eulspace, ptcl_field, symop, nptcls, pinds, fpls )
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
        real(dp), parameter   :: SELFW   = 0.1d0                        ! Weights atributed to self
        type :: ref_cache_type
            real    :: tR(3,3)
            real    :: normal(3)
            integer :: sym_self_proj
            logical :: sym_done
        end type ref_cache_type
        type :: ptcl_cache_type
            logical  :: active
            logical  :: even
            integer  :: proj
            integer  :: self_proj
            real(dp) :: w
            real(dp) :: pw_self  ! cached: w * SELFW / dkb02
            real(dp) :: pw_cl    ! cached: w / dkb03
            real     :: euls(3)
            real     :: R(3,3)
            real     :: Rsym(3,3)
        end type ptcl_cache_type
        type(kbinterpol)      :: kb
        type(ref_cache_type)  :: refs(self%ncls)
        type(ptcl_cache_type) :: ptcls(nptcls)
        complex(dp) :: rot_ptcl(self%pftsz, self%kfromto(1):self%interpklim), fcomp
        real(dp)    :: rot_ctfsq(self%pftsz, self%kfromto(1):self%interpklim)
        real(dp)    :: pw, ctfsq, w, dkb01 ,dkb02, dkb03
        real        :: proj_cl_addr(2,self%kfromto(1):self%interpklim), Rproj(3,3), tRproj(3,3)
        real        :: Rsym(3,3), R(3,3), Rptcl(3,3), R2d(2,2), pol2cart(2), hk(2), proj_euls(3)
        real        :: normal_ptcl(3), rhk(2), euls(3), phi, theta
        real        :: dang, sin_dang, dCL, dz, psi, sin_theta, dotp, best_dot
        integer     :: flims(2,3), nsym, iproj, i, iptcl, hh, kk, noris, nrefs, sh, irot, jrot, drot, isym
        integer     :: srcproj, best_proj
        logical     :: l_even, l_self, has_mirr
        ! Interpolation parameters.  The NN oversampled polar insertion path
        ! uses kb%apod_fast for both normalization and per-sample KB weights
        ! so that all factors come from the same degree-14 I0 power-series
        ! approximation documented in simple_kbinterpol.f90.
        kb          = kbinterpol(KBWINSZ, KBALPHA)
        flims       = transpose(fpls(1)%frlims)
        dkb01       = real(kb%apod_fast(0.),dp)
        dkb02       = dkb01*dkb01
        dkb03       = dkb01*dkb02
        ! Looping over the un-mirrored asymmetric unit
        noris       = eulspace%get_noris()
        nrefs       = noris/2
        has_mirr    = eulspace%isthere('mirr')
        ! Batch-level particle and reference geometry
        ptcls(:)%active    = .false.
        ptcls(:)%even      = .false.
        ptcls(:)%w         = 0.d0
        ptcls(:)%proj      = 0
        ptcls(:)%self_proj = 0
        !$omp parallel do default(shared) private(i,iptcl) schedule(static) proc_bind(close)
        do i = 1,nptcls
            iptcl = pinds(i)
            if( ptcl_field%get_state(iptcl) == 0 ) cycle
            ptcls(i)%w = real(ptcl_field%get(iptcl, 'w'), dp)
            if( ptcls(i)%w < 1.d-6 ) cycle
            ptcls(i)%active  = .true.
            ptcls(i)%pw_self = ptcls(i)%w * SELFW / dkb02
            ptcls(i)%pw_cl   = ptcls(i)%w / dkb03
            ptcls(i)%even   = ptcl_field%get_eo(iptcl) == 0
            ptcls(i)%proj   = ptcl_field%get_proj(iptcl)
            if( ptcls(i)%proj > nrefs .and. ptcls(i)%proj <= noris )then
                if( has_mirr )then
                    ptcls(i)%proj = eulspace%get_int(ptcls(i)%proj, 'mirr')
                else
                    ptcls(i)%proj = ptcls(i)%proj - nrefs
                endif
            endif
            ptcls(i)%euls = ptcl_field%get_euler(iptcl)
            ptcls(i)%R    = euler2m(ptcls(i)%euls)
        enddo
        !$omp end parallel do
        do iproj = 1,nrefs
            proj_euls          = eulspace%get_euler(iproj)
            Rproj              = euler2m(proj_euls)
            refs(iproj)%tR     = transpose(Rproj)
            refs(iproj)%normal = matmul(zvec, Rproj)
        enddo
        ! Main loop
        nsym = symop%get_nsym()
        do isym = 1,nsym
            call symop%get_sym_rmat(isym, Rsym)
            ptcls(:)%self_proj     = 0
            refs(:)%sym_self_proj  = 0
            refs(:)%sym_done       = .false.
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1,nptcls
                if( .not. ptcls(i)%active ) cycle
                if( isym == 1 )then
                    ptcls(i)%Rsym = ptcls(i)%R
                else
                    ptcls(i)%Rsym = matmul(ptcls(i)%R, Rsym)
                endif
            enddo
            !$omp end parallel do
            do i = 1,nptcls
                if( .not. ptcls(i)%active ) cycle
                srcproj = ptcls(i)%proj
                if( srcproj < 1 .or. srcproj > nrefs ) cycle
                if( isym == 1 )then
                    ptcls(i)%self_proj = srcproj
                else
                    if( .not. refs(srcproj)%sym_done )then
                        normal_ptcl = matmul(refs(srcproj)%normal, Rsym)
                        best_dot = -1.0
                        best_proj = 0
                        do iproj = 1,nrefs
                            dotp = abs(dot_product(refs(iproj)%normal, normal_ptcl))
                            if( dotp > best_dot )then
                                best_dot = dotp
                                best_proj = iproj
                            endif
                        enddo
                        refs(srcproj)%sym_self_proj = best_proj
                        refs(srcproj)%sym_done      = .true.
                    endif
                    ptcls(i)%self_proj = refs(srcproj)%sym_self_proj
                endif
            enddo
            !$omp parallel do default(shared) private(iproj,tRproj)&
            !$omp& private(i,pw,Rptcl,R,euls,psi,phi,pol2cart)&
            !$omp& private(sh,irot,jrot,R2d,hk,rhk,w,ctfsq,fcomp,hh,kk,theta,sin_theta,dang,dCL)&
            !$omp& private(sin_dang,dz,l_even,proj_cl_addr,rot_ptcl,rot_ctfsq,drot,l_self)&
            !$omp& proc_bind(close) schedule(static)
            do iproj = 1,nrefs
                ! Retrieves projection rotation matrix
                tRproj = refs(iproj)%tR
                ! loop over particles
                do i = 1,nptcls
                    if( .not. ptcls(i)%active ) cycle
                    ! e/o flag
                    l_even = ptcls(i)%even
                    ! particle euler angles & rotation matrix
                    Rptcl = ptcls(i)%Rsym
                    l_self = iproj == ptcls(i)%self_proj
                    if( l_self )then
                        ! PARTICLE INSERTION INTO SLICE (cached weight: w * SELFW / dkb02)
                        pw = ptcls(i)%pw_self
                        ! in-plane rotation index offset
                        if( isym == 1 )then
                            psi = ptcls(i)%euls(3)
                        else
                            euls = m2euler(Rptcl)
                            psi  = euls(3)
                        endif
                        drot = self%get_roind_fast(psi)-1
                        ! Loop over the PFT and interpolate
                        do irot = 1, self%pftsz
                            jrot = irot - drot
                            if( jrot < 1 ) jrot = jrot + self%nrots
                            do sh = self%kfromto(1), self%interpklim
                                ! Particle coordinate
                                rhk(:) = real(OSMPL_PAD_FAC) * self%polar(:,sh,jrot)
                                ! NN interpolation with fast KB weight
                                hh = nint(rhk(1))
                                kk = nint(rhk(2))
                                w  = real(kb%apod_fast(rhk(1)-real(hh)),dp) * real(kb%apod_fast(rhk(2)-real(kk)),dp)
                                w  = w * pw
                                if( w > DTINY )then
                                    hh = cyci_1d(flims(:,1), hh)
                                    kk = cyci_1d(flims(:,2), kk)
                                    rot_ptcl(irot,sh)  = PF2 * w * cmplx(fplane_get_cmplx(fpls(i), hh,kk), kind=dp)
                                    rot_ctfsq(irot,sh) =       w * real(fplane_get_ctfsq(fpls(i), hh,kk), dp)
                                else
                                    rot_ptcl(irot,sh)  = DCMPLX_ZERO
                                    rot_ctfsq(irot,sh) = 0.d0
                                endif
                            enddo
                        enddo
                        if( l_even )then
                            self%pfts_even(:,:,iproj) = self%pfts_even(:,:,iproj) + rot_ptcl
                            self%ctf2_even(:,:,iproj) = self%ctf2_even(:,:,iproj) + rot_ctfsq
                        else
                            self%pfts_odd(:,:,iproj)  = self%pfts_odd(:,:,iproj)  + rot_ptcl
                            self%ctf2_odd(:,:,iproj)  = self%ctf2_odd(:,:,iproj)  + rot_ctfsq
                        endif
                    else
                        ! COMMON LINES (cached weight: w / dkb03)
                        pw = ptcls(i)%pw_cl
                        ! Rotation of both planes by transpose of Rproj => the reference is on the (0,0,0))
                        R = matmul(Rptcl, tRproj)
                        ! Euler triplet identification
                        euls = m2euler(R)
                        ! In-plane angle of the particle CL
                        psi = 180.0 - euls(3)
                        ! in-plane angle of the reprojection CL
                        phi = euls(1)
                        if( phi > 180.0 )then
                            ! because we update only the [0;180[ range
                            phi = phi - 180.0
                            psi = psi + 180.0
                        else if( phi < 0.0 )then
                            phi = phi + 180.0
                            psi = psi - 180.0
                        endif
                        ! angle between the particle and reprojection slices
                        theta     = euls(2)
                        sin_theta = sin(deg2rad(theta))
                        ! CL logical coordinates on padded reprojection
                        pol2cart = [sin(deg2rad(phi)), -cos(deg2rad(phi))]
                        do sh = self%kfromto(1), self%interpklim
                            proj_cl_addr(:,sh) = real(sh*OSMPL_PAD_FAC) * pol2cart
                        enddo
                        ! In-plane rotation for mapping reprojection coordinates in particle-space
                        call rotmat2d(180.0+psi-phi, R2D)
                        ! Loop over the PFT and interpolate relevant components
                        do irot = 1, self%pftsz
                            ! angle betwen CL an current reprojection line
                            dang     = modulo(abs(self%angtab(irot) - phi), 180.0)
                            sin_dang = sin(deg2rad(dang))
                            ! dz increases with sh, so if the first shell is too far, all shells are too far
                            if( abs(sin_theta * sin_dang) * real(self%kfromto(1)*OSMPL_PAD_FAC) > DT ) cycle
                            do sh = self%kfromto(1), self%interpklim
                                ! rejects polar points to far from CL
                                hk(1) = real(OSMPL_PAD_FAC) * self%polar(1,sh,irot)
                                if( abs(proj_cl_addr(1,sh) - hk(1)) > DT ) exit
                                hk(2) = real(OSMPL_PAD_FAC) * self%polar(2,sh,irot)
                                if( abs(proj_cl_addr(2,sh) - hk(2)) > DT ) exit
                                ! distance from point to common line sin_dang * sqrt(sum(hk**2)) simplifies to:
                                dCL = sin_dang * real(sh*OSMPL_PAD_FAC)
                                ! relative altitude of the point to slice
                                dz  = sin_theta * dCL
                                ! rejects out-of-plane point
                                if( abs(dz) > DT ) exit
                                ! 2D mapping
                                rhk = matmul(hk, R2D)
                                ! NN interpolation with fast KB weight
                                hh    = nint(rhk(1))
                                kk    = nint(rhk(2))
                                w     = real(kb%apod_fast(rhk(1)-real(hh)),dp) * real(kb%apod_fast(rhk(2)-real(kk)),dp)
                                w     = w * pw * real(kb%apod_fast(dz),dp)
                                if( w < DTINY ) cycle
                                hh    = cyci_1d(flims(:,1), hh)
                                kk    = cyci_1d(flims(:,2), kk)
                                fcomp = PF2 * w * cmplx(fplane_get_cmplx(fpls(i), hh,kk), kind=dp)
                                ctfsq =       w * real(fplane_get_ctfsq(fpls(i),  hh,kk), dp)
                                if( l_even )then
                                    self%pfts_even(irot,sh,iproj) = self%pfts_even(irot,sh,iproj) + fcomp
                                    self%ctf2_even(irot,sh,iproj) = self%ctf2_even(irot,sh,iproj) + ctfsq
                                else
                                    self%pfts_odd(irot,sh,iproj)  = self%pfts_odd(irot,sh,iproj)  + fcomp
                                    self%ctf2_odd(irot,sh,iproj)  = self%ctf2_odd(irot,sh,iproj)  + ctfsq
                                endif
                            enddo   ! shell
                        enddo       ! rotation
                    endif           ! self
                enddo               ! particle
            enddo                   ! slice
            !$omp end parallel do
        enddo                       ! symmetry
    end subroutine polar_cavger_insert_ptcls_direct

    ! Observation-field restoration experiment. The particle planes are first
    ! accumulated into a dense expanded 3D Fourier numerator/density field using
    ! nearest-cell insertion, mirroring the polar=no "insert once, query later"
    ! split. Polar central sections are then gathered from that field and emitted
    ! as the usual even/odd partial-sum payload.
    module subroutine polar_cavger_insert_ptcls_obsfield( self, eulspace, ptcl_field, symop, nptcls, pinds, fpls )
        class(polarft_calc),        intent(inout) :: self
        class(oris),                intent(inout) :: eulspace
        class(oris), pointer,       intent(inout) :: ptcl_field
        class(sym),                 intent(inout) :: symop
        integer,                    intent(in)    :: nptcls, pinds(nptcls)
        class(fplane_type), target, intent(inout) :: fpls(nptcls)
        type(fgrid_obs_field_eo) :: obs
        type(ori) :: o
        complex(dp), allocatable :: pfts_even(:,:,:), pfts_odd(:,:,:)
        real(dp),    allocatable :: ctf2_even(:,:,:), ctf2_odd(:,:,:)
        real(sp) :: hcoords(self%pftsz,self%interpklim-self%kfromto(1)+1)
        real(sp) :: kcoords(self%pftsz,self%interpklim-self%kfromto(1)+1)
        real :: pw
        integer :: lims(3,2), kspan(2), kspan_len, nrefs, noris, i, iptcl, eo, box
        if( nptcls < 1 ) return
        box = self%p_ptr%box_crop
        if( is_even(box) )then
            lims(1,:) = [0, box/2]
            lims(2,:) = [-box/2, box/2-1]
            lims(3,:) = [-box/2, box/2-1]
        else
            lims(1,:) = [0, (box-1)/2]
            lims(2,:) = [-(box-1)/2, (box-1)/2]
            lims(3,:) = [-(box-1)/2, (box-1)/2]
        endif
        noris = eulspace%get_noris()
        if( eulspace%isthere('mirr') )then
            nrefs = noris / 2
        else
            nrefs = noris
        endif
        nrefs = min(nrefs, self%ncls)
        if( nrefs < 1 ) THROW_HARD('no references available; polar_cavger_insert_ptcls_obsfield')
        call obs%new(lims, self%interpklim)
        do i = 1, nptcls
            iptcl = pinds(i)
            if( ptcl_field%get_state(iptcl) == 0 ) cycle
            pw = 1.0
            if( ptcl_field%isthere(iptcl,'w') ) pw = ptcl_field%get(iptcl,'w')
            if( pw < TINY ) cycle
            eo = ptcl_field%get_eo(iptcl)
            call ptcl_field%get_ori(iptcl, o)
            call obs%insert_plane(symop, o, fpls(i), eo, pw)
        enddo
        kspan     = [self%kfromto(1), self%interpklim]
        kspan_len = kspan(2) - kspan(1) + 1
        hcoords   = transpose(self%polar(1,self%kfromto(1):self%interpklim,1:self%pftsz))
        kcoords   = transpose(self%polar(2,self%kfromto(1):self%interpklim,1:self%pftsz))
        allocate(pfts_even(self%pftsz,kspan_len,nrefs), pfts_odd(self%pftsz,kspan_len,nrefs))
        allocate(ctf2_even(self%pftsz,kspan_len,nrefs), ctf2_odd(self%pftsz,kspan_len,nrefs))
        call obs%even%extract_polar(eulspace, nrefs, kspan, hcoords, kcoords, pfts_even, ctf2_even)
        call obs%odd%extract_polar( eulspace, nrefs, kspan, hcoords, kcoords, pfts_odd,  ctf2_odd )
        self%pfts_even(:,kspan(1):kspan(2),1:nrefs) = self%pfts_even(:,kspan(1):kspan(2),1:nrefs) + pfts_even
        self%pfts_odd( :,kspan(1):kspan(2),1:nrefs) = self%pfts_odd( :,kspan(1):kspan(2),1:nrefs) + pfts_odd
        self%ctf2_even(:,kspan(1):kspan(2),1:nrefs) = self%ctf2_even(:,kspan(1):kspan(2),1:nrefs) + ctf2_even
        self%ctf2_odd( :,kspan(1):kspan(2),1:nrefs) = self%ctf2_odd( :,kspan(1):kspan(2),1:nrefs) + ctf2_odd
        deallocate(pfts_even, pfts_odd, ctf2_even, ctf2_odd)
        call obs%kill
        call o%kill
    end subroutine polar_cavger_insert_ptcls_obsfield

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

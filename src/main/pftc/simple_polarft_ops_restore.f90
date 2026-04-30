!@descr: submodule for class average restoration in the polar Fourier domain
submodule (simple_polarft_calc) simple_polarft_ops_restore
implicit none
#include "simple_local_flags.inc"
contains               

    !>  \brief  Filter references
    !! keep serial
    module subroutine polar_filterrefs( self, icls, filter )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: icls
        real,                intent(in)    :: filter(:)
        integer :: k
        real    :: fk
        do k = self%kfromto(1),self%interpklim
            fk = filter(k)
            self%pfts_merg(:,k,icls) = fk * self%pfts_merg(:,k,icls)
        enddo
        do k = self%kfromto(1),self%interpklim
            fk = filter(k)
            self%pfts_even(:,k,icls) = fk * self%pfts_even(:,k,icls)
        enddo
        do k = self%kfromto(1),self%interpklim
            fk = filter(k)
            self%pfts_odd(:,k,icls)  = fk * self%pfts_odd(:,k,icls)
        enddo
    end subroutine polar_filterrefs

    module subroutine polar_cavger_calc_frc( self, pft1, pft2, n, frc )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%interpklim), pft2(self%pftsz,self%kfromto(1):self%interpklim)
        integer,             intent(in)    :: n
        real(sp),            intent(inout) :: frc(1:n)
        real(dp) :: var1, var2, denom
        integer  :: k
        frc(1:self%kfromto(1)-1) = 0.999
        do k = self%kfromto(1), self%interpklim
            var1  = sum(csq_fast(pft1(:,k)))
            var2  = sum(csq_fast(pft2(:,k)))
            if( (var1>DTINY) .and. (var2>DTINY) )then
                denom  = sqrt(var1*var2)
                frc(k) = real(sum(pft1(:,k)*conjg(pft2(:,k))) / denom, sp)
            else
                frc(k) = 0.0
            endif
        enddo
        if( self%interpklim < n ) frc(self%interpklim+1:) = 0.0
    end subroutine polar_cavger_calc_frc

    ! 3D SECTION

    ! 3D SECTION - Polar references are generated from common lines of the polar sums of particles

    module subroutine polar_cavger_normalize_commonline_refs( self, reforis, symop, cline, update_frac )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        type(sym),           intent(in)    :: symop
        type(cmdline),       intent(in)    :: cline
        real,                intent(in)    :: update_frac
        type(class_frcs)         :: frcs
        complex(dp), allocatable :: prev_even(:,:,:), prev_odd(:,:,:)
        complex(dp) :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp) :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp)    :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp)    :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp)    :: fsc(self%kfromto(1):self%interpklim), fsc_state(self%kfromto(1):self%interpklim), ufrac_trec
        real        :: fsc_boxcrop(1:fdim(self%p_ptr%box_crop)-1)
        integer     :: i, state, base, nprojs
        ! Mirror Fourier & CTF2 slices
        call self%mirror_slices( reforis )
        ! Common-lines contribution
        call self%calc_comlin_contrib(reforis, symop, pfts_even, pfts_odd, ctf2_even, ctf2_odd)
        ! e/o trailing reconstruction part 1
        if( self%p_ptr%l_trail_rec )then
            ufrac_trec = real(merge(self%p_ptr%ufrac_trec ,update_frac , cline%defined('ufrac_trec')),dp)
            ! read previous polar references to be used first for FSC calculation
            call prepare_trail_rec_arrays( self, reforis, prev_even, prev_odd )
            ! The per-state FSCs below are written to disk and used for ML
            ! regularization. This full-range FSC is retained only for the
            ! legacy even/odd registration cutoff after restoration.
            call self%calc_fsc(pfts_even, pfts_odd, ctf2_even, ctf2_odd, fsc,&
                                &ufrac_trec=ufrac_trec, prev_even=prev_even, prev_odd=prev_odd)
        else
            call self%calc_fsc(pfts_even, pfts_odd, ctf2_even, ctf2_odd, fsc)
        endif
        nprojs = self%p_ptr%nspace
        call frcs%new(nprojs, self%p_ptr%box_crop, self%p_ptr%smpd_crop, self%p_ptr%nstates)
        do state = 1,self%p_ptr%nstates
            base = (state - 1) * nprojs
            if( self%p_ptr%l_trail_rec )then
                call calc_fsc_range(self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, base+1, base+nprojs, fsc_state, &
                    &ufrac_trec=ufrac_trec, prev_even=prev_even, prev_odd=prev_odd)
            else
                call calc_fsc_range(self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, base+1, base+nprojs, fsc_state)
            endif
            fsc_boxcrop(                 :self%kfromto(1)) = 1.0
            fsc_boxcrop(self%kfromto(1)  :self%interpklim) = real(fsc_state(self%kfromto(1):self%interpklim))
            if( self%interpklim < size(fsc_boxcrop) )then
                fsc_boxcrop(self%interpklim+1:)            = 0.0
            endif
            call arr2file(fsc_boxcrop, string(FSC_FBODY//int2str_pad(state,2)//BIN_EXT))
            do i = 1,nprojs
                ! FRCs are set to the state-local FSC. to check if we are using those
                call frcs%set_frc(i, fsc_boxcrop, state)
            enddo
            if( self%p_ptr%l_ml_reg ) call add_invtausq2rho_range(self, ctf2_even, ctf2_odd, base+1, base+nprojs, fsc_state)
        enddo
        call frcs%write(string(FRCS_FILE))
        call frcs%kill
        fsc_boxcrop(                 :self%kfromto(1)) = 1.0
        fsc_boxcrop(self%kfromto(1)  :self%interpklim) = real(fsc(self%kfromto(1):self%interpklim))
        if( self%interpklim < size(fsc_boxcrop) )then
            fsc_boxcrop(self%interpklim+1:)            = 0.0
        endif
        ! Wiener normalization
        call self%restore_references(reforis, pfts_even, pfts_odd, ctf2_even, ctf2_odd)
        call average_lowres_eo_for_docking(self, fsc_boxcrop)
        ! e/o trailing reconstruction part 2
        if( self%p_ptr%l_trail_rec )then
            call finalize_trail_rec( self, ufrac_trec, prev_even, prev_odd )
        endif
    end subroutine polar_cavger_normalize_commonline_refs

    ! Deals with summing slices and their mirror
    module subroutine mirror_slices( self, ref_space )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: ref_space
        complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
        real        :: psi
        integer     :: iref, m, iloc, mloc, iglob, base, istate, nprojs, nprojs_nomirr, nwork
        logical     :: l_rotm
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        nprojs        = self%p_ptr%nspace
        nprojs_nomirr = nprojs / 2
        if( nprojs * self%p_ptr%nstates /= self%ncls )then
            THROW_HARD('state-major polar reference count mismatch; mirror_slices')
        endif
        nwork = nprojs_nomirr * self%p_ptr%nstates
        !$omp parallel do default(shared) proc_bind(close) private(iref,m,iloc,mloc,iglob,base,istate,pft,ctf2,psi,l_rotm)
        do iglob = 1,nwork
            istate = (iglob - 1) / nprojs_nomirr + 1
            iloc   = mod(iglob - 1, nprojs_nomirr) + 1
            base   = (istate - 1) * nprojs
            iref   = base + iloc
            mloc   = ref_space%get_int(iloc,'mirr')
            if( mloc < 1 .or. mloc > nprojs ) THROW_HARD('mirror index out of state-local range; mirror_slices')
            m      = base + mloc
            psi    = abs(ref_space%get(mloc, 'psi'))
            l_rotm = (psi > 0.1) .and. (psi < 359.9)
            ! Fourier components
            if( l_rotm )then
                call self%mirror_pft(self%pfts_even(:,:,m), pft)
                self%pfts_even(:,:,iref) = self%pfts_even(:,:,iref) + conjg(pft)
                call self%mirror_pft(self%pfts_odd(:,:,m), pft)
                self%pfts_odd(:,:,iref)  = self%pfts_odd(:,:,iref)  + conjg(pft)
                call self%mirror_pft(conjg(self%pfts_even(:,:,iref)), self%pfts_even(:,:,m))
                call self%mirror_pft(conjg(self%pfts_odd(:,:,iref)),  self%pfts_odd(:,:,m))
            else
                call self%mirror_pft(self%pfts_even(:,:,m), pft)
                self%pfts_even(:,:,iref) = self%pfts_even(:,:,iref) + pft
                call self%mirror_pft(self%pfts_odd(:,:,m), pft)
                self%pfts_odd(:,:,iref)  = self%pfts_odd(:,:,iref)  + pft
                call self%mirror_pft(self%pfts_even(:,:,iref), self%pfts_even(:,:,m))
                call self%mirror_pft(self%pfts_odd(:,:,iref),  self%pfts_odd(:,:,m))
            endif
            ! CTF
            call self%mirror_ctf2(self%ctf2_even(:,:,m), ctf2)
            self%ctf2_even(:,:,iref) = self%ctf2_even(:,:,iref) + ctf2
            call self%mirror_ctf2(self%ctf2_odd(:,:,m),  ctf2)
            self%ctf2_odd(:,:,iref)  = self%ctf2_odd(:,:,iref)  + ctf2
            call self%mirror_ctf2(self%ctf2_even(:,:,iref), self%ctf2_even(:,:,m))
            call self%mirror_ctf2(self%ctf2_odd(:,:,iref),  self%ctf2_odd(:,:,m))
        enddo
        !$omp end parallel do
    end subroutine mirror_slices

    ! Calculate common-lines contributions from all the slices
    module subroutine calc_comlin_contrib( self, ref_space, symop, pfts_cl_even, pfts_cl_odd, ctf2_cl_even, ctf2_cl_odd )
        class(polarft_calc), intent(in)    :: self
        type(oris),          intent(in)    :: ref_space
        type(sym),           intent(in)    :: symop
        complex(kind=dp),    intent(inout) :: pfts_cl_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(kind=dp),    intent(inout) :: pfts_cl_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp), allocatable :: Rsym(:,:,:)
        complex(dp) :: pft_slice_e(self%kfromto(1):self%interpklim,self%pftsz)
        complex(dp) :: pft_slice_o(self%kfromto(1):self%interpklim,self%pftsz)
        complex(dp) :: cl_e(self%kfromto(1):self%interpklim), cl_o(self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2_slice_e(self%kfromto(1):self%interpklim,self%pftsz)
        real(dp)    :: ctf2_slice_o(self%kfromto(1):self%interpklim,self%pftsz)
        real(dp)    :: rl_e(self%kfromto(1):self%interpklim), rl_o(self%kfromto(1):self%interpklim)
        real(dp)    :: wl(self%kfromto(1):self%interpklim), wr(self%kfromto(1):self%interpklim)
        real(dp), allocatable :: R(:,:,:)
        real(dp)    :: Rj(3,3), tRi(3,3), eulers(3), psi
        real        :: Rtmp(3,3)
        integer     :: rotl, rotr, iref, jref, m, isym, nsym
        integer     :: iloc, jloc, mloc, iglob, base, istate, nprojs, nprojs_nomirr, nwork
        logical     :: l_rotm, conjgl, conjgr
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        ! Symmetry rotation matrices
        nsym          = symop%get_nsym()
        nprojs        = self%p_ptr%nspace
        nprojs_nomirr = nprojs / 2
        if( nprojs * self%p_ptr%nstates /= self%ncls )then
            THROW_HARD('state-major polar reference count mismatch; calc_comlin_contrib')
        endif
        nwork = nprojs_nomirr * self%p_ptr%nstates
        allocate(R(3,3,nprojs),source=0.d0)
        allocate(Rsym(3,3,nsym),source=0.d0)
        !$omp parallel default(shared) proc_bind(close)&
        !$omp private(iref,jref,m,isym,tRi,Rj,Rtmp,eulers,wl,wr,psi,l_rotm,cl_e,cl_o,rl_e,rl_o)&
        !$omp& private(iloc,jloc,mloc,iglob,base,istate)&
        !$omp& private(rotl,rotr,conjgl,conjgr,pft_slice_e,pft_slice_o,ctf2_slice_e,ctf2_slice_o)
        ! Init
        !$omp workshare
        pfts_cl_even = DCMPLX_ZERO
        pfts_cl_odd  = DCMPLX_ZERO
        ctf2_cl_even = 0.d0
        ctf2_cl_odd  = 0.d0
        !$omp end workshare
        ! Caching rotation matrices
        !$omp do schedule(static)
        do iref = 1,nprojs
            R(:,:,iref) = real(ref_space%get_mat(iref),dp)
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do isym = 1,nsym
            call symop%get_sym_rmat(isym, Rtmp)
            Rsym(:,:,isym) = real(Rtmp,dp)
        end do
        !$omp end do
        ! Common lines are restricted within a state; cross-state common lines
        ! are not physically meaningful.
        !$omp do schedule(static)
        do iglob = 1,nwork
            istate = (iglob - 1) / nprojs_nomirr + 1
            iloc   = mod(iglob - 1, nprojs_nomirr) + 1
            base   = (istate - 1) * nprojs
            iref   = base + iloc
            mloc   = ref_space%get_int(iloc,'mirr')
            if( mloc < 1 .or. mloc > nprojs ) THROW_HARD('mirror index out of state-local range; calc_comlin_contrib')
            m      = base + mloc
            tRi = transpose(R(:,:,iloc))
            pft_slice_e  = DCMPLX_ZERO
            pft_slice_o  = DCMPLX_ZERO
            ctf2_slice_e = 0.d0
            ctf2_slice_o = 0.d0
            do jloc = 1,nprojs_nomirr
                jref = base + jloc
                do isym = 1,nsym
                    if( isym == 1 )then
                        if( jref == iref ) cycle    ! self   exclusion
                        if( jref == m )    cycle    ! mirror exclusion
                        ! Rotation of both planes by transpose of Ri (tRixRi=I & Rsym=I)
                        Rj = matmul(R(:,:,jloc), tRi)
                    else
                        ! Symmetry operator
                        Rj = matmul(R(:,:,jloc), Rsym(:,:,isym))
                        ! Rotation of both planes by transpose of Ri (tRixRi -> I)
                        Rj = matmul(Rj, tRi)
                    endif
                    ! Euler angles identification
                    eulers = dm2euler(Rj)
                    ! in plane rotation angle of jref slice intersecting iref
                    psi = 360.d0 - eulers(3)
                    ! get the weights, rotation indices and compute the interpolated common line
                    call self%gen_clin_weights(psi, rotl, rotr, wl, wr)
                    call update_index(rotl, conjgl)
                    call update_index(rotr, conjgr)
                    if( conjgl )then
                        cl_e = wl * conjg(self%pfts_even(rotl,:,jref))
                        cl_o = wl * conjg(self%pfts_odd(rotl,:,jref))
                    else
                        cl_e = wl * self%pfts_even(rotl,:,jref)
                        cl_o = wl * self%pfts_odd(rotl,:,jref)
                    endif
                    if( conjgr )then
                        cl_e = cl_e + wr * conjg(self%pfts_even(rotr,:,jref))
                        cl_o = cl_o + wr * conjg(self%pfts_odd(rotr,:,jref))
                    else
                        cl_e = cl_e + wr * self%pfts_even(rotr,:,jref)
                        cl_o = cl_o + wr * self%pfts_odd(rotr,:,jref)
                    endif
                    rl_e = wl * self%ctf2_even(rotl,:,jref) + wr * self%ctf2_even(rotr,:,jref)
                    rl_o = wl * self%ctf2_odd(rotl,:,jref)  + wr * self%ctf2_odd(rotr,:,jref)
                    ! in plane rotation angle of iref slice
                    psi = eulers(1)
                    ! get the weights, rotation indices
                    call self%gen_clin_weights(psi, rotl, rotr, wl, wr)
                    call update_index(rotl, conjgl)
                    call update_index(rotr, conjgr)
                    ! interpolate intersecting line into iref plane
                    if( conjgl )then
                        pft_slice_e(:,rotl)  = pft_slice_e(:,rotl)  + wl * conjg(cl_e)
                        pft_slice_o(:,rotl)  = pft_slice_o(:,rotl)  + wl * conjg(cl_o)
                    else
                        pft_slice_e(:,rotl)  = pft_slice_e(:,rotl)  + wl * cl_e
                        pft_slice_o(:,rotl)  = pft_slice_o(:,rotl)  + wl * cl_o
                    endif
                    if( conjgr )then
                        pft_slice_e(:,rotr)  = pft_slice_e(:,rotr)  + wr * conjg(cl_e)
                        pft_slice_o(:,rotr)  = pft_slice_o(:,rotr)  + wr * conjg(cl_o)
                    else
                        pft_slice_e(:,rotr)  = pft_slice_e(:,rotr)  + wr * cl_e
                        pft_slice_o(:,rotr)  = pft_slice_o(:,rotr)  + wr * cl_o
                    endif
                    ctf2_slice_e(:,rotl) = ctf2_slice_e(:,rotl) + wl * rl_e
                    ctf2_slice_e(:,rotr) = ctf2_slice_e(:,rotr) + wr * rl_e
                    ctf2_slice_o(:,rotl) = ctf2_slice_o(:,rotl) + wl * rl_o
                    ctf2_slice_o(:,rotr) = ctf2_slice_o(:,rotr) + wr * rl_o
                enddo
            enddo
            ! memory layout
            pfts_cl_even(:,:,iref) = transpose(pft_slice_e)
            pfts_cl_odd(:,:,iref)  = transpose(pft_slice_o)
            ctf2_cl_even(:,:,iref) = transpose(ctf2_slice_e)
            ctf2_cl_odd(:,:,iref)  = transpose(ctf2_slice_o)
        enddo
        !$omp end do
        ! Mirroring contributions
        !$omp do schedule(static)
        do iglob = 1,nwork
            istate = (iglob - 1) / nprojs_nomirr + 1
            iloc   = mod(iglob - 1, nprojs_nomirr) + 1
            base   = (istate - 1) * nprojs
            iref   = base + iloc
            mloc   = ref_space%get_int(iloc,'mirr')
            if( mloc < 1 .or. mloc > nprojs ) THROW_HARD('mirror index out of state-local range; calc_comlin_contrib')
            m      = base + mloc
            psi    = abs(ref_space%get(mloc, 'psi'))
            l_rotm = (psi>0.1) .and. (psi<359.9)
            call self%mirror_pft(pfts_cl_even(:,:,iref), pfts_cl_even(:,:,m))
            call self%mirror_pft(pfts_cl_odd(:,:,iref),  pfts_cl_odd(:,:,m))
            if( l_rotm )then
                pfts_cl_even(:,:,m) = conjg(pfts_cl_even(:,:,m))
                pfts_cl_odd(:,:,m) = conjg(pfts_cl_odd(:,:,m))
            endif
            call self%mirror_ctf2(ctf2_cl_even(:,:,iref), ctf2_cl_even(:,:,m))
            call self%mirror_ctf2(ctf2_cl_odd(:,:,iref),  ctf2_cl_odd(:,:,m))
        enddo
        !$omp end do
        !$omp end parallel
        deallocate(R, Rsym)
    contains

        pure subroutine update_index( rot, conjugate )
            integer, intent(inout)  :: rot
            logical, intent(out) :: conjugate
            integer :: irot
            irot = merge(rot-self%nrots, rot, rot>self%nrots)
            if( irot < 1 )then
                irot      = irot + self%pftsz
                conjugate = .true.
            elseif( irot > self%pftsz )then
                irot      = irot - self%pftsz
                conjugate = .true.
            else
                conjugate = .false.
            endif
            rot = irot
        end subroutine update_index

    end subroutine calc_comlin_contrib

    ! 3D SECTION - Polar references are sampled from assembled Cartesian observation fields

    module subroutine polar_cavger_normalize_obsfield_refs( self, reforis, cline, update_frac )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        type(cmdline),       intent(in)    :: cline
        real,                intent(in)    :: update_frac
        type(class_frcs)         :: frcs
        complex(dp), allocatable :: prev_even(:,:,:), prev_odd(:,:,:)
        real(dp)    :: fsc(self%kfromto(1):self%interpklim), fsc_state(self%kfromto(1):self%interpklim)
        real(dp)    :: invtau2_even(self%kfromto(1):self%interpklim), invtau2_odd(self%kfromto(1):self%interpklim)
        real(sp)    :: hcoords(self%pftsz,self%interpklim-self%kfromto(1)+1)
        real(sp)    :: kcoords(self%pftsz,self%interpklim-self%kfromto(1)+1)
        real(dp)    :: ufrac_trec
        real        :: fsc_boxcrop(1:fdim(self%p_ptr%box_crop)-1)
        integer     :: i, state, base, nprojs, nrefs, noris, prior_start, kspan(2)
        integer     :: prev_pftsz, prev_nrefs
        logical     :: have_prev_refs, need_prev_refs
        if( .not. allocated(self%obsfields) ) THROW_HARD('obsfields are not allocated; polar_cavger_normalize_obsfield_refs')
        nprojs = self%p_ptr%nspace
        if( mod(nprojs,2) /= 0 )then
            THROW_HARD('obsfield polar reference extraction requires an even nspace for half-grid mirroring')
        endif
        if( .not. reforis%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space; polar_cavger_normalize_obsfield_refs')
        endif
        noris = reforis%get_noris()
        if( noris < nprojs )then
            THROW_HARD('not enough reference orientations for obsfield polar extraction')
        endif
        nrefs = nprojs / 2
        need_prev_refs = self%p_ptr%l_ml_reg .or. self%p_ptr%l_trail_rec
        have_prev_refs = need_prev_refs .and. file_exists(POLAR_REFS_FBODY//'_even'//BIN_EXT) .and. &
            &file_exists(POLAR_REFS_FBODY//'_odd'//BIN_EXT)
        if( have_prev_refs )then
            call self%get_pft_array_dims(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_pftsz, nrefs=prev_nrefs)
            have_prev_refs = (prev_pftsz == self%pftsz) .and. (prev_nrefs <= self%ncls) .and. &
                &(mod(self%ncls,self%p_ptr%nstates) == 0) .and. (mod(prev_nrefs,self%p_ptr%nstates) == 0)
        endif
        if( self%p_ptr%l_trail_rec .and. (.not. have_prev_refs) )then
            THROW_HARD('compatible previous polar references are required for obsfield trailing reconstruction')
        endif
        ! Previous restored references provide the FSC/SSNR estimate, as in the
        ! Cartesian route. Fractional updates trail the restored current refs
        ! against that previous pair; the current obsfield density still sets
        ! the ML denominator scale for this update.
        if( have_prev_refs )then
            call prepare_trail_rec_arrays( self, reforis, prev_even, prev_odd )
            ! Global FSC is retained only for the low-res even/odd docking cutoff;
            ! ML regularization uses the state-local FSC below.
            call calc_fsc_from_refs(self, prev_even, prev_odd, 1, self%ncls, fsc)
        else
            fsc = 0.d0
        endif
        if( self%p_ptr%l_trail_rec )then
            ufrac_trec = real(merge(self%p_ptr%ufrac_trec ,update_frac , cline%defined('ufrac_trec')),dp)
        endif
        ! write down FRCs
        kspan  = [self%kfromto(1), self%interpklim]
        if( self%p_ptr%l_ml_reg )then
            prior_start = ml_prior_start(self)
        else
            prior_start = self%interpklim + 1
        endif
        hcoords = transpose(self%polar(1,self%kfromto(1):self%interpklim,1:self%pftsz))
        kcoords = transpose(self%polar(2,self%kfromto(1):self%interpklim,1:self%pftsz))
        call frcs%new(nprojs, self%p_ptr%box_crop, self%p_ptr%smpd_crop, self%p_ptr%nstates)
        do state = 1,self%p_ptr%nstates
            base = (state - 1) * nprojs
            if( have_prev_refs )then
                call calc_fsc_from_refs(self, prev_even, prev_odd, base+1, base+nprojs, fsc_state)
            else
                fsc_state = 0.d0
            endif
            fsc_boxcrop(                 :self%kfromto(1)) = 1.0
            fsc_boxcrop(self%kfromto(1)  :self%interpklim) = real(fsc_state(self%kfromto(1):self%interpklim))
            if( self%interpklim < size(fsc_boxcrop) )then
                fsc_boxcrop(self%interpklim+1:)            = 0.0
            endif
            call arr2file(fsc_boxcrop, string(FSC_FBODY//int2str_pad(state,2)//BIN_EXT))
            do i = 1,nprojs
                ! FRCs are set to the state-local FSC. to check if we are using those
                call frcs%set_frc(i, fsc_boxcrop, state)
            enddo
            invtau2_even = 0.d0
            invtau2_odd  = 0.d0
            if( self%p_ptr%l_ml_reg )then
                call self%obsfields(state)%even%calc_invtau2(kspan, fsc_state, real(self%p_ptr%tau,dp), &
                    &prior_start, invtau2_even)
                call self%obsfields(state)%odd%calc_invtau2( kspan, fsc_state, real(self%p_ptr%tau,dp), &
                    &prior_start, invtau2_odd )
            endif
            call self%obsfields(state)%extract_polar(reforis, nrefs, kspan, hcoords, kcoords, &
                &invtau2_even, invtau2_odd, prior_start, &
                &self%pfts_even(:,kspan(1):kspan(2),base+1:base+nrefs), &
                &self%pfts_odd( :,kspan(1):kspan(2),base+1:base+nrefs), &
                &self%pfts_merg(:,kspan(1):kspan(2),base+1:base+nrefs))
        enddo
        call mirror_slices_obsfield( self, reforis )
        call frcs%write(string(FRCS_FILE))
        call frcs%kill
        fsc_boxcrop(                 :self%kfromto(1)) = 1.0
        fsc_boxcrop(self%kfromto(1)  :self%interpklim) = real(fsc(self%kfromto(1):self%interpklim))
        if( self%interpklim < size(fsc_boxcrop) )then
            fsc_boxcrop(self%interpklim+1:)            = 0.0
        endif
        call average_lowres_eo_for_docking(self, fsc_boxcrop)
        ! e/o trailing reconstruction part 2
        if( self%p_ptr%l_trail_rec .and. have_prev_refs )then
            call finalize_trail_rec( self, ufrac_trec, prev_even, prev_odd )
        endif
        if( allocated(prev_even) ) deallocate(prev_even)
        if( allocated(prev_odd)  ) deallocate(prev_odd)
    end subroutine polar_cavger_normalize_obsfield_refs

    subroutine average_lowres_eo_for_docking( self, fsc_boxcrop )
        class(polarft_calc), intent(inout) :: self
        real,                intent(in)    :: fsc_boxcrop(:)
        integer :: find4eoavg
        ! Only this low-resolution range is averaged into both half maps.
        ! Higher-frequency even/odd references remain independently restored.
        find4eoavg = max(K4EOAVGLB, calc_fourier_index(FREQ4EOAVG3D, self%p_ptr%box, self%p_ptr%smpd))
        find4eoavg = min(find4eoavg, get_find_at_crit(fsc_boxcrop, FSC4EOAVG3D))
        find4eoavg = min(self%interpklim, find4eoavg)
        if( find4eoavg < self%kfromto(1) ) return
        !$omp parallel workshare proc_bind(close)
        self%pfts_even(:,self%kfromto(1):find4eoavg,:) = self%pfts_merg(:,self%kfromto(1):find4eoavg,:)
        self%pfts_odd(:, self%kfromto(1):find4eoavg,:) = self%pfts_merg(:,self%kfromto(1):find4eoavg,:)
        !$omp end parallel workshare
    end subroutine average_lowres_eo_for_docking

    !>  \brief Generate the mirror slices for restored obsfield references
    !!         iref runs through the unique half (un-mirrored references), the
    !!         iref-th slice is final and after mirroring overwrites corresponding
    !!         mirror reference (m-th).
    module subroutine mirror_slices_obsfield( self, ref_space )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: ref_space
        complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
        real        :: psi
        integer     :: iref, m, iloc, mloc, iglob, base, istate, nprojs, nprojs_nomirr, nwork
        logical     :: l_rotm
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        nprojs        = self%p_ptr%nspace
        nprojs_nomirr = nprojs / 2
        if( nprojs * self%p_ptr%nstates /= self%ncls )then
            THROW_HARD('state-major polar reference count mismatch; mirror_slices_obsfield')
        endif
        nwork = nprojs_nomirr * self%p_ptr%nstates
        !$omp parallel do default(shared) proc_bind(close) private(iref,m,iloc,mloc,iglob,base,istate,pft,psi,l_rotm)
        do iglob = 1,nwork
            istate = (iglob - 1) / nprojs_nomirr + 1
            iloc   = mod(iglob - 1, nprojs_nomirr) + 1
            base   = (istate - 1) * nprojs
            iref   = base + iloc
            mloc   = ref_space%get_int(iloc,'mirr')
            if( mloc < 1 .or. mloc > nprojs ) THROW_HARD('mirror index out of state-local range; mirror_slices_obsfield')
            m      = base + mloc
            psi    = abs(ref_space%get(mloc, 'psi'))
            l_rotm = (psi > 0.1) .and. (psi < 359.9)
            ! Fourier components
            if( l_rotm )then
                call self%mirror_pft(self%pfts_merg(:,:,iref), pft)
                self%pfts_merg(:,:,m) = conjg(pft)
                call self%mirror_pft(self%pfts_even(:,:,iref), pft)
                self%pfts_even(:,:,m) = conjg(pft)
                call self%mirror_pft(self%pfts_odd(:,:,iref), pft)
                self%pfts_odd(:,:,m)  = conjg(pft)
            else
                call self%mirror_pft(self%pfts_merg(:,:,iref), self%pfts_merg(:,:,m))
                call self%mirror_pft(self%pfts_even(:,:,iref), self%pfts_even(:,:,m))
                call self%mirror_pft(self%pfts_odd(:,:,iref),  self%pfts_odd(:,:,m))
            endif
        enddo
        !$omp end parallel do
    end subroutine mirror_slices_obsfield

    !>  \brief local private routine for trailing reconstruction weighing
    !!         and dimension checking/editing for polar modes
    subroutine prepare_trail_rec_arrays( self, reforis, prev_even, prev_odd )
        class(polarft_calc),       intent(in) :: self
        type(oris),                intent(in) :: reforis
        complex(dp), allocatable, intent(out) :: prev_even(:,:,:), prev_odd(:,:,:)
        type(string) :: fname
        type(sym)    :: symop
        type(oris)   :: prev_eulspace
        integer      :: prev_pftsz, prev_klims(2), prev_nrefs, ks, ke
        integer      :: current_nspace, prev_nspace
        logical      :: l_pad
        call self%get_pft_array_dims(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_pftsz, prev_klims, prev_nrefs)
        ! PFT arrays are state-major blocks of per-projection sections.
        ! The angular shell size must match; frequency ranges can be padded.
        ! When per-state nspace grows, each current projection is sourced from
        ! the closest previous projection in the same state.
        if( (prev_nrefs == self%ncls) .and. (prev_pftsz == self%pftsz) )then
            call self%read_pft_array(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_even)
            call self%read_pft_array(string(POLAR_REFS_FBODY)//'_odd'//BIN_EXT,  prev_odd)
        else
            if( prev_pftsz /= self%pftsz )then
                THROW_HARD('Previous refs have different pftsz, cannot be used for trailing reconstruction')
            endif
            if( prev_nrefs > self%ncls )then
                THROW_HARD('Previous per-state nspace is larger than current, cannot remap trailing reconstruction refs')
            else
                ! Previous refs have different dimensions and need to be mapped to current space
                if( mod(self%ncls, self%p_ptr%nstates) /= 0 )then
                    THROW_HARD('Current state-major reference count is not divisible by nstates')
                endif
                if( mod(prev_nrefs, self%p_ptr%nstates) /= 0 )then
                    THROW_HARD('Previous state-major reference count is not divisible by nstates')
                endif
                current_nspace = self%ncls / self%p_ptr%nstates
                prev_nspace    = prev_nrefs / self%p_ptr%nstates
                call symop%new(self%p_ptr%pgrp)
                call prev_eulspace%new(prev_nspace, is_ptcl=.false.)
                call symop%build_refspiral(prev_eulspace)
                ks = max(self%kfromto(1), prev_klims(1))
                ke = min(self%interpklim, prev_klims(2))
                ! Compare full requested/stored ranges, not just intersection
                ! bounds, so missing low/high frequencies trigger zero padding.
                l_pad = (self%kfromto(1) /= prev_klims(1)) .or. (self%interpklim /= prev_klims(2))
                fname = string(POLAR_REFS_FBODY)//'_even'//BIN_EXT
                call map_current_refs_to_closest_previous_refs(fname, prev_even)
                fname = string(POLAR_REFS_FBODY)//'_odd'//BIN_EXT
                call map_current_refs_to_closest_previous_refs(fname, prev_odd)
                call symop%kill
                call prev_eulspace%kill
                call fname%kill
            endif
        endif
        contains

            subroutine map_current_refs_to_closest_previous_refs( fname, array )
                type(string),             intent(in)  :: fname
                complex(dp), allocatable, intent(out) :: array(:,:,:)
                complex(sp), allocatable :: src(:,:,:)
                type(ori) :: o
                integer   :: state, iproj, i, j, cur_base, prev_base, prev_ref
                call self%read_any_pft_array(fname, src)
                allocate(array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls), source=DCMPLX_ZERO)
                !$omp parallel do default(shared) proc_bind(close) private(i,o,j,state,iproj,cur_base,prev_base,prev_ref) &
                !$omp& schedule(static)
                do i = 1,self%ncls
                    state = (i - 1) / current_nspace + 1
                    iproj = i - (state - 1) * current_nspace
                    cur_base  = (state - 1) * current_nspace
                    prev_base = (state - 1) * prev_nspace
                    call reforis%get_ori(iproj, o)
                    j = prev_eulspace%find_closest_proj(o)
                    prev_ref = prev_base + j
                    if( l_pad )then
                        if( ks <= ke ) array(:,ks:ke,cur_base+iproj) = cmplx(src(:,ks:ke,prev_ref), kind=dp)
                    else
                        array(:,:,cur_base+iproj) = cmplx(src(:,:,prev_ref), kind=dp)
                    endif
                enddo
                !$omp end parallel do
                deallocate(src)
                call o%kill
            end subroutine map_current_refs_to_closest_previous_refs

    end subroutine prepare_trail_rec_arrays

    ! 3D SECTION - Common routines

    !>  \brief local private routine for trailing reconstruction  weighing and merging
    subroutine finalize_trail_rec( self, ufrac_trec, prev_even, prev_odd )
        class(polarft_calc),      intent(inout) :: self
        real(dp),                 intent(in)    :: ufrac_trec
        complex(dp), allocatable, intent(inout) :: prev_even(:,:,:)
        complex(dp), allocatable, intent(inout) :: prev_odd(:,:,:)
        type(string) :: fname
        integer      :: tmp_pftsz, tmp_kfromto(2), tmp_nrefs
        logical      :: have_merged_prev
        ! adding weighted previous refs after restoration of current refs
        !$omp parallel workshare proc_bind(close)
        self%pfts_even = ufrac_trec * self%pfts_even + (1.d0-ufrac_trec) * prev_even
        self%pfts_odd  = ufrac_trec * self%pfts_odd  + (1.d0-ufrac_trec) * prev_odd
        !$omp end parallel workshare
        fname = string(POLAR_REFS_FBODY)//BIN_EXT
        have_merged_prev = file_exists(fname)
        if( have_merged_prev )then
            call self%get_pft_array_dims(fname, tmp_pftsz, tmp_kfromto, tmp_nrefs)
            have_merged_prev = (tmp_nrefs == self%ncls) .and. (tmp_pftsz == self%pftsz)
        endif
        if( .not. have_merged_prev )then
            ! Previous handoff may be e/o-only. The halves have already been
            ! trailed independently above; derive only the merged companion.
            !$omp parallel workshare proc_bind(close)
            self%pfts_merg = 0.5d0 * (self%pfts_even + self%pfts_odd)
            !$omp end parallel workshare
            deallocate(prev_even, prev_odd)
        else
            deallocate(prev_even,prev_odd)
            call self%read_pft_array(fname, prev_even)
            !$omp parallel workshare proc_bind(close)
            self%pfts_merg = ufrac_trec * self%pfts_merg + (1.d0-ufrac_trec) * prev_even
            !$omp end parallel workshare
            deallocate(prev_even)
        endif
        call fname%kill
    end subroutine finalize_trail_rec

    ! Calculate global FSC within [self%kfromto(1);self%interpklim]
    module subroutine calc_fsc( self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, fsc, ufrac_trec, prev_even, prev_odd )
        class(polarft_calc),   intent(in)  :: self
        complex(dp),           intent(in)  :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp),           intent(in)  :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(in)  :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(in)  :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(out) :: fsc(self%kfromto(1):self%interpklim)
        real(dp),    optional, intent(in)  :: ufrac_trec
        complex(dp), optional, intent(in)  :: prev_even(:,:,:)
        complex(dp), optional, intent(in)  :: prev_odd(:,:,:)
        if( present(ufrac_trec) .and. present(prev_even) .and. present(prev_odd) )then
            call calc_fsc_range(self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, 1, self%ncls, fsc, &
                &ufrac_trec=ufrac_trec, prev_even=prev_even, prev_odd=prev_odd)
        else
            call calc_fsc_range(self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, 1, self%ncls, fsc)
        endif
    end subroutine calc_fsc

    subroutine calc_fsc_range( self, pfts_even, pfts_odd, ctf2_even, ctf2_odd, ref_first, ref_last, fsc, &
            &ufrac_trec, prev_even, prev_odd )
        class(polarft_calc),   intent(in)  :: self
        complex(dp),           intent(in)  :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp),           intent(in)  :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(in)  :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),              intent(in)  :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,               intent(in)  :: ref_first, ref_last
        real(dp),              intent(out) :: fsc(self%kfromto(1):self%interpklim)
        real(dp),    optional, intent(in)  :: ufrac_trec
        complex(dp), optional, intent(in)  :: prev_even(:,:,:)
        complex(dp), optional, intent(in)  :: prev_odd(:,:,:)
        complex(dp) :: even(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp) :: odd(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: vare(self%kfromto(1):self%interpklim)
        real(dp)    :: varo(self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
        integer     :: icls, k
        if( ref_first < 1 .or. ref_last > self%ncls .or. ref_first > ref_last )then
            THROW_HARD('invalid reference range in calc_fsc_range')
        endif
        if( self%p_ptr%l_trail_rec )then
            if( .not.present(ufrac_trec) .or. .not.present(prev_even) .or. .not.present(prev_odd) )then
                THROW_HARD('Trailing reconstruction requested but missing arguments in calc_fsc_range')
            endif
        endif
        fsc  = 0.d0; vare = 0.d0; varo = 0.d0
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(icls,even,odd,k,pft,ctf2) reduction(+:fsc,vare,varo)
        do icls = ref_first,ref_last
            pft  = pfts_even(:,:,icls)
            ctf2 = ctf2_even(:,:,icls)
            call self%shell_floor_norm(pft, ctf2, even)
            pft  = pfts_odd(:,:,icls)
            ctf2 = ctf2_odd(:,:,icls)
            call self%shell_floor_norm(pft, ctf2, odd)
            if( self%p_ptr%l_trail_rec )then
                even = ufrac_trec * even + (1.d0-ufrac_trec) * prev_even(:,:,icls)
                odd  = ufrac_trec * odd  + (1.d0-ufrac_trec) * prev_odd(:,:,icls)
            endif
            do k = self%kfromto(1),self%interpklim
                fsc(k)   = fsc(k)  + sum(real(even(:,k) * conjg(odd(:,k)), dp))
                vare(k)  = vare(k) + sum(real(even(:,k) * conjg(even(:,k)),dp))
                varo(k)  = varo(k) + sum(real(odd(:,k)  * conjg(odd(:,k)), dp))
            enddo
        enddo
        !$omp end parallel do
        vare = vare * varo
        where( vare > DTINY )
            fsc = fsc / sqrt(vare)
        elsewhere
            fsc = 0.d0
        end where
    end subroutine calc_fsc_range

    subroutine calc_fsc_from_refs( self, restored_even, restored_odd, ref_first, ref_last, fsc )
        class(polarft_calc), intent(in)  :: self
        complex(dp),         intent(in)  :: restored_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp),         intent(in)  :: restored_odd( self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,             intent(in)  :: ref_first, ref_last
        real(dp),            intent(out) :: fsc(self%kfromto(1):self%interpklim)
        real(dp) :: vare(self%kfromto(1):self%interpklim), varo(self%kfromto(1):self%interpklim)
        integer  :: icls, k
        if( ref_first < 1 .or. ref_last > self%ncls .or. ref_first > ref_last )then
            THROW_HARD('invalid reference range in calc_fsc_from_refs')
        endif
        fsc  = 0.d0
        vare = 0.d0
        varo = 0.d0
        ! Reference-domain FSC for already-restored polar refs read from POLAR_REFS_even/odd.
        !$omp parallel do default(shared) schedule(static) proc_bind(close) private(icls,k) reduction(+:fsc,vare,varo)
        do icls = ref_first,ref_last
            do k = self%kfromto(1),self%interpklim
                fsc(k)  = fsc(k)  + sum(real(restored_even(:,k,icls) * conjg(restored_odd(:,k,icls)), dp))
                vare(k) = vare(k) + sum(real(restored_even(:,k,icls) * conjg(restored_even(:,k,icls)), dp))
                varo(k) = varo(k) + sum(real(restored_odd( :,k,icls) * conjg(restored_odd( :,k,icls)), dp))
            enddo
        enddo
        !$omp end parallel do
        vare = vare * varo
        where( vare > DTINY )
            fsc = fsc / sqrt(vare)
        elsewhere
            fsc = 0.d0
        end where
    end subroutine calc_fsc_from_refs

    integer function ml_prior_start( self )
        class(polarft_calc), intent(in) :: self
        ml_prior_start = max(6, calc_fourier_index(self%p_ptr%hp, self%p_ptr%box_crop, self%p_ptr%smpd_crop))
    end function ml_prior_start

    ! Compute the ML regularization term to be added to the CTF2 in the denominator of the Wiener filter
    module subroutine add_invtausq2rho( self, ctf2_even, ctf2_odd, fsc )
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(inout) :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(inout) :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: fsc(self%kfromto(1):self%interpklim)
        call add_invtausq2rho_range(self, ctf2_even, ctf2_odd, 1, self%ncls, fsc)
    end subroutine add_invtausq2rho

    subroutine calc_invtausq2rho_range( self, ctf2_even, ctf2_odd, ref_first, ref_last, fsc, kstart, invtau2e, invtau2o )
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(in)    :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,             intent(in)    :: ref_first, ref_last
        real(dp),            intent(in)    :: fsc(self%kfromto(1):self%interpklim)
        integer,             intent(in)    :: kstart
        real(dp),            intent(out)   :: invtau2e(self%kfromto(1):self%interpklim)
        real(dp),            intent(out)   :: invtau2o(self%kfromto(1):self%interpklim)
        real(dp) :: sig2e(self%kfromto(1):self%interpklim), sig2o(self%kfromto(1):self%interpklim)
        real(dp) :: ssnr(self%kfromto(1):self%interpklim)
        real(dp) :: tau2(self%kfromto(1):self%interpklim)
        real(dp) :: cc, fudge, invtau2, shell_n
        integer  :: k
        if( ref_first < 1 .or. ref_last > self%ncls .or. ref_first > ref_last )then
            THROW_HARD('invalid reference range in calc_invtausq2rho_range')
        endif
        shell_n      = real(ref_last - ref_first + 1,dp) * real(self%pftsz,dp)
        invtau2e     = 0.d0
        invtau2o     = 0.d0
        !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k)
        do k = self%kfromto(1),self%interpklim
            sig2e(k) = sum(ctf2_even(:,k,ref_first:ref_last))
            sig2o(k) = sum(ctf2_odd(:, k,ref_first:ref_last))
        enddo
        !$omp end parallel do
        fudge = real(self%p_ptr%tau,dp)
        do k = self%kfromto(1),self%interpklim
            cc = max(0.001d0,min(0.999d0,fsc(k)))
            ssnr(k) = cc / (1.d0 - cc)
        enddo
        where( sig2e > DTINY )
            sig2e = shell_n / sig2e
        elsewhere
            sig2e = 0.d0
        end where
        tau2 = ssnr * sig2e
        do k = kstart,self%interpklim
            if( tau2(k) > DTINY )then
                invtau2 = 1.d0 / (fudge * tau2(k))
                invtau2e(k) = invtau2
            endif
        enddo
        where( sig2o > DTINY )
            sig2o = shell_n / sig2o
        elsewhere
            sig2o = 0.d0
        end where
        tau2 = ssnr * sig2o
        do k = kstart,self%interpklim
            if( tau2(k) > DTINY )then
                invtau2 = 1.d0 / (fudge * tau2(k))
                invtau2o(k) = invtau2
            endif
        enddo
    end subroutine calc_invtausq2rho_range

    subroutine add_invtausq2rho_range( self, ctf2_even, ctf2_odd, ref_first, ref_last, fsc )
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(inout) :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(inout) :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,             intent(in)    :: ref_first, ref_last
        real(dp),            intent(in)    :: fsc(self%kfromto(1):self%interpklim)
        real(dp) :: invtau2e(self%kfromto(1):self%interpklim), invtau2o(self%kfromto(1):self%interpklim)
        real(dp) :: invtau2
        integer  :: icls, k, kstart, p
        kstart = ml_prior_start(self)
        call calc_invtausq2rho_range(self, ctf2_even, ctf2_odd, ref_first, ref_last, fsc, kstart, invtau2e, invtau2o)
        !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
        do k = kstart,self%interpklim
            if( invtau2e(k) > 0.d0 )then
                invtau2 = invtau2e(k)
                ctf2_even(:,k,ref_first:ref_last) = ctf2_even(:,k,ref_first:ref_last) + invtau2
            else
                do icls = ref_first,ref_last
                    do p = 1,self%pftsz
                        invtau2 = unsampled_floor(ctf2_even(p,k,icls))
                        ctf2_even(p,k,icls) = ctf2_even(p,k,icls) + invtau2
                    enddo
                enddo
            endif
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
        do k = kstart,self%interpklim
            if( invtau2o(k) > 0.d0 )then
                invtau2 = invtau2o(k)
                ctf2_odd(:,k,ref_first:ref_last) = ctf2_odd(:,k,ref_first:ref_last) + invtau2
            else
                do icls = ref_first,ref_last
                    do p = 1,self%pftsz
                        invtau2 = unsampled_floor(ctf2_odd(p,k,icls))
                        ctf2_odd(p,k,icls) = ctf2_odd(p,k,icls) + invtau2
                    enddo
                enddo
            endif
        enddo
        !$omp end parallel do
    end subroutine add_invtausq2rho_range

    ! Performs the final normalization of references: CTF.I / (CTF2 + reg)
    module subroutine restore_references( self, reforis, pfts_even, pfts_odd, ctf2_even, ctf2_odd )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        complex(dp),         intent(in)    :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp),         intent(in)    :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp),            intent(in)    :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp) :: pfte(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp) :: pfto(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2e(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2o(self%pftsz,self%kfromto(1):self%interpklim)
        real        :: psi
        integer     :: icls, m, iloc, mloc, iglob, base, istate, nprojs, nprojs_nomirr, nwork
        logical     :: l_rotm
        nprojs        = self%p_ptr%nspace
        nprojs_nomirr = nprojs / 2
        if( nprojs * self%p_ptr%nstates /= self%ncls )then
            THROW_HARD('state-major polar reference count mismatch; restore_references')
        endif
        nwork = nprojs_nomirr * self%p_ptr%nstates
        ! Restoration
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(icls,m,iloc,mloc,iglob,base,istate,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
        do iglob = 1,nwork
            istate = (iglob - 1) / nprojs_nomirr + 1
            iloc   = mod(iglob - 1, nprojs_nomirr) + 1
            base   = (istate - 1) * nprojs
            icls   = base + iloc
            ! merged then e/o
            pfte  = pfts_even(:,:,icls)
            pfto  = pfts_odd(:,:,icls)
            ctf2e = ctf2_even(:,:,icls)
            ctf2o = ctf2_odd(:,:,icls)
            pft   = pfte  + pfto
            ctf2  = ctf2e + ctf2o
            call self%shell_floor_norm(pft,  ctf2,  self%pfts_merg(:,:,icls))
            call self%shell_floor_norm(pfte, ctf2e, self%pfts_even(:,:,icls))
            call self%shell_floor_norm(pfto, ctf2o, self%pfts_odd(:,:,icls))
            ! mirror the restored images (overwrites)
            mloc = reforis%get_int(iloc,'mirr')
            if( mloc < 1 .or. mloc > nprojs ) THROW_HARD('mirror index out of state-local range; restore_references')
            m    = base + mloc
            call self%mirror_pft(self%pfts_merg(:,:,icls), self%pfts_merg(:,:,m))
            call self%mirror_pft(self%pfts_even(:,:,icls), self%pfts_even(:,:,m))
            call self%mirror_pft(self%pfts_odd(:,:,icls),  self%pfts_odd(:,:,m))
            psi    = abs(reforis%get(mloc, 'psi'))
            l_rotm = (psi > 0.1) .and. (psi < 359.9)
            if( l_rotm )then
                self%pfts_merg(:,:,m) = conjg(self%pfts_merg(:,:,m))
                self%pfts_even(:,:,m) = conjg(self%pfts_even(:,:,m))
                self%pfts_odd(:,:,m)  = conjg(self%pfts_odd(:,:,m))
            endif
        enddo
        !$omp end parallel do
    end subroutine restore_references



    ! produces y-mirror of real (reciprocal) matrix 
    module pure subroutine mirror_ctf2( self, ctf2in, ctf2out )
        class(polarft_calc), intent(in)    :: self
        real(dp),            intent(in)    :: ctf2in(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp),            intent(inout) :: ctf2out(self%pftsz,self%kfromto(1):self%interpklim)
        integer :: i,j
        ctf2out(1,:) = ctf2in(1,:)
        do i = 2,self%pftsz/2
            j = self%pftsz-i+2
            ctf2out(i,:) = ctf2in(j,:)
            ctf2out(j,:) = ctf2in(i,:)
        enddo
        i = self%pftsz/2 + 1
        ctf2out(i,:) = ctf2in(i,:)
    end subroutine mirror_ctf2

    ! produces y-mirror of complex (reciprocal) matrix, for particles only
    module pure subroutine mirror_pft( self, pftin, pftout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pftin(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp),         intent(inout) :: pftout(self%pftsz,self%kfromto(1):self%interpklim)
        integer :: i,j
        pftout(1,:) = conjg(pftin(1,:))
        do i = 2,self%pftsz/2
            j = self%pftsz-i+2
            pftout(i,:) = pftin(j,:)
            pftout(j,:) = pftin(i,:)
        enddo
        i = self%pftsz/2 + 1
        pftout(i,:) = pftin(i,:)
    end subroutine mirror_pft

    ! Private utility to normalize PFTs by the CTF^2/sampling denominator.
    ! A per-shell floor is added to weak nonzero denominator entries before
    ! division. This is an empirical regularizer for sparse common-line support
    ! and CTF zeros, not only a guard against floating point zero division.
    module pure subroutine shell_floor_norm( self, Mnum, Mdenom, Mout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: Mnum(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp),            intent(inout) :: Mdenom(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp),         intent(inout) :: Mout(self%pftsz,self%kfromto(1):self%interpklim)
        logical  :: msk(self%pftsz)
        real(dp) :: avg, t
        integer  :: k, n
        do k = self%kfromto(1),self%interpklim
            msk = Mdenom(:,k) > DSMALL
            n   = count(msk)
            if( n == 0 ) cycle
            avg = sum(Mdenom(:,k),mask=msk) / real(n,dp)
            t   = avg/50.d0
            where((Mdenom(:,k) < t).and.msk) Mdenom(:,k) = Mdenom(:,k) + t
        enddo
        where( Mdenom > DSMALL)
            Mout = Mnum / Mdenom
        elsewhere
            Mout = DCMPLX_ZERO
        end where
    end subroutine shell_floor_norm

end submodule simple_polarft_ops_restore

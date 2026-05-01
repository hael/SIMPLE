!@descr: submodule for controlling various state-related things in the polarops module
submodule (simple_polarft_calc) simple_polarft_ops_state
use simple_memoize_ft_maps, only: memoize_ft_maps
use simple_refine3D_fnames, only: refine3D_obsfield_part_fname
use simple_fgrid_obsfield,  only: fgrid_obsfield_eo
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
                &self%pfts_merg(self%pftsz, self%kfromto(1):self%interpklim, self%ncls))
        if( trim(self%p_ptr%polar) /= 'obsfield' )then
            allocate(self%ctf2_even(self%pftsz, self%kfromto(1):self%interpklim, self%ncls),&
                    &self%ctf2_odd( self%pftsz, self%kfromto(1):self%interpklim, self%ncls))
        endif
        call self%polar_cavger_zero_pft_refs
        self%pfts_merg = DCMPLX_ZERO
    end subroutine polar_cavger_new

    module subroutine polar_cavger_zero_pft_refs( self )
        class(polarft_calc), intent(inout) :: self
        !$omp parallel workshare
        self%pfts_even = DCMPLX_ZERO
        self%pfts_odd  = DCMPLX_ZERO
        !$omp end parallel workshare
        if( allocated(self%ctf2_even) )then
            !$omp parallel workshare
            self%ctf2_even = 0.d0
            self%ctf2_odd  = 0.d0
            !$omp end parallel workshare
        endif
    end subroutine polar_cavger_zero_pft_refs

    subroutine obsfield_lims_from_params( self, lims )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(out) :: lims(3,2)
        integer :: box
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
    end subroutine obsfield_lims_from_params

    subroutine ensure_obsfields_allocated( self )
        class(polarft_calc), intent(inout) :: self
        integer :: lims(3,2), istate
        if( allocated(self%obsfields) ) return
        call obsfield_lims_from_params(self, lims)
        allocate(self%obsfields(self%p_ptr%nstates))
        do istate = 1, self%p_ptr%nstates
            call self%obsfields(istate)%new(lims, self%interpklim)
        enddo
    end subroutine ensure_obsfields_allocated

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
        integer :: i, icls, iptcl, eo, state, proj
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
        ! 3D polar references are state-major, while get_proj returns the
        ! state-local projection index in 1:nspace.
        !$omp parallel do schedule(guided) proc_bind(close) default(shared)&
        !$omp private(iptcl,eo,icls,state,proj)&
        !$omp reduction(+:eo_pops)
        do iptcl = 1,ptcl_field%get_noris()
            state = ptcl_field%get_state(iptcl)
            if( state == 0  ) cycle
            if( ptcl_field%get(iptcl,'w') < SMALL ) cycle
            eo = ptcl_field%get_eo(iptcl)+1
            if( l_3D )then
                proj = ptcl_field%get_proj(iptcl)
                if( proj < 1 .or. proj > self%p_ptr%nspace ) cycle
                icls = (state - 1) * self%p_ptr%nspace + proj
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
        integer     :: eopops(2,self%ncls), i, icls, iptcl, irot, k, state, proj
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
        !$omp private(i,iptcl,w,l_even,icls,irot,incr_shift,rptcl,rctf,rctf2,k,sigma2,state,proj)
        do i = 1,nptcls
            ! particles parameters
            iptcl = pinds(i)
            state = spproj_field%get_state(iptcl)
            if( state == 0  ) cycle
            w = real(spproj_field%get(iptcl,'w'),dp)
            if( w < DSMALL ) cycle
            l_even = spproj_field%get_eo(iptcl)==0
            if( l_3D )then
                proj = spproj_field%get_proj(iptcl)
                if( proj < 1 .or. proj > self%p_ptr%nspace ) cycle
                icls = (state - 1) * self%p_ptr%nspace + proj
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

    ! Observation-field restoration experiment. Workers only accumulate dense
    ! state-local observation fields here; assembly later extracts polar sections
    ! from the reduced fields, matching the polar=no insert-then-assemble shape.
    module subroutine polar_cavger_insert_ptcls_obsfield( self, eulspace, ptcl_field, symop, nptcls, pinds, fpls, &
        &reforis_in, nspace_out, incr_shifts )
        class(polarft_calc),           intent(inout) :: self
        class(oris), target,           intent(inout) :: eulspace
        class(oris), pointer,          intent(inout) :: ptcl_field
        class(sym),                    intent(inout) :: symop
        integer,                       intent(in)    :: nptcls, pinds(nptcls)
        class(fplane_type), target,    intent(inout) :: fpls(nptcls)
        ! Retained for source compatibility; assembly now owns the output
        ! reference space for obsfield.
        class(oris), target, optional, intent(inout) :: reforis_in
        integer,             optional, intent(in)    :: nspace_out
        real,                optional, intent(in)    :: incr_shifts(2,nptcls)
        type(ori) :: o
        real :: crop_factor, pw, shift_crop(2)
        integer :: i, iptcl, eo, pstate
        logical :: l_shift
        if( nptcls < 1 ) return
        call ensure_obsfields_allocated(self)
        l_shift = present(incr_shifts)
        if( l_shift ) crop_factor = real(self%p_ptr%box_crop) / real(self%p_ptr%box)
        do i = 1, nptcls
            iptcl  = pinds(i)
            pstate = ptcl_field%get_state(iptcl)
            if( pstate < 1 .or. pstate > self%p_ptr%nstates ) cycle
            pw = 1.0
            if( ptcl_field%isthere(iptcl,'w') ) pw = ptcl_field%get(iptcl,'w')
            if( pw < TINY ) cycle
            eo = ptcl_field%get_eo(iptcl)
            call ptcl_field%get_ori(iptcl, o)
            if( l_shift )then
                shift_crop = incr_shifts(:,i) * crop_factor
                call self%obsfields(pstate)%insert_plane(symop, o, fpls(i), eo, pw, shift_crop=shift_crop)
            else
                call self%obsfields(pstate)%insert_plane(symop, o, fpls(i), eo, pw)
            endif
        enddo
        call o%kill
    end subroutine polar_cavger_insert_ptcls_obsfield

    module subroutine polar_cavger_write_obsfield_parts( self )
        class(polarft_calc), intent(inout) :: self
        type(string) :: fname
        integer :: istate
        call ensure_obsfields_allocated(self)
        do istate = 1, self%p_ptr%nstates
            fname = refine3D_obsfield_part_fname(istate, self%p_ptr%part, self%p_ptr%numlen)
            call self%obsfields(istate)%write(fname)
        enddo
        call fname%kill
    end subroutine polar_cavger_write_obsfield_parts

    module subroutine polar_cavger_assemble_obsfields_from_parts( self, reforis, bench )
        class(polarft_calc), intent(inout) :: self
        class(oris), target, intent(inout) :: reforis
        real(timer_int_kind), optional, intent(out) :: bench(:)
        type(fgrid_obsfield_eo) :: obs_part
        type(string) :: fname
        integer(timer_int_kind) :: t_read, t_append
        integer :: istate, ipart
        real(timer_int_kind) :: rt_read, rt_append
        rt_read   = 0.
        rt_append = 0.
        if( present(bench) ) bench = 0.
        call ensure_obsfields_allocated(self)
        do istate = 1, self%p_ptr%nstates
            call self%obsfields(istate)%reset
            do ipart = 1, self%p_ptr%nparts
                fname = refine3D_obsfield_part_fname(istate, ipart, self%p_ptr%numlen)
                if( L_BENCH_GLOB ) t_read = tic()
                call obs_part%read(fname)
                if( L_BENCH_GLOB ) rt_read = rt_read + toc(t_read)
                if( L_BENCH_GLOB ) t_append = tic()
                call self%obsfields(istate)%append_field(obs_part)
                if( L_BENCH_GLOB ) rt_append = rt_append + toc(t_append)
                call obs_part%kill
            enddo
        enddo
        if( present(bench) )then
            if( size(bench) >= 1 ) bench(1) = rt_read
            if( size(bench) >= 2 ) bench(2) = rt_append
        endif
        call obs_part%kill
        call fname%kill
    end subroutine polar_cavger_assemble_obsfields_from_parts

    module subroutine polar_cavger_kill( self )
        class(polarft_calc), intent(inout) :: self
        integer :: istate
        if( allocated(self%pfts_even)    ) deallocate(self%pfts_even)
        if( allocated(self%pfts_odd)     ) deallocate(self%pfts_odd)
        if( allocated(self%pfts_merg)    ) deallocate(self%pfts_merg)
        if( allocated(self%ctf2_even)    ) deallocate(self%ctf2_even)
        if( allocated(self%ctf2_odd)     ) deallocate(self%ctf2_odd)
        if( allocated(self%obsfields) )then
            do istate = 1, size(self%obsfields)
                call self%obsfields(istate)%kill
            enddo
            deallocate(self%obsfields)
        endif
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

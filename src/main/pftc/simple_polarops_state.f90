!@descr: submodule for controlling various state-related things
submodule (simple_polarops) simple_polarops_state
implicit none
#include "simple_local_flags.inc"
contains

    !> Module initialization
    module subroutine polar_cavger_new( pftc, comlin, nrefs )
        class(polarft_calc),     intent(in) :: pftc
        logical,                 intent(in) :: comlin
        integer,       optional, intent(in) :: nrefs
        call polar_cavger_kill
        nrots    = pftc%get_nrots()
        pftsz    = pftc%get_pftsz()
        kfromto  = pftc%get_kfromto()
        l_comlin = comlin
        smpd     = params_glob%smpd
        if( present(nrefs) )then
            ncls = nrefs
        else
            ncls = pftc%get_nrefs()
        endif
        allocate(prev_eo_pops(2,ncls), eo_pops(2,ncls), source=0)
        allocate(pfts_even(pftsz,kfromto(1):kfromto(2),ncls),pfts_odd(pftsz,kfromto(1):kfromto(2),ncls),&
                &ctf2_even(pftsz,kfromto(1):kfromto(2),ncls),ctf2_odd(pftsz,kfromto(1):kfromto(2),ncls),&
                &pfts_merg(pftsz,kfromto(1):kfromto(2),ncls))
        call polar_cavger_zero_pft_refs
        pfts_merg = DCMPLX_ZERO
    end subroutine polar_cavger_new

    module subroutine polar_cavger_zero_pft_refs
        !$omp parallel workshare
        pfts_even = DCMPLX_ZERO
        pfts_odd  = DCMPLX_ZERO
        ctf2_even = 0.d0
        ctf2_odd  = 0.d0
        !$omp end parallel workshare
    end subroutine polar_cavger_zero_pft_refs

    module subroutine polar_cavger_set_ref_pftc( icls, which, pftc )
        integer,             intent(in)    :: icls
        character(len=*),    intent(in)    :: which
        class(polarft_calc), intent(inout) :: pftc
        select case(trim(which))
        case('merged')
            call pftc%set_ref_pft(icls, cmplx(pfts_merg(:,:,icls),kind=sp), .true.)
        case('even')
            call pftc%set_ref_pft(icls, cmplx(pfts_even(:,:,icls),kind=sp), .true.)
        case('odd')
            call pftc%set_ref_pft(icls, cmplx(pfts_odd(:,:,icls),kind=sp), .false.)
        end select
    end subroutine polar_cavger_set_ref_pftc

    module subroutine polar_cavger_calc_pops( spproj )
        class(sp_project), target, intent(in) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer :: i, icls, iptcl, eo
        logical :: l_3D
        l_3D = .false.
        select case(trim(params_glob%oritype))
        case('ptcl2D')
            ptcl_field => spproj%os_ptcl2D
            cls_field  => spproj%os_cls2D
        case('ptcl3D')
            ptcl_field => spproj%os_ptcl3D
            cls_field  => spproj%os_cls3D
            l_3D       = .true.
        case DEFAULT
            THROW_HARD('Unsupported ORITYPE: '//trim(params_glob%oritype))
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
            eo_pops(eo,icls) = eo_pops(eo,icls) + 1
        enddo
        !$omp end parallel do
        prev_eo_pops = 0
        if( cls_field%get_noris() == ncls )then
            do i = 1,ncls
                if( l_3D )then
                    icls = ptcl_field%get_proj(i)
                else
                    icls = ptcl_field%get_class(i)
                endif
                if( .not.cls_field%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(1,icls) = cls_field%get_int(i,'prev_pop_even')
                prev_eo_pops(2,icls) = cls_field%get_int(i,'prev_pop_odd')
            enddo
        endif
        eo_pops = eo_pops + prev_eo_pops
    end subroutine polar_cavger_calc_pops

    !>  \brief  Updates Fourier components and normalization matrices with new particles
    module subroutine polar_cavger_update_sums( nptcls, pinds, spproj, pftc, incr_shifts, is3D )
        use simple_euclid_sigma2, only: eucl_sigma2_glob
        integer,                         intent(in)    :: nptcls
        integer,                         intent(in)    :: pinds(nptcls)
        class(sp_project),               intent(inout) :: spproj
        class(polarft_calc), target,     intent(inout) :: pftc
        real,                  optional, intent(in)    :: incr_shifts(2,nptcls)
        logical,               optional, intent(in)    :: is3d
        class(oris), pointer :: spproj_field
        complex(sp), pointer :: pptcls(:,:,:), rptcl(:,:)
        real(sp),    pointer :: pctfmats(:,:,:), rctf(:,:)
        real(dp),    pointer :: rctf2(:,:)
        real(dp) :: w
        real     :: sigma2(kfromto(1):kfromto(2)), incr_shift(2)
        integer  :: eopops(2,ncls), i, icls, iptcl, irot, k
        logical  :: l_ctf, l_even, l_3D, l_shift
        l_3D = .false.
        if( present(is3D) ) l_3D = is3D
        l_shift = .false.
        if( present(incr_shifts) ) l_shift = .true.
        ! retrieve particle info & pointers
        call spproj%ptr2oritype(params_glob%oritype, spproj_field)
        call pftc%get_ptcls_ptr(pptcls)
        l_ctf = pftc%is_with_ctf()
        if( l_ctf ) call pftc%get_ctfmats_ptr(pctfmats)
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
            irot = pftc%get_roind_fast(spproj_field%e3get(iptcl))
            if( l_shift )then
                incr_shift = incr_shifts(:,i)
                ! weighted restoration
                if( any(abs(incr_shift) > 1.e-6) ) call pftc%shift_ptcl(iptcl, -incr_shift)
            endif
            call pftc%get_work_pft_ptr(rptcl)
            ! Particle rotation
            call pftc%rotate_pft(pptcls(:,:,i), irot, rptcl)
            ! Particle weight
            rptcl = real(w) * rptcl
            ! Particle ML regularization
            if( params_glob%l_ml_reg )then
                sigma2 = eucl_sigma2_glob%sigma2_noise(kfromto(1):kfromto(2),iptcl)
                do k = kfromto(1),kfromto(2)
                    rptcl(:,k) = rptcl(:,k) / sigma2(k)
                enddo
            endif
            ! Array updates
            if( l_ctf )then
                call pftc%get_work_rpft_ptr(rctf)
                call pftc%get_work_rpft8_ptr(rctf2)
                ! weighted CTF2
                call pftc%rotate_pft(pctfmats(:,:,i), irot, rctf)
                rctf2 = w * real(rctf,kind=dp)**2
                rptcl = rptcl * rctf    ! PhFlip(X).|CTF|
                ! CTF2 ML regularization
                if( params_glob%l_ml_reg )then
                    do k = kfromto(1),kfromto(2)
                        rctf2(:,k) = rctf2(:,k) / real(sigma2(k),dp)
                    enddo
                endif
                if( l_even )then
                    !$omp critical
                    pfts_even(:,:,icls) = pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                    ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + rctf2
                    !$omp end critical
                else
                    !$omp critical
                    pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                    ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + rctf2
                    !$omp end critical
                endif
            else
                if( params_glob%l_ml_reg )then
                    ! CTF2=1 & ML regularization
                    call pftc%get_work_rpft8_ptr(rctf2)
                    do k = kfromto(1),kfromto(2)
                        rctf2(:,k) = w / real(sigma2(k),dp)
                    enddo
                    if( l_even )then
                        !$omp critical
                        pfts_even(:,:,icls) = pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                        ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + rctf2
                        !$omp end critical
                    else
                        !$omp critical
                        pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                        ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + rctf2
                        !$omp end critical
                    endif
                else
                    if( l_even )then
                        !$omp critical
                        pfts_even(:,:,icls) = pfts_even(:,:,icls) + cmplx(rptcl,kind=dp)
                        ctf2_even(:,:,icls) = ctf2_even(:,:,icls) + w
                        !$omp end critical
                    else
                        !$omp critical
                        pfts_odd(:,:,icls)  = pfts_odd(:,:,icls)  + cmplx(rptcl,kind=dp)
                        ctf2_odd(:,:,icls)  = ctf2_odd(:,:,icls)  + w
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
        eo_pops = eo_pops + eopops
        ! cleanup
        nullify(spproj_field,rptcl,rctf,pptcls,pctfmats)
    end subroutine polar_cavger_update_sums

    module subroutine polar_cavger_kill
        if( allocated(pfts_even) )then
            deallocate(pfts_even,pfts_odd,ctf2_even,ctf2_odd,pfts_merg,eo_pops,prev_eo_pops)
        endif
        smpd       = 0.
        ncls       = 0
        nrots      = 0
        kfromto(2) = 0
        pftsz      = 0
        l_comlin   = .false.
    end subroutine polar_cavger_kill

    ! Estimates the center of the volume based on the distribution of
    ! the individual particles in-plane offsets and map the shifts to both
    ! the particles and the references stored in the pftc
    module subroutine center_3Dpolar_refs( pftc, algndoc, algnrefs )
        class(polarft_calc), intent(inout) :: pftc
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
        do iref = 1,ncls
            ! Projection direction rotation matrix
            R = euler2m([algnrefs%e1get(iref), algnrefs%e2get(iref), 0.0])
            ! 3D Shift rotated with respect to projection direction
            offset2D = matmul(R, offset3D)
            ! Apply offset to e/o references
            call pftc%shift_ref(iref, offset2D(1:2))
        enddo
        !$omp end parallel do
    end subroutine center_3Dpolar_refs

end submodule simple_polarops_state

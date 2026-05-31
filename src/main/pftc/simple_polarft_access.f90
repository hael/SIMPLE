!@descr: polarft class accessors submodule
submodule (simple_polarft_calc) simple_polarft_access
#include "simple_local_flags.inc"
implicit none

contains

    ! ===== GETTERS + POINTER ACCESSORS =====

    module pure function get_nrots(self) result(nrots)
        class(polarft_calc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots

    module pure function get_pdim_srch(self) result(pdim)
        class(polarft_calc), intent(in) :: self
        integer :: pdim(3)
        pdim = [self%pftsz,self%kfromto(1),self%kfromto(2)]
    end function get_pdim_srch

    module pure function get_pdim_interp(self) result(dims)
        class(polarft_calc), intent(in) :: self
        integer :: dims(3)
        dims = [self%pftsz,self%kfromto(1),self%interpklim]
    end function get_pdim_interp

    module pure function get_kfromto(self) result(kfromto)
        class(polarft_calc), intent(in) :: self
        integer :: kfromto(2)
        kfromto = [self%kfromto(1),self%kfromto(2)]
    end function get_kfromto

    module pure integer function get_pftsz(self)
        class(polarft_calc), intent(in) :: self
        get_pftsz = self%pftsz
    end function get_pftsz

    module function get_rot(self, roind) result(rot)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: roind
        real(sp) :: rot
        if( roind < 1 .or. roind > self%nrots )then
            write(logfhandle,*) 'roind: ', roind
            write(logfhandle,*) 'nrots: ', self%nrots
            THROW_HARD('roind is out of range; get_rot')
        endif
        rot = self%angtab(roind)
    end function get_rot

    module function get_roind(self, rot) result(ind)
        class(polarft_calc), intent(in) :: self
        real(sp),            intent(in) :: rot
        real(sp) :: dists(self%nrots)
        integer :: ind
        dists = abs(self%angtab-rot)
        where(dists>180.)dists = 360.-dists
        ind = minloc(dists, dim=1)
    end function get_roind

    module pure integer function get_roind_fast(self, psi)
        class(polarft_calc), intent(in) :: self
        real,                intent(in) :: psi
        get_roind_fast = modulo(nint(psi / self%dang), self%nrots) + 1
    end function get_roind_fast

    module pure real function get_dang(self)
        class(polarft_calc), intent(in) :: self
        get_dang = self%dang
    end function get_dang

    module pure function get_coord(self, rot, k) result(xy)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: rot, k
        real(sp) :: xy(2)
        xy(1) = self%polar(1,k,rot)
        xy(2) = self%polar(2,k,rot)
    end function get_coord

    module subroutine get_polar_coords(self, dim, coords)
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: dim
        real(sp),            intent(out) :: coords(self%kfromto(1):self%kfromto(2),self%pftsz)
        if( dim == 1 )then
            coords = self%polar(1,self%kfromto(1):self%kfromto(2),:self%pftsz)
        else if( dim == 2 )then
            coords = self%polar(2,self%kfromto(1):self%kfromto(2),:self%pftsz)
        else
            THROW_HARD('invalid coordinate dimension in get_polar_coords')
        endif
    end subroutine get_polar_coords

    module subroutine get_ref_pft(self, iref, iseven, pft)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iref
        logical,             intent(in) :: iseven
        complex(sp), intent(inout) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        if( iseven )then
            pft = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
    end subroutine get_ref_pft

    module subroutine get_ptcl_pft(self, iptcl, pft)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iptcl
        complex(sp),         intent(inout) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        if( iptcl < self%pfromto(1) .or. iptcl > self%pfromto(2) )then
            write(logfhandle,*) 'iptcl: ', iptcl
            write(logfhandle,*) 'pfromto: ', self%pfromto
            THROW_HARD('iptcl is out of range; get_ptcl_pft')
        endif
        pft = self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),self%pinds(iptcl))
    end subroutine get_ptcl_pft

    module subroutine get_ptcl_line(self, iptcl, irot, line)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iptcl, irot
        complex(sp),         intent(inout) :: line(self%kfromto(1):self%kfromto(2))
        if( iptcl < self%pfromto(1) .or. iptcl > self%pfromto(2) )then
            write(logfhandle,*) 'iptcl: ', iptcl
            write(logfhandle,*) 'pfromto: ', self%pfromto
            THROW_HARD('iptcl is out of range; get_ptcl_line')
        endif
        if( irot < 1 .or. irot > self%pftsz )then
            write(logfhandle,*) 'irot: ', irot
            write(logfhandle,*) 'pftsz: ', self%pftsz
            THROW_HARD('irot is out of range; get_ptcl_line')
        endif
        line = self%pfts_ptcls(irot,self%kfromto(1):self%kfromto(2),self%pinds(iptcl))
    end subroutine get_ptcl_line

    module pure integer function get_nrefs(self)
        class(polarft_calc), intent(in) :: self
        get_nrefs = self%nrefs
    end function get_nrefs

    module logical function exists(self)
        class(polarft_calc), intent(in) :: self
        exists = self%existence
    end function exists

    module logical function ptcl_iseven(self, iptcl)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iptcl
        ptcl_iseven = self%iseven(self%pinds(iptcl))
    end function ptcl_iseven

    module integer function get_nptcls(self)
        class(polarft_calc), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    module subroutine get_pinds(self, pinds)
        class(polarft_calc),  intent(in)  :: self
        integer, allocatable, intent(out) :: pinds(:)
        pinds = self%pinds
    end subroutine get_pinds

    module pure logical function is_with_ctf( self )
        class(polarft_calc), intent(in) :: self
        is_with_ctf = self%with_ctf
    end function is_with_ctf

    module pure integer function get_nctf_models_seen(self)
        class(polarft_calc), intent(in) :: self
        get_nctf_models_seen = self%nctf_models_seen
    end function get_nctf_models_seen

    module subroutine reset_ctf_model_audit(self)
        class(polarft_calc), intent(inout) :: self
        self%ctf_model_audit = .true.
        self%nctf_models_seen = 0
        if( allocated(self%ctf_models_seen) ) deallocate(self%ctf_models_seen)
        allocate(self%ctf_models_seen(max(1,self%nptcls)))
        if( allocated(self%ctfparams_ptcls) ) deallocate(self%ctfparams_ptcls)
        if( allocated(self%ctfparams_ptcls_set) ) deallocate(self%ctfparams_ptcls_set)
        allocate(self%ctfparams_ptcls(1:self%nptcls), self%ctfparams_ptcls_set(1:self%nptcls))
        self%ctfparams_ptcls_set = .false.
        if( allocated(self%ctfparams_scored) ) self%ctfparams_scored = .false.
    end subroutine reset_ctf_model_audit

    module subroutine disable_ctf_model_audit(self)
        class(polarft_calc), intent(inout) :: self
        self%ctf_model_audit = .false.
        self%nctf_models_seen = 0
        if( allocated(self%ctf_models_seen) ) deallocate(self%ctf_models_seen)
        if( allocated(self%ctfparams_ptcls) ) deallocate(self%ctfparams_ptcls)
        if( allocated(self%ctfparams_ptcls_set) ) deallocate(self%ctfparams_ptcls_set)
    end subroutine disable_ctf_model_audit

    module integer function get_nctf_models_scored(self)
        class(polarft_calc), intent(in) :: self
        get_nctf_models_scored = self%nctf_models_scored
    end function get_nctf_models_scored

    module subroutine reset_ctf_scoring_audit(self)
        class(polarft_calc), intent(inout) :: self
        self%ctf_scoring_audit = .true.
        self%nctf_models_scored = 0
        if( allocated(self%ctf_models_scored) ) deallocate(self%ctf_models_scored)
        allocate(self%ctf_models_scored(max(1,self%nptcls)))
        if( .not. allocated(self%ctfparams_scored) ) allocate(self%ctfparams_scored(1:self%nptcls))
        self%ctfparams_scored = .false.
    end subroutine reset_ctf_scoring_audit

    module subroutine disable_ctf_scoring_audit(self)
        class(polarft_calc), intent(inout) :: self
        self%ctf_scoring_audit = .false.
        self%nctf_models_scored = 0
        if( allocated(self%ctf_models_scored) ) deallocate(self%ctf_models_scored)
        if( allocated(self%ctfparams_scored) ) deallocate(self%ctfparams_scored)
    end subroutine disable_ctf_scoring_audit

    module function allocate_pft( self ) result( pft )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), allocatable :: pft(:,:)
        allocate(pft(self%pftsz,self%kfromto(1):self%interpklim), source=CMPLX_ZERO)
    end function allocate_pft

    module function allocate_ptcl_pft( self ) result( pft )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), allocatable :: pft(:,:)
        allocate(pft(self%pftsz,self%kfromto(1):self%interpklim), source=CMPLX_ZERO)
    end function allocate_ptcl_pft

    module pure subroutine get_precalc_objfun_vals(self, ind, ithr, vals)
        class(polarft_calc),  intent(in)  :: self
        integer,              intent(in)  :: ind, ithr
        real,                 intent(out) :: vals(self%nrots)
        vals(:) = self%crmat_many(ithr)%r(1:self%nrots, ind)
    end subroutine get_precalc_objfun_vals

    logical function same_ctf_model(lhs, rhs)
        type(ctfparams), intent(in) :: lhs, rhs
        real, parameter :: CTFTOL = 1.0e-5
        same_ctf_model = (lhs%ctfflag == rhs%ctfflag) .and. &
            (lhs%l_phaseplate .eqv. rhs%l_phaseplate) .and. &
            (abs(lhs%kv    - rhs%kv)    <= CTFTOL) .and. &
            (abs(lhs%cs    - rhs%cs)    <= CTFTOL) .and. &
            (abs(lhs%fraca - rhs%fraca) <= CTFTOL)
    end function same_ctf_model

end submodule simple_polarft_access

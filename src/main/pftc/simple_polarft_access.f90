submodule (simple_polarft_calc) simple_polarft_access
include 'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none

contains

    ! ===== GETTERS + POINTER ACCESSORS =====

    module pure function get_nrots(self) result(nrots)
        class(polarft_calc), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots

    module pure function get_pdim(self) result(pdim)
        class(polarft_calc), intent(in) :: self
        integer :: pdim(3)
        pdim = [self%pftsz,self%kfromto(1),self%kfromto(2)]
    end function get_pdim

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
        get_roind_fast = nint(psi / self%dang) + 1
        if( get_roind_fast <=         0 ) get_roind_fast = get_roind_fast + self%nrots
        if( get_roind_fast > self%nrots ) get_roind_fast = get_roind_fast - self%nrots
    end function get_roind_fast

    module pure real function get_dang(self)
        class(polarft_calc), intent(in) :: self
        get_dang = self%dang
    end function get_dang

    module pure function get_coord(self, rot, k) result(xy)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: rot, k
        real(sp) :: xy(2)
        xy(1) = self%polar(rot,k)
        xy(2) = self%polar(self%nrots+rot,k)
    end function get_coord

    module subroutine get_ref_pft(self, iref, iseven, pft)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iref
        logical,             intent(in) :: iseven
        complex(sp), intent(inout) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        if( iseven )then
            pft = self%pfts_refs_even(:,:,iref)
        else
            pft = self%pfts_refs_odd(:,:,iref)
        endif
    end subroutine get_ref_pft

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

    module integer function get_npix(self)
        class(polarft_calc), intent(in) :: self
        get_npix = sum(nint(self%npix_per_shell(self%kfromto(1):self%kfromto(2))))
    end function get_npix

    module pure logical function is_with_ctf( self )
        class(polarft_calc), intent(in) :: self
        is_with_ctf = self%with_ctf
    end function is_with_ctf

    module subroutine get_work_pft_ptr( self, ptr )
        class(polarft_calc),  intent(in)  :: self
        complex(sp), pointer, intent(out) :: ptr(:,:)
        integer :: ithr
        ithr = omp_get_thread_num()+1
        ptr => self%heap_vars(ithr)%pft_ref_tmp
    end subroutine get_work_pft_ptr

    module subroutine get_work_rpft_ptr( self, ptr )
        class(polarft_calc), intent(in)  :: self
        real(sp), pointer,   intent(out) :: ptr(:,:)
        integer :: ithr
        ithr = omp_get_thread_num()+1
        ptr => self%heap_vars(ithr)%pft_r
    end subroutine get_work_rpft_ptr

    module subroutine get_work_rpft8_ptr( self, ptr )
        class(polarft_calc), intent(in)  :: self
        real(dp), pointer,   intent(out) :: ptr(:,:)
        integer :: ithr
        ithr = omp_get_thread_num()+1
        ptr => self%heap_vars(ithr)%pft_r1_8
    end subroutine get_work_rpft8_ptr

    module subroutine get_ptcls_ptr( self, ptr )
        class(polarft_calc), target,  intent(in)  :: self
        complex(sp),         pointer, intent(out) :: ptr(:,:,:)
        ptr => self%pfts_ptcls
    end subroutine get_ptcls_ptr

    module subroutine get_ctfmats_ptr( self, ptr )
        class(polarft_calc), target,  intent(in)  :: self
        real(sp),            pointer, intent(out) :: ptr(:,:,:)
        ptr => self%ctfmats
    end subroutine get_ctfmats_ptr

    module subroutine get_refs_ptr( self, ptre, ptro )
        class(polarft_calc), target,  intent(in)  :: self
        complex(sp),         pointer, intent(out) :: ptre(:,:,:), ptro(:,:,:)
        ptre => self%pfts_refs_even
        ptro => self%pfts_refs_odd
    end subroutine get_refs_ptr

end submodule simple_polarft_access
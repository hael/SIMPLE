submodule (simple_polarft_calc) simple_polarft_core
include 'simple_lib.f08'
#include "simple_local_flags.inc"
use simple_parameters, only: params_glob
implicit none

contains

    ! ===== CORE (new, kill, setters, getters, pointer helpers) =====

    module subroutine new(self, nrefs, pfromto, kfromto, eoarr)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: nrefs
        integer,                     intent(in)    :: pfromto(2), kfromto(2)
        integer, optional,           intent(in)    :: eoarr(pfromto(1):pfromto(2))
        real(sp), allocatable :: polar_here(:)
        real(dp)              :: A(2)
        integer               :: irot, k, ithr, i
        logical               :: even_dims, test(2)
        call self%kill
        ! set particle index range
        self%pfromto = pfromto
        ! set band-pass Fourier index limits
        self%kfromto = kfromto
        self%nk      = self%kfromto(2) - self%kfromto(1) + 1
        ! error check
        if( self%pfromto(2) - self%pfromto(1) + 1 < 1 )then
            write(logfhandle,*) 'pfromto: ', self%pfromto(1), self%pfromto(2)
            THROW_HARD ('nptcls (# of particles) must be > 0; new')
        endif
        if( nrefs < 1 )then
            write(logfhandle,*) 'nrefs: ', nrefs
            THROW_HARD ('nrefs (# of reference sections) must be > 0; new')
        endif
        self%ldim = [params_glob%box,params_glob%box,1] !< logical dimensions of original cartesian image
        test      = .false.
        test(1)   = is_even(self%ldim(1))
        test(2)   = is_even(self%ldim(2))
        even_dims = all(test)
        if( .not. even_dims )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            THROW_HARD ('only even logical dims supported; new')
        endif
        ! set constants
        self%nptcls = self%pfromto(2) - self%pfromto(1) + 1   !< the total number of particles in partition
        self%nrefs  = nrefs                                   !< the number of references (logically indexded [1,nrefs])
        self%pftsz  = magic_pftsz(nint(params_glob%msk_crop)) !< size of reference (number of vectors used for matching,determined by radius of molecule)
        self%nrots  = 2 * self%pftsz                          !< number of in-plane rotations for one pft  (pftsz*2)
        ! generate polar coordinates
        allocate( self%polar(2*self%nrots,self%kfromto(1):self%kfromto(2)),&
                    &self%angtab(self%nrots), self%iseven(1:self%nptcls), polar_here(2*self%nrots))
        self%dang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = real(irot-1)*self%dang
            ! cycling over non-redundant logical dimensions
            do k=self%kfromto(1),self%kfromto(2)
                self%polar(irot,k)            =  sin(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,k) = -cos(self%angtab(irot))*real(k) ! y-coordinate
            end do
            ! for k = 1
            polar_here(irot)            =  sin(real(self%angtab(irot)))
            polar_here(irot+self%nrots) = -cos(real(self%angtab(irot)))
            ! angle (in degrees) from now
            self%angtab(irot) = rad2deg(self%angtab(irot))
        end do
        self%dang = rad2deg(self%dang)
        ! index translation table
        allocate( self%pinds(self%pfromto(1):self%pfromto(2)), source=0 )
        self%pinds = (/(i,i=1,self%nptcls)/)
        ! eo assignment
        if( present(eoarr) )then
            if( all(eoarr == - 1) )then
                self%iseven = .true.
            else
                do i=self%pfromto(1),self%pfromto(2)
                    if( self%pinds(i) > 0 )then
                        if( eoarr(i) == 0 )then
                            self%iseven(self%pinds(i)) = .true.
                        else
                            self%iseven(self%pinds(i)) = .false.
                        endif
                    endif
                end do
            endif
        else
            self%iseven = .true.
        endif
        ! generate the argument transfer constants for shifting reference polarfts
        allocate( self%argtransf(self%nrots,self%kfromto(1):self%kfromto(2)), self%argtransf_shellone(self%nrots) )
        A = DPI / real(self%ldim(1:2)/2,dp) ! argument transfer matrix normalization constant
        ! shell = 1
        self%argtransf_shellone(:self%pftsz  ) = real(polar_here(:self%pftsz),dp)                        * A(1) ! x-part
        self%argtransf_shellone(self%pftsz+1:) = real(polar_here(self%nrots+1:self%nrots+self%pftsz),dp) * A(2) ! y-part
        ! all shells in resolution range
        self%argtransf(:self%pftsz,:)     = real(self%polar(:self%pftsz,:),dp)                          * A(1)  ! x-part
        self%argtransf(self%pftsz + 1:,:) = real(self%polar(self%nrots + 1:self%nrots+self%pftsz,:),dp) * A(2)  ! y-part
        ! allocate others
        allocate(self%pfts_refs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%pfts_refs_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%pfts_drefs_even(self%pftsz,self%kfromto(1):self%kfromto(2),3,params_glob%nthr),&
                &self%pfts_drefs_odd (self%pftsz,self%kfromto(1):self%kfromto(2),3,params_glob%nthr),&
                &self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                &self%ctfmats(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                &self%sqsums_ptcls(1:self%nptcls),self%ksqsums_ptcls(1:self%nptcls),self%wsqsums_ptcls(1:self%nptcls),&
                &self%heap_vars(params_glob%nthr))
        do ithr=1,params_glob%nthr
            allocate(self%heap_vars(ithr)%pft_ref(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_ref_tmp(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%argvec(self%pftsz),&
                    &self%heap_vars(ithr)%shvec(self%pftsz),&
                    &self%heap_vars(ithr)%shmat(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%kcorrs(self%nrots),&
                    &self%heap_vars(ithr)%pft_ref_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_ref_tmp_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_dref_8(self%pftsz,self%kfromto(1):self%kfromto(2),3),&
                    &self%heap_vars(ithr)%shmat_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_r1_8(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                    &self%heap_vars(ithr)%pft_r(self%pftsz,self%kfromto(1):self%kfromto(2)))
        end do
        self%pfts_refs_even = zero
        self%pfts_refs_odd  = zero
        self%pfts_ptcls     = zero
        self%sqsums_ptcls   = 0.d0
        self%ksqsums_ptcls  = 0.d0
        self%wsqsums_ptcls  = 0.d0
        ! CTF
        self%ctfmats = 1.0
        self%with_ctf = .false.
        if( params_glob%ctf .ne. 'no' ) self%with_ctf = .true.
        ! setup npix_per_shell
        call self%setup_npix_per_shell
        ! allocation for memoization
        call self%allocate_ptcls_memoization
        ! flag existence
        self%existence = .true.
        ! set pointer to global instance
        pftc_glob => self
    end subroutine new

    module subroutine kill(self)
        class(polarft_calc), intent(inout) :: self
        integer :: ithr
        if( self%existence )then
            do ithr=1,params_glob%nthr
                deallocate(self%heap_vars(ithr)%pft_ref,self%heap_vars(ithr)%pft_ref_tmp,&
                    &self%heap_vars(ithr)%argvec, self%heap_vars(ithr)%shvec,&
                    &self%heap_vars(ithr)%shmat,self%heap_vars(ithr)%kcorrs,&
                    &self%heap_vars(ithr)%pft_ref_8,self%heap_vars(ithr)%pft_ref_tmp_8,&
                    &self%heap_vars(ithr)%pft_dref_8,self%heap_vars(ithr)%pft_r,&
                    &self%heap_vars(ithr)%shmat_8,self%heap_vars(ithr)%pft_r1_8)
            end do
            if( allocated(self%ctfmats)        ) deallocate(self%ctfmats)
            if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
            deallocate(self%sqsums_ptcls, self%ksqsums_ptcls, self%wsqsums_ptcls, self%angtab, self%argtransf,self%pfts_ptcls,&
                &self%polar, self%pfts_refs_even, self%pfts_refs_odd, self%pfts_drefs_even, self%pfts_drefs_odd,&
                &self%iseven, self%pinds, self%heap_vars, self%argtransf_shellone)
            call self%kill_memoized_ptcls
            call self%kill_memoized_refs
            nullify(self%sigma2_noise, pftc_glob)
            self%existence = .false.
        endif
    end subroutine kill

    module subroutine reallocate_ptcls(self, nptcls, pinds)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: nptcls
        integer,             intent(in)    :: pinds(nptcls)
        integer :: i,iptcl
        self%pfromto(1) = minval(pinds)
        self%pfromto(2) = maxval(pinds)
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( self%nptcls == nptcls )then
            ! just need to update particles indexing
        else
            ! re-index & reallocate
            self%nptcls = nptcls
            if( allocated(self%sqsums_ptcls) ) deallocate(self%sqsums_ptcls)
            if( allocated(self%ksqsums_ptcls)) deallocate(self%ksqsums_ptcls)
            if( allocated(self%wsqsums_ptcls)) deallocate(self%wsqsums_ptcls)
            if( allocated(self%iseven) )       deallocate(self%iseven)
            if( allocated(self%pfts_ptcls) )   deallocate(self%pfts_ptcls)
            if( allocated(self%ctfmats) )      deallocate(self%ctfmats)
            allocate( self%pfts_ptcls(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                     &self%ctfmats(self%pftsz,self%kfromto(1):self%kfromto(2),1:self%nptcls),&
                     &self%sqsums_ptcls(1:self%nptcls),self%ksqsums_ptcls(1:self%nptcls),&
                     &self%wsqsums_ptcls(1:self%nptcls),self%iseven(1:self%nptcls))
            call self%kill_memoized_ptcls
            call self%allocate_ptcls_memoization
        endif
        self%pfts_ptcls    = zero
        self%ctfmats       = 1.0
        self%sqsums_ptcls  = 0.d0
        self%ksqsums_ptcls = 0.d0
        self%wsqsums_ptcls = 0.d0
        self%iseven        = .true.
        allocate(self%pinds(self%pfromto(1):self%pfromto(2)), source=0)
        do i = 1,self%nptcls
            iptcl = pinds(i)
            self%pinds( iptcl ) = i
        enddo
    end subroutine reallocate_ptcls

    module subroutine set_ref_pft(self, iref, pft, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        complex(sp),         intent(in)    :: pft(:,:)
        logical,             intent(in)    :: iseven
        if( iseven )then
            self%pfts_refs_even(:,:,iref) = pft
        else
            self%pfts_refs_odd(:,:,iref)  = pft
        endif
    end subroutine set_ref_pft

    module subroutine set_ptcl_pft(self, iptcl, pft)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        complex(sp),         intent(in)    :: pft(:,:)
        self%pfts_ptcls(:,:,self%pinds(iptcl)) = pft
        call self%memoize_sqsum_ptcl(iptcl)
    end subroutine set_ptcl_pft

    module subroutine set_ref_fcomp(self, iref, irot, k, comp, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot, k
        complex(sp),         intent(in)    :: comp
        logical,             intent(in)    :: iseven
        if( iseven )then
            self%pfts_refs_even(irot,k,iref) = comp
        else
            self%pfts_refs_odd(irot,k,iref)  = comp
        endif
    end subroutine set_ref_fcomp

    module subroutine set_dref_fcomp(self, iref, irot, k, dcomp, iseven)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, irot, k
        complex(sp),         intent(in)    :: dcomp(3)
        logical,             intent(in)    :: iseven
        if( iseven )then
            self%pfts_drefs_even(irot,k,:,iref) = dcomp
        else
            self%pfts_drefs_odd(irot,k,:,iref)  = dcomp
        endif
    end subroutine set_dref_fcomp

    module subroutine set_ptcl_fcomp(self, iptcl, irot, k, comp)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl, irot, k
        complex(sp),         intent(in)    :: comp
        self%pfts_ptcls(irot,k,self%pinds(iptcl)) = comp
    end subroutine set_ptcl_fcomp

    module subroutine cp_even2odd_ref(self, iref)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        self%pfts_refs_odd(:,:,iref) = self%pfts_refs_even(:,:,iref)
    end subroutine cp_even2odd_ref

    module subroutine cp_odd2even_ref(self, iref)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref
        self%pfts_refs_even(:,:,iref) = self%pfts_refs_odd(:,:,iref)
    end subroutine cp_odd2even_ref

    module subroutine cp_refs(self, self2)
        class(polarft_calc), intent(inout) :: self, self2
        self%pfts_refs_odd  = self2%pfts_refs_odd
        self%pfts_refs_even = self2%pfts_refs_even
    end subroutine cp_refs

    module subroutine cp_even_ref2ptcl(self, iref, iptcl)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        self%pfts_ptcls(:,:,self%pinds(iptcl)) = self%pfts_refs_even(:,:,iref)
        call self%memoize_sqsum_ptcl(self%pinds(iptcl))
    end subroutine cp_even_ref2ptcl

    module subroutine swap_ptclsevenodd(self)
        class(polarft_calc), intent(inout) :: self
        self%iseven = .not.self%iseven
    end subroutine swap_ptclsevenodd

    module subroutine set_eo(self, iptcl, is_even)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        logical,             intent(in)    :: is_even
        self%iseven(self%pinds(iptcl)) = is_even
    end subroutine set_eo

    module subroutine set_eos(self, eoarr)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: eoarr(self%nptcls)
        integer :: i
        if( all(eoarr == - 1) )then
            self%iseven = .true.
        else
            do i=1,self%nptcls
                if( eoarr(i) == 0 )then
                    self%iseven(i) = .true.
                else
                    self%iseven(i) = .false.
                endif
            end do
        endif
    end subroutine set_eos

    module subroutine set_with_ctf(self, l_wctf)
        class(polarft_calc), intent(inout) :: self
        logical,             intent(in)    :: l_wctf
        self%with_ctf = l_wctf
    end subroutine set_with_ctf

    module subroutine assign_sigma2_noise(self, sigma2_noise)
        class(polarft_calc),       intent(inout) :: self
        real, allocatable, target, intent(inout) :: sigma2_noise(:,:)
        self%sigma2_noise => sigma2_noise
    end subroutine assign_sigma2_noise

end submodule simple_polarft_core

!>  \brief  SIMPLE polaracf_corrcalc class
module simple_polaracf_corrcalc
use simple_defs      ! singleton
use simple_params,   only: params
use simple_jiffys,   only: alloc_err
use simple_ran_tabu, only: ran_tabu
implicit none

public :: polaracf_corrcalc
private

type :: polaracf_corrcalc
    private
    integer           :: pfromto(2) = 1      !< from/to particle indices (in parallel execution)
    integer           :: nptcls     = 1      !< the total number of particles in partition (logically indexded [fromp,top])
    integer           :: nrefs      = 1      !< the number of references (logically indexded [1,nrefs])
    integer           :: nrots      = 0      !< number of in-plane rotations for one pacf (determined by radius of molecule)
    integer           :: refsz      = 0      !< size of reference (nrots/2) (number of vectors used for matching)
    integer           :: ptclsz     = 0      !< size of particle (2*nrots)
    integer           :: ldim(3)    = 0      !< logical dimensions of original cartesian image
    integer           :: rfromto(2) = 0      !< pixel ring range
    integer           :: nrings     = 0      !< number of rings
    real, allocatable :: angtab(:)           !< table of in-plane angles (in degrees)
    real, allocatable :: polar(:,:)          !< table of polar coordinates (in Cartesian coordinates)
    real, allocatable :: pacfs_refs(:,:,:)   !< 3D real matrix of polar reference sections (nrefs,refsz,nrings)
    real, allocatable :: pacfs_ptcls(:,:,:)  !< 3D real matrix of particle sections (phase-flipped)
    logical           :: existence=.false.   !< to indicate existence
  contains
    procedure          :: new
    procedure, private :: r_ind
    procedure, private :: ptcl_ind
    procedure          :: set_ref_pacf
    procedure          :: set_ptcl_pacf
    procedure          :: print
    procedure          :: exists
    procedure          :: vis_ptcl
    procedure          :: vis_ref
    procedure          :: gencorrs
    procedure, private :: corr
    procedure          :: kill
end type polaracf_corrcalc

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new( self, nrefs, pfromto, ldim, rfromto )
        use simple_math, only: rad2deg, is_even, round2even
        class(polaracf_corrcalc), intent(inout) :: self
        integer,                  intent(in)    :: nrefs, pfromto(2), ldim(3), rfromto(2)
        integer :: alloc_stat, irot, k, r_ind
        logical :: even_dims, test(3)
        real    :: ang
        ! kill possibly pre-existing object
        call self%kill
        ! error check
        if( rfromto(2)-rfromto(1) <= 2 )then
            write(*,*) 'rfromto: ', rfromto(1), rfromto(2)
            stop 'ring range too narrow; new; simple_polaracf_corrcalc'
        endif
        if( pfromto(2)-pfromto(1)+1 < 1 )then
            write(*,*) 'pfromto: ', pfromto(1), pfromto(2)
            stop 'nptcls (# of particles) must be > 0; new; simple_polaracf_corrcalc'
        endif
        if( nrefs < 1 )then
            write(*,*) 'nrefs: ', nrefs
            stop 'nrefs (# of reference sections) must be > 0; new; simple_polaracf_corrcalc'
        endif
        if( any(ldim == 0) )then
            write(*,*) 'ldim: ', ldim
            stop 'ldim is not conforming (is zero); new; simple_polaracf'
        endif
        if( ldim(3) > 1 )then
            write(*,*) 'ldim: ', ldim
            stop '3D polaracfs are not yet supported; new; simple_polaracf_corrcalc'
        endif
        test = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        if( .not. even_dims )then
            write(*,*) 'ldim: ', ldim
            stop 'only even logical dims supported; new; simple_polaracf_corrcalc'
        endif
        ! set constants
        self%pfromto = pfromto                      !< from/to particle indices (in parallel execution)
        self%nptcls  = pfromto(2)-pfromto(1)+1      !< the total number of particles in partition (logically indexded [fromp,top])
        self%nrefs   = nrefs                        !< the number of references (logically indexded [1,nrefs])
        self%rfromto = rfromto                      !< ring range
        self%nrots   = round2even(twopi*rfromto(2)) !< number of in-plane rotations for one pacf  (determined by radius of molecule)
        self%refsz   = self%nrots/2                 !< size of reference (nrots/2) (number of vectors used for matching)
        self%ptclsz  = self%nrots*2                 !< size of particle (2*nrots)
        self%ldim    = ldim                         !< logical dimensions of original cartesian image
        self%nrings  = rfromto(2)-rfromto(1)+1      !< number of rings
        ! generate polar coordinates
        allocate( self%polar(self%ptclsz,self%nrings), self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('polar coordinate arrays; new; simple_polaracf_corrcalc', alloc_stat)
        ang = twopi/real(self%nrots)
        do irot=1,self%nrots
            self%angtab(irot) = (irot-1)*ang
            do k=self%rfromto(1),self%rfromto(2)
                r_ind = self%r_ind(k) 
                self%polar(irot,r_ind)            = cos(self%angtab(irot))*real(k) ! x-coordinate
                self%polar(irot+self%nrots,r_ind) = sin(self%angtab(irot))*real(k) ! y-coordinate
            end do
            self%angtab(irot) = rad2deg(self%angtab(irot)) ! angle (in degrees)
        end do
        ! allocate polaracfs
        allocate(  self%pacfs_refs(self%nrefs,self%refsz,self%nrings),   &
                   self%pacfs_ptcls(self%nptcls,self%ptclsz,self%nrings), stat=alloc_stat)
        call alloc_err('polaracfs; new; simple_polaracf_corrcalc', alloc_stat)
        self%pacfs_refs  = 0.
        self%pacfs_ptcls = 0.
        self%existence   = .true.
    end subroutine new

    !>  \brief  is for generating rotational correlations for autocorrelation functions
    function gencorrs( self, iref, iptcl ) result( cc )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(polaracf_corrcalc), intent(inout) :: self        !< instance
        integer,                  intent(in)    :: iref, iptcl !< ref & ptcl indices
        real    :: cc(self%refsz)
        integer :: irot
        !$omp parallel do default(shared) private(irot)
        do irot=1,self%refsz
            cc(irot) = self%corr(iref, iptcl, irot)
        end do
        !$omp end parallel do
    end function gencorrs

    !>  \brief  for calculating the correlation between reference iref and particle iptcl in rotation irot
    function corr( self, iref, iptcl, irot ) result( cc )
        use simple_math, only: calc_corr
        class(polaracf_corrcalc), intent(inout) :: self              !< instance
        integer,                  intent(in)    :: iref, iptcl, irot !< reference, particle, rotation
        real, allocatable :: diffmat1(:,:), diffmat2(:,:)
        real              :: cc,ax,ay,sxx,syy,sxy,npix 
        integer           :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        allocate(diffmat1(self%refsz,self%nrings), diffmat2(self%refsz,self%nrings))
        npix     = real(self%refsz*self%nrings)
        ax       = sum(self%pacfs_refs(iref,:,:))/npix
        ay       = sum(self%pacfs_ptcls(ptcl_ind,irot:irot+self%refsz-1,:))/npix
        diffmat1 = self%pacfs_refs(iref,:,:)-ax
        diffmat2 = self%pacfs_ptcls(ptcl_ind,irot:irot+self%refsz-1,:)-ay
        sxx      = sum(diffmat1**2.)
        syy      = sum(diffmat2**2.)
        sxy      = sum(diffmat1*diffmat2)
        deallocate(diffmat1,diffmat2)
        cc = calc_corr(sxy,sxx*syy)
    end function corr

    !>  \brief  for mapping physical index to implemented
    function r_ind( self, r_in ) result( r_out )
        class(polaracf_corrcalc), intent(in) :: self
        integer,                  intent(in) :: r_in
        integer :: r_out
        r_out = r_in-self%rfromto(1)+1
    end function
    
    !>  \brief  for mapping physical particle index to implemented
    function ptcl_ind( self, iptcl_in ) result( iptcl_out )
        class(polaracf_corrcalc), intent(in) :: self
        integer,                  intent(in) :: iptcl_in
        integer :: iptcl_out
        iptcl_out = iptcl_in-self%pfromto(1)+1
    end function
    
    !>  \brief  sets reference pacf iref
    subroutine set_ref_pacf( self, iref, acfmat )
        use simple_math, only: quadri
        class(polaracf_corrcalc), intent(inout) :: self
        integer,                  intent(in)    :: iref
        real,                     intent(in)    :: acfmat(self%ldim(1),self%ldim(2))
        integer :: irot, k, rind
        do irot=1,self%refsz
            do k=self%rfromto(1),self%rfromto(2)
                rind = self%r_ind(k) 
                self%pacfs_refs(iref,irot,rind) =&
                quadri(self%polar(irot,rind),self%polar(irot+self%nrots,rind),acfmat,self%ldim(1),self%ldim(2))
            end do
        end do
    end subroutine
    
    !>  \brief  sets particle pacf iptcl
    subroutine set_ptcl_pacf( self, iptcl, acfmat )
        use simple_math, only: quadri
        class(polaracf_corrcalc), intent(inout) :: self
        integer,                  intent(in)    :: iptcl
        real,                     intent(in)    :: acfmat(self%ldim(1),self%ldim(2))
        integer :: ptcl_ind, irot, k, rind
        ptcl_ind = self%ptcl_ind(iptcl)
        do irot=1,self%nrots
            do k=self%rfromto(1),self%rfromto(2)
                rind = self%r_ind(k) 
                self%pacfs_ptcls(ptcl_ind,irot,rind) =&
                quadri(self%polar(irot,rind),self%polar(irot+self%nrots,rind),acfmat,self%ldim(1),self%ldim(2))
            end do
        end do
        ! expand rot dim
        self%pacfs_ptcls(ptcl_ind,self%nrots+1:,:) = self%pacfs_ptcls(ptcl_ind,:self%nrots,:) 
    end subroutine
    
    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ptcl( self, iptcl )
        use gnufor2
        class(polaracf_corrcalc), intent(in) :: self
        integer,                  intent(in) :: iptcl
        integer :: ptcl_ind
        ptcl_ind = self%ptcl_ind(iptcl)
        call gnufor_image(self%pacfs_ptcls(ptcl_ind,:self%refsz,:), palette='gray')
    end subroutine
    
    !>  \brief  is for plotting a particle polar FT
    subroutine vis_ref( self, iref )
        use gnufor2
        class(polaracf_corrcalc), intent(in) :: self
        integer,                  intent(in) :: iref
        call gnufor_image(self%pacfs_refs(iref,:,:),  palette='gray')
    end subroutine
      
    !>  \brief  for printing info about the object
    subroutine print( self )
        class(polaracf_corrcalc), intent(in) :: self
        write(*,*) "from/to particle indices               (self%pfromto): ", self%pfromto
        write(*,*) "total n particles in partition          (self%nptcls): ", self%nptcls
        write(*,*) "number of references                     (self%nrefs): ", self%nrefs
        write(*,*) "number of rotations                      (self%nrots): ", self%nrots
        write(*,*) "number of rings                         (self%nrings): ", self%nrings
        write(*,*) "nr of rots for ref (2nd dim of pacfmat)  (self%refsz): ", self%refsz
        write(*,*) "n rots for ptcl (2nd dim of pacfmat)    (self%ptclsz): ", self%ptclsz 
        write(*,*) "logical dim. of original Cartesian image  (self%ldim): ", self%ldim
    end subroutine

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(polaracf_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function

    !>  \brief  is a destructor
    subroutine kill( self )
        class(polaracf_corrcalc), intent(inout) :: self
        if( self%existence )then
            deallocate(self%angtab, self%polar, self%pacfs_refs, self%pacfs_ptcls)
            self%existence = .false.
        endif
    end subroutine
    
end module simple_polaracf_corrcalc
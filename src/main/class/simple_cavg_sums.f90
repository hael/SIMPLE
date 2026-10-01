!@descr: unregularized 2D class Fourier sums on disk: the partless carry-over set and per-worker contributions
! Per class: even/odd numerators and CTF^2 sums, captured before restoration. STATE: carried set, owner-written,
! records M(c). CONTRIBUTION: one worker's from-zero sums + centering offsets, e/o pops, l_frac.
! Written via temp name + rename, closed by a payload byte count. Ownership: doc/policies/2D/abinitio2D_policy.md
module simple_cavg_sums
use simple_core_module_api
use simple_ftiter, only: ftiter
implicit none

public :: cavg_sums
public :: CAVG_SUMS_STATE, CAVG_SUMS_CONTRIB
public :: CAVG_SUMS_OK, CAVG_SUMS_MISSING, CAVG_SUMS_CORRUPT, CAVG_SUMS_MISMATCH
public :: cavg_contrib_fname
private
#include "simple_local_flags.inc"

integer,          parameter :: CAVG_SUMS_STATE    = 1
integer,          parameter :: CAVG_SUMS_CONTRIB  = 2
integer,          parameter :: CAVG_SUMS_OK       = 0       ! read status
integer,          parameter :: CAVG_SUMS_MISSING  = 1
integer,          parameter :: CAVG_SUMS_CORRUPT  = 2
integer,          parameter :: CAVG_SUMS_MISMATCH = 3
integer,          parameter :: FORMAT_VERSION     = 1
character(len=16),parameter :: MAGIC              = 'SIMPLE_CAVGSUMS1'
character(len=*), parameter :: TMP_SUFFIX         = '.tmp'

type :: cavg_sums
    private
    integer                  :: kind      = 0
    integer                  :: ncls      = 0
    integer                  :: box_crop  = 0
    integer                  :: cshape(2) = 0
    integer                  :: part      = 0
    real                     :: smpd_crop = 0.
    logical                  :: l_frac    = .false.
    real,        allocatable :: mrep(:)          ! STATE: represented population M(c)
    real,        allocatable :: offsets(:,:)     ! CONTRIBUTION: centering offsets (2,ncls)
    integer,     allocatable :: eo_pops(:,:)     ! CONTRIBUTION: accumulated even/odd counts (2,ncls)
    complex(sp), allocatable :: cmat(:,:,:,:)    ! numerators (cshape(1),cshape(2),ncls,even:odd)
    real,        allocatable :: ctfsq(:,:,:,:)   ! CTF^2 sums, same shape
    logical                  :: exists    = .false.
  contains
    procedure          :: new
    procedure          :: zero
    procedure          :: set_sums
    procedure          :: get_sums
    procedure          :: set_mrep
    procedure          :: get_mrep
    procedure          :: set_contrib_meta
    procedure          :: get_offsets
    procedure          :: get_eo_pops
    procedure          :: get_l_frac
    procedure          :: get_part
    procedure          :: matches
    procedure          :: accumulate
    procedure          :: shift_class
    procedure          :: blend
    procedure          :: pad_to
    procedure          :: class_mass
    procedure          :: write
    procedure          :: read
    procedure          :: kill
end type cavg_sums

contains

    !> Name of the contribution file of one worker part (no dependence on numlen); the
    !! carried set is CAVG_STATE_FILE (simple_defs_fname)
    function cavg_contrib_fname( part ) result( fname )
        integer, intent(in) :: part
        type(string) :: fname
        fname = CAVG_CONTRIB_FBODY//int2str(part)//BIN_EXT
    end function cavg_contrib_fname

    !> Zeroed sums of a given kind for ncls classes of box_crop x box_crop at smpd_crop
    subroutine new( self, kind, ncls, box_crop, smpd_crop, part )
        class(cavg_sums),  intent(inout) :: self
        integer,           intent(in)    :: kind, ncls, box_crop
        real,              intent(in)    :: smpd_crop
        integer, optional, intent(in)    :: part
        call self%kill
        if( kind /= CAVG_SUMS_STATE .and. kind /= CAVG_SUMS_CONTRIB ) THROW_HARD('unknown cavg_sums kind')
        if( ncls < 1 .or. box_crop < 2 ) THROW_HARD('invalid cavg_sums dimensions')
        self%kind      = kind
        self%ncls      = ncls
        self%box_crop  = box_crop
        self%smpd_crop = smpd_crop
        self%cshape    = [fdim(box_crop), box_crop]
        self%part      = 0
        if( present(part) ) self%part = part
        allocate(self%cmat(self%cshape(1),self%cshape(2),ncls,2), source=cmplx(0.,0.,kind=sp))
        allocate(self%ctfsq(self%cshape(1),self%cshape(2),ncls,2), source=0.)
        allocate(self%mrep(ncls), source=0.)
        allocate(self%offsets(2,ncls), source=0.)
        allocate(self%eo_pops(2,ncls), source=0)
        self%exists = .true.
    end subroutine new

    subroutine zero( self )
        class(cavg_sums), intent(inout) :: self
        if( .not. self%exists ) return
        self%cmat    = cmplx(0.,0.,kind=sp)
        self%ctfsq   = 0.
        self%mrep    = 0.
        self%offsets = 0.
        self%eo_pops = 0
        self%l_frac  = .false.
    end subroutine zero

    !> Copy the four accumulators in (even/odd numerators and CTF^2 sums)
    subroutine set_sums( self, cmat_e, cmat_o, ctfsq_e, ctfsq_o )
        class(cavg_sums), intent(inout) :: self
        complex(sp),      intent(in)    :: cmat_e(:,:,:), cmat_o(:,:,:)
        real,             intent(in)    :: ctfsq_e(:,:,:), ctfsq_o(:,:,:)
        call check_shape(self, shape(cmat_e));  call check_shape(self, shape(cmat_o))
        call check_shape(self, shape(ctfsq_e)); call check_shape(self, shape(ctfsq_o))
        self%cmat(:,:,:,1)  = cmat_e
        self%cmat(:,:,:,2)  = cmat_o
        self%ctfsq(:,:,:,1) = ctfsq_e
        self%ctfsq(:,:,:,2) = ctfsq_o
    end subroutine set_sums

    !> Copy the four accumulators out
    subroutine get_sums( self, cmat_e, cmat_o, ctfsq_e, ctfsq_o )
        class(cavg_sums), intent(in)    :: self
        complex(sp),      intent(inout) :: cmat_e(:,:,:), cmat_o(:,:,:)
        real,             intent(inout) :: ctfsq_e(:,:,:), ctfsq_o(:,:,:)
        call check_shape(self, shape(cmat_e));  call check_shape(self, shape(cmat_o))
        call check_shape(self, shape(ctfsq_e)); call check_shape(self, shape(ctfsq_o))
        cmat_e  = self%cmat(:,:,:,1)
        cmat_o  = self%cmat(:,:,:,2)
        ctfsq_e = self%ctfsq(:,:,:,1)
        ctfsq_o = self%ctfsq(:,:,:,2)
    end subroutine get_sums

    subroutine set_mrep( self, mrep )
        class(cavg_sums), intent(inout) :: self
        real,             intent(in)    :: mrep(:)
        if( size(mrep) /= self%ncls ) THROW_HARD('mrep size does not match the class count')
        self%mrep = mrep
    end subroutine set_mrep

    subroutine get_mrep( self, mrep )
        class(cavg_sums),  intent(in)    :: self
        real, allocatable, intent(inout) :: mrep(:)
        if( allocated(mrep) ) deallocate(mrep)
        allocate(mrep, source=self%mrep)
    end subroutine get_mrep

    !> Metadata of a worker contribution: centering offsets, accumulated populations and
    !! whether the iteration blends carry-over
    subroutine set_contrib_meta( self, offsets, eo_pops, l_frac )
        class(cavg_sums), intent(inout) :: self
        real,             intent(in)    :: offsets(:,:)
        integer,          intent(in)    :: eo_pops(:,:)
        logical,          intent(in)    :: l_frac
        if( size(offsets,2) /= self%ncls .or. size(eo_pops,2) /= self%ncls )then
            THROW_HARD('contribution metadata does not match the class count')
        endif
        self%offsets = offsets
        self%eo_pops = eo_pops
        self%l_frac  = l_frac
    end subroutine set_contrib_meta

    subroutine get_offsets( self, offsets )
        class(cavg_sums),  intent(in)    :: self
        real, allocatable, intent(inout) :: offsets(:,:)
        if( allocated(offsets) ) deallocate(offsets)
        allocate(offsets, source=self%offsets)
    end subroutine get_offsets

    subroutine get_eo_pops( self, eo_pops )
        class(cavg_sums),     intent(in)    :: self
        integer, allocatable, intent(inout) :: eo_pops(:,:)
        if( allocated(eo_pops) ) deallocate(eo_pops)
        allocate(eo_pops, source=self%eo_pops)
    end subroutine get_eo_pops

    logical function get_l_frac( self )
        class(cavg_sums), intent(in) :: self
        get_l_frac = self%l_frac
    end function get_l_frac

    integer function get_part( self )
        class(cavg_sums), intent(in) :: self
        get_part = self%part
    end function get_part

    !> Whether the sums describe the run's classes: class count, box and sampling
    logical function matches( self, ncls, box_crop, smpd_crop )
        class(cavg_sums), intent(in) :: self
        integer,          intent(in) :: ncls, box_crop
        real,             intent(in) :: smpd_crop
        matches = self%exists
        if( .not. matches ) return
        matches = self%ncls == ncls .and. self%box_crop == box_crop .and. &
            &abs(self%smpd_crop - smpd_crop) <= 1.e-4 * max(smpd_crop, 1.e-6)
    end function matches

    !> Add another set's sums (and contribution populations): the reduction of worker parts
    subroutine accumulate( self, other )
        class(cavg_sums), intent(inout) :: self
        class(cavg_sums), intent(in)    :: other
        if( other%ncls /= self%ncls .or. any(other%cshape /= self%cshape) )then
            THROW_HARD('cannot accumulate class sums of different geometry')
        endif
        self%cmat    = self%cmat    + other%cmat
        self%ctfsq   = self%ctfsq   + other%ctfsq
        self%eo_pops = self%eo_pops + other%eo_pops
    end subroutine accumulate

    !> Fourier phase shift of both numerators of one class by offset (pixels at box_crop),
    !! the class-centering shift of the previous set; CTF^2 sums are shift invariant
    subroutine shift_class( self, icls, offset )
        class(cavg_sums), intent(inout) :: self
        integer,          intent(in)    :: icls
        real,             intent(in)    :: offset(2)
        type(ftiter) :: fit
        complex(dp), allocatable :: ph_h(:), ph_k(:)
        complex(dp) :: w1, w2
        real(dp)    :: sh(2)
        integer     :: lims(3,2), h, k, hphys, kphys, ieo
        if( icls < 1 .or. icls > self%ncls ) THROW_HARD('class index out of range; shift_class')
        if( all(abs(offset) < TINY) ) return
        fit   = ftiter([self%box_crop, self%box_crop, 1], 1.0)
        lims  = fit%loop_lims(2)
        sh(1) = real(offset(1) * PI / real(self%box_crop/2), dp)
        sh(2) = real(offset(2) * PI / real(self%box_crop/2), dp)
        w1    = cmplx(cos(sh(1)), sin(sh(1)), kind=dp)
        w2    = cmplx(cos(sh(2)), sin(sh(2)), kind=dp)
        allocate(ph_h(lims(1,1):lims(1,2)), ph_k(lims(2,1):lims(2,2)))
        ph_h(lims(1,1)) = cmplx(cos(real(lims(1,1),dp)*sh(1)), sin(real(lims(1,1),dp)*sh(1)), kind=dp)
        do h = lims(1,1)+1, lims(1,2)
            ph_h(h) = ph_h(h-1) * w1
        enddo
        ph_k(lims(2,1)) = cmplx(cos(real(lims(2,1),dp)*sh(2)), sin(real(lims(2,1),dp)*sh(2)), kind=dp)
        do k = lims(2,1)+1, lims(2,2)
            ph_k(k) = ph_k(k-1) * w2
        enddo
        do ieo = 1, 2
            do k = lims(2,1), lims(2,2)
                kphys = k + 1 + merge(self%box_crop, 0, k < 0)
                do h = lims(1,1), lims(1,2)
                    hphys = h + 1
                    self%cmat(hphys,kphys,icls,ieo) = self%cmat(hphys,kphys,icls,ieo) * cmplx(ph_k(k) * ph_h(h), kind=sp)
                enddo
            enddo
        enddo
        deallocate(ph_h, ph_k)
    end subroutine shift_class

    !> self = s(c) * self + w(c) * prev for every class, all four accumulators
    subroutine blend( self, prev, s, w )
        class(cavg_sums), intent(inout) :: self
        class(cavg_sums), intent(in)    :: prev
        real,             intent(in)    :: s(:), w(:)
        integer :: icls
        if( prev%ncls /= self%ncls .or. any(prev%cshape /= self%cshape) )then
            THROW_HARD('cannot blend class sums of different geometry')
        endif
        if( size(s) /= self%ncls .or. size(w) /= self%ncls ) THROW_HARD('blend weights do not match the class count')
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls = 1, self%ncls
            self%cmat(:,:,icls,:)  = s(icls) * self%cmat(:,:,icls,:)  + w(icls) * prev%cmat(:,:,icls,:)
            self%ctfsq(:,:,icls,:) = s(icls) * self%ctfsq(:,:,icls,:) + w(icls) * prev%ctfsq(:,:,icls,:)
        enddo
        !$omp end parallel do
    end subroutine blend

    !> Fourier-pad the sums to a larger box of the same physical extent (the streaming pool's
    !! crop-box upsample); M(c) is unchanged
    subroutine pad_to( self, box_crop, smpd_crop )
        class(cavg_sums), intent(inout) :: self
        integer,          intent(in)    :: box_crop
        real,             intent(in)    :: smpd_crop
        complex(sp), allocatable :: cmat(:,:,:,:)
        real,        allocatable :: ctfsq(:,:,:,:)
        type(ftiter) :: fit_old, fit_new
        integer      :: flims(3,2), h, k, phys(2), phys_pd(2), cshape_new(2)
        if( box_crop == self%box_crop ) return
        if( box_crop <  self%box_crop ) THROW_HARD('Fourier cropping of class sums is not supported; pad_to')
        fit_old    = ftiter([self%box_crop, self%box_crop, 1], 1.0)
        fit_new    = ftiter([box_crop, box_crop, 1], 1.0)
        flims      = fit_old%loop_lims(2)
        cshape_new = [fdim(box_crop), box_crop]
        allocate(cmat(cshape_new(1),cshape_new(2),self%ncls,2), source=cmplx(0.,0.,kind=sp))
        allocate(ctfsq(cshape_new(1),cshape_new(2),self%ncls,2), source=0.)
        do h = flims(1,1), flims(1,2)
            do k = flims(2,1), flims(2,2)
                phys    = fit_old%comp_addr_phys(h,k)
                phys_pd = fit_new%comp_addr_phys(h,k)
                cmat(phys_pd(1),phys_pd(2),:,:)  = self%cmat(phys(1),phys(2),:,:)
                ctfsq(phys_pd(1),phys_pd(2),:,:) = self%ctfsq(phys(1),phys(2),:,:)
            end do
        end do
        call move_alloc(cmat,  self%cmat)
        call move_alloc(ctfsq, self%ctfsq)
        self%box_crop  = box_crop
        self%smpd_crop = smpd_crop
        self%cshape    = cshape_new
    end subroutine pad_to

    !> Sampling mass of each class: the sum of its even and odd CTF^2 accumulators
    subroutine class_mass( self, mass )
        class(cavg_sums),      intent(in)    :: self
        real(dp), allocatable, intent(inout) :: mass(:)
        integer :: icls
        if( allocated(mass) ) deallocate(mass)
        allocate(mass(self%ncls), source=0._dp)
        do icls = 1, self%ncls
            mass(icls) = sum(real(self%ctfsq(:,:,icls,:), dp))
        enddo
    end subroutine class_mass

    !> Publish the set: write to a temporary name, then rename into place
    subroutine write( self, fname )
        class(cavg_sums), intent(in) :: self
        class(string),    intent(in) :: fname
        type(string)    :: tmpname
        integer(kind=8) :: nbytes
        integer         :: funit, ios
        if( .not. self%exists ) THROW_HARD('cannot write empty class sums')
        tmpname = fname//TMP_SUFFIX
        call del_file(tmpname)
        call fopen(funit, tmpname, status='REPLACE', action='WRITE', access='STREAM', iostat=ios)
        call fileiochk('cavg_sums write: '//tmpname%to_char(), ios)
        write(funit, iostat=ios) MAGIC, FORMAT_VERSION, self%kind, self%ncls, self%box_crop, self%cshape, &
            &self%part, merge(1, 0, self%l_frac), self%smpd_crop
        if( ios == 0 ) write(funit, iostat=ios) self%mrep, self%offsets, self%eo_pops
        if( ios == 0 ) write(funit, iostat=ios) self%cmat, self%ctfsq
        nbytes = payload_bytes(self%ncls, self%cshape)
        if( ios == 0 ) write(funit, iostat=ios) nbytes
        call fileiochk('cavg_sums write: '//tmpname%to_char(), ios)
        call fclose(funit)
        call simple_rename(tmpname, fname, overwrite=.false.)
        call tmpname%kill
    end subroutine write

    !> Read a set. status: CAVG_SUMS_OK; CAVG_SUMS_MISSING (no file); CAVG_SUMS_CORRUPT
    !! (unreadable, truncated, foreign or other format version); CAVG_SUMS_MISMATCH (another kind).
    !! A leftover temporary file is never read. On any status other than OK the object is empty.
    subroutine read( self, fname, kind, status )
        class(cavg_sums), intent(inout) :: self
        class(string),    intent(in)    :: fname
        integer,          intent(in)    :: kind
        integer,          intent(out)   :: status
        character(len=16) :: magic_read
        integer(kind=8)   :: nbytes, fsize
        integer :: funit, ios, version, kind_read, ncls, box_crop, cshape(2), part, ifrac
        real    :: smpd_crop
        call self%kill
        status = CAVG_SUMS_MISSING
        if( .not. file_exists(fname) ) return
        status = CAVG_SUMS_CORRUPT
        inquire(file=fname%to_char(), size=fsize)
        call fopen(funit, fname, status='OLD', action='READ', access='STREAM', iostat=ios)
        if( ios /= 0 ) return
        read(funit, iostat=ios) magic_read, version, kind_read, ncls, box_crop, cshape, part, ifrac, smpd_crop
        if( ios /= 0 .or. magic_read /= MAGIC .or. version /= FORMAT_VERSION )then
            call fclose(funit)
            return
        endif
        if( ncls < 1 .or. box_crop < 2 .or. any(cshape /= [fdim(box_crop), box_crop]) )then
            call fclose(funit)
            return
        endif
        if( fsize /= header_bytes() + payload_bytes(ncls, cshape) + 8_8 )then
            call fclose(funit)
            return
        endif
        if( kind_read /= kind )then
            call fclose(funit)
            status = CAVG_SUMS_MISMATCH
            return
        endif
        call self%new(kind_read, ncls, box_crop, smpd_crop, part)
        self%l_frac = ifrac == 1
        read(funit, iostat=ios) self%mrep, self%offsets, self%eo_pops
        if( ios == 0 ) read(funit, iostat=ios) self%cmat, self%ctfsq
        if( ios == 0 ) read(funit, iostat=ios) nbytes
        call fclose(funit)
        if( ios /= 0 .or. nbytes /= payload_bytes(ncls, cshape) )then
            call self%kill
            return
        endif
        status = CAVG_SUMS_OK
    end subroutine read

    subroutine kill( self )
        class(cavg_sums), intent(inout) :: self
        if( allocated(self%mrep)    ) deallocate(self%mrep)
        if( allocated(self%offsets) ) deallocate(self%offsets)
        if( allocated(self%eo_pops) ) deallocate(self%eo_pops)
        if( allocated(self%cmat)    ) deallocate(self%cmat)
        if( allocated(self%ctfsq)   ) deallocate(self%ctfsq)
        self%kind      = 0
        self%ncls      = 0
        self%box_crop  = 0
        self%cshape    = 0
        self%part      = 0
        self%smpd_crop = 0.
        self%l_frac    = .false.
        self%exists    = .false.
    end subroutine kill

    ! private helpers

    subroutine check_shape( self, shp )
        class(cavg_sums), intent(in) :: self
        integer,          intent(in) :: shp(:)
        if( size(shp) /= 3 ) THROW_HARD('class sums arrays are rank 3')
        if( any(shp /= [self%cshape(1), self%cshape(2), self%ncls]) ) THROW_HARD('class sums array shape mismatch')
    end subroutine check_shape

    ! bytes of the fixed header: magic, 8 integers, one real
    integer(kind=8) function header_bytes()
        header_bytes = 16_8 + 8_8 * 4_8 + 4_8
    end function header_bytes

    ! bytes of the per-class metadata and the four arrays
    integer(kind=8) function payload_bytes( ncls, cshape )
        integer, intent(in) :: ncls, cshape(2)
        integer(kind=8) :: nvox
        nvox = int(cshape(1),8) * int(cshape(2),8) * int(ncls,8)
        payload_bytes = 4_8 * int(ncls,8) * 5_8 &   ! mrep, offsets(2,:), eo_pops(2,:)
            &         + 8_8 * nvox * 2_8        &   ! cmat even/odd
            &         + 4_8 * nvox * 2_8            ! ctfsq even/odd
    end function payload_bytes

end module simple_cavg_sums

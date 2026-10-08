!@descr: flex_pca embedding artifact: the raw embedding and the state stage's trailing block, published atomically
!!
!! One binary file per run (flex_pca_embedding.bin): the raw MAP embedding with its posterior
!! moments and provenance (the working lattice), followed by a trailing block once the state stage
!! has calibrated the noise scale: the noise scale and, after the latent deconvolution, the
!! deconvolved coordinates, precisions and mixture labels. Every write goes to a temporary name that
!! is flushed and then atomically renamed over the artifact (the destination is never deleted first),
!! so a crash leaves the previous complete artifact, never a torn one. A states-only resume (infile=) reads the raw block and adopts the trailing block when present.
module simple_flex_pca_embedding_io
use simple_core_module_api, only: del_file, dp, file_exists, logfhandle, simple_copy_file, simple_exception, string
use simple_syslib,          only: simple_atomic_replace, simple_sync_file
use simple_flex_pca_records, only: flex_selection, flex_fit_model, flex_latent
implicit none
private
#include "simple_local_flags.inc"

public :: write_embedding_cache, read_embedding_cache
public :: write_noise_scale, read_noise_scale, write_deconv_block, read_deconv_block

character(len=8), parameter :: COV_CACHE_MAGIC   = 'SIMPLFXC'
! Bump whenever the artifact layout changes; read_embedding_cache rejects any other version.
integer,          parameter :: COV_CACHE_VERSION = 5   ! 5: trailing block = noise scale, then the optional deconvolution
!> trailing block of the SAME file: the noise scale, then (flag 1) the deconvolved coordinates and
!! precisions and the mixture labels; a resume adopts it from infile, so no second file and no K ladder
character(len=8), parameter :: COV_TRAIL_MAGIC   = 'SIMPLFXD'
character(len=*), parameter :: PART_SUFFIX       = '.part'

contains

    !> Binary cache of the state stages' inputs, so a different state count / bandwidth / placement can
    !! be tried without re-fitting the basis and re-embedding every particle.
    subroutine write_embedding_cache( fname, box_crop, smpd_crop, sel, model, latent )
        type(flex_selection), intent(in) :: sel
        type(flex_fit_model), intent(in) :: model
        type(flex_latent),    intent(in) :: latent
        character(len=*),     intent(in) :: fname
        integer,              intent(in) :: box_crop          !< the working lattice the embedding was made on (provenance)
        real,                 intent(in) :: smpd_crop
        character(len=:), allocatable :: tmp
        integer :: u
        tmp = trim(fname)//PART_SUFFIX
        call del_file(tmp)
        open(newunit=u, file=tmp, status='replace', action='write', access='stream', form='unformatted')
        write(u) COV_CACHE_MAGIC, COV_CACHE_VERSION
        write(u) sel%nptcls, model%ncomp
        write(u) box_crop, smpd_crop
        write(u) sel%pinds(1:sel%nptcls)
        write(u) latent%z
        write(u) model%eigvals
        write(u) latent%contrast
        write(u) latent%resid_energy
        write(u) latent%resid_mean_energy
        write(u) latent%precision
        write(u) model%sig2_eff
        close(u)
        call simple_sync_file(tmp)
        call simple_atomic_replace(tmp, fname)
        write(logfhandle,'(A,A,A,F8.1,A)') '>>> FLEX_PCA embedding cached to ',trim(fname), &
            &' (',real(8*(sel%nptcls*model%ncomp + model%ncomp*model%ncomp*sel%nptcls))/1048576.0,' MB); reuse with infile=<path>'
        call flush(logfhandle)
    end subroutine write_embedding_cache

    subroutine read_embedding_cache( fname, box_crop, smpd_crop, sel, model, latent )
        type(flex_selection), intent(in)    :: sel
        type(flex_fit_model), intent(inout) :: model     !< ncomp, eigvals, sig2_eff read
        type(flex_latent),    intent(inout) :: latent    !< z, contrast, residual energies, precision read
        character(len=*),     intent(in)    :: fname
        integer,              intent(in)    :: box_crop  !< the run's working lattice; 0 skips the lattice check
        real,                 intent(in)    :: smpd_crop
        integer, allocatable :: pinds_cached(:)
        character(len=len(COV_CACHE_MAGIC)) :: magic
        integer :: u, ver, nptcls_c, i, box_c
        real    :: smpd_c
        if( allocated(latent%z) ) deallocate(latent%z)
        if( allocated(model%eigvals) ) deallocate(model%eigvals)
        if( allocated(latent%contrast) ) deallocate(latent%contrast)
        if( allocated(latent%resid_energy) ) deallocate(latent%resid_energy)
        if( allocated(latent%resid_mean_energy) ) deallocate(latent%resid_mean_energy)
        if( allocated(latent%precision) ) deallocate(latent%precision)
        if( .not. file_exists(fname) ) THROW_HARD('flex_pca embedding cache not found: '//trim(fname))
        open(newunit=u, file=fname, status='old', action='read', access='stream', form='unformatted')
        read(u) magic, ver
        if( magic /= COV_CACHE_MAGIC ) THROW_HARD('not a flex_pca embedding cache: '//trim(fname))
        if( ver /= COV_CACHE_VERSION ) THROW_HARD('flex_pca embedding cache version mismatch; re-run the fit')
        read(u) nptcls_c, model%ncomp
        if( nptcls_c /= sel%nptcls ) THROW_HARD('flex_pca embedding cache particle count does not match the project')
        read(u) box_c, smpd_c
        if( box_crop > 0 .and. (box_c /= box_crop .or. &
            &abs(smpd_c - smpd_crop) > 1.e-3*smpd_crop) )then
            write(logfhandle,'(A,I0,A,F7.3,A,I0,A,F7.3,A)') '>>> FLEX_PCA embedding cache lattice box_crop=', &
                &box_c, ' @ ', smpd_c, ' A; this run: box_crop=', box_crop, ' @ ', smpd_crop, ' A'
            THROW_HARD('flex_pca embedding cache was built at a different working sampling (box_crop/lp); resume at the cached lattice or re-run the fit')
        endif
        allocate(pinds_cached(nptcls_c))
        read(u) pinds_cached
        do i = 1, sel%nptcls
            if( pinds_cached(i) /= sel%pinds(i) ) &
                &THROW_HARD('flex_pca embedding cache was built from a different particle selection')
        end do
        allocate(latent%z(sel%nptcls,model%ncomp), model%eigvals(model%ncomp), latent%contrast(sel%nptcls), latent%resid_energy(sel%nptcls), &
            &latent%resid_mean_energy(sel%nptcls), latent%precision(model%ncomp,model%ncomp,sel%nptcls))
        read(u) latent%z
        read(u) model%eigvals
        read(u) latent%contrast
        read(u) latent%resid_energy
        read(u) latent%resid_mean_energy
        read(u) latent%precision
        read(u) model%sig2_eff
        close(u)
        deallocate(pinds_cached)
    end subroutine read_embedding_cache

    !> Byte offset of the first byte after the raw block of an artifact written by write_embedding_cache
    !! (stream access: magic(8) ver(4) nptcls(4) ncomp(4) box(4) smpd(4) pinds(4n) z(8nk) eigvals(8k)
    !! contrast(8n) resid(8n) resid_mean(8n) precision(8kkn) sig2(8)).
    pure function raw_block_end( nptcls, ncomp ) result( pos )
        integer, intent(in) :: nptcls, ncomp
        integer(kind=8) :: pos, n, k
        n = int(nptcls, 8); k = int(ncomp, 8)
        pos = 8_8 + 4_8 + 4_8 + 4_8 + 4_8 + 4_8 + 4_8*n + 8_8*n*k + 8_8*k + 3_8*8_8*n + 8_8*k*k*n + 8_8 + 1_8
    end function raw_block_end

    !> The calibrated noise scale as the artifact's trailing block, replacing any earlier one; a resume
    !! that finds no deconvolution re-runs it at this scale
    subroutine write_noise_scale( fname, nptcls, ncomp, noise_scale )
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: nptcls, ncomp
        real(dp),         intent(in) :: noise_scale
        call write_trailing_block(fname, nptcls, ncomp, noise_scale)
    end subroutine write_noise_scale

    !> The deconvolved coordinates, precisions, mixture labels (0 when absent) and noise scale (1 when
    !! absent) as the artifact's trailing block, replacing any earlier one
    subroutine write_deconv_block( fname, nptcls, ncomp, z, precision, labels, noise_scale )
        character(len=*),   intent(in) :: fname
        integer,            intent(in) :: nptcls, ncomp
        real(dp),           intent(in) :: z(:,:), precision(:,:,:)
        integer,  optional, intent(in) :: labels(:)
        real(dp), optional, intent(in) :: noise_scale
        integer, allocatable :: lab(:)
        real(dp) :: ns
        allocate(lab(nptcls), source=0)
        if( present(labels) ) lab(1:nptcls) = labels(1:nptcls)
        ns = 1.d0
        if( present(noise_scale) ) ns = noise_scale
        call write_trailing_block(fname, nptcls, ncomp, ns, z, precision, lab)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA deconvolved coordinates, precisions and labels added to ', trim(fname)
        call flush(logfhandle)
        deallocate(lab)
    end subroutine write_deconv_block

    !> A copy of the artifact with the new trailing block after its raw block, renamed over the artifact
    subroutine write_trailing_block( fname, nptcls, ncomp, noise_scale, z, precision, labels )
        character(len=*),   intent(in) :: fname
        integer,            intent(in) :: nptcls, ncomp
        real(dp),           intent(in) :: noise_scale
        real(dp), optional, intent(in) :: z(:,:), precision(:,:,:)
        integer,  optional, intent(in) :: labels(:)
        character(len=:), allocatable :: tmp
        integer :: u, io, l_deconv
        if( .not. file_exists(fname) ) THROW_HARD('no flex_pca embedding artifact to extend: '//trim(fname))
        tmp = trim(fname)//PART_SUFFIX
        call simple_copy_file(string(fname), string(tmp))
        open(newunit=u, file=tmp, status='old', action='readwrite', access='stream', form='unformatted', iostat=io)
        if( io /= 0 ) THROW_HARD('cannot open the flex_pca embedding artifact copy '//tmp)
        l_deconv = 0
        if( present(z) .and. present(precision) .and. present(labels) ) l_deconv = 1
        write(u, pos=raw_block_end(nptcls, ncomp)) COV_TRAIL_MAGIC
        write(u) noise_scale, l_deconv
        if( l_deconv == 1 )then
            write(u) nptcls, ncomp
            write(u) z(1:nptcls,1:ncomp)
            write(u) precision(1:ncomp,1:ncomp,1:nptcls)
            write(u) labels(1:nptcls)
        endif
        endfile(u)
        close(u)
        call simple_sync_file(tmp)
        call simple_atomic_replace(tmp, fname)
    end subroutine write_trailing_block

    !> Open the artifact and position it at its trailing block: found=.false. when the file is absent,
    !! of another layout or (nptcls, ncomp), or has no trailing block yet
    subroutine open_trailing_block( fname, nptcls, ncomp, u, found )
        character(len=*), intent(in)  :: fname
        integer,          intent(in)  :: nptcls, ncomp
        integer,          intent(out) :: u
        logical,          intent(out) :: found
        character(len=len(COV_CACHE_MAGIC)) :: magic0
        character(len=len(COV_TRAIL_MAGIC)) :: magic
        integer :: io, ver, n_c, k_c
        found = .false.
        u     = 0
        ver   = 0
        n_c   = -1
        k_c   = -1
        if( .not. file_exists(fname) ) return
        open(newunit=u, file=fname, status='old', action='read', access='stream', form='unformatted', iostat=io)
        if( io /= 0 )then
            u = 0
            return
        endif
        read(u, iostat=io) magic0, ver
        if( io == 0 ) found = magic0 == COV_CACHE_MAGIC .and. ver == COV_CACHE_VERSION
        if( found ) read(u, iostat=io) n_c, k_c
        if( found ) found = io == 0 .and. n_c == nptcls .and. k_c == ncomp
        if( found ) read(u, pos=raw_block_end(nptcls, ncomp), iostat=io) magic
        if( found ) found = io == 0 .and. magic == COV_TRAIL_MAGIC
        if( .not. found )then
            close(u)
            u = 0
        endif
    end subroutine open_trailing_block

    !> The noise scale of the artifact's trailing block; found=.false. (and 1.0) when it has none
    subroutine read_noise_scale( fname, nptcls, ncomp, noise_scale, found )
        character(len=*), intent(in)  :: fname
        integer,          intent(in)  :: nptcls, ncomp
        real(dp),         intent(out) :: noise_scale
        logical,          intent(out) :: found
        integer :: u, io
        noise_scale = 1.d0
        call open_trailing_block(fname, nptcls, ncomp, u, found)
        if( .not. found ) return
        read(u, iostat=io) noise_scale
        close(u)
        found = io == 0
        if( .not. found ) noise_scale = 1.d0
    end subroutine read_noise_scale

    !> The deconvolution of the artifact's trailing block if present and consistent with (nptcls, ncomp);
    !! found=.false. otherwise (no trailing block, or the noise scale alone)
    subroutine read_deconv_block( fname, nptcls, ncomp, z, precision, labels, noise_scale, found )
        character(len=*),      intent(in)  :: fname
        integer,               intent(in)  :: nptcls, ncomp
        real(dp), allocatable, intent(out) :: z(:,:), precision(:,:,:)
        integer,  allocatable, intent(out) :: labels(:)
        real(dp),              intent(out) :: noise_scale
        logical,               intent(out) :: found
        integer :: u, io, l_deconv, n_b, k_b
        noise_scale = 1.d0
        l_deconv    = 0
        n_b         = -1
        k_b         = -1
        call open_trailing_block(fname, nptcls, ncomp, u, found)
        if( .not. found ) return
        found = .false.
        read(u, iostat=io) noise_scale, l_deconv
        if( io == 0 .and. l_deconv == 1 ) read(u, iostat=io) n_b, k_b
        if( io /= 0 .or. l_deconv /= 1 .or. n_b /= nptcls .or. k_b /= ncomp )then
            close(u)
            return
        endif
        allocate(z(nptcls,ncomp), precision(ncomp,ncomp,nptcls), labels(nptcls))
        read(u, iostat=io) z
        if( io == 0 ) read(u, iostat=io) precision
        if( io == 0 ) read(u, iostat=io) labels
        close(u)
        if( io /= 0 )then
            deallocate(z, precision, labels)
            return
        endif
        found = .true.
    end subroutine read_deconv_block

end module simple_flex_pca_embedding_io

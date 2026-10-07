!@descr: flex_pca embedding persistence: the resume cache and its deconvolved block
!!
!! One binary file per run (flex_pca_embedding.bin): the raw MAP embedding with its posterior
!! moments and provenance (the working lattice), optionally followed by the deconvolved
!! coordinates, precisions, mixture labels and noise scale as a trailing block. A states-only
!! resume (infile=) reads the raw block and adopts the trailing block when present.
module simple_flex_pca_embedding_io
use simple_core_module_api, only: del_file, dp, file_exists, logfhandle, simple_exception
use simple_flex_pca_records, only: flex_selection, flex_fit_model, flex_latent
implicit none
private
#include "simple_local_flags.inc"

public :: write_embedding_cache, read_embedding_cache
public :: append_deconv_block, read_deconv_block

character(len=8), parameter :: COV_CACHE_MAGIC   = 'SIMPLFXC'
! Bump whenever the cache layout changes; read_embedding_cache rejects any other version.
integer,          parameter :: COV_CACHE_VERSION = 4   ! 4: box_crop/smpd_crop provenance after (nptcls, ncomp)
!> optional trailing block of the SAME file: the deconvolved coordinates and precisions, the mixture
!! labels and the noise scale; a resume adopts it from infile, so no second cache file and no K ladder
character(len=8), parameter :: COV_DECONV_MAGIC  = 'SIMPLFXD'

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
        integer :: u
        call del_file(fname)
        open(newunit=u, file=fname, status='replace', action='write', access='stream', form='unformatted')
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

    !> Byte offset of the first byte after the raw block of a cache written by write_embedding_cache
    !! (stream access: magic(8) ver(4) nptcls(4) ncomp(4) box(4) smpd(4) pinds(4n) z(8nk) eigvals(8k)
    !! contrast(8n) resid(8n) resid_mean(8n) precision(8kkn) sig2(8)).
    pure function raw_block_end( nptcls, ncomp ) result( pos )
        integer, intent(in) :: nptcls, ncomp
        integer(kind=8) :: pos, n, k
        n = int(nptcls, 8); k = int(ncomp, 8)
        pos = 8_8 + 4_8 + 4_8 + 4_8 + 4_8 + 4_8 + 4_8*n + 8_8*n*k + 8_8*k + 3_8*8_8*n + 8_8*k*k*n + 8_8 + 1_8
    end function raw_block_end

    !> Append the deconvolved block to an existing raw cache. Anything already after the raw block
    !! (an older deconvolved block) is overwritten: the file is truncated at the raw block's end.
    subroutine append_deconv_block( fname, nptcls, ncomp, z, precision, labels, noise_scale )
        character(len=*),   intent(in) :: fname
        integer,            intent(in) :: nptcls, ncomp
        real(dp),           intent(in) :: z(:,:), precision(:,:,:)
        integer,  optional, intent(in) :: labels(:)
        real(dp), optional, intent(in) :: noise_scale
        integer, allocatable :: lab(:)
        real(dp) :: ns
        integer  :: u, io
        integer(kind=8) :: pos
        allocate(lab(nptcls), source=0)
        if( present(labels) ) lab(1:nptcls) = labels(1:nptcls)
        ns = 1.d0; if( present(noise_scale) ) ns = noise_scale
        pos = raw_block_end(nptcls, ncomp)
        open(newunit=u, file=fname, status='old', action='readwrite', access='stream', form='unformatted', iostat=io)
        if( io /= 0 ) THROW_HARD('append_deconv_block: cannot open '//trim(fname))
        write(u, pos=pos) COV_DECONV_MAGIC
        write(u) nptcls, ncomp
        write(u) z(1:nptcls,1:ncomp)
        write(u) precision(1:ncomp,1:ncomp,1:nptcls)
        write(u) lab
        write(u) ns
        endfile(u)
        close(u)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA deconvolved coordinates, precisions and labels appended to ', trim(fname)
        call flush(logfhandle)
        deallocate(lab)
    end subroutine append_deconv_block

    !> Read the deconvolved block of a cache if present and consistent with (nptcls, ncomp); found=.false.
    !! otherwise (older caches end after the raw block).
    subroutine read_deconv_block( fname, nptcls, ncomp, z, precision, labels, noise_scale, found )
        character(len=*),      intent(in)  :: fname
        integer,               intent(in)  :: nptcls, ncomp
        real(dp), allocatable, intent(out) :: z(:,:), precision(:,:,:)
        integer,  allocatable, intent(out) :: labels(:)
        real(dp),              intent(out) :: noise_scale
        logical,               intent(out) :: found
        character(len=len(COV_DECONV_MAGIC)) :: magic
        character(len=len(COV_CACHE_MAGIC))  :: magic0
        integer :: u, io, ver, n_c, k_c, n_b, k_b
        integer(kind=8) :: pos
        found = .false.; noise_scale = 1.d0
        if( .not. file_exists(fname) ) return
        open(newunit=u, file=fname, status='old', action='read', access='stream', form='unformatted', iostat=io)
        if( io /= 0 ) return
        read(u, iostat=io) magic0, ver
        if( io /= 0 .or. magic0 /= COV_CACHE_MAGIC .or. ver /= COV_CACHE_VERSION )then; close(u); return; endif
        read(u, iostat=io) n_c, k_c
        if( io /= 0 .or. n_c /= nptcls .or. k_c /= ncomp )then; close(u); return; endif
        pos = raw_block_end(nptcls, ncomp)
        read(u, pos=pos, iostat=io) magic
        if( io /= 0 .or. magic /= COV_DECONV_MAGIC )then; close(u); return; endif
        read(u, iostat=io) n_b, k_b
        if( io /= 0 .or. n_b /= nptcls .or. k_b /= ncomp )then; close(u); return; endif
        allocate(z(nptcls,ncomp), precision(ncomp,ncomp,nptcls), labels(nptcls))
        read(u, iostat=io) z
        if( io == 0 ) read(u, iostat=io) precision
        if( io == 0 ) read(u, iostat=io) labels
        if( io == 0 ) read(u, iostat=io) noise_scale
        close(u)
        if( io /= 0 )then
            deallocate(z, precision, labels); return
        endif
        found = .true.
    end subroutine read_deconv_block

end module simple_flex_pca_embedding_io

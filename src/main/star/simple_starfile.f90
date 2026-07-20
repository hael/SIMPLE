!@descr: STAR file I/O — writes optics, micrograph, and 2D-particle tables for RELION-compatible STAR files.
!        Tables are written to a staging file (<fname>.tmp) and atomically renamed on complete().
!        Parallel writers partition rows across OpenMP threads using an unequal-batch schedule; each
!        thread builds its sub-table independently and flushes in file order via !$omp ordered.
!        String fields are read inside !$omp critical regions to protect non-thread-safe oris accessors.
module simple_starfile
    use simple_core_module_api
    use simple_starfile_wrappers
    implicit none

    public :: starfile
    private
#include "simple_local_flags.inc"

    ! -------------------------------------------------------------------------
    ! starfile: atomic writer for RELION/STAR-format metadata files.
    !   fname   – final output path
    !   ftmp    – staging path (fname // '.tmp'); atomically renamed on complete()
    !   relroot – basename of the working-directory stem; used to relativise
    !             absolute file paths stored in oris string fields
    !   verbose – enable per-table wall-clock timing diagnostics
    ! -------------------------------------------------------------------------
    type starfile
        type(string) :: fname
        type(string) :: ftmp
        type(string) :: relroot
        integer      :: optics_offset = 0
        logical      :: verbose       = .false.
      contains
        procedure :: init
        procedure :: complete
        procedure :: write_optics_table
        procedure :: write_mics_table
        procedure :: write_ptcl2D_table
        procedure :: write_ptcl2D_table_parallel
    end type starfile

contains

    ! -------------------------------------------------------------------------
    ! init: initialise the starfile object.
    !   fname   – desired output file name
    !   verbose – (optional) enable timing diagnostics; default .false.
    ! Removes any stale .tmp file left by a previous run.
    ! -------------------------------------------------------------------------
    subroutine init( self, fname, verbose, optics_offset )
        class(starfile),   intent(inout) :: self
        class(string),     intent(in)    :: fname
        logical, optional, intent(in)    :: verbose
        integer, optional, intent(in)    :: optics_offset
        type(string) :: cwd
        call simple_getcwd(cwd)
        if( present(verbose)       ) self%verbose = verbose
        if( present(optics_offset) ) self%optics_offset = optics_offset
        self%fname   = fname
        self%ftmp    = fname%to_char() // '.tmp'
        self%relroot = basename(stemname(cwd))
        if( file_exists(self%ftmp) ) call del_file(self%ftmp)
    end subroutine init

    ! -------------------------------------------------------------------------
    ! complete: atomically promote the staging file to the final output path.
    ! No-op if the staging file does not exist (e.g. nothing was written).
    ! -------------------------------------------------------------------------
    subroutine complete( self )
        class(starfile), intent(inout) :: self
        if( file_exists(self%ftmp) ) then
            if( file_exists(self%fname) ) call del_file(self%fname)
            call simple_rename(self%ftmp, self%fname)
        end if
    end subroutine complete

    ! -------------------------------------------------------------------------
    ! calc_part_boundaries: divide n rows into nparts = omp_get_max_threads()
    ! non-overlapping contiguous sub-ranges stored in part_boundaries(:,1:2).
    ! Uses an unequal-batch schedule: later partitions receive slightly more
    ! rows so that total work is approximately load-balanced for operations
    ! whose per-row overhead grows with iteration index.
    ! When nparts <= 1 the single partition covers [1, n].
    ! -------------------------------------------------------------------------
    subroutine calc_part_boundaries( part_boundaries, n )
        integer, allocatable, intent(inout) :: part_boundaries(:,:)
        integer,              intent(in)    :: n
        integer :: nparts, batchdelta, batchbase, nstart, nend, part
        nparts = omp_get_max_threads()
        allocate( part_boundaries(nparts, 2) )
        if( nparts <= 1 ) then
            ! Single thread: one partition covers everything
            batchdelta = 0
            batchbase  = n
        else
            ! Unequal-batch: partition size grows linearly with part index
            batchdelta = ceiling( real(n) / ( real(nparts) * real(nparts - 1) ) )
            batchbase  = ceiling( real(n) / ( real(nparts) * 2 ) )
        end if
        nstart = 1
        nend   = 0
        do part = 1, nparts
            nend = nend + batchbase + (part - 1) * batchdelta
            if( nend > n ) nend = n
            part_boundaries(part, 1) = nstart
            part_boundaries(part, 2) = nend
            nstart = nend + 1
        end do
    end subroutine calc_part_boundaries

    ! -------------------------------------------------------------------------
    ! write_optics_table: write the 'optics' STAR block to ftmp.
    ! Rows with state == 0 are skipped.  Written in parallel ordered partitions;
    ! an empty table is emitted immediately when all rows are inactive.
    ! -------------------------------------------------------------------------
    subroutine write_optics_table( self, optics_oris )
        class(starfile),           intent(inout) :: self
        class(oris),               intent(in)    :: optics_oris
        type(starfile_table_type)                :: part_table
        integer, allocatable                     :: part_boundaries(:,:)
        integer                                  :: ipart, igrp, nobj
        integer(timer_int_kind)                  :: ms0 = 0_timer_int_kind
        real(timer_int_kind)                     :: ms_complete
        character(len=XLONGSTRLEN)               :: str_ogname
        logical                                  :: header_written
        if( self%verbose ) ms0 = tic()
        ! --- Fast path: empty input — emit an empty table and return ----------
        if( optics_oris%get_noris() == 0 ) then
            call starfile_table__new(part_table)
            call starfile_table__setIsList(part_table, .false.)
            call starfile_table__setname(part_table, 'optics')
            call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
            call starfile_table__write_ofile(part_table, tableend=.true.)
            call starfile_table__close_ofile(part_table)
            call starfile_table__delete(part_table)
            if( allocated(part_boundaries) ) deallocate(part_boundaries)
            if( self%verbose ) then
                ms_complete = toc(ms0)
                print *, 'wrote optics table in: ', ms_complete
                call flush(6)
            end if
            return
        end if
        call calc_part_boundaries( part_boundaries, optics_oris%get_noris() )
        ! --- Parallel ordered write -------------------------------------------
        header_written = .false.
        !$omp parallel default(shared) private(ipart, igrp, part_table, str_ogname, nobj)
        !$omp do ordered
        do ipart = 1, size(part_boundaries, 1)
            if( part_boundaries(ipart,1) > part_boundaries(ipart,2) ) cycle
            call starfile_table__new(part_table)
            call starfile_table__setIsList(part_table, .false.)
            call starfile_table__setname(part_table, 'optics')
            nobj = 0
            do igrp = part_boundaries(ipart,1), part_boundaries(ipart,2)
                if( optics_oris%get_state(igrp) == 0 ) cycle
                call starfile_table__addObject(part_table)
                nobj = nobj + 1
                ! ints
                if( optics_oris%isthere(igrp, 'ogid')  ) call starfile_table__setValue_int(part_table, EMDL_IMAGE_OPTICS_GROUP, optics_oris%get_int(igrp, 'ogid') + self%optics_offset)
                if( optics_oris%isthere(igrp, 'pop')   ) call starfile_table__setValue_int(part_table, SMPL_OPTICS_POPULATION,  optics_oris%get_int(igrp, 'pop' ))
                ! doubles
                if( optics_oris%isthere(igrp, 'kv')    ) call starfile_table__setValue_double(part_table, EMDL_CTF_VOLTAGE,      real(optics_oris%get(igrp, 'kv'),    dp))
                if( optics_oris%isthere(igrp, 'smpd')  ) call starfile_table__setValue_double(part_table, EMDL_IMAGE_PIXEL_SIZE, real(optics_oris%get(igrp, 'smpd'),  dp))
                if( optics_oris%isthere(igrp, 'cs')    ) call starfile_table__setValue_double(part_table, EMDL_CTF_CS,           real(optics_oris%get(igrp, 'cs'),    dp))
                if( optics_oris%isthere(igrp, 'fraca') ) call starfile_table__setValue_double(part_table, EMDL_CTF_Q0,           real(optics_oris%get(igrp, 'fraca'), dp))
                if( optics_oris%isthere(igrp, 'opcx')  ) call starfile_table__setValue_double(part_table, SMPL_OPTICS_CENTROIDX, real(optics_oris%get(igrp, 'opcx'),  dp))
                if( optics_oris%isthere(igrp, 'opcy')  ) call starfile_table__setValue_double(part_table, SMPL_OPTICS_CENTROIDY, real(optics_oris%get(igrp, 'opcy'),  dp))
                ! strings — get_static is not thread-safe; guard with critical
                if( optics_oris%isthere(igrp, 'ogname') ) then
                    !$omp critical
                    call optics_oris%get_static(igrp, 'ogname', str_ogname)
                    !$omp end critical
                    call starfile_table__setValue_string(part_table, EMDL_IMAGE_OPTICS_GROUP_NAME, trim(str_ogname))
                end if
            end do
            ! Flush this partition in file order; skip if no objects were added
            !$omp ordered
            if( nobj > 0 ) then
                call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
                if( .not. header_written ) then
                    call starfile_table__write_ofile(part_table, tableend=.false.)
                    header_written = .true.
                else
                    call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.false.)
                end if
                call starfile_table__close_ofile(part_table)
            end if
            call starfile_table__delete(part_table)
            !$omp end ordered
        end do
        !$omp end do
        !$omp end parallel
        ! Terminate the table, or emit an empty one when all rows were filtered out
        call starfile_table__new(part_table)
        call starfile_table__setIsList(part_table, .false.)
        call starfile_table__setname(part_table, 'optics')
        call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
        if( header_written ) then
            call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.true.)
        else
            call starfile_table__write_ofile(part_table, tableend=.true.)
        end if
        call starfile_table__close_ofile(part_table)
        call starfile_table__delete(part_table)
        if( allocated(part_boundaries) ) deallocate(part_boundaries)
        if( self%verbose ) then
            ms_complete = toc(ms0)
            print *, 'wrote optics table in: ', ms_complete
            call flush(6)
        end if
    end subroutine write_optics_table

    ! -------------------------------------------------------------------------
    ! write_mics_table: write the 'micrographs' STAR block to ftmp.
    ! Micrographs with state == 0 are skipped.  Scalar oris reads are
    ! unsynchronised; only get_static string reads are guarded with !$omp critical.
    ! -------------------------------------------------------------------------
    subroutine write_mics_table( self, mics_oris )
        class(starfile),           intent(inout) :: self
        class(oris),               intent(in)    :: mics_oris
        type(starfile_table_type)                :: part_table
        integer, allocatable                     :: part_boundaries(:,:)
        integer                                  :: ipart, imic, pathtrim, nobj
        integer(timer_int_kind)                  :: ms0 = 0_timer_int_kind
        real(timer_int_kind)                     :: ms_complete
        character(len=XLONGSTRLEN)               :: str_movie, str_intg, str_mc_starfile, str_boxfile, str_ctfjpg
        logical                                  :: header_written
        if( self%verbose ) ms0 = tic()
        ! --- Fast path: empty input — emit an empty table and return ----------
        if( mics_oris%get_noris() == 0 ) then
            call starfile_table__new(part_table)
            call starfile_table__setIsList(part_table, .false.)
            call starfile_table__setname(part_table, 'micrographs')
            call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
            call starfile_table__write_ofile(part_table, tableend=.true.)
            call starfile_table__close_ofile(part_table)
            call starfile_table__delete(part_table)
            if( allocated(part_boundaries) ) deallocate(part_boundaries)
            if( self%verbose ) then
                ms_complete = toc(ms0)
                print *, 'wrote micrographs table in: ', ms_complete
                call flush(6)
            end if
            return
        end if
        call calc_part_boundaries( part_boundaries, mics_oris%get_noris() )
        ! --- Parallel ordered write -------------------------------------------
        header_written = .false.
        !$omp parallel do ordered &
        !$omp   private(ipart, imic, part_table, pathtrim, nobj, str_movie, str_intg, str_mc_starfile, str_boxfile, str_ctfjpg) &
        !$omp   default(shared) proc_bind(close)
        do ipart = 1, size(part_boundaries, 1)
            if( part_boundaries(ipart,1) > part_boundaries(ipart,2) ) cycle
            call starfile_table__new(part_table)
            call starfile_table__setIsList(part_table, .false.)
            call starfile_table__setname(part_table, 'micrographs')
            pathtrim = 0  ! reset relative-path cache for each partition
            nobj = 0
            do imic = part_boundaries(ipart,1), part_boundaries(ipart,2)
                if( mics_oris%get_state(imic) == 0 ) cycle
                call starfile_table__addObject(part_table)
                nobj = nobj + 1
                ! ints
                if( mics_oris%isthere(imic, 'ogid'   ) ) call starfile_table__setValue_int(part_table, EMDL_IMAGE_OPTICS_GROUP, mics_oris%get_int(imic, 'ogid'   ) + self%optics_offset)
                if( mics_oris%isthere(imic, 'xdim'   ) ) call starfile_table__setValue_int(part_table, EMDL_IMAGE_SIZE_X,       mics_oris%get_int(imic, 'xdim'   ))
                if( mics_oris%isthere(imic, 'ydim'   ) ) call starfile_table__setValue_int(part_table, EMDL_IMAGE_SIZE_Y,       mics_oris%get_int(imic, 'ydim'   ))
                if( mics_oris%isthere(imic, 'nframes') ) call starfile_table__setValue_int(part_table, EMDL_IMAGE_SIZE_Z,       mics_oris%get_int(imic, 'nframes'))
                if( mics_oris%isthere(imic, 'nptcls' ) ) call starfile_table__setValue_int(part_table, SMPL_N_PTCLS,            mics_oris%get_int(imic, 'nptcls' ))
                if( mics_oris%isthere(imic, 'nmics'  ) ) call starfile_table__setValue_int(part_table, SMPL_N_MICS,             mics_oris%get_int(imic, 'nmics'  ))
                if( mics_oris%isthere(imic, 'micid'  ) ) call starfile_table__setValue_int(part_table, SMPL_MIC_ID,             mics_oris%get_int(imic, 'micid'  ))
                ! doubles (defocus stored in microns internally; STAR expects Angstroms → divide by 1e-4)
                if( mics_oris%isthere(imic, 'dfx'    ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_DEFOCUSU,      real(mics_oris%get(imic, 'dfx'    ) / 0.0001, dp))
                if( mics_oris%isthere(imic, 'dfy'    ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_DEFOCUSV,      real(mics_oris%get(imic, 'dfy'    ) / 0.0001, dp))
                if( mics_oris%isthere(imic, 'angast' ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_DEFOCUS_ANGLE, real(mics_oris%get(imic, 'angast' ),          dp))
                call starfile_table__setValue_double(part_table, EMDL_CTF_PHASESHIFT, &
                    &real(rad2deg(mics_oris%get(imic, 'phshift')), dp))
                if( mics_oris%isthere(imic, 'ctfres' ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_MAXRES,        real(mics_oris%get(imic, 'ctfres' ),          dp))
                if( mics_oris%isthere(imic, 'icefrac') ) call starfile_table__setValue_double(part_table, SMPL_ICE_FRAC,          real(mics_oris%get(imic, 'icefrac'),          dp))
                if( mics_oris%isthere(imic, 'astig'  ) ) call starfile_table__setValue_double(part_table, SMPL_ASTIGMATISM,       real(mics_oris%get(imic, 'astig'  ),          dp))
                ! strings — get_static is not thread-safe; guard with critical
                !$omp critical
                call mics_oris%get_static(imic, 'movie',       str_movie)
                call mics_oris%get_static(imic, 'intg',        str_intg)
                call mics_oris%get_static(imic, 'mc_starfile', str_mc_starfile)
                call mics_oris%get_static(imic, 'boxfile',     str_boxfile)
                call mics_oris%get_static(imic, 'ctfjpg',      str_ctfjpg)
                !$omp end critical
                str_intg        = get_relative_path_here(str_intg)
                str_mc_starfile = get_relative_path_here(str_mc_starfile)
                str_boxfile     = get_relative_path_here(str_boxfile)
                str_ctfjpg      = get_relative_path_here(str_ctfjpg)
                if( mics_oris%isthere(imic, 'movie'      ) ) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_MOVIE_NAME,    trim(str_movie))
                if( mics_oris%isthere(imic, 'intg'       ) ) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_NAME,          trim(str_intg))
                if( mics_oris%isthere(imic, 'mc_starfile') ) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_METADATA_NAME, trim(str_mc_starfile))
                if( mics_oris%isthere(imic, 'boxfile'    ) ) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_COORDINATES,   trim(str_boxfile))
                if( mics_oris%isthere(imic, 'ctfjpg'     ) ) call starfile_table__setValue_string(part_table, EMDL_CTF_PSPEC,                trim(str_ctfjpg))
            end do
            ! Flush this partition in file order; skip if no objects were added
            !$omp ordered
            if( nobj > 0 ) then
                call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
                if( .not. header_written ) then
                    call starfile_table__write_ofile(part_table, tableend=.false.)
                    header_written = .true.
                else
                    call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.false.)
                end if
                call starfile_table__close_ofile(part_table)
            end if
            call starfile_table__delete(part_table)
            !$omp end ordered
        end do
        !$omp end parallel do
        ! Terminate the table, or emit an empty one when all rows were filtered out
        call starfile_table__new(part_table)
        call starfile_table__setIsList(part_table, .false.)
        call starfile_table__setname(part_table, 'micrographs')
        call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
        if( header_written ) then
            call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.true.)
        else
            call starfile_table__write_ofile(part_table, tableend=.true.)
        end if
        call starfile_table__close_ofile(part_table)
        call starfile_table__delete(part_table)
        if( allocated(part_boundaries) ) deallocate(part_boundaries)
        if( self%verbose ) then
            ms_complete = toc(ms0)
            print *, 'wrote micrographs table in: ', ms_complete
            call flush(6)
        end if

        contains

            ! Trim path to a root-relative form using the cached pathtrim offset.
            ! On first call the offset of self%relroot in path is located and cached.
            function get_relative_path_here( path ) result( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if( pathtrim == 0 ) pathtrim = index(path, self%relroot%to_char())
                if( pathtrim > 0 ) then
                    newpath = trim(path(pathtrim:))
                else
                    newpath = trim(path)
                end if
            end function get_relative_path_here

    end subroutine write_mics_table

    ! -------------------------------------------------------------------------
    ! write_ptcl2D_table: write the 'particles' STAR block to ftmp (serial).
    ! Particles with state == 0 or invalid stkind are skipped.
    ! indstk is taken from ptcl2d_oris when present, otherwise computed from
    ! the stack's [fromp, top] range stored in stk_oris.
    ! When mics_oris is present and size matches stk_oris, the micrograph name
    ! is written for each particle.
    ! -------------------------------------------------------------------------
    subroutine write_ptcl2D_table( self, ptcl2d_oris, stk_oris, mics_oris )
        class(starfile),       intent(inout) :: self
        class(oris),           intent(in)    :: ptcl2d_oris, stk_oris
        class(oris), optional, intent(in)    :: mics_oris
        type(starfile_table_type)            :: ptcl_table
        integer                              :: iptcl, pathtrim, half_boxsize, stkind, indstk, fromp, top
        character(len=XLONGSTRLEN)           :: stkname, micname
        integer(timer_int_kind)              :: ms0 = 0_timer_int_kind
        real(timer_int_kind)                 :: ms_complete
        if( self%verbose ) ms0 = tic()
        call starfile_table__new(ptcl_table)
        call starfile_table__setIsList(ptcl_table, .false.)
        call starfile_table__setname(ptcl_table, 'particles')
        pathtrim = 0
        do iptcl = 1, ptcl2d_oris%get_noris()
            if( ptcl2d_oris%get_state(iptcl) == 0 ) cycle
            stkind = ptcl2d_oris%get_int(iptcl, 'stkind')
            if( stkind <= 0 ) cycle
            ! Resolve within-stack index
            if( ptcl2d_oris%isthere(iptcl, 'indstk') ) then
                indstk = ptcl2d_oris%get_int(iptcl, 'indstk')
            else
                fromp = stk_oris%get_fromp(stkind)
                top   = stk_oris%get_top(stkind)
                if( iptcl < fromp .or. iptcl > top ) cycle
                indstk = iptcl - fromp + 1
            end if
            call starfile_table__addObject(ptcl_table)
            half_boxsize = floor(stk_oris%get(stkind, 'box') / 2.0)
            ! ints
            if( ptcl2d_oris%isthere(iptcl, 'ogid'   ) ) call starfile_table__setValue_int(ptcl_table, EMDL_IMAGE_OPTICS_GROUP, ptcl2d_oris%get_int(iptcl, 'ogid') + self%optics_offset)
            if( ptcl2d_oris%isthere(iptcl, 'class'  ) ) call starfile_table__setValue_int(ptcl_table, EMDL_PARTICLE_CLASS,     ptcl2d_oris%get_class(iptcl))
            if( ptcl2d_oris%isthere(iptcl, 'gid'    ) ) call starfile_table__setValue_int(ptcl_table, EMDL_MLMODEL_GROUP_NO,   ptcl2d_oris%get_int(iptcl, 'gid'))
            ! doubles (defocus: microns → Angstroms; coords: pixel + half-box offset)
            if( ptcl2d_oris%isthere(iptcl, 'dfx'    ) ) call starfile_table__setValue_double(ptcl_table, EMDL_CTF_DEFOCUSU,             real(ptcl2d_oris%get(iptcl, 'dfx'    ) / 0.0001,        dp))
            if( ptcl2d_oris%isthere(iptcl, 'dfy'    ) ) call starfile_table__setValue_double(ptcl_table, EMDL_CTF_DEFOCUSV,             real(ptcl2d_oris%get(iptcl, 'dfy'    ) / 0.0001,        dp))
            if( ptcl2d_oris%isthere(iptcl, 'angast' ) ) call starfile_table__setValue_double(ptcl_table, EMDL_CTF_DEFOCUS_ANGLE,        real(ptcl2d_oris%get(iptcl, 'angast' ),                 dp))
            call starfile_table__setValue_double(ptcl_table, EMDL_CTF_PHASESHIFT, &
                &real(rad2deg(ptcl2d_oris%get(iptcl, 'phshift')), dp))
            if( ptcl2d_oris%isthere(iptcl, 'e3'     ) ) call starfile_table__setValue_double(ptcl_table, EMDL_ORIENT_PSI,               real(ptcl2d_oris%get(iptcl, 'e3'     ),                 dp))
            if( ptcl2d_oris%isthere(iptcl, 'xpos'   ) ) call starfile_table__setValue_double(ptcl_table, EMDL_IMAGE_COORD_X,            real(ptcl2d_oris%get(iptcl, 'xpos'   ) + half_boxsize,  dp))
            if( ptcl2d_oris%isthere(iptcl, 'ypos'   ) ) call starfile_table__setValue_double(ptcl_table, EMDL_IMAGE_COORD_Y,            real(ptcl2d_oris%get(iptcl, 'ypos'   ) + half_boxsize,  dp))
            if( ptcl2d_oris%isthere(iptcl, 'x'      ) ) call starfile_table__setValue_double(ptcl_table, EMDL_ORIENT_ORIGIN_X_ANGSTROM, real(ptcl2d_oris%get(iptcl, 'x'      ),                 dp))
            if( ptcl2d_oris%isthere(iptcl, 'y'      ) ) call starfile_table__setValue_double(ptcl_table, EMDL_ORIENT_ORIGIN_Y_ANGSTROM, real(ptcl2d_oris%get(iptcl, 'y'      ),                 dp))
            ! strings
            if( indstk > 0 ) then
                if( stk_oris%isthere(stkind, 'stk') ) then
                    call stk_oris%get_static(stkind, 'stk', stkname)
                    stkname = get_relative_path_here(stkname)
                    stkname = int2str(indstk) // '@' // trim(stkname)
                    call starfile_table__setValue_string(ptcl_table, EMDL_IMAGE_NAME, trim(stkname))
                    if( present(mics_oris) ) then
                        ! micrograph name available when stk 1:1 maps to mics
                        if( stk_oris%get_noris() == mics_oris%get_noris() ) then
                            if( mics_oris%isthere(stkind, 'intg') ) then
                                call mics_oris%get_static(stkind, 'intg', micname)
                                micname = get_relative_path_here(micname)
                                call starfile_table__setValue_string(ptcl_table, EMDL_MICROGRAPH_NAME, trim(micname))
                            end if
                        end if
                    end if
                end if
            end if
        end do
        call starfile_table__open_ofile(ptcl_table, self%ftmp%to_char(), 1)
        call starfile_table__write_ofile(ptcl_table, tableend=.true.)
        call starfile_table__close_ofile(ptcl_table)
        call starfile_table__delete(ptcl_table)
        if( self%verbose ) then
            ms_complete = toc(ms0)
            print *, 'wrote particles2D table in: ', ms_complete
            call flush(6)
        end if

        contains

            ! Trim path to a root-relative form using the cached pathtrim offset.
            ! On first call the offset of self%relroot in path is located and cached.
            function get_relative_path_here( path ) result( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if( pathtrim == 0 ) pathtrim = index(path, self%relroot%to_char())
                if( pathtrim > 0 ) then
                    newpath = trim(path(pathtrim:))
                else
                    newpath = trim(path)
                end if
            end function get_relative_path_here

    end subroutine write_ptcl2D_table

    ! -------------------------------------------------------------------------
    ! write_ptcl2D_table_parallel: parallel version of write_ptcl2D_table.
    ! Row range is split across omp_get_max_threads() partitions; each thread
    ! builds its own sub-table independently and flushes in file order via
    ! !$omp ordered.  get_static calls are protected by !$omp critical.
    ! -------------------------------------------------------------------------
    subroutine write_ptcl2D_table_parallel( self, ptcl2d_oris, stk_oris, mics_oris )
        class(starfile),       intent(inout) :: self
        class(oris),           intent(in)    :: ptcl2d_oris, stk_oris
        class(oris), optional, intent(in)    :: mics_oris
        type(starfile_table_type)            :: part_table
        integer, allocatable                 :: part_boundaries(:,:)
        character(len=XLONGSTRLEN)           :: micname, str_stk
        integer                              :: ipart, iptcl, pathtrim, half_boxsize, stkind, indstk, fromp, top, nobj
        integer(timer_int_kind)              :: ms0 = 0_timer_int_kind
        real(timer_int_kind)                 :: ms_complete
        logical                              :: header_written
        if( self%verbose ) ms0 = tic()
        ! --- Fast path: empty input — emit an empty table and return ----------
        if( ptcl2d_oris%get_noris() == 0 ) then
            call starfile_table__new(part_table)
            call starfile_table__setIsList(part_table, .false.)
            call starfile_table__setname(part_table, 'particles')
            call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
            call starfile_table__write_ofile(part_table, tableend=.true.)
            call starfile_table__close_ofile(part_table)
            call starfile_table__delete(part_table)
            if( allocated(part_boundaries) ) deallocate(part_boundaries)
            if( self%verbose ) then
                ms_complete = toc(ms0)
                print *, 'wrote particles2D table in: ', ms_complete
                call flush(6)
            end if
            return
        end if
        call calc_part_boundaries( part_boundaries, ptcl2d_oris%get_noris() )
        ! --- Parallel ordered write -------------------------------------------
        header_written = .false.
        !$omp parallel do ordered default(shared) &
        !$omp   private(ipart, iptcl, part_table, half_boxsize, stkind, indstk, fromp, top, pathtrim, nobj, str_stk, micname) &
        !$omp   proc_bind(close)
        do ipart = 1, size(part_boundaries, 1)
            if( part_boundaries(ipart,1) > part_boundaries(ipart,2) ) cycle
            call starfile_table__new(part_table)
            call starfile_table__setIsList(part_table, .false.)
            call starfile_table__setname(part_table, 'particles')
            pathtrim = 0  ! reset relative-path cache for each partition
            nobj = 0
            do iptcl = part_boundaries(ipart,1), part_boundaries(ipart,2)
                if( ptcl2d_oris%get_state(iptcl) == 0 ) cycle
                stkind = ptcl2d_oris%get_int(iptcl, 'stkind')
                if( stkind <= 0 ) cycle
                ! Resolve within-stack index
                if( ptcl2d_oris%isthere(iptcl, 'indstk') ) then
                    indstk = ptcl2d_oris%get_int(iptcl, 'indstk')
                else
                    fromp = stk_oris%get_fromp(stkind)
                    top   = stk_oris%get_top(stkind)
                    if( iptcl < fromp .or. iptcl > top ) cycle
                    indstk = iptcl - fromp + 1
                end if
                call starfile_table__addObject(part_table)
                nobj = nobj + 1
                half_boxsize = floor(stk_oris%get(stkind, 'box') / 2.0)
                ! ints
                if( ptcl2d_oris%isthere(iptcl, 'ogid'   ) ) call starfile_table__setValue_int(part_table, EMDL_IMAGE_OPTICS_GROUP, ptcl2d_oris%get_int(iptcl, 'ogid') + self%optics_offset)
                if( ptcl2d_oris%isthere(iptcl, 'class'  ) ) call starfile_table__setValue_int(part_table, EMDL_PARTICLE_CLASS,     ptcl2d_oris%get_class(iptcl))
                if( ptcl2d_oris%isthere(iptcl, 'gid'    ) ) call starfile_table__setValue_int(part_table, EMDL_MLMODEL_GROUP_NO,   ptcl2d_oris%get_int(iptcl, 'gid'))
                ! doubles (defocus: microns → Angstroms; coords: pixel + half-box offset)
                if( ptcl2d_oris%isthere(iptcl, 'dfx'    ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_DEFOCUSU,             real(ptcl2d_oris%get(iptcl, 'dfx'    ) / 0.0001,        dp))
                if( ptcl2d_oris%isthere(iptcl, 'dfy'    ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_DEFOCUSV,             real(ptcl2d_oris%get(iptcl, 'dfy'    ) / 0.0001,        dp))
                if( ptcl2d_oris%isthere(iptcl, 'angast' ) ) call starfile_table__setValue_double(part_table, EMDL_CTF_DEFOCUS_ANGLE,        real(ptcl2d_oris%get(iptcl, 'angast' ),                 dp))
                call starfile_table__setValue_double(part_table, EMDL_CTF_PHASESHIFT, &
                    &real(rad2deg(ptcl2d_oris%get(iptcl, 'phshift')), dp))
                if( ptcl2d_oris%isthere(iptcl, 'e3'     ) ) call starfile_table__setValue_double(part_table, EMDL_ORIENT_PSI,               real(ptcl2d_oris%get(iptcl, 'e3'     ),                 dp))
                if( ptcl2d_oris%isthere(iptcl, 'xpos'   ) ) call starfile_table__setValue_double(part_table, EMDL_IMAGE_COORD_X,            real(ptcl2d_oris%get(iptcl, 'xpos'   ) + half_boxsize,  dp))
                if( ptcl2d_oris%isthere(iptcl, 'ypos'   ) ) call starfile_table__setValue_double(part_table, EMDL_IMAGE_COORD_Y,            real(ptcl2d_oris%get(iptcl, 'ypos'   ) + half_boxsize,  dp))
                if( ptcl2d_oris%isthere(iptcl, 'x'      ) ) call starfile_table__setValue_double(part_table, EMDL_ORIENT_ORIGIN_X_ANGSTROM, real(ptcl2d_oris%get(iptcl, 'x'      ),                 dp))
                if( ptcl2d_oris%isthere(iptcl, 'y'      ) ) call starfile_table__setValue_double(part_table, EMDL_ORIENT_ORIGIN_Y_ANGSTROM, real(ptcl2d_oris%get(iptcl, 'y'      ),                 dp))
                ! strings — get_static is not thread-safe; guard with critical
                if( indstk > 0 ) then
                    if( stk_oris%isthere(stkind, 'stk') ) then
                        !$omp critical
                        call stk_oris%get_static(stkind, 'stk', str_stk)
                        str_stk = get_relative_path_here(str_stk)
                        !$omp end critical
                        call starfile_table__setValue_string(part_table, EMDL_IMAGE_NAME, int2str(indstk) // '@' // trim(str_stk))
                        if( present(mics_oris) ) then
                            ! micrograph name available when stk 1:1 maps to mics
                            if( stk_oris%get_noris() == mics_oris%get_noris() ) then
                                micname = ''
                                !$omp critical
                                if( mics_oris%isthere(stkind, 'intg') ) then
                                    call mics_oris%get_static(stkind, 'intg', micname)
                                end if
                                !$omp end critical
                                if( len_trim(micname) > 0 ) then
                                    micname = get_relative_path_here(micname)
                                    call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_NAME, trim(micname))
                                end if
                            end if
                        end if
                    end if
                end if
            end do
            ! Flush this partition in file order; skip if no objects were added
            !$omp ordered
            if( nobj > 0 ) then
                call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
                if( .not. header_written ) then
                    call starfile_table__write_ofile(part_table, tableend=.false.)
                    header_written = .true.
                else
                    call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.false.)
                end if
                call starfile_table__close_ofile(part_table)
            end if
            call starfile_table__delete(part_table)
            !$omp end ordered
        end do
        !$omp end parallel do
        ! Terminate the table, or emit an empty one when all rows were filtered out
        call starfile_table__new(part_table)
        call starfile_table__setIsList(part_table, .false.)
        call starfile_table__setname(part_table, 'particles')
        call starfile_table__open_ofile(part_table, self%ftmp%to_char(), 1)
        if( header_written ) then
            call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.true.)
        else
            call starfile_table__write_ofile(part_table, tableend=.true.)
        end if
        call starfile_table__close_ofile(part_table)
        call starfile_table__delete(part_table)
        if( allocated(part_boundaries) ) deallocate(part_boundaries)
        if( self%verbose ) then
            ms_complete = toc(ms0)
            print *, 'wrote particles2D table in: ', ms_complete
            call flush(6)
        end if

        contains

            ! Trim path to a root-relative form using the cached pathtrim offset.
            ! On first call the offset of self%relroot in path is located and cached.
            function get_relative_path_here( path ) result( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if( pathtrim == 0 ) pathtrim = index(path, self%relroot%to_char())
                if( pathtrim > 0 ) then
                    newpath = trim(path(pathtrim:))
                else
                    newpath = trim(path)
                end if
            end function get_relative_path_here

    end subroutine write_ptcl2D_table_parallel

end module simple_starfile

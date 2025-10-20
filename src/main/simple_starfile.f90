module simple_starfile
include 'simple_lib.f08'
!$ use omp_lib

use simple_starfile_wrappers
implicit none

public :: starfile
private
#include "simple_local_flags.inc"

type starfile

    type(starfile_table_type)     :: table
    character(len=LONGSTRLEN)     :: fname
    character(len=LONGSTRLEN)     :: ftmp
    character(len=:), allocatable :: relroot
    logical                       :: verbose = .false.

contains

    procedure :: init
    procedure :: complete
    procedure :: write_optics_table
    procedure :: write_mics_table
    procedure :: write_ptcl2D_table

end type starfile

contains

    subroutine init( self, fname, verbose)
        class(starfile),            intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        logical, optional,          intent(in)    :: verbose
        character(len=XLONGSTRLEN)                :: cwd
        call simple_getcwd(cwd)
        if(present(verbose)) self%verbose = verbose
        self%fname   = fname
        self%ftmp    = fname // '.tmp'
        self%relroot = basename(stemname(cwd))
        if(file_exists(trim(adjustl(self%ftmp)))) call del_file(trim(adjustl(self%ftmp)))
    end subroutine init

    subroutine complete( self )
        class(starfile),  intent(inout) :: self
        if(file_exists(trim(adjustl(self%ftmp)))) then
            if(file_exists(trim(adjustl(self%fname)))) call del_file(trim(adjustl(self%fname)))
            call simple_rename(trim(adjustl(self%ftmp)), trim(adjustl(self%fname)))
        end if
    end subroutine complete

    subroutine calc_part_boundaries( part_boundaries, n )
        integer, allocatable, intent(inout) :: part_boundaries(:,:)
        integer,              intent(in)    :: n
        integer                             :: nparts, batchdelta, batchbase, nstart, nend, part
        !$omp parallel
            !$omp master
                nparts = omp_get_num_threads()
            !$omp end master
        !$omp end parallel
        allocate(part_boundaries(nparts, 2))
        batchdelta = ceiling(real(n) / (real(nparts) * real(nparts - 1)))
        batchbase  = ceiling(real(n) / (real(nparts) * 2))
        nstart = 1
        nend   = 0
        do part=1, nparts
            nend = nend + batchbase + ((part - 1) * batchdelta)
            if(nend .gt. n) nend = n
            part_boundaries(part,1) = nstart
            part_boundaries(part,2) = nend
            nstart = nend + 1
        end do
    end subroutine calc_part_boundaries

    subroutine write_optics_table( self, optics_oris )
        class(starfile),            intent(inout) :: self
        class(oris),                intent(in)    :: optics_oris
        type(starfile_table_type)                 :: part_table
        integer,                    allocatable   :: part_boundaries(:,:)
        integer                                   :: ipart, igrp
        integer(timer_int_kind)                   :: ms0
        real(timer_int_kind)                      :: ms_complete
        ms0 = tic()
        call calc_part_boundaries(part_boundaries, optics_oris%get_noris())
        !$omp parallel private(ipart, igrp, part_table)
            !$omp do ordered
                do ipart=1, size(part_boundaries, 1)
                    call starfile_table__new(part_table)
                    call starfile_table__setIsList(part_table, .false.)
                    call starfile_table__setname(part_table, 'optics')
                    if(part_boundaries(ipart,1) .gt. part_boundaries(ipart,2)) cycle
                    do igrp=part_boundaries(ipart,1), part_boundaries(ipart,2)
                        if(optics_oris%get(igrp, 'state') .eq. 0.0 ) cycle
                        call starfile_table__addObject(part_table)
                        ! ints
                        if(optics_oris%isthere(igrp, 'ogid'))   call starfile_table__setValue_int(part_table, EMDL_IMAGE_OPTICS_GROUP, optics_oris%get_int(igrp, 'ogid'))
                        if(optics_oris%isthere(igrp, 'pop'))    call starfile_table__setValue_int(part_table, SMPL_OPTICS_POPULATION,  optics_oris%get_int(igrp, 'pop' ))
                        ! doubles
                        if(optics_oris%isthere(igrp, 'kv'))     call starfile_table__setValue_double(part_table, EMDL_CTF_VOLTAGE,      real(optics_oris%get(igrp, 'kv'),    dp))
                        if(optics_oris%isthere(igrp, 'smpd'))   call starfile_table__setValue_double(part_table, EMDL_IMAGE_PIXEL_SIZE, real(optics_oris%get(igrp, 'smpd'),  dp))
                        if(optics_oris%isthere(igrp, 'cs'))     call starfile_table__setValue_double(part_table, EMDL_CTF_CS,           real(optics_oris%get(igrp, 'cs'),    dp))
                        if(optics_oris%isthere(igrp, 'fraca'))  call starfile_table__setValue_double(part_table, EMDL_CTF_Q0,           real(optics_oris%get(igrp, 'fraca'), dp))
                        if(optics_oris%isthere(igrp, 'opcx'))   call starfile_table__setValue_double(part_table, SMPL_OPTICS_CENTROIDX, real(optics_oris%get(igrp, 'opcx'),  dp))
                        if(optics_oris%isthere(igrp, 'opcy'))   call starfile_table__setValue_double(part_table, SMPL_OPTICS_CENTROIDY, real(optics_oris%get(igrp, 'opcy'),  dp))
                        ! strings
                        if(optics_oris%isthere(igrp, 'ogname')) call starfile_table__setValue_string(part_table, EMDL_IMAGE_OPTICS_GROUP_NAME, trim(adjustl(optics_oris%get_static(igrp, 'ogname'))))
                    end do
                    !$omp ordered
                        call starfile_table__open_ofile(part_table, trim(self%ftmp), 1)
                        if(ipart == 1) then
                            call starfile_table__write_ofile(part_table, tableend=.false.)
                        else
                            call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.false.)
                        end if
                        call starfile_table__close_ofile(part_table)
                    !$omp end ordered
                end do
            !$omp end do
        !$omp end parallel
        if(allocated(part_boundaries)) deallocate(part_boundaries)
        ms_complete = toc(ms0)
        if (self%verbose) print *,'wrote optics table in: ', ms_complete; call flush(6)
    end subroutine write_optics_table

    subroutine write_mics_table( self, mics_oris )
        class(starfile),            intent(inout) :: self
        class(oris),                intent(in)    :: mics_oris
        type(starfile_table_type)                 :: part_table
        integer,                    allocatable   :: part_boundaries(:,:)
        integer                                   :: ipart, imic, pathtrim
        integer(timer_int_kind)                   :: ms0
        real(timer_int_kind)                      :: ms_complete
        ms0 = tic()
        call calc_part_boundaries(part_boundaries, mics_oris%get_noris())
        !$omp parallel private(ipart, imic, part_table, ms0, ms_complete, pathtrim)
            !$omp do ordered
                do ipart=1, size(part_boundaries, 1)
                    call starfile_table__new(part_table)
                    call starfile_table__setIsList(part_table, .false.)
                    call starfile_table__setname(part_table, 'micrographs')
                    if(part_boundaries(ipart,1) .gt. part_boundaries(ipart,2)) cycle
                    pathtrim = 0
                    do imic=part_boundaries(ipart,1), part_boundaries(ipart,2)
                        if(mics_oris%get_state(imic) .eq. 0 ) cycle
                        call starfile_table__addObject(part_table)
                        ! ints
                        if(mics_oris%isthere(imic, 'ogid'   )) call starfile_table__setValue_int(part_table, EMDL_IMAGE_OPTICS_GROUP, mics_oris%get_int(imic, 'ogid'   ))
                        if(mics_oris%isthere(imic, 'xdim'   )) call starfile_table__setValue_int(part_table, EMDL_IMAGE_SIZE_X,       mics_oris%get_int(imic, 'xdim'   ))
                        if(mics_oris%isthere(imic, 'ydim'   )) call starfile_table__setValue_int(part_table, EMDL_IMAGE_SIZE_Y,       mics_oris%get_int(imic, 'ydim'   ))
                        if(mics_oris%isthere(imic, 'nframes')) call starfile_table__setValue_int(part_table, EMDL_IMAGE_SIZE_Z,       mics_oris%get_int(imic, 'nframes'))
                        if(mics_oris%isthere(imic, 'nptcls' )) call starfile_table__setValue_int(part_table, SMPL_N_PTCLS,            mics_oris%get_int(imic, 'nptcls' ))
                        if(mics_oris%isthere(imic, 'nmics'  )) call starfile_table__setValue_int(part_table, SMPL_N_MICS,             mics_oris%get_int(imic, 'nmics'  ))
                        if(mics_oris%isthere(imic, 'micid'  )) call starfile_table__setValue_int(part_table, SMPL_MIC_ID,             mics_oris%get_int(imic, 'micid'  ))
                        ! doubles
                        if(mics_oris%isthere(imic, 'dfx'    )) call starfile_table__setValue_double(part_table,  EMDL_CTF_DEFOCUSU,      real(mics_oris%get(imic, 'dfx') / 0.0001, dp))
                        if(mics_oris%isthere(imic, 'dfy'    )) call starfile_table__setValue_double(part_table,  EMDL_CTF_DEFOCUSV,      real(mics_oris%get(imic, 'dfy') / 0.0001, dp))
                        if(mics_oris%isthere(imic, 'angast' )) call starfile_table__setValue_double(part_table,  EMDL_CTF_DEFOCUS_ANGLE, real(mics_oris%get(imic, 'angast'),       dp))
                        if(mics_oris%isthere(imic, 'phshift')) call starfile_table__setValue_double(part_table,  EMDL_CTF_PHASESHIFT,    real(mics_oris%get(imic, 'phshift'),      dp))
                        if(mics_oris%isthere(imic, 'ctfres' )) call starfile_table__setValue_double(part_table,  EMDL_CTF_MAXRES,        real(mics_oris%get(imic, 'ctfres'),       dp))
                        if(mics_oris%isthere(imic, 'icefrac')) call starfile_table__setValue_double(part_table,  SMPL_ICE_FRAC,          real(mics_oris%get(imic, 'icefrac'),      dp))
                        if(mics_oris%isthere(imic, 'astig'  )) call starfile_table__setValue_double(part_table,  SMPL_ASTIGMATISM,       real(mics_oris%get(imic, 'astig'),        dp))
                        ! strings
                        if(mics_oris%isthere(imic, 'movie'      )) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_MOVIE_NAME,    trim(mics_oris%get_static(imic, 'movie')))
                        if(mics_oris%isthere(imic, 'intg'       )) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_NAME,          trim(get_relative_path(trim(mics_oris%get_static(imic, 'intg')),        self%relroot, pathtrim)))
                        if(mics_oris%isthere(imic, 'mc_starfile')) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_METADATA_NAME, trim(get_relative_path(trim(mics_oris%get_static(imic, 'mc_starfile')), self%relroot, pathtrim)))
                        if(mics_oris%isthere(imic, 'boxfile'    )) call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_COORDINATES,   trim(get_relative_path(trim(mics_oris%get_static(imic, 'boxfile'    )), self%relroot, pathtrim)))
                        if(mics_oris%isthere(imic, 'ctfjpg'     )) call starfile_table__setValue_string(part_table, EMDL_CTF_PSPEC,                trim(get_relative_path(trim(mics_oris%get_static(imic, 'ctfjpg'     )), self%relroot, pathtrim)))
                    end do
                    !$omp ordered
                        call starfile_table__open_ofile(part_table, trim(self%ftmp), 1)
                        if(ipart == 1) then
                            call starfile_table__write_ofile(part_table, tableend=.false.)
                        else
                            call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.false.)
                        end if
                        call starfile_table__close_ofile(part_table)
                    !$omp end ordered
                end do
            !$omp end do
        !$omp end parallel
        if(allocated(part_boundaries)) deallocate(part_boundaries)
        ms_complete = toc(ms0)
        if (self%verbose) print *,'wrote micrographs table in: ', ms_complete; call flush(6)
    end subroutine write_mics_table

    subroutine write_ptcl2D_table(self, ptcl2d_oris, stk_oris, mics_oris)
        class(starfile),            intent(inout) :: self
        class(oris),                intent(in)    :: ptcl2d_oris, stk_oris
        class(oris), optional,      intent(in)    :: mics_oris
        type(starfile_table_type)                 :: part_table
        integer,                    allocatable   :: part_boundaries(:,:)
        character(len=XLONGSTRLEN)                :: stkname, micname
        integer                                   :: ipart, iptcl, pathtrim, half_boxsize, stkind, indstk, fromp, top
        integer(timer_int_kind)                   :: ms0
        real(timer_int_kind)                      :: ms_complete
        ms0 = tic()
        call calc_part_boundaries(part_boundaries, ptcl2d_oris%get_noris())
        !$omp parallel private(ipart, iptcl, part_table, half_boxsize, stkind, indstk, fromp, top, ms0, ms_complete, pathtrim, stkname, micname) 
            !$omp do ordered
                do ipart=1, size(part_boundaries, 1)
                    call starfile_table__new(part_table)
                    call starfile_table__setIsList(part_table, .false.)
                    call starfile_table__setname(part_table, 'particles')
                    if(part_boundaries(ipart,1) .gt. part_boundaries(ipart,2)) cycle
                    pathtrim = 0
                    do iptcl=part_boundaries(ipart,1), part_boundaries(ipart,2)
                        if(ptcl2d_oris%get_state(iptcl) .eq. 0 ) cycle
                        call starfile_table__addObject(part_table)
                    !     call spproj%get_stkname_and_ind('ptcl2D', i, stkname, ind_in_stk)
                        stkind = ptcl2d_oris%get_int(iptcl, 'stkind')
                        if( ptcl2d_oris%isthere(iptcl, 'indstk') ) then
                            indstk = ptcl2d_oris%get_int(iptcl, 'indstk')
                        else
                            fromp = stk_oris%get_fromp(stkind)
                            top   = stk_oris%get_top(stkind)
                            if( iptcl < fromp .or. iptcl > top ) cycle
                            indstk = iptcl - fromp + 1
                        endif
                        half_boxsize = floor(stk_oris%get(stkind, 'box') / 2.0)
                        ! ints
                        if(ptcl2d_oris%isthere(iptcl, 'ogid'   )) call starfile_table__setValue_int(part_table, EMDL_IMAGE_OPTICS_GROUP, ptcl2d_oris%get_int(iptcl, 'ogid'))
                        if(ptcl2d_oris%isthere(iptcl, 'class'  )) call starfile_table__setValue_int(part_table, EMDL_PARTICLE_CLASS,     ptcl2d_oris%get_class(iptcl))
                        if(ptcl2d_oris%isthere(iptcl, 'gid'    )) call starfile_table__setValue_int(part_table, EMDL_MLMODEL_GROUP_NO,   ptcl2d_oris%get_int(iptcl, 'gid'))
                        ! doubles
                        if(ptcl2d_oris%isthere(iptcl, 'dfx'    )) call starfile_table__setValue_double(part_table,  EMDL_CTF_DEFOCUSU,              real(ptcl2d_oris%get(iptcl, 'dfx') / 0.0001,        dp))
                        if(ptcl2d_oris%isthere(iptcl, 'dfy'    )) call starfile_table__setValue_double(part_table,  EMDL_CTF_DEFOCUSV,              real(ptcl2d_oris%get(iptcl, 'dfy') / 0.0001,        dp))
                        if(ptcl2d_oris%isthere(iptcl, 'angast' )) call starfile_table__setValue_double(part_table,  EMDL_CTF_DEFOCUS_ANGLE,         real(ptcl2d_oris%get(iptcl, 'angast'),              dp))
                        if(ptcl2d_oris%isthere(iptcl, 'phshift')) call starfile_table__setValue_double(part_table,  EMDL_CTF_PHASESHIFT,            real(ptcl2d_oris%get(iptcl, 'phshift'),             dp))
                        if(ptcl2d_oris%isthere(iptcl, 'e3'     )) call starfile_table__setValue_double(part_table,  EMDL_ORIENT_PSI,                real(ptcl2d_oris%get(iptcl, 'e3'),                  dp))
                        if(ptcl2d_oris%isthere(iptcl, 'xpos'   )) call starfile_table__setValue_double(part_table,  EMDL_IMAGE_COORD_X,             real(ptcl2d_oris%get(iptcl, 'xpos') + half_boxsize, dp))
                        if(ptcl2d_oris%isthere(iptcl, 'ypos'   )) call starfile_table__setValue_double(part_table,  EMDL_IMAGE_COORD_Y,             real(ptcl2d_oris%get(iptcl, 'ypos') + half_boxsize, dp))
                        if(ptcl2d_oris%isthere(iptcl, 'x'      )) call starfile_table__setValue_double(part_table,  EMDL_ORIENT_ORIGIN_X_ANGSTROM,  real(ptcl2d_oris%get(iptcl, 'x'),                   dp))
                        if(ptcl2d_oris%isthere(iptcl, 'y'      )) call starfile_table__setValue_double(part_table,  EMDL_ORIENT_ORIGIN_Y_ANGSTROM,  real(ptcl2d_oris%get(iptcl, 'y'),                   dp))
                        ! strings
                        if(stkind .gt. 0 .and. indstk .gt. 0) then
                            if(stk_oris%isthere(stkind, 'stk')) then
                                !$omp critical
                                stkname = int2str(indstk) // '@' // trim(get_relative_path(trim(stk_oris%get_static(stkind, 'stk')), self%relroot, pathtrim))
                                !$omp end critical
                                call starfile_table__setValue_string(part_table, EMDL_IMAGE_NAME, trim(stkname))
                                if(present(mics_oris)) then
                                    if(stk_oris%get_noris() .eq. mics_oris%get_noris()) then
                                        if(mics_oris%isthere(stkind, 'intg')) then
                                            micname = trim(get_relative_path(trim(mics_oris%get_static(stkind, 'intg')), self%relroot, pathtrim))
                                            call starfile_table__setValue_string(part_table, EMDL_MICROGRAPH_NAME, trim(micname))
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end do
                    !$omp ordered
                        call starfile_table__open_ofile(part_table, trim(self%ftmp), 1)
                        if(ipart == 1) then
                            call starfile_table__write_ofile(part_table, tableend=.false.)
                        else
                            call starfile_table__write_ofile(part_table, ignoreheader=.true., tableend=.false.)
                        end if
                        call starfile_table__close_ofile(part_table)
                    !$omp end ordered
                end do
            !$omp end do
        !$omp end parallel
        if(allocated(part_boundaries)) deallocate(part_boundaries)
        ms_complete = toc(ms0)
        if (self%verbose) print *,'wrote particles2D table in: ', ms_complete; call flush(6)
    end subroutine write_ptcl2D_table

end module simple_starfile
!@descr: for interoperability with RELION
module simple_relion
use, intrinsic :: iso_c_binding, only: C_int, C_double, C_NULL_CHAR
use simple_defs,       only: dp, logfhandle, XLONGSTRLEN, C_TRUE, C_FALSE
use simple_defs_fname, only: INTGMOV_SUFFIX, FORCTF_SUFFIX, MRC_EXT, STK_EXT, STAR_EXT, EXTRACT_STK_FBODY
use simple_error,      only: simple_exception
use simple_fileio,     only: basename, file_exists, fname_new_ext, simple_abspath, simple_rename
use simple_linalg,     only: rad2deg
use simple_oris,       only: oris
use simple_rnd,        only: seed_rnd, ran3
use simple_string,     only: string
use simple_string_utils, only: int2str
use simple_sp_project, only: sp_project
use simple_cmdline,    only: cmdline
use simple_syslib,     only: simple_rmfile, simple_mkdir, syslib_symlink
use simple_starfile_wrappers, only: starfile_table_type, starfile_table__new, starfile_table__delete, &
    starfile_table__addObject, starfile_table__setName, starfile_table__setValue_int, &
    starfile_table__setValue_double, starfile_table__setValue_string, starfile_table__open_ofile, &
    starfile_table__write_ofile, starfile_table__close_ofile, &
    EMDL_CTF_ASTIGMATISM, EMDL_CTF_CS, EMDL_CTF_DEFOCUSU, EMDL_CTF_DEFOCUSV, &
    EMDL_CTF_DEFOCUS_ANGLE, EMDL_CTF_IMAGE, EMDL_CTF_PHASESHIFT, EMDL_CTF_Q0, EMDL_CTF_VOLTAGE, &
    EMDL_IMAGE_COORD_X, EMDL_IMAGE_COORD_Y, EMDL_IMAGE_DIMENSIONALITY, EMDL_IMAGE_NAME, &
    EMDL_IMAGE_OPTICS_GROUP, EMDL_IMAGE_OPTICS_GROUP_NAME, EMDL_IMAGE_PIXEL_SIZE, EMDL_IMAGE_SIZE, &
    EMDL_MICROGRAPH_METADATA_NAME, EMDL_MICROGRAPH_NAME, EMDL_MICROGRAPH_NAME_WODOSE, &
    EMDL_MLMODEL_GROUP_NAME, EMDL_ORIENT_ORIGIN_X_ANGSTROM, EMDL_ORIENT_ORIGIN_Y_ANGSTROM, &
    EMDL_ORIENT_PSI, EMDL_ORIENT_ROT, EMDL_ORIENT_TILT, EMDL_PARTICLE_RANDOM_SUBSET
use simple_epu_xml,    only: read_epu_beamshift
implicit none
private
public :: relion_project
#include "simple_local_flags.inc"

integer,           parameter :: MIC_EXT_LEN    = len_trim(INTGMOV_SUFFIX) + len(MRC_EXT)
integer,           parameter :: FORCTF_EXT_LEN = len_trim(FORCTF_SUFFIX)  + len(MRC_EXT)
integer,           parameter :: STK_EXT_LEN    = len_trim(INTGMOV_SUFFIX) + len(STK_EXT)
integer,           parameter :: STK_PREFIX_LEN = len_trim(EXTRACT_STK_FBODY)
real,              parameter :: MICRON2ANGS    = 10000.0   ! 1 µm = 10000 Å
character(len=12), parameter :: MICPATH        = 'micrographs/'
character(len=10), parameter :: PTCLPATH       = 'particles/'
character(len=7),  parameter :: MOVIEPATH      = 'movies/'

type relion_project
    integer                   :: opticsgroups
    type(string), allocatable :: movienames(:)
    integer,      allocatable :: moviegroup(:)
contains
    procedure         :: create
    procedure         :: write_corrected_micrographs_star
    procedure         :: write_micrographs_star
    procedure         :: write_particles2D_star
    procedure         :: generate_epu_tiltgroups
    procedure         :: generate_xml_tiltgroups
    procedure         :: find_movienames
    procedure         :: allocate_opticsgroups
    procedure         :: generate_single_tiltgroup
    procedure, nopass :: h_clust ! hierarchical clustering
end type relion_project

contains

    subroutine write_corrected_micrographs_star(self, cline, spproj)
        class(relion_project),  intent(inout) :: self
        class(sp_project),      intent(inout) :: spproj
        class(cmdline),         intent(inout) :: cline
        type(string)              :: getstring, getstring_base
        type(starfile_table_type) :: table
        integer :: i,j
        logical :: exists
        logical :: state, movie, intg, smpd, kv, cs, fraca, opticsgroup
        if(spproj%os_mic%get_noris() == 0) then
            return
        endif
        write(logfhandle, *) 'Generating micrographs_corrected.star ... '
        call simple_mkdir('micrographs')
        call simple_mkdir('movies')
        state = .false.
        do i= 1,spproj%os_mic%get_noris()
            if(spproj%os_mic%isthere(i, 'state')) then
                state = .TRUE.
                if(spproj%os_mic%get(i,'state') .GT. 0) then
                    exit
                endif
            else
                exit
            endif
        end do
        i = min(i, spproj%os_mic%get_noris())
        movie       = spproj%os_mic%isthere(i, 'movie')
        intg        = spproj%os_mic%isthere(i, 'intg')
        smpd        = spproj%os_mic%isthere(i, 'smpd')
        kv          = spproj%os_mic%isthere(i, 'kv')
        cs          = spproj%os_mic%isthere(i, 'cs')
        fraca       = spproj%os_mic%isthere(i, 'fraca')
        opticsgroup = spproj%os_mic%isthere(i, 'opticsgroup')
        call starfile_table__new(table)
        call starfile_table__setName(table, "optics")
        do j=1, self%opticsgroups
            call starfile_table__addObject(table)
            call starfile_table__setValue_int(table, EMDL_IMAGE_OPTICS_GROUP, j + int(cline%get_rarg('optics_offset')))
            call starfile_table__setValue_string(table, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int((j + int(cline%get_rarg('optics_offset'))), 4))))
            if(fraca) call starfile_table__setValue_double(table, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
            if(cs) call starfile_table__setValue_double(table, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
            if(kv) call starfile_table__setValue_double(table, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
            if(smpd) call starfile_table__setValue_double(table, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
        end do
        call starfile_table__open_ofile(table, 'micrographs_corrected.star.tmp', 0)
        call starfile_table__write_ofile(table)
        call starfile_table__close_ofile(table)
        call starfile_table__delete(table)
        call starfile_table__new(table)
        call starfile_table__setName(table, "micrographs")
        !STAR data
        do i=1, spproj%os_mic%get_noris()
            if((state .AND. (spproj%os_mic%get(i,'state') .GT. 0)) .OR. (.NOT. state)) then
                call starfile_table__addObject(table)
                if(intg) then
                    call spproj%os_mic%getter(i, 'intg', getstring)
                    call append_table_with_string(getstring, table, MIC_EXT_LEN, MICPATH, MRC_EXT, EMDL_MICROGRAPH_NAME)
                    if( spproj%os_mic%isthere(i, 'mc_starfile') )then
                        call spproj%os_mic%getter(i, 'mc_starfile', getstring)
                        call append_table_with_string(getstring, table, len(STAR_EXT), MICPATH, STAR_EXT, EMDL_MICROGRAPH_METADATA_NAME)
                    endif
                endif
                if(movie) then
                    call spproj%os_mic%getter(i, 'movie', getstring)
                    getstring_base = basename(getstring)
                    inquire(file=MOVIEPATH // getstring_base%to_char() // ' ', exist=exists)
                    if(.NOT. exists) then
                        call simple_rmfile(MOVIEPATH // getstring_base%to_char())
                        getstring = simple_abspath(getstring, check_exists=.false.)
                        call syslib_symlink(getstring, string(MOVIEPATH // getstring_base%to_char()))
                    endif
                endif
                if(opticsgroup)then
                    call starfile_table__setValue_int(table, EMDL_IMAGE_OPTICS_GROUP, &
                        &int(spproj%os_mic%get(i, 'opticsgroup') + int(cline%get_rarg('optics_offset'))))
                endif
            endif
        end do
        call starfile_table__open_ofile(table, 'micrographs_corrected.star.tmp', 1)
        call starfile_table__write_ofile(table)
        call starfile_table__close_ofile(table)
        ! written as .tmp and renamed when complete: a reader never meets a partial file
        call simple_rename(string('micrographs_corrected.star.tmp'), string('micrographs_corrected.star'))
        call starfile_table__delete(table)
        call getstring%kill
        call getstring_base%kill
    end subroutine write_corrected_micrographs_star

    subroutine write_micrographs_star(self, cline, spproj)
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        type(string)              :: getstring, getstring_base
        type(starfile_table_type) :: table
        real(dp)                  :: dfx, dfy
        integer                   :: i,j
        logical                   :: exists
        logical :: state, intg, forctf, pspec, smpd, kv, cs, l_dfx, l_dfy, angast, fraca, opticsgroup
        if(spproj%os_mic%get_noris() == 0) then
            return
        endif
        write(logfhandle, *) 'Generating micrographs.star ... '
        call simple_mkdir('micrographs')
        state = .false.
        do i=1, spproj%os_mic%get_noris()
            if(spproj%os_mic%isthere(i, 'state')) then
                state = .TRUE.
                if(spproj%os_mic%get(i,'state') .GT. 0) then
                    exit
                endif
            else
                exit
            endif
        end do
        i = min(i, spproj%os_mic%get_noris())
        intg        = spproj%os_mic%isthere(i, 'intg')
        forctf      = spproj%os_mic%isthere(i, 'forctf')
        pspec       = spproj%os_mic%isthere(i, 'pspec')
        smpd        = spproj%os_mic%isthere(i, 'smpd')
        kv          = spproj%os_mic%isthere(i, 'kv')
        cs          = spproj%os_mic%isthere(i, 'cs')
        l_dfx       = spproj%os_mic%isthere(i, 'dfx')
        l_dfy       = spproj%os_mic%isthere(i, 'dfy')
        angast      = spproj%os_mic%isthere(i, 'angast')
        fraca       = spproj%os_mic%isthere(i, 'fraca')
        opticsgroup = spproj%os_mic%isthere(i, 'opticsgroup')
        call starfile_table__new(table)
        call starfile_table__setName(table, "optics")
        do j=1, self%opticsgroups
            call starfile_table__addObject(table)
            call starfile_table__setValue_int(table, EMDL_IMAGE_OPTICS_GROUP, j + int(cline%get_rarg('optics_offset')))
            call starfile_table__setValue_string(table, EMDL_IMAGE_OPTICS_GROUP_NAME, &
                &'opticsgroup' // trim(int2str(int((j + int(cline%get_rarg('optics_offset'))), 4))))
            if(fraca) call starfile_table__setValue_double(table, EMDL_CTF_Q0,           real(spproj%os_mic%get(i, 'fraca'), dp))
            if(   cs) call starfile_table__setValue_double(table, EMDL_CTF_CS,           real(spproj%os_mic%get(i, 'cs'), dp))
            if(   kv) call starfile_table__setValue_double(table, EMDL_CTF_VOLTAGE,      real(spproj%os_mic%get(i, 'kv'), dp))
            if( smpd) call starfile_table__setValue_double(table, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
        end do
        call starfile_table__open_ofile(table, 'micrographs.star.tmp', 0)
        call starfile_table__write_ofile(table)
        call starfile_table__close_ofile(table)
        call starfile_table__delete(table)
        call starfile_table__new(table)
        call starfile_table__setName(table, "micrographs")
        !STAR data
        do i=1, spproj%os_mic%get_noris()
            if((state .AND. (spproj%os_mic%get(i,'state') .GT. 0)) .OR. (.NOT. state)) then
                call starfile_table__addObject(table)
                if(intg) then
                    call spproj%os_mic%getter(i, 'intg', getstring)
                    call append_table_with_string(getstring, table, MIC_EXT_LEN, MICPATH, MRC_EXT, EMDL_MICROGRAPH_NAME )
                endif
                if(forctf) then
                    call spproj%os_mic%getter(i, 'forctf', getstring)
                    call append_table_with_string(getstring, table, FORCTF_EXT_LEN, MICPATH, MRC_EXT, EMDL_MICROGRAPH_NAME_WODOSE )
                endif
                if(pspec)then
                    call spproj%os_mic%getter(i, 'pspec', getstring)
                    getstring_base = basename(getstring)
                    call starfile_table__setValue_string(table, EMDL_CTF_IMAGE, MICPATH // getstring_base%to_char())
                    inquire(file=MICPATH // getstring_base%to_char(), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../'//getstring%to_char(), MICPATH//getstring_base%to_char())
                    endif
                endif
                if(l_dfx .AND. l_dfy) then
                    dfx = real(spproj%os_mic%get_dfx(i) * MICRON2ANGS, dp)
                    dfy = real(spproj%os_mic%get_dfy(i) * MICRON2ANGS, dp)
                    call starfile_table__setValue_double(table, EMDL_CTF_ASTIGMATISM, abs(dfx-dfy))
                    call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUSU, dfx)
                    call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUSV, dfy)
                endif
                if(angast) call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_mic%get(i, 'angast'), dp))
                call starfile_table__setValue_double(table, EMDL_CTF_PHASESHIFT, &
                    &real(rad2deg(spproj%os_mic%get(i, 'phshift')), dp))
                if(opticsgroup) call starfile_table__setValue_int(table, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'opticsgroup') + int(cline%get_rarg('optics_offset'))))
            endif
        end do
        call starfile_table__open_ofile(table, 'micrographs.star.tmp', 1)
        call starfile_table__write_ofile(table)
        call starfile_table__close_ofile(table)
        ! written as .tmp and renamed when complete: a reader never meets a partial file
        call simple_rename(string('micrographs.star.tmp'), string('micrographs.star'))
        call starfile_table__delete(table)
        call getstring%kill
        call getstring_base%kill
    end subroutine write_micrographs_star

    subroutine append_table_with_string( str, table, len_suffix, path, new_suffix, field )
        class(string),              intent(in)    :: str
        integer(C_int),             intent(in)    :: field
        integer,                    intent(in)    :: len_suffix
        character(len=*),           intent(in)    :: path, new_suffix
        class(starfile_table_type), intent(inout) :: table
        type(string) :: str_base, filename
        integer :: lgstr
        logical :: exists
        str_base = basename(str)
        lgstr          = len(str_base%to_char())
        filename       = str_base%to_char([1,lgstr-len_suffix]) // trim(new_suffix)
        call starfile_table__setValue_string(table, field, trim(path)//filename%to_char())
        inquire(file=trim(path) // filename%to_char(), exist=exists)
        if(.NOT. exists) then
            if( str%to_char([1,1]) == '/' )then
                call syslib_symlink(str, trim(path)//filename%to_char())
            else
                call syslib_symlink('../'//str%to_char(), trim(path)//filename%to_char())
            endif
        endif
        call str_base%kill
        call filename%kill
    end subroutine append_table_with_string

    subroutine write_particles2D_star(self, cline, spproj)
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        type(string)              :: getstring, getstring_base, getstring_mrcs
        character(len=1024)       :: groupname
        ! integer, allocatable      :: ptclcount(:)
        integer                   :: group
        type(starfile_table_type) :: table
        real(dp)                  :: dfx, dfy
        integer                   :: i, j, stkindex, nbase, ind_in_stk
        logical                   :: exists
        logical :: state, stkstate, stkind, smpd, kv, cs, l_dfx, l_dfy, angast, fraca
        logical :: ptcldfx, ptcldfy, ptclangast, coords, stk, box, opticsgroup, eo
        real                      :: dfxmin
        real                      :: dfxmax
        real                      :: dfxstep
        logical                   :: l_poses3D
        if(spproj%os_ptcl2D%get_noris() == 0) then
            return
        endif
        ! 3D poses when the 3D field holds them, else the 2D field's in-plane poses
        l_poses3D = .false.
        if( spproj%os_ptcl3D%get_noris() == spproj%os_ptcl2D%get_noris() )then
            l_poses3D = spproj%has_valid_ptcl3D_alignment()
        endif
        write(logfhandle, *) 'Generating particles2D.star ... '
        call simple_mkdir('particles')
        state    = .false.
        stkstate = .false.
        coords   = .false.
        do i=1, spproj%os_ptcl2D%get_noris()
            if(spproj%os_ptcl2D%isthere(i, 'state')) then
                state = .TRUE.
                if(spproj%os_ptcl2D%get(i,'state') .GT. 0) then
                    exit
                endif
            else
                exit
            endif
        end do
        i = min(i, spproj%os_ptcl2D%get_noris())
        stkind     = spproj%os_ptcl2D%isthere(i, 'stkind')
        ptcldfx    = spproj%os_ptcl2D%isthere(i, 'dfx')
        ptcldfy    = spproj%os_ptcl2D%isthere(i, 'dfy')
        ptclangast = spproj%os_ptcl2D%isthere(i, 'angast')
        eo         = spproj%os_ptcl2D%get_nodd() > 0
        !trial up to 1000 particles for non 0 values of xpos and ypos
        do j=i, min(i + 1000, spproj%os_ptcl2D%get_noris())
          if((spproj%os_ptcl2D%get(j, 'xpos') .GT. 0) .OR. (spproj%os_ptcl2D%get(j, 'ypos') .GT. 0)) then
            coords = .TRUE.
            exit
          endif
        end do
        do i=1, spproj%os_stk%get_noris()
            if(spproj%os_stk%isthere(i, 'state')) then
                stkstate = .TRUE.
                if(spproj%os_stk%get(i,'state') .GT. 0) then
                    exit
                endif
            else
                exit
            endif
        end do
        i = min(i, spproj%os_stk%get_noris())
        l_dfx       = spproj%os_stk%isthere(i, 'dfx')
        l_dfy       = spproj%os_stk%isthere(i, 'dfy')
        angast      = spproj%os_stk%isthere(i, 'angast')
        box         = spproj%os_stk%isthere(i, 'box')
        smpd        = spproj%os_stk%isthere(i, 'smpd')
        kv          = spproj%os_stk%isthere(i, 'kv')
        cs          = spproj%os_stk%isthere(i, 'cs')
        fraca       = spproj%os_stk%isthere(i, 'fraca')
        stk         = spproj%os_stk%isthere(i, 'stk')
        opticsgroup = spproj%os_stk%isthere(i, 'opticsgroup')
        call starfile_table__new(table)
        call starfile_table__setName(table, "optics")
        do j=1, self%opticsgroups
            call starfile_table__addObject(table)
            call starfile_table__setValue_int(table, EMDL_IMAGE_OPTICS_GROUP, j + int(cline%get_rarg('optics_offset')))
            call starfile_table__setValue_string(table, EMDL_IMAGE_OPTICS_GROUP_NAME, &
                &'opticsgroup' // trim(int2str(int((j + int(cline%get_rarg('optics_offset'))), 4))))
            if(fraca) call starfile_table__setValue_double(table, EMDL_CTF_Q0,           real(spproj%os_stk%get(i, 'fraca'), dp))
            if(   cs) call starfile_table__setValue_double(table, EMDL_CTF_CS,           real(spproj%os_stk%get(i, 'cs'), dp))
            if(   kv) call starfile_table__setValue_double(table, EMDL_CTF_VOLTAGE,      real(spproj%os_stk%get(i, 'kv'), dp))
            if( smpd) call starfile_table__setValue_double(table, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_stk%get(i, 'smpd'), dp))
            if(  box) call starfile_table__setValue_int(table, EMDL_IMAGE_SIZE,          int(spproj%os_stk%get(i, 'box')))
            call starfile_table__setValue_int(table, EMDL_IMAGE_DIMENSIONALITY, 2)
        end do
        call starfile_table__open_ofile(table, 'particles2D.star.tmp', 0)
        call starfile_table__write_ofile(table)
        call starfile_table__close_ofile(table)
        call starfile_table__delete(table)
        call starfile_table__new(table)
        call starfile_table__setName(table, "particles")
        dfxstep = 1.
        if(cline%get_iarg('reliongroups') > 0 .AND. l_dfx) then
            call spproj%os_stk%minmax('dfx', dfxmin, dfxmax)
            dfxstep = (dfxmax - dfxmin)/ real(cline%get_iarg('reliongroups') + 1)
        endif
        !STAR data
        do i = 1,spproj%os_ptcl2D%get_noris()
            if(stkind) then
                call spproj%map_ptcl_ind2stk_ind('ptcl2D', i, stkindex, ind_in_stk )
                if((state .AND. (spproj%os_ptcl2D%get_state(i) .GT. 0)) .OR. (.NOT. state)) then
                    if((stkstate .AND. (spproj%os_stk%get_state(stkindex) .GT. 0)) .OR. (.NOT. stkstate)) then
                        call starfile_table__addObject(table)
                        if(stk) then
                            call spproj%os_stk%getter(stkindex, 'stk', getstring)
                            getstring_base = basename(getstring)
                            getstring_mrcs = fname_new_ext(getstring_base, STK_EXT(2:))
                            call starfile_table__setValue_string(table, EMDL_IMAGE_NAME, &
                                &trim(int2str(int(ind_in_stk,4)))//'@'//PTCLPATH//getstring_mrcs%to_char())
                            inquire(file=PTCLPATH//getstring_mrcs%to_char(), exist=exists)
                            if(.NOT. exists) then
                                if( getstring%to_char([1,1]) == '/' )then
                                    call syslib_symlink( getstring%to_char(), PTCLPATH//getstring_mrcs%to_char())
                                else
                                    call syslib_symlink('../'//getstring%to_char(), PTCLPATH//getstring_mrcs%to_char())
                                endif
                            endif
                            nbase = len(getstring_base%to_char())
                            if( nbase > STK_PREFIX_LEN + STK_EXT_LEN )then
                                if( getstring_base%to_char([1,STK_PREFIX_LEN]) == EXTRACT_STK_FBODY )then
                                    if( getstring_base%to_char([nbase-STK_EXT_LEN+1,nbase]) == INTGMOV_SUFFIX//STK_EXT )then
                                        call starfile_table__setValue_string(table, EMDL_MICROGRAPH_NAME, &
                                            &MICPATH//getstring_base%to_char([STK_PREFIX_LEN+1,nbase-STK_EXT_LEN])//MRC_EXT)
                                    endif
                                endif
                            endif
                        endif
                        if(coords .AND. box) then
                            call starfile_table__setValue_double(table, EMDL_IMAGE_COORD_X, &
                                &real(int(spproj%os_ptcl2D%get(i,'xpos') + (spproj%os_stk%get(stkindex,'box') / 2)), dp))
                            call starfile_table__setValue_double(table, EMDL_IMAGE_COORD_Y, &
                                &real(int(spproj%os_ptcl2D%get(i,'ypos') + (spproj%os_stk%get(stkindex,'box') / 2)), dp))
                        endif
                        if(ptcldfx .AND. ptcldfy .AND. ptclangast) then
                            dfx = real(spproj%os_ptcl2D%get_dfx(i) * MICRON2ANGS, dp)
                            dfy = real(spproj%os_ptcl2D%get_dfy(i) * MICRON2ANGS, dp)
                            call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUSU, dfx)
                            call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUSV, dfy)
                            call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_ptcl2D%get(i,'angast'), dp))
                        else if(l_dfx .AND. l_dfy .AND. angast) then
                            dfx = real(spproj%os_stk%get_dfx(stkindex) * MICRON2ANGS, dp)
                            dfy = real(spproj%os_stk%get_dfy(stkindex) * MICRON2ANGS, dp)
                            call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUSU, dfx)
                            call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUSV, dfy)
                            call starfile_table__setValue_double(table, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_stk%get(stkindex,'angast'), dp))
                        endif
                        call starfile_table__setValue_double(table, EMDL_CTF_PHASESHIFT, &
                            &real(rad2deg(spproj%os_ptcl2D%get(i, 'phshift')), dp))
                        if( l_poses3D )then
                            call write_pose(spproj%os_ptcl3D)
                        else
                            call write_pose(spproj%os_ptcl2D)
                        endif
                        if(cline%get_rarg('reliongroups') > 0 .AND. l_dfx) then
                            group = ceiling((spproj%os_stk%get_dfx(stkindex) - dfxmin) / dfxstep)
                            write(groupname, *) group
                            call starfile_table__setValue_string(table, EMDL_MLMODEL_GROUP_NAME, trim(adjustl(groupname)))
                        endif
                        if( opticsgroup )then
                            call starfile_table__setValue_int(table, EMDL_IMAGE_OPTICS_GROUP, spproj%os_stk%get_int(stkindex, 'opticsgroup') + cline%get_iarg('optics_offset'))
                        endif
                        if( eo )then
                            call starfile_table__setValue_int(table, EMDL_PARTICLE_RANDOM_SUBSET, spproj%os_ptcl2D%get_eo(i)+1)
                        endif
                    endif
                endif
            endif
        end do
        call starfile_table__open_ofile(table, 'particles2D.star.tmp', 1)
        call starfile_table__write_ofile(table)
        call starfile_table__close_ofile(table)
        ! written as .tmp and renamed when complete: a reader never meets a partial file
        call simple_rename(string('particles2D.star.tmp'), string('particles2D.star'))
        call starfile_table__delete(table)
        call getstring%kill
        call getstring_base%kill
        call getstring_mrcs%kill
        contains

            ! the pose of particle i in os; origins from pixels of its stack to Angstrom
            subroutine write_pose( os )
                class(oris), intent(in) :: os
                character(len=XLONGSTRLEN) :: msg, stkname
                real :: smpd, x, y
                if( l_poses3D )then
                    call starfile_table__setValue_double(table, EMDL_ORIENT_ROT,  real(os%e1get(i), dp))
                    call starfile_table__setValue_double(table, EMDL_ORIENT_TILT, real(os%e2get(i), dp))
                endif
                call starfile_table__setValue_double(table, EMDL_ORIENT_PSI, real(os%e3get(i), dp))
                x = os%get(i, 'x')
                y = os%get(i, 'y')
                if( x /= 0. .or. y /= 0. )then
                    smpd = spproj%os_stk%get(stkindex, 'smpd')
                    if( smpd <= 0. )then
                        stkname = ''
                        if( spproj%os_stk%isthere(stkindex, 'stk') ) call spproj%os_stk%get_static(stkindex, 'stk', stkname)
                        write(msg,'(A,I0,A,I0,3A)') 'particle ', i, ' of stack ', stkindex, ' (', trim(stkname),&
                            &') has no pixel size to write its origin in Angstrom'
                        THROW_HARD(trim(msg))
                    endif
                    x = x * smpd
                    y = y * smpd
                endif
                call starfile_table__setValue_double(table, EMDL_ORIENT_ORIGIN_X_ANGSTROM, real(x, dp))
                call starfile_table__setValue_double(table, EMDL_ORIENT_ORIGIN_Y_ANGSTROM, real(y, dp))
            end subroutine write_pose

    end subroutine write_particles2D_star

    subroutine find_movienames(self, spproj)
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        integer      :: moviecount, i
        logical      :: micsource
        type(string) :: moviename, tmpname
        if (spproj%os_mic%get_noris() .gt. 0) then
            moviecount = spproj%os_mic%get_noris()
            micsource = .TRUE.
        else if (spproj%os_stk%get_noris() .gt. 0) then
            moviecount = spproj%os_stk%get_noris()
            micsource = .FALSE.
        else
            THROW_HARD('no micrographs or stacks in project file')
        endif
        if(allocated(self%movienames)) then
            call self%movienames%kill
            deallocate(self%movienames)
        endif
        if(allocated(self%moviegroup)) deallocate(self%moviegroup)
        allocate(self%movienames(moviecount), self%moviegroup(moviecount))
        do i=1, moviecount
            if(micsource) then
                moviename = basename(spproj%os_mic%get_str(i, 'intg'))
                tmpname   = moviename%to_char([1,moviename%strlen_trim()-MIC_EXT_LEN])
            else
                moviename = basename(spproj%os_stk%get_str(i, 'stk'))
                tmpname   = moviename%to_char([STK_PREFIX_LEN+1,moviename%strlen_trim()-STK_EXT_LEN])
            endif
            self%movienames(i) = tmpname
        end do
        call moviename%kill
        call tmpname%kill
    end subroutine find_movienames

    subroutine generate_epu_tiltgroups(self)
        class(relion_project),  intent(inout)   :: self
        integer                                 :: i,j,k
        character(len=20), allocatable :: tiltgroups(:)
        type(string) :: tiltname, tmpname
        if(.NOT. allocated(tiltgroups)) then
            allocate(tiltgroups(50))            ! Max 50 raw groups set here(ie per stage movement)
        end if
        do i=1, size(tiltgroups)
            tiltgroups(i) = ''
        end do
        j = 1
        do i=1, size(self%movienames)
            tiltname = self%movienames(i)
            tmpname  = tiltname%to_char([index(tiltname%to_char(),'Data_')+5,tiltname%strlen_trim()])
            tiltname = tmpname%to_char([1,index(tmpname%to_char(),'_')-1])
            if(.NOT. any(tiltgroups .eq. tiltname%to_char())) then
                tiltgroups(j) = tiltname%to_char()
                self%moviegroup(i) = j
                j = j+1
            else
                do k=1, size(tiltgroups)
                    if(tiltgroups(k) .eq. tiltname%to_char()) then
                        exit
                    endif
                end do
                self%moviegroup(i) = k
            endif
        end do
        self%opticsgroups = j - 1
        call tiltname%kill
        call tmpname%kill
        if(allocated(tiltgroups))deallocate(tiltgroups)
    end subroutine generate_epu_tiltgroups

    subroutine generate_single_tiltgroup(self)
        class(relion_project),  intent(inout)   :: self
        integer                                 :: i
        do i=1, size(self%movienames)
            self%moviegroup(i) = 1
        end do
        self%opticsgroups = 1
    end subroutine generate_single_tiltgroup

    subroutine generate_xml_tiltgroups(self, cline)
        use CPlot2D_wrapper_module, only: CPlot2D_type, CDataSet_type, CPlot2D__new, CPlot2D__delete, &
            CPlot2D__SetXAxisSize, CPlot2D__SetYAxisSize, CPlot2D__SetDrawXAxisGridLines, &
            CPlot2D__SetDrawYAxisGridLines, CPlot2D__SetDrawLegend, CPlot2D__SetFlipY, &
            CPlot2D__SetXAxisTitle, CPlot2D__SetYAxisTitle, CPlot2D__AddDataSet, &
            CPlot2D__OutputPostScriptPlot, CDataSet__new, CDataSet__delete, CDataSet__SetDrawMarker, &
            CDataSet__SetMarkerSize, CDataSet__SetDatasetColor, CDataSet_addpoint
        class(relion_project),  intent(inout)   :: self
        class(cmdline),         intent(inout)   :: cline
        character(len=20), allocatable :: tiltgroups(:)
        integer              :: i,j
        type(string)         :: tiltname, tmpname, fname, fname_eps, xmlloc_file
        type(string)         :: title
        type(CPlot2D_type)   :: plot2D
        type(CDataSet_type)  :: dataSet
        real                 :: beamtiltx, beamtilty, threshold
        logical              :: l_beamshift
        real,    allocatable :: tilts(:,:), centroids(:,:)
        integer, allocatable :: populations(:) 
        integer              :: pos
        call seed_rnd
        write(logfhandle, *) "Parsing movie XML files ... "
        if(cline%get_carg('xmlloc') .eq. '') then
            THROW_HARD('xmlloc is not set')
        endif
        if(cline%get_rarg('tilt_thres') <= 0) then
            threshold = 0.05
        else
            threshold = cline%get_rarg('tilt_thres')
        endif
        xmlloc_file = cline%get_carg('xmlloc')
        if(.not. file_exists(xmlloc_file)) THROW_HARD('xmlloc does not exist')
        if(.NOT. allocated(tilts)) then
           allocate(tilts(size(self%movienames), 2))
        end if
        if(.NOT. allocated(self%moviegroup)) then
           allocate(self%moviegroup(size(self%movienames)))
        end if
        if(.NOT. allocated(tiltgroups)) then
            allocate(tiltgroups(50))            ! Max 50 raw groups set here(ie per stage movement)
        end if
        do i=1, size(tiltgroups)
            tiltgroups(i) = ''
        end do
        j = 1
        do i=1, size(self%movienames)
            pos = index(self%movienames(i)%to_char(), '_fractions')
            if( pos .ne. 0 ) then
                tmpname = xmlloc_file%to_char() // "/" // self%movienames(i)%to_char([1,pos-1]) // ".xml"
            else
                tmpname = xmlloc_file%to_char() // "/" // self%movienames(i)%to_char() //".xml"
            endif
            fname = tmpname%substr_remove(string('_EER'))
            if(.NOT. file_exists(fname) )then
                write(logfhandle, *) fname%to_char() // ' does not exist. Ignoring'
                tilts(i,1) = 0.0
                tilts(i,2) = 0.0
            else
                call read_epu_beamshift(fname%to_char(), beamtiltx, beamtilty, l_beamshift)
                if( .not. l_beamshift ) write(logfhandle, *) fname%to_char() // ' holds no BeamShift. Ignoring'
                tilts(i,1) = beamtiltx
                tilts(i,2) = beamtilty
            endif  
            if(index(self%movienames(1)%to_char(), 'FoilHole') .ne. 0) then
                tiltname = self%movienames(i)
                tmpname  = tiltname%to_char([index(tiltname%to_char(),'Data_')+5,tiltname%strlen_trim()])
                tiltname = tmpname%to_char([1,index(tmpname%to_char(),'_')-1])
                if(.NOT. any(tiltgroups .eq. tiltname%to_char())) then
                    tiltgroups(j) = tiltname%to_char()
                    j = j+1
                end if
            end if
        end do
        call self%h_clust(tilts, threshold, self%moviegroup, centroids, populations) ! threshold needs to be variable!
        if(index(self%movienames(1)%to_char(), 'FoilHole') .ne. 0) then
            self%opticsgroups = int(size(populations) * (j - 1))
            do i=1, size(self%movienames)
                tiltname = self%movienames(i)
                tmpname  = tiltname%to_char([index(tiltname%to_char(),'Data_')+5,tiltname%strlen_trim()])
                tiltname = tmpname%to_char([1,index(tmpname%to_char(),'_')-1]) 
                do j=1, size(tiltgroups)
                    if(tiltgroups(j) .eq. tiltname%to_char()) then
                        exit
                    endif
                end do
                self%moviegroup(i) = self%moviegroup(i)  + (size(populations) * (j - 1))
            end do
        else
            self%opticsgroups = int(size(populations))
        end if
        call CPlot2D__new(plot2D, 'Optics Groups'//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.0_c_double)
        call CPlot2D__SetYAxisSize(plot2D, 400.0_c_double)
        call CPlot2D__SetDrawXAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetDrawYAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        title = 'Beamshift_X'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%to_char())
        title = 'Beamshift_Y'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%to_char())
        do i=1, self%opticsgroups
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_TRUE)
            call CDataSet__SetMarkerSize(dataSet, real(3.0, c_double))
            call CDataSet__SetDatasetColor(dataSet, real(ran3(), c_double), real(ran3(), c_double), real(ran3(), c_double))
            do j=1, size(self%moviegroup)
                if(self%moviegroup(j) .eq. i) then
                    call CDataSet_addpoint(dataSet, tilts(j,1), tilts(j,2))
                end if
            end do
            call CPlot2D__AddDataSet(plot2D, dataset)
            call CDataSet__delete(dataset)
        end do
        fname_eps = 'optics_groups.eps'//C_NULL_CHAR
        call CPlot2D__OutputPostScriptPlot(plot2D, fname_eps%to_char())
        call CPlot2D__delete(plot2D)
        call tiltname%kill
        call tmpname%kill
        call fname%kill
        call fname_eps%kill
        call xmlloc_file%kill
        call title%kill
        if(allocated(tilts)      )deallocate(tilts)
        if(allocated(centroids)  )deallocate(centroids)
        if(allocated(populations))deallocate(populations)
        if(allocated(tiltgroups) )deallocate(tiltgroups)
    end subroutine generate_xml_tiltgroups

    subroutine allocate_opticsgroups(self, cline, spproj)
        class(relion_project),  intent(inout) :: self
        class(sp_project),      intent(inout) :: spproj
        class(cmdline),         intent(inout) :: cline
        integer, allocatable :: opticsgroupmap(:,:)
        type(string) :: moviename
        integer      :: i, j
        real         :: tiltgroupmax
        tiltgroupmax = cline%get_rarg('tiltgroupmax')
        if(.NOT. allocated(opticsgroupmap)) then
            allocate(opticsgroupmap(2, self%opticsgroups)) !1:count, 2: current mapped cluster id
        end if
        do i=1, self%opticsgroups
            opticsgroupmap(1,i) = 0
            opticsgroupmap(2,i) = i
        end do
        do i=1, size(self%movienames)
            opticsgroupmap(1,self%moviegroup(i)) = opticsgroupmap(1,self%moviegroup(i)) + 1
            if(tiltgroupmax > 0 .AND. opticsgroupmap(1,self%moviegroup(i)) .gt. tiltgroupmax) then
                self%opticsgroups = self%opticsgroups + 1
                opticsgroupmap(2,self%moviegroup(i)) = self%opticsgroups
                opticsgroupmap(1,self%moviegroup(i)) = 1
            end if
            self%moviegroup(i) = opticsgroupmap(2,self%moviegroup(i))
        end do
        do i=1, spproj%os_mic%get_noris()
            moviename = basename(spproj%os_mic%get_str(i, 'intg'))
            moviename = moviename%to_char([1,moviename%strlen_trim()-MIC_EXT_LEN])
            do j=1, size(self%movienames)
                if(self%movienames(j) .eq. moviename) then
                    exit
                end if
            end do
            call spproj%os_mic%set(i, 'opticsgroup', self%moviegroup(j))
        end do
        do i=1, spproj%os_stk%get_noris()
            moviename = basename(spproj%os_stk%get_str(i,'stk'))
            moviename = moviename%to_char([STK_PREFIX_LEN+1,moviename%strlen_trim()-STK_EXT_LEN])
            do j=1, size(self%movienames)
                if(self%movienames(j) .eq. moviename) then
                    exit
                end if
            end do
            call spproj%os_stk%set(i, 'opticsgroup', self%moviegroup(j))
        end do
        call moviename%kill
    end subroutine allocate_opticsgroups

    subroutine create(self, cline, spproj)
        class(relion_project), intent(inout) :: self
        class(sp_project),     intent(inout) :: spproj
        class(cmdline),        intent(inout) :: cline
        write(logfhandle, *) '>>> Writing Relion 5.0 compatible STAR files'
        call self%find_movienames(spproj)
        if(cline%get_carg('xmlloc') .eq. '') then
            if(index(self%movienames(1)%to_char(), 'FoilHole') .ne. 0) then
                write(logfhandle, *) '>>> Using EPU filenames for beamtilt information'
                call self%generate_epu_tiltgroups()
            else
                write(logfhandle, *) '>>> Creating single beamtilt group'
                call self%generate_single_tiltgroup()
            end if
        else
            if(index(self%movienames(1)%to_char(), 'FoilHole') .ne. 0) then
                write(logfhandle, *) '>>> Using XML files and EPU filenames for beamtilt information'
            else
                write(logfhandle, *) '>>> Using XML files for beamtilt information'
            endif
            call self%generate_xml_tiltgroups(cline)
        endif
        call self%allocate_opticsgroups(cline, spproj)
        call self%write_corrected_micrographs_star(cline, spproj)
        call self%write_micrographs_star(cline, spproj)
        call self%write_particles2D_star(cline, spproj)
        if( allocated(self%moviegroup) ) deallocate(self%moviegroup)
        if( allocated(self%movienames) ) then
            call self%movienames%kill
            deallocate(self%movienames)
        endif
    end subroutine create

    ! distance threshold based yerarchical clustering
    ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
    subroutine h_clust(data_in,thresh,labels,centroids,populations)
        real,    intent(in)  :: data_in(:,:)   ! input data, point coords
        real,    intent(in)  :: thresh         ! threshold for class merging
        integer, intent(out) :: labels(:)      ! labels of the elements in vec
        real,    allocatable, intent(out) :: centroids(:,:) ! centroids of the classes
        integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
        real,    allocatable :: mat(:,:)                    ! pariwise distances matrix
        logical, allocatable :: mask(:), outliers(:)
        integer :: N, i, j, cnt, ncls
        integer :: index(1), loc1(1), loc2(1)
        real    :: d
        if( size(data_in, dim = 2) .ne. 2 ) THROW_HARD('Input data shouldbe two dimensional!; h_clust')
        N = size(data_in, dim = 1) ! number of data points
        ! 1) calc all the couples of distances, using euclid dist
        allocate(mat(N,N), source = 0.)
        do i = 1, N-1
          do j = i+1, N
            mat(i,j) = sqrt((data_in(i,1)-data_in(j,1))**2 + (data_in(i,2)-data_in(j,2))**2) ! pariwise euclidean distance
            mat(j,i) = mat(i,j)
          enddo
        enddo
        ! 2) Generate binary clusters
        allocate(mask(N),source = .true. )
        allocate(outliers(N), source = .false.)
        ncls = 0
        do i = 1, N
          if(mask(i)) then ! if it's not already been clustered
            mask(i) = .false.
            ! find the index of the couple
            d = minval(mat(i,:), mask)
            index(:) = minloc(mat(i,:), mask)
            ncls = ncls + 1
            ! assign lavels
            labels(i) = ncls
            if(d <= thresh) then ! if it's not an outlier (it has a couple)
              labels(index(1)) = labels(i)
              mask(index(1)) = .false. ! index(1) has already been clustered
            else
              outliers(i) = .true.
            endif
          endif
        enddo
        ! 3) Calculate centroids
        allocate(centroids(ncls,2), source = 0.)
        mask = .true. ! reset
        do i = 1, ncls
          ! find the first member of the class
          loc1(:) = minloc(abs(labels-i))
          if(.not. outliers(loc1(1))) then
            mask(loc1(1)) = .false.
            ! find the second member of the class
            loc2(:) = minloc(abs(labels-i), mask)
            mask(loc2(1)) = .false.
            centroids(i,1) = (data_in(loc1(1),1) + data_in(loc2(1),1))/2.
            centroids(i,2) = (data_in(loc1(1),2) + data_in(loc2(1),2))/2.
          else ! the class has just one element
            loc1(:) = minloc(abs(labels-i))
            centroids(i,1) = data_in(loc1(1),1)
            centroids(i,2) = data_in(loc1(1),2)
            mask(loc1(1)) = .false.
          endif
        enddo
        mask  = .true. ! reset
        ! 4) merge classes
        do i = 1, ncls-1
          do j = i+1, ncls
            if(sqrt((centroids(i,1)-centroids(j,1))**2+(centroids(i,2)-centroids(j,2))**2) <= thresh) then ! merge classes
              ! change label to class j
              where(labels == j) labels = i
            endif
          enddo
        enddo
        ! 5) Reoder labels
        cnt = 0
        do i = 1, ncls
           if(any(labels== i)) then !there is a class labelled i
              cnt = cnt + 1                  !increasin order cc
              where(labels == i) labels = cnt
            endif
        enddo
        ! 6) recalculate centroids
        deallocate(centroids)
        ncls = maxval(labels) ! the nb of classes is maxval(labels)
        allocate(centroids(ncls,2),   source = 0.)
        allocate(populations(ncls), source = 0 )
        mask = .true. ! reset
        do i = 1, ncls
          populations(i) = count(labels == i) ! counts the nb of elements in the class
          ! find the all the cnt member of the class and update the centroids
          do j = 1, populations(i)
            loc1(:) = minloc(abs(labels-i), mask)
            mask(loc1(1)) = .false. ! do not consider this member of the class anymore
            centroids(i,1) = centroids(i,1)+ data_in(loc1(1),1)
            centroids(i,2) = centroids(i,2)+ data_in(loc1(1),2)
          enddo
          centroids(i,:) = centroids(i,:)/real(populations(i))
        enddo
        ! output generation, Joseph please change it as you wish
        write(logfhandle, *) '>>> HIERARCHICAL CLUSTERING COMPLETED'
        write(logfhandle, *) '>>> Distance threshold selected: ', thresh
        do i = 1, ncls
            write(logfhandle, *) '>>>Class ', i, 'Population ', populations(i), ' Centroid [', centroids(i,1), ',', centroids(i,2), ']'
        enddo
    end subroutine h_clust

end module simple_relion

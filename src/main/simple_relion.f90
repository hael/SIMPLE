module simple_relion
include 'simple_lib.f08'
use simple_starfile_wrappers
use simple_sp_project,          only: sp_project
use simple_cmdline,             only: cmdline
use FoX_dom
implicit none
private
public :: relion_project
#include "simple_local_flags.inc"

type relion_project
    integer                                 :: opticsgroups
    character(len=128),     allocatable     :: movienames(:)
    integer,                allocatable     :: moviegroup(:)

contains

    procedure :: create
    procedure :: write_corrected_micrographs_star
    procedure :: write_micrographs_star
    procedure :: write_particles2D_star
    procedure :: generate_epu_tiltgroups
    procedure :: generate_xml_tiltgroups
    procedure :: find_movienames
    procedure :: allocate_opticsgroups
    procedure :: generate_single_tiltgroup
    procedure :: h_clust ! hierarchical clustering

end type relion_project

contains

    subroutine write_corrected_micrographs_star(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        character(len=:),       allocatable     :: getstring
        type(starfile_table_type)               :: startable
        integer                                 :: i,j
        logical                                 :: exists

        logical                                 :: state        = .FALSE.
        logical                                 :: movie        = .FALSE.
        logical                                 :: intg         = .FALSE.
        logical                                 :: smpd         = .FALSE.
        logical                                 :: kv           = .FALSE.
        logical                                 :: cs           = .FALSE.
        logical                                 :: fraca        = .FALSE.
        logical                                 :: opticsgroup  = .FALSE.

        if(spproj%os_mic%get_noris() == 0) then
            return
        endif

        write(logfhandle, *) 'Generating micrographs_corrected.star ... '

        call simple_mkdir('micrographs', errmsg= "simple_relion:: create micrographs directory")
        call simple_mkdir('movies', errmsg= "simple_relion:: create movies directory")

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

        movie       = spproj%os_mic%isthere(i, 'movie')
        intg        = spproj%os_mic%isthere(i, 'intg')
        smpd        = spproj%os_mic%isthere(i, 'smpd')
        kv          = spproj%os_mic%isthere(i, 'kv')
        cs          = spproj%os_mic%isthere(i, 'cs')
        fraca       = spproj%os_mic%isthere(i, 'fraca')
        opticsgroup = spproj%os_mic%isthere(i, 'opticsgroup')

        call starfile_table__new(startable)
        call starfile_table__setName(startable, "optics")

        do j=1, self%opticsgroups
            call starfile_table__addObject(startable)
            call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, j)
            call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int(j, 4))))
            if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
            if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
            if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
            if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
        end do

        call starfile_table__open_ofile(startable, 'micrographs_corrected.star', 0)
        call starfile_table__write_ofile(startable)
        call starfile_table__close_ofile(startable)
        call starfile_table__delete(startable)
        call starfile_table__new(startable)
        call starfile_table__setName(startable, "micrographs")

        !STAR data
        do i=1, spproj%os_mic%get_noris()
            if((state .AND. (spproj%os_mic%get(i,'state') .GT. 0)) .OR. (.NOT. state)) then
                call starfile_table__addObject(startable)

                if(intg) then
                    call spproj%os_mic%getter(i, 'intg', getstring)
                    call starfile_table__setValue_string(startable, EMDL_MICROGRAPH_NAME, 'micrographs/' // trim(adjustl(basename(getstring))))
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                    call starfile_table__setValue_string(startable, EMDL_MICROGRAPH_METADATA_NAME, 'micrographs/' // trim(adjustl(basename(getstring(1 : len(getstring) - 9)))) // '.star')
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring(1 : len(getstring) - 9)))) // '.star', exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring(1 : len(getstring) - 9))) // '.star ', 'micrographs/' // trim(adjustl(basename(getstring(1 : len(getstring) - 9)))) // '.star ', 'Failed to generate symlink')
                    endif
                endif

                if(movie) then
                    call spproj%os_mic%getter(i, 'movie', getstring)
                    inquire(file='movies/' // trim(adjustl(basename(getstring))) // ' ', exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink(trim(adjustl(getstring)), 'movies/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                endif

                if(opticsgroup) call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'opticsgroup')))

            endif
        end do

        call starfile_table__open_ofile(startable, 'micrographs_corrected.star', 1)
        call starfile_table__write_ofile(startable)
        call starfile_table__close_ofile(startable)

        call starfile_table__delete(startable)

        if(allocated(getstring))deallocate(getstring)

    end subroutine write_corrected_micrographs_star

    subroutine write_micrographs_star(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        character(len=:),       allocatable     :: getstring
        type(starfile_table_type)               :: startable
        integer                                 :: i,j
        logical                                 :: exists

        logical                                 :: state        = .FALSE.
        logical                                 :: movie        = .FALSE.
        logical                                 :: intg         = .FALSE.
        logical                                 :: forctf       = .FALSE.
        logical                                 :: pspec        = .FALSE.
        logical                                 :: smpd         = .FALSE.
        logical                                 :: kv           = .FALSE.
        logical                                 :: cs           = .FALSE.
        logical                                 :: dfx          = .FALSE.
        logical                                 :: dfy          = .FALSE.
        logical                                 :: angast       = .FALSE.
        logical                                 :: fraca        = .FALSE.
        logical                                 :: opticsgroup    = .FALSE.

        if(spproj%os_mic%get_noris() == 0) then
            return
        endif

        write(logfhandle, *) 'Generating micrographs.star ... '

        call simple_mkdir('micrographs', errmsg= "simple_relion:: create micrographs directory")

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

        movie       = spproj%os_mic%isthere(i, 'movie')
        intg        = spproj%os_mic%isthere(i, 'intg')
        forctf      = spproj%os_mic%isthere(i, 'forctf')
        pspec       = spproj%os_mic%isthere(i, 'pspec')
        smpd        = spproj%os_mic%isthere(i, 'smpd')
        kv          = spproj%os_mic%isthere(i, 'kv')
        cs          = spproj%os_mic%isthere(i, 'cs')
        dfx         = spproj%os_mic%isthere(i, 'dfx')
        dfy         = spproj%os_mic%isthere(i, 'dfy')
        angast      = spproj%os_mic%isthere(i, 'angast')
        fraca       = spproj%os_mic%isthere(i, 'fraca')
        opticsgroup = spproj%os_mic%isthere(i, 'opticsgroup')


        call starfile_table__new(startable)
        call starfile_table__setName(startable, "optics")

        do j=1, self%opticsgroups
            call starfile_table__addObject(startable)
            call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, j)
            call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int(j, 4))))
            if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
            if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
            if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
            if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
        end do

        call starfile_table__open_ofile(startable, 'micrographs.star', 0)
        call starfile_table__write_ofile(startable)
        call starfile_table__close_ofile(startable)
        call starfile_table__delete(startable)
        call starfile_table__new(startable)
        call starfile_table__setName(startable, "micrographs")

        !STAR data
        do i=1, spproj%os_mic%get_noris()
            if((state .AND. (spproj%os_mic%get(i,'state') .GT. 0)) .OR. (.NOT. state)) then
                call starfile_table__addObject(startable)
                if(intg) then
                    call spproj%os_mic%getter(i, 'intg', getstring)
                    call starfile_table__setValue_string(startable, EMDL_MICROGRAPH_NAME, 'micrographs/' // trim(adjustl(basename(getstring))))
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                endif

                if(forctf) then
                    call spproj%os_mic%getter(i, 'forctf', getstring)
                    call starfile_table__setValue_string(startable, EMDL_MICROGRAPH_NAME_WODOSE, 'micrographs/' // trim(adjustl(basename(getstring))))
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                endif

                if(pspec) then
                    call spproj%os_mic%getter(i, 'pspec', getstring)
                    call starfile_table__setValue_string(startable, EMDL_CTF_IMAGE, 'micrographs/' // trim(adjustl(basename(getstring))))
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                endif

                if(dfx .AND. dfy) then
                    call starfile_table__setValue_double(startable, EMDL_CTF_ASTIGMATISM, real(abs((spproj%os_mic%get(i, 'dfx') * 10000) - (spproj%os_mic%get(i, 'dfy') * 10000)), dp))
                    call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUSU, real(spproj%os_mic%get(i, 'dfx') * 10000, dp))
                    call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUSV, real(spproj%os_mic%get(i, 'dfy') * 10000, dp))
                endif

                if(angast) call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_mic%get(i, 'angast'), dp))

                if(opticsgroup) call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'opticsgroup')))

            endif
        end do

        call starfile_table__open_ofile(startable, 'micrographs.star', 1)
        call starfile_table__write_ofile(startable)
        call starfile_table__close_ofile(startable)

        call starfile_table__delete(startable)

        if(allocated(getstring))deallocate(getstring)

    end subroutine write_micrographs_star

    subroutine write_particles2D_star(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        character(len=:),       allocatable     :: getstring
        character(len=1024)                     :: groupname
        integer,                allocatable     :: ptclcount(:)
        integer                                 :: group
        type(starfile_table_type)               :: startable
        integer                                 :: i, j, stkindex
        logical                                 :: exists
        logical                                 :: state        = .FALSE.
        logical                                 :: stkstate     = .FALSE.
        logical                                 :: stkind       = .FALSE.
        logical                                 :: smpd         = .FALSE.
        logical                                 :: kv           = .FALSE.
        logical                                 :: cs           = .FALSE.
        logical                                 :: dfx          = .FALSE.
        logical                                 :: dfy          = .FALSE.
        logical                                 :: angast       = .FALSE.
        logical                                 :: fraca        = .FALSE.
        logical                                 :: ptcldfx      = .FALSE.
        logical                                 :: ptcldfy      = .FALSE.
        logical                                 :: ptclangast   = .FALSE.
        logical                                 :: xpos         = .FALSE.
        logical                                 :: ypos         = .FALSE.
        logical                                 :: stk          = .FALSE.
        logical                                 :: box          = .FALSE.
        logical                                 :: opticsgroup  = .FALSE.
        real                                    :: dfxmin
        real                                    :: dfxmax
        real                                    :: dfxstep

        if(spproj%os_ptcl2D%get_noris() == 0) then
            return
        endif

        write(logfhandle, *) 'Generating particles2D.star ... '

        call simple_mkdir('particles', errmsg= "simple_relion:: create particles directory")

        if(.NOT. allocated(ptclcount)) then
            allocate(ptclcount(spproj%os_stk%get_noris()))
            do i=1, spproj%os_stk%get_noris()
               ptclcount(i) = 0
            end do
        end if

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

        stkind      = spproj%os_ptcl2D%isthere(i, 'stkind')
        ptcldfx     = spproj%os_ptcl2D%isthere(i, 'dfx')
        ptcldfy     = spproj%os_ptcl2D%isthere(i, 'dfy')
        ptclangast  = spproj%os_ptcl2D%isthere(i, 'angast')
        xpos        = spproj%os_ptcl2D%isthere(i, 'xpos')
        ypos        = spproj%os_ptcl2D%isthere(i, 'ypos')

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

        dfx         = spproj%os_stk%isthere(i, 'dfx')
        dfy         = spproj%os_stk%isthere(i, 'dfy')
        angast      = spproj%os_stk%isthere(i, 'angast')
        box         = spproj%os_stk%isthere(i, 'box')
        smpd        = spproj%os_stk%isthere(i, 'smpd')
        kv          = spproj%os_stk%isthere(i, 'kv')
        cs          = spproj%os_stk%isthere(i, 'cs')
        fraca       = spproj%os_stk%isthere(i, 'fraca')
        stk         = spproj%os_stk%isthere(i, 'stk')
        opticsgroup = spproj%os_stk%isthere(i, 'opticsgroup')

        call starfile_table__new(startable)
        call starfile_table__setName(startable, "optics")

        do j=1, self%opticsgroups
            call starfile_table__addObject(startable)
            call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, j)
            call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int(j, 4))))
            if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_stk%get(i, 'fraca'), dp))
            if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_stk%get(i, 'cs'), dp))
            if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_stk%get(i, 'kv'), dp))
            if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_stk%get(i, 'smpd'), dp))
            if(box) call starfile_table__setValue_int(startable, EMDL_IMAGE_SIZE, int(spproj%os_stk%get(i, 'box')))
            call starfile_table__setValue_int(startable, EMDL_IMAGE_DIMENSIONALITY, 2)
        end do

        call starfile_table__open_ofile(startable, 'particles2D.star', 0)
        call starfile_table__write_ofile(startable)
        call starfile_table__close_ofile(startable)
        call starfile_table__delete(startable)
        call starfile_table__new(startable)
        call starfile_table__setName(startable, "particles")

        if(cline%get_rarg('reliongroups') > 0 .AND. dfx) then
            call spproj%os_stk%minmax('dfx', dfxmin, dfxmax)
            dfxstep = (dfxmax - dfxmin)/ (cline%get_rarg('reliongroups') + 1)
        endif

        !STAR data
        do i=1, spproj%os_ptcl2D%get_noris()
            if(stkind) then
                stkindex = int(spproj%os_ptcl2D%get(i, 'stkind'))
                ptclcount(stkindex) = ptclcount(stkindex) + 1
                if((state .AND. (spproj%os_ptcl2D%get(i,'state') .GT. 0)) .OR. (.NOT. state)) then
                    if((stkstate .AND. (spproj%os_stk%get(stkindex, 'state') .GT. 0)) .OR. (.NOT. stkstate)) then
                        call starfile_table__addObject(startable)

                        if(stk) then
                            call spproj%os_stk%getter(stkindex, 'stk', getstring)
                            call starfile_table__setValue_string(startable, EMDL_IMAGE_NAME,trim(int2str(int(ptclcount(stkindex), 4))) // '@particles/' // trim(adjustl(basename(getstring))) // 's')
                            inquire(file='particles/' // trim(adjustl(basename(getstring))) // 's', exist=exists)
                            if(.NOT. exists) then
                                call syslib_symlink('../' // trim(adjustl(getstring)), 'particles/' // trim(adjustl(basename(getstring))) // 's', 'Failed to generate symlink')
                            endif
                            getstring = trim(adjustl(basename(getstring)))
                            call starfile_table__setValue_string(startable, EMDL_MICROGRAPH_NAME, 'micrographs/' // getstring(12 : len(getstring)))
                        endif

                        if(xpos .AND. ypos .AND. box) then
                            call starfile_table__setValue_double(startable, EMDL_IMAGE_COORD_X, real(int(spproj%os_ptcl2D%get(i,'xpos') + (spproj%os_stk%get(stkindex,'box') / 2)), dp))
                            call starfile_table__setValue_double(startable, EMDL_IMAGE_COORD_Y, real(int(spproj%os_ptcl2D%get(i,'ypos') + (spproj%os_stk%get(stkindex,'box') / 2)), dp))
                        endif

                        if(ptcldfx .AND. ptcldfy .AND. ptclangast) then
                            call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUSU, real(spproj%os_ptcl2D%get(i,'dfx') * 10000, dp))
                            call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUSV, real(spproj%os_ptcl2D%get(i,'dfy') * 10000, dp))
                            call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_ptcl2D%get(i,'angast'), dp))
                        else if(dfx .AND. dfy .AND. angast) then
                            call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUSU, real(spproj%os_stk%get(stkindex,'dfx') * 10000, dp))
                            call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUSV, real(spproj%os_stk%get(stkindex,'dfy') * 10000, dp))
                            call starfile_table__setValue_double(startable, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_stk%get(stkindex,'angast'), dp))
                        endif

                        if(cline%get_rarg('reliongroups') > 0 .AND. dfx) then
                            group = ceiling((real(spproj%os_stk%get(stkindex,'dfx')) - dfxmin) / dfxstep)
                            write(groupname, *) group
                            call starfile_table__setValue_string(startable, EMDL_MLMODEL_GROUP_NAME, trim(adjustl(groupname)))
                        endif

                        if(opticsgroup) call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_stk%get(stkindex, 'opticsgroup')))

                    endif
                endif
            endif
        end do

        call starfile_table__open_ofile(startable, 'particles2D.star', 1)
        call starfile_table__write_ofile(startable)
        call starfile_table__close_ofile(startable)

        call starfile_table__delete(startable)

        if(allocated(getstring))deallocate(getstring)

    end subroutine write_particles2D_star

    subroutine find_movienames(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        integer                                 :: moviecount, i
        logical                                 :: micsource
        character (len=:),      allocatable     :: moviename, tmpname

        if (spproj%os_mic%get_noris() .gt. 0) then
            moviecount = spproj%os_mic%get_noris()
            micsource = .TRUE.
        else if (spproj%os_stk%get_noris() .gt. 0) then
            moviecount = spproj%os_stk%get_noris()
            micsource = .FALSE.
        else
            THROW_HARD('no micrographs or stacks in project file')
        endif

        if(.NOT. allocated(self%movienames)) then
            allocate(self%movienames(moviecount))
        end if

        if(.NOT. allocated(self%moviegroup)) then
            allocate(self%moviegroup(moviecount))
        end if

        do i=1, moviecount
            if(micsource) then
                moviename = trim(adjustl(basename(spproj%os_mic%get_static(i, 'intg'))))
                tmpname   = moviename(:len_trim(moviename)-9)
            else
                moviename = trim(adjustl(basename(spproj%os_stk%get_static(i, 'stk'))))
                tmpname   = moviename(12:len_trim(moviename)-9)
            endif
            self%movienames(i) = tmpname
        end do

    end subroutine find_movienames

    subroutine generate_epu_tiltgroups(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        integer                                 :: i,j,k
        character (len=:),      allocatable     :: tiltname, tmpname
        character(len=20),      allocatable     :: tiltgroups(:)

        if(.NOT. allocated(tiltgroups)) then
            allocate(tiltgroups(50))            ! Max 50 raw groups set here(ie per stage movement)
        end if

        do i=1, size(tiltgroups)
            tiltgroups(i) = ''
        end do

        j = 1

        do i=1, size(self%movienames)
            tiltname = trim(adjustl(self%movienames(i)))
            tmpname  = tiltname(index(tiltname,'Data_')+5:)
            tiltname = tmpname(:index(tmpname,'_')-1)

            if(.NOT. any(tiltgroups .eq. tiltname)) then
                tiltgroups(j) = tiltname
                self%moviegroup(i) = j
                j = j+1
            else
                do k=1, size(tiltgroups)
                    if(tiltgroups(k) .eq. tiltname) then
                        exit
                    endif
                end do
                self%moviegroup(i) = k
            endif
        end do

        self%opticsgroups = j - 1

        if(allocated(tiltname))deallocate(tiltname)
        if(allocated(tiltgroups))deallocate(tiltgroups)

    end subroutine generate_epu_tiltgroups

    subroutine generate_single_tiltgroup(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        integer                                 :: i

        do i=1, size(self%movienames)
            self%moviegroup(i) = 1
        end do

        self%opticsgroups = 1

    end subroutine generate_single_tiltgroup

    subroutine generate_xml_tiltgroups(self, cline, spproj)

        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        integer                                 :: i,j,k
        character (len=:),      allocatable     :: tiltname
        type(Node),             pointer         :: xmldoc, beamtiltnode, beamtiltnodex, beamtiltnodey
        real                                    :: beamtiltx, beamtilty
        real, allocatable                       :: tilts(:,:)
        real, allocatable                       :: clustercentres(:,:)
        real, allocatable                       :: work(:)
        logical                                 :: exists
        integer                                 :: fault

        write(logfhandle, *) "Parsing movie XML files ... "

        if(cline%get_carg('xmlloc') .eq. '') then
            THROW_HARD('xmlloc is not set')
        endif

        if(cline%get_rarg('tiltcount') .le. 0) then
            THROW_HARD('tiltcount is not set')
        endif

        inquire(file=trim(adjustl(cline%get_carg('xmlloc'))), exist=exists)

        if(.NOT. exists) then
            THROW_HARD('xmlloc does not exist')
        endif

        if(.NOT. allocated(tilts)) then
            allocate(tilts(2, size(self%movienames)))
        end if

        if(.NOT. allocated(clustercentres)) then
            allocate(clustercentres(2, int(cline%get_rarg('tiltcount'))))
        end if

        if(.NOT. allocated(work)) then
            allocate(work(size(self%movienames)))
        end if

        do i=1, size(self%movienames)

            inquire(file=trim(adjustl(cline%get_carg('xmlloc'))) // "/" // int2str(i) //".xml", exist=exists)
            if(.NOT. exists) then
                THROW_HARD(trim(adjustl(cline%get_carg('xmlloc'))) // "/" // int2str(i) //".xml" // ' does not exist')
            endif

            xmldoc => parseFile(trim(adjustl(cline%get_carg('xmlloc'))) // "/" // int2str(i) //".xml")
            beamtiltnode => item(getElementsByTagname(xmldoc, "BeamTilt"), 0)
            beamtiltnodex => item(getElementsByTagname(beamtiltnode, "a:_x"), 0)
            beamtiltnodey => item(getElementsByTagname(beamtiltnode, "a:_y"), 0)
            beamtiltx = str2real(getTextContent(beamtiltnodex))
            beamtilty = str2real(getTextContent(beamtiltnodey))
             write(logfhandle, *) i, size(self%movienames), beamtiltx, beamtilty
            tilts(1, i) = beamtiltx
            tilts(2, i) = beamtilty

            call destroy(xmldoc)
        end do

        call KMPP(tilts, 2, spproj%os_stk%get_noris(), int(cline%get_rarg('tiltcount')), clustercentres, self%moviegroup, work, fault)

        self%opticsgroups = int(cline%get_rarg('tiltcount'))

        if(allocated(tiltname))deallocate(tiltname)
        if(allocated(tilts))deallocate(tilts)
        if(allocated(clustercentres))deallocate(clustercentres)
        if(allocated(work))deallocate(work)

    end subroutine generate_xml_tiltgroups

    SUBROUTINE KMPP (X, P, N, K, C, Z, WORK, IFAULT)

       IMPLICIT NONE
       INTEGER P, N, K, Z, IFAULT
       REAL X, C, WORK
       DIMENSION X(P,N), C(P,K), Z(N), WORK(N)
       INTEGER ITER                 ! maximum iterations
       REAL BIG                     ! arbitrary large number
       PARAMETER (ITER = 1000, BIG = 1E33)

     INTEGER         H          ! count iterations
     INTEGER         I          ! count points
     INTEGER         I1         ! point marked as initial center
     INTEGER         J          ! count dimensions
     INTEGER         L          ! count clusters
     INTEGER         L0         ! present cluster ID
     INTEGER         L1          ! new cluster ID

     REAL      BEST                 ! shortest distance to a center
     REAL      D2                   ! squared distance
     REAL      TOT                  ! a total
     REAL      W                     ! temp scalar

       LOGICAL CHANGE             ! whether any points have been reassigned


       IFAULT = 0
       IF (K < 1 .OR. K > N) THEN       ! K out of bounds
         IFAULT = 3
         RETURN
       END IF
       DO I = 1, N                       ! clear Z
         Z(I) = 0
       END DO


       DO I = 1, N
         WORK(I) = BIG
       END DO

       CALL RANDOM_NUMBER (W)
       I1 = MIN(INT(W * FLOAT(N)) + 1, N)  ! choose first center at random
       DO J = 1, P
         C(J,1) = X(J,I1)
       END DO

       DO L = 2, K                    ! initialize other centers
         TOT = 0.
         DO I = 1, N                     ! measure from each point
           BEST = WORK(I)
           D2 = 0.                         ! to prior center
           DO J = 1, P
             D2 = D2 + (X(J,I) - C(J,L-1)) **2  ! Squared Euclidean distance
             IF (D2 .GE. BEST) GO TO 10               ! needless to add to D2
           END DO                          ! next J
           IF (D2 < BEST) BEST = D2          ! shortest squared distance
           WORK(I) = BEST
     10      TOT = TOT + BEST             ! cumulative squared distance
         END DO                      ! next data point


         CALL RANDOM_NUMBER (W)
         W = W * TOT    ! uniform at random over cumulative distance
         TOT = 0.
         DO I = 1, N
           I1 = I
           TOT = TOT + WORK(I)
           IF (TOT > W) GO TO 20
         END DO                ! next I
     20  CONTINUE
         DO J = 1, P         ! assign center
           C(J,L) = X(J,I1)
         END DO
       END DO               ! next center to initialize

       DO H = 1, ITER
         CHANGE = .FALSE.

         DO I = 1, N
           L0 = Z(I)
           L1 = 0
           BEST = BIG
           DO L = 1, K
             D2 = 0.
             DO J = 1, P
               D2 = D2 + (X(J,I) - C(J,L)) **2
               IF (D2 .GE. BEST) GO TO 30
             END DO
     30      CONTINUE
             IF (D2 < BEST) THEN           ! new nearest center
               BEST = D2
               L1 = L
             END IF
           END DO        ! next L

           IF (L0 .NE. L1) THEN
             Z(I) = L1                   !  reassign point
             CHANGE = .TRUE.
           END IF
         END DO         ! next I
         IF (.NOT. CHANGE) RETURN      ! success

         DO L = 1, K              ! zero population
           WORK(L) = 0.
         END DO
         DO L = 1, K               ! zero centers
           DO J = 1, P
             C(J,L) = 0.
           END DO
         END DO

         DO I = 1, N
           L = Z(I)
           WORK(L) = WORK(L) + 1.             ! count
           DO J = 1, P
             C(J,L) = C(J,L) + X(J,I)         ! add
           END DO
         END DO

         DO L = 1, K
           IF (WORK(L) < 0.5) THEN          ! empty cluster check
             IFAULT = 1                     ! fatal error
             RETURN
           END IF
           W = 1. / WORK(L)
           DO J = 1, P
             C(J,L) = C(J,L) * W     ! multiplication is faster than division
           END DO
         END DO

       END DO                   ! next H
       IFAULT = 2                ! too many iterations
       RETURN

      END SUBROUTINE KMPP ! of KMPP

    subroutine allocate_opticsgroups(self, cline, spproj)

        class(relion_project),  intent(inout) :: self
        class(sp_project),      intent(inout) :: spproj
        class(cmdline),         intent(inout) :: cline
        character (len=:),      allocatable   :: moviename
        integer                               :: i, j
        integer, allocatable                  :: opticsgroupmap(:,:)
        real                                  :: tiltgroupmax

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
            moviename = trim(adjustl(basename(spproj%os_mic%get_static(i, 'intg'))))
            moviename = moviename(:len_trim(moviename)-9)
            do j=1, size(self%movienames)
                if(self%movienames(j) .eq. moviename) then
                    exit
                end if
            end do
            call spproj%os_mic%set(i, 'opticsgroup', float(self%moviegroup(j)))
        end do

        do i=1, spproj%os_stk%get_noris()
            moviename = trim(adjustl(basename(spproj%os_stk%get_static(i, 'stk'))))
            moviename = moviename(12:len_trim(moviename)-9)
            do j=1, size(self%movienames)
                if(self%movienames(j) .eq. moviename) then
                    exit
                end if
            end do
            call spproj%os_stk%set(i, 'opticsgroup', float(self%moviegroup(j)))
        end do

        if(allocated(moviename)) deallocate(moviename)

    end subroutine allocate_opticsgroups

    subroutine create(self, cline, spproj)

        class(relion_project),  intent(inout) :: self
        class(sp_project),      intent(inout) :: spproj
        class(cmdline),         intent(inout) :: cline

        write(logfhandle, *) 'Writing Relion 3.1 compatible STAR files'

        call self%find_movienames(cline, spproj)

        if(cline%get_carg('tiltgroups') .eq. 'epu') then
            call self%generate_epu_tiltgroups(cline, spproj)
        else if(cline%get_carg('tiltgroups') .eq. 'xml') then
            call self%generate_xml_tiltgroups(cline, spproj)
        else
            call self%generate_single_tiltgroup(cline, spproj)
        endif

        call self%allocate_opticsgroups(cline, spproj)
        call self%write_corrected_micrographs_star(cline, spproj)
        call self%write_micrographs_star(cline, spproj)
        call self%write_particles2D_star(cline, spproj)

        if(allocated(self%movienames)) deallocate(self%movienames)
        if(allocated(self%moviegroup)) deallocate(self%moviegroup)

    end subroutine create


    ! distance threshold based yerarchical clustering
    ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
    subroutine h_clust(self,data_in,thresh,labels,centroids,populations)
        class(relion_project),  intent(inout) :: self ! unused
        real,    intent(in)  :: data_in(:,:)   ! input data, point coords
        real,    intent(in)  :: thresh         ! threshold for class merging
        integer, intent(out) :: labels(:)      ! labels of the elements in vec
        real,    allocatable, intent(out) :: centroids(:,:) ! centroids of the classes
        integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
        real,    allocatable :: mat(:,:)                    ! pariwise distances matrix
        logical, allocatable :: mask(:), outliers(:)
        integer :: N, i, j, cnt, ncls, filnum, io_stat
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
        write(logfhandle, *) 'HIERARCHICAL CLUSTERING COMPLETED'
        write(logfhandle, *) 'Distance threshold selected: ', thresh
        do i = 1, N
            write(logfhandle, *) 'Point [', data_in(i,1),',', data_in(i,2), '], class: ', labels(i)
        enddo
        do i = 1, ncls
            write(logfhandle, *) 'Class ', i, 'Population ', populations(i), ' Centroid [', centroids(i,1), ',', centroids(i,2), ']'
        enddo
    end subroutine h_clust

end module simple_relion

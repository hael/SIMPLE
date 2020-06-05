module simple_relion
include 'simple_lib.f08'
use simple_starfile_wrappers
use simple_sp_project,          only: sp_project
use simple_cmdline,             only: cmdline
implicit none
private
public :: relion_project
#include "simple_local_flags.inc"

type relion_project
    character(len=20),      allocatable     :: tiltgroups(:)


contains

    procedure :: create
    procedure :: write_corrected_micrographs_star
    procedure :: write_micrographs_star
    procedure :: write_particles2D_star
    procedure :: generate_epu_tiltgroups
    
    !procedure :: write_particles3D_star

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
        logical                                 :: tiltgroup    = .FALSE.
        
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
        tiltgroup   = spproj%os_mic%isthere(i, 'tiltgroup')
        
        if(cline%get_carg('relion31') .eq. 'yes') then
            call starfile_table__new(startable)
            call starfile_table__setName(startable, "optics")
            if(tiltgroup) then
                do j=1, size(self%tiltgroups)
                    if(self%tiltgroups(j) .ne. '') then
                        call starfile_table__addObject(startable)
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, j)
                        call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int(j, 4))))
                        if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                        if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                        if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                        if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
                    endif
                  
                end do
            else
                call starfile_table__addObject(startable)
                call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, 1)
                call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup1')
                if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
            endif
            call starfile_table__open_ofile(startable, 'micrographs_corrected.star', 0)
            call starfile_table__write_ofile(startable)
            call starfile_table__close_ofile(startable)
            call starfile_table__delete(startable)
            call starfile_table__new(startable)
            call starfile_table__setName(startable, "micrographs")
        else
            call starfile_table__new(startable)
        endif
        
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
                
                if(cline%get_carg('relion31') .eq. 'no') then
                    if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                    if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                    if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                    if(smpd) then
                        call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
                        call starfile_table__setValue_double(startable, EMDL_CTF_MAGNIFICATION, real(10000, dp))
                    endif
                else
                    if(tiltgroup) then
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'tiltgroup')))
                    else
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, 1)
                    endif
                endif
                
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
        logical                                 :: tiltgroup    = .FALSE.
        
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
        tiltgroup   = spproj%os_mic%isthere(i, 'tiltgroup')
        
        
        if(cline%get_carg('relion31') .eq. 'yes') then
            call starfile_table__new(startable)
            call starfile_table__setName(startable, "optics")
            if(tiltgroup) then
                do j=1, size(self%tiltgroups)
                    if(self%tiltgroups(j) .ne. '') then
                        call starfile_table__addObject(startable)
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, j)
                        call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int(j, 4))))
                        if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                        if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                        if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                        if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
                    endif
                  
                end do
            else
                call starfile_table__addObject(startable)
                call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, 1)
                call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup1')
                if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
            endif
            call starfile_table__open_ofile(startable, 'micrographs.star', 0)
            call starfile_table__write_ofile(startable)
            call starfile_table__close_ofile(startable)
            call starfile_table__delete(startable)
            call starfile_table__new(startable)
            call starfile_table__setName(startable, "micrographs")
        else
            call starfile_table__new(startable)
        endif
        
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
                
                if(cline%get_carg('relion31') .eq. 'no') then
                    if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                    if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                    if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                    if(smpd) then
                        call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
                        call starfile_table__setValue_double(startable, EMDL_CTF_MAGNIFICATION, real(10000, dp))
                    endif
                else
                    if(tiltgroup) then
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'tiltgroup')))
                    else
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, 1)
                    endif
                endif
                
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
        logical                                 :: tiltgroup    = .FALSE.
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
        tiltgroup   = spproj%os_stk%isthere(i, 'tiltgroup')
        
        if(cline%get_carg('relion31') .eq. 'yes') then
            call starfile_table__new(startable)
            call starfile_table__setName(startable, "optics")
            if(tiltgroup) then
                do j=1, size(self%tiltgroups)
                    if(self%tiltgroups(j) .ne. '') then
                        call starfile_table__addObject(startable)
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, j)
                        call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup' // trim(int2str(int(j, 4))))
                        if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_stk%get(i, 'fraca'), dp))
                        if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_stk%get(i, 'cs'), dp))
                        if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_stk%get(i, 'kv'), dp))
                        if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_stk%get(i, 'smpd'), dp))
                        if(box) call starfile_table__setValue_int(startable, EMDL_IMAGE_SIZE, int(spproj%os_stk%get(i, 'box')))
                        call starfile_table__setValue_int(startable, EMDL_IMAGE_DIMENSIONALITY, 2)
                        
                    endif
                  
                end do
            else
                call starfile_table__addObject(startable)
                call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, 1)
                call starfile_table__setValue_string(startable, EMDL_IMAGE_OPTICS_GROUP_NAME, 'opticsgroup1')
                if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_stk%get(i, 'fraca'), dp))
                if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_stk%get(i, 'cs'), dp))
                if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_stk%get(i, 'kv'), dp))
                if(smpd) call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_stk%get(i, 'smpd'), dp))
            endif
            call starfile_table__open_ofile(startable, 'particles2D.star', 0)
            call starfile_table__write_ofile(startable)
            call starfile_table__close_ofile(startable)
            call starfile_table__delete(startable)
            call starfile_table__new(startable)
            call starfile_table__setName(startable, "particles")
        else
            call starfile_table__new(startable)
        endif
        
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
                        
                        if(cline%get_carg('relion31') .eq. 'no') then
                            if(fraca) call starfile_table__setValue_double(startable, EMDL_CTF_Q0, real(spproj%os_stk%get(stkindex, 'fraca'), dp))
                            if(cs) call starfile_table__setValue_double(startable, EMDL_CTF_CS, real(spproj%os_stk%get(stkindex, 'cs'), dp))
                            if(kv) call starfile_table__setValue_double(startable, EMDL_CTF_VOLTAGE, real(spproj%os_stk%get(stkindex, 'kv'), dp))
                            if(smpd) then
                                call starfile_table__setValue_double(startable, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_stk%get(stkindex, 'smpd'), dp))
                                call starfile_table__setValue_double(startable, EMDL_CTF_MAGNIFICATION, real(10000, dp))
                            endif
                            if(tiltgroup) call starfile_table__setValue_int(startable, EMDL_PARTICLE_BEAM_TILT_CLASS, int(spproj%os_stk%get(stkindex, 'tiltgroup')))
                        else
                            if(tiltgroup) then
                                call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_stk%get(stkindex, 'tiltgroup')))
                            else
                                call starfile_table__setValue_int(startable, EMDL_IMAGE_OPTICS_GROUP, 1)
                            endif
                        endif
                                
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
    
    subroutine generate_epu_tiltgroups(self, cline, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        integer                                 :: i,j,k
        character (len=:),      allocatable     :: tiltname

        if(.NOT. allocated(self%tiltgroups)) then
            allocate(self%tiltgroups(50))            ! Max 50 groups set here
        end if
        
        do i=1, size(self%tiltgroups) 
            self%tiltgroups(i) = ''
        end do
        
        j = 1
        
        do i=1, spproj%os_stk%get_noris()
            tiltname = trim(adjustl(basename(spproj%os_stk%get_static(i, 'stk'))))
            tiltname = tiltname(index(tiltname,'Data_') + 5:)
            tiltname = tiltname(:index(tiltname,'_')-1)
            if(.NOT. any(self%tiltgroups .eq. tiltname)) then
                self%tiltgroups(j) = tiltname
                call spproj%os_stk%set(i, 'tiltgroup', float(j))
                j = j+1
            else
                do k=1, size(self%tiltgroups)
                    if(self%tiltgroups(k) .eq. tiltname) then
                        exit
                    endif
                end do
                call spproj%os_stk%set(i, 'tiltgroup', float(k))
            endif
        end do
       
        do i=1, spproj%os_mic%get_noris()
            tiltname = trim(adjustl(basename(spproj%os_mic%get_static(i, 'intg'))))
            tiltname = tiltname(index(tiltname,'Data_') + 5:)
            tiltname = tiltname(:index(tiltname,'_')-1)
            if(.NOT. any(self%tiltgroups .eq. tiltname)) then
                self%tiltgroups(j) = tiltname
                call spproj%os_mic%set(i, 'tiltgroup', float(j))
                j = j+1
            else
                do k=1, size(self%tiltgroups)
                    if(self%tiltgroups(k) .eq. tiltname) then
                        exit
                    endif
                end do
                call spproj%os_mic%set(i, 'tiltgroup', float(k))
            endif
        end do
        
        if(allocated(tiltname))deallocate(tiltname)
        
    end subroutine generate_epu_tiltgroups

    subroutine create(self, cline, spproj)

        class(relion_project),  intent(inout) :: self
        class(sp_project),      intent(inout) :: spproj
        class(cmdline),         intent(inout) :: cline

        if(cline%get_carg('relion31') .eq. 'yes') then
            write(logfhandle, *) 'Writing Relion 3.1 compatible STAR files'
        else
            write(logfhandle, *) 'Writing Relion 3.0 compatible STAR files'
        endif
        
        if(cline%get_carg('tiltgroups') .eq. 'epu') then
            call self%generate_epu_tiltgroups(cline, spproj)
        endif
        
        call self%write_corrected_micrographs_star(cline, spproj)
        call self%write_micrographs_star(cline, spproj)
        call self%write_particles2D_star(cline, spproj)
       ! call self%write_particles3D_star(spproj) !Disabled as not fully tested
        if(allocated(self%tiltgroups))deallocate(self%tiltgroups)
    end subroutine create

end module simple_relion

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

    character(len=8),   allocatable     :: mickey(:)
    logical,            allocatable     :: miclogical(:)
    character(len=8),   allocatable     :: stkkey(:)
    logical,            allocatable     :: stklogical(:)
    character(len=8),   allocatable     :: ptcl2Dkey(:)
    logical,            allocatable     :: ptcl2Dlogical(:)
    character(len=8),   allocatable     :: ptcl3Dkey(:)
    logical,            allocatable     :: ptcl3Dlogical(:)
    real                                :: dfxmin
    real                                :: dfxmax
    real                                :: dfxstep

contains

    procedure :: create
    procedure :: test_mic_params
    procedure :: test_stk_params
    procedure :: test_ptcl2D_params
    procedure :: test_ptcl3D_params
    procedure :: write_corrected_micrographs_star
    procedure :: write_micrographs_star
    procedure :: write_particles2D_star
    procedure :: get_epu_tiltgroup
    !procedure :: write_particles3D_star

end type relion_project

contains

    subroutine test_mic_params(self, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        integer                                 :: i
        integer                                 :: j
        
        if(.NOT. allocated(self%mickey)) then
            allocate(self%mickey(12))
        end if
        
        if(.NOT. allocated(self%miclogical)) then
            allocate(self%miclogical(size(self%mickey)))
        end if
        
        self%mickey(1) = 'movie'
        self%miclogical(1) = .TRUE.
        self%mickey(2) = 'intg'
        self%miclogical(2) = .TRUE.
        self%mickey(3) = 'forctf'
        self%miclogical(3) = .TRUE.
        self%mickey(4) = 'pspec'
        self%miclogical(4) = .TRUE.
        self%mickey(5) = 'smpd'
        self%miclogical(5) = .TRUE.
        self%mickey(6) = 'kv'
        self%miclogical(6) = .TRUE.
        self%mickey(7) = 'cs'
        self%miclogical(7) = .TRUE.
        self%mickey(8) = 'dfx'
        self%miclogical(8) = .TRUE.
        self%mickey(9) = 'dfy'
        self%miclogical(9) = .TRUE.
        self%mickey(10) = 'angast'
        self%miclogical(10) = .TRUE.
        self%mickey(11) = 'fraca'
        self%miclogical(11) = .TRUE.
        self%mickey(12) = 'state'
        self%miclogical(12) = .TRUE.

        do i=1, spproj%os_mic%get_noris()
            if(spproj%os_mic%isthere(i, 'state')) then
                if(spproj%os_mic%get(i,'state') .GT. 0) then
                    do j=1, size(self%mickey)
                        if(.NOT. spproj%os_mic%isthere(i,trim(self%mickey(j)))) then
                            self%miclogical(j) = .FALSE.
                        endif
                    end do
                endif
            else
                do j=1, size(self%mickey)
                    if(.NOT. spproj%os_mic%isthere(i,trim(self%mickey(j)))) then
                        self%miclogical(j) = .FALSE.
                    endif
                end do
            endif
        end do
        
    end subroutine test_mic_params
        
    subroutine test_stk_params(self, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        integer                                 :: i
        integer                                 :: j
        
        if(.NOT. allocated(self%stkkey)) then
            allocate(self%stkkey(9))
        end if
        
        if(.NOT. allocated(self%stklogical)) then
            allocate(self%stklogical(size(self%stkkey)))
        end if
        
        self%stkkey(1) = 'smpd'
        self%stklogical(1) = .TRUE.
        self%stkkey(2) = 'kv'
        self%stklogical(2) = .TRUE.
        self%stkkey(3) = 'cs'
        self%stklogical(3) = .TRUE.
        self%stkkey(4) = 'dfx'
        self%stklogical(4) = .TRUE.
        self%stkkey(5) = 'dfy'
        self%stklogical(5) = .TRUE.
        self%stkkey(6) = 'angast'
        self%stklogical(6) = .TRUE.
        self%stkkey(7) = 'fraca'
        self%stklogical(7) = .TRUE.
        self%stkkey(8) = 'box'
        self%stklogical(8) = .TRUE.
        self%stkkey(9) = 'stk'
        self%stklogical(9) = .TRUE.
        
        do i=1, spproj%os_stk%get_noris()
            if(spproj%os_stk%isthere(i, 'state')) then
                if(spproj%os_stk%get(i,'state') .GT. 0) then
                    do j=1, size(self%stkkey)
                        if(.NOT. spproj%os_stk%isthere(i,trim(self%stkkey(j)))) then
                            self%stklogical(j) = .FALSE.
                        endif
                    end do
                endif
            else
                do j=1, size(self%stkkey)
                    if(.NOT. spproj%os_stk%isthere(i,trim(self%stkkey(j)))) then
                        self%stklogical(j) = .FALSE.
                    endif
                end do
            endif
        end do
        
    end subroutine test_stk_params
    
    subroutine test_ptcl2D_params(self, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        integer                                 :: i
        integer                                 :: j
        
        if(.NOT. allocated(self%ptcl2Dkey)) then
            allocate(self%ptcl2Dkey(6))
        end if
        
        if(.NOT. allocated(self%ptcl2Dlogical)) then
            allocate(self%ptcl2Dlogical(size(self%ptcl2Dkey)))
        end if
        
        self%ptcl2Dkey(1) = 'xpos'
        self%ptcl2Dlogical(1) = .TRUE.
        self%ptcl2Dkey(2) = 'ypos'
        self%ptcl2Dlogical(2) = .TRUE.
        self%ptcl2Dkey(3) = 'stkind'
        self%ptcl2Dlogical(3) = .TRUE.
        self%ptcl2Dkey(4) = 'state'
        self%ptcl2Dlogical(4) = .TRUE.
        self%ptcl2Dkey(5) = 'dfx'
        self%ptcl2Dlogical(5) = .TRUE.
        self%ptcl2Dkey(6) = 'dfy'
        self%ptcl2Dlogical(6) = .TRUE.
        
        do i=1, spproj%os_ptcl2D%get_noris()
            if(spproj%os_ptcl2D%isthere(i, 'state')) then
                if(spproj%os_ptcl2D%get(i,'state') .GT. 0) then
                    do j=1, size(self%ptcl2Dkey)
                        if(.NOT. spproj%os_ptcl2D%isthere(i,trim(self%ptcl2Dkey(j)))) then
                            self%ptcl2Dlogical(j) = .FALSE.
                        endif
                    end do
                endif
            else
                do j=1, size(self%ptcl2Dkey)
                    if(.NOT. spproj%os_ptcl2D%isthere(i,trim(self%ptcl2Dkey(j)))) then
                        self%ptcl2Dlogical(j) = .FALSE.
                    endif
                end do
            endif
        end do

    end subroutine test_ptcl2D_params
    
    subroutine test_ptcl3D_params(self, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        integer                                 :: i
        integer                                 :: j
        
        if(.NOT. allocated(self%ptcl3Dkey)) then
            allocate(self%ptcl3Dkey(9))
        end if
        
        if(.NOT. allocated(self%ptcl3Dlogical)) then
            allocate(self%ptcl3Dlogical(size(self%ptcl3Dkey)))
        end if
        
        self%ptcl3Dkey(1) = 'state'
        self%ptcl3Dlogical(1) = .TRUE.
        self%ptcl3Dkey(2) = 'e1'
        self%ptcl3Dlogical(2) = .TRUE.
        self%ptcl3Dkey(3) = 'e2'
        self%ptcl3Dlogical(3) = .TRUE.
        self%ptcl3Dkey(4) = 'e3'
        self%ptcl3Dlogical(4) = .TRUE.
        self%ptcl3Dkey(5) = 'x'
        self%ptcl3Dlogical(5) = .TRUE.
        self%ptcl3Dkey(6) = 'y'
        self%ptcl3Dlogical(6) = .TRUE.
        self%ptcl3Dkey(7) = 'xpos'
        self%ptcl3Dlogical(7) = .TRUE.
        self%ptcl3Dkey(8) = 'ypos'
        self%ptcl3Dlogical(8) = .TRUE.
        self%ptcl3Dkey(9) = 'stkind'
        self%ptcl3Dlogical(9) = .TRUE.
              
        do i=1, spproj%os_ptcl3D%get_noris()
            if(spproj%os_ptcl3D%isthere(i, 'state')) then
                if(spproj%os_ptcl3D%get(i,'state') .GT. 0) then
                    do j=1, size(self%ptcl3Dkey)
                        if(.NOT. spproj%os_ptcl3D%isthere(i,trim(self%ptcl3Dkey(j)))) then
                            self%ptcl3Dlogical(j) = .FALSE.
                        endif
                    end do
                endif
            else
                do j=1, size(self%ptcl3Dkey)
                    if(.NOT. spproj%os_ptcl3D%isthere(i,trim(self%ptcl3Dkey(j)))) then
                        self%ptcl3Dlogical(j) = .FALSE.
                    endif
                end do
            endif
        end do

    end subroutine test_ptcl3D_params 
    
    subroutine write_corrected_micrographs_star(self, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        character(len=:),       allocatable     :: getstring
        type(starfile_table_type)               :: mic_starfile
        integer                                 :: i
        logical                                 :: exists
        
        if(spproj%os_mic%get_noris() == 0) then
            return
        endif
        
        if(.NOT. self%miclogical(1)) then
            return
        endif
        
        if(.NOT. self%miclogical(2)) then
            return
        endif
        
        write(logfhandle, *) 'Generating micrographs_corrected.star ... '

        call simple_mkdir('micrographs', errmsg= "simple_relion:: create micrographs directory")
        
        call simple_mkdir('movies', errmsg= "simple_relion:: create movies directory")
        
        call starfile_table__new(mic_starfile)
        call starfile_table__setcomment(mic_starfile, "SIMPLE 3.0; export_relion")
        
        !STAR data
        do i=1, spproj%os_mic%get_noris()
            call starfile_table__addObject(mic_starfile)
            call spproj%os_mic%getter(i, 'intg', getstring)
            call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_NAME, 'micrographs/' // trim(adjustl(basename(getstring))))
            inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
            if(.NOT. exists) then
                call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
            endif
            call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_METADATA_NAME, 'micrographs/' // trim(adjustl(basename(getstring(1 : len(getstring) - 9)))) // '.star')
            inquire(file='micrographs/' // trim(adjustl(basename(getstring(1 : len(getstring) - 9)))) // '.star', exist=exists)
            if(.NOT. exists) then
                call syslib_symlink('../' // trim(adjustl(getstring(1 : len(getstring) - 9))) // '.star ', 'micrographs/' // trim(adjustl(basename(getstring(1 : len(getstring) - 9)))) // '.star ', 'Failed to generate symlink')
            endif
            call spproj%os_mic%getter(i, 'movie', getstring)
            inquire(file='movies/' // trim(adjustl(basename(getstring))) // ' ', exist=exists)
            if(.NOT. exists) then
                call syslib_symlink(trim(adjustl(getstring)), 'movies/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
            endif
        end do

        call starfile_table__open_ofile(mic_starfile, 'micrographs_corrected.star')
        call starfile_table__write_ofile(mic_starfile)
        call starfile_table__close_ofile(mic_starfile)

        call starfile_table__delete(mic_starfile)
        
        if(allocated(getstring))deallocate(getstring)

    end subroutine write_corrected_micrographs_star
    
    subroutine write_micrographs_star(self, cline, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        character(len=:),       allocatable     :: getstring
        type(starfile_table_type)               :: mic_starfile
        integer                                 :: i
        logical                                 :: exists
        
        if(spproj%os_mic%get_noris() == 0) then
            return
        endif
        
        if(.NOT. self%miclogical(2)) then
            return
        endif
        
        write(logfhandle, *) 'Generating micrographs.star ... '

        call simple_mkdir('micrographs', errmsg= "simple_relion:: create micrographs directory")
        
        call starfile_table__new(mic_starfile)
        call starfile_table__setcomment(mic_starfile, "SIMPLE 3.0; export_relion")
        

        !STAR data
        do i=1, spproj%os_mic%get_noris()
            if((self%miclogical(12) .AND. (spproj%os_mic%get(i,'state') .GT. 0)) .OR. (.NOT. self%miclogical(12))) then
                call starfile_table__addObject(mic_starfile)
                call spproj%os_mic%getter(i, 'intg', getstring)
                call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_NAME, 'micrographs/' // trim(adjustl(basename(getstring))))
                inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                if(.NOT. exists) then
                    call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                endif
                if(cline%get_carg('eputiltgroups') .eq. 'yes') then
                    call starfile_table__setValue_int(mic_starfile, EMDL_PARTICLE_BEAM_TILT_CLASS, self%get_epu_tiltgroup(getstring))
                endif 
                if(self%miclogical(3)) then
                    call spproj%os_mic%getter(i, 'forctf', getstring)
                    call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_NAME_WODOSE, 'micrographs/' // trim(adjustl(basename(getstring))))
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                endif
                if(self%miclogical(4)) then
                    call spproj%os_mic%getter(i, 'pspec', getstring)
                    call starfile_table__setValue_string(mic_starfile, EMDL_CTF_IMAGE, 'micrographs/' // trim(adjustl(basename(getstring))))
                    inquire(file='micrographs/' // trim(adjustl(basename(getstring))), exist=exists)
                    if(.NOT. exists) then
                        call syslib_symlink('../' // trim(adjustl(getstring)), 'micrographs/' // trim(adjustl(basename(getstring))), 'Failed to generate symlink')
                    endif
                endif
                if(self%miclogical(5)) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DETECTOR_PIXEL_SIZE, real(spproj%os_mic%get(i, 'smpd'), dp))
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_MAGNIFICATION, real(10000, dp))
                endif
                if(self%miclogical(6)) call starfile_table__setValue_double(mic_starfile, EMDL_CTF_VOLTAGE, real(spproj%os_mic%get(i, 'kv'), dp))
                if(self%miclogical(7)) call starfile_table__setValue_double(mic_starfile, EMDL_CTF_CS, real(spproj%os_mic%get(i, 'cs'), dp))
                if(self%miclogical(8) .AND. self%miclogical(9)) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_ASTIGMATISM, real(abs((spproj%os_mic%get(i, 'dfx') * 10000) - (spproj%os_mic%get(i, 'dfy') * 10000)), dp))
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DEFOCUSU, real(spproj%os_mic%get(i, 'dfx') * 10000, dp))
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DEFOCUSV, real(spproj%os_mic%get(i, 'dfy') * 10000, dp))
                endif
                if(self%miclogical(10)) call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_mic%get(i, 'angast'), dp))
                if(self%miclogical(11)) call starfile_table__setValue_double(mic_starfile, EMDL_CTF_Q0, real(spproj%os_mic%get(i, 'fraca'), dp))
                   
            endif
        end do
        
        call starfile_table__open_ofile(mic_starfile, 'micrographs.star')
        call starfile_table__write_ofile(mic_starfile)
        call starfile_table__close_ofile(mic_starfile)

        call starfile_table__delete(mic_starfile)
        
        if(allocated(getstring))deallocate(getstring)
        
    end subroutine write_micrographs_star
  
    subroutine write_particles2D_star(self, cline, spproj)
    
        class(relion_project),  intent(inout)   :: self
        class(sp_project),      intent(inout)   :: spproj
        class(cmdline),         intent(inout)   :: cline
        character(len=:),       allocatable     :: getstring
        type(starfile_table_type)               :: ptcl_starfile
        integer,                allocatable     :: ptclcount(:)
        integer                                 :: i, stkind
        integer                                 :: group
        logical                                 :: exists
        character(len=1024)                     :: groupname
        
        if(spproj%os_ptcl2D%get_noris() == 0) then
            return
        endif
        
        if(.NOT. self%ptcl2Dlogical(3)) then
            return
        endif
        
        if(.NOT. allocated(ptclcount)) then
            allocate(ptclcount(spproj%os_stk%get_noris())) 
            do i=1, spproj%os_stk%get_noris()
               ptclcount(i) = 0 
            end do
        end if
        
        write(logfhandle, *) 'Generating particles2D.star ... '

        call simple_mkdir('particles', errmsg= "simple_relion:: create micrographs directory")
        
        if(self%ptcl2Dlogical(5) .AND. self%ptcl2Dlogical(6) .AND. cline%get_carg('reliongroups') .eq. 'yes') then
             call spproj%os_ptcl2D%minmax('dfx', self%dfxmin, self%dfxmax)
             self%dfxstep = (self%dfxmax -self%dfxmin)/(cline%get_rarg('reliongroups_count') + 1)
        endif
        
        call starfile_table__new(ptcl_starfile)
        call starfile_table__setcomment(ptcl_starfile, "SIMPLE 3.0; export_relion")
        
        !STAR data
        do i=1, spproj%os_ptcl2D%get_noris()
            stkind = int(spproj%os_ptcl2D%get(i, 'stkind'))
            ptclcount(stkind) = ptclcount(stkind) + 1
            if((self%ptcl2Dlogical(4) .AND. (spproj%os_ptcl2D%get(i,'state') .GT. 0)) .OR. (.NOT. self%ptcl2Dlogical(4))) then
                call starfile_table__addObject(ptcl_starfile)
                call spproj%os_stk%getter(stkind, 'stk', getstring)
                call starfile_table__setValue_string(ptcl_starfile, EMDL_IMAGE_NAME,trim(int2str(int(ptclcount(stkind), 4))) // '@particles/' // trim(adjustl(basename(getstring))) // 's')
                inquire(file='particles/' // trim(adjustl(basename(getstring))) // 's', exist=exists)
                if(.NOT. exists) then
                    call syslib_symlink('../' // trim(adjustl(getstring)), 'particles/' // trim(adjustl(basename(getstring))) // 's', 'Failed to generate symlink')
                endif
                if(cline%get_carg('eputiltgroups') .eq. 'yes') then
                    call starfile_table__setValue_int(ptcl_starfile, EMDL_PARTICLE_BEAM_TILT_CLASS, self%get_epu_tiltgroup(getstring))
                endif 
                getstring = trim(adjustl(basename(getstring)))
                call starfile_table__setValue_string(ptcl_starfile, EMDL_MICROGRAPH_NAME, 'micrographs/' // getstring(12 : len(getstring)))
                if(self%ptcl2Dlogical(1) .AND. self%ptcl2Dlogical(2) .AND. self%stklogical(8)) then
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_IMAGE_COORD_X, real(int(spproj%os_ptcl2D%get(i,'xpos') + (spproj%os_stk%get(stkind,'box') / 2)), dp))
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_IMAGE_COORD_Y, real(int(spproj%os_ptcl2D%get(i,'ypos') + (spproj%os_stk%get(stkind,'box') / 2)), dp))
                endif
                if(self%stklogical(1)) then
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_DETECTOR_PIXEL_SIZE, real(spproj%os_stk%get(stkind,'smpd'), dp))
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_MAGNIFICATION, real(10000, dp))
                endif
                if(self%stklogical(2)) call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_VOLTAGE, real(spproj%os_stk%get(stkind,'kv'), dp))
                if(self%stklogical(3)) call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_CS, real(spproj%os_stk%get(stkind,'cs'), dp))
                
                if(self%ptcl2Dlogical(5) .AND. self%ptcl2Dlogical(6)) then
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_DEFOCUSU, real(spproj%os_ptcl2D%get(i,'dfx') * 10000, dp))
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_DEFOCUSV, real(spproj%os_ptcl2D%get(i,'dfy') * 10000, dp))
                else if(self%stklogical(4) .AND. self%stklogical(5)) then
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_DEFOCUSU, real(spproj%os_stk%get(stkind,'dfx') * 10000, dp))
                    call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_DEFOCUSV, real(spproj%os_stk%get(stkind,'dfy') * 10000, dp))
                endif
                
                if(cline%get_carg('reliongroups') .eq. 'yes') then
                    group = ceiling((real(spproj%os_stk%get(stkind,'dfx')) - self%dfxmin) / self%dfxstep)
                    write(groupname, *) group
                    call starfile_table__setValue_string(ptcl_starfile, EMDL_MLMODEL_GROUP_NAME, trim(adjustl(groupname)))
                endif
                
                if(self%stklogical(6)) call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_stk%get(stkind,'angast'), dp))
                if(self%stklogical(7)) call starfile_table__setValue_double(ptcl_starfile, EMDL_CTF_Q0, real(spproj%os_stk%get(stkind,'fraca'), dp))
            endif
        end do

        call starfile_table__open_ofile(ptcl_starfile, 'particles2D.star')
        call starfile_table__write_ofile(ptcl_starfile)
        call starfile_table__close_ofile(ptcl_starfile)

        call starfile_table__delete(ptcl_starfile)
        
        if(allocated(getstring))deallocate(getstring)
        if(allocated(ptclcount))deallocate(ptclcount)
        
    end subroutine write_particles2D_star
  
    function get_epu_tiltgroup(self, fname) result (tiltint)
    
        class(relion_project),  intent(inout)              :: self
        character (len=:),      intent(in),allocatable     :: fname
        character (len=:),      allocatable                :: tiltname
        integer                                            :: i
        integer                                            :: tiltint

        tiltint=0
        tiltname = trim(adjustl(basename(fname)))
        i = index(tiltname,'Data_')
        tiltname = tiltname(i+5:)
        i = index(tiltname,'_')
        tiltname = tiltname(:i-1)
        read ( tiltname, * ) tiltint
        
        if(allocated(tiltname))deallocate(tiltname)

    end function get_epu_tiltgroup
  
    subroutine create(self, cline, spproj)
    
        class(relion_project),  intent(inout) :: self
        class(sp_project),      intent(inout) :: spproj
        class(cmdline),         intent(inout) :: cline
        
        call self%test_mic_params(spproj)
        call self%test_stk_params(spproj)
        call self%test_ptcl2D_params(spproj)
        call self%test_ptcl3D_params(spproj)
        call self%write_corrected_micrographs_star(spproj)
        call self%write_micrographs_star(cline, spproj)
        call self%write_particles2D_star(cline, spproj)
       ! call self%write_particles3D_star(spproj) !Disabled as not fully tested
        
        if(allocated(self%mickey))deallocate(self%mickey)
        if(allocated(self%miclogical))deallocate(self%miclogical)
        if(allocated(self%stkkey))deallocate(self%stkkey)
        if(allocated(self%stklogical))deallocate(self%stklogical)
        if(allocated(self%ptcl2Dkey))deallocate(self%ptcl2Dkey)
        if(allocated(self%ptcl2Dlogical))deallocate(self%ptcl2Dlogical)
        if(allocated(self%ptcl3Dkey))deallocate(self%ptcl3Dkey)
        if(allocated(self%ptcl3Dlogical))deallocate(self%ptcl3Dlogical)
    
    end subroutine create
    
end module simple_relion

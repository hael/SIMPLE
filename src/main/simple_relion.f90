module simple_relion
include 'simple_lib.f08'
use simple_starfile_wrappers
use simple_sp_project, only: sp_project

implicit none
private
public :: relion_project
#include "simple_local_flags.inc"

type relion_project

    character(len=LONGSTRLEN), allocatable :: filestab(:,:) ! 1-movie, 2-intg, 3-forctf, 4-pspec
    character(len=LONGSTRLEN), allocatable :: stktab(:,:) ! 1-stk, 2-written
    real, allocatable :: filesparams(:,:) ! 1-smpd, 2-kv, 3-cs, 4-dfx, 5-dfy, 6-angast, 7-fraca, 8-state
    real, allocatable :: stkparams(:,:) ! 1-smpd, 2-kv, 3-cs, 4-dfx, 5-dfy, 6-angast, 7-fraca, 8-count, 9-box
    real, allocatable :: ptclparams(:,:) ! 1-xpos, 2-ypos, 3-stkind, 4-state, 5-e1, 6-e2, 7-e3, 8-x, 9-y
    logical :: fmovie, fintg, fforctf, fpspec, fsmpd, fkv, fcs, fdf, fangast, ffraca, fstate
    logical :: sstk, ssmpd, skv, scs, sdf, sangast, sfraca, sbox
    logical :: pxpos, pypos, pstkind, pstate, p3d

contains

    procedure :: create
    procedure :: get_files
    procedure :: get_file_parameters
    procedure :: get_stacks
    procedure :: get_stack_parameters
    procedure :: get_particles2D_parameters
    procedure :: get_particles3D_parameters
    procedure :: write_micrographs_star
    procedure :: write_particles2D_star
    procedure :: find_intg

end type relion_project

contains

    subroutine create(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        self%fmovie = .TRUE.
        self%fintg = .TRUE.
        self%fforctf = .TRUE.
        self%fpspec = .TRUE.
        self%fsmpd = .TRUE.
        self%fkv = .TRUE.
        self%fcs = .TRUE.
        self%fdf = .TRUE.
        self%fangast = .TRUE.
        self%ffraca = .TRUE.
        self%fstate = .TRUE.
        self%sstk = .TRUE.
        self%ssmpd = .TRUE.
        self%skv = .TRUE.
        self%scs = .TRUE.
        self%sdf = .TRUE.
        self%sangast = .TRUE.
        self%sfraca = .TRUE.
        self%sbox = .TRUE.
        self%pxpos = .TRUE.
        self%pypos = .TRUE.
        self%pstkind = .TRUE.
        self%pstate = .TRUE.
        self%p3d = .FALSE.

        call self%get_files(spproj)
        call self%get_file_parameters(spproj)
        call self%get_stacks(spproj)
        call self%get_stack_parameters(spproj)
        call self%get_particles2D_parameters(spproj)
        call self%get_particles3D_parameters(spproj)
        call self%write_micrographs_star(spproj)
        call self%write_particles2D_star(spproj)

        if(allocated(self%filesparams))deallocate(self%filesparams)
        if(allocated(self%filestab))deallocate(self%filestab)
        if(allocated(self%stkparams))deallocate(self%stkparams)
        if(allocated(self%stktab))deallocate(self%stktab)
        if(allocated(self%ptclparams))deallocate(self%ptclparams)

    end subroutine create

    subroutine get_files(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        character(len=:), allocatable :: getstring
        integer :: i, lines

        lines = spproj%os_mic%get_noris()

        if(lines > 0 .AND. (.NOT. allocated(self%filestab))) then
            allocate(self%filestab(lines, 4))
        endif

        do i=1, lines
            call spproj%os_mic%getter(i,'movie',getstring)
            if(len(getstring) == 0)then
                self%fmovie = .FALSE.
            else
                self%filestab(i,1) = getstring
            endif

            call spproj%os_mic%getter(i,'intg',getstring)
            if(len(getstring) == 0)then
                self%fintg = .FALSE.
            else
                self%filestab(i,2) = getstring
            endif

            call spproj%os_mic%getter(i,'forctf',getstring)
            if(len(getstring) == 0)then
                self%fforctf = .FALSE.
            else
                self%filestab(i,3) = getstring
            endif

            call spproj%os_mic%getter(i,'pspec',getstring)
            if(len(getstring) == 0)then
                self%fpspec = .FALSE.
            else
                self%filestab(i,4) = getstring
            endif
       end do

       if(allocated(getstring))deallocate(getstring)

    end subroutine get_files

    subroutine get_file_parameters(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        real :: getreal
        integer :: i, lines

        lines = spproj%os_mic%get_noris()

        if(lines > 0 .AND. (.NOT. allocated(self%filesparams))) then
            allocate(self%filesparams(lines, 8))
        endif

        do i=1, lines
            getreal = spproj%os_mic%get(i,'smpd')
            if(getreal == 0)then
                self%fsmpd = .FALSE.
            else
                self%filesparams(i,1) = getreal
            endif

            getreal = spproj%os_mic%get(i,'kv')
            if(getreal == 0)then
                self%fkv = .FALSE.
            else
                self%filesparams(i,2) = getreal
            endif

            getreal = spproj%os_mic%get(i,'cs')
            if(getreal == 0)then
                self%fcs = .FALSE.
            else
                self%filesparams(i,3) = getreal
            endif

            getreal = spproj%os_mic%get(i,'dfx')
            if(getreal == 0)then
                self%fdf = .FALSE.
            else
                self%filesparams(i,4) = getreal
            endif

            getreal = spproj%os_mic%get(i,'dfy')
            if(getreal == 0)then
                self%fdf = .FALSE.
            else
                self%filesparams(i,5) = getreal
            endif

            getreal = spproj%os_mic%get(i,'angast')
            if(getreal == 0)then
                self%fangast = .FALSE.
            else
                self%filesparams(i,6) = getreal
            endif

            getreal = spproj%os_mic%get(i,'fraca')
            if(getreal == 0)then
                self%ffraca = .FALSE.
            else
                self%filesparams(i,7) = getreal
            endif

            if(spproj%os_mic%isthere(i,'state') ) then
                getreal = spproj%os_mic%get(i,'state')
                self%filesparams(i,8) = getreal
            else
                self%filesparams(i,8) = 1
            endif

        end do

    end subroutine get_file_parameters

    subroutine get_stacks(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        character(len=:), allocatable :: getstring
        integer :: i, lines

        lines = spproj%os_stk%get_noris()

        if(lines > 0 .AND. (.NOT. allocated(self%stktab))) then
            allocate(self%stktab(lines, 2))
        endif

        do i=1, lines
            call spproj%os_stk%getter(i,'stk',getstring)
            if(len(getstring) == 0)then
                self%sstk = .FALSE.
            else
                self%stktab(i,1) = getstring
                self%stktab(i,2) = 'false'
            endif
       end do

       if(allocated(getstring))deallocate(getstring)

    end subroutine get_stacks

    subroutine get_stack_parameters(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        real :: getreal
        integer :: i, lines

        lines = spproj%os_stk%get_noris()

        if(lines > 0 .AND. (.NOT. allocated(self%stkparams))) then
            allocate(self%stkparams(lines, 9))
        endif

        do i=1, lines
            getreal = spproj%os_stk%get(i,'smpd')
            if(getreal == 0)then
                self%ssmpd = .FALSE.
            else
                self%stkparams(i,1) = getreal
            endif

            getreal = spproj%os_stk%get(i,'kv')
            if(getreal == 0)then
                self%skv = .FALSE.
            else
                self%stkparams(i,2) = getreal
            endif

            getreal = spproj%os_stk%get(i,'cs')
            if(getreal == 0)then
                self%scs = .FALSE.
            else
                self%stkparams(i,3) = getreal
            endif

            getreal = spproj%os_stk%get(i,'dfx')
            if(getreal == 0)then
                self%sdf = .FALSE.
            else
                self%stkparams(i,4) = getreal
            endif

            getreal = spproj%os_stk%get(i,'dfy')
            if(getreal == 0)then
                self%sdf = .FALSE.
            else
                self%stkparams(i,5) = getreal
            endif

            getreal = spproj%os_stk%get(i,'angast')
            if(getreal == 0)then
                self%sangast = .FALSE.
            else
                self%stkparams(i,6) = getreal
            endif

            getreal = spproj%os_stk%get(i,'fraca')
            if(getreal == 0)then
                self%sfraca = .FALSE.
            else
                self%stkparams(i,7) = getreal
            endif

            self%stkparams(i,8) = 1

            getreal = spproj%os_stk%get(i,'box')
            if(getreal == 0)then
                self%sbox = .FALSE.
            else
                self%stkparams(i,9) = getreal
            endif


        end do

    end subroutine get_stack_parameters

    subroutine get_particles2D_parameters(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        real :: getreal
        integer :: i, lines

        lines = spproj%os_ptcl2D%get_noris()

        if(lines > 0 .AND. (.NOT. allocated(self%ptclparams))) then
            allocate(self%ptclparams(lines, 9))
        endif

        do i=1, lines
            getreal = spproj%os_ptcl2D%get(i,'xpos')
            self%ptclparams(i,1) = getreal

            getreal = spproj%os_ptcl2D%get(i,'ypos')
            self%ptclparams(i,2) = getreal

            getreal = spproj%os_ptcl2D%get(i,'stkind')
            if(getreal == 0)then
                self%pstkind = .FALSE.
            else
                self%ptclparams(i,3) = getreal
            endif

            getreal = spproj%os_ptcl2D%get(i,'state')
            self%ptclparams(i,4) = getreal

        end do

    end subroutine get_particles2D_parameters

    subroutine get_particles3D_parameters(self, spproj)

        class(relion_project), intent(inout) :: self
        class(sp_project), intent(inout) :: spproj

        real :: getreal
        integer :: i, lines

        lines = spproj%os_ptcl3D%get_noris()

        if(lines > 0 .AND. (.NOT. allocated(self%ptclparams))) then
            allocate(self%ptclparams(lines, 9))
        endif

        do i=1, lines
            self%p3d = .TRUE.
            getreal = spproj%os_ptcl3D%get(i,'e1')
            self%ptclparams(i,5) = getreal

            getreal = spproj%os_ptcl3D%get(i,'e2')
            self%ptclparams(i,6) = getreal

            getreal = spproj%os_ptcl3D%get(i,'e3')
            self%ptclparams(i,7) = getreal

            getreal = spproj%os_ptcl3D%get(i,'x')
            self%ptclparams(i,8) = getreal

            getreal = spproj%os_ptcl3D%get(i,'y')
            self%ptclparams(i,9) = getreal

        end do

    end subroutine get_particles3D_parameters

    subroutine write_micrographs_star(self, spproj)
        class(relion_project), intent(inout) :: self
        class(sp_project),     intent(inout) :: spproj
        type(starfile_table_type) :: mic_starfile
        integer :: i, lines

        lines = spproj%os_mic%get_noris()

        if(lines == 0)then
            return
        endif

        if(self%fintg) write(logfhandle,*) 'Exporting micrographs with the following fields'
        if(self%fintg) write(logfhandle,*) ' - intg'
        if(self%fmovie) write(logfhandle,*) ' - movie'
        if(self%fforctf) write(logfhandle,*) ' - forctf'
        if(self%fpspec) write(logfhandle,*) ' - pspec'
        if(self%fsmpd) write(logfhandle,*) ' - smpd'
        if(self%fkv) write(logfhandle,*) ' - kv'
        if(self%fcs) write(logfhandle,*) ' - cs'
        if(self%fdf) write(logfhandle,*) ' - dfx'
        if(self%fdf) write(logfhandle,*) ' - dfy'
        if(self%fangast) write(logfhandle,*) ' - angast'
        if(self%ffraca) write(logfhandle,*) ' - fraca'

        if(self%fintg) call simple_mkdir('micrographs', errmsg= "simple_relion:: create micrographs directory")
        if(self%fmovie) call simple_mkdir('movies', errmsg= "simple_relion:: create movies directory")

        call starfile_table__new(mic_starfile)
        call starfile_table__setcomment(mic_starfile, "SIMPLE 3.0; export_relion")

        do i=1, lines
            if(self%fstate .AND. (self%filesparams(i,8) > 0)) then
                call starfile_table__addObject(mic_starfile)
                if(self%fintg) then
                    call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_NAME, &
                        'micrographs/' // basename(self%filestab(i, 2)))
                    call syslib_symlink('../' // self%filestab(i, 2), 'micrographs/' // basename(self%filestab(i, 2)), 'Failed to generate symlink')
                endif
                if(self%fmovie) then
                    call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_MOVIE_NAME, &
                        'movies/' // basename(self%filestab(i, 1)))
                    call syslib_symlink(self%filestab(i, 1), 'movies/' // basename(self%filestab(i, 1)), 'Failed to generate symlink')
                endif
                if(self%fforctf) then
                    call starfile_table__setValue_string(mic_starfile, EMDL_MICROGRAPH_NAME_WODOSE, &
                        'micrographs/' // basename(self%filestab(i, 3)))
                    call syslib_symlink('../' // self%filestab(i, 3), 'micrographs/' // basename(self%filestab(i, 3)), 'Failed to generate symlink')
                endif

                if(self%fpspec) then
                    call starfile_table__setValue_string(mic_starfile, EMDL_CTF_IMAGE, &
                        'micrographs/' // basename(self%filestab(i, 4)))
                    call syslib_symlink('../' // self%filestab(i, 4), 'micrographs/' // basename(self%filestab(i, 4)), 'Failed to generate symlink')
                endif

                if(self%fsmpd) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DETECTOR_PIXEL_SIZE, &
                        real(self%filesparams(i, 1), dp))
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_MAGNIFICATION, real(10000, dp))
                end if
                if(self%fkv) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_VOLTAGE, &
                        real(int(self%filesparams(i, 2)), dp))
                end if
                if(self%fcs) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_CS, &
                        real(self%filesparams(i, 3), dp))
                end if
                if(self%fdf) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_ASTIGMATISM, &
                        real(abs(self%filesparams(i, 4) * 10000) - (self%filesparams(i, 5) * 10000),dp))
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DEFOCUSU, &
                        real(self%filesparams(i, 4) * 10000, dp))
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DEFOCUSV, &
                        real(self%filesparams(i, 5) * 10000, dp))
                end if
                if(self%fangast) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_DEFOCUS_ANGLE, &
                        real(self%filesparams(i, 6), dp))
                end if
                if(self%ffraca) then
                    call starfile_table__setValue_double(mic_starfile, EMDL_CTF_Q0, &
                        real(self%filesparams(i, 7), dp))
                end if
            end if
        end do

        call starfile_table__open(mic_starfile, 'micrographs.star')
        call starfile_table__write(mic_starfile)
        call starfile_table__close(mic_starfile)

        call starfile_table__delete(mic_starfile)

    end subroutine write_micrographs_star

    subroutine write_particles2D_star(self, spproj)
        class(relion_project), intent(inout) :: self
        class(sp_project),     intent(inout) :: spproj
        type(starfile_table_type) :: ptcl2D_starfile
        integer :: i, lines
        character(len=LONGSTRLEN) :: imgfilename
        character(len=LONGSTRLEN) :: intgname

        lines = spproj%os_ptcl2D%get_noris()

        if(lines == 0)then
            return
        endif

        if(self%sstk) write(logfhandle,*) 'Exporting particles with the following fields'
        if(self%sstk) write(logfhandle,*) ' - stk'
        if(self%ssmpd) write(logfhandle,*) ' - smpd'
        if(self%skv) write(logfhandle,*) ' - kv'
        if(self%scs) write(logfhandle,*) ' - cs'
        if(self%sdf) write(logfhandle,*) ' - dfx'
        if(self%sdf) write(logfhandle,*) ' - dfy'
        if(self%sangast) write(logfhandle,*) ' - angast'
        if(self%sfraca) write(logfhandle,*) ' - fraca'
        if(self%pxpos) write(logfhandle,*) ' - xpos'
        if(self%pypos) write(logfhandle,*) ' - ypos'
        if(self%p3d) then
            write(logfhandle,*) ' - e1'
            write(logfhandle,*) ' - e2'
            write(logfhandle,*) ' - e3'
            write(logfhandle,*) ' - x'
            write(logfhandle,*) ' - y'
        endif

        if(self%sstk) call simple_mkdir('particles', errmsg= "simple_relion:: create particles directory")

        call starfile_table__new(ptcl2D_starfile)
        call starfile_table__setcomment(ptcl2D_starfile, "SIMPLE 3.0; export_relion")

        !STAR data
        do i=1, lines
            if(self%pstate .AND. (self%ptclparams(i,4) > 0)) then
                call starfile_table__addObject(ptcl2D_starfile)
                if(self%sstk) then
                    imgfilename = trim(int2str(int(self%stkparams(int(self%ptclparams(i, 3)), 8))))
                    imgfilename = imgfilename // '@'
                    imgfilename = imgfilename // 'particles/' // basename(self%stktab(int(self%ptclparams(i, 3)), 1)) // 's   '
                    call starfile_table__setValue_string(ptcl2D_starfile, EMDL_IMAGE_NAME, &
                        imgfilename)
                    if(self%stktab(int(self%ptclparams(i, 3)), 2) == 'false') then
                        call syslib_symlink('../' // self%stktab(int(self%ptclparams(i, 3)), 1), 'particles/' // basename(self%stktab(int(self%ptclparams(i, 3)), 1)) // 's', 'Failed to generate symlink')
                        self%stktab(int(self%ptclparams(i, 3)), 2) = 'true'
                    endif
                    if(self%fintg) then
                      intgname = basename(self%stktab(int(self%ptclparams(i, 3)), 1))
                      call self%find_intg(intgname)
                      call starfile_table__setCValue_string(ptcl2D_starfile, EMDL_MICROGRAPH_NAME, &
                          'micrographs/' // trim(intgname))
                    endif
                endif
                if(self%pxpos .AND. self%pypos .AND. self%sbox ) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_IMAGE_COORD_X, &
                        real(int(self%ptclparams(i, 1) + (self%stkparams(int(self%ptclparams(i, 3)), 9) / 2)), dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_IMAGE_COORD_Y, &
                        real(int(self%ptclparams(i, 2) + (self%stkparams(int(self%ptclparams(i, 3)), 9) / 2)), dp))
                endif
                if(self%ssmpd) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_DETECTOR_PIXEL_SIZE, &
                        real(self%stkparams(int(self%ptclparams(i, 3)), 1), dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_MAGNIFICATION, &
                        real(10000, dp))
                end if
                if(self%skv) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_VOLTAGE, &
                        real(int(self%stkparams(int(self%ptclparams(i, 3)), 2)), dp))
                end if
                if(self%scs) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_CS, &
                        real(self%stkparams(int(self%ptclparams(i, 3)), 3), dp))
                end if
                if(self%sdf) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_DEFOCUSU, &
                        real(self%stkparams(int(self%ptclparams(i, 3)), 4) * 10000, dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_DEFOCUSV, &
                        real(self%stkparams(int(self%ptclparams(i, 3)), 5) * 10000, dp))
                end if
                if(self%sangast) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_DEFOCUS_ANGLE, &
                        real(self%stkparams(int(self%ptclparams(i, 3)), 6), dp))
                end if
                if(self%sfraca) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_CTF_Q0, &
                        real(self%stkparams(int(self%ptclparams(i, 3)), 7), dp))
                end if
                if(self%p3d) then
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_ORIENT_ROT, &
                        real(self%ptclparams(i, 5), dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_ORIENT_TILT, &
                        real(self%ptclparams(i, 6),dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_ORIENT_PSI, &
                        real(self%ptclparams(i, 7), dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_ORIENT_ORIGIN_X, &
                        real(self%ptclparams(i, 8), dp))
                    call starfile_table__setValue_double(ptcl2D_starfile, EMDL_ORIENT_ORIGIN_Y, &
                        real(self%ptclparams(i, 9), dp))
                endif
            endif
            self%stkparams(int(self%ptclparams(i, 3)), 8) = self%stkparams(int(self%ptclparams(i, 3)), 8) + 1
        end do

        call starfile_table__open(ptcl2D_starfile, 'particles.star')
        call starfile_table__write(ptcl2D_starfile)
        call starfile_table__close(ptcl2D_starfile)

        call starfile_table__delete(ptcl2D_starfile)

    end subroutine write_particles2D_star

    subroutine find_intg(self, intgname)
        class(relion_project), intent(inout) :: self

        character(*) :: intgname

        intgname = intgname(12:len(intgname))
    end subroutine find_intg

end module simple_relion

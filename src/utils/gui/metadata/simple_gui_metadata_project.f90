!@descr: GUI metadata for the top-level SIMPLE project — populated from a project file
!==============================================================================
! MODULE: simple_gui_metadata_project
!
! PURPOSE:
!   Extends gui_metadata_base with fields identifying the SIMPLE project a
!   GUI session is attached to.  set() reads the projinfo segment and the
!   segment record counts of a *.simple project file to populate the project
!   name, project file path, and micrograph/stack/particle counts, without
!   loading the full particle/stack/micrograph segments.  The Unix timestamp
!   at which the metadata was first assigned is also recorded.
!
! TYPES:
!   gui_metadata_project — extends gui_metadata_base
!     set()     — read projfile and populate project name and record counts,
!                 or populate the same fields from an in-memory sp_project
!     get()     — retrieve project name, project file path, record counts,
!                 and created timestamp; returns the l_assigned flag
!     jsonise() — serialise all fields to a json_value tree (base override)
!
! DEPENDENCIES:
!   unix, json_kinds, json_module, simple_defs, simple_defs_fname, simple_fileio,
!   simple_string, simple_error, simple_string_utils, simple_sp_project, simple_gui_metadata_base,
!   simple_gui_metadata_types, simple_gui_metadata_micrograph, simple_gui_metadata_cavg2D
!==============================================================================
module simple_gui_metadata_project
  use unix,                           only: c_long, c_time
  use json_kinds
  use json_module,                    only: json_core, json_value
  use simple_defs,                    only: LONGSTRLEN
  use simple_defs_fname,              only: MRC_EXT, JPG_EXT
  use simple_fileio,                  only: swap_suffix
  use simple_string,                  only: string
  use simple_error,                   only: simple_exception
  use simple_string_utils,            only: int2str
  use simple_sp_project,              only: sp_project
  use simple_gui_metadata_base,       only: gui_metadata_base
  use simple_gui_metadata_types,      only: GUI_METADATA_MICROGRAPH_TYPE, GUI_METADATA_CAVG2D_TYPE
  use simple_gui_metadata_micrograph, only: gui_metadata_micrograph
  use simple_gui_metadata_cavg2D,     only: gui_metadata_cavg2D, sprite_sheet_pos

  implicit none

  public :: gui_metadata_project
  private
#include "simple_local_flags.inc"

  ! one classification stage's set of 2D class-average metadata entries
  type :: gui_metadata_cavg2D_stage
    private
    logical                                :: is_final = .false.
    type(gui_metadata_cavg2D), allocatable :: cavgs(:)
  end type gui_metadata_cavg2D_stage

  type, extends(gui_metadata_base) :: gui_metadata_project
    private
    character(len=LONGSTRLEN)     :: projname = '' ! SIMPLE project name
    character(len=LONGSTRLEN)     :: projfile = '' ! path to the *.simple project file
    integer                       :: nmics    = 0  ! number of records in the mic segment
    integer                       :: nstks    = 0  ! number of records in the stk segment
    integer                       :: nptcls   = 0  ! number of records in the ptcl2D segment
    integer                       :: ncls2D   = 0  ! number of records in the cls2D segment
    integer                       :: created  = 0  ! Unix timestamp of first assignment
    type(gui_metadata_micrograph), allocatable :: meta_micrographs(:)
    type(gui_metadata_cavg2D_stage), allocatable :: meta_cavg2D(:)
  contains
    procedure :: set_1
    procedure :: set_2
    generic   :: set => set_1, set_2
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_project

contains

  ! Read projfile and populate the project name and segment record counts.
  ! Only the projinfo segment and segment headers are read, so this is cheap
  ! even for projects with large particle/stack/micrograph segments.
  subroutine set_1( self, projfile )
    class(gui_metadata_project), intent(inout) :: self
    type(string),                intent(in)    :: projfile
    type(sp_project)                           :: proj
    type(string)                               :: projname
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    call proj%read_segment('projinfo', projfile)
    call proj%read_data_info(projfile, self%nmics, self%nstks, self%nptcls)
    call proj%projinfo%getter(1, 'projname', projname)
    if( .not. self%l_assigned ) self%created = int(c_time(0_c_long))
    self%l_assigned = .true.
    self%projname   = projname%to_char()
    self%projfile   = projfile%to_char()
    call proj%kill()
  end subroutine set_1

  ! Populate the project name and segment record counts from an already
  ! in-memory project, without touching disk.
  subroutine set_2( self, spproj, stage2D )
    class(gui_metadata_project), intent(inout) :: self
    type(sp_project),            intent(inout) :: spproj
    integer,          optional,  intent(in)    :: stage2D
    type(string)                               :: projname, projfile, cavgsstk, cavgsjpg
    integer                                    :: i, nmeta_micrographs, ncls_stk
    integer                                    :: xtiles, ytiles, xtile, ytile
    integer                                    :: nstage2D, array_idx
    logical                                    :: l_final
    real                                       :: smpd_cavgs
    type(gui_metadata_cavg2D_stage), allocatable :: meta_cavg2D_tmp(:)
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_final = present(stage2D)
    if( l_final ) l_final = stage2D == 0
    if( l_final ) then
        nstage2D = 0
    else
        nstage2D = 1
        if( present(stage2D) ) nstage2D = stage2D
        if( nstage2D < 1 ) THROW_HARD('stage2D must be >= 1, or 0 for final output')
    end if
    call spproj%projinfo%getter(1, 'projname', projname)
    call spproj%projinfo%getter(1, 'projfile', projfile)
    if( .not. self%l_assigned ) self%created = int(c_time(0_c_long))
    self%l_assigned = .true.
    self%projname   = projname%to_char()
    self%projfile   = projfile%to_char()
    self%nmics      = spproj%os_mic%get_noris()
    self%nstks      = spproj%os_stk%get_noris()
    self%nptcls     = spproj%os_ptcl2D%get_noris()
    self%ncls2D     = spproj%os_cls2D%get_noris()
    ! add micrographs (max 50)
    if( allocated(self%meta_micrographs) ) deallocate(self%meta_micrographs)
    if(spproj%os_mic%isthere('thumb')) then
        nmeta_micrographs = min(50, self%nmics)
        allocate(self%meta_micrographs(nmeta_micrographs))
        do i = 1, nmeta_micrographs
            call self%meta_micrographs(i)%new(GUI_METADATA_MICROGRAPH_TYPE)
            call self%meta_micrographs(i)%set(path  =spproj%os_mic%get_str(i, "thumb")  , &
                                              dfx   =spproj%os_mic%get(i,     "dfx")    , &
                                              dfy   =spproj%os_mic%get(i,     "dfy")    , &
                                              ctfres=spproj%os_mic%get(i,      "ctfres"), &
                                              i_max =nmeta_micrographs                  , &
                                              i     =i                                    )
        end do
    end if
    ! add 2D classes
    if( nstage2D == 1 ) then
        if( allocated(self%meta_cavg2D) ) deallocate(self%meta_cavg2D)
    end if
    if( self%ncls2D > 0 ) then
        if( l_final ) then
            ! reuse an existing final slot, or append a new one at the end of the stage array
            array_idx = 0
            if( allocated(self%meta_cavg2D) ) then
                do i = 1, size(self%meta_cavg2D)
                    if( self%meta_cavg2D(i)%is_final ) then
                        array_idx = i
                        exit
                    end if
                end do
            end if
            if( array_idx == 0 ) then
                if( .not. allocated(self%meta_cavg2D) ) then
                    allocate(self%meta_cavg2D(1))
                else
                    allocate(meta_cavg2D_tmp(size(self%meta_cavg2D) + 1))
                    meta_cavg2D_tmp(1:size(self%meta_cavg2D)) = self%meta_cavg2D
                    call move_alloc(meta_cavg2D_tmp, self%meta_cavg2D)
                end if
                array_idx = size(self%meta_cavg2D)
            end if
            self%meta_cavg2D(array_idx)%is_final = .true.
        else
            if( .not. allocated(self%meta_cavg2D) ) then
                allocate(self%meta_cavg2D(nstage2D))
            else if( size(self%meta_cavg2D) < nstage2D ) then
                ! grow the stage array, preserving previously recorded stage containers
                allocate(meta_cavg2D_tmp(nstage2D))
                meta_cavg2D_tmp(1:size(self%meta_cavg2D)) = self%meta_cavg2D
                call move_alloc(meta_cavg2D_tmp, self%meta_cavg2D)
            end if
            array_idx = nstage2D
            self%meta_cavg2D(array_idx)%is_final = .false.
        end if
        if( allocated(self%meta_cavg2D(array_idx)%cavgs) ) deallocate(self%meta_cavg2D(array_idx)%cavgs)
        allocate(self%meta_cavg2D(array_idx)%cavgs(self%ncls2D))
        call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd_cavgs, fail=.false.)
        if( ncls_stk /= self%ncls2D ) THROW_HARD('cavgs stack ncls does not match os_cls2D record count')
        cavgsjpg = swap_suffix(cavgsstk, JPG_EXT, MRC_EXT)
        xtiles   = floor(sqrt(real(self%ncls2D)))
        ytiles   = ceiling(real(self%ncls2D) / real(xtiles))
        do i = 1, self%ncls2D
            xtile = mod(i-1, xtiles)
            ytile = (i-1) / xtiles
            call self%meta_cavg2D(array_idx)%cavgs(i)%new(GUI_METADATA_CAVG2D_TYPE)
            call self%meta_cavg2D(array_idx)%cavgs(i)%set(path    = cavgsjpg,               &
                                         mrcpath = cavgsstk,                           &
                                         i       = i,                                  &
                                         i_max   = self%ncls2D,                        &
                                         res     = spproj%os_cls2D%get(i, 'res'),      &
                                         pop     = spproj%os_cls2D%get_int(i, 'pop'),  &
                                         idx     = i,                                  &
                                         sprite  = sprite_sheet_pos(                   &
                                             x = xtile * (100.0 / max(1, xtiles - 1)), &
                                             y = ytile * (100.0 / max(1, ytiles - 1)), &
                                             h = 100 * ytiles,                         &
                                             w = 100 * xtiles)                         )
        end do
    end if

  end subroutine set_2

  ! Retrieve the project name, project file path, segment record counts, and
  ! created timestamp. Returns .true. if the object has been assigned.
  function get( self, projname, projfile, nmics, nstks, nptcls, created ) result( l_assigned )
    class(gui_metadata_project), intent(in)  :: self
    type(string),                 intent(out) :: projname, projfile
    integer,                      intent(out) :: nmics, nstks, nptcls, created
    logical                                   :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    projname   = trim(self%projname)
    projfile   = trim(self%projfile)
    nmics      = self%nmics
    nstks      = self%nstks
    nptcls     = self%nptcls
    created    = self%created
  end function get

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_project), intent(inout) :: self
    type(json_core)                            :: json
    type(json_value),             pointer      :: json_ptr, json_mics_ptr, json_cls2D_ptr, json_stage_ptr
    type(string)                               :: stage_key
    integer                                    :: i_mic, i_stage2D
    logical                                    :: l_add
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'projname', trim(self%projname))
      call json%add(json_ptr, 'projfile', trim(self%projfile))
      call json%add(json_ptr, 'nmics',    self%nmics         )
      call json%add(json_ptr, 'nstks',    self%nstks         )
      call json%add(json_ptr, 'nptcls',   self%nptcls        )
      call json%add(json_ptr, 'ncls2D',   self%ncls2D        )
      call json%add(json_ptr, 'created',  self%created       )
      ! Add micrographs section if available
      if( allocated(self%meta_micrographs) ) then
        l_add = .false.
        call json%create_array(json_mics_ptr, 'micrographs')
        do i_mic=1, size(self%meta_micrographs)
          if( self%meta_micrographs(i_mic)%assigned() ) then
            l_add = .true.
            call json%add(json_mics_ptr, self%meta_micrographs(i_mic)%jsonise())
          endif
        enddo
        if( l_add ) then
          call json%add(json_ptr, json_mics_ptr)
        else
          call json%destroy(json_mics_ptr)
        endif
      endif
      ! Add cls2D section if available
      if( allocated(self%meta_cavg2D) ) then
        l_add = .false.
        call json%create_object(json_cls2D_ptr, 'cls2D')
        do i_stage2D=1, size(self%meta_cavg2D)
          if( self%meta_cavg2D(i_stage2D)%is_final ) then
            stage_key = 'final'
          else
            stage_key = 'stage' // int2str(i_stage2D)
          end if
          call json%create_array(json_stage_ptr, stage_key%to_char())
          if( allocated(self%meta_cavg2D(i_stage2D)%cavgs) ) then
            do i_mic=1, size(self%meta_cavg2D(i_stage2D)%cavgs)
              if( self%meta_cavg2D(i_stage2D)%cavgs(i_mic)%assigned() ) then
                l_add = .true.
                call json%add(json_stage_ptr, self%meta_cavg2D(i_stage2D)%cavgs(i_mic)%jsonise())
              endif
            enddo
          endif
          call json%add(json_cls2D_ptr, json_stage_ptr)
        enddo
        if( l_add ) then
          call json%add(json_ptr, json_cls2D_ptr)
        else
          call json%destroy(json_cls2D_ptr)
        endif
      endif
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_project

!@descr: GUI metadata for the top-level SIMPLE project — a record of its sections, filled by setters
! simple_gui_project_builder reads the sp_project, writes the previews and fills this record:
! set_project first, then start_update, then the setters of the sections an update touches. The
! 2D and 3D sections keep one slot per stage (stage s >= 1 in slot s) and one final slot (stage 0).
! Holds allocatable components, so it is jsonised in-process and never serialised over IPC.
module simple_gui_metadata_project
  use unix,                           only: c_long, c_time
  use json_module,                    only: json_core, json_value
  use simple_defs,                    only: LONGSTRLEN
  use simple_string,                  only: string
  use simple_error,                   only: simple_exception
  use simple_string_utils,            only: int2str
  use simple_gui_metadata_base,       only: gui_metadata_base
  use simple_gui_metadata_micrograph, only: gui_metadata_micrograph
  use simple_gui_metadata_ptcl,       only: gui_metadata_ptcl
  use simple_gui_metadata_cavg2D,     only: gui_metadata_cavg2D
  use simple_gui_metadata_vol3D,      only: gui_metadata_vol3D

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

  ! one refinement stage's set of per-state 3D volume metadata entries
  type :: gui_metadata_vol3D_stage
    private
    logical                               :: is_final = .false.
    type(gui_metadata_vol3D), allocatable :: states(:)
  end type gui_metadata_vol3D_stage

  type, extends(gui_metadata_base) :: gui_metadata_project
    private
    character(len=LONGSTRLEN)     :: projname = '' ! SIMPLE project name
    character(len=LONGSTRLEN)     :: projfile = '' ! path to the *.simple project file
    integer                       :: nmics    = 0  ! number of records in the mic segment
    integer                       :: nmics_selected = 0  ! number of selected records in the mic segment
    integer                       :: nstks    = 0  ! number of records in the stk segment
    integer                       :: nptcls   = 0  ! number of records in the ptcl2D segment
    integer                       :: nptcls_selected = 0  ! number of selected records in the ptcl2D segment
    integer                       :: ncls2D   = 0  ! number of records in the cls2D segment
    integer                       :: ncls2D_selected = 0  ! number of selected records in the cls2D segment
    integer                       :: nstates3D = 0 ! number of states in the ptcl3D 'state' label
    integer                       :: created  = 0  ! Unix timestamp of first assignment
    real                          :: mskdiam   = 0. ! mask diameter (in A) used for the cls2D run
    real                          :: mskscale  = 0. ! cavgs box size in A (box * smpd), for overlay scaling
    integer                       :: dim_cavgs = 0  ! cavgs box size (in pixels)
    integer                       :: xdim_mic   = 0 ! micrograph width in pixels
    integer                       :: ydim_mic   = 0 ! micrograph height in pixels
    real                          :: smpd_mic   = 0. ! micrograph pixel size (in A)
    integer                       :: pspec_size = 0 ! power spectrum thumbnail size (in pixels)
    character(len=LONGSTRLEN)     :: ptcls_jpg = '' ! path to a JPEG montage of a random particle sample
    integer                       :: nptcls_shown = 0 ! number of particles included in ptcls_jpg
    type(gui_metadata_micrograph),   allocatable :: meta_movies(:)
    type(gui_metadata_micrograph),   allocatable :: meta_micrographs(:)
    type(gui_metadata_cavg2D_stage), allocatable :: meta_cavg2D(:)
    type(gui_metadata_vol3D_stage),  allocatable :: meta_vol3D(:)
    type(gui_metadata_ptcl),         allocatable :: meta_ptcls(:)
  contains
    procedure :: kill => kill_override
    procedure :: set_project
    procedure :: start_update
    procedure :: set_summary
    procedure :: set_movies
    procedure :: set_micrographs
    procedure :: set_particles
    procedure :: set_cavgs2D
    procedure :: set_vols3D
    procedure :: get
    procedure :: serialise => serialise_override
    procedure :: jsonise => jsonise_override
  end type gui_metadata_project

contains

  !---------------- setters ----------------

  ! Name the project and mark the record assigned; the first call stamps the created time.
  subroutine set_project( self, projname, projfile )
    class(gui_metadata_project), intent(inout) :: self
    type(string),                intent(in)    :: projname, projfile
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not. self%l_assigned ) self%created = int(c_time(0_c_long))
    self%l_assigned = .true.
    self%projname   = projname%to_char()
    self%projfile   = projfile%to_char()
  end subroutine set_project

  ! Start an update: the movie, micrograph and particle lists are dropped, and so are the 2D and
  ! 3D stages when @p stage is 1, the first of a new run. The counts and sizes are kept.
  subroutine start_update( self, stage )
    class(gui_metadata_project), intent(inout) :: self
    integer,                     intent(in)    :: stage
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( allocated(self%meta_movies)      ) deallocate(self%meta_movies)
    if( allocated(self%meta_micrographs) ) deallocate(self%meta_micrographs)
    if( allocated(self%meta_ptcls)       ) deallocate(self%meta_ptcls)
    if( stage == 1 )then
      if( allocated(self%meta_cavg2D) ) deallocate(self%meta_cavg2D)
      if( allocated(self%meta_vol3D)  ) deallocate(self%meta_vol3D)
    endif
  end subroutine start_update

  ! Set the counts and sizes given; the others keep their values.
  subroutine set_summary( self, nmics, nmics_selected, nstks, nptcls, nptcls_selected, ncls2D, ncls2D_selected, nstates3D, pspec_size )
    class(gui_metadata_project), intent(inout) :: self
    integer, optional,           intent(in)    :: nmics, nmics_selected, nstks, nptcls, nptcls_selected
    integer, optional,           intent(in)    :: ncls2D, ncls2D_selected, nstates3D, pspec_size
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( present(nmics)           ) self%nmics           = nmics
    if( present(nmics_selected)  ) self%nmics_selected  = nmics_selected
    if( present(nstks)           ) self%nstks           = nstks
    if( present(nptcls)          ) self%nptcls          = nptcls
    if( present(nptcls_selected) ) self%nptcls_selected = nptcls_selected
    if( present(ncls2D)          ) self%ncls2D          = ncls2D
    if( present(ncls2D_selected) ) self%ncls2D_selected = ncls2D_selected
    if( present(nstates3D)       ) self%nstates3D       = nstates3D
    if( present(pspec_size)      ) self%pspec_size      = pspec_size
  end subroutine set_summary

  ! The movie thumbnails.
  subroutine set_movies( self, movies )
    class(gui_metadata_project),   intent(inout) :: self
    type(gui_metadata_micrograph), intent(in)    :: movies(:)
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%meta_movies = movies
  end subroutine set_movies

  ! The micrographs shown, with the dimensions (pixels) and pixel size (A) of the micrographs.
  subroutine set_micrographs( self, micrographs, xdim, ydim, smpd )
    class(gui_metadata_project),   intent(inout) :: self
    type(gui_metadata_micrograph), intent(in)    :: micrographs(:)
    integer,                       intent(in)    :: xdim, ydim
    real,                          intent(in)    :: smpd
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%meta_micrographs = micrographs
    self%xdim_mic         = xdim
    self%ydim_mic         = ydim
    self%smpd_mic         = smpd
  end subroutine set_micrographs

  ! The particle sample: the montage @p ptcls_jpg and one entry per tile.
  subroutine set_particles( self, ptcls_jpg, ptcls )
    class(gui_metadata_project), intent(inout) :: self
    type(string),                intent(in)    :: ptcls_jpg
    type(gui_metadata_ptcl),     intent(in)    :: ptcls(:)
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%ptcls_jpg    = ptcls_jpg%to_char()
    self%nptcls_shown = size(ptcls)
    self%meta_ptcls   = ptcls
  end subroutine set_particles

  ! The class averages of a 2D stage (0 the final one), with the class-average box (pixels), the
  ! mask diameter (A) and the box in A.
  subroutine set_cavgs2D( self, stage, cavgs, dim_cavgs, mskdiam, mskscale )
    class(gui_metadata_project), intent(inout) :: self
    integer,                     intent(in)    :: stage
    type(gui_metadata_cavg2D),   intent(in)    :: cavgs(:)
    integer,                     intent(in)    :: dim_cavgs
    real,                        intent(in)    :: mskdiam, mskscale
    type(gui_metadata_cavg2D_stage), allocatable :: tmp(:)
    integer :: islot, nslots
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( stage < 0 ) THROW_HARD('stage must be >= 1, or 0 for the final one')
    if( allocated(self%meta_cavg2D) )then
      call stage_slot(self%meta_cavg2D%is_final, stage, islot, nslots)
      if( nslots > size(self%meta_cavg2D) )then
        allocate(tmp(nslots))
        tmp(1:size(self%meta_cavg2D)) = self%meta_cavg2D
        call move_alloc(tmp, self%meta_cavg2D)
      endif
    else
      call stage_slot([logical ::], stage, islot, nslots)
      allocate(self%meta_cavg2D(nslots))
    endif
    self%meta_cavg2D(islot)%is_final = stage == 0
    self%meta_cavg2D(islot)%cavgs    = cavgs
    self%dim_cavgs = dim_cavgs
    self%mskdiam   = mskdiam
    self%mskscale  = mskscale
  end subroutine set_cavgs2D

  ! The state volumes of a 3D stage (0 the final one).
  subroutine set_vols3D( self, stage, vols )
    class(gui_metadata_project), intent(inout) :: self
    integer,                     intent(in)    :: stage
    type(gui_metadata_vol3D),    intent(in)    :: vols(:)
    type(gui_metadata_vol3D_stage), allocatable :: tmp(:)
    integer :: islot, nslots
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( stage < 0 ) THROW_HARD('stage must be >= 1, or 0 for the final one')
    if( allocated(self%meta_vol3D) )then
      call stage_slot(self%meta_vol3D%is_final, stage, islot, nslots)
      if( nslots > size(self%meta_vol3D) )then
        allocate(tmp(nslots))
        tmp(1:size(self%meta_vol3D)) = self%meta_vol3D
        call move_alloc(tmp, self%meta_vol3D)
      endif
    else
      call stage_slot([logical ::], stage, islot, nslots)
      allocate(self%meta_vol3D(nslots))
    endif
    self%meta_vol3D(islot)%is_final = stage == 0
    self%meta_vol3D(islot)%states   = vols
  end subroutine set_vols3D

  ! The slot of @p stage in a stage array whose slots' final flags are @p is_final: stage s >= 1
  ! is slot s; stage 0 is the final slot, reused when there is one and appended otherwise.
  ! @p nslots is the size the array needs.
  subroutine stage_slot( is_final, stage, islot, nslots )
    logical, intent(in)  :: is_final(:)
    integer, intent(in)  :: stage
    integer, intent(out) :: islot, nslots
    nslots = size(is_final)
    if( stage == 0 )then
      islot = findloc(is_final, .true., 1)
      if( islot == 0 )then
        nslots = nslots + 1
        islot  = nslots
      endif
    else
      islot  = stage
      nslots = max(nslots, stage)
    endif
  end subroutine stage_slot

  !---------------- getters ----------------

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

  !---------------- serialisation ----------------

  ! The record's allocatable components would travel as descriptors of this process's memory:
  ! it is jsonised in-process only, and never sent over a pipe.
  subroutine serialise_override( self, buffer )
    class(gui_metadata_project),           intent(in)    :: self
    character(len=:),         allocatable, intent(inout) :: buffer
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( allocated(buffer) ) deallocate(buffer)
    THROW_HARD('gui_metadata_project holds allocatable components and is jsonised in-process only')
  end subroutine serialise_override

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_project), intent(in)    :: self
    type(json_core)                            :: json
    type(json_value),             pointer      :: json_ptr, json_mics_ptr, json_cls2D_ptr, json_stage_ptr, json_ptcls_ptr
    type(json_value),             pointer      :: json_cls3D_ptr
    type(string)                               :: stage_key
    integer                                    :: i_mic, i_stage2D, i_stage3D, i_state3D, i_ptcl
    logical                                    :: l_add
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'projname', trim(self%projname))
      call json%add(json_ptr, 'projfile', trim(self%projfile))
      if(self%nmics  > 0) call json%add(json_ptr, 'nmics',    self%nmics  )
      if(self%nstks  > 0) call json%add(json_ptr, 'nstks',    self%nstks  )
      if(self%nptcls > 0) call json%add(json_ptr, 'nptcls',   self%nptcls )
      if(self%ncls2D > 0) call json%add(json_ptr, 'ncls2D',   self%ncls2D )
      if(self%nmics_selected >  0) call json%add(json_ptr, 'nmics_selected', self%nmics_selected)
      if(self%nptcls_selected > 0) call json%add(json_ptr, 'nptcls_selected', self%nptcls_selected)
      if(self%ncls2D_selected > 0) call json%add(json_ptr, 'ncls2D_selected', self%ncls2D_selected)
      call json%add(json_ptr, 'created',  self%created       )
      if( self%dim_cavgs > 0 ) call json%add(json_ptr, 'dim_cavgs', self%dim_cavgs)
      if( self%mskscale > 0. ) then
        call json%add(json_ptr, 'mskdiam',  dble(self%mskdiam) )
        call json%add(json_ptr, 'mskscale', dble(self%mskscale))
      end if
      if( self%pspec_size > 0 ) call json%add(json_ptr, 'pspec_size', self%pspec_size)
      if( self%smpd_mic > 0.  ) call json%add(json_ptr, 'smpd_mic', dble(self%smpd_mic))
      if( self%xdim_mic > 0 .and. self%ydim_mic > 0 ) then
        call json%add(json_ptr, 'xdim_mic', self%xdim_mic)
        call json%add(json_ptr, 'ydim_mic', self%ydim_mic)
      end if
      if( len_trim(self%ptcls_jpg) > 0 ) then
        call json%add(json_ptr, 'ptcls_jpg', trim(self%ptcls_jpg))
        call json%add(json_ptr, 'nptcls_shown', self%nptcls_shown)
      end if
      ! Add movies section if available
      if( allocated(self%meta_movies) ) then
        l_add = .false.
        call json%create_array(json_mics_ptr, 'movies')
        do i_mic=1, size(self%meta_movies)
          if( self%meta_movies(i_mic)%assigned() ) then
            l_add = .true.
            call json%add(json_mics_ptr, self%meta_movies(i_mic)%jsonise())
          endif
        enddo
        if( l_add ) then
          call json%add(json_ptr, json_mics_ptr)
        else
          call json%destroy(json_mics_ptr)
        endif
      endif
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
      ! Add particles section if available
      if( allocated(self%meta_ptcls) ) then
        l_add = .false.
        call json%create_array(json_ptcls_ptr, 'particles')
        do i_ptcl=1, size(self%meta_ptcls)
          if( self%meta_ptcls(i_ptcl)%assigned() ) then
            l_add = .true.
            call json%add(json_ptcls_ptr, self%meta_ptcls(i_ptcl)%jsonise())
          endif
        enddo
        if( l_add ) then
          call json%add(json_ptr, json_ptcls_ptr)
        else
          call json%destroy(json_ptcls_ptr)
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
      ! Add cls3D section if available
      if( allocated(self%meta_vol3D) ) then
        l_add = .false.
        call json%create_object(json_cls3D_ptr, 'cls3D')
        do i_stage3D=1, size(self%meta_vol3D)
          if( self%meta_vol3D(i_stage3D)%is_final ) then
            stage_key = 'final'
          else
            stage_key = 'stage' // int2str(i_stage3D)
          end if
          call json%create_array(json_stage_ptr, stage_key%to_char())
          if( allocated(self%meta_vol3D(i_stage3D)%states) ) then
            do i_state3D=1, size(self%meta_vol3D(i_stage3D)%states)
              if( self%meta_vol3D(i_stage3D)%states(i_state3D)%assigned() ) then
                l_add = .true.
                call json%add(json_stage_ptr, self%meta_vol3D(i_stage3D)%states(i_state3D)%jsonise())
              endif
            enddo
          endif
          call json%add(json_cls3D_ptr, json_stage_ptr)
        enddo
        if( l_add ) then
          call json%add(json_ptr, json_cls3D_ptr)
        else
          call json%destroy(json_cls3D_ptr)
        endif
      endif
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_project), intent(inout) :: self
    select type( self )
      type is( gui_metadata_project )
        self = gui_metadata_project()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_project

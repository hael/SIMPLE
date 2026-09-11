!@descr: GUI metadata for the top-level SIMPLE project — populated from a project file
!==============================================================================
! MODULE: simple_gui_metadata_project
!
! PURPOSE:
!   Extends gui_metadata_base with fields identifying the SIMPLE project a
!   GUI session is attached to.  set() populates the project name, project
!   file path, and micrograph/stack/particle counts from an in-memory
!   sp_project.  The Unix timestamp at which the metadata was first assigned
!   is also recorded.
!
! TYPES:
!   gui_metadata_project — extends gui_metadata_base
!     set()     — populate project name and record counts from an in-memory sp_project
!     get()     — retrieve project name, project file path, record counts,
!                 and created timestamp; returns the l_assigned flag
!     jsonise() — serialise all fields to a json_value tree (base override)
!
! DEPENDENCIES:
!   unix, json_kinds, json_module, simple_defs, simple_defs_fname, simple_fileio,
!   simple_string, simple_error, simple_string_utils, simple_sp_project, simple_gui_metadata_base,
!   simple_gui_metadata_types, simple_gui_metadata_micrograph, simple_gui_metadata_cavg2D,
!   simple_image, simple_math, simple_motion_gain_helpers,
!   simple_procimgstk, simple_gui_utils, simple_syslib
!==============================================================================
module simple_gui_metadata_project
  use unix,                           only: c_long, c_time
  use json_kinds
  use json_module,                    only: json_core, json_value
  use simple_defs,                    only: LONGSTRLEN, GUI_PSPECSZ, SHORTSTRLEN
  use simple_defs_fname,              only: MRC_EXT, JPG_EXT
  use simple_fileio,                  only: swap_suffix, file_exists
  use simple_string,                  only: string
  use simple_error,                   only: simple_exception
  use simple_string_utils,            only: int2str
  use simple_sp_project,              only: sp_project
  use simple_gui_metadata_base,       only: gui_metadata_base
  use simple_gui_metadata_types,      only: GUI_METADATA_MICROGRAPH_TYPE, GUI_METADATA_CAVG2D_TYPE, GUI_METADATA_PTCL_TYPE
  use simple_gui_metadata_micrograph, only: gui_metadata_micrograph
  use simple_gui_metadata_ptcl,       only: gui_metadata_ptcl
  use simple_gui_metadata_cavg2D,     only: gui_metadata_cavg2D, sprite_sheet_pos
  use simple_nrtxtfile,               only: nrtxtfile
  use simple_image,                   only: image
  use simple_math,                    only: round2even
  use simple_motion_gain_helpers,     only: read_movies_and_sum_frames
  use simple_procimgstk,              only: random_selection_from_imgfile, bp_imgfile
  use simple_gui_utils,               only: mrc2jpeg_tiled
  use simple_syslib,                  only: del_file, simple_abspath

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
    integer                       :: nmics_selected = 0  ! number of selected records in the mic segment
    integer                       :: nstks    = 0  ! number of records in the stk segment
    integer                       :: nptcls   = 0  ! number of records in the ptcl2D segment
    integer                       :: nptcls_selected = 0  ! number of selected records in the ptcl2D segment
    integer                       :: ncls2D   = 0  ! number of records in the cls2D segment
    integer                       :: ncls2D_selected = 0  ! number of selected records in the cls2D segment
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
    type(gui_metadata_ptcl),         allocatable :: meta_ptcls(:)
  contains
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_project

contains

  ! Populate the project name and segment record counts from an already
  ! in-memory project, without touching disk.
  subroutine set( self, spproj, oritype, stage2D, selection )
    class(gui_metadata_project),     intent(inout) :: self
    type(sp_project),                intent(inout) :: spproj
    character(len=*),    optional,   intent(in)    :: oritype
    integer,             optional,   intent(in)    :: stage2D
    logical,             optional,   intent(in)    :: selection
    type(gui_metadata_cavg2D_stage), allocatable   :: meta_cavg2D_tmp(:)
    type(gui_metadata_micrograph),   allocatable   :: meta_micrographs_tmp(:), meta_movies_tmp(:)
    type(gui_metadata_cavg2D),       allocatable   :: cavgs_tmp(:)
    real,                            allocatable   :: boxdata(:)
    integer,                         allocatable   :: micrograph_indices(:)
    type(nrtxtfile)                                :: boxfile
    type(string)                                   :: projname, projfile, cavgsstk, cavgsjpg, boxpath
    character(len=SHORTSTRLEN)                     :: md_oritype
    integer                                        :: i, j, x, y, nmeta_micrographs, ncls_stk, n_valid_micrographs, n_valid_cavgs
    integer                                        :: xtiles, ytiles, xtile, ytile, nrecs, nlines
    integer                                        :: nstage2D, array_idx, out_ind
    integer                                        :: nptcls_all, nptcls_valid, nptcls_sample, box_ptcls, n_valid_ptcls
    logical                                        :: l_final, l_selection, l_ctf
    real                                           :: smpd_cavgs, box_cavgs, mskdiam_cavgs, smpd_ptcls
    type(string)                                   :: ptclsstk, ptclsjpg, ptclslpstk, ptclsjpglp
    integer,                            parameter  :: N_PTCLS_SAMPLE = 100
    integer,                            parameter  :: N_MOV_THUMBS = 10
    type(image)                                    :: movsum, movthumb
    type(string)                                   :: movfname, movthumbfname
    integer                                        :: n_movthumbs, n_movies_sum, n_frames_sum, ldim_mov(3), ldim_thumb(3)
    real                                            :: scale_thumb
                
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_final = present(stage2D)
    if( l_final ) l_final = stage2D == 0
    l_selection = .false.
    if( present(selection) ) l_selection = selection
    md_oritype = 'all'
    if( present(oritype) ) md_oritype = oritype
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
    
   ! self%nstks      = spproj%os_stk%get_noris() ! do we really need this?
   ! self%nptcls     = spproj%os_ptcl2D%get_noris()
   ! self%ncls2D     = spproj%os_cls2D%get_noris()
   ! self%pspec_size = GUI_PSPECSZ
    ! add movies (max 10)
    if( allocated(self%meta_movies) ) deallocate(self%meta_movies)
    if( md_oritype == 'mov' .or. md_oritype == 'all' ) then
        self%nmics  = spproj%os_mic%get_noris()
        n_movthumbs = 0
        if( self%nmics > 0 ) allocate(self%meta_movies(N_MOV_THUMBS)) 
        do i = 1, self%nmics
            if( n_movthumbs >= N_MOV_THUMBS ) exit
            if( .not. spproj%os_mic%isthere(i, 'imgkind') ) cycle
            if( spproj%os_mic%get_str(i, 'imgkind') /= 'movie' ) cycle
            if( spproj%os_mic%isthere(i, 'movthumb') ) cycle ! already generated
            n_movthumbs = n_movthumbs + 1
            movfname = spproj%os_mic%get_str(i, 'movie')
            call read_movies_and_sum_frames([movfname], spproj%os_mic%get(i, 'smpd'), movsum, n_movies_sum, n_frames_sum)
            ldim_mov      = movsum%get_ldim()
            scale_thumb   = real(GUI_PSPECSZ) / real(ldim_mov(1))
            ldim_thumb(1) = round2even(real(ldim_mov(1)) * scale_thumb)
            ldim_thumb(2) = round2even(real(ldim_mov(2)) * scale_thumb)
            ldim_thumb(3) = 1
            call movthumb%new(ldim_thumb, spproj%os_mic%get(i, 'smpd'))
            call movsum%fft()
            call movsum%clip(movthumb)
            call movthumb%ifft()
            movthumbfname = 'movthumb' // int2str(i) // JPG_EXT
            call movthumb%write_jpg(movthumbfname, norm=.true., quality=90)
            movthumbfname = simple_abspath(movthumbfname)
            call self%meta_movies(n_movthumbs)%new(GUI_METADATA_MICROGRAPH_TYPE)
            call self%meta_movies(n_movthumbs)%set(path  = movthumbfname  , &
                                                  i_max  = N_MOV_THUMBS   , &
                                                  i      = i                )                                 
            call movsum%kill()
            call movthumb%kill()
        end do
        ! trim unused (unassigned) slots left by skipped micrographs
        if( n_movthumbs < N_MOV_THUMBS ) then
            allocate(meta_movies_tmp(n_movthumbs))
            meta_movies_tmp = self%meta_movies(1:n_movthumbs)
            call move_alloc(meta_movies_tmp, self%meta_movies)
        end if
    end if
    ! add micrographs (max 50)
    if( allocated(self%meta_micrographs) ) deallocate(self%meta_micrographs)
    if( md_oritype == 'mic' .or. md_oritype == 'ptcl' .or. md_oritype == 'all' ) then
        self%nmics      = spproj%os_mic%get_noris()
        self%pspec_size = GUI_PSPECSZ
        if( md_oritype == 'ptcl' ) then
            self%nstks  = spproj%os_stk%get_noris()
            self%nptcls = spproj%os_ptcl2D%get_noris()
        end if
        if( l_selection) self%nmics_selected = spproj%os_mic%count_state_gt_zero()
        if(spproj%os_mic%isthere('thumb')) then
            nmeta_micrographs = min(50, self%nmics)
            allocate(self%meta_micrographs(nmeta_micrographs))
            self%xdim_mic       = nint(spproj%os_mic%get(1, "xdim"))
            self%ydim_mic       = nint(spproj%os_mic%get(1, "ydim"))
            self%smpd_mic       = spproj%os_mic%get(1, "smpd")
            n_valid_micrographs = 0
            l_ctf = .false.
            if( spproj%os_mic%isthere('ctfjpg') ) l_ctf = .true.
            do i = 1, nmeta_micrographs
                if( spproj%os_mic%get_state(i) == 0 ) cycle ! needs improvement to work with pagination
                n_valid_micrographs = n_valid_micrographs + 1
                call self%meta_micrographs(n_valid_micrographs)%new(GUI_METADATA_MICROGRAPH_TYPE)
                if( l_ctf ) then
                    call self%meta_micrographs(n_valid_micrographs)%set(path  =spproj%os_mic%get_str(i, "thumb")  , &
                                                      dfx   =spproj%os_mic%get(i,     "dfx")    , &
                                                      dfy   =spproj%os_mic%get(i,     "dfy")    , &
                                                      ctfres=spproj%os_mic%get(i,      "ctfres"), &
                                                      ctfimg=spproj%os_mic%get_str(i, "ctfjpg") , &
                                                      i_max =nmeta_micrographs                  , &
                                                      i     =i                                    )
                else
                    call self%meta_micrographs(n_valid_micrographs)%set(path  =spproj%os_mic%get_str(i, "thumb")  , &
                                                      dfx   =spproj%os_mic%get(i,     "dfx")    , &
                                                      dfy   =spproj%os_mic%get(i,     "dfy")    , &
                                                      ctfres=spproj%os_mic%get(i,      "ctfres"), &
                                                      i_max =nmeta_micrographs                  , &
                                                      i     =i                                    )
                end if
                call self%meta_micrographs(n_valid_micrographs)%clear_coordinates()
                boxpath = spproj%os_mic%get_str(i, "boxfile")
                if( boxpath%strlen() > 0 .and. file_exists(boxpath) ) then
                    call boxfile%new(boxpath, 1)
                    nrecs  = boxfile%get_nrecs_per_line()
                    nlines = boxfile%get_ndatalines()
                    if( nrecs >= 4 ) then
                        allocate(boxdata(nrecs))
                        do j = 1, nlines
                            call boxfile%readNextDataLine(boxdata)
                            x = nint(boxdata(1) + boxdata(3)/2)
                            y = nint(boxdata(2) + boxdata(4)/2)
                            call self%meta_micrographs(n_valid_micrographs)%set_coordinate(j, x, y, self%xdim_mic, self%ydim_mic)
                        enddo
                        deallocate(boxdata)
                    endif
                    call boxfile%kill()
                endif
            end do
            ! trim unused (unassigned) slots left by skipped micrographs
            if( n_valid_micrographs < nmeta_micrographs ) then
                allocate(meta_micrographs_tmp(n_valid_micrographs))
                meta_micrographs_tmp = self%meta_micrographs(1:n_valid_micrographs)
                call move_alloc(meta_micrographs_tmp, self%meta_micrographs)
            end if
        end if
    end if
    ! add particles: JPEG montage of a random sample of (selected) particles
    if( allocated(self%meta_ptcls) ) deallocate(self%meta_ptcls)
    if( md_oritype == 'ptcl' .or. md_oritype == 'all' ) then
        nptcls_all = spproj%os_ptcl2D%get_noris()
        if( nptcls_all > 0 ) then
            if( spproj%os_ptcl2D%isthere('state') ) then
                nptcls_valid = spproj%os_ptcl2D%count_state_gt_zero()
            else
                nptcls_valid = nptcls_all
            end if
            nptcls_sample = min(N_PTCLS_SAMPLE, nptcls_valid)
            if( nptcls_sample > 0 ) then
                box_ptcls  = nint(spproj%os_stk%get(1, 'box'))
                smpd_ptcls = spproj%os_stk%get(1, 'smpd')
                ptclsstk   = 'ptcls_sample' // MRC_EXT
                ptclsjpg   = 'ptcls_sample' // JPG_EXT
                ptclslpstk = 'ptcls_sample_lp' // MRC_EXT
                ptclsjpglp = 'ptcls_sample_lp' // JPG_EXT
                call random_selection_from_imgfile(spproj, ptclsstk, box_ptcls, nptcls_sample, pinds=micrograph_indices)
                call mrc2jpeg_tiled(ptclsstk, ptclsjpg, ntiles=n_valid_ptcls)
                call bp_imgfile(ptclsstk, ptclslpstk, smpd_ptcls, 0., 10.)
                call mrc2jpeg_tiled(ptclslpstk, ptclsjpglp, ntiles=n_valid_ptcls)
                call del_file(ptclslpstk)
                call del_file(ptclsstk)
                ptclsjpg          = simple_abspath(ptclsjpg)
                ptclsjpglp        = simple_abspath(ptclsjpglp)
                self%ptcls_jpg    = ptclsjpg%to_char()
                self%nptcls_shown = n_valid_ptcls
                allocate(self%meta_ptcls(n_valid_ptcls))
                xtiles         = floor(sqrt(real(n_valid_ptcls)))
                ytiles         = ceiling(real(n_valid_ptcls) / real(xtiles))
                n_valid_cavgs  = 0
                do i = 1, n_valid_ptcls
                    n_valid_cavgs = n_valid_cavgs + 1
                    xtile = mod(i-1, xtiles)
                    ytile = (i-1) / xtiles
                    call self%meta_ptcls(i)%new(GUI_METADATA_PTCL_TYPE)
                    call self%meta_ptcls(i)%set(path    = ptclsjpg,                           &
                                                pathlp  = ptclsjpglp,                         &
                                                i       = i,                                  &
                                                i_max   = n_valid_ptcls,                      &
                                                df      = (spproj%os_ptcl2D%get(micrograph_indices(i), 'dfx') + spproj%os_ptcl2D%get(micrograph_indices(i), 'dfy')) / 2.0, &
                                                box     = box_ptcls,                          &
                                                idx     = micrograph_indices(i),              &
                                                sprite  = sprite_sheet_pos(                   &
                                                    x = xtile * (100.0 / max(1, xtiles - 1)), &
                                                    y = ytile * (100.0 / max(1, ytiles - 1)), &
                                                    h = 100 * ytiles,                         &
                                                    w = 100 * xtiles)                         )
                end do
                if( allocated(micrograph_indices) ) deallocate(micrograph_indices)
            end if
        end if
    end if
    ! add 2D classes
    if( nstage2D == 1 ) then
        if( allocated(self%meta_cavg2D) ) deallocate(self%meta_cavg2D)
    end if
    if( md_oritype == 'cls2D' .or. md_oritype == 'all' ) then
        self%nstks  = spproj%os_stk%get_noris()
        self%nptcls = spproj%os_ptcl2D%get_noris()
        self%ncls2D = spproj%os_cls2D%get_noris()
        if( l_selection ) then
            self%nptcls_selected = nint(spproj%os_cls2D%get_sum('pop'))
            self%ncls2D_selected = spproj%os_cls2D%count_state_gt_zero()
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
            box_cavgs = 0.
            out_ind   = 0
            call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd_cavgs, fail=.false., out_ind=out_ind, box=box_cavgs)
            if( ncls_stk /= self%ncls2D ) THROW_HARD('cavgs stack ncls does not match os_cls2D record count')
            mskdiam_cavgs = 0.
            if( out_ind > 0 .and. spproj%os_out%isthere(out_ind, 'mskdiam') ) mskdiam_cavgs = spproj%os_out%get(out_ind, 'mskdiam')
            self%dim_cavgs = nint(box_cavgs)
            self%mskdiam   = mskdiam_cavgs
            self%mskscale  = box_cavgs * smpd_cavgs
            cavgsjpg       = swap_suffix(cavgsstk, JPG_EXT, MRC_EXT)
            xtiles         = floor(sqrt(real(self%ncls2D)))
            ytiles         = ceiling(real(self%ncls2D) / real(xtiles))
            n_valid_cavgs  = 0
            do i = 1, self%ncls2D
                if( spproj%os_cls2D%get_state(i) == 0 ) cycle
                n_valid_cavgs = n_valid_cavgs + 1
                xtile = mod(i-1, xtiles)
                ytile = (i-1) / xtiles
                call self%meta_cavg2D(array_idx)%cavgs(n_valid_cavgs)%new(GUI_METADATA_CAVG2D_TYPE)
                call self%meta_cavg2D(array_idx)%cavgs(n_valid_cavgs)%set(path    = cavgsjpg,               &
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
            ! trim unused (unassigned) slots left by skipped classes
            if( n_valid_cavgs < self%ncls2D ) then
                allocate(cavgs_tmp(n_valid_cavgs))
                cavgs_tmp = self%meta_cavg2D(array_idx)%cavgs(1:n_valid_cavgs)
                call move_alloc(cavgs_tmp, self%meta_cavg2D(array_idx)%cavgs)
            end if
        end if
    end if

  end subroutine set

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
    type(json_value),             pointer      :: json_ptr, json_mics_ptr, json_cls2D_ptr, json_stage_ptr, json_ptcls_ptr
    type(string)                               :: stage_key
    integer                                    :: i_mic, i_stage2D, i_ptcl
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
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_project

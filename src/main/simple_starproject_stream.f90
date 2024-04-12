module simple_starproject_stream
include 'simple_lib.f08'
!$ use omp_lib
use simple_sp_project, only: sp_project
use simple_cmdline,    only: cmdline
use simple_parameters, only: params_glob
use simple_starfile_wrappers
use simple_math
use CPlot2D_wrapper_module
use simple_rnd
use FoX_dom
implicit none

public :: starproject_stream
private
#include "simple_local_flags.inc"

type starproject_stream
    type(starfile_table_type)  :: starfile
    character(len=LONGSTRLEN)  :: starfile_name
    character(len=LONGSTRLEN)  :: starfile_tmp
    character(len=LONGSTRLEN)  :: rootpath
    real                       :: shift_threshold = 0.05
    integer                    :: group_offset = 0
    logical                    :: nicestream   = .false.
    logical                    :: use_beamtilt = .false.
contains
    ! export
    procedure          :: stream_export_micrographs
    procedure          :: stream_export_optics
    !starfile
    procedure, private :: starfile_init
    procedure, private :: starfile_deinit
    procedure, private :: starfile_write_table
    procedure, private :: starfile_set_optics_table
    procedure, private :: starfile_set_optics_group_table
    procedure, private :: starfile_set_micrographs_table
    !optics
    procedure, private :: assign_optics_single
    procedure, private :: assign_optics


end type starproject_stream

contains

    ! starfiles

    subroutine starfile_init( self, fname, outdir)
        class(starproject_stream),  intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        character(len=LONGSTRLEN),  intent(in)    :: outdir
        character(len=XLONGSTRLEN)                :: cwd
        character(len=:),           allocatable   :: stem
        self%starfile_name = fname
        self%starfile_tmp  = fname // '.tmp'
        if(len_trim(outdir) > 0) then
            call simple_getcwd(cwd)
            stem = basename(stemname(cwd))
            self%rootpath = trim(stem)
            self%nicestream = .true.
        end if
        call starfile_table__new(self%starfile)
        if(allocated(stem)) deallocate(stem)
    end subroutine starfile_init

    subroutine starfile_deinit( self )
        class(starproject_stream),  intent(inout) :: self
        call starfile_table__delete(self%starfile)
        if(file_exists(trim(adjustl(self%starfile_tmp)))) then
            if(file_exists(trim(adjustl(self%starfile_name)))) call del_file(trim(adjustl(self%starfile_name)))
            call simple_rename(trim(adjustl(self%starfile_tmp)), trim(adjustl(self%starfile_name)))
        end if
    end subroutine starfile_deinit

    subroutine starfile_write_table( self, append )
        class(starproject_stream),  intent(inout) :: self
        logical,                    intent(in)    :: append
        integer                                   :: append_int
        append_int = 0
        if(append) append_int = 1
        call starfile_table__open_ofile(self%starfile, trim(self%starfile_tmp), append_int)
        call starfile_table__write_ofile(self%starfile)
        call starfile_table__close_ofile(self%starfile)
    end subroutine starfile_write_table

    subroutine starfile_set_optics_table( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer                                   :: i
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'optics')
        do i=1,spproj%os_optics%get_noris()
            call starfile_table__addObject(self%starfile)
            if(spproj%os_optics%get(i, 'state') .eq. 0.0 ) cycle
            ! ints
            if(spproj%os_optics%isthere(i, 'ogid'))   call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_optics%get(i, 'ogid')))
            if(spproj%os_optics%isthere(i, 'pop'))    call starfile_table__setValue_int(self%starfile, SMPL_OPTICS_POPULATION,  int(spproj%os_optics%get(i, 'pop' )))
            ! doubles
            if(spproj%os_optics%isthere(i, 'kv'))     call starfile_table__setValue_double(self%starfile, EMDL_CTF_VOLTAGE,      real(spproj%os_optics%get(i, 'kv'),    dp))
            if(spproj%os_optics%isthere(i, 'smpd'))   call starfile_table__setValue_double(self%starfile, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_optics%get(i, 'smpd'),  dp))
            if(spproj%os_optics%isthere(i, 'cs'))     call starfile_table__setValue_double(self%starfile, EMDL_CTF_CS,           real(spproj%os_optics%get(i, 'cs'),    dp))
            if(spproj%os_optics%isthere(i, 'fraca'))  call starfile_table__setValue_double(self%starfile, EMDL_CTF_Q0,           real(spproj%os_optics%get(i, 'fraca'), dp))
            if(spproj%os_optics%isthere(i, 'opcx'))   call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_CENTROIDX, real(spproj%os_optics%get(i, 'opcx'),  dp))
            if(spproj%os_optics%isthere(i, 'opcy'))   call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_CENTROIDY, real(spproj%os_optics%get(i, 'opcy'),  dp))
            ! strings
            if(spproj%os_optics%isthere(i, 'ogname')) call starfile_table__setValue_string(self%starfile, EMDL_IMAGE_OPTICS_GROUP_NAME, trim(adjustl(spproj%os_optics%get_static(i, 'ogname'))))
        end do
    end subroutine starfile_set_optics_table

    subroutine starfile_set_optics_group_table( self, spproj, ogid)
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer,                    intent(in)    :: ogid
        integer                                   :: i
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'opticsgroup_' // int2str(ogid))
        do i=1,spproj%os_mic%get_noris()
            if(spproj%os_mic%isthere(i, 'ogid') .and. int(spproj%os_mic%get(i, 'ogid')) == ogid) then
                call starfile_table__addObject(self%starfile)
                if(spproj%os_mic%get(i, 'state') .eq. 0.0 ) cycle
                ! ints
                call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'ogid')))
                ! doubles
                if(spproj%os_mic%isthere(i, 'shiftx')) call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_SHIFTX, real(spproj%os_mic%get(i, 'shiftx'), dp))
                if(spproj%os_mic%isthere(i, 'shifty')) call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_SHIFTY, real(spproj%os_mic%get(i, 'shifty'), dp))
            end if
        end do
    end subroutine starfile_set_optics_group_table

    subroutine starfile_set_micrographs_table( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer                                   :: i
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'micrographs')
        do i=1,spproj%os_mic%get_noris()
            call starfile_table__addObject(self%starfile)
            if(spproj%os_mic%get(i, 'state') .eq. 0.0 ) cycle
            ! ints
            if(spproj%os_mic%isthere(i, 'ogid'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, int(spproj%os_mic%get(i, 'ogid'   )))
            if(spproj%os_mic%isthere(i, 'xdim'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_SIZE_X,       int(spproj%os_mic%get(i, 'xdim'   )))           
            if(spproj%os_mic%isthere(i, 'ydim'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_SIZE_Y,       int(spproj%os_mic%get(i, 'ydim'   )))
            if(spproj%os_mic%isthere(i, 'nframes')) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_SIZE_Z,       int(spproj%os_mic%get(i, 'nframes')))
            if(spproj%os_mic%isthere(i, 'nptcls' )) call starfile_table__setValue_int(self%starfile, SMPL_N_PTCLS,            int(spproj%os_mic%get(i, 'nptcls' )))
            if(spproj%os_mic%isthere(i, 'nmics'  )) call starfile_table__setValue_int(self%starfile, SMPL_N_MICS,             int(spproj%os_mic%get(i, 'nmics'  )))
            if(spproj%os_mic%isthere(i, 'micid'  )) call starfile_table__setValue_int(self%starfile, SMPL_MIC_ID,             int(spproj%os_mic%get(i, 'micid'  )))
            ! doubles
            if(spproj%os_mic%isthere(i, 'dfx'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSU,      real(spproj%os_mic%get(i, 'dfx') / 0.0001, dp))
            if(spproj%os_mic%isthere(i, 'dfy'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSV,      real(spproj%os_mic%get(i, 'dfy') / 0.0001, dp))
            if(spproj%os_mic%isthere(i, 'angast' )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_mic%get(i, 'angast'),       dp))
            if(spproj%os_mic%isthere(i, 'phshift')) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_PHASESHIFT,    real(spproj%os_mic%get(i, 'phshift'),      dp))
            if(spproj%os_mic%isthere(i, 'ctfres' )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_MAXRES,        real(spproj%os_mic%get(i, 'ctfres'),       dp))
            if(spproj%os_mic%isthere(i, 'icefrac')) call starfile_table__setValue_double(self%starfile,  SMPL_ICE_FRAC,          real(spproj%os_mic%get(i, 'icefrac'),      dp))
            if(spproj%os_mic%isthere(i, 'astig'  )) call starfile_table__setValue_double(self%starfile,  SMPL_ASTIGMATISM,       real(spproj%os_mic%get(i, 'astig'),        dp))
            ! strings
            if(spproj%os_mic%isthere(i, 'movie'      )) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_MOVIE_NAME,    trim(get_root_path(trim(spproj%os_mic%get_static(i, 'movie'      )))))
            if(spproj%os_mic%isthere(i, 'intg'       )) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_NAME,          trim(get_root_path(trim(spproj%os_mic%get_static(i, 'intg'       )))))
            if(spproj%os_mic%isthere(i, 'mc_starfile')) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_METADATA_NAME, trim(get_root_path(trim(spproj%os_mic%get_static(i, 'mc_starfile')))))
            if(spproj%os_mic%isthere(i, 'boxfile'    )) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_COORDINATES,   trim(get_root_path(trim(spproj%os_mic%get_static(i, 'boxfile')), 'reference_pick_extract')))
        end do

        contains

            function get_root_path ( path, addpath) result ( newpath )
                character(len=*),           intent(in) :: path
                character(len=*), optional, intent(in) :: addpath
                character(len=STDLEN)                  :: newpath, l_addpath
                newpath = ''
                l_addpath = ''
                if(present(addpath)) l_addpath = addpath
                if( len_trim(path) > 3 ) then
                    if( path(1:3) == '../' ) then
                        if(self%nicestream) then
                            ! nice stream
                            newpath = trim(self%rootpath) // '/' // trim(l_addpath) // '/' // trim(path(3:))
                        else
                            newpath = trim(path(4:))
                        end if
                    else
                       newpath = trim(path) 
                    end if
                end if
            end function get_root_path

    end subroutine starfile_set_micrographs_table

    ! export

    subroutine stream_export_micrographs( self, spproj, outdir)
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        character(len=LONGSTRLEN),  intent(in)    :: outdir
        call self%assign_optics_single(spproj)
        call self%starfile_init('micrographs.star', outdir)
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        call self%starfile_set_micrographs_table(spproj)
        call self%starfile_write_table(append = .true.)
        call self%starfile_deinit()
    end subroutine stream_export_micrographs
 
    subroutine stream_export_optics( self, spproj, outdir)
        class(starproject_stream),            intent(inout) :: self
        class(sp_project),                    intent(inout) :: spproj
        character(len=LONGSTRLEN),            intent(in)    :: outdir
        integer                                             :: ioptics
        if(params_glob%beamtilt .eq. 'yes') then
            self%use_beamtilt = .true.
        else
            self%use_beamtilt = .false.
        end if
        call self%assign_optics(spproj)
        call self%starfile_init('optics.star', outdir)
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        do ioptics = 1, spproj%os_optics%get_noris()
            call self%starfile_set_optics_group_table(spproj, ioptics)
            call self%starfile_write_table(append = .true.)
        end do
        call self%starfile_deinit()
    end subroutine stream_export_optics

    ! optics

    subroutine assign_optics_single( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer                                   :: i
        call spproj%os_mic%set_all2single('ogid', 1.0)
        call spproj%os_stk%set_all2single('ogid', 1.0)
        call spproj%os_ptcl2D%set_all2single('ogid', 1.0)
        call spproj%os_ptcl3D%set_all2single('ogid', 1.0)
        call spproj%os_optics%new(1, is_ptcl=.false.)
        call spproj%os_optics%set(1, "ogid",   1.0)
        call spproj%os_optics%set(1, "ogname", "opticsgroup1")
        do i=1,spproj%os_mic%get_noris()
            if(spproj%os_mic%get(i, 'state') .gt. 0.0 ) exit      
        end do
        call spproj%os_optics%set(1, "smpd",  spproj%os_mic%get(i, "smpd")   )
        call spproj%os_optics%set(1, "cs",    spproj%os_mic%get(i, "cs")     )
        call spproj%os_optics%set(1, "kv",    spproj%os_mic%get(i, "kv")     )
        call spproj%os_optics%set(1, "fraca", spproj%os_mic%get(i, "fraca")  )
        call spproj%os_optics%set(1, "state", 1.0)
        call spproj%os_optics%set(1, "opcx",  0.0)
        call spproj%os_optics%set(1, "opcy",  0.0)
        call spproj%os_optics%set(1, "pop",   real(spproj%os_mic%get_noris()))
    end subroutine assign_optics_single

    subroutine assign_optics( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        type(ori)                                 :: template_ori
        real                                      :: grp_info(10000, 3) ! max 10000 optics groups ! 1: centroid x, 2: centroid y, 3: population
        integer                                   :: i, ntilt, nshift
        self%shift_threshold = params_glob%tilt_thres
        call assign_tiltgroups()
        call assign_shiftgroups()
        do i=1, spproj%os_mic%get_noris()
            if(spproj%os_mic%get(i, 'state') .gt. 0.0 ) exit      
        end do
        call template_ori%new(.false.)
        call template_ori%set("smpd",  spproj%os_mic%get(i, "smpd")   )
        call template_ori%set("cs",    spproj%os_mic%get(i, "cs")     )
        call template_ori%set("kv",    spproj%os_mic%get(i, "kv")     )
        call template_ori%set("fraca", spproj%os_mic%get(i, "fraca")  )
        call template_ori%set("state", 1.0)
        call template_ori%set("pop",   0.0)
        call template_ori%set("ogid",  0.0)
        call template_ori%set("opcx",  0.0)
        call template_ori%set("opcy",  0.0)
        call template_ori%set("ogname", "opticsgroup")
        call spproj%os_optics%new(nshift, is_ptcl=.false.)
        do i = 1, nshift
            call spproj%os_optics%append(i, template_ori)
            call spproj%os_optics%set(i, "ogid", real(i + self%group_offset))
            call spproj%os_optics%set(i, "pop",  grp_info(i, 3))
            call spproj%os_optics%set(i, "opcx", grp_info(i, 1))
            call spproj%os_optics%set(i, "opcy", grp_info(i, 2))
            call spproj%os_optics%set(i, "ogname", "opticsgroup" // int2str(i + self%group_offset))
        end do
        call template_ori%kill()
        call spproj%write()

        contains

            subroutine assign_tiltgroups()
                real, allocatable :: tiltuniq(:)
                integer           :: itilt, imic
                if(self%use_beamtilt) then
                    call elim_dup(spproj%os_mic%get_all('tiltgrp'), tiltuniq)
                    ntilt = size(tiltuniq)
                    do itilt=1, ntilt 
                        do imic=1, spproj%os_mic%get_noris()
                            if(spproj%os_mic%get(imic, 'tiltgrp') == tiltuniq(itilt)) call spproj%os_mic%set(imic, 'tmpgrp', real(itilt))
                        end do
                    end do
                else
                    ntilt = 1
                    call spproj%os_mic%set_all2single('tmpgrp', 1.0)
                end if
                if(allocated(tiltuniq)) deallocate(tiltuniq)
                write(logfhandle,'(A,I8)') '>>> # TILT GROUPS ASSIGNED : ', ntilt
            end subroutine

            subroutine assign_shiftgroups()
                real,    allocatable :: tiltgrps(:), shiftxs(:), shiftys(:)
                real,    allocatable :: shifts(:,:), centroids(:,:)
                integer, allocatable :: populations(:), labels(:), indices(:)
                integer              :: itilt, imic, ishift, tiltgrppop
                write(logfhandle,'(A,F8.2)') '>>> CLUSTERING TILT GROUPS USING SHIFTS AND THRESHOLD : ', self%shift_threshold
                tiltgrps = spproj%os_mic%get_all('tmpgrp')
                shiftxs  = spproj%os_mic%get_all('shiftx')
                shiftys  = spproj%os_mic%get_all('shifty')
                nshift   = 0
                do itilt = 1, ntilt
                    write(logfhandle,'(A,I8)') '      CLUSTERING TILT GROUP : ', itilt
                    tiltgrppop = count(tiltgrps == itilt)
                    allocate(shifts(tiltgrppop, 2))
                    allocate(labels(tiltgrppop   ))
                    allocate(indices(tiltgrppop  ))
                    ishift = 1
                    do imic = 1, spproj%os_mic%get_noris()
                        if(tiltgrps(imic) == itilt) then
                            shifts(ishift, 1) = shiftxs(imic)
                            shifts(ishift, 2) = shiftys(imic)
                            indices(ishift) = imic 
                            ishift = ishift + 1
                        end if
                    end do
                    call h_clust(shifts, self%shift_threshold, labels, centroids, populations)
                    do ishift = 1, size(labels)
                        call spproj%os_mic%set(indices(ishift), 'ogid', real(labels(ishift) + nshift + self%group_offset))
                    end do
                    do ishift = 1, size(populations)
                        grp_info(ishift + nshift, 1) = centroids(ishift, 1)
                        grp_info(ishift + nshift, 2) = centroids(ishift, 2)
                        grp_info(ishift + nshift, 3) = populations(ishift)
                    end do
                    nshift = nshift + size(populations)
                    deallocate(shifts, labels, indices)
                    write(logfhandle,'(A,I8)') '        # SHIFT GROUPS ASSIGNED : ', size(populations)
                end do
                call spproj%os_mic%delete_entry('tmpgrp')
                if(allocated(populations)) deallocate(populations)
                if(allocated(labels))      deallocate(labels)
                if(allocated(indices))     deallocate(indices)
                if(allocated(shifts))      deallocate(shifts)
                if(allocated(centroids))   deallocate(centroids)
                if(allocated(tiltgrps))    deallocate(tiltgrps)
                if(allocated(shiftxs))     deallocate(shiftxs)
                if(allocated(shiftys))     deallocate(shiftys)
            end subroutine assign_shiftgroups

            ! distance threshold based yerarchical clustering
            ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
            subroutine h_clust(data_in, thresh, labels, centroids, populations)
                real,                 intent(in)  :: data_in(:,:)   ! input data, point coords
                real,                 intent(in)  :: thresh         ! threshold for class merging
                integer,              intent(out) :: labels(:)      ! labels of the elements in vec
                real,    allocatable, intent(out) :: centroids(:,:) ! centroids of the classes
                integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
                real,    allocatable :: mat(:,:)                    ! pariwise distances matrix
                logical, allocatable :: mask(:), outliers(:)
                integer :: N, i, j, cnt, ncls
                integer :: index(1), loc1(1), loc2(1)
                real    :: d
                if( size(data_in, dim = 2) .ne. 2 )then
                    THROW_HARD('Input data should be two dimensional!; h_clust')
                endif
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
                allocate(mask(N),     source = .true. )
                allocate(outliers(N), source = .false.)
                ncls = 0
                do i = 1, N
                    if( mask(i) )then ! if it's not already been clustered
                        mask(i) = .false.
                        ! find the index of the couple
                        d = minval(mat(i,:), mask)
                        index(:) = minloc(mat(i,:), mask)
                        ncls = ncls + 1
                        ! assign labels
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
                    if( any(labels== i) )then !there is a class labelled i
                        cnt = cnt + 1
                        where(labels == i) labels = cnt
                    endif
                enddo
                ! 6) recalculate centroids
                deallocate(centroids)
                ncls = maxval(labels) ! the nr of classes is maxval(labels)
                allocate(centroids(ncls,2), source = 0.)
                allocate(populations(ncls), source = 0 )
                mask = .true. ! reset
                do i = 1, ncls
                    populations(i) = count(labels == i) ! counts the nr of elements in the class
                    ! find the all the cnt member of the class and update the centroids
                    do j = 1, populations(i)
                        loc1(:) = minloc(abs(labels-i), mask)
                        mask(loc1(1)) = .false. ! do not consider this member of the class anymore
                        centroids(i,1) = centroids(i,1)+ data_in(loc1(1),1)
                        centroids(i,2) = centroids(i,2)+ data_in(loc1(1),2)
                    enddo
                    centroids(i,:) = centroids(i,:)/real(populations(i))
                enddo
            end subroutine h_clust

    end subroutine assign_optics

end module simple_starproject_stream

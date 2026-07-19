!@descr: fixed-record prior support artifact for probabilistic 3D refinement
module simple_prob_prior3D
use, intrinsic :: iso_fortran_env, only: int64
use simple_error, only: simple_exception
implicit none

public :: prior3d_writer, prior3d_reader
public :: PRIOR3D_FNAME, PRIOR3D_MAGIC, PRIOR3D_VERSION, PRIOR3D_REMAP_VERSION, PRIOR3D_NREMAP
private
#include "simple_local_flags.inc"

character(len=*), parameter :: PRIOR3D_FNAME = 'posterior_support3d.dat'
integer(int64),   parameter :: PRIOR3D_MAGIC = 731984201_int64
integer,          parameter :: PRIOR3D_VERSION = 1
integer,          parameter :: PRIOR3D_REMAP_VERSION = 1
integer,          parameter :: PRIOR3D_NREMAP = 3
integer,          parameter :: PRIOR3D_PGRP_LEN = 64
integer,          parameter :: PRIOR3D_LABEL_LEN = 32

type :: prior3d_writer
    integer :: unit = -1
    integer :: recl = 0
    integer :: kmax = 0
    integer :: nptcls = 0
    integer :: nsel = 0
    real :: support_mass = 0.
    integer, allocatable :: state(:), proj(:), inpl(:), has_sh(:)
    real,    allocatable :: dist(:), weight(:), x(:), y(:), euls(:,:)
contains
    procedure :: open  => open_writer
    procedure :: write_row
    procedure :: close => close_writer
    procedure :: kill  => kill_writer
end type prior3d_writer

type :: prior3d_reader
    integer :: unit = -1
    integer :: recl = 0
    integer :: kmax = 0
    integer :: nptcls = 0
    integer :: nstates = 0
    integer :: source_nspace = 0
    integer :: source_nspace_sub = 0
    integer :: remap_version = 0
    real :: support_mass = 0.
    character(len=PRIOR3D_PGRP_LEN) :: pgrp = ''
    character(len=PRIOR3D_LABEL_LEN) :: likelihood_transform = ''
    character(len=PRIOR3D_LABEL_LEN) :: support_convention = ''
    integer, allocatable :: pinds(:)
    integer :: nsel = 0
    integer, allocatable :: state(:), proj(:), inpl(:), has_sh(:)
    real,    allocatable :: dist(:), weight(:), x(:), y(:), euls(:,:)
contains
    procedure :: open  => open_reader
    procedure :: read_row
    procedure :: close => close_reader
    procedure :: kill  => kill_reader
end type prior3d_reader

contains

subroutine open_writer(self, binfname, nptcls, nstates, source_nspace, source_nspace_sub, kmax, pgrp, pinds)
    class(prior3d_writer), intent(inout) :: self
    character(len=*),      intent(in)    :: binfname
    integer,               intent(in)    :: nptcls, nstates, source_nspace, source_nspace_sub, kmax
    character(len=*),      intent(in)    :: pgrp
    integer,               intent(in)    :: pinds(:)
    integer :: ios, meta_unit
    character(len=512) :: metafname
    integer(int64) :: magic
    integer :: version, remap_version
    character(len=PRIOR3D_PGRP_LEN) :: pgrp_buf
    character(len=PRIOR3D_LABEL_LEN) :: likelihood_buf, convention_buf
    call self%kill
    if( nptcls < 1 .or. size(pinds) /= nptcls ) THROW_HARD('prior3d_writer: invalid particle count')
    if( kmax < 1 ) THROW_HARD('prior3d_writer: invalid support width')
    self%nptcls = nptcls
    self%kmax   = kmax
    self%nsel = 0
    allocate(self%state(kmax), self%proj(kmax), self%inpl(kmax), self%has_sh(kmax))
    allocate(self%dist(kmax), self%weight(kmax), self%x(kmax), self%y(kmax), self%euls(3,kmax))
    inquire(iolength=self%recl) self%nsel, self%support_mass, self%state, self%proj, self%inpl, self%has_sh,&
        &self%dist, self%weight, self%x, self%y, self%euls
    open(newunit=self%unit, file=trim(binfname), status='REPLACE', action='WRITE', access='DIRECT',&
        &form='UNFORMATTED', recl=self%recl, iostat=ios)
    if( ios /= 0 ) THROW_HARD('prior3d_writer: failed to open data file '//trim(binfname))
    metafname = trim(binfname)//'.meta'
    open(newunit=meta_unit, file=trim(metafname), status='REPLACE', action='WRITE', access='STREAM',&
        &form='UNFORMATTED', iostat=ios)
    if( ios /= 0 )then
        close(self%unit)
        self%unit = -1
        THROW_HARD('prior3d_writer: failed to open metadata file '//trim(metafname))
    endif
    magic = PRIOR3D_MAGIC
    version = PRIOR3D_VERSION
    remap_version = PRIOR3D_REMAP_VERSION
    pgrp_buf = ''
    if( len_trim(pgrp) > 0 ) pgrp_buf(1:min(len_trim(pgrp),len(pgrp_buf))) = &
        &pgrp(1:min(len_trim(pgrp),len(pgrp_buf)))
    likelihood_buf = 'profile_exp_minus_dist'
    convention_buf = 'dense_uniform_support'
    write(meta_unit, iostat=ios) magic, version, remap_version, nptcls, nstates, source_nspace, source_nspace_sub,&
        &kmax, pgrp_buf, likelihood_buf, convention_buf, pinds
    close(meta_unit)
    if( ios /= 0 ) THROW_HARD('prior3d_writer: failed writing metadata file '//trim(metafname))
end subroutine open_writer

subroutine write_row(self, pind, nsel, support_mass, state, proj, inpl, has_sh, dist, weight, x, y, euls)
    class(prior3d_writer), intent(inout) :: self
    integer, intent(in) :: pind, nsel, state(:), proj(:), inpl(:), has_sh(:)
    real,    intent(in) :: support_mass
    real,    intent(in) :: dist(:), weight(:), x(:), y(:), euls(:,:)
    integer :: ios, n
    if( self%unit < 0 ) THROW_HARD('prior3d_writer: data file is not open')
    if( pind < 1 .or. pind > huge(1) ) return
    n = max(0, min(self%kmax, nsel))
    self%nsel = 0
    self%support_mass = 0.
    self%state = 0
    self%proj = 0
    self%inpl = 0
    self%has_sh = 0
    self%dist = 0.
    self%weight = 0.
    self%x = 0.
    self%y = 0.
    self%euls = 0.
    self%nsel = n
    self%support_mass = max(0., min(1., support_mass))
    if( n > 0 )then
        self%state(1:n) = state(1:n)
        self%proj(1:n) = proj(1:n)
        self%inpl(1:n) = inpl(1:n)
        self%has_sh(1:n) = has_sh(1:n)
        self%dist(1:n) = dist(1:n)
        self%weight(1:n) = weight(1:n)
        self%x(1:n) = x(1:n)
        self%y(1:n) = y(1:n)
        self%euls(:,1:n) = euls(:,1:n)
    endif
    write(self%unit, rec=pind, iostat=ios) self%nsel, self%support_mass, self%state, self%proj, self%inpl, self%has_sh,&
        &self%dist, self%weight, self%x, self%y, self%euls
    if( ios /= 0 ) THROW_HARD('prior3d_writer: failed writing particle record')
end subroutine write_row

subroutine close_writer(self)
    class(prior3d_writer), intent(inout) :: self
    if( self%unit >= 0 ) close(self%unit)
    self%unit = -1
end subroutine close_writer

subroutine kill_writer(self)
    class(prior3d_writer), intent(inout) :: self
    call self%close
    if( allocated(self%state) ) deallocate(self%state)
    if( allocated(self%proj) ) deallocate(self%proj)
    if( allocated(self%inpl) ) deallocate(self%inpl)
    if( allocated(self%has_sh) ) deallocate(self%has_sh)
    if( allocated(self%dist) ) deallocate(self%dist)
    if( allocated(self%weight) ) deallocate(self%weight)
    if( allocated(self%x) ) deallocate(self%x)
    if( allocated(self%y) ) deallocate(self%y)
    if( allocated(self%euls) ) deallocate(self%euls)
    self%recl = 0
    self%kmax = 0
    self%nptcls = 0
    self%support_mass = 0.
end subroutine kill_writer

subroutine open_reader(self, binfname, ok)
    class(prior3d_reader), intent(inout) :: self
    character(len=*),      intent(in)    :: binfname
    logical,               intent(out)   :: ok
    integer :: ios, meta_unit, version
    integer(int64) :: magic
    character(len=512) :: metafname
    character(len=PRIOR3D_PGRP_LEN) :: pgrp_buf
    character(len=PRIOR3D_LABEL_LEN) :: likelihood_buf, convention_buf
    logical :: exists
    call self%kill
    ok = .false.
    metafname = trim(binfname)//'.meta'
    inquire(file=trim(binfname), exist=exists)
    if( .not. exists ) return
    inquire(file=trim(metafname), exist=exists)
    if( .not. exists ) return
    open(newunit=meta_unit, file=trim(metafname), status='OLD', action='READ', access='STREAM',&
        &form='UNFORMATTED', iostat=ios)
    if( ios /= 0 ) return
    read(meta_unit, iostat=ios) magic, version, self%remap_version, self%nptcls, self%nstates, self%source_nspace,&
        &self%source_nspace_sub, self%kmax, pgrp_buf, likelihood_buf, convention_buf
    if( ios /= 0 .or. magic /= PRIOR3D_MAGIC .or. version /= PRIOR3D_VERSION .or.&
        &self%remap_version /= PRIOR3D_REMAP_VERSION .or. trim(likelihood_buf) /= 'profile_exp_minus_dist' .or.&
        &trim(convention_buf) /= 'dense_uniform_support' .or.&
        &self%nptcls < 1 .or. self%kmax < 1 )then
        close(meta_unit)
        call self%kill
        return
    endif
    allocate(self%pinds(self%nptcls))
    read(meta_unit, iostat=ios) self%pinds
    close(meta_unit)
    if( ios /= 0 )then
        call self%kill
        return
    endif
    self%pgrp = pgrp_buf
    self%likelihood_transform = likelihood_buf
    self%support_convention = convention_buf
    allocate(self%state(self%kmax), self%proj(self%kmax), self%inpl(self%kmax), self%has_sh(self%kmax))
    allocate(self%dist(self%kmax), self%weight(self%kmax), self%x(self%kmax), self%y(self%kmax), self%euls(3,self%kmax))
    inquire(iolength=self%recl) self%nsel, self%support_mass, self%state, self%proj, self%inpl, self%has_sh,&
        &self%dist, self%weight, self%x, self%y, self%euls
    open(newunit=self%unit, file=trim(binfname), status='OLD', action='READ', access='DIRECT',&
        &form='UNFORMATTED', recl=self%recl, iostat=ios)
    if( ios /= 0 )then
        call self%kill
        return
    endif
    ok = .true.
end subroutine open_reader

subroutine read_row(self, pind, ok)
    class(prior3d_reader), intent(inout) :: self
    integer, intent(in) :: pind
    logical, intent(out) :: ok
    integer :: ios
    ok = .false.
    if( self%unit < 0 .or. pind < 1 ) return
    read(self%unit, rec=pind, iostat=ios) self%nsel, self%support_mass, self%state, self%proj, self%inpl, self%has_sh,&
        &self%dist, self%weight, self%x, self%y, self%euls
    if( ios /= 0 .or. self%nsel < 0 .or. self%nsel > self%kmax .or.&
        &self%support_mass < 0. .or. self%support_mass > 1.0001 )then
        self%nsel = 0
        return
    endif
    ok = .true.
end subroutine read_row

subroutine close_reader(self)
    class(prior3d_reader), intent(inout) :: self
    if( self%unit >= 0 ) close(self%unit)
    self%unit = -1
end subroutine close_reader

subroutine kill_reader(self)
    class(prior3d_reader), intent(inout) :: self
    call self%close
    if( allocated(self%pinds) ) deallocate(self%pinds)
    if( allocated(self%state) ) deallocate(self%state)
    if( allocated(self%proj) ) deallocate(self%proj)
    if( allocated(self%inpl) ) deallocate(self%inpl)
    if( allocated(self%has_sh) ) deallocate(self%has_sh)
    if( allocated(self%dist) ) deallocate(self%dist)
    if( allocated(self%weight) ) deallocate(self%weight)
    if( allocated(self%x) ) deallocate(self%x)
    if( allocated(self%y) ) deallocate(self%y)
    if( allocated(self%euls) ) deallocate(self%euls)
    self%recl = 0
    self%kmax = 0
    self%nptcls = 0
    self%nstates = 0
    self%source_nspace = 0
    self%source_nspace_sub = 0
    self%remap_version = 0
    self%pgrp = ''
    self%likelihood_transform = ''
    self%support_convention = ''
    self%support_mass = 0.
end subroutine kill_reader

end module simple_prob_prior3D

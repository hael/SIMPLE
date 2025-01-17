! orientation eul_prob_tab2D, used in refine3D
module simple_eul_prob_tab2D
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_polarft_corrcalc,  only: pftcc_glob
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad
use simple_eul_prob_tab,      only: eulprob_dist_switch!, eulprob_corr_switch
implicit none

public :: eul_prob_tab2D
private
#include "simple_local_flags.inc"

type :: eul_prob_tab2D
    type(ptcl_ref), allocatable :: loc_tab(:,:)   !< 2D search table (ncls x nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)   !< assignment map  (nptcls)
    integer,        allocatable :: pinds(:)       !< particle indices
    logical,        allocatable :: populated(:)   !< nonempty classes mask
    integer                     :: nptcls         !< size of pinds array
    integer                     :: ncls           !< # of classes
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! TABLE
    procedure          :: fill_tab_greedy_inpl
    procedure, private :: normalize_cls
    procedure          :: assign_cls_greedy
    ! I/O
    procedure          :: write_tab
    procedure          :: read_tab_to_glob
    procedure          :: write_assignment
    procedure          :: read_assignment
    ! DESTRUCTOR
    procedure          :: kill
end type eul_prob_tab2D

contains

    ! CONSTRUCTORS

    subroutine new( self, pinds )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: pinds(:)
        integer, allocatable :: pops(:)
        integer :: i, iptcl, icls
        real    :: x
        call self%kill
        call seed_rnd
        self%nptcls = size(pinds)
        self%ncls   = params_glob%ncls
        allocate(self%loc_tab(self%ncls,self%nptcls), self%assgn_map(self%nptcls),self%pinds(self%nptcls))
        ! Particles
        !$omp parallel do default(shared) private(i,iptcl,icls) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            self%pinds(i) = pinds(i)
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%istate = 1
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
            do icls = 1,self%ncls
                self%loc_tab(icls,i)%pind   = iptcl
                self%loc_tab(icls,i)%istate = 1
                self%loc_tab(icls,i)%iproj  = icls
                self%loc_tab(icls,i)%inpl   = 0
                self%loc_tab(icls,i)%dist   = huge(x)
                self%loc_tab(icls,i)%x      = 0.
                self%loc_tab(icls,i)%y      = 0.
                self%loc_tab(icls,i)%has_sh = .false.
            end do
        end do
        !$omp end parallel do
        ! Classes
        if( build_glob%spproj%os_cls2D%get_noris() == 0 )then
            if( build_glob%spproj_field%isthere('class') )then
                call build_glob%spproj%os_ptcl2D%get_pops(pops, 'class', maxn=self%ncls)
            else
                allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
            endif
        else
            if( build_glob%spproj_field%isthere('class') )then
                if( build_glob%spproj%os_cls2D%get_noris() /= self%ncls )then
                    ! to be able to restart after having run cleanup with fewer classes
                    allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
                else
                    pops = nint(build_glob%spproj%os_cls2D%get_all('pop'))
                    where( pops < 2 ) pops = 0 ! ignoring classes with one particle
                endif
            else
                allocate(pops(self%ncls), source=MINCLSPOPLIM+1)
            endif
        endif
        if( all(pops == 0) ) THROW_HARD('All class pops cannot be zero!')
        self%populated = pops > 0
    end subroutine new

    ! TABLE

    ! fill the probability table
    subroutine fill_tab_greedy_inpl( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(pftcc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        real    :: scores(pftcc_glob%get_nrots())
        real    :: lims(2,2), lims_init(2,2), cxy(3), best_score
        integer :: i, iptcl, ithr, irot, icls, best_rot
        if( params_glob%l_doshift )then
            ! search objects
            lims(:,1)      = -params_glob%trs
            lims(:,2)      =  params_glob%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
                    &maxits=params_glob%maxits_sh, opt_angle=.true.)
            end do
            ! search
            !$omp parallel do default(shared) private(i,iptcl,ithr,icls,irot,best_rot,best_score,scores,cxy)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                do icls = 1, self%ncls
                    if( .not.self%populated(icls) ) cycle
                    call pftcc_glob%gencorrs(icls, iptcl, scores)
                    irot     = maxloc(scores, dim=1)
                    best_rot = irot
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if( irot == 0 )then
                        best_score = scores(best_rot)
                        cxy(2:3)   = 0.
                    else
                        best_rot   = irot
                        best_score = cxy(1)
                    endif
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(best_score)
                    self%loc_tab(icls,i)%inpl   = best_rot
                    self%loc_tab(icls,i)%x      = cxy(2)
                    self%loc_tab(icls,i)%y      = cxy(3)
                    self%loc_tab(icls,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,icls,irot,scores)&
            !$omp proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                do icls = 1, self%ncls
                    if( .not.self%populated(icls) ) cycle
                    call pftcc_glob%gencorrs(icls, iptcl, scores)
                    irot = maxloc(scores, dim=1)
                    self%loc_tab(icls,i)%dist = eulprob_dist_switch(scores(irot))
                    self%loc_tab(icls,i)%inpl = irot
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_tab_greedy_inpl

    subroutine normalize_cls( self )
        class(eul_prob_tab2D), intent(inout) :: self
        ! PLACEHOLDER
    end subroutine normalize_cls

    subroutine assign_cls_greedy( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: i, icls
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,icls)
        do i = 1,self%nptcls
            icls = minloc(self%loc_tab(:,i)%dist, dim=1, mask=self%populated)
            self%assgn_map(i) = self%loc_tab(icls,i)
        enddo
        !$omp end parallel do
    end subroutine assign_cls_greedy

    ! FILE I/O

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_tab( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        character(len=*),      intent(in) :: binfname
        integer :: funit, addr, io_stat, file_header(2)
        file_header(1) = self%ncls
        file_header(2) = self%nptcls
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit,pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_tab

    ! read the partition-wise dist value binary file to global object's table
    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        character(len=*),      intent(in)    :: binfname
        type(ptcl_ref), allocatable :: mat_loc(:,:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, ncls_loc, i_loc, i_glob
        if( file_exists(trim(binfname)) )then
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab2D; read_tab_to_glob; file: '//trim(binfname), io_stat)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        ! reading header and the ncls/nptcls in this partition file
        read(unit=funit,pos=1) file_header
        ncls_loc = file_header(1)
        nptcls_loc = file_header(2)
        if( ncls_loc .ne. params_glob%ncls ) THROW_HARD( 'NCLS should be the same in this partition file!' )
        allocate(mat_loc(ncls_loc, nptcls_loc))
        ! read partition information
        addr = sizeof(file_header) + 1
        read(unit=funit,pos=addr) mat_loc
        call fclose(funit)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_glob = 1, self%nptcls
            do i_loc = 1, nptcls_loc
                if( mat_loc(1,i_loc)%pind == self%loc_tab(1,i_glob)%pind )then
                    self%loc_tab(:,i_glob) = mat_loc(:,i_loc)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine read_tab_to_glob

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        character(len=*),      intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,trim(binfname),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        character(len=*),      intent(in)    :: binfname
        type(ptcl_ref), allocatable :: assgn_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(trim(binfname)) )then
            THROW_HARD('file '//trim(binfname)//' does not exists!')
        else
            call fopen(funit,trim(binfname),access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab2D; read_assignment; file: '//trim(binfname), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit,pos=headsz + 1) assgn_glob
        call fclose(funit)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( self%assgn_map(i_loc)%pind == assgn_glob(i_glob)%pind )then
                    self%assgn_map(i_loc) = assgn_glob(i_glob)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine read_assignment

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab2D), intent(inout) :: self
        self%ncls   = 0
        self%nptcls = 0
        if( allocated(self%loc_tab)   ) deallocate(self%loc_tab)
        if( allocated(self%assgn_map) ) deallocate(self%assgn_map)
        if( allocated(self%pinds)     ) deallocate(self%pinds)
        if( allocated(self%populated) ) deallocate(self%populated)
    end subroutine kill

end module simple_eul_prob_tab2D
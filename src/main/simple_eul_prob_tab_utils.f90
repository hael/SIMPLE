!@descr: shared utility routines for probabilistic alignment tables
module simple_eul_prob_tab_utils
use simple_math,      only: rotmat2d
use simple_type_defs, only: ptcl_ref
implicit none

public :: build_pind_lookup, materialize_seed_shift, read_seed_shift_table, write_seed_shift_table
private

contains

    subroutine build_pind_lookup( glob_pinds, loc_pinds, pind2glob, max_pind )
        integer,              intent(in)  :: glob_pinds(:), loc_pinds(:)
        integer, allocatable, intent(out) :: pind2glob(:)
        integer,              intent(out) :: max_pind
        integer :: i, pind
        max_pind = max(maxval(glob_pinds), maxval(loc_pinds))
        if( max_pind < 1 )then
            allocate(pind2glob(0))
            return
        endif
        allocate(pind2glob(max_pind), source=0)
        do i = 1, size(glob_pinds)
            pind = glob_pinds(i)
            if( pind > 0 .and. pind <= max_pind ) pind2glob(pind) = i
        enddo
    end subroutine build_pind_lookup

    subroutine materialize_seed_shift( assgn, seed_shift, seed_has_sh, l_doshift, seed_nrots )
        type(ptcl_ref),       intent(inout) :: assgn
        real,                 intent(in)    :: seed_shift(2)
        logical,              intent(in)    :: seed_has_sh, l_doshift
        integer,              intent(in)    :: seed_nrots
        real    :: rotmat(2,2), rot_xy(2)
        integer :: irot
        if( .not. l_doshift ) return
        if( assgn%has_sh    ) return
        if( .not. seed_has_sh ) return
        irot = assgn%inpl
        if( irot < 1 .or. irot > seed_nrots .or. seed_nrots < 1 )then
            assgn%x      = 0.
            assgn%y      = 0.
            assgn%has_sh = .true.
            return
        endif
        call rotmat2d(real(irot - 1) * 360. / real(seed_nrots), rotmat)
        rot_xy(1) = seed_shift(1) * rotmat(1,1) + seed_shift(2) * rotmat(2,1)
        rot_xy(2) = seed_shift(1) * rotmat(1,2) + seed_shift(2) * rotmat(2,2)
        assgn%x      = rot_xy(1)
        assgn%y      = rot_xy(2)
        assgn%has_sh = .true.
    end subroutine materialize_seed_shift

    subroutine write_seed_shift_table( funit, addr, seed_nrots, seed_shifts, seed_has_sh )
        integer, intent(in)    :: funit
        integer, intent(inout) :: addr
        integer, intent(in)    :: seed_nrots
        real,    intent(in)    :: seed_shifts(:,:)
        logical, intent(in)    :: seed_has_sh(:)
        write(funit, pos=addr) seed_nrots
        addr = addr + sizeof(seed_nrots)
        write(funit, pos=addr) seed_shifts
        addr = addr + sizeof(seed_shifts)
        write(funit, pos=addr) seed_has_sh
        addr = addr + sizeof(seed_has_sh)
    end subroutine write_seed_shift_table

    subroutine read_seed_shift_table( funit, addr, seed_nrots, seed_shifts, seed_has_sh )
        integer, intent(in)    :: funit
        integer, intent(inout) :: addr
        integer, intent(out)   :: seed_nrots
        real,    intent(out)   :: seed_shifts(:,:)
        logical, intent(out)   :: seed_has_sh(:)
        read(funit, pos=addr) seed_nrots
        addr = addr + sizeof(seed_nrots)
        read(funit, pos=addr) seed_shifts
        addr = addr + sizeof(seed_shifts)
        read(funit, pos=addr) seed_has_sh
        addr = addr + sizeof(seed_has_sh)
    end subroutine read_seed_shift_table

end module simple_eul_prob_tab_utils

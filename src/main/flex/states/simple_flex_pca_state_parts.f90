!@descr: flex_pca state-reconstruction parts codec.
!! Owns the per-round weight table and both backends' part-artifact names.
module simple_flex_pca_state_parts
use simple_core_module_api, only: fclose, file_exists, fileiochk, fopen, int2str, int2str_pad, mrc_ext, real2str, &
    &simple_exception, simple_rename, string
use simple_parameters,         only: parameters
use simple_flex_pca_artifacts, only: FLEX_PCA_PART_MAGIC
use simple_flex_pca_rounds,    only: flex_pca_rounds
implicit none

public :: write_state_weights_round, read_state_weights_round
public :: flex_state_part_fbody, flex_pca_rho_part_name
public :: flex_pcg_state_raw_fname, flex_pcg_state_provenance
private
#include "simple_local_flags.inc"

character(len=*), parameter :: WEIGHTS_FNAME = 'flex_pca_round_weights.bin'

contains

    !> Per-round state weights, written with the global pinds so a worker matches its rows by index.
    !! Rewritten before every state round since each round uses a different weight table.
    subroutine write_state_weights_round( pinds, weights, nptcls, nstates, split_eo )
        integer, intent(in) :: pinds(:), nptcls, nstates
        real,    intent(in) :: weights(:,:)
        !! .true. when this round accumulates even and odd into separate reconstructors; a worker
        !! cannot infer it from params%stage (every state round arrives as PCA_STAGE_STATES)
        logical, intent(in) :: split_eo
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat, eo_flag
        fname     = string(WEIGHTS_FNAME)
        tmp_fname = fname//'.tmp'
        eo_flag   = merge(1, 0, split_eo)
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_state_weights_round; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, nptcls, nstates, eo_flag
        call fileiochk('write_state_weights_round; header', io_stat)
        write(funit, iostat=io_stat) pinds(1:nptcls)
        call fileiochk('write_state_weights_round; pinds', io_stat)
        write(funit, iostat=io_stat) weights(1:nptcls,1:nstates)
        call fileiochk('write_state_weights_round; weights', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
    end subroutine write_state_weights_round

    !> Worker side: return the weight rows for this part's pinds, in this part's order.
    subroutine read_state_weights_round( my_pinds, my_nptcls, weights_out, nstates, split_eo )
        integer,           intent(in)  :: my_pinds(:), my_nptcls
        real, allocatable, intent(out) :: weights_out(:,:)
        integer,           intent(out) :: nstates
        logical,           intent(out) :: split_eo !< see write_state_weights_round
        type(string) :: fname
        integer, allocatable :: gpinds(:)
        real,    allocatable :: gw(:,:)
        integer :: funit, io_stat, magic, gn, i, j, hit, eo_flag
        logical :: l_sorted
        fname = string(WEIGHTS_FNAME)
        if( .not. file_exists(fname) ) THROW_HARD('flex_pca worker found no '//WEIGHTS_FNAME)
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_state_weights_round; open', io_stat)
        read(funit, iostat=io_stat) magic, gn, nstates, eo_flag
        call fileiochk('read_state_weights_round; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad round-weights magic')
        split_eo = eo_flag == 1
        allocate(gpinds(gn), gw(gn,nstates))
        read(funit, iostat=io_stat) gpinds
        call fileiochk('read_state_weights_round; pinds', io_stat)
        read(funit, iostat=io_stat) gw
        call fileiochk('read_state_weights_round; weights', io_stat)
        call fclose(funit)
        allocate(weights_out(my_nptcls,nstates), source=0.)
        ! both lists are ascending project rows, so a merge walk replaces the O(N_local x N_global)
        ! scan; falls back to the scan if either list is unordered
        l_sorted = .true.
        do i = 2, my_nptcls
            if( my_pinds(i) <= my_pinds(i-1) ) l_sorted = .false.
        end do
        do j = 2, gn
            if( gpinds(j) <= gpinds(j-1) ) l_sorted = .false.
        end do
        j = 1
        do i = 1, my_nptcls
            hit = 0
            if( l_sorted )then
                do while( j <= gn )
                    if( gpinds(j) >= my_pinds(i) ) exit
                    j = j + 1
                end do
                if( j <= gn )then
                    if( gpinds(j) == my_pinds(i) ) hit = j
                endif
            else
                do j = 1, gn
                    if( gpinds(j) == my_pinds(i) )then
                        hit = j
                        exit
                    endif
                end do
            endif
            if( hit == 0 ) THROW_HARD('flex_pca worker particle absent from the master weight table')
            weights_out(i,:) = gw(hit,:)
        end do
        deallocate(gpinds, gw)
        call fname%kill
    end subroutine read_state_weights_round

    !> Part-file body for one worker's partial reconstruction of one state; eo=0 is the even/single
    !! accumulator, eo=1 the odd one (a split round writes both per state).
    function flex_state_part_fbody( params, rounds, part, state, eo ) result( fbody )
        class(parameters), intent(in) :: params
        class(flex_pca_rounds), intent(in) :: rounds
        integer,           intent(in) :: part, state, eo
        type(string) :: fbody
        fbody = string('flex_pca_statepart')//int2str_pad(part,max(1,params%numlen))// &
            &'_'//int2str_pad(state,2)
        if( eo == 1 ) fbody = fbody//'_o'
        fbody = rounds%part_path(fbody%to_char())
    end function flex_state_part_fbody

    !> rho companion of a state part: same directory, 'rho_' on the file name only
    function flex_pca_rho_part_name( pf ) result( fn )
        type(string), intent(in) :: pf
        type(string) :: fn
        character(len=:), allocatable :: c
        integer :: k
        c = pf%to_char()
        k = index(c, '/', back=.true.)
        fn = string(c(1:k))//'rho_'//c(k+1:)//MRC_EXT
    end function flex_pca_rho_part_name

    !> raw artifact of one worker's (state, half) accumulation, under the part directory
    function flex_pcg_state_raw_fname( params, rounds, part, state, eo ) result( fname )
        class(parameters), intent(in) :: params
        class(flex_pca_rounds), intent(in) :: rounds
        integer,           intent(in) :: part, state, eo
        type(string) :: fname
        character(len=2) :: half
        half = '_e'
        if( eo == 1 ) half = '_o'
        fname = rounds%part_path('flex_pca_pcgraw_part'//int2str_pad(part,max(1,params%numlen))// &
            &'_'//int2str_pad(state,2)//half//'.bin')
    end function flex_pcg_state_raw_fname

    !> identity every raw artifact of one run carries; the master refuses a part that disagrees
    function flex_pcg_state_provenance( params, box_rec, smpd_rec ) result( provenance )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: box_rec
        real,              intent(in) :: smpd_rec
        character(len=256) :: provenance
        provenance = 'flexpcg-v1|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|box='//trim(int2str(params%box))// &
            &'|smpd='//trim(real2str(params%smpd))//'|box_rec='//trim(int2str(box_rec))// &
            &'|smpd_rec='//trim(real2str(smpd_rec))//'|mskdiam='//trim(real2str(params%mskdiam))// &
            &'|ctf='//trim(params%ctf)
    end function flex_pcg_state_provenance

end module simple_flex_pca_state_parts

module simple_commanders_prob
use simple_commanders_api
use simple_eul_prob_tab_utils, only: read_seed_shift_table
use simple_pftc_srch_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_prob_tab
  contains
    procedure :: execute      => exec_prob_tab
end type commander_prob_tab

type, extends(commander_base) :: commander_prob_tab_neigh
    contains
        procedure :: execute      => exec_prob_tab_neigh
end type commander_prob_tab_neigh

type, extends(commander_base) :: commander_prob_align
  contains
    procedure :: execute      => exec_prob_align
end type commander_prob_align

type, extends(commander_base) :: commander_prob_align_neigh
    contains
        procedure :: execute      => exec_prob_align_neigh
end type commander_prob_align_neigh

type, extends(commander_base) :: commander_prob_tab2D
  contains
    procedure :: execute      => exec_prob_tab2D
end type commander_prob_tab2D

type, extends(commander_base) :: commander_prob_align2D
  contains
    procedure :: execute      => exec_prob_align2D
end type commander_prob_align2D

contains

    subroutine exec_prob_tab( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils,    only: read_reprojection_model
        use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_eul_prob_tab,            only: eul_prob_tab
        class(commander_prob_tab), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab)       :: eulprob_obj_part
        integer :: nptcls, batchsz_max, nbatches, ibatch, batch_start, batch_end, batchsz
        integer, allocatable :: batches(:,:)
        logical :: l_state_only
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        if( nptcls < 1 ) THROW_HARD('exec_prob_tab selected no particles')
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call read_reprojection_model(params, build, batchsz_max)
        call prep_sigmas_objfun(params, build, .false.)
        call alloc_ptcl_imgs( params, build, tmp_imgs, tmp_imgs_pad, batchsz_max )
        call build%pftc%memoize_refs(eulspace=build%eulspace)
        ! Fill the partition table in matcher-sized batches to cap particle PFT memo memory.
        l_state_only = str_has_substr(params%refine, 'prob_state')
        call eulprob_obj_part%new(params, build, pinds, state_only=l_state_only)
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles3D(params, build, batchsz, pinds(batch_start:batch_end), tmp_imgs, tmp_imgs_pad)
            if( l_state_only )then
                call eulprob_obj_part%fill_tab_state_only_range(batch_start, batch_end)
            else
                call eulprob_obj_part%fill_tab_range(batch_start, batch_end)
            endif
        end do
        if( l_state_only )then
            call eulprob_obj_part%write_state_tab(fname)
        else
            call eulprob_obj_part%write_tab(fname)
        endif
        call eulprob_obj_part%kill
        if( allocated(batches) ) deallocate(batches)
        call build%pftc%kill
        call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_tab_neigh( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils,    only: read_reprojection_model
        use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        class(commander_prob_tab_neigh), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname, fname_batch
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab_neigh) :: eulprob_obj_part_neigh
        integer :: nptcls, batchsz_shadow, nbatches_shadow
        logical :: l_batch_prob_neigh, l_compare_batch_prob_neigh
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_neigh requires prior particle sampling (in exec_prob_align)')
        endif
        if( nptcls < 1 ) THROW_HARD('exec_prob_tab_neigh selected no particles')
        l_batch_prob_neigh = trim(params%prob_neigh_mode) == 'shc' .or. trim(params%prob_neigh_mode) == 'snhc'
        l_compare_batch_prob_neigh = .not. l_batch_prob_neigh
        fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(params%part,params%numlen)//'.dat'
        if( l_compare_batch_prob_neigh )then
            fname_batch = fname//'_batch_probe'
            call run_prob_tab_neigh_full(fname)
            call run_prob_tab_neigh_batch(fname_batch, batchsz_shadow, nbatches_shadow)
            call compare_prob_neigh_sparse_tables(fname, fname_batch, params, batchsz_shadow, nbatches_shadow)
            if( file_exists(fname_batch) ) call del_file(fname_batch)
            call fname_batch%kill
        else
            call run_prob_tab_neigh_batch(fname, batchsz_shadow, nbatches_shadow)
        endif
        call fname%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab_neigh'))
        call simple_end('**** SIMPLE_PROB_TAB_NEIGH NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine prepare_prob_neigh_workspace(batchsz_here)
            integer, intent(in) :: batchsz_here
            call read_reprojection_model(params, build, batchsz_here)
            call prep_sigmas_objfun(params, build, .false.)
            call alloc_ptcl_imgs(params, build, tmp_imgs, tmp_imgs_pad, batchsz_here)
            call build%pftc%memoize_refs(eulspace=build%eulspace)
        end subroutine prepare_prob_neigh_workspace

        subroutine cleanup_prob_neigh_workspace
            call build%pftc%kill
            call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        end subroutine cleanup_prob_neigh_workspace

        subroutine run_prob_tab_neigh_full(outfname)
            class(string), intent(in) :: outfname
            call prepare_prob_neigh_workspace(nptcls)
            call build_batch_particles3D(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
            call eulprob_obj_part_neigh%new_neigh(params, build, pinds)
            call eulprob_obj_part_neigh%fill_tab
            call eulprob_obj_part_neigh%write_tab(outfname)
            call eulprob_obj_part_neigh%kill
            call cleanup_prob_neigh_workspace
        end subroutine run_prob_tab_neigh_full

        subroutine run_prob_tab_neigh_batch(outfname, batchsz_used, nbatches_used)
            class(string), intent(in)  :: outfname
            integer,       intent(out) :: batchsz_used, nbatches_used
            integer, allocatable :: batches(:,:)
            integer :: ibatch, batch_start, batch_end, batchsz
            batchsz_used   = min(nptcls, max(1, params%nthr * BATCHTHRSZ))
            nbatches_used  = ceiling(real(nptcls) / real(batchsz_used))
            batches        = split_nobjs_even(nptcls, nbatches_used)
            batchsz_used   = maxval(batches(:,2) - batches(:,1) + 1)
            call prepare_prob_neigh_workspace(batchsz_used)
            call eulprob_obj_part_neigh%new_neigh(params, build, pinds)
            do ibatch = 1, nbatches_used
                batch_start = batches(ibatch,1)
                batch_end   = batches(ibatch,2)
                batchsz     = batch_end - batch_start + 1
                call build_batch_particles3D(params, build, batchsz, pinds(batch_start:batch_end), tmp_imgs, tmp_imgs_pad)
                call eulprob_obj_part_neigh%fill_tab_range(batch_start, batch_end)
            enddo
            call eulprob_obj_part_neigh%write_tab(outfname)
            call eulprob_obj_part_neigh%kill
            if( allocated(batches) ) deallocate(batches)
            call cleanup_prob_neigh_workspace
        end subroutine run_prob_tab_neigh_batch

    end subroutine exec_prob_tab_neigh

    subroutine compare_prob_neigh_sparse_tables(fname_full, fname_batch, params, batchsz_used, nbatches_used)
        class(string),     intent(in) :: fname_full, fname_batch
        class(parameters), intent(in) :: params
        integer,           intent(in) :: batchsz_used, nbatches_used

        type :: neigh_sparse_table
            integer :: nrefs = 0, nptcls = 0, nnz = 0, seed_nrots = 0
            integer,        allocatable :: pinds(:), counts(:), refs(:), offsets(:)
            real,           allocatable :: seed_shifts(:,:)
            logical,        allocatable :: seed_has_sh(:)
            type(ptcl_ref), allocatable :: tab(:)
        end type neigh_sparse_table

        type(neigh_sparse_table) :: full_tab, batch_tab
        real, parameter :: DIST_TOL = 1.e-5, SHIFT_TOL = 1.e-4
        integer :: i, pos_full, pos_batch, best_full, best_batch
        integer :: pinds_mismatches, count_mismatches, particles_refset_diff, particles_any_diff
        integer :: only_full, only_batch, common_refs, common_dist_mismatches
        integer :: common_inpl_mismatches, common_shift_mismatches, common_meta_mismatches
        integer :: best_ref_mismatches, best_inpl_mismatches, best_shift_mismatches, best_meta_mismatches
        integer :: seed_has_mismatches, seed_shift_mismatches
        integer :: first_pind, first_full_count, first_batch_count, first_full_ref, first_batch_ref
        integer :: first_full_inpl, first_batch_inpl
        real    :: dist_abs, dist_abs_sum, dist_abs_sq_sum, dist_abs_max
        real    :: shift_abs, shift_abs_sum, shift_abs_max
        real    :: best_dist_abs, best_dist_abs_sum, best_dist_abs_max
        real    :: seed_shift_abs, seed_shift_abs_sum, seed_shift_abs_max
        real    :: first_full_dist, first_batch_dist
        logical :: particle_diff

        call read_neigh_sparse_table(fname_full, full_tab)
        call read_neigh_sparse_table(fname_batch, batch_tab)
        call init_stats()
        write(logfhandle,'(A)') '>>> PROB_NEIGH BATCH_COMPARE BEGIN'
        write(logfhandle,'(A,I0,A,A,A,I0,A,I0)') &
            &'>>> PROB_NEIGH BATCH_COMPARE part ', params%part, ' mode ', trim(params%prob_neigh_mode), &
            &' batchsz ', batchsz_used, ' nbatches ', nbatches_used
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
            &'>>> PROB_NEIGH BATCH_COMPARE headers full(nrefs,nptcls,nnz)=(', full_tab%nrefs, ',', &
            &full_tab%nptcls, ',', full_tab%nnz, ') batch=(', batch_tab%nrefs, ',', batch_tab%nptcls, ',', &
            &batch_tab%nnz, ')'
        if( full_tab%nrefs /= batch_tab%nrefs .or. full_tab%nptcls /= batch_tab%nptcls )then
            write(logfhandle,'(A)') '>>> PROB_NEIGH BATCH_COMPARE ABORT header mismatch'
            call cleanup_tables()
            call flush(logfhandle)
            return
        endif
        do i = 1, full_tab%nptcls
            particle_diff = .false.
            if( full_tab%pinds(i) /= batch_tab%pinds(i) )then
                pinds_mismatches = pinds_mismatches + 1
                particle_diff = .true.
            endif
            if( full_tab%counts(i) /= batch_tab%counts(i) )then
                count_mismatches = count_mismatches + 1
                particle_diff = .true.
            endif
            if( full_tab%seed_has_sh(i) .neqv. batch_tab%seed_has_sh(i) )then
                seed_has_mismatches = seed_has_mismatches + 1
                particle_diff = .true.
            endif
            seed_shift_abs = maxval(abs(full_tab%seed_shifts(:,i) - batch_tab%seed_shifts(:,i)))
            seed_shift_abs_max = max(seed_shift_abs_max, seed_shift_abs)
            seed_shift_abs_sum = seed_shift_abs_sum + seed_shift_abs
            if( seed_shift_abs > SHIFT_TOL )then
                seed_shift_mismatches = seed_shift_mismatches + 1
                particle_diff = .true.
            endif
            call compare_ref_sets_for_particle(i, particle_diff)
            call compare_best_for_particle(i, particle_diff)
            if( particle_diff )then
                particles_any_diff = particles_any_diff + 1
                if( first_pind == 0 ) call capture_first_mismatch(i)
            endif
        enddo
        call print_stats()
        call cleanup_tables()
        call flush(logfhandle)

    contains

        subroutine init_stats()
            pinds_mismatches       = 0
            count_mismatches       = 0
            particles_refset_diff  = 0
            particles_any_diff     = 0
            only_full              = 0
            only_batch             = 0
            common_refs            = 0
            common_dist_mismatches = 0
            common_inpl_mismatches = 0
            common_shift_mismatches = 0
            common_meta_mismatches = 0
            best_ref_mismatches    = 0
            best_inpl_mismatches   = 0
            best_shift_mismatches  = 0
            best_meta_mismatches   = 0
            seed_has_mismatches    = 0
            seed_shift_mismatches  = 0
            dist_abs_sum           = 0.
            dist_abs_sq_sum        = 0.
            dist_abs_max           = 0.
            shift_abs_sum          = 0.
            shift_abs_max          = 0.
            best_dist_abs_sum      = 0.
            best_dist_abs_max      = 0.
            seed_shift_abs_sum     = 0.
            seed_shift_abs_max     = 0.
            first_pind             = 0
            first_full_count       = 0
            first_batch_count      = 0
            first_full_ref         = 0
            first_batch_ref        = 0
            first_full_inpl        = 0
            first_batch_inpl       = 0
            first_full_dist        = 0.
            first_batch_dist       = 0.
        end subroutine init_stats

        subroutine read_neigh_sparse_table(fname, table)
            class(string),             intent(in)    :: fname
            type(neigh_sparse_table),  intent(inout) :: table
            integer :: funit, addr, io_stat, file_header(3), pos
            if( .not. file_exists(fname) ) THROW_HARD('missing prob_neigh sparse table: '//fname%to_char())
            call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('compare_prob_neigh_sparse_tables; file: '//fname%to_char(), io_stat)
            read(unit=funit,pos=1) file_header
            table%nrefs  = file_header(1)
            table%nptcls = file_header(2)
            table%nnz    = file_header(3)
            allocate(table%pinds(table%nptcls), table%counts(table%nptcls))
            allocate(table%seed_shifts(2,table%nptcls), table%seed_has_sh(table%nptcls))
            allocate(table%refs(table%nnz), table%tab(table%nnz), table%offsets(table%nptcls+1))
            addr = sizeof(file_header) + 1
            read(funit, pos=addr) table%pinds
            addr = addr + int(sizeof(table%pinds))
            call read_seed_shift_table(funit, addr, table%seed_nrots, table%seed_shifts, table%seed_has_sh)
            read(funit, pos=addr) table%counts
            addr = addr + int(sizeof(table%counts))
            read(funit, pos=addr) table%refs
            addr = addr + int(sizeof(table%refs))
            read(funit, pos=addr) table%tab
            call fclose(funit)
            if( any(table%counts < 0) .or. sum(table%counts) /= table%nnz )&
                &THROW_HARD('bad sparse counts in compare_prob_neigh_sparse_tables')
            table%offsets(1) = 1
            do pos = 1, table%nptcls
                table%offsets(pos+1) = table%offsets(pos) + table%counts(pos)
            enddo
        end subroutine read_neigh_sparse_table

        subroutine compare_ref_sets_for_particle(iptcl_loc, particle_diff)
            integer, intent(in)    :: iptcl_loc
            logical, intent(inout) :: particle_diff
            integer :: k_loc, ref_loc
            logical :: refset_diff
            refset_diff = full_tab%counts(iptcl_loc) /= batch_tab%counts(iptcl_loc)
            do k_loc = full_tab%offsets(iptcl_loc), full_tab%offsets(iptcl_loc+1)-1
                ref_loc = full_tab%refs(k_loc)
                pos_batch = find_ref_pos(batch_tab, iptcl_loc, ref_loc)
                if( pos_batch == 0 )then
                    only_full = only_full + 1
                    refset_diff = .true.
                else
                    common_refs = common_refs + 1
                    call compare_common_entry(k_loc, pos_batch, particle_diff)
                endif
            enddo
            do k_loc = batch_tab%offsets(iptcl_loc), batch_tab%offsets(iptcl_loc+1)-1
                ref_loc = batch_tab%refs(k_loc)
                pos_full = find_ref_pos(full_tab, iptcl_loc, ref_loc)
                if( pos_full == 0 )then
                    only_batch = only_batch + 1
                    refset_diff = .true.
                endif
            enddo
            if( refset_diff )then
                particles_refset_diff = particles_refset_diff + 1
                particle_diff = .true.
            endif
        end subroutine compare_ref_sets_for_particle

        subroutine compare_common_entry(pos_full_loc, pos_batch_loc, particle_diff)
            integer, intent(in)    :: pos_full_loc, pos_batch_loc
            logical, intent(inout) :: particle_diff
            dist_abs = abs(full_tab%tab(pos_full_loc)%dist - batch_tab%tab(pos_batch_loc)%dist)
            dist_abs_sum    = dist_abs_sum + dist_abs
            dist_abs_sq_sum = dist_abs_sq_sum + dist_abs * dist_abs
            dist_abs_max    = max(dist_abs_max, dist_abs)
            shift_abs = max(abs(full_tab%tab(pos_full_loc)%x - batch_tab%tab(pos_batch_loc)%x),&
                &abs(full_tab%tab(pos_full_loc)%y - batch_tab%tab(pos_batch_loc)%y))
            shift_abs_sum = shift_abs_sum + shift_abs
            shift_abs_max = max(shift_abs_max, shift_abs)
            if( dist_abs > DIST_TOL )then
                common_dist_mismatches = common_dist_mismatches + 1
                particle_diff = .true.
            endif
            if( full_tab%tab(pos_full_loc)%inpl /= batch_tab%tab(pos_batch_loc)%inpl )then
                common_inpl_mismatches = common_inpl_mismatches + 1
                particle_diff = .true.
            endif
            if( shift_abs > SHIFT_TOL )then
                common_shift_mismatches = common_shift_mismatches + 1
                particle_diff = .true.
            endif
            if( full_tab%tab(pos_full_loc)%pind /= batch_tab%tab(pos_batch_loc)%pind .or.&
                &full_tab%tab(pos_full_loc)%istate /= batch_tab%tab(pos_batch_loc)%istate .or.&
                &full_tab%tab(pos_full_loc)%iproj /= batch_tab%tab(pos_batch_loc)%iproj .or.&
                &full_tab%tab(pos_full_loc)%has_sh .neqv. batch_tab%tab(pos_batch_loc)%has_sh )then
                common_meta_mismatches = common_meta_mismatches + 1
                particle_diff = .true.
            endif
        end subroutine compare_common_entry

        subroutine compare_best_for_particle(iptcl_loc, particle_diff)
            integer, intent(in)    :: iptcl_loc
            logical, intent(inout) :: particle_diff
            best_full  = best_ref_pos(full_tab,  iptcl_loc)
            best_batch = best_ref_pos(batch_tab, iptcl_loc)
            if( best_full == 0 .or. best_batch == 0 )then
                if( best_full /= best_batch )then
                    best_ref_mismatches = best_ref_mismatches + 1
                    particle_diff = .true.
                endif
                return
            endif
            if( full_tab%refs(best_full) /= batch_tab%refs(best_batch) )then
                best_ref_mismatches = best_ref_mismatches + 1
                particle_diff = .true.
            endif
            best_dist_abs = abs(full_tab%tab(best_full)%dist - batch_tab%tab(best_batch)%dist)
            best_dist_abs_sum = best_dist_abs_sum + best_dist_abs
            best_dist_abs_max = max(best_dist_abs_max, best_dist_abs)
            if( full_tab%tab(best_full)%inpl /= batch_tab%tab(best_batch)%inpl )then
                best_inpl_mismatches = best_inpl_mismatches + 1
                particle_diff = .true.
            endif
            shift_abs = max(abs(full_tab%tab(best_full)%x - batch_tab%tab(best_batch)%x),&
                &abs(full_tab%tab(best_full)%y - batch_tab%tab(best_batch)%y))
            if( shift_abs > SHIFT_TOL )then
                best_shift_mismatches = best_shift_mismatches + 1
                particle_diff = .true.
            endif
            if( full_tab%tab(best_full)%istate /= batch_tab%tab(best_batch)%istate .or.&
                &full_tab%tab(best_full)%iproj /= batch_tab%tab(best_batch)%iproj .or.&
                &full_tab%tab(best_full)%has_sh .neqv. batch_tab%tab(best_batch)%has_sh )then
                best_meta_mismatches = best_meta_mismatches + 1
                particle_diff = .true.
            endif
        end subroutine compare_best_for_particle

        integer function find_ref_pos(table, iptcl_loc, ref_loc) result(pos_found)
            type(neigh_sparse_table), intent(in) :: table
            integer,                  intent(in) :: iptcl_loc, ref_loc
            integer :: pos_loc
            pos_found = 0
            do pos_loc = table%offsets(iptcl_loc), table%offsets(iptcl_loc+1)-1
                if( table%refs(pos_loc) == ref_loc )then
                    pos_found = pos_loc
                    return
                endif
            enddo
        end function find_ref_pos

        integer function best_ref_pos(table, iptcl_loc) result(pos_best)
            type(neigh_sparse_table), intent(in) :: table
            integer,                  intent(in) :: iptcl_loc
            integer :: pos_loc
            real    :: best_dist
            pos_best = 0
            best_dist = huge(best_dist)
            do pos_loc = table%offsets(iptcl_loc), table%offsets(iptcl_loc+1)-1
                if( table%tab(pos_loc)%dist < best_dist )then
                    best_dist = table%tab(pos_loc)%dist
                    pos_best  = pos_loc
                endif
            enddo
        end function best_ref_pos

        subroutine capture_first_mismatch(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            first_pind        = full_tab%pinds(iptcl_loc)
            first_full_count  = full_tab%counts(iptcl_loc)
            first_batch_count = batch_tab%counts(iptcl_loc)
            best_full         = best_ref_pos(full_tab,  iptcl_loc)
            best_batch        = best_ref_pos(batch_tab, iptcl_loc)
            if( best_full > 0 )then
                first_full_ref  = full_tab%refs(best_full)
                first_full_inpl = full_tab%tab(best_full)%inpl
                first_full_dist = full_tab%tab(best_full)%dist
            endif
            if( best_batch > 0 )then
                first_batch_ref  = batch_tab%refs(best_batch)
                first_batch_inpl = batch_tab%tab(best_batch)%inpl
                first_batch_dist = batch_tab%tab(best_batch)%dist
            endif
        end subroutine capture_first_mismatch

        subroutine print_stats()
            real :: common_avg, common_rms, shift_avg, best_avg, seed_shift_avg
            common_avg = 0.
            common_rms = 0.
            shift_avg  = 0.
            if( common_refs > 0 )then
                common_avg = dist_abs_sum / real(common_refs)
                common_rms = sqrt(dist_abs_sq_sum / real(common_refs))
                shift_avg  = shift_abs_sum / real(common_refs)
            endif
            best_avg = 0.
            if( full_tab%nptcls > 0 ) best_avg = best_dist_abs_sum / real(full_tab%nptcls)
            seed_shift_avg = 0.
            if( full_tab%nptcls > 0 ) seed_shift_avg = seed_shift_abs_sum / real(full_tab%nptcls)
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') &
                &'>>> PROB_NEIGH BATCH_COMPARE particles any_diff=', particles_any_diff, &
                &' pind_mismatch=', pinds_mismatches, ' count_mismatch=', count_mismatches, &
                &' refset_diff=', particles_refset_diff
            write(logfhandle,'(A,I0,A,I0,A,I0)') &
                &'>>> PROB_NEIGH BATCH_COMPARE refs common=', common_refs, &
                &' only_full=', only_full, ' only_batch=', only_batch
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
                &'>>> PROB_NEIGH BATCH_COMPARE common dist_mismatch=', common_dist_mismatches, &
                &' inpl_mismatch=', common_inpl_mismatches, ' shift_mismatch=', common_shift_mismatches, &
                &' meta_mismatch=', common_meta_mismatches, ' dist_abs_avg=', common_avg, &
                &' dist_abs_rms=', common_rms, ' dist_abs_max=', dist_abs_max
            write(logfhandle,'(A,ES12.5,A,ES12.5)') &
                &'>>> PROB_NEIGH BATCH_COMPARE common shift_abs_avg=', shift_avg, &
                &' shift_abs_max=', shift_abs_max
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,ES12.5,A,ES12.5)') &
                &'>>> PROB_NEIGH BATCH_COMPARE best ref_mismatch=', best_ref_mismatches, &
                &' inpl_mismatch=', best_inpl_mismatches, ' shift_mismatch=', best_shift_mismatches, &
                &' meta_mismatch=', best_meta_mismatches, ' dist_abs_avg=', best_avg, &
                &' dist_abs_max=', best_dist_abs_max
            write(logfhandle,'(A,I0,A,I0,A,ES12.5,A,ES12.5)') &
                &'>>> PROB_NEIGH BATCH_COMPARE seed has_mismatch=', seed_has_mismatches, &
                &' shift_mismatch=', seed_shift_mismatches, ' shift_abs_avg=', seed_shift_avg, &
                &' shift_abs_max=', seed_shift_abs_max
            if( first_pind > 0 )then
                write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                    &'>>> PROB_NEIGH BATCH_COMPARE first_mismatch pind=', first_pind, &
                    &' counts(full,batch)=', first_full_count, ',', first_batch_count, &
                    &' best_ref(full,batch)=', first_full_ref, ',', first_batch_ref, &
                    &' best_inpl(full,batch)=', first_full_inpl, ',', first_batch_inpl
                write(logfhandle,'(A,ES12.5,A,ES12.5)') &
                    &'>>> PROB_NEIGH BATCH_COMPARE first_mismatch best_dist(full,batch)=', &
                    &first_full_dist, ',', first_batch_dist
            endif
            write(logfhandle,'(A)') '>>> PROB_NEIGH BATCH_COMPARE END'
        end subroutine print_stats

        subroutine cleanup_tables()
            call cleanup_table(full_tab)
            call cleanup_table(batch_tab)
        end subroutine cleanup_tables

        subroutine cleanup_table(table)
            type(neigh_sparse_table), intent(inout) :: table
            if( allocated(table%pinds)        ) deallocate(table%pinds)
            if( allocated(table%counts)       ) deallocate(table%counts)
            if( allocated(table%refs)         ) deallocate(table%refs)
            if( allocated(table%offsets)      ) deallocate(table%offsets)
            if( allocated(table%seed_shifts)  ) deallocate(table%seed_shifts)
            if( allocated(table%seed_has_sh)  ) deallocate(table%seed_has_sh)
            if( allocated(table%tab)          ) deallocate(table%tab)
        end subroutine cleanup_table

    end subroutine compare_prob_neigh_sparse_tables

    subroutine exec_prob_align( self, cline )
        use simple_eul_prob_tab,            only: eul_prob_tab
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(commander_prob_tab) :: xprob_tab
        type(eul_prob_tab)       :: eulprob_obj_glob
        type(cmdline)            :: cline_prob_tab
        type(qsys_env)           :: qenv
        type(chash)              :: job_descr
        integer :: nptcls, ipart
        logical :: l_state_only
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sampled incremented
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype)
        call cleanup_prob_align_outputs(params, .false.)
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute(cline_prob_tab)
        else
            ! setup the environment for distributed execution
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! Build the global table only after worker tables are complete.  Keeping it
        ! live while workers build dense partition tables roughly doubles peak RSS.
        l_state_only = str_has_substr(params%refine, 'prob_state')
        call eulprob_obj_glob%new(params, build, pinds, state_only=l_state_only)
        ! reading corrs from all parts
        if( l_state_only )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_state_tab(fname)
            enddo
            call eulprob_obj_glob%state_assign
        else
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_tab_to_glob(fname)
            enddo
            call eulprob_obj_glob%ref_assign
        endif
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        call eulprob_obj_glob%kill
        ! cleanup
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align

    subroutine exec_prob_align_neigh( self, cline )
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(commander_prob_align_neigh), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,           allocatable :: pinds(:)
        type(string)                   :: fname
        type(builder)                  :: build
        type(parameters)               :: params
        type(commander_prob_tab_neigh) :: xprob_tab_neigh
        type(eul_prob_tab_neigh)       :: eulprob_obj_glob_neigh
        type(cmdline)                  :: cline_prob_tab
        type(qsys_env)                 :: qenv
        type(chash)                    :: job_descr
        integer :: nptcls
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call build%spproj%write_segment_inside(params%oritype)
        call cleanup_prob_align_outputs(params, .true.)
        ! Global object only needs sampled-set maps before reading partition sparse tables.
        ! The neighborhood scoring itself is performed in each prob_tab_neigh partition job.
        call eulprob_obj_glob_neigh%new_neigh(params, build, pinds)
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab_neigh')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab_neigh%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        call eulprob_obj_glob_neigh%read_tabs_to_glob(string(DIST_FBODY)//'_neigh_', params%nparts, params%numlen)
        call eulprob_obj_glob_neigh%ref_assign
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob_neigh%write_assignment(fname)
        call eulprob_obj_glob_neigh%kill
        ! cleanup
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align_neigh'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align_neigh

    subroutine exec_prob_tab2D( self, cline )
        use simple_matcher_smpl_and_lplims, only: set_bp_range2D
        use simple_strategy2D_matcher,  only: set_b_p_ptrs2D, &
                                              ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
        use simple_matcher_pftc_prep,      only: prep_pftc4align2D
        use simple_matcher_ptcl_batch,  only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D
        use simple_imgarr_utils,        only: alloc_imgarr
        use simple_classaverager,       only: cavger_new, cavger_read_all, cavger_kill
        use simple_eul_prob_tab2D,      only: eul_prob_tab2D
        class(commander_prob_tab2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab2D)     :: eulprob_obj_part
        real    :: frac_srch_space
        integer :: nptcls, batchsz_max, nbatches, ibatch, batch_start, batch_end, batchsz
        integer, allocatable :: batches(:,:)
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        frac_srch_space  = build%spproj_field%get_avg('frac')
        call set_bp_range2D(params, build, cline, params%which_iter, frac_srch_space)
        ! reproduce particle sampling from exec_prob_align2D
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab2D requires prior particle sampling (in exec_prob_align2D)')
        endif
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        call set_b_p_ptrs2D(params, build)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        call alloc_imgarr(batchsz_max, [params%box, params%box, 1], params%smpd, ptcl_imgs)
        ! mirror cluster2D_exec reference setup
        call cavger_new(params, build)
        if( .not. cline%defined('refs') ) THROW_HARD('exec_prob_tab2D requires refs on the command line')
        call cavger_read_all
        call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, batchsz_max, params%which_iter, .false.)
        ! Fill the partition table in matcher-sized batches to cap polar FT memo memory.
        call eulprob_obj_part%new(params, build, pinds)
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles2D(params, build, batchsz, pinds(batch_start:batch_end),&
                &ptcl_imgs(1:batchsz), ptcl_match_imgs, ptcl_match_imgs_pad)
            call eulprob_obj_part%fill_tab_range(batch_start, batch_end)
        end do
        ! write the 2D probability table
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part%write_tab(fname)
        call eulprob_obj_part%kill
        if( allocated(batches) ) deallocate(batches)
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        call cavger_kill
        call build%pftc%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_tab2D'))
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

    subroutine exec_prob_align2D( self, cline )
        use simple_eul_prob_tab2D,          only: eul_prob_tab2D
        use simple_strategy2D_matcher,      only: set_b_p_ptrs2D
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4update2D
        use simple_builder,                 only: builder
        class(commander_prob_align2D), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(string)               :: fname
        type(builder)              :: build
        type(parameters)           :: params
        type(commander_prob_tab2D) :: xprob_tab2D
        type(eul_prob_tab2D)       :: eulprob_obj_glob
        type(cmdline)              :: cline_prob_tab
        type(qsys_env)             :: qenv
        type(chash)                :: job_descr
        integer :: nptcls, ipart
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call set_b_p_ptrs2D(params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! Mirror the 3D workflow: sampled-update is active from the first stage onward.
        ! In probabilistic mode the sampled subset is reused within the current iteration
        ! by prob_tab2D/cluster2D_exec, but it is redrawn on later iterations.
        call sample_ptcls4update2D(params, build, [params%fromp,params%top], params%l_update_frac, nptcls, pinds)
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> PROB_ALIGN2D: sampled ', nptcls, ' particles over ', params%nparts, ' part(s)'
        call flush(logfhandle)
        ! write sampling to project
        call build%spproj%write_segment_inside(params%oritype)
        ! build the global prob table (nclasses x nptcls)
        call eulprob_obj_glob%new(params, build, pinds)
        ! generate partition-wise dist tables
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab2D')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab2D%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: prob_tab2D workers completed; merging partition tables'
        call flush(logfhandle)
        ! merge all partition tables into global
        do ipart = 1, params%nparts
            fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
            call eulprob_obj_glob%read_tab_to_glob(fname)
        end do
        ! global probabilistic class assignment
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: running global probabilistic assignment'
        call flush(logfhandle)
        call eulprob_obj_glob%ref_assign
        ! write assignment to file
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        write(logfhandle,'(A,A)') '>>> PROB_ALIGN2D: writing assignment ', fname%to_char()
        call flush(logfhandle)
        call eulprob_obj_glob%write_assignment(fname)
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: assignment written'
        call flush(logfhandle)
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_align2D'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align2D

    subroutine cleanup_prob_align_outputs( params, neigh )
        class(parameters), intent(in) :: params
        logical,           intent(in) :: neigh
        type(string) :: fname
        integer :: ipart
        if( neigh )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(ipart,params%numlen)//'.dat'
                call del_file(fname)
            end do
        else
            call del_files(DIST_FBODY, params%nparts, ext='.dat', numlen=params%numlen)
        endif
        call del_file(string(ASSIGNMENT_FBODY)//'.dat')
        call fname%kill
    end subroutine cleanup_prob_align_outputs

end module simple_commanders_prob

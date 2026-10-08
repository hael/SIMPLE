!@descr: flex_pca project gateway: the one owner of what the run reads from and writes into the project
!!
!! Selection and validation of the run's particles (a worker's partition list, the master's state>0
!! rows), the even/odd repair the master persists for its workers, the canonical sigma2 load with
!! the master/worker pin, and the deliveries that mutate the run's project copy: the state weight set
!! registered in the out segment, the hard state labels written into ptcl3D, and the state maps
!! reconstructed from that set and registered as vol and FSC entries.
module simple_flex_pca_project_gateway
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api, only: fclose, file_exists, fileiochk, fopen, logfhandle, nlines, simple_exception, string, tiny
use simple_builder,            only: builder
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_sp_project,         only: sp_project
use simple_sigma2_files,       only: load_sigma2_groups
use simple_flex_pca_rounds,    only: flex_pca_rounds
use simple_state_weight_set,   only: state_weight_set, STATE_WEIGHTS_KIND_PARTITION
implicit none
private
#include "simple_local_flags.inc"

public :: validate_covariance_inputs, ensure_canonical_sigma_state, load_and_validate_sigma
public :: write_state_weight_set, write_discrete_state_project, register_embedding_artifact
public :: deliver_state_maps

character(len=*), parameter :: SIGMA_STATE_FNAME = 'flex_pca_sigma_state.txt'

contains

    !> Read a part's particle-index list (one integer per line, as written by the master's
    !! flex_pca_plan_partitions through arr2txtfile).
    subroutine read_pind_list( fname, pinds, nptcls )
        character(len=*),     intent(in)  :: fname
        integer, allocatable, intent(out) :: pinds(:)
        integer,              intent(out) :: nptcls
        type(string) :: fn
        integer :: funit, io_stat, i, n
        fn = trim(fname)
        if( .not. file_exists(fn) ) THROW_HARD('flex_pca worker: particle list not found: '//trim(fname))
        n = nlines(fn)
        if( n < 1 ) THROW_HARD('flex_pca worker: empty particle list: '//trim(fname))
        allocate(pinds(n))
        call fopen(funit, file=fn, action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pind_list', io_stat)
        do i = 1, n
            read(funit,*) pinds(i)
        end do
        call fclose(funit)
        call fn%kill
        nptcls = n
    end subroutine read_pind_list

    subroutine validate_covariance_inputs( params, l_vol1_explicit, l_pindfile, build, pinds, nptcls, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters),       intent(inout) :: params
        logical,                intent(in)    :: l_vol1_explicit, l_pindfile
        type(builder),          intent(inout) :: build
        integer, allocatable,   intent(out)   :: pinds(:)
        integer,                intent(out)   :: nptcls
        integer :: q, i, cnt
        integer, allocatable :: sel(:)
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('flex_pca requires oritype=ptcl3D')
        if( .not. l_vol1_explicit )then
            THROW_HARD('flex_pca requires a consensus mean map: pass vol1 or register one in the project out segment')
        endif
        ! a run writes its state labels into its project copy, so a rerun must start from the original
        ! project (a subset is selected with pindfile=, never with the labels)
        if( build%spproj_field%get_n('state') /= 1 ) THROW_HARD('flex_pca works on one population but the project carries several state labels (a delivered copy): rerun from the original project, selecting particles with pindfile= if needed')
        if( l_pindfile )then
            ! a distributed worker takes the master's partition as it was planned: its own
            ! particle-index list, never a re-derivation of the selection from fromp/top
            call read_pind_list(params%pindfile%to_char(), pinds, nptcls)
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_PCA (WORKER) ',nptcls,' particles from ', &
                &params%pindfile%to_char()
            call flush(logfhandle)
        else
            ! not sample4rec: its updatecnt > 0 condition belongs to trailing reconstruction
            allocate(sel(max(0, params%top - params%fromp + 1)))
            cnt = 0
            do i = params%fromp, params%top
                if( build%spproj_field%get_state(i) > 0 )then
                    cnt = cnt + 1; sel(cnt) = i
                endif
            end do
            nptcls = cnt
            if( allocated(pinds) ) deallocate(pinds)
            allocate(pinds(nptcls), source=sel(1:nptcls))
            deallocate(sel)
        endif
        if( nptcls < 100 ) THROW_HARD('flex_pca requires at least 100 active particles')
        if( build%spproj%os_ptcl2D%get_noris() > 0 .and. &
            &build%spproj%os_ptcl2D%get_noris() /= build%spproj%os_ptcl3D%get_noris() ) &
            &THROW_HARD('flex_pca requires matching ptcl2D and ptcl3D rows')
        ! eo split by index parity when the project's is degenerate. MASTER-ONLY and must be persisted: a
        ! worker counts over its own range, so it would silently build a different -- equally valid -- split.
        if( rounds%is_worker() )then
            if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 1 .or. &
                &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 1 )then
                THROW_HARD('flex_pca worker sees a degenerate eo split; the master did not persist its repair')
            endif
        else if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA assigning alternating even/odd halfsets (project eo was degenerate)'
            call build%spproj_field%partition_eo
            if( rounds%is_master() )then
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                write(logfhandle,'(A)') '>>> FLEX_PCA persisted the repaired eo split for the workers'
            endif
        endif
        if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 ) &
            &THROW_HARD('flex_pca requires populated even and odd halfsets')
    end subroutine validate_covariance_inputs

    !> The canonical sigma2 state, validated (canonical_sigma2_consumable) or seeded from particle
    !! power (simple_sigma2_bootstrap) before any part command line is generated; the master and
    !! shared-memory roles only, workers consume what it committed. Per-stack (group) sigma2 needs
    !! particles in both halves of every stack: a 2D selection deselects whole stacks, and the
    !! canonical reduce then refuses ("empty even/odd half"), so a selection with such a stack falls
    !! back to the pooled spectrum here and sigma_est=global never has to be typed for a subset.
    subroutine ensure_canonical_sigma_state( params, build, cline )
        use simple_core_module_api,  only: OBJFUN_EUCLID
        use simple_sigma2_bootstrap, only: ensure_sigma2_for_iteration
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer, allocatable :: cnt(:,:)
        integer :: iptcl, ngroups, g, e, nempty
        logical :: l_bootstrapped
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        if( .not. params%l_sigma_glob )then
            ngroups = 0
            do iptcl = 1, params%nptcls
                if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                ngroups = max(ngroups, build%spproj_field%get_int(iptcl, 'stkind'))
            enddo
            if( ngroups > 0 )then
                allocate(cnt(0:1,ngroups), source=0)
                do iptcl = 1, params%nptcls
                    if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                    g = build%spproj_field%get_int(iptcl, 'stkind')
                    e = build%spproj_field%get_eo(iptcl)
                    if( g < 1 .or. g > ngroups .or. e < 0 .or. e > 1 ) cycle
                    cnt(e,g) = cnt(e,g) + 1
                enddo
                nempty = 0
                do g = 1, ngroups
                    if( cnt(0,g) == 0 .or. cnt(1,g) == 0 ) nempty = nempty + 1
                enddo
                if( nempty > 0 )then
                    if( sum(cnt(0,:)) == 0 .or. sum(cnt(1,:)) == 0 )then
                        write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: the project carries no even/odd assignment, so &
                            &per-stack sigma2 halves cannot be formed; using the pooled (global) noise spectrum &
                            &(sigma_est=global)'
                    else
                        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA SIGMA: ', nempty, ' of ', ngroups, &
                            &' sigma2 groups (stacks) have no particles in one half after the selection; &
                            &using the pooled (global) noise spectrum (sigma_est=global)'
                    endif
                    params%sigma_est   = 'global'
                    params%l_sigma_glob = .true.
                    call cline%set('sigma_est', 'global')
                endif
                deallocate(cnt)
            endif
        endif
        call ensure_sigma2_for_iteration(cline, params%projfile, 1, params%box, params%smpd, params%l_sigma_glob, &
            &'FLEX_PCA SIGMA', l_bootstrapped, sigma_est=merge('global', 'group ', params%l_sigma_glob))
        if( l_bootstrapped ) call build%spproj%read_segment('projinfo', params%projfile)
    end subroutine ensure_canonical_sigma_state

    subroutine load_and_validate_sigma( params, build, cline, pinds, loaded , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters),       intent(inout) :: params
        type(builder),          intent(inout) :: build
        class(cmdline),         intent(inout) :: cline
        integer,                intent(in)    :: pinds(:)
        logical,                intent(out)   :: loaded
        integer :: i, k, iptcl, noris, fromp_save, top_save
        ! The sigma table is allocated on params%fromp:top. A distributed worker's particle list is the
        ! master's partition (pindfile), not the fromp/top split of the project, so widen the range to
        ! every particle for the load (the canonical reader fills rows for state>0 within the range) and
        ! restore after. The unit fallback below has always done the same.
        noris      = build%spproj_field%get_noris()
        fromp_save = params%fromp
        top_save   = params%top
        params%fromp = 1
        params%top   = noris
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, build%spproj_field, loaded)
        params%fromp = fromp_save
        params%top   = top_save
        if( .not. loaded )then
            ! gen_fplane4rec follows norm_noise_taper_edge_pad_fft, so a unit spectrum is the correct fallback
            ! Construct the sigma object through its own API (the legacy loader's idiom: widen
            ! fromp/top to every particle so the table covers the whole project, restore after):
            ! esig%new also registers the table with the polar calculator; a bare allocate of the
            ! component inside a never-constructed object leaves the consumers with an
            ! inconsistent object
            params%fromp = 1
            params%top   = noris
            call build%esig%new(params, build%pftc, string('flex_pca_unit_sigma2.dat'), params%box)
            params%fromp = fromp_save
            params%top   = top_save
            build%esig%sigma2_noise = 1.0
            loaded = .true.
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: using unit white-noise spectrum after background normalization'
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: provide sigma2_group*.star for coloured experimental noise'
        endif
        do i = 1, size(pinds)
            iptcl = pinds(i)
            do k = lbound(build%esig%sigma2_noise,1), ubound(build%esig%sigma2_noise,1)
                if( .not. ieee_is_finite(build%esig%sigma2_noise(k,iptcl)) .or. &
                    &build%esig%sigma2_noise(k,iptcl) <= TINY )then
                    THROW_HARD('flex_pca found a nonpositive or nonfinite sigma2 value')
                endif
            end do
        end do
        params%ml_reg   = 'yes'
        params%l_ml_reg = .true.
        call cline%set('ml_reg','yes')
        ! Pin the sigma2 decision across the master/worker boundary (see write_sigma_state).
        if( rounds%is_master() )then
            call write_sigma_state(loaded)
        else if( rounds%is_worker() )then
            call check_sigma_state(loaded)
        endif
    end subroutine load_and_validate_sigma

    !> Whether the master resolved real sigma2 spectra or fell back to the unit white-noise spectrum.
    !! Discovery is re-run independently in every worker, and a worker that resolves it differently
    !! whitens against a different noise model -- which changes every inner product without any error
    !! being raised, because both outcomes are individually legitimate. Pin it and make workers assert.
    subroutine write_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_sigma_state; open', io_stat)
        write(funit,'(L1)') loaded
        call fclose(funit)
    end subroutine write_sigma_state

    subroutine check_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        logical :: loaded_master
        if( .not. file_exists(string(SIGMA_STATE_FNAME)) )then
            THROW_HARD('flex_pca worker found no flex_pca_sigma_state.txt from the master')
        endif
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('check_sigma_state; open', io_stat)
        read(funit,'(L1)') loaded_master
        call fclose(funit)
        if( loaded_master .neqv. loaded )then
            THROW_HARD('flex_pca worker resolved sigma2 differently from the master; the parts would be whitened against different noise models')
        endif
    end subroutine check_sigma_state

    !>  The delivered weight table as the project's state weight set (simple_state_weight_set): one
    !!  file per state over every project row (zero outside the selection) with the hard labels, a
    !!  manifest published last, registered in the out segment of the run's own project copy.
    subroutine write_state_weight_set( params, build, pinds, weights, labels, l_merged )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        integer,          intent(in)    :: pinds(:), labels(:)
        real,             intent(in)    :: weights(:,:)
        logical,          intent(in)    :: l_merged
        type(state_weight_set) :: weight_set
        call weight_set%publish(build%spproj, build%spproj_field, params%projfile, pinds, weights, labels, &
            &merge('flex_pca_merged', 'flex_pca       ', l_merged))
        write(logfhandle,'(A,I0,A,I0,A,I0,A,A)') '>>> FLEX_PCA STATE WEIGHT SET PUBLISHED: generation ', &
            &weight_set%get_generation(), ', ', weight_set%get_nptcls(), ' particles x ', weight_set%get_nstates(), &
            &' states, kind ', merge('partition', 'kernel   ', weight_set%get_kind() == STATE_WEIGHTS_KIND_PARTITION)
        call flush(logfhandle)
        call weight_set%kill
    end subroutine write_state_weight_set

    !>  The run's embedding artifact as the project's out-segment entry flex_embedding (absolute path), so
    !!  consumers find it through the project; written with the out segment. No artifact, no entry.
    subroutine register_embedding_artifact( params, build, fname )
        use simple_syslib, only: simple_abspath
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        character(len=*),  intent(in)    :: fname
        type(string) :: path
        integer :: ind
        if( .not. file_exists(fname) ) return
        path = simple_abspath(fname)
        call build%spproj%add_entry2os_out('flex_embedding', ind)
        call build%spproj%os_out%set(ind, 'imgkind',        'flex_embedding')
        call build%spproj%os_out%set(ind, 'flex_embedding', path%to_char())
        call build%spproj%write_segment_inside('out', params%projfile)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA EMBEDDING ARTIFACT REGISTERED: ', path%to_char()
        call path%kill
    end subroutine register_embedding_artifact

    !>  Write the hard state assignment INTO the run's own project: ptcl3D/state carries each embedded
    !!  particle's label, 0 elsewhere. Judge the clusters independently of the kernel-weighted backend with
    !!      simple_exec prg=reconstruct3D projfile=<projfile> nstates=<nstates>
    !!  mkdir=yes already gave the master a private copy of the project, so this rewrites that copy and
    !!  never the project the user pointed at. MUTATES the live field, so it must run after every stage
    !!  that reads the input particle selection, and on the master only.
    subroutine write_discrete_state_project( spproj, pinds, labels, nstates, projfile )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:), labels(:), nstates
        type(string),     intent(in)    :: projfile
        logical, allocatable :: assigned(:)
        integer :: i, iptcl, state, nptcls, nexcluded
        if( size(pinds) < 1 .or. size(labels) /= size(pinds) .or. nstates < 2 ) &
            &THROW_HARD('invalid flex_pca discrete-state assignment')
        if( len_trim(projfile%to_char()) == 0 ) THROW_HARD('flex_pca discrete-state project file is empty')
        nptcls = spproj%os_ptcl3D%get_noris()
        ! validate BEFORE mutating: this overwrites the live project field rather than a private copy,
        ! so a mid-loop abort would leave the input selection half-replaced
        allocate(assigned(nptcls), source=.false.)
        nexcluded = 0
        do i = 1,size(pinds)
            iptcl = pinds(i)
            state = labels(i)
            if( iptcl < 1 .or. iptcl > nptcls ) THROW_HARD('flex_pca discrete-state particle index outside project')
            if( assigned(iptcl) ) THROW_HARD('duplicate particle in flex_pca discrete-state assignment')
            if( state > nstates ) THROW_HARD('flex_pca discrete-state label outside state range')
            if( state < 1 ) nexcluded = nexcluded + 1
            assigned(iptcl) = .true.
        end do
        do iptcl = 1,nptcls
            call spproj%os_ptcl3D%set_state(iptcl,0)
        end do
        do i = 1,size(pinds)
            if( labels(i) >= 1 ) call spproj%os_ptcl3D%set_state(pinds(i),labels(i))
        end do
        ! ptcl3D only, as refine3D does for nstates>1: a state INDEX carries no ptcl2D meaning, and only
        ! selection (0/1) is mirrored across the two segments
        call spproj%write_segment_inside('ptcl3D', projfile)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA HARD STATES WRITTEN TO: ',projfile%to_char()
        do state = 1,nstates
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA DISCRETE-STATE state=',state, &
                &' population=',count(labels==state)
        end do
        if( nexcluded > 0 ) write(logfhandle,'(A,I0)') &
            &'>>> FLEX_PCA DISCRETE-STATE unassigned particles left at state=0: ',nexcluded
        write(logfhandle,'(A,A,A,I0)') '>>> RECONSTRUCT WITH: simple_exec prg=reconstruct3D projfile=', &
            &projfile%to_char(),' nstates=',nstates
        call flush(logfhandle)
        deallocate(assigned)
    end subroutine write_discrete_state_project

    !> The delivered state maps (box decision D9: the native box, as reconstruct3D does): reconstruct3D
    !! with m_estimator=flex from the weight set just published, in this process and directory, then
    !! registered as the project's ordinary vol and FSC entries of every state (ruling 5); the vol and FSC
    !! entries of states it does not deliver, and an earlier release's vol_flex entries, are removed
    subroutine deliver_state_maps( params, build, nstates )
        use simple_commanders_rec,  only: commander_rec3D
        use simple_refine3D_fnames, only: refine3D_state_vol_fname, refine3D_fsc_fname
        use simple_syslib,          only: simple_abspath
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nstates
        type(commander_rec3D) :: xrec3D
        type(cmdline)         :: cline_rec
        type(string)          :: vol_fname, fsc_fname
        integer :: state, nstates_out
        call cline_rec%set('prg',         'reconstruct3D')
        call cline_rec%set('projfile',    params%projfile%to_char())
        call cline_rec%set('mkdir',       'no')
        call cline_rec%set('oritype',     'ptcl3D')
        call cline_rec%set('nstates',     nstates)
        call cline_rec%set('m_estimator', 'flex')
        call cline_rec%set('rec_backend', trim(params%rec_states_backend))
        call cline_rec%set('mskdiam',     params%mskdiam)
        call cline_rec%set('pgrp',        trim(params%pgrp))
        call cline_rec%set('objfun',      trim(params%objfun))
        call cline_rec%set('sigma_est',   trim(params%sigma_est))
        call cline_rec%set('nthr',        params%nthr)
        ! the run's solver controls carry over to its delivered maps
        call cline_rec%set('maxits_pcg',  params%maxits_pcg)
        call cline_rec%set('rtol',        params%rtol)
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA DELIVERED STATE MAPS: reconstruct3D m_estimator=flex nstates=', nstates, &
            &' from the published state weight set'
        call flush(logfhandle)
        call xrec3D%execute(cline_rec)
        call cline_rec%kill
        ! the reconstruction updated the project file: register on top of what it holds
        call build%spproj%read_segment('out', params%projfile)
        nstates_out = 0
        if( build%spproj%os_out%get_noris() > 0 ) nstates_out = build%spproj%os_out%get_n('state')
        do state = 1, max(nstates, nstates_out)
            call build%spproj%remove_entry_from_osout('vol_flex', state)
            vol_fname = refine3D_state_vol_fname(state)
            fsc_fname = refine3D_fsc_fname(state)
            if( state > nstates .or. .not. file_exists(vol_fname) .or. .not. file_exists(fsc_fname) )then
                call build%spproj%remove_entry_from_osout('vol', state)
                call build%spproj%remove_entry_from_osout('fsc', state)
                if( state <= nstates ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA state ', state, &
                    &' has no delivered map; not registered'
                cycle
            endif
            vol_fname = simple_abspath(vol_fname)
            call build%spproj%add_vol2os_out(vol_fname, params%smpd, state, 'vol', box=params%box)
            call build%spproj%add_fsc2os_out(fsc_fname, state, params%box)
        end do
        call build%spproj%write_segment_inside('out', params%projfile)
        call vol_fname%kill
        call fsc_fname%kill
    end subroutine deliver_state_maps

end module simple_flex_pca_project_gateway

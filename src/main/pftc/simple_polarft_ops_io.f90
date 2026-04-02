!@descr: submodule for parallel I/O and polar->Cartesian conversion
submodule (simple_polarft_calc) simple_polarft_ops_io
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine vol_pad2ref_pfts_write_range(self, vol_pad, eulspace, state, iproj_from, iproj_to, mask, tmpl_fname)
        use simple_projector, only: projector
        class(polarft_calc), intent(inout) :: self
        class(projector),    intent(in)    :: vol_pad
        class(oris),         intent(inout) :: eulspace
        integer,             intent(in)    :: state, iproj_from, iproj_to
        logical,             intent(in)    :: mask(:)
        class(string),       intent(in)    :: tmpl_fname
        integer :: funit_e, funit_o, nprojs_write, iref_from, iref_to, i
        if( .not. self%existence ) THROW_HARD('polarft_calc does not exist; vol_pad2ref_pfts_write_range')
        if( state < 1 .or. state > self%p_ptr%nstates )then
            write(logfhandle,*) 'state, nstates: ', state, self%p_ptr%nstates
            THROW_HARD('state out of range; vol_pad2ref_pfts_write_range')
        endif
        if( iproj_from < 1 .or. iproj_to > self%p_ptr%nspace .or. iproj_from > iproj_to )then
            write(logfhandle,*) 'iproj_from, iproj_to, nspace: ', iproj_from, iproj_to, self%p_ptr%nspace
            THROW_HARD('projection range out of bounds; vol_pad2ref_pfts_write_range')
        endif
        if( size(mask) < self%kfromto(2) )then
            write(logfhandle,*) 'size(mask), kfromto(2): ', size(mask), self%kfromto(2)
            THROW_HARD('mask too short for requested k-range; vol_pad2ref_pfts_write_range')
        endif
        iref_from = (state - 1) * self%p_ptr%nspace + iproj_from
        iref_to   = (state - 1) * self%p_ptr%nspace + iproj_to
        call vol_pad2ref_pfts(self, vol_pad, eulspace, state, .true.,  mask)
        call vol_pad2ref_pfts(self, vol_pad, eulspace, state, .false., mask)
        nprojs_write = iproj_to - iproj_from + 1
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_even'//BIN_EXT, funit_e)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_odd'//BIN_EXT,  funit_o)
        !$omp parallel do default(shared) private(i) num_threads(2) schedule(static)
        do i = 1, 2
            select case(i)
                case(1)
                    call write_ref_pft_range_local(funit_e, self%pfts_refs_even, nprojs_write)
                case(2)
                    call write_ref_pft_range_local(funit_o, self%pfts_refs_odd,  nprojs_write)
            end select
        end do
        !$omp end parallel do
        call fclose(funit_e)
        call fclose(funit_o)

    contains

        subroutine write_ref_pft_range_local(funit, array, nrefs_out)
            integer,     intent(in) :: funit, nrefs_out
            complex(sp), intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%nrefs)
            write(unit=funit,pos=1) [self%pftsz, self%kfromto(1), self%interpklim, nrefs_out]
            write(unit=funit,pos=(4*sizeof(funit)+1)) array(:,:,iref_from:iref_to)
        end subroutine write_ref_pft_range_local
        
    end subroutine vol_pad2ref_pfts_write_range

    !>  \brief  Converts the polar references to a cartesian grid
    module subroutine polar_cavger_refs2cartesian( self, cavgs, which )
        class(polarft_calc),     intent(in)    :: self
        type(image),             intent(inout) :: cavgs(self%ncls)
        character(len=*),        intent(in)    :: which
        complex, allocatable :: cmat(:,:)
        real,    allocatable :: norm(:,:)
        integer, parameter :: EVEN_CASE = 0, ODD_CASE = 1, MERGED_CASE = 2  
        complex :: pft(1:self%pftsz,self%kfromto(1):self%interpklim), fc
        real    :: phys(2), dh,dk,mdk,mdh
        integer :: k,c,irot,physh,physk,box,icls,case_sel
        box = self%p_ptr%box_crop
        c   = box/2+1
        select case(trim(which))
            case('even')
                case_sel = EVEN_CASE
            case('odd')
                case_sel = ODD_CASE
            case('merged')
                case_sel = MERGED_CASE
        end select
        allocate(cmat(c,box),norm(c,box))
        !$omp parallel do schedule(guided) proc_bind(close) default(shared)&
        !$omp private(icls,pft,cmat,norm,irot,k,phys,fc,physh,physk,dh,dk,mdh,mdk)
        do icls = 1,self%ncls
            select case(case_sel)
                case(EVEN_CASE)
                    pft = cmplx(self%pfts_even(1:self%pftsz,self%kfromto(1):self%interpklim,icls), kind=sp)
                case(ODD_CASE)
                    pft = cmplx(self%pfts_odd(1:self%pftsz,self%kfromto(1):self%interpklim,icls), kind=sp)
                case(MERGED_CASE)
                    pft = cmplx(self%pfts_merg(1:self%pftsz,self%kfromto(1):self%interpklim,icls), kind=sp)
            end select
            ! Bi-linear interpolation
            cmat = CMPLX_ZERO
            norm = 0.0
            do irot = 1,self%pftsz
                do k = self%kfromto(1),self%interpklim
                    phys  = self%get_coord(irot,k) + [1.,real(c)]
                    fc    = pft(irot,k)
                    physh = floor(phys(1))
                    physk = floor(phys(2))
                    dh    = phys(1) - real(physh)
                    dk    = phys(2) - real(physk)
                    mdh   = 1.0 - dh
                    mdk   = 1.0 - dk
                    if( physh > 0 .and. physh <= c )then
                        if( physk <= box .and. physk >= 1 )then
                            cmat(physh,physk) = cmat(physh,physk) + mdh*mdk*fc
                            norm(physh,physk) = norm(physh,physk) + mdh*mdk
                            if( physk+1 <= box .and. physk+1 >= 1 )then
                                cmat(physh,physk+1) = cmat(physh,physk+1) + mdh*dk*fc
                                norm(physh,physk+1) = norm(physh,physk+1) + mdh*dk
                            endif
                        endif
                    endif
                    physh = physh + 1
                    if( physh > 0 .and. physh <= c )then
                        if( physk <= box .and. physk >= 1 )then
                            cmat(physh,physk) = cmat(physh,physk) + dh*mdk*fc
                            norm(physh,physk) = norm(physh,physk) + dh*mdk
                            if( physk+1 <= box .and. physk+1 >= 1 )then
                                cmat(physh,physk+1) = cmat(physh,physk+1) + dh*dk*fc
                                norm(physh,physk+1) = norm(physh,physk+1) + dh*dk
                            endif
                        endif
                    endif
                end do
            end do
            where( norm > TINY )
                cmat = cmat / norm
            elsewhere
                cmat = 0.0
            end where
            ! irot = self%pftsz+1, eg. angle=180.
            do k = 1,box/2-1
                cmat(1,k+c) = conjg(cmat(1,c-k))
            enddo
            ! arbitrary magnitude
            cmat(1,c) = CMPLX_ZERO
            ! set image
            call cavgs(icls)%new([box,box,1], self%p_ptr%smpd_crop, wthreads=.false.)
            call cavgs(icls)%set_cmat(cmat)
            call cavgs(icls)%shift_phorig()
            call cavgs(icls)%ifft
        enddo
        !$omp end parallel do
    end subroutine polar_cavger_refs2cartesian

    ! Read cavgs PFT array, poorly optimized but only used in unit test
    module subroutine polar_cavger_read( self, fname, which )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        character(len=*),    intent(in)    :: which
        select case(which)
            case('even')
                call self%read_pft_array(fname, self%pfts_even)
            case('odd')
                call self%read_pft_array(fname, self%pfts_odd)
            case('merged')
                call self%read_pft_array(fname, self%pfts_merg)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_read

    ! Writes cavgs PFT array, poorly optimized but only used in unit test
    module subroutine polar_cavger_write( self, fname, which )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: fname
        character(len=*),    intent(in) :: which
        select case(which)
            case('even')
                call self%write_pft_array(self%pfts_even, fname)
            case('odd')
                call self%write_pft_array(self%pfts_odd,  fname)
            case('merged')
                call self%write_pft_array(self%pfts_merg, fname)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_write

    ! Writes all cavgs PFT arrays
    !! performance critical code
    module subroutine polar_cavger_writeall( self, tmpl_fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: tmpl_fname
        integer :: funit_e, funit_o, funit_m
        integer :: i
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_even'//BIN_EXT, funit_e)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_odd'//BIN_EXT, funit_o)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//BIN_EXT, funit_m)
        !$omp parallel do default(shared) private(i) num_threads(3) schedule(static)
        do i = 1, 3
            select case(i)
                case(1)
                    call self%write_pft_array_local(funit_e, self%pfts_even)
                case(2)
                    call self%write_pft_array_local(funit_o, self%pfts_odd)
                case(3)
                    call self%write_pft_array_local(funit_m, self%pfts_merg)
            end select
        end do
        !$omp end parallel do
        call fclose(funit_e)
        call fclose(funit_o)
        call fclose(funit_m)
    end subroutine polar_cavger_writeall

    ! Write references contained in pftc
    !! performance critical code
    module subroutine polar_cavger_write_eo_pftcrefs( self, tmpl_fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: tmpl_fname
        integer :: funit_e, funit_o, funit_m, i
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_even'//BIN_EXT, funit_e)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_odd'//BIN_EXT, funit_o)
        !$omp parallel do default(shared) private(i) num_threads(2) schedule(static)
        do i = 1, 2
            select case(i)
                case(1)
                    self%pfts_even = cmplx(self%pfts_refs_even,kind=dp)
                    call self%write_pft_array_local(funit_e, self%pfts_even)
                case(2)
                    self%pfts_odd  = cmplx(self%pfts_refs_odd,kind=dp)
                    call self%write_pft_array_local(funit_o, self%pfts_odd)
            end select
        end do
        !$omp end parallel do
        call fclose(funit_e)
        call fclose(funit_o)
        call self%polar_cavger_zero_pft_refs !! removed after writing
    end subroutine polar_cavger_write_eo_pftcrefs

    ! Converts cavgs PFTS to cartesian grids and writes them
    !! this should be removed in performance critical code beacause it is costly
    !! and the outputted averages suffer from severe interpolation artefacts
    !! if we need Cartesian class averages we should generate them properly periodically
    module subroutine polar_cavger_write_cartrefs( self, tmpl_fname, which )
        class(polarft_calc), intent(in) :: self
        class(string),       intent(in) :: tmpl_fname
        character(len=*),    intent(in) :: which
        type(image), allocatable :: imgs(:)
        call alloc_imgarr(self%ncls, [self%p_ptr%box_crop, self%p_ptr%box_crop,1], self%p_ptr%smpd_crop, imgs)
        select case(trim(which))
            case('even','odd')
                call self%polar_cavger_refs2cartesian(imgs, trim(which) )
                call write_imgarr(imgs, tmpl_fname//'_'//trim(which)//MRC_EXT)
            case('merged')
                call self%polar_cavger_refs2cartesian(imgs, 'merged' )
                call write_imgarr(imgs, tmpl_fname//MRC_EXT)
        end select
        call dealloc_imgarr(imgs)
    end subroutine polar_cavger_write_cartrefs

    ! Reads all cavgs PFT arrays
    !! performance critical code
    module subroutine polar_cavger_read_all( self, fname )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        complex(sp), allocatable :: buf_e(:,:,:), buf_o(:,:,:), buf_m(:,:,:)
        type(string) :: refs, refs_even, refs_odd, ext
        integer :: funit_e, funit_o, funit_m
        integer :: dims_e(4), dims_o(4), dims_m(4)
        integer :: i
        ext = string('.')//fname2ext(fname)
        if( ext == MRC_EXT )then
            refs = get_fbody(fname, MRC_EXT, separator=.false.)//BIN_EXT
        elseif( ext == BIN_EXT )then
            refs = fname
        else
            THROW_HARD('Unsupported file format: '//ext%to_char())
        endif
        refs_even = get_fbody(refs,BIN_EXT,separator=.false.)//'_even'//BIN_EXT
        refs_odd  = get_fbody(refs,BIN_EXT,separator=.false.)//'_odd'//BIN_EXT
        if( .not. file_exists(refs) )then
            THROW_HARD('Polar references do not exist in cwd: '//refs%to_char())
        endif
        if( file_exists(refs_even) )then ! assume all files are there
            call self%open_pft_array_for_read(refs,      self%pfts_merg, funit_m, dims_m, buf_m)
            call self%open_pft_array_for_read(refs_even, self%pfts_even, funit_e, dims_e, buf_e)
            call self%open_pft_array_for_read(refs_odd,  self%pfts_odd,  funit_o, dims_o, buf_o)
            !$omp parallel do default(shared) private(i) num_threads(3) schedule(static)
            do i = 1, 3
                select case(i)
                    case(1)
                        call self%transfer_pft_array_buffer(self%pfts_merg, funit_m, dims_m, buf_m)
                    case(2)
                        call self%transfer_pft_array_buffer(self%pfts_even, funit_e, dims_e, buf_e)
                    case(3)
                        call self%transfer_pft_array_buffer(self%pfts_odd,  funit_o, dims_o, buf_o)
                end select
            end do
            !$omp end parallel do
            call fclose(funit_m)
            call fclose(funit_e)
            call fclose(funit_o)
            deallocate(buf_m, buf_e, buf_o)
        else
            call self%read_pft_array(refs, self%pfts_merg)
            !$omp parallel workshare
            self%pfts_even = self%pfts_merg
            self%pfts_odd  = self%pfts_merg
            !$omp end parallel workshare
        endif
    end subroutine polar_cavger_read_all

    !>  \brief  writes partial class averages (PFTS + CTF2) to disk (distributed execution)
    !! performance critical code
    module subroutine polar_cavger_readwrite_partial_sums( self, which )
        class(polarft_calc), intent(inout) :: self
        character(len=*),    intent(in)    :: which
        complex(sp), allocatable :: pfte_buf(:,:,:),  pfto_buf(:,:,:)
        real(sp),    allocatable :: ctf2e_buf(:,:,:), ctf2o_buf(:,:,:)
        type(string) :: cae, cao, cte, cto
        integer :: dims_cae(4), dims_cao(4), dims_cte(4), dims_cto(4)
        integer :: funit_cae, funit_cao, funit_cte, funit_cto
        integer :: i
        cae = 'cavgs_even_part'//int2str_pad(self%p_ptr%part,self%p_ptr%numlen)//BIN_EXT
        cao = 'cavgs_odd_part'//int2str_pad(self%p_ptr%part,self%p_ptr%numlen)//BIN_EXT
        cte = 'ctfsqsums_even_part'//int2str_pad(self%p_ptr%part,self%p_ptr%numlen)//BIN_EXT
        cto = 'ctfsqsums_odd_part'//int2str_pad(self%p_ptr%part,self%p_ptr%numlen)//BIN_EXT
        select case(trim(which))
            case('read')
                call self%open_pft_array_for_read(cae,  self%pfts_even, funit_cae, dims_cae, pfte_buf)
                call self%open_pft_array_for_read(cao,  self%pfts_odd,  funit_cao, dims_cao, pfto_buf)
                call self%open_ctf2_array_for_read(cte, self%ctf2_even, funit_cte, dims_cte, ctf2e_buf)
                call self%open_ctf2_array_for_read(cto, self%ctf2_odd,  funit_cto, dims_cto, ctf2o_buf)
                !$omp parallel do default(shared) private(i) num_threads(4) schedule(static)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call self%transfer_pft_array_buffer(self%pfts_even,  funit_cae, dims_cae, pfte_buf)
                        case(2)
                            call self%transfer_pft_array_buffer(self%pfts_odd,   funit_cao, dims_cao, pfto_buf)
                        case(3)
                            call self%transfer_ctf2_array_buffer(self%ctf2_even, funit_cte, dims_cte, ctf2e_buf)
                        case(4)
                            call self%transfer_ctf2_array_buffer(self%ctf2_odd,  funit_cto, dims_cto, ctf2o_buf)
                    end select
                end do
                !$omp end parallel do
                call fclose(funit_cae)
                call fclose(funit_cao)
                call fclose(funit_cte)
                call fclose(funit_cto)
                deallocate(pfte_buf, pfto_buf, ctf2e_buf, ctf2o_buf)
            case('write')
                call open_pft_or_ctf2_array_for_write(cae, funit_cae)
                call open_pft_or_ctf2_array_for_write(cao, funit_cao)
                call open_pft_or_ctf2_array_for_write(cte, funit_cte)
                call open_pft_or_ctf2_array_for_write(cto, funit_cto)
                !$omp parallel do default(shared) private(i) num_threads(4) schedule(static)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call self%write_pft_array_local(funit_cae, self%pfts_even)
                        case(2)
                            call self%write_pft_array_local(funit_cao, self%pfts_odd)
                        case(3)
                            call self%write_ctf2_array_local(funit_cte, self%ctf2_even)
                        case(4)
                            call self%write_ctf2_array_local(funit_cto, self%ctf2_odd)
                    end select
                end do
                !$omp end parallel do
                call fclose(funit_cae)
                call fclose(funit_cao)
                call fclose(funit_cte)
                call fclose(funit_cto)
            case DEFAULT
                THROW_HARD('unknown which flag; only read & write supported; cavger_readwrite_partial_sums')
        end select
        call cae%kill
        call cao%kill
        call cte%kill
        call cto%kill
    end subroutine polar_cavger_readwrite_partial_sums

    !>  \brief  Reads in and reduces partial matrices prior to restoration
    !! performance critical code
    module subroutine polar_cavger_assemble_sums_from_parts( self )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), allocatable :: pfte(:,:,:),      pfto(:,:,:)
        complex(sp), allocatable :: pfte_buf(:,:,:),  pfto_buf(:,:,:)
        real(dp),    allocatable :: ctf2e(:,:,:),     ctf2o(:,:,:)
        real(sp),    allocatable :: ctf2e_buf(:,:,:), ctf2o_buf(:,:,:)
        type(string) :: cae, cao, cte, cto
        integer :: dims_cae(4), dims_cao(4), dims_cte(4), dims_cto(4)
        integer :: funit_cae, funit_cao, funit_cte, funit_cto
        integer :: ipart, i 
        allocate(pfte(self%pftsz,self%kfromto(1):self%interpklim,self%ncls),  pfto(self%pftsz,self%kfromto(1):self%interpklim,self%ncls),&
               &ctf2e(self%pftsz,self%kfromto(1):self%interpklim,self%ncls), ctf2o(self%pftsz,self%kfromto(1):self%interpklim,self%ncls))
        call self%polar_cavger_zero_pft_refs
        do ipart = 1,self%p_ptr%nparts
            cae = 'cavgs_even_part'    //int2str_pad(ipart,self%p_ptr%numlen)//BIN_EXT
            cao = 'cavgs_odd_part'     //int2str_pad(ipart,self%p_ptr%numlen)//BIN_EXT
            cte = 'ctfsqsums_even_part'//int2str_pad(ipart,self%p_ptr%numlen)//BIN_EXT
            cto = 'ctfsqsums_odd_part' //int2str_pad(ipart,self%p_ptr%numlen)//BIN_EXT
            call self%open_pft_array_for_read(cae, pfte,  funit_cae, dims_cae, pfte_buf)
            call self%open_pft_array_for_read(cao, pfto,  funit_cao, dims_cao, pfto_buf)
            call self%open_ctf2_array_for_read(cte, ctf2e, funit_cte, dims_cte, ctf2e_buf)
            call self%open_ctf2_array_for_read(cto, ctf2o, funit_cto, dims_cto, ctf2o_buf)
            !$omp parallel do default(shared) private(i) num_threads(4) schedule(static)
            do i = 1, 4
                select case(i)
                    case(1)
                        call self%transfer_pft_array_buffer(pfte,  funit_cae, dims_cae, pfte_buf)
                        self%pfts_even = self%pfts_even + pfte
                    case(2)
                        call self%transfer_pft_array_buffer(pfto,  funit_cao, dims_cao, pfto_buf)
                        self%pfts_odd  = self%pfts_odd  + pfto
                    case(3)
                        call self%transfer_ctf2_array_buffer(ctf2e, funit_cte, dims_cte, ctf2e_buf)
                        self%ctf2_even = self%ctf2_even + ctf2e
                    case(4)
                        call self%transfer_ctf2_array_buffer(ctf2o, funit_cto, dims_cto, ctf2o_buf)
                        self%ctf2_odd  = self%ctf2_odd  + ctf2o
                end select
            end do
            !$omp end parallel do
            call fclose(funit_cae)
            call fclose(funit_cao)
            call fclose(funit_cte)
            call fclose(funit_cto)
        enddo
        deallocate(pfte_buf, pfto_buf, ctf2e_buf, ctf2o_buf) ! needs to be explicit, the others are in local scope
    end subroutine polar_cavger_assemble_sums_from_parts

    ! PFT IO HELPERS THAT ENABLE PARALLEL IO
    ! Format for PFT I/O
    ! First  integer: PFTSZ
    ! Second integer: KFROMTO(1)
    ! Third  integer: INTERPKLIM (upper bound of available on-disk k-range)
    ! Fourth integer: NCLS
    ! input/ouput in kind=dp but read/written in kind=sp

    ! private helper
    subroutine open_pft_or_ctf2_array_for_write( fname, funit )
        class(string), intent(in)  :: fname
        integer,       intent(out) :: funit
        integer :: io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("open_pft_or_ctf2_array_for_write: "//fname%to_char(),io_stat)
    end subroutine open_pft_or_ctf2_array_for_write

    ! private method
    module subroutine write_pft_array_local( self, funit, array )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: funit
        complex(dp),         intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        ! Header contract: dims = [pftsz, kfromto(1), interpklim, ncls].
        ! Therefore dims(2:3) describes the available on-disk k-range, not a search window.
        write(unit=funit,pos=1) [self%pftsz, self%kfromto(1), self%interpklim, self%ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) cmplx(array,kind=sp)
    end subroutine write_pft_array_local

    ! private method
    module subroutine write_pft_array( self, array, fname )
        class(polarft_calc), intent(in) :: self
        complex(dp),         intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        class(string),       intent(in) :: fname
        integer :: funit
        call open_pft_or_ctf2_array_for_write(fname, funit)
        call self%write_pft_array_local(funit, array)
        call fclose(funit)
    end subroutine write_pft_array

    ! private method
    module subroutine write_ctf2_array_local( self, funit, array )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: funit
        real(dp),            intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        ! Same header contract as PFT arrays: dims(2:3) = [kfromto(1), interpklim].
        write(unit=funit,pos=1) [self%pftsz, self%kfromto(1), self%interpklim, self%ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) real(array,kind=sp)
    end subroutine write_ctf2_array_local

    ! private method
    module subroutine write_ctf2_array( self, array, fname )
        class(polarft_calc), intent(in) :: self
        real(dp),            intent(in) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        class(string),       intent(in) :: fname
        integer :: funit
        call open_pft_or_ctf2_array_for_write(fname, funit)
        call self%write_ctf2_array_local(funit, array)
        call fclose(funit)
    end subroutine write_ctf2_array

    ! private method
    module subroutine get_pft_array_dims( self, fname, pftsz, kfromto, nrefs )
        class(polarft_calc), intent(in)  :: self
        class(string),       intent(in)  :: fname
        integer,             intent(out) :: pftsz, kfromto(2), nrefs
        integer :: io_stat, dims(4), funit
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        call fclose(funit)
        pftsz   = dims(1)
        ! Returned kfromto is the on-disk available range [kfromto(1), interpklim].
        kfromto = dims(2:3)
        nrefs   = dims(4)
    end subroutine get_pft_array_dims

    ! private method
    module subroutine open_pft_array_for_read( self, fname, array, funit, dims, buffer )
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,                  intent(out)   :: funit, dims(4)
        complex(sp), allocatable, intent(inout) :: buffer(:,:,:)
        integer :: io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pft_array; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        if( .not.allocated(array) )then
            allocate(array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls))
        endif
        if( dims(1) /= self%pftsz )then
            write(logfhandle,*) 'PFTZ header mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%pftsz
            write(logfhandle,*) 'found:    ', dims(1)
            THROW_HARD('Incompatible PFT header; open_pft_array_for_read')
        endif
        if( dims(4) /= self%ncls )then
            write(logfhandle,*) 'NCLS header mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%ncls
            write(logfhandle,*) 'found:    ', dims(4)
            THROW_HARD('Incompatible PFT header; open_pft_array_for_read')
        endif
        if( self%interpklim >= dims(3) )then
            ! identical size and padding allowed
        else if( self%interpklim < dims(3) )then
            write(logfhandle,*) 'INTERPKLIM header mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%interpklim
            write(logfhandle,*) 'found:    ', dims(3)
            THROW_HARD('Incompatible PFT header; open_pft_array_for_read')
        endif
        if( .not. allocated(buffer) ) allocate(buffer(dims(1),dims(2):dims(3),dims(4)))
    end subroutine open_pft_array_for_read
    
    ! private method
    module subroutine transfer_pft_array_buffer( self, array, funit, dims, buffer )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(inout) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,             intent(in)    :: funit, dims(4)
        complex(sp),         intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
        integer :: klo, khi
        ! Read stored (single-precision) array payload
        read(unit=funit, pos=(sizeof(dims)+1)) buffer
        ! Default to zero padding everywhere
        array(:,:,:) = (0.0_dp, 0.0_dp)
        ! Copy only the overlap in k between requested kfromto and stored dims(2:3)
        klo = max(self%kfromto(1), dims(2))
        khi = min(self%interpklim, dims(3))
        if( klo <= khi ) array(:,klo:khi,:) = cmplx(buffer(:,klo:khi,:), kind=dp)
    end subroutine transfer_pft_array_buffer

    ! private helper
    module subroutine read_pft_array( self, fname, array)
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        complex(sp), allocatable :: buffer(:,:,:)
        integer :: dims(4), funit
        call self%open_pft_array_for_read(fname, array, funit, dims, buffer)
        call self%transfer_pft_array_buffer(     array, funit, dims, buffer )
        deallocate(buffer)
        call fclose(funit)
    end subroutine read_pft_array

    ! private helper, no checks are performed on dimensions
    module subroutine read_any_pft_array( self, fname, array )
        class(polarft_calc),      intent(in)    :: self
        class(string),            intent(in)    :: fname
        complex(sp), allocatable, intent(inout) :: array(:,:,:)
        integer :: dims(4), funit
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call open_pft_or_ctf2_array_for_write(fname, funit)
        read(unit=funit,pos=1) dims
        if( allocated(array) ) deallocate(array)
        allocate(array(dims(1),dims(2):dims(3),dims(4)))
        read(unit=funit, pos=(sizeof(dims)+1)) array
        call fclose(funit)
    end subroutine read_any_pft_array

    ! private helper
    module subroutine open_ctf2_array_for_read( self, fname, array, funit, dims, buffer )
        class(polarft_calc),   intent(in)    :: self
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,               intent(out)   :: funit, dims(4)
        real(sp), allocatable, intent(inout) :: buffer(:,:,:)
        integer :: io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('open_ctf2_array_for_read; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        if( .not.allocated(array) )then
            allocate(array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls))
        endif
        if( dims(1) /= self%pftsz )then
            write(logfhandle,*) 'PFTZ header mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%pftsz
            write(logfhandle,*) 'found:    ', dims(1)
            THROW_HARD('Incompatible CTF2 header; open_pft_array_for_read')
        endif
        if( dims(4) /= self%ncls )then
            write(logfhandle,*) 'NCLS header mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%ncls
            write(logfhandle,*) 'found:    ', dims(4)
            THROW_HARD('Incompatible CTF2 header; open_pft_array_for_read')
        endif
        if( self%interpklim >= dims(3) )then
            ! identical size and padding allowed
        else if( self%interpklim < dims(3) )then
            write(logfhandle,*) 'INTERPKLIM header mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%interpklim
            write(logfhandle,*) 'found:    ', dims(3)
            THROW_HARD('Incompatible CTF2 header; open_pft_array_for_read')
        endif
        if( .not. allocated(buffer) ) allocate(buffer(dims(1),dims(2):dims(3),dims(4)))
    end subroutine open_ctf2_array_for_read

    ! private helper
    module subroutine transfer_ctf2_array_buffer( self, array, funit, dims, buffer )
        class(polarft_calc), intent(in) :: self
        real(dp), intent(inout) :: array(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        integer,  intent(in)    :: funit, dims(4)
        real(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
        integer :: klo, khi
        ! Read stored (single-precision) array payload
        read(unit=funit, pos=(sizeof(dims)+1)) buffer
        ! Default to zero padding everywhere
        array(:,:,:) = 0.0_dp
        ! Copy only the overlap in k between requested kfromto and stored dims(2:3)
        klo = max(self%kfromto(1), dims(2))
        khi = min(self%interpklim, dims(3))
        if( klo <= khi )then
            array(:,klo:khi,:) = real(buffer(:,klo:khi,:), dp)
        endif
    end subroutine transfer_ctf2_array_buffer

    ! Reads an even or odd projection-range binary written by vol_pad2ref_pfts_write_range
    ! and places the payload at the given global reference positions in pfts_refs_even/odd.
    module subroutine read_ref_pfts_range( self, fname, iseven, iref_from, iref_to )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        logical,             intent(in)    :: iseven
        integer,             intent(in)    :: iref_from, iref_to
        complex(sp), allocatable :: buf(:,:,:)
        integer :: funit, io_stat, dims(4), klo, khi
        if( .not. file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist; read_ref_pfts_range')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_ref_pfts_range; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit, pos=1) dims
        if( dims(1) /= self%pftsz )then
            write(logfhandle,*) 'pftsz mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', self%pftsz, ' found: ', dims(1)
            THROW_HARD('Incompatible header; read_ref_pfts_range')
        endif
        if( dims(4) /= iref_to - iref_from + 1 )then
            write(logfhandle,*) 'nrefs mismatch in: ', trim(fname%to_char())
            write(logfhandle,*) 'expected: ', iref_to-iref_from+1, ' found: ', dims(4)
            THROW_HARD('Incompatible header; read_ref_pfts_range')
        endif
        allocate(buf(dims(1), dims(2):dims(3), dims(4)))
        read(unit=funit, pos=(sizeof(dims)+1)) buf
        call fclose(funit)
        klo = max(self%kfromto(1), dims(2))
        khi = min(self%interpklim, dims(3))
        if( klo <= khi )then
            if( iseven )then
                self%pfts_refs_even(:, klo:khi, iref_from:iref_to) = buf(:, klo:khi, :)
            else
                self%pfts_refs_odd( :, klo:khi, iref_from:iref_to) = buf(:, klo:khi, :)
            endif
        endif
        deallocate(buf)
    end subroutine read_ref_pfts_range

    ! Reassembles pfts_refs_even/odd from per-state/per-part files produced by workers
    ! that called vol_pad2ref_pfts_write_range.  The partition over nspace is re-derived
    ! with split_nobjs_even so that it exactly mirrors the worker split.
    module subroutine assemble_projected_refs_from_parts( self, nparts, numlen )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: nparts, numlen
        integer, allocatable :: parts(:,:)
        type(string) :: fname
        integer :: ipart, s, iref_from, iref_to
        if( .not. self%existence ) THROW_HARD('polarft_calc does not exist; assemble_projected_refs_from_parts')
        parts = split_nobjs_even(self%p_ptr%nspace, nparts)
        do s = 1, self%p_ptr%nstates
            do ipart = 1, nparts
                iref_from = (s - 1) * self%p_ptr%nspace + parts(ipart, 1)
                iref_to   = (s - 1) * self%p_ptr%nspace + parts(ipart, 2)
                fname = string(POLAR_REFS_FBODY)//'_s'//int2str_pad(s,2)  &
                      & //'_part'//int2str_pad(ipart, numlen)//'_even'//BIN_EXT
                call self%read_ref_pfts_range(fname, .true.,  iref_from, iref_to)
                fname = string(POLAR_REFS_FBODY)//'_s'//int2str_pad(s,2)  &
                      & //'_part'//int2str_pad(ipart, numlen)//'_odd'//BIN_EXT
                call self%read_ref_pfts_range(fname, .false., iref_from, iref_to)
            end do
        end do
        call fname%kill
        deallocate(parts)
    end subroutine assemble_projected_refs_from_parts

end submodule simple_polarft_ops_io

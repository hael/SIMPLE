submodule (simple_polarops) simple_polarops_io
implicit none
#include "simple_local_flags.inc"
contains

    !>  \brief  Converts the polar references to a cartesian grid
    module subroutine polar_cavger_refs2cartesian( pftc, cavgs, which, pfts_in )
        use simple_image
        class(polarft_calc),     intent(in)    :: pftc
        type(image),             intent(inout) :: cavgs(ncls)
        character(len=*),        intent(in)    :: which
        complex(dp),   optional, intent(in)    :: pfts_in(1:pftsz,kfromto(1):kfromto(2),1:ncls)
        complex, allocatable :: cmat(:,:)
        real,    allocatable :: norm(:,:)
        integer, parameter :: EVEN_CASE = 0, ODD_CASE = 1, MERGED_CASE = 2  
        complex :: pft(1:pftsz,kfromto(1):kfromto(2)), fc
        real    :: phys(2), dh,dk,mdk,mdh
        integer :: k,c,irot,physh,physk,box,icls,case_sel
        logical :: pfts_in_present
        pfts_in_present = present(pfts_in)
        box = params_glob%box_crop
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
        do icls = 1, ncls
            if( pfts_in_present )then
                pft = cmplx(pfts_in(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
            else
                select case(case_sel)
                    case(EVEN_CASE)
                        pft = cmplx(pfts_even(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
                    case(ODD_CASE)
                        pft = cmplx(pfts_odd(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
                    case(MERGED_CASE)
                        pft = cmplx(pfts_merg(1:pftsz,kfromto(1):kfromto(2),icls), kind=sp)
                end select
            endif
            ! Bi-linear interpolation
            cmat = CMPLX_ZERO
            norm = 0.0
            do irot = 1,pftsz
                do k = kfromto(1),kfromto(2)
                    phys  = pftc%get_coord(irot,k) + [1.,real(c)]
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
            call cavgs(icls)%new([box,box,1], smpd, wthreads=.false.)
            call cavgs(icls)%set_cmat(cmat)
            call cavgs(icls)%shift_phorig()
            call cavgs(icls)%ifft
        enddo
        !$omp end parallel do
    end subroutine polar_cavger_refs2cartesian

    ! Read cavgs PFT array, poorly optimized but only used in unit test
    module subroutine polar_cavger_read( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        select case(which)
            case('even')
                call read_pft_array(fname, pfts_even)
            case('odd')
                call read_pft_array(fname, pfts_odd)
            case('merged')
                call read_pft_array(fname, pfts_merg)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_read

    ! Writes cavgs PFT array, poorly optimized but only used in unit test
    module subroutine polar_cavger_write( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        select case(which)
            case('even')
                call write_pft_array(pfts_even, fname)
            case('odd')
                call write_pft_array(pfts_odd,  fname)
            case('merged')
                call write_pft_array(pfts_merg, fname)
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine polar_cavger_write

    ! Writes all cavgs PFT arrays
    !! performance critical code
    module subroutine polar_cavger_writeall( tmpl_fname )
        class(string), intent(in) :: tmpl_fname
        integer :: funit_e, funit_o, funit_m
        integer :: i
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_even'//BIN_EXT, funit_e)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_odd'//BIN_EXT, funit_o)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//BIN_EXT, funit_m)
        !$omp parallel do default(shared) private(i) num_threads(3) schedule(static)
        do i = 1, 3
            select case(i)
                case(1)
                    call write_pft_array_local(funit_e, pfts_even)
                case(2)
                    call write_pft_array_local(funit_o, pfts_odd)
                case(3)
                    call write_pft_array_local(funit_m, pfts_merg)
            end select
        end do
        !$omp end parallel do
        call fclose(funit_e)
        call fclose(funit_o)
        call fclose(funit_m)
    end subroutine polar_cavger_writeall

    ! Write references contained in pftc
    !! performance critical code
    module subroutine polar_cavger_writeall_pftcrefs( tmpl_fname )
        class(string),  intent(in) :: tmpl_fname
        complex(sp),  pointer :: ptre(:,:,:), ptro(:,:,:)
        integer :: funit_e, funit_o, i
        call pftc_glob%get_refs_ptr(ptre, ptro)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_even'//BIN_EXT, funit_e)
        call open_pft_or_ctf2_array_for_write(tmpl_fname//'_odd'//BIN_EXT, funit_o)
        !$omp parallel do default(shared) private(i) num_threads(2) schedule(static)
        do i = 1, 2
            select case(i)
                case(1)
                    pfts_even = cmplx(ptre,kind=dp)
                    call write_pft_array_local(funit_e, pfts_even)
                case(2)
                    pfts_odd  = cmplx(ptro,kind=dp)
                    call write_pft_array_local(funit_o, pfts_odd)
            end select
        end do
        !$omp end parallel do
        call fclose(funit_e)
        call fclose(funit_o)
        call polar_cavger_zero_pft_refs !! removed after writing
        nullify(ptre, ptro)
    end subroutine polar_cavger_writeall_pftcrefs

    ! Converts cavgs PFTS to cartesian grids and writes them
    !! this should be removed in performance critical code beacause it is costly
    !! and the outputted averages suffer from severe interpolation artefacts
    !! if we need Cartesian class averages we should generate them properly periodically
    module subroutine polar_cavger_write_cartrefs( pftc, tmpl_fname, which )
        class(polarft_calc), intent(in) :: pftc
        class(string),       intent(in) :: tmpl_fname
        character(len=*),    intent(in) :: which
        type(image), allocatable :: imgs(:)
        call alloc_imgarr(ncls, [params_glob%box_crop, params_glob%box_crop,1], smpd, imgs)
        select case(trim(which))
            case('even','odd')
                call polar_cavger_refs2cartesian( pftc, imgs, trim(which) )
                call write_imgarr(imgs, tmpl_fname//'_'//trim(which)//params_glob%ext%to_char())
            case('merged')
                call polar_cavger_refs2cartesian( pftc, imgs, 'merged' )
                call write_imgarr(imgs, tmpl_fname//params_glob%ext)
        end select
        call dealloc_imgarr(imgs)
    end subroutine polar_cavger_write_cartrefs

    ! Reads all cavgs PFT arrays
    !! performance critical code
    module subroutine polar_cavger_read_all( fname )
        class(string),  intent(in) :: fname
        complex(sp), allocatable :: buf_e(:,:,:), buf_o(:,:,:), buf_m(:,:,:)
        type(string) :: refs, refs_even, refs_odd, ext
        integer :: funit_e, funit_o, funit_m
        integer :: dims_e(4), dims_o(4), dims_m(4)
        integer :: i
        ext = string('.')//fname2ext(fname)
        if( ext == params_glob%ext )then
            refs = get_fbody(fname, params_glob%ext, separator=.false.)//BIN_EXT
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
            call open_pft_array_for_read(refs,      pfts_merg, funit_m, dims_m, buf_m)
            call open_pft_array_for_read(refs_even, pfts_even, funit_e, dims_e, buf_e)
            call open_pft_array_for_read(refs_odd,  pfts_odd,  funit_o, dims_o, buf_o)
            !$omp parallel do default(shared) private(i) num_threads(3) schedule(static)
            do i = 1, 3
                select case(i)
                    case(1)
                        call transfer_pft_array_buffer(pfts_merg, funit_m, dims_m, buf_m)
                    case(2)
                        call transfer_pft_array_buffer(pfts_even, funit_e, dims_e, buf_e)
                    case(3)
                        call transfer_pft_array_buffer(pfts_odd,  funit_o, dims_o, buf_o)
                end select
            end do
            !$omp end parallel do
            call fclose(funit_m)
            call fclose(funit_e)
            call fclose(funit_o)
            deallocate(buf_m, buf_e, buf_o)
        else
            call read_pft_array(refs, pfts_merg)
            !$omp parallel workshare
            pfts_even = pfts_merg
            pfts_odd  = pfts_merg
            !$omp end parallel workshare
        endif
    end subroutine polar_cavger_read_all

    !>  \brief  writes partial class averages (PFTS + CTF2) to disk (distributed execution)
    !! performance critical code
    subroutine polar_cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        complex(dp), allocatable :: pfte(:,:,:),      pfto(:,:,:)
        complex(sp), allocatable :: pfte_buf(:,:,:),  pfto_buf(:,:,:)
        real(dp),    allocatable :: ctf2e(:,:,:),     ctf2o(:,:,:)
        real(sp),    allocatable :: ctf2e_buf(:,:,:), ctf2o_buf(:,:,:)
        type(string) :: cae, cao, cte, cto
        integer :: dims_cae(4), dims_cao(4), dims_cte(4), dims_cto(4)
        integer :: funit_cae, funit_cao, funit_cte, funit_cto
        integer :: i
        cae = 'cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        cao = 'cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        cte = 'ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        cto = 'ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//BIN_EXT
        select case(trim(which))
            case('read')
                call open_pft_array_for_read(cae,  pfts_even, funit_cae, dims_cae, pfte_buf)
                call open_pft_array_for_read(cao,  pfts_odd,  funit_cao, dims_cao, pfto_buf)
                call open_ctf2_array_for_read(cte, ctf2_even, funit_cte, dims_cte, ctf2e_buf)
                call open_ctf2_array_for_read(cto, ctf2_odd,  funit_cto, dims_cto, ctf2o_buf)
                !$omp parallel do default(shared) private(i) num_threads(4) schedule(static)
                do i = 1, 4
                    select case(i)
                        case(1)
                            call transfer_pft_array_buffer(pfts_even,  funit_cae, dims_cae, pfte_buf)
                        case(2)
                            call transfer_pft_array_buffer(pfts_odd,   funit_cao, dims_cao, pfto_buf)
                        case(3)
                            call transfer_ctf2_array_buffer(ctf2_even, funit_cte, dims_cte, ctf2e_buf)
                        case(4)
                            call transfer_ctf2_array_buffer(ctf2_odd,  funit_cto, dims_cto, ctf2o_buf)
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
                            call write_pft_array_local(funit_cae, pfts_even)
                        case(2)
                            call write_pft_array_local(funit_cao, pfts_odd)
                        case(3)
                            call write_ctf2_array_local(funit_cte, ctf2_even)
                        case(4)
                            call write_ctf2_array_local(funit_cto, ctf2_odd)
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
    subroutine polar_cavger_assemble_sums_from_parts( reforis, clin_anneal )
        type(oris), optional, intent(in) :: reforis
        real,       optional, intent(in) :: clin_anneal
        complex(dp), allocatable :: pfte(:,:,:),      pfto(:,:,:)
        complex(sp), allocatable :: pfte_buf(:,:,:),  pfto_buf(:,:,:)
        real(dp),    allocatable :: ctf2e(:,:,:),     ctf2o(:,:,:)
        real(sp),    allocatable :: ctf2e_buf(:,:,:), ctf2o_buf(:,:,:)
        type(string) :: cae, cao, cte, cto
        integer :: dims_cae(4), dims_cao(4), dims_cte(4), dims_cto(4)
        integer :: funit_cae, funit_cao, funit_cte, funit_cto
        integer :: ipart, i 
        allocate(pfte(pftsz,kfromto(1):kfromto(2),ncls),  pfto(pftsz,kfromto(1):kfromto(2),ncls),&
               &ctf2e(pftsz,kfromto(1):kfromto(2),ncls), ctf2o(pftsz,kfromto(1):kfromto(2),ncls))
        call polar_cavger_zero_pft_refs
        do ipart = 1,params_glob%nparts
            cae = 'cavgs_even_part'    //int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            cao = 'cavgs_odd_part'     //int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            cte = 'ctfsqsums_even_part'//int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            cto = 'ctfsqsums_odd_part' //int2str_pad(ipart,params_glob%numlen)//BIN_EXT
            call open_pft_array_for_read(cae, pfte,  funit_cae, dims_cae, pfte_buf)
            call open_pft_array_for_read(cao, pfto,  funit_cao, dims_cao, pfto_buf)
            call open_ctf2_array_for_read(cte, ctf2e, funit_cte, dims_cte, ctf2e_buf)
            call open_ctf2_array_for_read(cto, ctf2o, funit_cto, dims_cto, ctf2o_buf)
            !$omp parallel do default(shared) private(i) num_threads(4) schedule(static)
            do i = 1, 4
                select case(i)
                    case(1)
                        call transfer_pft_array_buffer(pfte,  funit_cae, dims_cae, pfte_buf)
                        pfts_even = pfts_even + pfte
                    case(2)
                        call transfer_pft_array_buffer(pfto,  funit_cao, dims_cao, pfto_buf)
                        pfts_odd  = pfts_odd  + pfto
                    case(3)
                        call transfer_ctf2_array_buffer(ctf2e, funit_cte, dims_cte, ctf2e_buf)
                        ctf2_even = ctf2_even + ctf2e
                    case(4)
                        call transfer_ctf2_array_buffer(ctf2o, funit_cto, dims_cto, ctf2o_buf)
                        ctf2_odd  = ctf2_odd  + ctf2o
                end select
            end do
            !$omp end parallel do
            call fclose(funit_cae)
            call fclose(funit_cao)
            call fclose(funit_cte)
            call fclose(funit_cto)
        enddo
        deallocate(pfte_buf, pfto_buf, ctf2e_buf, ctf2o_buf) ! needs to be explicit, the others are in local scope
        ! merge eo-pairs and normalize
        select case(trim(params_glob%ref_type))
            case('polar_cavg')
                call polar_cavger_merge_eos_and_norm2D
            case DEFAULT
                call polar_cavger_merge_eos_and_norm(reforis=reforis, cl_weight=clin_anneal)
        end select
    end subroutine polar_cavger_assemble_sums_from_parts

    subroutine polar_cavger_dims_from_header( fname, pftsz_here, kfromto_here, ncls_here )
        class(string), intent(in)    :: fname
        integer,       intent(inout) :: pftsz_here, kfromto_here(2), ncls_here
        integer :: dims(4), funit, io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('dims_from_header; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        call fclose(funit)
        pftsz_here   = dims(1)
        kfromto_here = dims(2:3)
        ncls_here    = dims(4)
    end subroutine polar_cavger_dims_from_header

    ! PFT IO HELPERS THAT ENABLE PARALLEL IO
    ! Format for PFT I/O
    ! First  integer: PFTSZ
    ! Second integer: KFROMTO(1)
    ! Third  integer: KFROMTO(2)
    ! Fourth integer: NCLS
    ! input/ouput in kind=dp but read/written in kind=sp

    module subroutine open_pft_or_ctf2_array_for_write( fname, funit )
        class(string), intent(in)  :: fname
        integer,       intent(out) :: funit
        integer :: io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("open_pft_or_ctf2_array_for_write: "//fname%to_char(),io_stat)
    end subroutine open_pft_or_ctf2_array_for_write

    module subroutine write_pft_array_local( funit, array )
        integer,     intent(in) :: funit
        complex(dp), intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) cmplx(array,kind=sp)
    end subroutine write_pft_array_local

    module subroutine write_pft_array( array, fname )
        complex(dp),   intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        class(string), intent(in) :: fname
        integer :: funit
        call open_pft_or_ctf2_array_for_write(fname, funit)
        call write_pft_array_local(funit, array)
        call fclose(funit)
    end subroutine write_pft_array

    module subroutine write_ctf2_array_local( funit, array )
        integer,  intent(in) :: funit
        real(dp), intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        write(unit=funit,pos=1) [pftsz, kfromto(1), kfromto(2), ncls]
        write(unit=funit,pos=(4*sizeof(funit)+1)) real(array,kind=sp)
    end subroutine write_ctf2_array_local

    module subroutine write_ctf2_array( array, fname )
        real(dp),      intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        class(string), intent(in) :: fname
        integer :: funit
        call open_pft_or_ctf2_array_for_write(fname, funit)
        call write_ctf2_array_local(funit, array)
        call fclose(funit)
    end subroutine write_ctf2_array

    module subroutine open_pft_array_for_read( fname, array, funit, dims, buffer )
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
            allocate(array(pftsz,kfromto(1):kfromto(2),ncls))
        endif
        if( .not. all(dims == [pftsz, kfromto(1), kfromto(2), ncls]) )then
            if( pftsz /= dims(1) )then
                THROW_HARD('Incompatible PFT size in '//fname%to_char()//': '//int2str(pftsz)//' vs '//int2str(dims(1)))
            endif
            if( ncls /= dims(4) )then
                THROW_HARD('Incompatible NCLS in '//fname%to_char()//': '//int2str(ncls)//' vs '//int2str(dims(4)))
            endif
        endif
        if( .not. allocated(buffer) ) allocate(buffer(dims(1),dims(2):dims(3),dims(4)))
    end subroutine open_pft_array_for_read

    module subroutine transfer_pft_array_buffer( array, funit, dims, buffer )
        complex(dp), intent(inout) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        integer,     intent(in)    :: funit, dims(4)
        complex(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
        integer :: k, klo, khi
        ! Read stored (single-precision) array payload
        read(unit=funit, pos=(sizeof(dims)+1)) buffer
        ! Default to zero padding everywhere
        array(:,:,:) = (0.0_dp, 0.0_dp)
        ! Copy only the overlap in k between requested kfromto and stored dims(2:3)
        klo = max(kfromto(1), dims(2))
        khi = min(kfromto(2), dims(3))
        if( klo <= khi )then
            do k = klo, khi
                array(:,k,:) = cmplx(buffer(:,k,:), kind=dp)
            enddo
        endif
    end subroutine transfer_pft_array_buffer

    module subroutine read_pft_array( fname, array)
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        complex(sp), allocatable :: buffer(:,:,:)
        integer :: dims(4), funit
        call open_pft_array_for_read(fname, array, funit, dims, buffer)
        call transfer_pft_array_buffer(     array, funit, dims, buffer )
        deallocate(buffer)
        call fclose(funit)
    end subroutine read_pft_array

    module subroutine open_ctf2_array_for_read( fname, array, funit, dims, buffer )
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,               intent(out)   :: funit, dims(4)
        real(sp), allocatable, intent(inout) :: buffer(:,:,:)
        integer :: io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_ctf2_array; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        if( .not.allocated(array) )then
            allocate(array(pftsz,kfromto(1):kfromto(2),ncls))
        endif
        if( .not. all(dims == [pftsz, kfromto(1), kfromto(2), ncls]) )then
            if( pftsz /= dims(1) )then
                THROW_HARD('Incompatible real array size in '//fname%to_char()//': '//int2str(pftsz)//' vs '//int2str(dims(1)))
            endif
            if( ncls /= dims(4) )then
                THROW_HARD('Incompatible NCLS in '//fname%to_char()//': '//int2str(ncls)//' vs '//int2str(dims(4)))
            endif
        endif
        if( .not. allocated(buffer) ) allocate(buffer(dims(1),dims(2):dims(3),dims(4)))
    end subroutine open_ctf2_array_for_read

    module subroutine transfer_ctf2_array_buffer( array, funit, dims, buffer )
        real(dp), intent(inout) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        integer,  intent(in)    :: funit, dims(4)
        real(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
        integer :: k, klo, khi
        ! Read stored (single-precision) array payload
        read(unit=funit, pos=(sizeof(dims)+1)) buffer
        ! Default to zero padding everywhere
        array(:,:,:) = 0.0_dp
        ! Copy only the overlap in k between requested kfromto and stored dims(2:3)
        klo = max(kfromto(1), dims(2))
        khi = min(kfromto(2), dims(3))
        if( klo <= khi )then
            do k = klo, khi
                array(:,k,:) = real(buffer(:,k,:), dp)
            enddo
        endif
    end subroutine transfer_ctf2_array_buffer

    module subroutine read_ctf2_array( fname, array )
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        real(sp), allocatable :: buffer(:,:,:)
        integer :: dims(4), funit
        call open_ctf2_array_for_read(fname, array, funit, dims, buffer)
        call transfer_ctf2_array_buffer(     array, funit, dims, buffer)
        deallocate(buffer)
        call fclose(funit)
    end subroutine read_ctf2_array

    ! writes out |PFTs_MERG| as mrcs
    module subroutine pft2img()
        type(image) :: img
        integer :: nk,i,k,icls
        nk = kfromto(2)
        if( .not.is_even(nk) ) nk = nk+1
        call img%new([pftsz,nk,1],1.0)
        do icls = 1,ncls
            img = 0.0
            do i = 1,pftsz
                do k = kfromto(1),kfromto(2)
                    call img%set([i,k,1], real(abs(pfts_merg(i,k,icls))))
                enddo
            enddo
            call img%write(string('pfts_it'//int2str(params_glob%which_iter)//'.mrc'),icls)
        enddo
        call img%kill
    end subroutine pft2img

end submodule simple_polarops_io

module simple_comlin
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_oris,             only: oris
use simple_polarft_corrcalc, only: polarft_corrcalc
implicit none

public :: comlin_map, polar_comlin_pfts, gen_polar_comlins, read_write_comlin

contains

    subroutine comlin_map( lims, e_ind, eall, coord_map, all_coords )
        integer,         intent(in)    :: lims(3,2)
        integer,         intent(in)    :: e_ind
        class(oris),     intent(in)    :: eall
        type(fplan_map), intent(inout) :: coord_map
        type(fplan_map), intent(inout) :: all_coords
        type(ori) :: e, e2
        integer   :: i, h, k, sqlp, sqarg, xy2(2), cnt, ori_phys(3), phys(3)
        logical   :: good_coord
        call eall%get_ori(e_ind, e)
        sqlp = (maxval(lims(:,2)))**2
        cnt  = 0
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                sqarg = dot_product([h,k],[h,k])
                if( sqarg > sqlp ) cycle
                if (h .ge. 0) then
                    ori_phys(1) = h + 1
                    ori_phys(2) = k + 1 + MERGE(params_glob%box,0,k < 0)
                    ori_phys(3) = 1
                else
                    ori_phys(1) = -h + 1
                    ori_phys(2) = -k + 1 + MERGE(params_glob%box,0,-k < 0)
                    ori_phys(3) = 1
                endif
                do i = 1, eall%get_noris()
                    if( i == e_ind ) cycle
                    call eall%get_ori(i, e2)
                    call comlin_coord(lims, [h,k], e, e2, xy2, good_coord)
                    if( good_coord )then
                        cnt = cnt + 1
                        all_coords%tar_find(  cnt) = i
                        all_coords%ori_phys(:,cnt) = ori_phys
                        all_coords%ori_four(:,cnt) = [h,k]
                        if (xy2(1) .ge. 0) then
                            phys(1) = xy2(1) + 1
                            phys(2) = xy2(2) + 1 + MERGE(params_glob%box,0,xy2(2) < 0)
                            phys(3) = 1
                        else
                            phys(1) = -xy2(1) + 1
                            phys(2) = -xy2(2) + 1 + MERGE(params_glob%box,0,-xy2(2) < 0)
                            phys(3) = 1
                        endif
                        all_coords%tar_phys(:,cnt) = phys
                        all_coords%tar_four(:,cnt) = xy2
                    endif
                enddo
            enddo
        enddo
        allocate(coord_map%tar_find(  cnt), source=all_coords%tar_find(  1:cnt))
        allocate(coord_map%ori_phys(3,cnt), source=all_coords%ori_phys(:,1:cnt))
        allocate(coord_map%ori_four(2,cnt), source=all_coords%ori_four(:,1:cnt))
        allocate(coord_map%tar_phys(3,cnt), source=all_coords%tar_phys(:,1:cnt))
        allocate(coord_map%tar_four(2,cnt), source=all_coords%tar_four(:,1:cnt))
        coord_map%n_points = cnt
    end subroutine comlin_map

    ! projecting coordinates xy1 of the plane at orientation e1 to the coordinates xy1 of the plane at e2
    subroutine comlin_coord(lims, xy1, e1, e2, xy2, good_coord)
        integer,      intent(in)    :: lims(3,2)
        integer,      intent(in)    :: xy1(2)
        type(ori),    intent(in)    :: e1
        type(ori),    intent(in)    :: e2
        integer,      intent(inout) :: xy2(2)
        logical,      intent(out)   :: good_coord
        real    :: e1_rotmat(3,3), e2_rotmat(3,3), e2_inv(3,3), loc1_3D(3), loc2_3D(3), e1_inv(3,3)
        integer :: xy(3), sqlp
        good_coord = .false.
        if( all(xy1 == 0) )return
        sqlp      = (maxval(lims(:,2)))**2
        e1_rotmat = e1%get_mat()
        e2_rotmat = e2%get_mat()
        ! rotating coordinates xy1 of the first plane
        loc1_3D   = matmul(real([xy1(1),xy1(2),0]), e1_rotmat)
        ! inversely rotating rotated xy1 to the second plane
        e2_inv    = transpose(e2_rotmat)
        xy        = nint(matmul(loc1_3D, e2_inv))
        xy2       = xy(1:2)
        if( all(xy2 == 0) ) return
        ! checking the inversely rotated coords make sense
        if( xy(1) >= lims(1,1) .and. xy(1) <= lims(1,2) .and. &
           &xy(2) >= lims(2,1) .and. xy(2) <= lims(2,2) .and. &
           &xy(3) == 0         .and. dot_product(xy2,xy2) <= sqlp )then
            ! double check after this
        else
            return
        endif
        ! reverse checking
        loc2_3D = matmul(real([xy2(1),xy2(2),0]), e2_rotmat)
        e1_inv  = transpose(e1_rotmat)
        xy      = nint(matmul(loc2_3D, e1_inv))
        ! checking the inversely rotated coords make sense
        if( xy(1) >= lims(1,1) .and. xy(1) <= lims(1,2) .and. &
           &xy(2) >= lims(2,1) .and. xy(2) <= lims(2,2) .and. &
           &xy(3) == 0         .and. dot_product(xy(1:2),xy(1:2)) <= sqlp .and. all(xy(1:2) == xy1) )then
            good_coord = .true.
        endif
    end subroutine comlin_coord

    subroutine gen_polar_comlins( pftcc, ref_space, pcomlines)
        use simple_oris
        class(polarft_corrcalc), intent(in)    :: pftcc
        type(oris),              intent(in)    :: ref_space
        type(polar_fmap),        intent(inout) :: pcomlines(:,:)
        integer :: iref, jref, irot, kind, irot_l, irot_r, nrefs, pftsz
        real    :: loc1_3D(3), loc2_3D(3), denom, a1, a2, b1, b2, line2D(3), irot_real, k_real, w, hk1(2), hk2(2),&
                  &rotmat(3,3),invmats(3,3,pftcc%get_nrefs()), loc1s(3,pftcc%get_nrefs()), loc2s(3,pftcc%get_nrefs()), line3D(3)
        nrefs = pftcc%get_nrefs()
        pftsz = pftcc%get_pftsz()
        ! randomly chosing two sets of (irot, kind) to generate the polar common lines
        irot  = 5
        kind  = params_glob%kfromto(1) + 5
        if( kind >= params_glob%kfromto(2) ) kind = params_glob%kfromto(1) ! using kfromto(1) as kind in the first set of (irot,kind)
        hk1   = pftcc%get_coord(irot,kind)
        irot  = 16
        kind  = params_glob%kfromto(1) + 15
        if( kind >= params_glob%kfromto(2) ) kind = params_glob%kfromto(2) ! using kfromto(2) as kind in the second set of (irot,kind)
        hk2   = pftcc%get_coord(irot,kind)
        ! caching rotation matrices and their corresponding inverse matrices
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,rotmat)
        do iref = 1, nrefs
            rotmat            = ref_space%get_mat(iref)
            invmats(:,:,iref) = transpose(rotmat)
            loc1s(:,iref)     = matmul([hk1(1), hk1(2), 0.], rotmat)
            loc2s(:,iref)     = matmul([hk2(1), hk2(2), 0.], rotmat)
        enddo
        !$omp end parallel do
        ! constructing polar common lines
        pcomlines%legit = .false.
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(iref,loc1_3D,loc2_3D,denom,a1,b1,jref,a2,b2,line3D,line2D,irot_real,k_real,irot_l,irot_r,w)
        do iref = 1, nrefs
            loc1_3D = loc1s(:,iref)
            loc2_3D = loc2s(:,iref)
            denom   = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
            a1      = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
            b1      = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
            do jref = 1, nrefs
                if( jref == iref )cycle
                ! getting the 3D common line
                loc1_3D     = loc1s(:,jref)
                loc2_3D     = loc2s(:,jref)
                denom       = (loc1_3D(1) * loc2_3D(2) - loc1_3D(2) * loc2_3D(1))
                if( abs(denom) < TINY )cycle
                a2          = (loc1_3D(3) * loc2_3D(2) - loc1_3D(2) * loc2_3D(3)) / denom
                b2          = (loc1_3D(1) * loc2_3D(3) - loc1_3D(3) * loc2_3D(1)) / denom
                if( abs(b1-b2) < TINY )cycle
                line3D(1:2) = [1., -(a1-a2)/(b1-b2)]
                line3D(3)   = a2*line3D(1) + b2*line3D(2)
                ! projecting the 3D common line to a polar line on the jref-th reference
                line2D      = matmul(line3D, invmats(:,:,jref))
                call pftcc%get_polar_coord(line2D(1:2), irot_real, k_real)
                if( irot_real < 1. ) irot_real = irot_real + real(pftsz)
                ! caching the indices irot_j and irot_j+1 and the corresponding linear weight
                irot_l = floor(irot_real)
                irot_r = irot_l + 1
                w      = irot_real - real(irot_l)
                if( irot_l > pftsz ) irot_l = irot_l - pftsz
                if( irot_r > pftsz ) irot_r = irot_r - pftsz
                pcomlines(jref,iref)%targ_irot_l = irot_l
                pcomlines(jref,iref)%targ_irot_r = irot_r
                pcomlines(jref,iref)%targ_w      = w
                ! projecting the 3D common line to a polar line on the iref-th reference
                line2D = matmul(line3D, invmats(:,:,iref))
                call pftcc%get_polar_coord(line2D(1:2), irot_real, k_real)
                if( irot_real < 1. ) irot_real = irot_real + real(pftsz)
                ! caching the indices irot_i and irot_i+1 and the corresponding linear weight
                irot_l = floor(irot_real)
                irot_r = irot_l + 1
                w      = irot_real - real(irot_l)
                if( irot_l > pftsz ) irot_l = irot_l - pftsz
                if( irot_r > pftsz ) irot_r = irot_r - pftsz
                pcomlines(jref,iref)%self_irot_l = irot_l
                pcomlines(jref,iref)%self_irot_r = irot_r
                pcomlines(jref,iref)%self_w      = w
                pcomlines(jref,iref)%legit       = .true.
            enddo
        enddo
        !$omp end parallel do
    end subroutine gen_polar_comlins

    subroutine polar_comlin_pfts( pcomlines, pfts_in, pfts)
        type(polar_fmap), intent(in)    :: pcomlines(:,:)
        complex,          intent(in)    :: pfts_in(:,:,:)
        complex,          intent(inout) :: pfts(:,:,:)
        complex :: pft_line(params_glob%kfromto(1):params_glob%kfromto(2))
        integer :: iref, jref, irot_l, irot_r, nrefs
        real    :: w
        nrefs = size(pcomlines,1)
        pfts  = complex(0.,0.)
        !$omp parallel do default(shared) private(iref,jref,irot_l,irot_r,w,pft_line)&
        !$omp proc_bind(close) schedule(static)
        do iref = 1, nrefs
            do jref = 1, nrefs
                if( .not. pcomlines(jref,iref)%legit )cycle
                ! compute the interpolated polar common line, between irot_j and irot_j+1
                irot_l   = pcomlines(jref,iref)%targ_irot_l
                irot_r   = pcomlines(jref,iref)%targ_irot_r
                w        = pcomlines(jref,iref)%targ_w
                pft_line = (1.-w) * pfts_in(irot_l,:,jref) + w * pfts_in(irot_r,:,jref)
                ! extrapolate the interpolated polar common line to irot_i and irot_i+1 of iref-th reference
                irot_l   = pcomlines(jref,iref)%self_irot_l
                irot_r   = pcomlines(jref,iref)%self_irot_r
                w        = pcomlines(jref,iref)%self_w
                pfts(irot_l,:,iref) = pfts(irot_l,:,iref) + (1.-w) * pft_line
                pfts(irot_r,:,iref) = pfts(irot_r,:,iref) +     w  * pft_line
            enddo
        enddo
        !$omp end parallel do
    end subroutine polar_comlin_pfts

    subroutine read_write_comlin( pcomlines, pftcc, eulspace )
        type(polar_fmap), allocatable, intent(inout) :: pcomlines(:,:)
        class(polarft_corrcalc),       intent(in)    :: pftcc
        type(oris),                    intent(in)    :: eulspace
        integer :: io_stat, funit
        if( .not. allocated(pcomlines) ) allocate(pcomlines(params_glob%nspace,params_glob%nspace))
        if( file_exists(trim(POLAR_COMLIN)) )then
            call fopen(funit,trim(POLAR_COMLIN),access='STREAM',action='READ',status='OLD', iostat=io_stat)
            read(unit=funit,pos=1) pcomlines
        else
            call gen_polar_comlins(pftcc, eulspace, pcomlines)
            call fopen(funit,trim(POLAR_COMLIN),access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
            write(unit=funit,pos=1) pcomlines
        endif
        call fclose(funit)
    end subroutine read_write_comlin

end module simple_comlin

module simple_comlin
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_oris,       only: oris
implicit none

public :: comlin_map

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
                        all_coords%target_find(cnt) = i
                        all_coords%ori_phys(:,cnt)  = ori_phys
                        all_coords%ori_four(:,cnt)  = [h,k]
                        if (xy2(1) .ge. 0) then
                            phys(1) = xy2(1) + 1
                            phys(2) = xy2(2) + 1 + MERGE(params_glob%box,0,xy2(2) < 0)
                            phys(3) = 1
                        else
                            phys(1) = -xy2(1) + 1
                            phys(2) = -xy2(2) + 1 + MERGE(params_glob%box,0,-xy2(2) < 0)
                            phys(3) = 1
                        endif
                        all_coords%target_phys(:,cnt) = phys
                        all_coords%target_four(:,cnt) = xy2
                    endif
                enddo
            enddo
        enddo
        allocate(coord_map%target_find(  cnt), source=all_coords%target_find(  1:cnt))
        allocate(coord_map%ori_phys(   3,cnt), source=all_coords%ori_phys(:,   1:cnt))
        allocate(coord_map%ori_four(   3,cnt), source=all_coords%ori_four(:,   1:cnt))
        allocate(coord_map%target_phys(3,cnt), source=all_coords%target_phys(:,1:cnt))
        allocate(coord_map%target_four(3,cnt), source=all_coords%target_four(:,1:cnt))
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
        integer :: errflg, xy(3), sqlp
        good_coord = .false.
        if( all(xy1 == 0) )return
        sqlp      = (maxval(lims(:,2)))**2
        e1_rotmat = e1%get_mat()
        e2_rotmat = e2%get_mat()
        ! rotating coordinates xy1 of the first plane
        loc1_3D   = matmul(real([xy1(1),xy1(2),0]), e1_rotmat)
        ! inversely rotating rotated xy1 to the second plane
        call matinv(e2_rotmat, e2_inv, 3, errflg)
        if( errflg < 0 ) return
        xy  = nint(matmul(loc1_3D, e2_inv))
        xy2 = xy(1:2)
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
        call matinv(e1_rotmat, e1_inv, 3, errflg)
        if( errflg < 0 ) return
        xy = nint(matmul(loc2_3D, e1_inv))
        ! checking the inversely rotated coords make sense
        if( xy(1) >= lims(1,1) .and. xy(1) <= lims(1,2) .and. &
           &xy(2) >= lims(2,1) .and. xy(2) <= lims(2,2) .and. &
           &xy(3) == 0         .and. dot_product(xy(1:2),xy(1:2)) <= sqlp .and. all(xy(1:2) == xy1) )then
            good_coord = .true.
        endif
    end subroutine comlin_coord

end module simple_comlin
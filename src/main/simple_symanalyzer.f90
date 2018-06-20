module simple_symanalyzer


contains

    subroutine symmetrise_map( vol_in, pgrp, hp, lp  )
        use simple_volpft_symsrch
        use simple_sym,       only: sym
        use simple_ori,       only: ori, m2euler
        use simple_projector, only: projector
        class(projector), intent(inout) :: vol_in
        character(len=*), intent(in)    :: pgrp
        real,             intent(in)    :: hp, lp
        type(ori)         :: symaxis, o
        type(sym)         :: symobj
        integer           :: isym, nsym
        real              :: rmat_symaxis(3,3), rmat(3,3)
        real, allocatable :: sym_rmats(:,:,:)


        ! make point-group object
        call symobj%new(pgrp)
        nsym = symobj%get_nsym()


        ! extract the rotation matrices for the symops
        allocate(sym_rmats(isym,3,3))
        do isym=1,nsym
            o = symobj%get_symori(isym)
            sym_rmats(isym,:,:) = o%get_mat()
        end do



        ! init search object
        call volpft_symsrch_init(vol_in, pgrp, hp, lp)
        ! search
        call volpft_srch4symaxis(symaxis)
        ! get the rotation matrix of the symaxis
        rmat_symaxis = symaxis%get_mat()

        do isym=1,nsym
            rmat = matmul(sym_rmats(isym,:,:), rmat_symaxis)
            call o%set_euler(m2euler(rmat))



        end do
        ! destruct
        call o%kill
        call symaxis%kill
        call symobj%kill
    end subroutine symmetrise_map

end module simple_symanalyzer

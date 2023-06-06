! Subtracts stack 2 from stack 1
program simple_test_signal_subtr

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none
    
    type(image),    allocatable :: cavgs(:), reprojs(:), out(:)
    character(*), parameter    :: stks='cavgs_vs_reprojections_rec_and_sim.mrc'

    ! Read in stacks

    ! Subtract reprojs from cavgs

    ! Write output

end program simple_test_signal_subtr


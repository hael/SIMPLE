module simple_polarftcc_online



complex(sp), allocatable         :: sums_trailing(:,:), sums_leading(:,:)
complex(sp), allocatable         :: refs(:,:), refs_ctf(:,:)


! express all rotations as moving windows
! average in polar Fourier space
! update the model with a learning rate:
!        sums_trailing * lambda + sums_leading * (1 - lambda)
! randomly assign parameters for all imgs
! shellweights in polar coordinates
! particles with no assignment search the current fraction of search space






end module simple_polarftcc_online
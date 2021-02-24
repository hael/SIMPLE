module simple_np_coordination_number
  use simple_math
  use simple_fileio

  implicit none

  public :: run_coord_number_analysis
  private

contains

  ! This function calculates the coordination number for each atom
  ! in the input model and prints it on a txt file.
  ! ATTENTION: input coords of model have to be in ANGSTROMS.
  subroutine run_coord_number_analysis (model, d, contact_scores )
    real, allocatable,    intent(in)    :: model(:,:)
    real,                 intent(in)    :: d   ! radius for coordination number analyis
    integer, optional,    intent(inout) :: contact_scores(size(model,2))
    integer :: cn(size(model,2))
    integer :: filnum, io_stat, natoms, iatom
    natoms = size(model, 2)
    call fopen(filnum, file='CN.txt', iostat=io_stat)
    call calc_cn(cn)
    do iatom  = 1, natoms
      write(filnum,*) cn(iatom)
    enddo
    call fclose(filnum)
    if(present(contact_scores)) contact_scores = cn

  contains
    subroutine calc_cn(coordination_number)
      integer,intent(inout) :: coordination_number(natoms)
      real    :: dist
      integer :: jatom,cnt ! counter
      coordination_number = 0
      do iatom = 1, natoms
        cnt = 0 ! reset counter
        do jatom = 1, natoms
          if(iatom /= jatom) then
            dist = euclid(model(:,iatom), model(:,jatom))
            if(dist < d) cnt = cnt + 1
          endif
        enddo
        coordination_number(iatom) = cnt
      enddo
    end subroutine calc_cn
  end subroutine run_coord_number_analysis


end module simple_np_coordination_number

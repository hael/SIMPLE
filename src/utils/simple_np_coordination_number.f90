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
  subroutine run_coord_number_analysis (model, d, coord_numbers, coord_numbers_gen)
    real, allocatable,    intent(in)    :: model(:,:)
    real,                 intent(in)    :: d   ! radius for coordination number analyis
    integer,              intent(inout) :: coord_numbers(size(model,2))
    real,                 intent(inout) :: coord_numbers_gen(size(model,2))
    integer :: filnum, io_stat, natoms, iatom
    natoms = size(model, 2)
    call fopen(filnum, file='CN.txt', iostat=io_stat)
    call calc_cn(coord_numbers,coord_numbers_gen)
    do iatom  = 1, natoms
      write(filnum,*) coord_numbers(iatom)
    enddo
    call fclose(filnum)

  contains
    subroutine calc_cn(coordination_numbers,coordination_numbers_gen)
      integer,intent(inout) :: coordination_numbers(natoms)
      real,   intent(inout) :: coordination_numbers_gen(natoms)
      real    :: dist
      integer :: jatom,cnt ! counter
      integer :: cn_max(natoms)
      coordination_numbers     = 0
      coordination_numbers_gen = 0
      cn_max(:)                = 0
      !cn
      do iatom = 1, natoms
        cnt = 0 ! reset counter, nb of neighbours
        do jatom = 1, natoms
          if(iatom /= jatom) then
            dist = euclid(model(:,iatom), model(:,jatom))
            if(dist < d) cnt = cnt + 1
          endif
        enddo
        coordination_numbers(iatom) = cnt
      enddo
      !cn_gen
      do iatom = 1,natoms
        cnt    = 0  ! sum of coordination numbers
        do jatom = 1, natoms
          if(iatom /= jatom) then
            dist = euclid(model(:,iatom), model(:,jatom))
            if(dist < d) then ! if they are neighbours
              cnt = cnt + coordination_numbers(jatom)
              if(coordination_numbers(jatom) > cn_max(iatom)) then
                cn_max(iatom) = coordination_numbers(jatom)
              endif
            endif
          endif
        enddo
        if(cn_max(iatom)>0) then
          coord_numbers_gen(iatom) = real(cnt)/real(cn_max(iatom))
        else
          coord_numbers_gen(iatom) = 0.
        endif
      enddo
    end subroutine calc_cn
  end subroutine run_coord_number_analysis


end module simple_np_coordination_number

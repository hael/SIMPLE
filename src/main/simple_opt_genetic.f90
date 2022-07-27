! simple probabilistic genetic algorithm
module simple_opt_genetic
include 'simple_lib.f08'
implicit none

contains
    subroutine crossover(parent1, parent2, cross_rate, child1, child2)
        integer, intent(in)    :: parent1(:), parent2(:)
        real   , intent(in)    :: cross_rate
        integer, intent(inout) :: child1(:), child2(:)
        integer :: cross_pt, N
        child1 = parent1
        child2 = parent2
        N      = size(parent1)
        if( rand() < cross_rate )then
            cross_pt = 1 + floor(rand()*(N-1))
            child1(         1:cross_pt) = parent1(         1:cross_pt)
            child1(cross_pt+1:N)        = parent2(cross_pt+1:N)
            child2(         1:cross_pt) = parent2(         1:cross_pt)
            child2(cross_pt+1:N)        = parent1(cross_pt+1:N)
        endif
    end subroutine crossover

    subroutine mutation(bitstring, mut_rate)
        integer, intent(inout) :: bitstring(:)
        real   , intent(in)    :: mut_rate
        integer :: k
        do k = 1, size(bitstring)
            if( rand() < mut_rate ) bitstring(k) = 1 - bitstring(k)
        enddo
    end subroutine mutation

    function select_pop(pop, costs, times) result(ind)
        integer, intent(in)  :: pop(:,:)
        real,    intent(in)  :: costs(:)
        integer, intent(in)  :: times
        integer :: ind
        integer :: k, tmp_ind
        ind = 1+floor(rand()*size(pop,1))
        do k = 1, times
            tmp_ind = 1+floor(rand()*size(pop,1))
            if( costs(tmp_ind) > costs(ind) ) ind = tmp_ind
        enddo
    end function select_pop

    subroutine genetic_opt(num_bits, max_iter, pop_size, cross_rate, mut_rate, best, best_val)
        integer, intent(in)    :: num_bits
        integer, intent(in)    :: max_iter
        integer, intent(in)    :: pop_size
        real   , intent(in)    :: cross_rate
        real   , intent(in)    :: mut_rate
        integer, intent(inout) :: best(:)
        real   , intent(inout) :: best_val
        integer, allocatable :: pop(:,:), children(:,:), selected(:,:), p1(:), p2(:), c1(:), c2(:)
        real   , allocatable :: costs(:)
        integer :: k,l
        call srand(time())
        allocate(pop(pop_size, num_bits), children(pop_size, num_bits), selected(pop_size, num_bits),&
                &p1(num_bits), p2(num_bits), c1(num_bits), c2(num_bits), source = 0)
        allocate(costs(pop_size), source = 0.)
        do k = 1, pop_size
            do l = 1, num_bits
                pop(k,l) = floor(rand()*2)
            enddo
        enddo
        best     = pop(1,:)
        best_val = objective(best)
        do k = 1, max_iter
            do l = 1, pop_size
                costs(l) = objective(pop(l,:))
                if( costs(l) > best_val )then
                    best     = pop(l,:)
                    best_val = costs(l)
                    print *, 'new best = ', best
                    print *, 'current cost = ', best_val
                endif
            enddo
            do l = 1, pop_size
                selected(l,:) = pop(select_pop(pop, costs, times = 3), :)
            enddo
            do l = 1, pop_size, 2
                p1 = selected(l  , :)
                p2 = selected(l+1, :)
                call crossover(p1, p2, cross_rate, c1, c2)
                call mutation(c1, mut_rate)
                call mutation(c2, mut_rate)
                children(l  ,:) = c1
                children(l+1,:) = c2
            enddo
            pop = children
        enddo
    end subroutine genetic_opt

    function objective(bitstring) result(val)
        integer, intent(in)  :: bitstring(:)
        real :: val
        integer, parameter :: GOAL(20) = [0,1,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,1,1,0]
        integer :: k
        val = 0
        do k = 1, size(GOAL)
            if( bitstring(k) == GOAL(k) ) val = val + 1
        enddo
    end function objective
end module simple_opt_genetic

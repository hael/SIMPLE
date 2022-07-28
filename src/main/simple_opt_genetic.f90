! simple probabilistic genetic algorithm
module simple_opt_genetic
    include 'simple_lib.f08'
    implicit none

    abstract interface
        function objective_func(bitstring) result(val)
            integer, intent(in)  :: bitstring(:)
            real                 :: val
        end function objective_func
    end interface
contains
    subroutine crossover(bitstr1, bitstr2, cross_rate)
        integer, intent(inout) :: bitstr1(:), bitstr2(:)
        real   , intent(in)    :: cross_rate
        integer :: cross_pt, N
        N = size(bitstr1)
        if( rand() < cross_rate )then
            cross_pt = 1 + floor(rand()*(N-1))
            bitstr2(         1:cross_pt) = bitstr1(         1:cross_pt)
            bitstr1(cross_pt+1:N)        = bitstr2(cross_pt+1:N)
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

    subroutine genetic_opt(obj_func, num_bits, max_iter, pop_size, cross_rate, mut_rate, best, best_val)
        procedure(objective_func), pointer :: obj_func
        integer, intent(in)    :: num_bits
        integer, intent(in)    :: max_iter
        integer, intent(in)    :: pop_size
        real   , intent(in)    :: cross_rate
        real   , intent(in)    :: mut_rate
        integer, intent(inout) :: best(:)
        real   , intent(inout) :: best_val
        integer, allocatable :: pop(:,:), selected(:,:)
        real   , allocatable :: costs(:)
        integer :: k,l
        logical :: got_new_best
        call srand(time())
        allocate(pop(pop_size, num_bits), selected(pop_size, num_bits), source = 0)
        allocate(costs(pop_size), source = 0.)
        do k = 1, pop_size
            do l = 1, num_bits
                pop(k,l) = floor(rand()*2)
            enddo
        enddo
        best     = pop(1,:)
        best_val = obj_func(best)
        do k = 1, max_iter
            got_new_best = .false.
            do l = 1, pop_size
                costs(l) = obj_func(pop(l,:))
                if( costs(l) > best_val )then
                    best         = pop(l,:)
                    best_val     = costs(l)
                    got_new_best = .true.
                endif
            enddo
            if( got_new_best )then
                print *, 'new best = ', best
                print *, 'current cost = ', best_val
            endif
            do l = 1, pop_size, 2
                selected(l  ,:) = pop(select_pop(pop, costs, times = 3), :)
                selected(l+1,:) = pop(select_pop(pop, costs, times = 3), :)
                call crossover(selected(l,:), selected(l+1,:), cross_rate)
                call mutation(selected(l  ,:), mut_rate)
                call mutation(selected(l+1,:), mut_rate)
            enddo
            pop = selected
        enddo
    end subroutine genetic_opt
end module simple_opt_genetic

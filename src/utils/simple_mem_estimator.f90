! job ram usage estimation
module simple_mem_estimator
include 'simple_lib.f08'
use simple_parameters, only: parameters
implicit none

contains

    subroutine estimate_mem_usage(job_descr, q_descr, extra_params)
        class(chash),               intent(inout) :: q_descr
        class(chash),               intent(in)    :: job_descr
        type(parameters), optional, intent(in)    :: extra_params
        type(string) :: prg
        integer      :: top, fromp, np, part, io_stat      
        if(job_descr%isthere("prg")) then
            prg   = job_descr%get("prg")
            top   = 0
            fromp = 0
            part  = 1
            if(job_descr%isthere("top"))   top   = str2int(job_descr%get("top"),   io_stat)
            if(job_descr%isthere("fromp")) fromp = str2int(job_descr%get("fromp"), io_stat)
            if(job_descr%isthere("part"))  part  = str2int(job_descr%get("part"),  io_stat)
            np = top - fromp + 1
            select case(prg%to_char())
                case("calc_pspec")
                    write(logfhandle,'(A)') "MEMPSPEC " // int2str(extra_params%box)
                    if(present(extra_params) .and. extra_params%box > 0) then
                        call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_pspec(part, extra_params%box, np, extra_params%nthr)))
                    endif
                case("cluster2D")
                    if(present(extra_params) .and. extra_params%box > 0) then
                        call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_2D(part, extra_params%box, np, extra_params%ncls, extra_params%nthr)))
                    endif
                case("ctf_estimate")
                    if(present(extra_params) .and. extra_params%box > 0) then
                     !   call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_2D(part, extra_params%box, np, extra_params%ncls, extra_params%nthr)))
                    endif
                case("motion_correct")
                    write(logfhandle,'(A)') 'XDIM ' // trim(int2str(extra_params%xdim)) // ' YDIM ' // trim(int2str(extra_params%ydim)) 
                    if(present(extra_params) .and. extra_params%box > 0) then
                      !  call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_motion_corr(part, extra_params%box, extra_params%nthr)))
                    endif
                case("pick")
                    if(present(extra_params) .and. extra_params%box > 0) then
                       ! call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_pick(part, extra_params%box, extra_params%nthr)))
                    endif    
                case("preprocess")
                    if(present(extra_params) .and. extra_params%box > 0) then
                      !  call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_preprocess(part, extra_params%box, extra_params%nthr)))
                    endif   
                case("scale")
                    if(present(extra_params) .and. extra_params%box > 0) then
                        call q_descr%set('job_memory_per_task', int2str(estimate_mem_usage_scale(part, extra_params%box, extra_params%nthr)))
                    endif
            end select
            call prg%kill
        endif
    end subroutine estimate_mem_usage
    
    function estimate_mem_usage_2D(part, boxsize, nptcls, ncls, nthr) result (memusagerounded)
        integer, intent(in) :: part, boxsize, nptcls, ncls, nthr
        integer             :: memusagerounded
        real                :: memusage, baseline, clsgradient, threadcoefficient, ptcl1000coefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        ptcl1000coefficient = 0.02 ! % increase in mem per 1000 particles
        baseline    = 450000 + 23 * boxsize * boxsize
        clsgradient = (0.1 * boxsize * boxsize) * (1 + (nptcls/1000.0) * ptcl1000coefficient)
        memusage = (baseline + clsgradient * ncls) * (1 + threadcoefficient * nthr)
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (CLS2D) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_2D
    
    function estimate_mem_usage_pspec(part, boxsize, nptcls, nthr) result (memusagerounded)
        integer, intent(in) :: part, boxsize, nptcls, nthr
        integer             :: memusagerounded
        real                :: memusage, baseline, threadcoefficient, ptclcoefficient
        threadcoefficient   = 0.125 ! % increase in mem per thread
        ptclcoefficient = 0.02 ! % increase in mem per 1000 particles
        baseline = (50000 + 1.2 * boxsize * boxsize) + ( nptcls * boxsize * ptclcoefficient)
        memusage = baseline * (1 + threadcoefficient * nthr)
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (PSPEC) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_pspec
    
    function estimate_mem_usage_scale(part, boxsize, nthr) result (memusagerounded)
        integer, intent(in) :: part, boxsize, nthr
        integer             :: memusagerounded
        real                :: memusage, baseline, threadcoefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        baseline    = 250000 + 1100 * boxsize
        memusage = baseline * (1 + threadcoefficient * nthr)
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (SCALE) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_scale
   
    function estimate_mem_usage_motion_corr(part, micsizex, micsizey, nthr) result (memusagerounded)
        integer, intent(in) :: part, micsizex, micsizey, nthr
        integer             :: memusagerounded
        real                :: memusage, threadcoefficient, ptcl1000coefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        ptcl1000coefficient = 0.02 ! % increase in mem per 1000 particles
       ! baseline    = 450000 + 23 * boxsize * boxsize
        !clsgradient = (0.1 * boxsize * boxsize) * (1 + (nptcls/1000.0) * ptcl1000coefficient)
        !memusage = (baseline + clsgradient * ncls) * (1 + threadcoefficient * nthr)
        memusage = 100000000
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (MOTION CORRECTION) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_motion_corr
    
    function estimate_mem_usage_preprocess(part, micsizex, micsizey, nthr) result (memusagerounded)
        integer, intent(in) :: part, micsizex, micsizey, nthr
        integer             :: memusagerounded
        real                :: memusage, threadcoefficient, ptcl1000coefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        ptcl1000coefficient = 0.02 ! % increase in mem per 1000 particles
       ! baseline    = 450000 + 23 * boxsize * boxsize
        !clsgradient = (0.1 * boxsize * boxsize) * (1 + (nptcls/1000.0) * ptcl1000coefficient)
        !memusage = (baseline + clsgradient * ncls) * (1 + threadcoefficient * nthr)
        memusage = 100000000
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (PREPROCESS) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_preprocess
    
    function estimate_mem_usage_ctf_estimation(part, micsizex, micsizey, nthr) result (memusagerounded)
        integer, intent(in) :: part, micsizex, micsizey, nthr
        integer             :: memusagerounded
        real                :: memusage, threadcoefficient, ptcl1000coefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        ptcl1000coefficient = 0.02 ! % increase in mem per 1000 particles
       ! baseline    = 450000 + 23 * boxsize * boxsize
        !clsgradient = (0.1 * boxsize * boxsize) * (1 + (nptcls/1000.0) * ptcl1000coefficient)
        !memusage = (baseline + clsgradient * ncls) * (1 + threadcoefficient * nthr)
        memusage = 100000000
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (CTF ESTIMATION) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_ctf_estimation  
    
    function estimate_mem_usage_pick(part, micsizex, micsizey, nthr) result (memusagerounded)
        integer, intent(in) :: part, micsizex, micsizey, nthr
        integer             :: memusagerounded
        real                :: memusage, threadcoefficient, ptcl1000coefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        ptcl1000coefficient = 0.02 ! % increase in mem per 1000 particles
       ! baseline    = 450000 + 23 * boxsize * boxsize
        !clsgradient = (0.1 * boxsize * boxsize) * (1 + (nptcls/1000.0) * ptcl1000coefficient)
        !memusage = (baseline + clsgradient * ncls) * (1 + threadcoefficient * nthr)
        memusage = 100000000
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (PICK) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_pick  
    
    function estimate_mem_usage_extract(part, micsizex, micsizey, nthr) result (memusagerounded)
        integer, intent(in) :: part, micsizex, micsizey, nthr
        integer             :: memusagerounded
        real                :: memusage, threadcoefficient, ptcl1000coefficient
        threadcoefficient   = 0.01 ! % increase in mem per thread
        ptcl1000coefficient = 0.02 ! % increase in mem per 1000 particles
       ! baseline    = 450000 + 23 * boxsize * boxsize
        !clsgradient = (0.1 * boxsize * boxsize) * (1 + (nptcls/1000.0) * ptcl1000coefficient)
        !memusage = (baseline + clsgradient * ncls) * (1 + threadcoefficient * nthr)
        memusage = 100000000
        memusagerounded = ceiling(memusage / 1000000) * 1000
        if(part == 1) write(logfhandle,'(A)') '>>> ESTIMATED MEMORY REQUIREMENT (EXTRACT) : ' // trim(int2str(memusagerounded)) // ' (' // trim(int2str(ceiling(memusage / 1000))) // ') MB'
    end function estimate_mem_usage_extract

    
end module simple_mem_estimator

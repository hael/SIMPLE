module simple_environment
include 'simple_lib.f08'

implicit none

contains

    subroutine print_environment()
            
            CHARACTER(len = 255) :: environment_value
            CHARACTER(len = 63), Allocatable :: environment_variables(:);
            INTEGER :: i, stat, len
                
            len = 63

            Allocate (environment_variables(7))
            
            environment_variables = [CHARACTER(len = 63) :: "SLURM_JOBID", "SLURM_JOB_USER", "SLURM_JOB_CPUS_PER_NODE", "SLURM_MEM_PER_CPU", "SLURMD_NODENAME", "SLURM_JOB_ACCOUNT", "SLURM_SUBMIT_DIR"]
            
            CALL get_environment_variable("SLURM_JOBID", environment_value, len, stat)
            
            if(stat .eq. 0) then
                    
                    write(logfhandle,*) ""
                    write(logfhandle,*) "##### SIMPLE SLURM Environment #####"
                    write(logfhandle,*) ""

                    do i=1,7
                
                        CALL get_environment_variable(TRIM(environment_variables(i)), environment_value, len, stat)
            
                        if(stat .eq. 0) then
            
                                write(logfhandle,*) TRIM(environment_variables(i)), achar(9), " : ", achar(9), TRIM(environment_value)
            
                        endif
            
                    end do
           

                    write(logfhandle,*) ""
                    write(logfhandle,*) "####################################"
                    write(logfhandle,*) ""

            endif

            Deallocate (environment_variables)

    end subroutine print_environment

end module simple_environment

module continuous_3D_pcg_refinement_noise_test_support
use iso_fortran_env, only: real64
implicit none
private

integer, parameter, public :: NOISE_DP = real64
public :: centered_correlation
public :: population_mean
public :: population_variance
public :: relative_error

contains

pure real(NOISE_DP) function population_mean(values) result(mean_value)
    real, intent(in) :: values(:,:,:)

    mean_value = sum(real(values,NOISE_DP)) / real(size(values),NOISE_DP)
end function population_mean

pure real(NOISE_DP) function population_variance(values) result(variance)
    real, intent(in) :: values(:,:,:)
    real(NOISE_DP) :: mean_value

    mean_value = population_mean(values)
    variance = sum((real(values,NOISE_DP) - mean_value)**2) / real(size(values),NOISE_DP)
end function population_variance

pure real(NOISE_DP) function centered_correlation(first, second) result(correlation)
    real, intent(in) :: first(:,:,:), second(:,:,:)
    real(NOISE_DP) :: denominator, first_mean, second_mean

    first_mean = population_mean(first)
    second_mean = population_mean(second)
    denominator = sqrt(sum((real(first,NOISE_DP) - first_mean)**2) * &
        &sum((real(second,NOISE_DP) - second_mean)**2))
    if( denominator <= tiny(denominator) )then
        correlation = huge(correlation)
    else
        correlation = sum((real(first,NOISE_DP) - first_mean) * &
            &(real(second,NOISE_DP) - second_mean)) / denominator
    endif
end function centered_correlation

pure real(NOISE_DP) function relative_error(actual, expected) result(error)
    real(NOISE_DP), intent(in) :: actual, expected

    error = abs(actual - expected) / max(abs(expected), tiny(expected))
end function relative_error

end module continuous_3D_pcg_refinement_noise_test_support

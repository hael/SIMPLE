#include <cstddef>
#include <cstdint>
#include <cstring>

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__) && !defined(__clang__)
#define SIMPLE_HAS_F16C_DISPATCH 1
#include <immintrin.h>
#else
#define SIMPLE_HAS_F16C_DISPATCH 0
#endif

namespace
{
using converter_type = void (*)(const std::int16_t *, float *, std::size_t);
using encoder_type = void (*)(const float *, std::int16_t *, std::size_t);
static_assert(sizeof(std::int16_t) == 2, "binary16 storage must be 16 bits");
static_assert(sizeof(float) == sizeof(std::uint32_t), "float32 storage must be 32 bits");

inline float bits_to_float(std::uint32_t bits) noexcept
{
    float value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

inline float half_to_float_exact(std::int16_t input) noexcept
{
    const std::uint16_t half = static_cast<std::uint16_t>(input);
    const std::uint32_t sign = static_cast<std::uint32_t>(half & 0x8000u) << 16;
    const std::uint32_t exponent = (half >> 10) & 0x1fu;
    std::uint32_t mantissa = half & 0x03ffu;
    std::uint32_t bits;

    if (exponent == 0u)
    {
        if (mantissa == 0u)
        {
            bits = sign;
        }
        else
        {
            std::uint32_t exponent32 = 113u;
            while ((mantissa & 0x0400u) == 0u)
            {
                mantissa <<= 1;
                --exponent32;
            }
            mantissa &= 0x03ffu;
            bits = sign | (exponent32 << 23) | (mantissa << 13);
        }
    }
    else if (exponent == 31u)
    {
        bits = sign | 0x7f800000u | (mantissa << 13);
    }
    else
    {
        bits = sign | ((exponent + 112u) << 23) | (mantissa << 13);
    }
    return bits_to_float(bits);
}

void convert_exact(const std::int16_t *input, float *output, std::size_t count) noexcept
{
    for (std::size_t index = 0; index < count; ++index)
    {
        output[index] = half_to_float_exact(input[index]);
    }
}

#if SIMPLE_HAS_F16C_DISPATCH
__attribute__((target("avx,f16c")))
void convert_f16c(const std::int16_t *input, float *output, std::size_t count) noexcept
{
    std::size_t index = 0;
    for (; index + 8 <= count; index += 8)
    {
        __m128i packed;
        std::memcpy(&packed, input + index, sizeof(packed));
        const __m256 expanded = _mm256_cvtph_ps(packed);
        _mm256_storeu_ps(output + index, expanded);
    }
    _mm256_zeroupper();
    convert_exact(input + index, output + index, count - index);
}

__attribute__((target("avx,f16c")))
void encode_f16c(const float *input, std::int16_t *output, std::size_t count) noexcept
{
    for (std::size_t index = 0; index < count; index += 8)
    {
        const __m256 expanded = _mm256_loadu_ps(input + index);
        const __m128i packed = _mm256_cvtps_ph(
            expanded, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
        std::memcpy(output + index, &packed, sizeof(packed));
    }
    _mm256_zeroupper();
}

converter_type select_converter() noexcept
{
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx") && __builtin_cpu_supports("f16c"))
    {
        return convert_f16c;
    }
    return convert_exact;
}

encoder_type select_encoder() noexcept
{
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx") && __builtin_cpu_supports("f16c"))
    {
        return encode_f16c;
    }
    return nullptr;
}
#else
converter_type select_converter() noexcept
{
    return convert_exact;
}

encoder_type select_encoder() noexcept
{
    return nullptr;
}
#endif

const converter_type selected_converter = select_converter();
const encoder_type selected_encoder = select_encoder();
}

extern "C" void simple_float16_to_float32(const std::int16_t *input, float *output,
                                           std::size_t count) noexcept
{
    selected_converter(input, output, count);
}

extern "C" std::size_t simple_float32_to_float16_f16c(const float *input,
                                                        std::int16_t *output,
                                                        std::size_t count) noexcept
{
    const std::size_t vector_count = count - count % 8;
    if (selected_encoder == nullptr || vector_count == 0)
    {
        return 0;
    }
    selected_encoder(input, output, vector_count);
    return vector_count;
}

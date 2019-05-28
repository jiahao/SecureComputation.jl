"""
CPU Random Number Generators

Provided by
- Intel's Secure Key Digital RNG
- AMD's Cryptographic Co-Processor

References:
https://software.intel.com/en-us/articles/intel-digital-random-number-generator-drng-software-implementation-guide
https://www.amd.com/system/files/TechDocs/amd-random-number-generator.pdf
"""
module CPURNGs

export RdRand, RdSeed

using Random
import Base: rand

include("cpuid.jl")
using .CPUID

abstract type CPURNG <: AbstractRNG end
struct RdRand <: CPURNG end
struct RdSeed <: CPURNG end
NativeBitsTypes = Int64, Int32, Int16, UInt64, UInt32, UInt16

Random.rng_native_52(::CPURNG) = UInt64

for (issupported, rng, RNGType, stepfnname) in (
    (CPUID.hasRdRand(), :rdrand, RdRand, :rdrand_step),
    (CPUID.hasRdSeed(), :rdseed, RdSeed, :rdseed_step))

if issupported

@eval begin
"""
$($stepfnname)(T) -> v::T, cf::UInt32

Return a value v from the hardware assembly instruction.

The second return value, cf, encodes the carry flag value set by the
instruction. It is 0xffffffff if the carry flag is set, indicating that the
first return value is a valid random number and 0x00000000 if the carry
flag is not set, indicating that the first return value is 0 and there was
not enough hardware entropy to return a random number.
"""
function $stepfnname end

"""
rand($($RNGType)(), T) -> v::T

Return a value from the hardware random number generator.
"""
function rand(::$RNGType, ::Random.SamplerType{T}) where {T <: Union{NativeBitsTypes...}}
    for i=1:100
        x, c = $stepfnname(T)
        if c==0xffffffff
            return x
        end
    end
end
end #eval

for T in NativeBitsTypes
    nb = 8 * sizeof(T)
    t = "i$nb"
    @eval begin
    function $stepfnname(::Type{$T})
        # HACK: if we pass a homogeneous tuple, Julia uses [2 x iXX] instead of
        # {iXX, iXX}, which is incompatible with the InlineAsm signature.
        # Add in another argument to force the tuple to be heterogeneous.
        a, b, _ = Base.llvmcall((
            $"declare {$t, i32} @llvm.x86.$rng.$nb()",
            $"""%call = call { $t, i32, i64 } asm sideeffect "$rng \$0; sbb \$1,\$1; xor \$2, \$2;", "=r,=r,=r"()
            ret { $t, i32, i64 } %call"""),
        Tuple{$T, UInt32, Int64}, Tuple{})
        return $T(a), UInt32(b)
    end
    end #eval
end
else
    @eval rand(::$RNGType, ::Random.SamplerType{T}) where T =
        error("Your CPU does not support ", $RNGType)
end
end #rng

end #module

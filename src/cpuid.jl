"""
Functions for querying the CPUID assembly instruction
"""
module CPUID

"""
cpuid_llvm(eax, ecx) -> eax, ebx, ecx, edx

Call CPUID assembly instruction, returning the values set in CPU registers.

Reference: See, for example,
https://c9x.me/x86/html/file_module_x86_id_45.html

Source for code: https://github.com/m-j-w/CpuId.jl
"""
@noinline cpuid_llvm(eax::UInt32, ecx::UInt32) =
    Base.llvmcall("""
        ; eax = %0, ecx = %1, %2 is some label
        ; call 'cpuid' with arguments loaded into registers EAX = eax, ECX = ecx
        %3 = tail call { i32, i32, i32, i32 } asm sideeffect "cpuid",
            "={ax},={bx},={cx},={dx},{ax},{cx},~{dirflag},~{fpsr},~{flags}"
            (i32 %0, i32 %1) #2
        ; retrieve the result values and convert to vector [4 x i32]
        %4 = extractvalue { i32, i32, i32, i32 } %3, 0
        %5 = extractvalue { i32, i32, i32, i32 } %3, 1
        %6 = extractvalue { i32, i32, i32, i32 } %3, 2
        %7 = extractvalue { i32, i32, i32, i32 } %3, 3
        ; return the values as a new tuple
        %8  = insertvalue [4 x i32] undef, i32 %4, 0
        %9  = insertvalue [4 x i32]   %8 , i32 %5, 1
        %10 = insertvalue [4 x i32]   %9 , i32 %6, 2
        %11 = insertvalue [4 x i32]  %10 , i32 %7, 3
        ret [4 x i32] %11
    """
    # llvmcall requires actual types, rather than the usual (...) tuple
    , NTuple{4,UInt32}, Tuple{UInt32,UInt32}
    , eax, ecx)
cpuid_llvm(eax=UInt32(0), ecx=UInt32(0)) =
    cpuid_llvm(UInt32(eax), UInt32(ecx))

"""
Check if CPU manufacturer is Intel or AMD
"""
function isIntelOrAMD()
    eax, ebx, ecx, edx = cpuid_llvm()
    ManufacturerID = String(reinterpret(UInt8, [ebx, edx, ecx]))
    (ManufacturerID == "GenuineIntel") | (ManufacturerID == "AuthenticAMD")
end

"""
Check if the RdRand instruction is supported by current CPU
"""
function hasRdRand()
    eax, ebx, ecx, edx = cpuid_llvm(1)
    return isIntelOrAMD() & Bool((ecx >> 30) & 1)
end

"""
Check if the RdSeed instruction is supported by current CPU
"""
function hasRdSeed()
    eax, ebx, ecx, edx = cpuid_llvm(7)
    return isIntelOrAMD() & Bool((ebx >> 18) & 1)
end

end #module

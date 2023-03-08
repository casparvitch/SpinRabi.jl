module Run

using FromFile: @from
@from "./SpinRabi.jl" import SpinRabi

if true
    a = 1e-3
    b = 9e-4
    𝔹₁ = (0, b, sqrt(abs(a^2 - b^2))) # T 
    𝔹₂ = (1, 0, 0) # T (keep it arbitrary) 
    D = 0
    E = 0
    S = 1 // 2

    rabi_freqs, esr_freqs, rabi_ids =
        SpinRabi.rabi_freqs(𝔹₁, 𝔹₂, D, E, S::Rational)
    println(abs.(rabi_freqs))
    println(esr_freqs)
    println(rabi_ids)
end

println()
println("---")
println()

if true
    a = 1e-3
    b = 9e-4
    𝔹₁ = (0, b, sqrt(abs(a^2 - b^2))) # T 
    𝔹₂ = (1, 0, 0) # T (keep it arbitrary) 
    D = 10e6
    E = 0
    S = 1 // 1
    
    rabi_freqs, esr_freqs, rabi_ids =
        SpinRabi.rabi_freqs(𝔹₁, 𝔹₂, D, E, S::Rational)
    println(abs.(rabi_freqs))
    println(esr_freqs)
    println(rabi_ids)
end

end # module

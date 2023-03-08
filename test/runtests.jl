using FromFile: @from
@from "../src/SpinRabi.jl" import SpinRabi
using Test
import LinearAlgebra: eigen

@testset "spinop" begin
    @test SpinRabi.spinop(1 // 2)[1] ≈ (1 / 2) * [0 1; 1 0] atol = 1e-5
    @test SpinRabi.spinop(1 // 2)[2] ≈ (1 / 2) * [0 -1im; 1im 0] atol = 1e-5
    @test SpinRabi.spinop(1 // 2)[3] ≈ (1 / 2) * [1 0; 0 -1] atol = 1e-5
    @test SpinRabi.spinop(1 // 1)[1] ≈ (1 / √(2)) * [0 1 0; 1 0 1; 0 1 0] atol =
        1e-5
    @test SpinRabi.spinop(1 // 1)[2] ≈
          (1 / √(2)) * [0 -1im 0; 1im 0 -1im; 0 1im 0] atol = 1e-5
    @test SpinRabi.spinop(1 // 1)[3] ≈ [1 0 0; 0 0 0; 0 0 -1] atol = 1e-5
end

@testset "calc_esr_freqs" begin
    esr_freqs = SpinRabi.calc_esr_freqs(-2e-3, 2.87e9, 0, 1 // 1)[1]
    sesr = sort(esr_freqs, rev = true)
    @test (sesr[1] - sesr[2]) / (2 * 28e9) ≈ 2e-3 atol = 1e-2

    esr_freqs2 = SpinRabi.calc_esr_freqs(2e-3, 2.87e9, 0, 1 // 1)[1]
    sesr2 = sort(esr_freqs2, rev = true)
    @test (sesr2[1] - sesr2[2]) / (2 * 28e9) ≈ 2e-3 atol = 1e-2

    esr_freqs3 = SpinRabi.calc_esr_freqs(-10e-3, 2.87e9, 0, 1 // 1)[1]
    sesr3 = sort(esr_freqs3, rev = true)
    @test (sesr3[1] - sesr3[2]) / (2 * 28e9) ≈ 10e-3 atol = 1e-2

    esr_freqs4 = SpinRabi.calc_esr_freqs(0.2e-3, 10e6, 50e6, 3 // 2)[1]
    sesr4 = sort(esr_freqs4, rev = true)
    @test (sesr4[1] - sesr4[2]) / (2 * 28e9) ≈ 0.2e-3 atol = 1e-2

    esr_freqs5 = SpinRabi.calc_esr_freqs(-0.2e-3, 10e6, 50e6, 3 // 2)[1]
    sesr5 = sort(esr_freqs5, rev = true)
    @test (sesr5[1] - sesr5[2]) / (2 * 28e9) ≈ 0.2e-3 atol = 1e-2
end

@testset "simple_x-drive_rabi" begin
    𝔹₁ = (0, 0, 1e-3) # T 
    𝔹₂ = (1, 0, 0) # T (keep it arbitrary) 
    D = 2.87e9
    E = 0
    S = 1 // 1

    rabi_freqs, _, _ = SpinRabi.rabi_freqs(𝔹₁, 𝔹₂, D, E, S::Rational)
    # 2 * gamma [Hz/T] * √2S / 2 [S-factor] * (1/2π) [convert to linear freq]
    @test rabi_freqs[1] ≈ 2 * 14e9 * (√2 / 2) / (2 * π) atol = 1e-2
end

@testset "simple_y-drive_rabi" begin
    𝔹₁ = (0, 0, 1e-3) # T 
    𝔹₂ = (0, 1, 0) # T (keep it arbitrary) 
    D = 2.87e9
    E = 0
    S = 1 // 1

    rabi_freqs, _, _ = SpinRabi.rabi_freqs(𝔹₁, 𝔹₂, D, E, S::Rational)
    # 2 * gamma [Hz/T] * √2S / 2 [S-factor] * (1/2π) [convert to linear freq]
    @test rabi_freqs[1] ≈ 2 * 14e9 * (√2 / 2) / (2 * π) atol = 1e-2
end

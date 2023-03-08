module SpinRabi

using LinearAlgebra

# ----------------------------------- ESR Stuff -------------------------------
# s_idxs is probably poor naming here, I mean ms_idxs ?
const 𝐠 = 2
# swap sign if you want ms=+1 to be higher energy that ms=-1
const 𝛄 = -1 * 𝐠 * 14e9 # Hz/T

# https://en.wikipedia.org/wiki/Spin_(physics)#Higher_spins
# basis ordering is mₛ = (1, 0, -1) etc.
function spinop(S::Rational)
    δ(x, y) = ==(x, y)
    σx = Matrix{ComplexF64}(undef, (Int(2 * S + 1), Int(2 * S + 1)))
    σy = Matrix{ComplexF64}(undef, (Int(2 * S + 1), Int(2 * S + 1)))
    σz = Matrix{ComplexF64}(undef, (Int(2 * S + 1), Int(2 * S + 1)))

    @inbounds for i in 1:Int(2 * S + 1)
        for j in 1:Int(2 * S + 1)
            sqrtfac = √(Complex((S + 1) * (i + j - 1) - i * j))
            σx[i, j] = (1 / 2) * (δ(i, j + 1) + δ(i + 1, j)) * sqrtfac
            σy[i, j] = (1im / 2) * (δ(i, j + 1) - δ(i + 1, j)) * sqrtfac
            σz[i, j] = (S + 1 - i) * δ(i, j)
        end
    end
    return [σx, σy, σz]
end

# ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
function rad2cart(r, θ, ϕ)
    return r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)
end

# t-dep hamiltonian method here: 
# https://www.nature.com/articles/s41467-019-09429-x    

function _ham_sd(𝔹, D, E, S::Rational)
    sx, sy, sz = spinop(S)
    # hmmm... any reason to not have a negative before γ there?
    return D * sz * sz + E * (sx * sx - sy * sy) + 𝛄 * sum(𝔹 .* [sx, sy, sz])
end

# new estates are returned ordered by increasing energy
function _esr_calculator(S, evals, evecs; selrules)
    ms_charac = _define_ms_charac(S, evecs)

    new_charac = _define_new_charac()
    s_idxs = collect(1:Int(2 * S + 1))

    esr_freqs = Array{Float64}(undef, 0)
    # vec of pairs identifying old/new estates for above
    esr_old_ids = Array{Set{Rational}}(undef, 0)
    esr_new_ids = Array{Set{Symbol}}(undef, 0)
    for s_idx in s_idxs
        for other_s_idx in s_idxs
            this_energy = evals[s_idx]
            other_energy = evals[other_s_idx]
            this_ms = -1 * (ms_charac[s_idx] - 1 - S)
            other_ms = -1 * (ms_charac[other_s_idx] - 1 - S)
            this_new = new_charac[s_idx]
            other_new = new_charac[other_s_idx]
            if other_s_idx == s_idx || Set((this_new, other_new)) in esr_new_ids
                continue
            else
                if !selrules # simple mode... (no selection rules)
                    push!(esr_freqs, abs(this_energy - other_energy))
                    push!(esr_evecs)
                    push!(esr_old_ids, Set((this_ms, other_ms)))
                    push!(esr_new_ids, Set((this_new, other_new)))
                else
                    # check for angular momentum selection rules
                    # on the 'old' identities only
                    # only valid when there's minimal state-mixing 
                    #   (e.g. low off-axis 𝔹)
                    if abs(this_ms - other_ms) != 1
                        continue
                    else
                        push!(esr_freqs, abs(this_energy - other_energy))
                        push!(esr_old_ids, Set((this_ms, other_ms)))
                        push!(esr_new_ids, Set((this_new, other_new)))
                    end
                end
            end
        end
    end
    return esr_freqs, esr_old_ids, esr_new_ids
end

function calc_esr_freqs(𝔹, D, E, S::Rational; selrules::Bool = true)
    ham = _ham_sd(𝔹, D, E, S)
    evals, evecs = eigen(ham)

    return _esr_calculator(S, evals, evecs, selrules = selrules)
end

function calc_esr_sys(𝔹, D, E, S::Rational; selrules::Bool = true)
    ham = _ham_sd(𝔹, D, E, S)
    evals, evecs = eigen(ham)

    esr_freqs, esr_old_ids, esr_new_ids =
        _esr_calculator(S, evals, evecs, selrules = selrules)
    esr_evec_map = Dict(
        _define_new_charac()[i] => evec for
        (i, evec) in enumerate(eachcol(evecs))
    )
    return esr_freqs, esr_evec_map, esr_old_ids, esr_new_ids
end

function _define_new_charac()
    return [
        :α,
        :β,
        :γ,
        :δ,
        :ϵ,
        :ζ,
        :η,
        :θ,
        :ι,
        :κ,
        :λ,
        :μ,
        :ν,
        :ξ,
        :ο,
        :π,
        :ρ,
        :σ,
        :τ,
        :υ,
        :χ,
        :ψ,
        :ω,
    ]
end

function _define_ms_charac(S, evecs)
    s_idxs = collect(1:Int(2 * S + 1))
    # idx in this array represents new_ms, value is idx in old_ms array
    what_old_ms_is_this_new_s_most_like =
        Array{Union{Nothing, Number}}(nothing, Int(2 * S + 1))
    # print("+++ "); println(-1 .* (collect(1:Int(2 * S + 1)) .- S .- 1))

    # very messy algorithm...
    # iterate over new estate indexes, find candidate it is most like (by
    # projection onto evectors). Assume 1:1 old-new estate transform.

    for new_idx in s_idxs
        projs = abs.(evecs[:, new_idx]) # projection onto evec bases
        old_s_candidates = collect(1:Int(2 * S + 1)) # candidates (idx in s_idxs)
        while true
            _, most_similar_ms_idx = findmax(projs)
            old_s = old_s_candidates[most_similar_ms_idx]
            if isnothing(what_old_ms_is_this_new_s_most_like[new_idx]) &&
               !(old_s in what_old_ms_is_this_new_s_most_like)
                what_old_ms_is_this_new_s_most_like[new_idx] = old_s
                break # ok we found 'it', break out *to* for loop
            else
                # println("haha I'm in danger")
                # this new_s basis already used, retry after deleting this optn
                deleteat!(projs, most_similar_ms_idx)
                deleteat!(old_s_candidates, most_similar_ms_idx)
            end
        end
    end
    return what_old_ms_is_this_new_s_most_like
end

# ----------------------------------- Rabi stuff -------------------------------

# return in Hz
function rabi_freqs(B₁, B₂, D, E, S::Rational; selrules::Bool = true)
    esr_freqs, esr_evec_map, esr_old_ids, esr_new_ids =
        calc_esr_sys(B₁, D, E, S, selrules = selrules)

    H₂ = _ham_sd(B₂, 0, 0, S::Rational)

    rb_freqs = Array{ComplexF64}(undef, 0)
    rabi_ids = Array{Set{Symbol}}(undef, 0)
    for (i, esr_freq) in enumerate(esr_freqs)
        ids = collect(esr_new_ids[i])
        evec_g = esr_evec_map[ids[1]]
        evec_e = esr_evec_map[ids[2]]
        esr_freq = evec_g' * H₂ * evec_e / (2 * π)
        push!(rb_freqs, esr_freq)
        push!(rabi_ids, Set((ids[1], ids[2])))
    end

    return rb_freqs, esr_freqs, rabi_ids
end

end # module

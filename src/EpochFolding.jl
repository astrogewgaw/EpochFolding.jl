module EpochFolding

mutable struct FoldStatistics
    N::Int64
    n::Int64
    μ::Float64
    σ²::Float64
    μₚ::Float64
    σₚ²::Float64
    χᵣ²::Float64
end

FoldStatistics() = FoldStatistics(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

function isone(x)
    abs((x) - 1.0) < 1e-12
end

function sprinkle(Y, Ybuf, N, ϕlow, Δϕ, x, ϕadded)
    Δξ = 1 / N
    valperϕ = x / Δϕ
    ξcurr = ϕlow * N
    while (Δϕ > 1e-12)
        ξᵢ = floor(Int, ξcurr + 1e-12)
        δϕ = ((ξᵢ + 1.0) - ξcurr) * Δξ
        ξᵢ = mod(ξᵢ, N)
        if (δϕ > Δϕ)
            δϕ = Δϕ
        end
        ϕleft = 1 - ϕadded
        fracδϕ = (ϕleft > δϕ) ? δϕ : ϕleft
        Ybuf[ξᵢ+1] += fracδϕ * valperϕ
        ϕadded += fracδϕ
        if isone(ϕadded)
            for j ∈ eachindex(Y)
                Y[j] += Ybuf[j]
            end
            for j ∈ eachindex(Ybuf)
                Ybuf[j] = 0.0
            end
            fracδϕ = δϕ - ϕleft
            Ybuf[ξᵢ+1] += fracδϕ * valperϕ
            ϕadded = fracδϕ
        end
        Δϕ -= δϕ
        ξcurr += δϕ * N
    end
    ϕadded
end

function delta(Y, Ybuf, N, ϕlow, Δϕ, x)
    ξ = (Int(floor((ϕlow + 0.5 * Δϕ) * N + 1e-12)) % N) + 1
    Ybuf[ξ] += 1.0
    Y[ξ] += x
end

function fold!(
    Y,
    X;
    δt,
    f₀,
    t₀=0.0,
    ϕ₀=0.0,
    fd=0.0,
    fdd=0.0,
    standard=true,
)
    n = size(Y, 1)

    stats = FoldStatistics()
    Ybuf = zeros(eltype(Y), n)

    ϕ = 0.0
    ϕlow = 0.0
    ϕnext = 0.0
    ϕadded = 0.0

    fd /= 2.0
    fdd /= 6.0
    stats.n = n
    stats.σ² *= (stats.N - 1.0)

    ϕ = t₀ * (t₀ * (t₀ * fdd + fd) + f₀) + ϕ₀
    ϕlow = (ϕ < 0.0) ? 1.0 + modf(ϕ)[1] : modf(ϕ)[1]
    for i ∈ eachindex(X)
        tnext = t₀ + i * δt
        ϕnext = tnext * (tnext * (tnext * fdd + fd) + f₀) + ϕ₀
        Δϕ = ϕnext - ϕ
        if standard
            ϕadded = sprinkle(Y, Ybuf, n, ϕlow, Δϕ, X[i], ϕadded)
        else
            delta(Y, Ybuf, n, ϕlow, Δϕ, X[i])
        end
        ϕhigh = ϕlow + Δϕ
        ϕlow = modf(ϕhigh)[1]
        ϕ = ϕnext

        stats.N += 1.0
        σ = X[i] - stats.μ
        stats.μ += σ / stats.N
        stats.σ² += σ * (X[i] - stats.μ)
    end
    for i ∈ eachindex(Y)
        stats.μₚ += Y[i]
    end
    stats.μₚ /= n

    stats.χᵣ² = 0.0
    for i ∈ eachindex(Y)
        stats.χᵣ² += (Y[i] - stats.μₚ)^2
    end
    stats.σ² /= (stats.N - 1.0)
    stats.σₚ² = stats.σ² * stats.N * (1 / n)
    stats.χᵣ² /= (stats.σₚ² * (n - 1))
    ϕnext = (ϕnext < 0.0) ? 1.0 + modf(ϕnext)[1] : modf(ϕnext)[1]
    stats, ϕnext
end

function fold(
    X;
    δt,
    f₀,
    n=64,
    t₀=0.0,
    ϕ₀=0.0,
    fd=0.0,
    fdd=0.0,
    standard=true,
)
    Y = zeros(Float64, n)
    stats, ϕnext = fold!(Y, X; δt, f₀, t₀, ϕ₀, fd, fdd, standard)
    Y, stats, ϕnext
end

end

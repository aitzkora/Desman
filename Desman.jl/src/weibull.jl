"""
computes F₀(t_{ij}^r/w*exp(σ.z); θ)-F₀(t_{ij}^d/w*exp(σ.z); θ)
where μ is the variable vectors
"""
function weibull_diff(σ,μ, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=exp(sum(σ.*μ[4:end])) # kontuz σ must be a linear form 1xn
  θ = w*μ[1]*σ_γ
  α = μ[2]
  fd =exp(-(xd/θ)^α)
  if (xg <= zero(T))
    fg = zero(T)
  else
    fg = exp(-(xg/θ)^α)
  end
  return fd-fg
end

"""
computes ∇ᵤG = [∂μ1G, ∂αG, 0, ∂zG]ᵀ
"""
function grad_weibull_diff(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false},
                           μ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=exp(sum(σ.*μ[4:end]))
  ∂₁θ = w*σ_γ
  θ = μ[1]*∂₁θ
  α = μ[2]
  fd = exp(-(xd/θ)^α)
  g = zeros(T, length(μ))
  g[1] = fd*α*(xd/θ)^α*∂₁θ/θ
  g[4:end]= fd*α*(xd/θ)^α*σ[:]
  if (xd != zero(T))
    g[2] = -fd*(xd/θ)^α*log(xd/θ)
  else
    g[2] = zero(T)
  end
  if (xg != Inf)
    fg = exp(-(xg/θ)^α) 
    g[1] += -fg*α*(xg/θ)^α*∂₁θ/θ
    if (xg != zero(T))
      g[2] += fg*(xg/θ)^α*log(xg/θ)
    else
      g[2] = zero(T)
    end
    g[4:end] -= fg*α*(xg/θ)^α*σ[:]
  end
  g[3] = zero(T)
  return g
end

"""
computes ∇²ᵤG = 
|∂μ1(∂μ1G) | ∂α(∂μ1)G | 0|   ∂z(∂μ1G) | 
|  ∂μ1(∂αG)|      ∂α²G| 0|   ∂z(∂αG)  |
|         0|         0| 0|        0   |
|  ∂μ1(∂zG)|   ∂α(∂zG)| 0|     ∂z²G   |
with analytical formulas : some bench
- @benchmark hessian_weibull_diff(σ, λ, w, xg, xd)
  - Time  (mean ± σ):   2.441 μs ±  16.341 μs  
  - Memory estimate: 2.62 KiB, allocs estimate: 41.

- @benchmark FiniteDiff.finite_difference_hessian(g,λ)
  - Time  (mean ± σ):   50.454 μs ± 200.699 μs  
  - Memory estimate: 32.45 KiB, allocs estimate: 911.

- @benchmark ForwardDiff.hessian(g,λ)
  - Time  (mean ± σ):   10.723 μs ±  38.284 μs  
  - Memory estimate: 19.32 KiB, allocs estimate: 18.
- => analytical method is faster than other
- FIXME : why it has more allocation than ForwardDiff
"""
function hessian_weibull_diff(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false},
                      μ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  n = size(μ, 1)

  # alloc and easy part filling
  H = zeros(T, n, n)
  H[3,:] .= zero(T)
  H[:,3] .= zero(T)

  # reused values
  σ_γ=exp(sum(σ.*μ[4:end]))
  ∂₁θ = w*σ_γ
  θ = μ[1]*∂₁θ 
  α = μ[2]
  ∂zθ = θ.*σ

  fd =exp(-(xd/θ)^α)
  if (xg <= zero(T))
    fg = zero(T)
  else
    fg = exp(-(xg/θ)^α)
  end

  g1d = fd*α*(xd/θ)^α*∂₁θ/θ
  g1g = fg*α*(xg/θ)^α*∂₁θ/θ

  v_block = zero(T)

  if (xd != zero(T))
    # upper left value
    H[1,1]     = g1d * α*(xd/θ)^α               * ∂₁θ/θ - 
                  fd * α*α*(xd/θ)^α*∂₁θ/θ       * ∂₁θ/θ - 
                  fd * α*(xd/θ)^α               * (∂₁θ/θ)^2
    # second off diagonal values
    H[1,2]  = fd*∂₁θ/θ*(xd/θ)^α*(T(1)+α*log(xd/θ)*(T(1)-(xd/θ)^α))

    # second diagonal value
    H[2,2]  = -fd * log(xd/θ)^2 * (xd/θ)^α * (T(1)-(xd/θ)^α)

    # last first row value
    H[1,4:end]  = -α^2 * fd * ∂₁θ/θ*(xd/θ)^α *(T(1) - (xd/θ)^α).*σ[:]

    # last second row value
    H[2,4:end]  = fd * (xd/θ)^α * (T(1) + α * log(xd/θ) * (T(1) - (xd/θ)^α)) .* σ[:]

    # last block, rank one tensor
    v_block  = - α^2* fd * (xd/θ)^α * (T(1) - (xd/θ)^α)
  end
  if (xg!=Inf) && (xg != zero(T))
    # second off diagonal values
    H[1,1]    -=(g1g * α*(xg/θ)^α               * ∂₁θ/θ - 
                  fg * α*α*(xg/θ)^α*∂₁θ/θ       * ∂₁θ/θ - 
                  fg * α*(xg/θ)^α               * (∂₁θ/θ)^2)
    # upper left value
    H[1,2] -= fg*∂₁θ/θ*(xg/θ)^α*(T(1)+α*log(xg/θ)*(T(1)-(xg/θ)^α))

    # second diagonal value
    H[2,2] -= -fg * log(xg/θ)^2 * (xg/θ)^α * (T(1)-(xg/θ)^α)
             
    # last first row value
    H[1,4:end] -= -α^2 * fg * ∂₁θ/θ*(xg/θ)^α *(T(1) - (xg/θ)^α).*σ[:]

    # last second row value
    H[2,4:end] -= fg * (xg/θ)^α * (T(1) + α * log(xg/θ) * (T(1) - (xg/θ)^α)) .* σ[:]

    # last block, rank one tensor
    v_block += α^2* fg * (xg/θ)^α * (T(1) - (xg/θ)^α)
  end

  H[4:end,4:end] = v_block * σ[:] * σ[:]'
  # symmetrization
  H[2,1]  = H[1,2]
  H[4:end, 1] = view(H, 1, 4:n)
  H[4:end, 2] = view(H, 2, 4:n)

  return H
end


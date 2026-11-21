"""
computes F₀(t_{ij}^r/w*exp(σ.z); θ)-F₀(t_{ij}^d/w*exp(σ.z); θ)
where μ is the variable vectors
"""
function weibull_diff(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false}
                     ,μ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
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
computes ∇²ᵤG = [∂μ1(∂μ1G), ∂α(∂μ1)G, 0,   ∂z(∂μ1G)
                  ∂μ1(∂αG),      ∂α²G, 0,   ∂z(∂αG)
                         0,         0, 0,        0
                  ∂μ1(∂zG),   ∂α(∂zG), 0,     ∂z²G]
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

  # upper left value
  h[1,1]     = g1d * α*(xd/θ)^α               * ∂₁θ/θ - 
                fd * α*α*(xd/θ)^(α-1)*∂₁θ/θ^2 * ∂₁θ/θ - 
                fd * α*(xd/θ)^α               * ∂₁θ^2/θ^2
  h[1,1]    -=(g1g * α*(xg/θ)^α               * ∂₁θ/θ - 
                fg * α*α*(xg/θ)^(α-1)*∂₁θ/θ^2 * ∂₁θ/θ - 
                fg * α*(xg/θ)^α               * ∂₁θ^2/θ^2)


  # second off diagonal values
  h[1,2]     = α* ∂₁θ/θ * (-fd * (xd/θ)^α * log(xd/θ) +
                            fd * (xd/θ)^α * log(xd/θ))

  h[1,2]    -= α* ∂₁θ/θ * (-fg * (xg/θ)^α * log(xg/θ) +
                            fg * (xg/θ)^α * log(xg/θ))

  h[2,1] = h[1,2]

  # second diagonal value
  h[2,2]  = -log(xd/θ)*(-fd*(xd/θ)^α * log(xd/θ)
                        -fd*(xd/θ)^α * log(xd/θ))
  h[2,2] -= -log(xd/θ)*(-fg*(xg/θ)^α * log(xg/θ)
                        -fg*(xg/θ)^α * log(xg/θ))

  # last first row value
  h[1,4:end] = α*(fd*(xd/θ)^α*∂₁θ/θ
                  fd*(xd/θ)^α*∂₁θ/θ
                  fd*(xd/θ)^α*∂₁θ/θ)

  return H
end


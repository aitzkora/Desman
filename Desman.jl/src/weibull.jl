function weibull_diff(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false}
                     ,λ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=exp(sum(σ.*λ[4:end])) # kontuz σ must be a linear form 1xn
  θ = w*λ[1]*σ_γ
  α = λ[2]
  fd =exp(-(xd/θ)^α)
  if (xg <= 0.0)
    fg = 0.0
  else
    fg = exp(-(xg/θ)^α)
  end
  return fd-fg
end

function grad_weibull(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false},
                      λ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=exp(sum(σ.*λ[4:end]))
  θ = w*λ[1]*σ_γ
  α = λ[2]
  g = zeros(3+length(λ[4:end]))
  g[1] = (exp(-(xd/θ)^α)*α*(xd/θ)^α/θ)*w*σ_γ
  g[4:end]= exp(-(xd/θ)^α)*α*(xd/θ)^α*σ[:]
  if (xd != 0.0)
    g[2] = (-exp(-(xd/θ)^α)*(xd/θ)^α*log(xd/θ))
  else
    g[2] = 0.0
  end
  if (xg!=Inf)
    g[1] += -w*exp(-(xg/θ)^α)*α*(xg/θ)^α/θ*σ_γ
    if (xg != 0.0)
      g[2] += exp(-(xg/θ)^α)*(xg/θ)^α*log(xg/θ)
    else
      g[2] = 0.0
    end
    g[4:end] -= exp(-(xg/θ)^α)*α*(xg/θ)^α*σ[:]
  end
  g[3] = zero(T)
  return g
end


using QuadGK
using FastGaussQuadrature
using DataFrames
using Distributions
using Integrals
using SpecialFunctions

function weibull_diff(σ::Matrix{T}, λ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=sum(σ*λ[4:end]) # kontuz σ must be a linear form 1xn
  θ = w*(λ[1]+σ_γ)
  α = λ[2]
  return exp(-(xd/θ)^α)-exp(-(xg/θ)^α)
end

function grad_weibull(σ::Matrix{T}, λ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=sum(σ*λ[4:end])
  θ = w*(λ[1]+σ_γ)
  α = λ[2]
  g = zeros(3+length(λ[4:end]))
  g[1] = (exp(-(xd/θ)^α)*α*(xd/θ)^α/θ)*w
  g[2] = (-exp(-(xd/θ)^α)*(xd/θ)^α*log(xd/θ))
  if (xg!=Inf)
    g[1] += -w*exp(-(xg/θ)^α)*α*(xg/θ)^α/θ
    g[2] += exp(-(xg/θ)^α)*(xg/θ)^α*log(xg/θ)
  end
  g[3] = zero(T)
  g[4:end]= g[1] .* σ[:]
  return g
end

function prod_weibull(bio::Biotope, i::Int64, selVar::Vector{Int64}, w::T, λ::Vector{T}) where {T<:Real}
  idx = collect([size(bio.df,2)-size(bio.Σ,2)+2:size(bio.df,2);])
  γ = λ[4:end]
  z = zeros(T,length(bio.idx[i]))
  for (jj, j) in enumerate(bio.idx[i])
      z[jj] = weibull_diff(Matrix{T}(bio.df[j:j,idx[selVar]]), λ, w, T(bio.df.durationd[j]), T(bio.df.durationg[j]))
  end 
  return prod(z)
end

"""
pdf_frailty(bio, i, selVar)
computes the density of the `i`-th latrine against the `selVar` collection of covariates
"""
function pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}, ::Type{T}) where {T<:Real}
  function proxyDistCov(w::T, λ::Vector{T})
    return prod_weibull(bio,i, selVar,w,λ) * Γh(λ[3], w)
  end
  return proxyDistCov
end

# default type is Float64
pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}) = pdf_frailty(bio, i, selVar, Float64)


function grad_pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}, ::Type{T}) where {T<:Real}
  function proxyGrad(w::T, λ::Vector{T})
    idx = collect([size(bio.df,2)-size(bio.Σ,2)+2:size(bio.df,2);])
    p_wei = prod_weibull(bio, i, selVar, w, λ)
    ∂pdf = zeros(size(λ))
    for (jj, j) in enumerate(bio.idx[i])
        σ = Matrix{T}(bio.df[j:j,idx[selVar]])
        p_exclu_j = p_wei ./ weibull_diff(σ, λ, w, T(bio.df.durationd[j]), T(bio.df.durationg[j]))
        ∂pdf += p_exclu_j .* grad_weibull(σ, λ, w, T(bio.df.durationd[j]), T(bio.df.durationg[j]))
    end 
    ∂pdf *= Γh(λ[3], w)
    ∂pdf[3] = p_wei * ∂ᵤΓh(λ[3], w)
    return ∂pdf
  end
  return proxyGrad
end

grad_pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}) =  grad_pdf_frailty(bio, i, selVar, Float64)

"""
Γh is a optimized version of x->pdf(Gamma(u,1/u),x)
Γh(3.4,5.2) -> 25.503 ns (0 allocations: 0 bytes)
pdf(Gamma(3.4,1. / 3.4), 5.2) ->  235.772 ns (2 allocations: 48 bytes)
"""
function Γh(u::T, w::T) where {T<:Real}
  return w^(u-1)*exp(-w*u)*u^u/gamma(u)
end 

function ∂ᵤΓh(u::T, w::T) where {T<:Real}
  return Γh(u,w)*(log(w*u)-w-digamma(u)+1)
end

"""
logLikelihood(b::Biotope, selVar::Vector{Int64}, methodInt = GaussLegendre)
computes the logLikelihood conditioned by a subset of covariates given by `selVar`
"""
function logLikelihood(b::Biotope, selVar::Vector{Int64}, methodInt = GaussLegendre)
  return λ->-sum(log.([Integrals.solve(IntegralProblem(pdf_frailty(b,i,selVar), (0.0, Inf), λ), methodInt())[1] for i=1:b.N]))
end

function grad_logLikelihood(b::Biotope, selVar::Vector{Int64}, methodInt = QuadGKJL)
    return λ->sum([solve(IntegralProblem(grad_pdf_frailty(b, i, selVar), (0.0, Inf), λ), methodInt()).u ./
                   solve(IntegralProblem(pdf_frailty(b, i, selVar), (0.0, Inf), λ), methodInt())[1] for i=1:b.N])
end

function logLikelihood2(b::Biotope, selVar::Vector{Int64}, rtol =1e-8)
  return λ->-sum(log.([quadgk(x->pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]))
end

function logLikelihood(b::Biotope)
  logLikelihood(b, Int64[])
end


using QuadGK
using FastGaussQuadrature
using DataFrames
using Distributions
using Integrals
using SpecialFunctions

function weibull_diff(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false}
                      ,λ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=sum(σ.*λ[4:end]) # kontuz σ must be a linear form 1xn
  θ = w*(λ[1]+σ_γ)
  α = λ[2]
  return exp(-(xd/θ)^α)-exp(-(xg/θ)^α)
end

function grad_weibull(σ::SubArray{T, 1, Matrix{T},Tuple{Int64,Vector{Int64}}, false},
                      λ::Vector{T}, w::T, xd::T, xg::T) where {T<:Real}
  σ_γ=sum(σ.*λ[4:end])
  θ = w*(λ[1]+σ_γ)
  α = λ[2]
  g = zeros(3+length(λ[4:end]))
  g[1] = (exp(-(xd/θ)^α)*α*(xd/θ)^α/θ)*w
  if (xd != 0.0)
    g[2] = (-exp(-(xd/θ)^α)*(xd/θ)^α*log(xd/θ))
  else
    g[2] = 0.0
  end
  if (xg!=Inf)
    g[1] += -w*exp(-(xg/θ)^α)*α*(xg/θ)^α/θ
    if (xg != 0.0)
      g[2] += exp(-(xg/θ)^α)*(xg/θ)^α*log(xg/θ)
    else
      g[2] = 0.0
    end
  end
  g[3] = zero(T)
  g[4:end]= g[1] .* σ[:]
  return g
end

function prod_weibull(bio::Biotope, i::Int64, selVar::Vector{Int64}, w::T, λ::Vector{T}) where {T<:Real}
  z = one(T)
  for j in bio.idxLat[i]
      z *= weibull_diff(view(bio.matCov,j,selVar), λ, w, T(bio.durd[j]), T(bio.durg[j]))
  end 
  return z
end

function prod_weibull_exclude(bio::Biotope, i::Int64, selVar::Vector{Int64}, w::T, λ::Vector{T}, k::Int64) where {T<:Real}
  z = one(T)
  for j in bio.idxLat[i]
    if (j != k)
        z *= weibull_diff(view(bio.matCov,j,selVar), λ, w, T(bio.durd[j]), T(bio.durg[j]))
      end
  end 
  return z
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
    idx = collect([bio.covStart:bio.covEnd;])
    p_wei = prod_weibull(bio, i, selVar, w, λ)
    ∂pdf = zeros(size(λ))
    for (jj, j) in enumerate(bio.idxLat[i])
        σ = view(bio.matCov, j, selVar)
        p_exclu_j = p_wei ./ weibull_diff(σ, λ, w, T(bio.durd[j]), T(bio.durg[j]))
        ∂pdf += p_exclu_j .* grad_weibull(σ, λ, w, T(bio.durd[j]), T(bio.durg[j]))
    end 
    ∂pdf *= Γh(λ[3], w)
    ∂pdf[3] = p_wei * ∂ᵤΓh(λ[3], w)
    return ∂pdf
  end
  return proxyGrad
end

grad_pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}) =  grad_pdf_frailty(bio, i, selVar, Float64)

function grad_pdf_frailty_check(bio::Biotope, i::Int64, selVar::Vector{Int64}, ::Type{T}) where {T<:Real}
  function proxyGrad(w::T, λ::Vector{T})
    p_wei = prod_weibull(bio, i, selVar, w, λ)
    ∂pdf = zeros(size(λ))
    for (jj, j) in enumerate(bio.idxLat[i])
        σ = view(bio.matCov, j, selVar)
        ∂pdf += prod_weibull_exclude(bio, i, selVar, w, λ, j) .* grad_weibull(σ, λ, w, T(bio.durd[j]), T(bio.durg[j]))
    end 
    ∂pdf *= Γh(λ[3], w)
    ∂pdf[3] = p_wei * ∂ᵤΓh(λ[3], w)
    return ∂pdf
  end
  return proxyGrad
end

grad_pdf_frailty_check(bio::Biotope, i::Int64, selVar::Vector{Int64}) =  grad_pdf_frailty_check(bio, i, selVar, Float64)

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

"""
in place function to compute f and/or ∇f if there are needeed
"""
function getfg!(b::Biotope, selVar::Vector{Int64})
  function fg!(F, G, λ::Vector{Float64})
    rtol=1.e-8  
    Li = [quadgk(x->pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
    if G !== nothing
      ∇Li =[quadgk(x->grad_pdf_frailty_check(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
      G[:] = -sum(∇Li./Li)
    end
    if F !== nothing
      F = -sum(log.(Li))
      return F
    end
  end
  return fg!
end

function getg!(b::Biotope, selVar::Vector{Int64})
  function g!(G, λ::Vector{Float64})
    rtol=1.e-8  
    Li = [quadgk(x->pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
    ∇Li =[quadgk(x->grad_pdf_frailty_check(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
    G[:] = -sum(∇Li./Li)
  end
  return g!
end

function grad_logLikelihood(b::Biotope, selVar::Vector{Int64}, methodInt = QuadGKJL)
    return λ->-sum([solve(IntegralProblem(grad_pdf_frailty_check(b, i, selVar), (0.0, Inf), λ), methodInt()).u ./
                   solve(IntegralProblem(pdf_frailty(b, i, selVar), (0.0, Inf), λ), methodInt())[1] for i=1:b.N])
end

function grad_logLikelihood_div(b::Biotope, selVar::Vector{Int64}, methodInt = QuadGKJL)
  return λ->-sum([solve(IntegralProblem(grad_pdf_frailty(b, i, selVar), (0.0, Inf), λ), methodInt()).u ./
                  solve(IntegralProblem(pdf_frailty(b, i, selVar), (0.0, Inf), λ), methodInt())[1] for i=1:b.N])
end


function logLikelihood(b::Biotope)
  logLikelihood(b, Int64[])
end


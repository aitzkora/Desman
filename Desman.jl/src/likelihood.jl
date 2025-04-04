using QuadGK
using FastGaussQuadrature
using DataFrames
using Distributions
using Integrals


function prod_weibull(bio::Biotope, i::Int64, selVar::Vector{Int64}, w::Float64, λ::Vector{Float64})
  idx = collect([size(bio.df,2)-size(bio.Σ,2)+2:size(bio.df,2);])
  γ = λ[4:end]
  z = zeros(length(bio.idx[i]))
  for (jj, j) in enumerate(bio.idx[i])
      sum_γ_j = sum(Matrix{Float64}(bio.df[j:j,idx[selVar]])*γ) 
      wei_j = x->ccdf(Weibull(exp(λ[2]), w.*exp(λ[1]+sum_γ_j)), x)
      z[jj] = wei_j(bio.df.durationd[j]) - wei_j(bio.df.durationg[j])
  end 
  return prod(z)
end

"""
pdf_frailty(bio, i, selVar)
computes the density of the `i`-th latrine against the `selVar` collection of covariables
"""
function pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64})
  function proxyDistCov(w::Float64, λ::Vector{Float64})
    return prod_weibull(bio,i, selVar,w,λ) * pdf(Gamma(exp(λ[3]), exp(-λ[3])),  w)
  end
  return proxyDistCov
end

"""
logLikelihood(b::Biotope, selVar::Vector{Int64}, methodInt = GaussLegendre)
computes the logLikelihood conditionned by a subset of covariables given by `selVar`
"""
function logLikelihood(b::Biotope, selVar::Vector{Int64}, methodInt = GaussLegendre)
  return λ->-sum(log.([Integrals.solve(IntegralProblem(pdf_frailty(b,i,selVar), (0.0, Inf), λ), methodInt())[1] for i=1:b.N]))
end

function logLikelihood2(b::Biotope, selVar::Vector{Int64}, rtol =1e-8)
  return λ->-sum(log.([quadgk(x->pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]))
end

function logLikelihood(b::Biotope)
  logLikelihood(b, Int64[])
end


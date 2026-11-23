using QuadGK
using DataFrames
using Distributions
using SpecialFunctions


"""
logLikelihood(b::Biotope, selVar::Vector{Int64})
computes the logLikelihood conditioned by a subset of covariates given by `selVar`
"""
function logLikelihood(b::Biotope, selVar::Vector{Int64})
  rtol = 1.e-8
  return λ->-sum(log.([quadgk(x->pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]))
end

"""
in-place function to compute f and/or ∇f if there are needed
"""
function getfg!(b::Biotope, selVar::Vector{Int64})
  function fg!(F, G, λ::Vector{Float64})
    rtol=1.e-8  
    Li = [quadgk(x->pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
    if G !== nothing
      ∇Li =[quadgk(x->grad_pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
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
    ∇Li =[quadgk(x->grad_pdf_frailty(b,i,selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] for i=1:b.N]
    G[:] = -sum(∇Li./Li)
  end
  return g!
end

function grad_logLikelihood(b::Biotope, selVar::Vector{Int64})
 rtol=1.e-8
 return λ->-sum([quadgk(x->grad_pdf_frailty(b, i, selVar)(x,λ), 0.0, Inf, rtol=rtol)[1] ./
                 quadgk(x->pdf_frailty(b, i, selVar)(x,λ), 0.0, Inf,rtol=rtol)[1] for i=1:b.N])
end


function hessian_logLikelihood(b::Biotope, selVar::Vector{Int64})
  function proxyHessianLog(λ)
    rtol=1.e-8
    n = size(λ, 1)
    H = zeros(Float64, n,n)
    for i=1:b.N
      ∇Lᵢ = quadgk(x->grad_pdf_frailty(b, i, selVar)(x,λ), 0.0, Inf, rtol=rtol)[1]
      Lᵢ = quadgk(x->pdf_frailty(b, i, selVar)(x,λ), 0.0, Inf,rtol=rtol)[1]
      ∇²Lᵢ = quadgk(x->hessian_pdf_frailty(b, i, selVar)(x,λ), 0.0, Inf, rtol=rtol)[1]
      H += ∇Lᵢ * ∇Lᵢ' ./ Lᵢ^2 - ∇²Lᵢ ./ Lᵢ
    end 
    return H
  end 
end

function logLikelihood(b::Biotope)
  logLikelihood(b, Int64[])
end


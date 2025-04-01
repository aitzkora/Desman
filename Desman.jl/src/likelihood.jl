using QuadGK
using FastGaussQuadrature
using DataFrames
using Distributions
using Integrals

# Without covariances
function distNoCov(df::DataFrame)
    function proxyDist(w, λ) 
        pweibull = x->ccdf(Weibull(exp(λ[2]), w.*exp(λ[1])), x) # complementary distribution fonction 1 - F(x)
      prod(pweibull.(df.durationd) .- pweibull.(df.durationg)) .* pdf(Gamma(exp(λ[3]), exp(-λ[3])),  w)
  end    
  return proxyDist
end

function logLikelihoodNoCov(b::Biotope, methodInt = GaussLegendre)
  return λ->-sum(log.([Integrals.solve(IntegralProblem(distNoCov(b.df[b.idx[i],:]), (0.0, Inf), λ), methodInt())[1] for i=1:b.N]))
end

function distCov(bio::Biotope, i::Int64, selVar::Vector{Int64})
  function proxyDistCov(w, λ)
    idx = collect([size(bio.df,2)-size(bio.Σ,2)+3:size(bio.df,2);]) # Kontuz! 3  pas bon a terme
    γ = λ[4:end]
    z = zeros(bio.N)
    for j in bio.idx[i]
        z[j] = ccdf(Weibull(exp(λ[2]), w.*exp(λ[1]+sum(bio.df[j,idx[selVar]].*γ))), bio.df.durationd[j]) 
             - ccdf(Weibull(exp(λ[2]), w.*exp(λ[1]+sum(bio.df[j,idx[selVar]].*γ))), bio.df.durationg[j])
    end 
    return prod(z).* pdf(Gamma(exp(λ[3]), exp(-λ[3])),  w)
  end
  return proxyDistCov
end

function logLikelihoodCov(b::Biotope, selVar::Vector{Int64}, methodInt = GaussLegendre)
  return λ->-sum(log.([Integrals.solve(IntegralProblem(distCov(b,i,selVar), (0.0, Inf), λ), methodInt())[1] for i=1:b.N]))
end


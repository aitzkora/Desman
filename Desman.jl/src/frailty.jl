"""
pdf_frailty(bio, i, selVar)
computes the density of the `i`-th latrine against the `selVar` collection of covariates
"""
function pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}, ::Type{T}) where {T<:Real}
  function proxyDistCov(w::T, λ::Vector{T})
    return prod_diff(weibull_diff, bio,i, selVar,w,λ) * Γh(λ[3], w)
  end
  return proxyDistCov
end

# default type is Float64
pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}) = pdf_frailty(bio, i, selVar, Float64)

"""
grad_pdf_frailty(bio, i, selVar)
computes the derivativ of the density of the `i`-th latrine against the `selVar` collection of covariates
"""

function grad_pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}, ::Type{T}) where {T<:Real}
  function proxyGrad(w::T, λ::Vector{T})
    p_wei = prod_diff(weibull_diff, bio, i, selVar, w, λ)
    ∂pdf = zeros(size(λ))
    for (jj, j) in enumerate(bio.idxLat[i])
        σ = view(bio.matCov, j, selVar)
        ∂pdf += prod_diff_exclude(weibull_diff, bio, i, selVar, w, λ, j) .* grad_weibull(σ, λ, w, T(bio.durd[j]), T(bio.durg[j]))
    end 
    ∂pdf *= Γh(λ[3], w)
    ∂pdf[3] = p_wei * ∂ᵤΓh(λ[3], w)
    return ∂pdf
  end
  return proxyGrad
end

grad_pdf_frailty(bio::Biotope, i::Int64, selVar::Vector{Int64}) =  grad_pdf_frailty(bio, i, selVar, Float64)

function hessian_pdf(bio::Biotope,  i::Int64, selVar::Vector{Int64}, ::Type{T}) where {T<:Real}

end

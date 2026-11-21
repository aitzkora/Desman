"""
Γh is a optimized version of x->pdf(Gamma(u,1/u),x)
Γh(3.4,5.2) -> 25.503 ns (0 allocations: 0 bytes)
pdf(Gamma(3.4,1. / 3.4), 5.2) ->  235.772 ns (2 allocations: 48 bytes)
"""
function Γh(u::T, w::T) where {T<:Real}
  return w^(u-1)*exp(-w*u)*u^u/gamma(u)
end 

"""
∂ᵤΓh : computes the derivative respect to u of the pdf of Γ(u,1/u)
"""
function ∂ᵤΓh(u::T, w::T) where {T<:Real}
  return Γh(u,w)*(log(w*u)-w-digamma(u)+1)
end

"""
∂²ᵤᵤΓh : computes the second derivative respect to u of the pdf of Γ(u,1/u)
"""
function ∂²ᵤᵤΓh(u::T, w::T) where {T<:Real}
  return ∂ᵤΓh(u::T, w::T)*(log(w*u)-w-digamma(u)+1) +Γh(u,w)*(1/u-trigamma(u))
end

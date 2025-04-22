using LBFGSB
function optimizeEskolZaharra(bio, selVar, m::Int64, λ, print::Int64 = -1, factr::Float64=1e7, lbval::Float64 = 1e-6)
  f = logLikelihood(bio,selVar)
  nVar = length(selVar)
  g! = getg!(bio,selVar)
  μ₀ = [λ; ones(nVar)]
  lb = lbval*[ones(3); ones(nVar)];
  ub = Inf*[ones(3) ; ones(nVar)];
  fout, xout = lbfgsb(f, g!, μ₀, lb=lb, ub=ub, m=m, factr=1e7, pgtol=1e-5, iprint=print, maxfun=15000, maxiter=15000)
  gout = zeros(length(μ₀))
  g!(gout, xout)
  return fout, xout, gout
end 

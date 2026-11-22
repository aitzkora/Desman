module Desman

using Infiltrator
using Revise

# Desmanen komunei buruzko informazioa 😉
export Biotope
include("biotope.jl")

# weibull functions
export weibull_diff, grad_weibull_diff, hessian_weibull_diff
include("weibull.jl")

# product functions
export prod_diff_exclude, prod_diff
include("products.jl")

# Γ related functions
export Γh, ∂ᵤΓh, ∂²ᵤᵤΓh
include("Γ-functions.jl")

# frailty functions
export pdf_frailty, grad_pdf_frailty, hessian_pdf_frailty
include("frailty.jl")

# likelihood related stuff
export logLikelihood, grad_logLikelihood, getfg!, getg!
include("likelihood.jl")

# optimization algos
export optimizeINRIA
include("optimize.jl")

# Miscellanous
export blanksPad
include("misc.jl")

# code called at the beginning of the import
function __init__()
  init_N2QN1()
end

export datapath 
datapath = joinpath(dirname(pathof(Desman)), "..", "data")  # ☣️ non-relocatable

end # module Desman

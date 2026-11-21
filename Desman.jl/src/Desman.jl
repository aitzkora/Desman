module Desman

using DataFrames
                
export Biotope
include("biotope.jl")

# Likelihood related stuff
export logLikelihood, pdf_frailty, logLikelihood2, prod_weibull, weibull_diff, grad_pdf_frailty, grad_logLikelihood, grad_weibull, 
	   Γh, ∂ᵤΓh, ∂²ᵤᵤΓh, prod_weibull_exclude, prod_weibull_view, pdf_frailty_view, grad_logLikelihood_div, getfg!, getg!
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

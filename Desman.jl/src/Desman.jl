module Desman

using DataFrames
export Biotope, logLikelihood, pdf_frailty, logLikelihood2, prod_weibull, weibull_diff,
       grad_pdf_frailty, Γh, ∂ᵤΓh, grad_logLikelihood, grad_weibull, grad_pdf_frailty_check, prod_weibull_exclude,
       prod_weibull_view, pdf_frailty_view, grad_logLikelihood_div, getfg!, getg!, blanksPad, optimizeINRIA
  

include("biotope.jl")
include("likelihood.jl")
include("optimize.jl")
include("misc.jl")

# code called at the beginning of the import
function __init__()
  init_N2QN1()
end

export datapath 
datapath = joinpath(dirname(pathof(Desman)), "..", "data")  # ☣️ non-relocatable

end # module Desman

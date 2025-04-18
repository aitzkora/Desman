module Desman

using DataFrames
export Biotope, logLikelihood, pdf_frailty, logLikelihood2, prod_weibull, weibull_diff,
       grad_pdf_frailty, Γh, ∂ᵤΓh, grad_logLikelihood, grad_weibull, grad_pdf_frailty_check, prod_weibull_exclude  

include("biotope.jl")
include("likelihood.jl")

export datapath 
datapath = joinpath(dirname(pathof(Desman)), "..", "data")  # ☣️ non-relocatable

end # module Desman

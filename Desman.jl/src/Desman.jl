module Desman

using DataFrames
export Biotope, logLikelihood, pdf_frailty, logLikelihood2, prod_weibull

include("biotope.jl")
include("likelihood.jl")

export datapath 
datapath = joinpath(dirname(pathof(Desman)), "..", "data")  # ☣️ non-relocatable

end # module Desman

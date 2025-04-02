module Desman

using DataFrames

export Biotope, logLikelihoodNoCov, distNoCov, distCov, logLikelihoodCov 

include("biotope.jl")
include("likelihood.jl")

export datapath 
datapath = joinpath(dirname(pathof(Desman)), "..", "data")  # ☣️ non-relocatable

end # module Desman

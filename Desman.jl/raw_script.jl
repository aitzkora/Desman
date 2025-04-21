using DataFrames
using CSV
using Desman
using FiniteDiff
using Optim
using Dates
using LineSearches
# read raw data
df1 = CSV.read(joinpath(datapath,"survcalibfeb2024_1.csv"), DataFrame; delim=';', decimal=',');
df2 = CSV.read(joinpath(datapath,"survcalibaug2024_1.csv"), DataFrame; delim=';', decimal=',');
df  = vcat(df1, df2, cols=:union)

# convert string to float and NA to Inf
df.durationg =(x->(x=="NA") ? Inf : parse(Float64,replace(x, "," => "."))).(df[:,:durationg])
finite_idx=findall(df.durationg .< Inf)

## covariance matrix
Σ = CSV.read(joinpath(datapath,"Latrine_cov.csv"), DataFrame; delim=";", decimal=',')

# join data and covariance
df = outerjoin(df, Σ; on=:siteID, matchmissing=:equal, makeunique=true)

## initial λ
λ₀ = sum((df.durationg[finite_idx]+df.durationd[finite_idx])./2.)/length(finite_idx)
λ = [ λ₀, 1, 1]

bio = Biotope(df, Σ)

using BenchmarkTools
using LinearAlgebra
using Test

function optimSelVar(selVar, m)
  fun! = getfg!(bio,selVar)
  nVar= length(selVar)
  lower = 1e-6.* [ones(3) ; ones(nVar)]
  upper = Inf*[ones(3) ; ones(nVar)]
  #linesearch = LineSearches.HagerZhang(linesearchmax = 20)
  inner_optimize =  LBFGS(; m = m)
  μ₀ = [λ; ones(size(selVar,1))]
  sol = optimize(Optim.only_fg!(fun!), lower, upper, μ₀, Fminbox(inner_optimize), 
                 Optim.Options(show_trace = true, show_every = 3, iterations=50, g_tol=1e-3, f_tol=1e-3); inplace=true)
end 

# nCov = size(bio.matCov,2)
# for i=0:2^nCov-1
# selVar = findall(digits(i, base=2, pad=nCov).!=0)
#  selVar = [1,2]
#  nVar = length(selVar)
#
#  println("ₘᵢₙ ℒ(α,β,θ) -> ", sol.minimum)
#end

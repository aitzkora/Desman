using Pkg
Pkg.activate("/home/fux/iturriak/Desman_hub/Desman.jl");
using Desman
using DataFrames
using CSV
using FiniteDiff
using Dates
using Printf


include("optim_bidon.jl")
init_qnb()

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

nCov = size(bio.matCov,2)
sols = [Float64[] for _=1:2^nCov]
selVars = [Int64[] for _=1:2^nCov]
fsols = zeros(2^nCov)
gsols = zeros(2^nCov)
#errors = [ [1,2,3], [3,4] , [1, 2, 3, 4, 5], [2,3,4], [1,2,3,4], [5], [3,5], [4,5], [2,3,4,5]]
errors = []
@printf("| covariates |   f(x)    | |∇f(x)| |\n")

for i=0:2^nCov-1
  selVars[i+1] = findall(digits(i, base=2, pad=nCov).!=0)
  if (!(selVars[i+1] in errors))
   nVar = length(selVars[i+1])
   #fsol, sol, gsol = optimizeEskolZaharra(bio, selVars[i+1], 1+nVar, λ, -1, 1e7, 1e-6) ;
   fsol, sol, gsol = optimizeINRIA(bio, selVars[i+1], 1+nVar, λ, -1, 1e7, 1e-6) ;
   @printf("| %s | %03.5f | %03.5f |\n", blanksPad(selVars[i+1],nCov) ,fsol, norm(gsol))
   sols[i+1]= sol 
   fsols[i+1]= fsol
   gsols[i+1]= norm(gsol)
  else
   sols[i+1] = Float64[]
   fsols[i+1] = Inf
   gsols[i+1] = Inf
  end
end

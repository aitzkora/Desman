#using Pkg
#Pkg.activate("/home/fux/iturriak/Desman_hub/Desman.jl");
using Desman
using DataFrames
using CSV
using FiniteDiff
using Dates
using Printf
using LinearAlgebra

# read raw data
# df = CSV.read(joinpath(datapath,"Simulation-Case-1.csv"), DataFrame; delim=';', decimal='.');
df = DataFrame(CSV.File("Desman.jl/data/Simulation-Case-1.csv"))

# convert string to float and NA to Inf
df.durationg =(x->(x=="NA") ? Inf : parse(Float64,replace(x, "," => "."))).(df[:,:durationg])
finite_idx=findall(df.durationg .< Inf)

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
AIC = zeros(2^nCov)
#errors = [ [1,2,3], [3,4] , [1, 2, 3, 4, 5], [2,3,4], [1,2,3,4], [5], [3,5], [4,5], [2,3,4,5]]
errors = [ [4,5]]
@printf("| covariates |   f(x)       | |∇f(x)| |  #it  |  #sim |\n")

for i=0:2^nCov-1
  selVars[i+1] = findall(digits(i, base=2, pad=nCov).!=0)
  if (!(selVars[i+1] in errors))
   nVar = length(selVars[i+1])
   fsol, sol, gsol, it, nsim = optimizeINRIA(bio, selVars[i+1], λ; lbval= 1e-6, ϵ = 5e-5, print_iter=false)
   @printf("| %s | %03.8f | %03.5f | %05d | %05d |\n", blanksPad(selVars[i+1],nCov) ,fsol, norm(gsol), it, nsim)
   sols[i+1]= sol 
   fsols[i+1]= fsol
   gsols[i+1]= norm(gsol)
   AIC[i+1]= 2*fsol + 2*nVar
  else
   sols[i+1] = Float64[]
   fsols[i+1] = Inf
   gsols[i+1] = Inf
   AIC[i+1] = Inf
  end
end

using DataFrames
using CSV
using Desman
using FiniteDiff
using Optim
using Dates

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
# generates selections of covariates 
nCov = size(bio.Σ,2) - 1
#for i=2^nCov-1:-1:0
#for i=0:2^nCov-1
#  selVar = findall(digits(i, base=2, pad=nCov).!=0)
#  nVar = length(selVar)
#  ## optim problem data
#  f = logLikelihood(bio, selVar)
#  #∇f(x) = FiniteDiff.finite_difference_gradient(f,x)
#  ∇f(x) = grad_logLikelihood(bio, selVar)(x)
#  #
#  ## 0 ⩽ x ⩽ ∞ 
#  lower = 1e-6.* [ones(3) ; ones(nVar)]
#  upper = Inf*[ones(3) ; ones(nVar)]
#  μ₀ = [λ; ones(size(selVar,1))]
#  #
#  ## we choose bfgs with box constrained optimization
#  println("covariates = ", selVar)
#  sol = optimize(f, ∇f, lower, upper, μ₀, Fminbox(LBFGS()), Optim.Options(show_trace = true, show_every = 3, iterations=50, g_tol=1e-3); inplace=false)
#  println("ₘᵢₙ ℒ(α,β,θ) -> ", sol.minimum)
#end

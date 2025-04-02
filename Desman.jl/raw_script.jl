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

# optim problem data
f = logLikelihood(bio)
∇f(x) = FiniteDiff.finite_difference_gradient(f,x)

# 0 ⩽ x ⩽ ∞ 
lower = zeros(3) 
upper = Inf*ones(3) 

# we choose bfgs with box constrained optimization
#sol = optimize(f, ∇f, lower, upper, λ, Fminbox(LBFGS()); inplace=false)

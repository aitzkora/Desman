using Desman
using Test
using CSV
using DataFrames
using Distributions
using SpecialFunctions 
using FiniteDiff
using LinearAlgebra

# internal function to build data in common between tests
function fixture()
  df1 = CSV.read(joinpath(datapath,"survcalibfeb2024_1.csv"), DataFrame; delim=';', decimal=',');
  df2 = CSV.read(joinpath(datapath,"survcalibaug2024_1.csv"), DataFrame; delim=';', decimal=',');
  df  = vcat(df1, df2, cols=:union)
  df.durationg =(x->(x=="NA") ? Inf : parse(Float64,replace(x, "," => "."))).(df[:,:durationg])
  finite_idx=findall(df.durationg .< Inf)
  Σ = CSV.read(joinpath(datapath,"Latrine_cov.csv"), DataFrame; delim=";", decimal=',')
  df = outerjoin(df, Σ; on=:siteID, matchmissing=:equal, makeunique=true)
  # initial λ
  λ₀ = sum((df.durationg[finite_idx]+df.durationd[finite_idx])./2.)/length(finite_idx)
  λ = [ λ₀, 1, 1]
  bio = Biotope(df, Σ)
  return bio, λ
end


function exec_tests(args)
  nArgs = length(args)
  test_files = ["likelihood-λ",
                "Γ-functions",
                "products",
                "gradients" ]

  if nArgs != 0  
     if ARGS[1] in test_files
       test_files = [ARGS[1]]
     else
       @error "test suite must ∈ $(test_files)"
       return nothing # manage correctly aborting execution 
     end
  end
  for tst in test_files
     @testset "$tst" begin
       include("$(tst).jl")
     end    
  end
end

# do the exec trough a function
exec_tests(ARGS)

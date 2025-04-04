using Desman
using Test
using CSV
using DataFrames
using Distributions
using SpecialFunctions 


@testset "likelihood-λ" begin
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
    
   @test logLikelihood(bio)(λ) ≈  936.14681223 atol = 1e-5
end

@testset "Gamma functions" begin
   @test pdf(Gamma(3.33, 1. /3.33), 54.) ≈ Γh(3.33, 54.) atol=1e-8
   @test ∂ᵤΓh(3.2,5.1) ≈ (-0.00011584186052754) atol=1e-15
end

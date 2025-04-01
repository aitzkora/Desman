using Desman
using Test
using CSV
using DataFrames
using QuadGK


@testset "likelihood" begin
    dataroot =     # read raw data
    df1 = CSV.read(joinpath(datapath, "survcalibfeb2024_1.csv"), DataFrame; delim=';', decimal=',');
    df2 = CSV.read(joinpath(datapath, "survcalibaug2024_1.csv"), DataFrame; delim=';', decimal=',');
    df  = vcat(df1, df2, cols=:union)
    
    # convert string to float and NA to Inf
    df.durationg =(x->(x=="NA") ? Inf : parse(Float64,replace(x, "," => "."))).(df[:,:durationg])
    finite_idx=findall(df.durationg .< Inf)
    #
    ## covariance matrix
    Σ = CSV.read(joinpath(datapath, "Latrine_cov.csv"), DataFrame)
    #
    ## initial λ
    λ₀ = sum((df.durationg[finite_idx]+df.durationd[finite_idx])./2.)/length(finite_idx)
    λ = [ λ₀, 1, 1]
    #
    bio = Biotope(df, Σ, λ)
    
    l1 = logVraisemblanceNoCov(bio, λ, quadgk)

    @info  l1

   @test l1 ≈  -936.2474 atol = 1e-5

end


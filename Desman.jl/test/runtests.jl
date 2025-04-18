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

@testset "likelihood-λ" begin
    # compare logLikelihood to R result
    bio, λ = fixture()
    @test logLikelihood(bio, Int64[])(exp.(λ)) ≈  936.14681223 atol = 1e-5

    # gradient check up to 2%
    selVar = [1,2]
    μ₀ = [λ; ones(length(selVar))]
    g_ana =  grad_logLikelihood(bio, selVar)(μ₀) 
    g_diff = FiniteDiff.finite_difference_gradient(logLikelihood(bio, selVar), μ₀)
    @test norm(g_ana-g_diff)/norm(g_diff) < 0.02
end

@testset "Gamma functions" begin
   @test pdf(Gamma(3.33, 1. /3.33), 54.) ≈ Γh(3.33, 54.) atol=1e-8
   @test ∂ᵤΓh(3.2,5.1) ≈ (-0.00011584186052754) atol=1e-15
end

@testset "Products" begin

   bio, _ = fixture()
   selVar = [2,3,4]
   nvar = length(selVar)
   i = 15
   j = bio.idxLat[i][end]
   λ = rand(3+nvar)
   w = (rand(1)*20 .+ 15)[1]
   idx = collect([bio.covStart:bio.covEnd;])
   σ = Matrix{Float64}(bio.df[j:j,idx[selVar]])
   @test prod_weibull_exclude(bio, i, selVar, w, λ, j) ≈ prod_weibull(bio, i, selVar, w, λ) ./ weibull_diff(σ, λ, w, bio.df.durationd[j], bio.df.durationg[j]) atol=1e-10
end


@testset "Gradients" begin
   # weibull
   nvar = 3
   σ = rand(1,nvar)
   w = (rand(1)*10 .+ 10)[1]
   λ = rand(3+nvar)
   xd = rand(1)[1]
   xg = rand(1)[1]
   @test FiniteDiff.finite_difference_gradient(x->weibull_diff(σ,x, w, xg, xd), λ) ≈ grad_weibull(σ, λ, w, xg, xd) atol=1e-10


   # take into account some special values for xd or xg
   @test isfinite(norm(grad_weibull(σ, λ, w, 0.0, xg)))
   @test isfinite(norm(grad_weibull(σ, λ, w, xd, 0.0)))
   @test isfinite(norm(grad_weibull(σ, λ, w, xd, Inf)))

   # pdf_frailty
  
   bio, λ = fixture()
   i = 12
   selVar = [3,4];
   γ = rand(2)
   g = pdf_frailty(bio, i, selVar)
   w = w
   μ = [λ;γ]
   @test FiniteDiff.finite_difference_gradient(x->g(w,x), μ) ≈ grad_pdf_frailty(bio,i, selVar)(w, μ) atol=1e-10

   # pdf_frailty_check
   bio, λ = fixture()
   i = 12
   selVar = [3,4];
   γ = rand(2)
   g = pdf_frailty(bio, i, selVar)
   w = (rand(1)*10 .+ 10)[1]
   w = w
   μ = [λ;γ]
   @test grad_pdf_frailty(bio,i, selVar)(w, μ) ≈ grad_pdf_frailty(bio,i, selVar)(w, μ) atol=1e-10

end  

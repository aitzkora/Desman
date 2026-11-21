# compare logLikelihood to R result
bio, λ = fixture()
@test logLikelihood(bio, Int64[])(exp.(λ)) ≈  936.14681223 atol = 1e-5

# gradient check up to 2%
selVar = [1,2]
μ₀ = [λ; ones(length(selVar))]
g_ana =  grad_logLikelihood(bio, selVar)(μ₀) 
g_diff = FiniteDiff.finite_difference_gradient(logLikelihood(bio, selVar), μ₀)
@test norm(g_ana-g_diff)/norm(g_diff) < 0.02

#test inplace fun and gradient computation
fg! = getfg!(bio, selVar)
lf = 0.
lg = zeros(length(μ₀))
lf = fg!(lf, lg, μ₀)
of = logLikelihood(bio, selVar)(μ₀)

@test lf ≈ of atol = 1e-10
@test norm(lg-g_ana) / norm(g_ana)  < 0.02

#test with just gradient in place computation
g_in_place! = getg!(bio, selVar)
lg = zeros(length(μ₀))
g_in_place!(lg,μ₀)
@test norm(lg-g_ana) / norm(g_ana)  < 0.02


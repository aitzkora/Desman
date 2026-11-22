using ForwardDiff
# weibull

nvar = 5
σ = view(collect(rand(nvar,1)'),1,collect([1:nvar;]))
w = (rand(1)*10 .+ 10)[1]
λ = rand(3+nvar)
xg = rand(1)[1]
xd = xg + rand(1)[1]

@test FiniteDiff.finite_difference_gradient(x->weibull_diff(σ,x, w, xg, xd), λ) ≈ grad_weibull_diff(σ, λ, w, xg, xd) atol=1e-10

@test ForwardDiff.hessian(x->weibull_diff(σ,x, w, xg, xd),λ) ≈ hessian_weibull_diff(σ, λ, w, xg, xd) atol=1e-10

# take into account some special values for xd or xg
@test isfinite(norm(grad_weibull_diff(σ, λ, w, 0.0, xg)))
@test isfinite(norm(grad_weibull_diff(σ, λ, w, xd, 0.0)))
@test isfinite(norm(grad_weibull_diff(σ, λ, w, xd, Inf)))

# pdf_frailty

bio, λ = fixture()
i = 12
selVar = [3,4];
γ = rand(2)
g = pdf_frailty(bio, i, selVar)
w = w
μ = [λ;γ]
@test FiniteDiff.finite_difference_gradient(x->g(w,x), μ) ≈ grad_pdf_frailty(bio,i, selVar)(w, μ) atol=1e-10

# pdf_frailty
bio, λ = fixture()
i = 12
selVar = [3,4];
γ = rand(2)
g = pdf_frailty(bio, i, selVar)
w = (rand(1)*10 .+ 10)[1]
w = w
μ = [λ;γ]
@test grad_pdf_frailty(bio,i, selVar)(w, μ) ≈ grad_pdf_frailty(bio,i, selVar)(w, μ) atol=1e-10


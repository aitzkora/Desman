@test pdf(Gamma(3.33, 1. /3.33), 54.) ≈ Γh(3.33, 54.) atol=1e-6
@test ∂ᵤΓh(3.2,5.1) ≈ (-0.00011584186052754) atol=1e-15

let 
  g(x)=∂ᵤΓh(x[1], 0.4)
  @test ∂²ᵤᵤΓh(0.3646, 0.4) ≈ FiniteDiff.finite_difference_gradient(g, [0.3646])[1] atol=1e-8
end

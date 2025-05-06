using N2QN1_jll
using LinearAlgebra
using Libdl:dlopen, dlsym
using Printf


function init_N2QN1()
   global n2qn1_ptr
   qnb_lib = dlopen(libn2qn1r)
   n2qn1_ptr = dlsym(qnb_lib, :n2qn1_)
end

function call_N2QN1R!(x::Vector{Float64}, f::Float64, g::Vector{Float64},
                   dxmin::Vector{Float64}, df1::Float64, epsabs::Vector{Float64},
                   imp::Int64, io::Int64, mode::Vector{Cint}, iter::Vector{Cint}, nsim::Vector{Cint},
                   binf::Vector{Float64}, bsup::Vector{Float64}, iz::Vector{Cint}, rz::Vector{Float64},
                   reverse::Bool)
    @assert length(x) == length(g) == length(bsup) == length(binf)
    n = length(x)
    ccall(n2qn1_ptr, Cvoid, (Ref{Cint},    Ref{Float64}, Ref{Float64}, Ref{Float64}, # f, f, g
                             Ptr{Float64}, Ref{Float64}, Ptr{Float64}, # dxmin, df1, epsabs
                             Ref{Cint},    Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, # imp, io, mode, iter, nsim   
                             Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Ptr{Float64}, Ref{Cint}), # binf, bsup, iz, rz, reverse
                             n, x, f, g,
                             dxmin, df1, epsabs,
                             imp, io, mode, iter, nsim,
                             binf, bsup, iz, rz, reverse) 
end


function bfgsb(f,g!, x0::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64};
               dxmin::Vector{Float64}=1e-12*ones(size(x0,1)), ϵ::Float64=1e-4, max_iter::Int64=300, print_iter::Bool=false)
  # allocate 
  xk = copy(x0)
  gk = zeros(length(x0))
  # init parameters
  df1 = f(x0)
  mode = ones(Cint, 1)
  imp = 0
  it = 0
  iter = Cint[max_iter]
  nsim = Cint[iter[1] + iter[1]÷2] # nsim ≈ max_iter * 1.5
  iz = zeros(Cint, 1024) # why 1024
  rz = zeros(Float64, 1024)
  if (print_iter) 
    imp = 1
    io = 6 
  else
    imp = 0
    io = 22 # !FIXME remove that file at end
  end
  fk=df1
  g!(gk,x0)
  epsabs = ϵ*ones(1)*norm(gk)
  call_N2QN1R!(xk, fk, gk, dxmin, df1, epsabs, imp, io, mode, iter, nsim, lb, ub, iz, rz, true)
  while (mode[1]>7 && it < max_iter)
    it += 1
    fk = f(xk)
    g!(gk, xk)
    if (print_iter)
       @printf("| %05d | %.3e | %.3e |\n", it, fk, norm(gk))
    end
    call_N2QN1R!(xk, fk, gk, dxmin, df1, epsabs, imp,
              io, mode, iter, nsim, lb, ub, iz, rz, true)
  end
  return fk, xk, iter[1], nsim[1]
end

# proxy function calling the bfgs with variables of the project
function optimizeINRIA(bio, selVar, λ; lbval::Float64 = 1e-6, ϵ::Float64=5e-5)
    f = logLikelihood(bio,selVar)
    nVar = length(selVar)
    g! = getg!(bio,selVar)
    μ₀ = [λ; ones(nVar)]
    lb = lbval*[ones(3); ones(nVar)];
    ub = Inf*[ones(3) ; ones(nVar)];
    fout, xout, it, sim = bfgsb(f,g!,μ₀, lb, ub; print_iter=false, ϵ = ϵ) 
    gout = zeros(length(μ₀))
    g!(gout, xout)
    return fout, xout, gout, it, sim
end 

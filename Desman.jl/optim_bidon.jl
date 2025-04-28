using LinearAlgebra

using Libdl:dlopen, dlsym

using Printf

function init_qnb()
   global n2qn1_ptr

   qnb_lib = dlopen("/home/fux/iturriak/Desman_hub/Desman.jl/qnb.so")

   n2qn1_ptr = dlsym(qnb_lib, :n2qn1_)

end

function call_qnb!(x::Vector{Float64},
                   f::Float64,
                   g::Vector{Float64},
                   dxmin::Vector{Float64},
                   df1::Float64,
                   epsabs::Vector{Float64},
                   imp::Int64,
                   io::Int64,
                   mode::Vector{Cint},
                   iter::Int64,
                   nsim::Int64,
                   binf::Vector{Float64},
                   bsup::Vector{Float64},
                   iz::Vector{Cint},
                   rz::Vector{Float64},
                   reverse::Bool)
    @assert length(x) == length(g) == length(bsup) == length(binf)
    n = length(x)
    ccall(n2qn1_ptr, Cvoid, (Ref{Cint},    # n
                             Ref{Float64}, # x(1:n)
                             Ref{Float64}, # f
                             Ref{Float64}, # g
                             Ptr{Float64}, # dxmin
                             Ref{Float64}, # df1
                             Ptr{Float64}, # epsabs
                             Ref{Cint},    # imp
                             Ref{Cint},    # io
                             Ptr{Cint},    # mode
                             Ref{Cint},    # iter
                             Ref{Cint},    # nsim
                             Ptr{Float64}, # binf
                             Ptr{Float64}, # bsup
                             Ptr{Cint},    # iz
                             Ptr{Float64}, # rz
                             Ref{Cint}), # reverse
                             n,
                             x, 
                             f,
                             g,
                             dxmin,
                             df1,
                             epsabs,
                             imp,
                             io,
                             mode,
                             iter,
                             nsim,
                             binf,
                             bsup,
                             iz,
                             rz,
                             reverse) 
end


function bfgsb(f,g!, x0::Vector{Float64}, 
                     lb::Vector{Float64}, 
                     ub::Vector{Float64};
                        dxmin::Vector{Float64}=1e-12*ones(size(x0,1)),
                        print_iter::Bool=false,
                        max_iter::Int64=300)
  # allocate 
  xk = copy(x0)
  gk = zeros(length(x0))
  # init parameters
  df1 = f(x0)
  mode = ones(Cint, 1)
  imp = 0
  it = 0
  iter =  100
  nsim=3*iter
  iz = zeros(Cint, 1024) # why 1024
  rz = zeros(Float64, 1024)
  io = 06  # descripteur de fichier du journal 06-> std output
  fk=df1
  g!(gk,x0)
  epsabs = 5e-5*ones(1)*norm(gk)
  call_qnb!(xk, fk, gk, dxmin, df1, epsabs, imp, io, mode, iter, nsim, lb, ub, iz, rz, true)
  while (mode[1]>7 && it < max_iter)
    it += 1
    fk = f(xk)
    g!(gk, xk)
    if (print_iter)
       @printf("| %05d | %.3e | %.3e |\n", it, fk, norm(gk))
    end
    nsim=3*iter
    call_qnb!(xk, fk, gk, dxmin, df1, epsabs, imp,
              io, mode, iter, nsim, lb, ub, iz, rz, true)
  end
  return fk, xk
end


function optimizeINRIA(bio, selVar, m::Int64, λ, print::Int64 = -1, factr::Float64=1e7, lbval::Float64 = 1e-6)

    f = logLikelihood(bio,selVar)
    nVar = length(selVar)
    g! = getg!(bio,selVar)
    μ₀ = [λ; ones(nVar)]
    lb = lbval*[ones(3); ones(nVar)];
    ub = Inf*[ones(3) ; ones(nVar)];
    fout, xout = bfgsb(f,g!,μ₀, lb, ub; print_iter=false)
    gout = zeros(length(μ₀))
    g!(gout, xout)
    return fout, xout, gout
end 



using LBFGSB
function optimizeEskolZaharra(bio, selVar, m::Int64, λ, print::Int64 = -1, factr::Float64=1e7, lbval::Float64 = 1e-6)
  f = logLikelihood(bio,selVar)
  nVar = length(selVar)
  g! = getg!(bio,selVar)
  μ₀ = [λ; ones(nVar)]
  lb = lbval*[ones(3); ones(nVar)];
  ub = Inf*[ones(3) ; ones(nVar)];
  fout, xout = lbfgsb(f, g!, μ₀, lb=lb, ub=ub, m=m, factr=1e7, pgtol=1e-5, iprint=print, maxfun=15000, maxiter=15000)
  gout = zeros(length(μ₀))
  g!(gout, xout)
  return fout, xout, gout
end



using N2QN1_jll
using LinearAlgebra
using Libdl:dlopen, dlsym
function init_n2qn1()
   global n2qn1_ptr

   n2qn1_lib = dlopen(libn2qn1)
   n2qn1_ptr = dlsym(n2qn1_lib, :n2qn1_)

end

function call_n2qn1!(simu,
                     x::Vector{Float64},
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
                     bsup::Vector{Float64})
    @assert length(x) == length(g) == length(bsup) == length(binf)
    n = length(x)
    iz = zeros(Cint, 2*n+1) # why 1024
    rz = zeros(Float64, n*(n+9)÷2)
    ccall(n2qn1_ptr, Cvoid, (Ptr{Cvoid}, Ref{Cint},  Ref{Float64}, Ref{Float64}, Ref{Float64},  # simu, n , x, f , g
                             Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Cint}, Ref{Cint},   # dxmin, df1, epsabs, imp, io
                             Ptr{Cint},    Ref{Cint},    Ref{Cint},                            # mode, iter, nsim
                             Ptr{Float64}, Ptr{Float64}, Ptr{Cint},  Ptr{Float64},             # binf, bsup, iz, rz   
                             Ptr{Cint}, Ptr{Float32}, Ptr{Float64}),                           # izs, rzs, dzs 
                             simu, n, x, f, g,
                             dxmin, df1, epsabs, imp, io,
                             mode, iter, nsim,
                             binf, bsup, iz, rz, 
                             C_NULL, C_NULL, C_NULL) 
end


function getSimul(n, obj, grad!)
  function proxy!(indic::Ptr{Cint},n::Cint,x::Ptr{Float64},f::Ptr{Float64},g::Ptr{Float64},izs::Ptr{Cint},rzs::Ptr{Float32},dzs::Ptr{Float64})
    if (unsafe_load(indic, 1)==4)
       xv = unsafe_wrap(Array, x, n)
       gv = unsafe_wrap(Array, g, n)
       grad!(gv,xv)
       fv = obj(xv)
       unsafe_store!(f,fv)
    end     
    nothing
  end
  return proxy!
end 

function bfgsb_direct(f,g!, 
                      x0::Vector{Float64}, 
                      lb::Vector{Float64}, 
                      ub::Vector{Float64};
                      dxmin::Vector{Float64}=1e-12*ones(size(x0,1)),
                      print_iter::Bool=false,
                      max_iter::Int64=300)
  # allocate 
  n = length(x0)
  xk = copy(x0)
  gk = zeros(length(x0))
  # init parameters
  df1 = f(x0)
  mode = ones(Cint, 1)
  imp = 0 
  iter = 300
  nsim=3*iter
  io = 22  
  fk=df1
  g!(gk,x0)
  epsabs = 5e-5*ones(1)*norm(gk)
  simu = @cfunction($(getSimul(n, f, g!)), Cvoid, (Ptr{Cint}, Ref{Cint}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Ptr{Float32}, Ptr{Float64}))
  call_n2qn1!(simu, xk, fk, gk, dxmin, df1, epsabs, imp, io, mode, iter, nsim, lb, ub)
  return fk, xk, iter, nsim
end

function optimizeINRIA_direct(bio, selVar, m::Int64, λ, lbval::Float64 = 1e-6)
    f = logLikelihood(bio,selVar)
    nVar = length(selVar)
    g! = getg!(bio,selVar)
    μ₀ = [λ; ones(nVar)]
    lb = lbval*[ones(3); ones(nVar)];
    ub = Inf*[ones(3) ; ones(nVar)];
    fout, xout, it , nsim= bfgsb_direct(f,g!,μ₀, lb, ub; print_iter=false)
    gout = zeros(length(μ₀))
    g!(gout, xout)
    return fout, xout, gout, it, nsim
end 


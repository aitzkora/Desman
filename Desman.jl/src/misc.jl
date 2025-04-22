function blanksPad(s::String, n::Int64)
  @assert length(s) <= n
  reduce(*,[string(c)*" " for c in s])*(" "^(2*(n-length(s))))
end

function blanksPad(l::Vector{Int64}, n::Int64)
    blanksPad(reduce(*,string.(l)), n)
end

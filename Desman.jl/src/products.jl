function prod_diff(base_diff, bio::Biotope, i::Int64, selVar::Vector{Int64}, w::T, λ::Vector{T}) where {T<:Real}
  z = one(T)
  for j in bio.idxLat[i]
      z *= base_diff(view(bio.matCov,j,selVar), λ, w, T(bio.durd[j]), T(bio.durg[j]))
  end 
  return z
end

function prod_diff_exclude(base_diff, bio::Biotope, i::Int64, selVar::Vector{Int64}, w::T, λ::Vector{T}, k::Int64) where {T<:Real}
  z = one(T)
  for j in bio.idxLat[i]
    if (j != k)
        z *= base_diff(view(bio.matCov,j,selVar), λ, w, T(bio.durd[j]), T(bio.durg[j]))
      end
  end 
  return z
end

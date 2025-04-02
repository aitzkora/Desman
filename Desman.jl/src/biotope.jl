"""
structure for storing a Desman biotope place
"""
struct Biotope 
  df::DataFrame 
  Σ::DataFrame # covariance dataframe
  latrine::Vector{String}
  N::Int64 # number of distinct latrines
  freq::Vector{Int64}   # to store frequencies
  idx::Vector{Vector{Int64}} #
  # constructor
  function Biotope(df::DataFrame, cov::DataFrame)
    lat = sort(unique(df.siteID))
    N = length(lat)
    freq = combine(groupby(df, :siteID), nrow => :Freq).Freq
    idx = [zeros(Int64,0) for _ in 1:N]
    for i=1:N
      idx[i] = findall(df.siteID .== lat[i])
    end
    new(df, cov, lat, N, freq, idx)
  end 
end

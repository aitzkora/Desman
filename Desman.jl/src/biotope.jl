"""
structure for storing a Desman biotope place
"""
struct Biotope 
  df::DataFrame 
  Σ::DataFrame # covariance dataframe
  latrine::Vector{String}
  N::Int64 # number of distinct latrines
  freq::Vector{Int64}   # to store frequencies
  idxLat::Vector{Vector{Int64}} # sites indexes corresponding to a Latrin
  covStart::Int64
  covEnd::Int64
  # constructor
  function Biotope(df::DataFrame, cov::DataFrame)
    lat = sort(unique(df.siteID))
    N = length(lat)
    freq = combine(groupby(df, :siteID), nrow => :Freq).Freq
    idxLat = [zeros(Int64,0) for _ in 1:N]
    for i=1:N
      idxLat[i] = findall(df.siteID .== lat[i])
    end
    covStart = size(df,2)-size(cov,2)+2
    covEnd = size(df,2)
    new(df, cov, lat, N, freq, idxLat, covStart, covEnd)
  end 
end

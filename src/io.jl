function make_results_df(times::AbstractVector{T}, rvs_clean::AbstractVector{T}, σ_rvs_clean::AbstractVector{T},
              rvs_dirty::AbstractVector{T} = Vector{T}(undef,0), σ_rvs_dirty::AbstractVector{T} = Vector{T}(undef,0),
              scores::AbstractArray{T,2} = Array{T,2}(undef,0,0), σ_scores::AbstractArray{T,2} = Array{T,2}(undef,0,0),
              #unit_scores::AbstractVector{String} = fill("nondim",size(scores,2))
              ) where { T<:Real }
  @assert length(times) == length(rvs_clean) == length(σ_rvs_clean)
  @assert length(rvs_dirty) == length(σ_rvs_dirty)
  output_df = DataFrame("Time [MJD]"=>times, "RV_C [m/s]"=>rvs_clean, "e_RV_C [m/s]"=>σ_rvs_clean,
                    "RV_D [m/s]"=>rvs_dirty, "e_RV_D [m/s]"=>σ_rvs_dirty )
  if length(rvs_dirty) == length(times)
    output_df[!,"RV_D [m/s]"] = rvs_dirty
    output_df[!,"e_RV_D [m/s]"] = σ_rvs_dirty
  end
  if size(scores,1) >= 1
    @assert length(times) == size(scores,1)
    if size(σ_scores,2) == 0
      σ_scores = fill(NaN,length(times),size(scores,2)  )
    end
    @assert length(times) == size(σ_scores,1)
    for i in 1:size(scores,2)
      output_df[!,"Score_" * string(i) * " [nondim]"] = vec(scores[:,i])
      output_df[!,"e_Score_" * string(i) * " [nondim]"] = vec(σ_scores[:,i])
    end
  end
  return output_df
end

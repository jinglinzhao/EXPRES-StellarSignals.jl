 verbose = true
 make_plots = false
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
if verbose   println("# Loading other packages")    end
 using CSV, DataFrames, Query, StatsBase, Statistics, Dates
 if verbose   println("# Loading unpackaged code from this repo")    end
 include(joinpath("..","..","src","io.jl"))

 function calc_epoch_offset_vector(rvs_ccf, σ_rvs, tidx_epoch4, tidx_epoch5)
     vcat(rvs_ccf[tidx_epoch4].-mean(rvs_ccf[tidx_epoch4], weights(1.0 ./σ_rvs[tidx_epoch4].^2)), (rvs_ccf[tidx_epoch5].-mean(rvs_ccf[tidx_epoch5],weights(1.0 ./σ_rvs[tidx_epoch5].^2))))
 end

 using Scalpels
 function compute_scalpels_output(mnb::Integer; sort_by_responce::Bool=true)
   (rv_scalpels, scores, ) = calc_clean_rvs_scores_basis_scalpels(rvs_ccf, ccfs, num_basis=mnb, σ_rvs=σ_rvs, sort_by_responce=sort_by_responce)
   rv_scalpels_offsets = (rv_scalpels .- vcat(fill(mean(rv_scalpels[tidx_epoch4], weights(1.0 ./σ_rvs[tidx_epoch4].^2)), sum(tidx_epoch4)), fill(mean(rv_scalpels[tidx_epoch5],weights(1.0 ./σ_rvs[tidx_epoch5].^2)), sum(tidx_epoch5)) ) )

   println("max_num_basis = ", mnb, "  std RV = ", std(rv_scalpels_offsets) )

   output_df = DataFrame("Time [MJD]"=>order_list_timeseries.times, "RV_C [m/s]"=>rv_scalpels_offsets, "e_RV_C [m/s]"=>σ_rvs,
                       "RV_D [m/s]"=>rvs_ccf.-rv_scalpels_offsets, "e_RV_D [m/s]"=>σ_rvs )
   for i in 1:mnb
      output_df[!,"Score_" * string(i) * " [nondim]"] = vec(scores[:,i])
      output_df[!,"e_Score_" * string(i) * " [nondim]"] = fill(NaN,length(rv_scalpels_offsets))
   end
   return output_df
 end


all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_expres_data_" * string(starid) * ".jl"))
 GC.gc()
 ref_obs_idx = RvSpectMLBase.choose_obs_idx_for_init_guess(df_files_use,get_inst(all_spectra))
 tidx_epoch4 = df_files_use.expres_epoch.==4
 tidx_epoch5 = df_files_use.expres_epoch.==5

# Something buggy if include order 1.  Order 86 essentially ignored due to tellurics. First and last orders have NaN issues
lfc_orders = 43:72
max_orders = 12:83
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=lfc_orders, recalc=true )
  #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )
  order_list_timeseries_ref_obs_only = extract_chunklist_timeseries_with_subset_obs(order_list_timeseries, ref_obs_idx)

linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",espresso_mask_filename)
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=false)
 ((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries_ref_obs_only, line_list_espresso, pipeline_plan,  mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, v_step=100, recalc=true, calc_ccf_var=true)
 line_width_50 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.5)
 line_sigma = line_width_50/sqrt(2*log(2))
 line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)
 σ_rvs_formal_espresso = map(obsid->calc_formal_rv_precission(all_spectra[obsid],RvSpectMLBase.make_chunk_list_from_loc_df(all_spectra[obsid],get_inst(all_spectra), line_list_espresso),smooth_factor=2.0),1:length(all_spectra) )
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf_ref = RvSpectML.calc_rvs_from_ccf_total( reshape(ccfs,size(ccfs,1),1), reshape(ccf_vars,size(ccf_vars,1),1), pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, bin_nightly=false, recalc=true, alg_fit_rv=alg_fit_rv)
 mean(σ_rvs_formal_espresso), median(σ_rvs_formal_espresso), extrema(σ_rvs_formal_espresso), std(rvs_ccf_ref)

((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries_ref_obs_only, line_list_espresso, pipeline_plan,  mask_scale_factor=4.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, recalc=true, calc_ccf_var=true)
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf_ref = RvSpectML.calc_rvs_from_ccf_total( reshape(ccfs,size(ccfs,1),1), reshape(ccf_vars,size(ccf_vars,1),1), pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, bin_nightly=false, recalc=true, alg_fit_rv=alg_fit_rv)
 σ_rv_fit_ref = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 σ_rvs_formal_espresso
 ccf_var_scale_factor = σ_rvs_formal_espresso[ref_obs_idx]/σ_rv_fit_ref[1]

((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan,  mask_scale_factor=4.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, ccf_var_scale=ccf_var_scale_factor, calc_ccf_var=true, recalc=true)
 # ccf_vars_espresso .*= ccf_var_scale_factor  # If wanted to scale ccf_var manually, rather than passing ccf_var_scale=ccf_var_scale_factor above
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, ccf_vars_espresso, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true,  bin_nightly=false, alg_fit_rv=alg_fit_rv)
 σ_rvs_ccf_espresso = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 rvs_ccf_espresso, σ_rvs_ccf_espresso
 σ_rvs_ccf_espresso./σ_rvs_formal_espresso
 println("CCF ESPRESSO, LFC: ")
 rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_espresso.-mean(rvs_ccf_espresso)))
 rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_espresso.-mean(rvs_ccf_espresso))
 println("# RMS of RVs: ", std(rvs_ccf_espresso), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night, " σ_rvs_formal = ", median(σ_rvs_ccf_espresso))

#=
# If you wanted to compute the full covariance matrix for the CCFs.  VERY slow.  And not that different.
((ccfs_espresso2, ccf_covars_espresso2), v_grid2) = ccf_total(order_list_timeseries_ref_obs_only, line_list_espresso, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, calc_ccf_covar=true, recalc=true)
alg_fit_rv2 =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
rvs_ccf_espresso2 = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso2, ccf_covars_espresso2, pipeline_plan, v_grid=v_grid2, times = order_list_timeseries.times, recalc=true, bin_nightly=false, alg_fit_rv=alg_fit_rv2)
rvs_ccf_espresso2, σ_rvs_ccf_espresso2
=#


num_basis_scalpels = 2
  println("# ESPRESSO mask, LFC, Scalpels, num_basis_vectors = ", num_basis_scalpels)
  ccfs = ccfs_espresso; ccf_vars = ccf_vars_espresso; rvs_ccf = rvs_ccf_espresso; σ_rvs = σ_rvs_ccf_espresso; println("# Setting ccfs to ESPRESSO mask ", num_basis_scalpels, " basis vectors");  @time df_scalpels0 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
  CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), df_scalpels0)
  rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels0[!,Symbol("RV_C [m/s]")].-mean(df_scalpels0[!,Symbol("RV_C [m/s]")])))
  rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels0[!,Symbol("RV_C [m/s]")].-mean(df_scalpels0[!,Symbol("RV_C [m/s]")]))
  println("# RMS of RVs: ", std(df_scalpels0[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

#= =#

chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_telluric_free(all_spectra, line_list_espresso, min_pixels_in_chunk=5, verbose=false)
 dont_need_to!(pipeline_plan,:extract_orders)
 rv_offset_epoch4_minus_epoch5 = calc_epoch_offset_vector(rvs_ccf_espresso, σ_rvs_formal_espresso, tidx_epoch4, tidx_epoch5 )
 set_rv_est!(chunk_list_timeseries, rvs_ccf_espresso .+ tidx_epoch5.*rv_offset_epoch4_minus_epoch5 )
 chunk_grids = RvSpectML.make_grids_for_chunklist_timeseries(chunk_list_timeseries)
 λmatrix = reduce(vcat,chunk_grids)
 idx_use_pca = trues(length(all_spectra));
 if starid == 101501 idx_use_pca[2] = false  end
 size(λmatrix)

 #using Plots
 #plt = scatter(chunk_list_timeseries.times,rvs_ccf_espresso.-mean(rvs_ccf_espresso),yerr=σ_rvs_ccf_espresso, label=:none, xlabel="Time (d)", ylabel="RV (m/s)", title="RVs CCF, ESPRESSO, LFC")
 #plt = scatter(rvs_ccf_espresso.-mean(rvs_ccf_espresso), label="RVs CCF, ESPRESSO LFC")
 #savefig("rvs_" * string(starid) * "_espresso_lfc.png")

 num_basis_dcpca = 2
 println("# ESPRESSO mask, LFC, DCPCA, num_basis_vectors = ", num_basis_dcpca)
 #for num_basis_dcpca in 1:6
   dcpca_result_espresso = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_espresso, num_basis_vectors_pca=num_basis_dcpca, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=40)
   #print("DCPCA ", num_basis_dcpca, ": " ) #std(rvs) = ", std(dcpca_result_espresso.rvs_clean))
   rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=dcpca_result_espresso.rvs_clean.-mean(dcpca_result_espresso.rvs_clean)))
   rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=dcpca_result_espresso.rvs_clean.-mean(dcpca_result_espresso.rvs_clean))
   println("# RMS of RVs: ", std(dcpca_result_espresso.rvs_clean), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

   #scatter!(plt,dcpca_result_espresso.rvs_dirty, label=string(num_basis_dcpca))
 #end
 #display(plt)
 σ_rvs_clean_espresso = σ_rvs_formal_espresso # sqrt.(σ_rvs_formal_espresso.^2 .+ dcpca_result_espresso.σ_rvs_dirty.^2 )
 df_espresso_dcpca0 = make_results_df(order_list_timeseries.times,dcpca_result_espresso.rvs_clean,σ_rvs_clean_espresso,dcpca_result_espresso.rvs_dirty,NaN.+dcpca_result_espresso.σ_rvs_dirty,dcpca_result_espresso.scores,NaN.+dcpca_result_espresso.σ_scores)
 CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_dcpca0_results.csv"), df_espresso_dcpca0)



#dcpca_result_espresso = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_espresso, num_basis_vectors_pca=6, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=10)
#dcpca_result_espresso.σ_scores

#=
using Plots
scatter(chunk_list_timeseries.times)
plt = scatter(dcpca_result_espresso.rvs_dirty)
plt = scatter(dcpca_result_espresso.scores[:,1])
scatter!(plt,dcpca_result_espresso.scores[:,1])
scatter!(plt,dcpca_result_espresso.scores[:,2])
scatter!(plt,dcpca_result_espresso.scores[:,3])
scatter!(plt,dcpca_result_espresso.scores[:,4])
scatter!(plt,dcpca_result_espresso.scores[:,5])
=#

#scatter(chunk_list_timeseries.times,dcpca_result_espresso.rvs_clean.-mean(dcpca_result_espresso.rvs_clean))
#scatter!(chunk_list_timeseries.times,dcpca_result_espresso.rvs_dirty)

#=
plt = scatter(λmatrix,dcpca_result_espresso.model.mean,ms=1)
xlims!(5200,5400.0)
scatter!(plt,λmatrix,dcpca_result_espresso.model.mean,ms=2)
scatter!(plt,λmatrix,1.0.+dcpca_result_espresso.model.proj[:,1],ms=2)
plt =scatter(λmatrix,1.0.+dcpca_result_espresso.model.proj[:,1],ms=1)
xlims!(5200,5210.0)
scatter!(plt,λmatrix,1.0.+dcpca_result_espresso.model.proj[:,2],ms=1.5)
plot!(plt,λmatrix,1.0.+dcpca_result_espresso.model.proj[:,3],ms=1.5)
xlims!(5202.5,5204.16)
xlims!(5298,5302)
xlims!(5590,5593)
xlims!(5299,5300)
xlims!(5340,5400)
xlims!(5347,5348)
xlims!(5393,5395)
xlims!(6159,6162)
xlims!(6158,6165)
xlims!(5211.0,5213.)
xlims!(5445.0,5455.)
xlims!(5540.0,5555.)
xlims!(5546.0,5550.)
xlims!(5250.0,5450.)
xlims!(5440.0,5450.)
xlims!(5300.0,5310.)
xlims!(5430.0,5450.)

idxstart = sum(length.(chunk_grids)[1:94]); idxstop = idxstart + 30; scatter!(plt,λmatrix[idxstart:idxstop],1.0.+dcpca_result_espresso.model.proj[idxstart:idxstop,3],ms=2)

plt = scatter(λmatrix,1.0.+dcpca_result_espresso.model.proj[:,2],ms=1.5)
xlims!(5400,5580)
plt = plot(λmatrix,dcpca_result_espresso.model.mean,ms=1)
xlims!(5500,5580)
scatter!(plt,λmatrix,1.0.+10.0 .* dcpca_result_espresso.model.proj[:,2],ms=1.5)
ylims!(0.9,1.1)
xlims!(5500,5510)
=#



println("CCF custom mask common, LFC: ")
linelist_for_ccf_fn_w_path = joinpath("..","..","data",string(starid),psu_mask3_filename)
 line_list_clean1 = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=0*ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=false)
 σ_rvs_formal_clean1 = map(obsid->calc_formal_rv_precission(all_spectra[obsid],RvSpectMLBase.make_chunk_list_from_loc_df(all_spectra[obsid],get_inst(all_spectra), line_list_clean1),smooth_factor=2.0),1:length(all_spectra) )
 mean(σ_rvs_formal_clean1), median(σ_rvs_formal_clean1), extrema(σ_rvs_formal_clean1)
t_start_1 = now()
 ((ccfs_clean1, ccf_vars_clean1), v_grid) = ccf_total(order_list_timeseries, line_list_clean1, pipeline_plan,  mask_scale_factor=4.0, range_no_mask_change=line_width_05, ccf_mid_velocity=0*ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)

 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2.0)
 rvs_ccf_clean1 = calc_rvs_from_ccf_total(ccfs_clean1, ccf_vars_clean1, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
 σ_rvs_clean1 = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 #std(rvs_ccf_clean1)
 rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_clean1.-mean(rvs_ccf_clean1)))
 rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_clean1.-mean(rvs_ccf_clean1))
 println("# RMS of RVs: ", std(rvs_ccf_clean1), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night, " σ_rvs_formal = ", median(σ_rvs_formal_clean1))


 num_basis_scalpels = 2
   println("# CCF custom mask common, LFC, Scalpels, num_basis_vectors = ", num_basis_scalpels)
   ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
   @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
   CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels3_results.csv"), df_scalpels1)
   #plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
   rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
   rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
   println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

println("# CCF custom mask common, LFC, DCPCA, num_basis_vectors = ", num_basis_dcpca)
chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_telluric_free(all_spectra, line_list_clean1, min_pixels_in_chunk=5, verbose=false)
 dont_need_to!(pipeline_plan,:extract_orders)
 rv_offset_epoch4_minus_epoch5 = calc_epoch_offset_vector(rvs_ccf_clean1, σ_rvs_formal_clean1, tidx_epoch4, tidx_epoch5 )
 set_rv_est!(chunk_list_timeseries, rvs_ccf_clean1 .+ tidx_epoch5.*rv_offset_epoch4_minus_epoch5 )
 chunk_grids = RvSpectML.make_grids_for_chunklist_timeseries(chunk_list_timeseries)
 idx_use_pca = trues(length(all_spectra)); if starid == 101501 idx_use_pca[2] = false  end
 dcpca_result_clean1 = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_clean1, num_basis_vectors_pca=2, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=40)
 σ_rvs_clean_clean1 = σ_rvs_formal_clean1 # sqrt.(σ_rvs_formal_clean1.^2 .+ dcpca_result_clean1.σ_rvs_dirty.^2 )
 df_clean1 = make_results_df(order_list_timeseries.times,dcpca_result_clean1.rvs_clean,σ_rvs_clean_clean1,dcpca_result_clean1.rvs_dirty,NaN.+dcpca_result_clean1.σ_rvs_dirty,dcpca_result_clean1.scores,NaN.+dcpca_result_clean1.σ_scores)
 t_stop_1 = now()
 CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_dcpca3_results.csv"), df_clean1)
 t_stop_1-t_start_1
 std(rvs_ccf_clean1), std(dcpca_result_clean1.rvs_clean)
 rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_clean1[!,Symbol("RV_C [m/s]")].-mean(df_clean1[!,Symbol("RV_C [m/s]")])))
 rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_clean1[!,Symbol("RV_C [m/s]")].-mean(df_clean1[!,Symbol("RV_C [m/s]")]))
 println("# RMS of RVs: ", std(df_clean1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)


#=
println("CCF custom mask common: std(rvs) = ", std(rvs_ccf_clean1))
 plt = scatter(rvs_ccf_clean1.-mean(rvs_ccf_clean1), label="RV CCF")
  #for num_basis_dcpca in 1:6
  println("# Custom mask common, LFC, DCPCA, num_basis_vectors = ", num_basis_dcpca)
    dcpca_result_clean1 = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_clean1, num_basis_vectors_pca=num_basis_dcpca, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=0)
    #print("DCPCA ", num_basis_dcpca, ": " ) #std(rvs) = ", std(dcpca_result_espresso.rvs_clean))
    rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=dcpca_result_clean1.rvs_clean.-mean(dcpca_result_clean1.rvs_clean)))
    rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=dcpca_result_clean1.rvs_clean.-mean(dcpca_result_clean1.rvs_clean))
    println("# RMS of RVs: ", std(dcpca_result_clean1.rvs_clean), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

    scatter!(plt,dcpca_result_clean1.rvs_dirty, label=string(num_basis_dcpca))
  #end
  σ_rvs_clean_espresso = σ_rvs_formal_espresso # sqrt.(σ_rvs_formal_espresso.^2 .+ dcpca_result_espresso.σ_rvs_dirty.^2 )
  df_clean1_dcpca1 = make_results_df(order_list_timeseries.times,dcpca_result_espresso.rvs_clean,σ_rvs_clean_espresso,dcpca_result_espresso.rvs_dirty,NaN.+dcpca_result_espresso.σ_rvs_dirty,dcpca_result_espresso.scores,NaN.+dcpca_result_espresso.σ_scores)
  CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_dcpca3_results.csv"), df_clean1_dcpca1)
  display(plt)
=#

#=
using Plots
scatter(order_list_timeseries.times,rvs_ccf_espresso.-mean(rvs_ccf_espresso),yerr=σ_rvs_ccf_espresso)
scatter!(order_list_timeseries.times,rvs_ccf_clean1.-mean(rvs_ccf_clean1),yerr=σ_rvs_clean_clean1)
std(rvs_ccf_clean1)
scatter(order_list_timeseries.times,rvs_ccf_clean1.-mean(rvs_ccf_clean1) .- (rvs_ccf_espresso.-mean(rvs_ccf_espresso)),yerr=σ_rvs_ccf_espresso, xlabel="Time (d)", ylabel="ΔRV (m/s)")
=#

#=
num_basis_scalpels = 1
 println("# Custom mask common, LFC, Scalpels, num_basis_vectors = ", num_basis_scalpels)
 ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
 @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
 #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), df_scalpels1)
 plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
 rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
 rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
 println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

num_basis_scalpels = 2
  ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
  @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
  CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels3_results.csv"), df_scalpels1)
  plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
  rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
  rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
  println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 3
   ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
   @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
   #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), df_scalpels1)
   plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
   rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
   rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
   println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 4
    ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
    @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
    #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), df_scalpels1)
    plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
    rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
    rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
    println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)
 num_basis_scalpels = 5

     ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
     @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
     #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), df_scalpels1)
     plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
     rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
     rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
     println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 6
      ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_clean_clean1; println("# Setting ccfs to clean1 mask ", num_basis_scalpels, " basis vectors");
      @time df_scalpels1 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
      #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), df_scalpels1)
      plt = scatter(df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels1[!,Symbol("e_RV_C [m/s]")], label="1")
      rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")])))
      rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels1[!,Symbol("RV_C [m/s]")].-mean(df_scalpels1[!,Symbol("RV_C [m/s]")]))
      println("# RMS of RVs: ", std(df_scalpels1[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)
=#

println("CCF custom mask best, LFC: ")
linelist_for_ccf_fn_w_path = joinpath("..","..","data",string(starid),psu_mask4_filename)
  line_list_clean2 = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=0*ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=false)
  σ_rvs_formal_clean2 = map(obsid->calc_formal_rv_precission(all_spectra[obsid],RvSpectMLBase.make_chunk_list_from_loc_df(all_spectra[obsid],get_inst(all_spectra), line_list_clean2),smooth_factor=2.0),1:length(all_spectra) )
  mean(σ_rvs_formal_clean2), median(σ_rvs_formal_clean2), extrema(σ_rvs_formal_clean2)
  t_start_2 = now()
  ((ccfs_clean2, ccf_vars_clean2), v_grid) = ccf_total(order_list_timeseries, line_list_clean2, pipeline_plan,  mask_scale_factor=4.0, range_no_mask_change=line_width_05, ccf_mid_velocity=0*ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
  alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
  rvs_ccf_clean2 = calc_rvs_from_ccf_total(ccfs_clean2, ccf_vars_clean2, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
  σ_rvs_clean2 = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
  rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=rvs_ccf_clean2.-mean(rvs_ccf_clean2)))
  rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=rvs_ccf_clean2.-mean(rvs_ccf_clean2))
  println("# RMS of RVs: ", std(rvs_ccf_clean2), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night, " σ_rvs_formal = ", median(σ_rvs_formal_clean2))


  num_basis_scalpels = 2
    println("# CCF custom mask best, LFC, Scalpels, num_basis_vectors = ", num_basis_scalpels)
    ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean2; #println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
    @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
    CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels1)
    #plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
    rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
    rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
    println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)


println("# Custom mask best, DCPCA: ")
  chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_telluric_free(all_spectra, line_list_clean2, min_pixels_in_chunk=5, verbose=false)
  dont_need_to!(pipeline_plan,:extract_orders)
  rv_offset_epoch4_minus_epoch5 = calc_epoch_offset_vector(rvs_ccf_clean2, σ_rvs_formal_clean2, tidx_epoch4, tidx_epoch5 )
  set_rv_est!(chunk_list_timeseries, rvs_ccf_clean2 .+ tidx_epoch5.*rv_offset_epoch4_minus_epoch5 )
  chunk_grids = RvSpectML.make_grids_for_chunklist_timeseries(chunk_list_timeseries)
  idx_use_pca = trues(length(all_spectra)); if starid == 101501 idx_use_pca[2] = false  end
  dcpca_result_clean2 = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_clean2, num_basis_vectors_pca=2, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=40)
  σ_rvs_clean_clean2 = σ_rvs_formal_clean2 # sqrt.(σ_rvs_formal_clean2.^2 .+ dcpca_result_clean2.σ_rvs_dirty.^2 )
  df_clean2 = make_results_df(order_list_timeseries.times,dcpca_result_clean2.rvs_clean,σ_rvs_clean_clean2,dcpca_result_clean2.rvs_dirty,NaN.+dcpca_result_clean2.σ_rvs_dirty,dcpca_result_clean2.scores,NaN.+dcpca_result_clean2.σ_scores)
  t_stop_2 = now()
  rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_clean2[!,Symbol("RV_C [m/s]")].-mean(df_clean2[!,Symbol("RV_C [m/s]")])))
  rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_clean2[!,Symbol("RV_C [m/s]")].-mean(df_clean2[!,Symbol("RV_C [m/s]")]))
  println("# RMS of RVs: ", std(df_clean2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)
  CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_dcpca4_results.csv"), df_clean2)

  t_stop_2-t_start_2
#std(dcpca_result_espresso.rvs_clean)
#std(dcpca_result_clean1.rvs_clean), std(dcpca_result_clean2.rvs_clean)
#df_clean2

#=
num_basis_scalpels = 1
 ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean_clean2; println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
 @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
 #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels2)
 plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
 rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
 rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
 println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 2
  ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean_clean2; println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
  @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
  #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels2)
  plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
  rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
  rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
  println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 3
   ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean_clean2; println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
   @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
   #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels2)
   plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
   rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
   rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
   println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 4
    ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean_clean2; println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
    @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
    #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels2)
    plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
    rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
    rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
    println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 5
     ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean_clean2; println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
     @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
     #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels2)
     plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
     rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
     rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
     println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)

 num_basis_scalpels = 6
      ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_clean_clean2; println("# Setting ccfs to clean2 mask ", num_basis_scalpels, " basis vectors");
      @time df_scalpels2 = compute_scalpels_output(num_basis_scalpels, sort_by_responce=true)
      #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), df_scalpels2)
      plt = scatter(df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]), yerr=df_scalpels2[!,Symbol("e_RV_C [m/s]")], label="1")
      rms_rv_nightly = std(bin_rvs_nightly(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")])))
      rms_rv_within_night = rms_rvs_within_night(times=order_list_timeseries.times,rvs=df_scalpels2[!,Symbol("RV_C [m/s]")].-mean(df_scalpels2[!,Symbol("RV_C [m/s]")]))
      println("# RMS of RVs: ", std(df_scalpels2[!,Symbol("RV_C [m/s]")]), "  nightly RVs: ", rms_rv_nightly, "  within night: ",rms_rv_within_night)
=#



#=
println("std rvs = ", std(rvs_ccf_espresso) )
num_basis_scalpels = 1
 ccfs = ccfs_espresso; ccf_vars = ccf_vars_espresso; rvs_ccf = rvs_ccf_espresso; σ_rvs = σ_rvs_formal_espresso; println("# Setting ccfs to ESPRESSO mask ", num_basis_scalpels, " basis vectors");
 @time output_df1 = compute_scalpels_output(num_basis_scalpels,sort_by_responce=true)
 #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
 plt = scatter(output_df1[!,Symbol("RV_C [m/s]")].-mean(output_df1[!,Symbol("RV_C [m/s]")]), yerr=output_df1[!,Symbol("e_RV_C [m/s]")], label="1")

num_basis_scalpels = 2
 @time output_df = compute_scalpels_output(num_basis_scalpels)
 #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
 plt = scatter!(output_df[!,Symbol("RV_C [m/s]")].-mean(output_df[!,Symbol("RV_C [m/s]")])#= , yerr=output_df[!,Symbol("e_RV_C [m/s]")] =#, label="2")

num_basis_scalpels = 3
  @time output_df3 = compute_scalpels_output(num_basis_scalpels)
  #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
  plt = scatter(output_df3[!,Symbol("RV_C [m/s]")].-mean(output_df3[!,Symbol("RV_C [m/s]")]),  #=yerr=output_df3[!,Symbol("e_RV_C [m/s]")],=# label="3")

num_basis_scalpels = 4
  @time output_df4 = compute_scalpels_output(num_basis_scalpels)
  #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
  plt = scatter!(output_df4[!,Symbol("RV_C [m/s]")].-mean(output_df4[!,Symbol("RV_C [m/s]")]),  #=yerr=output_df4[!,Symbol("e_RV_C [m/s]")], =# label="4")

num_basis_scalpels = 5
  @time output_df5 = compute_scalpels_output(num_basis_scalpels)
  #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
  plt = scatter!(output_df5[!,Symbol("RV_C [m/s]")].-mean(output_df5[!,Symbol("RV_C [m/s]")]), #=yerr=output_df5[!,Symbol("e_RV_C [m/s]")], =# label="5")

num_basis_scalpels = 6
   @time output_df6 = compute_scalpels_output(num_basis_scalpels)
   #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
   plt = scatter!(output_df6[!,Symbol("RV_C [m/s]")].-mean(output_df6[!,Symbol("RV_C [m/s]")]), #=yerr=output_df6[!,Symbol("e_RV_C [m/s]")], =# label="6")

num_basis_scalpels = 7
     @time output_df6 = compute_scalpels_output(num_basis_scalpels)
     #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
     plt = scatter!(output_df6[!,Symbol("RV_C [m/s]")].-mean(output_df6[!,Symbol("RV_C [m/s]")]), #=yerr=output_df6[!,Symbol("e_RV_C [m/s]")], =# label="6")

num_basis_scalpels = 8
          @time output_df6 = compute_scalpels_output(num_basis_scalpels)
          #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
          plt = scatter!(output_df6[!,Symbol("RV_C [m/s]")].-mean(output_df6[!,Symbol("RV_C [m/s]")]), #=yerr=output_df6[!,Symbol("e_RV_C [m/s]")], =# label="6")

num_basis_scalpels = 9
                    @time output_df6 = compute_scalpels_output(num_basis_scalpels)
                    #CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels0_results.csv"), output_df)
                    plt = scatter!(output_df6[!,Symbol("RV_C [m/s]")].-mean(output_df6[!,Symbol("RV_C [m/s]")]), #=yerr=output_df6[!,Symbol("e_RV_C [m/s]")], =# label="6")


scatter!(dcpca_result_espresso.rvs_clean.-mean(dcpca_result_espresso.rvs_clean), yerr=σ_rvs_clean_espresso)
# scatter!(df_espresso[!,Symbol("RV_C [m/s]")], yerr=df_espresso[!,Symbol("e_RV_C [m/s]")])

(scalpels_score,scalpels_basis) = Scalpels.calc_basis_scores_scalpels(rvs_ccf, ccfs, σ_rvs = σ_rvs, num_basis=12, assume_centered=false, sort_by_responce=false )

scatter(scalpels_score[:,1])
scatter!(scalpes_score[:,2])
scatter!(scalpels_score[:,3])
scatter(scalpels_score[:,4])
scatter!(scalpels_score[:,5])
scatter(scalpels_score[:,6])

scatter(chunk_list_timeseries.times,scalpels_score[:,1])
scatter!(chunk_list_timeseries.times,scalpels_score[:,2])
scatter!(chunk_list_timeseries.times,scalpels_score[:,3])
scatter(chunk_list_timeseries.times,scalpels_score[:,4])
scatter!(chunk_list_timeseries.times,scalpels_score[:,5])
scatter(chunk_list_timeseries.times,scalpels_score[:,6])

plot(scalpels_basis[:,1])
plot!(scalpels_basis[:,2])
plot!(scalpels_basis[:,3])
plot!(scalpels_basis[:,4])
plot!(scalpels_basis[:,5])
plot!(scalpels_basis[:,6])
scatter!(scalpels_basis[:,7])
scatter!(scalpels_basis[:,8])
scatter!(scalpels_basis[:,9])
scatter!(scalpels_basis[:,10])
scatter!(scalpels_basis[:,11])
scatter!(scalpels_basis[:,12])

scatter(scalpels_score[:,4])

output_df[tidx_epoch4,Symbol("RV_C [m/s]")]
=#

 #=
 ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_formal_clean1;  println("# Setting ccfs to clean1 mask");
 output_df = compute_scalpels_output(num_basis_scalpels)
 CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels3_results.csv"), output_df)

 ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_formal_clean2;  println("# Setting ccfs to clean2 mask");
 output_df = compute_scalpels_output(num_basis_scalpels)
 CSV.write(joinpath("..","..","pennstate",string(starid),string(starid) * "_pennstate_scalpels4_results.csv"), output_df)
 =#

println("# Done")
#= =#

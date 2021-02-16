if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("EXPRES-StellarSignals")   end
using Pkg
Pkg.activate(".")

verbose = true
 make_plots = false
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectML
if verbose   println("# Loading other packages")    end
 using CSV, DataFrames, Query, StatsBase, Statistics, Dates
 if verbose   println("# Loading unpackaged code from this repo")    end
 include(joinpath(pkgdir(RvSpectML),"..","EXPRES-StellarSignals","src","io.jl"))

all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_expres_data_101501.jl"))
 if typeof(get_inst(all_spectra)) <: AnyEXPRES   continuum_normalize_spectra!(all_spectra)   end
 #if typeof(get_inst(all_spectra)) <: AnyEXPRES   EXPRES.blaze_normalize_spectra!(all_spectra)   end
 #RvSpectMLBase.discard_tellurics(all_spectra)
 RvSpectMLBase.discard_blaze(all_spectra)
 RvSpectMLBase.discard_continuum(all_spectra)
 GC.gc()
 ref_obs_idx = RvSpectMLBase.choose_obs_idx_for_init_guess(df_files_use,get_inst(all_spectra))
 tidx_epoch4 = df_files_use.expres_epoch.==4
 tidx_epoch5 = df_files_use.expres_epoch.==5

# Something buggy if include order 1.  Order 86 essentially ignored due to tellurics. First and last orders have NaN issues
lfc_orders = 43:72
  #max_orders = 12:83
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=lfc_orders, recalc=true )
  #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )
  order_list_timeseries_ref_obs_only = extract_chunklist_timeseries_with_subset_obs(order_list_timeseries, ref_obs_idx)

function calc_epoch_offset_vector(rvs_ccf, σ_rvs, tidx_epoch4, tidx_epoch5)
    vcat(rvs_ccf[tidx_epoch4].-mean(rvs_ccf[tidx_epoch4], weights(1.0 ./σ_rvs[tidx_epoch4].^2)), (rvs_ccf[tidx_epoch5].-mean(rvs_ccf[tidx_epoch5],weights(1.0 ./σ_rvs[tidx_epoch5].^2))))
end

linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks","G9.espresso.mas")
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=false)
 ((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries_ref_obs_only, line_list_espresso, pipeline_plan,  mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, v_step=100, recalc=true, calc_ccf_var=true)
 line_width_50 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.5)
 line_sigma = line_width_50/sqrt(2*log(2))
 line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)
 σ_rvs_formal_espresso = map(obsid->calc_formal_rv_precission(all_spectra[obsid],RvSpectMLBase.make_chunk_list_from_loc_df(all_spectra[obsid],get_inst(all_spectra), line_list_espresso),smooth_factor=2.0),1:length(all_spectra) )
 mean(σ_rvs_formal_espresso), median(σ_rvs_formal_espresso), extrema(σ_rvs_formal_espresso)

((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries_ref_obs_only, line_list_espresso, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, recalc=true, calc_ccf_var=true)
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf_ref = RvSpectML.calc_rvs_from_ccf_total( reshape(ccfs,size(ccfs,1),1), reshape(ccf_vars,size(ccf_vars,1),1), pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, bin_nightly=false, recalc=true, alg_fit_rv=alg_fit_rv)
 σ_rv_fit_ref = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 σ_rvs_formal_espresso
 ccf_var_scale_factor = σ_rvs_formal_espresso[ref_obs_idx]/σ_rv_fit_ref[1]

((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, ccf_var_scale=ccf_var_scale_factor, calc_ccf_var=true, recalc=true)
 # ccf_vars_espresso .*= ccf_var_scale_factor  # If wanted to scale ccf_var manually, rather than passing ccf_var_scale=ccf_var_scale_factor above
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, ccf_vars_espresso, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true,  bin_nightly=false, alg_fit_rv=alg_fit_rv)
 σ_rvs_ccf_espresso = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 rvs_ccf_espresso, σ_rvs_ccf_espresso
 σ_rvs_ccf_espresso./σ_rvs_formal_espresso

#=
# If you wanted to compute the full covariance matrix for the CCFs.  VERY slow.  And not that different.
((ccfs_espresso2, ccf_covars_espresso2), v_grid2) = ccf_total(order_list_timeseries_ref_obs_only, line_list_espresso, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, calc_ccf_covar=true, recalc=true)
alg_fit_rv2 =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
rvs_ccf_espresso2 = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso2, ccf_covars_espresso2, pipeline_plan, v_grid=v_grid2, times = order_list_timeseries.times, recalc=true, bin_nightly=false, alg_fit_rv=alg_fit_rv2)
rvs_ccf_espresso2, σ_rvs_ccf_espresso2
=#

#= =#

chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_telluric_free(all_spectra, line_list_espresso, min_pixels_in_chunk=5, verbose=false)
 dont_need_to!(pipeline_plan,:extract_orders)
 rv_offset_epoch4_minus_epoch5 = calc_epoch_offset_vector(rvs_ccf_espresso, σ_rvs_formal_espresso, tidx_epoch4, tidx_epoch5 )
 set_rv_est!(chunk_list_timeseries, rvs_ccf_espresso .+ tidx_epoch5.*rv_offset_epoch4_minus_epoch5 )
 chunk_grids = RvSpectML.make_grids_for_chunklist_timeseries(chunk_list_timeseries)
 idx_use_pca = trues(length(all_spectra)); idx_use_pca[2] = false
 dcpca_result_espresso = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_espresso, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=0)

σ_rvs_clean_espresso = σ_rvs_formal_espresso # sqrt.(σ_rvs_formal_espresso.^2 .+ dcpca_result_espresso.σ_rvs_dirty.^2 )
 df_espresso = make_results_df(order_list_timeseries.times,dcpca_result_espresso.rvs_clean,σ_rvs_clean_espresso,dcpca_result_espresso.rvs_dirty,NaN.+dcpca_result_espresso.σ_rvs_dirty,dcpca_result_espresso.scores,NaN.+dcpca_result_espresso.σ_scores)
 CSV.write(joinpath("pennstate","101501","101501_pennstate_dcpca0_results.csv"), df_espresso)


linelist_for_ccf_fn_w_path = joinpath("data","101501","HD101501_f75s75v3_mask.csv")
 line_list_clean1 = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=false)
 σ_rvs_formal_clean1 = map(obsid->calc_formal_rv_precission(all_spectra[obsid],RvSpectMLBase.make_chunk_list_from_loc_df(all_spectra[obsid],get_inst(all_spectra), line_list_clean1),smooth_factor=2.0),1:length(all_spectra) )
 mean(σ_rvs_formal_clean1), median(σ_rvs_formal_clean1), extrema(σ_rvs_formal_clean1)
 t_start_1 = now()
 ((ccfs_clean1, ccf_vars_clean1), v_grid) = ccf_total(order_list_timeseries, line_list_clean1, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
 rvs_ccf_clean1 = calc_rvs_from_ccf_total(ccfs_clean1, ccf_vars_clean1, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
 #σ_rvs_clean1 = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_telluric_free(all_spectra, line_list_clean1, min_pixels_in_chunk=5, verbose=false)
 dont_need_to!(pipeline_plan,:extract_orders)
 rv_offset_epoch4_minus_epoch5 = calc_epoch_offset_vector(rvs_ccf_clean1, σ_rvs_formal_clean1, tidx_epoch4, tidx_epoch5 )
 set_rv_est!(chunk_list_timeseries, rvs_ccf_clean1 .+ tidx_epoch5.*rv_offset_epoch4_minus_epoch5 )
 chunk_grids = RvSpectML.make_grids_for_chunklist_timeseries(chunk_list_timeseries)
 idx_use_pca = trues(length(all_spectra)); idx_use_pca[2] = false
 dcpca_result_clean1 = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_clean1, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=0)
 σ_rvs_clean_clean1 = σ_rvs_formal_clean1 # sqrt.(σ_rvs_formal_clean1.^2 .+ dcpca_result_clean1.σ_rvs_dirty.^2 )
 df_clean1 = make_results_df(order_list_timeseries.times,dcpca_result_clean1.rvs_clean,σ_rvs_clean_clean1,dcpca_result_clean1.rvs_dirty,NaN.+dcpca_result_clean1.σ_rvs_dirty,dcpca_result_clean1.scores,NaN.+dcpca_result_clean1.σ_scores)
 t_stop_1 = now()
 CSV.write(joinpath("pennstate","101501","101501_pennstate_dcpca1_results.csv"), df_clean1)
 t_stop_1-t_start_1

linelist_for_ccf_fn_w_path = joinpath("data","101501","HD101501_f80s85v2_mask.csv")
  line_list_clean2 = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=false)
  σ_rvs_formal_clean2 = map(obsid->calc_formal_rv_precission(all_spectra[obsid],RvSpectMLBase.make_chunk_list_from_loc_df(all_spectra[obsid],get_inst(all_spectra), line_list_clean2),smooth_factor=2.0),1:length(all_spectra) )
  mean(σ_rvs_formal_clean2), median(σ_rvs_formal_clean2), extrema(σ_rvs_formal_clean2)
  t_start_2 = now()
  ((ccfs_clean2, ccf_vars_clean2), v_grid) = ccf_total(order_list_timeseries, line_list_clean2, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
  alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
  rvs_ccf_clean2 = calc_rvs_from_ccf_total(ccfs_clean2, ccf_vars_clean2, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
  #σ_rvs_clean2 = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
  chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_telluric_free(all_spectra, line_list_clean2, min_pixels_in_chunk=5, verbose=false)
  dont_need_to!(pipeline_plan,:extract_orders)
  rv_offset_epoch4_minus_epoch5 = calc_epoch_offset_vector(rvs_ccf_clean2, σ_rvs_formal_clean2, tidx_epoch4, tidx_epoch5 )
  set_rv_est!(chunk_list_timeseries, rvs_ccf_clean2 .+ tidx_epoch5.*rv_offset_epoch4_minus_epoch5 )
  chunk_grids = RvSpectML.make_grids_for_chunklist_timeseries(chunk_list_timeseries)
  idx_use_pca = trues(length(all_spectra)); idx_use_pca[2] = false
  dcpca_result_clean2 = RvSpectML.DCPCA.clean_rvs_dcpca(chunk_list_timeseries, rvs_ccf_clean2, idx_use_pca=idx_use_pca, chunk_grids=chunk_grids, num_samples_σ_scores=0)
  σ_rvs_clean_clean2 = σ_rvs_formal_clean2 # sqrt.(σ_rvs_formal_clean2.^2 .+ dcpca_result_clean2.σ_rvs_dirty.^2 )
  df_clean2 = make_results_df(order_list_timeseries.times,dcpca_result_clean2.rvs_clean,σ_rvs_clean_clean2,dcpca_result_clean2.rvs_dirty,NaN.+dcpca_result_clean2.σ_rvs_dirty,dcpca_result_clean2.scores,NaN.+dcpca_result_clean2.σ_scores)
  t_stop_2 = now()
  CSV.write(joinpath("pennstate","101501","101501_pennstate_dcpca2_results.csv"), df_clean2)
  t_stop_2-t_start_2

std(dcpca_result_espresso.rvs_clean), std(dcpca_result_clean1.rvs_clean), std(dcpca_result_clean2.rvs_clean)


using Scalpels
function compute_scalpels_output(mnb::Integer)
  (rv_scalpels, scores, ) = calc_clean_rvs_scores_basis_scalpels(rvs_ccf, ccfs, num_basis=mnb, σ_rvs=σ_rvs)
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

num_basis_scalpels = 2
 ccfs = ccfs_espresso; ccf_vars = ccf_vars_espresso; rvs_ccf = rvs_ccf_espresso; σ_rvs = σ_rvs_formal_espresso; println("# Setting ccfs to ESPRESSO mask ", num_basis_scalpels, " basis vectors");
 @time output_df = compute_scalpels_output(num_basis_scalpels)
 CSV.write(joinpath("pennstate","101501","101501_pennstate_scalpels0_results.csv"), output_df)

 ccfs = ccfs_clean1; ccf_vars = ccf_vars_clean1; rvs_ccf = rvs_ccf_clean1; σ_rvs = σ_rvs_formal_clean1;  println("# Setting ccfs to clean1 mask");
 output_df = compute_scalpels_output(num_basis_scalpels)
 CSV.write(joinpath("pennstate","101501","101501_pennstate_scalpels1_results.csv"), output_df)

 ccfs = ccfs_clean2; ccf_vars = ccf_vars_clean2; rvs_ccf = rvs_ccf_clean2; σ_rvs = σ_rvs_formal_clean2;  println("# Setting ccfs to clean2 mask");
 output_df = compute_scalpels_output(num_basis_scalpels)
 CSV.write(joinpath("pennstate","101501","101501_pennstate_scalpels2_results.csv"), output_df)

println("# Done")
#= =#

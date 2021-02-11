if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("EXPRES-StellarSignals")   end
using Pkg
Pkg.activate(".")

verbose = true
 make_plots = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments, EchelleInstruments.EXPRES
 using EchelleCCFs
 using RvSpectML
 if verbose   println("# Loading other packages")    end
 using CSV, DataFrames, Query, StatsBase, Statistics, Dates

 target_subdir = "Sim1"   # USER: Replace with directory of your choice
  fits_target_str = "Sim1"
  output_dir = "examples/output/"
  paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/group/ebf11/default/ebf11/expres"]
  # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
  pipeline_plan = PipelinePlan()
  dont_make_plot!(pipeline_plan, :movie)

 reset_all_needs!(pipeline_plan)
 if need_to(pipeline_plan,:read_spectra)
    if verbose println("# Finding what data files are avaliable.")  end
    eval(read_data_paths(paths_to_search=paths_to_search_for_param))
    @assert isdefined(Main,:expres_data_path)
    df_files = make_manifest(expres_data_path, target_subdir, EXPRES )

    if verbose println("# Reading in customized parameters from param.jl.")  end
    eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))

    if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
    @time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, normalization=:continuum),eachrow(df_files_use))
    #@time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, store_blaze=true, store_continuum=true, store_pixel_mask=true,normalization=:continuum),eachrow(df_files_use))
    GC.gc()
    dont_need_to!(pipeline_plan,:read_spectra)
 end

 #@time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, store_blaze=true, store_continuum=true, store_pixel_mask=true, normalization=:continuum),eachrow(df_files_use))
  #GC.gc()

# Something buggy if include order 1.  Order 86 essentially ignored due to tellurics. First and last orders have NaN issues
#max_orders = 12:83
lfc_orders = 43:72
  #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=lfc_orders, recalc=true )

line_list_filename = "G2.espresso.mas"
linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)
 line_width_50_default = 10.0e3
 lsf_width = 2.0e3
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 2*RvSpectMLBase.max_bc+3*lsf_width,
   orders_to_use=lfc_orders, recalc=true, verbose=true)
 #CSV.write(joinpath("pennstate","Sim1","Sim1_line_list_espresso.csv"), line_list_espresso)

(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan,  mask_scale_factor=4.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true)

using Plots
ccfs_norm = EchelleCCFs.calc_normalized_ccfs( v_grid, ccfs, v_mid=ccf_mid_velocity, dv_min=RvSpectMLBase.max_bc,dv_max= 2*RvSpectMLBase.max_bc)
ccf_template = EchelleCCFs.calc_ccf_template(ccfs_norm, assume_normalized=true)
plot(v_grid,ccfs, labels=:none)
plot(v_grid,ccfs_norm , labels=:none)
plot(v_grid,(ccfs_norm.-ccf_template)[:,range(1,size(ccfs,2),step=10)], labels=:none)
heatmap(v_grid,1:size(ccfs_norm,2),ccfs_norm'.-ccf_template')
xlabel!("v (m/s)")
ylabel!("Observation ID")
title!("Sim1: CCF - <CCF>, ESPRESSO mask + SNR weights")
#savefig("ccf_heatmaps_Sim1.png")
line_width_50 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.5)
line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)

linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)
 line_list_espresso_old = copy(line_list_espresso)
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
   Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()), orders_to_use=lfc_orders, recalc=true, verbose=true)
 #CSV.write(joinpath("pennstate","Sim1","Sim1_line_list_espresso.csv"), line_list_espresso)
 #line_list_espresso.weight .= copy(line_list_espresso.weight_input)

#1.5 0.35,
#1.4, 0.4,
msf = 1.4 ; fwtf = 0.4
msf = 4.0; fwtf = 2.0
 println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
 ((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan,  mask_scale_factor=msf, range_no_mask_change=0*line_width_05+3*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=100,
   #v_max=RvSpectMLBase.max_bc,
   v_max= 0*RvSpectMLBase.max_bc+5*line_width_50,
   calc_ccf_var=true, recalc=true)
 #fwtf = 2.0; println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
  alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
  rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, ccf_vars_espresso, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
  σ_rvs_espresso = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))

ccfs_norm = EchelleCCFs.calc_normalized_ccfs( v_grid, ccfs_espresso, v_mid=ccf_mid_velocity, dv_min=2.5*line_width_50,dv_max= 2*RvSpectMLBase.max_bc)
ccf_template = EchelleCCFs.calc_ccf_template(ccfs_norm, assume_normalized=true)
plot(v_grid,ccfs_espresso, labels=:none)
plot(v_grid,ccfs_norm , labels=:none)
plot(v_grid,(ccfs_norm.-ccf_template)[:,range(1,size(ccfs,2),step=10)], labels=:none)
xlims!(ccf_mid_velocity-3*line_width_50,ccf_mid_velocity+3*line_width_50)
heatmap(v_grid,1:size(ccfs_norm,2),(ccfs_norm.-ccf_template)')
cols_to_fit = findlast(x->x<ccf_mid_velocity-3*line_width_50,v_grid):findfirst(x->x>ccf_mid_velocity+3*line_width_50,v_grid)
heatmap(view(v_grid,cols_to_fit),1:size(ccfs_norm,2),view(ccfs_norm,cols_to_fit,:)'.-view(ccf_template,cols_to_fit,:)')
xlabel!("v (m/s)")
ylabel!("Observation ID")
title!("Sim1: CCF - <CCF>, ESPRESSO mask + SNR weights")
#xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
#xlims!(ccf_mid_velocity-3*line_width_50,ccf_mid_velocity+3*line_width_50)
savefig("ccf_heatmaps_Sim1_espresso_msf=4.0.png")

#write_each_ccf_fits(map(s->s.metadata,all_spectra), v_grid, ccfs_espresso, ccf_vars_espresso,   fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
#            total_vel=rvs_ccf_espresso, sigma_total_vel=σ_rvs_espresso )




max_orders = 12:83
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )

line_list_filename = "G2.espresso.mas"
 linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
      Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),
      orders_to_use=max_orders, recalc=true, verbose=true)

(order_ccfs, order_ccf_vars, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_espresso, pipeline_plan, mask_scale_factor=4.0, ccf_mid_velocity=ccf_mid_velocity,
 v_step=100, v_max= 0*RvSpectMLBase.max_bc+5*line_width_50, orders_to_use=max_orders, calc_ccf_var=true,recalc=true)

obs_id = 1
heatmap(v_grid, max_orders, (order_ccfs[:,:,obs_id]./mean(order_ccfs[:,:,obs_id],dims=1))')
xlabel!("v (m/s)")
ylabel!("Order index")
title!("Sim1: CCF_" *string(obs_id) * ", ESPRESSO mask + SNR weights")
#xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
xlims!(ccf_mid_velocity-3*line_width_50,ccf_mid_velocity+3*line_width_50)

alg_fit_rv2 =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
order_rvs_g = zeros(length(all_spectra), length(max_orders))
order_rv_std_g = zeros(length(all_spectra), length(max_orders))
order_rvs_t = zeros(length(all_spectra), length(max_orders))
order_rv_std_t = zeros(length(all_spectra), length(max_orders))
for (i,ord) in enumerate(max_orders)
   println("# Order: ", i)
   if sum(view(order_ccfs,:,i,1)) > 0
      rvs_order_ccf = RvSpectML.calc_rvs_from_ccf_total(view(order_ccfs,:,i,:), view(order_ccf_vars,:,i,:), pipeline_plan, v_grid=v_grid_order_ccfs, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv2)
      alg_fit_rv3 =  EchelleCCFs.MeasureRvFromCCFTemplate(v_grid=v_grid_order_ccfs, frac_of_width_to_fit=fwtf, template = vec(mean(view(order_ccfs,:,i,:),dims=2)) )
      order_rvs_g[:,i] .= read_cache(pipeline_plan, :rvs_ccf_total )
      order_rv_std_g[:,i] .= read_cache(pipeline_plan, :σ_rvs_ccf_total)
      rvs_order_ccf = RvSpectML.calc_rvs_from_ccf_total(view(order_ccfs,:,i,:), view(order_ccf_vars,:,i,:), pipeline_plan, v_grid=v_grid_order_ccfs, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv3)
      order_rvs_t[:,i] .= read_cache(pipeline_plan, :rvs_ccf_total )
      order_rv_std_t[:,i] .= read_cache(pipeline_plan, :σ_rvs_ccf_total)
   end
end


#EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), v_grid, ccfs_espresso, ccf_vars_espresso,   fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
#                  orders=161.-max_orders,ccf_orders=order_ccfs,ccf_orders_var=order_ccf_vars,
#                  total_vel=rvs_ccf_espresso, sigma_total_vel=σ_rvs_espresso,
#                  order_vels = order_rvs_g, sigma_order_vels = order_rv_std_g )

EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), 100.0 .* v_grid, ccfs_espresso, ccf_vars_espresso,   fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
                  line_list_filename=line_list_filename, orders=161 .-max_orders,ccf_orders=order_ccfs,ccf_orders_var=order_ccf_vars,
                  total_vel=100.0 .* rvs_ccf_espresso, sigma_total_vel= 100.0 .* σ_rvs_espresso,
                  order_vels = 100.0 .* order_rvs_g, sigma_order_vels = 100.0 .* order_rv_std_g )

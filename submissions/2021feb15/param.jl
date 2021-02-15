global max_spectra_to_use = 200
if max_spectra_to_use < 200
   @warn "param.in setting max_spectra_to_use to " * string(max_spectra_to_use)
end

global tophap_ccf_mask_scale_factor=1.6

global fits_target_str
if fits_target_str == "Solar"
   global linelist_for_ccf_filename = "G2.espresso.mas"
   hostname = gethostname()
   if occursin("aci.ics.psu.edu",hostname)
      global ancilary_solar_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar"
   elseif occursin("nuc8",hostname)  # Eric's home machine :)
      global ancilary_solar_data_path = "/home/eford/Data/SolarSpectra/NEID_solar/"
   end
   global ccf_mid_velocity = 0
   global bjd_first_good = 2458745.1296134139
   global bjd_last_good = 2458745.283
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      @take(max_spectra_to_use) |>
      DataFrame
end





if fits_target_str == "101501"
   global espresso_mask_filename = "G9.espresso.mas"
   global psu_mask3_filename = "HD101501q90f95n8e=false_mask.csv"
   global psu_mask4_filename = "HD101501q45f95n8e=false_mask.csv"

   global ccf_mid_velocity = -5.0e3
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @take(max_spectra_to_use) |>
      DataFrame
end

if fits_target_str == "10700"
   global espresso_mask_filename = "G9.espresso.mas"
   global psu_mask3_filename = "HD10700q90f80n7e=false_mask.csv"
   global psu_mask4_filename = "HD10700q50f80n7e=false_mask.csv"
   global ccf_mid_velocity = -16640.0
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @take(max_spectra_to_use) |>
      DataFrame
end

if fits_target_str == "34411"  #G1?
   global espresso_mask_filename = "G2.espresso.mas"
   global psu_mask3_filename = "HD34411q90f80n0e=false_mask.csv"
   global psu_mask4_filename = "HD34411q65f80n0e=false_mask.csv"
   global ccf_mid_velocity = 66500.0
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @take(max_spectra_to_use) |>
      DataFrame
end

if fits_target_str == "26965" # K0?
   global espresso_mask_filename = "G9.espresso.mas"
   global psu_mask3_filename = "HD26965q90f95n3e=false_mask.csv"
   global psu_mask4_filename = "HD26965q55f95n3e=false_mask.csv"
   global ccf_mid_velocity = -40320.0
   global df_files
   global df_files_use = df_files |>
      @filter( _.target == fits_target_str ) |>
      @take(max_spectra_to_use) |>
      DataFrame
end

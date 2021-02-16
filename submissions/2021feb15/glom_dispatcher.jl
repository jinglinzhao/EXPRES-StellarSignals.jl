using Pkg
Pkg.activate("submissions/2021feb15")
Pkg.instantiate()

# for std()
using Statistics

# for det()
using LinearAlgebra

# for importing the data from CSV
using DataFrames
using CSV

# For GLOM
import GPLinearODEMaker; GLOM = GPLinearODEMaker

# For this module
using GLOM_RV_Example; GLOM_RV = GLOM_RV_Example

# For units in orbit fitting functions
using UnitfulAstro, Unitful

using Plots


rv_indicators_type = "dcpca0"
kernel_choice = 3
#starid = 101501
#starid = 10700
starid = 26965
#starid = 34411
if starid == 101501
    star_rot_rate = 16.28  # days
elseif starid == 10700
    star_rot_rate = 34.0  # days  # wikipedia, how embarassing
elseif starid == 26965
    star_rot_rate = 37.1  # days  #  Saar & Osten (1997) via https://iopscience.iop.org/article/10.3847/1538-3881/aaa896
elseif starid == 34411
    star_rot_rate = 27.9  # days # https://onlinelibrary.wiley.com/doi/full/10.1002/asna.202013790
end
 init_timescale_1 = star_rot_rate/4
 init_timescale_2 = 1.0/24
 tighten_lengthscale_priors = 3

#=
if starid == 101501
    star_rot_rate = 16.28  # days
elseif starid == 10700
    star_rot_rate = 34.0  # days  # wikipedia, how embarassing
elseif starid == 26965
    star_rot_rate = 37.1  # days  #  Saar & Osten (1997) via https://iopscience.iop.org/article/10.3847/1538-3881/aaa896
elseif starid == 34411
    star_rot_rate = 27.9  # days # https://onlinelibrary.wiley.com/doi/full/10.1002/asna.202013790
end

init_timescale_1 = star_rot_rate/4
init_timescale_2 = 1.0/24
tighten_lengthscale_priors = 3

include("glom_generic.jl")
=#

starid_list = [101501, 26965, 34411, 10700]
rv_indicators_type_list = ["dcpca0", "scalpels0", "dcpca3", "scalpels3", "dcpca4", "scalpels4" ]
kernel_choice_list = [3, 5, 4]

#for kernel_choice in kernel_choice_list
kernel_choice = 3
    for starid in starid_list
#starid = 26965
        if starid == 101501
            star_rot_rate = 16.28  # days
        elseif starid == 10700
            star_rot_rate = 34.0  # days  # wikipedia, how embarassing
        elseif starid == 26965
            star_rot_rate = 37.1  # days  #  Saar & Osten (1997) via https://iopscience.iop.org/article/10.3847/1538-3881/aaa896
        elseif starid == 34411
            star_rot_rate = 27.9  # days # https://onlinelibrary.wiley.com/doi/full/10.1002/asna.202013790
        end
        init_timescale_1 = star_rot_rate/4
        init_timescale_2 = 1.0/24
        if kernel_choice == 4 ||  kernel_choice == 6
            tighten_lengthscale_priors = 3
        else
            tighten_lengthscale_priors = 1
        end

        for rv_indicators_type in rv_indicators_type_list
            #rv_indicators_type = "dcpca0"

            include("glom_generic.jl")
            run_glom_rv(starid, rv_indicators_type, kernel_choice, star_rot_rate, init_timescale_1, init_timescale_2, tighten_lengthscale_priors)
        end
    end
#end

starid = 26965
rv_indicators_type = "dcpca0"
 include("glom_generic.jl")
rv_indicators_type = "scalpels0"
 include("glom_generic.jl")
rv_indicators_type = "scalpels3"
 include("glom_generic.jl")
rv_indicators_type = "scalpels4"
 include("glom_generic.jl")
rv_indicators_type = "dcpca3"
 include("glom_generic.jl")
rv_indicators_type = "dcpca4"
 include("glom_generic.jl")

starid = 34411
 rv_indicators_type = "dcpca0"
  include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
  include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
  include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
  include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
  include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
  include("glom_generic.jl")

starid = 101501
 rv_indicators_type = "dcpca0"
   include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
    include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
    include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
    include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
    include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
    include("glom_generic.jl")

starid = 10700
 rv_indicators_type = "dcpca0"
       include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
        include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
        include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
        include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
        include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
        include("glom_generic.jl")

# NEW KERNEL
kernel_choice = 5
starid = 26965
rv_indicators_type = "dcpca0"
 include("glom_generic.jl")
rv_indicators_type = "scalpels0"
 include("glom_generic.jl")
rv_indicators_type = "scalpels3"
 include("glom_generic.jl")
rv_indicators_type = "scalpels4"
 include("glom_generic.jl")
rv_indicators_type = "dcpca3"
 include("glom_generic.jl")
rv_indicators_type = "dcpca4"
 include("glom_generic.jl")

starid = 34411
 rv_indicators_type = "dcpca0"
  include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
  include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
  include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
  include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
  include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
  include("glom_generic.jl")

starid = 101501
 rv_indicators_type = "dcpca0"
   include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
    include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
    include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
    include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
    include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
    include("glom_generic.jl")

starid = 10700
 rv_indicators_type = "dcpca0"
       include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
        include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
        include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
        include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
        include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
        include("glom_generic.jl")


        # NEW KERNEL
kernel_choice = 4
 starid = 26965
 rv_indicators_type = "dcpca0"
 include("glom_generic.jl")
 rv_indicators_type = "scalpels0"
 include("glom_generic.jl")
 rv_indicators_type = "scalpels3"
 include("glom_generic.jl")
 rv_indicators_type = "scalpels4"
 include("glom_generic.jl")
 rv_indicators_type = "dcpca3"
 include("glom_generic.jl")
 rv_indicators_type = "dcpca4"
 include("glom_generic.jl")

starid = 34411
         rv_indicators_type = "dcpca0"
          include("glom_generic.jl")
         rv_indicators_type = "scalpels0"
          include("glom_generic.jl")
         rv_indicators_type = "scalpels3"
          include("glom_generic.jl")
         rv_indicators_type = "scalpels4"
          include("glom_generic.jl")
         rv_indicators_type = "dcpca3"
          include("glom_generic.jl")
         rv_indicators_type = "dcpca4"
          include("glom_generic.jl")

        starid = 101501
         rv_indicators_type = "dcpca0"
           include("glom_generic.jl")
         rv_indicators_type = "scalpels0"
            include("glom_generic.jl")
         rv_indicators_type = "scalpels3"
            include("glom_generic.jl")
         rv_indicators_type = "scalpels4"
            include("glom_generic.jl")
         rv_indicators_type = "dcpca3"
            include("glom_generic.jl")
         rv_indicators_type = "dcpca4"
            include("glom_generic.jl")

        starid = 10700
         rv_indicators_type = "dcpca0"
               include("glom_generic.jl")
         rv_indicators_type = "scalpels0"
                include("glom_generic.jl")
         rv_indicators_type = "scalpels3"
                include("glom_generic.jl")
         rv_indicators_type = "scalpels4"
                include("glom_generic.jl")
         rv_indicators_type = "dcpca3"
                include("glom_generic.jl")
         rv_indicators_type = "dcpca4"
                include("glom_generic.jl")

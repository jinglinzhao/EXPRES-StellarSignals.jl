using Pkg
if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("EXPRES-StellarSignals")   end
cd("submissions/2021feb15")
Pkg.activate(".")
starid = 101501
include("scr_eford_round2_generic.jl")

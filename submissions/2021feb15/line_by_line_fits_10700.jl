using Pkg
cd("EXPRES-StellarSignals")
Pkg.activate(".")
cd("submissions/2021feb15")
starid = 10700
include("line_by_line_fits_generic.jl")

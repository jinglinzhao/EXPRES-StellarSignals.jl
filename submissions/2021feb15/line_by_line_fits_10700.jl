using Pkg
cd("EXPRES-StellarSignals")
cd("submissions/2021feb15")
Pkg.activate(".")
starid = 10700
include("line_by_line_fits_generic.jl")

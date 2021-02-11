using Pkg
cd("EXPRES-StellarSignals")
cd("submissions/2021feb15")
Pkg.activate(".")
starid = 34411
include("line_by_line_fits_generic.jl")

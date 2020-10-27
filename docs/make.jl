using EXPRES-StellarSignals
using Documenter

makedocs(;
    modules=[EXPRES-StellarSignals],
    authors="Eric Ford",
    repo="https://github.com/eford/EXPRES-StellarSignals.jl/blob/{commit}{path}#L{line}",
    sitename="EXPRES-StellarSignals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://eford.github.io/EXPRES-StellarSignals.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/eford/EXPRES-StellarSignals.jl",
)

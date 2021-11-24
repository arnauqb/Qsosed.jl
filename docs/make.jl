using qsosed
using Documenter

DocMeta.setdocmeta!(qsosed, :DocTestSetup, :(using qsosed); recursive=true)

makedocs(;
    modules=[qsosed],
    authors="Arnau Quera-Bofarull",
    repo="https://github.com/arnauqb/qsosed.jl/blob/{commit}{path}#{line}",
    sitename="qsosed.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

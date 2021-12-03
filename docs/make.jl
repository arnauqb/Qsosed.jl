using Qsosed
using Documenter

DocMeta.setdocmeta!(Qsosed, :DocTestSetup, :(using Qsosed); recursive=true)

makedocs(;
    modules=[Qsosed],
    authors="Arnau Quera-Bofarull",
    repo="https://github.com/arnauqb/qsosed.jl/blob/{commit}{path}#{line}",
    sitename="Qsosed.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

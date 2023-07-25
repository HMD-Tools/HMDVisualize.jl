using HMDVisualize
using Documenter

DocMeta.setdocmeta!(HMDVisualize, :DocTestSetup, :(using HMDVisualize); recursive=true)

makedocs(;
    modules=[HMDVisualize],
    authors="Atsushi Yoshida <atsushi.yoshida.github@gmail.com>",
    repo="https://github.com/yoshiatsu163/HMDVisualize.jl/blob/{commit}{path}#{line}",
    sitename="HMDVisualize.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yoshiatsu163.github.io/HMDVisualize.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yoshiatsu163/HMDVisualize.jl",
    devbranch="master",
)

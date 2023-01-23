using PhysicsTools
using Documenter

DocMeta.setdocmeta!(PhysicsTools, :DocTestSetup, :(using PhysicsTools); recursive=true)

makedocs(;
    modules=[PhysicsTools],
    authors="Christian Haack",
    repo="https://github.com/PLEnuM-group/PhysicsTools.jl/blob/{commit}{path}#{line}",
    sitename="PhysicsTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://PLEnuM-group.github.io/PhysicsTools.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/PLEnuM-group/PhysicsTools.jl",
    devbranch="main"
)

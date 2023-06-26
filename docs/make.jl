using EduMaterialModels
using Documenter

DocMeta.setdocmeta!(EduMaterialModels, :DocTestSetup, :(using EduMaterialModels); recursive=true)

makedocs(;
    modules=[EduMaterialModels],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/EduMaterialModels.jl/blob/{commit}{path}#{line}",
    sitename="EduMaterialModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KnutAM.github.io/EduMaterialModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/EduMaterialModels.jl",
    devbranch="main",
)

using EduMaterialModels
using PlutoSliderServer
using Documenter

DocMeta.setdocmeta!(EduMaterialModels, :DocTestSetup, :(using EduMaterialModels); recursive=true)

const is_CI = get(ENV, "CI", "false") == "true"

const nb_dir = joinpath(@__DIR__, "src", "pluto_notebooks")

notebooks = [
    "Plasticity" => "pluto_notebooks/plasticity_md.md",
    "Viscoplasticity" => "pluto_notebooks/viscoplasticity_md.md",
    "Viscoelasticity" => "pluto_notebooks/viscoelasticity_md.md"
]
if Sys.iswindows()
    @warn("Cannot create html notebooks on windows")
else
    # Create static html versions
    PlutoSliderServer.export_directory(nb_dir; Export_create_index=false)
end

makedocs(;
    modules=[EduMaterialModels],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/EduMaterialModels.jl/blob/{commit}{path}#{line}",
    sitename="EduMaterialModels.jl",
    format=Documenter.HTML(;
        prettyurls=is_CI,
        canonical="https://KnutAM.github.io/EduMaterialModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Notebooks" => notebooks,
        "Material Models" => "models.md",
        "Utility functions" => "utils.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/EduMaterialModels.jl",
    devbranch="main",
)

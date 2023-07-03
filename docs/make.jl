using EduMaterialModels
using PlutoStaticHTML
using Documenter

DocMeta.setdocmeta!(EduMaterialModels, :DocTestSetup, :(using EduMaterialModels); recursive=true)

const is_CI = get(ENV, "CI", "false") == "true"

const nb_dir = joinpath(@__DIR__, "src", "pluto_notebooks")
notebooks = [file[1:end-3] for file in readdir(nb_dir) if endswith(file, ".jl")]
nb_md(nb) = joinpath("pluto_notebooks", nb*".md")

for nb in notebooks
    try
        build_notebooks(BuildOptions(nb_dir; output_format=documenter_output), [nb*".jl"])
    catch e
        @info "building notebook failed with error $e"
        @info "Outputting an empty markdown file"
        open(nb_md(nb), "w") do md
            println(md, "notebook is missing")
        end
    end
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
        "Notebooks" => [uppercasefirst(nb) => nb_md(nb) for nb in notebooks],
        "Material Models" => "models.md",
        "Utility functions" => "utils.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/EduMaterialModels.jl",
    devbranch="main",
)

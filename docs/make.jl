using Documenter
using BusInjectionModel

makedocs(
    sitename = "BusInjectionModel",
    format = Documenter.HTML(),
    modules = [BusInjectionModel],
    workdir = joinpath(@__DIR__, ".."),
    pages = [
        "User Documentation" => "index.md",
        "Methods" => "methods.md",
        "Math" => "math.md"
    ],
    warnonly = true,  # TODO rm this and fix all the docs
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/NLaws/BusInjectionModel.jl"
)

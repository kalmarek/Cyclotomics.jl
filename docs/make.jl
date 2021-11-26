push!(LOAD_PATH, "../src/")
using Cyclotomics
using Documenter

makedocs(
    sitename = "Cyclotomics.jl",
    modules = [Cyclotomics],
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/kalmarek/Cyclotomics.jl")

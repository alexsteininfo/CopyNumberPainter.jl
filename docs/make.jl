using Documenter
using CopyNumberEvolution

DocMeta.setdocmeta!(CopyNumberEvolution, :DocTestSetup,
                    :(using CopyNumberEvolution); recursive = true)

makedocs(
    sitename = "CopyNumberEvolution.jl",
    modules = [CopyNumberEvolution],
    authors = "Alexander Stein",
    pages = [
        "Home" => "index.md",
        "Concepts" => "concepts.md",
        "Input: trees" => "trees.md",
        "The alteration model" => "model.md",
        "Output" => "output.md",
        "MEDICC2 interoperability" => "interop.md",
        "MEDICC2 benchmark" => "benchmark.md",
        "Limitations and open questions" => "limitations.md",
        "API reference" => "api.md",
    ],
    doctest = true,
    checkdocs = :exports,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true",
                            size_threshold = 400 * 2^10, size_threshold_warn = 250 * 2^10),
)

deploydocs(repo = "github.com/alexander-stein/CopyNumberEvolution.jl.git",
           devbranch = "main", push_preview = false)

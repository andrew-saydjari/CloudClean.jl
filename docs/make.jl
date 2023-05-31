using CloudClean
using Documenter

makedocs(
    modules = [CloudClean],
    clean=true,
    highlightsig = true,
    sitename= "CloudClean.jl",
    format = Documenter.HTML(
        # assets = ["assets/favicon.ico"],
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages    = [
        "Introduction" => "index.md",
        "API Reference" => "api.md",
        # "Contributing" => "contrib.md"
    ]
)

deploydocs(
    repo = "github.com/andrew-saydjari/CloudClean.jl.git",
    branch = "gh-pages",
    devbranch = "main"
)

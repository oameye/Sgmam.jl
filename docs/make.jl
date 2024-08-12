CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using Sgmam
using Documenter
using DocumenterVitepress
using DocumenterCitations

DocMeta.setdocmeta!(Sgmam, :DocTestSetup, :(using Sgmam); recursive=true)

# using Downloads: Downloads
# Downloads.download(
#   "https://raw.githubusercontent.com/oameye/sgmam.jl/master/README.md",
#   joinpath(@__DIR__, "src/index.md"),
# )

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric  # default
)

makedocs(;
  modules=[Sgmam],
  repo=Remotes.GitHub("oameye", "Sgmam.jl"),
  authors="Orjan Ameye <orjan.ameye@hotmail.com>",
  sitename="Sgmam.jl",
  format=DocumenterVitepress.MarkdownVitepress(;
    repo="github.com/oameye/Sgmam.jl",
    devbranch = "master",
    devurl = "dev",
  ),
  pages=["Home" => "index.md", "API" => "API.md"],
  plugins = [bib,],
)

if CI
  deploydocs(;
   repo="github.com/oameye/Sgmam.jl",
   devbranch="master",
   target = "build",
   branch = "gh-pages",
   push_preview = true,)
end

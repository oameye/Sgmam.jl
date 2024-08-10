CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using Sgmam
using Documenter

DocMeta.setdocmeta!(Sgmam, :DocTestSetup, :(using Sgmam); recursive=true)

using Downloads: Downloads
Downloads.download(
  "https://raw.githubusercontent.com/oameye/sgmam.jl/master/README.md",
  joinpath(@__DIR__, "src/index.md"),
)

makedocs(;
  modules=[Sgmam],
  authors="Orjan Ameye <orjan.ameye@hotmail.com>",
  sitename="Sgmam.jl",
  format=Documenter.HTML(;
    canonical = "https://oameye.github.io/Sgmam.jl",
    edit_link = "master",
    assets    = String[],
  ),
  pages=["Home" => "index.md", "API" => "api.md"],
)

if CI
  deploydocs(; repo="github.com/oameye/Sgmam.jl", devbranch="master")
end

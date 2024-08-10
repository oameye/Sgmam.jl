using Sgmam
using Documenter

DocMeta.setdocmeta!(Sgmam, :DocTestSetup, :(using Sgmam); recursive=true)

makedocs(;
  modules=[Sgmam],
  authors="Orjan Ameye <orjan.ameye@hotmail.com>",
  sitename="Sgmam.jl",
  format=Documenter.HTML(;
    canonical = "https://oameye.github.io/Sgmam.jl",
    edit_link = "master",
    assets    = String[],
  ),
  pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/oameye/Sgmam.jl", devbranch="master")

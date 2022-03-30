using Ascender
using Documenter

DocMeta.setdocmeta!(Ascender, :DocTestSetup, :(using Ascender); recursive=true)

makedocs(;
    modules=[Ascender],
    authors="Yaozhenghang Ma <yaozhenghang.ma@gmail.com> and contributors",
    repo="https://github.com/yaozhenghangma/Ascender.jl/blob/{commit}{path}#{line}",
    sitename="Ascender.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yaozhenghangma.github.io/Ascender.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yaozhenghangma/Ascender.jl",
    devbranch="main",
)

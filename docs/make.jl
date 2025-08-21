using Documenter, SMGReader

makedocs(sitename="SMGReader.jl",
    pages =["index.md",
    "Chromatin Stencilling" => "chromstencil.md",
    "Direct RNA" => "directrna.md",
    "API Reference" => "api.md"])
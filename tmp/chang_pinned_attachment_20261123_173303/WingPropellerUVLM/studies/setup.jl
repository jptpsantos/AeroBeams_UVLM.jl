# Run once from any directory: julia --project=<studies> <studies>/setup.jl
using Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__) do
    Pkg.develop(path="..")
end
Pkg.instantiate()

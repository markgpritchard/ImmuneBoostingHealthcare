
using Pkg
Pkg.add(Pkg.PackageSpec(; name="CSV", version="0.10.14"))
Pkg.add(Pkg.PackageSpec(; name="CairoMakie", version="0.12.3"))
Pkg.add(Pkg.PackageSpec(; name="CategoricalArrays", version="0.10.8"))
Pkg.add(Pkg.PackageSpec(; name="DataFrames", version="1.6.1"))
Pkg.add(Pkg.PackageSpec(; name="DifferentialEquations", version="7.10.0"))
Pkg.add(Pkg.PackageSpec(; name="Distributions", version="0.25.109"))
Pkg.add(Pkg.PackageSpec(; name="DrWatson", version="2.15.0"))
Pkg.add(Pkg.PackageSpec(; name="Pigeons", version="0.4.2"))
Pkg.add(Pkg.PackageSpec(; name="StaticArrays", version="1.9.7"))
Pkg.add(Pkg.PackageSpec(; name="Turing", version="0.33.0"))

using DrWatson
@quickactivate "ImmuneBoostingHealthcare"
using Pkg
Pkg.instantiate()

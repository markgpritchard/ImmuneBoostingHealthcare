
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CairoMakie, MakieTeX 

healthcaredag = let 
    td1 = TeXDocument(read(scriptsdir("dag.tex"), String))
    fig = Figure(; size=( 500, 400 ))
    ga = GridLayout(fig[1, 1:2])
    lt1 = LTeX(ga[1, 1], td1; tellwidth=false)

    fig
end

safesave(plotsdir("healthcaredag.pdf"), healthcaredag)


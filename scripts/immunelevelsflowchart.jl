
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CairoMakie, MakieTeX 

immunelevelsflowchart = let 
    td1 = TeXDocument(read(scriptsdir("immunelevelsflowchart.tex"), String))
    fig = Figure(; size=( 300, 300 ))
    ga = GridLayout(fig[1, 1:2])
    lt1 = LTeX(ga[1, 1], td1; tellwidth=false)
    fig
end

safesave(plotsdir("immunelevelsflowchart.svg"), immunelevelsflowchart)


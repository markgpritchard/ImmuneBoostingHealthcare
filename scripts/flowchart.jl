
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CairoMakie, MakieTeX 

ch4flowchart = let 
    td1 = TeXDocument(read(scriptsdir("tikz_chapter4.tex"), String))
    fig = Figure(; size=( 500, 800 ))
    ga = GridLayout(fig[1, 1:2])
    lt1 = LTeX(ga[1, 1], td1; tellwidth=false)

    fig
end

safesave(plotsdir("ch4flowchart.pdf"), ch4flowchart)



using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))

using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions

# ~~~~~~~~
#
# 

dates = [ Date("2020-03-18") + Day(t) for t ∈ 1:831 ]
s = zeros(831)
s[1] = 1.0
v = zeros(831)
for t ∈ 2:831 
    if t < 500 
        s[t] = s[t-1]
    else
        vacc = vaccinate(dates[t]) * s[t-1]
        v[t] = v[t-1] + vacc
        s[t] = s[t-1] - vacc
    end
end 

lines(dates, v)

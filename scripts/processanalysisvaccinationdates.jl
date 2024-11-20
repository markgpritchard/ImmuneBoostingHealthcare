
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include("analysesimssetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alternativevaccinations = let 
    minus2months = [ vaccinatestaff(t; boostdateoffset=-62) for t ∈ 1:ndates ] 
    minus1month = [ vaccinatestaff(t; boostdateoffset=-31) for t ∈ 1:ndates ] 
    plus1month = [ vaccinatestaff(t; boostdateoffset=30) for t ∈ 1:ndates ] 
    plus2months = [ vaccinatestaff(t; boostdateoffset=61) for t ∈ 1:ndates ] 
    @strdict minus2months minus1month plus1month plus2months
end

safesave(datadir("sims", "alternativevaccinations.jld2"), alternativevaccinations)

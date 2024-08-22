
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CategoricalArrays, CSV, DataFrames, Dates


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load UK data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hospitaldata = CSV.read(datadir("exp_raw", "SiteData.csv"), DataFrame)
rename!(hospitaldata, Dict(Symbol("Trust Code") => "TrustCode"))
rename!(hospitaldata, Dict(Symbol("Trust Type") => "TrustType"))
rename!(hospitaldata, Dict(Symbol("Site Code") => "SiteCode"))
rename!(hospitaldata, Dict(Symbol("Site Type") => "SiteType"))
rename!(hospitaldata, Dict(Symbol("Site heated volume (m³)") => "HeatedVolumeString"))
rename!(
    hospitaldata, 
    Dict(
        Symbol("Single bedrooms for patients with en-suite facilities (No.)") => 
            "SingleBedsString"
    )
)
filter!(:SiteType => x -> x[1] == '1' || x[1] == '2', hospitaldata)
insertcols!(
    hospitaldata,
    :HeatedVolume => [ 
        parse(Int, replace(x, ',' => "")) 
        for x ∈ hospitaldata.HeatedVolumeString
    ], 
    :SingleBedsEnsuite => [ 
        x == "Not Applicable" ? 
            missing :
            parse(Int, replace(x, ',' => "")) 
        for x ∈ hospitaldata.SingleBedsString
    ]
)

hospitalbeds = CSV.read(datadir("exp_raw", "GeneralAcuteOccupiedBedsbyTrust.csv"), DataFrame)
filter!(:TotalBedsAvailable => x -> !ismissing(x) && x != "" && x != " -   ", hospitalbeds)
for i ∈ axes(hospitalbeds, 1)
    hospitalbeds.OrgCode[i] = replace(hospitalbeds.OrgCode[i], ' ' => "")
end
insertcols!(
    hospitalbeds,
    :TotalBeds => [ 
        parse(Int, replace(hospitalbeds.TotalBedsAvailable[i], ',' => "")) 
        for i ∈ axes(hospitalbeds, 1)
    ]
)

leftjoin!(hospitaldata, hospitalbeds; on= :TrustCode => :OrgCode )

insertcols!(
    hospitaldata,
    :VolumePerBed => Vector{Union{Missing, Float64}}(missing, size(hospitaldata, 1)),
    :ProportionSingleBeds => Vector{Union{Missing, Float64}}(missing, size(hospitaldata, 1)),
)

select!(hospitaldata, :TrustCode, :TrustType, :SiteCode, :SiteType, :TotalBeds, :HeatedVolume, :SingleBedsEnsuite, :VolumePerBed, :ProportionSingleBeds)

for trust ∈ unique(hospitaldata.TrustCode)
    totalsinglebeds = sum(hospitaldata.SingleBedsEnsuite .* (hospitaldata.TrustCode .== trust))
    totalvolume = sum(hospitaldata.HeatedVolume .* (hospitaldata.TrustCode .== trust))
    inds = findall(x -> x == trust, hospitaldata.TrustCode)
    for i ∈ inds 
        hospitaldata.VolumePerBed[i] = totalvolume / hospitaldata.TotalBeds[i]
        hospitaldata.ProportionSingleBeds[i] = min(
            totalsinglebeds / hospitaldata.TotalBeds[i],
            one(totalsinglebeds / hospitaldata.TotalBeds[i])
        )
    end
end

select!(hospitaldata, :TrustCode, :VolumePerBed, :ProportionSingleBeds)
unique!(hospitaldata)

coviddata = CSV.read(datadir("exp_raw", "dataset.csv"), DataFrame)

# make hospital codes categorical variables 
rename!(coviddata, :Codes => "StringCodes")
insertcols!(coviddata, 1, :Code => CategoricalArray(coviddata.StringCodes))

# calculate values for the analysis
insertproportions!(
    coviddata; 
    I=:StaffAbsences, N=:StaffTotal, Y=:CovidBeds, M=:AllBeds
)

# remove data where values are missing or NaN
for name ∈ names(coviddata)
    name ∈ [ "Code", "StringCodes", "Date" ] && continue
    filter!(name => x -> !ismissing(x) && !isnan(x), coviddata)
end

# remove outlying values where proportions >= 1 
for c ∈ [ :PatientsProportion, :StaffProportion ]
    for i ∈ axes(coviddata, 1)
        if ismissing(getproperty(coviddata, c)[i]) || getproperty(coviddata, c)[i] < 0
            getproperty(coviddata, c)[i] = 0  # most of these are filtered out later 
        elseif getproperty(coviddata, c)[i] > 1
            getproperty(coviddata, c)[i] = 1
        end
    end
    #filter!(c => x -> 0 <= x < 1, coviddata)
end

leftjoin!(coviddata, hospitaldata; on= :Code => :TrustCode )

communitydata = CSV.read(datadir("exp_raw", "OxCGRT_compact_subnational_v1.csv"), DataFrame)

#filter!(:CountryCode => x-> x == "GBR", communitydata)
filter!(:RegionCode => x-> x == "UK_ENG", communitydata)
insertcols!(
    communitydata, 
    :newcases => [ 
        i == 1 ? 
            0 : 
            max(communitydata.ConfirmedCases[i] - communitydata.ConfirmedCases[i-1], 0) 
        for i ∈ axes(communitydata, 1) 
    ],
    :FormattedDate => [ Date("$d", "yyyymmdd") for d ∈ communitydata.Date ]
)
insertcols!(communitydata, :weeklycases => [ i <= 7 ? 0 : maximum(@view communitydata.newcases[i-6:i]) for i ∈ axes(communitydata, 1) ])
select!(communitydata, :FormattedDate, :weeklycases, :StringencyIndex_Average)

leftjoin!(coviddata, communitydata; on= :Date => :FormattedDate )

filter!(:VolumePerBed => x -> !ismissing(x), coviddata)

# need to find the first and last date for each hospital
maxstartdate, minenddate = let 
    hospcodes = unique(coviddata.Code)
    startdate = zeros(Date, length(hospcodes))
    enddate = zeros(Date, length(hospcodes))
    for (i, c) ∈ enumerate(hospcodes)
        _tdf = filter(:Code => x -> x == c, coviddata)
        startdate[i] = minimum(_tdf.Date)
        enddate[i] = maximum(_tdf.Date)
    end
    ( maximum(startdate), maximum(enddate) )
end
# (Date("2020-04-04"), Date("2022-06-08"))

filter!(:Date => x -> Date("2020-04-04") <= x <= Date("2022-06-08"), coviddata)
insertcols!(coviddata, :t => Dates.value.(coviddata.Date - Date("2020-04-03")))

# remove hospitals with very little data 
removecodes = String7[ ]
for (i, c) ∈ enumerate(unique(coviddata.Code))
    if sum(coviddata.Code .== c) < 780 
        push!(removecodes, String(c))
    end
#    _tdf = filter(:Code => x -> x == c, coviddata)
 #   println(size(_tdf))
 #   patients[:, i] .= _tdf.PatientsProportion
  #  staff[:, i] .= _tdf.StaffProportion
end

filter!(:StringCodes => x -> x ∉ removecodes, coviddata)
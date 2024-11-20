
using DrWatson
@quickactivate :ImmuneBoostingHealthcare
using CategoricalArrays, CSV, DataFrames, Dates


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load UK data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

finaldata = let 
    hospitaldata = CSV.read(datadir("exp_raw", "SiteData.csv"), DataFrame)
    rename!(hospitaldata, Dict(Symbol("Trust Code") => "TrustCode"))
    rename!(hospitaldata, Dict(Symbol("Commissioning Region") => "CommissioningRegion"))
    rename!(hospitaldata, Dict(Symbol("Trust Type") => "TrustType"))
    rename!(hospitaldata, Dict(Symbol("Site Code") => "SiteCode"))
    rename!(hospitaldata, Dict(Symbol("Site Type") => "SiteType"))
    rename!(hospitaldata, Dict(Symbol("Site heated volume (m³)") => "HeatedVolumeString"))
    _sbsymbol = Symbol("Single bedrooms for patients with en-suite facilities (No.)")
    rename!(hospitaldata, Dict(_sbsymbol => "SingleBedsString"))
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

    hospitalbeds = CSV.read(
        datadir("exp_raw", "GeneralAcuteOccupiedBedsbyTrust.csv"), DataFrame
    )
    filter!(
        :TotalBedsAvailable => x -> !ismissing(x) && x != "" && x != " -   ", hospitalbeds
    )
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
        :ProportionSingleBeds => Vector{Union{Missing, Float64}}(
            missing, size(hospitaldata, 1)
        ),
    )

    select!(
        hospitaldata, 
        :TrustCode, :CommissioningRegion, :TrustType, :SiteCode, :SiteType, 
        :TotalBeds, :HeatedVolume, :SingleBedsEnsuite, :VolumePerBed, :ProportionSingleBeds
    )

    for trust ∈ unique(hospitaldata.TrustCode)
        totalsinglebeds = sum(
            hospitaldata.SingleBedsEnsuite .* (hospitaldata.TrustCode .== trust)
        )
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

    select!(
        hospitaldata, 
        :TrustCode, :CommissioningRegion, :VolumePerBed, :ProportionSingleBeds
    )
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
    end

    for i ∈ axes(coviddata, 1)
        c = :PatientsProportion
        if ismissing(getproperty(coviddata, c)[i]) || getproperty(coviddata, c)[i] < 0
            getproperty(coviddata, c)[i] = 0  # most of these are filtered out later 
        elseif getproperty(coviddata, c)[i] > 1
            getproperty(coviddata, c)[i] = 1
        end
        c = :StaffProportion
        if ismissing(getproperty(coviddata, c)[i]) || getproperty(coviddata, c)[i] < 0
            getproperty(coviddata, c)[i] = missing  # most of these are filtered out later 
        elseif getproperty(coviddata, c)[i] > 1
            getproperty(coviddata, c)[i] = missing
        end
    end

    coviddata = innerjoin(coviddata, hospitaldata; on= :Code => :TrustCode )

    communitydata = CSV.read(
        datadir("exp_raw", "OxCGRT_compact_subnational_v1.csv"), DataFrame
    )

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
    insertcols!(communitydata, :weeklycases => [ 
        i <= 7 ? 
            0 :
            maximum(@view communitydata.newcases[i-6:i]) 
        for i ∈ axes(communitydata, 1) 
    ])
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

    filter!(:Date => x -> Date("2020-04-04") <= x <= Date("2022-06-08"), coviddata)
    insertcols!(coviddata, :t => Dates.value.(coviddata.Date - Date("2020-03-19")))

    # remove hospitals with very little data 
    removecodes = String7[ ]
    for (i, c) ∈ enumerate(unique(coviddata.Code))
        if sum(coviddata.Code .== c) < 780 
            push!(removecodes, String(c))
        end
    end

    filter!(:StringCodes => x -> x ∉ removecodes, coviddata)

    # generate final dataset 

    finaldata = DataFrame(
        :Code => Any[ ],
        :StringCodes => String7[ ],
        :CommissioningRegion => String[ ],
        :Date => Date[ ],
        :CovidBeds => Int[ ],
        :AllBeds => Int[ ],
        :StaffAbsences => Int[ ],
        :StaffTotal => Int[ ],
        :CovidPatients => Float64[ ],
        :CovidAbsences => Float64[ ],
        :PatientsProportion => Float64[ ],
        :StaffProportion => Float64[ ],
        :VolumePerBed => Float64[ ],
        :ProportionSingleBeds => Float64[ ],
        :weeklycases => Int[ ],
        :StringencyIndex_Average => Float64[ ],
        :t => Float64[ ],
    )
    for c ∈ unique(coviddata.Code)
        ldf = DataFrame(
            :Code => c,
            :StringCodes => String(c),
            :Date => Date("2020-03-19"):Day(1):Date("2022-06-28")
        )
        leftjoin!(ldf, coviddata; on=[ :Code, :StringCodes, :Date ])
        maximum(describe(ldf).nmissing) == 832 && continue  # all missing
        for i ∈ axes(ldf, 1)
            for v ∈ [ 
                :CovidBeds, 
                :StaffAbsences, 
                :CovidPatients, 
                :CovidAbsences, 
                :PatientsProportion, 
                :StaffProportion, 
                :weeklycases, 
                :StringencyIndex_Average,
            ]
                if ismissing(getproperty(ldf, v)[i])
                    getproperty(ldf, v)[i] = 0 
                end
            end
            for v ∈ [ 
                :CommissioningRegion, 
                :AllBeds, 
                :StaffTotal, 
                :VolumePerBed, 
                :ProportionSingleBeds, 
            ]
                if ismissing(getproperty(ldf, v)[i])
                    getproperty(ldf, v)[i] = maximum(skipmissing(getproperty(ldf, v)))
                end
            end
            if ismissing(ldf.t[i])
                ldf.t[i] = Dates.value.(ldf.Date[i] - Date("2020-03-19"))
            end
        end
        append!(finaldata, ldf)
    end
    finaldata.Code = CategoricalArray(finaldata.Code)
        
    Dict("finaldata" => finaldata)
end

safesave(datadir("exp_pro", "finaldata.jld2"), finaldata)

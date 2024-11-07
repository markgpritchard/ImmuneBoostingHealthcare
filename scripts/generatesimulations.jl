

using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CategoricalArrays, CSV, DataFrames, Dates, Distributions, Random, StaticArrays
using StatsBase

#using BenchmarkTools

if isfile(datadir("exp_pro", "finaldata.jld2"))
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
else 
    include("loaddata.jl")
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
end

## Simulate community

simulations = let 
    function modeltransmission(t) 
        if t < 80 
            return 0.5 * (1 + 0.05 * cos(2π * t / 365))
        elseif t < 130 
            return 0.4
        else 
            return 0.475 * (1 + 0.05 * cos(2π * t / 365))
        end
    end
    
    function modelvaccination(t)
        if 264 <= t < 469  # 8 December 2020 to 1 July 2021
            return 0.004
        elseif 531 <= t <= 621  # 1 September to 30 November 2021
            return 0.004
        elseif t >= 712  # from 1 March 2021
            return 0.002
        else 
            return 0.0 
        end
    end
    
    Random.seed!(1729)
    
    communityvalues = let 
        u0_community = [ 55_999_900, 100, 0, 0, 0, 0, 0, 0 ]
        communityp = SEIIRRRSp(
            modeltransmission,  # infection rate 
            0.5,  # rate of leaving exposed compartments 
            0.2,  # rate of leaving infectious compartment
            0.0,  # strength of "force of boosting" relative to λ
            modelvaccination,  # vaccination rate
            0.01  # rate of immune waning 
        )
        stochasticseiirrrs(u0_community, 1:831, communityp)
    end
    
    global const COMMUNITYVALUES = deepcopy(communityvalues)
    
    simweeklycases = [
        (
            start = max(1, 7 * round(Int, (t - 0.5 ) / 7, RoundDown));
            lst = min(831, 7 * round(Int, (t - 0.5 ) / 7, RoundUp));
            sum(@view communityvalues[start:lst, 3:4])
        )
        for t ∈ 1:831
    ]
    
    function modelcommunitytransmissiontohcw(t)
        num = sum(@view COMMUNITYVALUES[round(Int, t, RoundUp), 3:4]) 
        denom = sum(@view COMMUNITYVALUES[round(Int, t, RoundUp), 1:7])
        return num * modeltransmission(t) / (1.5 * denom)
    end
    
    ## Simulate hospitals
    
    simdata_noboost = DataFrame(
        :Code => String7[ ],
        :StringCodes => String7[ ],
        :Date => Missing[ ],
        :CovidBeds => Int64[ ],
        :AllBeds => Int64[ ],
        :StaffAbsences => Int64[ ],
        :StaffTotal => Int64[ ],
        :CovidPatients => Float64[ ],
        :CovidAbsences => Float64[ ],
        :PatientsProportion => Float64[ ],
        :StaffProportion => Float64[ ],
        :VolumePerBed => Float64[ ],
        :ProportionSingleBeds => Float64[ ],
        :weeklycases => Int64[ ],
        :StringencyIndex_Average => Float64[ ],
        :t => Float64[ ],
        :betahh => Float64[ ],
        :betahp => Float64[ ],
        :betaph => Float64[ ],
        :betapp => Float64[ ],
    )
    
    for h ∈ 1:100
        Random.seed!(h)
        patienttotal = max(10, sample(finaldata.AllBeds)) 
        stafftotal = sample(finaldata.StaffTotal)  
        vpd = sample(finaldata.VolumePerBed) 
        psb = sample(finaldata.ProportionSingleBeds) 
        
        u0 = zeros(Int, 17)
        u0[1] = patienttotal 
        u0[7] = stafftotal
    
        p = WXYYZSEIIRRRSp(
            rand(Uniform(0.1, 0.3)) - vpd / 20000,  # rate of infection from healthcare worker to healthcare worker 
            rand(Uniform(0.05, 0.15)) - vpd / 40000,  # rate of infection from patient to healthcare worker  
            rand(Uniform(0.05, 0.15)) - vpd / 40000,  # rate of infection from healthcare worker to patient  
            rand(Uniform(0.05, 0.15)) - psb / 20,  # rate of infection from patient to patient  
            0.5,  # rate of leaving exposed compartments 
            0.2,  # rate of leaving infectious compartment
            0.0,  # strength of "force of boosting" relative to λ
            vaccinatestaff,  # healthcare worker vaccination rate
            0.01,  # rate of immune waning 
            rand(Uniform(0.2, 0.25)),  # discharge rate of non-infected 
            rand(Uniform(0.1, 0.17)),  # discharge rate of infected  
            modelcommunitytransmissiontohcw,  # community force of infection  
            2/7  # daily proportion diagnosed
        )
        
        @unpack hospitaldiagnoses, hospitalvalues = stochasticwxyyzseiirrrs(
            u0, 1.0:831.0, p, communityvalues
        )
    
        df = DataFrame(
            :Code => [ "$h" for _ ∈ 1:831 ],
            :StringCodes => [ "$h" for _ ∈ 1:831 ],
            :Date => [ missing for _ ∈ 1:831 ],
            :CovidBeds => [ sum(@view hospitalvalues[t, 3:5]) for t ∈ 1:831 ],
            :AllBeds => [ sum(@view hospitalvalues[t, 1:6]) for t ∈ 1:831 ],
            :StaffAbsences => hospitaldiagnoses,
            :StaffTotal => [ sum(@view hospitalvalues[t, 7:14]) for t ∈ 1:831 ],
            :CovidPatients => [ 
                sum(@view hospitalvalues[t, 3:5]) / sum(@view hospitalvalues[t, 1:6]) 
                for t ∈ 1:831 
            ],
            :CovidAbsences => [ 
                hospitaldiagnoses[t] / sum(@view hospitalvalues[t, 7:14]) 
                for t ∈ 1:831 
            ],
            :PatientsProportion => [ 
                sum(@view hospitalvalues[t, 3:5]) / sum(@view hospitalvalues[t, 1:6]) 
                for t ∈ 1:831 
            ],
            :StaffProportion => [ 
                hospitaldiagnoses[t] / sum(@view hospitalvalues[t, 7:14]) 
                for t ∈ 1:831 
            ],
            :VolumePerBed => [ vpd for _ ∈ 1:831 ],
            :ProportionSingleBeds => [ psb for _ ∈ 1:831 ],
            :weeklycases => simweeklycases,
            :StringencyIndex_Average => finaldata.StringencyIndex_Average[2:832],
            :t => 1.0:831.0,
            :betahh => [ betahh(p, t) for t ∈ 1:831 ],
            :betahp => [ betahp(p, t) for t ∈ 1:831 ],
            :betaph => [ betaph(p, t) for t ∈ 1:831 ],
            :betapp => [ betapp(p, t) for t ∈ 1:831 ],
        )
        append!(simdata_noboost, df)
    end
    
    simdata_boost = DataFrame(
        :Code => String7[ ],
        :StringCodes => String7[ ],
        :Date => Missing[ ],
        :CovidBeds => Int64[ ],
        :AllBeds => Int64[ ],
        :StaffAbsences => Int64[ ],
        :StaffTotal => Int64[ ],
        :CovidPatients => Float64[ ],
        :CovidAbsences => Float64[ ],
        :PatientsProportion => Float64[ ],
        :StaffProportion => Float64[ ],
        :VolumePerBed => Float64[ ],
        :ProportionSingleBeds => Float64[ ],
        :weeklycases => Int64[ ],
        :StringencyIndex_Average => Float64[ ],
        :t => Float64[ ],
        :betahh => Float64[ ],
        :betahp => Float64[ ],
        :betaph => Float64[ ],
        :betapp => Float64[ ],
    )
    
    for h ∈ 1:100
        Random.seed!(h)
        patienttotal = max(10, sample(finaldata.AllBeds)) 
        stafftotal = sample(finaldata.StaffTotal)  
        vpd = sample(finaldata.VolumePerBed) 
        psb = sample(finaldata.ProportionSingleBeds) 
        
        u0 = zeros(Int, 17)
        u0[1] = patienttotal 
        u0[7] = stafftotal
        
        p = WXYYZSEIIRRRSp(
            rand(Uniform(0.1, 0.3)) - vpd / 20000,  # rate of infection from healthcare worker to healthcare worker 
            rand(Uniform(0.05, 0.15)) - vpd / 40000,  # rate of infection from patient to healthcare worker  
            rand(Uniform(0.05, 0.15)) - vpd / 40000,  # rate of infection from healthcare worker to patient  
            rand(Uniform(0.05, 0.15)) - psb / 20,  # rate of infection from patient to patient  
            0.5,  # rate of leaving exposed compartments 
            0.2,  # rate of leaving infectious compartment
            2.0,  # strength of "force of boosting" relative to λ
            modelvaccination, # healthcare worker vaccination rate
            0.01,  # rate of immune waning 
            rand(Uniform(0.2, 0.25)),  # discharge rate of non-infected 
            rand(Uniform(0.1, 0.17)),  # discharge rate of infected  
            modelcommunitytransmissiontohcw,  # community force of infection  
            2/7  # daily proportion diagnosed
        )
        
        @unpack hospitaldiagnoses, hospitalvalues = stochasticwxyyzseiirrrs(
            u0, 1.0:831.0, p, communityvalues
        )
    
        df = DataFrame(
            :Code => [ "$h" for _ ∈ 1:831 ],
            :StringCodes => [ "$h" for _ ∈ 1:831 ],
            :Date => [ missing for _ ∈ 1:831 ],
            :CovidBeds => [ sum(@view hospitalvalues[t, 3:5]) for t ∈ 1:831 ],
            :AllBeds => [ sum(@view hospitalvalues[t, 1:6]) for t ∈ 1:831 ],
            :StaffAbsences => hospitaldiagnoses,
            :StaffTotal => [ sum(@view hospitalvalues[t, 7:14]) for t ∈ 1:831 ],
            :CovidPatients => [
                sum(@view hospitalvalues[t, 3:5]) / sum(@view hospitalvalues[t, 1:6]) 
                for t ∈ 1:831 
            ],
            :CovidAbsences => [ 
                hospitaldiagnoses[t] / sum(@view hospitalvalues[t, 7:14]) 
                for t ∈ 1:831 
            ],
            :PatientsProportion => [ 
                sum(@view hospitalvalues[t, 3:5]) / sum(@view hospitalvalues[t, 1:6]) 
                for t ∈ 1:831 
            ],
            :StaffProportion => [ 
                hospitaldiagnoses[t] / sum(@view hospitalvalues[t, 7:14]) 
                for t ∈ 1:831 
            ],
            :VolumePerBed => [ vpd for _ ∈ 1:831 ],
            :ProportionSingleBeds => [ psb for _ ∈ 1:831 ],
            :weeklycases => simweeklycases,
            :StringencyIndex_Average => finaldata.StringencyIndex_Average[2:832],
            :t => 1.0:831.0,
            :betahh => [ betahh(p, t) for t ∈ 1:831 ],
            :betahp => [ betahp(p, t) for t ∈ 1:831 ],
            :betaph => [ betaph(p, t) for t ∈ 1:831 ],
            :betapp => [ betapp(p, t) for t ∈ 1:831 ],
        )
        append!(simdata_boost, df)
    end

    Dict( 
        "boostedsimulation" => simdata_boost, 
        "unboostedsimulation" => simdata_noboost,
        "communitycases" => communityvalues, 
        "simweeklycases" => simweeklycases,
    )
end

safesave(datadir("sims", "simulations.jld2"), simulations)

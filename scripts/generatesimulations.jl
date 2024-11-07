

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

filter!(:CommissioningRegion => x -> x == "MIDLANDS COMMISSIONING REGION", finaldata)

## Simulate number of isolating healthcare workers in each hospital

### Scenario with no immune boosting

Random.seed!(1729)

unboostedsimulation = let 
    sim = deepcopy(finaldata)
    bc = deepcopy(sim)  # betacoefficients
    select!(bc, :StringCodes, :VolumePerBed, :ProportionSingleBeds) 
    unique!(bc)
    insertcols!(
        bc,
        :betahh => Vector{Float64}(undef, size(bc, 1)),
        :betahp => Vector{Float64}(undef, size(bc, 1)),
    )

    for i ∈ axes(bc, 1)
        bc.betahh[i] = rand(
            truncated(
                Normal(
                    1 - 0.0004 * bc.VolumePerBed[i] - 0.05 * bc.ProportionSingleBeds[i], 
                    1
                ), 
                0.0, 
                10.0
            )
        ) / 7
        bc.betahp[i] = rand(
            truncated(
                Normal(
                    1.4 - 0.001 * bc.VolumePerBed[i] - 0.25 * bc.ProportionSingleBeds[i], 
                    1
                ), 
                0.0, 
                10.0
            )
        ) / 7
    end

    select!(bc, :StringCodes, :betahh, :betahp) 
    leftjoin!(sim, bc; on=:StringCodes)
    insertcols!(
        sim, 
        :λc => [ 
            (0.23 + 0.0006 * (100 - s)) * c / 56_550_138 
            for (s, c) ∈ zip(finaldata.StringencyIndex_Average, finaldata.weeklycases) 
        ],
        :λp => sim.betahp .* sim.PatientsProportion,
        :λh => Vector{Float64}(undef, size(sim, 1)),
        :StaffNewAbsences => Vector{Int}(undef, size(sim, 1)),
        :StaffSusceptible => Vector{Int}(undef, size(sim, 1)),
        :UndiagnosedStaffInfected => Vector{Int}(undef, size(sim, 1)),
        :StaffR1 => Vector{Int}(undef, size(sim, 1)),
        :StaffR2 => Vector{Int}(undef, size(sim, 1)),
        :StaffR3 => Vector{Int}(undef, size(sim, 1)),
    )

    for i ∈ axes(sim, 1)
        if i == 1 || sim.Code[i] != sim.Code[i-1]
            sim.λh[i] = 0.0 
            sim.StaffAbsences[i] = 0 
            sim.CovidAbsences[i] = 0.0 
            sim.StaffProportion[i] = 0.0 
            sim.StaffSusceptible[i] = sim.StaffTotal[i]
            sim.UndiagnosedStaffInfected[i] = 0 
            sim.StaffR1[i] = 0
            sim.StaffR2[i] = 0
            sim.StaffR3[i] = 0
        else
            sim.λh[i] = sim.StaffProportion[i-1] * sim.betahh[i]
            i_s = rand(
                truncated(
                    Poisson((sim.λc[i] + sim.λh[i] + sim.λp[i]) * sim.StaffSusceptible[i-1]), 
                    0, 
                    sim.StaffSusceptible[i-1]
                )
            )
            d_i = rand(
                truncated(
                    Poisson(0.1 * sim.UndiagnosedStaffInfected[i-1]), 
                    0,
                    sim.UndiagnosedStaffInfected[i-1]
                )
            )
            r_i = rand(
                truncated(
                    Poisson(0.15 * (sim.UndiagnosedStaffInfected[i-1] - d_i)), 
                    0, 
                    sim.UndiagnosedStaffInfected[i-1] - d_i
                )
            )
            v_s = rand(
                truncated(
                    Poisson(vaccinatestaff(sim.t[i-1]) * (sim.StaffSusceptible[i-1] - i_s)), 
                    0, 
                    sim.StaffSusceptible[i-1] - i_s
                )
            )
            b_2 = 0 
            b_3 = 0
            v_r2 = rand(
                truncated(
                    Poisson(vaccinatestaff(sim.t[i-1]) * sim.StaffR2[i-1]), 
                    0, 
                    sim.StaffR2[i-1] - b_2
                )
            )
            v_r3 = rand(
                truncated(
                    Poisson(vaccinatestaff(sim.t[i-1]) * sim.StaffR3[i-1]), 
                    0, 
                    sim.StaffR3[i-1] - b_3
                )
            )
            w_1 = rand(truncated(Poisson(0.01 * sim.StaffR1[i-1]), 0, sim.StaffR1[i-1]))
            w_2 = rand(
                truncated(
                    Poisson(0.01 * (sim.StaffR2[i-1] - v_r2)), 
                    0, 
                    sim.StaffR2[i-1] - v_r2 - b_2
                )
            ) 
            w_3 = rand(
                truncated(
                    Poisson(0.01 * (sim.StaffR3[i-1] - v_r3)), 
                    0, 
                    sim.StaffR3[i-1] - v_r3 - b_3
                )
            ) 

            sim.StaffNewAbsences[i] = d_i
            if i <= 10 || sim.Code[i] != sim.Code[i-10]
                sim.StaffAbsences[i] = sim.StaffAbsences[i-1] + sim.StaffNewAbsences[i]
                sim.StaffR1[i] = sim.StaffR1[i-1] - w_1 + v_s + v_r2 + v_r3 + r_i + b_2 + b_3
            else
                sim.StaffAbsences[i] = sum(@view sim.StaffNewAbsences[i-9:i])
                sim.StaffR1[i] = (
                    sim.StaffR1[i-1] - 
                    w_1 + 
                    sim.StaffNewAbsences[i-10] + 
                    v_s + v_r2 + v_r3 + r_i + b_2 + b_3
                )
            end
            sim.CovidAbsences[i] = sim.StaffAbsences[i] / sim.StaffTotal[i]
            sim.StaffProportion[i] = sim.StaffAbsences[i] / sim.StaffTotal[i]
            sim.StaffSusceptible[i] = sim.StaffSusceptible[i-1] - v_s - i_s + w_3
            sim.UndiagnosedStaffInfected[i] = sim.UndiagnosedStaffInfected[i-1] - d_i - r_i + i_s
            sim.StaffR2[i] = sim.StaffR2[i-1] - w_2 - v_r2 - b_2 + w_1
            sim.StaffR3[i] = sim.StaffR3[i-1] - w_3 - v_r3 - b_3  + w_2
        end
    end
    Dict("unboostedsimulation" => sim)
end

safesave(datadir("sims", "unboostedsimulation.jld2"), unboostedsimulation)


### Scenario with immune boosting

boostedsimulation = let  # psi = 2 
    ψ = 2
    sim = deepcopy(finaldata)
    bc = deepcopy(sim)  # betacoefficients
    select!(bc, :StringCodes, :VolumePerBed, :ProportionSingleBeds) 
    unique!(bc)
    insertcols!(
        bc,
        :betahh => Vector{Float64}(undef, size(bc, 1)),
        :betahp => Vector{Float64}(undef, size(bc, 1)),
    )

    for i ∈ axes(bc, 1)
        bc.betahh[i] = rand(
            truncated(
                Normal(
                    1 - 0.0004 * bc.VolumePerBed[i] - 0.05 * bc.ProportionSingleBeds[i], 
                    1
                ), 
                0.0, 
                10.0
            )
        ) / 7
        bc.betahp[i] = rand(
            truncated(
                Normal(
                    1.4 - 0.001 * bc.VolumePerBed[i] - 0.25 * bc.ProportionSingleBeds[i], 
                    1
                ), 
                0.0, 
                10.0
            )
        ) / 7
    end

    select!(bc, :StringCodes, :betahh, :betahp) 
    leftjoin!(sim, bc; on=:StringCodes)
    insertcols!(
        sim, 
        :λc => [ 
            (0.23 + 0.0006 * (100 - s)) * c / 56_550_138 
            for (s, c) ∈ zip(finaldata.StringencyIndex_Average, finaldata.weeklycases) 
        ],
        :λp => sim.betahp .* sim.PatientsProportion,
        :λh => Vector{Float64}(undef, size(sim, 1)),
        :StaffNewAbsences => Vector{Int}(undef, size(sim, 1)),
        :StaffSusceptible => Vector{Int}(undef, size(sim, 1)),
        :UndiagnosedStaffInfected => Vector{Int}(undef, size(sim, 1)),
        :StaffR1 => Vector{Int}(undef, size(sim, 1)),
        :StaffR2 => Vector{Int}(undef, size(sim, 1)),
        :StaffR3 => Vector{Int}(undef, size(sim, 1)),
    )

    for i ∈ axes(sim, 1)
        if i == 1 || sim.Code[i] != sim.Code[i-1]
            sim.λh[i] = 0.0 
            sim.StaffAbsences[i] = 0 
            sim.CovidAbsences[i] = 0.0 
            sim.StaffProportion[i] = 0.0 
            sim.StaffSusceptible[i] = sim.StaffTotal[i]
            sim.UndiagnosedStaffInfected[i] = 0 
            sim.StaffR1[i] = 0
            sim.StaffR2[i] = 0
            sim.StaffR3[i] = 0
        else
            sim.λh[i] = sim.StaffProportion[i-1] * sim.betahh[i]
            i_s = rand(
                truncated(
                    Poisson((sim.λc[i] + sim.λh[i] + sim.λp[i]) * sim.StaffSusceptible[i-1]), 
                    0, 
                    sim.StaffSusceptible[i-1]
                )
            )
            d_i = rand(
                truncated(
                    Poisson(0.1 * sim.UndiagnosedStaffInfected[i-1]), 
                    0,
                    sim.UndiagnosedStaffInfected[i-1]
                )
            )
            r_i = rand(
                truncated(
                    Poisson(0.15 * (sim.UndiagnosedStaffInfected[i-1] - d_i)), 
                    0, 
                    sim.UndiagnosedStaffInfected[i-1] - d_i
                )
            )
            v_s = rand(
                truncated(
                    Poisson(vaccinatestaff(sim.t[i-1]) * (sim.StaffSusceptible[i-1] - i_s)), 
                    0, 
                    sim.StaffSusceptible[i-1] - i_s
                )
            )
            b_2 = rand(
                truncated(
                    Poisson((sim.λc[i] + sim.λh[i] + sim.λp[i]) * ψ * sim.StaffR2[i-1]), 
                    0, 
                    sim.StaffR2[i-1]
                )
            )
            b_3 = rand(
                truncated(
                    Poisson((sim.λc[i] + sim.λh[i] + sim.λp[i]) * ψ * sim.StaffR3[i-1]), 
                    0, 
                    sim.StaffR3[i-1]
                )
            )
            v_r2 = rand(
                truncated(
                    Poisson(vaccinatestaff(sim.t[i-1]) * sim.StaffR2[i-1]), 
                    0, 
                    sim.StaffR2[i-1] - b_2
                )
            )
            v_r3 = rand(
                truncated(
                    Poisson(vaccinatestaff(sim.t[i-1]) * sim.StaffR3[i-1]), 
                    0, 
                    sim.StaffR3[i-1] - b_3
                )
            )
            w_1 = rand(truncated(Poisson(0.01 * sim.StaffR1[i-1]), 0, sim.StaffR1[i-1]))
            w_2 = rand(
                truncated(
                    Poisson(0.01 * (sim.StaffR2[i-1] - v_r2)), 
                    0, 
                    sim.StaffR2[i-1] - v_r2 - b_2
                )
            ) 
            w_3 = rand(
                truncated(
                    Poisson(0.01 * (sim.StaffR3[i-1] - v_r3)), 
                    0, 
                    sim.StaffR3[i-1] - v_r3 - b_3
                )
            ) 

            sim.StaffNewAbsences[i] = d_i
            if i <= 10 || sim.Code[i] != sim.Code[i-10]
                sim.StaffAbsences[i] = sim.StaffAbsences[i-1] + sim.StaffNewAbsences[i]
                sim.StaffR1[i] = sim.StaffR1[i-1] - w_1 + v_s + v_r2 + v_r3 + r_i + b_2 + b_3
            else
                sim.StaffAbsences[i] = sum(@view sim.StaffNewAbsences[i-9:i])
                sim.StaffR1[i] = (
                    sim.StaffR1[i-1] - 
                    w_1 + 
                    sim.StaffNewAbsences[i-10] + 
                    v_s + v_r2 + v_r3 + r_i + b_2 + b_3
                )
            end
            sim.CovidAbsences[i] = sim.StaffAbsences[i] / sim.StaffTotal[i]
            sim.StaffProportion[i] = sim.StaffAbsences[i] / sim.StaffTotal[i]
            sim.StaffSusceptible[i] = sim.StaffSusceptible[i-1] - v_s - i_s + w_3
            sim.UndiagnosedStaffInfected[i] = sim.UndiagnosedStaffInfected[i-1] - d_i - r_i + i_s
            sim.StaffR2[i] = sim.StaffR2[i-1] - w_2 - v_r2 - b_2 + w_1
            sim.StaffR3[i] = sim.StaffR3[i-1] - w_3 - v_r3 - b_3  + w_2
        end
    end
    Dict("boostedsimulation" => sim)
end

safesave(datadir("sims", "boostedsimulation.jld2"), boostedsimulation)


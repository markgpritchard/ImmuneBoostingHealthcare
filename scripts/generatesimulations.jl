
include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare

using DifferentialEquations, Distributions, Random, StaticArrays

## Simulate community

const COMMUNITYSOL = let
    u0_community = seiirrrs_u0(; S=55_990_000, E=10_000)

    p0_community = SEIIRRRSp( ; 
        beta=0.4, 
        gamma=0.2, 
        epsilon=0.0,
        rho=0.5, 
        omega=0.01
    ) 
    
    communitylockdown!(integrator) = integrator.p = modifyp(integrator.p; beta=0.2)
    communitylockdowncb = PresetTimeCallback(
        80, communitylockdown!; 
        save_positions=( false, false )
    )
    communityendlockdown!(integrator) = integrator.p = modifyp(integrator.p; beta=0.3)
    communityendlockdowncb = PresetTimeCallback(
        130, communityendlockdown!; 
        save_positions=( false, false )
    )
    
    communitycbs = CallbackSet(communitylockdowncb, communityendlockdowncb)
    
    communityprob = ODEProblem(seiirrrs!, u0_community, ( 0.0, 800.0 ), p0_community)
    communitysol = solve(communityprob, Vern9(; lazy=false); callback=communitycbs, saveat=10)

    communitysol
end

const COMMUNITYCASES = [ sum(@view COMMUNITYSOL[i][3:4]) for i ∈ eachindex(COMMUNITYSOL) ]

simulations = let
    function hospitalaffect!(integrator)
        i = round(Int, integrator.t / 10) + 1
        n = sum(@view integrator.u[1:6])
        adjustv = 450 / n
        alpha1 = 0.2 * adjustv * COMMUNITYSOL[i][1] / sum(@view COMMUNITYSOL[i][1:8])
        alpha2 = 0.2 * adjustv * COMMUNITYSOL[i][2] / sum(@view COMMUNITYSOL[i][1:8])
        alpha3 = 0.2 * adjustv * sum(@view COMMUNITYSOL[i][3:4]) / sum(@view COMMUNITYSOL[i][1:8])
        alpha4 = 0.2 * adjustv * sum(@view COMMUNITYSOL[i][5:8]) / sum(@view COMMUNITYSOL[i][1:8])
        lambdac = COMMUNITYSOL[i][11] / 10
        integrator.p = modifyp(
            integrator.p; 
            alpha=SA[ alpha1, alpha2, alpha3, alpha4 ], lambdac
        )
    end
    hospitalaffecttimes = collect(10:10:800)
    cb = PresetTimeCallback(hospitalaffecttimes, hospitalaffect!; save_positions=( false, false ))
    
    function vaccinate!(integrator) 
        vaccn = 0.0
        for i ∈ [ 7, 11, 12, 13 ]
            v = 0.9 * integrator.u[i]
            integrator.u[i] += -v
            vaccn += v
        end
        integrator.u[14] += vaccn 
    end
    vcb = PresetTimeCallback(300, vaccinate!; save_positions=( false, false ))
    
    cbs = CallbackSet(cb, vcb)
    
    # beta parameters for each hospital -- kept constant across different boosting conditions 
    Random.seed!(1729)
    vpds = [ rand(truncated(Normal(510, 290_000), 0, 824_000)) for _ ∈ 1:135 ]
    psbs = [ rand(truncated(Normal(0.221, 0.0205), 0, 1)) for _ ∈ 1:135 ]
    
    allbetas = [ 
        [ rand(truncated(Normal(x - rand(Beta(4, 6)) * vpds[i] / 220_000 - rand(Beta(4, 6)) * psbs[i], 0.1), 0, 1)) for x ∈ [ 0.8, 0.8, 0.6, 0.405 ] ] 
        for i ∈ 1:135 
    ]
    
    # Simulate 135 hospitals with no immune boosting.
    
    u0 = wxyyzseiirrrs_u0(; W=450, S=900)
    
    unboostedp0 = WXYYZSEIIRRRSp( ; 
        alpha=SA[ 0.2, 0.0, 0.0, 0.0 ], 
        beta=SA[ 0.4, 0.4, 0.2, 0.05 ],  # (βhh, βhp, βph, βpp)
        gamma=0.2, 
        delta=SA[ 0.2, 0.1], 
        epsilon=0.0,
        lambdac=0.01, 
        rho=0.5, 
        omega=0.01
    ) 
    
    unboostedsimulation = simulatehospitals(
        135, u0, 800, unboostedp0; 
        abstol=1e-15, callback=cbs, maxiters=5e4, allbetas, saveat=1
    )
    
    boostedp0 = WXYYZSEIIRRRSp( ; 
        alpha=SA[ 0.2, 0.0, 0.0, 0.0 ], 
        beta=SA[ 0.4, 0.2, 0.2, 0.05 ],  # (βhh, βhp, βph, βpp)
        gamma=0.2, 
        delta=SA[ 0.2, 0.1], 
        epsilon=1.0,
        lambdac=0.01, 
        rho=0.5, 
        omega=0.01
    ) 
    
    boostedsimulation = simulatehospitals(
        135, u0, 800, boostedp0; 
        abstol=1e-15, callback=cbs, maxiters=5e4, allbetas, saveat=1
    )
    
    for sim ∈ [ unboostedsimulation, boostedsimulation ]
        insertcols!(
            sim,
            :CommunityCases => [ COMMUNITYCASES[sim.t[i] >= 800 ? 80 : round(Int, sim.t[i] / 10, RoundDown)+1] for i ∈ axes(sim, 1) ],
            :VolumePerBed => [ vpds[levelcode(sim.Code[i])] for i ∈ axes(sim, 1) ],
            :ProportionSingleBeds => [ psbs[levelcode(sim.Code[i])] for i ∈ axes(sim, 1) ],
            :DiagnosedY => [ rand(Binomial(round(Int, y), 0.9)) for y ∈ sim.Y ],
            :DiagnosedI => [ rand(Binomial(round(Int, y), 0.9)) for y ∈ sim.I ],
        )
        insertcols!(
            sim,
            :PatientsProportion => sim.DiagnosedY ./ sim.M,
            :StaffProportion => sim.DiagnosedI ./ sim.N,
        )
    end

    Dict( 
        "boostedsimulation" => boostedsimulation, 
        "unboostedsimulation" => unboostedsimulation 
    )
end

safesave(datadir("sims", "simulations.jld2"), simulations)

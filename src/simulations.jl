
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model with force of infection from patients, other staff and community
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function seiirrrs!(du, u, p, t)
    # Proportion infectious
    p_i = sum(@view u[3:4]) / sum(@view u[1:8])  
    
    # Force of infection 
    λ = p.β[1] * p_i

    # Equations 
    du[1] = -λ * u[1] + 3 * p.ω * u[7]                                  # dS
    du[2] = λ * u[1] - p.rho * u[2]                                     # dE
    du[3] = p.rho * u[2] - 2 * p.γ * u[3]                               # dI1
    du[4] = 2 * p.γ * (u[3] - u[4])                                     # dI2
    du[5] = 2 * p.γ * u[4] - (3 * p.ω + p.ε * λ) * u[5] + p.rho * u[8]  # dR1
    du[6:7] = [                                                         # dR2:dR3
        3 * p.ω * u[i-1] - (3 * p.ω + p.ε * λ) * u[i] 
        for i ∈ 6:7 
    ]
    du[8] = p.ε * λ * sum(u[5:7]) - p.rho * u[8]                        # dE*
    du[9] = λ * u[1]  # cumulative number of infections 
    du[10] = 0.0  # median time in resistant compartment (filled by a Callback)
    u[11] = λ  # force of infection (not a differential equation)
end

function wxyyzseiirrrs!(du, u, p, t)
    np = sum(@view u[1:6])
    # Proportions infectious
    pp = sum(@view u[3:4]) / np  # proportion of patients infectious
    ph = sum(@view u[9:10]) / sum(@view u[7:14])  # proportion healthcare workers infectious
    
    # Force of infection 
    λh = p.β[1] * ph + p.β[2] * pp + p.λc
    λp = p.β[3] * ph + p.β[4] * pp

    # Equations 
    du[1] = p.α[1] * np - (p.δ[1] + λp) * u[1]                                       # dW
    du[2] = p.α[2] * np + λp * u[1] - (p.rho + p.δ[1]) * u[2]                        # dX
    du[3] = p.α[3] * np / 2 + p.rho * u[2] - (2 * p.γ + p.δ[2]) * u[3]               # dY1
    du[4] = p.α[3] * np / 2 + 2 * p.γ * u[3] - (2 * p.γ + p.δ[2]) * u[4]             # dY2
    du[5] = 2 * p.γ * u[4] - p.δ[1] * u[5]  # (recovered in hospital)           # dZ1
    du[6] = p.α[4] * np - p.δ[1] * u[6]  # (admitted immune)                         # dZ2
    du[7] = -λh * u[7] + 3 * p.ω * u[13]                                        # dS
    du[8] = λh * u[7] - p.rho * u[8]                                            # dE
    du[9] = p.rho * u[8] - 2 * p.γ * u[9]                                       # dI1
    du[10] = 2 * p.γ * (u[9] - u[10])                                           # dI2
    du[11] = 2 * p.γ * u[10] - (3 * p.ω + p.ε * λh) * u[11] + p.rho * u[14]     # dR1    
    du[12:13] = [                                                               # dR2:dR3
        3 * p.ω * u[i-1] - (3 * p.ω + p.ε * λh) * u[i] 
        for i ∈ 12:13 
    ]
    du[14] = p.ε * λh * sum(u[11:13]) - p.rho * u[14]                           # dE*
    du[15] = λh * u[7]  # cumulative number of infections 
    du[16] = 0.0  # median time in resistant compartment (filled by a Callback)
    u[17] = λh  # force of infection (not a differential equation)
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Constant and discrete forces of infection 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function simulatehospitals(hospitals, args...; kwargs...) 
    df = simulatehospital(args...; hospitalid=1, kwargs...)
    for i ∈ 2:hospitals 
        append!(df, simulatehospital(args...; hospitalid=i, kwargs...))
    end 
    insertcols!(df, 1, :Code => CategoricalArray(df.IntCode))
    return df 
end 

function simulatehospital(
    u0, duration, p0; 
    hospitalid, 
    infectionrate=true, lambdacvalues=nothing, changelambdafrequency=nothing, 
    allbetas=nothing, z1=:resistant, kwargs...
)
    p0 = _updatebetas(p0, allbetas, hospitalid)
    lambdac_callbacks = _changelambdac(p0.λc, lambdacvalues, changelambdafrequency, duration) 
    sol = _simulatehospital(u0, duration, p0, lambdac_callbacks; kwargs...)
    df = _sol_data(sol; infectionrate, z1)
    insertcols!(df, :IntCode => hospitalid)
    return df
end

function _simulatehospital(u0, duration::T, p0, callbackset; kwargs...) where T <: Real
    tspan = ( zero(T), duration )
    return _simulatehospital(u0, tspan, p0, callbackset; kwargs...)
end

function _simulatehospital(u0, tspan::Tuple, p0, callbackset; kwargs...)
    @assert minimum(u0) >= 0 "Cannot have negative numbers in compartments"
    @assert minimum(p0) >= 0 "Cannot have negative parameters"
    prob = ODEProblem(wxyyzseiirrrs!, u0, tspan, p0)
    return _simulatehospital(prob, callbackset; kwargs...)
end

function _simulatehospital(prob, callbackset::CallbackSet; kwargs...)
    return _simulatehospital(prob; callback=callbackset, kwargs...)
end

_simulatehospital(prob, ::Nothing; kwargs...) = _simulatehospital(prob; kwargs...)

function _simulatehospital(prob; alg=Vern9(; lazy=false), saveat=1, kwargs...)
    return solve(prob, alg; saveat, kwargs...)
end

_changelambdac(::Any, ::Nothing, ::Nothing, ::Any) = nothing

function _changelambdac(initiallambda, lambdacvalues, changelambdafrequency, duration) 
    # When does lambda change 
    changetimes = _changelambdac_changetimes(changelambdafrequency, duration)
    # What does it change to 
    changevalues = _changelambdac_changevalues(initiallambda, lambdacvalues, changetimes)
    changecallbacks = [ 
        _changelambdac_callback(changetimes[i], changevalues[i]) 
        for i ∈ axes(changetimes, 1) 
    ]
    return CallbackSet(changecallbacks...)
end 

function _changelambdac_changetimes(changelambdafrequency, duration; theta=1)
    changetimes = Float64[ ]
    t = 0.0 
    while t < duration 
        t += rand(Gamma(changelambdafrequency, theta))
        push!(changetimes, t)
    end 
    pop!(changetimes)
    return changetimes 
end 

function _changelambdac_changevalues(initiallambda, lambdacvalues, changetimes)
    lambdapossibilities = size(lambdacvalues, 1)
    if initiallambda ∈ lambdacvalues 
        v = findall(x -> x == initiallambda, lambdacvalues)[1]
    else 
        v = 1 
    end 
    changevalues = Float64[ ] 
    for _ ∈ axes(changetimes, 1)
        if v == lambdapossibilities 
            v = 1 
        else 
            v += 1 
        end 
        push!(changevalues, lambdacvalues[v])
    end 
    return changevalues
end 

function _changelambdac_callback(changetime, newlambda)
    affect!(integrator) = integrator.p = modifyp(integrator.p; lambdac=newlambda)
    return PresetTimeCallback(changetime, affect!; save_positions=( false, false ))
end 

function _simulateabsences(data) 
    return [ 
        data.t[i] <= 10 ? 
            data.cumulativeinfections[i] : 
            data.cumulativeinfections[i] - data.cumulativeinfections[i-10] 
        for i ∈ axes(data, 1) 
    ] 
end

# Simulate vaccination by moving susceptibles and resistants to E*
function _vaccination!(u, proportionvaccinated)
    @assert 0 <= proportionvaccinated <= 1
    vaccinated = proportionvaccinated * (u[7] + sum(u[11:13]))  # proportion of S and R1:R3
    u[14] += vaccinated  # to E*
    for i ∈ [ 7; collect(11:13) ]   # from S and R1:R3
        u[i] = (1 - proportionvaccinated) * u[i] 
    end 
end 

function _sol_data(sol; infectionrate, z1)
    labels = [ :W, :X, :Y, :Z, :M, :S, :E, :I, :R, :N, :cumulativeinfections, :λh ]
    if z1 == :resistant
        indices = [ 1, 2, 3:4, 5:6, 1:6, 7, 8, 9:10, 11:14, 7:14, 15, 17 ]
    else 
        indices = [ 1, 2, 3:5, 6, 1:6, 7, 8, 9:10, 11:14, 7:14, 15, 17 ]
    end 
    return _sol_data(sol, labels, indices; infectionrate)
end 

function _sol_data(sol, labels, indices; infectionrate) 
    data = DataFrame(t = sol.t)
    _insertdatavector!(data, labels, sol, indices) 
    if infectionrate
        insertcols!(
            data, 
            :infectionrate => _calculateinfectionrate(data.cumulativeinfections, size(sol, 2))
        )
    end 
    return data
end 

function _insertdatavector!(data, labels::AbstractVector{Symbol}, sol, indices) 
    for (i, label) ∈ enumerate(labels)
        _insertdatavector!(data, label, sol, indices[i])
    end 
end 

function _insertdatavector!(data, label::Symbol, sol, index)
    insertcols!(data, label => _datavector(sol, index))
end

_datavector(sol, index::Real) = [ sol.u[i][index] for i ∈ eachindex(sol) ]

function _datavector(sol, index::Vector{<:Real})
    return [ sum(@view sol.u[i][index]) for i ∈ eachindex(sol) ]
end

_datavector(sol, index::UnitRange) = _datavector(sol, collect(index)) 

function _calculateinfectionrate(cumulativeinfections, nrows) 
    infectionrate = deepcopy(cumulativeinfections)
    for i ∈ 2:nrows
        infectionrate[i] = cumulativeinfections[i] - cumulativeinfections[i-1]
    end
    return infectionrate
end 

function _makerandombetas(p0, randombetas::Bool)
    @assert randombetas = false 
    return p0 
end

function _makerandombetas(p0, randombetas::AbstractVector{<:Real})
    betas = [ rand(truncated(Normal(b, 0.1), 0, 1)) for b ∈ randombetas ]
    return modifyp(p0; beta=betas)
end

function _makerandombetas(p0, randombetas::AbstractVector{<:Distribution})
    betas = [ rand(b) for b ∈ randombetas ]
    return modifyp(p0; beta=betas)
end

_updatebetas(p0, ::Nothing, ::Any) = p0
_updatebetas(p0, allbetas, hospitalid) = modifyp(p0; beta=allbetas[hospitalid])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to make u0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function seiirrrs_u0(; S=0.0, E=0.0, I1=0.0, I2=0.0, R1=0.0, R2=0.0, R3=0.0, Estar=0.0)
    return seiirrrs_u0(S, E, I1, I2, R1, R2, R3, Estar)
end 

function seiirrrs_u0(S, E, I1, I2, R1, R2, R3, Estar)
    return seiirrrs_u0([ S, E, I1, I2, R1, R2, R3, Estar ])
end

seiirrrs_u0(compartments) = makeu0!(compartments, 11)

function wxyyzseiirrrs_u0( ; 
    W=0.0, X=0.0, Y1=0.0, Y2=0.0, Z1=0.0, Z2=0.0, 
    S=0.0, E=0.0, I1=0.0, I2=0.0, R1=0.0, R2=0.0, R3=0.0, Estar=0.0
)
    return wxyyzseiirrrs_u0(W, X, Y1, Y2, Z1, Z2, S, E, I1, I2, R1, R2, R3, Estar)
end 

function wxyyzseiirrrs_u0(W, X, Y1, Y2, Z1, Z2, S, E, I1, I2, R1, R2, R3, Estar)
    return wxyyzseiirrrs_u0([ W, X, Y1, Y2, Z1, Z2, S, E, I1, I2, R1, R2, R3, Estar ])
end

wxyyzseiirrrs_u0(compartments) = makeu0!(compartments, 17)

function makeu0!(compartments, n)
    # add compartments for cumulativeinfections, mediantime, λ
    makeu0!(compartments)
    @assert length(compartments) == n
    return compartments
end 

function makeu0!(compartments::Vector{T}) where T <: Real
    # add compartments for cumulativeinfections, mediantime, λ
    append!(compartments, zeros(T, 3)) 
    @assert minimum(compartments) >= 0 "Cannot initialize u0 with negative values"
    return compartments
end 

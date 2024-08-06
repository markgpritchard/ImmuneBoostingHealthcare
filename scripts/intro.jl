
using DrWatson
@quickactivate "ImmuneBoostingHealthcare"

using Bootstrap, CairoMakie, CSV, DataFrames, Random, Turing, Pigeons

###########################################################################################
# Functions 
###########################################################################################

using CategoricalArrays, CSV, DataFrames, Dates, DifferentialEquations, Distributions, GLM, StaticArrays
import Base: minimum

###########################################################################################
# Structs 
###########################################################################################


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters for ODE models 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abstract type AbstractParameters{T} end 

struct SEIIRRRSp{T} <: AbstractParameters{T} 
    β       :: T                # infection rate 
    γ       :: T                # rate of leaving infectious compartment
    ε       :: T                # strength of "force of boosting" relative to λ
    rho     :: T                # rate of leaving exposed compartments 
    ω       :: T                # rate of immune waning 
end

struct WXYYZSEIIRRRSp{T} <: AbstractParameters{T}
    α       :: SVector{4, T}    # rate of patients being admitted 
    β       :: SVector{4, T}    # rates of infections between patients and staff (βhh, βhp, βph, βpp)
    γ       :: T                # rate of leaving infectious compartment
    δ       :: SVector{2, T}    # rate of patients leaving hospital 
    ε       :: T                # strength of "force of boosting" relative to λ
    λc      :: T                # community force of infection  
    rho     :: T                # rate of leaving exposed compartments 
    ω       :: T                # rate of immune waning 
end

# produce parameters with keyword arguments or non-static arrays
SEIIRRRSp(; beta, gamma, epsilon, rho, omega) = SEIIRRRSp(beta, gamma, epsilon, rho, omega)

function WXYYZSEIIRRRSp(; alpha, beta, gamma, delta, epsilon, lambdac, rho, omega) 
    return WXYYZSEIIRRRSp(alpha, beta, gamma, delta, epsilon, lambdac, rho, omega)
end 

function WXYYZSEIIRRRSp(
    alpha::Vector{T}, beta, gamma::T, delta, epsilon::T, lambdac::T, rho::T, omega::T
) where T
    a = SVector{4}(alpha)
    return WXYYZSEIIRRRSp(a, beta, gamma, delta, epsilon, lambdac, rho, omega)
end 

function WXYYZSEIIRRRSp(
    alpha::SVector{4, T}, beta::Vector{T}, gamma::T, delta, 
    epsilon::T, lambdac::T, rho::T, omega::T
) where T
    b = SVector{4}(beta)
    return WXYYZSEIIRRRSp(alpha, b, gamma, delta, epsilon, lambdac, rho, omega)
end 

function WXYYZSEIIRRRSp(
    alpha::SVector{4, T}, beta::SVector{4, T}, gamma::T, delta::Vector{T}, 
    epsilon::T, lambdac::T, rho::T, omega::T
) where T
    d = SVector{2}(delta)
    return WXYYZSEIIRRRSp(alpha, beta, gamma, d, epsilon, lambdac, rho, omega)
end 

# parameter structs are non-mutating, so function to "modify" then by creating a new version 

function modifyp(
    p::SEIIRRRSp; 
    beta=nothing, gamma=nothing, epsilon=nothing, rho=nothing, omega=nothing
)
return SEIIRRRSp( 
    _modifyp(p, beta, :β), 
    _modifyp(p, gamma, :γ), 
    _modifyp(p, epsilon, :ε), 
    _modifyp(p, rho, :rho), 
    _modifyp(p, omega, :ω), 
)
end 

function modifyp(
    p::WXYYZSEIIRRRSp; 
    alpha=nothing, beta=nothing, gamma=nothing, delta=nothing, epsilon=nothing, 
    lambdac=nothing, rho=nothing, omega=nothing
)
return WXYYZSEIIRRRSp( 
    _modifyp(p, alpha, :α), 
    _modifyp(p, beta, :β), 
    _modifyp(p, gamma, :γ), 
    _modifyp(p, delta, :δ), 
    _modifyp(p, epsilon, :ε), 
    _modifyp(p, lambdac, :λc), 
    _modifyp(p, rho, :rho), 
    _modifyp(p, omega, :ω), 
)
end 

_modifyp(p, ::Nothing, symbol) = getfield(p, symbol)
_modifyp(::AbstractParameters{T}, v::T, ::Symbol) where T = v
_modifyp(::AbstractParameters{T}, v::AbstractVector{T}, ::Symbol) where T = v



###########################################################################################
# Simulations 
###########################################################################################

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

function simulationproportions!(data; I=:I, Y=:Y, M=:M, N=:N)
    insertcols!(data, :PatientsProportion => getproperty(data, Y) ./ getproperty(data, M))
    insertcols!(data, :StaffProportion => getproperty(data, I) ./ getproperty(data, N))
    #absences = _simulateabsences(data)
    #insertcols!(data, :AbsencesProportion => absences ./ data.N)
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

seiirrrs_u0(compartments) = makeu0(compartments, 11)

function wxyyzseiirrrs_u0( ; 
    W=0.0, X=0.0, Y1=0.0, Y2=0.0, Z1=0.0, Z2=0.0, 
    S=0.0, E=0.0, I1=0.0, I2=0.0, R1=0.0, R2=0.0, R3=0.0, Estar=0.0
)
    return wxyyzseiirrrs_u0(W, X, Y1, Y2, Z1, Z2, S, E, I1, I2, R1, R2, R3, Estar)
end 

function wxyyzseiirrrs_u0(W, X, Y1, Y2, Z1, Z2, S, E, I1, I2, R1, R2, R3, Estar)
    return wxyyzseiirrrs_u0([ W, X, Y1, Y2, Z1, Z2, S, E, I1, I2, R1, R2, R3, Estar ])
end

wxyyzseiirrrs_u0(compartments) = makeu0(compartments, 17)

function makeu0(compartments::Vector{T}, n) where T <: Real
    # add compartments for cumulativeinfections, mediantime, λ
    append!(compartments, zeros(T, 3)) 
    @assert minimum(compartments) >= 0 "Cannot initialize u0 with negative values"
    @assert length(compartments) == n
    return compartments
end 


###########################################################################################
# Regression 
###########################################################################################

function insertpreviousvalues!(data, column; kwargs...)
    insertpreviousvalues!(data, column, 1:40; kwargs...)
end

function insertpreviousvalues!(data, column, trange::Integer; kwargs...) 
    insertpreviousvalues!(data, column, 1:trange; kwargs...)
end

function insertpreviousvalues!(
    data, column::AbstractVector{<:Symbol}, trange::AbstractArray; 
    label=column
)
    for (i, c) ∈ enumerate(column) 
        insertpreviousvalues!(data, c, trange; label=label[i]) 
    end
end 

function insertpreviousvalues!(data, column::Symbol, trange::AbstractArray; kwargs...)
    for i ∈ trange
        insertpreviousvalues!(data, column, i, maximum(trange); kwargs...) 
    end
end 

function insertpreviousvalues!(data, column, t, maxt; label=column)
    insertcols!(data, "$(label)_t$t" => _previousvalues(data, column, maxt - t, maxt))
end 

function _previousvalues(data, column, deltas, maxt)
    d = data[:, column]
    return [ data.t[i] <= maxt ? zero(d[1]) : d[i-deltas] for i ∈ axes(data, 1) ]
end 

function calculategroupedvalues(data, prefix, numberrange::AbstractArray)
    tot = sum([ getproperty(data, Symbol("$prefix$i")) for i ∈ numberrange ])
    return tot / length(numberrange)
end

function calculategroupedvalues(data, prefix, numb::Integer)
    return getproperty(data, Symbol("$prefix$numb"))
end

function insertgroupedvalues!(data, range1, range2, range3, range4)
    insertcols!(
        data,
        :i4 => calculategroupedvalues(data, "I_t", range4),
        :i3 => calculategroupedvalues(data, "I_t", range3),
        :i2 => calculategroupedvalues(data, "I_t", range2),
        :i1 => calculategroupedvalues(data, "I_t", range1),
        :y4 => calculategroupedvalues(data, "Y_t", range4),
        :y3 => calculategroupedvalues(data, "Y_t", range3),
        :y2 => calculategroupedvalues(data, "Y_t", range2),
        :y1 => calculategroupedvalues(data, "Y_t", range1),
    )
end
#=
function calculatedirecteffect(data)
    Ei2nformula = @formula(i2 ~ t^0.5 + log(t) + Code)
    Ei2nregr = fit(LinearModel, Ei2nformula, data)
    Ei2n = predict(Ei2nregr)
    di2n = Ei2n .- data.i2
    sigman2 = var(di2n)
    
    Ei2dformula = @formula(i2 ~ y3 + i3 + i1 + y1 + t^0.5 + log(t) + Code)
    Ei2dregr = fit(LinearModel, Ei2dformula, data)
    Ei2d = predict(Ei2dregr)
    di2d = Ei2d .- data.i2
    sigmad2 = var(di2d)
    
    wn = @. exp(-(data.i2 - Ei2n)^2 / (2π * sigman2)) / (2π * sigman2)
    wd = @. exp(-(data.i2 - Ei2d)^2 / (2π * sigmad2)) / (2π * sigmad2)
    w_unlimited = wn ./ wd
    w98 = quantile(w_unlimited, 0.98)
    w = [ x > w98 ? w98 : x for x ∈ w_unlimited ]
    
    directeffectformula = @formula(i4 ~ y2 + t^0.5 + log(t) + Code)
    Ei2dregr = fit(LinearModel, directeffectformula, data; wts=w)
    return @ntuple Ei2dregr w wd wn w_unlimited w98
end
=#
#=
function calculatedirecteffect(data)
    Ey2nformula = @formula(y2 ~ t^0.5 + log(t) + Code)
    Ey2nregr = fit(LinearModel, Ey2nformula, data)
    Ey2n = predict(Ey2nregr)
    dy2n = Ey2n .- data.y2
    sigman2 = var(dy2n)
    
    Ey2dformula = @formula(y2 ~ y1 + i1 + y3 + i3 + t^0.5 + log(t) + Code)

    Ey2dregr = fit(LinearModel, Ey2dformula, data)
    Ey2d = predict(Ey2dregr)
    dy2d = Ey2d .- data.y2
    sigmad2 = var(dy2d)
    
    wn = @. exp(-(data.y2 - Ey2n)^2 / (2π * sigman2)) / (2π * sigman2)
    wd = @. exp(-(data.y2 - Ey2d)^2 / (2π * sigmad2)) / (2π * sigmad2)
    w_unlimited = wn ./ wd
    w98 = quantile(w_unlimited, 0.98)
    w = [ x > w98 ? w98 : x for x ∈ w_unlimited ]
    
    directeffectformula = @formula(i4 ~ y2 + t^0.5 + log(t) + Code)
    directeffectregr = fit(LinearModel, directeffectformula, data; wts=w)
    return @ntuple directeffectregr w wd wn w_unlimited w98
end

function calculatetotaleffect(data)
    Ey2nformula = @formula(y2 ~ t^0.5 + log(t) + Code)
    Ey2nregr = fit(LinearModel, Ey2nformula, data)
    Ey2n = predict(Ey2nregr)
    dy2n = Ey2n .- data.y2
    sigman2 = var(dy2n)
    
    Ey2dformula = @formula(y2 ~ y1 + i1 + t^0.5 + log(t) + Code)

    Ey2dregr = fit(LinearModel, Ey2dformula, data)
    Ey2d = predict(Ey2dregr)
    dy2d = Ey2d .- data.y2
    sigmad2 = var(dy2d)
    
    wn = @. exp(-(data.y2 - Ey2n)^2 / (2π * sigman2)) / (2π * sigman2)
    wd = @. exp(-(data.y2 - Ey2d)^2 / (2π * sigmad2)) / (2π * sigmad2)
    w_unlimited = wn ./ wd
    w98 = quantile(w_unlimited, 0.98)
    w = [ x > w98 ? w98 : x for x ∈ w_unlimited ]
    
    totaleffectformula = @formula(i4 ~ y2 + t^0.5 + log(t) + Code)
    totaleffectregr = fit(LinearModel, totaleffectformula, data; wts=w)
    return @ntuple totaleffectregr w wd wn w_unlimited w98
end
=#
function calculatedirecteffect(data)
    Ey2nformula = @formula(Y_t11 ~ t^0.5 + log(t) + Code)
    Ey2nregr = fit(LinearModel, Ey2nformula, data)
    Ey2n = predict(Ey2nregr)
    dy2n = Ey2n .- data.Y_t11
    sigman2 = var(dy2n)
    
    Ey2dformula = @formula(
        Y_t11 ~ 
        Y_t1 + Y_t2 + Y_t3 + Y_t4 + Y_t5 + Y_t6 + Y_t7 + Y_t8 + Y_t9 + Y_t10 + 
        Y_t12 + Y_t13 + Y_t14 + Y_t15 + Y_t16 + Y_t17 + Y_t18 + Y_t19 + Y_t20 + Y_t21 + 
        I_t1 + I_t2 + I_t3 + I_t4 + I_t5 + I_t6 + I_t7 + I_t8 + I_t9 + I_t10 + 
        I_t12 + I_t13 + I_t14 + I_t15 + I_t16 + I_t17 + I_t18 + I_t19 + I_t20 + I_t21 + 
        t^0.5 + log(t) + Code
    )

    Ey2dregr = fit(LinearModel, Ey2dformula, data)
    Ey2d = predict(Ey2dregr)
    dy2d = Ey2d .- data.Y_t11
    sigmad2 = var(dy2d)
    
    wn = @. exp(-(data.Y_t11 - Ey2n)^2 / (2π * sigman2)) / (2π * sigman2)
    wd = @. exp(-(data.Y_t11 - Ey2d)^2 / (2π * sigmad2)) / (2π * sigmad2)
    println("wn=$wn")

    println("wd=$wd")
    w_unlimited = wn ./ wd
    w98 = quantile(w_unlimited, 0.98)
    w = [ x > w98 ? w98 : x for x ∈ w_unlimited ]
    
    directeffectformula = @formula(I_t22 ~ Y_t11 + t^0.5 + log(t) + Code)
    directeffectregr = fit(LinearModel, directeffectformula, data; wts=w)
    return @ntuple directeffectregr w wd wn w_unlimited w98
end

function calculatetotaleffect(data)
    Ey2nformula = @formula(Y_t11 ~ t^0.5 + log(t) + Code)
    Ey2nregr = fit(LinearModel, Ey2nformula, data)
    Ey2n = predict(Ey2nregr)
    dy2n = Ey2n .- data.Y_t11
    sigman2 = var(dy2n)
    
    Ey2dformula = @formula(
        Y_t11 ~ 
        Y_t1 + Y_t2 + Y_t3 + Y_t4 + Y_t5 + Y_t6 + Y_t7 + Y_t8 + Y_t9 + Y_t10 + 
        I_t1 + I_t2 + I_t3 + I_t4 + I_t5 + I_t6 + I_t7 + I_t8 + I_t9 + I_t10 + 
        t^0.5 + log(t) + Code
    )

    Ey2dregr = fit(LinearModel, Ey2dformula, data)
    Ey2d = predict(Ey2dregr)
    dy2d = Ey2d .- data.Y_t11
    sigmad2 = var(dy2d)
    
    wn = @. exp(-(data.Y_t11 - Ey2n)^2 / (2π * sigman2)) / (2π * sigman2)
    wd = @. exp(-(data.Y_t11 - Ey2d)^2 / (2π * sigmad2)) / (2π * sigmad2)
    w_unlimited = wn ./ wd
    w98 = quantile(w_unlimited, 0.98)
    w = [ x > w98 ? w98 : x for x ∈ w_unlimited ]
    
    totaleffectformula = @formula(I_t22 ~ Y_t11 + t^0.5 + log(t) + Code)
    totaleffectregr = fit(LinearModel, totaleffectformula, data; wts=w)
    return @ntuple totaleffectregr w wd wn w_unlimited w98
end

function pointestimatedirecteffect(data)
    @unpack directeffectregr = calculatedirecteffect(data)
    return coef(directeffectregr)[2]
end

function pointestimatetotaleffect(data)
    @unpack totaleffectregr = calculatetotaleffect(data)
    return coef(totaleffectregr)[2]
end

###########################################################################################
# Extra functions 
###########################################################################################

minimum(p::SEIIRRRSp) = minimum([ p.β, p.γ, p.ε, p.rho, p.ω ])

function minimum(p::WXYYZSEIIRRRSp)
    a = minimum(p.α)
    b = minimum(p.β)
    d = minimum(p.δ)
    return minimum([ a, b, p.γ, d, p.ε, p.λc, p.rho, p.ω ])
end 

###########################################################################################
# scripts 
###########################################################################################

# Simulate community

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

const COMMUNITYSOL = deepcopy(communitysol)

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
cb = PresetTimeCallback(hospitalaffecttimes, hospitalaffect!)

function vaccinate!(integrator) 
    vaccn = 0.0
    for i ∈ [ 7, 11, 12, 13 ]
        v = 0.9 * integrator.u[i]
        integrator.u[i] += -v
        vaccn += v
    end
    integrator.u[14] += vaccn 
end
vcb = PresetTimeCallback(300, vaccinate!)

cbs = CallbackSet(cb, vcb)

# beta parameters for each hospital -- kept constant across different boosting conditions 
Random.seed!(1729)
allbetas = [ 
    [ rand(truncated(Normal(x, 0.1), 0, 1)) for x ∈ [ 0.4, 0.4, 0.2, 0.005 ] ] 
    for _ ∈ 1:300 
]

# Simulate 300 hospitals with no immune boosting.

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
    300, u0, 800, unboostedp0; 
    abstol=1e-15, callback=cbs, maxiters=5e4, allbetas, saveat=1
)

simulationproportions!(unboostedsimulation)
insertcols!(
    unboostedsimulation, 
    :Rproportion => unboostedsimulation.R ./ unboostedsimulation.N
)
insertpreviousvalues!(
    unboostedsimulation, [ :PatientsProportion, :StaffProportion ], 43; 
    label=[ :Y, :I ]
)

function immunewaningvector(
    t::AbstractVector, locationcodes::AbstractVector, 
    beta2, beta1, betahalf, betazero, betaminus1
)
    singlesitet = t[findall(x -> x == locationcodes[1], locationcodes)]
    logchanges = [ 
        x == 0 ? 
            zero(beta2) : 
            beta2 * x^2 + beta1 * x + betahalf * sqrt(x) + betazero * log(x) + betaminus1 / x  
        for x ∈ singlesitet 
    ]
    f = ones(length(logchanges))
    for i ∈ eachindex(f)
        i == 1 && continue 
        if exp(logchanges[i]) == Inf
            f[i] = 0.0 
        else
            f[i] = f[(i - 1)] * (1 - exp(logchanges[i]) / (1 + exp(logchanges[i])))
        end
    end
    return f
end

function immunewaningvector(
    data::DataFrame, t, locationcodes, beta2, beta1, betahalf, betazero, betaminus1
)
    return immunewaningvector(
        getproperty(data, t), getproperty(data,locationcodes), 
        beta2, beta1, betahalf, betazero, betaminus1
    )
end

#
#
#
#
#
#
#
#beta2 = 0.02; beta1 = 0.3; betahalf =-0.2; betazero = 0.2; betaminus1 = -1.2

function immunewaning(t::AbstractVector, locationcodes::AbstractVector, infections::AbstractVector, beta2, beta1, betahalf, betazero, betaminus1)
    f = immunewaningvector(t, locationcodes, beta2, beta1, betahalf, betazero, betaminus1)
    immune = zeros(length(infections))
    for code ∈ unique(locationcodes)
        inds = findall(x -> x == code, locationcodes)
        inds0 = inds[1] - 1
        t_loc = t[inds]
        infect_loc = infections[inds]
        @assert length(t_loc) == length(infect_loc)
        for i ∈ eachindex(infect_loc)
            for j ∈ i+1:length(infect_loc)
                immune[(j + inds0)] += infect_loc[i] * f[(j - i)]
            end
        end
    end
    return immune
end

function immunewaning(data::DataFrame, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1)
    return immunewaning(getproperty(data, t), getproperty(data, locationcodes), getproperty(data, infections), beta2, beta1, betahalf, betazero, betaminus1)
end

function addimmunewaning!(data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1)
    insertcols!(data, :immune => immunewaning(data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1))
end

function replaceimmunewaning!(data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1)
    data.immune = immunewaning(data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1)
end

@model bcodevalue(prior) = b ~ prior

@model function q3fitmodel(data)
    # function for immune waning 
    ℓ_codes = length(unique(data.Code))
    beta2  ~ Normal(0, 1)
    beta1  ~ Normal(0, 1)
    betahalf ~ Normal(0, 1)
    betazero ~ Normal(0, 1) 
    betaminus1 ~ Normal(0, 1)

    b1 ~ Normal(0, 1)
    b2 ~ Normal(0, 1)
    b3 ~ Normal(0, 1)
    b4 ~ Normal(0, 1)
    b5 ~ Normal(0, 1)
    b6 ~ Normal(0, 1)
    b_code = Vector{typeof(b1)}(undef, ℓ_codes)
    for i ∈ eachindex(b_code)
        b_code[i] = @submodel prefix="b_code$i" b = bcodevalue(Normal(0, 1))
    end

    sigma2 ~ Exponential(1)

    replaceimmunewaning!(data, :t, :Code, :i4, beta2, beta1, betahalf, betazero, betaminus1)

    q3regr = (
        b1 .* data.dichy2int .+ 
        b2 .* data.y3 .+ 
        b3 .* data.i3 .+ 
        b4 .* data.i1 .+ 
        b5 .* data.y1 .+ 
        b6 .* data.immune .+ 
        sum([
            b_code[i] .* (data.Code .== code) for (i, code) ∈ enumerate(unique(data.Code))
        ])
    )

    data.i4 ~ MvNormal(q3regr, sigma2)
end

function pol_q3fitmodel(config)
    @unpack model, modelname, n_rounds, n_chains, seed = config
    roundconfig = @ntuple modelname model n_rounds=1 n_chains seed
    round1 = produce_or_load(pol_q3fitmodel1, roundconfig, datadir("sims"))
    for n ∈ 2:(n_rounds-1)
        roundconfig = @ntuple modelname model n_rounds=n n_chains seed
        roundn = produce_or_load(pol_q3fitmodelelement, roundconfig, datadir("sims"))
    end
    oldfilename = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds - 1)_seed=$(seed).jld2"
    pt = load(datadir("sims", oldfilename))["pt"]
    pt = increment_n_rounds!(pt, 1)
    new_pt = pigeons(pt)
    new_chains = Chains(new_pt)

    return Dict(
        "chain" => new_chains, 
        "pt" => new_pt, 
        "modelname" => modelname, 
        "n_rounds" => n_rounds, 
        "n_chains" => n_chains,
    )
end

function pol_q3fitmodel1(config)
    @unpack model, modelname, n_chains, seed = config
    fitted_pt = pigeons( ;
        target=TuringLogPotential(model),
        n_rounds=1,
        n_chains,
        multithreaded=true,
        record=[ traces; record_default() ],
        seed,
        variational=GaussianReference(),
    )
    fitted_chains = Chains(fitted_pt)
    return Dict(
        "chain" => fitted_chains, 
        "pt" => fitted_pt, 
        "modelname" => modelname, 
        "n_rounds" => n_rounds, 
        "n_chains" => n_chains,
    )
end

function pol_q3fitmodelelement()
    @unpack model, modelname, n_rounds, n_chains, seed = config
    oldfilename = "modelname=$(modelname)_n_chains=$(n_chains)_n_rounds=$(n_rounds - 1)_seed=$(seed).jld2"
    pt = load(datadir("sims", oldfilename))["pt"]
    pt = increment_n_rounds!(pt, 1)
    new_pt = pigeons(pt)
    new_chains = Chains(new_pt)

    return Dict(
        "chain" => new_chains, 
        "pt" => new_pt, 
        "modelname" => modelname, 
        "n_rounds" => n_rounds, 
        "n_chains" => n_chains,
    )

end


#=
fitted_pt = pigeons( ;
    target=TuringLogPotential(model),
    n_rounds=1,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=1,
    variational=GaussianReference(),
)

fitted_chains = Chains(fitted_pt)


insertcols!(unboostedsimulation, :immune => immune)
=#

dc = filter(:t => x -> x > 40, unboostedsimulation)
insertgroupedvalues!(dc, 1:7, 8:14, 15:21, 22)
select!(dc, :Code, :t, :i4, :i3, :i2, :i1, :y4, :y3, :y2, :y1, :S, :E, :I, :R)
#filter!(:t => x -> x > 340, dc)

# dichotomize exposure `y2` 
insertcols!(
    dc, 
    :dichy2 => cut(dc.y2, quantile(dc.y2, [ 0, 0.5, 1 ]); extend=true, labels=[ 0, 1 ])
)
# make an integer version that will be accepted by `glm`
insertcols!(dc, :dichy2int => levelcode.(dc.dichy2) .- 1)

# create column `immune` which is used in fitting values 
addimmunewaning!(dc, :t, :Code, :i4, 0, 0, 0, 0, 0)

unboostedq3model = q3fitmodel(dc)

# temp 
n_rounds = 2; id = 1

unboostedq3config = @ntuple modelname="unboostedq3model" model=unboostedq3model n_rounds n_chains=4 seed=100+id
unboostedq3dict = produce_or_load(pol_q3fitmodel, unboostedq3config, datadir("sims"))


# Q1 
# P(dichy2 == 1)
w1n = sum(dc.dichy2int .== 1) / (sum(dc.dichy2int .== 0) + sum(dc.dichy2int .== 1))
# P(dichy2 == 1 | I1 == i1, Y1 == y1, Code == code)
w1dformula = @formula(dichy2int ~ i1 + y1 + Code) 
w1dregr = glm(w1dformula, dc, Binomial(), LogitLink())
w1d = predict(w1dregr)
w1 = w1n ./ w1d
w1_98 = quantile(w1, 0.98)
w1 = [ x > w1_98 ? w1_98 : x for x ∈ w1 ]
q1 = sum(dc.i4 .* (dc.dichy2int .== 1) .* w1) / sum((dc.dichy2int .== 1) .* w1)

# Q2
insertcols!(dc, :notdichy2int => 1 .- dc.dichy2int)
# P(dichy2 == 0)
w2n = sum(dc.dichy2int .== 0) / (sum(dc.dichy2int .== 0) + sum( dc.dichy2int .== 1))
# P(dichy2 == 0 | I1 == i1, Y1 == y1, Code == code)
w2dformula = @formula(notdichy2int ~ i1 + y1 + Code) 
w2dregr = glm(w2dformula, dc, Binomial(), LogitLink())
w2d = predict(w2dregr)
w2 = w2n ./ w2d
w2_98 = quantile(w2, 0.98)
w2 = [ x > w2_98 ? w2_98 : x for x ∈ w2 ]
q2 = sum(dc.i4 .* (dc.notdichy2int .== 1) .* w2) / sum((dc.notdichy2int .== 1) .* w2)

# Q3
q3formula = @formula(i4 ~ dichy2int + y3 + i3 + i1 + y1 + immune + Code)
q3dregr = lm(q3formula, dc)
q3pred_a = predict(q3dregr)
# reverse effect of dichy2int 
_effectdichy2int = coef(q3dregr)[2]
q3pred = [ x == 1 ? -_effectdichy2int : _effectdichy2int for x ∈ dc.dichy2int] .+ q3pred_a
q3 = sum(q3pred .* (dc.dichy2int .== 1) .* w1) / sum((dc.dichy2int .== 1) .* w1)

nde = q3 - q1
nie = q2 - q3
#=
# correlation between `i4` and `y2` adjusting only for site and time
correlationformula = @formula(i4 ~ y2 + t^0.5 + log(t) + Code)
unboostedcorrelationregr = fit(LinearModel, correlationformula, dc)

w_unboosted, directeffectunboosted = let 
    @unpack w, directeffectregr = calculatedirecteffect(dc)
    #@unpack w, Ei2dregr = calculatedirecteffect(unboostedsimulation)
    ( w, directeffectregr )
end

pointestimatetotaleffect(dc)
pointestimatedirecteffect(dc)

pointestimatetotaleffect(unboostedsimulation)
pointestimatedirecteffect(unboostedsimulation)

#bs1 = bootstrap(pointestimatedirecteffect, dc, BasicSampling(1000))
#bci1 = confint(bs1, BasicConfInt(0.95))
=#
# Simulation with immune boostings, ε = 1

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
    300, u0, 800, boostedp0; 
    abstol=1e-15, callback=cbs, maxiters=5e4, allbetas, saveat=1
)

simulationproportions!(boostedsimulation)
insertcols!(boostedsimulation, :Rproportion => boostedsimulation.R ./ boostedsimulation.N)
insertpreviousvalues!(
    boostedsimulation, [ :λh, :PatientsProportion, :StaffProportion, :Rproportion ], 43; 
    label=[ :λh, :Y, :I, :R ]
)

dc2 = filter(:t => x -> x > 40, boostedsimulation)
insertgroupedvalues!(dc2, 1:7, 8:14, 15:21, 22)
#insertgroupedvalues!(dc2, 1, 8, 15, 22)
select!(dc2, :Code, :t, :i4, :i3, :i2, :i1, :y4, :y3, :y2, :y1, :S, :E, :I, :R)
#filter!(:t => x -> x > 340, dc2)

# dichotomize exposure `y2` 
insertcols!(
    dc2, 
    :dichy2 => cut(dc2.y2, quantile(dc2.y2, [ 0, 0.5, 1 ]); extend=true, labels=[ 0, 1 ])
)
# make an integer version that will be accepted by `glm`
insertcols!(dc2, :dichy2int => levelcode.(dc2.dichy2) .- 1)

boostedq3model = q3fitmodel(dc2)

boostedq3config = @ntuple modelname="boostedq3model" model=boostedq3model n_rounds n_chains=8 seed=110+id
boostedq3dict = produce_or_load(pol_q3fitmodel, boostedq3config, datadir("sims"))

# Q1 
# P(dichy2 == 1)
w1n = sum(dc2.dichy2int .== 1) / (sum(dc2.dichy2int .== 0) + sum( dc2.dichy2int .== 1))
# P(dichy2 == 1 | I1 = i1, Y1 = y1, Code == code)
w1dformula = @formula(dichy2int ~ i1 + y1 + Code) 
w1dregr = glm(w1dformula, dc2, Binomial(), LogitLink())
w1d = predict(w1dregr)
w1 = w1n ./ w1d
w1_98 = quantile(w1, 0.98)
w1 = [ x > w1_98 ? w1_98 : x for x ∈ w1 ]
q1 = sum(dc2.i4 .* (dc2.dichy2int .== 1) .* w1) / sum((dc2.dichy2int .== 1) .* w1)

# Q2
insertcols!(dc2, :notdichy2int => 1 .- dc2.dichy2int)
# P(dichy2 == 0)
w2n = sum(dc2.dichy2int .== 0) / (sum(dc2.dichy2int .== 0) + sum( dc2.dichy2int .== 1))
# P(dichy2 == 0 | I1 = i1, Y1 = y1, Code == code)
w2dformula = @formula(notdichy2int ~ i1 + y1 + Code) 
w2dregr = glm(w2dformula, dc2, Binomial(), LogitLink())
w2d = predict(w2dregr)
w2 = w2n ./ w2d
w2_98 = quantile(w2, 0.98)
w2 = [ x > w2_98 ? w2_98 : x for x ∈ w2 ]
q2 = sum(dc2.i4 .* (dc2.notdichy2int .== 1) .* w2) / sum((dc2.notdichy2int .== 1) .* w2)

# Q3
q3formula = @formula(i4 ~ dichy2 + y3 + i3 + i1 + i2 + y1 + Code)
q3dregr = lm(q3formula, dc2)
q3pred_a = predict(q3dregr)
# reverse effect of dichy2int 
_effectdichy2int = coef(q3dregr)[2]
q3pred = [ x == 1 ? -_effectdichy2int : _effectdichy2int for x ∈ dc2.dichy2int] .+ q3pred_a
q3 = sum(q3pred .* (dc2.dichy2int .== 1) .* w1) / sum((dc2.dichy2int .== 1) .* w1)

nde = q3 - q1
nie = q2 - q3



boostedcorrelationregr = fit(LinearModel, correlationformula, dc2)

w_boosted, directeffectboosted = let 
    @unpack w, directeffectregr = calculatedirecteffect(dc2)
    ( w, directeffectregr )
end

pointestimatetotaleffect(dc2)
pointestimatedirecteffect(dc2)

pointestimatetotaleffect(boostedsimulation)
pointestimatedirecteffect(boostedsimulation)

#bs2 = bootstrap(pointestimatedirecteffect, dc2, BasicSampling(1000))
#bci2 = confint(bs2, BasicConfInt(0.95))


###########################################################################################
# UK data
###########################################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load UK data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coviddata = CSV.read(datadir("exp_raw", "dataset.csv"), DataFrame)

# make hospital codes categorical variables 
rename!(coviddata, :Codes => "StringCodes")
insertcols!(coviddata, 1, :Code => CategoricalArray(coviddata.StringCodes))

# need to find the first date for each hospital
maxstartdate = let 
    hospcodes = unique(coviddata.Code)
    startdate = zeros(Date, length(hospcodes))
    for (i, c) ∈ enumerate(hospcodes)
        _tdf = filter(:Code => x -> x == c, coviddata)
        startdate[i] = minimum(_tdf.Date)
    end
    maximum(startdate)
end
# 2020-04-02

filter!(:Date => x -> x >= Date("2020-04-02"), coviddata)
insertcols!(coviddata, :t => Dates.value.(coviddata.Date - Date("2020-04-02")))

# calculate values for the analysis
simulationproportions!(
    coviddata; 
    I=:StaffAbsences, N=:StaffTotal, Y=:CovidBeds, M=:AllBeds
)

insertpreviousvalues!(
    coviddata, [ :PatientsProportion, :StaffProportion ], 43; 
    label=[ :Y, :I ]
)

# remove all data before 2021-01-09, when staff largely unvaccinated, and for the next 40 days
Date("2021-01-08") + Day(40)
#2021-02-17
filter!(:Date => x -> x >= Date("2021-02-17"), coviddata)
# and for the next 0

# remove data where values are missing or NaN
for name ∈ names(coviddata)
    name ∈ [ "Code", "StringCodes", "Date" ] && continue
    filter!(name => x -> !ismissing(x) && !isnan(x), coviddata)
end

#insertgroupedvalues!(coviddata, 1:14, 15:28, 29:42, 23)
insertgroupedvalues!(coviddata, 1, 7, 14, 21)

# remove outlying values where proportions == 1 
for name ∈ [ :i4, :i3, :i2, :i1, :y4, :y3, :y2, :y1 ]
    filter!(name => x -> 0 <= x < 1, coviddata)
end

# correlation between `i4` and `y2` adjusting only for site and time

covidcorrelationregr = fit(LinearModel, correlationformula, coviddata)

w_covid, directeffectcovid = let 
    @unpack w, directeffectregr = calculatedirecteffect(coviddata)
    ( w, directeffectregr )
end

pointestimatetotaleffect(coviddata)
pointestimatedirecteffect(coviddata)

#bscovid = bootstrap(pointestimatedirecteffect, coviddata, BasicSampling(1000))
#bcicovid = confint(bscovid, BasicConfInt(0.95))

fig, ax = scatter(coviddata.t, coviddata.i4)
xs = 321:1:797
ys = [ -0.00635744 * t^0.5 + 0.0804562 * log(t) - 0.343375 for t ∈ xs ]
lines!(ax, xs, ys; color=:red)

fig

scatter(coviddata.y2, coviddata.i4)
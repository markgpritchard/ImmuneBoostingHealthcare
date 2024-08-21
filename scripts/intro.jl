
using DrWatson
@quickactivate "ImmuneBoostingHealthcare"
#=
#using Bootstrap, CairoMakie, CSV, DataFrames, Random, Turing, Pigeons
using CSV, DataFrames, Random, Turing#, Pigeons
using SparseArrays 

using CairoMakie

###########################################################################################
# Functions 
###########################################################################################

#using CategoricalArrays, CSV, DataFrames, Dates, DifferentialEquations, Distributions, GLM, StaticArrays
#using CategoricalArrays, CSV, DataFrames, Dates, DifferentialEquations, Distributions, StaticArrays
using CategoricalArrays, CSV, DataFrames, Dates, Distributions, StaticArrays
#using LinearAlgebra: I
import Base: minimum
using CSV, DataFrames, Dates, Distributions, StaticArrays


###########################################################################################
# Structs 
###########################################################################################

makechangematrix(compartments, pairs) = _makechangematrix(zeros, compartments, pairs)
spmakechangematrix(compartments, pairs) = _makechangematrix(spzeros, compartments, pairs)

function _makechangematrix(f, compartments, pairs::Vector{Pair{Int64, Int64}})
    mat = f(Int, length(pairs), compartments)
    for (i, p) ∈ enumerate(pairs)
        mat[i, p.first] = -1
        mat[i, p.second] = 1
    end
    return mat
end

HOSPITALSISCHANGEMATRIX = makechangematrix(
    19,
    [
        # hospital patients: SIS with a "diagnosed" compartment and "diagnosed and recovered"
        1 => 2,  # Sp -> Ip, patient infections
        2 => 3,  # Ip -> Ip*, patient diagnosis
        2 => 1,  # Ip -> Sp, patient recovery 
        3 => 4,  # Ip* -> Sp*, patient diagnosed recovery (not susceptible to reinfection while still in hospital)
        # healthcare workers: SIS with ten "isolated" (non-infectious) compartments  
        5 => 6,  # Sh -> Ih, healthcare worker infections
        6 => 10,  # Ih -> X1h*, healthcare worker diagnosis 1 
        6 => 5,  # Ih -> Sh, healthcare worker recovery 
        # isolated compartment transitions come later
        # discharges
        1 => 7,  # Sp -> Sc, susceptible discharges
        2 => 8,  # Ip -> Ic, infectious discharges 
        3 => 9,  # Ip* -> Ic*, diagnosed infectious discharges
        4 => 7,  # Sp* -> Sc, recovered discharges
        # community (will be approximately deterministic): SEIIRRRS with two "diagnosed" compartments
        7 => 8,  # Sc -> Ic, community infections
        8 => 9,  # Ic -> Ic*, community diagnosis
        8 => 7,  # Ic -> Sc, community recovery 
        9 => 7,  # Ic* -> Sc, community diagnosed recovery 
        # admissions (will be calculated to replace discharges)
        7 => 1,  # Sp -> Sc, susceptible admissions
        8 => 2,  # Ic -> Ip, infectious admissions
        9 => 5,  # Ic* -> Ip*, diagnosed infectious admissions
        # progression of isolated healthcare workers (all occur at rate = 1)
        10 => 11,  # X1h* -> X2h*
        11 => 12,  # X2h* -> X3h*
        12 => 13,  # X3h* -> X4h*
        13 => 14,  # X4h* -> X5h*
        14 => 15,  # X5h* -> X6h*
        15 => 16,  # X6h* -> X7h*
        16 => 17,  # X7h* -> X8h*
        17 => 18,  # X8h* -> X9h*
        18 => 19,  # X9h* -> X10h*
        19 => 5,  # X10h* -> Sh
    ]
)

HOSPITALSEIIRRRSCHANGEMATRIX = spmakechangematrix(
    34,
    [
        # hospital patients: SEIIR with two "diagnosed" compartments plus "diagnosed and recovered"
        1 => 2,  # Sp -> Ep, patient infections
        2 => 3,  # Ep -> I1p, patient progression 1
        3 => 4,  # I1p -> I2p, patient progression 2
        3 => 5,  # I1p -> I1p*, patient diagnosis 1 
        4 => 6,  # I2p -> I2p*, patient diagnosis 2
        5 => 6,  # I1p* -> I2p*, patient diagnosed progression
        4 => 7,  # I2p -> Rp, patient recovery 
        6 => 8,  # I2p* -> Rp*, patient diagnosed recovery 
        # healthcare workers: SEIIRRRS with ten "isolated" (non-infectious) compartments  
        9 => 10,  # Sh -> Eh, healthcare worker infections
        10 => 11,  # Eh -> I1h, healthcare worker progression 1
        11 => 12,  # I1h -> I2h, healthcare worker progression 2
        11 => 25,  # I1h -> X1h*, healthcare worker diagnosis 1 
        12 => 25,  # I2h -> X1h*, healthcare worker diagnosis 2
        12 => 13,  # I2h -> R1h, healthcare worker recovery 
        13 => 14,  # R1h -> R2h, healthcare worker waning 1 
        14 => 15,  # R2h -> R3h, healthcare worker waning 2 
        15 => 9,  # R3h -> Sh, healthcare worker waning 3
        # isolated compartment transitions come later
        # discharges
        1 => 16,  # Sp -> Sc, susceptible discharges
        2 => 17,  # Ep -> Ec, exposed discharges
        3 => 18,  # I1p -> I1c, infectious discharges 1
        4 => 19,  # I2p -> I2c, infectious discharges 2
        5 => 20,  # I1p* -> I1c*, diagnosed infectious discharges 1
        6 => 21,  # I2p* -> I2c*, diagnosed infectious discharges 2
        7 => 23,  # Rp -> R2c, immune discharges
        8 => 23,  # Rp* -> R2c, recovered discharges
        # community (will be approximately deterministic): SEIIRRRS with two "diagnosed" compartments
        16 => 17,  # Sc -> Ec, community infections
        17 => 18,  # Ec -> I1c, community progression 1
        18 => 19,  # I1c -> I2c, community progression 2
        18 => 20,  # I1c -> I1c*, community diagnosis 1 
        19 => 21,  # I2c -> I2c*, community diagnosis 2
        20 => 21,  # I1c* -> I2c*, community diagnosed progression
        19 => 22,  # I2c -> R1c, community recovery 
        21 => 22,  # I2c* -> R1c, community diagnosed recovery 
        22 => 23,  # R1c -> R2c, community waning 1 
        23 => 24,  # R2c -> R3c, community waning 2 
        24 => 16,  # R3c -> Sc, community waning 3
        # admissions (will be calculated to replace discharges)
        16 => 1,  # Sp -> Sc, susceptible admissions
        17 => 2,  # Ec -> Ep, exposed admissions
        18 => 3,  # I1c -> I1p, infectious admissions 1
        19 => 4,  # I2c -> I2p, infectious admissions 2
        20 => 5,  # I1c* -> I1p*, diagnosed infectious admissions 1
        21 => 6,  # I2c* -> I2p*, diagnosed infectious admissions 2
        22 => 7,  # R1c -> Rp, immune admissions 1
        23 => 7,  # R2c -> Rp, immune admissions 2
        24 => 7,  # R3c -> Rp, immune admissions 3
        # progression of isolated healthcare workers (all occur at rate = 1)
        25 => 26,  # X1h* -> X2h*
        26 => 27,  # X2h* -> X3h*
        27 => 28,  # X3h* -> X4h*
        28 => 29,  # X4h* -> X5h*
        29 => 30,  # X5h* -> X6h*
        30 => 31,  # X6h* -> X7h*
        31 => 32,  # X7h* -> X8h*
        32 => 33,  # X8h* -> X9h*
        33 => 34,  # X9h* -> X10h*
        34 => 22,  # X10h* -> R1h
    ]
)

#=
function _changematrixsources(mat)
    ℓ = size(mat, 1)
    sources = Vector{Int}(undef, ℓ)
    for i ∈ axes(mat, 1)
        @assert sum(mat[i, :] .== -1) == 1  # exactly one source per row 
        sources[i] = findfirst(x -> x == -1, mat[i, :])
    end
    return sources
end

WXYYZSEIIRRRSCHANGEMATRIX = sparse


XYXSISCHANGEMATRIX = [
    -1   1   0   0  # patient infections 
     1  -1   0   0  # patient recoveries
     0   0  -1   1  # healthcare worker infections 
     0   0   1  -1  # healthcare worker recoveries
    -1   0   0   0  # susceptible discharges
     0  -1   0   0  # infected discharges
]

WXYYZSEIIRRRSCHANGEMATRIX = [
    #Sp   Ep   I1p  I2p  Rp   I1p* I2p* Rp*  Sh   Eh   I1h  I2h  R1h  R2h  R3h  I1h* 
    # community infections 
    -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # patient infections 
     0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0  # patient infection progressions 1
     0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0  # patient infection progressions 2
     0   0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0  # patient recoveries
     0   0  -1   0   0   1   0   0   0   0   0   0   0   0   0   0  # patient diagnosis 1 
     0   0   0  -1   0   0   1   0   0   0   0   0   0   0   0   0  # patient diagnosis 2 
     0   0   0   0   0  -1   1   0   0   0   0   0   0   0   0   0  # diagnosed patient infection progressions
     0   0   0   0   0   0  -1   1   0   0   0   0   0   0   0   0  # diagnosed patient infection recoveries
     0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0   0  # healthcare worker infections 
     0   0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0  # healthcare infection progressions 1
     0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0   0  # healthcare infection progressions 2
     0   0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0  # healthcare recoveries
     0   0   0   0   0   0   0   0   0   0   0   0  -1   1   0   0  # healthcare waning 1
     0   0   0   0   0   0   0   0   0   0   0   0   0  -1   1   0  # healthcare waning 2
     0   0   0   0   0   0   0   0   1   0   0   0   0   0  -1   0  # healthcare waning 3
     0   0   0   0   0   0   0   0   0   0   0   0   1  -1   0   0  # healthcare boosting 2
     0   0   0   0   0   0   0   0   0   0   0   0   1   0  -1   0  # healthcare boosting 3
     0   0   0   0   0   0   0   0   0   0  -1   0   0   0   0   1  # healthcare diagnosis 1 
     0   0   0   0   0   0   0   0   0   0   0  -1   0   0   0   1  # healthcare diagnosis 2      
    -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # susceptible discharges
     0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # exposed discharges
     0   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0  # infectious 1 discharges
     0   0   0  -1   0   0   0   0   0   0   0   0   0   0   0   0  # infectious 2 discharges
     0   0   0   0  -1   0   0   0   0   0   0   0   0   0   0   0  # immune discharges
     0   0   0   0   0  -1   0   0   0   0   0   0   0   0   0   0  # diagnosed infectious 1 discharges
     0   0   0   0   0   0  -1   0   0   0   0   0   0   0   0   0  # diagnosed infectious 2 discharges
     0   0   0   0   0   0   0  -1   0   0   0   0   0   0   0   0  # diagnosed immune discharges
]

WXYYZSEIIRRRSSOURCES = _changematrixsources(WXYYZSEIIRRRSCHANGEMATRIX)
=#


function wxyyzseiirrrs!(u, p)
    # forces of infection 
    λh = p.βhh * sum(@view u[8:9]) / sum(@view u[6:12]) + 
        p.βhp * sum(@view u[3:4]) / sum(@view u[1:5]) + 
        p.λc
    λp = p.βph * sum(@view u[8:9])  / sum(@view u[6:12]) + 
        p.βpp * sum(@view u[3:4]) / sum(@view u[1:5])

    rates = [
        [
            λp * u[1],  # patient infections 
            p.ρ * u[2]  # patient infection progressions 1
        ];
        2 .* p.γ .* u[3:4];  # patient infection progressions / recoveries
        p.ηp .* u[3:4];  # patient diagnosis 
        2 .* p.γ .* u[6:7];  # diagnosed patient infection progressions / recoveries
        [
            λh * u[9],  # healthcare worker infections 
            p.ρ * u[10],  # healthcare infection progressions 1
        ];
        2 .* p.γ .* u[11:12];  # healthcare infection progressions  / recoveries
        3 .* p.ω .* u[13:15];  # healthcare waning 
        λh * p.ψ .* u[14:15];  # healthcare boosting 
        p.ηh .* u[11:12];  # healthcare diagnosis
        p.δ1 .* u[1:2];  # susceptible / exposed discharges
        p.δ2 .* u[3:4];  # infectious discharges
        [ p.δ1 .* u[5] ]; # immune discharges
        p.δ2 .* u[6:7];  # diagnosed infectious discharges
        [ p.δ1 .* u[8] ]  # diagnosed immune discharges
    ]

    events = [ taustepvalue(u, rate, i) for (rate, i) ∈ zip(rates, WXYYZSEIIRRRSSOURCES) ]
        
    return events
end

function taustepvalue(source::T, rate) where T
    estimate::T = rand(Poisson(rate))
    if estimate > source
        return source 
    else 
        return estimate 
    end
end

taustepvalue(sources::AbstractVector, rate, ind) = taustepvalue(sources[ind], rate)

function taustep!(ratesfunction, u, p, changematrix)
    rates = ratesfunction(u, p)
    events = [ rand(Poisson(rate)) for rate ∈ rates ]

    @assert size(changematrix, 1) == length(events)
    @assert size(changematrix, 2) == length(u)

    for (i, event) ∈ enumerate(events) 
        u += event * changematrix[i, :] 
        if minimum(u) < 0  # then this cycle has made more events than it could have 
            reduction = -minimum(u)  # how many events were too many 
            u += -event * changematrix[i, :]  # undo what you just did 
            u += (event - reduction) * changematrix[i, :]  # re-do it with reduced number of events
        end 
    end 
    return u
end

function xyxsisrates(u, p)
    # forces of infection 
    λh = p.βhh * u[4] / sum(@view u[3:4]) + p.βhp * u[2] / sum(@view u[1:2])
    λp = p.βph * u[4]  / sum(@view u[3:4]) + p.βpp * u[2] / sum(@view u[1:2])

    rates = [
        λp * u[1],  # patient infections 
        p.γ * u[2],  # patient recoveries
        λh * u[3],  # healthcare worker infections 
        p.γ * u[4],  # healthcare worker recoveries
        p.δ1 * u[1],  # susceptible discharges
        p.δ2 * u[2]  # infected discharges
    ]

    return rates
end

function wxyyzseiirrrsrates(u, p)
    # forces of infection 
    λh = p.βhh * sum(@view u[8:9]) / sum(@view u[6:12]) + 
        p.βhp * sum(@view u[3:4]) / sum(@view u[1:5]) + 
        p.λc
    λp = p.βph * sum(@view u[8:9])  / sum(@view u[6:12]) + 
        p.βpp * sum(@view u[3:4]) / sum(@view u[1:5])

    rates = [
        [
            λp * u[1],  # patient infections 
            p.ρ * u[2],  # patient infection progressions 1
        ];
        2 .* p.γ .* u[3:4];  # patient infection progressions / recoveries
        [
            λh * u[6],  # healthcare worker infections 
            p.ρ * u[7],  # healthcare infection progressions 1
        ];
        2 .* p.γ .* u[8:9];  # healthcare infection progressions  / recoveries
        3 .* p.ω .* u[10:12];  # healthcare waning 
        p.δ1 .* u[1:2];  # susceptible discharges
        p.δ2 .* u[3:4];  # infectious discharges
        [ p.δ1 * u[5] ]  # immune discharges
    ]

    return rates
end

wxyyzseiirrrs!(u, p) = taustep!(wxyyzseiirrrsrates, u, p, WXYYZSEIIRRRSCHANGEMATRIX)

function wxyyzseiirrrs(u0, p, t, αvec, rvec)
    umat = zeros(Int, t, 12)
    u = u0
    poph = sum(@view u0[1:5])
    for i ∈ 1:t
        umat[i, :] .= u
        u = wxyyzseiirrrs!(u, p)
        upopi = sum(@view u[1:5])
        u1 = round(Int, αvec[i] * (poph - upopi))
        u2 = poph - upopi - u1 
        u[1] += u1 
        u[2] += u2
    end
    return umat
end

function xyxsisdf(u0, p, t, αvec, id::Int=1)
    umat = xyxsis(u0, p, t, αvec)
    df = DataFrame(
        :t => 1:t,
        :Hospital => id,
        :X => umat[:, 1],
        :Y => umat[:, 2],
        :M => [ sum(@view umat[i, 1:2]) for i ∈ 1:t ],
        :S => umat[:, 1],
        :I => umat[:, 1],
        :N => [ sum(@view umat[i, 3:4]) for i ∈ 1:t ],
    )
    return df 
end

function xyxsisdf(u0, p, t, αvec, ids::AbstractVector{T}) where T
    df = DataFrame(
        :t => Int[ ],
        :Hospital => T[ ],
        :X => Int[ ],
        :Y => Int[ ],
        :M => Int[ ],
        :S => Int[ ],
        :I => Int[ ],
        :N => Int[ ],
    )

    for id ∈ ids
        append!(df, xyxsisdf(u0, p, t, αvec, id))
    end 

    return df
end




ic = let 
    i = zeros(800)
    i[1] = 0.005
    for j ∈ 2:800
        i[j] = 0.8 * i[j-1] + 0.4 * (1 + 0.75 * cos(2π * j / 365)) * i[j-1] * (1 - i[j-1])
    end
    i
end

asdf = xyxsis([ 450, 0, 900, 0 ], ( βhh=0.2, βhp=0.1, βph=0.4, βpp=0.05, γ=0.2, δ1=0.2, δ2=0.16 ), 800, ic)

u0 = [ 450, 0, 900, 0 ]

uvector = let 
    uv = zeros(Int, 800, 4)
    u = u0 
    for i ∈ 1:800
        println("i=$i, u=$u")
        uv[i, :] .= u
        p = ( α=(ic[i]), βhh=0.2, βhp=0.1, βph=0.4, βpp=0.05, γ=0.2, δ1=0.2, δ2=0.16 )
        u = xyxsis!(u, p)
        un = u[1] + u[2]
        u1 = round(Int, ic[i] * (450 - un))
        u2 = 450 - un - u1 
        u[1] += u1 
        u[2] += u2
    end
    uv
end



#=
function model(u::AbstractVector{T}, p, t) where T 
    # compartments 
    Sp, I1p, I2p, Rp, XI1p, XI2p, XRp, Sh, I1h, I2h, R1h, R2h, R3h, Vh, X1h, X2h, X3h, X4h, X5h, X6h, X7h, X8h, X9h, X10h, Sc, Ic, Rc = u 

    # parameters (some parameters are hard-coded)
    α1, α2, α3, βhc, βhh, βhp, βph, βpp, δp, δh, η1, η2, ν, ψ, ω

    u2 = Vector{T}(undef, 27)
    
    u2[1] = Sp * 5//6 
=#


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
# Susceptible--infectious--susceptible model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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
=#
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

## Download hospital data 

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
#=
insertcols!(
    hospitaldata,
    :VolumePerBed => hospitaldata.HeatedVolume ./ hospitaldata.TotalBeds,
    :ProportionSingleBeds => min.(hospitaldata.SingleBedsEnsuite ./ hospitaldata.TotalBeds, 1.0)
)
=#
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

communitycases = [ sum(@view COMMUNITYSOL[i][3:4]) for i ∈ eachindex(COMMUNITYSOL) ]

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
        :CommunityCases => [ communitycases[sim.t[i] >= 800 ? 80 : round(Int, sim.t[i] / 10, RoundDown)+1] for i ∈ axes(sim, 1) ],
        :HeatedVolumePerBed => [ vpds[levelcode(sim.Code[i])] for i ∈ axes(sim, 1) ],
        :ProportionSingleBeds => [ psbs[levelcode(sim.Code[i])] for i ∈ axes(sim, 1) ],
        :DiagnosedY => [ rand(Binomial(round(Int, y), 0.9)) for y ∈ sim.Y ],
        :DiagnosedI => [ rand(Binomial(round(Int, y), 0.9)) for y ∈ sim.I ],
    )
    insertcols!(
        sim,
        :PatientsProportion => sim.DiagnosedY ./ sim.M,
        :StaffProportion => sim.DiagnosedI ./ sim.N,
    )
    insertpreviousvalues!(
        sim, [ :PatientsProportion, :StaffProportion ], 31; 
        label=[ :Y, :I ]
    )
    insertgroupedvalues!(sim, 1:10, 11:20, 21:30, 31)
end

function f(beta0, beta1, beta2, beta3, beta4, beta5)
    g = [ beta0 + beta1 * x + beta2 * log(x) + beta3 / x + beta4 * sqrt(x) + beta5 * x^2 for x ∈ 1:800 ]
    gproportion = exp.(g) ./ (exp.(g) .+ 1)
    h = [ prod(@view gproportion[1:i]) for i ∈ eachindex(gproportion) ]
    return h
end

#=

julia> @benchmark f(0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  142.500 μs … 135.339 ms  ┊ GC (min … max): 0.00% … 99.83%
 Time  (median):     166.800 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   201.093 μs ±   1.353 ms  ┊ GC (mean ± σ):  6.72% ±  1.00%

     █       
  ▃▃▅█▆▄▅▃▂▂▂▁▂▂▃▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  142 μs           Histogram: frequency by time          497 μs <

 Memory estimate: 19.12 KiB, allocs estimate: 3.



=#

#=
w_vec = f(0.2, 0.1, 0.3, -0.2, 0.2, 0.1)
r_vec = zeros(size(unboostedsimulation, 1))
for i ∈ axes(unboostedsimulation, 1)
    unboostedsimulation.t[i] == 0 && continue 
    r_vec[i] = sum([ unboostedsimulation.StaffProportion[i-x] * w_vec[x] for x ∈ 1:(Int(unboostedsimulation.t[i])) ])
end
=#

function makervec(data, beta0, beta1, beta2, beta3, beta4, beta5)
    w_vec = f(beta0, beta1, beta2, beta3, beta4, beta5)
    #r_vec = zeros(size(data, 1))
    #r_vec = [ data.t[i] == 0 ? 0.0 : sum(reverse(@view data.StaffProportion[i-(Int(data.t[i])):i-1]) .* @view w_vec[1:(Int(data.t[i]))]) for i ∈ axes(data, 1) ]
    tvec = Int.(data.t)
    #r_vec = [ t == 0 ? 0.0 : sum(view(data.StaffProportion, i-1:-1:i-t) .* view(w_vec, 1:t)) for (i, t) ∈ enumerate(tvec) ]
    sp = data.StaffProportion
    @fastmath r_vec = [ t == 0 ? 0.0 : sum(view(sp, i-1:-1:i-t) .* view(w_vec, 1:t)) for (i, t) ∈ enumerate(tvec) ]
    #@fastmath r_vec = [ data.t[i] == 0 ? 0.0 : sum(view(data.StaffProportion, i-1:-1:i-Int(data.t[i])) .* view(w_vec, 1:Int(data.t[i]))) for i ∈ axes(data, 1) ]
    #for i ∈ axes(data, 1)
    #    data.t[i] == 0 && continue
    #    @fastmath r_vec[i] = sum(reverse(@view data.StaffProportion[i-(Int(data.t[i])):i-1]) .* @view w_vec[1:(Int(data.t[i]))])
        #data.t[i] == 0 && continue
        #r_vec[i] = sum(reverse(@view data.StaffProportion[i-(Int(data.t[i])):i-1]) .* @view w_vec[1:(Int(data.t[i]))])
    #end
    return r_vec
end

#=
using BenchmarkTools


#=
julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 14.864 s (9.11% GC) to evaluate,
 with a memory estimate of 7.30 GiB, over 428917482 allocations.

julia> @benchmark makervec!(g, gproportion, w_vec, unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 14.894 s (9.02% GC) to evaluate,
 with a memory estimate of 7.30 GiB, over 428917479 allocations.

julia> @benchmark makervec!(r_vec, g, gproportion, w_vec, unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)     
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 15.105 s (8.97% GC) to evaluate,
 with a memory estimate of 7.29 GiB, over 428917476 allocations.

julia> @benchmark makervec!(r_vec, g, gproportion, w_vec, unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)     
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 41.677 s (8.10% GC) to evaluate,
 with a memory estimate of 16.23 GiB, over 982812460 allocations.

julia> @benchmark makervec!(r_vec, g, gproportion, w_vec, unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.239 s …    1.448 s  ┊ GC (min … max): 12.06% … 12.83%
 Time  (median):     1.348 s               ┊ GC (median):    12.32%
 Time  (mean ± σ):   1.346 s ± 111.521 ms  ┊ GC (mean ± σ):  12.40% ±  0.40%

  █    █                                                █  █  
  █▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁█ ▁
  1.24 s         Histogram: frequency by time         1.45 s <

 Memory estimate: 1.77 GiB, allocs estimate: 5214857.
 
julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.060 s …   1.116 s  ┊ GC (min … max): 11.99% … 11.37%
 Time  (median):     1.074 s              ┊ GC (median):    11.84%
 Time  (mean ± σ):   1.082 s ± 22.150 ms  ┊ GC (mean ± σ):  12.22% ±  1.00%

  █        █    █                 █                       █
  █▁▁▁▁▁▁▁▁█▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.06 s         Histogram: frequency by time        1.12 s <

 Memory estimate: 1.77 GiB, allocs estimate: 5214863.
 
 #####


julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  1.684 s …    2.588 s  ┊ GC (min … max): 42.11% … 43.71%
 Time  (median):     2.566 s               ┊ GC (median):    44.09%
 Time  (mean ± σ):   2.279 s ± 515.842 ms  ┊ GC (mean ± σ):  45.89% ±  4.49%

  █                                                       ██
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁██ ▁
  1.68 s         Histogram: frequency by time         2.59 s <

 Memory estimate: 2.37 GiB, allocs estimate: 14824760.
 
julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  1.645 s …    2.293 s  ┊ GC (min … max): 37.77% … 41.43%
 Time  (median):     1.821 s               ┊ GC (median):    40.01%
 Time  (mean ± σ):   1.919 s ± 335.051 ms  ┊ GC (mean ± σ):  39.93% ±  1.85%

  █              █                                         █
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.64 s         Histogram: frequency by time         2.29 s <

 Memory estimate: 1.81 GiB, allocs estimate: 9079994. 

julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  773.415 ms …    1.362 s  ┊ GC (min … max): 32.90% … 47.17%
 Time  (median):        1.086 s               ┊ GC (median):    37.53%
 Time  (mean ± σ):      1.045 s ± 239.990 ms  ┊ GC (mean ± σ):  39.97% ±  8.91%

  █      █                        █       █                   █
  █▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  773 ms           Histogram: frequency by time          1.36 s <

 Memory estimate: 1.18 GiB, allocs estimate: 7327420.


julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  761.423 ms …    1.393 s  ┊ GC (min … max): 31.45% … 45.20%
 Time  (median):     848.242 ms               ┊ GC (median):    33.21%
 Time  (mean ± σ):      1.028 s ± 297.130 ms  ┊ GC (mean ± σ):  39.32% ±  7.12%

  █     █ █                                           █       █
  █▁▁▁▁▁█▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█ ▁
  761 ms           Histogram: frequency by time          1.39 s <

 Memory estimate: 1.13 GiB, allocs estimate: 6117525.

julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 6 samples with 1 evaluation.
 Range (min … max):  734.703 ms …    1.374 s  ┊ GC (min … max): 31.67% … 47.45%
 Time  (median):     899.307 ms               ┊ GC (median):    31.33%
 Time  (mean ± σ):   990.718 ms ± 240.868 ms  ┊ GC (mean ± σ):  38.27% ±  8.94%

  ▁          █       ▁                       ▁                ▁
  █▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  735 ms           Histogram: frequency by time          1.37 s <

 Memory estimate: 1.10 GiB, allocs estimate: 6273229.

julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 6 samples with 1 evaluation.
 Range (min … max):  667.358 ms …    1.070 s  ┊ GC (min … max): 31.65% … 33.36%
 Time  (median):     897.182 ms               ┊ GC (median):    38.41%
 Time  (mean ± σ):   882.302 ms ± 191.369 ms  ┊ GC (mean ± σ):  37.13% ±  6.20%

  ▁  ▁            ▁                                    ▁      █
  █▁▁█▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█ ▁
  667 ms           Histogram: frequency by time          1.07 s <

 Memory estimate: 882.63 MiB, allocs estimate: 3087999.

julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 12 samples with 1 evaluation.
 Range (min … max):  335.745 ms … 658.931 ms  ┊ GC (min … max): 42.81% … 39.47%
 Time  (median):     372.866 ms               ┊ GC (median):    44.18%
 Time  (mean ± σ):   418.665 ms ± 105.730 ms  ┊ GC (mean ± σ):  43.49% ±  2.27%

  ▁█ ▁▁ ▁▁▁    ▁                       ▁    ▁                 ▁
  ██▁██▁███▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  336 ms           Histogram: frequency by time          659 ms <

 Memory estimate: 603.12 MiB, allocs estimate: 1219732.

julia> @benchmark makervec(unboostedsimulation, 0.2, -0.3, 0.1, -0.1, 0.4, -0.03)
BenchmarkTools.Trial: 14 samples with 1 evaluation.
 Range (min … max):  288.043 ms … 441.850 ms  ┊ GC (min … max): 40.58% … 42.34%
 Time  (median):     362.122 ms               ┊ GC (median):    45.36%
 Time  (mean ± σ):   371.542 ms ±  41.678 ms  ┊ GC (mean ± σ):  42.65% ±  2.69%

  █                 █   ███ ██   █   ██           ██       █  █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁███▁██▁▁▁█▁▁▁██▁▁▁▁▁▁▁▁▁▁▁██▁▁▁▁▁▁▁█▁▁█ ▁
  288 ms           Histogram: frequency by time          442 ms <

 Memory estimate: 605.29 MiB, allocs estimate: 1336027.
=#
=#

@model function fitmodel(data, recordedy2)
    #beta0 ~ Normal(0, 10)
    #beta1 ~ Normal(0, 10)
    #beta2 ~ Normal(0, 10)
    #beta3 ~ Normal(0, 10)

    gamma0 ~ Normal(0, 1)
    gamma1 ~ Normal(0, 1)
    gamma2 ~ Normal(0, 1)
    gamma3 ~ Normal(0, 1)
    gamma4 ~ Normal(0, 1)
    gamma5 ~ Normal(0, 1)
    gamma6 ~ Normal(0, 1)
    gamma7 ~ Normal(0, 1)
    gamma8 ~ Normal(0, 1)
    gamma9 ~ Normal(0, 1)
    gamma10 ~ Normal(0, 1)

    sigma2 ~ Exponential(1)
    
    #r_vec = makervec(data, beta0, beta1, beta2, beta3, 0, 0)

    tvec = data.t

    #y2predod = [ gamma0 + gamma1 * data.y1[i] + gamma2 * data.i1[i] + gamma3 * data.i2[i] + gamma4 * data.y3[i] + gamma5 * data.i3[i] + gamma6 * r_vec[i] + gamma7 * data.ProportionSingleBeds[i] + gamma8 * data.HeatedVolumePerBed[i] + gamma9 * data.CommunityCases[i] for i ∈ axes(data, 1) ]
    y2predod = [ gamma0 + gamma1 * data.y1[i] + gamma2 * data.i1[i] + gamma3 * data.i2[i] + gamma4 * data.y3[i] + gamma5 * data.i3[i] + gamma6 * data.ProportionSingleBeds[i] + gamma7 * data.HeatedVolumePerBed[i] + gamma8 * data.CommunityCases[i] + gamma9 * tvec[i] + gamma10 * (tvec[i])^2 for i ∈ axes(data, 1) ]
    y2pred = exp.(y2predod) ./ (1 .+ exp.(y2predod))

    

    for (i, t) ∈ enumerate(tvec)
        t <= 40 && continue 
        Turing.@addlogprob! logpdf(Normal(recordedy2[i], sigma2), y2predod[i])
    end
end

recordedy2 = log.(1e-10 .+ unboostedsimulation.y2 ./ (1 .- unboostedsimulation.y2))


model = fitmodel(unboostedsimulation, recordedy2)
chain = sample(model, NUTS(0.65), 1_000)





function fitmodel_target(data, recordedy2)
    return Pigeons.TuringLogPotential(fitmodel(data, recordedy2))
end

recordedy2 = log.(1e-10 .+ unboostedsimulation.y2 ./ (1 .- unboostedsimulation.y2))

const FitmodelType_unboosted = typeof(fitmodel_target(unboostedsimulation, recordedy2))

function Pigeons.initialization(target::FitmodelType_unboosted, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    #Pigeons.update_state!(result, :beta0, 1, 0.0)
    #Pigeons.update_state!(result, :beta1, 1, 0.0)
    #Pigeons.update_state!(result, :beta2, 1, 0.0)
    #Pigeons.update_state!(result, :beta3 , 1, 0.0)

    Pigeons.update_state!(result, :gamma0, 1, 0.0)
    Pigeons.update_state!(result, :gamma1, 1, 0.0)
    Pigeons.update_state!(result, :gamma2, 1, 0.0)
    Pigeons.update_state!(result, :gamma3, 1, 0.0)
    Pigeons.update_state!(result, :gamma4, 1, 0.0)
    Pigeons.update_state!(result, :gamma5, 1, 0.0)
    Pigeons.update_state!(result, :gamma6, 1, 0.0)
    Pigeons.update_state!(result, :gamma7, 1, 0.0)
    Pigeons.update_state!(result, :gamma8, 1, 0.0)
    Pigeons.update_state!(result, :gamma9, 1, 0.0)
    Pigeons.update_state!(result, :gamma10, 1, 0.0)

    Pigeons.update_state!(result, :sigma2 , 1, 0.5)

    return result
end

unboostedfitted_pt = pigeons( ;
    target=fitmodel_target(unboostedsimulation, recordedy2), 
    n_rounds=2,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(1),
    variational=GaussianReference(),
)

unboostedfitted_chain = Chains(unboostedfitted_pt)


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

function immunewaning(
    t::AbstractVector, locationcodes::AbstractVector, infections::AbstractVector, 
    beta2, beta1, betahalf, betazero, betaminus1
)
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

function immunewaning(
    data::DataFrame, t, locationcodes, infections, 
    beta2, beta1, betahalf, betazero, betaminus1
)
    return immunewaning(
        getproperty(data, t), getproperty(data, locationcodes), getproperty(data, infections), 
        beta2, beta1, betahalf, betazero, betaminus1
    )
end

function addimmunewaning!(
    data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1
)
    insertcols!(
        data, 
        :immune => immunewaning(
            data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1
        )
    )
end

function replaceimmunewaning!(
    data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1
)
    data.immune = immunewaning(
        data, t, locationcodes, infections, beta2, beta1, betahalf, betazero, betaminus1
    )
end

#using FillArrays

@model bcodevalue(prior) = b ~ prior

@model function q3fitmodel(prevalence, data)
    # function for immune waning 
    println("check1")
    ℓ_codes = length(unique(data.Code))
    beta2  ~ Normal(0, 1)
    beta1  ~ Normal(0, 1)
    betahalf ~ Normal(0, 1)
    betazero ~ Normal(0, 1) 
    betaminus1 ~ Normal(0, 1)
    println("check2")

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
    println("check3")

    replaceimmunewaning!(data, :t, :Code, :i4, beta2, beta1, betahalf, betazero, betaminus1)
    println("check4")

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
    println("check5")

    println("length(prevalence) = $(length(prevalence)) == length(q3regr) = $(length(q3regr)) : $(length(data.i4) == length(q3regr))")

    #println(q3regr)

    println(describe(q3regr))

    #datai4 = data.i4

    println("sigma2=$sigma2")
    #datai4 ~ MvNormal(q3regr, sigma2 * I)
    for i ∈ eachindex(prevalence)
        prevalence[i] ~ Normal(q3regr[i], sigma2)
    end
end

function pol_q3fitmodel(config)
    println("using pol_q3fitmodel")
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
    println("using pol_q3fitmodel1")

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

function pol_q3fitmodelelement(config)
    println("using pol_q3fitmodelelement")

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




#insertcols!(unboostedsimulation, :immune => immune)


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

#unboostedq3model = q3fitmodel(dc)

# temp 
n_rounds = 3; id = 1

unboostedq3model = q3fitmodel(dc.i4, dc)

fitted_pt = pigeons( ;
    target=TuringLogPotential(unboostedq3model),
    n_rounds=1,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=1,
    variational=GaussianReference(),
)

fitted_chains = Chains(fitted_pt)

asdf_pt = pigeons( ;
    target=TuringLogPotential(unboostedq3model),
    n_rounds=2,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=1,
    variational=GaussianReference(),
)

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

=#
###########################################################################################
# UK data
###########################################################################################

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
#=
insertcols!(
    hospitaldata,
    :VolumePerBed => hospitaldata.HeatedVolume ./ hospitaldata.TotalBeds,
    :ProportionSingleBeds => min.(hospitaldata.SingleBedsEnsuite ./ hospitaldata.TotalBeds, 1.0)
)
=#
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



#=
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
#=
# remove all data before 2021-01-09, when staff largely unvaccinated, and for the next 40 days
Date("2021-01-08") + Day(40)
#2021-02-17
filter!(:Date => x -> x >= Date("2021-02-17"), coviddata)
# and for the next 0
=#
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


=#


###

function simulationproportions!(data; I=:I, Y=:Y, M=:M, N=:N)
    insertcols!(data, :PatientsProportion => getproperty(data, Y) ./ getproperty(data, M))
    insertcols!(data, :StaffProportion => getproperty(data, I) ./ getproperty(data, N))
    #absences = _simulateabsences(data)
    #insertcols!(data, :AbsencesProportion => absences ./ data.N)
end 



coviddata = CSV.read(datadir("exp_raw", "dataset.csv"), DataFrame)

# make hospital codes categorical variables 
rename!(coviddata, :Codes => "StringCodes")
insertcols!(coviddata, 1, :Code => CategoricalArray(coviddata.StringCodes))

# calculate values for the analysis
simulationproportions!(
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


nhospitals = length(unique(coviddata.Code))
ndates = length(unique(coviddata.Date))

patients = zeros(ndates, nhospitals)
staff = zeros(ndates, nhospitals) 
newstaff = zeros(ndates, nhospitals) 
for (i, c) ∈ enumerate(unique(coviddata.Code))
    _tdf = filter(:Code => x -> x == c, coviddata)
    for t ∈ 1:796
        ind = findfirst(x -> x == t, _tdf.t) 
        if !isnothing(ind)
            patients[t, i] = _tdf.PatientsProportion[ind]
            staff[t, i] = _tdf.StaffProportion[ind]
            if t == 1 
                newstaff[t, i] = staff[t, i] / 10
            elseif t <= 10 
                newstaff[t, i] =  max(
                    0,
                    min(
                        1,
                        staff[t, i] * t / 10 - sum(@view newstaff[1:(t - 1), i])
                    )

                )
            else
                newstaff[t, i] =  max(
                    0,
                    min(
                        1,
                        staff[t, i] - sum(@view newstaff[(t - 10):(t - 1), i])
                    )

                )
            end
        end
#    println(size(_tdf))
 #   patients[:, i] .= _tdf.PatientsProportion
  #  staff[:, i] .= _tdf.StaffProportion
    end
end

vpd, psb = let 
    _tdf = select(coviddata, :Code, :VolumePerBed, :ProportionSingleBeds)
    unique!(_tdf)
    ( _tdf.VolumePerBed, _tdf.ProportionSingleBeds )
end

stringency = coviddata.StringencyIndex_Average[1:ndates]

weeklycases = coviddata.weeklycases[1:ndates]

vaccinated = zeros(ndates, nhospitals)
for d ∈ axes(vaccinated, 1), h ∈ axes(vaccinated, 2)
    if d ∈ [ 300, 450, 650, 800 ]
        vaccinated[d, h] = 0.8
    end
end
#=
α1 = 0.02
α2 = 0.05
α3 = 0.02
α4 = 0.6
α5 = 0.02
α6 = 0.005
α7 = 0.02
α8 = 0.005
α9 = 2 
α10 = 2.5
α11 = 0.3 
α12 = 500 
α13 = 1
α14 = 0.3 
α15 = 500
α16 = 0.5
α17 = 0.5

βp = @. α1 + α2 * vpd + α3 * psb
βh = @. α4 + α5 * vpd + α6 * psb
βc = @. α7 + α8 * (100 - stringency)

foi = zeros(ndates, nhospitals)
for d ∈ axes(foi, 1), h ∈ axes(foi, 2)
    foi[d, h] = βp[h]* patients[d, h] + βh[h]* staff[d, h] + βc[d] * weeklycases[d] / 56_000_000
end

wane_noboost = zeros(ndates, nhospitals)
for d ∈ axes(wane_noboost, 1), h ∈ axes(wane_noboost, 2)
    d == 1 && continue
    wane_noboost[d, h] = sum([ (staff[x, h] + vaccinated[x, h]) * pdf(Weibull(α11, α12), d - x) for x ∈ 1:(d - 1) ])
end

delaystaff = zeros(ndates, nhospitals)
immunestaff = zeros(ndates, nhospitals)
wane_boost = zeros(ndates, nhospitals)
for d ∈ axes(wane_boost, 1), h ∈ axes(wane_boost, 2)
    delaystaff[d, h] = staff[d, h] + vaccinated[d, h]
    immunestaff[d, h] = staff[d, h] + vaccinated[d, h]
    d == 1 && continue

    for d_2 ∈ 1:(d-1)
        de = min(α13 * foi[d_2, h], 1) * immunestaff[d_2, h]
        delaystaff[d, h] += de
        immunestaff[d, h] += de 
        immunestaff[d_2, h] += -de 
    end
    wane = 0 
    for d_2 ∈ 1:(d-1)
        newwane = immunestaff[d_2, h] * pdf(Weibull(α14, α15), d - d_2)
        wane += newwane
        immunestaff[d_2, h] += -newwane 
    end
    wane_boost[d, h] = wane 
end
    
susceptible = ones(ndates, nhospitals)
for d ∈ axes(susceptible, 1), h ∈ axes(susceptible, 2)
    d == 1 && continue
    susceptible[d, h] = max(
        0,
        min(
            1,
            susceptible[(d - 1), h] - staff[(d - 1), h] + α16 * wane_noboost[(d - 1), h] + α17 * wane_boost[(d - 1), h]
        )
    )
end

pred = foi .* susceptible
=#

@model function fitmodel(newstaff, patients, staff, vaccinated, weeklycases, vpd, psb, stringency, ndates, nhospitals)
    α1 ~ Beta(1, 1)
    α2 ~ Exponential(1)
    α3 ~ Exponential(1)
    α4 ~ Beta(1, 1)
    α5 ~ Exponential(1)
    α6 ~ Exponential(1)
    α7 ~ Beta(1, 1)
    α8 ~ Exponential(1)

    a1 ~ Exponential(0.3) 
    θ1 ~ Exponential(500) 
    ψ ~ Exponential(1)
    a2 ~ Exponential(0.3)
    θ2 ~ Exponential(500)
    α ~ Beta(1, 1)
    β ~ Beta(1, 1)

    sigma2 ~ Exponential(1)

    βp = [ α1 + α2 * v + α3 * p for(v, p) ∈ zip(vpd, psb) ]
    βh = [ α4 + α5 * v + α6 * p for(v, p) ∈ zip(vpd, psb) ]
    βc = @. α7 + α8 * (100 - stringency)

    foi = zeros(ndates, nhospitals)
    for d ∈ axes(foi, 1), h ∈ axes(foi, 2)
        foi[d, h] = βp[h]* patients[d, h] + βh[h]* staff[d, h] + βc[d] * weeklycases[d] / 56_000_000
    end

    wane_noboost = zeros(ndates, nhospitals)
    for d ∈ axes(wane_noboost, 1), h ∈ axes(wane_noboost, 2)
        d == 1 && continue
        wane_noboost[d, h] = sum([ (staff[x, h] + vaccinated[x, h]) * pdf(Weibull(a1, θ1), d - x) for x ∈ 1:(d - 1) ])
    end

    delaystaff = zeros(ndates, nhospitals)
    immunestaff = zeros(ndates, nhospitals)
    wane_boost = zeros(ndates, nhospitals)
    for d ∈ axes(wane_boost, 1), h ∈ axes(wane_boost, 2)
        delaystaff[d, h] = staff[d, h] + vaccinated[d, h]
        immunestaff[d, h] = staff[d, h] + vaccinated[d, h]
        d == 1 && continue

        for d_2 ∈ 1:(d-1)
            de = min(ψ * foi[d_2, h], 1) * immunestaff[d_2, h]
            delaystaff[d, h] += de
            immunestaff[d, h] += de 
            immunestaff[d_2, h] += -de 
        end
        wane = 0 
        for d_2 ∈ 1:(d-1)
            newwane = immunestaff[d_2, h] * pdf(Weibull(a2, θ2), d - d_2)
            wane += newwane
            immunestaff[d_2, h] += -newwane 
        end
        wane_boost[d, h] = wane 
    end
        
    susceptible = ones(ndates, nhospitals)
    for d ∈ axes(susceptible, 1), h ∈ axes(susceptible, 2)
        d == 1 && continue
        susceptible[d, h] = max(
            0,
            min(
                1,
                susceptible[(d - 1), h] - staff[(d - 1), h] + α * wane_noboost[(d - 1), h] + β * wane_boost[(d - 1), h]
            )
        )
    end

    pred = min.(foi .* susceptible, 1)

    for d ∈ axes(pred, 1), h ∈ axes(pred, 2)
        d <= 14 && continue 
        Turing.@addlogprob! logpdf(Normal(newstaff[d, h], sigma2), pred[d, h])
    end
end

model = fitmodel(newstaff, patients, staff, vaccinated, weeklycases, vpd, psb, stringency, ndates, nhospitals)
chain = sample(model, NUTS(0.65), 3)

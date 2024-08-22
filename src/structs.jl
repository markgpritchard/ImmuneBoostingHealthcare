
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extra functions 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

minimum(p::SEIIRRRSp) = minimum([ p.β, p.γ, p.ε, p.rho, p.ω ])

function minimum(p::WXYYZSEIIRRRSp)
    a = minimum(p.α)
    b = minimum(p.β)
    d = minimum(p.δ)
    return minimum([ a, b, p.γ, d, p.ε, p.λc, p.rho, p.ω ])
end 




# for general page 
import Base: ==, HasLength, iterate, length, showerror




################################

abstract type AbstractModelOutputs{T} end 

mutable struct HCWSEIRRRVOutput{T} <: AbstractModelOutputs{T}
    t       :: Int  # time
    S       :: T
    V       :: T
    E       :: T
    I       :: T
    I′1     :: T
    I′2     :: T
    I′3     :: T
    I′4     :: T
    I′5     :: T
    I′6     :: T
    I′7     :: T
    I′8     :: T
    I′9     :: T
    I′10    :: T
    R1      :: T
    R2      :: T
    R3      :: T
end

mutable struct HCWSEIRRRVCPOutput{T} <: AbstractModelOutputs{T}
    t       :: Int  # time
    S       :: T
    V       :: T
    E       :: T
    I       :: T
    I′1     :: T
    I′2     :: T
    I′3     :: T
    I′4     :: T
    I′5     :: T
    I′6     :: T
    I′7     :: T
    I′8     :: T
    I′9     :: T
    I′10    :: T
    R1      :: T
    R2      :: T
    R3      :: T
    I_c     :: Int 
    N_c     :: Int 
    I_p     :: Int 
    N_p     :: Int
end

HCWSEIRRRVOutput() = HCWSEIRRRVOutput(0, zeros(17)...)
HCWSEIRRRVOutput(S::T) where T = HCWSEIRRRVOutput(0, S, zeros(T, 16)...)

abstract type AbstractParameters{S, T} end 

struct HCWSEIRRRParameters{S, T} <: AbstractParameters{S, T}
    beta_c  :: S 
    beta_h  :: S 
    beta_p  :: S 
    gamma   :: S 
    sigma   :: S 
    theta   :: S 
    psi     :: S 
    nu      :: T 
    eta     :: S 
    omega   :: S
end

function HCWSEIRRRParameters(; beta_c=0.0, beta_h=0.0, beta_p=0.0, gamma=0.0, sigma=0.0, theta=0.0, psi=0.0, nu=0.0, eta=0.0, omega=0.0)
    return HCWSEIRRRParameters(beta_c, beta_h, beta_p, gamma, sigma, theta, psi, nu, eta, omega)
end

gettime(a::AbstractModelOutputs) = a.t
get_n(a::AbstractModelOutputs) = sum([ x for x ∈ a ])
get_S(a::AbstractModelOutputs) = a.S
get_E(a::AbstractModelOutputs) = a.E
get_I(a::AbstractModelOutputs) = a.I
get_v(a::AbstractModelOutputs) = a.V
get_R1(a::AbstractModelOutputs) = a.R1
get_R2(a::AbstractModelOutputs) = a.R2
get_R3(a::AbstractModelOutputs) = a.R3

countnewlydiagnosed(a::AbstractModelOutputs) = a.I′1
gettotaldiagnosed(a::AbstractModelOutputs) = sum(_diagnosedvector(a))
gettotalimmune(a::AbstractModelOutputs) = sum([ getproperty(a, Symbol("R$i")) for i ∈ 1:3 ])

function ==(a::AbstractModelOutputs, b::AbstractModelOutputs)
    for (x, y) ∈ zip(a, b)
        x != y && return false
    end
    return true
end

function iterate(a::HCWSEIRRRVOutput, state=0)
    state >= 17 && return nothing 
    return Base.getfield(a, state + 2), state + 1 
end

Base.HasLength(::HCWSEIRRRVOutput) = 17
length(::HCWSEIRRRVOutput) = 17

function iterate(a::HCWSEIRRRParameters, state=0)
    state >= 10 && return nothing 
    return Base.getfield(a, state + 1), state + 1 
end

Base.HasLength(::HCWSEIRRRParameters) = 10
length(::HCWSEIRRRParameters) = 10


function advancetime!(a::AbstractModelOutputs) 
    a.t += 1
    dv = _diagnosedvector(a)
    a.I′1 = 0 
    for i ∈ 1:9 
        setproperty!(a, Symbol("I′$(i + 1)"), dv[i])
    end   
    a.R1 += dv[10] 
end

function _moveindividuals!(a::AbstractModelOutputs, from, to, x)
    setproperty!(a, from, getproperty(a, from) -x)
    setproperty!(a, to, getproperty(a, to) +x)
    getproperty(a, from) < 0 && _negativecompartmentwarning(a, from)
end

expose!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :S, :E, x)
progress!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :E, :I, x)

function progress!(a::AbstractModelOutputs, p::HCWSEIRRRParameters)
    _parameterdomain(p::HCWSEIRRRParameters)
    progress!(a, a.E * p.sigma)
end

diagnose!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :I, :I′1, x) 

function diagnose!(a::AbstractModelOutputs, p::HCWSEIRRRParameters)
    _parameterdomain(p)
    diagnose!(a, a.I * p.theta)
end 

recover!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :I, :R1, x) 
recover!(a::AbstractModelOutputs, p::HCWSEIRRRParameters) = recover!(a, a.I * p.gamma * (1 - p.theta)) 
wane1!(a::AbstractModelOutputs, x) = _wane!(a, :R1, :R2, x) 
wane2!(a::AbstractModelOutputs, x) = _wane!(a, :R2, :R3, x) 
wane3!(a::AbstractModelOutputs, x) = _wane!(a, :R3, :S, x) 
vaccinateS!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :S, :V, x) 

function vaccinateS!(a::HCWSEIRRRVOutput, p::HCWSEIRRRParameters{S,T}) where {S <: Any, T <: Number}
    _parameterdomain(p)
    vaccinateS!(a, a.S * p.nu)
end 

function vaccinateS!(a::HCWSEIRRRVOutput, p::HCWSEIRRRParameters{S,T}, t) where {S <: Any, T <: Function}
    _parameterdomain(p, t)
    vaccinateS!(a, a.S * p.nu(t))
end 

becomeimmunefromvaccine!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :V, :R1, x) 
exposevaccinated!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :V, :E, x) 
vaccinateR2!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R2, :R1, x) 
vaccinateR3!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R3, :R1, x) 
boostR2!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R2, :R1, x) 
boostR3!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R3, :R1, x) 

_wane!(a::AbstractModelOutputs, from, to, x::Number) = _moveindividuals!(a, from, to, x)

function _wane!(a::AbstractModelOutputs, from, to, p::HCWSEIRRRParameters)
    _parameterdomain(p)
    _wane!(a, from, to, getproperty(a, from) * 3 * p.omega)
end

function _negativecompartmentwarning(a, from) 
    @warn "Negative value in compartment $from, $(getproperty(a, from)), at time $(gettime(a))"
end

_diagnosedvector(a::AbstractModelOutputs) = [ getproperty(a, Symbol("I′$i")) for i ∈ 1:10 ]

function _parametermin(p::HCWSEIRRRParameters{S, T}) where {S <: Any, T <: Number}
    return minimum([ x for x ∈ p ])
end

function _parametermin(p::HCWSEIRRRParameters{S, T}) where {S <: Any, T <: Function}
    return minimum(
        [ 
            getproperty(p, a) 
            for a ∈ [ :beta_c, :beta_h, :beta_p, :gamma, :sigma, :theta, :psi, :eta, :omega ] 
        ]
    )
end

function _parameterdomain(p::HCWSEIRRRParameters, t=nothing) 
    _parametermin(p) < 0 && throw(DomainError(p, "all parameters must be non-negative"))
    p.gamma > 1 && throw(DomainError(p.gamma, "gamma must be between 0 and 1"))
    p.sigma > 1 && throw(DomainError(p.sigma, "sigma must be between 0 and 1"))
    p.theta > 1 && throw(DomainError(p.theta, "theta must be between 0 and 1"))
    p.eta > 1 && throw(DomainError(p.eta, "eta must be between 0 and 1"))
    p.omega > 1 / 3 && throw(DomainError(p.omega, "omega must be between 0 and 1/3"))
    _nuparameterdomain(p, t)
end

function _nuparameterdomain(p::HCWSEIRRRParameters{S, T}, t) where {S <: Any, T <: Number}
    p.nu > 1 && throw(DomainError(p.nu, "nu must be between 0 and 1"))
end

function _nuparameterdomain(p::HCWSEIRRRParameters{S, T}, t) where {S <: Any, T <: Function}
    p.nu(t) < 0 && throw(DomainError(p.nu(t), "nu must be between 0 and 1 (error at time $t)"))
    p.nu(t) > 1 && throw(DomainError(p.nu(t), "nu must be between 0 and 1 (error at time $t)"))
end

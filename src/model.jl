


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

HCWSEIRRRVCPOutput(; I_c=0, N_c=0, I_p=0, N_p=0) = HCWSEIRRRVCPOutput(0, zeros(17)..., I_c, N_c, I_p, N_p)
HCWSEIRRRVCPOutput(S::T; I_c=0, N_c=0, I_p=0, N_p=0)  where T = HCWSEIRRRVCPOutput(0, S, zeros(T, 16)..., I_c, N_c, I_p, N_p)


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
get_n(a::AbstractModelOutputs) = sum(@view [ x for x ∈ a ][2:18])
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
    state >= 18 && return nothing 
    return Base.getfield(a, state + 1), state + 1 
end

Base.HasLength(::HCWSEIRRRVOutput) = 18
length(::HCWSEIRRRVOutput) = 18

function iterate(a::HCWSEIRRRVCPOutput, state=0)
    state >= 22 && return nothing 
    return Base.getfield(a, state + 1), state + 1 
end

Base.HasLength(::HCWSEIRRRVCPOutput) = 22
length(::HCWSEIRRRVCPOutput) = 22

function iterate(a::HCWSEIRRRParameters, state=0)
    state >= 10 && return nothing 
    return Base.getfield(a, state + 1), state + 1 
end

Base.HasLength(::HCWSEIRRRParameters) = 10
length(::HCWSEIRRRParameters) = 10

function runmodel(a, p)
    nextgen = deepcopy(a)
    runmodel!(nextgen, a, p)
    return nextgen
end

function runmodel!(nextgen, a, p)
    advancetime!(a)
end


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

function expose!(a::AbstractModelOutputs, p::HCWSEIRRRParameters)
    _parameterdomain(p::HCWSEIRRRParameters)
    _moveindividuals!(a, :S, :E, _numbertoinfect(a, p))
end

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
wane1!(a::AbstractModelOutputs, x, t=nothing) = _wane!(a, :R1, :R2, x, t) 
wane2!(a::AbstractModelOutputs, x, t=nothing) = _wane!(a, :R2, :R3, x, t) 
wane3!(a::AbstractModelOutputs, x, t=nothing) = _wane!(a, :R3, :S, x, t) 
vaccinateS!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :S, :V, x) 

function vaccinateS!(a::AbstractModelOutputs, p::HCWSEIRRRParameters{S,T}) where {S <: Any, T <: Number}
    _parameterdomain(p)
    lambda = calculatelambda(a, p)
    vaccinateS!(a, a.S * p.nu * (1 - lambda))
end 

function vaccinateS!(a::AbstractModelOutputs, p::HCWSEIRRRParameters{S,T}, t) where {S <: Any, T <: Function}
    _parameterdomain(p, t)
    lambda = calculatelambda(a, p)
    vaccinateS!(a, a.S * p.nu(t))
end 

becomeimmunefromvaccine!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :V, :R1, x) 
exposevaccinated!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :V, :E, x) 
vaccinateR2!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R2, :R1, x) 
vaccinateR3!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R3, :R1, x) 
boostR2!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R2, :R1, x) 
boostR3!(a::AbstractModelOutputs, x) = _moveindividuals!(a, :R3, :R1, x) 

function calculatelambdaprime_c(a::HCWSEIRRRVCPOutput, p::HCWSEIRRRParameters)
    a.I_c < 0 && throw(DomainError(a.I_c, "number of infections in the community cannot be negative"))
    a.N_c < 0 && throw(DomainError(a.I_c, "community population cannot be negative"))
    return p.beta_c * a.I_c / a.N_c
end

function calculatelambdaprime_p(a::HCWSEIRRRVCPOutput, p::HCWSEIRRRParameters)
    a.I_p < 0 && throw(DomainError(a.I_p, "number of infections among patients cannot be negative"))
    a.N_p < 0 && throw(DomainError(a.I_p, "number of patients cannot be negative"))
    return p.beta_p * a.I_p / a.N_p
end

function calculatelambdaprime_h(a::AbstractModelOutputs, p::HCWSEIRRRParameters)
    return p.beta_h * get_I(a) / get_n(a)
end

calculatelambdaprime(a::HCWSEIRRRVOutput, p) = calculatelambdaprime_h(a, p)

function calculatelambdaprime(a::HCWSEIRRRVCPOutput, p)
    return calculatelambdaprime_c(a, p) + calculatelambdaprime_p(a, p) + calculatelambdaprime_h(a, p)
end

function calculatelambda(a, p)
    return 1 - exp(-calculatelambdaprime(a, p))
end

function calculateboosting(a, p)
    p.psi == 0 && return zero(1 - exp(-p.psi * calculatelambdaprime(a, p)))
    return 1 - exp(-p.psi * calculatelambdaprime(a, p))
end

_wane!(a::AbstractModelOutputs, from, to, x::Number, ::Nothing) = _moveindividuals!(a, from, to, x)

function _wane!(a::AbstractModelOutputs, from, to, p::HCWSEIRRRParameters{S, T}, ::Nothing) where {S <: Any, T <: Number}
    _parameterdomain(p)
    boost = calculateboosting(a, p)
    _wane!(a, from, to, getproperty(a, from) * 3 * p.omega * (1 - p.nu) * (1 - boost), nothing)
end

function _wane!(a::AbstractModelOutputs, from, to, p::HCWSEIRRRParameters{S, T}, t) where {S <: Any, T <: Function}
    _parameterdomain(p)
    boost = calculateboosting(a, p)
    _wane!(a, from, to, getproperty(a, from) * 3 * p.omega * (1 - p.nu(t)) * (1 - boost), nothing)
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

function _numbertoinfect(a, p)
    lambda = calculatelambda(a, p)
    infect = a.S * lambda 
    return infect 
end
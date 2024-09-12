
struct Automatic end  # not exported
automatic = Automatic()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters for simulation models 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abstract type AbstractParameters{S, T, U} end 

struct SEIIRRRSp{S, T, U} <: AbstractParameters{S, T, U} 
    β       :: S                # infection rate 
    rho     :: T                # rate of leaving exposed compartments 
    γ       :: T                # rate of leaving infectious compartment
    ψ       :: T                # strength of "force of boosting" relative to λ
    nu      :: U                # vaccination rate
    ω       :: T                # rate of immune waning 
end

struct WXYYZSEIIRRRSp{S, T, U} <: AbstractParameters{S, T, U} 
    βhh     :: T                # rate of infection from healthcare worker to healthcare worker 
    βhp     :: T                # rate of infection from patient to healthcare worker    
    βph     :: T                # rate of infection from healthcare worker to patient  
    βpp     :: T                # rate of infection from patient to patient  
    rho     :: T                # rate of leaving exposed compartments 
    γ       :: T                # rate of leaving infectious compartment
    ψ       :: T                # strength of "force of boosting" relative to λ
    nu      :: U                # healthcare worker vaccination rate
    ω       :: T                # rate of immune waning 
    δn      :: T                # discharge rate of non-infected 
    δi      :: T                # discharge rate of infected  
    λc      :: S                # community force of infection  
end

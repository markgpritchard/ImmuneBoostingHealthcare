
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CategoricalArrays, CSV, DataFrames, Dates, Distributions, Random, Turing, Pigeons




#############################################


##############################################


#chain = sample(model, NUTS(0.65), 8)


include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare

using Turing, Pigeons

if length(ARGS) == 2 
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
else
    id = 1 
    n_rounds = 10
end

nhospitals = counthospitals(coviddata)
ndates = countdates(coviddata)

@unpack newstaff, patients, staff = datamatrices(coviddata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(coviddata)

stringency = coviddata.StringencyIndex_Average[1:ndates]
weeklycases = coviddata.weeklycases[1:ndates]

# numbers vaccinated currently simulated 
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for d ∈ axes(vaccinated, 1), h ∈ axes(vaccinated, 2)
        if d ∈ [ 300, 450, 650, 800 ]
            vaccinated[d, h] = 0.8
        end
    end
    vaccinated
end

rjt = zeros(Float64, ndates)  # levels of immunity


βp = [ .4 + .3 * v + .9 * p for(v, p) ∈ zip(vpd, psb) ]
βh = [ .6 + .33 * v + 1 * p for(v, p) ∈ zip(vpd, psb) ]
βc = @. .9 + .04 * (100 - stringency)
foi = zeros(Float64, ndates, nhospitals)
for d ∈ axes(foi, 1), h ∈ axes(foi, 2)
    foi[d, h] = βp[h] * patients[d, h] + βh[h]* staff[d, h] + βc[d] * weeklycases[d] / 56_000_000
end
a = 0.04; θ = 500.7

function predictinfections!(T, rjt, ndates, nhospitals, foi, vaccinated, ψ, a, θ)
    predictedinfections = zeros(T, ndates, nhospitals)
    for h ∈ 1:nhospitals
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  # reset the vector rjt 
            rjt1 = predictedinfections[(d - 1), h] + 
                vaccinated[(d - 1), h] * (1 - sum(rjt)) + 
                (1 - exp(-ψ * foi[(d - 1), h])) * sum(rjt)
            for x ∈ d:-1:2
                rjt[x] = rjt[x-1] * 
                    exp(-ψ * foi[(d - 1), h]) * 
                    (1 - pdf(Weibull(a, θ), x) / (ccdf(Weibull(a, θ), x - 1)))
            end
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
    return predictedinfections
end

predictedinfections = zeros(Float64, ndates, nhospitals)


function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    for h ∈ 1:nhospitals
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        predictedinfections[1, h] = zero(T)
        for d ∈ 2:ndates  # reset the vector rjt 
            rjt1 = predictedinfections[(d - 1), h] + 
                vaccinated[(d - 1), h] * (1 - sum(rjt)) + 
                (1 - exp(-ψ * foi[(d - 1), h])) * sum(rjt)
            for x ∈ d:-1:2
                rjt[x] = rjt[x-1] * 
                    exp(-ψ * foi[(d - 1), h]) * 
                    (1 - pdf(Weibull(a, θ), x) / (ccdf(Weibull(a, θ), x - 1)))
            end
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
    return predictedinfections
end

##

function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    #predictedinfections[1, :] .= zeros(T, nhospitals)
    for h ∈ 1:nhospitals
        @assert predictedinfections[1, h] == 0
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  
            v = vaccinated[(d - 1), h]
            rjt1 = predictedinfections[(d - 1), h] + 
                #vaccinated[(d - 1), h] * (1 - sum(rjt)) #+ 
                v * (1 - sum(rjt)) #+ 
                #(1 - exp(-ψ * foi[(d - 1), h])) * sum(rjt)
            for x ∈ d:-1:2
                b = (1 - exp(-ψ * foi[(d - 1), h]))
                #b′ = b + (1 - b) * vaccinated[(d - 1), h]
                b′ = b + (1 - b) * v
                rjt1 += b′ * rjt[x-1]
                rjt[x] = (1 - b′) * rjt[x-1] * 
                   # exp(-ψ * foi[(d - 1), h]) * 
                    (1 - pdf(Weibull(a, θ), x) / (ccdf(Weibull(a, θ), x - 1)))
            end
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
    return predictedinfections
end

function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    for h ∈ 1:nhospitals
        @assert predictedinfections[1, h] == 0
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  
            v = vaccinated[(d - 1), h]
            b = (1 - exp(-ψ * foi[(d - 1), h])) * (1 - v) + v
            rjt1 = predictedinfections[(d - 1), h] + v * (1 - sum(rjt)) + b * sum(rjt)
            for x ∈ d:-1:2
                #rjt1 += b * rjt[x-1]
                rjt[x] = (1 - b) * rjt[x-1] * 
                    (1 - pdf(Weibull(a, θ), x) / ccdf(Weibull(a, θ), x - 1))
            end
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
end

function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    for h ∈ 1:nhospitals
        @assert predictedinfections[1, h] == 0
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  
            v = vaccinated[(d - 1), h]
            b = (1 - exp(-ψ * foi[(d - 1), h])) * (1 - v) + v
            rjt1 = predictedinfections[(d - 1), h] + v * (1 - sum(rjt)) + b * sum(rjt)
            for x ∈ d:-1:2
                #rjt1 += b * rjt[x-1]
                rjt[x] = (1 - b) * rjt[x-1] * 
                    (1 - pdf(Weibull(a, θ), x) / ccdf(Weibull(a, θ), x - 1))
            end
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
end
#=
function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    for h ∈ 1:nhospitals
        _predictinfections!(@view predictedinfections[:, h], rjt, ndates, @view foi[:, h], @view vaccinated[:, h], ψ, a, θ)
    end
end

function _predictinfections!(predictedinfections::SubArray, rjt::Vector{T}, ndates, foi::SubArray, vaccinated::SubArray, ψ, a, θ) where T
    @assert predictedinfections[1] == 0
    for d ∈ 1:ndates  # reset the vector rjt 
        rjt[d] = zero(T)
    end
    for d ∈ 2:ndates  
        v = vaccinated[d-1]
        b = (1 - exp(-ψ * foi[d-1])) * (1 - v) + v
        rjt1 = predictedinfections[d-1] + v * (1 - sum(rjt)) + b * sum(rjt)
        for x ∈ d:-1:2
            #rjt1 += b * rjt[x-1]
            rjt[x] = (1 - b) * rjt[x-1] * 
                (1 - pdf(Weibull(a, θ), x) / ccdf(Weibull(a, θ), x - 1))
        end
        rjt[1] = rjt1 
        predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
    end
end
=#



#=
function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    for h ∈ 1:nhospitals
        @assert predictedinfections[1, h] == 0
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  
            v = vaccinated[(d - 1), h]
            b = (1 - exp(-ψ * foi[(d - 1), h])) * (1 - v) + v
            carry0 = zero(T)
            carry1 = rjt[1]  # held for the next step 
            rjt[1] = predictedinfections[(d - 1), h] + v * (1 - sum(rjt)) + b * sum(rjt)
            for x ∈ 2:d 
                carry0 = carry1  # the previous value of rjt[x-1]
                carry1 = rjt[x]  # held for the next step
                rjt[x] = (1 - b) * carry0 * 
                    (1 - pdf(Weibull(a, θ), x) / ccdf(Weibull(a, θ), x - 1))
            end
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
end
=#
@benchmark predictinfections!(Float64, rjt, ndates, nhospitals, foi, vaccinated, 2.0, 0.04, 500.7)

function predictinfections!(predictedinfections::Matrix{T}, rjt::Vector{T}, ndates, nhospitals, foi, vaccinated, ψ, a, θ) where T
    for h ∈ 1:nhospitals
        @assert predictedinfections[1, h] == 0
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  
            v = vaccinated[(d - 1), h]
            b = (1 - exp(-ψ * foi[(d - 1), h])) * (1 - v) + v
            rjt1 = predictedinfections[(d - 1), h] + v * (1 - sum(rjt)) + b * sum(rjt)
            @view(rjt[2:d]) .+= [ 
                -rjt[x] + (1 - b) * rjt[x-1] * (1 - pdf(Weibull(a, θ), x) / ccdf(Weibull(a, θ), x - 1))
                for x ∈ 2:d
            ]
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end
end




asfd = [ 1., 2, 3, 4, 5, 6, 7 ]

carry0 = zero(asfd[1])
carry1 = asfd[1] 
asfd[1] = 0 
for i ∈ 2:7 
    carry0 = carry1 
    carry1 = asfd[i] 
    asfd[i] = carry0 
end

asfd

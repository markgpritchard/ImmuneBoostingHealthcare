
counthospitals(data) = length(unique(data.Code))
countdates(data; dateid=:Date) = length(unique(getproperty(data, dateid))) 

function datamatrices(data, ndates, nhospitals)
    patients = zeros(ndates, nhospitals)
    staff = zeros(ndates, nhospitals) 
    newstaff = zeros(ndates, nhospitals) 
    for (i, c) ∈ enumerate(unique(data.Code))
        _tdf = filter(:Code => x -> x == c, data)
        for t ∈ 1:ndates
            patients[t, i] = _tdf.PatientsProportion[t]
            staff[t, i] = _tdf.StaffProportion[t]
            if t == 1 
                newstaff[t, i] = staff[t, i] / 10
            elseif t <= 10 
                newstaff[t, i] = max(
                    0,
                    min(1, staff[t, i] * t / 10 - sum(@view newstaff[1:(t - 1), i]))
                )
            else
                newstaff[t, i] = max(
                    0,
                    min(1, staff[t, i] - sum(@view newstaff[(t - 10):(t - 1), i]))
                )
            end
        end
    end
    return @ntuple newstaff patients staff
end

function hospitalconditionmatrices(data)
    _tdf = select(data, :Code, :VolumePerBed, :ProportionSingleBeds)
    unique!(_tdf)
    return @ntuple vpd=(_tdf.VolumePerBed).^(1/3) psb=_tdf.ProportionSingleBeds 
end

function hcwseiirrr(
    u, p::HCWSEIIRRRp, t::Integer, λc::Number, patients::Number, vaccinated::Number
)
    S, E, I, I′1, I′2, I′3, I′4, I′5, I′6, I′7, I′8, I′9, I′10, R1, R2, R3 = u 

    λ = 1 - exp(-(λc + p.βp * patients + p.βh * I))
    λψ = vaccinated + (1 - exp(-p.ψ * (λc + p.βp * patients + p.βh * I))) * (1 - vaccinated)

    new_u = [
        S * (1 - λ - vaccinated * (1 - λ)) + R3 * 3 * p.ω * (1 - λψ),  # S 
        E * (1 - p.η) + λ * S,  # E 
        I * (1 - p.θ - p.γ * (1 - p.θ)) + p.η * E,  # I 
        p.θ * I,  # I′1
        I′1,  # I′2
        I′2,  # I′3
        I′3,  # I′4
        I′4,  # I′5
        I′5,  # I′6
        I′6,  # I′7
        I′7,  # I′8
        I′8,  # I′9
        I′9,  # I′10
        # R1:
        R1 * (1 - 3 * p.ω * (1 - λψ)) +  # previous R1 that has not waned
            S * (vaccinated * (1 - λ)) +   # vaccinated from S 
            p.γ * (1 - p.θ) * I +  # recovered 
            I′10 +  # ended isolation 
            λψ * (R2 + R3),  # boosted by exposure or vaccination from R2 and R3 
        R2 * (1 - 3 * p.ω * (1 - λψ) - λψ) + R1 * 3 * p.ω * (1 - λψ),  # R3
        R3 * (1 - 3 * p.ω * (1 - λψ) - λψ) + R2 * 3 * p.ω * (1 - λψ),  # R3
    ]
    return new_u
end

function hcwseiirrr_isolating(
    u0, p, tspan::AbstractVector{<:Integer}, λc, patients, vaccinated, j=1
)
    isolating = zeros(Float64, length(tspan))
    hcwseiirrr_isolating!(isolating, u0, p, tspan, λc, patients, vaccinated, j)
    return isolating
end

function hcwseiirrr_isolating!(
    isolating, u0::AbstractVector, p::HCWSEIIRRRp, tspan::AbstractVector{<:Integer}, 
    λc::AbstractVector, patients::AbstractVector, vaccinated::AbstractVector, j
)
    u = u0 
    isolating[1] = sum(@view u[4:13])
    for t ∈ tspan 
        u = hcwseiirrr(u, p, t, λc[t], patients[t], vaccinated[t])
        @assert sum(u) ≈ 1 "sum(u)=$(sum(u))≠1, u=$u, p=$p, t=$t"
        isolating[t] = sum(@view u[4:13])
    end
end

function hcwseiirrr_isolating!(
    isolating, ::Automatic, p::HCWSEIIRRRp, tspan::AbstractVector{<:Integer}, 
    λc::AbstractVector, patients::AbstractVector, vaccinated::AbstractVector, j
)
    u = zeros(16)
    u[1] = 1.0
    isolating[1] = sum(@view u[4:13])
    for t ∈ tspan 
        u = hcwseiirrr(u, p, t, λc[t], patients[t], vaccinated[t])
        @assert sum(u) ≈ 1 "sum(u)=$(sum(u))≠1, u=$u, p=$p, t=$t"
        isolating[t] = sum(@view u[4:13])
    end
end

function hcwseiirrr_isolating!(
    isolating, u0, p::HCWSEIIRRRp, tspan::AbstractVector{<:Integer}, 
    λc, patients, vaccinated::AbstractMatrix, j
)
    hcwseiirrr_isolating!(isolating, u0, p, tspan, λc, patients, view(vaccinated, :, j), j)
end

function hcwseiirrr_isolating!(
    isolating, u0, p::HCWSEIIRRRp, tspan::AbstractVector{<:Integer}, 
    λc, patients::AbstractMatrix, vaccinated::AbstractVector, j
)
    hcwseiirrr_isolating!(isolating, u0, p, tspan, λc, view(patients, :, j), vaccinated, j)
end

function hcwseiirrr_isolating!(
    isolating, u0, p::HCWSEIIRRRp, tspan::AbstractVector{<:Integer}, 
    λc::AbstractMatrix, patients::AbstractVector, vaccinated::AbstractVector, j
)
    hcwseiirrr_isolating!(isolating, u0, p, tspan, view(λc, :, j), patients, vaccinated, j)
end

function hcwseiirrr_isolating!(isolating, u0, p, tspan, λc, patients, vaccinated)
    hcwseiirrr_isolating!(isolating, u0, p, tspan, λc, patients, vaccinated, 1)
end

@model function fitmodel( 
    patients, staff, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals;
    alpha1prior=truncated(Normal(0.2, 1), -1, 10),
    alpha2prior=Normal(0, 0.1),
    alpha3prior=Normal(0, 1),
    alpha4prior=truncated(Normal(0.2, 1), -1, 10),
    alpha5prior=Normal(0, 0.1),
    alpha6prior=Normal(0, 1),
    alpha7prior=truncated(Normal(0.2, 1), -1, 10),
    alpha8prior=truncated(Normal(0.1, 0.1), 0, 10),  # require greater stringency leads to less transmission
    omegaprior=Uniform(0, 0.33),
    psiprior=truncated(Exponential(1), 0, 1000),
    thetaprior=Beta(1, 1),
    sigma2prior=Exponential(1),
)
    α1 ~ alpha1prior
    α2 ~ alpha2prior
    α3 ~ alpha3prior
    α4 ~ alpha4prior
    α5 ~ alpha5prior
    α6 ~ alpha6prior
    α7 ~ alpha7prior
    α8 ~ alpha8prior
    ω ~ omegaprior 

    if psiprior isa Number 
        ψ = psiprior 
    else
        ψ ~ psiprior
    end

    θ ~ thetaprior
    sigma2 ~ sigma2prior

    T = typeof(α1)
   # u0 = zeros(16)
   # u0[1] = 1.0
    isolating = zeros(Float64, ndates)

    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community
    
    for j ∈ 1:nhospitals
        p = HCWSEIIRRRp(
            max(zero(T), α4 + α5 * vpd[j] + α6 * psb[j]),  # βh 
            max(zero(T), α1 + α2 * vpd[j] + α3 * psb[j]),  # βp 
            0.5,  # η 
            0.2,  # γ 
            ψ,
            ω,  
            θ
        )
        hcwseiirrr_isolating!(isolating, automatic, p, 1:ndates, λc, patients, vaccinated, j)

        for t ∈ 1:ndates
            Turing.@addlogprob! logpdf(Normal(isolating[t], sigma2), staff[t, j])
        end
    end
end

function loadchainsdf(filenamestart; psiprior=Exponential(1), kwargs...)
    return loadchainsdf(filenamestart, psiprior; kwargs...)
end

function loadchainsdf(filenamestart, ::Sampleable; ids=1:5, maxrounds=12,)
    df = DataFrame(
        :iteration => Int[ ],
        :chain => Int[ ],
        :α1 => Float64[ ],
        :α2 => Float64[ ],
        :α3 => Float64[ ],
        :α4 => Float64[ ],
        :α5 => Float64[ ],
        :α6 => Float64[ ],
        :α7 => Float64[ ],
        :α8 => Float64[ ],
        :ω => Float64[ ],
        :ψ => Float64[ ],
        :θ => Float64[ ],
        :sigma2 => Float64[ ],
        :log_density => Float64[ ],
    )
    for i ∈ ids
        _loaded = false 
        j = maxrounds
        while !_loaded && j >= 1
            filename = "$(filenamestart)_id_$(i)_round_$(j).jld2"
            if isfile(datadir("sims", filename))
                _df = DataFrame(load(datadir("sims", filename))["chain"])
                _df.chain = [ i for _ ∈ axes(_df, 1) ]
                df = vcat(df, _df)
                _loaded = true
            end
            j += -1
        end
    end
    return df
end

function loadchainsdf(filenamestart, psiprior::Number; ids=1:5, maxrounds=12,)
    df = DataFrame(
        :iteration => Int[ ],
        :chain => Int[ ],
        :α1 => Float64[ ],
        :α2 => Float64[ ],
        :α3 => Float64[ ],
        :α4 => Float64[ ],
        :α5 => Float64[ ],
        :α6 => Float64[ ],
        :α7 => Float64[ ],
        :α8 => Float64[ ],
        :ω => Float64[ ],
        :θ => Float64[ ],
        :sigma2 => Float64[ ],
        :log_density => Float64[ ],
    )
    for i ∈ ids
        _loaded = false 
        j = maxrounds
        while !_loaded && j >= 1
            filename = "$(filenamestart)_id_$(i)_round_$(j).jld2"
            if isfile(datadir("sims", filename))
                _df = DataFrame(load(datadir("sims", filename))["chain"])
                _df.chain = [ i for _ ∈ axes(_df, 1) ]
                df = vcat(df, _df)
                _loaded = true
            end
            j += -1
        end
    end
    insertcols!(df, 12, :ψ => psiprior)
    return df
end

function predictdiagnoses(
    df::DataFrame, psi, 
    patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
)
    diagnoses = zeros(Float64, ndates, nhospitals, size(df, 1))
    for i ∈ axes(df, 1)
        λc = [ 
            max(zero(Float64), df.α7[i] + df.α8[i] * (100 - s)) 
            for s ∈ stringency 
        ] .* community
        for j ∈ 1:nhospitals
            p = HCWSEIIRRRp(
                max(zero(Float64), df.α4[i] + df.α5[i] * vpd[j] + df.α6[i] * psb[j]),  # βh 
                max(zero(Float64), df.α1[i] + df.α2[i] * vpd[j] + df.α3[i] * psb[j]),  # βp 
                0.5,  # η 
                0.2,  # γ 
                _predictdiagnosespsi(psi, df, i),
                df.ω[i],  
                df.θ[i]
            )
            hcwseiirrr_isolating!(
                view(diagnoses, :, j, i), automatic, p, 1:ndates, 
                λc, patients, vaccinated, j
            )
        end
    end
    return diagnoses
end

function predictdiagnoses(
    df::DataFrame, patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
)
    return predictdiagnoses(
        df, automatic, 
        patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
end

_predictdiagnosespsi(psi::Number, ::Any, ::Any) = psi
_predictdiagnosespsi(::Automatic, df, i) = df.ψ[i]

function predicttotaldiagnoses(args...; cri=( 0.05, 0.95 ))
    predicteddiagnoses = predictdiagnoses(args...)
    totaldiagnoses = zeros(Float64, size(predicteddiagnoses, 2), size(predicteddiagnoses, 3))
    mediantotaldiagnoses = zeros(Float64, size(predicteddiagnoses, 2))
    lcitotaldiagnoses = zeros(Float64, size(predicteddiagnoses, 2))
    ucitotaldiagnoses = zeros(Float64, size(predicteddiagnoses, 2))
    for i ∈ axes(predicteddiagnoses, 2)
        for j ∈ axes(predicteddiagnoses, 3)
            totaldiagnoses[i, j] = sum(@view predicteddiagnoses[:, i, j]) / 10
        end
        mediantotaldiagnoses[i] = quantile(totaldiagnoses[i, :], 0.5)
        lcitotaldiagnoses[i] = quantile(totaldiagnoses[i, :], cri[1])
        ucitotaldiagnoses[i] = quantile(totaldiagnoses[i, :], cri[2])
    end
    return (
        totaldiagnoses=totaldiagnoses, 
        mediantotaldiagnoses=mediantotaldiagnoses, 
        lcitotaldiagnoses=lcitotaldiagnoses, 
        predicteddiagnoses=predicteddiagnoses, 
        ucitotaldiagnoses=ucitotaldiagnoses,
    )
end

function calculatebetah(df, vpd, psb, i, j) 
    return max(zero(Float64), df.α4[i] + df.α5[i] * vpd[j] + df.α6[i] * psb[j])
end

function calculatebetap(df, vpd, psb, i, j) 
    return max(zero(Float64), df.α1[i] + df.α2[i] * vpd[j] + df.α3[i] * psb[j])
end

function calculatelambdac(df, stringency, community, i)
    λc = [ 
        max(zero(Float64), df.α7[i] + df.α8[i] * (100 - s)) 
        for s ∈ stringency 
    ] .* community
    return λc
end

function calculatebetahs(df, vpd, psb, nhospitals; cri=( 0.05, 0.95 )) 
    medianbetah = zeros(Float64, nhospitals)
    lcibetah = zeros(Float64, nhospitals)
    ucibetah = zeros(Float64, nhospitals)

    for j ∈ 1:nhospitals
        betas = [ calculatebetah(df, vpd, psb, i, j) for i ∈ axes(df, 1) ]
        medianbetah[j] = quantile(betas, 0.5)
        lcibetah[j] = quantile(betas, cri[1])
        ucibetah[j] = quantile(betas, cri[2])
    end

    return @ntuple medianbetah lcibetah ucibetah
end

function calculatebetaps(df, vpd, psb, nhospitals; cri=( 0.05, 0.95 )) 
    medianbetap = zeros(Float64, nhospitals)
    lcibetap = zeros(Float64, nhospitals)
    ucibetap = zeros(Float64, nhospitals)

    for j ∈ 1:nhospitals
        betas = [ calculatebetap(df, vpd, psb, i, j) for i ∈ axes(df, 1) ]
        medianbetap[j] = quantile(betas, 0.5)
        lcibetap[j] = quantile(betas, cri[1])
        ucibetap[j] = quantile(betas, cri[2])
    end

    return @ntuple medianbetap lcibetap ucibetap
end

function calculatelambdacs(df, stringency, community; cri=( 0.05, 0.95 )) 
    medianlambdac = zeros(Float64, length(community))
    lcilambdac = zeros(Float64, length(community)) 
    ucilambdac = zeros(Float64, length(community))
    alllambddas = [ calculatelambdac(df, stringency, community, i) for i ∈ axes(df, 1) ]

    for t ∈ axes(community, 1)
        λ = [ alllambddas[i][t] for i ∈ axes(df, 1) ]
        medianlambdac[t] = quantile(λ, 0.5)
        lcilambdac[t] = quantile(λ, cri[1])
        ucilambdac[t] = quantile(λ, cri[2])
    end

    return @ntuple medianlambdac lcilambdac ucilambdac
end

function calculatebetas(df, vpd, psb, nhospitals; kwargs...)
    @unpack medianbetah, lcibetah, ucibetah = calculatebetahs(
        df, vpd, psb, nhospitals; 
        kwargs...
    ) 
    @unpack medianbetap, lcibetap, ucibetap = calculatebetaps(
        df, vpd, psb, nhospitals; 
        kwargs...
    ) 
    @ntuple medianbetah lcibetah ucibetah medianbetap lcibetap ucibetap 
end

function processoutputs(
    data::DataFrame, coviddata::DataFrame, 
    chaindf::DataFrame, chaindf_psi0::DataFrame, 
    vaccinated::Vector; 
    dateid=:Date
)
    # `data` can be the covid data or simulated data. `coviddata` must be the covid data
    nhospitals = counthospitals(data)
    ndates = countdates(data; dateid)
    @unpack newstaff, patients, staff = datamatrices(data, ndates, nhospitals)
    totalinfections = [ sum(@view newstaff[:, i]) for i ∈ 1:nhospitals ]
    @unpack vpd, psb = hospitalconditionmatrices(data)
    stringency = coviddata.StringencyIndex_Average[1:ndates]
    community = data.weeklycases[1:ndates] ./ 56_000_000
    boostpredictions = predicttotaldiagnoses(
        chaindf, patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
    noboostpredictions = predicttotaldiagnoses(
        chaindf_psi0, 0, 
        patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
    @unpack medianbetah, lcibetah, ucibetah, medianbetap, lcibetap, ucibetap = calculatebetas(
        chaindf, vpd, psb, nhospitals
    )
    @unpack medianlambdac, lcilambdac, ucilambdac = calculatelambdacs(
        chaindf, stringency, community
    ) 

    return Dict(
        :chaindf => chaindf,
        :chaindf_psi0 => chaindf_psi0,
        :community => community,
        :data => data,
        :lcibetah => lcibetah,
        :lcibetap => lcibetap,
        :lcilambdac => lcilambdac,
        :lcitotaldiagnoses => boostpredictions.lcitotaldiagnoses,
        :lcitotaldiagnosesnoboost => noboostpredictions.lcitotaldiagnoses,
        :medianbetah => medianbetah,
        :medianbetap => medianbetap,
        :medianlambdac => medianlambdac,
        :mediantotaldiagnoses => boostpredictions.mediantotaldiagnoses,
        :mediantotaldiagnosesnoboost => noboostpredictions.mediantotaldiagnoses,
        :ndates => ndates,
        :nhospitals => nhospitals,
        :patients => patients, 
        :predictdiagnoses => boostpredictions.predicteddiagnoses,
        :predictdiagnosesnoboost => noboostpredictions.predicteddiagnoses,
        :psb => psb,
        :staff => staff,
        :stringency => stringency,
        :totaldiagnoses => boostpredictions.totaldiagnoses,
        :totaldiagnosesnoboost => noboostpredictions.totaldiagnoses,
        :totalinfections => totalinfections,
        :ucibetah => ucibetah,
        :ucibetap => ucibetap,
        :ucilambdac => ucilambdac,
        :ucitotaldiagnoses => boostpredictions.ucitotaldiagnoses,
        :ucitotaldiagnosesnoboost => noboostpredictions.ucitotaldiagnoses,
        :vpd => vpd,
    )
end

function processoutputs(
    data::DataFrame, chaindf::DataFrame, chaindf_psi0::DataFrame, vaccinated::Vector; 
    kwargs...
)
    return processoutputs(data, data, chaindf, chaindf_psi0, vaccinated; kwargs...)
end

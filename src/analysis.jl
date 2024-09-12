
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
    return @ntuple vpd=_tdf.VolumePerBed psb=_tdf.ProportionSingleBeds 
end

function predictinfections(
    βp::T, βh::T, βc::T, ψ::T, ω::T, patients, staff, community, vaccinated; 
    immunevectorlength=10,
) where T 
    ndates = length(βc)
    nhospitals = length(βp)
    @assert nhospitals == length(βh)

    predictedinfections = zeros(T, ndates, nhospitals)
    immunevector = SizedVector{immunevectorlength}(zeros(T, immunevectorlength))  

    predictinfections!(
        predictedinfections, immunevector, ndates, nhospitals, 
        βp, βh, βc, ψ, ω, patients, staff, community, vaccinated
    )

    return predictedinfections
end

function predictinfections(
    df::DataFrame, patients, staff, community, vaccinated, vpd, psb, stringency; 
    immunevectorlength=10, inds=automatic, ψ=automatic,
)
    return _predictinfections(
        df, patients, staff, community, vaccinated, vpd, psb, stringency, ψ, inds; 
        immunevectorlength
    )
end

function _predictinfections(
    df, patients, staff, community, vaccinated, vpd, psb, stringency, ψ, ::Automatic; 
    immunevectorlength
)
    allinds = axes(df, 1)
    return _predictinfections(
        df, patients, staff, community, vaccinated, vpd, psb, stringency, ψ, allinds; 
        immunevectorlength
    )
end

function _predictinfections(
    df, patients, staff, community, vaccinated, vpd, psb, stringency, ψ, ind::Number; 
    immunevectorlength
)
    allinds = [ ind ]
    return _predictinfections(
        df, patients, staff, community, vaccinated, vpd, psb, stringency, ψ, allinds; 
        immunevectorlength
    )
end

function _predictinfections(
    df, patients, staff, community, vaccinated, vpd, psb, stringency, 
    ::Automatic, inds::AbstractVector; 
    immunevectorlength
)
    ndates, nhospitals = size(patients)
    nsamples = length(inds)
    predictedinfections = zeros(ndates, nhospitals, nsamples)
    immunevector = SizedVector{immunevectorlength}(zeros(immunevectorlength))

    @unpack βp, βh, βc = calculatebetas(df, vpd, psb, stringency)

    for (i, ind) ∈ enumerate(inds)
        predictinfections!(
            @view(predictedinfections[:, :, i]), immunevector, 
            ndates, nhospitals, 
            βp[ind], βh[ind], βc[ind], getproperty(df, :ψ)[ind], getproperty(df, :ω)[ind], 
            patients, staff, community, vaccinated
        )
    end

    return predictedinfections
end

function _predictinfections(
    df, patients, staff, community, vaccinated, vpd, psb, stringency, 
    ψ::Number, inds::AbstractVector; 
    immunevectorlength
)
    ndates, nhospitals = size(patients)
    nsamples = length(inds)
    predictedinfections = zeros(ndates, nhospitals, nsamples)
    immunevector = SizedVector{immunevectorlength}(zeros(immunevectorlength))

    @unpack βp, βh, βc = calculatebetas(df, vpd, psb, stringency)

    for (i, ind) ∈ enumerate(inds)
        predictinfections!(
            @view(predictedinfections[:, :, i]), immunevector, 
            ndates, nhospitals, 
            βp[ind], βh[ind], βc[ind], ψ, getproperty(df, :ω)[ind], 
            patients, staff, community, vaccinated
        )
    end

    return predictedinfections
end

function predictinfections!(
    predictedinfections::AbstractMatrix{T}, immunevector::SizedVector{N, T}, 
    βp, βh, βc, ψ, ω, patients, staff, community, vaccinated
) where {N, T}
    ndates = length(βc)
    nhospitals = length(βp)
    @assert nhospitals == length(βh)

    predictinfections!(
        predictedinfections, immunevector, ndates, nhospitals, 
        βp, βh, βc, ψ, ω, patients, staff, community, vaccinated
    )
end

function predictinfections!(
    predictedinfections::AbstractMatrix{T}, immunevector::SizedVector{N, T}, 
    ndates, nhospitals, βp, βh, βc, ψ, ω, patients, staff, community, vaccinated
) where {N, T}
    foi = zeros(T, ndates, nhospitals)
    for t ∈ axes(foi, 1), j ∈ axes(foi, 2)
        foi[t, j] = βp[j] * patients[t, j] + βh[j]* staff[t, j] + βc[t] * community[t]
    end

    predictinfections!(
        predictedinfections, immunevector, ndates, nhospitals, foi, ψ, ω, vaccinated
    )
end

function predictinfections!(
    predictedinfections::AbstractMatrix{T}, immunevector::SizedVector{N, T}, 
    ndates, nhospitals, foi, ψ, ω, vaccinated
) where {N, T}
    for j ∈ 1:nhospitals
        @assert predictedinfections[1, j] == 0
        for x ∈ 1:N  # reset immunevector
            immunevector[x] = zero(T)
        end
        for t ∈ 2:ndates  
            v = vaccinated[(t - 1), j]
            immune10 = predictedinfections[(t - 1), j] + 
                v * (1 - sum(immunevector) - predictedinfections[(t - 1), j])
            # probability of boosting from natural immune boosting plus vaccination
            pb = (1 - exp(-ψ * foi[(t - 1), j])) * (1 - v) + v  
            for x ∈ 1:N-1
                immune10 += pb * immunevector[x]
                immunevector[x] += -(pb + N * ω) * immunevector[x] + 
                    N * ω * (1 - pb) * immunevector[x+1]
            end
            immunevector[N] += -(N * ω * (1 - pb)) * immunevector[N] + immune10
            predictedinfections[t, j] = (1 - sum(immunevector)) * (1 - exp(-foi[(t - 1), j]))
        end
    end
end

#= (version 1)
@model function fitmodel(
    newstaff, patients, staff, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals;
    alpha1prior=truncated(Normal(0, 1), -1, 10),
    alpha2prior=Normal(0, 0.001),
    alpha3prior=Normal(0, 1),
    alpha4prior=truncated(Normal(0, 1), -1, 10),
    alpha5prior=Normal(0, 0.001),
    alpha6prior=Normal(0, 1),
    alpha7prior=truncated(Normal(0, 1), -1, 10),
    alpha8prior=Normal(0, 0.1),
    omegaprior=Uniform(0, 0.1),
    psiprior=Exponential(1),
    sigma2prior=Exponential(1),
    immunevectorlength=10,
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
    ψ ~ psiprior
    sigma2 ~ sigma2prior

    T = typeof(α1)

    # levels of immunity
    immunevector = SizedVector{immunevectorlength}(zeros(T, immunevectorlength))  

    βc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community

    for j ∈ 1:nhospitals
        predictedinfections = 0
        for x ∈ 1:immunevectorlength  # reset immunevector
            immunevector[x] = zero(T)
        end
        βp = max(zero(T), α1 + α2 * vpd[j] + α3 * psb[j]) .* patients[:, j]
        βh = max(zero(T), α4 + α5 * vpd[j] + α6 * psb[j]) .* staff[:, j]
        foi = βp .+ βh .+ βc
        for t ∈ 2:ndates 
            v = vaccinated[(t - 1), j]
            immune10 = predictedinfections + v * (1 - sum(immunevector) - predictedinfections)
            # probability of boosting from natural immune boosting plus vaccination
            pb = (1 - exp(-ψ * foi[t])) * (1 - v) + v  
            for x ∈ 1:immunevectorlength-1
                immune10 += pb * immunevector[x]
                immunevector[x] += -(pb + immunevectorlength * ω) * immunevector[x] + 
                    immunevectorlength * ω * (1 - pb) * immunevector[x+1]
            end
            immunevector[immunevectorlength] += -(immunevectorlength * ω * (1 - pb)) * 
                immunevector[immunevectorlength] + 
                immune10
            predictedinfections = (1 - sum(immunevector)) * (1 - exp(-foi[t]))
            Turing.@addlogprob! logpdf(Normal(newstaff[t, j], sigma2), predictedinfections)
        end
    end
end
=#


@model function fitmodel(
    newstaff, patients, staff, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals;
    alpha1prior=truncated(Normal(0, 1), -1, 10),
    alpha2prior=Normal(0, 0.001),
    alpha3prior=Normal(0, 1),
    alpha4prior=truncated(Normal(0, 1), -1, 10),
    alpha5prior=Normal(0, 0.001),
    alpha6prior=Normal(0, 1),
    alpha7prior=truncated(Normal(0, 1), -1, 10),
    alpha8prior=truncated(Normal(0, 0.1), 0, 10),  # require greater stringency leads to less transmission
    omegaprior=Uniform(0, 0.1),
    psiprior=Exponential(1),
    sigma2prior=Exponential(1),
    immunevectorlength=3,
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

    sigma2 ~ sigma2prior

    T = typeof(α1)

    # levels of immunity
    immunevector = SizedVector{immunevectorlength}(zeros(T, immunevectorlength))  

    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community
    
    for j ∈ 1:nhospitals
        predictedinfections = 0
        for x ∈ 1:immunevectorlength  # reset immunevector
            immunevector[x] = zero(T)
        end
        λp = max(zero(T), α1 + α2 * vpd[j] + α3 * psb[j]) .* patients[:, j]
        βh = max(zero(T), α4 + α5 * vpd[j] + α6 * psb[j]) 
        for t ∈ 2:ndates 
            foi = λc[t] + λp[j] + βh * predictedinfections
            v = vaccinated[(t - 1), j]
            immune10 = predictedinfections + v * (1 - sum(immunevector) - predictedinfections)
            # probability of boosting from natural immune boosting plus vaccination
            pb = (1 - exp(-ψ * foi)) * (1 - v) + v  
            for x ∈ 1:immunevectorlength-1
                immune10 += pb * immunevector[x]
                immunevector[x] += -(pb + immunevectorlength * ω) * immunevector[x] + 
                    immunevectorlength * ω * (1 - pb) * immunevector[x+1]
            end
            immunevector[immunevectorlength] += -(immunevectorlength * ω * (1 - pb)) * 
                immunevector[immunevectorlength] + 
                immune10
            predictedinfections = (1 - sum(immunevector)) * (1 - exp(-foi))
            Turing.@addlogprob! logpdf(Normal(newstaff[t, j], sigma2), predictedinfections)
        end
    end
end

#=
@model function fitmodel(
    newstaff, patients, staff, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals;
    alpha1prior=truncated(Normal(0, 1), -1, 10),
    alpha2prior=Normal(0, 0.001),
    alpha3prior=Normal(0, 1),
    alpha4prior=truncated(Normal(0, 1), -1, 10),
    alpha5prior=Normal(0, 0.001),
    alpha6prior=Normal(0, 1),
    alpha7prior=truncated(Normal(0, 1), -1, 10),
    alpha8prior=truncated(Normal(0, 0.1), 0, 10),  # require greater stringency leads to less transmission
    omegaprior=Uniform(0, 0.1),
    psiprior=Exponential(1),
    sigma2prior=Exponential(1),
    immunevectorlength=3,
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

    sigma2_bp ~ sigma2prior
    sigma2_bh ~ sigma2prior
    sigma2 ~ sigma2prior

    T = typeof(α1)

    # levels of immunity
    immunevector = SizedVector{immunevectorlength}(zeros(T, immunevectorlength))  

    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community
    βp = zeros(T, nhospitals)
    βh = zeros(T, nhospitals)
    
    for j ∈ 1:nhospitals
        predictedinfections = 0
        for x ∈ 1:immunevectorlength  # reset immunevector
            immunevector[x] = zero(T)
        end
        βp[j] ~ truncated(Normal(α1 + α2 * vpd[j] + α3 * psb[j], sigma2_bp), 0, 10)
        βh[j] ~ truncated(Normal(α4 + α5 * vpd[j] + α6 * psb[j], sigma2_bh), 0, 10)
        for t ∈ 2:ndates 
            foi = λc[t] + βp[j] * patients[t, j] + βh[j] * predictedinfections
            v = vaccinated[(t - 1), j]
            immune10 = predictedinfections + v * (1 - sum(immunevector) - predictedinfections)
            # probability of boosting from natural immune boosting plus vaccination
            pb = (1 - exp(-ψ * foi)) * (1 - v) + v  
            for x ∈ 1:immunevectorlength-1
                immune10 += pb * immunevector[x]
                immunevector[x] += -(pb + immunevectorlength * ω) * immunevector[x] + 
                    immunevectorlength * ω * (1 - pb) * immunevector[x+1]
            end
            immunevector[immunevectorlength] += -(immunevectorlength * ω * (1 - pb)) * 
                immunevector[immunevectorlength] + 
                immune10
            predictedinfections = (1 - sum(immunevector)) * (1 - exp(-foi))
            Turing.@addlogprob! logpdf(Normal(newstaff[t, j], sigma2), predictedinfections)
        end
    end
end
=#

#=
@memoize function _memofitfunction!(
    immunevector::SizedVector{immunevectorlength, T, Vector{T}}, 
    α1::T, α2, α3, α4, α5, α6, α7, α8, ω, ψ, sigma2, 
    newstaff, patients, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals
) where {immunevectorlength, T}
    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community

    for j ∈ 1:nhospitals
        predictedinfections = 0
        for x ∈ 1:immunevectorlength  # reset immunevector
            immunevector[x] = zero(T)
        end
        λp = max(zero(T), α1 + α2 * vpd[j] + α3 * psb[j]) .* patients[:, j]
        βh = max(zero(T), α4 + α5 * vpd[j] + α6 * psb[j]) 
        for t ∈ 2:ndates 
            foi = λc[t] + λp[j] + βh * predictedinfections
            v = vaccinated[(t - 1), j]
            immune10 = predictedinfections + v * (1 - sum(immunevector) - predictedinfections)
            # probability of boosting from natural immune boosting plus vaccination
            pb = (1 - exp(-ψ * foi)) * (1 - v) + v  
            for x ∈ 1:immunevectorlength-1
                immune10 += pb * immunevector[x]
                immunevector[x] += -(pb + immunevectorlength * ω) * immunevector[x] + 
                    immunevectorlength * ω * (1 - pb) * immunevector[x+1]
            end
            immunevector[immunevectorlength] += -(immunevectorlength * ω * (1 - pb)) * 
                immunevector[immunevectorlength] + 
                immune10
            predictedinfections = (1 - sum(immunevector)) * (1 - exp(-foi))
            Turing.@addlogprob! logpdf(Normal(newstaff[t, j], sigma2), predictedinfections)
        end
    end
end
=#
function loadchainsdf(filenamestart; psiprior=Exponential(1), kwargs...)
    return _loadchainsdf(filenamestart, psiprior; kwargs...)
end

function _loadchainsdf(filenamestart, ::Sampleable; ids=1:5, maxrounds=12,)
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

function _loadchainsdf(filenamestart, psiprior::Number; ids=1:5, maxrounds=12,)
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

function _calculatebetap(α1::T, α2::T, α3::T, vpd, psb) where T
    return [ max(zero(T), α1 + α2 * v + α3 * p) for (v, p) ∈ zip(vpd, psb) ]
end

function _calculatebetap(df::DataFrame, vpd, psb) 
    return [
        _calculatebetap(
            getproperty(df, :α1)[i], getproperty(df, :α2)[i], getproperty(df, :α3)[i], 
            vpd, psb
        )
        for i ∈ axes(df, 1)
    ]
end

function _calculatebetah(α4::T, α5::T, α6::T, vpd, psb) where T
    return [ max(zero(T), α4 + α5 * v + α6 * p) for (v, p) ∈ zip(vpd, psb) ]
end

function _calculatebetah(df::DataFrame, vpd, psb) 
    return [
        _calculatebetah(
            getproperty(df, :α4)[i], getproperty(df, :α5)[i], getproperty(df, :α6)[i], 
            vpd, psb
        )
        for i ∈ axes(df, 1)
    ]
end

function _calculatebetac(α7::T, α8::T, stringency) where T
    return [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ]
end

function _calculatebetac(df::DataFrame, stringency) 
    return [
        _calculatebetac(getproperty(df, :α7)[i], getproperty(df, :α8)[i], stringency)
        for i ∈ axes(df, 1)
    ]
end

function calculatebetas(df, vpd, psb, stringency)
    βp = _calculatebetap(df, vpd, psb)
    βh = _calculatebetah(df, vpd, psb)
    βc = _calculatebetac(df, stringency)
    return @ntuple βp βh βc
end

function summarizepredictedinfections(
    predictedinfections::Array{<:Real, 3}; 
    cri=[ 0.05, 0.95 ],
)
    return _summarizepredictedinfections(predictedinfections, cri)
end

function summarizepredictedinfections(
    df::DataFrame, args...; 
    cri=[ 0.05, 0.5, 0.95 ], kwargs...
)
    predictedinfections = predictinfections(df, args...; kwargs...) 
    return _summarizepredictedinfections(predictedinfections, cri)
end

function _summarizepredictedinfections(predictedinfections, cri::Number)
    @assert cri <= 1 
    lcri = (1 - cri) / 2
    ucri = 1 - lcri 
    return _summarizepredictedinfections(predictedinfections, [ lcr, 0.5, ucri ])
end

function _summarizepredictedinfections(predictedinfections, cri)
    ndates, nhospitals, nsamples = size(predictedinfections)      
    totals = zeros(nsamples, nhospitals)
    means = zeros(nhospitals)
    lcris = zeros(nhospitals)
    medians = zeros(nhospitals)
    ucris = zeros(nhospitals)

    for i ∈ axes(totals, 2), j ∈ axes(totals, 1)
        totals[j, i] = sum(@view predictedinfections[:, i, j])
    end

    for i ∈ 1:nhospitals
        means[i] = mean(totals[:, i])
        lcri, median, ucri = quantile(totals[:, i], cri)
        lcris[i] = lcri
        medians[i] = median
        ucris[i] = ucri
    end

    return @ntuple totals means lcris medians ucris
end

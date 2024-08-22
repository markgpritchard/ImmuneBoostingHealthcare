
counthospitals(data) = length(unique(data.Code))
countdates(data) = length(unique(data.Date)) 

function datamatrices(data, ndates, nhospitals)
    patients = zeros(ndates, nhospitals)
    staff = zeros(ndates, nhospitals) 
    newstaff = zeros(ndates, nhospitals) 
    for (i, c) ∈ enumerate(unique(data.Code))
        _tdf = filter(:Code => x -> x == c, data)
        for t ∈ 1:ndates
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
        end
    end
    return @ntuple newstaff patients staff
end

function hospitalconditionmatrices(data)
    _tdf = select(data, :Code, :VolumePerBed, :ProportionSingleBeds)
    unique!(_tdf)
    return @ntuple vpd=_tdf.VolumePerBed psb=_tdf.ProportionSingleBeds 
end

@model function fitmodel(
    newstaff, patients, staff, vaccinated, weeklycases, 
    vpd, psb, stringency, ndates, nhospitals;
    alpha1prior=Beta(1, 1),
    alpha2prior=Exponential(1),
    alpha3prior=Exponential(1),
    alpha4prior=Beta(1, 1),
    alpha5prior=Exponential(1),
    alpha6prior=Exponential(1),
    alpha7prior=Beta(1, 1),
    alpha8prior=Exponential(1),
    aprior=Exponential(0.3),
    thetaprior=Exponential(500),
    psiprior=Exponential(1),
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
    a ~ aprior
    θ ~ thetaprior 
    ψ ~ psiprior
    sigma2 ~ sigma2prior

    T = typeof(α1)

    rjt = zeros(T, ndates)  # levels of immunity

    βp = [ α1 + α2 * v + α3 * p for(v, p) ∈ zip(vpd, psb) ]
    βh = [ α4 + α5 * v + α6 * p for(v, p) ∈ zip(vpd, psb) ]
    βc = @. α7 + α8 * (100 - stringency)

    foi = zeros(T, ndates, nhospitals)
    for d ∈ axes(foi, 1), h ∈ axes(foi, 2)
        foi[d, h] = βp[h] * patients[d, h] + βh[h]* staff[d, h] + βc[d] * weeklycases[d] / 56_000_000
    end

    predictedinfections = zeros(T, ndates, nhospitals)
    for h ∈ 1:nhospitals
        for d ∈ 1:ndates  # reset the vector rjt 
            rjt[d] = zero(T)
        end
        for d ∈ 2:ndates  # reset the vector rjt 
            rjt1 = predictedinfections[(d - 1), h] + (1 - exp(-ψ * foi[(d - 1), h])) * sum(rjt)
            for x ∈ d:-1:2
                rjt[x] = rjt[x-1] * 
                    exp(-ψ * foi[(d - 1), h]) * 
                    (1 - pdf(Weibull(a, θ), x) / (ccdf(Weibull(a, θ), x - 1)))
            end
            rjt[1] = rjt1 
            predictedinfections[d, h] = (1 - sum(rjt)) * (1 - exp(-foi[(d - 1), h]))
        end
    end

    for d ∈ axes(predictedinfections, 1), h ∈ axes(predictedinfections, 2)
        d <= 14 && continue 
        Turing.@addlogprob! logpdf(Normal(newstaff[d, h], sigma2), predictedinfections[d, h])
    end
end

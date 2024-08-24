
counthospitals(data) = length(unique(data.Code))
countdates(data; dateid=:Date) = length(unique(getproperty(data, dateid))) 

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
    end
    return @ntuple newstaff patients staff
end

function hospitalconditionmatrices(data)
    _tdf = select(data, :Code, :VolumePerBed, :ProportionSingleBeds)
    unique!(_tdf)
    return @ntuple vpd=_tdf.VolumePerBed psb=_tdf.ProportionSingleBeds 
end

function _predictinfections!(
    predictedinfections::AbstractMatrix{T}, immunevector::SizedVector{N, T}, 
    ndates, nhospitals, foi, vaccinated, ψ, ω
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
                immunevector[x] += -(pb + N * ω) * immunevector[x] + N * ω * (1 - pb) * immunevector[x+1]
            end
            immunevector[N] += -(N * ω * (1 - pb)) * immunevector[N] + immune10
            predictedinfections[t, j] = (1 - sum(immunevector)) * (1 - exp(-foi[(t - 1), j]))
        end
    end
end

@model function fitmodel(
    newstaff, patients, staff, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals;
    alpha1prior=Beta(1, 1),
    alpha2prior=Exponential(0.1),
    alpha3prior=Exponential(0.1),
    alpha4prior=Beta(1, 1),
    alpha5prior=Exponential(0.1),
    alpha6prior=Exponential(0.1),
    alpha7prior=Beta(1, 1),
    alpha8prior=Exponential(0.1),
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

    βp = [ α1 + α2 * v + α3 * p for(v, p) ∈ zip(vpd, psb) ]
    βh = [ α4 + α5 * v + α6 * p for(v, p) ∈ zip(vpd, psb) ]
    βc = @. α7 + α8 * (100 - stringency)

    foi = zeros(T, ndates, nhospitals)
    for t ∈ axes(foi, 1), j ∈ axes(foi, 2)
        foi[t, j] = βp[j] * patients[t, j] + βh[j]* staff[t, j] + βc[t] * community[t]
    end

    predictedinfections = zeros(T, ndates, nhospitals)
    _predictinfections!(
        predictedinfections, immunevector, ndates, nhospitals, foi, vaccinated, ψ, ω
    )

    for t ∈ axes(predictedinfections, 1), j ∈ axes(predictedinfections, 2)
        t <= 14 && continue 
        Turing.@addlogprob! logpdf(Normal(newstaff[t, j], sigma2), predictedinfections[t, j])
    end
end

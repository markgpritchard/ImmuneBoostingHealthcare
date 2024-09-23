
_beta(p::SEIIRRRSp{<:Number, <:Any, <:Any}, ::Any) = p.β
_beta(p::SEIIRRRSp{<:Function, <:Any, <:Any}, t::Number) = p.β(t)
_betahh(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.βhh
_betahp(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.βhp
_betaph(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.βph
_betapp(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.βpp
betahh(args...) = _betahh(args...)
betahp(args...) = _betahp(args...)
betaph(args...) = _betaph(args...)
betapp(args...) = _betapp(args...)
_eta(p::AbstractParameters{<:Any, <:Number, <:Any}, ::Any) = p.η
_gamma(p::AbstractParameters{<:Any, <:Number, <:Any}, ::Any) = p.γ
_psi(p::AbstractParameters{<:Any, <:Number, <:Any}, ::Any) = p.ψ
_nu(p::AbstractParameters{<:Any, <:Any, <:Number}, ::Any) = p.nu
_nu(p::AbstractParameters{<:Any, <:Any, <:Function}, t::Number) = p.nu(t)
_omega(p::AbstractParameters{<:Any, <:Number, <:Any}, ::Any) = p.ω
_lambdac(p::WXYYZSEIIRRRSp{<:Number, <:Any, <:Any}, ::Any) = p.λc
_lambdac(p::WXYYZSEIIRRRSp{<:Function, <:Any, <:Any}, t::Number) = p.λc(t)
_deltan(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.δn
_deltai(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.δi
_theta(p::WXYYZSEIIRRRSp{<:Any, <:Number, <:Any}, ::Any) = p.θ

function _stochasticseiirrrs!(u, t, p)
    S, E, I1, I2, R1, R2, R3 = u
    pI = (I1 + I2) / sum(u)

    rates = [
        _beta(p, t) * S * pI,  # infection 
        _eta(p, t) * E,  # transition E -> I1
        _gamma(p, t) * I1 * 2,  # transition I1 -> I2
        _gamma(p, t) * I2 * 2,  # transition I2 -> R1
        _omega(p, t) * R1 * 3,  # transition R1 -> R2
        _omega(p, t) * R2 * 3,  # transition R2 -> R3
        _omega(p, t) * R3 * 3,  # transition R3 -> S
        (_beta(p, t) * _psi(p, t) * pI + _nu(p, t)) * R2,  # boosting plus vaccination from R2 
        (_beta(p, t) * _psi(p, t) * pI + _nu(p, t)) * R3,  # boosting plus vaccination from R3 
        _nu(p, t) * S,  # vaccination from S
    ]

    # what happens?
    eventid = sample(1:10, Weights(rates))
    u += SEIIRRRSTRANSITIONMATRIX[eventid, :]

    # when does it happen?
    rnum = rand() 
    timestep = -log(rnum) / sum(rates)
    t += timestep
    return ( u, t )
end

function stochasticseiirrrs(u0, trange, p)
    communityvalues = zeros(Int, length(trange), 8)
    u = deepcopy(u0)
    communityvalues[1, :] = u
    t = trange[1] 
    tind = 1
    while t < maximum(trange)
        u, t = _stochasticseiirrrs!(u, t, p)
        if t >= trange[tind+1]
            communityvalues[tind+1, :] = u
            tind += 1
        end
    end
    return communityvalues 
end

function _stochasticwxyyzseiirrrs!(u, t, p, communityvalues)
    W, X, Y1, Y2, Zr, Za, S, E, I1, I2, I′, R1, R2, R3, = u
    pY = (Y1 + Y2) / sum(@view u[1:6])
    pI = (I1 + I2) / sum(@view u[7:13])

    round(Int, t, RoundUp)
    # probability of admission in each group 
    denom = (
        sum(@view communityvalues[round(Int, t, RoundUp), 1:7]) + 
        sum(@view communityvalues[round(Int, t, RoundUp), 3:4])
    )
    probads = [
        communityvalues[round(Int, t, RoundUp), 1] / denom,
        communityvalues[round(Int, t, RoundUp), 2] / denom,
        communityvalues[round(Int, t, RoundUp), 3] * 2 / denom,
        communityvalues[round(Int, t, RoundUp), 4] * 2 / denom,
        sum(@view communityvalues[round(Int, t, RoundUp), 5:7]) / denom,
    ]

    rates = [
        _deltan(p, t) * W * probads[2],  # discharge from W and admission to X 
        _deltan(p, t) * W * probads[3],  # discharge from W and admission to Y1 
        _deltan(p, t) * W * probads[4],  # discharge from W and admission to Y2 
        _deltan(p, t) * W * probads[5],  # discharge from W and admission to Za 
        _deltan(p, t) * X * probads[1],  # discharge from X and admission to W 
        (_deltan(p, t) * probads[3] + _eta(p, t)) * X,  # discharge from X and admission to Y1 plus transition X -> Y1 
        _deltan(p, t) * X * probads[4],  # discharge from X and admission to Y2 
        _deltan(p, t) * X * probads[5],  # discharge from X and admission to Za 
        _deltai(p, t) * Y1 * probads[1],  # discharge from Y1 and admission to W 
        _deltai(p, t) * Y1 * probads[2],  # discharge from Y1 and admission to X 
        (_deltai(p, t) * probads[4] + _gamma(p, t) * 2) * Y1,  # discharge from Y1 and admission to Y2 plus transition Y1 -> Y2 
        _deltai(p, t) * Y1 * probads[5],  # discharge from Y1 and admission to Za 
        _deltai(p, t) * Y2 * probads[1],  # discharge from Y2 and admission to W 
        _deltai(p, t) * Y2 * probads[2],  # discharge from Y2 and admission to X 
        _deltai(p, t) * Y2 * probads[3],  # discharge from Y2 and admission to Y1 
        _deltai(p, t) * Y2 * probads[5],  # discharge from Y2 and admission to Za 
        _deltan(p, t) * Zr * probads[1],  # discharge from Zr and admission to W 
        _deltan(p, t) * Zr * probads[2],  # discharge from Zr and admission to X 
        _deltan(p, t) * Zr * probads[3],  # discharge from Zr and admission to Y1 
        _deltan(p, t) * Zr * probads[4],  # discharge from Zr and admission to Y2 
        _deltan(p, t) * Zr * probads[5],  # discharge from Zr and admission to Za 
        _deltan(p, t) * Za * probads[1],  # discharge from Za and admission to W 
        _deltan(p, t) * Za * probads[2],  # discharge from Za and admission to X 
        _deltan(p, t) * Za * probads[3],  # discharge from Za and admission to Y1 
        _deltan(p, t) * Za * probads[4],  # discharge from Za and admission to Y2 
        W * (_betaph(p, t) * pI + _betapp(p, t) * pY),  # patient infection 
        _gamma(p, t) * Y2 * 2,  # transition Y2 -> Zr
        S * (_betahh(p, t) * pI + _betahp(p, t) * pY + _lambdac(p, t)),  # healthcare worker infection 
        _eta(p, t) * E,  # transition E -> I1
        _gamma(p, t) * I1 * 2,  # transition I1 -> I2
        _theta(p, t) * I1,  # diagnosis from I1
        _theta(p, t) * I2,  # diagnosis from I2
        _gamma(p, t) * I2 * 2,  # transition I2 -> R1
        _gamma(p, t) * I′,  # transition I′ -> R1
        _omega(p, t) * R1 * 3,  # transition R1 -> R2
        _omega(p, t) * R2 * 3,  # transition R2 -> R3
        _omega(p, t) * R3 * 3,  # transition R3 -> S
        ((_betaph(p, t) * pI + _betapp(p, t) * pY) * _psi(p, t) + _nu(p, t)) * R2,  # boosting from R2 plus vaccination from R2 
        ((_betaph(p, t) * pI + _betapp(p, t) * pY) * _psi(p, t) + _nu(p, t)) * R3,  # boosting from R3 plus vaccination from R3 
        _nu(p, t) * S,  # vaccination from S
    ]

    # what happens?
    eventid = sample(1:40, Weights(rates))
    u += WXYYZSEIIRRRSTRANSITIONMATRIX[eventid, :]

    # when does it happen?
    rnum = rand() 
    timestep = -log(rnum) / sum(rates)
    t += timestep
    return ( u, t )
end

function stochasticwxyyzseiirrrs(u0, trange, p, communityvalues)
    hospitalvalues = zeros(Int, length(trange), 17)
    u = deepcopy(u0)
    hospitalvalues[1, :] = u
    t = trange[1] 
    tind = 1
    while t <= maximum(trange)
        u, t = _stochasticwxyyzseiirrrs!(u, t, p, communityvalues)
        if tind + 2 <= length(trange) && t >= trange[tind+2]
            hospitalvalues[tind+1, :] = hospitalvalues[tind, :]
            u = deepcopy(hospitalvalues[tind+1, :])
            tind += 1 
            t = trange[tind]
        elseif t >= trange[tind+1]
            hospitalvalues[tind+1, :] = u
            tind += 1
        end
    end
    hospitaldiagnoses = [  # NB hospitalvalues[:, 17] is cumulative diagnoses
        t <= 10 ? 
            hospitalvalues[t, 17] : 
            hospitalvalues[t, 17] - hospitalvalues[t-10, 17]
        for t ∈ round.(Int, trange)
    ]
    return @ntuple hospitaldiagnoses hospitalvalues 
end

function vaccinatestaff(t::Number)
    if 264 <= t <= 280  # 8 to 24 December 2020
        return 0.035
    elseif 285 <= t < 305  # 29 December 2020 to 18 January 2021
        return 0.055
    elseif 305 <= t <= 323  # 18 January to 5 February 2021
        return 0.025
        # values above give 89% vaccinated by 5 February 2021
    elseif 531 <= t <= 621  # 1 September to 30 November 2021
        return 0.0125  # 68% boosted over 90 days
    else 
        return 0.0
    end
end

vaccinatestaff(date::Date) = vaccinatestaff(datetot(date))

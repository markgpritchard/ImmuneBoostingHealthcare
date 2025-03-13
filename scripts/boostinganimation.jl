
using DrWatson
@quickactivate "ImmuneBoostingHealthcare"

using GLMakie, Distributions 

videosdir(args...) = projectdir("videos", args...)

#fixedimmunity(t_vaccinate) = [ t_vaccinate <= t <= t_vaccinate + 100 ? 1 : 0 for t ∈ 1:1:500 ]
fixedimmunity(t, t_vaccinate) = t_vaccinate <= t <= t_vaccinate + 100 ? 1 : 0
fixedimmunity(::Any, t, t_vaccinate) = fixedimmunity(t, t_vaccinate)

function boostedimmunity(epidemic, t, t_vaccinate)
    if t < t_vaccinate 
        return 0.0 
    elseif t <= t_vaccinate + 100 
        return 1.0 
    end
    vec = zeros(100); vec[1] = 1.0 
    for τ ∈ t_vaccinate:1:t 
        _boostedimmunity!(epidemic, vec, τ) 
    end
    return sum(vec) * maximum([ epidemic(x) for x ∈ t:500 ])
end

function _boostedimmunity!(epidemic, vec, τ)
    ep = epidemic(τ)
    _advanceboostedimmunity!(vec)
    _boostboostedimmunity!(vec, ep)
end

function _advanceboostedimmunity!(vec)
    for i ∈ 100:-1:2 
        vec[i] = vec[i-1]
    end    
    vec[1] = 0.0
end

function _boostboostedimmunity!(vec, ep)
    for i ∈ 2:1:100 
        v = ep * vec[i] 
        vec[1] += v 
        vec[i] += -v 
    end    
end

function minepiimmunity(epidemic, immunefunction, t, t_vaccinate; multiplier=1)
    epi = epidemic(t)
    rawimm = immunefunction(epidemic, t, t_vaccinate) 
    if rawimm == 1 
        imm = multiplier
    else 
        imm = rawimm 
    end
    return min(epi, imm)
end

function auc(epidemic, immunefunction, t_vaccinate; multiplier=1, trange=1:500)
    return sum(
        [ 
            minepiimmunity(epidemic, immunefunction, t, t_vaccinate; multiplier) 
            for t ∈ trange 
        ]
    )
end

normalepidemic(t) = 100 < t <= 400 ? 5 * pdf(Normal(250, 70), t) - 0.00286 : 0


# no immune boosting 

let 
    aucdata = [ auc(normalepidemic, fixedimmunity, x; multiplier=0.02565) for x ∈ 1:410 ]

    ep = [ normalepidemic(t) for t ∈ 1:410 ]
    data = [ 
        (
            v = fixedimmunity(normalepidemic, t, 1);
            v == 1 ? 0.02565 : v
        )
        for t ∈ 1:410 
    ]
    data2 = [ minepiimmunity(normalepidemic, fixedimmunity, t, 1; multiplier=0.02565) for t ∈ 1:410 ]
    ind = 1
    plotaucxs = [ 1:ind; ind * ones(410 - ind) ]
    plotaucys = [ aucdata[1:ind]; aucdata[ind] * ones(410 - ind) ]
    values = Observable(data)
    values2 = Observable(data2)
    indvalues = Observable(ind)
    valueaucxs = Observable(plotaucxs)
    valueaucys = Observable(plotaucys)
    
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    lines!(ax1, 1:410, ep; color=:black, linewidth=1)
    lines!(ax1, 1:410, values; color=:blue, linewidth=1)
    band!(ax1, 1:410, zeros(410), values2; color=( :blue, 0.5 ))
    lines!(ax2, valueaucxs, valueaucys; color=:blue, linewidth=1)
    vlines!(ax2, indvalues; color=:blue, linewidth=1)
    scatter!(ax2, 410, 4; markersize=0)
    Label(
        fig[1, 0], "Infections and Immunity"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[2, 0], "Infections averted"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[3, 1], "Time"; 
        fontsize=11.84, tellwidth=false,
    )
    
    record(fig, videosdir("staticimm.mp4")) do io
        for τ ∈ 1:2:410
            values[ ] = [ fixedimmunity(normalepidemic, t, τ) * 0.02565 for t ∈ 1:410 ]
            values2[ ] = [ minepiimmunity(normalepidemic, fixedimmunity, t, τ; multiplier=0.02565) for t ∈ 1:410 ]
            indvalues[ ] = τ
            valueaucxs[ ] = [ 1:τ; τ * ones(410 - τ) ]
            valueaucys[ ] = [ aucdata[1:τ]; aucdata[τ] * ones(410 - τ) ]
            recordframe!(io)                               
        end
    end
end

staticfig = let 
    aucdata = [ auc(normalepidemic, fixedimmunity, x; multiplier=0.02565) for x ∈ 1:410 ]

    ep = [ normalepidemic(t) for t ∈ 1:410 ]
    data = [ 
        (
            v = fixedimmunity(normalepidemic, t, 200);
            v == 1 ? 0.02565 : v
        )
        for t ∈ 1:410 
    ]
    data2 = [ minepiimmunity(normalepidemic, fixedimmunity, t, 200; multiplier=0.02565) for t ∈ 1:410 ]
    ind = 200
    values = Observable(data)
    values2 = Observable(data2)
    indvalues = Observable(ind)
    
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    lines!(ax1, 1:410, ep; color=:black, linewidth=1)
    lines!(ax1, 1:410, values; color=:blue, linewidth=1)
    band!(ax1, 1:410, zeros(410), values2; color=( :blue, 0.5 ))
    lines!(ax2, 1:410, aucdata; color=:blue, linewidth=1)
    vlines!(ax2, indvalues; color=:blue, linewidth=1)
    scatter!(ax2, 410, 4; markersize=0)

    Label(
        fig[1, 0], "Infections and Immunity"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[2, 0], "Infections averted"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[3, 1], "Time"; 
        fontsize=11.84, tellwidth=false,
    )

    fig
end
safesave(plotsdir("staticfig.jpg"), staticfig)

let 
    aucdata = [ auc(normalepidemic, boostedimmunity, x; multiplier=0.02565) for x ∈ 1:410 ]

    ep = [ normalepidemic(t) for t ∈ 1:410 ]
    data = [ 
        (
            v = boostedimmunity(normalepidemic, t, 1);
            v == 1 ? 0.02565 : v
        )
        for t ∈ 1:410 
    ]
    data2 = [ minepiimmunity(normalepidemic, boostedimmunity, t, 1; multiplier=0.02565) for t ∈ 1:410 ]
    ind = 1
    plotaucxs = [ 1:ind; ind * ones(410 - ind) ]
    plotaucys = [ aucdata[1:ind]; aucdata[ind] * ones(410 - ind) ]
    values = Observable(data)
    values2 = Observable(data2)
    indvalues = Observable(ind)
    valueaucxs = Observable(plotaucxs)
    valueaucys = Observable(plotaucys)
    
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    lines!(ax1, 1:410, ep; color=:black, linewidth=1)
    lines!(ax1, 1:410, values; color=:blue, linewidth=1)
    band!(ax1, 1:410, zeros(410), values2; color=( :blue, 0.5 ))
    lines!(ax2, valueaucxs, valueaucys; color=:blue, linewidth=1)
    vlines!(ax2, indvalues; color=:blue, linewidth=1)
    scatter!(ax2, 410, 4; markersize=0)
    Label(
        fig[1, 0], "Infections and Immunity"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[2, 0], "Infections averted"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[3, 1], "Time"; 
        fontsize=11.84, tellwidth=false,
    )
    
    record(fig, videosdir("boostedimm.mp4")) do io
        for τ ∈ 1:2:410
            values[ ] = [ 
                (
                    v = boostedimmunity(normalepidemic, t, τ);
                    v == 1 ? 0.02565 : v
                )
                for t ∈ 1:410 
            ]
            values2[ ] = [ minepiimmunity(normalepidemic, boostedimmunity, t, τ; multiplier=0.02565) for t ∈ 1:410 ]
            indvalues[ ] = τ
            valueaucxs[ ] = [ 1:τ; τ * ones(410 - τ) ]
            valueaucys[ ] = [ aucdata[1:τ]; aucdata[τ] * ones(410 - τ) ]
            recordframe!(io)                               
        end
    end
end

boostedfig = let 
    aucdata = [ auc(normalepidemic, boostedimmunity, x; multiplier=0.02565) for x ∈ 1:410 ]

    ep = [ normalepidemic(t) for t ∈ 1:410 ]
    data = [ 
        (
            v = boostedimmunity(normalepidemic, t, 160);
            v == 1 ? 0.02565 : v
        )
        for t ∈ 1:410 
    ]
    data2 = [ minepiimmunity(normalepidemic, boostedimmunity, t, 160; multiplier=0.02565) for t ∈ 1:410 ]
    ind = 160
    values = Observable(data)
    values2 = Observable(data2)
    indvalues = Observable(ind)
    
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    lines!(ax1, 1:410, ep; color=:black, linewidth=1)
    lines!(ax1, 1:410, values; color=:blue, linewidth=1)
    band!(ax1, 1:410, zeros(410), values2; color=( :blue, 0.5 ))
    lines!(ax2, 1:410, aucdata; color=:blue, linewidth=1)
    vlines!(ax2, indvalues; color=:blue, linewidth=1)
    scatter!(ax2, 410, 4; markersize=0)

    Label(
        fig[1, 0], "Infections and Immunity"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[2, 0], "Infections averted"; 
        fontsize=11.84, rotation=π/2, tellheight=false,
    )
    Label(
        fig[3, 1], "Time"; 
        fontsize=11.84, tellwidth=false,
    )

    fig
end
safesave(plotsdir("boostedfig.jpg"), boostedfig)

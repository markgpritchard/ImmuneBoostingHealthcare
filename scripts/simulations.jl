
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CairoMakie
using PlotFormatting

include("loaddata.jl")
include("analysesimssetup.jl")

alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sinusoidal force of infection 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function coslambda(t, m=0.015, n=0.005)
    h = (m + n) / 2 
    d = (m - n) / 2 
    return h + d * cos(2π * t)
end

simfigurescoslambda = with_theme(theme_latexfonts()) do
    λc = [ coslambda((t - 308) / 365) for t ∈ 1:832]
    patients = zeros(832)
    u0 = zeros(16)
    u0[1] = 1

    fig = Figure(; size=( 500, 275 ))
    axs = [ 
        Axis(
            fig[1, j]; 
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for j ∈ 1:6
    ]

    vac = vaccinated
    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])
        p = HCWSEIIRRRp(
            0.075,  # βh 
            0.075,  # βp 
            0.5,  # η 
            0.2,  # γ 
            ψ,  # ψ
            0.01,  
            0.5
        )
        u0 = zeros(16)
        u0[1] = 1.0
        output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
        s = output[:, 1]
        ei = [ sum(@view output[t, 2:13]) for t ∈ 1:832 ]
        r = [ sum(@view output[t, 14:16]) for t ∈ 1:832 ]

        lines!(axs[j], 1:832, s[1:832]; color=COLOURVECTOR[1], linewidth=1, label="Susceptible")
        lines!(axs[j], 1:832, ei[1:832]; color=COLOURVECTOR[2], linewidth=1, label="Infected")
        lines!(axs[j], 1:832, r[1:832]; color=COLOURVECTOR[3], linewidth=1, label="Resistant")

        formataxis!(
            axs[j]; 
            hidex=true, hidexticks=true, hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t, :b )
        )
        if j != 1 hidespines!(axs[j], :l) end
    end 

    formataxis!(
        axs[6]; 
        hidex=true, hidexticks=true, hidey=true, hideyticks=true,
        trimspines=true, hidespines=( :l, :r, :t, :b )
    )

    linkaxes!(axs...)
    #linkxaxes!(axs..., laxs...)
    #for ax ∈ axs 
    #    for x ∈ [ 104, 288, 469, 653, 834 ]
    #        vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
    #    end
    #end
    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    #Label(fig[2, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    #Label(fig[1, 0], "Prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)

    #for r ∈ [ 1, 2 ] rowgap!(fig.layout, r, 5) end
    #colgap!(fig.layout, 1, 5)

    #colsize!(fig.layout, 6, Auto(0.1))

    #leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)


    #formataxis!(leg; horizontal=true)
        
    fig
end

#safesave(plotsdir("simfigurescoslambda.pdf"), simfigurescoslambda)


simdifffigurescoslambda = with_theme(theme_latexfonts()) do
    λc = [ coslambda((t - 308) / 365) for t ∈ 1:832]
    patients = zeros(832)
    u0 = zeros(16)
    u0[1] = 1

    fig2 = Figure(; size=( 500, 400 ))
    axs = [ 
        Axis(
            fig2[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:6 
    ]
    laxs = [ Axis(fig2[1:4, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:5 ]
    haxs = [ Axis(fig2[i, 1:5]) for i ∈ 1:4 ]
    
    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])
        p = HCWSEIIRRRp(
            0.075,  # βh 
            0.075,  # βp 
            0.5,  # η 
            0.2,  # γ 
            ψ,  # ψ
            0.01,  
            0.5
        )
    
        u0 = zeros(16)
        u0[1] = 1
        mainoutput = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vaccinated)
        mainip = [ mainoutput[t, 2] * 0.5 for t ∈ 1:832 ]
    
        for (i, vac) ∈ enumerate([ 
            alternativevaccinations["minus2months"],
            alternativevaccinations["minus1month"],
            alternativevaccinations["plus1month"],
            alternativevaccinations["plus2months"]
        ])
    
            u0 = zeros(16)
            u0[1] = 1
            output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
            ip = [ output[t, 2] * 0.5 for t ∈ 1:832 ]
    
            ipdiff = cumsum(ip[288:832] .- mainip[288:832]) .* 4572
            lines!(axs[i, j], 288:832, ipdiff; color=COLOURVECTOR[1], linewidth=1,)
    
            formataxis!(
                axs[i, j]; 
                hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i, j], :b) end
            if j != 1 hidespines!(axs[i, j], :l) end
        end
    end 

    for i ∈ 1:4 
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end
    
    text!(haxs[1, 1], 325, 500; text="2 months earlier", align=( :center, :center ), fontsize=10)    
    text!(haxs[2, 1], 325, 500; text="1 month earlier", align=( :center, :center ), fontsize=10)    
    text!(haxs[3, 1], 325, 600; text="1 month later", align=( :center, :center ), fontsize=10)    
    text!(haxs[3, 1], 325, -650; text="2 months later", align=( :center, :center ), fontsize=10)    
    
    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    linkxaxes!(haxs..., laxs[1])
    linkyaxes!(axs..., haxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    for ax ∈ haxs 
        hlines!(ax, 1; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
        setvalue!(ax, 288, 0)
        setvalue!(ax, 834, 0)
    end
    Label(fig2[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig2[5, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig2[1:4, 0], "Cumulative difference"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 5 ] rowgap!(fig2.layout, r, 5) end
    for r ∈ 2:4 rowgap!(fig2.layout, r, 7) end
    colgap!(fig2.layout, 1, 5)

    colsize!(fig2.layout, 6, Auto(0.1))
        
    fig2
end

safesave(plotsdir("simdifffigurescoslambda.pdf"), simdifffigurescoslambda)

##


datarbk = filter(:StringCodes => x -> x == "RBK", finaldata)
insertcols!(datarbk, :Community => community)

λc = 0.02 .* community

patients = let 
    pts = zeros(832)
    prev = 0.0
    for (t, p) ∈ enumerate(datarbk.PatientsProportion) 
        if prev == 0 
            pts[t] = p 
            prev = p
        elseif p > 5 * prev || p < 0.2 * prev
            pts[t] = pts[t-1]
        else 
            pts[t] = p 
            prev = p
        end 
    end
    pts
end

staff = let 
    stf = zeros(832)
    prev = 0.0
    for (t, p) ∈ enumerate(datarbk.StaffProportion) 
        if prev == 0 
            stf[t] = p 
            prev = p
        elseif p > 5 * prev || p < 0.2 * prev
            stf[t] = stf[t-1]
        else 
            stf[t] = p 
            prev = p
        end 
    end
    stf
end

u0 = zeros(16)
u0[1] = 1


simfigures = with_theme(theme_latexfonts()) do
    fig = simfigurescoslambda
    axs = [ 
        Axis(
            fig[2, j]; 
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for j ∈ 1:6
    ]
    laxs = [ Axis(fig[1:2, j]; xticks=[ 104, 288, 469, 653, 834 ]) for j ∈ 1:5 ]

    vac = vaccinated
    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])
        p = HCWSEIIRRRp(
            0.075,  # βh 
            0.075,  # βp 
            0.5,  # η 
            0.2,  # γ 
            ψ,  # ψ
            0.01,  
            0.5
        )
        u0 = zeros(16)
        u0[1] = 1.0
        output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
        s = output[:, 1]
        ei = [ sum(@view output[t, 2:13]) for t ∈ 1:832 ]
        r = [ sum(@view output[t, 14:16]) for t ∈ 1:832 ]

        lines!(axs[j], 1:832, s[1:832]; color=COLOURVECTOR[1], linewidth=1, label="Susceptible")
        lines!(axs[j], 1:832, ei[1:832]; color=COLOURVECTOR[2], linewidth=1, label="Infected")
        lines!(axs[j], 1:832, r[1:832]; color=COLOURVECTOR[3], linewidth=1, label="Resistant")

        formataxis!(
            axs[j]; 
            hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if j != 1 hidespines!(axs[j], :l) end
    end 

    formataxis!(
        axs[6]; 
        hidex=true, hidexticks=true, hidey=true, hideyticks=true,
        trimspines=true, hidespines=( :l, :r, :t, :b )
    )

    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    
    Label(fig[3, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig[1:2, 0], "Prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 3 ] rowgap!(fig.layout, r, 5) end
    colgap!(fig.layout, 1, 5)

    colsize!(fig.layout, 6, Auto(0.1))

    leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)


    formataxis!(leg; horizontal=true)
    
    

    fig
end

simdifffigures = with_theme(theme_latexfonts()) do
    fig2 = Figure(; size=( 500, 400 ))
    axs = [ 
        Axis(
            fig2[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:6 
    ]
    laxs = [ Axis(fig2[1:4, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:5 ]
    haxs = [ Axis(fig2[i, 1:5]) for i ∈ 1:4 ]
    
    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])
        p = HCWSEIIRRRp(
            0.075,  # βh 
            0.075,  # βp 
            0.5,  # η 
            0.2,  # γ 
            ψ,  # ψ
            0.01,  
            0.5
        )
    
        u0 = zeros(16)
        u0[1] = 1
        mainoutput = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vaccinated)
        mainip = [ mainoutput[t, 2] * 0.5 for t ∈ 1:832 ]
    
        for (i, vac) ∈ enumerate([ 
            alternativevaccinations["minus2months"],
            alternativevaccinations["minus1month"],
            alternativevaccinations["plus1month"],
            alternativevaccinations["plus2months"]
        ])
    
            u0 = zeros(16)
            u0[1] = 1
            output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
            ip = [ output[t, 2] * 0.5 for t ∈ 1:832 ]
    
            ipdiff = cumsum(ip[288:832] .- mainip[288:832]) .* 4572
            lines!(axs[i, j], 288:832, ipdiff; color=COLOURVECTOR[1], linewidth=1,)
    
            formataxis!(
                axs[i, j]; 
                hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i, j], :b) end
            if j != 1 hidespines!(axs[i, j], :l) end
        end
    end 

    for i ∈ 1:4 
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end
    
    text!(axs[1, 1], 298, 500; text="2 months earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[2, 1], 298, 500; text="1 month earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[3, 1], 298, 500; text="1 month later", align=( :left, :center ), fontsize=10)    
    text!(axs[4, 1], 298, 650; text="2 months later", align=( :left, :center ), fontsize=10)    
    
    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    linkyaxes!(axs..., haxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    for ax ∈ haxs 
        hlines!(ax, 1; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    Label(fig2[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig2[5, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig2[1:4, 0], "Cumulative difference"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 5 ] rowgap!(fig2.layout, r, 5) end
    for r ∈ 2:4 rowgap!(fig2.layout, r, 7) end
    colgap!(fig2.layout, 1, 5)

    colsize!(fig2.layout, 6, Auto(0.1))
    
    #[ 0, 1, 2, 5, 10 ]
    
    
    fig2
    
    

end

using Distributions, Pigeons, Random, Turing

@model function fitmodelonehospital( 
    patients, staff, vaccinated, community, stringency, ndates;
    betahprior=truncated(Exponential(1), 0, 10),
    betapprior=truncated(Exponential(1), 0, 10),
    alpha7prior=truncated(Normal(0.2, 1), -1, 10),
    alpha8prior=truncated(Normal(0.1, 0.1), 0, 10),  # require greater stringency leads to less transmission
    omegaprior=truncated(Exponential(0.02), 0, 0.33),
    psiprior=truncated(Exponential(1), 0, 1000),
    thetaprior=Beta(1, 1),
    sigma2prior=Exponential(1),
)
    βh ~ betahprior
    βp ~ betapprior
    α7 ~ alpha7prior
    α8 ~ alpha8prior
     
    if omegaprior isa Number 
        ω = omegaprior
    else
        ω ~ omegaprior
    end

    if psiprior isa Number 
        ψ = psiprior 
    else
        ψ ~ psiprior
    end

    θ ~ thetaprior
    sigma2 ~ sigma2prior

    T = typeof(α7)
   # u0 = zeros(16)
   # u0[1] = 1.0
    isolating = zeros(Float64, ndates)

    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community
    
    p = HCWSEIIRRRp(
        βh, 
        βp, 
        0.5, 
        0.2, 
        ψ, 
        ω, 
        θ
    )
        #println("$p")
    hcwseiirrr_isolating!(isolating, ImmuneBoostingHealthcare.automatic, p, 1:ndates, λc, patients, vaccinated, 1)

    for t ∈ 1:ndates
        if isnan(isolating[t])
            Turing.@addlogprob! -Inf
        else
            Turing.@addlogprob! logpdf(Normal(isolating[t], sigma2), staff[t])
        end
    end
end


function fitdatamodel_target(
    patients=patients, 
    staff=staff, 
    vaccinated=vaccinated, 
    community=community, 
    stringency=stringency, 
    ndates=832;
    omegaprior=truncated(Exponential(0.02), 0, 0.33),
)
    return Pigeons.TuringLogPotential(
        fitmodelonehospital(
            patients, staff, vaccinated, community, stringency, ndates;
            omegaprior
        )
    )
end

const FitdatamodelType = typeof(fitdatamodel_target())

function Pigeons.initialization(target::FitdatamodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :βh, 1, 0.075)
    Pigeons.update_state!(result, :βp, 1, 0.075)
    Pigeons.update_state!(result, :α7, 1, 0.1)
    Pigeons.update_state!(result, :α8, 1, 0.0)
    Pigeons.update_state!(result, :ψ, 1, 1.0)
    Pigeons.update_state!(result, :θ, 1, 2/7)
    Pigeons.update_state!(result, :sigma2, 1, 1.0)
    return result
end

fitted_pt = pigeons( ;
    target=fitdatamodel_target(
        patients, staff, vaccinated, community, stringency, ndates
    ),
    n_rounds=0,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(id),
    variational=GaussianReference(),
)

new_pt = fitted_pt

for i ∈ 1:5
    filename = "fittedvalues_onehospital_id_$(id)_round_$(i).jld2"
    nextfilename = "fittedvalues_onehospital_id_$(id)_round_$(i + 1).jld2"
    isfile(datadir("sims", nextfilename)) && continue
    if isfile(datadir("sims", filename))
        global new_pt = load(datadir("sims", filename))["pt"]
    else
        pt = increment_n_rounds!(new_pt, 1)
        global new_pt = pigeons(pt)
        new_chains = Chains(new_pt)
        resultdict = Dict(
            "chain" => new_chains, 
            "pt" => new_pt, 
            "n_rounds" => i, 
            "n_chains" => 4,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end

fittedvalues_onehospital = load(datadir("sims", "fittedvalues_onehospital_id_1_round_5.jld2"))
fittedvalues_onehospitaldf = DataFrame(fittedvalues_onehospital["chain"])


fitted_pt_omega001 = pigeons( ;
    target=fitdatamodel_target(
        patients, staff, vaccinated, community, stringency, ndates; omegaprior=0.01
    ),
    n_rounds=0,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(id),
    variational=GaussianReference(),
)

new_pt = fitted_pt_omega001

for i ∈ 1:5
    filename = "fittedvalues_onehospital_omega_001_id_$(id)_round_$(i).jld2"
    nextfilename = "fittedvalues_onehospital_omega_001_id_$(id)_round_$(i + 1).jld2"
    isfile(datadir("sims", nextfilename)) && continue
    if isfile(datadir("sims", filename))
        global new_pt = load(datadir("sims", filename))["pt"]
    else
        pt = increment_n_rounds!(new_pt, 1)
        global new_pt = pigeons(pt)
        new_chains = Chains(new_pt)
        resultdict = Dict(
            "chain" => new_chains, 
            "pt" => new_pt, 
            "n_rounds" => i, 
            "n_chains" => 4,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end

fittedvalues_onehospital_omega001 = load(datadir("sims", "fittedvalues_onehospital_omega_001_id_1_round_5.jld2"))
fittedvalues_onehospitaldf_omega001 = DataFrame(fittedvalues_onehospital_omega001["chain"])




fitdifffigures = with_theme(theme_latexfonts()) do
    fig2 = Figure(; size=( 500, 400 ))
    axs = [ 
        Axis(
            fig2[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:7
    ]
    laxs = [ Axis(fig2[1:4, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:6 ]
    haxs = [ Axis(fig2[i, 1:6]) for i ∈ 1:4 ]
    
    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10, 100 ])
        for k ∈ axes(fittedvalues_onehospitaldf, 1)
            if ψ == 100 
                _ψ = fittedvalues_onehospitaldf.ψ[k] 
            else 
                _ψ = ψ
            end
            p = HCWSEIIRRRp(
                fittedvalues_onehospitaldf.βh[k],  # βh 
                fittedvalues_onehospitaldf.βp[k],  # βp 
                0.5,  # η 
                0.2,  # γ 
                _ψ,  # ψ
                fittedvalues_onehospitaldf.ω[k],  
                fittedvalues_onehospitaldf.θ[k]
            )
    
            u0 = zeros(16)
            u0[1] = 1
            mainoutput = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vaccinated)
            mainip = [ mainoutput[t, 2] * 0.5 for t ∈ 1:832 ]
        
            for (i, vac) ∈ enumerate([ 
                alternativevaccinations["minus2months"],
                alternativevaccinations["minus1month"],
                alternativevaccinations["plus1month"],
                alternativevaccinations["plus2months"]
            ])
        
                u0 = zeros(16)
                u0[1] = 1
                output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
                ip = [ output[t, 2] * 0.5 for t ∈ 1:832 ]
        
                ipdiff = cumsum(ip[288:832] .- mainip[288:832]) .* 4572
                lines!(axs[i, j], 288:832, ipdiff; color=( COLOURVECTOR[1], 0.1 ), linewidth=1,)
            end
        end
        for i ∈ 1:4
            formataxis!(
                axs[i, j]; 
                hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i, j], :b) end
            if j != 1 hidespines!(axs[i, j], :l) end
        end
    end 

    for i ∈ 1:4 
        formataxis!(
            axs[i, 7]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end
    
    text!(axs[1, 1], 298, 500; text="2 months earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[2, 1], 298, 500; text="1 month earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[3, 1], 298, 500; text="1 month later", align=( :left, :center ), fontsize=10)    
    text!(axs[4, 1], 298, 650; text="2 months later", align=( :left, :center ), fontsize=10)    
    
    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    linkyaxes!(axs..., haxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    for ax ∈ haxs 
        hlines!(ax, 1; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    Label(fig2[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig2[5, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig2[1:4, 0], "Cumulative difference"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 5 ] rowgap!(fig2.layout, r, 5) end
    for r ∈ 2:4 rowgap!(fig2.layout, r, 7) end
    colgap!(fig2.layout, 1, 5)

    colsize!(fig2.layout, 7, Auto(0.1))
    
    #[ 0, 1, 2, 5, 10 ]
    
    
    fig2
    
    

end

fitdifffigures_omega001 = with_theme(theme_latexfonts()) do
    fig2 = Figure(; size=( 500, 400 ))
    axs = [ 
        Axis(
            fig2[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:7
    ]
    laxs = [ Axis(fig2[1:4, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:6 ]
    haxs = [ Axis(fig2[i, 1:6]) for i ∈ 1:4 ]
    
    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10, 100 ])
        for k ∈ axes(fittedvalues_onehospitaldf_omega001, 1)
            if ψ == 100 
                _ψ = fittedvalues_onehospitaldf_omega001.ψ[k] 
            else 
                _ψ = ψ
            end
            p = HCWSEIIRRRp(
                fittedvalues_onehospitaldf_omega001.βh[k],  # βh 
                fittedvalues_onehospitaldf_omega001.βp[k],  # βp 
                0.5,  # η 
                0.2,  # γ 
                _ψ,  # ψ
                0.01,  
                fittedvalues_onehospitaldf_omega001.θ[k]
            )
    
            u0 = zeros(16)
            u0[1] = 1
            mainoutput = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vaccinated)
            mainip = [ mainoutput[t, 2] * 0.5 for t ∈ 1:832 ]
        
            for (i, vac) ∈ enumerate([ 
                alternativevaccinations["minus2months"],
                alternativevaccinations["minus1month"],
                alternativevaccinations["plus1month"],
                alternativevaccinations["plus2months"]
            ])
        
                u0 = zeros(16)
                u0[1] = 1
                output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
                ip = [ output[t, 2] * 0.5 for t ∈ 1:832 ]
        
                ipdiff = cumsum(ip[288:832] .- mainip[288:832]) .* 4572
                lines!(axs[i, j], 288:832, ipdiff; color=( COLOURVECTOR[1], 0.1 ), linewidth=1,)
            end
        end
        for i ∈ 1:4
            formataxis!(
                axs[i, j]; 
                hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i, j], :b) end
            if j != 1 hidespines!(axs[i, j], :l) end
        end
    end 

    for i ∈ 1:4 
        formataxis!(
            axs[i, 7]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end
    
    text!(axs[1, 1], 298, 500; text="2 months earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[2, 1], 298, 500; text="1 month earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[3, 1], 298, 500; text="1 month later", align=( :left, :center ), fontsize=10)    
    text!(axs[4, 1], 298, 650; text="2 months later", align=( :left, :center ), fontsize=10)    
    
    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    linkyaxes!(axs..., haxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    for ax ∈ haxs 
        hlines!(ax, 1; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    Label(fig2[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig2[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig2[5, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig2[1:4, 0], "Cumulative difference"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 5 ] rowgap!(fig2.layout, r, 5) end
    for r ∈ 2:4 rowgap!(fig2.layout, r, 7) end
    colgap!(fig2.layout, 1, 5)

    colsize!(fig2.layout, 7, Auto(0.1))
    
    #[ 0, 1, 2, 5, 10 ]
    
    
    fig2
    
    

end

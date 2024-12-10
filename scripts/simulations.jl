
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using CairoMakie
using PlotFormatting

include("loaddata.jl")
include("analysesimssetup.jl")

alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))

include("loaddata.jl")
include("analysesimssetup.jl")

alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))

datarbk = filter(:StringCodes => x -> x == "RBK", finaldata)
insertcols!(datarbk, :Community => community)

#λc = 0.02 .* community

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

fittedvalues_onehospitaldf = let 
    fv = load(datadir("sims", "fittedvalues_onehospital_omega_001_id_1_round_12nn.jld2")) 
    df = DataFrame(fv["chain"])
    for i ∈ 2:4 
        fv = load(
            datadir("sims", "fittedvalues_onehospital_omega_001_id_$(i)_round_12nn.jld2")
        ) 
        tdf = DataFrame(fv["chain"]) 
        for j ∈ axes(tdf, 1)
            tdf.chain[j] = i 
        end 
        append!(df, tdf)
    end 
    df 
end

dataoutputs100 = load(datadir("sims", "dataoutputs100.jld2"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function coslambda(t, m=0.015, n=0.005)
    h = (m + n) / 2 
    d = (m - n) / 2 
    return h + d * cos(2π * t)
end

cosmodels = let
    d = Dict{String, Matrix{Float64}}()
    λc = [ coslambda((t - 319) / 365) for t ∈ 1:832 ]
    for ψ ∈ [ 0, 1, 2, 5, 10 ]
        for (vac, lbl) ∈ zip(
            [
                vaccinated,
                alternativevaccinations["minus2months"],
                alternativevaccinations["minus1month"],
                alternativevaccinations["plus1month"],
                alternativevaccinations["plus2months"],
            ], 
            [ "sept", "minus2", "minus1", "plus1", "plus2" ]
        )  
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
            output = runhcwseiirrr(u0, p, 1:1:832, λc, zeros(832), vac)
            push!(d, "model_ψ=$(ψ)_vac=$(lbl)" => output)
        end
    end
    d
end

betahmodels = let
    d = Dict{String, Matrix{Float64}}()
    λc = [ (t < 10 ? 1.0 : 0.1) * coslambda((t - 319) / 365) for t ∈ 1:832 ] 
    for ψ ∈ [ 0, 1, 2, 5, 10 ]
        for (vac, lbl) ∈ zip(
            [
                vaccinated,
                alternativevaccinations["minus2months"],
                alternativevaccinations["minus1month"],
                alternativevaccinations["plus1month"],
                alternativevaccinations["plus2months"],
            ], 
            [ "sept", "minus2", "minus1", "plus1", "plus2" ]
        )  
            p = HCWSEIIRRRp(
                0.4,  # βh 
                0.0,  # βp 
                0.5,  # η 
                0.2,  # γ 
                ψ,  # ψ
                0.01,  
                0.5
            )
            u0 = zeros(16)
            u0[1] = 1.0
            output = runhcwseiirrr(u0, p, 1:1:832, λc, zeros(832), vac)
            push!(d, "model_ψ=$(ψ)_vac=$(lbl)" => output)
        end
    end
    d
end

hospmodels = let
    d = Dict{String, Matrix{Float64}}()
    λc = 0.02 .* community
    for ψ ∈ [ 0, 1, 2, 5, 10 ]
        for (vac, lbl) ∈ zip(
            [
                vaccinated,
                alternativevaccinations["minus2months"],
                alternativevaccinations["minus1month"],
                alternativevaccinations["plus1month"],
                alternativevaccinations["plus2months"],
            ], 
            [ "sept", "minus2", "minus1", "plus1", "plus2" ]
        )  
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
            push!(d, "model_ψ=$(ψ)_vac=$(lbl)" => output)
        end
    end
    d
end

fitmodels = let
    d = Dict{String, Matrix{Float64}}()
    λc = 0.02 .* community
    for (vac, lbl) ∈ zip(
        [
            vaccinated,
            alternativevaccinations["minus2months"],
            alternativevaccinations["minus1month"],
            alternativevaccinations["plus1month"],
            alternativevaccinations["plus2months"],
        ], 
        [ "sept", "minus2", "minus1", "plus1", "plus2" ]
    )  
        modelleds = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        modelledei = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        modelledi0 = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        modelledr = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        modelledn = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        lambdacs = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        lambdahs = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        lambdaps = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        totallambdas = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        diagnosedprevalence = zeros(size(fittedvalues_onehospitaldf, 1), 832)
        for i ∈ axes(fittedvalues_onehospitaldf, 1)
            λc = [ 
                max(
                    0.0, 
                    (
                        fittedvalues_onehospitaldf.α7[i] + 
                        fittedvalues_onehospitaldf.α8[i] * (100 - s)
                    )
                ) 
                for s ∈ stringency 
            ] .* community
            p = HCWSEIIRRRp(
                fittedvalues_onehospitaldf.βh[i],  # βh 
                fittedvalues_onehospitaldf.βp[i],  # βp 
                0.5,  # η 
                0.2,  # γ 
                fittedvalues_onehospitaldf.ψ[i] ,  # ψ
                0.01,  
                fittedvalues_onehospitaldf.θ[i]
            )
            u0 = zeros(16)
            u0[1] = 1.0
            output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
            s = output[:, 1]
            iv = output[:, 3]
            ei = [ sum(@view output[t, 2:13]) for t ∈ 1:832 ]
            r = [ sum(@view output[t, 14:16]) for t ∈ 1:832 ]
            n = [ sum(@view output[t, 1:16]) for t ∈ 1:832 ]
            diagnosedprev = [ sum(@view output[t, 4:13]) for t ∈ 1:832 ]
            λh = [ fittedvalues_onehospitaldf.βh[i] * output[t, 3] for t ∈ 1:832 ]
            λp = fittedvalues_onehospitaldf.βp[i] .* patients
    
            modelleds[i, :] .= s
            modelledei[i, :] .= ei
            modelledi0[i, :] .= iv
            modelledr[i, :] .= r
            modelledn[i, :] .= n
            lambdacs[i, :] .= λc
            lambdahs[i, :] .= λh
            lambdaps[i, :] .= λp
            totallambdas[i, :] .= λc .+ λh .+ λp
            diagnosedprevalence[i, :] .= diagnosedprev
        end 
        push!(d, "model_vac=$(lbl)_modelleds" => modelleds)
        push!(d, "model_vac=$(lbl)_modelledei" => modelledei)
        push!(d, "model_vac=$(lbl)_modelledi0" => modelledi0)
        push!(d, "model_vac=$(lbl)_modelledr" => modelledr)
        push!(d, "model_vac=$(lbl)_modelledn" => modelledn)
        push!(d, "model_vac=$(lbl)_lambdacs" => lambdacs)
        push!(d, "model_vac=$(lbl)_lambdahs" => lambdahs)
        push!(d, "model_vac=$(lbl)_lambdaps" => lambdaps)
        push!(d, "model_vac=$(lbl)_totallambdas" => totallambdas)
        push!(d, "model_vac=$(lbl)_diagnosedprevalence" => diagnosedprevalence)
    end
    d
end

fitmodelsmultihospital = let
    d = Dict{String, Matrix{Float64}}()
    λc = 0.02 .* community
    for 
        (vac, lbl) ∈ zip(
            [
                vaccinated,
                alternativevaccinations["minus2months"],
                alternativevaccinations["minus1month"],
                alternativevaccinations["plus1month"],
                alternativevaccinations["plus2months"],
            ], 
            [ "sept", "minus2", "minus1", "plus1", "plus2" ]
        ),
        hosp ∈ 1:23  

        #modelleds = zeros(size(dataoutputs100["df"], 1), 832)
        #modelledei = zeros(size(dataoutputs100["df"], 1), 832)
        modelledi0 = zeros(size(dataoutputs100["df"], 1), 832)
        #modelledr = zeros(size(dataoutputs100["df"], 1), 832)
       # modelledn = zeros(size(dataoutputs100["df"], 1), 832)
        #lambdacs = zeros(size(dataoutputs100["df"], 1), 832)
        #lambdahs = zeros(size(dataoutputs100["df"], 1), 832)
        #lambdaps = zeros(size(dataoutputs100["df"], 1), 832)
        #totallambdas = zeros(size(dataoutputs100["df"], 1), 832)
        #diagnosedprevalence = zeros(size(dataoutputs100["df"], 1), 832)
        @info "hosp=$hosp"
        for i ∈ axes(dataoutputs100["df"], 1)
            λc = [ 
                max(
                    0.0, 
                    (
                        dataoutputs100["df"].α7[i] + 
                        dataoutputs100["df"].α8[i] * (100 - s)
                    )
                ) 
                for s ∈ stringency 
            ] .* community
            p = HCWSEIIRRRp(
                getproperty(dataoutputs100["df"], "betahs[$(hosp)]")[i],  # βh 
                getproperty(dataoutputs100["df"], "betaps[$(hosp)]")[i],  # βp 
                0.5,  # η 
                0.2,  # γ 
                dataoutputs100["df"].ψ[i] ,  # ψ
                0.01,  
                dataoutputs100["df"].θ[i]
            )
            u0 = zeros(16)
            u0[1] = 1.0
            output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
            #s = output[:, 1]
            iv = output[:, 3]
            #ei = [ sum(@view output[t, 2:13]) for t ∈ 1:832 ]
            #r = [ sum(@view output[t, 14:16]) for t ∈ 1:832 ]
            #n = [ sum(@view output[t, 1:16]) for t ∈ 1:832 ]
            #diagnosedprev = [ sum(@view output[t, 4:13]) for t ∈ 1:832 ]
            #λh = [ fittedvalues_onehospitaldf.βh[i] * output[t, 3] for t ∈ 1:832 ]
            #λp = fittedvalues_onehospitaldf.βp[i] .* patients
    
            #modelleds[i, :] .= s
            #modelledei[i, :] .= ei
            modelledi0[i, :] .= iv
            #modelledr[i, :] .= r
            #modelledn[i, :] .= n
            #lambdacs[i, :] .= λc
            #lambdahs[i, :] .= λh
            #lambdaps[i, :] .= λp
            #totallambdas[i, :] .= λc .+ λh .+ λp
            #diagnosedprevalence[i, :] .= diagnosedprev
            if i / 1000 == round(Int, i / 1000) print("$i ") end
        end 
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_modelleds" => modelleds)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_modelledei" => modelledei)
        push!(d, "model_hosp=$(hosp)_vac=$(lbl)_modelledi0" => modelledi0)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_modelledr" => modelledr)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_modelledn" => modelledn)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_lambdacs" => lambdacs)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_lambdahs" => lambdahs)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_lambdaps" => lambdaps)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_totallambdas" => totallambdas)
        #push!(d, "model_hosp=$(hosp)_vac=$(lbl)_diagnosedprevalence" => diagnosedprevalence)
    end
    d
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot compartments 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

compartmentfigures = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 700 ))
    axs = [ 
        Axis(
            fig[i, j]; 
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:9, j ∈ 1:6
    ]
    laxs = [ Axis(fig[1:9, j]; xticks=[ 104, 288, 469, 653, 834 ]) for j ∈ 1:5 ]

    i = 0
    for (m, model) ∈ enumerate([ cosmodels, betahmodels, hospmodels ])
        i += 1
        for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])
            lines!(
                axs[i, j], 1:832, model["model_ψ=$(ψ)_vac=sept"][:, 1]; 
                color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
            )
            lines!(
                axs[i, j], 
                1:832, 
                [ sum(@view model["model_ψ=$(ψ)_vac=sept"][k, 2:13]) for k ∈ 1:832 ]; 
                color=COLOURVECTOR[2], linewidth=1, label="Infected"
            )
            lines!(
                axs[i, j], 
                1:832, 
                [ sum(@view model["model_ψ=$(ψ)_vac=sept"][k, 14:16]) for k ∈ 1:832 ]; 
                color=COLOURVECTOR[3], linewidth=1, label="Resistant"
            )
        end
        i += 1 
        if m == 1 
            for j ∈ 1:5
                lines!(
                    axs[i, j], 1:832, [ coslambda((t - 319) / 365) for t ∈ 1:832 ]; 
                    color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
                )
            end
            i += 1 
            for j ∈ 1:5
                lines!(
                    axs[i, j], 1:832, zeros(832); 
                    color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
                )
                #axs[i, j].yticks = [ 0.0, 0.5, 1.0 ]
                setvalue!(axs[i, j], 104, 0.5)
            end
        else
            if m == 2 
                #λc = [ t < 350 ? 0.08 : 0.002 for t ∈ 1:832 ] .* community
                #λp = [ t < 180 ? 0.075 : 0.0075 for t ∈ 1:832 ] .* patients
                #λc = [ t < 30 ? 0.015 : 0.0 for t ∈ 1:832 ] 
                #λc = [ t < 10 ? coslambda((t - 319) / 365) : 0.0 for t ∈ 1:832 ] 
                λc = [ (t < 10 ? 1.0 : 0.1) * coslambda((t - 319) / 365) for t ∈ 1:832 ] 
                λp = 0.0#1 .* patients
                βh = 0.4
            else
                λc = 0.02 .* community
                λp = 0.075 .* patients
                βh = 0.075
            end
            for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])
                λh = βh .* model["model_ψ=$(ψ)_vac=sept"][:, 3]
                λt = λc .+ λp .+ λh
                lines!(
                    axs[i, j], 1:832, λt; 
                    color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
                )
                lines!(
                    axs[i+1, j], 1:832, λh ./ λt; 
                    color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
                )
            end
            i += 1
        end
    end
    setvalue!(axs[2, 1], 104, 0)


    for i ∈ 1:9, j ∈ 1:5
        formataxis!(
            axs[i, j];
            hidex=(i != 9), hidexticks=(i != 9), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 9 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end

    linkxaxes!(axs..., laxs...)
    for i ∈ 1:9 linkyaxes!(axs[i, :]...) end

    for ax ∈ laxs 
        for x ∈ [ 104, 288, 469, 653, 834 ]
            vlines!(
                ax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
            )
        end
        formataxis!(
            ax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end
    
    Label(fig[10, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    for i ∈ 1:3 
        Label(
            fig[(3 * i - 2), 0], "Prevalence"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[(3 * i - 1), 0], L"$\lambda$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[(3 * i), 0], L"$\lambda_h/\lambda$"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
    end

    colsize!(fig.layout, 6, Auto(0.1))
    for i ∈ 1:9 
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end

    leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)
    formataxis!(leg; horizontal=true)
    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    colgap!(fig.layout, 1, 5)
    for r ∈ [ 1, 2, 11 ] rowgap!(fig.layout, r, 5) end
    labelplots!([ "A", "B", "C" ], fig; rows=[ 1, 4, 7 ], padding=( 0, 5, 0, 0 ))

    fig
end

safesave(plotsdir("compartmentfigures.pdf"), compartmentfigures)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Effect of changing vaccination date 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simdifffigurescoslambda = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 350 ))
    axs = [ 
        Axis(
            fig[(2 * i), j]; 
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:6
    ]
    laxs = [ Axis(fig[2:8, j]; xticks=[ 104, 288, 469, 653, 834 ]) for j ∈ 1:5 ]
    haxs = [ Axis(fig[(2 * i), 1:5]; xticks=[ 104, 288, 469, 653, 834 ]) for i ∈ 1:4 ]

    for 
        (i, vac) ∈ enumerate([ "minus2", "minus1", "plus1", "plus2" ]),
        (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])

        lines!(
            axs[i, j], 
            469:832, 
            (
                cumsum(cosmodels["model_ψ=$(ψ)_vac=$(vac)"][469:832, 3]) .- 
                cumsum(cosmodels["model_ψ=$(ψ)_vac=sept"][469:832, 3])
            ); 
            color=COLOURVECTOR[1], linewidth=1, 
        )

        formataxis!(
            axs[i, j];
            hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 4 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end

    for i ∈ 1:4 
        hlines!(
            haxs[i], 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
        )
    end

    linkxaxes!(axs..., laxs...)
    for i ∈ 1:4
        linkyaxes!(axs[i, :]..., haxs[i])
    end

    for ax ∈ laxs 
        #for x ∈ [ 104, 288, 469, 653, 834 ]
        for x ∈ [ 469, 653, 834 ]
            vlines!(
                ax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
            )
        end
        formataxis!(
            ax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end
    
    Label(fig[9, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(
        fig[2:8, 0], "Difference in cumulative incidence"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )

    colsize!(fig.layout, 6, Auto(0.1))
    for i ∈ 1:4 
        formataxis!(
            haxs[i]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end

    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[1, 1], "2 months earlier"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[3, 1], "1 month earlier"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[5, 1], "1 month later"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[7, 1], "2 months later"; fontsize=11.84, halign=:left, tellwidth=false)
    
    colgap!(fig.layout, 1, 5)
    for r ∈ 1:9 rowgap!(fig.layout, r, 5) end

    fig
end

safesave(plotsdir("simdifffigurescoslambda.pdf"), simdifffigurescoslambda)

simdifffiguresbetah = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 350 ))
    axs = [ 
        Axis(
            fig[(2 * i), j]; 
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:6
    ]
    laxs = [ Axis(fig[2:8, j]; xticks=[ 104, 288, 469, 653, 834 ]) for j ∈ 1:5 ]
    haxs = [ Axis(fig[(2 * i), 1:5]; xticks=[ 104, 288, 469, 653, 834 ]) for i ∈ 1:4 ]

    for 
        (i, vac) ∈ enumerate([ "minus2", "minus1", "plus1", "plus2" ]),
        (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])

        lines!(
            axs[i, j], 
            469:832, 
            (
                cumsum(betahmodels["model_ψ=$(ψ)_vac=$(vac)"][469:832, 3]) .- 
                cumsum(betahmodels["model_ψ=$(ψ)_vac=sept"][469:832, 3])
            ); 
            color=COLOURVECTOR[1], linewidth=1, 
        )

        formataxis!(
            axs[i, j];
            hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 4 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end

    for i ∈ 1:4 
        hlines!(
            haxs[i], 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
        )
    end

    linkxaxes!(axs..., laxs...)
    for i ∈ 1:4
        linkyaxes!(axs[i, :]..., haxs[i])
    end

    for ax ∈ laxs 
        #for x ∈ [ 104, 288, 469, 653, 834 ]
        for x ∈ [ 469, 653, 834 ]
            vlines!(
                ax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
            )
        end
        formataxis!(
            ax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end
    
    Label(fig[9, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(
        fig[2:8, 0], "Difference in cumulative incidence"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )

    colsize!(fig.layout, 6, Auto(0.1))
    for i ∈ 1:4 
        formataxis!(
            haxs[i]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end

    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[1, 1], "2 months earlier"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[3, 1], "1 month earlier"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[5, 1], "1 month later"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[7, 1], "2 months later"; fontsize=11.84, halign=:left, tellwidth=false)
    
    colgap!(fig.layout, 1, 5)
    for r ∈ 1:9 rowgap!(fig.layout, r, 5) end

    fig
end

safesave(plotsdir("simdifffiguresbetah.pdf"), simdifffiguresbetah)

simdifffigureshosp = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 350 ))
    axs = [ 
        Axis(
            fig[(2 * i), j]; 
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:4, j ∈ 1:6
    ]
    laxs = [ Axis(fig[2:8, j]; xticks=[ 104, 288, 469, 653, 834 ]) for j ∈ 1:5 ]
    haxs = [ Axis(fig[(2 * i), 1:5]; xticks=[ 104, 288, 469, 653, 834 ]) for i ∈ 1:4 ]

    for 
        (i, vac) ∈ enumerate([ "minus2", "minus1", "plus1", "plus2" ]),
        (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ])

        lines!(
            axs[i, j], 
            469:832, 
            (
                cumsum(hospmodels["model_ψ=$(ψ)_vac=$(vac)"][469:832, 3]) .- 
                cumsum(hospmodels["model_ψ=$(ψ)_vac=sept"][469:832, 3])
            ); 
            color=COLOURVECTOR[1], linewidth=1, 
        )

        formataxis!(
            axs[i, j];
            hidex=(i != 4), hidexticks=(i != 4), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 4 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end

    for i ∈ 1:4 
        hlines!(
            haxs[i], 0; 
            color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
        )
    end

    linkxaxes!(axs..., laxs...)
    for i ∈ 1:4
        linkyaxes!(axs[i, :]..., haxs[i])
    end

    for ax ∈ laxs 
        #for x ∈ [ 104, 288, 469, 653, 834 ]
        for x ∈ [ 469, 653, 834 ]
            vlines!(
                ax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
            )
        end
        formataxis!(
            ax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end
    
    Label(fig[9, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(
        fig[2:8, 0], "Difference in cumulative incidence"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )

    colsize!(fig.layout, 6, Auto(0.1))
    for i ∈ 1:4 
        formataxis!(
            haxs[i]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
    end

    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[1, 1], "2 months earlier"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[3, 1], "1 month earlier"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[5, 1], "1 month later"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[7, 1], "2 months later"; fontsize=11.84, halign=:left, tellwidth=false)
    
    colgap!(fig.layout, 1, 5)
    for r ∈ 1:9 rowgap!(fig.layout, r, 5) end

    fig
end

safesave(plotsdir("simdifffigureshosp.pdf"), simdifffigureshosp)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitted to one hospital
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fittedcompartmentfigures = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 350 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])
    gc = GridLayout(fig[1, 3])
    gd = GridLayout(fig[1, 4])

    let 
        axs = [ 
            Axis(
                ga[i, 1]; 
                xticks=( 
                    [ 104, 288, 469, 653, 834 ], 
                    [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
                ), 
                xticklabelrotation=-π/4,
            ) 
            for i ∈ 1:4 
        ]
        lax = Axis(ga[1:4, 1]; xticks=[ 104, 288, 469, 653, 834 ])

        sl = [ quantile(fitmodels["model_vac=sept_modelleds"][:, t], 0.05) for t ∈ 1:832 ]
        sm = [ quantile(fitmodels["model_vac=sept_modelleds"][:, t], 0.5) for t ∈ 1:832 ]
        su = [ quantile(fitmodels["model_vac=sept_modelleds"][:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[1], 1:832, sm; 
            color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
        )
        band!(
            axs[1], 1:832, sl, su; 
            color=( COLOURVECTOR[1], 0.5 ),
        )
        il = [ quantile(fitmodels["model_vac=sept_modelledei"][:, t], 0.05) for t ∈ 1:832 ]
        im = [ quantile(fitmodels["model_vac=sept_modelledei"][:, t], 0.5) for t ∈ 1:832 ]
        iu = [ quantile(fitmodels["model_vac=sept_modelledei"][:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[2], 1:832, im; 
            color=COLOURVECTOR[2], linewidth=1, label="Susceptible"
        )
        band!(
            axs[2], 1:832, il, iu; 
            color=( COLOURVECTOR[2], 0.5 ),
        )
        rl = [ quantile(fitmodels["model_vac=sept_modelledr"][:, t], 0.05) for t ∈ 1:832 ]
        rm = [ quantile(fitmodels["model_vac=sept_modelledr"][:, t], 0.5) for t ∈ 1:832 ]
        ru = [ quantile(fitmodels["model_vac=sept_modelledr"][:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[3], 1:832, rm; 
            color=COLOURVECTOR[3], linewidth=1, label="Susceptible"
        )
        band!(
            axs[3], 1:832, rl, ru; 
            color=( COLOURVECTOR[3], 0.5 ),
        )
        dl = [ 
            quantile(fitmodels["model_vac=sept_diagnosedprevalence"][:, t], 0.05) 
            for t ∈ 1:832 
        ]
        dm = [ 
            quantile(fitmodels["model_vac=sept_diagnosedprevalence"][:, t], 0.5) 
            for t ∈ 1:832 
        ]
        du = [ 
            quantile(fitmodels["model_vac=sept_diagnosedprevalence"][:, t], 0.95) 
            for t ∈ 1:832 
        ]
        lines!(
            axs[4], 1:832, dm; 
            color=COLOURVECTOR[2], linewidth=1, label="Susceptible"
        )
        band!(
            axs[4], 1:832, dl, du; 
            color=( COLOURVECTOR[2], 0.5 ),
        )
        scatter!(
            axs[4], 1:832, staff;
            markersize=1, color=:black
        )

        linkxaxes!(axs..., lax)
        for x ∈ [ 104, 288, 469, 653, 834 ]
            vlines!(
                lax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
            )
        end
        for i ∈ 1:4 
            formataxis!(
                axs[i];
                hidex=(i != 4), hidexticks=(i != 4),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i], :b) end
        end

        formataxis!(
            lax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        
        Label(ga[5, 1], "Date"; fontsize=11.84, tellwidth=false)
        Label(ga[1, 0], "Susceptible"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(ga[2, 0], "Infected"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(ga[3, 0], "Resistant"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(ga[4, 0], "Diagnosed"; fontsize=11.84, rotation=π/2, tellheight=false)
        colgap!(ga, 1, 5)
        rowgap!(ga, 4, 5)
    end
    let
        axs = [ 
            Axis(
                gb[i, 1]; 
                xticks=( 
                    [ 104, 288, 469, 653, 834 ], 
                    [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
                ), 
                xticklabelrotation=-π/4,
            ) 
            for i ∈ 1:4 
        ]
        lax = Axis(gb[1:4, 1]; xticks=[ 104, 288, 469, 653, 834 ])

        tl = [ quantile(fitmodels["model_vac=sept_totallambdas"][:, t], 0.05) for t ∈ 1:832 ]
        tm = [ quantile(fitmodels["model_vac=sept_totallambdas"][:, t], 0.5) for t ∈ 1:832 ]
        tu = [ quantile(fitmodels["model_vac=sept_totallambdas"][:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[1], 1:832, tm; 
            color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
        )
        band!(
            axs[1], 1:832, tl, tu; 
            color=( COLOURVECTOR[1], 0.5 ),
        )
        lc = [ 
            fitmodels["model_vac=sept_lambdacs"][i, j] == 0 ? 
                0.0 :
                (
                    fitmodels["model_vac=sept_lambdacs"][i, j] / 
                    fitmodels["model_vac=sept_totallambdas"][i, j] 
                )
            for 
                i ∈ axes(fitmodels["model_vac=sept_lambdacs"], 1), 
                j ∈ axes(fitmodels["model_vac=sept_lambdacs"], 2)
        ] 
        cl = [ quantile(lc[:, t], 0.05) for t ∈ 1:832 ]
        cm = [ quantile(lc[:, t], 0.5) for t ∈ 1:832 ]
        cu = [ quantile(lc[:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[2], 1:832, cm; 
            color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
        )
        band!(
            axs[2], 1:832, cl, cu; 
            color=( COLOURVECTOR[1], 0.5 ),
        )
        lp = [ 
            fitmodels["model_vac=sept_lambdaps"][i, j] == 0 ? 
                0.0 :
                (
                    fitmodels["model_vac=sept_lambdaps"][i, j] / 
                    fitmodels["model_vac=sept_totallambdas"][i, j] 
                )
            for 
                i ∈ axes(fitmodels["model_vac=sept_lambdaps"], 1), 
                j ∈ axes(fitmodels["model_vac=sept_lambdaps"], 2)
        ] 
        pl = [ quantile(lp[:, t], 0.05) for t ∈ 1:832 ]
        pm = [ quantile(lp[:, t], 0.5) for t ∈ 1:832 ]
        pu = [ quantile(lp[:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[3], 1:832, pm; 
            color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
        )
        band!(
            axs[3], 1:832, pl, pu; 
            color=( COLOURVECTOR[1], 0.5 ),
        )
        lh = [ 
            fitmodels["model_vac=sept_lambdahs"][i, j] == 0 ? 
                0.0 :
                (
                    fitmodels["model_vac=sept_lambdahs"][i, j] / 
                    fitmodels["model_vac=sept_totallambdas"][i, j] 
                )
            for 
                i ∈ axes(fitmodels["model_vac=sept_lambdahs"], 1), 
                j ∈ axes(fitmodels["model_vac=sept_lambdahs"], 2)
        ] 
        hl = [ quantile(lh[:, t], 0.05) for t ∈ 1:832 ]
        hm = [ quantile(lh[:, t], 0.5) for t ∈ 1:832 ]
        hu = [ quantile(lh[:, t], 0.95) for t ∈ 1:832 ]
        lines!(
            axs[4], 1:832, hm; 
            color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
        )
        band!(
            axs[4], 1:832, hl, hu; 
            color=( COLOURVECTOR[1], 0.5 ),
        )

        linkxaxes!(axs..., lax)
        for x ∈ [ 104, 288, 469, 653, 834 ]
            vlines!(
                lax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,
            )
        end
        for i ∈ 1:4 
            formataxis!(
                axs[i];
                hidex=(i != 4), hidexticks=(i != 4),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i], :b) end
        end

        formataxis!(
            lax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )
        
        Label(gb[5, 1], "Date"; fontsize=11.84, tellwidth=false)
        Label(gb[1, 0], L"Total $\lambda$"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(gb[2, 0], L"$\lambda_c/\lambda$"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(gb[3, 0], L"$\lambda_p/\lambda$"; fontsize=11.84, rotation=π/2, tellheight=false)
        Label(gb[4, 0], L"$\lambda_h/\lambda$"; fontsize=11.84, rotation=π/2, tellheight=false)
        colgap!(gb, 1, 5)
        rowgap!(gb, 4, 5)
    end
    let
        axs = [ 
            Axis(
                gc[(2 * i), 1]; 
                xticks=( 
                    [ 104, 288, 469, 653, 834 ], 
                    [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
                ), 
                xticklabelrotation=-π/4,
            ) 
            for i ∈ 1:4 
        ]
        lax = Axis(gc[2:8, 1]; xticks=[ 104, 288, 469, 653, 834 ])

        for (i, mod) ∈ enumerate([ "minus2", "minus1", "plus1", "plus2" ])
            dif = fitmodels["model_vac=$(mod)_modelledi0"] - fitmodels["model_vac=sept_modelledi0"]
            for j ∈ axes(dif, 1)
                for t ∈ 2:832
                    dif[j, t] += dif[j, t-1]
                end
            end
            dl = [ quantile(dif[:, t], 0.05) for t ∈ 1:832 ]
            dm = [ quantile(dif[:, t], 0.5) for t ∈ 1:832 ]
            du = [ quantile(dif[:, t], 0.95) for t ∈ 1:832 ]
            lines!(
                axs[i], 469:832, dm[469:832]; 
                color=COLOURVECTOR[1], linewidth=1, label="Susceptible"
            )
            band!(
                axs[i], 469:832, dl[469:832], du[469:832]; 
                color=( COLOURVECTOR[1], 0.5 ),
            )
            hlines!(
                axs[i], 0; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )
            @info "Model $mod, final distribution $(dm[832]) ($(dl[832]), $(du[832]))"
        end

        for x ∈ [ 469, 653, 834 ]
            vlines!(
                lax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )
        end
        for i ∈ 1:4 
            formataxis!(
                axs[i];
                hidex=(i != 4), hidexticks=(i != 4),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i], :b) end
        end

        formataxis!(
            lax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )

        Label(gc[1, 1], "2 months earlier"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(gc[3, 1], "1 month earlier"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(gc[5, 1], "1 month later"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(gc[7, 1], "2 months later"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(gc[9, 1], "Date"; fontsize=11.84, tellwidth=false)
        Label(
            gc[2:8, 0], "Difference in cumulative incidence"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        colgap!(gc, 1, 5)
        for r ∈ 1:8 rowgap!(gc, r, 5) end
    end

    colsize!(fig.layout, 4, Auto(0.04))

    labelplots!([ "A", "B", "C" ], [ ga, gb, gc ]; rows=1,)

    fig
end

safesave(plotsdir("fittedcompartmentfigures.pdf"), fittedcompartmentfigures)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitted to multiple hospitals 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

multiplehospitalfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 350 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    let
        axs = [ 
            Axis(
                ga[(2 * i), 1]; 
                xticks=( 
                    [ 104, 288, 469, 653, 834 ], 
                    [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
                ), 
                xticklabelrotation=-π/4,
            ) 
            for i ∈ 1:4 
        ]
        lax = Axis(ga[2:8, 1]; xticks=[ 104, 288, 469, 653, 834 ])

        for (i, mod) ∈ enumerate([ "minus2", "minus1", "plus1", "plus2" ])
            for k ∈ 1:23
                dif = (
                    fitmodelsmultihospital["model_hosp=$(k)_vac=$(mod)_modelledi0"] - 
                    fitmodelsmultihospital["model_hosp=$(k)_vac=sept_modelledi0"]
                )
                for j ∈ axes(dif, 1)
                    for t ∈ 2:832
                        dif[j, t] += dif[j, t-1]
                    end
                end
                #dl = [ quantile(dif[:, t], 0.05) for t ∈ 1:832 ]
                dm = [ quantile(dif[:, t], 0.5) for t ∈ 1:832 ]
                #du = [ quantile(dif[:, t], 0.95) for t ∈ 1:832 ]
                lines!(
                    axs[i], 469:832, dm[469:832]; 
                    #color=COLOURVECTOR[1], linewidth=1, label="Susceptible",
                    color=dataoutputs100["observationssincejuly"][i],
                    colorrange=extrema(dataoutputs100["observationssincejuly"]),
                    linewidth=1, 
                )
            end

            hlines!(
                axs[i], 0; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )
        end

        for x ∈ [ 469, 653, 834 ]
            vlines!(
                lax, x; 
                color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth=1,
            )
        end
        for i ∈ 1:4 
            formataxis!(
                axs[i];
                hidex=(i != 4), hidexticks=(i != 4),
                trimspines=true, hidespines=( :r, :t )
            )
            if i != 4 hidespines!(axs[i], :b) end
        end

        formataxis!(
            lax; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
            hidespines=( :l, :r, :t, :b )
        )

        Label(ga[1, 1], "2 months earlier"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(ga[3, 1], "1 month earlier"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(ga[5, 1], "1 month later"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(ga[7, 1], "2 months later"; fontsize=11.84, halign=:left, tellwidth=false)
        Label(ga[9, 1], "Date"; fontsize=11.84, tellwidth=false)
        Label(
            ga[2:8, 0], "Difference in cumulative incidence"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        colgap!(ga, 1, 5)
        for r ∈ 1:8 rowgap!(ga, r, 5) end
    end

    colsize!(fig.layout, 2, Auto(0.04))

    fig
end

safesave(plotsdir("multiplehospitalfig.pdf"), multiplehospitalfig)


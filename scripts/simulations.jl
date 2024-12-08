
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sinusoidal force of infection 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function coslambda(t, m=0.015, n=0.005)
    h = (m + n) / 2 
    d = (m - n) / 2 
    return h + d * cos(2π * t)
end

simfigurescoslambda = with_theme(theme_latexfonts()) do
    #λc = [ coslambda((t - 308) / 365) for t ∈ 1:832]
    λc = [ coslambda((t - 288) / 365) for t ∈ 1:832]
    patients = zeros(832)
    u0 = zeros(16)
    u0[1] = 1

    fig = Figure(; size=( 500, 300 ))
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
    #λc = [ coslambda((t - 308) / 365) for t ∈ 1:832]
    λc = [ coslambda((t - 319) / 365) for t ∈ 1:832]
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
        for x ∈ [ 104, 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    
    Label(fig[3, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig[1:2, 0], "Prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)

    colgap!(fig.layout, 1, 5)

    colsize!(fig.layout, 6, Auto(0.1))

    leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)
    for r ∈ [ 1, 2, 4 ] rowgap!(fig.layout, r, 5) end


    formataxis!(leg; horizontal=true)
    
    labelplots!([ "A", "B" ], fig; rows=[ 1, 2 ], padding=( 0, 5, 0, 0 ))

    fig
end

safesave(plotsdir("simfigures.pdf"), simfigures)



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


fittedvalues_onehospitaldf = let 
    fv = load(datadir("sims", "fittedvalues_onehospital_omega_001_id_1_round_12nn.jld2")) 
    df = DataFrame(fv["chain"])
    for i ∈ 2:4 
        fv = load(datadir("sims", "fittedvalues_onehospital_omega_001_id_$(i)_round_12nn.jld2")) 
        tdf = DataFrame(fv["chain"]) 
        for j ∈ axes(tdf, 1)
            tdf.chain[j] = i 
        end 
        append!(df, tdf)
    end 
    df 
end



fittedvalues_onehospitalchainfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 700 ))    
    chaxs = [ Axis(fig[(2 * i - 1), 1]) for i ∈ 1:6 ]
    deaxs = [ Axis(fig[(2 * i - 1), 3]) for i ∈ 1:6 ]
    for (i, v) ∈ enumerate([ :βh, :βp, :α7, :α8, :ψ, :θ ])
        for ch ∈ 3:-1:1
            d = filter(:chain => x -> x == ch, fittedvalues_onehospitaldf)
            lines!(chaxs[i], d.iteration, getproperty(d, v); color=COLOURVECTOR[ch], linewidth =1,)
            density!(
                deaxs[i], getproperty(d, v); 
                color=( :white, 0 ), strokecolor=COLOURVECTOR[ch], strokewidth=1
            )

        end
        Label(
            fig[(2 * i - 1), 0], 
            [ 
                L"$\beta_h$", 
                L"$\beta_p$", 
                L"$\alpha_1$", 
                L"$\alpha_2$", 
                L"$\psi$", 
                L"$\theta$" 
            ][i]; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[(2 * i - 1), 2], "Density"; 
            fontsize=11.84, rotation=π/2, tellheight=false
        )
        Label(
            fig[(2 * i), 3], 
            [ 
                L"$\beta_h$", 
                L"$\beta_p$", 
                L"$\alpha_1$", 
                L"$\alpha_2$", 
                L"$\psi$", 
                L"$\theta$" 
            ][i]; 
            fontsize=11.84, tellwidth=false
        )
        formataxis!(
            chaxs[i];
            hidex=(i != 6), hidexticks=(i !=6 ), trimspines=true, hidespines=( :t, :r )
        )
        if i != 6 hidespines!(chaxs[i], :b) end
        formataxis!(
            deaxs[i];
            trimspines=true, hidespines=( :t, :r )
        )

    end
    Label(fig[12, 1], "iteration"; fontsize=11.84, tellwidth=false)
    for c ∈ [ 1, 3 ] colgap!(fig.layout, c, 5) end 
    for ri ∈ 1:6 rowgap!(fig.layout, 2 * ri - 1, 5) end
    fig
end

safesave(plotsdir("fittedvalues_onehospitalchainfig.pdf"), fittedvalues_onehospitalchainfig)

fittedvalues_dict = let 
    modelleds = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    modelledei = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    modelledr = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    modelledn = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    lambdacs = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    lambdahs = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    lambdaps = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    totallambdas = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    diagnosedprevelence = zeros(size(fittedvalues_onehospitaldf, 1), 832)
    vac = vaccinated
    for i ∈ axes(fittedvalues_onehospitaldf, 1)
        λc = [ max(0.0, fittedvalues_onehospitaldf.α7[i] + fittedvalues_onehospitaldf.α8[i] * (100 - s)) for s ∈ stringency ] .* community
        p = HCWSEIIRRRp(
            fittedvalues_onehospitaldf.βh[i],  # βh 
            fittedvalues_onehospitaldf.βp[i],  # βp 
            0.5,  # η 
            0.2,  # γ 
            fittedvalues_onehospitaldf.ψ[i],  # ψ
            0.01,  
            fittedvalues_onehospitaldf.θ[i]
        )
        u0 = zeros(16)
        u0[1] = 1.0
        output = runhcwseiirrr(u0, p, 1:1:832, λc, patients, vac)
        s = output[:, 1]
        ei = [ sum(@view output[t, 2:13]) for t ∈ 1:832 ]
        r = [ sum(@view output[t, 14:16]) for t ∈ 1:832 ]
        n = [ sum(@view output[t, 1:16]) for t ∈ 1:832 ]
        diagnosedprev = [ sum(@view output[t, 4:13]) for t ∈ 1:832 ]
        λh = [ fittedvalues_onehospitaldf.βh[i] * output[t, 3] for t ∈ 1:832 ]
        λp = fittedvalues_onehospitaldf.βp[i] .* patients

        modelleds[i, :] .= s
        modelledei[i, :] .= ei
        modelledr[i, :] .= r
        modelledn[i, :] .= n
        lambdacs[i, :] .= λc
        lambdahs[i, :] .= λh
        lambdaps[i, :] .= λp
        totallambdas[i, :] .= λc .+ λh .+ λp
        diagnosedprevelence[i, :] .= diagnosedprev
    end 

    @dict modelleds modelledei modelledr modelledn lambdacs lambdahs lambdaps totallambdas diagnosedprevelence
end

fittedvaluesfig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 500 ))    
    axs = [ 
        Axis(fig[i,1];
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                    [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        )
        for i ∈ 1:6
    ]
    lax = Axis(fig[1:6, 1]; xticks=[ 104, 288, 469, 653, 834 ])

    for (i, v) ∈ enumerate([ 
        fittedvalues_dict[:modelleds], 
        fittedvalues_dict[:diagnosedprevelence], 
        fittedvalues_dict[:totallambdas], 
        fittedvalues_dict[:lambdahs], 
        fittedvalues_dict[:lambdacs], 
        fittedvalues_dict[:lambdaps] 
    ])
        medv = zeros(832)
        lv = zeros(832)
        uv = zeros(832) 

        for t ∈ 1:832 
            l, m, u = quantile(v[:, t], [ 0.05, 0.5, 0.95 ])
            medv[t] = m 
            lv[t] = l
            uv[t] = u 
        end 

        if i <= 2 
            lines!(axs[i], 1:832, medv; color=COLOURVECTOR[1], linewidth=1)
            band!(axs[i], 1:832, lv, uv; color=( COLOURVECTOR[1], 0.5 ))
        else
            lines!(axs[i], 1:832, 1 .- exp.(.-medv); color=COLOURVECTOR[1], linewidth=1)
            band!(axs[i], 1:832, 1 .- exp.(.-lv), 1 .- exp.(.-uv); color=( COLOURVECTOR[1], 0.5 ))
        end


        formataxis!(
            axs[i];
            hidex=(i != 6), hidexticks=(i != 6), trimspines=true, hidespines=( :t, :r )
        )
        if i != 6 hidespines!(axs[i], :b) end
    end

    scatter!(axs[2], 1:832, staff; markersize=2, color=:black)

    for x ∈ [ 104, 288, 469, 653, 834 ]
        vlines!(lax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
    end
    formataxis!(
        lax;
        hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b )
    )

    linkxaxes!(axs..., lax)

    Label(fig[1, 0], "Susceptible"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[2, 0], "Diagnosed"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[3, 0], L"Total $\lambda$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[4, 0], L"$\lambda_h$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[5, 0], L"$\lambda_c$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[6, 0], L"$\lambda_p$"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[7, 1], "Date"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 6, 5)

#=
    medinc = zeros(832)
    linc = zeros(832)
    uinc = zeros(832) 
    for t ∈ 1:832 
        l, m, u = quantile(fittedvalues_dict[:diagnosedprevelence][:, t], [ 0.05, 0.5, 0.95 ])
        medinc[t] = m 
        linc[t] = l
        uinc[t] = u 
    end 

    fig = Figure(; size=( 500, 250 ))    
    ax = Axis(fig[1, 1])
    lines!(ax, 1:832, medinc; color=COLOURVECTOR[1], linewidth=1)
    band!(ax, 1:832, linc, uinc; color=( COLOURVECTOR[1], 0.5 ))
    scatter!(ax, 1:832, staff; markersize=3, color=:black)

    formataxis!(
        ax;
        trimspines=true, hidespines=( :t, :r )
    )=#
    fig
end



fittedcommunitylambdafig = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 250 ))    
    axs = [ 
        Axis(fig[i,1];
            xticks=( 
                [ 104, 288, 469, 653, 834 ], 
                    [ "July 2020", "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        )
        for i ∈ 1:4
    ]

    for (i, λ) ∈ enumerate([ fittedvalues_dict[:lambdahs], fittedvalues_dict[:lambdacs], fittedvalues_dict[:lambdaps] ])
        medλc = zeros(832)
        lλc = zeros(832)
        uλc = zeros(832) 

        for t ∈ 1:832 
            l, m, u = quantile(λ[:, t], [ 0.05, 0.5, 0.95 ])
            medλc[t] = m 
            lλc[t] = l
            uλc[t] = u 
        end 

        if i == 1 
            lines!(axs[i], 1:832, 1 .- exp.(.-medλc); color=COLOURVECTOR[1], linewidth=1)
            band!(axs[i], 1:832, 1 .- exp.(.-lλc), 1 .- exp.(.-uλc); color=( COLOURVECTOR[1], 0.5 ))
        else
            lines!(axs[i], 1:832, 1 .- exp.(.-medλc); color=COLOURVECTOR[1], linewidth=1)
            band!(axs[i], 1:832, 1 .- exp.(.-lλc), 1 .- exp.(.-uλc); color=( COLOURVECTOR[1], 0.5 ))
        end
        

        formataxis!(
            axs[i];
            hidex=(i != 4), hidexticks=(i != 4), trimspines=true, hidespines=( :t, :r )
        )
    end




   #= Label(
        fig[1, 0], L"$\lambda_c$"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        fig[2, 1], "Date";
        fontsize=11.84, tellwidth=false
    )
    formataxis!(
        ax;
        trimspines=true, hidespines=( :t, :r )
    )

    colgap!(fig.layout, 1, 5)  
    rowgap!(fig.layout, 1, 5) =#
    fig
end

safesave(plotsdir("fittedcommunitylambdafig.pdf"), fittedcommunitylambdafig)



prevalencefig = with_theme(theme_latexfonts()) do
    vac = vaccinated

    fig = Figure(; size=( 500, 275 ))

    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)

    let 
        λc = [ coslambda((t - 308) / 365) for t ∈ 1:832]
        patients = zeros(832)
        u0 = zeros(16)
        u0[1] = 1

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
    end
    
    let
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

        laxs = [ Axis(fig[1:2, j]; xticks=[ 104, 288, 469, 653, 834 ]) for j ∈ 1:5 ]
        for lax ∈ laxs
            for x ∈ [ 104, 288, 469, 653, 834 ]
                vlines!(lax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
            end
            formataxis!(
                lax;
                hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b )
            )   
        end

        formataxis!(
            axs[6];
            hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b )
        ) 

        leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)
        formataxis!(leg; horizontal=true)

        linkxaxes!(axs..., laxs...)
        
    end

    colsize!(fig.layout, 6, Auto(0.1))

    Label(fig[1:2, 0], "Compartment sizes"; fontsize=11.84, rotation=π/2, tellheight=false)
    Label(fig[3, 1:5], "Date"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    for r ∈ [ 1, 2, 4 ] rowgap!(fig.layout, r, 5) end
        
    fig
end













#############################################################



using DrWatson
using CairoMakie
using PlotFormatting

@quickactivate :ImmuneBoostingHealthcare

include("loaddata.jl")
include("analysesimssetup.jl")

alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fixed force of infection 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function fixedlambda(t, m=0.015, n=0.005)
    if t < 365 
        return m 
    elseif t >= 550 && t < 713 
        return m
    else 
        return n
    end 
end

## High 

λc = [ fixedlambda(t) for t ∈ 1:832]
patients = zeros(832)
u0 = zeros(16)
u0[1] = 1

simfigures = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 500 ))
    axs = [ 
        Axis(
            fig[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:5, j ∈ 1:6
    ]
    laxs = [ Axis(fig[1:5, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:5 ]

    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ]), (i, vac) ∈ enumerate(
        [ 
            alternativevaccinations["minus2months"],
            alternativevaccinations["minus1month"],
            vaccinated,
            alternativevaccinations["plus1month"],
            alternativevaccinations["plus2months"]
        ]
    )
        p = HCWSEIIRRRp(
            0.0,  # βh 
            0.0,  # βp 
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

        lines!(axs[i, j], 288:832, s[288:832]; color=COLOURVECTOR[1], linewidth=1, label="Susceptible")
        lines!(axs[i, j], 288:832, ei[288:832]; color=COLOURVECTOR[2], linewidth=1, label="Infected")
        lines!(axs[i, j], 288:832, r[288:832]; color=COLOURVECTOR[3], linewidth=1, label="Resistant")

        formataxis!(
            axs[i, j]; 
            hidex=(i != 5), hidexticks=(i != 5), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 5 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end 

    for i ∈ 1:5 
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end

    text!(axs[1, 1], 288, 1; text="2 months earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[2, 1], 288, 1; text="1 month earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[4, 1], 288, 1; text="1 month later", align=( :left, :center ), fontsize=10)    
    text!(axs[5, 1], 288, 1; text="2 months later", align=( :left, :center ), fontsize=10)    

    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig[6, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig[1:5, 0], "Prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 6 ] rowgap!(fig.layout, r, 5) end
    colgap!(fig.layout, 1, 5)

    colsize!(fig.layout, 6, Auto(0.1))

    leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)


    formataxis!(leg; horizontal=true)
    #[ 0, 1, 2, 5, 10 ]
    
    

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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sinusoidal force of infection 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function coslambda(t, m=0.015, n=0.005)
    h = (m + n) / 2 
    d = (m - n) / 2 
    return h + d * cos(2π * t)
end

## High 

λc = [ coslambda((t - 288) / 365) for t ∈ 1:832]
patients = zeros(832)
u0 = zeros(16)
u0[1] = 1

simfigures = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 500 ))
    axs = [ 
        Axis(
            fig[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:5, j ∈ 1:6
    ]
    laxs = [ Axis(fig[1:5, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:5 ]

    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ]), (i, vac) ∈ enumerate(
        [ 
            alternativevaccinations["minus2months"],
            alternativevaccinations["minus1month"],
            vaccinated,
            alternativevaccinations["plus1month"],
            alternativevaccinations["plus2months"]
        ]
    )
        p = HCWSEIIRRRp(
            0.0,  # βh 
            0.0,  # βp 
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

        lines!(axs[i, j], 288:832, s[288:832]; color=COLOURVECTOR[1], linewidth=1, label="Susceptible")
        lines!(axs[i, j], 288:832, ei[288:832]; color=COLOURVECTOR[2], linewidth=1, label="Infected")
        lines!(axs[i, j], 288:832, r[288:832]; color=COLOURVECTOR[3], linewidth=1, label="Resistant")

        formataxis!(
            axs[i, j]; 
            hidex=(i != 5), hidexticks=(i != 5), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 5 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end 

    for i ∈ 1:5 
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end

    text!(axs[1, 1], 288, 1; text="2 months earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[2, 1], 288, 1; text="1 month earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[4, 1], 288, 1; text="1 month later", align=( :left, :center ), fontsize=10)    
    text!(axs[5, 1], 288, 1; text="2 months later", align=( :left, :center ), fontsize=10)    

    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig[6, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig[1:5, 0], "Prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 6 ] rowgap!(fig.layout, r, 5) end
    colgap!(fig.layout, 1, 5)

    colsize!(fig.layout, 6, Auto(0.1))

    leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)


    formataxis!(leg; horizontal=true)
    #[ 0, 1, 2, 5, 10 ]
    
    

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
    
    #[ 0, 1, 2, 5, 10 ]
    
    
    fig2
    
    

end


##


datarbk = filter(:StringCodes => x -> x == "RBK", finaldata["finaldata"])
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

u0 = zeros(16)
u0[1] = 1


simfigures = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 500, 500 ))
    axs = [ 
        Axis(
            fig[i, j]; 
            xticks=( 
                [ 288, 469, 653, 834 ], 
                [ "Jan. 2021", "July 2021", "Jan. 2022", "July 2022" ] 
            ), 
            xticklabelrotation=-π/4,
        ) 
        for i ∈ 1:5, j ∈ 1:6
    ]
    laxs = [ Axis(fig[1:5, j]; xticks=[ 288, 469, 653, 834 ]) for j ∈ 1:5 ]

    for (j, ψ) ∈ enumerate([ 0, 1, 2, 5, 10 ]), (i, vac) ∈ enumerate(
        [ 
            alternativevaccinations["minus2months"],
            alternativevaccinations["minus1month"],
            vaccinated,
            alternativevaccinations["plus1month"],
            alternativevaccinations["plus2months"]
        ]
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
        s = output[:, 1]
        ei = [ sum(@view output[t, 2:13]) for t ∈ 1:832 ]
        r = [ sum(@view output[t, 14:16]) for t ∈ 1:832 ]

        lines!(axs[i, j], 288:832, s[288:832]; color=COLOURVECTOR[1], linewidth=1, label="Susceptible")
        lines!(axs[i, j], 288:832, ei[288:832]; color=COLOURVECTOR[2], linewidth=1, label="Infected")
        lines!(axs[i, j], 288:832, r[288:832]; color=COLOURVECTOR[3], linewidth=1, label="Resistant")

        formataxis!(
            axs[i, j]; 
            hidex=(i != 5), hidexticks=(i != 5), hidey=(j != 1), hideyticks=(j != 1),
            trimspines=true, hidespines=( :r, :t )
        )
        if i != 5 hidespines!(axs[i, j], :b) end
        if j != 1 hidespines!(axs[i, j], :l) end
    end 

    for i ∈ 1:5 
        formataxis!(
            axs[i, 6]; 
            hidex=true, hidexticks=true, hidey=true, hideyticks=true,
            trimspines=true, hidespines=( :l, :r, :t, :b )
        )
    end

    text!(axs[1, 1], 288, 1; text="2 months earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[2, 1], 288, 1; text="1 month earlier", align=( :left, :center ), fontsize=10)    
    text!(axs[4, 1], 288, 1; text="1 month later", align=( :left, :center ), fontsize=10)    
    text!(axs[5, 1], 288, 1; text="2 months later", align=( :left, :center ), fontsize=10)    

    linkaxes!(axs...)
    linkxaxes!(axs..., laxs...)
    for ax ∈ laxs 
        for x ∈ [ 288, 469, 653, 834 ]
            vlines!(ax, x; color=RGBAf(0, 0, 0, 0.12), linestyle=( :dot, :dense ), linewidth =1,)
        end
        formataxis!(ax; hidex=true, hidexticks=true, hidey=true, hideyticks=true, hidespines=( :l, :r, :t, :b ))
    end
    Label(fig[0, 1], L"$\psi=0$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 2], L"$\psi=1$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 3], L"$\psi=2$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 4], L"$\psi=5$"; fontsize=11.84, halign=:left, tellwidth=false)
    Label(fig[0, 5], L"$\psi=10$"; fontsize=11.84, halign=:left, tellwidth=false)
    
    Label(fig[6, 1:5], "Date"; fontsize=11.84, tellwidth=false)
    Label(fig[1:5, 0], "Prevalence"; fontsize=11.84, rotation=π/2, tellheight=false)

    for r ∈ [ 1, 6 ] rowgap!(fig.layout, r, 5) end
    colgap!(fig.layout, 1, 5)

    colsize!(fig.layout, 6, Auto(0.1))

    leg = fig[-1, 1:5] = Legend(fig, axs[1, 1]; orientation=:horizontal)


    formataxis!(leg; horizontal=true)
    #[ 0, 1, 2, 5, 10 ]
    
    

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









@recipe function f(ms::MassSpectrum; max_labels=15, font_size=8)
    mass = ms.mass
    intensity = ms.intensity
    ylim --> (0.0, 1.05maximum(intensity))
    legend --> false
    @series begin
        seriestype := :sticks
        ms.mass, ms.intensity
    end
    @series begin
        max_labels = min(length(mass), max_labels)
        ordering = sortperm(intensity; rev=true)[1:max_labels]
        mass = ms.mass[ordering]
        intensity = ms.intensity[ordering]
        seriestype := :scatter
        markercolor --> [:black]
        markersize --> 0
        markeralpha --> 0.0
        series_annotation --> Main.Plots.text.([mass_format(m) for m in mass], font_size)
        mass, intensity
    end
end

@recipe function f(mss::AbstractArray{<:MassSpectrum})
    for (i, ms) in enumerate(mss)
        @series begin
            label --> str(i)
            ms
        end
    end
    ylim := maximum(intensity(ms) for ms in mss)
end

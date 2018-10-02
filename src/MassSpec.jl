module MassSpec

using CSV
using DataFrames
using Printf
using RecipesBase
using StatsBase

import Statistics

struct MassSpec{T<:Number,S<:Number} <: AbstractArray{T,1}
    mass::Array{T,1}
    intensity::Array{S,1}
end

Base.size(ms::MassSpec) = (size(ms.mass)[1],)
Base.getindex(ms::MassSpec, i::Integer) = (ms.mass[i], ms.intensity[i])
Base.getindex(ms::MassSpec, i) = MassSpec(ms.mass[i], ms.intensity[i])

function Base.setindex!(ms::MassSpec, val, i::Integer)
    setindex!(ms.mass, val[1], i)
    setindex!(ms.intensity, val[2], i)
end

function Base.setindex!(ms::MassSpec, val, i::Union{AbstractRange, Colon})
    setindex!(ms.mass, val[:,1], i)
    setindex!(ms.intensity, val[:,2], i)
end

function Base.deleteat!(ms::MassSpec, x)
    deleteat!(ms.mass, x)
    deleteat!(ms.intensity, x)
end


function MassSpec(path::AbstractString) 
    df = CSV.read(path, header=["mass", "mV"], datarow=2, delim=',', allowmissing=:none)
    MassSpec(df.mass, df.mV)
end

Statistics.mean(ms::MassSpec) = Statistics.mean(ms.intensity)

@recipe function f(ms::MassSpec; max_labels=15, font_size=8)
    mass = ms.mass
    intensity = ms.intensity
    ylim := (0.0, 1.05maximum(intensity))
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

function discretize(ms::MassSpec, rng)
    wv = weights(ms.intensity)
    fit(Histogram, ms.mass, wv, rng; closed=:left).weights
end

mass_format(m) = Printf.@sprintf "%.1f" m

export MassSpec, discretize

end

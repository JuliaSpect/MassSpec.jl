module MassSpec

import CSV
import Statistics

using DataFrames
using Printf
using RecipesBase
using StatsBase

struct MassSpectrum{T<:Number} <: AbstractArray{Tuple{T,T},1}
    mass::Array{T,1}
    intensity::Array{T,1}
end

intensity(ms::MassSpectrum) = maximum(ms.intensity)
Base.sum(ms::MassSpectrum) = sum(ms.intensity)
Base.size(ms::MassSpectrum) = (size(ms.mass)[1],)
Base.getindex(ms::MassSpectrum, i::Union{Integer, CartesianIndex}) = (ms.mass[i], ms.intensity[i])
Base.getindex(ms::MassSpectrum, i) = MassSpectrum(ms.mass[i], ms.intensity[i])

function Base.setindex!(ms::MassSpectrum, val, i::Integer)
    setindex!(ms.mass, val[1], i)
    setindex!(ms.intensity, val[2], i)
end

function Base.setindex!(ms::MassSpectrum, val, i::Union{AbstractRange, Colon})
    setindex!(ms.mass, val[:,1], i)
    setindex!(ms.intensity, val[:,2], i)
end

function Base.deleteat!(ms::MassSpectrum, x)
    deleteat!(ms.mass, x)
    deleteat!(ms.intensity, x)
end


function MassSpectrum(path::AbstractString)
    df = CSV.read(path, header=["mass", "mV"], datarow=2, delim=',', allowmissing=:none)
    MassSpectrum(df.mass, df.mV)
end

mass_format(m) = Printf.@sprintf "%.1f" m

include("operations.jl")
include("plotting.jl")

export MassSpectrum, Mass, discretize, subtract, distance

end

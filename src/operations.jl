import LinearAlgebra

"""Remove peaks near `mass` from spectrum `ms`.
Peaks within `tol` of `mass` are removed."""
function subtract(ms::MassSpectrum, mass, tol)
    ms[findall(peak -> abs(peak[1] - mass) > tol, ms)]
end

function subtract(m1::MassSpectrum, m2::MassSpectrum, tol)
    reduce(ms -> subtract(ms, m, tol), m2.mass, m1)
end

"""
Calculated the distance between two mass spectra.

The return value is between 0.0 and 1.0 and is calculated based
on the modelling each spectrum's peaks two-dimensional Gaussian in the
mass-intensity plane and subtracting the two. The square of the resulting
difference signal has a nice closed-form integral over the mass-intensity plane.
"""
function distance(ms1::MassSpectrum, ms2::MassSpectrum, σmass, σintensity)
    mass = [ms1.mass; ms2.mass]
    intensity = [ms1.intensity; ms2.intensity]
    α = [ms1.intensity; -1 .* ms2.intensity]
    d = 0.0
    max_d = 0.0
    for i in eachindex(mass)
        for j in eachindex(mass)
            mass_contrib = -(mass[i] - mass[j])^2 / (4σmass^2)
            intensity_contrib = -(intensity[i] - intensity[j])^2 / (4σintensity^2)
            term = α[i] * α[j] * exp(mass_contrib + intensity_contrib)
            d += term
            if term > 0.0
                max_d += term
            end
        end
    end
    max(0.0, d) / max_d
end

function discretize(ms::MassSpectrum, rng)
    wv = weights(ms.intensity)
    fit(Histogram, ms.mass, wv, rng; closed=:left).weights
end

"""
Convert binned intensity data to a mass spectrum.
"""
function undiscretize(mass_bins::AbstractArray, intensities::AbstractArray, intensity_threshold = 0.05maximum(intensities))
    indices = intensities .> intensity_threshold
    MassSpectrum(mass_bins[indices], intensities[indices])
end

function LinearAlgebra.normalize!(ms::MassSpectrum, p::Real=1)
    LinearAlgebra.normalize!(ms.intensity, p)
    ms
end

LinearAlgebra.normalize(ms::MassSpectrum, p::Real=1) = LinearAlgebra.normalize!(deepcopy(ms), p)

Statistics.mean(ms::MassSpectrum) = Statistics.mean(ms.intensity)

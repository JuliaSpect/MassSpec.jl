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
function distance(ms1::MassSpectrum{T}, ms2::MassSpectrum{T}, σmass, σintensity) where T
    mass = [ms1.mass; ms2.mass]
    intensity = [ms1.intensity; ms2.intensity]
    α = [ms1.intensity; -1 .* ms2.intensity]
    d = zero(T)
    max_d = zero(T)
    for i in eachindex(mass)
        for j in eachindex(mass)
            mass_contrib = -(mass[i] - mass[j])^2 / (4σmass^2)
            intensity_contrib = -(intensity[i] - intensity[j])^2 / (4σintensity^2)
            term = α[i] * α[j] * exp(mass_contrib + intensity_contrib)
            d += term
            if term > zero(T)
                max_d += term
            end
        end
    end
    max(zero(T), d) / max_d
end

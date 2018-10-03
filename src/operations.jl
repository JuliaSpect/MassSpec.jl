"""Remove peaks near `mass` from spectrum `ms`.
Peaks within `tol` of `mass` are removed."""
function subtract(ms::MassSpectrum, mass, tol)
    ms[findall(peak -> abs(peak[1] - mass) > tol, ms)]
end

function subtract(m1::MassSpectrum, m2::MassSpectrum, tol)
    reduce(ms -> subtract(ms, m, tol), m2.mass, m1)
end

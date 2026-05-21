"""
    hellwarth_b_scheme(modes)

Hellwarth et al. 1999 PRB B-scheme effective frequency for a matrix whose first
column is LO phonon frequency and second column is infrared activity.
"""
function hellwarth_b_scheme(modes::AbstractMatrix{<:Real})
    size(modes, 2) >= 2 || throw(ArgumentError("modes must have frequency and activity columns."))
    activity2 = modes[:, 2] .^ 2
    return sqrt(sum(activity2) / sum(activity2 ./ (modes[:, 1] .^ 2)))
end

"""
    hellwarth_a_scheme(modes; temperature=295, convergence=1e-6)

Temperature-dependent Hellwarth et al. 1999 PRB A-scheme effective frequency.
"""
function hellwarth_a_scheme(modes::AbstractMatrix{<:Real}; temperature::Real = 295, convergence::Real = 1e-6)
    frequencies = modes[:, 1]
    activities = modes[:, 2]
    target = sum(activities .* coth.(π .* frequencies .* 1e12 .* hbar ./ (kB * temperature)) ./ frequencies) / sum(activities)
    condition(f) = coth(π * f * 1e12 * hbar / (kB * temperature)) / f - target
    lo = minimum(frequencies)
    hi = maximum(frequencies)
    mid = (lo + hi) / 2
    while (hi - lo) / 2 > convergence
        if sign(condition(mid)) == sign(condition(lo))
            lo = mid
        else
            hi = mid
        end
        mid = (lo + hi) / 2
    end
    return mid
end

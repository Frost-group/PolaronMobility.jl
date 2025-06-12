# Wrap PlotPolaron (and it's dependency on Plots) within it's own (sub)module.
module PlotPolaron

export plotpolaron # Should unify case here c.f. Module 

using PolaronMobility
using Plots

function plotpolaron(fileprefix, p::FrohlichPolaron; extension="png")
    println("Plotting polaron to $fileprefix...")

    #####
    ## Mass vs. Temperature plot
    plot(p.T, p.M, label="Phonon effective-mass",
        markersize=3, marker=:rect,
        xlab="Temperature (K)", ylab="Phonon effective-mass",
        ylim=(0, 1.2))

    savefig("$fileprefix-mass.$extension")

    #####
    ## Relaxationtime vs. Temperature plot
    plot(p.T, p.τ, label="Kadanoff relaxation time (ps)", markersize=3, marker=:rect, xlab="Temperature (K)", ylab="Relaxation time (ps)", ylim=(0, 1.2))

    savefig("$fileprefix-tau.$extension")

    ## Mass + relaxation time vs. Temperature plot
    plot(p.T, p.M, label="Phonon effective-mass (m\$_b\$)", markersize=3, marker=:rect,
        xlab="Temperature (K)", ylab="Effective-mass / relaxation time", ylim=(0, 1.2))
    plot!(p.T, p.τ, label="Kadanoff relaxation time (ps)", markersize=3, marker=:diamond,
        xlab="Temperature (K)", ylab="Relaxation time (ps)", ylim=(0, 1.2))

    savefig("$fileprefix-mass-tau.$extension")

    ####
    ## Variational parameters, v and w vs. Temperature plot
    plot(p.T, p.v, label="v", markersize=3, marker=:rect, xlab="Temperature (K)", ylab="hbar-omega")
    plot!(p.T, p.w, label="w", markersize=3, marker=:diamond)

    savefig("$fileprefix-vw.$extension")

    #####
    ## Spring Constants vs. Temperature plot
    plot(p.T, p.κ, label="Polaron spring-constant", markersize=3, marker=:uptriangle, xlab="Temperature (K)", ylab="Spring-constant",)

    savefig("$fileprefix-spring.$extension")

    #####
    ## Variation Energy vs. Temperature plots
    plot(p.T, p.A, label="A", markersize=3, marker=:downtriangle, xlab="Temperature (K)", ylab="Polaron free-energy")
    plot!(p.T, p.B, label="B", markersize=3, marker=:diamond)
    plot!(p.T, p.C, label="C", markersize=3, marker=:uptriangle)
    plot!(p.T, p.F, label="F", markersize=3, marker=:rect)
    #plot!(Ts,Fs,label="F=-(A+B+C)",markersize=3,marker=:rect)

    savefig("$fileprefix-variational.$extension")

    #####
    ## Polaron radius vs. Temperature
    plot(p.T, p.R .* 10^10, markersize=3, marker=:rect,
        label="Polaron radius", xlab="Temperature (K)", ylab="Polaron Radius (Angstrom)", ylims=(0, Inf))
    #    plot!(p.T,p.rfsmallalpha.*10^10,label="T=0 Schultz small alpha polaron radius") # obsolete
    savefig("$fileprefix-radius.$extension")

    #####
    ## Calculated mobility comparison plot
    plot(p.T, p.μK, label="Kadanoff", markersize=3, marker=:rect, xlab="Temperature (K)", ylab="Mobility (cm\$^2\$/Vs)", ylims=(0, 1000))
    plot!(p.T, p.μFHIP, label="FHIP", markersize=3, marker=:diamond)
    plot!(p.T, p.μH, label="Hellwarth1999", markersize=3, marker=:uptriangle)

    savefig("$fileprefix-mobility-calculated.$extension")
end

end


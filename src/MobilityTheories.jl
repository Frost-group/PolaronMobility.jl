# MobilityTheories.jl

"""
    polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)

Solves the Feynman polaron problem variationally with finite temperature Osaka energies. From the resulting v, and w parameters, calculates polaron structure (wave function size, etc.).  Uses FHIP, Kadanoff (Boltzmann relaxation time) and Hellwarth direct contour integration to predict a temperature-dependent mobility for the material system.

# Arguments  
- `Trange::range`: temperature range.
- `ε_Inf`: reduced optical dielectric constant.
- `ε_S`: reduced static dielectric constant.
- `freq`: characteristic dielectric phonon frequency (THz).
- `effectivemass`: bare-band effective mass (mₑ).

Returns a structure of type `Polaron`, containing arrays of useful
information.  Also prints a lot of information to the standard out - which
may be more useful if you're just inquiring as to a particular data point,
rather than plotting a temperature-dependent parameter.

As an example, to calculate the electron polaron in MAPI at 300 K:
# Examples
```jldoctest
polaronmobility(300, 4.5, 24.1, 2.25, 0.12)
```
"""
function polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)
    println("\n\nPolaron mobility for system ε_Inf=$ε_Inf, ε_S=$ε_S, freq=$freq,
                 effectivemass=$effectivemass; with Trange $Trange ...")

    # Internally we we use 'mb' for the 'band mass' in SI units, of the effecitve-mass of the electron
    mb = effectivemass * MassElectron
    ω = (2 * pi) * freq * 1e12 # angular-frequency, of the phonon mode

    α = frohlichalpha(ε_Inf, ε_S, freq, effectivemass)
    #α=2.395939683378253 # Hard coded; from MAPI params, 4.5, 24.1, 2.25THz, 0.12me

    v = 7.1 # starting guess for v,w variational parameters
    w = 6.5

    @printf("Polaron mobility input parameters: ε_Inf=%f ε_S=%f freq=%g α=%f \n", ε_Inf, ε_S, freq, α)
    @printf("Derived params in SI: ω =%g mb=%g \n", ω, mb)

    # Empty struct for storing data
    # A slightly better way of doing this ф_ф ...
    p = OldPolaron()
    # populate data structure with (athermal) parameters supplied...
    append!(p.α, α) # appending so as not to mess with type immutability
    append!(p.mb, mb)
    append!(p.ω, ω)

    # We define βred as the subsuming the energy of the phonon; i.e. kbT c.f. ħω
    for T in Trange
        β = 1 / (kB * T)
        βred = ħ * ω * β
        append!(p.βred, βred)
        @printf("T: %f β: %.3g βred: %.3g ħω  = %.3g meV\t", T, β, βred, 1E3 * ħ * ω / q)

        if T == 0
            v, w = feynmanvw(v, w, α, 1.0)
        else
            v, w = feynmanvw(v, w, α, 1.0, βred)
        end

        @printf("\n Polaron Parameters:  v= %.4f  w= %.4f", v, w)

        # From 1962 Feynman, definition of v and w in terms of the coupled Mass and spring-constant
        # See Page 1007, just after equation (18)
        # Units of M appear to be 'electron masses'
        # Unsure of units for k, spring coupling constant
        k = (v^2 - w^2)
        M = (v^2 - w^2) / w^2
        append!(p.k, k)
        append!(p.M, M)

        @printf("  ||   M=%f  k=%f\t", M, k)
        @printf("\n Bare-band effective mass: %f Polaron effective mass: %f Polaron mass enhancement: %f%%", effectivemass, effectivemass * (1 + M), M * 100)

        @printf("\n Polaron frequency (SI) v= %.2g Hz  w= %.2g Hz",
            v * ω / (2 * pi), w * ω / (2 * pi))

        # (46) in Feynman1955
        meSmallAlpha(α) = α / 6 + 0.025 * α^2
        # (47) In Feynman1955
        meLargeAlpha(α) = 16 * α^4 / (81 * π^4)
        #meLargeAlpha(α )=202*(α /10)^4
        if (verbose) # asymptotic solutions - not that interesting when you have the actual ones!
            @printf("\n Feynman1955(46,47): meSmallAlpha(α)= %.3f meLargeAlpha(α)= %.3f",
                meSmallAlpha(α), meLargeAlpha(α))
            @printf("\n Feynman1962: Approximate ~ Large alpha limit, v/w = %.2f  =~approx~= alpha^2 = %.2f ",
                v / w, α^2)
        end

        # POLARON SIZE
        @printf("\n Polaron size (rf), following Schultz1959. (s.d. of Gaussian polaron ψ )")
        # Schultz1959 - rather nicely he actually specifies everything down into units!
        # just before (2.4) in Schultz1959
        mu = ((v^2 - w^2) / v^2)
        # (2.4)
        rf = sqrt(3 / (2 * mu * v))
        # (2.4) SI scaling inferred from units in (2.5a) and Table II
        rfsi = rf * sqrt(2 * me * ω)
        @printf("\n\t Schultz1959(2.4): rf= %g (int units) = %g m [SI]", rf, rfsi)
        append!(p.rfsi, rfsi)

        if (verbose)
            scale = sqrt(2 * mb * ω) # Note we're using mb;
            #band effective-mass in SI units (i.e. meff*melectron)

            rfa = (3 / (0.44 * α))^0.5 # As given in Schultz1959(2.5a), but that 0.44 is actually 4/9
            @printf("\n\t Schultz1959(2.5a) with 0.44: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]", rfa, scale * rfa)
            rfa = (3 / ((4 / 9) * α))^0.5 # Rederived from Feynman1955, 8-8-2017; Yellow 2017.B Notebook pp.33-34
            @printf("\n\t Schultz1959(2.5a) with 4/9 re-derivation: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]", rfa, scale * rfa)
            append!(p.rfsmallalpha, scale * rfa)

            rfb = 3 * (pi / 2)^0.5 * α
            @printf("\n\t Schultz1959(2.5b): Feynman α→∞ expansion: rf= %g (int units) = %g m [SI]", rfb, scale * rfb)

            # Schultz1959 - Between (5.7) and (5.8) - resonance of Feynman SHM system
            phononfreq = sqrt(k / M)
            @printf("\n\t Schultz1959(5.7-5.8): fixed-e: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                phononfreq, phononfreq * ω / (2 * pi), phononfreq * hbar * ω * 1000 / q)

            phononfreq = sqrt(k / mu) # reduced mass
            @printf("\n\t Schultz1959(5.7-5.8): reducd mass: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                phononfreq, phononfreq * ω / (2 * pi), phononfreq * hbar * ω * 1000 / q)

            @printf("\n\t Schultz1959: electronfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                sqrt(k / 1), sqrt(k / 1) * ω / (2 * pi), sqrt(k / 1) * hbar * ω * 1000 / q)
            @printf("\n\t Schultz1959: combinedfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
                sqrt(k / (1 + M)), sqrt(k / (1 + M)) * ω / (2 * pi), sqrt(k / (1 + M)) * hbar * ω * 1000 / q)

            # Devreese1972: 10.1103/PhysRevB.5.2367
            # p.2371, RHS.
            @printf("\n Devreese1972: (Large Alpha) Franck-Condon frequency = %.2f", 4 * α^2 / (9 * pi))
        end

        # F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) #(62a) - Hellwarth 1999
        if (T > 0)
            @printf("\n Polaron Free Energy: A= %f B= %f C= %f F= %f", A(v, w, βred), B(v, w, βred, α), C(v, w, βred), F(v, w, βred, α)[1])
            @printf("\t = %f meV", 1000.0 * F(v, w, βred, α)[1] * ħ * ω / q)
            append!(p.A, A(v, w, βred))
            append!(p.B, B(v, w, βred, α))
            append!(p.C, C(v, w, βred))
            append!(p.F, F(v, w, βred, α))
        else # Athermal case; Enthalpy
            @printf("\n Polaron Enthalpy: F= %f = %f meV \n", F(v, w, α, 1.0), 1_000 * F(v, w, α, 1.0) * ħ * ω / q)

            return # return early, as if T=0, all mobility theories = infinite / fall over
        end

        # FHIP
        #    - low-T mobility, final result of Feynman1962
        # [1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv)
        # I believe here β is in SI (expanded) units
        @printf("\n Polaron Mobility theories:")
        μ = (w / v)^3 * (3 * q) / (4 * mb * ħ * ω^2 * α * β) * exp(ħ * ω * β) * exp((v^2 - w^2) / (w^2 * v))
        @printf("\n\tμ(FHIP)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
        append!(p.T, T)
        append!(p.FHIPμ, μ * 100^2)


        # Kadanoff
        #     - low-T mobility, constructed around Boltzmann equation.
        #     - Adds factor of 3/(2*beta) c.f. FHIP, correcting phonon emission behaviour
        # [1.61] in Devreese2016 - Kadanoff's Boltzmann eqn derived mob
        μ = (w / v)^3 * (q) / (2 * mb * ω * α) * exp(ħ * ω * β) * exp((v^2 - w^2) / (w^2 * v))
        @printf("\n\tμ(Kadanoff,via Devreese2016)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)

        append!(p.Kμ, μ * 100^2)

        ######
        # OK, now deep-diving into Kadanoff1963 itself to extract Boltzmann equation components
        # Particularly right-hand-side of page 1367
        #
        # From 1963 Kadanoff, (9), define eqm. number of phonons (just from T and phonon omega)
        Nbar = (exp(βred) - 1)^-1
        @printf("\n\t\tEqm. Phonon. pop. Nbar: %f ", Nbar)
        if (verbose)
            @printf("\n\texp(Bred): %f exp(-Bred): %f exp(Bred)-1: %f", exp(βred), exp(-βred), exp(βred) - 1)
        end
        Nbar = exp(-βred)
        #Note - this is only way to get Kadanoff1963 to be self-consistent with
        #FHIP, and later statements (Devreese) of the Kadanoff mobility.
        #It suggests that Kadanoff used the wrong identy for Nbar in 23(b) for
        #the Gamma0 function, and should have used a version with the -1 to
        #account for Bose / phonon statistics

        #myv=sqrt(k*(1+M)/M) # cross-check maths between different papers
        #@printf("\nv: %f myv: %f\n",v,myv)

        # Between 23 and 24 in Kadanoff 1963, for small momenta skip intergration --> Gamma0
        Gamma0 = 2 * α * Nbar * (M + 1)^(1 / 2) * exp(-M / v)
        Gamma0 *= ω  #* (ω *hbar) # Kadanoff 1963 uses hbar=omega=mb=1 units
        # Factor of omega to get it as a rate relative to phonon frequency
        # Factor of omega*hbar to get it as a rate per energy window
        μ = q / (mb * (M + 1) * Gamma0) #(25) Kadanoff 1963, with SI effective mass

        if (verbose) # these are cross-checks
            @printf("\n\tμ(Kadanoff1963 [Eqn. 25]) = %f m^2/Vs \t = %.2f cm^2/Vs", μ, μ * 100^2)
            @printf("\n\t\t Eqm. Phonon. pop. Nbar: %f ", Nbar)
        end

        @printf("\n\t\tGamma0 = %g rad/s = %g /s ",
            Gamma0, Gamma0 / (2 * pi))
        @printf(" \n\t\tTau=1/Gamma0 = %g s = %f ps",
            2 * pi / Gamma0, 2 * pi * 1E12 / Gamma0)
        Eloss = hbar * ω * Gamma0 / (2 * pi) # Simply Energy * Rate
        @printf("\n\t\tEnergy Loss = %g J/s = %g meV/ps", Eloss, Eloss * 1E3 / (q * 1E12))
        append!(p.Tau, 2 * pi * 1E12 / Gamma0) # Boosted into ps ?


        # Hellwarth1999 - directly do contour integration in Feynman1962, for
        # finite temperature DC mobility
        # Hellwarth1999 Eqn (2) and (1) - These are going back to the general
        # (pre low-T limit) formulas in Feynman1962.  to evaluate these, you
        # need to do the explicit contour integration to get the polaron
        # self-energy
        R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)

        #b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
        # OK! I now understand that there is a typo in Hellwarth1999 and
        # Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
        # the original, Feynman1962:
        b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)

        a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
        k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
        K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
        μ = RHS^-1 * (q) / (ω * mb)
        @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
        append!(p.Hμ, μ * 100^2)

        if (verbose)
            # Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
            # So let's test this
            R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
            b = 0
            a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
            #k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
            K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

            #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
            RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
            μ = RHS^-1 * (q) / (ω * mb)
            @printf("\n\tμ(Hellwarth1999,b=0)= %f m^2/Vs \t= %.2f cm^2/Vs", μ, μ * 100^2)
            @printf("\n\tError due to b=0; %f", (100^2 * μ - p.Hμ[length(p.Hμ)]) / (100^2 * μ))
            #append!(Hμs,μ*100^2)
        end
        @printf("\n") # blank line at end of spiel.

        # Recycle previous variation results (v,w) as next guess
        initial = [v, w] # Caution! Might cause weird sticking in local minima

        append!(p.v, v)
        append!(p.w, w)

    end

    return (p)
end

function polaron(αrange, Trange, Ωrange, ω, v_guesses, w_guesses; verbose = false)

    num_α = size(αrange, 1)
    num_T = length(Trange)
    num_Ω = length(Ωrange)
    num_ω = length(ω)

    function reduce_array(a) 
        if length(a) == 1
            only(a)
        else 
            dropdims(a, dims = tuple(findall(size(a) .== 1)...))
        end
    end

    @assert length(v_guesses) == length(w_guesses) "v and w guesses must be the same length."
    num_vw = length(v_guesses)

    # Instantiate 
    p = Dict(
        "α" => αrange,                                           # alphas
        "αeff" => sum(αrange, dims=2),                           # alphas sums
        "T" => Trange,                                           # temperatures
        "ω" => ω,                                                # phonon frequencies
        "β" => Matrix{Float64}(undef, num_T, num_ω),             # betas
        "Ω" => Ωrange,                                           # photon frequencies
        "v0" => Matrix{Float64}(undef, num_α, num_vw),           # v ground state params
        "w0" => Matrix{Float64}(undef, num_α, num_vw),           # w ground state params
        "F0" => Vector{Float64}(undef, num_α),                   # ground state energies
        "A0" => Vector{Float64}(undef, num_α),                   # A ground state parameter
        "B0" => Vector{Float64}(undef, num_α),                   # B ground state parameter
        "C0" => Vector{Float64}(undef, num_α),                   # C ground state parameter
        "v" => Array{Float64, 3}(undef, num_T, num_α, num_vw),   # v params
        "w" => Array{Float64, 3}(undef, num_T, num_α, num_vw),   # w params
        "F" => Matrix{Float64}(undef, num_T, num_α),             # energies
        "A" => Matrix{Float64}(undef, num_T, num_α),             # A parameter
        "B" => Matrix{Float64}(undef, num_T, num_α),             # B parameter
        "C" => Matrix{Float64}(undef, num_T, num_α),             # C parameter
        "κ" => Array{Float64, 3}(undef, num_T, num_α, num_vw),   # spring constants
        "M" => Array{Float64, 3}(undef, num_T, num_α, num_vw),   # fictitious masses
        "R" => Array{Float64, 3}(undef, num_T, num_α, num_vw),   # polaron radii
        "z" => Array{ComplexF64, 3}(undef, num_Ω, num_T, num_α), # complex impedences
        "σ" => Array{ComplexF64, 3}(undef, num_Ω, num_T, num_α), # complex conductivities
        "μ" => Matrix{Float64}(undef, num_T, num_α)              # mobilities
    )

    if verbose
        println("\e[?25l\e[K-----------------------------------------------------------------------")
        println("\e[K                         Polaron Information:                          ")
        println("\e[K-----------------------------------------------------------------------")
        if num_ω == 1
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KPhonon frequencies         | ω = ", ω, " ω₀")
        else
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KPhonon frequencies         | ω = ", join(round.(first(ω, 2), digits = 1), ", ")..., " ... ", join(round.(last(ω, 2), digits = 1), ", ")..., " ω₀")
        end
        process = 1
        αprocess = 1
    end

    for j in axes(αrange, 1)

        α = reduce_array(αrange[j, :])

        if verbose
            if num_ω == 1
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KFröhlich coupling          | α = ", join(α, ", ")...)
            else
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KFröhlich coupling          | α = ", join(round.(first(α, 2), digits = 3), ", ")..., " ... ", join(round.(last(α, 2), digits = 3), ", ")...)
            end
        end

        v0, w0, F0, A0, B0, C0 = feynmanvw(v_guesses, w_guesses, α, ω)
        v_guesses, w_guesses = v0, w0
        p["v0"][j, :] .= v0
        p["w0"][j, :] .= w0
        p["F0"][j] = F0
        p["A0"][j] = A0
        p["B0"][j] = B0
        p["C0"][j] = C0

        if verbose
            println("\e[K-----------------------------------------------------------------------") 
            println("\e[K              Ground State Information: [$(αprocess[]) / $(num_α) ($(round(αprocess[] / (num_α) * 100, digits=1)) %)]")
            println("\e[K-----------------------------------------------------------------------")
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KGS variational parameter   | v₀ = ", v0, " ω₀")
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KGS variational parameter   | w₀ = ", w0, " ω₀")
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KGS Energy                  | E₀ = ", F0 ./ 2π, " ħω₀")
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KGS Electron energy         | A₀ = ", A0, " ħω₀")
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KGS Interaction energy      | B₀ = ", B0, " ħω₀")
            println(IOContext(stdout, :compact => true, :limit => true), "\e[KGS Trial energy            | C₀ = ", C0, " ħω₀")
            Tprocess = 1
            αprocess += 1
        end

        for i in eachindex(Trange)
            T = Trange[i]

            if verbose
                println("\e[K-----------------------------------------------------------------------") 
                println("\e[K         Finite Temperature Information: [$(Tprocess[]) / $(num_T) ($(round(Tprocess[] / (num_T) * 100, digits=1)) %)]")
                println("\e[K-----------------------------------------------------------------------")
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KTemperatures               | T = ", T, " K")
            end

            # Calculate reduced thermodynamic betas for each phonon mode.
            # If temperature is zero, return Inf.
            β = ω ./ T
            p["β"][i, :] .= β  

            if verbose
                if num_ω == 1
                    println(IOContext(stdout, :compact => true, :limit => true), "\e[KReduced thermodynamic      | β = ", β)
                else
                    println(IOContext(stdout, :compact => true, :limit => true), "\e[KReduced thermodynamic      | β = ", join(round.(first(β, 2), digits = 3), ", ")..., " ... ", join(round.(last(β, 2), digits = 3), ", ")...)
                end
            end

            # Calculate variational parameters for each alpha parameter and temperature. Returns a Matrix of tuples.
            v, w, F, A, B, C = feynmanvw(v_guesses, w_guesses, α, ω, β)
            v_guesses, w_guesses = v, w
            p["v"][i, j, :] .= v
            p["w"][i, j, :] .= w
            p["F"][i, j] = F
            p["A"][i, j] = A
            p["B"][i, j] = B
            p["C"][i, j] = C

            if verbose
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KVariational parameter      | v = ", v, " ω₀")
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KVariational parameter      | w = ", w, " ω₀")
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KFree energy                | F = ", F ./ 2π, " ħω₀")
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KElectron energy            | A = ", A ./ 2π, " ħω₀")
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KInteraction energy         | B = ", B ./ 2π, " ħω₀")
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KTrial energy               | C = ", C ./ 2π, " ħω₀")
            end

            # Calculate fictitious spring constants for each alpha parameter and temperature. Returns a Matrix.
            κ = v .^2 .- w .^2
            p["κ"][i, j, :] .= κ

            if verbose
                println("\e[K-----------------------------------------------------------------------") 
                println("\e[K                      Trial System Information:                        ")
                println("\e[K-----------------------------------------------------------------------") 
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KFictitious spring constant | κ = ", κ, " m₀/ω₀²")
            end

            # Calculate fictitious masses for each alpha parameter and temperature. Returns a Matrix.
            M = κ ./ w .^2
            p["M"][i, j, :] .= M

            if verbose
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KFictitious mass            | M = ", M, " m₀")
            end

            R = sqrt.(3 .* v ./ (v .^2 .- w .^2) .^2)
            p["R"][i, j, :] .= R

            if verbose
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KPolaron radius             | R = ", R, " √(ħ/2m₀ω₀)")
            end

            # Calculates the dc mobility for each alpha parameter and each temperature.
            μ = polaron_mobility(v, w, α, ω, β)
            p["μ"][i, j] = μ
            
            if verbose
                println(IOContext(stdout, :compact => true, :limit => true), "\e[KMobility                   | μ = ", μ, " q/m₀ω₀")
                Ωprocess = 1
                Tprocess += 1
            end

            for k in eachindex(Ωrange)
                Ω = Ωrange[k]

                if verbose
                    println("\e[K-----------------------------------------------------------------------") 
                    println("\e[K           Linear Reponse Information: [$(Ωprocess[]) / $(num_Ω) ($(round(Ωprocess[] / (num_Ω) * 100, digits=1)) %)]")
                    println("\e[K-----------------------------------------------------------------------")
                    println(IOContext(stdout, :compact => true, :limit => true), "\e[KElectric field frequency   | Ω = ", Ω, " ω₀")
                end

                χ = polaron_memory_function(v, w, α, ω, β, Ω)
                if verbose
                    println(IOContext(stdout, :compact => true, :limit => true), "\e[KMemory function            | χ = ", χ .|> y -> ComplexF64.(y), " ω₀m₀V₀/q²")
                end

                # Calculate complex impedances for each alpha parameter, frequency and temperature. Returns a 3D Array.
                z = -im * Ω + im * χ
                p["z"][k, i, j] = z

                if verbose
                    println(IOContext(stdout, :compact => true, :limit => true), "\e[KComplex impedance          | z = ", z .|> y -> ComplexF64.(y), " ω₀m₀V₀/q²")
                end

                # Calculate complex conductivities for each alpha parameter, frequency and temperature. Returns a 3D array.
                σ = 1 / z
                p["σ"][k, i, j] = σ

                if verbose
                    println(IOContext(stdout, :compact => true, :limit => true), "\e[KComplex conductivity       | σ = ", σ .|> y -> ComplexF64.(y), " q²/ω₀m₀V₀")
                    println("\e[K-----------------------------------------------------------------------") 
                    println("\e[K[Total Progress: $(process[]) / $(num_α * num_T * num_Ω) ($(round(process[] / (num_α * num_T * num_Ω) * 100, digits=1)) %)]")
                    print("\e[9F")
                    Ωprocess += 1
                    process += 1
                end
            end
            if verbose print("\e[18F") end
        end
        if verbose 
            print("\e[10F") 
        end
    end
    if verbose print("\e[5F\e[?25h") end

    polaron_data = [
        p["α"],     # alphas
        p["αeff"],  # alphas sums
        p["T"],     # temperatures
        p["ω"],     # phonon frequencies
        p["β"],     # betas
        p["Ω"],     # photon frequencies
        p["v0"],    # v params ground state
        p["w0"],    # w params ground state
        p["F0"],    # energies ground state
        p["A0"],    # A parameter, bare electron energy ground state
        p["B0"],    # B parameter, interaction energy ground state
        p["C0"],    # C parameter, trial system free energy ground state
        p["v"],     # v params
        p["w"],     # w params
        p["F"],     # energies
        p["A"],     # A parameter, bare electron energy
        p["B"],     # B parameter, interaction energy
        p["C"],     # C parameter, trial system free energy
        p["κ"],     # spring constants
        p["M"],     # fictitious masses
        p["R"],     # polaron radii
        p["z"],     # complex impedences
        p["σ"],     # complex conductivities
        p["μ"]      # mobilities
    ]
    return Polaron(polaron_data...)
end

# Default guesses for variational paramaters.
polaron(αrange, Trange, Ωrange, ω; verbose = false) = polaron(αrange, Trange, Ωrange, ω, 3.11, 2.87; verbose = verbose)

# Single alpha parameter.
polaron(α::Real, Trange, Ωrange, ω; verbose = false) = polaron([α], Trange, Ωrange, ω; verbose = verbose)

# DC limit, zero frequency.
polaron(αrange, Trange, ω; verbose = false) = polaron(αrange, Trange, 0, ω; verbose = verbose)

polaron(αrange, Trange; verbose = false) = polaron(αrange, Trange, 1; verbose = verbose)

polaron(αrange; verbose = false) = polaron(αrange, 298, 1; verbose = verbose)

# Material specific
function polaron(material::Material, Trange, Ωrange; verbose = false)

    if verbose
        display(material)
    end

    ω₀ = 2π * 1e12 
    m₀ = me * 0.12
    r₀ = sqrt(ħ / (m₀ * ω₀))
    T₀ = ħ * ω₀ / kB

    m_eff = material.mb
    volume = material.volume
    phonon_freqs = material.freqs

    p = polaron(material.α', Trange ./ T₀, Ωrange, phonon_freqs, verbose = verbose)

    # Add Units
    F0_unit = p.F0 .* ħ * ω₀ / q * 1e3                                               # meV
    A0_unit = p.A0 .* ħ * ω₀ / q * 1e3                                             # meV
    B0_unit = p.B0 .* ħ * ω₀ / q * 1e3                                           # meV
    C0_unit = p.C0 .* ħ * ω₀ / q * 1e3                                       # meV
    F_unit = p.F .* ħ * ω₀ / q * 1e3                                         # meV
    A_unit = p.A .* ħ * ω₀ / q * 1e3                                         # meV
    B_unit = p.B .* ħ * ω₀ / q * 1e3                                        # meV
    C_unit = p.C .* ħ * ω₀ / q * 1e3                                            # meV
    M_units = p.M .* m_eff                                                           # electron mass units
    R_unit = [R * sqrt(ħ / 2 / m_eff / me / ω / 1e12) * 1e10 for R in p.R, ω in p.ω] # Angstroms
    z_unit = p.z .* (m₀ * r₀^2 * ω₀ / q^2)                    # Ohms
    σ_unit = p.σ ./ (m₀ * r₀^2 * ω₀ / q^2)                  # Siemens
    μ_unit = p.μ ./ m₀ / ω₀ * q * 1e4                # cm^2/Vs

    return Polaron(p.α, p.αeff, p.T, phonon_freqs, p.β, p.Ω, p.v0, p.w0, F0_unit, A0_unit, B0_unit, C0_unit, p.v, p.w, F_unit, A_unit, B_unit, C_unit, p.κ, M_units, R_unit, z_unit, σ_unit, μ_unit)
end

# DC limit for material specific code, zero frequency.
polaron(material::Material, Trange; verbose = false) = polaron(material::Material, Trange, 0; verbose = verbose)


"""
    save_polaron(p::NewPolaron, prefix)

Saves data from 'polaron' into file "prefix".
This is a .jdl file for storing the polaron data whilst preserving types. Allows for saving multidimensional arrays that sometimes arise in the polaron data.
Each parameter in the NewPolaron type is saved as a dictionary entry. E.g. NewPolaron.α is saved under JLD.load("prefix.jld")["alpha"].
"""
function save_polaron(polaron::Polaron, prefix)

    println("Saving polaron data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "alpha", polaron.α,
        "temperature", polaron.T,
        "beta", polaron.β,
        "phonon freq", polaron.ω,
        "v", polaron.v,
        "w", polaron.w,
        "spring", polaron.κ,
        "mass", polaron.M,
        "energy", polaron.F,
        "efield freq", polaron.Ω,
        "impedance", polaron.Z,
        "conductivity", polaron.σ,
        "mobility", polaron.μ
    )

    println("... Polaron data saved.")
end

"""
    load_polaron(p::NewPolaron, prefix)

Loads data from file "polaron_file_path" into a NewPolaron type.
"""
function load_polaron(polaron_file_path)

    println("Loading polaron data from $polaron_file_path ...")

    data = JLD.load("$polaron_file_path")

    polaron = NewPolaron(
        data["alpha"],
        data["temperature"],
        data["beta"],
        data["phonon freq"],
        data["v"],
        data["w"],
        data["spring"],
        data["mass"],
        data["energy"],
        data["efield freq"],
        data["impedance"],
        data["conductivity"],
        data["mobility"]
    )

    println("... Polaron loaded.")

    return polaron
end

"""
    savepolaron(fileprefix, p::Polaron)

Saves data from polaron 'p' into file "fileprefix".
This is a simple space-delimited text file, with each entry a separate temperature, for plotting with Gnuplot or similar.

Structure of file is written to the header:
# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis
# 1    2     3    4     5      6   7   8   9  10  11  12  13    14     15
"""
function savepolaron(fileprefix, p::Polaron)
    println("Saving data to $fileprefix.dat ...")
    f = open("$fileprefix.dat", "w")

    @printf(f, "# %s\n", fileprefix) # put name / material at header
    @printf(f, "# Params in SI: ω =%g mb=%g \n", p.ω[1], p.mb[1])
    @printf(f, "# Alpha parameter: α = %f  \n", p.α[1])

    @printf(f, "# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis\n")
    @printf(f, "#  1    2     3    4     5      6   7   8   9  10  11  12  13    14     15\n") # columns for GNUPLOT etc.

    for i in 1:length(p.T)
        @printf(f, "%d %03f %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
            p.T[i], p.βred[i], p.Kμ[i], p.Hμ[i], p.FHIPμ[i],
            p.v[i], p.w[i],
            p.k[i], p.M[i], p.A[i], p.B[i], p.C[i], p.F[i],
            p.Tau[i], p.rfsi[i])
    end
    close(f)
end

"""
    Hellwarth1999mobilityRHS((α, (v, w) ,f), effectivemass, T)

Calculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).

See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
"""
function Hellwarth1999mobilityRHS((α, (v, w), f), effectivemass, T)
    mb = effectivemass * MassElectron
    ω = f * 1e12 * 2π
    βred = ħ * ω / (kB * T)

    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
    b = R * βred / sinh(βred * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
    a = sqrt((βred / 2)^2 + R * βred * coth(βred * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

    # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * βred^(5 / 2) / sinh(βred / 2) * (v^3 / w^3) * K
    μ = RHS^(-1) * q / (ω * mb)

    return 1 / μ
end

function Hellwarth_mobility(β, α, v, w; ω = ω)
    R = (v^2 - w^2) / (w^2 * v) # inline, page 300 just after Eqn (2)
    b = R * β / sinh(β * v / 2) # Feynman1962 version; page 1010, Eqn (47b)
    a = sqrt((β / 2)^2 + R * β * coth(β * v / 2))
    k(u, a, b, v) = (u^2 + a^2 - b * cos(v * u))^(-3 / 2) * cos(u) # integrand in (2)
    K = quadgk(u -> k(u, a, b, v), 0, Inf)[1] # numerical quadrature integration of (2)

    # Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
    RHS = α / (3 * sqrt(π)) * β^(5 / 2) / sinh(β / 2) * (v^3 / w^3) * K
    μ = RHS^(-1)
    return μ
end

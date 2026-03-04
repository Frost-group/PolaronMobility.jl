# FrohlichPolaron.jl
# A do-it-all function for calculating properties of the Frohlich polaron. Data is stored as a struct.

"""
    Polaron(x...)

Type for storing the polaron information, `x...`.
"""
mutable struct FrohlichPolaron
    α       # Fröhlich coupling.
    αeff    # Effective Fröhlich coupling summed for multiple modes.
    mb      # Band mass
    T       # Temperature.
    ω       # Phonon mode frequencies.
    ωeff    # Effective phonon mode frequency.
    β       # Reduced thermodynamic beta ħω₀/kBT.
    Ω       # Electric field frequency.
    d       # Number of spatial dimensions.
    v0      # Athermal variational parameter v.
    w0      # Athermal variational parameter w.
    F0      # Polaron athermal energy (enthalpy).
    A0      # Bare electron enthalpy.
    B0      # ⟨S⟩ₜ interaction enthalpy.
    C0      # ⟨Sₜ⟩ₜ enthalpy of trial system.
    Fs      # Small alpha (α→0) approximate energy.
    Fl      # Large alpha (α→∞) approximate energy.
    κ0      # Fictitious spring constant.
    M0      # Athermal fictitious mass.
    Ms      # Small alpha (α→0) approximate fictitious mass.
    Ml      # Large alpha (α→∞) approximate fictitious mass.
    M0a     # Athermal asymptotic approximate fictitious mass (v0/w0).
    M0r     # Athermal reduced mass of particle + fictitious particle system.
    R0      # Athermal polaron radius (s. d. of a Gaussian wavefunction).
    Rs      # Small alpha (α→0) approximate polaron radius.
    Rl      # Large alpha (α→∞) approximate polaron radius.
    ΩFC     # Approximate Franck-Condon peak frequency.
    v       # Thermal variational parameter v.
    w       # Thermal variational parameter w.
    F       # Polaron free energy.
    A       # Bare electron free energy.
    B       # ⟨S⟩ₜ interaction energy.
    C       # ⟨Sₜ⟩ₜ free energy of trial system.
    κ       # Fictitious spring constant.
    M       # Thermal fictitious mass.
    Ma      # Thermal asymptotic approximate fictitious mass (v/w).
    Mr      # Thermal reduced mass of particle + fictitious particle system.
    R       # Schultz polaron radius.
    μ       # DC mobility.
    μFHIP   # FHIP low-temperature DC mobiltiy approximation.
    μD      # Kadanoff DC mobility, low-temperature approximation using Devreese2016.
    μK      # Kadanoff DC mobility, low-temperature approximation using Kadanoff1963.
    τ       # Relaxation time from Kadanoff Boltzmann transport equation.
    μH      # Hellwarth DC mobility.
    μH0     # Hellwarth DC mobility with b=0.
    χ       # Memory function or polaron self energy.
    z       # Complex impedence.
    σ       # Complex conductivity.

    function FrohlichPolaron(x...)
        new(reduce_array.(x)...)
    end
end

"""
    polaron(αrange, Trange, Ωrange, ω, v_guesses, w_guesses; verbose=false)

Outer constructor for the Polaron type. This function evaluates model data for the polaron, i.e. unitless and not material specific. 

# Examples
```jldoctest
julia> polaron(6, 300, 3, 1.0, 3.6, 2.8)
```
"""
function frohlichpolaron(αrange, Trange, Ωrange; ω=1, ωeff=1, mb=1, v_guesses=false, w_guesses=false, dims=3, kspace=false, cycleguesses=false, verbose=false)

    # v_guesses and w_guesses are initial values for v and w (including many v and w parameters).
    # These guesses are generally not needed unless instabilities are found in the minimisation and better initial values improve stability.

    # Get the length of any arrays etc.
    num_α = size(αrange, 1)
    num_T = length(Trange)
    num_Ω = length(Ωrange)
    num_ω = length(ω)
    num_d = length(dims)

    ω = reduce_array(ω)
    
    # For multiple variational modes, ensure that the number of v and w parameters is the same.
    @assert length(v_guesses) == length(w_guesses) "v and w guesses must be the same length."
    num_vw = length(v_guesses)

    # Instantiate all the polaron data that will go into the Polaron type.
    p = Dict(
        "α"     => αrange, # Alphas.
        "αeff"  => sum(αrange, dims=2), # Alphas sums.
        "mb"    => mb, # Band mass.
        "T"     => Trange, # Temperatures.
        "ω"     => ω, # Phonon frequencies.
        "ωeff"  => ωeff, # Effective Phonon frequency.
        "β"     => Vector{Float64}(undef, num_T), # Betas.
        "Ω"     => Ωrange, # Photon frequencies.
        "d"     => dims, # Number of spatial dimensions.
        "v0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal v params.
        "w0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal w params.
        "F0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal energies.
        "A0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal A parameter.
        "B0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal B parameter.
        "C0"    => Matrix{Float64}(undef, num_d, num_α), # Athermal C parameter.
        "Fs"    => Vector{Float64}(undef, num_α), # Small alpha (α→0) approximate energy.
        "Fl"    => Vector{Float64}(undef, num_α), # Large alpha (α→∞) approximate energy.
        "κ0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Fictitious spring constant.
        "M0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal fictitious mass.
        "Ms"    => Vector{Float64}(undef, num_α), # Small alpha (α→0) approximate fictitious mass.
        "Ml"    => Vector{Float64}(undef, num_α), # Large alpha (α→∞) approximate fictitious mass.
        "M0a"   => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal asymptotic approximate fictitious mass (v0/w0).
        "M0r"   => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal reduced mass of particle + fictitious particle system.
        "R0"    => Array{Float64,3}(undef, num_d, num_α, num_vw), # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        "Rs"    => Vector{Float64}(undef, num_α), # Small alpha (α→0) approximate polaron radius.
        "Rl"    => Vector{Float64}(undef, num_α), # Large alpha (α→∞) approximate polaron radius.
        "ΩFC"   => Vector{Float64}(undef, num_α), # Large alpha (α → ∞)approximate Franck-Condon peak frequency.
        "v"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # v params.
        "w"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # w params.
        "F"     => Array{Float64,3}(undef, num_T, num_d, num_α), # Energies.
        "A"     => Array{Float64,3}(undef, num_T, num_d, num_α), # A parameter.
        "B"     => Array{Float64,3}(undef, num_T, num_d, num_α), # B parameter.
        "C"     => Array{Float64,3}(undef, num_T, num_d, num_α), # C parameter.
        "κ"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # Spring constants.
        "M"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # Fictitious masses.
        "Ma"    => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw),  # Thermal asymptotic approximate fictitious mass (v/w).
        "Mr"    => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw),  # Thermal reduced mass of particle + fictitious particle system.
        "R"     => Array{Float64,4}(undef, num_T, num_d, num_α, num_vw), # Polaron radii.
        "μ"     => Array{Float64,3}(undef, num_T, num_d, num_α), # Mobilities.
        "μFHIP" => Matrix{Float64}(undef, num_T, num_α), # FHIP low-temperature DC mobiltiy approximation.
        "μD"    => Matrix{Float64}(undef, num_T, num_α), # Kadanoff DC mobility, low-temperature approximation using Devreese2016.
        "μK"    => Matrix{Float64}(undef, num_T, num_α), # Kadanoff DC mobility, low-temperature approximation.
        "τ"     => Matrix{Float64}(undef, num_T, num_α), # relaxation time from Kadanoff Boltzmann transport equation.
        "μH"    => Matrix{Float64}(undef, num_T, num_α), # Hellwarth DC mobility.
        "μH0"   => Matrix{Float64}(undef, num_T, num_α), # Hellwarth DC mobility with b=0.
        "χ"     => Array{ComplexF64,4}(undef, num_Ω, num_T, num_d, num_α), # 'Memory function' or self energy.
        "z"     => Array{ComplexF64,4}(undef, num_Ω, num_T, num_d, num_α), # Complex impedences.
        "σ"     => Array{ComplexF64,4}(undef, num_Ω, num_T, num_d, num_α), # Complex conductivities.
    )

    # Print the phonon frequencies. 
    if verbose
        # Instantiate IO. Compact and limit reduce printed info.
        io = IOContext(stdout, :compact => true, :limit => true) 
        process = 1     # Counter for total iterations.
        αprocess = 1    # Counter for αrange iterations. 
    end

    for j in axes(αrange, 1)    # αrange loop.

        dprocess = 1

        # Reduce αrange to either a Vector of values per phonon mode, or a single scalar. 
        α = reduce_array(αrange[j, :])

        # Effective alpha parameter is the sum over all contributions per phonon mode.
        αeff = sum(α)

        if verbose
            # Hide cursor whilst printing live data.
            print(io, "\e[?25l") 

            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K               Polaron Information: [$(αprocess[]) / $(num_α) ($(round(αprocess[] / (num_α) * 100, digits=1)) %)]")
            println(io, "\e[K-----------------------------------------------------------------------")

            # Different formatting for single vs multiple frequencies (limits an array to head and tail to limit prints).
            if num_ω == 1
                println(io, "\e[KPhonon frequencies             | ωeff = ", ωeff, " | ω = ", ω)
                println(io, "\e[KFröhlich coupling              | αeff = ", αeff, " | α = ", join(α, ", ")...)
            else
                println(io, "\e[KPhonon frequencies             | ωeff = ", ωeff, " | ω = ", join(round.(first(ω, 2), digits=2), ", ")..., " ... ", join(round.(last(ω, 2), digits=2), ", ")...)
                println(io, "\e[KFröhlich coupling              | αeff = ", αeff, " | α = ", join(round.(first(α, 2), digits=3), ", ")..., " ... ", join(round.(last(α, 2), digits=3), ", ")...)
            end

            αprocess += 1   # Increment αrange iteration.
        end

        # Small alpha (α → 0) approximate energy.
        F_small = (-αeff - αeff^2 / 81)
        p["Fs"][j] = F_small

        # Print small alpha energy.
        if verbose
            println(io, "\e[KSmall α→0 energy               | Fs = ", F_small)
        end

        # Large alpha (α → ∞) approximate energy.
        F_large = (-αeff^2 / 3π - 3 * log(2) - 3 / 4)
        p["Fl"][j] = F_large

        # Print large alpha energy.
        if verbose
            println(io, "\e[KLarge α→∞ energy               | Fl = ", F_large)
        end

        # Small coupling (α → 0) polaron mass approximation. Eqn. (46) in Feynman1955.
        M_small = (αeff / 6 + 0.025 * αeff^2)
        p["Ms"][j] = M_small

        # Print small alpha fictitious mass.
        if verbose
            println(io, "\e[KSmall α→0 fictitious mass      | Ms = ", M_small)
        end

        # Large coupling (α → ∞) polaron mass approximation. Eqn. (47) In Feynman1955.
        M_large = 16 * αeff^4 / (81 * π^4)
        p["Ml"][j] = M_large

        # Print large alpha fictitious mass.
        if verbose
            println(io, "\e[KLarge α→∞ fictitious mass      | Ml = ", M_large)
        end

        # Small coupling (α → 0) polaron radii approximiation. Eqn. (2.5a) in Schultz1959.
        R_small = sqrt(3 / (4 / 9 * αeff))
        p["Rs"][j] = R_small

        # Print small alpha polaron radius.
        if verbose
            println(io, "\e[KSmall α→0 polaron radius       | Rs = ", R_small)
        end

        # Large coupling (α → ∞) polaron radii approximiation. Eqn. (2.5b) in Schultz1959.
        R_large = 3 * √(π / 2) * αeff
        p["Rl"][j] = R_large

        # Print large alpha polaron radius.
        if verbose
            println(io, "\e[KLarge α→∞ polaron radius       | Rl = ", R_large)
        end

        # Franck-Condon (FC) frequency in large coupling (α → ∞) limit. RHS of pg. 2371 in Devreese1972.
        Ω_FC = 4 / 9π * αeff ^2
        p["ΩFC"][j] = Ω_FC

        # Print large alpha Franck-Condon peak frequency.
        if verbose
            println(io, "\e[KLarge α→∞ FC peak freq.        | ΩFC = ", Ω_FC)
        end

        for d in 1:num_d
            if verbose
                println(io, "\e[KNumber of dimensions", " [", dprocess, " / ", num_d, "]", "   | d = ", dims[d])
                dprocess += 1
            end

        # Extract the ground-state, athermal polaron properties (energy (enthalpy) and variational parameters v and w).
        # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
        # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.


        if !cycleguesses && v_guesses == false
            w_guesses = 2 .+ tanh.((6 .- α) ./ 3)
            v_guesses = αeff < 7 ? 3 .+ α ./ 4 : 4 .* α .^2 / 9π .- 3/2 * (2 * log(2) + 0.5772) .- 3/4
        end
        
        athermal_energy(v, w) = !kspace ? frohlich_energy(v, w, α, ω; dims = dims[d], mb = mb) : frohlich_energy_k_space(v, w, α, ω; dims = dims[d])
        
        v_gs, w_gs, F_gs, A_gs, B_gs, C_gs = vw_variation(athermal_energy, v_guesses, w_guesses)

        # Update the guesses to keep them close-ish to the true solutions during loops over alphas.
        if cycleguesses && (v_guesses == false)
            v_guesses, w_guesses = v_gs, w_gs
        end

        # Store the athermal data.
        p["v0"][d, j, :] .= v_gs 
        p["w0"][d, j, :] .= w_gs
        p["F0"][d, j] = F_gs
        p["A0"][d, j] = A_gs
        p["B0"][d, j] = B_gs
        p["C0"][d, j] = C_gs

        # Print athermal variational parameter and energy data.
        if verbose
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[K                   Zero Temperature Information:                       ")
            println(io, "\e[K-----------------------------------------------------------------------")
            println(io, "\e[KVariational parameter          | v0 = ", v_gs)
            println(io, "\e[KVariational parameter          | w0 = ", w_gs)
            println(io, "\e[KEnergy                         | E0 = ", F_gs)
            println(io, "\e[KElectron energy                | A0 = ", A_gs)
            println(io, "\e[KInteraction energy             | B0 = ", B_gs)
            println(io, "\e[KTrial energy                   | C0 = ", C_gs)
            Tprocess = 1    # Counter for Trange.
        end

        # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
        κ_gs = (v_gs .^ 2 .- w_gs .^ 2)
        p["κ0"][d, j, :] .= κ_gs

        # Print athermal fictitious spring constant.
        if verbose
            println(io, "\e[KFictitious spring constant     | κ0 = ", κ_gs)
        end

        # Calculate and store fictitious masses. Athermal.
        M_gs = κ_gs ./ w_gs .^ 2
        p["M0"][d, j, :] .= M_gs

        # Print athermal fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass                | M0 = ", M_gs)
        end

        # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Athermal
        M_asymp_gs = v_gs ./ w_gs
        p["M0a"][d, j, :] .= M_asymp_gs

        # Print athermal asymptotic fictitious mass.
        if verbose
            println(io, "\e[KFictitious mass (asymptotic)   | M0a = ", M_asymp_gs)
        end

        # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
        M_reduced_gs = (v_gs .^2 .- w_gs .^2) ./ v_gs .^ 2
        p["M0r"][d, j, :] .= M_reduced_gs

        # Print athermal reduced mass.
        if verbose
            println(io, "\e[KReduced mass                   | M0r = ", M_reduced_gs)
        end

        # Calculate and store polaron radii. Approximates the polaron wavefunction as a Gaussian and relates the size to the standard deviation. Eqn. (2.4) in Schultz1959. Athermal.
        R_gs = sqrt.(3 ./ (2 .* M_reduced_gs .* v_gs))
        p["R0"][d, j, :] .= R_gs ./ sqrt(ωeff)

        # Print athermal polaron radius.
        if verbose
            println(io, "\e[KPolaron radius                 | R0 = ", R_gs)
        end

        for i in eachindex(Trange)  # Temperatures loop.
            T = Trange[i]

            if !iszero(T)

                # Print temperature.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K         Finite Temperature Information: [$(Tprocess[]) / $(num_T) ($(round(Tprocess[] / (num_T) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KTemperatures                   | T = ", T)
                end

                # Calculate the reduced (unitless) thermodynamic betas for each phonon mode.
                β = 1 ./ T
                p["β"][i] = β

                # Print thermodynamic betas.
                if verbose
                    println(io, "\e[KReduced thermodynamic          | β = ", β)
                end

                # Calculate thermal polaron properties (energy (Gibbs free energy) and variational parameters v and w).
                # w is also the frequency of oscillation of the SHM trial system composed of the bare particle and fictitous mass.
                # A, B, C are components of the total energy: A is the bare electron energy, B the electron-phonon interaction energy, C is the energy of the harmonic trial system.

                if !cycleguesses && (v_guesses == false)
                    w_guesses = w_gs
                    v_guesses = v_gs
                end

                thermal_energy(v, w) = !kspace ? frohlich_energy(v, w, α, ω, β; dims = dims[d], mb = mb) : frohlich_energy_k_space(v, w, α, ω, β; dims = dims[d])
                v, w, F, A, B, C = vw_variation(thermal_energy, v_guesses, w_guesses)

                # Update the guesses to keep them close-ish to the true solutions during loops over temperatures.
                if cycleguesses && (v_guesses == false)
                    println("true")
                    v_guesses, w_guesses = v, w
                end

                # Store thermal data.
                p["v"][i, d, j, :] .= v
                p["w"][i, d, j, :] .= w
                p["F"][i, d, j] = F
                p["A"][i, d, j] = A
                p["B"][i, d, j] = B
                p["C"][i, d, j] = C

                # Print thermal data.
                if verbose
                    println(io, "\e[KVariational parameter          | v = ", v)
                    println(io, "\e[KVariational parameter          | w = ", w)
                    println(io, "\e[KFree energy                    | F = ", F)
                    println(io, "\e[KElectron energy                | A = ", A)
                    println(io, "\e[KInteraction energy             | B = ", B)
                    println(io, "\e[KTrial energy                   | C = ", C)
                end

                # Calculate and store fictitious spring constants. See after Eqn. (18), pg. 1007 of Feynman1962. Thermal
                κ = (v .^ 2 .- w .^ 2)
                p["κ"][i, d, j, :] .= κ

                # Print spring constants.
                if verbose
                    println(io, "\e[KFictitious spring constant     | κ = ", κ)
                end

                # Calculate and store fictitious masses. Thermal.
                M = κ ./ w .^ 2
                p["M"][i, d, j, :] .= M

                # Print masses.
                if verbose
                    println(io, "\e[KFictitious mass                | M = ", M)
                end

                # Approximate large coupling asymptotic limit for the polaron mass. Feynman1962. Thermal
                M_asymp = v ./ w
                p["Ma"][i, d, j, :] .= M_asymp

                # Print asymptotic masses.
                if verbose
                    println(io, "\e[KFictitious mass (asymptotic)   | Ma = ", M_asymp)
                end

                # Reduced mass of particle and fictitious mass system. Before eqn. (2.4) in Schultz1959. Athermal.
                M_reduced = (v .^2 .- w .^2) ./ v .^ 2
                p["Mr"][i, d, j, :] .= M_reduced
            
                # Print redcued masses.
                if verbose
                    println(io, "\e[KReduced mass                   | Mr = ", M_reduced)
                end

                # Calculate and store polaron radii.
                R = sqrt.(3 ./ (2 .* M_reduced .* v))
                p["R"][i, d, j, :] .= R ./ sqrt(ωeff)
            
                # Print polaron radius.
                if verbose
                    println(io, "\e[KPolaron radius                 | R = ", R)
                end

                if verbose
                    println("\e[K-----------------------------------------------------------------------")
                    println("\e[K                      DC Mobility Information:                         ")
                    println("\e[K-----------------------------------------------------------------------")
                end

                # Calculate and store the DC mobiliies.
                μ = !kspace ? frohlich_mobility(v, w, α, ω, β; dims = dims[d]) : frohlich_mobility_k_space(v, w, α, ω, β; dims = dims[d])
                p["μ"][i, d, j] = μ 
                
                # Print DC mobilities.
                if verbose
                    println(io, "\e[KFinite temperature mobility    | μ = ", μ)
                end

                # FHIP low-temperature mobility, final result of Feynman1962.
                # Eqn. (1.60) in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv).
                μ_FHIP = FHIP_mobility_lowT(v[1] ./ ωeff, w[1] ./ ωeff, α, ω, β)
                p["μFHIP"][i, j] = μ_FHIP
                
                # Print low-temperature FHIP mobility.
                if verbose
                    println(io, "\e[KFHIP low-temp. mobility        | μFHIP = ", μ_FHIP)
                end

                # Kadanoff low-temperaure mobility, constructed around Boltzmann equation.
                # Adds factor of 3 / (2 * β) c.f. FHIP, correcting phonon emission behaviour.
                # Provide also the Kadanoff mobility that is consistent with the
                # FHIP, and later statements (Devreese) of the Kadanoff mobility.
                # It suggests that Kadanoff used the wrong identy for Nbar in Eqn. (23b) for
                # the Γ₀ function, and should have used a version with the -1 to
                # account for Bose / phonon statistics!
                μ_Kadanoff_Devreese, μ_Kadanoff, rel_time = Kadanoff_mobility_lowT(v[1] ./ ωeff, w[1] ./ ωeff, α, ω, β)
                p["μD"][i, j] = μ_Kadanoff_Devreese
                p["μK"][i, j] = μ_Kadanoff
                p["τ"][i, j] = rel_time
            
                # Print low-temperature Kadanoff mobility (both original and Devreese-corrected versions).
                if verbose
                    println(io, "\e[KDevreese low-temp. mobility    | μD = ", μ_Kadanoff_Devreese)
                    println(io, "\e[KKadanoff low-temp. mobility    | μK = ", μ_Kadanoff)
                end

                #C alculates the DC mobility using Hellwarth et al. 1999 Eqn. (2).
                # Directly performs contour integration in Feynman1962, for finite temperature DC mobility.
                # Eqns. (2) and (1) are going back to the general (pre low-T limit) formulas in Feynman1962.  
                # To evaluate these, you need to do the explicit contour integration to get the polaron self-energy.
                # See Hellwarth et a. 1999: https://doi.org/10.1103/PhysRevB.60.299.
                μ_Hellwarth, μ_Hellwarth_b0 = Hellwarth_mobility(v[1], w[1], α, ω, β)
                p["μH"][i, j] = μ_Hellwarth
                p["μH0"][i, j] = μ_Hellwarth_b0
            
                # Print Hellwarth mobility (both with b and b=0) and Kadanoff relaxation time.
                if verbose
                    println(io, "\e[KHellwarth mobility             | μH = ", μ_Hellwarth)
                    println(io, "\e[KHellwarth mobility (b=0)       | μH0 = ", μ_Hellwarth_b0)
                    println(io, "\e[KKadanoff relaxation time       | τ = ", rel_time)
                end

            else
                # If zero temperature.
                v, w, β = v_gs, w_gs, Inf
                p["β"][i] = β
                p["v"][i, d, j, :] .= v_gs
                p["w"][i, d, j, :] .= w_gs 
                p["F"][i, d, j] = F_gs
                p["A"][i, d, j] = A_gs
                p["B"][i, d, j] = B_gs
                p["C"][i, d, j] = C_gs
                p["κ"][i, d, j, :] .= κ_gs
                p["M"][i, d, j, :] .= M_gs
                p["μH"][i, j] = Inf
                p["μH0"][i, j] = Inf
                p["μD"][i, j] = Inf
                p["μK"][i, j] = Inf
                p["τ"][i, j] = 0
                p["μFHIP"][i, j] = Inf
                p["μ"][i, d, j] = Inf
                p["R"][i, d, j, :] .= R_gs
                p["Mr"][i, d, j, :] .= M_reduced_gs
                p["Ma"][i, d, j, :] .= M_asymp_gs
            end # End of zero temperature if statement.

            if verbose
                Ωprocess = 1    # Counter for Ωrange.
                Tprocess += 1   # Increment Trange iterator.
            end

            for k in eachindex(Ωrange)  # E-field frequencies loop. 
                Ω = Ωrange[k] 

                if !iszero(T) || !iszero(Ω)

                # Print E-field frequency.
                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K         Frequency Response Information: [$(Ωprocess[]) / $(num_Ω) ($(round(Ωprocess[] / (num_Ω) * 100, digits=1)) %)]")
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[KElectric field frequency       | Ω = ", Ω)
                end

                # Calculate and store polaron memory functions (akin to self energy).
                χ = !kspace ? frohlich_memory_function(Ω, v, w, α, ω, β; dims = dims[d]) : frohlich_memory_function_k_space(Ω, v, w, α, ω, β; dims = dims[d])
                p["χ"][k, i, d, j] = χ

                # Print memory function.
                if verbose
                    println(io, "\e[KMemory function                | χ = ", χ)
                end

                # Calculate and store polaron complex impedances.

                z = -(im * Ω + im * χ)
                p["z"][k, i, d, j] = z 

                # Print complex impedances.
                if verbose
                    println(io, "\e[KComplex impedance              | z = ", z)
                end

                # Calculate and store polaron complex conductivities.
                σ = 1 / z
                p["σ"][k, i, d, j] = σ 

                # Print complex conductivities and show total algorithm progress.
                if verbose
                    println(io, "\e[KComplex conductivity           | σ = ", σ)
                end

                else

                    # If zero frequency.
                    p["χ"][k, i, d, j] = Inf + 0 * im
                    p["z"][k, i, d, j] = 0 + im * Inf
                    p["σ"][k, i, d, j] = 0 + 0 * im

                end # End of zero temp if statement.

                if verbose
                    println(io, "\e[K-----------------------------------------------------------------------")
                    println(io, "\e[K[Total Progress: $(process[]) / $(num_α * num_T * num_Ω * num_d) ($(round(process[] / (num_α * num_T * num_Ω * num_d) * 100, digits=1)) %)]")
                    Ωprocess += 1   # Increment Ωrange iterator.
                    process += 1    # Increment total iterator.
                    print(io, "\e[2F")
                end

                if verbose && (!iszero(T) || !iszero(Ω)) print(io, "\e[7F") end
            end
            if verbose && !iszero(T) print(io, "\e[26F") end   # Move up 26 lines and erase.
        end 
        if verbose print(io, "\e[15F") end   # Move up 15 lines and erase. 
        end # End dimensions loop.

        if verbose print(io, "\e[12F") end   # Move up 12 lines and erase. 
    end
    if verbose print("\e[?25h") end # Show cursor again.

    # Extract polaron data as Vector ready for Polaron type.
    polaron_data = [
        p["α"], # Alphas.
        p["αeff"], # Alphas sums.
        p["mb"], # Band mass.
        p["T"], # Temperatures.
        p["ω"], # Phonon frequencies.
        p["ωeff"], # Hellwarth eff phonon frequency.
        p["β"], # Betas.
        p["Ω"], # Photon frequencies.
        p["d"], # Number of spatial dimensions.
        p["v0"], # Athermal v params.
        p["w0"], # Athermal w params.
        p["F0"], # Athermal energies.
        p["A0"], # Athermal A parameter.
        p["B0"], # Athermal B parameter.
        p["C0"], # Athermal C parameter.
        p["Fs"], # Small alpha (α→0) approximate energy.
        p["Fl"], # Large alpha (α→∞) approximate energy.
        p["κ0"], # Fictitious spring constant.
        p["M0"], # Athermal fictitious mass.
        p["Ms"], # Small alpha (α→0) approximate fictitious mass.
        p["Ml"], # Large alpha (α→∞) approximate fictitious mass.
        p["M0a"], # Athermal asymptotic approximate fictitious mass (v0/w0).
        p["M0r"], # Athermal reduced mass of particle + fictitious particle system.
        p["R0"], # Athermal polaron radius (s. d. of a Gaussian wavefunction).
        p["Rs"], # Small alpha (α→0) approximate polaron radius.
        p["Rl"], # Large alpha (α→∞) approximate polaron radius.
        p["ΩFC"], # Large alpha (α → ∞)approximate Franck-Condon peak frequency.
        p["v"], # v params.
        p["w"], # w params.
        p["F"], # Energies.
        p["A"], # A parameter.
        p["B"], # B parameter.
        p["C"], # C parameter.
        p["κ"], # Spring constants.
        p["M"], # Fictitious masses.
        p["Ma"],  # Thermal asymptotic approximate fictitious mass (v/w).
        p["Mr"],  # Thermal reduced mass of particle + fictitious particle system.
        p["R"], # Polaron radii.
        p["μ"], # Mobilities.
        p["μFHIP"], # FHIP low-temperature DC mobiltiy approximation.
        p["μD"], # Kadanoff DC mobility, low-temperature approximation using Devreese2016.
        p["μK"], # Kadanoff DC mobility, low-temperature approximation.
        p["τ"], # relaxation time from Kadanoff Boltzmann transport equation.
        p["μH"], # Hellwarth DC mobility.
        p["μH0"], # Hellwarth DC mobility with b=0.
        p["χ"], # Memory function or self energy.
        p["z"], # Complex impedences.
        p["σ"], # Complex conductivities.
    ]

    # Return Polaron type containing all generated data over any coupling strengths, temperatures and frequencies.
    return FrohlichPolaron(polaron_data...)
end

"""
Single alpha parameter. polaron() expects alpha parameters to be in a Vector.
"""
frohlichpolaron(α::Real, Trange, Ωrange; ω=1, ωeff=1, mb=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, cycleguesses=false, verbose=false) = frohlichpolaron([α], Trange, Ωrange; ω=ω, ωeff=ωeff, mb=mb, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, cycleguesses=cycleguesses, verbose=verbose)

"""
No frequency input.
"""
frohlichpolaron(αrange, Trange; ω=1, ωeff=1, mb=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, cycleguesses=false, verbose=false) = frohlichpolaron(αrange, Trange, 0; ω=ω, ωeff=ωeff, mb=mb, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, verbose=verbose)

"""
No temperature input => 300 K.
"""
frohlichpolaron(αrange; ω=1, ωeff=1, mb=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, cycleguesses=false, verbose=false) = frohlichpolaron(αrange, 0, 0; ω=ω, ωeff=ωeff, mb=mb,v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, cycleguesses=cycleguesses, verbose=verbose)

"""
No input => α = 1
"""
frohlichpolaron(; ω=1, ωeff=1, mb=1, v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, cycleguesses=false, verbose=false) = frohlichpolaron(1, 0, 0; ω=ω, ωeff=ωeff, mb=mb, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, cycleguesses=cycleguesses, verbose=verbose)

"""
    polaron(material::Material, TΩrange...; v_guesses=3.11, w_guesses=2.87, verbose=false)
Material specific constructors that use material specific parameters to parameterise the polaron.
Material data is inputted through the `Material` type.
Returns all data in either SI units or other common, suitable units otherwise.
"""
function frohlichpolaron(material::Material, TΩrange...; v_guesses=3.11, w_guesses=2.87, dims=3, kspace=false, cycleguesses=false, verbose=false)

    # Show material data.
    if verbose
        display(material)
    end
    
    # Extract material data from Material type.
    phonon_freqs = material.f .* pustrip(1u"THz2π")
    phonon_eff_freq = material.feff .* pustrip(1u"THz2π")
    mb = material.mb

    TΩrange = length(TΩrange) == 1 ? TΩrange .* pustrip(1u"K") : TΩrange .* (pustrip(1u"K"),1)

    # Generate polaron data from the arbitrary model constructor.
    p = frohlichpolaron(material.α', TΩrange...; ω=phonon_freqs, ωeff=phonon_eff_freq, mb=mb, v_guesses=v_guesses, w_guesses=w_guesses, dims=dims, kspace=kspace, cycleguesses=cycleguesses, verbose=verbose)

    # Return material-specific, unitful Polaron type.
    return p
end

# Broadcast Polaron data.
function Base.show(io::IO, ::MIME"text/plain", x::FrohlichPolaron)

    io_lim = IOContext(io, :limit => true, :compact => true)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                         Polaron Information:                          ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KPhonon frequencies             | ωeff = ", x.ωeff, " | ω = ", x.ω)
    println(io_lim, "\e[KFröhlich coupling              | αeff = ", x.αeff, " | α = ", x.α)
    println(io_lim, "\e[KNumber of spatial dimensions   | d = ", x.d)
    println(io_lim, "\e[KSmall α→0 energy               | Fs = ", x.Fs)
    println(io_lim, "\e[KLarge α→∞ energy               | Fl = ", x.Fl)
    println(io_lim, "\e[KSmall α→0 fictitious mass      | Ms = ", x.Ms)
    println(io_lim, "\e[KLarge α→∞ fictitious mass      | Ml = ", x.Ml)
    println(io_lim, "\e[KSmall α→0 polaron radius       | Rs = ", x.Rs)
    println(io_lim, "\e[KLarge α→∞ polaron radius       | Rl = ", x.Rl)
    println(io_lim, "\e[KLarge α→∞ FC peak freq.        | ΩFC = ", x.ΩFC)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                     Zero Temperature Information:                     ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KVariational parameter          | v0 = ", x.v0)
    println(io_lim, "\e[KVariational parameter          | w0 = ", x.w0)
    println(io_lim, "\e[KEnergy                         | E0 = ", x.F0)
    println(io_lim, "\e[KElectron energy                | A0 = ", x.A0)
    println(io_lim, "\e[KInteraction energy             | B0 = ", x.B0)
    println(io_lim, "\e[KTrial energy                   | C0 = ", x.C0)
    println(io_lim, "\e[KFictitious spring constant     | κ0 = ", x.κ0)
    println(io_lim, "\e[KFictitious mass                | M0 = ", x.M0)
    println(io_lim, "\e[KFictitious mass (asymptotic)   | M0a = ", x.M0a)
    println(io_lim, "\e[KReduced mass                   | M0r = ", x.M0r)
    println(io_lim, "\e[KPolaron radius                 | R0 = ", x.R0)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                    Finite Temperature Information:                    ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KTemperatures                   | T = ", x.T)
    println(io_lim, "\e[KReduced thermodynamic          | β = ", x.β)
    println(io_lim, "\e[KVariational parameter          | v = ", x.v)
    println(io_lim, "\e[KVariational parameter          | w = ", x.w)
    println(io_lim, "\e[KFree energy                    | F = ", x.F)
    println(io_lim, "\e[KElectron energy                | A = ", x.A)
    println(io_lim, "\e[KInteraction energy             | B = ", x.B)
    println(io_lim, "\e[KTrial energy                   | C = ", x.C)
    println(io_lim, "\e[KFictitious spring constant     | κ = ", x.κ)
    println(io_lim, "\e[KFictitious mass                | M = ", x.M)
    println(io_lim, "\e[KFictitious mass (asymptotic)   | Ma = ", x.Ma)
    println(io_lim, "\e[KReduced mass                   | Mr = ", x.Mr)
    println(io_lim, "\e[KPolaron radius                 | R = ", x.R)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                      DC Mobility Information:                         ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KFinite temperature mobility    | μ = ", x.μ)
    println(io_lim, "\e[KFHIP low-temp. mobility        | μFHIP = ", x.μFHIP)
    println(io_lim, "\e[KDevreese low-temp. mobility    | μD = ", x.μD)
    println(io_lim, "\e[KKadanoff low-temp. mobility    | μK = ", x.μK)
    println(io_lim, "\e[KHellwarth mobility             | μH = ", x.μH)
    println(io_lim, "\e[KHellwarth mobility (b=0)       | μH0 = ", x.μH0)
    println(io_lim, "\e[KKadanoff relaxation time       | τ = ", x.τ)

    println("\e[K-----------------------------------------------------------------------")
    println("\e[K                  Frequency Response Information:                      ")
    println("\e[K-----------------------------------------------------------------------")

    println(io_lim, "\e[KElectric field frequency       | Ω = ", x.Ω)
    println(io_lim, "\e[KMemory function                | χ = ", x.χ)
    println(io_lim, "\e[KComplex impedance              | z = ", x.z)
    println(io_lim, "\e[KComplex conductivity           | σ = ", x.σ)

    println("\e[K-----------------------------------------------------------------------\n")
end

"""
    save_polaron(p::FrohlichPolaron, prefix)

Saves data from 'polaron' into file "prefix".
This is a .jdl file for storing the polaron data whilst preserving types. Allows for saving multidimensional arrays that sometimes arise in the polaron data.
Each parameter in the NewPolaron type is saved as a dictionary entry. E.g. NewPolaron.α is saved under JLD.load("prefix.jld")["alpha"].
"""
function save_frohlich_polaron(polaron::FrohlichPolaron, prefix)

    println("Saving polaron data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "alpha", polaron.α,
        "alpha eff", polaron.αeff,
        "band mass", polaron.mb,
        "temperature", polaron.T,
        "phonon freq", polaron.ω,
        "phonon freq eff", polaron.ωeff,
        "beta", polaron.β,
        "E-field freq", polaron.Ω,
        "dims", polaron.d,
        "v0", polaron.v0,
        "w0", polaron.w0,
        "athermal energy", polaron.F0,
        "athermal A", polaron.A0,
        "athermal B", polaron.B0,
        "athermal C", polaron.C0,
        "small alpha energy", polaron.Fs,
        "large alpha energy", polaron.Fl,
        "athermal spring", polaron.κ0,
        "athermal mass", polaron.M0,
        "small alpha mass", polaron.Ms,
        "large alpha mass", polaron.Ml,
        "athermal asympt mass", polaron.M0a,
        "athermal reduced mass", polaron.M0r,
        "athermal size", polaron.R0,
        "small alpha size", polaron.Rs,
        "large alpha size", polaron.Rl,
        "FC freq", polaron.ΩFC,
        "v", polaron.v,
        "w", polaron.w,
        "thermal energy", polaron.F,
        "thermal A", polaron.A,
        "thermal B", polaron.B,
        "thermal C", polaron.C,
        "thermal spring", polaron.κ,
        "thermal mass", polaron.M,
        "thermal asympt mass", polaron.Ma,
        "thermal reduced mass", polaron.Mr,
        "thermal size", polaron.R,
        "mobility", polaron.μ,
        "FHIP mobility", polaron.μFHIP,
        "Devreese mobility", polaron.μD,
        "Kadanoff mobility", polaron.μK,
        "relaxation time", polaron.τ,
        "Hellwarth mobility", polaron.μH,
        "b=0 mobility", polaron.μH0,
        "memory function", polaron.χ,
        "impedance", polaron.z,
        "conductivity", polaron.σ
    )

    println("... Polaron data saved.")
end

"""
    load_polaron(p::FrohlichPolaron, prefix)

Loads data from file "polaron_file_path" into a NewPolaron type.
"""
function load_frohlich_polaron(polaron_file_path)

    println("Loading polaron data from $polaron_file_path ...")

    data = JLD.load("$polaron_file_path")

    polaron = FrohlichPolaron(
        data["alpha"],
        data["alpha eff"],
        data["band mass"],
        data["temperature"],
        data["phonon freq"],
        data["phonon freq eff"],
        data["beta"],
        data["E-field freq"],
        data["dims"],
        data["v0"],
        data["w0"],
        data["athermal energy"],
        data["athermal A"],
        data["athermal B"],
        data["athermal C"],
        data["small alpha energy"],
        data["large alpha energy"],
        data["athermal spring"],
        data["athermal mass"],
        data["small alpha mass"],
        data["large alpha mass"],
        data["athermal asympt mass"],
        data["athermal reduced mass"],
        data["athermal size"],
        data["small alpha size"],
        data["large alpha size"],
        data["FC freq"],
        data["v"],
        data["w"],
        data["thermal energy"],
        data["thermal A"],
        data["thermal B"],
        data["thermal C"],
        data["thermal spring"],
        data["thermal mass"],
        data["thermal asympt mass"],
        data["thermal reduced mass"],
        data["thermal size"],
        data["mobility"],
        data["FHIP mobility"],
        data["Devreese mobility"],
        data["Kadanoff mobility"],
        data["relaxation time"],
        data["Hellwarth mobility"],
        data["b=0 mobility"],
        data["memory function"],
        data["impedance"],
        data["conductivity"]
    )
    println("... Polaron loaded.")

    return polaron
end


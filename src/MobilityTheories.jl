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
polaronmobility(300, 4.5, 24.1, 2.25E12, 0.12)
```
"""
function polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)
    println("\n\nPolaron mobility for system ε_Inf=$ε_Inf, ε_S=$ε_S, freq=$freq,
                 effectivemass=$effectivemass; with Trange $Trange ...")

    # Internally we we use 'mb' for the 'band mass' in SI units, of the effecitve-mass of the electron
    mb=effectivemass*MassElectron
    ω = (2*pi)*freq*1e12 # angular-frequency, of the phonon mode

    α=frohlichalpha(ε_Inf, ε_S,  freq,    effectivemass)
    #α=2.395939683378253 # Hard coded; from MAPI params, 4.5, 24.1, 2.25THz, 0.12me

    v=7.1 # starting guess for v,w variational parameters
    w=6.5

    @printf("Polaron mobility input parameters: ε_Inf=%f ε_S=%f freq=%g α=%f \n",ε_Inf, ε_S, freq, α )
    @printf("Derived params in SI: ω =%g mb=%g \n",ω ,mb)

    # Empty struct for storing data
    # A slightly better way of doing this ф_ф ...
    p=Polaron()
    # populate data structure with (athermal) parameters supplied...
    append!(p.α,α) # appending so as not to mess with type immutability
    append!(p.mb,mb)
    append!(p.ω,ω)

    # We define βred as the subsuming the energy of the phonon; i.e. kbT c.f. ħω
    for T in Trange
        β=1/(kB*T)
        βred=ħ*ω*β
        append!(p.βred,βred)
        @printf("T: %f β: %.3g βred: %.3g ħω  = %.3g meV\t",T,β,βred, 1E3* ħ*ω  / q)

        if T==0
            v,w=feynmanvw(α, v=v,w=w)
        else
            v,w=feynmanvw(α,βred, v=v,w=w)
        end

        @printf("\n Polaron Parameters:  v= %.4f  w= %.4f",v,w)

        # From 1962 Feynman, definition of v and w in terms of the coupled Mass and spring-constant
        # See Page 1007, just after equation (18)
        # Units of M appear to be 'electron masses'
        # Unsure of units for k, spring coupling constant
        k=(v^2-w^2)
        M=(v^2-w^2)/w^2
        append!(p.k,k)
        append!(p.M,M)

        @printf("  ||   M=%f  k=%f\t",M,k)
        @printf("\n Bare-band effective mass: %f Polaron effective mass: %f Polaron mass enhancement: %f%%",effectivemass,effectivemass*(1+M),M*100)

        @printf("\n Polaron frequency (SI) v= %.2g Hz  w= %.2g Hz",
                v*ω /(2*pi),w*ω /(2*pi))

		# (46) in Feynman1955
		meSmallAlpha(α )=α /6 + 0.025*α ^2
		# (47) In Feynman1955
		meLargeAlpha(α )=16*α ^4 / (81*π ^4)
		#meLargeAlpha(α )=202*(α /10)^4
        if (verbose) # asymptotic solutions - not that interesting when you have the actual ones!
		    @printf("\n Feynman1955(46,47): meSmallAlpha(α)= %.3f meLargeAlpha(α)= %.3f",
                    meSmallAlpha(α),meLargeAlpha(α))
            @printf("\n Feynman1962: Approximate ~ Large alpha limit, v/w = %.2f  =~approx~= alpha^2 = %.2f ",
                    v/w,α ^2)
        end

        # POLARON SIZE
        @printf("\n Polaron size (rf), following Schultz1959. (s.d. of Gaussian polaron ψ )")
        # Schultz1959 - rather nicely he actually specifies everything down into units!
        # just before (2.4) in Schultz1959
        mu=((v^2-w^2)/v^2)
        # (2.4)
        rf=sqrt(3/(2*mu*v))
        # (2.4) SI scaling inferred from units in (2.5a) and Table II
        rfsi=rf*sqrt(2*me*ω )
        @printf("\n\t Schultz1959(2.4): rf= %g (int units) = %g m [SI]",rf,rfsi )
        append!(p.rfsi,rfsi)

        if (verbose)
        scale=sqrt(2*mb*ω) # Note we're using mb;
        #band effective-mass in SI units (i.e. meff*melectron)

        rfa=(3/(0.44*α ))^0.5 # As given in Schultz1959(2.5a), but that 0.44 is actually 4/9
        @printf("\n\t Schultz1959(2.5a) with 0.44: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]",rfa,scale*rfa )
        rfa=(3/((4/9)*α ))^0.5 # Rederived from Feynman1955, 8-8-2017; Yellow 2017.B Notebook pp.33-34
        @printf("\n\t Schultz1959(2.5a) with 4/9 re-derivation: Feynman α→0 expansion: rfa= %g (int units) = %g m [SI]",rfa,scale*rfa )
        append!(p.rfsmallalpha,scale*rfa)

        rfb=3*(pi/2)^0.5 * α
        @printf("\n\t Schultz1959(2.5b): Feynman α→∞ expansion: rf= %g (int units) = %g m [SI]",rfb,scale*rfb )

        # Schultz1959 - Between (5.7) and (5.8) - resonance of Feynman SHM system
        phononfreq=sqrt(k/M)
        @printf("\n\t Schultz1959(5.7-5.8): fixed-e: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)

        phononfreq=sqrt(k/mu) # reduced mass
        @printf("\n\t Schultz1959(5.7-5.8): reducd mass: phononfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            phononfreq,phononfreq*ω /(2*pi), phononfreq*hbar*ω *1000/q)

        @printf("\n\t Schultz1959: electronfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            sqrt(k/1),sqrt(k/1)*ω /(2*pi), sqrt(k/1)*hbar*ω *1000/q)
        @printf("\n\t Schultz1959: combinedfreq= %g (int units) = %g [SI, Hz] = %g [meV]",
            sqrt(k/(1+M)),sqrt(k/(1+M))*ω /(2*pi), sqrt(k/(1+M))*hbar*ω *1000/q)

        # Devreese1972: 10.1103/PhysRevB.5.2367
        # p.2371, RHS.
        @printf("\n Devreese1972: (Large Alpha) Franck-Condon frequency = %.2f", 4*α ^2/(9*pi))
        end

        # F(v,w,β,α)=-(A(v,w,β)+B(v,w,β,α)+C(v,w,β)) #(62a) - Hellwarth 1999
        if (T>0)
            @printf("\n Polaron Free Energy: A= %f B= %f C= %f F= %f",A(v,w,βred),B(v,w,βred,α),C(v,w,βred),F(v,w,βred,α))
            @printf("\t = %f meV",1000.0 * F(v,w,βred,α) * ħ*ω  / q)
            append!(p.A,A(v,w,βred))
            append!(p.B,B(v,w,βred,α))
            append!(p.C,C(v,w,βred))
            append!(p.F,F(v,w,βred,α))
        else # Athermal case; Enthalpy
            @printf("\n Polaron Enthalpy: F= %f = %f meV \n",F(v,w,α), 1_000*F(v,w,α) * ħ*ω/q)
        
            return # return early, as if T=0, all mobility theories = infinite / fall over
        end
        
        # FHIP
        #    - low-T mobility, final result of Feynman1962
        # [1.60] in Devreese2016 page 36; 6th Edition of Frohlich polaron notes (ArXiv)
        # I believe here β is in SI (expanded) units
        @printf("\n Polaron Mobility theories:")
        μ=(w/v)^3 * (3*q)/(4*mb*ħ*ω^2*α*β) * exp(ħ*ω*β) * exp((v^2-w^2)/(w^2*v))
        @printf("\n\tμ(FHIP)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        append!(p.T,T)
        append!(p.FHIPμ,μ*100^2)


        # Kadanoff
        #     - low-T mobility, constructed around Boltzmann equation.
        #     - Adds factor of 3/(2*beta) c.f. FHIP, correcting phonon emission behaviour
        # [1.61] in Devreese2016 - Kadanoff's Boltzmann eqn derived mob
        μ=(w/v)^3 * (q)/(2*mb*ω*α) * exp(ħ*ω*β) * exp((v^2-w^2)/(w^2*v))
        @printf("\n\tμ(Kadanoff,via Devreese2016)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)

        append!(p.Kμ,μ*100^2)

        ######
        # OK, now deep-diving into Kadanoff1963 itself to extract Boltzmann equation components
        # Particularly right-hand-side of page 1367
        #
        # From 1963 Kadanoff, (9), define eqm. number of phonons (just from T and phonon omega)
        Nbar=(exp(βred) - 1)^-1
        @printf("\n\t\tEqm. Phonon. pop. Nbar: %f ",Nbar)
        if (verbose)
            @printf("\n\texp(Bred): %f exp(-Bred): %f exp(Bred)-1: %f",exp(βred),exp(-βred),exp(βred)-1)
        end
        Nbar=exp(-βred)
        #Note - this is only way to get Kadanoff1963 to be self-consistent with
        #FHIP, and later statements (Devreese) of the Kadanoff mobility.
        #It suggests that Kadanoff used the wrong identy for Nbar in 23(b) for
        #the Gamma0 function, and should have used a version with the -1 to
        #account for Bose / phonon statistics

        #myv=sqrt(k*(1+M)/M) # cross-check maths between different papers
        #@printf("\nv: %f myv: %f\n",v,myv)

        # Between 23 and 24 in Kadanoff 1963, for small momenta skip intergration --> Gamma0
        Gamma0=2*α * Nbar * (M+1)^(1/2) * exp(-M/v)
        Gamma0*=ω  #* (ω *hbar) # Kadanoff 1963 uses hbar=omega=mb=1 units
            # Factor of omega to get it as a rate relative to phonon frequency
            # Factor of omega*hbar to get it as a rate per energy window
        μ=q/( mb*(M+1) * Gamma0 ) #(25) Kadanoff 1963, with SI effective mass

        if (verbose) # these are cross-checks
            @printf("\n\tμ(Kadanoff1963 [Eqn. 25]) = %f m^2/Vs \t = %.2f cm^2/Vs",μ,μ*100^2)
            @printf("\n\t\t Eqm. Phonon. pop. Nbar: %f ",Nbar)
        end

        @printf("\n\t\tGamma0 = %g rad/s = %g /s ",
                Gamma0, Gamma0/(2*pi))
        @printf(" \n\t\tTau=1/Gamma0 = %g s = %f ps",
            2*pi/Gamma0, 2*pi*1E12/Gamma0)
        Eloss=hbar*ω * Gamma0/(2*pi) # Simply Energy * Rate
        @printf("\n\t\tEnergy Loss = %g J/s = %g meV/ps",Eloss,Eloss * 1E3 / (q*1E12) )
        append!(p.Tau, 2*pi*1E12/Gamma0) # Boosted into ps ?


        # Hellwarth1999 - directly do contour integration in Feynman1962, for
        # finite temperature DC mobility
        # Hellwarth1999 Eqn (2) and (1) - These are going back to the general
        # (pre low-T limit) formulas in Feynman1962.  to evaluate these, you
        # need to do the explicit contour integration to get the polaron
        # self-energy
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)

        #b=R*βred/sinh(b*βred*v/2) # This self-references b! What on Earth?
        # OK! I now understand that there is a typo in Hellwarth1999 and
        # Biaggio1997. They've introduced a spurious b on the R.H.S. compared to
        # the original, Feynman1962:
        b=R*βred/sinh(βred*v/2) # Feynman1962 version; page 1010, Eqn (47b)

        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        append!(p.Hμ,μ*100^2)

        if (verbose)
        # Hellwarth1999/Biaggio1997, b=0 version... 'Setting b=0 makes less than 0.1% error'
        # So let's test this
        R=(v^2-w^2)/(w^2*v) # inline, page 300 just after Eqn (2)
        b=0
        a=sqrt( (βred/2)^2 + R*βred*coth(βred*v/2))
        #k(u,a,b,v) = (u^2+a^2-b*cos(v*u))^(-3/2)*cos(u) # integrand in (2)
        K=quadgk(u->k(u,a,b,v),0,Inf)[1] # numerical quadrature integration of (2)

        #Right-hand-side of Eqn 1 in Hellwarth 1999 // Eqn (4) in Baggio1997
        RHS= α/(3*sqrt(π)) * βred^(5/2) / sinh(βred/2) * (v^3/w^3) * K
        μ=RHS^-1 * (q)/(ω*mb)
        @printf("\n\tμ(Hellwarth1999,b=0)= %f m^2/Vs \t= %.2f cm^2/Vs",μ,μ*100^2)
        @printf("\n\tError due to b=0; %f",(100^2*μ-p.Hμ[length(p.Hμ)])/(100^2*μ))
        #append!(Hμs,μ*100^2)
        end
        @printf("\n") # blank line at end of spiel.

        # Recycle previous variation results (v,w) as next guess
        initial=[v,w] # Caution! Might cause weird sticking in local minima

        append!(p.v,v)
        append!(p.w,w)

    end

    return(p)
end

"""
    make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff; temp = 300.0, efield_freq = 0.0, volume = nothing, ir_activity = nothing, N_params = 1, rtol = 1e-3, verbose = false, threads = false)

Solves the Feynman polaron problem variationally with finite temperature
Osaka energies. From the resulting v and w parameters, calculates polaron
structure (wave function size, etc.). 

Uses FHIP to predict a temperature- and frequency-dependent mobility, complex conductivity and impedance for the material system.

Can evaluate polaron properties for materials with multiple phonon branches using infrared activities and phonon branch frequencies.

# Arguments
- `ϵ_optic`: reduced optical dielectric constant.
- `ϵ_static`: reduced static dielectric constant.
- `phonon_freq`: vector of characteristic dielectric phonon frequencies (THz).
- `m_eff`: bare-band effective-mass (mₑ).
- `Trange`: temperature value or range.
- `Ω`: electric field frequency value or range (THz).
- `volume`: unit cell volume for the material (m³). `nothing` for one phonon mode.
- `ir_activity`: vector of infrared activities. `nothing` for one phonon  mode.
- `N_params`: number of variational parameters to minimise the polaron energy.
- `rtol`: relative tolerance for the accuracy of any quadrature integrations.
- 'verbose': `true` to print progress. `false` to print nothing.
- 'threads'; `true` activates multithreading through use of @tullio macro. Tullio utilises einsum notation.

Returns a structure of type NewPolaron, containing arrays of useful
information.  Also prints a lot of information to the standard out - which
may be more useful if you're just inquiring as to a particular data point,
rather than plotting a temperature-dependent parameter.

As an example, to calculate the electron polaron in MAPI, at temperatures 0:100:400 K and electric field frequencies 0.0:0.1:5.0 THz, and inclusive of 15 optical phonon modes:

# Examples
```jldoctest
make_polaron(
    4.5, 
    24.1, 
    [4.016471586720514, 3.887605410774121, 3.5313112232401513, 2.755392921480459, 2.4380741812443247, 2.2490917637719408, 2.079632190634424, 2.0336707697261187, 1.5673011873879714, 1.0188379384951798, 1.0022960504442775, 0.9970130778462072, 0.9201781906386209, 0.800604081794174, 0.5738689505255512], 
    0.12; 
    temp = 0.0:100.0:400.0, 
    efield_freq = 0.0:0.1:5.0, 
    volume = (6.29e-10)^3, 
    ir_activity = [0.08168931020200264, 0.006311654262282101, 0.05353548710183397, 0.021303020776321225, 0.23162784335484837, 0.2622203718355982, 0.23382298607799906, 0.0623239656843172, 0.0367465760261409, 0.0126328938653956, 0.006817361620021601, 0.0103757951973341, 0.01095811116040592, 0.0016830270365341532, 0.00646428491253749], 
    N_params = 1)
```
    
"""
function make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, Trange, Ω; volume = nothing, ir_activity = nothing, rtol = 1e-4, N_params = 1, verbose = false, threads = false)
   
    # Number of phonon modes.
    N_modes = length(phonon_freq)

    if verbose
        println("Calculating Fröhlich alpha...")
    end

    # Convert THz phonon frequencies to 2π THz.
    ω = 2π .* phonon_freq
    
    if N_modes == 1 
        # One phonon mode.

        # Calculate Frohlich alpha coupling parameters.
        α = frohlichalpha(ϵ_optic, ϵ_static, phonon_freq, m_eff)
    else 
        # Multiple phonon modes

        # Calculate contribution to the ionic delectric constant for each phonon mode.
        ϵ_ionic = [ϵ_ionic_mode(i, j, volume) for (i, j) in zip(phonon_freq, ir_activity)]

        # Calculate contribution to Frohlich alpha for each phonon mode.
        α = [multi_frohlichalpha(ϵ_optic, i, sum(ϵ_ionic), j, m_eff) for (i, j) in zip(ϵ_ionic, phonon_freq)]
    end

    if verbose
        show(IOContext(stdout, :limit => true), round.(α, digits = 3))
        print("\n\n")
        println("Calculating thermodynamic betas...")
    end

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    @tullio threads = threads betas[m, i] := Trange[i] == 0.0 ? Inf64 : ħ * ω[m] / (kB * Trange[i]) * 1e12 (i in eachindex(Trange))
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(betas, digits = 3))
        print("\n\n")
        println("Calculating variational parameters...")
        global count = 0
        global processes = length(Trange)
    end

    # Calculate variational parameters for each temperature from multiple phonon frequencies.
    @tullio threads = threads params[i] := Trange[i] == 0.0 ? var_params(α; v = 5.6, w = 3.4, ω = ω, rtol = rtol, N = N_params, verbose = verbose) : var_params(α, betas[:, i]; v = 5.6, w = 3.4, ω = ω, rtol = rtol, T = Trange[i], N = N_params, verbose = verbose) (i in eachindex(Trange))

    # Separate tuples of variational parameters into a list of 'v' and 'w' parameters for each temperature.
    @tullio threads = threads v_params[i] := params[i][1] (i in eachindex(Trange))
    @tullio threads = threads w_params[i] := params[i][2] (i in eachindex(Trange))
    
    if verbose
        println("\nv: ")
        show(IOContext(stdout, :limit => true), round.(v_params, digits = 3))
        print("\n\n")
        println("w: ")
        show(IOContext(stdout, :limit => true), round.(w_params, digits = 3))
        print("\n\n")
        println("Calculating spring constants...")
    end

    # Calculate fictitious spring constants for each temperature.
    @tullio threads = threads spring_constants[i] := v_params[i]^2 - w_params[i]^2

    if verbose
        show(IOContext(stdout, :limit => true), round.(spring_constants, digits = 3))
        print("\n\n")
        println("Calculating fictitious masses...")
    end

    # Calculate fictitious masses for each temperature.
    @tullio threads = threads masses[i] := spring_constants[i] / w_params[i]^2
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(masses, digits = 3))
        print("\n\n")
        println("Calculating free energies...")
        global count = 0
        global processes = length(Trange)
    end

    # Calculate ground-state free energies for each temperature.
    @tullio threads = threads energies[i] := Trange[i] == 0.0 ? multi_F(v_params[i], w_params[i], α; ω = ω, rtol = rtol, verbose = verbose) * 1000 * ħ / eV * 1e12 : multi_F(v_params[i], w_params[i], α, betas[:, i]; ω = ω, rtol = rtol, T = Trange[i],  verbose = verbose) * 1000 * ħ / eV * 1e12 (i in eachindex(Trange))
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(energies, digits = 3))
        print("\n\n")
        println("Calculating complex impedances...")
        global count = 0
        global processes = length(Trange) * length(Ω)
    end

    # Calculate complex impedances for each frequency and temperature. Returns a matrix.
    @tullio threads = threads impedances[f, i] := Ω[f] == Trange[i] == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ω[f], betas[:, i], α, v_params[i], w_params[i]; ω = ω, rtol = rtol, T = Trange[i], verbose = verbose) / eV^2 * 1e12 * me * m_eff * volume * 100^3 (i in eachindex(Trange), f in eachindex(Ω))

    if verbose
        show(IOContext(stdout, :limit => true), round.(impedances, digits = 3))
        print("\n\n")
        println("Calculating complex conductivities...")
    end

    # Calculate complex conductivities for each frequency and temperature. Returns a matrix.
    @tullio threads = threads conductivities[f, i] := Ω[f] == Trange[i] == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedances[f, i] (i in eachindex(Trange), f in eachindex(Ω))

    if verbose
        show(IOContext(stdout, :limit => true), round.(conductivities, digits = 3))
        print("\n\n")
        println("Calculating mobilities...")
        global count = 0
        global processes = length(Trange)
    end

    # Calculates the dc mobility for each temperature.
    @tullio threads = threads mobilities[i] := Trange[i] == 0.0 ? Inf64 : polaron_mobility(betas[:, i], α, v_params[i], w_params[i]; ω = ω, rtol = rtol, T = Trange[i], verbose = verbose) * eV / (1e12 * me * m_eff) * 100^2 (i in eachindex(Trange))
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(mobilities, digits = 3))
        print("\n")
    end

    # Return Polaron mutable struct with evaluated data.
    return NewPolaron(α, Trange, betas, phonon_freq, v_params, w_params, spring_constants, masses, energies, Ω, impedances, conductivities, mobilities)
end

"""
make_polaron(α, Trange, Ω; ω = 1.0, rtol = 1e-4, verbose = false, threads = false)

Same as above but from a model system with specified alpha values rather than from material properties. Here we only have one phonon mode with a normalised frequency `ω = 1.0`.
"""
function make_polaron(α, Trange, Ω; ω = 1.0, rtol = 1e-4, verbose = false, threads = false)

    if verbose
        println("\nCalculating thermodynamic betas...")
    end

    # Calculate reduced thermodynamic betas for each phonon mode.
    # If temperature is zero, return Inf.
    @tullio threads = threads betas[i] := Trange[i] == 0.0 ? Inf64 : ω / Trange[i]
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(betas, digits = 3))
        print("\n\n")
        println("Calculating variational parameters...")
        global count = 0
        global processes = length(Trange) * length(α)
    end

    # Calculate variational parameters for each alpha parameter and temperature. Returns a Matrix of tuples.
    @tullio threads = threads params[j, i] := Trange[i] == 0.0 ? var_params(α[j]; v = 5.6, w = 3.4, ω = ω, rtol = rtol, verbose = verbose) : var_params(α[j], betas[i]; v = 5.6, w = 3.4, ω = ω, rtol = rtol, T = Trange[i], verbose = verbose) (i in eachindex(Trange), j in eachindex(α))

    # Separate tuples of variational parameters into a Matrices of 'v' and 'w' parameters for each alpha parameter and temperature.
    @tullio threads = threads v_params[j, i] := params[j, i][1]
    @tullio threads = threads w_params[j, i] := params[j, i][2]
    
    if verbose
        println("\nv: ")
        show(IOContext(stdout, :limit => true), round.(v_params, digits = 3))
        print("\n\n")
        println("w: ")
        show(IOContext(stdout, :limit => true), round.(w_params, digits = 3))
        print("\n\n")
        println("Calculating spring constants...")
    end

    # Calculate fictitious spring constants for each alpha parameter and temperature. Returns a Matrix.
    @tullio threads = threads spring_constants[j, i] := v_params[j, i]^2 - w_params[j, i]^2

    if verbose
        show(IOContext(stdout, :limit => true), round.(spring_constants, digits = 3))
        print("\n\n")
        println("Calculating fictitious masses...")
    end

    # Calculate fictitious masses for each alpha parameter and temperature. Returns a Matrix.
    @tullio threads = threads masses[j, i] := spring_constants[j, i] / w_params[j, i]^2
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(masses, digits = 3))
        print("\n\n")
        println("Calculating free energies...")
        global count = 0
        global processes = length(T) * length(α)
    end

    # Calculate ground-state free energies for each alpha parameter and temperature. Returns a Matrix.
    @tullio threads = threads energies[j, i] := Trange[i] == 0.0 ? multi_F(v_params[j, i], w_params[j, i], α[j]; ω = ω, rtol = rtol, verbose = verbose) : multi_F(v_params[j, i], w_params[j, i], α[j], betas[i]; ω = ω, rtol = rtol, T = Trange[i],  verbose = verbose) (i in eachindex(Trange), j in eachindex(α))
    
    if verbose
        show(IOContext(stdout, :limit => true), round.(energies, digits = 3))
        print("\n\n")
        println("Calculating complex impedances...")
        global count = 0
        global processes = length(T) * length(Ω) * length(α)
    end

    # Calculate complex impedances for each alpha parameter, frequency and temperature. Returns a 3D Array.
    @tullio threads = threads impedances[k, j, i] := Ω[j] == Trange[i] == 0.0 ? 0.0 + 1im * 0.0 : polaron_complex_impedence(Ω[j], [betas[i]], α[k], v_params[k, i], w_params[k, i]; ω = ω, rtol = rtol, T = Trange[i],  verbose = verbose) (i in eachindex(Trange), j in eachindex(Ω), k in eachindex(α))

    if verbose
        show(IOContext(stdout, :limit => true), round.(impedances, digits = 3))
        print("\n\n")
        println("Calculating complex conductivities...")
    end

    # Calculate complex conductivities for each alpha parameter, frequency and temperature. Returns a 3D array.
    @tullio threads = threads conductivities[k, j, i] := Ω[j] == Trange[i] == 0.0 ? Inf64 + 1im * 0.0 : 1 / impedances[k, j, i] (i in eachindex(Trange), j in eachindex(Ω), k in eachindex(α))

    if verbose
        show(IOContext(stdout, :limit => true), round.(conductivities, digits = 3))
        print("\n\n")
        println("Calculating mobilities...")
        global count = 0
        global processes = length(T) * length(α)
    end

    # Calculates the dc mobility for each alpha parameter and each temperature.
    @tullio threads = threads mobilities[j, i] := Trange[i] == 0.0 ? Inf64 : polaron_mobility(betas[i], α[j], v_params[j, i], w_params[j, i]; ω = ω, rtol = rtol, T = Trange[i], verbose = verbose) (i in eachindex(Trange), j in eachindex(α))

    if verbose
        show(IOContext(stdout, :limit => true), round.(mobilities, digits = 3))
        print("\n")
    end

    # Return Polaron mutable struct with evaluated data.
    return NewPolaron(α, Trange, betas, ω, v_params, w_params, spring_constants, masses, energies, Ω, impedances, conductivities, mobilities)
end

"""
    save_polaron(p::NewPolaron, prefix)

Saves data from 'polaron' into file "prefix".
This is a .jdl file for storing the polaron data whilst preserving types. Allows for saving multidimensional arrays that sometimes arise in the polaron data.
Each parameter in the NewPolaron type is saved as a dictionary entry. E.g. NewPolaron.α is saved under JLD.load("prefix.jld")["alpha"].
"""
function save_polaron(polaron::NewPolaron, prefix)

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
    f=open("$fileprefix.dat","w")

    @printf(f,"# %s\n",fileprefix) # put name / material at header
    @printf(f,"# Params in SI: ω =%g mb=%g \n",p.ω[1] ,p.mb[1])
    @printf(f,"# Alpha parameter: α = %f  \n",p.α[1] )

    @printf(f,"# Ts, βreds, Kμs, Hμs, FHIPμs, vs, ws, ks, Ms, As, Bs, Cs, Fs, Taus, rfsis\n")
    @printf(f,"#  1    2     3    4     5      6   7   8   9  10  11  12  13    14     15\n") # columns for GNUPLOT etc.

    for i in 1:length(p.T)
        @printf(f,"%d %03f %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
        p.T[i], p.βred[i], p.Kμ[i], p.Hμ[i], p.FHIPμ[i],
        p.v[i], p.w[i],
        p.k[i], p.M[i], p.A[i], p.B[i], p.C[i], p.F[i],
        p.Tau[i], p.rfsi[i])
    end
    close(f)
end

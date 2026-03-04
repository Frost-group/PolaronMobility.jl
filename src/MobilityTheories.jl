# MobilityTheories.jl

"""
    polaronmobility(Trange,
                    ε_Inf, ε_S, freq, effectivemass;
                    verbose::Bool=false)

    Solves the Feynman polaron problem variationally with finite temperature
    Osaka energies.  From the resulting v, and w parameters, calculates polaron
    structure (wave function size, etc.).  Uses FHIP, Kadanoff (Boltzmann
    relaxation time) and Hellwarth direct contour integration to predict
    a temperature-dependent mobility for the material system.
    Input is a temperature range (e.g. 10:50:1000),
    reduced dielectric constants (e.g. 5, 20),
    characteristic dielectric phonon frequency (e.g. 2.25E12) - units Hertz
    bare-band effective-mass (e.g. 012) - units electron mass.

    Returns a structure of type Polaron, containing arrays of useful
    information.  Also prints a lot of information to the standard out - which
    may be more useful if you're just inquiring as to a particular data point,
    rather than plotting a temperature-dependent parameter.

    As an example, to calculate the electron polaron in MAPI at 300 K:
    polaronmobility(300, 4.5, 24.1, 2.25E12, 0.12)
"""
function polaronmobility(Trange, ε_Inf, ε_S, freq, effectivemass; verbose::Bool=false)
    println("\n\nPolaron mobility for system ε_Inf=$ε_Inf, ε_S=$ε_S, freq=$freq,
                 effectivemass=$effectivemass; with Trange $Trange ...")

    # Internally we we use 'mb' for the 'band mass' in SI units, of the effecitve-mass of the electron
    mb=effectivemass*MassElectron
    ω = (2*pi)*freq # angular-frequency, of the phonon mode

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

        mott = 1e-6 * ( 2 * (2 * rfsi)^3)^-1
        @printf("\n\t Mott critereon (cm^-3): %g",mott)

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


# PolaronMobility.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![made-with-julia](https://img.shields.io/badge/Made%20with-Julia-ff69bf.svg)](https://julialang.org)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00566/status.svg)](https://doi.org/10.21105/joss.00566)
[![docs-latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://Frost-group.github.io/PolaronMobility.jl/)


[![Build status](https://github.com/Frost-group/PolaronMobility.jl/workflows/CI/badge.svg)](https://github.com/Frost-group/PolaronMobility.jl/actions)
[![codecov.io](http://codecov.io/github/Frost-group/PolaronMobility.jl/coverage.svg?branch=main)](http://codecov.io/github/Frost-group/PolaronMobility.jl?branch=main)

`PolaronMobility.jl` is a Julia package which calculates the
temperature-dependent polaron mobility for a material. 

This is based on the Feynman variational solution to the Polaron problem. 
The electron-phonon coupling is treated as an effective α (alpha) Frohlich
Hamiltonian dimensionless parameter. 
The band structure is treated with an effective mass theory. 
The variational problem is solved numerically for finite-temperature free
energies. 
(The original 1960s work, and thus textbook solutions, often use asymptotic approximations to the integrals, with a more simple athermal action.)   
The mobility is calculated in three ways:
1) numerically by integrating the polaron self-energy along the imaginary axis (`Hellwarth1999`)
2) using Kadanoff's Boltzmann equation approximation (`Kadanoff1963`)
3) using the FHIP low-temperature asymptotic solution (`FHIP`)

These three methods are in approximately descending order of accuracy. 

We provide parameters for various metal-halide Perovskites, and other
interesting systems.

The motivation for developing these codes was to enable polaron mobility
calculations on arbitrary materials. 
They also provide the only extant implementation of Feynman's variational
method.  
They offer a convenient basis for writing codes that build on these variational
solutions. 

More [extensive documentation](https://Frost-group.github.io/PolaronMobility.jl/),
is perhaps easiest to read and understand alongside the first paper:
[ArXiv:1704.05404](https://arxiv.org/abs/1704.05404)
/ [Frost2017PRB](https://doi.org/10.1103/PhysRevB.96.195202).


## Installation

To install, type the following at the Julia (>1.0) REPL:

```
julia> import Pkg; Pkg.add("PolaronMobility")
```

## Cloud notebook

There is an [example notebook](JuliaBox-Example.ipynb) which can be run interactively on the (free) MyBinder notebook server. This is the fastest way to calculate a few polaron parameters, if you do not have Julia installed locally.

1) Click on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Frost-group/PolaronMobility.jl/master?filepath=JuliaBox-Example.ipynb)
2) That's it!

(Currently plotting does not work, as the Docker image is not built with the (heavy weight) Plots dependency, and I'm not sure how I can do this just for MyBinder, without requiring it generally for PolaronMobility.jl. If this is problematic for you, please open an issue and I'll try to fix it!)

## Using

As an example:

```
julia> using PolaronMobility
julia> MAPIe = material(4.5, 24.1, 0.12, 2.25)
julia> MAPIe_polaron = frohlichpolaron(MAPIe, 300; verbose = true)
```

Will calculate the polaron mobility for methyl-ammonium lead halide perovskite
(ϵ_optical = 4.5, ϵ_static = 24.1, effective_mass = 0.12 electron-masses, f = 2.25 THz) at 300 K. 

An abbreviated output should look like:
```
------------------------------------------
           Material Information           
------------------------------------------
Optic dielectric   | ϵo = 4.5
Static dielectric  | ϵs = 24.1
Ionic dielectric   | ϵi = 19.6
Band mass          | mb = 0.12
Fröhlich coupling  | α = 2.39394
Phonon frequencies | f = 2.25
Eff Phonon freq    | feff = 2.25
IR activities      | ir = 1
Unit cell volume   | V = 1
-------------------------------------------

-----------------------------------------------------------------------
                         Polaron Information:                          
-----------------------------------------------------------------------
Phonon frequencies             | ωeff = 2.25 | ω = 2.25
Fröhlich coupling              | αeff = 2.39394 | α = 2.39394
Number of spatial dimensions   | d = 3
Small α→0 energy               | Fs = -2.46469
Large α→∞ energy               | Fl = -3.43751
Small α→0 fictitious mass      | Ms = 0.542264
Large α→∞ fictitious mass      | Ml = 0.0666023
Small α→0 polaron radius       | Rs = 1.67917
Large α→∞ polaron radius       | Rl = 9.00108
Large α→∞ FC peak freq.        | ΩFC = 0.810764
-----------------------------------------------------------------------
                     Zero Temperature Information:                     
-----------------------------------------------------------------------
Variational parameter          | v0 = 3.30877
Variational parameter          | w0 = 2.66327
Energy                         | E0 = -5.56947
Electron energy                | A0 = -2.17859
Interaction energy             | B0 = 5.78198
Trial energy                   | C0 = 1.96608
Fictitious spring constant     | κ0 = 3.855
Fictitious mass                | M0 = 0.543495
Fictitious mass (asymptotic)   | M0a = 1.24237
Reduced mass                   | M0r = 0.35212
Polaron radius                 | R0 = 0.817278
-----------------------------------------------------------------------
                    Finite Temperature Information:                    
-----------------------------------------------------------------------
Temperatures                   | T = 300
Reduced thermodynamic          | β = 0.159975
Variational parameter          | v = 19.8612
Variational parameter          | w = 16.9599
Free energy                    | F = -8.59191
Electron energy                | A = -14.5096
Interaction energy             | B = 16.5498
Trial energy                   | C = 6.55175
Fictitious spring constant     | κ = 106.831
Fictitious mass                | M = 0.371408
Fictitious mass (asymptotic)   | Ma = 1.17107
Reduced mass                   | Mr = 0.270822
Polaron radius                 | R = 0.0722548
-----------------------------------------------------------------------
                      DC Mobility Information:                         
-----------------------------------------------------------------------
Finite temperature mobility    | μ = 0.487357
FHIP low-temp. mobility        | μFHIP = 2.93118
Devreese low-temp. mobility    | μD = 0.703373
Kadanoff low-temp. mobility    | μK = 0.703373
Hellwarth mobility             | μH = 1.09655
Hellwarth mobility (b=0)       | μH0 = 1.09824
Kadanoff relaxation time       | τ = 0.964611
-----------------------------------------------------------------------
```

The `frohlichpolaron` method returns a `Frohlich` struct type. To access a value we can call the variable from the table above.
For example, to get the mobility for MAPIe we do: 
```
julia> MAPIe_polaron.μ
0.48735686547867546
```

This package makes use of Unitful.jl to convert these values into unitful polaron types.
For example, to add units to the above MAPIe polaron:
```
julia> using Unitful
julia> addunits!(MAPIe_polaron)
julia> MAPIe_polaron.μ
0.48735686547867546 μ₀
```

All the values above are given in "polaron" units. We can convert to a specific unit using Unitful.
For example, to obtain the MAPIe polaron ground-state energy (meV), room-temperature free energy (meV) and mobility (cm^2/V/s):
```
julia> MAPIe_polaron.F0 |> u"meV"
-23.03349215713356 meV

julia> MAPIe_polaron.F |> u"meV"
-35.5332946810453 meV

julia> MAPIe_polaron.μ |> u"cm^2/V/s"
136.42332082324646 cm² s⁻¹ V⁻¹
```

The `frohlichpolaron` method can generally accept a range of α electron-phonon coupling parameters, a range of temperatures for finite temperature properties such as free energy and mobility, and range of Electric Field frequencies for calculating dynamical properties such as the complex conductivity.
For example:
```
julia> frohlich_polaron = frohlichpolaron(1:12, 1:400, 0:0.1:10)
```
will calculate all properties for α's 1 to 12, temperature 1K to 400K and Electric field frequencies 0 to 10 THz. 

It is also possible to change the number of spatial dimensions for the polaron using `dims = #` or as an array/range `dims = [2, 3]`. 
For example:
```
julia> frohlich_polaron = frohlichpolaron(1:12, 1:400, 0:0.1:10; dims = 2)
julia> frohlich_polaron = frohlichpolaron(1:12, 1:400, 0:0.1:10; dims = [2, 3])
```

Note that the Frohlich polaron does not exist in 1D.

TIP: If the `frohlichpolaron` errors or you interrupt the calculation whilst you have `verbose = true`, use `print("\e[?25h")` to return your cursor.

## Multiple phonon modes

To calculate polaron properties for a material with multiple phonon modes the code requires a vector of phonon frequencies (at the gamma point), a vector of corresponding infrared activities for each phonon mode, and the unitcell volume of the material. The code uses this additional information to calculate the ionic dielectric contributions to the static dielectric function.
For example, for MAPIe:
```
julia> volume = (6.29e-10)^3
julia> m_eff = 0.12
julia> ϵ_static = 24.1
julia> ϵ_optic = 4.5
julia> phonon_freqs = [4.016471586720514, 3.887605410774121, 3.5313112232401513, 2.755392921480459, 2.4380741812443247, 2.2490917637719408, 2.079632190634424, 2.0336707697261187, 1.5673011873879714, 1.0188379384951798, 1.0022960504442775, 0.9970130778462072, 0.9201781906386209, 0.800604081794174, 0.5738689505255512]
julia> ir_activities = [0.08168931020200264, 0.006311654262282101, 0.05353548710183397, 0.021303020776321225, 0.23162784335484837, 0.2622203718355982, 0.23382298607799906, 0.0623239656843172, 0.0367465760261409, 0.0126328938653956, 0.006817361620021601, 0.0103757951973341, 0.01095811116040592, 0.0016830270365341532, 0.00646428491253749]

julia> MAPIe = material(ϵ_optic, ϵ_static, m_eff, phonon_freqs, ir_activities, volume)

Hellwarth B Scheme... (athermal)
Hellwarth (58) summation: 0.03803673767058733
Hellwarth (59) summation (total ir activity ^2): 0.19283002835623678
Hellwarth (59) W_e (total ir activity ): 0.4391241605243747
Hellwarth (61) Omega (freq): 2.251571287857919
------------------------------------------
           Material Information           
------------------------------------------
Optic dielectric   | ϵo = 4.5
Static dielectric  | ϵs = 24.1
Ionic dielectric   | ϵi = [0.299962, 0.0247382, 0.254308, 0.166213, 2.30827, 3.07073, 3.20261, 0.892655, 0.886139, 0.720913, 0.401989, 0.618315, 0.766623, 0.155541, 1.16275]
Band mass          | mb = 0.12
Fröhlich coupling  | α = [0.0340093, 0.0028509, 0.03075, 0.0227523, 0.335905, 0.465256, 0.50462, 0.142232, 0.160835, 0.162287, 0.0912368, 0.140706, 0.181593, 0.0394994, 0.348765]
Phonon frequencies | f = [4.01647, 3.88761, 3.53131, 2.75539, 2.43807, 2.24909, 2.07963, 2.03367, 1.5673, 1.01884, 1.0023, 0.997013, 0.920178, 0.800604, 0.573869]
Eff Phonon freq    | feff = 2.25157
IR activities      | ir = [0.0816893, 0.00631165, 0.0535355, 0.021303, 0.231628, 0.26222, 0.233823, 0.062324, 0.0367466, 0.0126329, 0.00681736, 0.0103758, 0.0109581, 0.00168303, 0.00646428]
Unit cell volume   | V = 2.48858e-28
-------------------------------------------
```

We can then put this into the `frohlichpolaron` function as before.
For example here we would get:
```
julia> MAPIe_polaron = frohlichpolaron(MAPIe, 300; verbose = true)

-----------------------------------------------------------------------
                         Polaron Information:                          
-----------------------------------------------------------------------
Phonon frequencies             | ωeff = 2.25157 | ω = [4.01647, 3.88761, 3.53131, 2.75539, 2.43807, 2.24909, 2.07963, 2.03367, 1.5673, 1.01884, 1.0023, 0.997013, 0.920178, 0.800604, 0.573869]
Fröhlich coupling              | αeff = 2.6633 | α = [0.0340093, 0.0028509, 0.03075, 0.0227523, 0.335905, 0.465256, 0.50462, 0.142232, 0.160835, 0.162287, 0.0912368, 0.140706, 0.181593, 0.0394994, 0.348765]
Number of spatial dimensions   | d = 3
Small α→0 energy               | Fs = -2.75087
Large α→∞ energy               | Fl = -3.58205
Small α→0 fictitious mass      | Ms = 0.621212
Large α→∞ fictitious mass      | Ml = 0.102027
Small α→0 polaron radius       | Rs = 1.592
Large α→∞ polaron radius       | Rl = 10.0138
Large α→∞ FC peak freq.        | ΩFC = 1.00348
-----------------------------------------------------------------------
                     Zero Temperature Information:                     
-----------------------------------------------------------------------
Variational parameter          | v0 = 3.29227
Variational parameter          | w0 = 2.67919
Energy                         | E0 = -4.71904
Electron energy                | A0 = -1.83134
Interaction energy             | B0 = 4.88955
Trial energy                   | C0 = 1.66083
Fictitious spring constant     | κ0 = 3.66096
Fictitious mass                | M0 = 0.510019
Fictitious mass (asymptotic)   | M0a = 1.22883
Reduced mass                   | M0r = 0.337757
Polaron radius                 | R0 = 0.858447
-----------------------------------------------------------------------
                    Finite Temperature Information:                    
-----------------------------------------------------------------------
Temperatures                   | T = 300
Reduced thermodynamic          | β = [0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975, 0.159975]
Variational parameter          | v = 35.1921
Variational parameter          | w = 32.4542
Free energy                    | F = -10.3578
Electron energy                | A = -11.5987
Interaction energy             | B = 15.4728
Trial energy                   | C = 6.48363
Fictitious spring constant     | κ = 185.207
Fictitious mass                | M = 0.175839
Fictitious mass (asymptotic)   | Ma = 1.08436
Reduced mass                   | Mr = 0.149544
Polaron radius                 | R = 0.0554786
-----------------------------------------------------------------------
                      DC Mobility Information:                         
-----------------------------------------------------------------------
Finite temperature mobility    | μ = 0.572962
FHIP low-temp. mobility        | μFHIP = 4.71425
Devreese low-temp. mobility    | μD = 0.0673372
Kadanoff low-temp. mobility    | μK = 0.0673372
Hellwarth mobility             | μH = 0.0360092
Hellwarth mobility (b=0)       | μH0 = 0.0359907
Kadanoff relaxation time       | τ = 0.0572674
-----------------------------------------------------------------------
```

Using Unitful we then get the following energies and mobility:
```
julia> MAPIe_polaron.F0 |> u"meV"
-19.516365044323575 meV

julia> MAPIe_polaron.F |> u"meV"
-42.83622939936868 meV

julia> MAPIe_polaron.μ |> u"cm^2/V/s"
160.38638000423998 cm² s⁻¹ V⁻¹
```

## Multiple variational parameters

The base trial action used in the Feynman-Jensen variational path integral approximation has only two variational parameters called `v` and `w`.
This can be generalised to more variational parameters `v1`, `w1`, `v2`, `w2` etc where `v1 > w1 > v2 > w2 > ...`.
Including more variational paramaters can produce better variational results, improving the upper-bound on the polaron free energy.

To use this functionality, you have to provide an initial vector of guesses for the variational parameters.
For example for the single effective phonon mode MAPIe:
```
julia> MAPIe_polaron = frohlichpolaron(MAPIe; v_guesses = [4, 3], w_guesses = [3.5, 2.5], verbose = true)

-----------------------------------------------------------------------
               Polaron Information: [1 / 1 (100.0 %)]
-----------------------------------------------------------------------
Phonon frequencies             | ωeff = 2.25 | ω = 2.25
Fröhlich coupling              | αeff = 2.39394 | α = 2.3939410167951287
Small α→0 energy               | Fs = -2.46469
Large α→∞ energy               | Fl = -3.43751
Small α→0 fictitious mass      | Ms = 0.542264
Large α→∞ fictitious mass      | Ml = 0.0666023
Small α→0 polaron radius       | Rs = 1.67917
Large α→∞ polaron radius       | Rl = 9.00108
Large α→∞ FC peak freq.        | ΩFC = 0.810764
Number of dimensions [1 / 1]   | d = 3
-----------------------------------------------------------------------
                   Zero Temperature Information:                       
-----------------------------------------------------------------------
Variational parameter          | v0 = [2.28596, 12.4981]
Variational parameter          | w0 = [1.86037, 12.1284]
Energy                         | E0 = -5.57292
Electron energy                | A0 = -2.684
Interaction energy             | B0 = 5.78926
Trial energy                   | C0 = 2.46766
Fictitious spring constant     | κ0 = [1.76463, 9.1037]
Fictitious mass                | M0 = [0.509867, 0.0618885]
Fictitious mass (asymptotic)   | M0a = [1.22877, 1.03048]
Reduced mass                   | M0r = [0.33769, 0.0582815]
Polaron radius                 | R0 = [1.48402, 0.672612]
-----------------------------------------------------------------------
```

which gives a (polaron units) ground-state energy of `-5.57292` compared to the previous prediction `-5.56947`.

Note that whilst the code is capable of finite temperatures, more variational parameters and multiple modes all at once - the resultant energy landscape is difficult to optimise on and find converged results and requires well informed initial guesses to arrive at a correct solution.

## The Holstein polaron model

To calculate properties for a small Holstein-like polaron, use the `holstenpolaron()` function. It accepts all the same arguements as for the previous case with Frohlich, but will use the Holstein model instead and store the information in the `Holstein` type.

For example:
```
julia> holstein_polaron = holsteinpolaron(1:3, [0, 300], [0, 1]; dims = [1, 2, 3], verbose = true)

-----------------------------------------------------------------------
                         Polaron Information:                          
-----------------------------------------------------------------------
Phonon frequencies             | ωeff = 1 | ω = 1
Holstein coupling              | αeff = [1, 2, 3] | α = [1, 2, 3]
Number of spatial dimensions   | d = 3
-----------------------------------------------------------------------
                     Zero Temperature Information:                     
-----------------------------------------------------------------------
Variational parameter          | v0 = [4.06018, 4.80557, 7.07573]
Variational parameter          | w0 = [3.30811, 1.51507, 1.16092]
Total energy                   | F0 = [-7.3982, -9.21291, -12.1558]
Electron energy                | A0 = [-1.1281, -4.93576, -8.87222]
Interaction energy             | B0 = [1.50268, 4.90273, 9.86404]
Trial energy                   | C0 = [1.02362, 3.24593, 5.16394]
Fictitious spring constant     | κ0 = [5.54143, 20.7981, 48.7182]
Fictitious mass                | M0 = [0.506361, 9.06064, 36.1484]
Fictitious mass (asymptotic)   | M0a = [1.22734, 3.17185, 6.09495]
Reduced mass                   | M0r = [0.336149, 0.900603, 0.973081]
Polaron radius                 | R0 = [0.629813, 0.182562, 0.0945703]
-----------------------------------------------------------------------
                    Finite Temperature Information:                    
-----------------------------------------------------------------------
Temperatures                   | T = [0, 300]
Reduced thermodynamic          | β = [Inf, 0.00333333]
Variational parameter          | v = [4.06018 4.80557 7.07573; 6.76201 8.98156 10.7774]
Variational parameter          | w = [3.30811 1.51507 1.16092; 3.75187 3.22804 2.85246]
Free energy                    | F = [-7.3982 -9.21291 -12.1558; -1751.64 -1757.63 -1763.61]
Electron energy                | A = [-1.1281 -4.93576 -8.87222; 1739.64 1739.63 1739.61]
Interaction energy             | B = [1.50268 4.90273 9.86404; 5.98483 11.9697 17.9545]
Trial energy                   | C = [1.02362 3.24593 5.16394; 0.0131867 0.0292696 0.0450057]
Fictitious spring constant     | κ = [5.54143 20.7981 48.7182; 31.6483 70.2482 108.016]
Fictitious mass                | M = [0.506361 9.06064 36.1484; 2.2483 6.7415 13.2754]
Fictitious mass (asymptotic)   | Ma = [1.22734 3.17185 6.09495; 1.80231 2.78236 3.77828]
Reduced mass                   | Mr = [0.336149 0.900603 0.973081; 0.692147 0.870826 0.929949]
Polaron radius                 | R = [0.629813 0.182562 0.0945703; 0.142314 0.0738927 0.0526417]
-----------------------------------------------------------------------
                      DC Mobility Information:                         
-----------------------------------------------------------------------
Finite temperature mobility    | μ = [Inf Inf Inf; 440.224 199.198 95.5733]
-----------------------------------------------------------------------
                  Frequency Response Information:                      
-----------------------------------------------------------------------
Electric field frequency       | Ω = [0, 1]
Memory function                | χ = [Inf+0.0im Inf+0.0im; -2.90365+5.53403im -3.85738e-5+0.00226794im;;; Inf+0.0im Inf+0.0im; -185.241+292.499im -0.000229113+0.00496805im;;; Inf+0.0im Inf+0.0im; -767.896+2281.86im -0.00127186+0.0101249im]
Complex impedance              | z = [0.0+Inf*im 0.0+Inf*im; 5.53403+1.90365im 0.00226794-0.999961im;;; 0.0+Inf*im 0.0+Inf*im; 292.499+184.241im 0.00496805-0.999771im;;; 0.0+Inf*im 0.0+Inf*im; 2281.86+766.896im 0.0101249-0.998728im]
Complex conductivity           | σ = [0.0+0.0im 0.0+0.0im; 0.16158-0.0555821im 0.0022681+1.00003im;;; 0.0+0.0im 0.0+0.0im; 0.00244768-0.00154176im 0.00497021+1.0002im;;; 0.0+0.0im 0.0+0.0im; 0.000393763-0.000132337im 0.0101496+1.00117im]
-----------------------------------------------------------------------
```

Using `addunits!(holstein_polaron)` will then add the specific Holstien polaron units (energy in terms of the electron hopping energy instead of phonon energy as for the Frohlich model) which can then be transformed into other SI units using Unitful.

## General electron-phonon matrices and k-space integration

TBC

## Saving and loading polaron data

The `Frohlich` and `Holstein` polaron types can be saved as `.jld` files and loaded back into Julia.
For example:

```
julia> save_frohlich_polaron(MAPIe_polaron, "MAPIe_polaron")
Saving polaron data to MAPIe_polaron.jld ...
... Polaron data saved.
```

```
julia> MAPIe_polaron = load_frohlich_polaron("MAPIe_polaron.jld")
Loading polaron data from MAPIe_polaron.jld ...
... Polaron loaded.
-----------------------------------------------------------------------
                         Polaron Information:                          
-----------------------------------------------------------------------
Phonon frequencies             | ωeff = 2.25 | ω = 2.25
Fröhlich coupling              | αeff = 2.39394 | α = 2.39394
Number of spatial dimensions   | d = 3
Small α→0 energy               | Fs = -2.46469
Large α→∞ energy               | Fl = -3.43751
Small α→0 fictitious mass      | Ms = 0.542264
Large α→∞ fictitious mass      | Ml = 0.0666023
Small α→0 polaron radius       | Rs = 1.67917
Large α→∞ polaron radius       | Rl = 9.00108
Large α→∞ FC peak freq.        | ΩFC = 0.810764
-----------------------------------------------------------------------
                     Zero Temperature Information:                     
-----------------------------------------------------------------------
Variational parameter          | v0 = 3.30877
Variational parameter          | w0 = 2.66327
Energy                         | E0 = -5.56947
Electron energy                | A0 = -2.17859
Interaction energy             | B0 = 5.78198
Trial energy                   | C0 = 1.96608
Fictitious spring constant     | κ0 = 3.855
Fictitious mass                | M0 = 0.543495
Fictitious mass (asymptotic)   | M0a = 1.24237
Reduced mass                   | M0r = 0.35212
Polaron radius                 | R0 = 0.817278
-----------------------------------------------------------------------
                    Finite Temperature Information:                    
-----------------------------------------------------------------------
Temperatures                   | T = 300
Reduced thermodynamic          | β = 0.159975
Variational parameter          | v = 19.8612
Variational parameter          | w = 16.9599
Free energy                    | F = -8.59191
Electron energy                | A = -14.5096
Interaction energy             | B = 16.5498
Trial energy                   | C = 6.55175
Fictitious spring constant     | κ = 106.831
Fictitious mass                | M = 0.371408
Fictitious mass (asymptotic)   | Ma = 1.17107
Reduced mass                   | Mr = 0.270822
Polaron radius                 | R = 0.0722548
-----------------------------------------------------------------------
                      DC Mobility Information:                         
-----------------------------------------------------------------------
Finite temperature mobility    | μ = 0.487357
FHIP low-temp. mobility        | μFHIP = 2.93118
Devreese low-temp. mobility    | μD = 0.703373
Kadanoff low-temp. mobility    | μK = 0.703373
Hellwarth mobility             | μH = 1.09655
Hellwarth mobility (b=0)       | μH0 = 1.09824
Kadanoff relaxation time       | τ = 0.964611
-----------------------------------------------------------------------
                  Frequency Response Information:                      
-----------------------------------------------------------------------
Electric field frequency       | Ω = 0
Memory function                | χ = Inf+0.0im
Complex impedance              | z = 0.0+Inf*im
Complex conductivity           | σ = 0.0+0.0im
-----------------------------------------------------------------------
```

To do the same for a Holstein polaron use `save_holstein_polaron` and `load_holstein_polaron`.

Further details in the
[documentation](https://Frost-group.github.io/PolaronMobility.jl/).

## Research outputs

The central output of this model are temperature-dependent polaron mobilities: 

![MAPI Polaron mobility, plotted vs experimental data](mobility-calculated-experimental.png)

From the variational solution, you have characterised the polarons in your
system. 
This gives access to the effective mass renormalisations (phonon drag), polaron
binding energies, effective electron-phonon coupling parameters, etc.

## Community guidelines

Contributions to the code (extending that which is calculated), or additional
physical systems / examples, are very welcome. 

If you have questions about the software, scientific questions, or find errors,
please create a [GitHub issue](https://github.com/Frost-group/PolaronMobility.jl/issues). 

## Reference

If you find this package (or snippets, such as the entered and tested
free-energy expressions) useful for your work, please cite the paper 
[Frost2017PRB](https://doi.org/10.1103/PhysRevB.96.195202). 

```
@article{Frost2017,
  doi = {10.1103/physrevb.96.195202},
  url = {https://doi.org/10.1103/physrevb.96.195202},
  year  = {2017},
  month = {nov},
  publisher = {American Physical Society ({APS})},
  volume = {96},
  number = {19},
  author = {Jarvist Moore Frost},
  title = {Calculating polaron mobility in halide perovskites},
  journal = {Physical Review B}
}
```

These codes use the `Optim.jl` optimisation library to do the essential calculation of the Feynman variational theory. 
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00615/status.svg)](https://doi.org/10.21105/joss.00615)


# Scientific Discussion

## Variational Structure

The package implements the Feynman-Jensen variational idea in a model-agnostic
form. After any exactly solvable degrees of freedom have been integrated out,
the partition function is written as a path integral over a carrier path
measure. A tractable trial process with action `S0` gives

```math
\mathcal F
\le
\mathcal F_0
+
\frac{1}{\beta}
\left\langle
S-S_0
\right\rangle_0 .
```

In code the objective is decomposed as

```julia
free_energy(trial, parameters, beta) +
entropy_cost(trial, parameters, beta) +
interaction_free_energy(model, trial, parameters, beta)
```

This is the central abstraction. The solver only manipulates parameters and
bounds; the physics lives in methods attached to each model/trial pair.

## Fröhlich Hamiltonian And Coupling

The continuum Fröhlich model describes a carrier coupled to longitudinal
optical phonons through the macroscopic polar dielectric field. In the
single-mode convention used by Feynman, the material enters through

```math
\alpha =
\frac{1}{2}
\frac{1}{4\pi\epsilon_0}
\left(
\frac{1}{\epsilon_\infty}
-
\frac{1}{\epsilon_s}
\right)
\frac{e^2}{\hbar\omega}
\sqrt{\frac{2m_b\omega}{\hbar}} .
```

For multimode materials, `FrohlichMaterial` stores mode-resolved couplings
`\alpha_j` and reduced phonon frequencies `\omega_j/\omega_{\mathrm{eff}}`.
The Hellwarth B effective frequency is used for Kelvin, THz, mobility, and
length conversions, while the variational interaction itself still sums over
the individual modes.

## Feynman Gaussian Process

Feynman's 1955 insight was to integrate out the phonons exactly and replace
the resulting retarded self-interacting electron path by a solvable Gaussian
process: the electron is harmonically coupled to a fictitious coordinate. The
public parameters are

```math
w,\delta,
\qquad
v = w + \delta .
```

At zero temperature the Feynman mean-square-displacement kernel is

```math
D_F(\tau)
=
\frac{w^2}{v^2}\tau
+
\frac{v^2-w^2}{v^3}\left(1-e^{-v\tau}\right).
```

The trial contribution is

```math
A = -\frac{d}{2}(v-w),
\qquad
C = \frac{d}{4v}(v^2-w^2),
```

and the interaction contribution is evaluated as

```math
B =
\sum_j
\alpha_j K_d(\omega_j)
\int_0^\infty
\frac{e^{-\omega_j\tau}}
{\sqrt{\omega_j D_F(\tau)}}\,d\tau .
```

The minimized zero-temperature energy is stored as

```math
E_{\mathrm{var}} = -(A+B+C).
```

At finite temperature, `D_F` is replaced by the corresponding periodic
Brownian-bridge kernel and the phonon propagator becomes thermal. These are
the Osaka/Hellwarth finite-temperature formulas used in the original
`PolaronMobility.jl` workflow.

## Multiple Gaussian Modes

`MultiGaussianTrial` is the finite-mode extension used for Martin-Frost-style
multimode variational work. The parameter order is

```math
(w_1,\delta_1,w_2,\delta_2,\ldots),
\qquad
v_i = w_i+\delta_i .
```

The zero-temperature displacement is

```math
D(\tau)
=
\tau
+
\sum_i
\frac{h_i}{v_i^2}
\left[
\frac{1-e^{-v_i\tau}}{v_i}
-
\tau
\right],
```

with algebraic weights `h_i` obtained from the coupled fictitious oscillator
spectrum. For one mode this reduces exactly to Feynman's kernel, and the test
suite enforces that reduction.

For more than one mode, the variational energy should be no higher than the
one-mode value when the optimizer finds the same branch. Small improvements
are expected; dramatic drops at weak or intermediate coupling should be
treated as a warning sign.

## Profile-Functional Gaussian Trial

`ProfileGaussianTrial` is the package's current implementation of the older
general quadratic-memory functional idea associated with Adamowski, Gerlach,
and Leschke. It represents a positive profile

```math
\Gamma(\omega)
=
\sum_i
a_i
\frac{\nu_i^2}{\omega^2+\nu_i^2},
\qquad
a_i \ge 0 .
```

The profile defines the Gaussian displacement through

```math
D_\Gamma(\tau)
=
\frac{2}{\pi}
\int_0^\infty
\frac{1-\cos(\omega\tau)}
{\omega^2[1+\Gamma(\omega)]}\,d\omega ,
```

and the zero-temperature entropy/log-determinant contribution through

```math
E_{\mathrm{det}}
=
\frac{d}{2\pi}
\int_0^\infty
\left[
\log(1+\Gamma(\omega))
-
\frac{\Gamma(\omega)}{1+\Gamma(\omega)}
\right]d\omega .
```

For numerical stability, the package does not evaluate the displacement
integral by direct quadrature. Since `Γ` is a rational positive profile, it
decomposes

```math
\frac{1}{1+\Gamma(\omega)}
=
1 + \sum_j \frac{c_j}{\omega^2+p_j^2},
```

so the displacement is evaluated analytically as

```math
D_\Gamma(\tau)
=
\tau
+
\sum_j
\frac{c_j}{p_j^2}
\left[
\tau
-
\frac{1-e^{-p_j\tau}}{p_j}
\right] .
```

At finite temperature the same decomposition is used with the periodic
oscillator kernel

```math
O_p(\tau,\beta)
=
\frac{
1+e^{-p\beta}-e^{-p\tau}-e^{p(\tau-\beta)}
}
{p(1-e^{-p\beta})}
```

and the Brownian bridge `\tau(1-\tau/\beta)`.

Feynman's one-mode trial is embedded exactly by choosing

```math
\Gamma_F(\omega)
=
\frac{v^2-w^2}{\omega^2+w^2}.
```

In the package basis this means `basis_frequencies = [w]` and
`a1 = (v^2-w^2)/w^2`. The test suite checks that this profile reproduces the
Feynman displacement, entropy correction, interaction, and total objective.

## Experimental Nonlocal Gaussian Kernel

`NonlocalGaussianTrial` is deliberately more provisional than
`ProfileGaussianTrial`. It uses a finite nonlocal amplitude basis for the
displacement kernel and a configurable quadratic regularizer as its entropy
cost. Because that regularizer is not the Gaussian profile log determinant,
the returned objective is a useful kernel-exploration score but not a
validated variational upper bound.

If a nonlocal calculation returns an apparently much lower "energy" than the
Feynman or profile result, do not interpret that as a physical variational
improvement. It means the trial is outside the currently validated bound
machinery. Use `ProfileGaussianTrial` for energy comparisons.

## Fröhlich Mobility And Frequency Response

The one-mode Feynman solution feeds the historical mobility formulas:

- FHIP low-temperature mobility from Feynman, Hellwarth, Iddings, and
  Platzman.
- Kadanoff/Devreese low-temperature variants.
- Hellwarth finite-temperature mobility and multimode effective-frequency
  handling.

The package also evaluates a frequency-dependent memory function
`\Sigma(\Omega)`. It stores

```math
Z(\Omega) = -i[\Omega+\Sigma(\Omega)],
\qquad
\sigma(\Omega) = Z(\Omega)^{-1}.
```

The scalar FHIP, Kadanoff, and Hellwarth mobility formulas assume the
two-parameter Feynman action. For profile and nonlocal Gaussian trials those
literature-reference fields are reported as `NaN`, while the direct
memory-function response is still evaluated from the trial displacement kernel.

## Lattice Models As Poisson CTMC Variational Problems

The lattice models use a different solvable path measure from the continuum Fröhlich Gaussian trials. The electron lives on a hypercubic lattice and the trial path is a symmetric continuous-time nearest-neighbor Markov chain with hopping rate `\kappa`. The same `PoissonTrial` supports Holstein site coupling, Peierls bond coupling, and compatible composites.

The free CTMC contribution and relative-entropy rate for changing the bare hopping `J` to the variational rate `\kappa` are

```math
F_0(\kappa)=-2d\kappa,
\qquad
S_{\rm rel}(\kappa)=2d\kappa\log\frac{\kappa}{J}.
```

The hypercubic CTMC propagator factorizes:

```math
P_{\mathbf r}^{(d)}(t)=
e^{-2d\kappa t}\prod_{\mu=1}^d I_{r_\mu}(2\kappa t).
```

The package names the return and nearest-neighbor components

```math
q_0(t)=\left[e^{-2\kappa t}I_0(2\kappa t)\right]^d,
```

```math
q_1(t)=e^{-2d\kappa t}I_1(2\kappa t)I_0(2\kappa t)^{d-1}.
```

In code, `lattice_q0` and `lattice_q1` evaluate these with exponentially scaled Bessel functions, so the formulas are stable for `d = 1, 2, 3` and large arguments.

## Holstein CTMC Free Energy

The reduced Holstein Hamiltonian is

```math
H_H=
-J\sum_{\langle ij\rangle}(c_i^\dagger c_j+c_j^\dagger c_i)
+\omega_H\sum_i b_i^\dagger b_i
+g\sum_i n_i(b_i+b_i^\dagger).
```

After integrating out the local Einstein phonons, the lattice path acquires a retarded attraction between times at which it occupies the same site. At finite inverse temperature the package uses

```math
D_\beta(u;\omega)=
\frac{e^{-\omega u}+e^{-\omega(\beta-u)}}{1-e^{-\beta\omega}},
```

and the periodic CTMC bridge

```math
C_{\rm site}(u;\beta)=
\frac{q_0(u)q_0(\beta-u)}{q_0(\beta)}.
```

The variational objective is

```math
F_{\rm var}^{H}(\kappa,\beta)=
-2d\kappa+2d\kappa\log\frac{\kappa}{J}
-
\frac{g^2}{2}
\int_0^\beta D_\beta(u;\omega_H)C_{\rm site}(u;\beta)\,du.
```

At zero temperature this becomes

```math
F_{\rm var}^{H}(\kappa,\infty)=
-2d\kappa+2d\kappa\log\frac{\kappa}{J}
-
g^2\int_0^\infty e^{-\omega_Ht}\operatorname{ive}_0(2\kappa t)^d\,dt.
```

For `d = 1`, the integral reduces to

```math
\int_0^\infty e^{-\omega_Ht}e^{-2\kappa t}I_0(2\kappa t)\,dt
=
\frac{1}{\sqrt{\omega_H(\omega_H+4\kappa)}}.
```

## Peierls Bond-Coupled Lattice Model

The standalone Peierls model moves the electron-phonon coupling from site density to bond hopping modulation:

```math
H_P=
-J\sum_{\langle ij\rangle}B_{ij}
+\omega_P\sum_b a_b^\dagger a_b
+g_P\sum_b B_b(a_b+a_b^\dagger),
\qquad
B_{ij}=c_i^\dagger c_j+c_j^\dagger c_i.
```

`PeierlsMaterial` follows the same boundary convention as `HolsteinMaterial`: a relaxation energy and phonon energy define the reduced coupling,

```math
g_P=\sqrt{\frac{E_P}{h\nu_P}}.
```

The Peierls influence depends on bond-order correlations. The periodic bond bridge used in the package is

```math
C_{\rm bond}^{+}(u;\beta)=
2d\frac{q_0(u)q_0(\beta-u)+q_1(u)q_1(\beta-u)}{q_0(\beta)}.
```

The finite-temperature Peierls interaction term is

```math
F_P(\kappa,\beta)=
-\frac{g_P^2}{2}
\int_0^\beta D_\beta(u;\omega_P)C_{\rm bond}^{+}(u;\beta)\,du.
```

At zero temperature,

```math
F_P(\kappa,\infty)=
-\frac{g_P^2}{2}
\int_0^\infty e^{-\omega_Pt}
2d[\operatorname{ive}_0(2\kappa t)+\operatorname{ive}_1(2\kappa t)]
\operatorname{ive}_0(2\kappa t)^{d-1}\,dt.
```

The factor `1/2` is part of the Gaussian influence prefactor. The helper `peierls_integral_d` returns the full `2d[...]` bond integral; `interaction_free_energy(::PeierlsModel, ...)` multiplies it by `-g_P^2/2`.

## Lattice Mobility And Frequency Response

The lattice response implemented in the package is a CTMC first-return
current-blip projection. It is analogous in spirit to building response on top
of an optimized solvable trial process, but it is not the continuum FHIP
memory-function construction. Any older lattice "memory-function" wording in
the project should therefore be read historically rather than as the active
runtime theory.

Conceptually, the lattice calculation has two pieces:

1. the variational free energy chooses the effective hopping rate `\kappa`
   from periodic CTMC bridge correlations; and
2. the transport calculation inserts a current blip into that optimized CTMC
   path measure and dresses it with the exact phonon factor appropriate to the
   model.

The electron dynamics is encoded in the first-return kernel, while the phonons
enter through sideband weights multiplying that kernel.

Current insertions select hop-like path segments. The duration distribution is approximated by the CTMC first-return kernel from a nearest-neighbor site back to the origin:

```math
\widehat f_{a\to0}^{(d)}(s)=\frac{G_a^{(d)}(s)}{G_0^{(d)}(s)}.
```

Using the lattice resolvent identity, the package evaluates this as

```math
\widehat f_{a\to0}^{(d)}(s)=
\frac{s+2d\kappa-[G_0^{(d)}(s)]^{-1}}{2d\kappa},
```

where

```math
G_0^{(d)}(s)=\int_0^\infty e^{-st}\operatorname{ive}_0(2\kappa t)^d\,dt,
\qquad \Re s>0.
```

For `d = 1`, this gives the analytic expression

```math
\widehat f_{1\to0}^{(1)}(s)=
\frac{s+2\kappa-\sqrt{s(s+4\kappa)}}{2\kappa}.
```

For `d = 2`, the return Green's function is evaluated through the
square-lattice elliptic-integral form

```math
G_0^{(2)}(s)=
\frac{2}{\pi(s+4\kappa)}
K\!\left[\left(\frac{4\kappa}{s+4\kappa}\right)^2\right].
```

For `d = 3`, the package uses the cubic-lattice reduction

```math
G_0^{(3)}(s)=
\frac{1}{\pi}\int_0^\pi
G_0^{(2)}\!\left(s+2\kappa(1-\cos q)\right)\,dq.
```

So the standard `d = 2, 3` Holstein and Peierls transport calculations do not
repeatedly integrate the oscillatory `0..∞` Bessel/Laplace kernel. The small
positive broadening `\epsilon` fixes the retarded sheet through
`s = \epsilon - i\nu`.

The Einstein scale is retained as a diagnostic,

```math
\mu_E=\beta\kappa.
```

Given transport sidebands `(\nu_m,w_m)`, the reduced kernel is

```math
K(\Omega)=
\sum_m w_m\widehat f_{a\to0}^{(d)}(\epsilon-i(\Omega+\nu_m)).
```

The package reports

```math
\mu(\Omega)=\mu_EK(\Omega),
\qquad
\sigma(\Omega)=\mu(\Omega),
\qquad
Z(\Omega)=\sigma(\Omega)^{-1}.
```

At zero frequency,

```math
\mu_{\rm DC}=\beta\kappa\operatorname{Re}K(0).
```

At zero temperature this Einstein prefactor diverges, so the package treats DC
and optical response separately. The DC mobility remains infinite, but for
finite frequency the reported reduced optical response is

```math
\mu(\Omega)=\sigma(\Omega)=K(\Omega), \qquad \Omega \neq 0,\; T=0,
```

with `Z(\Omega)=K(\Omega)^{-1}`. This is the natural reduced current-current
normalization of the CTMC first-return/blip ansatz and keeps the Holstein and
Peierls optical spectra finite.

Thus the lattice mobility is not obtained by first defining a separate
force-memory function and then inverting a Drude-like impedance. Instead, the
package evaluates the retarded CTMC first-return kernel directly at the
phonon-shifted frequencies generated by the blip cloud or current vertex.

## Holstein, Peierls, And Composite Transport Sidebands

Holstein transport uses the exact current-blip cloud factor

```math
C_H(t)=
\exp[-S_H C_\beta(1-\cos\omega_H t)]
\exp[iS_H\sin\omega_Ht],
\qquad
S_H=2(g/\omega_H)^2,
```

where `C_\beta = coth(\beta\omega_H/2)`. At zero temperature this is a
Poisson ladder. At finite temperature the package evaluates the equivalent
Bessel sideband weights, so the runtime response is a direct exact-cloud
sideband sum over the first-return kernel rather than a separate memory
closure.

Peierls transport is different. Since the Peierls phonon modulates the current vertex itself, the normalized vertex cloud is

```math
\mathcal C_P(t)=
\frac{J_{\rm resp}^2+g_P^2D_P^>(t)}{J_{\rm resp}^2+g_P^2D_P^>(0)}.
```

Here `J_resp` is selected by `response_hopping`, with default `:bare`. At finite temperature,

```math
D_P^>(t)=(N_P+1)e^{-i\omega_Pt}+N_Pe^{i\omega_Pt}.
```

Thus Peierls has one zero-phonon current channel and `\pm\omega_P`
phonon-assisted channels. It is not an exponentiated Franck-Condon ladder
unless the hopping depends nonlinearly on displacement. In the code, this
means the Peierls response is assembled from a small exact sideband list
attached directly to the CTMC first-return kernel; it is not a fallback
renewal approximation and not a lattice memory-function ansatz.

For independent Holstein and Peierls baths, the composite transport cloud is
the product

```math
C_{HP}(t)=C_H(t)\mathcal C_P(t),
```

so the sideband list is the convolution of Holstein and Peierls sidebands.
This is the response-side analogue of adding independent influence functionals
in the free energy. If `\omega_H \neq \omega_P`, the implementation keeps the
actual shifted frequencies rather than forcing the composite transport onto one
integer phonon ladder. Shared-mode Holstein-Peierls cross terms are not
included.

## Combining Lattice Influence Functionals

Holstein and Peierls models live on the same lattice path space when they use the same dimension, hopping scale, and `PoissonTrial`. The composite model therefore adds independent influence functionals on the same path measure:

```math
F_{\rm int}^{\rm composite}(\kappa,\beta)=\sum_mF_{\rm int}^{(m)}(\kappa,\beta).
```

The corresponding code is:

```julia
holstein = holstein_poisson_problem(hopping = 1.0, phonon_frequency = 1.0, coupling = 0.4)
peierls = peierls_poisson_problem(hopping = 1.0, phonon_frequency = 0.8, coupling = 0.35)

problem = VariationalProblem(
    combine_models(holstein.model, peierls.model),
    PoissonTrial(dimension = 1, bare_hopping = 1.0),
)
result = solve(problem; temperatures = [1.0], frequencies = [0.0, 0.5])
```

For material-derived rubrene calculations, `rubrene_holstein_peierls_problem` builds both Ordejón Table II material models and reduces them to a shared Holstein phonon reference frequency before composition.

For transport-focused studies, the package also exposes guide-style helpers
that can reuse a frozen zero-temperature optimum:

```julia
rows = holstein_peierls_transport_sweep(
    problem;
    temperatures = 0.05:0.05:0.5,
    frequencies = [0.0, 0.25, 0.5],
    kappa_source = :zero_temperature,
)
```

These rows expose `mobility`, `mobility_einstein`, `mobility_factor`,
`conductivity_real`, `conductivity_imag`, sideband normalization diagnostics,
and whether the rate came from `:zero_temperature` or
`:per_temperature`.

## Limiting Regimes

Useful analytical limits are encoded in examples and tests because they are
excellent guardrails for future model additions.

For Fröhlich polarons, weak coupling gives the perturbative energy
`E \sim -\alpha` and the Feynman parameters remain close to the bare oscillator
scale. In the package convention this is often summarized as `v,w -> 3ω` for
the weak-coupling initial branch. At strong coupling, the optimized fictitious
oscillator softens toward the phonon scale, `w -> \omega`, while

```math
v \sim \frac{4\alpha^2}{9\pi}.
```

Hellwarth finite-temperature mobility decreases with thermal phonon
occupation. Adiabatic versus antiadiabatic behavior is read through the
relative carrier and phonon time scales: slow phonons produce long memory,
while fast phonons produce short retarded kernels.

For Holstein polarons, zero coupling removes the interaction and the
Donsker-Varadhan objective is minimized at the bare rate,

```math
\kappa \to J .
```

Weak coupling is obtained by evaluating the retarded interaction at the bare
rate. In the strong antiadiabatic/localized limit, the leading energy is the
Lang-Firsov shift

```math
-E_p = -\frac{g^2}{\omega_0},
```

with exponentially narrowed hopping, schematically

```math
-E_p - 2dJ\exp[-(g/\omega_0)^2].
```

High temperature shrinks imaginary-time memory through
`\coth(\beta\omega_0/2)` and the finite `\beta` integration window. The
adiabatic limit has `J/\omega_0 >> 1`; the antiadiabatic limit has
`J/\omega_0 << 1`.

For Peierls polarons, zero Peierls coupling recovers the bare Poisson walk or,
in a composite model, the remaining component. Weak Peierls corrections scale
as `-g_P^2`. Antiadiabatic Peierls phonons make the bond memory short-ranged;
adiabatic Peierls phonons create long-lived bond modulation, the variational
analogue of dynamic disorder.

## Tables, Sweeps, And Plots

`solve` returns typed result structs. Table helpers flatten those into
`NamedTuple` rows:

```julia
solution_table(result)
mobility_table(result)
response_table(result)
```

Continuation sweeps warm-start forward and backward branches and select the
lower-free-energy row. Frequency sweeps optimize once per model/coupling and
temperature point, then evaluate the response on the requested frequency grid
without reoptimizing per frequency.

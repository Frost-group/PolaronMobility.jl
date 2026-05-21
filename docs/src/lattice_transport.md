# Lattice Transport Theory

This page records the active lattice implementation used by `HolsteinModel`,
`PeierlsModel`, and compatible `CompositePolaronModel` calculations.

There are two related pieces:

1. a CTMC/Feynman-Jensen variational free energy built from periodic bridge
   correlations; and
2. a CTMC first-return transport projection dressed by Holstein and Peierls
   phonon sidebands.

The same optimized Poisson hopping rate `κ` can be used for both, but the two
approximations are conceptually distinct. The free energy is a normalized
imaginary-time variational calculation. The transport calculation is a
current-vertex/blip approximation using the first-return kernel of the CTMC.

## CTMC Trial Process

The lattice trial is a symmetric continuous-time nearest-neighbor Markov chain
with hopping rate `κ` on a `d`-dimensional hypercubic lattice. The free trial
contribution and relative-entropy mismatch are

```math
F_0(\kappa)=-2d\kappa,
\qquad
S_{\mathrm{rel}}(\kappa)=2d\kappa\log\frac{\kappa}{J}.
```

The CTMC propagator factorizes as

```math
P_{\mathbf r}^{(d)}(t)=
e^{-2d\kappa t}\prod_{\mu=1}^{d}I_{r_\mu}(2\kappa t).
```

The package uses

```math
q_0(t)=\left[e^{-2\kappa t}I_0(2\kappa t)\right]^d,
```

```math
q_1(t)=e^{-2d\kappa t}I_1(2\kappa t)I_0(2\kappa t)^{d-1}.
```

These are implemented as `lattice_q0` and `lattice_q1` using exponentially
scaled Bessel functions.

## Periodic Bridge Free Energies

For finite inverse temperature `β`, the Einstein phonon kernel is

```math
D_\beta(u;\omega)=
\frac{e^{-\omega u}+e^{-\omega(\beta-u)}}{1-e^{-\beta\omega}}.
```

Holstein uses the site bridge

```math
C_{\mathrm{site}}(u;\beta)=
\frac{q_0(u)q_0(\beta-u)}{q_0(\beta)}.
```

The finite-temperature Holstein objective is

```math
F_{\rm var}^{H}(\kappa,\beta)=
-2d\kappa
+2d\kappa\log\frac{\kappa}{J}
-
\frac{g^2}{2}
\int_0^\beta D_\beta(u;\omega_H)C_{\rm site}(u;\beta)\,du.
```

Peierls uses the bond-order bridge

```math
C_{\mathrm{bond}}^{+}(u;\beta)=
2d\frac{q_0(u)q_0(\beta-u)+q_1(u)q_1(\beta-u)}{q_0(\beta)}.
```

The finite-temperature Peierls interaction contribution is

```math
F_P(\kappa,\beta)=
-\frac{g_P^2}{2}
\int_0^\beta D_\beta(u;\omega_P)C_{\rm bond}^{+}(u;\beta)\,du.
```

At zero temperature, the package exposes Laguerre-integral kernels

```math
I_H^{(d)}(\kappa)=
\int_0^\infty dt\,e^{-\omega_Ht}\operatorname{ive}_0(2\kappa t)^d,
```

```math
I_P^{(d)}(\kappa)=
\int_0^\infty dt\,e^{-\omega_Pt}
2d[\operatorname{ive}_0(2\kappa t)+\operatorname{ive}_1(2\kappa t)]
\operatorname{ive}_0(2\kappa t)^{d-1}.
```

Thus

```math
F_H(\kappa,\infty)=-g^2I_H^{(d)}(\kappa),
\qquad
F_P(\kappa,\infty)=-\frac{g_P^2}{2}I_P^{(d)}(\kappa).
```

For a composite Holstein-Peierls model, independent influence functionals are
summed on the same Poisson path measure.

## General-d First-Return Kernel

Transport uses the CTMC first-return kernel from a nearest neighbor back to the
origin. Let

```math
G_0^{(d)}(s)=
\int_0^\infty dt\,e^{-st}\operatorname{ive}_0(2\kappa t)^d,
\qquad \Re s>0.
```

Then the lattice resolvent identity gives

```math
\widehat f_{a\to0}^{(d)}(s)=
\frac{s+2d\kappa-G_0^{(d)}(s)^{-1}}{2d\kappa}.
```

For `d = 1`, this reduces to

```math
\widehat f_{1\to0}^{(1)}(s)=
\frac{s+2\kappa-\sqrt{s(s+4\kappa)}}{2\kappa}.
```

For `d = 2`, the package uses the square-lattice elliptic-integral Green's
function

```math
G_0^{(2)}(s)=
\frac{2}{\pi(s+4\kappa)}
K\!\left[\left(\frac{4\kappa}{s+4\kappa}\right)^2\right].
```

For `d = 3`, it uses the cubic-lattice reduction

```math
G_0^{(3)}(s)=
\frac{1}{\pi}\int_0^\pi
G_0^{(2)}\!\left(s+2\kappa(1-\cos q)\right)\,dq.
```

Only dimensions beyond the standard hypercubic `d = 1, 2, 3` path fall back
to the generic Bessel/Laplace quadrature. The retarded response uses
`s = broadening - im * frequency_shift` with positive `broadening`.

## Holstein, Peierls, And Composite Sidebands

Holstein transport uses the exact-cloud current-blip factor

```math
C_H(t)=
\exp[-S_H C_\beta(1-\cos\omega_Ht)]
\exp[iS_H\sin\omega_Ht],
\qquad
S_H=2(g/\omega_H)^2.
```

At `T = 0`, this is a Poisson sideband ladder. At finite temperature, the
sidebands are equivalent to the usual Bessel-weighted Holstein cloud expansion.

Peierls transport is a current-vertex factor, not a Franck-Condon ladder:

```math
\mathcal C_P(t)=
\frac{J_{\rm resp}^2+g_P^2D_P^>(t)}{J_{\rm resp}^2+g_P^2D_P^>(0)}.
```

At finite temperature,

```math
D_P^>(t)=(N_P+1)e^{-i\omega_Pt}+N_Pe^{i\omega_Pt}.
```

So Peierls contributes a zero-phonon current channel and `±ω_P` assisted
channels. `J_resp` comes from `response_hopping`, whose default is the bare
model hopping.

For independent Holstein and Peierls phonons, the composite cloud is

```math
C_{HP}(t)=C_H(t)\mathcal C_P(t),
```

so the sidebands are a convolution. Shared-mode cross terms are not included.

## Mobility, Conductivity, And Impedance

The reduced Einstein scale is

```math
\mu_E=\beta\kappa.
```

Given model sidebands `(ν_m, w_m)`, the reduced transport kernel is

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

For `T = 0`, the DC mobility still diverges, but the finite-frequency optical
response is reported in the reduced current-current normalization,

```math
\mu(\Omega)=\sigma(\Omega)=K(\Omega), \qquad \Omega \neq 0,\; T=0,
```

so the Holstein and Peierls optical spectra remain finite while
`mobility_factor` continues to store the same transport kernel.

`LatticeMobilityResult` stores the DC value `mobility`, the Einstein estimate
`mobility_einstein`, the reduced kernel `mobility_factor`, and diagnostics such
as sideband count, sideband normalization, broadening, Laguerre points, and
rate provenance. The same metadata is carried by `LatticeResponseResult`.

## Guide-Style Sweeps

Besides the main `solve(...)` API, the package exposes

```julia
lattice_transport_sweep(problem; temperatures, frequencies, kappa_source = :zero_temperature)
holstein_transport_sweep(...)
peierls_transport_sweep(...)
holstein_peierls_transport_sweep(...)
```

These helpers are intended for theory-driven mobility studies:

- `kappa_source = :zero_temperature` reuses one frozen `β = ∞` optimized rate,
- `kappa_source = :per_temperature` reoptimizes at every requested temperature,
- sidebands are built over the requested temperature grid,
- repeated complex first-return evaluations are cached over reused frequency
  shifts,
- returned rows include `mobility`, `mobility_factor`, `conductivity_real`,
  `conductivity_imag`, `conductivity_abs`, and transport diagnostics.

A typical DC mobility sweep is

```julia
rows = holstein_peierls_transport_sweep(
    holstein_coupling = 1.5,
    peierls_coupling = 0.5,
    dimension = 3,
    temperatures = 0.1:0.1:2.0,
    frequencies = [0.0],
    kappa_source = :zero_temperature,
)
```

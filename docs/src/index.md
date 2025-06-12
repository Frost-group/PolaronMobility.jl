# PolaronMobility.jl

These codes calculate the temperature-dependent polaron mobility for
a material. 
We use the Feynman path-integral variational approach. 
We have parameters for various metal-halide Perovskites; as well as other
systems. 
These codes implement methods described across a wide range of now quite old
literature. 
The methods have been tested against literature values, and the calculated
units more well understood. 
They enable relative 'turn key' calculation of polaron parameters (most
particularly the finite temperature charge carrier mobility) for an
arbitrary material system, based on parameters that are standard to calculate
with modern ab-initio electronic structure methods. 

This documentation is intended to be read alongside the paper which was the
first application (and motivation) for these codes: 
[Frost2017PRB](https://doi.org/10.1103/PhysRevB.96.195202)
([ArXiv:1704.05404](https://arxiv.org/abs/1704.05404) ).

The required inputs are the dielectric constants (ϵ-static and ϵ-optic)
, a characteristic phonon frequency (ω), and the bare-electron band
effective-mass (me).  These values can be relatively easily calculated in the
ab-initio electronic structure package of your choosing, or measured directly.

From these four values, the code solves a temperature-dependent polaron model.
This is done by variationally optimising the temperature-dependent
free-energies for the coupled
electron-phonon system.
These optimised parameters describe the polaron with the infinite quantum field
of lattice vibrations 'integrated through', and replaced with a phonon-drag
term.
From this the polaron features such as effective-mass, size of the
wavefunction, frequency of energy oscillation etc. can be calculated.

This polaron state can then be used as an input to further models for polaron
mobility.

The codes are designed to produce a set of temperature-dependent mobilities and
other data, for direct incorporation into a scientific publication.

May your phonons drag in a manner truly sublime.

## Installation

These codes require Julia >1.7 . 

To install, type the following at the Julia REPL:

```
julia> Pkg.add("PolaronMobility")
```

## Community guidelines

Contributions to the code (extending that which is calculated), or additional
physical systems / examples, are very welcome. 

If you have questions about the software, scientific questions, or find errors,
please create a [GitHub issue](https://github.com/Frost-group/PolaronMobility.jl/issues). 

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


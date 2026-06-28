# QMParticleCircle

Monte Carlo study of the topological properties of a quantum particle 
constrained on a circle in the presence of a confined magnetic field.

## Theoretical Background

A magnetic field confined inside the circle, although not directly coupled 
to the particle, induces a topological term in the action proportional to 
the winding number Q of the trajectory, which weights the different 
homotopy classes differently in the path integral. 
The topological susceptibility,

  χ(β) = (1/β) ⟨Q²⟩ |θ=0,

quantifies the sensitivity of the free energy to this topology and is 
estimated here from the winding number sampled on a 1D Euclidean lattice 
with periodic boundary conditions.

## Contents

- `windingMC.c` — lattice Monte Carlo simulation and winding number measurement
- `observable.c` — Jackknife analysis with optional reweighting
- `lib/` — random number generator utilities

## Method

- Metropolis updates on a 1D periodic lattice
- Parallel Tempering across replicas at different lattice spacings to 
  overcome topological freezing near the continuum limit
- In the high temperature regime a reweighting potential is added to the 
  action to enhance sampling of rare topological sectors. The 
  correct statistical weight is restored  in post processing analysis.
- Jackknife error estimation with Kahan summation for numerical stability

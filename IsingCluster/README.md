# IsingCluster

Monte Carlo study of the critical properties of the 2D Ising model, with 
critical exponents extracted via finite-size scaling.

## Theoretical Background

Near a second order phase transition the magnetic susceptibility χ and 
correlation length ξ diverge as χ ~ |t|^(-γ), ξ ~ |t|^(-ν), with 
t = (β-βc)/βc. On a finite lattice these singularities are regularized 
by the lattice size L and the critical exponents are instead extracted 
from the L dependence of the observables through finite-size-scaling 
ansatz such as

  χ(β,L) = L^(γ/ν) χ₂[(β-βc)L^(1/ν)] + corr.

## Contents

- `isingCluster.c` — Wolff cluster Monte Carlo simulation
- `observable.c` — Jackknife analysis of susceptibility and Binder cumulant
- `lib/` — random number generator utilities

## Method

- Wolff cluster updates to mitigate critical slowing down near βc
- Magnetic susceptibility χ' = βL^D(⟨m²⟩-⟨|m|⟩²) and Binder cumulant 
  U = ⟨m⁴⟩/⟨m²⟩² computed from the sampled magnetization
- Jackknife error estimation with Kahan summation for numerical stability

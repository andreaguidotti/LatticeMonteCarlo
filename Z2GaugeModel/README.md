# Z2GaugeModel3D

Monte Carlo study of the static potential and string tension in the 
confined and deconfined phases of the 3D Z2 lattice gauge theory.

## Theoretical Background

The Z2 lattice gauge theory exhibits both a confined and a 
deconfined phase. Since gauge symmetry is local the 
Elitzur theorem forbids the use of local order parameters to distinguish 
the two phases. The static potential between two static sources,

  V(R) = − lim(T→∞) (1/T) log⟨W(R,T)⟩,

extracted from the large time behavior of the Wilson loop W(R,T), is used 
instead: a linear potential (area law) signals confinement, a constant 
potential at large distance (perimeter law) signals deconfinement.

## Contents

- `src/z2model.c` — lattice Monte Carlo simulation and Wilson loop measurement
- `src/observables.c` — Jackknife analysis of Wilson loop data
- `lib/` — random number generator and lattice geometry utilities

## Method

- Ordered Metropolis updates with precomputed Boltzmann factors
- Multihit variance reduction via the exact closed form link expectation, 
  tanh(β·sumStaples)
- For each configuration, Wilson loops of all extensions are measured 
  at all lattice sites and spatial directions and averaged, providing a 
  variance-reduced measurement per sweep and limiting simulation time. 
  This introduces correlations between loops of different extensions, 
  neglected in the (uncorrelated) fits used in the analysis
- Jackknife error estimation with Kahan summation for numerical stability

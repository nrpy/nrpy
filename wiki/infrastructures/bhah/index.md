# BHaH

> Router for standalone BHaH generated-application infrastructure.

| Page | Go here when... |
| --- | --- |
| [Project Assembly, main.c, And Makefile Emission](lifecycle-and-project-assembly.md) | You need executable/library generation order, `main.c` assembly, `Makefile_helpers.py`, or generated-output boundaries. |
| [Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md) | You need generated structs, parameters, `BHaH_defines.h`, parfiles, CLI overrides, or CUDA host/device structs. |
| [Grids, Coordinates, And Reference Metrics](grids-coordinates-and-reference-metrics.md) | You need numerical grid allocation, `ReferenceMetric`/fisheye setup, coordinate conversions, or CFL `dt`/`dsmin` evaluation. |
| [Curvilinear Boundaries And Stencils](curvilinear-boundaries-and-stencils.md) | You need `CurviBoundaryConditions`, inner parity/outer radiation BCs, angular grid-count parity rules, or stencil ghost bounds. |
| [MoL Time Integration](mol-time-integration.md) | You need Method of Lines storage, Runge-Kutta method tables, stage stepping, or time-integration hooks. |
| [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md) | You need diagnostics dispatch, 0D/1D/2D/volume output, progress indicators, or checkpoint/restart readers and writers. |
| [GR Evolution And Constraints](gr-evolution-and-constraints.md) | You need BHaH BSSN/fCCZ4 RHS, Ricci, YBS/CAHD options, algebraic repairs, `H`/`M^2`/`Lambda^2` constraints, or tile helpers. |
| [GR Initial Data And Psi4](gr-initial-data-and-psi4.md) | You need Cartesian/spherical ADM, TwoPunctures (`W` lapse), NRPyElliptic, TOVola initial data, basis transforms, or `Psi4` decomposition. |
| [SO(3) Rotation Helpers](so3-rotation-helpers.md) | You need BHaH generated C helper registration, conventions, or string validation under `nrpy/infrastructures/BHaH/rotation/`. |
| [BHaHAHA](bhahaha/index.md) | You need the BHaHAHA apparent-horizon C API, spherical grid interpolation, or multiresolution hyperbolic relaxation. |
| [GRoovy GRHD Runtime](groovy-grhd-runtime.md) | You need GRoovy GRHD gridfunctions, GRHayL/EOS calls, reconstruction, HLL fluxes, or primitive recovery. |
| [NRPyElliptic Workflow](nrpyelliptic-workflow.md) | You need NRPyElliptic initial guesses, hyperbolic-relaxation RHSs, residual norms, stop conditions, or diagnostics. |
| [Geodesics And Raytracing Runtime](geodesics-and-raytracing-runtime.md) | You need geodesic executables, evolution raytracing export, RKF45 integration, Lagrange interpolation, or output files. |
| [SEOBNR BOB Generated Library](seobnr-bob-generated-library.md) | You need SEOBNRv5/BOB/SEBOBv2 generated C library assembly, `commondata`, GSL ODE integration, or strain resampling. |

Back to [Infrastructures](../index.md).

# How to Cite NRPy

If NRPy contributes to published scientific results, please cite the papers below that are relevant to your use of NRPy.

> **Not sure which papers apply?** Cite the **general NRPy paper** plus each paper relevant to the numerical methods and generated projects you used.

## Always cite

- **General NRPy / curvilinear-coordinate numerical relativity foundation**  
  Ian Ruchlin, Zachariah B. Etienne, and Thomas W. Baumgarte,  
  ["SENR/NRPy+: Numerical relativity in singular curvilinear coordinate systems"](https://doi.org/10.1103/PhysRevD.97.064036),  
  *Physical Review D* **97**, 064036 (2018). ([arXiv](https://arxiv.org/abs/1712.07658))

## Quick reference

| If you used... | Also cite... |
|---|---|
| `nrpy.examples.nrpyelliptic_conformally_flat` | Assumpção et al. (2022) |
| GPU / CUDA-enabled elliptic projects | Tootle et al. (2025) |
| `superB` / Charm++ distributed-memory projects | Jadoo et al. (2025) |
| `nrpy.examples.bhahaha` or `BHaHAHA` | Etienne et al. (2026) |
| OpenMP builds of `nrpy.examples.blackhole_spectroscopy`, `nrpy.examples.spinning_blackhole`, or `nrpy.examples.two_blackholes_collide` with BHaHAHA horizon finding | Etienne et al. (2026) |
| `nrpy.examples.groovy_TOV_BSSN` or other `GRoovy` GRHD projects | Pierre Jacques et al. (2025) |
| `nrpy.examples.seobnrv5_aligned_spin_inspiral -seobnrv5_nrpy` | Mahesh et al. (2025) |
| `nrpy.examples.sebobv2` or `nrpy.examples.sebobv1_jax` | Mahesh et al. (2025) |
| `nrpy.examples.seobnrv5_aligned_spin_inspiral -seobnrv5_bob` | Mahesh et al. (2025) |
| `nrpy.examples.seobnrv5_aligned_spin_inspiral -seobnrv5_nrnqc_bob` | Mahesh et al. (2025) |

## NRPy project citations

| NRPy Project | Citation |
|---|---|
| **`NRPyElliptic` projects** | Thiago Assumpção, Leonardo R. Werneck, Terrence Pierre Jacques, and Zachariah B. Etienne, ["Fast hyperbolic relaxation elliptic solver for numerical relativity: Conformally flat, binary puncture initial data"](https://doi.org/10.1103/PhysRevD.105.104037), *Physical Review D* **105**, 104037 (2022). ([arXiv](https://arxiv.org/abs/2111.02424)) |
| **GPU / CUDA-enabled projects** | Samuel D. Tootle, Leonardo R. Werneck, Thiago Assumpção, Terrence Pierre Jacques, and Zachariah B. Etienne, ["Accelerating numerical relativity with code generation: CUDA-enabled hyperbolic relaxation"](https://doi.org/10.1088/1361-6382/add63e), *Classical and Quantum Gravity* **42**, 115005 (2025). ([arXiv](https://arxiv.org/abs/2501.14030)) |
| **`superB` / Charm++ distributed-memory projects** | Nishita Jadoo, Terrence Pierre Jacques, and Zachariah B. Etienne, ["superB/NRPy: Scalable, Task-Based Numerical Relativity for 3G Gravitational Wave Science"](https://doi.org/10.1088/1361-6382/adee71), *Classical and Quantum Gravity* **42**, 155006 (2025). ([arXiv](https://arxiv.org/abs/2505.00097)) |
| **`BHaHAHA` apparent-horizon finding** | Zachariah B. Etienne, Thiago Assumpção, Leonardo Rosa Werneck, and Samuel D. Tootle, ["BHaHAHA: A Fast, Robust Apparent Horizon Finder Library for Numerical Relativity"](https://doi.org/10.1088/1361-6382/ae09e9), *Classical and Quantum Gravity* **43**, 075007 (2026). ([arXiv](https://arxiv.org/abs/2505.15912)) |
| **`GRoovy` general relativistic hydrodynamics** | Terrence Pierre Jacques, Samuel Cupp, Leonardo R. Werneck, Samuel D. Tootle, Maria C. Babiuc Hamilton, and Zachariah B. Etienne, ["General relativistic hydrodynamics code for dynamical spacetimes with curvilinear coordinates, tabulated equations of state, and neutrino physics"](https://doi.org/10.1103/hc9l-1thx), *Physical Review D* **112**, 084044 (2025). ([arXiv](https://arxiv.org/abs/2412.03659)) |
| **NRPy `SEOBNR` / `SEBOB` waveform projects** | Siddharth Mahesh, Sean T. McWilliams, and Zachariah Etienne, ["Combining effective one-body inspirals and backwards one-body merger-ringdowns for aligned spin black hole binaries"](https://doi.org/10.1103/67c8-1rv3), *Physical Review D* **112**, 124065 (2025). ([arXiv](https://arxiv.org/abs/2508.20418)) |

## Minimal default recommendation

- Cite the **general NRPy paper**.
- Add **each project paper relevant to the numerical methods and generated projects you used**; a calculation may use more than one project, such as GRoovy evolution with BHaHAHA horizon finding.

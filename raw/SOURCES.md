# Source Manifest

> Pointer manifest for KB source material. Preserved historical documents live under `raw/source-docs/` so `AGENTS.md` is the sole root KB router; code, configuration, tests, and build inputs stay in place. Status is `living` for mutable sources, `frozen` for version-pinned specifications, and `frozen-historical` for superseded documents preserved verbatim.

## Aggregate Sources

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `core-top-level-package-modules` | Top-level package files matching `nrpy/*.py`, `nrpy/*.txt`, and `nrpy/py.typed`. | living | partial |
| `helpers-validation-and-reference-metric-tests` | Helper modules, validation utilities, and reference-metric tests under `nrpy/helpers`, `nrpy/validate_expressions`, and `nrpy/tests`. | living | partial |
| `equation-modules-and-trusted-values` | Equation modules, generated trusted-value dictionaries, and BOB fit metadata under `nrpy/equations`. | living | partial |
| `infrastructure-modules-and-embedded-headers` | Infrastructure Python modules, embedded C/C++/CUDA headers, C/C++ test drivers, generated C/CUDA trusted files, and package notes under `nrpy/infrastructures`. | living | partial |
| `example-generators-and-companion-scripts` | Example generators and companion scripts under `nrpy/examples`. | living | partial |
| `ci-and-local-automation` | CI workflows and local static-analysis/validation scripts under `.github`. | living | partial |

## Source Documents Moved Below Root

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `raw/source-docs/original-agents.md` | Pre-KB root instructions preserved verbatim; `coding_style.md` and `AGENTS.md` govern current rules. | frozen-historical | ingested |
| `raw/source-docs/kb-instructions.md` | Bootstrap KB template moved from repository root; `wiki/SCHEMA.md` and `wiki/workflows.md` govern current rules. | frozen-historical | ingested |
| `raw/source-docs/ybs-momentum-damping-spec.md` | Superseded YBS-MOM auxiliary-relaxation proposal preserved verbatim; `BSSN_RHSs.py` governs current behavior. | frozen-historical | ingested |

## Cited Code And Config Sources

Add an exact local path here only when no aggregate row above covers it.

| Source | Status |
| --- | --- |
| `AGENTS.md` | living |
| `README.md` | living |
| `CITATION.md` | living |
| `setup.py` | living |
| `requirements.txt` | living |
| `requirements-dev.txt` | living |
| `pyproject.toml` | living |
| `.mypy.ini` | living |
| `.pydocstyle` | living |
| `.pylintrc` | living |
| `.pylintrc_python36` | living |
| `.darglint` | living |
| `bin/nrpyinline.py` | living |
| `coding_style.md` | living |
| `tools/kb_lint.py` | living |

## External Sources

| Source | Tier | Provenance | Status | Ingest | Notes |
| --- | --- | --- | --- | --- | --- |
| `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md` | background | Andrej Karpathy `LLM Wiki` gist note. | frozen | partial | Source, compiled-page, and format-rule separation. |
| `https://github.com/paralab/Dendro-5.01` | external-spec | Dendrolib upstream host library repository (URL prefix). | living | partial | Cited octree, block, zip/unzip, remesh, AH, and CMake interfaces. |
| `https://github.com/paralab/Dendro-GR` | external-spec | Dendro-GR upstream application repository (URL prefix). | living | partial | Cited block padding, AMR, KO, puncture, checkpoint, and Psi4 conventions. |
| `https://arxiv.org/pdf/gr-qc/0104063v3` | external-spec | Baker, Campanelli, and Lousto, arXiv:gr-qc/0104063v3. | frozen | ingested | Psi4 tetrad in Sec. V.A, Eqs. (5.6)-(5.9). |
| `https://doi.org/10.1103/PhysRevD.65.044001` | background | Baker, Campanelli, and Lousto, Phys. Rev. D 65, 044001 (2002). | frozen | ingested | Published journal metadata for arXiv:gr-qc/0104063v3. |
| `https://arxiv.org/pdf/astro-ph/0503420v2` | external-spec | Duez et al., arXiv:astro-ph/0503420v2. | frozen | ingested | GRMHD four-magnetic field, fluxes, sources, and HLL speeds. |
| `https://arxiv.org/pdf/2412.03659v2` | external-spec | Jacques et al., arXiv:2412.03659v2. | frozen | ingested | Reference-metric GRHD conserved variables, fluxes, sources, and rescaling. |
| `https://arxiv.org/pdf/0902.3652v2` | external-spec | Brown, arXiv:0902.3652v2. | frozen | ingested | Conformal connection constraint and RHS adjustment, Eqs. (12), (15), (21)-(22). |
| `https://doi.org/10.1103/PhysRevD.79.104029` | background | Brown, Phys. Rev. D 79, 104029 (2009). | frozen | ingested | Published journal metadata for arXiv:0902.3652v2. |
| `https://arxiv.org/pdf/gr-qc/0209066v2` | external-spec | Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066v2. | frozen | ingested | Gamma-constraint adjustment coefficient in Eq. (45). |
| `https://arxiv.org/pdf/1205.5111v2` | external-spec | Yo, Lin, and Cao, arXiv:1205.5111v2. | frozen | ingested | `xi` parameterization in Eq. (47) and STF momentum gradient in Eq. (56). |
| `https://doi.org/10.1103/PhysRevD.110.064045` | external-spec | Etienne, Phys. Rev. D 110, 064045 (2024). | frozen | ingested | CAHD CFL/local grid-spacing scaling pattern. |
| `https://arxiv.org/pdf/1106.2254v2` | external-spec | Alic et al., arXiv:1106.2254v2. | frozen | ingested | Covariant CCZ4 (`kappa3=1`), connection stretch, and Eq. (20) 1+log lapse. |
| `https://arxiv.org/pdf/1403.3653v1` | external-spec | Sanchis-Gual et al., arXiv:1403.3653v1. | frozen | ingested | Reference-metric fCCZ4 variables, constraints, gauge, and Eq. (3.14). |
| `https://arxiv.org/pdf/2002.06225v2` | external-spec | Mewes et al., arXiv:2002.06225v2. | frozen | ingested | Reference-metric fCCZ4 system, Ricci, gauge, and matter projections. |
| `https://arxiv.org/pdf/2201.08857v1` | external-spec | Baumgarte and de Oliveira, arXiv:2201.08857v1. | frozen | ingested | Bona-Masso lapse equation and `f(alpha)=1` harmonic slicing. |
| `https://arxiv.org/pdf/gr-qc/9810065v1` | external-spec | Baumgarte and Shapiro, arXiv:gr-qc/9810065v1. | frozen | ingested | Zero-shift harmonic lapse relation, Eqs. (30)-(32). |
| `https://arxiv.org/pdf/gr-qc/9902024v1` | external-spec | Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024v1. | frozen | ingested | Zero-shift harmonic slicing `partial_t(alpha)=partial_t(exp(6*phi))`. |
| `https://arxiv.org/pdf/1712.07658v2` | external-spec | Ruchlin, Etienne, and Baumgarte, arXiv:1712.07658v2. | frozen | ingested | StaticTrumpet gauge pairing in Eqs. (65), (67), and (69). |
| `https://arxiv.org/pdf/2303.18143v2` | external-spec | Khalil et al., arXiv:2303.18143v2. | frozen | ingested | Precessing spin-evolution Eqs. (65)-(66), (69), and (71)-(72). |
| `https://github.com/zachetienne/nrpytutorial` | background | Historical NRPy+ tutorial notebooks (URL prefix). | living | partial | Computational-physics style, `HarmonicSlicing` name, and Step 2.b chain rule. |
| `https://arxiv.org/abs/1605.01938` | background | HBR2016 final-spin paper landing page. | living | partial | Final-spin fit background. |
| `https://arxiv.org/abs/1611.00332` | background | UIB2016 final-state paper landing page. | living | partial | Final-mass/spin fit background. |
| `https://arxiv.org/abs/2111.02424` | background | arXiv:2111.02424 landing page. | living | partial | Waveform model background. |
| `https://arxiv.org/abs/2303.18039` | background | SEOBNRv5HM aligned-spin paper landing page. | living | partial | Aligned-spin waveform background. |
| `https://arxiv.org/abs/2303.18046` | background | SEOBNRv5PHM precessing paper landing page. | living | partial | Precessing waveform background. |
| `https://arxiv.org/abs/2303.18143` | background | SEOBNRv5 precessing dynamics landing page. | living | partial | Precessing dynamics background. |
| `https://arxiv.org/abs/2508.20418` | background | SEBOB aligned-spin attachment landing page. | living | partial | SEBOB attachment background. |
| `https://arxiv.org/abs/2510.25012` | background | BOBv2 waveform paper landing page. | living | partial | BOBv2 model background. |
| `https://pypi.org/pypi/black/26.5.1/json` | external-spec | PyPI metadata for Black 26.5.1. | frozen | partial | Black Python `>=3.10` runtime requirement. |
| `https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html` | background | Intel Intrinsics Guide. | living | partial | SSE/AVX/AVX-512 intrinsic terminology. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/index.html` | background | NVIDIA CUDA Programming Guide. | living | partial | CUDA kernel/memory terminology. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html` | background | NVIDIA CUDA NVCC guide. | living | partial | `nvcc` compiler-driver terminology. |
| `https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case` | background | Einstein Toolkit test-case guide. | living | partial | ET test-suite conventions. |
| `https://docs.jax.dev/en/latest/installation.html` | background | JAX installation guide. | living | partial | JAX platform and install context. |
| `https://www.gnu.org/software/gsl/doc/html/usage.html` | background | GNU GSL 2.8 usage manual. | living | partial | GSL compilation and linking flags. |
| `https://link.aps.org/doi/10.1103/PhysRev.55.374` | background | Oppenheimer and Volkoff, Phys. Rev. 55, 374. | living | partial | TOV stellar-equilibrium background. |
| `https://web2.ph.utexas.edu/~gsudama/pub/1967_008.pdf` | background | Goldberg et al. spin-weighted spherical harmonics PDF. | living | partial | Spin-weighted harmonic formula background. |
| `https://pubs.aip.org/aip/jmp/article/57/9/092504/648118/How-should-spin-weighted-spherical-functions-be` | background | J. Math. Phys. 57, 092504 landing page. | living | partial | Quaternion spin-weighted function background. |
| `https://charm.readthedocs.io/en/v8.0.0/` | external-spec | Charm++ v8.0.0 quickstart, language, and libraries manuals (URL prefix). | frozen | partial | `.ci` files, chares, SDAG, PUP, reductions, `charmc`, `charmrun`, and CkIO. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/` | background | Charm++ v8.0.0 upstream compiler, `charmxi`, SDAG, and `ck-core` sources (URL prefix). | frozen | partial | Pinned Charm++ 8.0.0 chare-array, callback, location, and reduction code. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech9.html` | background | Cactus Users Guide Ch. C1 (`Application thorns`). | living | partial | ETLegacy thorn and CCL terminology. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech12.html` | background | Cactus Users Guide Ch. D2 (CCL reference). | living | partial | CarpetX CCL file terminology. |
| `https://einsteintoolkit.org/thornguide/CarpetX/CarpetX/documentation.html` | background | Einstein Toolkit `CarpetX` thorn guide. | living | partial | CarpetX scheduling and gridfunction tags. |

## Exclusions

- `project/`, `build/`, `dist/`, and tool caches (`__pycache__/`, `.mypy_cache/`, `.pytest_cache/`).
- Object files, executables, PDFs, archives, images, and non-Markdown outputs.
- Unregistered generated C/CUDA/Charm++/JAX projects and thorns.
- Scratch logs, transcripts, audit/scan/token/latest reports, and uncommissioned planning or ranking files.

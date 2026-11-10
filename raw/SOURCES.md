# Source Manifest

> Pointer manifest for source material. Root-level documentation sources live
> under `raw/source-docs/` so `AGENTS.md` is the only root KB document. Code,
> config, fixtures, selected logs, and build inputs stay in place. Status is
> `frozen` when a source is meant not to change and `living` when drift must
> trigger re-ingest.

## Aggregate Sources

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `core-top-level-package-modules` | Core top-level package modules from `find nrpy -maxdepth 1 -type f \( -name '*.py' -o -name '*.txt' -o -name 'py.typed' \)`. | living | partial |
| `helpers-package-modules` | Helper package files from `find nrpy/helpers -type f \( -name '*.py' -o -name '*.h' \)`. | living | ingested |
| `helpers-validation-and-reference-metric-tests` | Helpers, validation helpers, and reference metric tests from `nrpy/helpers`, `nrpy/validate_expressions`, and `nrpy/tests`. | living | partial |
| `equation-modules-and-trusted-values` | Equation modules, generated trusted-value files, and BOB test metadata from `nrpy/equations`, including `nrpy/equations/seobnr/tests/BOB_v2_fit_sim_list.md`. | living | partial |
| `infrastructure-modules-and-embedded-headers` | Python modules and embedded headers from `find nrpy/infrastructures -type f \( -name '*.py' -o -name '*.h' \)`. | living | partial |
| `carpetx-package-inventory` | CarpetX Python package inventory from `find nrpy/infrastructures/CarpetX -type f -name '*.py'`. | living | ingested |
| `example-generators-and-companion-scripts` | Example generators and companion scripts from `nrpy/examples` inventoried by `wiki/examples/example-generator-catalog.md`. | living | partial |
| `ci-and-local-automation` | CI and local automation files from `.github`. | living | partial |

## Source Documents Moved Below Root

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `raw/source-docs/original-agents.md` | Previous root agent instructions, preserved byte-for-byte before replacing `AGENTS.md`. Its root style-guide reference is superseded; see the source map row for current routing. | frozen | ingested |
| `raw/source-docs/kb-instructions.md` | KB schema and governance source moved from the repository root. | frozen | ingested |
| `raw/source-docs/ybs-momentum-damping-spec.md` | Superseded commissioned YBS-MOM auxiliary-relaxation proposal preserved as historical input. See CONTR-0007 and the current YBS-MOM page for governing behavior. | frozen | ingested |

## Cited Code And Config Sources

Add an exact cited-file row here only for a path that no aggregate row above
covers; `tools/kb_lint.py` `_source_registered` lists the aggregate path rules.
Rows here abbreviate to source and status per `wiki/SCHEMA.md`: the repository
path is the provenance, and any ingest state belongs in `wiki/source-map.md`.

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
| `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_vector_so3_apply_R_to_vector__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_RT_to_vector_so3_apply_RT_to_vector__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_tensorDD_so3_apply_R_to_tensorDD__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_axis_angle_to_R_so3_axis_angle_to_R__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_build_R_from_hats_so3_build_R_from_hats__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_find_nU_and_dphi_from_unit_vectors_so3_find_nU_and_dphi_from_unit_vectors__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_left_multiply_hats_with_R_so3_left_multiply_hats_with_R__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_relative_R_dst_from_src_so3_relative_R_dst_from_src__openmp.c` | living |
| `nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp` | living |
| `nrpy/infrastructures/Dendro/tests_infra/README.md` | living |

## External Sources

| Source | Provenance | Status | Ingest | Notes |
| --- | --- | --- | --- | --- |
| `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md` | Andrej Karpathy gist raw note, `LLM Wiki`; the full URL is its stable source locator. | frozen | partial | Background approach source for persistent LLM-maintained wiki governance: raw/wiki/schema layers and ingest/query/lint workflows. |
| `https://github.com/paralab/Dendro-5.01` | Dendrolib upstream repository used by Dendro-GR as its host-library dependency. | living | partial | Deciding source for cited balanced-octree input, octant refinement flags, regular-block decomposition, block axis spacings, zip/unzip, remeshing, intergrid transfer, the apparent-horizon finder checkpoint interface, and the CMake link targets of the default BLAS and LAPACK build. Only cited interfaces are ingested. |
| `https://github.com/paralab/Dendro-GR` | Dendro-GR upstream application repository and its Dendrolib integration. | living | partial | Source for cited host conventions that the generated solvers follow: padded-origin arithmetic, block reads, refinement modes, remesh/intergrid-transfer sequencing, constraint interfaces, derivative-order selection, high-order centered advection functions, KO selection, namespace declarations, GridInfo CSV output, mesh-level output-frequency scaling, the radial RIT `eta` prescription, lapse replacement, the puncture tracker and relative-position history, integrator selection and CFL spacing, checkpoint slot selection on restore, Psi4 extraction, and parameter defaults. Only cited interfaces are ingested. |
| `https://arxiv.org/pdf/gr-qc/0104063v3` | Version-pinned v3 PDF for Baker, Campanelli, and Lousto, arXiv:gr-qc/0104063. | frozen | ingested | Primary deciding source for the Psi4 tetrad contract in Sec. V.A: Eqs. (5.6)-(5.7), the following unnumbered Gram-Schmidt procedure, and the later Eq. (5.9) rotation. |
| `https://arxiv.org/pdf/astro-ph/0503420v2` | Version-pinned v2 paper for Duez et al., arXiv:astro-ph/0503420. | frozen | ingested | Magnetic four-vector, stress-energy, conserved fluid variables, fluxes, sources, HLL flux, and approximate MHD dispersion in Eqs. (16), (23)-(24), (31)-(39), and (48)-(50). NRPy reference-densitizes the Cartesian fluid equations. |
| `https://arxiv.org/pdf/2412.03659v2` | Version-pinned v2 paper for Jacques et al., arXiv:2412.03659. | frozen | ingested | Conformal metric split, reference-metric volume factor, GRHD conserved variables, fluxes, sources, and rescaling in Eqs. (2)-(5), (11), (13)-(24). Eq. (11) rewrites the divergence of a chosen current; the paper does not define or assert conservation of NRPy's entropy current. |
| `https://arxiv.org/pdf/0902.3652v2` | Version-pinned v2 PDF for Brown, arXiv:0902.3652. | frozen | ingested | Deciding source for the conformal connection constraint in Eqs. (12a), (12b), and (15), and its RHS adjustment in Eqs. (21e) and (22e); it does not prescribe the local norm diagnostic. |
| `https://arxiv.org/pdf/gr-qc/0209066v2` | Version-pinned v2 PDF for Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066. | frozen | ingested | Deciding source for the full Gamma-constraint adjustment coefficient in Eq. (45); current NRPy already owns the baseline `2/3` contribution and therefore adds only the incremental `chi` term. |
| `https://arxiv.org/pdf/1205.5111v2` | Version-pinned v2 PDF for Yo, Lin, and Cao, arXiv:1205.5111. | frozen | ingested | Deciding source for the Eq. (47) `xi` parameterization and Eq. (56) direct covariant STF momentum-gradient adjustment; comparison with the earlier full Gamma coefficient gives NRPy's incremental mapping `chi=2*xi/3`. |
| `https://doi.org/10.1103/PhysRevD.110.064045` | Published article for Etienne, Phys. Rev. D 110, 064045 (2024). | frozen | ingested | Deciding source for the CAHD level-local CFL/grid-spacing scaling pattern transferred to YBS-MOM; CAHD itself remains parabolic. |
| `https://arxiv.org/pdf/1106.2254v2` | Version-pinned v2 PDF for Alic et al., arXiv:1106.2254. | frozen | ingested | Corroborating source for covariant CCZ4, the `kappa3=1` covariance choice, the one-stretch Cartesian connection composition, and the exact advective 1+log formula in Eq. (20). |
| `https://arxiv.org/pdf/1403.3653v1` | Version-pinned v1 PDF for Sanchis-Gual et al., arXiv:1403.3653. | frozen | ingested | Corroborating source for reference-metric fCCZ4 variables, constraints, gauge equations, reduction to BSSN, and one connection stretch. Its general display retains a geometric divergence coefficient while spherical Eq. (3.14) contains the compatible Z4 divergence promotion, so the two are not an unqualified full off-constraint cross-check. |
| `https://arxiv.org/pdf/2002.06225v2` | Version-pinned v2 PDF for Mewes et al., arXiv:2002.06225. | frozen | ingested | Principal scientific specification for NRPy's three-dimensional reference-metric fCCZ4 variables, Ricci tensor, evolution equations, gauge, and matter projections in Eqs. (3)-(41), including the unnumbered evolution system after Eq. (32). Eqs. (5)-(7) and (10) also support the ADM metric and BSSN-to-ADM conversion used by GRHD/GRMHD. Its printed connection shift bracket contains an apparent convention error under its own full-Lie time-operator definition; NRPy documents and corrects the duplicated constraint-vector stretch. |
| `https://arxiv.org/pdf/2201.08857v1` | Version-pinned v1 PDF for Baumgarte and de Oliveira, arXiv:2201.08857. | frozen | ingested | Deciding source for the general Bona-Masso lapse equation and its `f(alpha)=1` harmonic specialization. |
| `https://arxiv.org/pdf/gr-qc/9810065v1` | Version-pinned v1 PDF for Baumgarte and Shapiro, arXiv:gr-qc/9810065. | frozen | ingested | Deciding source for the zero-shift harmonic lapse relation in Eqs. (30)-(32), including the special `C(x)=1` choice. |
| `https://arxiv.org/pdf/gr-qc/9902024v1` | Version-pinned v1 PDF for Baumgarte, Hughes, and Shapiro, arXiv:gr-qc/9902024. | frozen | ingested | Deciding source for the zero-shift statement that their harmonic slicing reduces to `partial_t(alpha)=partial_t(exp(6*phi))`. |
| `https://arxiv.org/pdf/1712.07658v2` | Version-pinned v2 PDF for Ruchlin, Etienne, and Baumgarte, arXiv:1712.07658. | frozen | ingested | Deciding source for the StaticTrumpet gauge pairing in Eqs. (65), (67), and (69). |
| `https://arxiv.org/pdf/2303.18143v2` | Version-pinned v2 PDF for Khalil et al., arXiv:2303.18143. | frozen | ingested | Corroborating source for the spin-evolution equations (65), (66), (69), (71), and (72). |
| `https://github.com/zachetienne/nrpytutorial` | Historical NRPy tutorial repository; selected overview, scalar-wave, and BSSN notebooks. | living | partial | Background for NRPy's computational-physics vocabulary and explanatory style, plus provenance for the `HarmonicSlicing` name and the `W`/`phi` chain-rule implementation in Step 2.b. |
| `https://arxiv.org/abs/1605.01938` | arXiv abstract page for HBR2016 final-spin paper. | living | partial | Background landing page; exact revision-to-fit mapping is outside this ingest. |
| `https://arxiv.org/abs/1611.00332` | arXiv abstract page for UIB2016 final-state paper and ancillary implementation. | living | partial | Background landing page; exact revision-to-fit/ancillary mapping is outside this ingest. |
| `https://arxiv.org/abs/2111.02424` | arXiv abstract page for arXiv:2111.02424. | living | partial | Background landing page; exact revision-to-claim mapping is outside this ingest. |
| `https://arxiv.org/abs/2303.18039` | arXiv abstract page for SEOBNRv5HM aligned-spin waveform model. | living | partial | Background landing page; exact revision-to-equation mapping is outside this ingest. |
| `https://arxiv.org/abs/2303.18046` | arXiv abstract page for SEOBNRv5PHM precessing waveform model. | living | partial | Background landing page; exact revision-to-equation mapping is outside this ingest. |
| `https://arxiv.org/abs/2303.18143` | arXiv abstract page for SEOBNRv5 precessing dynamics groundwork. | living | partial | Background landing page; exact revision-to-equation mapping is outside this ingest. |
| `https://arxiv.org/abs/2508.20418` | arXiv abstract page for SEBOB aligned-spin attachment model. | living | partial | Background landing page; exact revision-to-equation mapping is outside this ingest. |
| `https://arxiv.org/abs/2510.25012` | arXiv abstract page for BOBv2 waveform model. | living | partial | Background landing page; exact revision-to-equation/validation mapping is outside this ingest. |
| `https://pypi.org/pypi/black/26.5.1/json` | PyPI release metadata for Black 26.5.1. | frozen | partial | Official release metadata for Black's Python 3.10+ runtime requirement and installation surface. |
| `https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html` | Intel Intrinsics Guide landing page. | living | partial | Official SSE, AVX, AVX512, and intrinsic-family terminology for SIMD helper documentation; not authority for NRPy macro behavior. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/index.html` | NVIDIA CUDA Programming Guide. | living | partial | Official CUDA programming-model terminology for helper documentation; not authority for NRPy kernel behavior. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html` | NVIDIA CUDA Programming Guide, NVCC chapter. | living | partial | Official NVCC compiler-driver terminology used by CUDA example prerequisites. |
| `https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case` | Einstein Toolkit documentation, `Adding a test case`. | living | partial | Official test-suite terminology; local workflow configuration decides NRPy CI behavior. |
| `https://docs.jax.dev/en/latest/installation.html` | JAX installation documentation. | living | partial | Official supported-platform and installation context; local requirements and generator code decide NRPy behavior. |
| `https://www.gnu.org/software/gsl/doc/html/usage.html` | GNU GSL 2.8 documentation, `Using the Library`. | living | partial | Official compiling/linking context for GSL-dependent generated examples. |
| `https://link.aps.org/doi/10.1103/PhysRev.55.374` | DOI landing page for `10.1103/PhysRev.55.374`. | living | partial | Original Oppenheimer-Volkoff stellar-equilibrium background. |
| `https://web2.ph.utexas.edu/~gsudama/pub/1967_008.pdf` | PDF URL for Goldberg et al. spin-weighted spherical-harmonic formula reference. | living | partial | Goldberg-formula background for spin-weighted spherical harmonics. |
| `https://pubs.aip.org/aip/jmp/article/57/9/092504/648118/How-should-spin-weighted-spherical-functions-be` | Journal of Mathematical Physics article landing page. | living | partial | Background on spin-weighted functions and quaternion viewpoints. |
| `https://charm.readthedocs.io/en/v8.0.0/quickstart.html` | Charm++ quickstart from the Charm++ repository. | frozen | partial | External spec context for `.ci` files, generated `.decl.h` and `.def.h` files, `charmc`, and `charmrun`. |
| `https://charm.readthedocs.io/en/v8.0.0/charm%2B%2B/manual.html` | Charm++ language manual from the Charm++ repository. | frozen | partial | External spec context for chares, entry methods, proxies, SDAG, PUP, checkpoint/restart, reductions, chare arrays, and generic message-transport caveats. |
| `https://charm.readthedocs.io/en/v8.0.0/libraries/manual.html` | Charm++ and Converse libraries manual from the Charm++ repository. | frozen | partial | External spec context for CkIO and Charm++ library APIs. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/scripts/charmc` | Charm++ pinned Charm++ 8.0.0 compiler-driver script. | frozen | partial | Background implementation context for `charmc`; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Chare.C` | Charm++ pinned Charm++ 8.0.0 charmxi chare translation source. | frozen | partial | Background implementation context for chare and chare-array generation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Parameter.C` | Charm++ pinned Charm++ 8.0.0 charmxi parameter-marshalling source. | frozen | partial | Background implementation context for parameter-marshaled entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarray.h` | Charm++ pinned Charm++ 8.0.0 chare-array header. | frozen | partial | Background implementation context for chare arrays, proxies, and maps; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarrayindex.h` | Charm++ pinned Charm++ 8.0.0 chare-array index header. | frozen | partial | Background implementation context for `CkArrayIndex3D` and index objects; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckcallback.C` | Charm++ pinned Charm++ 8.0.0 callback implementation. | frozen | partial | Background implementation context for callback dispatch; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/CParsedFile.C` | Charm++ pinned Charm++ 8.0.0 SDAG parsed-file source. | frozen | partial | Background implementation context for SDAG translation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/CEntry.C` | Charm++ pinned Charm++ 8.0.0 SDAG entry translation source. | frozen | partial | Background implementation context for SDAG entry behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/constructs/When.C` | Charm++ pinned Charm++ 8.0.0 SDAG `when` construct source. | frozen | partial | Background implementation context for SDAG `when` matching; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/sdag.h` | Charm++ pinned Charm++ 8.0.0 SDAG runtime header. | frozen | partial | Background implementation context for SDAG runtime buffering; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/cklocation.h` | Charm++ pinned Charm++ 8.0.0 location-manager header. | frozen | partial | Background implementation context for location, migration, and placement caveats; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/cklocation.C` | Charm++ pinned Charm++ 8.0.0 location-manager implementation. | frozen | partial | Background implementation context for location-manager PUP and migration behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckreduction.h` | Charm++ pinned Charm++ 8.0.0 reduction header. | frozen | partial | Background implementation context for `CkReduction` and `CkReductionMsg`; pinned to Charm++ 8.0.0. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech9.html` | Cactus users guide chapter C1, `Application thorns`. | living | partial | Background for CCL and thorn-file terminology emitted by ETLegacy code. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech12.html` | Cactus users guide chapter D2, CCL reference. | living | partial | Background for CCL configuration-file terminology emitted by CarpetX thorn-assembly code. |
| `https://einsteintoolkit.org/thornguide/CarpetX/CarpetX/documentation.html` | Einstein Toolkit thorn guide page, `CarpetX`. | living | partial | Background for CarpetX CCL and dependency terminology emitted by CarpetX thorn-assembly code. |

## Exclusions

- `project/`, `build/`, and `dist/`.
- Caches such as `__pycache__/`, `.mypy_cache/`, and `.pytest_cache/`.
- Compiled files, object files, executables, generated PDFs, archives, images,
  rendered artifacts, and other non-text outputs.
- Generated C/CUDA projects and generated thorns unless registered as selected
  frozen evidence.
- `optimal-plan-runs/`, `plan*.md`, `synth_plan*.md`, and `rank*.md`.
- Scratch logs, prompt transcripts, token-count reports, and latest-snapshot
  reports.

# Source Manifest

> Pointer manifest for source material. Root-level documentation sources live
> under `raw/source-docs/` so `AGENTS.md` is the only root KB document. Code,
> config, fixtures, selected logs, and build inputs stay in place. Status is
> `frozen` when a source is meant not to change and `living` when drift must
> trigger re-ingest. Last audited: 2026-06-29.

## Aggregate Sources

| Source | Provenance | Status | Mtime | Hash | Ingest |
| --- | --- | --- | --- | --- | --- |
| `core-top-level-package-modules` | Core top-level package modules from `find nrpy -maxdepth 1 -type f \( -name '*.py' -o -name '*.txt' -o -name 'py.typed' \)`, 12 files. | living | 2026-06-29 11:58:33 | sha256:4898e5f7d33e03dd3f5b30bc668746bd6043c20a29f68556aff22197ee69c6b0 | partial |
| `helpers-validation-and-reference-metric-tests` | Helpers, validation helpers, and reference metric tests from `nrpy/helpers`, `nrpy/validate_expressions`, and `nrpy/tests`, 38 files. | living | 2026-06-29 11:58:33 | sha256:039963e2bbdd73a4a7a19f62d360a3a5f5c755ca1f7d37f86841134027e80269 | partial |
| `equation-modules-and-trusted-values` | Equation modules and generated trusted-value files from `nrpy/equations`, 311 files. | living | 2026-06-29 11:58:33 | sha256:8355737f318de5fa31ba3a56f4260565a6b646e61514c6569464329feff44090 | partial |
| `infrastructure-modules-and-embedded-headers` | Infrastructure modules and embedded headers from `nrpy/infrastructures`, 364 files. | living | 2026-06-29 11:58:33 | sha256:502ba4f40c287b8f0748deee00f5fe4e2b71b2a4771efb89fa696cc57d66ec65 | partial |
| `example-generators-and-companion-scripts` | Example generators and companion scripts from `nrpy/examples`, 35 files. | living | 2026-06-29 11:58:33 | sha256:744144d1ce74d36bd2b64021e7e1e2ccd9e94e38f35a5b56fa2a90c8d672eba0 | partial |
| `ci-and-local-automation` | CI and local automation files from `.github`, 3 files. | living | 2026-06-09 17:07:33 | sha256:7a0c354fec49d1874eb1bfe456017e61e868b96127b7b0f4dd6deb52ee518ae4 | partial |

## Source Documents Moved Below Root

| Source | Provenance | Status | Mtime | Hash | Ingest |
| --- | --- | --- | --- | --- | --- |
| `raw/source-docs/original-agents.md` | Previous root agent instructions, preserved byte-for-byte before replacing `AGENTS.md`. | frozen | 2026-06-29 12:24:28 | sha256:17e4452a4449d909844ba4cd1068929460a7111f19b679643c3b3dcd1f61b693 | ingested |
| `raw/source-docs/kb-instructions.md` | KB schema and governance source moved from the repository root. | frozen | 2026-06-29 11:59:39 | sha256:5e99693989331789ddadb4fa3dafe8aa0e275fec7521e0781fd7d5017c81f73f | ingested |

## Cited Code And Config Sources

Exact cited files are registered below; they may also be covered by an aggregate
row above.

| Source | Status | Mtime | Hash |
| --- | --- | --- | --- |
| `README.md` | living | 2026-06-29 11:58:33 | sha256:624146e60406f834d15541c7af89cfd48cb4b105268c06be4c94f2325108c6ab |
| `CITATION.md` | living | 2026-06-29 11:58:33 | sha256:5b86ec909e9a85a937cc14919eb6821f15a475af05b354839460147040b62844 |
| `coding_style.md` | living | 2026-06-29 11:58:33 | sha256:4bac95dced40591cf75fd2aa3363638617e1dbd5e54b582b129382d38129ab09 |
| `setup.py` | living | 2026-06-29 11:58:33 | sha256:4ed649ad92b8daa9ff6a34e8f7a4f0e06a512984dfdb98be1b9184f273ffbb7b |
| `requirements.txt` | living | 2026-06-29 11:58:33 | sha256:6920ff9765b68b37297b36bfa2094c1346f00e6aea40797b91f70c8a9d838744 |
| `requirements-dev.txt` | living | 2026-06-29 11:58:33 | sha256:3768d3d19513df055b4063b3d196a0c8416882f69f8c681d455cede4efcb33b3 |
| `pyproject.toml` | living | 2026-06-29 11:58:33 | sha256:e11ec40a4cc26a2d97cfa833bd498cc34970a6599628ac6c769eced1e6c3317e |
| `.mypy.ini` | living | 2026-06-09 17:07:33 | sha256:076478905ec6ca0d7859954b278a1ed65c2d18429509df6f8861043ec80cba5a |
| `.pydocstyle` | living | 2026-06-09 17:07:33 | sha256:d5aa547aa92bd7a4ccd36dbf214a8a6efb359546da2b0f8d66b07fbd8eab79e3 |
| `.pylintrc` | living | 2026-06-09 17:07:33 | sha256:3421615d39e00ceb3e5467a9a50fb4969696579f80c9786cada5fd5698ed0361 |
| `.pylintrc_python36` | living | 2026-06-09 17:07:33 | sha256:f1ac03512cc879e326d638e29a6128af5cbe45cff0aa7c49b969a4e7dbf7cb5a |
| `.darglint` | living | 2026-06-09 17:07:33 | sha256:2b5da712393c7d02d88056590ab93567aeb250d44b3b96d3aa22c3c1d22dfeab |
| `.github/workflows/main.yml` | living | 2026-06-09 17:07:33 | sha256:49f831ad6976f1bf8b55719a46cf61636dcb00a7deb74d8d5a533776b7d18ab6 |
| `.github/single_file_static_analysis.sh` | living | 2026-06-09 17:07:33 | sha256:ac927452e0ada5a0f1fe9a141b414dff2898896d6c8956b1c3ce4f2acbb6df96 |
| `.github/full_nrpy_local_ci.sh` | living | 2026-06-09 17:07:33 | sha256:a14053488856de143b3302d304c2ad2b94f0c7e2c0af728ac70c11922bd1007b |
| `nrpy/c_codegen.py` | living | 2026-06-29 11:58:33 | sha256:7b2be17cea0769db3379925a79646536e026dfeabbd06b48037dbea364be08d5 |
| `nrpy/c_function.py` | living | 2026-06-29 11:58:33 | sha256:336ebca412f844a67257d560d07423397d18b4cd95ede90fd18b2e5a011c5a5f |
| `nrpy/grid.py` | living | 2026-06-29 11:58:33 | sha256:99c8844226246b0c41a1327f547cc22f321aa4bdf1f7cb8657360efc40e86dbe |
| `nrpy/params.py` | living | 2026-06-29 11:58:33 | sha256:576a2eb79c4d709fbe111ee11f6be6e383be3b2da401c1301c1f01befb800079 |
| `nrpy/indexedexp.py` | living | 2026-06-29 11:58:33 | sha256:c1b2a0021aecd24a8270c3cb1436c8f48aaca06395a73fe6a9e585025060a5a1 |
| `nrpy/reference_metric.py` | living | 2026-06-29 11:58:33 | sha256:d1ebe9753c165e14d9c8973a92e1a8b98e744a50cdc5585cbea5073a553771a7 |
| `nrpy/finite_difference.py` | living | 2026-06-29 11:58:33 | sha256:d6e4e5a216718b310c9612491fcfcca1107aca9b8d7a758cd67a9f7768422e53 |
| `nrpy/validate_expressions/validate_expressions.py` | living | 2026-06-29 11:58:33 | sha256:25c98a2f11ba153ee471317b8f51de7603b0063e1b808671d1b882e3040f76b3 |
| `nrpy/equations/general_relativity/BSSN_RHSs.py` | living | 2026-06-29 11:58:33 | sha256:d59ff6c2902060f5564159785bc4d21190621284ff2a07475a4de90642279dc3 |
| `nrpy/equations/general_relativity/BSSN_quantities.py` | living | 2026-06-29 11:58:33 | sha256:65464c8af6ba2eaddcb7c4a814f35eb8980f24d41ace4c612b38818f64906c9d |
| `nrpy/equations/general_relativity/BSSN_gauge_RHSs.py` | living | 2026-06-29 11:58:33 | sha256:b1869773d07808d3303e88c8479ec562e331324c0c68358b8714eefdf0ca3d52 |
| `nrpy/equations/general_relativity/BSSN_constraints.py` | living | 2026-06-29 11:58:33 | sha256:062030036a49afe4fb0791b7b0d8194b4ae99f3e3d759803d4804076cefe95c1 |
| `nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:e6f37670fd739a83d937dc6d6824e8019d1df1fd75242fb3f02e3534253132eb |
| `nrpy/equations/general_relativity/tests/BSSN_quantities_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:3d31524a9cc222af152cdd4276133808c9df231b2d78a94f5c5239b64293882f |
| `nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:4a4926596b65245a8f198a2c2a62a9ae75c511d16142e51647faaf774512e107 |
| `nrpy/infrastructures/BHaH/main_c.py` | living | 2026-06-29 11:58:33 | sha256:83754bad86f6449f08080852bafac056d3b7a0a0b70e73be85e28a7991c8a9fa |
| `nrpy/infrastructures/BHaH/bhah_lib.py` | living | 2026-06-29 11:58:33 | sha256:0531d413da1eb74ad48f99d8428d786939339ba214f8b607e99a2563105b5a23 |
| `nrpy/infrastructures/BHaH/simple_loop.py` | living | 2026-06-29 11:58:33 | sha256:f99c0ab022817a03b02e476f0b071ee9eee79d8bf7dd12877bd35d161bfa795d |
| `nrpy/infrastructures/BHaH/griddata_commondata.py` | living | 2026-06-29 11:58:33 | sha256:d65600e26c665d7af4c267868ee96744c9502bcea20a11d25f9ff6e8e789a836 |
| `nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py` | living | 2026-06-29 11:58:33 | sha256:7bf1bc7a473cf0a3629991249d8d9264cbbd373c354949ae32d8403efa6068e9 |
| `nrpy/infrastructures/BHaH/rfm_precompute.py` | living | 2026-06-29 11:58:33 | sha256:02b8ad1a43a65d3000a9e923a03f2292052d2f9486860962855078e148615b5b |
| `nrpy/tests/reference_metric_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:8dfc73ed8d53b1c56f2721eb085702390e04b351ef035ad167016e1acc07f476 |
| `nrpy/tests/reference_metric_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:fa562e7e0dd7ba64b95129f4c4a26b565eccad9b96b783b12e3573d093efa474 |
| `nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py` | living | 2026-06-29 11:58:33 | sha256:4dd315d1777e16c6cf03efe76713f377979e9baef55da8ea8206d71c186b51b5 |
| `nrpy/examples/wave_equation_cartesian.py` | living | 2026-06-29 11:58:33 | sha256:d3337a4f2aec6ee2d6e44967d3b33d8397116acc4fdffcd606960c0250f5d986 |
| `nrpy/examples/two_blackholes_collide.py` | living | 2026-06-29 11:58:33 | sha256:47751a1bf77ce2dee2873b437bae41e492f2e807987c52ee3fe8a58ea10ecf2c |
| `nrpy/examples/tests/sebob_consistency_check.py` | living | 2026-06-29 11:58:33 | sha256:af7a994e86157842efd5ad93459c19c8fa811a4b523ea3585a63bb2172aa67d2 |
| `nrpy/examples/tests/sebobv2_consistency_check.py` | living | 2026-06-29 11:58:33 | sha256:cfe3b83c41d9006940c3f5388dc6479abeeffb1da0791db4e9850a2bdaa394f8 |

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

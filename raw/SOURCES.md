# Source Manifest

> Pointer manifest for source material. Root-level documentation sources live
> under `raw/source-docs/` so `AGENTS.md` is the only root KB document. Code,
> config, fixtures, selected logs, and build inputs stay in place. Status is
> `frozen` when a source is meant not to change and `living` when drift must
> trigger re-ingest. Last audited: 2026-06-30.

## Aggregate Sources

| Source | Provenance | Status | Mtime | Hash | Ingest |
| --- | --- | --- | --- | --- | --- |
| `core-top-level-package-modules` | Core top-level package modules from `find nrpy -maxdepth 1 -type f \( -name '*.py' -o -name '*.txt' -o -name 'py.typed' \)`, 12 files. | living | 2026-06-29 11:58:33 | sha256:4898e5f7d33e03dd3f5b30bc668746bd6043c20a29f68556aff22197ee69c6b0 | partial |
| `helpers-package-modules` | Helper package files from `find nrpy/helpers -type f \( -name '*.py' -o -name '*.h' \)`, 21 files. | living | 2026-06-29 13:58:41 | sha256:a7f4734cc55bd7c571cd5879db516de03fb2a6551c454774ad04a94b99e567f6 | ingested |
| `helpers-validation-and-reference-metric-tests` | Helpers, validation helpers, and reference metric tests from `nrpy/helpers`, `nrpy/validate_expressions`, and `nrpy/tests`, 38 files. | living | 2026-06-29 11:58:33 | sha256:039963e2bbdd73a4a7a19f62d360a3a5f5c755ca1f7d37f86841134027e80269 | partial |
| `equation-modules-and-trusted-values` | Equation modules and generated trusted-value files from `nrpy/equations`, 311 files. | living | 2026-06-29 11:58:33 | sha256:8355737f318de5fa31ba3a56f4260565a6b646e61514c6569464329feff44090 | partial |
| `infrastructure-modules-and-embedded-headers` | Infrastructure modules and embedded headers from `nrpy/infrastructures`, 364 files. | living | 2026-06-29 15:34:39 | sha256:ae23973e29b30242334f4a1e6b9ea96884e7467edc5151805242fbde3c8b5bdc | partial |
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
| `nrpy/py_codegen.py` | living | 2026-06-29 13:58:41 | sha256:924430dce028f9c51e413e70cb8aa0bc1ca45cb6e1c9510d78a828e89caeb226 |
| `nrpy/py_function.py` | living | 2026-06-29 13:58:41 | sha256:dfe29d4cf4959e7a6f45beaf785492a49790448470ddd26a46668004fd4bc22a |
| `nrpy/c_function.py` | living | 2026-06-29 11:58:33 | sha256:336ebca412f844a67257d560d07423397d18b4cd95ede90fd18b2e5a011c5a5f |
| `nrpy/grid.py` | living | 2026-06-29 11:58:33 | sha256:99c8844226246b0c41a1327f547cc22f321aa4bdf1f7cb8657360efc40e86dbe |
| `nrpy/params.py` | living | 2026-06-29 11:58:33 | sha256:576a2eb79c4d709fbe111ee11f6be6e383be3b2da401c1301c1f01befb800079 |
| `nrpy/indexedexp.py` | living | 2026-06-29 11:58:33 | sha256:c1b2a0021aecd24a8270c3cb1436c8f48aaca06395a73fe6a9e585025060a5a1 |
| `nrpy/reference_metric.py` | living | 2026-06-29 11:58:33 | sha256:d1ebe9753c165e14d9c8973a92e1a8b98e744a50cdc5585cbea5073a553771a7 |
| `nrpy/finite_difference.py` | living | 2026-06-29 11:58:33 | sha256:d6e4e5a216718b310c9612491fcfcca1107aca9b8d7a758cd67a9f7768422e53 |
| `nrpy/helpers/expr_tree.py` | living | 2026-06-29 13:58:41 | sha256:8dc68d2981a486c869aaa46e43acf1a127e681a93fa1bfc2fe9c02bcc1fedd8e |
| `nrpy/helpers/expression_utils.py` | living | 2026-06-29 13:58:41 | sha256:d068450d37e8f2c5117a3a909b8a04559bd21a632715d8e2bda89907d141cefa |
| `nrpy/helpers/float_to_rational.py` | living | 2026-06-29 13:58:41 | sha256:5fcb30e62cd1bdbea04dd67c074bc68631ac151ba655edcd9ce5adea1ce4dda0 |
| `nrpy/helpers/functional.py` | living | 2026-06-29 13:58:41 | sha256:60914a6df754c4a0998f5f9d5a46fb35f27802e9de4bee6066fa742d1ed8a52c |
| `nrpy/helpers/cse_preprocess_postprocess.py` | living | 2026-06-29 13:58:41 | sha256:db76a8c01651365de5ab00f69ea4aff631a0e9020f84b7411b8d872549aaea71 |
| `nrpy/helpers/custom_c_codegen_functions.py` | living | 2026-06-29 13:58:41 | sha256:c0da4d28d247f60b107e3832cfd36bbe8fa16c7d12a143f96a2ee20682f4d8e9 |
| `nrpy/helpers/jax_printer.py` | living | 2026-06-29 13:58:41 | sha256:3371ea1989ae5a0fe872f8fe4916d3c875494a856865101b5c2caa5c312e6974 |
| `nrpy/helpers/simd.py` | living | 2026-06-29 13:58:41 | sha256:39afa854bd3fb780bf35cf0bc05f3a21bd74d6797bfc63d5eb31761505ade8dc |
| `nrpy/helpers/simd_intrinsics.h` | living | 2026-06-29 13:58:41 | sha256:449aaf17b8af5e71d19a9d1983e256e73479d4f9081f938c4420102a7086efca |
| `nrpy/helpers/cuda_intrinsics.h` | living | 2026-06-29 13:58:41 | sha256:1099aa36f9b9151ea73656e6dca9eb33c2f1a13e8cda230b8477e71e6e12b34f |
| `nrpy/helpers/loop.py` | living | 2026-06-29 13:58:41 | sha256:4271134f8c60895d3243778e356b9add7d034260ac2a9315cf2ec3e1cfc5221f |
| `nrpy/helpers/parallelization/gpu_kernel.py` | living | 2026-06-29 13:58:41 | sha256:eddc36f080048b6efa676ec1a3ce72aa81bdf0dd3c51a0411fa894c41788de12 |
| `nrpy/helpers/parallelization/utilities.py` | living | 2026-06-29 13:58:41 | sha256:6078e43f54d509312105cf85a5b19cb91a15c28ad63996e79cbe1ca3ea17e73c |
| `nrpy/helpers/parallel_codegen.py` | living | 2026-06-29 13:58:41 | sha256:ddadeae7617cce81e8d3d57410237579e6a8b07aba2cd26379261930f9b1b214 |
| `nrpy/helpers/parallelization/__init__.py` | living | 2026-06-29 13:58:41 | sha256:5776ba86b2f78cdae3487b55c1a954aa968e5671a326726874dd3a9c164afaa5 |
| `nrpy/helpers/generic.py` | living | 2026-06-29 13:58:41 | sha256:b795553804c68b267da49bee17b785f697f352ebd5d8a3ae7e9e0528535af7cc |
| `nrpy/helpers/cached_functions.py` | living | 2026-06-29 13:58:41 | sha256:30d914c4be69d113e1ba810d24fb47522d5420e55eabfaab457678e0deb597ca |
| `nrpy/helpers/conditional_file_updater.py` | living | 2026-06-29 13:58:41 | sha256:b63696b30756281d91a95e62964bb865aaa5515776e2238831253a9da1e8fc42 |
| `nrpy/helpers/colorize_text.py` | living | 2026-06-29 13:58:41 | sha256:17de9bfab45dcc32cc753901034426ffc79b3ecb65e2ffe767d77e6ad944726e |
| `nrpy/helpers/type_annotation_utilities.py` | living | 2026-06-29 13:58:41 | sha256:3020b9613f0f276bda86010e73e6b5ea8ec5fde21bfde35557e63e4e7a525f61 |
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
| `nrpy/infrastructures/ETLegacy/interface_ccl.py` | living | 2026-06-29 17:38:36 | sha256:8d3f47d455f421d148f38dcdc10a1cf12208e31841e0e2e5ac5c55e9c090abcc |
| `nrpy/infrastructures/ETLegacy/param_ccl.py` | living | 2026-06-29 17:38:36 | sha256:a072c2ac4b1bcad76318b6d90f06dd3924f387243adc8c2071b1a6d58a3b035e |
| `nrpy/infrastructures/ETLegacy/schedule_ccl.py` | living | 2026-06-29 17:38:36 | sha256:2cd7e09ac76f3c3e9f1ad4f86e713d124c47294b0f18d41e8e4e127be4ece090 |
| `nrpy/infrastructures/ETLegacy/make_code_defn.py` | living | 2026-06-29 17:38:36 | sha256:a198287be5e894a8afed6a3cd11b5f47e95ee7089fae1b39686e4659ea23bf9c |
| `nrpy/infrastructures/ETLegacy/CodeParameters.py` | living | 2026-06-29 17:38:36 | sha256:35f506a0bd67f27ab52c7579a8c5b4c36c03d2ecb186551d2bec2a20c414bcff |
| `nrpy/infrastructures/ETLegacy/ETLegacy_include_header.py` | living | 2026-06-29 17:38:36 | sha256:ad2c756d21847a20b62e02bbe51120989dfb051d39d4f31b073a2f2211e32056 |
| `nrpy/infrastructures/ETLegacy/simple_loop.py` | living | 2026-06-29 17:38:36 | sha256:53f36746d16264a62211fe270d357bd9b50e9c45ea737cb97469e88e4bed09b4 |
| `nrpy/infrastructures/ETLegacy/MoL_registration.py` | living | 2026-06-29 17:38:36 | sha256:854eff17bc909e506b34dfa3df4cc2589a0e0fec65244bd658bf2c6009052507 |
| `nrpy/infrastructures/ETLegacy/Symmetry_registration.py` | living | 2026-06-29 17:38:36 | sha256:3644897ef2aedc2f146a15a240594f8c32c1530a4d8182d515f0aca19f3dca1e |
| `nrpy/infrastructures/ETLegacy/boundary_conditions.py` | living | 2026-06-29 17:38:36 | sha256:d52878df89de7c1e3f58f6d730e58bbdcef45ed200174a9d0a877499ba580571 |
| `nrpy/infrastructures/ETLegacy/zero_rhss.py` | living | 2026-06-29 17:38:36 | sha256:bd354fda74ba940cc2d8cae89e52923d88e61ed3e11fa6dd638c395035a815b4 |
| `nrpy/infrastructures/ETLegacy/general_relativity/rhs_eval.py` | living | 2026-06-29 17:38:36 | sha256:db8bf50c1760ac141133659e9c8c115a878bd6ea19af6e404562df219eab7751 |
| `nrpy/infrastructures/ETLegacy/general_relativity/Ricci_eval.py` | living | 2026-06-29 17:38:36 | sha256:c3f713ea241523a0dcfcf262fcd9515a58ba5a4d68d834ed91b33ef82a7a65ff |
| `nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py` | living | 2026-06-29 17:38:36 | sha256:4e8d5129718b875e0ebe72ae9af952927e75c2ba3d4b606944c21083202b99e5 |
| `nrpy/infrastructures/ETLegacy/general_relativity/ADM_to_BSSN.py` | living | 2026-06-29 17:38:36 | sha256:d336f496ef65051c4752a1abafcc86d3f419fc3f49c0399ed3ccb0dcb3141cd6 |
| `nrpy/infrastructures/ETLegacy/general_relativity/BSSN_to_ADM.py` | living | 2026-06-29 17:38:36 | sha256:e9f092005b3f6f1c0cd1d394e4b6b4985410b9175d30e9c9723a0a3a8beab951 |
| `nrpy/infrastructures/ETLegacy/general_relativity/T4DD_to_T4UU.py` | living | 2026-06-29 17:38:36 | sha256:1199ca1ef19df9c909e0f1e693d13f0a0c32033464fae72765446c03aedce33a |
| `nrpy/infrastructures/ETLegacy/general_relativity/RegisterSlicing.py` | living | 2026-06-29 17:38:36 | sha256:2b7af0cb13282ef1c487a879c530803704198449e3f1715acdb4d8d2e622a0d7 |
| `nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py` | living | 2026-06-29 17:38:36 | sha256:31cfc957c1d0a24ea1e88b10e65b9b6ecd9c9384dd891eaa89f35b2bc5bf640e |
| `nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgammahat_constraint.py` | living | 2026-06-29 17:38:36 | sha256:be033ba2762217b552f265ea9e1674ea2199eb219445f6ca0c0b06f8f4940c0a |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py` | living | 2026-06-29 17:38:36 | sha256:dc4b0d899b21ea8545bd67fd67f27c4b3aab68e26eab279707905753c23b331f |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py` | living | 2026-06-29 17:38:36 | sha256:f69abf66f34bbdc67104f1b7b3fc2debdda88314b684dcf42d4ac4cc51288d60 |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py` | living | 2026-06-29 17:38:36 | sha256:15b73c0db1fafb6feb99e3cf0ad4af888060abe929a5b2f8b9e61c53baff808f |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py` | living | 2026-06-29 17:38:36 | sha256:ddc3d11fcd53c5a651935e655d658efde2cef54c221f5376197d69cc398de5a4 |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py` | living | 2026-06-29 17:38:36 | sha256:6f5e8144fb6e343038310f238e29c72389c8533528b0f58a56ca7b61718d6e48 |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py` | living | 2026-06-29 17:38:36 | sha256:b58fb9cb25725c9480d840fd2dee99db4b500210a36ed3b33df3759a7f5dbf73 |
| `nrpy/infrastructures/superB/main_chare.py` | living | 2026-06-29 14:47:32 | sha256:bf9e686ebf58f86cfa5b4d23aca1eab5523dde80001d87b377d6bd49b6d2620d |
| `nrpy/infrastructures/superB/timestepping_chare.py` | living | 2026-06-29 14:47:32 | sha256:4a88fadf7223f4ef5bcdd3757792bfdf8dd105139c8126ae1bbcbb02d5392d6b |
| `nrpy/infrastructures/superB/Makefile_helpers.py` | living | 2026-06-29 14:47:32 | sha256:f52184fc4cfd001bafdc48a838cf8f5fc36799459b27f580a382bd3ad6d014f4 |
| `nrpy/infrastructures/superB/numerical_grids.py` | living | 2026-06-29 14:47:32 | sha256:ffcae7d68c857cb412ea4b1f9f7203df6d8c346b1dedc27b722c0305546774b8 |
| `nrpy/infrastructures/superB/chare_communication_maps.py` | living | 2026-06-29 14:47:32 | sha256:0038c5502115bba130967c50c24a57e0b783a15e02d1158f2154ca4f0079ab1e |
| `nrpy/infrastructures/superB/CurviBoundaryConditions.py` | living | 2026-06-29 14:47:32 | sha256:84eb1899c81cf74fae6fae9bd7ce6d6dd700ebbf0914d5a3aa4fc7c6d0a24692 |
| `nrpy/infrastructures/superB/MoL.py` | living | 2026-06-29 14:47:32 | sha256:f4e69bd2f074924bda00d8c296cc5e9db1b58f4dc770fb989b109932f86f87e5 |
| `nrpy/infrastructures/superB/initial_data.py` | living | 2026-06-29 14:47:32 | sha256:397c6d6c16867be2a0f655ca2454e375b04652dc4a0b624ca6f93b5bf351ede6 |
| `nrpy/infrastructures/superB/BHaH_implementation.py` | living | 2026-06-29 14:47:32 | sha256:3373dfedefc3a55909952652522e01c9d4bcaf928174a9df81f06d6dbc62331d |
| `nrpy/infrastructures/superB/horizon_finder_chare.py` | living | 2026-06-29 14:47:32 | sha256:602b2ca252444b635ffc6d4fc05e3e2770a580b26f7f49a95b9b27146803d19c |
| `nrpy/infrastructures/superB/interpolator3d_chare.py` | living | 2026-06-29 14:47:32 | sha256:b6ee8d36c45c4bde3357020973e2ed9661d1abcf7407b94f6b0ddc264cc911eb |
| `nrpy/infrastructures/superB/diagnostics/diagnostics.py` | living | 2026-06-29 14:47:32 | sha256:904530c81902199a1ffc2426d49b36ce794d0434a6da826d77fe28a72914e5d2 |
| `nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_grid_center.py` | living | 2026-06-29 14:47:32 | sha256:625d55f6e3650d27d8661f8657cf2afdf7c14a9be2f7971d89211f2e52f2a762 |
| `nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_1d_y_and_z_axes.py` | living | 2026-06-29 14:47:32 | sha256:131276098bf002023e1b5cd6a910d84dc937cf8fe6b92478680901d6f7a64e8b |
| `nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_2d_xy_and_yz_planes.py` | living | 2026-06-29 14:47:32 | sha256:4345d6e169e9ca2b23b7d3d06b32233890a70f19074e354447a502e76dea3b64 |
| `nrpy/infrastructures/superB/general_relativity/diagnostics_nearest.py` | living | 2026-06-29 14:47:32 | sha256:1fb42b27715f846406135bcfe6da9cda18a4edb5622d5604deb8187f0e821ea7 |
| `nrpy/infrastructures/superB/general_relativity/psi4_spinweightm2_decomposition.py` | living | 2026-06-29 14:47:32 | sha256:a9f1a844019ccb8a8ec7d6707d18cea5c580b38fc5cacebadce2d889a6363915 |
| `nrpy/infrastructures/superB/nrpyelliptic/diagnostics_nearest.py` | living | 2026-06-29 14:47:32 | sha256:cdd19f27b851f3d7098a8d3e344d3cd0dd0c3e46c8f0d8f277f9c9e0822ab998 |
| `nrpy/infrastructures/superB/superB/superB_pup.py` | living | 2026-06-29 14:47:32 | sha256:e3136e5f892d543ebfdfac98ead3c568704dc93c7d72efb7acf9378e42e3dd01 |
| `nrpy/infrastructures/superB/superB/superB.h` | living | 2026-06-29 14:47:32 | sha256:47859ee15811ce98fb367ae99bb3f1acb24abb5748d8d8748c27f46c72cb7813 |
| `nrpy/infrastructures/superB/superB/superB_pup_function_prototypes.h` | living | 2026-06-29 15:34:39 | sha256:51572f03d5121d99174be9ccca70ec1d7ca3526a0b2ca263b4acfb1a112f9053 |
| `nrpy/tests/reference_metric_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:8dfc73ed8d53b1c56f2721eb085702390e04b351ef035ad167016e1acc07f476 |
| `nrpy/tests/reference_metric_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:fa562e7e0dd7ba64b95129f4c4a26b565eccad9b96b783b12e3573d093efa474 |
| `nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py` | living | 2026-06-29 11:58:33 | sha256:4dd315d1777e16c6cf03efe76713f377979e9baef55da8ea8206d71c186b51b5 |
| `nrpy/examples/wave_equation_cartesian.py` | living | 2026-06-29 11:58:33 | sha256:d3337a4f2aec6ee2d6e44967d3b33d8397116acc4fdffcd606960c0250f5d986 |
| `nrpy/examples/two_blackholes_collide.py` | living | 2026-06-29 11:58:33 | sha256:47751a1bf77ce2dee2873b437bae41e492f2e807987c52ee3fe8a58ea10ecf2c |
| `nrpy/examples/superB_two_blackholes_collide.py` | living | 2026-06-29 14:47:32 | sha256:7423369f397c72b0f61eb479ceeea3e5e613ce05b2b5ee0a1acae8ad7a84751f |
| `nrpy/examples/superB_blackhole_spectroscopy.py` | living | 2026-06-29 14:47:32 | sha256:fc69305f0472d889e4f354cb352250f38a99b68642938fcade0a0e6287a4cdf6 |
| `nrpy/examples/superB_nrpyelliptic_conformally_flat.py` | living | 2026-06-29 14:47:32 | sha256:af0ba41012b4d27258e1f99a20789d458cdbae21a179b4edcc42f7745b4a4fe2 |
| `nrpy/examples/tests/sebob_consistency_check.py` | living | 2026-06-29 11:58:33 | sha256:af7a994e86157842efd5ad93459c19c8fa811a4b523ea3585a63bb2172aa67d2 |
| `nrpy/examples/tests/sebobv2_consistency_check.py` | living | 2026-06-29 11:58:33 | sha256:cfe3b83c41d9006940c3f5388dc6479abeeffb1da0791db4e9850a2bdaa394f8 |
| `nrpy/equations/basis_transforms/jacobians.py` | living | 2026-06-29 11:58:33 | sha256:108739db4e1643d56a42fe867c17544be49858effef7be759bee798ac36c03f7 |
| `nrpy/equations/basis_transforms/tests/jacobians_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:a27cad9211a86a5bb6ed0f242fa06fd71c9300ffe774cecfdf88df5f234b4bf4 |
| `nrpy/equations/basis_transforms/tests/jacobians_GeneralRFM_fisheyeN2.py` | living | 2026-06-29 11:58:33 | sha256:68b3256a2a2fa054e4fac84ead5b0f126965334f6af3395d5344675d66644a35 |
| `nrpy/equations/basis_transforms/tests/jacobians_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:c6e76a7bdcdfa6fc9791948642fe0b1e18f919dbd06ac3ef953acacab868e219 |
| `nrpy/equations/general_relativity/ADM_to_BSSN.py` | living | 2026-06-29 11:58:33 | sha256:be03d74dce18099d8adefc1464bab872575a27661dcb4f44ce4cddfd6e408c5f |
| `nrpy/equations/general_relativity/BSSN_to_ADM.py` | living | 2026-06-29 11:58:33 | sha256:ee8e576434ffac843596f2a773c651d998e83eb7cb6fe813c582bc49428afd18 |
| `nrpy/equations/general_relativity/BSSN_to_g4Christoffel.py` | living | 2026-06-29 11:58:33 | sha256:7494813506c1ba6428b38da993816f81a92ac6bc4c9ee1013184292fba30242c |
| `nrpy/equations/general_relativity/InitialData_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:e1ace04212eca3d6f08b2b4862b3233a4a391570b05ce304120e28465648a72a |
| `nrpy/equations/general_relativity/InitialData_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:89f045217d279dfbe991c5421b7b4cec20d6b4fda622e983aeee5e6edcc146d1 |
| `nrpy/equations/general_relativity/LorentzBoost.py` | living | 2026-06-29 11:58:33 | sha256:9af347801fecffd0d41050511280812970b1af510ebc2cc4af43f879bcc79621 |
| `nrpy/equations/general_relativity/T4munu.py` | living | 2026-06-29 11:58:33 | sha256:e19e05804adeb676ad8ce0c0c2a497d3f2e6c13985682d1450e70c23b1727622 |
| `nrpy/equations/general_relativity/bhahaha/ExpansionFunctionTheta.py` | living | 2026-06-29 11:58:33 | sha256:e3178d3245e1ba8afc7933b4be4e12c83c0ce2f2eebc4d0638424025ee2cf88b |
| `nrpy/equations/general_relativity/bhahaha/HorizonSpinVorticityDipole.py` | living | 2026-06-29 11:58:33 | sha256:e545a46b7782d05847f2cbbba9d3c59849046e808a0f00cff1488e10fb98b6f2 |
| `nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py` | living | 2026-06-29 11:58:33 | sha256:18debe2904a85d7f7e1271ac5df2ccc9bda123568ea945986333e7d0b7bb31cc |
| `nrpy/equations/general_relativity/bhahaha/approx_killing_vector_spin.py` | living | 2026-06-29 11:58:33 | sha256:8fd1a37b0f9715b24a7b430f295b60ce729d63370bcecd72d529539abab00bf9 |
| `nrpy/equations/general_relativity/bhahaha/area.py` | living | 2026-06-29 11:58:33 | sha256:3ce7c6c6bac8da7b512397ac722e1b21869fd401a75dc9292a0391a221aea47c |
| `nrpy/equations/general_relativity/bhahaha/tests/ExpansionFunctionTheta_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:9cc3d8f7c6e24c24cd281b53ec5a0ebad801aa3ceaefdb00cc25b278bba388d5 |
| `nrpy/equations/general_relativity/bhahaha/tests/HorizonSpinVorticityDipole_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:f4e7bbd387d46147119b9c1b04ad7db6d02c234682d9fc12855b9a06644819ba |
| `nrpy/equations/general_relativity/bhahaha/tests/SpECTRESpinEstimate_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:6d275bad64150d0c41370dd294647722189faff3dea83ede63f0c54b2bba592e |
| `nrpy/equations/general_relativity/bhahaha/tests/approx_killing_vector_spin_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:9bf3769ed65ba9e551e4784b71b35ce09437e231646769628c697f6175faf501 |
| `nrpy/equations/general_relativity/bhahaha/tests/area_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:69a0a7412fa361e8e6af0444a06d400966eda28d1766a988b5c827185829e0cf |
| `nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py` | living | 2026-06-29 11:58:33 | sha256:cfb1b5899341eb232b331431f93f88c468f7b70fb0083dfb7fcb63bec279fbe1 |
| `nrpy/equations/general_relativity/fishbone_moncrief/tests/fishbone_moncrief.py` | living | 2026-06-29 11:58:33 | sha256:3ce59437f68bc0a119f7e27ed499c7b9c184b9b7e8aafa2f07562e196e909440 |
| `nrpy/equations/general_relativity/g4munu_conversions.py` | living | 2026-06-29 11:58:33 | sha256:44e6b856e88ab3e21968ea028aeb44751a4465eca5d37704ed2c55c3a849fd2d |
| `nrpy/equations/general_relativity/geodesics/analytic_spacetimes.py` | living | 2026-06-29 11:58:33 | sha256:e58102a859ec0809ce26037f7c8fbb878f384f614f19f2799c86c1703b2753cf |
| `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/conserved_quantities.py` | living | 2026-06-29 11:58:33 | sha256:1666ef7c79725237dd2f00d18dc17ed8a2af2c8c8d6db6b39467684b4d8396e6 |
| `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_massive.py` | living | 2026-06-29 11:58:33 | sha256:41a2639d33c5fd9af356bb5f4282883b44b95fb7848c27ddab1edc4c539d4730 |
| `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_photon.py` | living | 2026-06-29 11:58:33 | sha256:02589a7d01047eb040e370f58a3014da8b4f630a0c45600ab8f90107fd8f3bf5 |
| `nrpy/equations/general_relativity/geodesics/geodesics.py` | living | 2026-06-29 11:58:33 | sha256:9b55f35c675ca02a4f523c7221a9637a4c11b800af89acac4594d54358ef2f5e |
| `nrpy/equations/general_relativity/geodesics/tests/analytic_spacetimes_KerrSchild_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:5e6541942912c40ce2042dd81828e55d43c01158ab0d5b94eb6eba6eaa785a0b |
| `nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_massive.py` | living | 2026-06-29 11:58:33 | sha256:1e691b037f3a810df306efc1aa6168f169bb44924b861abdd726b1b5b7e9550f |
| `nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_photon.py` | living | 2026-06-29 11:58:33 | sha256:d3896c74f2899d706ecb6e15eaba61f4fc61f34370edfc418065fe6cd943255d |
| `nrpy/equations/general_relativity/psi4.py` | living | 2026-06-29 11:58:33 | sha256:d2409b60157f7205da01e2ec915f6739bf4d2dea8952c4713a0d8088b77e4940 |
| `nrpy/equations/general_relativity/psi4_tetrads.py` | living | 2026-06-29 11:58:33 | sha256:2098929f0a2a841048625da007e494979d6deee7ed95d317f839738710819112 |
| `nrpy/equations/general_relativity/tests/ADM_to_BSSN_StaticTrumpet.py` | living | 2026-06-29 11:58:33 | sha256:a25d860e19c2a65b168b383cded56f89c4531f7cc987edc0a030562aeb2688af |
| `nrpy/equations/general_relativity/tests/BSSN_to_ADM_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:790327f5cddb211acbffd806e8b92266402ad70c66bc17c7964a1ee276f7a1ac |
| `nrpy/equations/general_relativity/tests/BSSN_to_g4Christoffel_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:475f6c4858196ea4d0df4d33571aa6577ee1c892805630229493d63b040b4b0e |
| `nrpy/equations/general_relativity/tests/InitialData_Cartesian_BrillLindquist.py` | living | 2026-06-29 11:58:33 | sha256:6f82928cfc1b4a84442d40bb6a26e9bcaf2e29e47b6737b45756a93f3c348c4f |
| `nrpy/equations/general_relativity/tests/InitialData_Cartesian_Kasner.py` | living | 2026-06-29 11:58:33 | sha256:d465454fd8268be282ed2ff8d9c54ebd5ed3246fe23d2b89b96c4614a89c0e6f |
| `nrpy/equations/general_relativity/tests/InitialData_Spherical_OffsetKerrSchild.py` | living | 2026-06-29 11:58:33 | sha256:50b6f67f8f8b4a1b6f6d76737d2970292620b99f41b7a0c8551a6cbf0b452cf1 |
| `nrpy/equations/general_relativity/tests/InitialData_Spherical_StaticTrumpet.py` | living | 2026-06-29 11:58:33 | sha256:99115d6091f3627564eafccb02b951fada14497943098226c1ac265bcae0c472 |
| `nrpy/equations/general_relativity/tests/InitialData_Spherical_UIUCBlackHole.py` | living | 2026-06-29 11:58:33 | sha256:875b6649f2c89d7e5d48684670b3188984e60e89c01604f366942f665f1bf9ef |
| `nrpy/equations/general_relativity/tests/LorentzBoost.py` | living | 2026-06-29 11:58:33 | sha256:475af2a54edc28d0209d3e68852a75d9ecc2b75d5a1979e855d9729cf1314467 |
| `nrpy/equations/general_relativity/tests/T4munu.py` | living | 2026-06-29 11:58:33 | sha256:99a88e7ab763db2bc34f2945bf5c181b48a87814b96c9c7e342438fd2a040a4d |
| `nrpy/equations/general_relativity/tests/g4munu_conversions.py` | living | 2026-06-29 11:58:33 | sha256:3c1d2a0d4966f6bc61ed4b8f27d4359da5e2b29648905848d61bea3afab82757 |
| `nrpy/equations/general_relativity/tests/psi4_leave_symbolic_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:11223713bafeae40de76ceeb8a3340550b1ad33d57560f2971b8ab0a258b2dea |
| `nrpy/equations/general_relativity/tests/psi4_quasiKinnersley_SinhSpherical_rfm_precompute.py` | living | 2026-06-29 11:58:33 | sha256:6275e0ad7bc5cb32cd2ccbbf408c8f213f0af8b5ce7d576fc3ded69acecd4952 |
| `nrpy/equations/general_relativity/tests/psi4_quasiKinnersley_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:a0a5e4afef5a20fd63da2ba39ebf20defa96b7ecd0325db3eef02081a0cbebfd |
| `nrpy/equations/general_relativity/tests/psi4_tetrads_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:8b602a4173a53a1948489caecf67532c527b0f3054b97e8cf5db2b1028f74e05 |
| `nrpy/equations/generalrfm/fisheye.py` | living | 2026-06-29 11:58:33 | sha256:b50aae8404d3c30ea93f9f6b6875e184859f19a38555c2f936027e8071f63a2d |
| `nrpy/equations/generalrfm/tests/fisheye_N1.py` | living | 2026-06-29 11:58:33 | sha256:267764855eabcc037274cc092d4cfe969d77e26f352a7693f106e421c5dcf01b |
| `nrpy/equations/generalrfm/tests/fisheye_N2.py` | living | 2026-06-29 11:58:33 | sha256:8e742ffe9cdd23e2918296fd59d0e7a56fc308f868cd98f9c8834d6762ca0dd8 |
| `nrpy/equations/grhd/GRHD_equations.py` | living | 2026-06-29 11:58:33 | sha256:ad5ef0653959988c7ab57aa05c41481db63db1ee57fab3d8f59f1ee7aca6cd34 |
| `nrpy/equations/grhd/HLL_fluxes.py` | living | 2026-06-29 11:58:33 | sha256:9475b47d2e89a30df051fe09434e4010a6fcf3012c16476b5c7b199b0a258160 |
| `nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py` | living | 2026-06-29 11:58:33 | sha256:3f49914a9432dcbdfd15200ca748603c140efc1dd3fbc048a8b845fb2f093ea3 |
| `nrpy/equations/grhd/characteristic_speeds.py` | living | 2026-06-29 11:58:33 | sha256:7ad4a9cb14e9c784b3d78ed27b11c001385aebabc2099859f0041e73c0deb112 |
| `nrpy/equations/grhd/tests/GRHD_equations_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:70554a5c0ce467195677071c5e900fe8257726aa0c0f56f2951349d821fc433b |
| `nrpy/equations/grhd/tests/GRHD_equations_SinhSpherical_rfm_precompute.py` | living | 2026-06-29 11:58:33 | sha256:e43dac1e5a7564119a5eb5ea8b2bedd9e0ac3e7b517effa68b18756d15725bf9 |
| `nrpy/equations/grhd/tests/GRHD_equations_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:3f7f3b466e911e8548de65841c8f312af8c67782d919fb9d5f97bc5a11850e90 |
| `nrpy/equations/grhd/tests/HLL_fluxes.py` | living | 2026-06-29 11:58:33 | sha256:728cfd971cab2ff43035e5df96d37d1651dd2bd413ce915365d83cae061eb4e7 |
| `nrpy/equations/grhd/tests/Min_Max_and_Piecewise_Expressions.py` | living | 2026-06-29 11:58:33 | sha256:0dd22eb12a79fbaa51987d22fa398ecbb5506800cf888f9d25b36738ae97bec9 |
| `nrpy/equations/grhd/tests/characteristic_speeds.py` | living | 2026-06-29 11:58:33 | sha256:85cbacb1a01a46d601c88f20c6a20f437c3143cd77d71ecfe956888e6536c454 |
| `nrpy/equations/nrpyelliptic/ConformallyFlat_RHSs.py` | living | 2026-06-29 11:58:33 | sha256:0b8f10e30ab438540660ee32be4a8ca9d2506a026ed2da6efad627e202e56994 |
| `nrpy/equations/nrpyelliptic/ConformallyFlat_SourceTerms.py` | living | 2026-06-29 11:58:33 | sha256:248f20f6069ba794e627b56b6f69da21aef17fcd2b1d02657d6ed18f4bab0050 |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_RHSs_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:9820a7fe42e77922b1f6cc929cec79d98dce92d7e78420c238069c1306b796e7 |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_RHSs_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:1ac440559ff2384a185db78ffe0edd4613912ea6f310f8cdea8fa8ec0534f8b7 |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_SourceTerms_Cartesian.py` | living | 2026-06-29 11:58:33 | sha256:7ef6d5fcf7a99fc3e8d38c48273ad8d229d93d98adcec6eb3905a6bc76f479f9 |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_SourceTerms_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:549ed4b7ef4cb7a9697faa5260f02a7331382bfc50f2d691472fb0547e26d546 |
| `nrpy/equations/quaternion_rotations/tensor_rotation.py` | living | 2026-06-29 11:58:33 | sha256:fe986ec2b0bc79d0b1583f11ba743d4a1e04e39c98c26dc5cc3f32fbfb23255e |
| `nrpy/equations/rotation/SO3_rotations.py` | living | 2026-06-29 11:58:33 | sha256:d3c8d704b7bd8ad3e8a60ffb2d53d9ee060ec24e9895381c5d74fb94856ec4dc |
| `nrpy/equations/rotation/tests/SO3_rotations.py` | living | 2026-06-29 11:58:33 | sha256:d45bf7c986fed5beee5e9e12ff353803cd3e08387d6a5df42a726efc4cce7f2a |
| `nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities.py` | living | 2026-06-29 11:58:33 | sha256:fa3475fb420396d3ca4682073917e547d86de25c1a8e152a62c71d25a2a3a852 |
| `nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities_higher_modes.py` | living | 2026-06-29 11:58:33 | sha256:cda5d979a91be71b0a9fa74739a736fe843555a92634b4d0b69894e1c481607b |
| `nrpy/equations/seobnr/BOB_v2_waveform_quantities_kankani_etal.py` | living | 2026-06-29 11:58:33 | sha256:7b89a0e37e9f91b46903175d9cf7b1cf7a0accfe7866a57fa949cf229b128d0a |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py` | living | 2026-06-29 11:58:33 | sha256:0481b48e82f65de18356c0b7064d3a7fd7380bd734e6d7b4c05f6472e8ebc2ce |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py` | living | 2026-06-29 11:58:33 | sha256:f889447cc79a4bb8d0ae495a654c73e761fd63ed41018402cc91f657de1f9511 |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_merger_quantities.py` | living | 2026-06-29 11:58:33 | sha256:2506da2e0b3df7d3717b5083eaf5aca1532f656440ad8c773b465502d29e4ebc |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_waveform_quantities.py` | living | 2026-06-29 11:58:33 | sha256:baf87f8f8b15fbcb0ded28521f7b9fef463c18a243aaa0f98dd9b6424fcac9d2 |
| `nrpy/equations/seobnr/SEOBNRv5_coprecessing_rotations_quantities.py` | living | 2026-06-29 11:58:33 | sha256:043fb6dafce3b527a942e525bbaef11f96b24db601cdda4f148edda628fae31b |
| `nrpy/equations/seobnr/SEOBNRv5_merger_ringdown.py` | living | 2026-06-29 11:58:33 | sha256:388be5783217ee5c6c1cac1980f16f4ec5175b0b24582fca5482b228d69d1fcb |
| `nrpy/equations/seobnr/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py` | living | 2026-06-29 11:58:33 | sha256:be42c42e954932434e068c6bde3f79513da31be1e1b778ab64d4fa0f392dd0e9 |
| `nrpy/equations/seobnr/SEOBNRv5_spin_evolution_equations.py` | living | 2026-06-29 11:58:33 | sha256:2160909a521e51366a93fd3aab1798281e34e5388462f84bccd9216a73606f6f |
| `nrpy/equations/seobnr/tests/BOB_aligned_spin_waveform_quantities.py` | living | 2026-06-29 11:58:33 | sha256:094ba094b8b2aceb6cfd177c9eee910adf844407c6cd816db9f71fda79c3e682 |
| `nrpy/equations/seobnr/tests/BOB_aligned_spin_waveform_quantities_higher_modes.py` | living | 2026-06-29 11:58:33 | sha256:0ff57a4cd1f6e02fe36cd238fb87ea21a66802ed437f2d7ee60c5e8062a24e41 |
| `nrpy/equations/seobnr/tests/BOB_v2_waveform_quantities_kankani_etal.py` | living | 2026-06-29 11:58:33 | sha256:f4bebc8622bc9847756adb7d0e6e4913b6ab638b9cdf071df5107a75acad5b5c |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_Hamiltonian.py` | living | 2026-06-29 11:58:33 | sha256:1e185579b712f99ed06bcf9c3435b24df3182f64d4a2fcf029897e265c433db4 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_constants.py` | living | 2026-06-29 11:58:33 | sha256:a78f81f0294be7e88102211b6577b1260848b34c459777747fa51891671c5b35 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_merger_quantities.py` | living | 2026-06-29 11:58:33 | sha256:04c641edf9b5deec19a39a9b06fa274b257875d38132e5e27026665936244b30 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_waveform_quantities.py` | living | 2026-06-29 11:58:33 | sha256:be2ab43bbc813701b44c19f4df01a0e71abef8c449134a62458dce10129e0da4 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_coprecessing_rotations_quantities.py` | living | 2026-06-29 11:58:33 | sha256:4c804f9920e852dcf87ae5b24ccf55e33ebc64ef2be70d9a11f67592207869b4 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_merger_ringdown.py` | living | 2026-06-29 11:58:33 | sha256:aa067642b93f324be8f98af3e49f0710f9585776c760a7b490cbf8deae42e404 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py` | living | 2026-06-29 11:58:33 | sha256:9a351352b8e40f1620cb92ca5bcf99595d1a44e0512bff867ca4e61575d149f9 |
| `nrpy/equations/seobnr/tests/SEOBNRv5_spin_evolution_equations.py` | living | 2026-06-29 11:58:33 | sha256:6515f3497031d17176531b8b066e04652f685da84c63257e8e134264a205617d |
| `nrpy/equations/special_functions/spin_weighted_spherical_harmonics.py` | living | 2026-06-29 11:58:33 | sha256:a96367dda36fa0063d59c8a07881091af7cc02843041ea3029e13e78c7bf0c90 |
| `nrpy/equations/special_functions/tests/spin_weighted_spherical_harmonics.py` | living | 2026-06-29 11:58:33 | sha256:dba274310d12250690789b919af87190fc6c5d262de6295651f7a63462194fe4 |
| `nrpy/equations/tov/TOV_equations.py` | living | 2026-06-29 11:58:33 | sha256:cadc2a34fddc08a2592f956d4406764872609d51f19936d64afa0514e97f4f0d |
| `nrpy/equations/tov/tests/TOV_equations.py` | living | 2026-06-29 11:58:33 | sha256:c2386fb5d9c5b4a79e32f1adffdd1fc0cbe8260287dd8041621e08380f13834c |
| `nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py` | living | 2026-06-29 11:58:33 | sha256:8d989c28d2499916c9955b90e0f70d8249c91e9e1c1bee0758063d9e811f9c2c |
| `nrpy/equations/wave_equation/WaveEquation_RHSs.py` | living | 2026-06-29 11:58:33 | sha256:83a10adcca206cb842bc95d2cf32daaa21049b0775ab799f17e7f714e30108bd |
| `nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py` | living | 2026-06-29 11:58:33 | sha256:ada9e9f00086c92de05fd4bb9afc585f62dfd5d32ab0cceb224061ab77cfbee7 |
| `nrpy/equations/wave_equation/tests/WaveEquationCurvilinear_RHSs_Spherical.py` | living | 2026-06-29 11:58:33 | sha256:0d09ebc02ed9a718151be8894494ca039df01978a6df2229b33e18f2bf91f8c3 |
| `nrpy/equations/wave_equation/tests/WaveEquation_RHSs_WaveEquation.py` | living | 2026-06-29 11:58:33 | sha256:8589c9a62a996f54c6ea728f5a7b6332d3ec6f51ae870746028491723ad2c0ff |
| `nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_PlaneWave.py` | living | 2026-06-29 11:58:33 | sha256:d53c4153d65ce029b3a23d2cf000ba3c238fa1f222ae6031bc0bb3d433dcffb0 |
| `nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_SphericalGaussian.py` | living | 2026-06-29 11:58:33 | sha256:fedc95c05700150d150614b1e924f97824d444a71d2cfe7ece611df1f6527037 |

## External Background Sources

| Source | Provenance | Status | Accessed | Ingest | Notes |
| --- | --- | --- | --- | --- | --- |
| `https://arxiv.org/abs/2111.02424` | arXiv abstract page for arXiv:2111.02424. | living | 2026-06-29 | partial | NRPyElliptic background for hyperbolic relaxation and conformally flat binary-puncture initial data. |
| `https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html` | Intel Intrinsics Guide landing page. | living | 2026-06-29 | partial | SSE, AVX, AVX512, and intrinsic-family terminology for SIMD helper documentation. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/index.html` | NVIDIA CUDA C++ Programming Guide. | living | 2026-06-29 | partial | CUDA C++ terminology for SIMD, CUDA header, and GPU kernel helper documentation. |
| `https://link.aps.org/doi/10.1103/PhysRev.55.374` | DOI landing page for `10.1103/PhysRev.55.374`. | living | 2026-06-29 | partial | Original Oppenheimer-Volkoff stellar-equilibrium background. |
| `https://web2.ph.utexas.edu/~gsudama/pub/1967_008.pdf` | PDF URL for Goldberg et al. spin-weighted spherical-harmonic formula reference. | living | 2026-06-29 | partial | Goldberg-formula background for spin-weighted spherical harmonics. |
| `https://pubs.aip.org/aip/jmp/article/57/9/092504/648118/How-should-spin-weighted-spherical-functions-be` | Journal of Mathematical Physics article landing page. | living | 2026-06-29 | partial | Background on spin-weighted functions and quaternion viewpoints. |
| `https://rotations.berkeley.edu/geodesics-of-the-rotation-group-so3/` | Berkeley rotations course page. | living | 2026-06-29 | partial | Background on SO(3) and quaternion rotation geometry. |
| `https://github.com/charmplusplus/charm/blob/main/doc/quickstart.rst` | Charm++ quickstart from the Charm++ repository. | living | 2026-06-29 | partial | Background for `.ci` files, generated `.decl.h` and `.def.h` files, `charmc`, and `charmrun`. |
| `https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst` | Charm++ language manual from the Charm++ repository. | living | 2026-06-29 | partial | Background for chares, entry methods, proxies, SDAG, PUP, checkpoint/restart, reductions, and chare arrays. |
| `https://github.com/charmplusplus/charm/blob/main/doc/libraries/manual.rst` | Charm++ and Converse libraries manual from the Charm++ repository. | living | 2026-06-29 | partial | Background for CkIO. |
| `https://einsteintoolkit.org/usersguide/UsersGuide.html` | Cactus 4.20 Users Guide page from the Einstein Toolkit site. | living | 2026-06-30 | partial | Background for Cactus thorn-writing and Cactus terminology emitted by ETLegacy code. |
| `https://einsteintoolkit.org/referencemanual/ReferenceManual.html` | Cactus 4.20 Reference Manual page from the Einstein Toolkit site. | living | 2026-06-30 | partial | Background for `CCTK_*` thorn-writer function terminology emitted by ETLegacy code. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech9.html` | Cactus users guide chapter C1, `Application thorns`. | living | 2026-06-30 | partial | Background for CCL and thorn-file terminology emitted by ETLegacy code. |
| `https://einsteintoolkit.org/thornguide/CactusNumerical/MoL/documentation.html` | Einstein Toolkit thorn guide page, `Method of Lines`. | living | 2026-06-30 | partial | Background for MoL terminology emitted by ETLegacy registration code. |
| `https://einsteintoolkit.org/thornguide/CactusBase/Boundary/documentation.html` | Einstein Toolkit thorn guide page, `Boundary Conditions`. | living | 2026-06-30 | partial | Background for Boundary terminology emitted by ETLegacy boundary-condition code. |
| `https://einsteintoolkit.org/thornguide/EinsteinEvolve/NewRad/documentation.html` | Einstein Toolkit thorn guide page, `NewRad`. | living | 2026-06-30 | partial | Background for NewRad terminology emitted by ETLegacy boundary-condition code. |
| `https://einsteintoolkit.org/thornguide/CactusBase/CartGrid3D/documentation.html` | Einstein Toolkit thorn guide page, `CartGrid3D`. | living | 2026-06-30 | partial | Background for CartGrid3D terminology emitted by ETLegacy symmetry code. |
| `https://einsteintoolkit.org/thornguide/EinsteinBase/ADMBase/documentation.html` | Einstein Toolkit thorn guide page, `ADMBase`. | living | 2026-06-30 | partial | Background for ADMBase terminology emitted by ETLegacy GR coupling code. |
| `https://einsteintoolkit.org/thornguide/EinsteinBase/TmunuBase/documentation.html` | Einstein Toolkit thorn guide page, `TmunuBase`. | living | 2026-06-30 | partial | Background for TmunuBase terminology emitted by ETLegacy matter-coupling code. |

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

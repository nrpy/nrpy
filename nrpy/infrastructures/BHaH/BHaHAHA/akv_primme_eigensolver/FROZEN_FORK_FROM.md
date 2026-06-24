Upstream: https://github.com/primme/primme
Imported version: 3.2.3
Imported tag: v3.2.3

This directory is a frozen internal fork of the PRIMME eigensolver core used
by BHaHAHA's AKV implementation. BHaHAHA does not track or update this code as
an upstream vendor subtree.

Retained public headers:
- `primme.h`
- `primme_eigs.h`
- `akv_primme.h`
- `akv_primme_prefix_symbols.h`

Retained eigensolver sources:
- `auxiliary_eigs.c`
- `auxiliary_eigs_normal.c`
- `convergence.c`
- `correction.c`
- `factorize.c`
- `init.c`
- `inner_solve.c`
- `main_iter.c`
- `ortho.c`
- `primme_c.c`
- `primme_interface.c`
- `restart.c`
- `solve_projection.c`
- `update_W.c`
- `update_projection.c`
- `auxiliary.c`
- `blaslapack.c`
- `akv_internal_blaslapack.c`
- `memman.c`
- `wtime.c`

Retained internal headers:
- `*.h`
- `*.h`
- `blaslapack_private.h`

Deliberately removed:
- `svds/`
- Fortran interfaces
- GPU wrapper sources
- upstream makefiles, dependency generators, and tools
- generated object files

Local adaptations:
- the PRIMME symbol surface is prefixed through `akv_primme_prefix_symbols.h`
  to avoid symbol collisions when BHaHAHA is linked into larger applications
  that may also link another PRIMME build
- the public entry point for BHaHAHA code is `akv_primme.h`
- `primme.h` was reduced to the eigensolver-only interface and no longer pulls
  in the SVDS API
- BHaHAHA builds this fork as a double-only AKV eigensolver and routes
  BLAS/LAPACK calls to the local `akv_internal_blaslapack.c` fallback, so
  generated AKV projects do not require system BLAS/LAPACK libraries

License:
- BSD 3-Clause, reproduced in `COPYING.txt`

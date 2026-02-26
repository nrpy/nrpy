import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc

def p0_reverse(p0_expr: sp.Expr) -> None:
    """
    Generate and register the C function to compute p^0 for a photon.
    Fixed to correctly map symbolic metric names to flat SoA array indices.
    """
    includes = ["BHaH_defines.h"]
    desc = """@brief Computes the initial time-component of the 4-momentum (p^0).
        Solves the quadratic Hamiltonian constraint for the negative root.
        Updated for SoA architecture using global and batch indexing."""
    name = "p0_reverse"

    # Signature matching your batch integrator requirements
    params = (
        "const double *restrict metric_g4DD, "
        "double *restrict all_photons_f, "
        "const long int num_rays, "
        "const int batch_size, "
        "const long int photon_idx, "
        "const int batch_id, "
        "double *restrict p0_out"
    )

    # 1. Map symbolic names (metric_g4DD01) to local batch array indices
    # We loop exactly as g4DD_metric.py does to ensure index 'k' matches.
    metric_map = {}
    k = 0
    for i in range(4):
        for j in range(i, 4):
            # The generator outputs symbols named 'metric_g4DDij'
            old_symbol = f"metric_g4DD{i}{j}"
            new_access = f"metric_g4DD[IDX_LOCAL({k}, batch_id, batch_size)]"
            metric_map[old_symbol] = new_access
            k += 1

    # 2. Generate the Math Body (using CSE)
    # We pass the expression for the negative root of p0
    body_math = ccg.c_codegen(
        [p0_expr], ["*p0_out"], enable_cse=True, verbose=False, include_braces=False
    )
    
    # 3. Perform the substitution
    # Replace the symbolic variables with the specific SoA array indexing
    for old, new in metric_map.items():
        body_math = body_math.replace(old, new)

    # 4. Generate the Preamble and Postamble
    # Unpack spatial momentum from the global state vector (indices 5, 6, 7)
    preamble = f"""
  const double pU1 = all_photons_f[IDX_GLOBAL(5, photon_idx, num_rays)];
  const double pU2 = all_photons_f[IDX_GLOBAL(6, photon_idx, num_rays)];
  const double pU3 = all_photons_f[IDX_GLOBAL(7, photon_idx, num_rays)];
"""

    # Final assignment to the global state at index 4 (p^t)
    postamble = f"\n  all_photons_f[IDX_GLOBAL(4, photon_idx, num_rays)] = *p0_out;"

    # Project Singularity-Axiom: Portable Body Wrapper
    body = f"""
        {preamble}
        {body_math}
        {postamble}
        """
    prefunc = """
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    """
    
    postfunc = """
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """
    
    # Step 6: Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,      
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
        postfunc=postfunc  
    )
#ifndef AKV_PRIMME_H
#define AKV_PRIMME_H

#include "akv_primme_prefix_symbols.h"

#ifndef BHAHAHA_AKV_PRIMME_NAMESPACE
#define BHAHAHA_AKV_PRIMME_NAMESPACE
#endif
#include "primme.h"

/* Public frontend entrypoints exported by primme_c.c. */
#define cprimme BHAHAHA_AKV_PRIMME_PUBLIC(cprimme)
#define cprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(cprimme_normal)
#define cublas_cprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_cprimme)
#define cublas_cprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(cublas_cprimme_normal)
#define cublas_dprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_dprimme)
#define cublas_hprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_hprimme)
#define cublas_hsprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_hsprimme)
#define cublas_kcprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(cublas_kcprimme_normal)
#define cublas_kprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_kprimme)
#define cublas_kprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(cublas_kprimme_normal)
#define cublas_ksprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_ksprimme)
#define cublas_sprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_sprimme)
#define cublas_zprimme BHAHAHA_AKV_PRIMME_PUBLIC(cublas_zprimme)
#define cublas_zprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(cublas_zprimme_normal)
#define dprimme BHAHAHA_AKV_PRIMME_PUBLIC(dprimme)
#define hprimme BHAHAHA_AKV_PRIMME_PUBLIC(hprimme)
#define hsprimme BHAHAHA_AKV_PRIMME_PUBLIC(hsprimme)
#define kcprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(kcprimme_normal)
#define kprimme BHAHAHA_AKV_PRIMME_PUBLIC(kprimme)
#define kprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(kprimme_normal)
#define ksprimme BHAHAHA_AKV_PRIMME_PUBLIC(ksprimme)
#define magma_cprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_cprimme)
#define magma_cprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(magma_cprimme_normal)
#define magma_dprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_dprimme)
#define magma_hprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_hprimme)
#define magma_hsprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_hsprimme)
#define magma_kcprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(magma_kcprimme_normal)
#define magma_kprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_kprimme)
#define magma_kprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(magma_kprimme_normal)
#define magma_ksprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_ksprimme)
#define magma_sprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_sprimme)
#define magma_zprimme BHAHAHA_AKV_PRIMME_PUBLIC(magma_zprimme)
#define magma_zprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(magma_zprimme_normal)
#define sprimme BHAHAHA_AKV_PRIMME_PUBLIC(sprimme)
#define zprimme BHAHAHA_AKV_PRIMME_PUBLIC(zprimme)
#define zprimme_normal BHAHAHA_AKV_PRIMME_PUBLIC(zprimme_normal)

#endif

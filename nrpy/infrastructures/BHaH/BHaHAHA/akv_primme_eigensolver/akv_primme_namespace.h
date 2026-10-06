#ifndef AKV_PRIMME_NAMESPACE_H
#define AKV_PRIMME_NAMESPACE_H

#define BHAHAHA_AKV_JOIN2_INNER(a, b) a##b
#define BHAHAHA_AKV_JOIN2(a, b) BHAHAHA_AKV_JOIN2_INNER(a, b)

#ifdef BHAHAHA_AKV_PRIMME_NAMESPACE
#  define BHAHAHA_AKV_PRIMME_PUBLIC(name) bah_akv_##name
#  define BHAHAHA_AKV_PRIMME_PUBLIC_EXPAND(name) BHAHAHA_AKV_JOIN2(bah_akv_, name)
#else
#  define BHAHAHA_AKV_PRIMME_PUBLIC(name) name
#  define BHAHAHA_AKV_PRIMME_PUBLIC_EXPAND(name) name
#endif

#endif

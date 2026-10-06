#ifndef __CUDA_INTRINSICS_H__
#define __CUDA_INTRINSICS_H__

// Check for CUDA
#if defined(__CUDACC__)

#if defined(UPWIND_ALG)
#undef UPWIND_ALG
#endif // UPWIND_ALG

// If CUDA instructions are unavailable:
#define REAL_CUDA_ARRAY REAL
#define CUDA_WIDTH 1 // 1 double per loop iteration

// Basic Operations (Scalar)
// Plain operators, not __dadd_rn/__dsub_rn/__dmul_rn: under nvcc's default
// -fmad=true the compiler contracts a*b+c written with operators into a single
// fma.rn.f64, but it never contracts the opaque round-to-nearest intrinsics, so
// every add/multiply pair the symbolic lowering did not fuse itself would stay two
// FP64 instructions. Each operator alone is still IEEE round-to-nearest. This matches
// the non-CUDA branch below.
#define ConstCUDA(a) ((a))
#define AbsCUDA(a) (fabs((a)))
#define AddCUDA(a, b) ((a) + (b))
#define SubCUDA(a, b) ((a) - (b))
#define MulCUDA(a, b) ((a) * (b))
#define DivCUDA(a, b) __ddiv_rn((a), (b))

// Fused Multiply-Add/Subtract Operations (Scalar)
// The negations are exact and become free source modifiers on the FMA, so the
// subtract and negated forms are one instruction each rather than two.
#define FusedMulAddCUDA(a, b, c) __fma_rn((a), (b), (c))
#define FusedMulSubCUDA(a, b, c) __fma_rn((a), (b), -(c))
#define NegFusedMulAddCUDA(a, b, c) __fma_rn(-(a), (b), (c))
#define NegFusedMulSubCUDA(a, b, c) (-__fma_rn((a), (b), (c)))

// Mathematical Functions (Scalar)
#define SqrtCUDA(a) (__dsqrt_rn((a)))
#define ExpCUDA(a) (exp((a)))
#define SinCUDA(a) (sin((a)))
#define CosCUDA(a) (cos((a)))

// Load and Store Operations (Scalar)
#define WriteCUDA(a, b) (*(a) = (b))
#define ReadCUDA(a) __ldg((a))

// Upwind Algorithm (Scalar Version)
// *NOTE*: This upwinding is reversed from usual upwinding algorithms,
// because the upwinding control vector in BSSN (the shift)
// acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) ((UpwindVecU) > 0.0 ? 1.0 : 0.0)

// Initialize vector (of size one) to zero (output is REAL_CUDA_ARRAY)
#define SetZeroCUDA 0.0
// Horizontal addition (output is a double)
#define HorizAddCUDA(a) ((a)) // For scalar fallback, no horizontal addition needed

// If compiled with AVX512F CUDA instructions enabled:
#else
// If CUDA instructions are unavailable:
#define REAL_CUDA_ARRAY REAL
#define CUDA_WIDTH 1 // 1 double per loop iteration

// Basic Operations (Scalar)
#define ConstCUDA(a) ((a))
#define AbsCUDA(a) (fabs((a)))
#define AddCUDA(a, b) ((a) + (b))
#define SubCUDA(a, b) ((a) - (b))
#define MulCUDA(a, b) ((a) * (b))
#define DivCUDA(a, b) ((a) / (b))

// Fused Multiply-Add/Subtract Operations (Scalar)
#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))
#define FusedMulSubCUDA(a, b, c) ((a) * (b) - (c))
#define NegFusedMulAddCUDA(a, b, c) ((c) - (a) * (b))
#define NegFusedMulSubCUDA(a, b, c) (-((a) * (b) + (c))) // -a*b - c = -(a*b + c)

// Mathematical Functions (Scalar)
#define SqrtCUDA(a) (sqrt((a)))
#define ExpCUDA(a) (exp((a)))
#define SinCUDA(a) (sin((a)))
#define CosCUDA(a) (cos((a)))

// Load and Store Operations (Scalar)
#define WriteCUDA(a, b) (*(a) = (b))
#define ReadCUDA(a) (*(a))

// Upwind Algorithm (Scalar Version)
// *NOTE*: This upwinding is reversed from usual upwinding algorithms,
// because the upwinding control vector in BSSN (the shift)
// acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) ((UpwindVecU) > 0.0 ? 1.0 : 0.0)

// Initialize vector (of size one) to zero (output is REAL_CUDA_ARRAY)
#define SetZeroCUDA 0.0
// Horizontal addition (output is a double)
#define HorizAddCUDA(a) ((a)) // For scalar fallback, no horizontal addition needed

#endif // defined(__CUDACC__)

#endif // __CUDA_INTRINSICS_H__

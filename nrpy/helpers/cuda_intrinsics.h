// Check for CUDA
#ifdef __NVCC__
// If CUDA instructions are unavailable:
#define REAL_CUDA_ARRAY REAL
#define CUDA_width 1 // 1 double per loop iteration

// Basic Operations (Scalar)
#define ConstCUDA(a) ((a))
#define AbsCUDA(a) (fabs(a))
#define AddCUDA(a, b) __dadd_rn((a), (b))
#define SubCUDA(a, b) __dsub_rn((a), (b))
#define MulCUDA(a, b) __dmul_rn((a), (b))
#define DivCUDA(a, b) __ddiv_rn((a), (b))

// Fused Multiply-Add/Subtract Operations (Scalar)
#define FusedMulAddCUDA(a, b, c) __fma_rn((a), (b), (c))
#define FusedMulSubCUDA(a, b, c) FusedMulAddCUDA((a), (b), MulCUDA((-1.0), c))
#define NegFusedMulAddCUDA(a, b, c) SubCUDA((c), MulCUDA((a), (b)))
#define NegFusedMulSubCUDA(a, b, c) MulCUDA((-1.0),(FusedMulAddCUDA((a), (b), (c))))

// Mathematical Functions (Scalar)
#define SqrtCUDA(a) (__dsqrt_rn((a)))
#define ExpCUDA(a) (exp(a))
#define SinCUDA(a) (sin(a))
#define CosCUDA(a) (cos(a))

// Load and Store Operations (Scalar)
#define WriteCUDA(a, b) *(a) = (b)
#define ReadCUDA(a) __ldg(a)

// Upwind Algorithm (Scalar Version)
// *NOTE*: This upwinding is reversed from usual upwinding algorithms,
// because the upwinding control vector in BSSN (the shift)
// acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) ((UpwindVecU) > 0.0 ? 1.0 : 0.0)

// Initialize vector (of size one) to zero (output is REAL_CUDA_ARRAY)
#define SetZeroCUDA 0.0
// Horizontal addition (output is a double)
#define HorizAddCUDA(a) (a) // For scalar fallback, no horizontal addition needed

// If compiled with AVX512F CUDA instructions enabled:
#else
// If CUDA instructions are unavailable:
#define REAL_CUDA_ARRAY REAL
#define CUDA_width 1 // 1 double per loop iteration

// Basic Operations (Scalar)
#define ConstCUDA(a) (a)
#define AbsCUDA(a) (fabs(a))
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
#define SqrtCUDA(a) (sqrt(a))
#define ExpCUDA(a) (exp(a))
#define SinCUDA(a) (sin(a))
#define CosCUDA(a) (cos(a))

// Load and Store Operations (Scalar)
#define WriteCUDA(a, b) *(a) = (b)
#define ReadCUDA(a) *(a)

// Upwind Algorithm (Scalar Version)
// *NOTE*: This upwinding is reversed from usual upwinding algorithms,
// because the upwinding control vector in BSSN (the shift)
// acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) ((UpwindVecU) > 0.0 ? 1.0 : 0.0)

// Initialize vector (of size one) to zero (output is REAL_CUDA_ARRAY)
#define SetZeroCUDA 0.0
// Horizontal addition (output is a double)
#define HorizAddCUDA(a) (a) // For scalar fallback, no horizontal addition needed

#endif

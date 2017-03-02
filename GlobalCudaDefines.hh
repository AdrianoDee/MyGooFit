#ifndef __GLOBAL_CUDA_HH__
#define __GLOBAL_CUDA_HH__

#include <thrust/functional.h> // Needed for Thrust constants
#include <cmath>
#include <string>
using namespace std;
extern int host_callnumber;

#ifdef OMP_ON
#include "omp.h"
#define MAX_THREADS 8
#pragma omp threadprivate (host_callnumber)
#endif

// Thrust 1.7 will make the use of THRUST_DEVICE_BACKEND an error
#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_BACKEND_OMP
#if THRUST_VERSION < 100699
// Ensure backwards compatibility with Thrust 1.5
#undef THRUST_DEVICE_BACKEND
#define THRUST_DEVICE_BACKEND THRUST_DEVICE_BACKEND_OMP
#endif
// OMP target - all 'device' memory is actually on host.
#define MEM_DEVICE
#define MEM_SHARED
#define MEM_CONSTANT
#define EXEC_TARGET __host__
#define THREAD_SYNCH _Pragma("omp barrier") // valid in C99 and C++11, but probably not C++93
#define DEVICE_VECTOR thrust::host_vector
// Use char* here because I need +1 to mean "offset by one byte", not "by one sizeof(whatever)".
// Can't use void* because then the compiler doesn't know how to do pointer arithmetic.
// This will fail if sizeof(char) is more than 1. But that should never happen, right?
#define MEMCPY(target, source, count, direction) memcpy((char*) target, source, count)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) memcpy(((char*) target)+offset, source, count)
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) memcpy((char*) target, ((char*) source)+offset, count)
#define GET_FUNCTION_ADDR(fname) host_fcn_ptr = (void*) fname
#define GET_INTEGRAL_ADDR(fname) host_int_ptr = (void*) fname
#define GET_ATPOINTS_ADDR(fname) host_pnt_ptr = (void*) fname
#define SYNCH dummySynch
#define THREADIDX (omp_get_thread_num())
#define BLOCKDIM (omp_get_num_threads())
#define BLOCKIDX (1)
void dummySynch ();
// Create my own error type to avoid __host__ redefinition
// conflict in Thrust from including driver_types.h
enum gooError {gooSuccess = 0, gooErrorMemoryAllocation};
#else
// CUDA target - defaults
#define MEM_DEVICE __device__
#define MEM_SHARED __shared__
#define MEM_CONSTANT __constant__
#define EXEC_TARGET __device__
#define SYNCH cudaDeviceSynchronize
#define THREAD_SYNCH __syncthreads();
#define DEVICE_VECTOR thrust::device_vector
#define MEMCPY(target, source, count, direction) cudaMemcpy(target, source, count, direction)
#define MEMCPY_TO_SYMBOL(target, source, count, offset, direction) cudaMemcpyToSymbol(target, source, count, offset, direction)
#define GET_FUNCTION_ADDR(fname) cudaMemcpyFromSymbol((void**) &host_fcn_ptr, fname, sizeof(void*))
#define GET_INTEGRAL_ADDR(fname) cudaMemcpyFromSymbol((void**) &host_int_ptr, fname, sizeof(void*))
#define GET_ATPOINTS_ADDR(fname) cudaMemcpyFromSymbol((void**) &host_pnt_ptr, fname, sizeof(void*))
#define MEMCPY_FROM_SYMBOL(target, source, count, offset, direction) cudaMemcpyFromSymbol(target, source, count, offset, direction)
// For CUDA case, just use existing errors, renamed
#include <driver_types.h>      // Needed for cudaError_t
enum gooError {gooSuccess = cudaSuccess,
	       gooErrorMemoryAllocation = cudaErrorMemoryAllocation};
#define THREADIDX (threadIdx.x)
#define BLOCKDIM (blockDim.x)
#define BLOCKIDX (blockIdx.x)
#endif

gooError gooMalloc (void** target, size_t bytes);
gooError gooFree (void* ptr);

#define DOUBLES 1

void abortWithCudaPrintFlush (std::string file, int line);

#ifdef DOUBLES
#define root2 1.4142135623730951
#define invRootPi 0.5641895835477563
#define devPi 3.14159265358979323846

#define MLb  5.61951
#define MBd  5.27961
#define MPsi2S  3.686109
#define MJpsi  3.096916
#define MProton  0.938272046
#define MKaon  0.493677
#define MPion  0.13957018

// From PDG neutral only K*(892)
#define M892  0.89581
#define G892  0.0474

// From EvtGen
//#define M892  0.8961
//#define G892  0.0507

// From PDG
#define M800  0.682
#define G800  0.547

// K*800 Belle values: M  0.946, G  0.736 ?
//#define M800  0.931

// From Belle
//#define G800  0.578
// From PDG neutral only
#define M1410  1.414
#define G1410  0.232
#define M1430_0  1.425
#define G1430_0  0.270
#define M1430_2  1.4324
#define G1430_2  0.109
// From PDG neutral only
#define M1780_3  1.776
#define G1780_3  0.159
// From PDG
#define M2380_5  2.382
#define G2380_5  0.178

typedef double fptype;

const fptype MLb2=MLb*MLb;
const fptype MLb4=MLb2*MLb2;
const fptype MBd2=MBd*MBd;
const fptype MBd4=MBd2*MBd2;
const fptype MPsi2S2=MPsi2S*MPsi2S;
const fptype MPsi2S4=MPsi2S2*MPsi2S2;
const fptype MJpsi2=MJpsi*MJpsi;
const fptype MJpsi4=MJpsi2*MJpsi2;
const fptype MProton2=MProton*MProton;
const fptype MProton4=MProton2*MProton2;
const fptype MKaon2=MKaon*MKaon;
const fptype MKaon4=MKaon2*MKaon2;
const fptype MPion2=MPion*MPion;
const fptype MPion4=MPion2*MPion2;

// Double math functions
#define ACOS acos
#define ATAN2 atan2
#define COS cos
#define COSH cosh
#define SINH sinh
#define ERF erf
#define ERFC erfc
#define EXP exp
#define FABS fabs
#define FMOD fmod
#define LOG log
#define MODF modf
#define SIN sin
#define SQRT sqrt
#define FLOOR floor
#define POW pow
#else
typedef float fptype;

#define root2 1.4142135623730951f
#define invRootPi 0.5641895835477563f
#define devPi 3.14159265358979323846f

#define MLb  5.61951f
#define MBd  5.27961f
#define MPsi2S  3.686109f
#define MJpsi  3.096916f
#define MProton  0.938272046f
#define MKaon  0.493677f
#define MPion  0.13957018f

// From PDG neutral only K*(892)
#define M892  0.89581f
#define G892  0.0474f

// From EvtGen
//#define M892  0.8961
//#define G892  0.0507f

// From PDG
#define M800  0.682f
#define G800  0.547f

// K*800 Belle values: M  0.946, G  0.736 ?
//#define M800  0.931
// From Belle
//#define G800  0.578f
// From PDG neutral only
#define M1410  1.414f
#define G1410  0.232f
#define M1430_0  1.425f
#define G1430_0  0.270f
#define M1430_2  1.4324f
#define G1430_2  0.109f
// From PDG neutral only
#define M1780_3  1.776f
#define G1780_3  0.159f
// From PDG
#define M2380_5  2.382f
#define G2380_5  0.178f

fptype MLb2 = MLb*MLb;
fptype MLb4 = MLb2*MLb2;
fptype MBd2 = MBd*MBd;
fptype MBd4=MBd2*MBd2;
fptype MPsi2S2=MPsi2S*MPsi2S;
fptype MPsi2S4=MPsi2S2*MPsi2S2;
fptype MJpsi2=MJpsi*MJpsi;
fptype MJpsi4=MJpsi2*MJpsi2;
fptype MProton2=MProton*MProton;
fptype MProton4=MProton2*MProton2;
fptype MKaon2=MKaon*MKaon;
fptype MKaon4=MKaon2*MKaon2;
fptype MPion2=MPion*MPion;
fptype MPion4=MPion2*MPion2;

// Float math functions
#define ACOS acosf
#define ATAN2 atan2f
#define COS cosf
#define COSH coshf
#define SINH sinhf
#define ERF erff
#define ERFC erfcf
#define EXP expf
#define FABS fabsf
#define FMOD fmodf
#define LOG logf
#define MODF modff
#define SIN sinf
#define SQRT sqrtf
#define FLOOR floorf
#define POW powf
#endif



#endif

/**
 *
 * @file config.h
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Config file of the library
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-30
 *
 */
#ifndef _config_h_
#define _config_h_

/* #undef ENABLE_MPI */
#define ENABLE_STARPU
/* #undef ENABLE_CUDA */

#if defined(ENABLE_MPI)
#include <mpi.h>
#endif

#if defined(ENABLE_STARPU)
#include <starpu.h>
#if defined(ENABLE_MPI)
#include <starpu_mpi.h>
#endif
#endif

#if defined(ENABLE_CUDA)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#endif /* _config_h_ */

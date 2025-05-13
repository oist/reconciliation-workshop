#ifndef CORAX_MODEL_EIGEN_H_
#define CORAX_MODEL_EIGEN_H_

#include "corax/corax_core.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in eigen.c */

  CORAX_EXPORT int corax_update_eigen(corax_partition_t *partition,
                                      unsigned int       params_index);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_MODEL_EIGEN_H_ */

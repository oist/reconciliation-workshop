#ifndef CORAX_MODEL_GAMMA_H_
#define CORAX_MODEL_GAMMA_H_

#include "corax/corax_core.h"

/* GAMMA discretization modes */
#define CORAX_GAMMA_RATES_MEAN 0
#define CORAX_GAMMA_RATES_MEDIAN 1

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in gamma.c */

  CORAX_EXPORT int corax_compute_gamma_cats(double       alpha,
                                            unsigned int categories,
                                            double *     output_rates,
                                            int          rates_mode);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_MODEL_GAMMA_H_ */

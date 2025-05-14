#ifndef CORAX_IO_OUTPUT_H_
#define CORAX_IO_OUTPUT_H_

#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in output.c */

  CORAX_EXPORT void corax_show_pmatrix(const corax_partition_t *partition,
                                       unsigned int             index,
                                       unsigned int float_precision);

  CORAX_EXPORT void corax_show_clv(const corax_partition_t *partition,
                                   unsigned int             clv_index,
                                   int                      scaler_index,
                                   unsigned int             float_precision);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_IO_OUTPUT_H_ */

#ifndef CORAX_UTIL_COMPRESS_H_
#define CORAX_UTIL_COMPRESS_H_

#include "corax/core/common.h"
#include "corax/util/msa.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in compress.c */

  /**
   * Compresses the MSA in place. This is to say, the buffer `sequence` is
   * changed to store the compressed alignment.
   *
   * @param[in,out] sequence The alignment to compress, should be the one from a
   * `corax_msa_t`.
   *
   * @param map The sequence encoding map. For example, `corax_map_nt`.
   *
   * @param count The number of sequences, also the number of tips, also the
   * number of taxa.
   *
   * @param[out] length The length of the compressed alignment.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT unsigned int *corax_compress_site_patterns(
      char **sequence, const corax_state_t *map, int count, int *length);

  CORAX_EXPORT
  unsigned int *
  corax_compress_site_patterns_msa(corax_msa_t *        msa,
                                   const corax_state_t *map,
                                   unsigned int *       site_pattern_map);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_UTIL_COMPRESS_H_ */

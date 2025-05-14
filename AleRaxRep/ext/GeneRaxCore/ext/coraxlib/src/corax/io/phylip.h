#ifndef CORAX_IO_PHYLIP_H_
#define CORAX_IO_PHYLIP_H_

#include "corax/core/common.h"
#include "corax/util/msa.h"

/* Simple structure for handling PHYLIP parsing */
typedef struct corax_phylip_s
{
  FILE *              fp;
  char *              line;
  size_t              line_size;
  size_t              line_maxsize;
  char                buffer[CORAX_LINEALLOC];
  const unsigned int *chrstatus;
  long                no;
  long                filesize;
  long                lineno;
  long                stripped_count;
  long                stripped[256];
} corax_phylip_t;

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in phylip.c */

  CORAX_EXPORT void corax_msa_destroy(corax_msa_t *msa);

  CORAX_EXPORT corax_phylip_t *corax_phylip_open(const char *        filename,
                                                 const unsigned int *map);

  CORAX_EXPORT int corax_phylip_rewind(corax_phylip_t *fd);

  CORAX_EXPORT void corax_phylip_close(corax_phylip_t *fd);

  CORAX_EXPORT corax_msa_t *corax_phylip_parse_interleaved(corax_phylip_t *fd);

  CORAX_EXPORT corax_msa_t *corax_phylip_parse_sequential(corax_phylip_t *fd);

  CORAX_EXPORT corax_msa_t *corax_phylip_load(const char * fname,
                                              corax_bool_t interleaved);

  CORAX_EXPORT int corax_phylip_save(const char *       out_fname,
                                     const corax_msa_t *msa);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_IO_PHYLIP_H_ */

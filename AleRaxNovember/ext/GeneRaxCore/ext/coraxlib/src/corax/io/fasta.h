#ifndef CORAX_IO_FASTA_H_
#define CORAX_IO_FASTA_H_

#include "corax/core/common.h"
#include "corax/util/msa.h"

/* Simple structure for handling FASTA parsing */

typedef struct corax_fasta
{
  FILE *              fp;
  char                line[CORAX_LINEALLOC];
  const unsigned int *chrstatus;
  long                no;
  long                filesize;
  long                lineno;
  long                stripped_count;
  long                stripped[256];
} corax_fasta_t;

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in fasta.c */

  CORAX_EXPORT corax_fasta_t *corax_fasta_open(const char *        filename,
                                               const unsigned int *map);

  CORAX_EXPORT int corax_fasta_getnext(corax_fasta_t *fd,
                                       char **        head,
                                       long *         head_len,
                                       char **        seq,
                                       long *         seq_len,
                                       long *         seqno);

  CORAX_EXPORT void corax_fasta_close(corax_fasta_t *fd);

  CORAX_EXPORT long corax_fasta_getfilesize(const corax_fasta_t *fd);

  CORAX_EXPORT long corax_fasta_getfilepos(const corax_fasta_t *fd);

  CORAX_EXPORT int corax_fasta_rewind(corax_fasta_t *fd);

  corax_msa_t *corax_fasta_load(const char *fname);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_IO_FASTA_H_ */

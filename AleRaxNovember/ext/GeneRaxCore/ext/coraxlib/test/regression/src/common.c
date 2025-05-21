#include "common.h"
#include <search.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const corax_state_t odd5_map[256] = {
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0x1f, 0, 0, 0x1f, 0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0x1f,
    0, 0x01, 0x02, 0x04, 0x08, 0x0c, 0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0x01, 0x02, 0x04, 0x08, 0x0c, 0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0,
    0, 0,    0,    0,    0,    0,    0, 0, 0, 0, 0,    0, 0, 0,    0, 0};

unsigned int get_attributes(int argc, char **argv)
{
  int          i;
  unsigned int attributes = CORAX_ATTRIB_ARCH_CPU;

  for (i = 1; i < argc; ++i)
  {
    if (!strcmp(argv[i], "tv"))
    {
      /* tipvector */
      attributes |= CORAX_ATTRIB_PATTERN_TIP;
    }
    else if (!strcmp(argv[i], "sr"))
    {
      /* avx vectorization */
      attributes |= CORAX_ATTRIB_SITE_REPEATS;
    }
    else if (!strcmp(argv[i], "avx"))
    {
      /* avx vectorization */
      attributes |= CORAX_ATTRIB_ARCH_AVX;
    }
    else if (!strcmp(argv[i], "sse"))
    {
      /* sse3 vectorization */
      attributes |= CORAX_ATTRIB_ARCH_SSE;
    }
    else if (!strcmp(argv[i], "avx2"))
    {
      /* avx2 vectorization */
      attributes |= CORAX_ATTRIB_ARCH_AVX2;
    }
    else
    {
      printf("Unrecognised attribute: %s\n", argv[i]);
      exit(1);
    }
  }
  return attributes;
}

void skip_test()
{
  printf("Skip\n");
  exit(0);
}

/* note: There is no exhaustive error checking in this function.
         It is intended to use with the test datasets that were
         validated in advance. */
corax_partition_t *parse_msa(const char * filename,
                           unsigned int states,
                           unsigned int rate_cats,
                           unsigned int rate_matrices,
                           corax_utree_t *tree,
                           unsigned int attributes)
{
  return parse_msa_reduced(
      filename, states, rate_cats, rate_matrices, tree, attributes, -1);
}

corax_partition_t *parse_msa_reduced(const char * filename,
                                   unsigned int states,
                                   unsigned int rate_cats,
                                   unsigned int rate_matrices,
                                   corax_utree_t *tree,
                                   unsigned int attributes,
                                   int max_sites)
{
  unsigned int     i;
  unsigned int     taxa_count = tree->tip_count;
  corax_partition_t *partition;
  long             hdrlen, seqlen, seqno;
  char *           seq = NULL, *hdr = NULL;

  /* open FASTA file */
  corax_fasta_t *fp = corax_fasta_open(filename, corax_map_fasta);
  if (!fp)
  {
    printf("Error opening file %s", filename);
    return NULL;
  }

  /* allocate arrays to store FASTA headers and sequences */
  char **headers = (char **)calloc(taxa_count, sizeof(char *));
  char **seqdata = (char **)calloc(taxa_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; corax_fasta_getnext(fp, &hdr, &hdrlen, &seq, &seqlen, &seqno); ++i)
  {
    if (sites == -1) sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  if (corax_errno != CORAX_ERROR_FILE_EOF)
  {
    printf("Error while reading file %s", filename);
    free(headers);
    free(seqdata);
    corax_fasta_close(fp);
    return NULL;
  }

  /* close FASTA file */
  corax_fasta_close(fp);

  if (sites == -1)
  {
    printf("Unable to read alignment");
    free(headers);
    free(seqdata);
    return NULL;
  }

  if (max_sites != -1) sites = max_sites;

  partition = corax_partition_create(taxa_count,     /* tip nodes */
                                   taxa_count - 2, /* inner nodes */
                                   states,
                                   (unsigned int)sites,
                                   rate_matrices,      /* rate matrices */
                                   2 * taxa_count - 3, /* prob matrices */
                                   rate_cats,          /* rate categories */
                                   taxa_count - 2,     /* scale buffers */
                                   attributes);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(taxa_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int *data =
      (unsigned int *)malloc(taxa_count * sizeof(unsigned int));
  for (i = 0; i < taxa_count; ++i)
  {
    data[i] = tree->nodes[i]->clv_index;
    ENTRY entry;
#ifdef __APPLE__
    entry.key = xstrdup(tree->nodes[i]->label);
#else
    entry.key = tree->nodes[i]->label;
#endif
    entry.data = (void *)(data + i);

    hsearch(entry, ENTER);
  }

  for (i = 0; i < taxa_count; ++i)
  {
    ENTRY query;
    query.key    = headers[i];
    ENTRY *found = NULL;

    found = hsearch(query, FIND);

    if (!found)
    {
      printf("Sequence with header %s does not appear in the tree", headers[i]);
      free(headers);
      free(seqdata);
      return NULL;
    }

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    const corax_state_t *map = states == 4 ? corax_map_nt : corax_map_aa;
    corax_set_tip_states(partition, tip_clv_index, map, seqdata[i]);

    free(headers[i]);
    free(seqdata[i]);
  }

  /* clean */
  hdestroy();
  free(data);
  free(headers);
  free(seqdata);

  return partition;
}

int cb_full_traversal(corax_unode_t *node) { return 1; }

void show_tree (corax_unode_t * tree, int SHOW_ASCII_TREE)
{
  if(SHOW_ASCII_TREE)
  {
    printf ("\n");
    corax_utree_show_ascii (
        tree,
        CORAX_UTREE_SHOW_LABEL |
        CORAX_UTREE_SHOW_BRANCH_LENGTH |
        CORAX_UTREE_SHOW_CLV_INDEX | CORAX_UTREE_SHOW_PMATRIX_INDEX
            | CORAX_UTREE_SHOW_SCALER_INDEX);
    char * newick = corax_utree_export_newick (tree, NULL);
    printf ("%s\n\n", newick);
    free (newick);
  }
  else
  {
    printf ("ASCII tree not shown (SHOW_ASCII_TREE flag)\n");
    return;
  }
}

__attribute__((format(printf, 1, 2))) void fatal(const char *format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

void *xmalloc(size_t size)
{
  void *t;
  t = malloc(size);
  if (!t) fatal("Unable to allocate enough memory.");

  return t;
}

char *xstrdup(const char *s)
{
  size_t len = strlen(s);
  char * p   = (char *)xmalloc(len + 1);
  return strcpy(p, s);
}

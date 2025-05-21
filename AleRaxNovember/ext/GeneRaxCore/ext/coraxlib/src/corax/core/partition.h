#ifndef CORAX_CORE_PARTITION_H_
#define CORAX_CORE_PARTITION_H_

#include "corax/core/common.h"

/** @defgroup corax_attributes Attributes
 * These are flags which are used to control the behavior of a partition.
 * @{
 */

/* attribute flags */

/** Flag specifying no SIMD operations */
#define CORAX_ATTRIB_ARCH_CPU 0
/** Flag specifying only SSE3 SIMD operations */
#define CORAX_ATTRIB_ARCH_SSE (1 << 0)
/** Flag specifying only AVX SIMD operations */
#define CORAX_ATTRIB_ARCH_AVX (1 << 1)
/** Flag specifying only AVX2 SIMD operations */
#define CORAX_ATTRIB_ARCH_AVX2 (1 << 2)
/** Flag specifying only AVX512 SIMD operations */
#define CORAX_ATTRIB_ARCH_AVX512 (1 << 3)
/** Mask for the CPU architecture attributes */
#define CORAX_ATTRIB_ARCH_MASK 0xF

/**
 * Flag which indicates the use of the pattern tip optimization. Mutually
 * exclusive with the `CORAX_ATTRIB_SITE_REPEATS` flag.
 */
#define CORAX_ATTRIB_PATTERN_TIP (1 << 4)

/* ascertainment bias correction */
#define CORAX_ATTRIB_AB_LEWIS (1 << 5)
#define CORAX_ATTRIB_AB_FELSENSTEIN (2 << 5)
#define CORAX_ATTRIB_AB_STAMATAKIS (3 << 5)
#define CORAX_ATTRIB_AB_MASK (7 << 5)
#define CORAX_ATTRIB_AB_FLAG (1 << 8)

#define CORAX_ATTRIB_RATE_SCALERS (1 << 9)

/* site repeats */

/**
 * Flag indicating the use of the site repeats optimization. Mutually exclusive
 * with the `CORAX_ATTRIB_PATTERN_TIP` flag.
 */
#define CORAX_ATTRIB_SITE_REPEATS (1 << 10)

#define CORAX_ATTRIB_NONREV (1 << 11)

/** Mask for all the attributes currently defined */
#define CORAX_ATTRIB_MASK ((1 << 12) - 1)

/** @} */

struct corax_repeats;

/** @defgroup corax_partition_t corax_partition_t
 * Module concerning the corax_partition_t
 */

/**
 * A partition is a section of the genome for which all sites evolved under the
 * same model. Informally, one can think of a partition being a gene, but
 * understand that this is not always the case. In particular partitions of a
 * genome might not code for something, and there are no requirements that the
 * sites which make up a partition are even contiguous.
 *
 * Here is a checklist for creating and using a new partition:
 * - Create the partition.
 * - Set the substitution parameters,
 * - Set the tip states:
 *   - Compress the MSA first.
 *   - Set the pattern weights as well.
 * - Set the frequencies.
 * - Set the proportion of invariant sites.
 * - Set the category rates.
 * - Set the category weights.
 *
 * @ingroup corax_partition_t
 */
typedef struct corax_partition
{
  /**
   * Number of tips present in this partition. Also, the column length in the
   * MSA.
   */
  unsigned int tips;

  /**
   * Number of CLV buffers. Typically, this is the number of edges in the tree
   */
  unsigned int clv_buffers;

  /**
   * The number of "conceptual" nodes in the tree. This is to say, the number of
   * nodes in the tree, and not the number of `corax_unode_t` present in the
   * tree structure. Includes the tips.
   */
  unsigned int nodes; // tips + clv_buffer

  /**
   * Number of states the for the current partition model. For example, a DNA
   * model will use 4 states.
   */
  unsigned int states;

  /**
   * Number of sites in the MSA. Typically, this is the number of _unique_
   * sites, as there are functions to compress the MSA into only unique sites.
   */
  unsigned int sites;

  /**
   * The sum of the pattern weights. When an MSA is compressed, it is compressed
   * into a list of unique sites with each associated with a pattern weight.
   * Typically, this is the length of the uncompressed MSA. The exception is
   * when a weighted MSA is being used. In this case, it will be the sum of the
   * MSA weights.
   */
  unsigned int pattern_weight_sum;

  /**
   * How many rate matrices are present in the partition. This is different than
   * the number of rate _categories_. This is instead to be able to specify a
   * mixture model. Typically though, this is 1
   */
  unsigned int rate_matrices;

  /**
   * The number of probability matrices that will be used for computation of a
   * likelihood. Practically, this is going to be the number of edges in the
   * tree.
   */
  unsigned int prob_matrices;

  /**
   * Number of rate categories for the partition. When specifying other sizes,
   * the number of rate categories does not need to be accounted for, as it will
   * be tracked by the partition.
   */
  unsigned int rate_cats;

  /**
   * Number of scale buffers.
   */
  unsigned int scale_buffers;

  /**
   * Bitvector of the flags used for computation in the `corax_partition_t`.
   *
   * @ingroup corax_attributes
   */
  unsigned int attributes;

  /* vectorization options */

  /**
   * One of three constants depending on what architecture is being used. At the
   * time of writing these are:
   * - `CORAX_ALIGNMENT_CPU`
   * - `CORAX_ALIGNMENT_SSE`
   * - `CORAX_ALIGNMENT_AVX`
   */
  size_t alignment;

  /**
   * How many states are used, after padding. This is also the size of an
   * individual CLV buffer.
   */
  unsigned int states_padded;

  double **      clv;
  double **      pmatrix;
  double *       rates;
  double *       rate_weights;
  double **      subst_params;
  unsigned int **scale_buffer;
  double **      frequencies;
  double *       prop_invar;
  int *          invariant;
  unsigned int * pattern_weights;

  int *    eigen_decomp_valid;
  double **eigenvecs;
  double **inv_eigenvecs;
  double **eigenvals;

  /* tip-tip precomputation data */
  unsigned int    maxstates;
  unsigned char **tipchars;
  unsigned char * charmap;
  double *        ttlookup;
  corax_state_t * tipmap;

  /* ascertainment bias correction */
  int asc_bias_alloc;
  int asc_additional_sites; // partition->asc_bias_alloc ? states : 0

  /* site repeats */
  /**
   * If repeats are disabled, repeats is set to NULL. Otherwise, it points to a
   * structure holding all information required to use the site repeats
   * optimization.
   */
  struct corax_repeats *repeats;
} corax_partition_t;

/** @defgroup corax_operation_t corax_operation_t
 * Module containing structures and functions related to operations.
 */

/**
 * Structure for driving likelihood operations. In general should only be
 * created using `corax_utree_create_operations`.
 *
 * @ingroup corax_operation_t
 */
typedef struct corax_operation
{
  unsigned int parent_clv_index;
  int          parent_scaler_index;
  unsigned int child1_clv_index;
  unsigned int child1_matrix_index;
  int          child1_scaler_index;
  unsigned int child2_clv_index;
  unsigned int child2_matrix_index;
  int          child2_scaler_index;
} corax_operation_t;


#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in partition.c */

  /**
   * Creates a partition. The checklist for creating a new partition is:
   *
   * @param tips The number of tips of the tree. In phylogenetic terms, this is
   * the number of taxa.
   *
   * @param clv_buffers This is the number of CLVs that will be required to
   * compute the tree. Practically, this is the number of edges, or the number
   * of inner nodes. The number of rate categories is automatically accounted
   * for, so no need to add it.
   *
   * @param states The number of states that the model has, I.E. the type of
   * sequence data that is being worked on. Practically, this is:
   * - 2 for binary data,
   * - 4 for nucleotide data,
   * - 20 for amino acid data,
   * - 61 for codon data,
   *
   * @param sites How long is the alignment. Note that, this is going to be the
   * post compressed length of the sequence, i.e. the length that is from
   * corax_compress_site_patterns, or the number of unique site patterns.
   *
   * @param rate_matrices The number of rate matrices that are allocated. In a
   * simple and standard model, this is 1. In the case of a mixture model, this
   * should be equal to the number of classes of models. Note that having rate
   * categories still only requires 1 rate matrix, as those are computed based
   * on the single rate matrix.
   *
   * @param prob_matrices  The number of probability matrices need for
   * calculation. This will almost always be equal to the number of branches in
   * the tree. **IMPORTANT**: the number of rate categories is automatically
   * accounted for.
   *
   * @param rate_cats Number of different rate categories to consider.
   *
   * @param scale_buffers Number of scaling buffers to allocate. Practically,
   * this is equal to the number of inner nodes in the tree.
   *
   * @param attributes Please see the Attributes module.
   *
   * @return The created partition.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT corax_partition_t *
               corax_partition_create(unsigned int tips,
                                      unsigned int clv_buffers,
                                      unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_matrices,
                                      unsigned int prob_matrices,
                                      unsigned int rate_cats,
                                      unsigned int scale_buffers,
                                      unsigned int attributes);

  /**
   * Destroys the partition, deallocating the memory.
   *
   * @param partition The partition to be destroyed.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT void corax_partition_destroy(corax_partition_t *partition);

  /**
   * Set the tip states based on an MSA for a partition. This will initialize
   * the tip CLVs based on the MSA that is passed.
   *
   * @param partition Partition to set the tips for.
   *
   * @param tip_index Index of the tip to be initialized.
   *
   * @param map A predefined map from `char` to `int` The choices are:
   * - corax_map_bin: For binary data with an alphabet of 0 and 1.
   * - corax_map_nt: For nucleotide data.
   * - corax_map_aa: For amino acid data.
   *
   * @param sequence The sequence associated with that tip.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT int corax_set_tip_states(corax_partition_t *  partition,
                                        unsigned int         tip_index,
                                        const corax_state_t *map,
                                        const char *         sequence);

  CORAX_EXPORT int corax_set_tip_clv(corax_partition_t *partition,
                                     unsigned int       tip_index,
                                     const double *     clv,
                                     int                padding);

  /**
   * Sets the pattern weights for a partition.
   *
   * @param partition The partition to set the weights on.
   *
   * @param pattern_weights The array of weights. While you could set this
   * yourself, it should typically be the output of
   * `corax_compress_site_patterns`
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT void
  corax_set_pattern_weights(corax_partition_t * partition,
                            const unsigned int *pattern_weights);

  /**
   * Sets a substitution matrix for a partition.
   *
   * @param partition Partition for which the substitution matrix will be set.
   *
   * @param params_index Index of which rate matrix to use.
   *
   * @param params An array of substitution parameters. If we wanted to use the
   * following matrix
   * ```
   * *  a  b  c
   * a  *  d  e
   * b  d  *  f
   * c  e  f  *
   * ```
   * Then we would use the array
   * ```
   * double subst_params[] = {a, b, c, d, e, f}
   * ```
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT void corax_set_subst_params(corax_partition_t *partition,
                                           unsigned int       params_index,
                                           const double *     params);

  /**
   * Sets the based distribution frequencies for a partition. This needs to be
   * done before a likelihood can be computed.
   *
   * @param partition The partition for which the frequencies will be set for.
   *
   * @param params_index The model index to set the frequencies for.
   *
   * @params frequencies The array of frequencies which will be used to compute
   * a likelihood.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT void corax_set_frequencies(corax_partition_t *partition,
                                          unsigned int       params_index,
                                          const double *     frequencies);

  CORAX_EXPORT void corax_set_category_rates(corax_partition_t *partition,
                                             const double *     rates);

  CORAX_EXPORT void corax_set_category_weights(corax_partition_t *partition,
                                               const double *     rate_weights);

  CORAX_EXPORT int corax_set_asc_bias_type(corax_partition_t *partition,
                                           int                asc_bias_type);

  CORAX_EXPORT void
  corax_set_asc_state_weights(corax_partition_t * partition,
                              const unsigned int *state_weights);

  CORAX_EXPORT void corax_fill_parent_scaler(unsigned int        scaler_size,
                                             unsigned int *      parent_scaler,
                                             const unsigned int *left_scaler,
                                             const unsigned int *right_scaler);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_CORE_PARTITION_H_ */

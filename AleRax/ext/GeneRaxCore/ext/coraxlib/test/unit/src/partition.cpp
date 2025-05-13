#include <gtest/gtest.h>

#include "corax/util/hardware.h"
#include "corax/core/partition.h"

// CORAX_EXPORT corax_partition_t *
//            corax_partition_create(unsigned int tips,
//                                   unsigned int clv_buffers,
//                                   unsigned int states,
//                                   unsigned int sites,
//                                   unsigned int rate_matrices,
//                                   unsigned int prob_matrices,
//                                   unsigned int rate_cats,
//                                   unsigned int scale_buffers,
//                                   unsigned int attributes);
//
//  Creates a partition. The checklist for creating a new partition is:
//
//  @param tips The number of tips of the tree. In phylogenetic terms, this is
//  the number of taxa.
//
//  @param clv_buffers This is the number of CLVs that will be required to
//  compute the tree. Practically, this is the number of edges, or the number
//  of inner nodes. The number of rate categories is automatically accounted
//  for, so no need to add it.
//
//  @param states The number of states that the model has, I.E. the type of
//  sequence data that is being worked on. Practically, this is:
//  - 2 for binary data,
//  - 4 for nucleotide data,
//  - 20 for amino acid data,
//  - 61 for codon data,
//
//  @param sites How long is the alignment. Note that, this is going to be the
//  post compressed length of the sequence, i.e. the length that is from
//  corax_compress_site_patterns, or the number of unique site patterns.
//
//  @param rate_matrices The number of rate matrices that are allocated. In a
//  simple and standard model, this is 1. In the case of a mixture model, this
//  should be equal to the number of classes of models. Note that having rate
//  categories still only requires 1 rate matrix, as those are computed based
//  on the single rate matrix.
//
//  @param prob_matrices  The number of probability matrices need for
//  calculation. This will almost always be equal to the number of branches in
//  the tree. **IMPORTANT**: the number of rate categories is automatically
//  accounted for.
//
//  @param rate_cats Number of different rate categories to consider.
//
//  @param scale_buffers Number of scaling buffers to allocate. Practically,
//  this is equal to the number of inner nodes in the tree.
//
// The `attributes` field is a bitset that has the is combination of the
// following flags.
// - Architecture attributes: Only one may be set
//   - `CORAX_ATTRIB_ARCH_CPU`
//   - `CORAX_ATTRIB_ARCH_SSE`
//   - `CORAX_ATTRIB_ARCH_AVX`
//   - `CORAX_ATTRIB_ARCH_AVX2`
// - Ascertainment Bias: Only one of the "types" may be set, and if they are,
//   `CORAX_ATTRIB_AB_FLAG` must also be set
//   - `CORAX_ATTRIB_AB_LEWIS`
//   - `CORAX_ATTRIB_AB_FELSENSTEIN`
//   - `CORAX_ATTRIB_AB_STAMATAKIS`
//   - `CORAX_ATTRIB_AB_FLAG`
// - Scalers
//   - `CORAX_ATTRIB_RATE_SCALERS`
// - Optimizations: Only one may be set
//   - `CORAX_ATTRIB_PATTERN_TIP`
//   - `CORAX_ATTRIB_SITE_REPEATS`
//
//  @return The created partition.
//

const unsigned int tips_default          = 64;
const unsigned int clv_buffers_default   = tips_default - 1;
const unsigned int states_binary         = 2;
const unsigned int states_nucleotide     = 4;
const unsigned int states_amino_acid     = 20;
const unsigned int states_codon          = 61;
const unsigned int states_default        = states_binary;
const unsigned int sites_default         = 1337;
const unsigned int rate_matrices_default = 1;
const unsigned int rate_matrices_mixture = 2;
const unsigned int prob_matrices_default = tips_default - 1;
const unsigned int rate_cats_default     = 1;
const unsigned int scale_buffers_default = tips_default - 1;
const unsigned int attributes_default =
    CORAX_ATTRIB_ARCH_CPU | CORAX_ATTRIB_SITE_REPEATS;

void cpu_feature_detection_override()
{
  corax_hardware_probe();
  corax_hardware.sse3_present = true;
  corax_hardware.avx_present  = true;
  corax_hardware.avx2_present = false;
}

void cpu_feature_detection_reset() { corax_hardware_probe(); }

TEST(Partition, Basic)
{
  corax_partition_t *partition = nullptr;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes_default));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  corax_partition_destroy(partition);
}

TEST(Partition, MultipleArchitectureAttributesSet_none)
{
  cpu_feature_detection_override();
  assert(CORAX_HAS_CPU_FEATURE(sse3_present));
  assert(CORAX_HAS_CPU_FEATURE(avx_present));
  assert(!CORAX_HAS_CPU_FEATURE(avx2_present));

  corax_partition_t *partition  = nullptr;
  const unsigned int attributes = 0;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  // The best one should be autoselected.
  EXPECT_TRUE(partition->alignment & CORAX_ALIGNMENT_CPU);
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_SSE);
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_AVX);
  EXPECT_EQ(partition->states_padded, states_default);
  corax_partition_destroy(partition);
  cpu_feature_detection_reset();
}

TEST(Partition, MultipleArchitectureAttributesSet_cpu)
{
  cpu_feature_detection_override();
  assert(CORAX_HAS_CPU_FEATURE(sse3_present));
  assert(CORAX_HAS_CPU_FEATURE(avx_present));
  assert(!CORAX_HAS_CPU_FEATURE(avx2_present));

  corax_partition_t *partition  = nullptr;
  const unsigned int attributes = CORAX_ATTRIB_ARCH_CPU;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  // The best one should be autoselected.
  EXPECT_TRUE(partition->alignment & CORAX_ALIGNMENT_CPU);
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_SSE);
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_AVX);
  EXPECT_EQ(partition->states_padded, states_default);
  corax_partition_destroy(partition);
  cpu_feature_detection_reset();
}

TEST(Partition, MultipleArchitectureAttributesSet_sse3)
{
  cpu_feature_detection_override();
  assert(CORAX_HAS_CPU_FEATURE(sse3_present));
  assert(CORAX_HAS_CPU_FEATURE(avx_present));
  assert(!CORAX_HAS_CPU_FEATURE(avx2_present));

  corax_partition_t *partition = nullptr;
  const unsigned int attributes =
      CORAX_ATTRIB_ARCH_CPU | CORAX_ATTRIB_ARCH_SSE | CORAX_ATTRIB_ARCH_AVX2;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  // The best one should be autoselected.
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_CPU);
  EXPECT_TRUE(partition->alignment & CORAX_ALIGNMENT_SSE);
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_AVX);
  EXPECT_EQ(partition->states_padded, (states_default + 1) & 0xFFFFFFFE);
  corax_partition_destroy(partition);
  cpu_feature_detection_reset();
}

TEST(Partition, MultipleArchitectureAttributesSet_avx)
{
  cpu_feature_detection_override();
  assert(CORAX_HAS_CPU_FEATURE(sse3_present));
  assert(CORAX_HAS_CPU_FEATURE(avx_present));
  assert(!CORAX_HAS_CPU_FEATURE(avx2_present));

  corax_partition_t *partition  = nullptr;
  const unsigned int attributes = CORAX_ATTRIB_ARCH_SSE | CORAX_ATTRIB_ARCH_AVX;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  // The best one should be autoselected.
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_CPU);
  EXPECT_FALSE(partition->alignment & CORAX_ALIGNMENT_SSE);
  EXPECT_TRUE(partition->alignment & CORAX_ALIGNMENT_AVX);
  EXPECT_EQ(partition->states_padded, (states_default + 3) & 0xFFFFFFFC);
  corax_partition_destroy(partition);
  cpu_feature_detection_reset();
}

TEST(Partition, CheckFields)
{

  corax_partition_t *partition  = nullptr;
  const unsigned int attributes = CORAX_ATTRIB_ARCH_CPU | CORAX_ATTRIB_ARCH_SSE
                                  | CORAX_ATTRIB_ARCH_AVX
                                  | CORAX_ATTRIB_ARCH_AVX2;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  /* Check clv fields */

  EXPECT_EQ(partition->clv_buffers, clv_buffers_default);
  EXPECT_EQ(partition->tips, tips_default);
  EXPECT_EQ(partition->nodes, tips_default + clv_buffers_default);

  EXPECT_EQ(partition->sites, sites_default);

  EXPECT_EQ(partition->clv_buffers, clv_buffers_default);
  EXPECT_NE(partition->clv, nullptr);

  for (size_t i = 0; i < partition->clv_buffers; ++i)
  {
    EXPECT_NE(partition->clv[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; ++j)
    {
      EXPECT_EQ(partition->clv[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->frequencies, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; ++i)
  {
    EXPECT_NE(partition->frequencies[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; ++j)
    {
      EXPECT_EQ(partition->frequencies[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->frequencies, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; ++i)
  {
    EXPECT_NE(partition->frequencies[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; ++j)
    {

      EXPECT_EQ(partition->frequencies[i][j], 0.0);
    }
  }

  EXPECT_EQ(partition->rate_cats, rate_cats_default);
  EXPECT_NE(partition->rates, nullptr);
  EXPECT_NE(partition->rate_weights, nullptr);

  for (size_t i = 0; i < partition->rate_cats; i++)
  {
    EXPECT_EQ(partition->rates[i], 0.0);
    EXPECT_EQ(partition->rate_weights[i], 1.0);
  }

  EXPECT_EQ(partition->scale_buffers, scale_buffers_default);
  EXPECT_NE(partition->scale_buffer, nullptr);
  size_t scalar_count = partition->sites * partition->rate_cats;
  for (size_t i = 0; i < partition->scale_buffers; i++)
  {
    EXPECT_NE(partition->scale_buffer[i], nullptr);
    for (size_t j = 0; j < scalar_count; j++)
    {
      EXPECT_EQ(partition->scale_buffer[i][j], 0);
    }
  }

  /* Check rate matrix fields */

  EXPECT_EQ(partition->rate_matrices, rate_matrices_default);
  EXPECT_NE(partition->subst_params, nullptr);

  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->subst_params[i], nullptr);
  }

  /* Check prob matrix fields */

  EXPECT_EQ(partition->prob_matrices, prob_matrices_default);
  EXPECT_NE(partition->pmatrix, nullptr);
  size_t matrix_size = partition->states * partition->states_padded;
  for (size_t i = 0; i < partition->prob_matrices; i++)
  {
    EXPECT_NE(partition->pmatrix[i], nullptr);
    for (size_t j = 0; j < matrix_size; j++)
    {
      EXPECT_EQ(partition->pmatrix[i][j], 0.0);
    }
  }

  /* Check eigen fields */

  EXPECT_NE(partition->eigenvecs, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->eigenvecs[i], nullptr);
    for (size_t j = 0; j < matrix_size; j++)
    {
      EXPECT_EQ(partition->eigenvecs[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->inv_eigenvecs, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->inv_eigenvecs[i], nullptr);
    for (size_t j = 0; j < matrix_size; j++)
    {
      EXPECT_EQ(partition->inv_eigenvecs[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->eigenvals, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->eigenvals[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; j++)
    {
      EXPECT_EQ(partition->eigenvals[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->eigen_decomp_valid, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_EQ(partition->eigen_decomp_valid[i], 0);
  }

  /* Check site repeats fields */

  EXPECT_EQ(partition->repeats, nullptr);

  corax_partition_destroy(partition);
}

TEST(Partition, CheckFieldsSiteRepeats)
{

  corax_partition_t *partition = nullptr;
  const unsigned int attributes =
      CORAX_ATTRIB_ARCH_CPU | CORAX_ATTRIB_ARCH_SSE | CORAX_ATTRIB_ARCH_AVX
      | CORAX_ATTRIB_ARCH_AVX2 | CORAX_ATTRIB_SITE_REPEATS;
  ASSERT_NO_THROW(partition = corax_partition_create(tips_default,
                                                     clv_buffers_default,
                                                     states_default,
                                                     sites_default,
                                                     rate_matrices_default,
                                                     prob_matrices_default,
                                                     rate_cats_default,
                                                     scale_buffers_default,
                                                     attributes));
  ASSERT_NE(partition, nullptr);
  ASSERT_NE(partition, static_cast<corax_partition_t *>(CORAX_FAILURE));
  ASSERT_FALSE(corax_errno);

  /* Check clv fields */

  EXPECT_EQ(partition->clv_buffers, clv_buffers_default);
  EXPECT_EQ(partition->tips, tips_default);
  EXPECT_EQ(partition->nodes, tips_default + clv_buffers_default);

  EXPECT_EQ(partition->sites, sites_default);

  EXPECT_EQ(partition->clv_buffers, clv_buffers_default);
  EXPECT_NE(partition->clv, nullptr);

  for (size_t i = 0; i < partition->clv_buffers; ++i)
  {
    EXPECT_EQ(partition->clv[i], nullptr);
  }

  EXPECT_NE(partition->frequencies, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; ++i)
  {
    EXPECT_NE(partition->frequencies[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; ++j)
    {
      EXPECT_EQ(partition->frequencies[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->frequencies, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; ++i)
  {
    EXPECT_NE(partition->frequencies[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; ++j)
    {

      EXPECT_EQ(partition->frequencies[i][j], 0.0);
    }
  }

  EXPECT_EQ(partition->rate_cats, rate_cats_default);
  EXPECT_NE(partition->rates, nullptr);
  EXPECT_NE(partition->rate_weights, nullptr);

  for (size_t i = 0; i < partition->rate_cats; i++)
  {
    EXPECT_EQ(partition->rates[i], 0.0);
    EXPECT_EQ(partition->rate_weights[i], 1.0);
  }

  EXPECT_EQ(partition->scale_buffers, scale_buffers_default);
  EXPECT_NE(partition->scale_buffer, nullptr);
  for (size_t i = 0; i < partition->scale_buffers; i++)
  {
    EXPECT_EQ(partition->scale_buffer[i], nullptr);
  }

  /* Check rate matrix fields */

  EXPECT_EQ(partition->rate_matrices, rate_matrices_default);
  EXPECT_NE(partition->subst_params, nullptr);

  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->subst_params[i], nullptr);
  }

  /* Check prob matrix fields */

  EXPECT_EQ(partition->prob_matrices, prob_matrices_default);
  EXPECT_NE(partition->pmatrix, nullptr);
  size_t matrix_size = partition->states * partition->states_padded;
  for (size_t i = 0; i < partition->prob_matrices; i++)
  {
    EXPECT_NE(partition->pmatrix[i], nullptr);
    for (size_t j = 0; j < matrix_size; j++)
    {
      EXPECT_EQ(partition->pmatrix[i][j], 0.0);
    }
  }

  /* Check eigen fields */

  EXPECT_NE(partition->eigenvecs, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->eigenvecs[i], nullptr);
    for (size_t j = 0; j < matrix_size; j++)
    {
      EXPECT_EQ(partition->eigenvecs[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->inv_eigenvecs, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->inv_eigenvecs[i], nullptr);
    for (size_t j = 0; j < matrix_size; j++)
    {
      EXPECT_EQ(partition->inv_eigenvecs[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->eigenvals, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_NE(partition->eigenvals[i], nullptr);
    for (size_t j = 0; j < partition->states_padded; j++)
    {
      EXPECT_EQ(partition->eigenvals[i][j], 0.0);
    }
  }

  EXPECT_NE(partition->eigen_decomp_valid, nullptr);
  for (size_t i = 0; i < partition->rate_matrices; i++)
  {
    EXPECT_EQ(partition->eigen_decomp_valid[i], 0);
  }

  /* Check site repeats fields */

  EXPECT_NE(partition->repeats, nullptr);

  corax_partition_destroy(partition);
}

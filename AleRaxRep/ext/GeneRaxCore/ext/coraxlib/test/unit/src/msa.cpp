#include "corax/corax.h"
#include "environment.hpp"
#include <gtest/gtest.h>


TEST(MSA, column_entropies)
{
    std::string filename = env->msa_filename();
    const char* c_filename = filename.c_str();
    corax_msa_t* msa = corax_phylip_load(c_filename, CORAX_TRUE);

    EXPECT_TRUE(msa != NULL);

    int n_taxa = msa->count;
    EXPECT_EQ(n_taxa, 68);

    int n_sites = msa->length;
    EXPECT_EQ(n_sites, 766);

    double *column_entropies = corax_msa_column_entropies(msa, 4, corax_map_nt);

    EXPECT_NEAR(column_entropies[0], 0.0, 0.01);
    EXPECT_NEAR(column_entropies[5], 0.322757, 0.01);

    corax_msa_destroy(msa);
    free(column_entropies);
}

TEST(MSA, entropy)
{
    std::string filename = env->msa_filename();
    const char* c_filename = filename.c_str();

    corax_msa_t* msa = corax_phylip_load(c_filename, CORAX_TRUE);

    double entropy = corax_msa_entropy(msa, 4, corax_map_nt);
    EXPECT_NEAR(entropy, 0.19863, 0.01);

    corax_msa_destroy(msa);
}

TEST(MSA, weighted_entropy)
{
    std::string filename = env->msa_filename();
    const char* c_filename = filename.c_str();

    corax_msa_t* msa = corax_phylip_load(c_filename, CORAX_TRUE);

    ASSERT_NE(msa, nullptr);

    unsigned int * weights = corax_compress_site_patterns_msa(msa,
                                                              corax_map_nt,
                                                              NULL);

    ASSERT_NE(weights, nullptr);

    corax_msa_stats_t *msa_stats = corax_msa_compute_stats(msa,
                                                           4, corax_map_nt,
                                                           weights,
                                                           CORAX_MSA_STATS_ENTROPY);

    ASSERT_NE(msa_stats, nullptr);

    double entropy = msa_stats->entropy;
    EXPECT_NEAR(entropy, 0.19863, 0.01);

    corax_msa_destroy_stats(msa_stats);
    corax_msa_destroy(msa);
    free(weights);
}


TEST(MSA, pattern_entropy)
{
  std::string filename = env->small_msa_filename();
  const char* c_filename = filename.c_str();
  corax_msa_t* msa = corax_fasta_load(c_filename);
  unsigned int* site_pattern_map = (unsigned int *)calloc(msa->length, sizeof(unsigned int));
  unsigned int* site_weights = corax_compress_site_patterns_msa(msa, corax_map_nt, site_pattern_map);

  double pattern_entropy = corax_msa_pattern_entropy(msa, site_weights, corax_map_nt);

  /* The MSA "small.fasta" has 10 taxa and 5 different patterns with the following number of occurrences
   *  'TTTTTTTTTT': 26,
      'CCCCCCCCCC': 16,
      'AAAAAAAAAA': 28,
      'GGGGGGGGGG': 21,
      '----------': 431
   * => The pattern entropy computes as:
   * mult = (26 * log(26) + 16 * log(16) + ... + 431 * log(431))
   * In this case this is ~ 2900.801214
   */

  EXPECT_NEAR(pattern_entropy, 2900.801214, 0.001);

  free(site_weights);
  free(site_pattern_map);
  corax_msa_destroy(msa);
}

TEST(MSA, bollback)
{
    std::string filename = env->small_msa_filename();
    const char* c_filename = filename.c_str();
    corax_msa_t* msa = corax_fasta_load(c_filename);
    unsigned int* site_pattern_map = (unsigned int *)calloc(msa->length, sizeof(unsigned int));
    unsigned int* site_weights = corax_compress_site_patterns_msa(msa, corax_map_nt, site_pattern_map);

    double bollback = corax_msa_bollback_multinomial(msa, site_weights, corax_map_nt);

    /* The MSA "small.fasta" has 10 taxa and 5 different patterns with the following number of occurrences
     *  'TTTTTTTTTT': 26,
        'CCCCCCCCCC': 16,
        'AAAAAAAAAA': 28,
        'GGGGGGGGGG': 21,
        '----------': 431
     * => The bollback computes as:
     * mult = (26 * log(26) + 16 * log(16) + ... + 431 * log(431)) - num_sites * log(num_sites)
     * In this case this is ~ -365.70
     */

    EXPECT_NEAR(bollback, -365.701267, 0.001);

    free(site_weights);
    free(site_pattern_map);
    corax_msa_destroy(msa);
}

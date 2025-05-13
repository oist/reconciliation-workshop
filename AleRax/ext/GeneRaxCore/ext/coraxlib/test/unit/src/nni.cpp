/* 
 * TESTS INFO:
 * In TESTS simple0 (testing corax_algo_nni_local() ) and simple1 (testing corax_algo_nni_round() ), 
 * the topology of the tree expected to change. 
 * 
 * In TESTS simple2 (testing corax_algo_nni_local() ) and simple3 (testing corax_algo_nni_round() ), 
 * the topology of the tree expected to remain the same. 
 * 
 * In test cases simple 0,1,2,3 , we apply the corresponding functions to the initial tree topologies 
 * (either corax_algo_nni_local() or corax_algo_nni_round() ), and, afterwads, we check whether
 * the rf distance between the output topologies and the expected ones is zero.
 * 
 * In simple4 we check whether the CORAX_NNI_ROUND_TRIPLET_ERROR
 */

#include "corax/corax.h"
#include "environment.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <unistd.h>

#define check_triplet_error(nni_logl)                               \
    {                                                               \
        EXPECT_DOUBLE_EQ(nni_logl, CORAX_FAILURE);                  \
        EXPECT_EQ(corax_errno, CORAX_NNI_ROUND_TRIPLET_ERROR);      \
    }

#define check_likelihoods(logl1, logl2)                             \
    {                                                               \
        EXPECT_DOUBLE_EQ(logl1, logl2);                             \
    }

#define check_nni_logl(nni_logl, logl1)                             \
    {                                                               \
        EXPECT_TRUE(nni_logl>=logl1);                               \
    }

#define check_rf_dist(split1, split2, tip_nodes_count)                              \
    {                                                                               \
        auto rf = corax_utree_split_rf_distance(splits[0],                          \
                                                splits[1],                          \
                                                (unsigned int) tip_nodes_count);    \
        EXPECT_EQ(rf, 0);                                                           \
    }                                                                               \


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

void check_trees(corax_utree_t* tree, const char* _check_tree, int root_index, int tip_nodes_count){
    
    // check
    corax_utree_t* check_tree = corax_utree_parse_newick_string_unroot(_check_tree);
    
    // make taxon a root 
    check_tree->vroot = check_tree->nodes[root_index]->back;
    //show_tree(check_tree->vroot, 1);
    
    // Create splits for calculating rf distance
    size_t _num_trees = 2;
    std::vector<corax_split_t *> splits(_num_trees);
    std::vector<corax_utree_t *> pars_trees(_num_trees);
    std::vector<corax_unode_t *> pars_roots(_num_trees);

    pars_trees[0] = tree;
    pars_roots[0] = tree->vroot;

    pars_trees[1] = check_tree;
    pars_roots[1] = check_tree->vroot;


    std::map<std::string, unsigned int> labelToId;
    for (int i = 0; i < tip_nodes_count; ++i) {
      labelToId.insert({std::string(pars_trees[0]->nodes[i]->label), i});
    }
    for (auto _tree: pars_trees) {
      for (int i = 0; i < tip_nodes_count; ++i) {
        auto leaf = _tree->nodes[i];
        auto id = labelToId.at(std::string(leaf->label));
        leaf->node_index = leaf->clv_index = id;
      }
    }

    for (size_t i = 0; i < _num_trees; ++i)
    {
        splits[i] = corax_utree_split_create(pars_roots[i],
                                              (unsigned int) tip_nodes_count,
                                              nullptr);
    }

    check_rf_dist(splits[0], splits[1],  (unsigned int) tip_nodes_count);

    for (auto s: splits)
        corax_utree_split_destroy(s);
    
    corax_utree_destroy(check_tree, NULL);
}

double calculate_likelihood(corax_partition_t *partition,
                            corax_unode_t * node,
                            corax_unode_t ** travbuffer,
                            corax_operation_t * operations,
                            unsigned int * matrix_indices,
                            double * branch_lengths,
                            unsigned int* params_indices)
{
    unsigned int traversal_size,
                    ops_count,
                    matrix_count;

    corax_utree_traverse(node,
                        CORAX_TREE_TRAVERSE_POSTORDER,
                        [](corax_unode_t*){return 1;},
                        travbuffer,
                        &traversal_size);
    

    corax_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                                matrix_indices, operations, &matrix_count,
                                &ops_count);

    corax_update_prob_matrices (partition, params_indices, matrix_indices, branch_lengths,
                                    matrix_count);

    corax_update_clvs (partition, operations, ops_count);

    return corax_compute_edge_loglikelihood (partition,
                                        node->clv_index,
                                        node->scaler_index,
                                        node->back->clv_index,
                                        node->back->scaler_index,
                                        node->pmatrix_index,
                                        params_indices,
                                        NULL);
}


corax_treeinfo_t* initialize_treeinfo(corax_unode_t* root,
                                        corax_partition_t* partition,
                                        double alpha,
                                        unsigned int* params_indices,
                                        int tip_nodes_count)
{
    corax_treeinfo_t* treeinfo;
    treeinfo = corax_treeinfo_create(root, (unsigned int) tip_nodes_count, 1, CORAX_BRLEN_UNLINKED);
    corax_treeinfo_init_partition(treeinfo, 
                                    0,
                                    partition, 
                                    CORAX_OPT_PARAM_ALL, 
                                    CORAX_GAMMA_RATES_MEAN, 
                                    alpha, 
                                    params_indices, 
                                    NULL);
    
    return treeinfo;
}

TEST(coraxlib_nni, simple0)
{  
    corax_utree_t* tree = corax_utree_parse_newick_string_unroot("((a:0.1,b:0.1)root,(c:0.1,d:0.1):0.1);");
    corax_partition_t * partition;
    corax_operation_t * operations;
    corax_unode_t ** trav_buffer;
    corax_treeinfo_t * treeinfo;

    double alpha = 1.0;
    double logl1, logl2;
    
    int tip_nodes_count = 4;
    int inner_nodes_count = 2;

    unsigned int params_indices[4] = {0,0,0,0};

    partition = corax_partition_create(4,       /* Tip CLVs */
                                        2,       /* Inner CLVs */
                                        4,       /* States */
                                        6,       /* Sequence length */
                                        1,       /* Models (sets of subst. params)*/
                                        5,       /* P matrices */
                                        4,       /* Rate categories */
                                        2,       /* Scale buffers */
                                        CORAX_ATTRIB_ARCH_AVX);
    
    double branch_lengths[5] = { 0.1, 0.1, 0.1, 0.1, 0.1};
    unsigned int matrix_indices[5] = { 0, 1, 2, 3, 4};
    double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
    double subst_params[6] = {1,1,1,1,1,1};
    double rate_cats[4];
    corax_compute_gamma_cats(alpha, 4, rate_cats, CORAX_GAMMA_RATES_MEAN);
    
    /* set */
    corax_set_frequencies(partition, 0, frequencies);
    corax_set_subst_params(partition, 0, subst_params);
    corax_set_category_rates(partition, rate_cats);

    corax_set_tip_states(partition, 0, corax_map_nt, "AGGCTT"); // a
    corax_set_tip_states(partition, 1, corax_map_nt, "ACGATC"); // b
    corax_set_tip_states(partition, 2, corax_map_nt, "AGGCTT"); // c
    corax_set_tip_states(partition, 3, corax_map_nt, "ACGATC"); // d

    corax_unode_t* final_root = tree->nodes[2];
    
    // calculating lolg 
    trav_buffer = new corax_unode_t *[2*tip_nodes_count-2];
    operations = new corax_operation_t[inner_nodes_count];

    corax_unode_t* node = tree->vroot;
    logl1 = calculate_likelihood(partition,
                                node,
                                trav_buffer,
                                operations,
                                matrix_indices,
                                branch_lengths,
                                params_indices);
    
    
    treeinfo = initialize_treeinfo(node, partition, alpha, params_indices, tip_nodes_count);

    logl2 = corax_treeinfo_compute_loglh(treeinfo, 0);
    
    // Checking if two likelihoods are equal
    check_likelihoods(logl1, logl2);
    
    // NNI local
    double nni_logl = corax_algo_nni_local(treeinfo, 
                                            CORAX_OPT_BLO_NEWTON_FAST, 
                                            CORAX_OPT_MIN_BRANCH_LEN, 
                                            CORAX_OPT_MAX_BRANCH_LEN, 
                                            CORAX_OPT_DEFAULT_SMOOTHINGS, 
                                            CORAX_OPT_DEFAULT_EPSILON);
    tree->vroot = final_root->back;
    // checking

    // Checking for improved likelihood + topology
    check_nni_logl(nni_logl,logl1);

    const char * _check_tree = "((a:0.1,c:0.1)root,(b:0.1,d:0.1):0.1);";
    check_trees(tree, _check_tree, 2, tip_nodes_count);

    delete[] trav_buffer;
    delete[] operations;

    corax_partition_destroy(partition);
    corax_treeinfo_destroy(treeinfo);
    corax_utree_destroy(tree, NULL);
}


TEST(coraxlib_nni, simple1)
{  
    corax_utree_t* tree = corax_utree_parse_newick_string_unroot("(((d:0.1,e:0.1),c:0.1),(a:0.1,b:0.1));");
    corax_partition_t * partition;
    corax_operation_t * operations;
    corax_unode_t ** trav_buffer;
    corax_treeinfo_t * treeinfo;

    double alpha = 1.0;
    double logl1, logl2;
    
    int tip_nodes_count = 5;
    int inner_nodes_count = 3;

    unsigned int params_indices[4] = {0,0,0,0};

    partition = corax_partition_create((unsigned int) tip_nodes_count,  /* Tip CLVs */
                                        (unsigned int) inner_nodes_count,       /* Inner CLVs */
                                        4,       /* States */
                                        6,       /* Sequence length */
                                        1,       /* Models (sets of subst. params)*/
                                        7,       /* P matrices */
                                        4,       /* Rate categories */
                                        (unsigned int) inner_nodes_count,       /* Scale buffers */
                                        CORAX_ATTRIB_ARCH_AVX);
    
    double branch_lengths[7] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    unsigned int matrix_indices[7] = { 0, 1, 2, 3, 4, 5, 6};
    double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
    double subst_params[6] = {1,1,1,1,1,1};
    double rate_cats[4];
    corax_compute_gamma_cats(alpha, 4, rate_cats, CORAX_GAMMA_RATES_MEAN);
    
    /* set */
    corax_set_frequencies(partition, 0, frequencies);
    corax_set_subst_params(partition, 0, subst_params);
    corax_set_category_rates(partition, rate_cats);

    corax_set_tip_states(partition, 0, corax_map_nt, "CAATCC"); // a
    corax_set_tip_states(partition, 1, corax_map_nt, "CAAACG"); // b
    corax_set_tip_states(partition, 2, corax_map_nt, "CAATCC"); // c
    corax_set_tip_states(partition, 3, corax_map_nt, "CAATCG"); // d
    corax_set_tip_states(partition, 4, corax_map_nt, "CAAACG"); // e

    corax_unode_t* final_root = tree->nodes[0];
    
    // calculating lolg 
    trav_buffer = new corax_unode_t *[2*tip_nodes_count-2];
    operations = new corax_operation_t[inner_nodes_count];

    corax_unode_t* node = tree->vroot;
    logl1 = calculate_likelihood(partition,
                                node,
                                trav_buffer,
                                operations,
                                matrix_indices,
                                branch_lengths,
                                params_indices);
    
    treeinfo = initialize_treeinfo(node, partition, alpha, params_indices, tip_nodes_count);

    logl2 = corax_treeinfo_compute_loglh(treeinfo, 0);

    // Checking if two likelihoods are equal
    check_likelihoods(logl1, logl2);

    // --------------- NNI -----------------------------
    // NNI round
    double tolerance = 0.1;
    double nni_logl = corax_algo_nni_round(treeinfo, tolerance, 
                                            CORAX_OPT_BLO_NEWTON_FAST, 
                                            CORAX_OPT_MIN_BRANCH_LEN, 
                                            CORAX_OPT_MAX_BRANCH_LEN, 
                                            CORAX_OPT_DEFAULT_SMOOTHINGS, 
                                            CORAX_OPT_DEFAULT_EPSILON,
                                            false);
    // --------------------------------------------------

    // logl test
    check_nni_logl(nni_logl,logl1);

    // make taxon a root
    tree->vroot = final_root->back;
    
    //show_tree(tree->vroot, 1);
    const char* _check_tree = "(((a:0.1,d:0.1),e:0.1),(b:0.1,c:0.1));";
    check_trees(tree, _check_tree, 2, tip_nodes_count);
    
    delete[] trav_buffer;
    delete[] operations;

    corax_partition_destroy(partition);
    corax_treeinfo_destroy(treeinfo);
    corax_utree_destroy(tree, NULL);
}

TEST(coraxlib_nni, simple2)
{  
    corax_utree_t* tree = corax_utree_parse_newick_string_unroot("((a:0.1,b:0.1)root,(c:0.1,d:0.1):0.1);");
    corax_partition_t * partition;
    corax_operation_t * operations;
    corax_unode_t ** trav_buffer;
    corax_treeinfo_t * treeinfo;

    double alpha = 1.0;
    double logl1, logl2;
    
    int tip_nodes_count = 4;
    int inner_nodes_count = 2;

    unsigned int params_indices[4] = {0,0,0,0};

    partition = corax_partition_create(4,       /* Tip CLVs */
                                        2,       /* Inner CLVs */
                                        4,       /* States */
                                        5,       /* Sequence length */
                                        1,       /* Models (sets of subst. params)*/
                                        5,       /* P matrices */
                                        4,       /* Rate categories */
                                        2,       /* Scale buffers */
                                        CORAX_ATTRIB_ARCH_AVX);
    
    double branch_lengths[5] = { 0.1, 0.1, 0.1, 0.1, 0.1};
    unsigned int matrix_indices[5] = { 0, 1, 2, 3, 4};
    double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
    double subst_params[6] = {1,1,1,1,1,1};
    double rate_cats[4];
    corax_compute_gamma_cats(alpha, 4, rate_cats, CORAX_GAMMA_RATES_MEAN);
    
    /* set */
    corax_set_frequencies(partition, 0, frequencies);
    corax_set_subst_params(partition, 0, subst_params);
    corax_set_category_rates(partition, rate_cats);

    corax_set_tip_states(partition, 0, corax_map_nt, "ACGTA"); // a
    corax_set_tip_states(partition, 1, corax_map_nt, "ACGTA"); // b
    corax_set_tip_states(partition, 2, corax_map_nt, "CCGTA"); // c
    corax_set_tip_states(partition, 3, corax_map_nt, "CCGTA"); // d

    corax_unode_t* final_root = tree->nodes[2];
    
    // calculating lolg 
    trav_buffer = new corax_unode_t *[2*tip_nodes_count-2];
    operations = new corax_operation_t[inner_nodes_count];

    corax_unode_t* node = tree->vroot;
    logl1 = calculate_likelihood(partition,
                                node,
                                trav_buffer,
                                operations,
                                matrix_indices,
                                branch_lengths,
                                params_indices);
    
    
    treeinfo = initialize_treeinfo(node, partition, alpha, params_indices, tip_nodes_count);

    logl2 = corax_treeinfo_compute_loglh(treeinfo, 0);
    
    // Checking if two likelihoods are equal
    check_likelihoods(logl1, logl2);
    
    // NNI local
    double nni_logl = corax_algo_nni_local(treeinfo, 
                                            CORAX_OPT_BLO_NEWTON_FAST, 
                                            CORAX_OPT_MIN_BRANCH_LEN, 
                                            CORAX_OPT_MAX_BRANCH_LEN, 
                                            CORAX_OPT_DEFAULT_SMOOTHINGS, 
                                            CORAX_OPT_DEFAULT_EPSILON);
    tree->vroot = final_root->back;

    // checking

    // Checking for improved likelihood + topology
    check_nni_logl(nni_logl,logl1);
    
    const char * _check_tree = "((a:0.1,b:0.1)root,(c:0.1,d:0.1):0.1);";
    check_trees(tree, _check_tree, 2, tip_nodes_count);

    delete[] trav_buffer;
    delete[] operations;

    corax_partition_destroy(partition);
    corax_treeinfo_destroy(treeinfo);
    corax_utree_destroy(tree, NULL);
}

TEST(coraxlib_nni, simple3)
{  
    corax_utree_t* tree = corax_utree_parse_newick_string_unroot("((a:0.1,b:0.1),(c:0.1, d:0.1)root,((e:0.1,f:0.1),(g:0.1,h:0.1)):0.1);");
    corax_partition_t * partition;
    corax_operation_t * operations;
    corax_unode_t ** trav_buffer;
    corax_treeinfo_t * treeinfo;

    double alpha = 1.0;
    double logl1, logl2;
    
    int tip_nodes_count = 8;
    int inner_nodes_count = tip_nodes_count-2;

    unsigned int params_indices[4] = {0,0,0,0};
    unsigned int num_of_p_matrices = 2*((unsigned int)tip_nodes_count) - 3;

    partition = corax_partition_create((unsigned int) tip_nodes_count,       /* Tip CLVs */
                                        (unsigned int) inner_nodes_count,       /* Inner CLVs */
                                        4,       /* States */
                                        8,       /* Sequence length */
                                        1,       /* Models (sets of subst. params)*/
                                        num_of_p_matrices,       /* P matrices */
                                        4,       /* Rate categories */
                                        (unsigned int) inner_nodes_count,       /* Scale buffers */
                                        CORAX_ATTRIB_ARCH_AVX);
    
    double branch_lengths[15];
    unsigned int matrix_indices[15];
    for(int i = 0; i<15; i++){
        branch_lengths[i] = 0.1;
        matrix_indices[i] = (unsigned int) i;
    }

    double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
    double subst_params[6] = {1,1,1,1,1,1};
    double rate_cats[4];
    corax_compute_gamma_cats(alpha, 4, rate_cats, CORAX_GAMMA_RATES_MEAN);
    
    /* set */
    corax_set_frequencies(partition, 0, frequencies);
    corax_set_subst_params(partition, 0, subst_params);
    corax_set_category_rates(partition, rate_cats);

    corax_set_tip_states(partition, 0, corax_map_nt, "AGCTACGT"); // a
    corax_set_tip_states(partition, 1, corax_map_nt, "AGCTACGT"); // b
    corax_set_tip_states(partition, 2, corax_map_nt, "AGGAACGT"); // c
    corax_set_tip_states(partition, 3, corax_map_nt, "AGGAACGT"); // d
    corax_set_tip_states(partition, 4, corax_map_nt, "ACGTCCCT"); // e
    corax_set_tip_states(partition, 5, corax_map_nt, "ACGTCCCT"); // b
    corax_set_tip_states(partition, 6, corax_map_nt, "ACGTAGCT"); // c
    corax_set_tip_states(partition, 7, corax_map_nt, "ACGTAGCT"); // d

    corax_unode_t* final_root = tree->nodes[0];

    // calculating lolg 
    trav_buffer = new corax_unode_t *[2*tip_nodes_count-2];
    operations = new corax_operation_t[inner_nodes_count];

    corax_unode_t* node = tree->vroot;
    logl1 = calculate_likelihood(partition,
                                node,
                                trav_buffer,
                                operations,
                                matrix_indices,
                                branch_lengths,
                                params_indices);
    
    treeinfo = initialize_treeinfo(node, partition, alpha, params_indices, tip_nodes_count);

    logl2 = corax_treeinfo_compute_loglh(treeinfo, 0);

    // Checking if two likelihoods are equal
    check_likelihoods(logl1, logl2);
    
    // --------------- NNI -----------------------------
    // NNI round
    double tolerance = 0.1;
    double nni_logl = corax_algo_nni_round(treeinfo, tolerance, 
                                            CORAX_OPT_BLO_NEWTON_FAST, 
                                            CORAX_OPT_MIN_BRANCH_LEN, 
                                            CORAX_OPT_MAX_BRANCH_LEN, 
                                            CORAX_OPT_DEFAULT_SMOOTHINGS, 
                                            CORAX_OPT_DEFAULT_EPSILON,
                                            false);
    // --------------------------------------------------

    tree->vroot = final_root->back;

    // logl check
    check_nni_logl(nni_logl,logl1);

    // checking 
    const char* _check_tree = "((a:0.1,b:0.1),(c:0.1, d:0.1)root,((e:0.1,f:0.1),(g:0.1,h:0.1)):0.1);";
    check_trees(tree, _check_tree, 0, tip_nodes_count);
    
    delete[] trav_buffer;
    delete[] operations;

    corax_partition_destroy(partition);
    corax_treeinfo_destroy(treeinfo);
    corax_utree_destroy(tree, NULL);
}

TEST(coraxlib_nni, simple4)
{  
    corax_utree_t* tree = corax_utree_parse_newick_string_unroot("((a:0.1,b:0.1)root,c:0.1);");
    corax_partition_t * partition;
    corax_operation_t * operations;
    corax_unode_t ** trav_buffer;
    corax_treeinfo_t * treeinfo;

    double alpha = 1.0;
    double logl1, logl2;
    
    int tip_nodes_count = 3;
    int inner_nodes_count = 1;

    unsigned int params_indices[4] = {0,0,0,0};

    partition = corax_partition_create(3,       /* Tip CLVs */
                                        1,       /* Inner CLVs */
                                        4,       /* States */
                                        5,       /* Sequence length */
                                        1,       /* Models (sets of subst. params)*/
                                        3,       /* P matrices */
                                        4,       /* Rate categories */
                                        1,       /* Scale buffers */
                                        CORAX_ATTRIB_ARCH_AVX);
    
    double branch_lengths[3] = { 0.1, 0.1, 0.1};
    unsigned int matrix_indices[3] = { 0, 1, 2};
    double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
    double subst_params[6] = {1,1,1,1,1,1};
    double rate_cats[4];
    corax_compute_gamma_cats(alpha, 4, rate_cats, CORAX_GAMMA_RATES_MEAN);
    
    /* set */
    corax_set_frequencies(partition, 0, frequencies);
    corax_set_subst_params(partition, 0, subst_params);
    corax_set_category_rates(partition, rate_cats);

    corax_set_tip_states(partition, 0, corax_map_nt, "ACGTA"); // a
    corax_set_tip_states(partition, 1, corax_map_nt, "ACGTA"); // b
    corax_set_tip_states(partition, 2, corax_map_nt, "CCGTA"); // c

    // calculating lolg 
    trav_buffer = new corax_unode_t *[2*tip_nodes_count-2];
    operations = new corax_operation_t[inner_nodes_count];

    corax_unode_t* node = tree->vroot;
    logl1 = calculate_likelihood(partition,
                                node,
                                trav_buffer,
                                operations,
                                matrix_indices,
                                branch_lengths,
                                params_indices);
    
    
    treeinfo = initialize_treeinfo(node, partition, alpha, params_indices, tip_nodes_count);

    logl2 = corax_treeinfo_compute_loglh(treeinfo, 0);
    
    // Checking if two likelihoods are equal
    check_likelihoods(logl1, logl2);
    
    // --------------- NNI -----------------------------
    // NNI round
    double tolerance = 0.1;
    double nni_logl = corax_algo_nni_round(treeinfo, tolerance, 
                                            CORAX_OPT_BLO_NEWTON_FAST, 
                                            CORAX_OPT_MIN_BRANCH_LEN, 
                                            CORAX_OPT_MAX_BRANCH_LEN, 
                                            CORAX_OPT_DEFAULT_SMOOTHINGS, 
                                            CORAX_OPT_DEFAULT_EPSILON,
                                            false);
    // --------------------------------------------------

    check_triplet_error(nni_logl);
    
    delete[] trav_buffer;
    delete[] operations;

    corax_partition_destroy(partition);
    corax_treeinfo_destroy(treeinfo);
    corax_utree_destroy(tree, NULL);
}
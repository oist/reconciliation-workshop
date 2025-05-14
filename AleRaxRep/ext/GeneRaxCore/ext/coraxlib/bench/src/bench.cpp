#include <benchmark/benchmark.h>

#include <corax/corax.h>
#include <iostream>
#include "PLLTreeInfo.hpp"
#include <string>
#include <vector>

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GET_DATA(s) STRINGIFY(BENCH_DATA/s)


template <class ...Args>
void BM_Parse_Newick(benchmark::State& state, Args&&... args) {
  auto argsTuple = std::make_tuple(std::move(args)...);
  auto treePath = std::get<0>(argsTuple);
  std::vector<corax_utree_t *> toDestroy;
  for (auto _ : state) {
    auto tree = corax_utree_parse_newick_rooted(treePath);
    corax_utree_destroy(tree, nullptr);
  }
}

template <class T>
std::unique_ptr<PLLTreeInfo> createTreeInfo(T &argsTuple)
{
  auto treePath = std::get<0>(argsTuple);
  auto msaPath= std::get<1>(argsTuple);
  auto modelStr = std::get<2>(argsTuple);
  auto repeatsEnabled = std::get<3>(argsTuple);
  auto vectorizationAttribute = std::get<4>(argsTuple);
  bool isNewickAFile = false;
  return std::make_unique<PLLTreeInfo>(treePath,
      isNewickAFile,
      msaPath,
      modelStr,
      repeatsEnabled,
      vectorizationAttribute);
}

template <class ...Args>
void BM_kernel_likelihood(benchmark::State& state, Args&&... args) {
  auto argsTuple = std::make_tuple(std::move(args)...);
  auto treeInfoWrapper = createTreeInfo(argsTuple);
  auto treeinfo = treeInfoWrapper->getTreeInfo();
  double ll = corax_treeinfo_compute_loglh(treeinfo, false);
  for (auto _ : state) {
      auto v = corax_compute_edge_loglikelihood(treeinfo->partitions[0],
                                       treeinfo->root->clv_index,
                                       treeinfo->root->scaler_index,
                                       treeinfo->root->back->clv_index,
                                       treeinfo->root->back->scaler_index,
                                       treeinfo->root->pmatrix_index,
                                       treeinfo->param_indices[0],
                                       nullptr);
      assert(v == ll);
  }
}

static int cb_full_traversal(corax_unode_t *node)
{
  CORAX_UNUSED(node);
  return CORAX_SUCCESS;
}

template <class ...Args>
void BM_kernel_partial(benchmark::State& state, Args&&... args) {
  auto argsTuple = std::make_tuple(std::move(args)...);
  auto treeInfoWrapper = createTreeInfo(argsTuple);
  auto treeinfo = treeInfoWrapper->getTreeInfo();
  double ll = corax_treeinfo_compute_loglh(treeinfo, false);
  unsigned int traversal_size = 0;
  unsigned int ops_count = 0;
  corax_utree_traverse(treeinfo->root,
                          CORAX_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          treeinfo->travbuffer,
                          &traversal_size);
  corax_utree_create_operations(
      (const corax_unode_t *const *)treeinfo->travbuffer,
      traversal_size,
      NULL,
      NULL,
      treeinfo->operations,
      NULL,
      &ops_count);
  for (auto _ : state) {
    corax_update_clvs(treeinfo->partitions[0], 
        treeinfo->operations, 
        ops_count);
  }
}

template <class ...Args>
void BM_search(benchmark::State& state, Args&&... args) {
  auto argsTuple = std::make_tuple(std::move(args)...);
  for (auto _ : state) {
    auto treeInfoWrapper = createTreeInfo(argsTuple);
    auto treeinfo = treeInfoWrapper->getTreeInfo();
    double ll1 = corax_treeinfo_compute_loglh(treeinfo, false);
    
    double cutoff = 0.0; //wtf is that
    cutoff_info_t cutoff_info;
    int minRadius = 1;
    int maxRadius = 1;
    int toKeep = 20;
    int thorough = true;
    const double RAXML_BRLEN_MIN = 0.000001;
    const double RAXML_BRLEN_MAX = 100.0;
    const int RAXML_BRLEN_SMOOTHINGS = 8;
    unsigned int traversal_size = 0;
    unsigned int ops_count = 0;
    corax_utree_traverse(treeinfo->root,
                            CORAX_TREE_TRAVERSE_POSTORDER,
                            cb_full_traversal,
                            treeinfo->travbuffer,
                            &traversal_size);
    corax_utree_create_operations(
        (const corax_unode_t *const *)treeinfo->travbuffer,
        traversal_size,
        NULL,
        NULL,
        treeinfo->operations,
        NULL,
        &ops_count);
    corax_update_clvs(treeinfo->partitions[0], 
        treeinfo->operations, 
        ops_count);
    double ll2 = corax_treeinfo_compute_loglh(treeinfo, false);
    auto ll3 = corax_algo_spr_round(treeinfo,
      static_cast<int>(minRadius),
      static_cast<int>(maxRadius),
      static_cast<int>(toKeep), // params.ntopol_keep
      static_cast<int>(thorough), // THOROUGH
      0, //int brlen_opt_method,
      RAXML_BRLEN_MIN,
      RAXML_BRLEN_MAX,
      RAXML_BRLEN_SMOOTHINGS,
      0.1,
      nullptr,
      cutoff,
      0.1); //double subtree_cutoff);
  }
}


#define GET_TREE(DATA) GET_DATA(trees/DATA ## .newick)
#define GET_MSA(DATA) GET_DATA(msas/DATA ## .phy)
#define GET_MSA_FASTA(DATA) GET_DATA(msas/DATA ## .fasta)

#define BENCH_TREE(TREE) BENCHMARK_CAPTURE(BM_Parse_Newick, \
    tree_ ## TREE ## _taxa,  \
    GET_TREE(TREE));
 

#define BENCH_KERNEL(KER, VEC, REPEATS, MODEL, DATA) \
  BENCHMARK_CAPTURE(BM_kernel_ ## KER, \
    Dataset_ ## DATA  ## _ ## MODEL ## _ ## VEC ## _Repeat ## REPEATS ,\
    GET_TREE(DATA), \
    GET_MSA(DATA),\
    STRINGIFY(MODEL),\
    REPEATS,\
    CORAX_ATTRIB_ARCH_ ## VEC \
    );

#define BENCH_KERNEL_ALL_REPEATS(KER, VEC, MODEL, DATA) \
  BENCH_KERNEL(KER, VEC, true, MODEL, DATA) \
  BENCH_KERNEL(KER, VEC, false, MODEL, DATA) \

#define BENCH_KERNEL_ALL_VEC(KER, MODEL, DATA) \
   BENCH_KERNEL_ALL_REPEATS(KER, AVX2, MODEL, DATA) \
   BENCH_KERNEL_ALL_REPEATS(KER, AVX, MODEL, DATA) \
   BENCH_KERNEL_ALL_REPEATS(KER, SSE, MODEL, DATA)
 
#define BENCH_KERNEL_ALL_DATA(KER) \
  BENCH_KERNEL_ALL_VEC(KER, LG+G, 140) \
  BENCH_KERNEL_ALL_VEC(KER, GTR+G, 128)

#define BENCH_KERNEL_ALL() \
  BENCH_KERNEL_ALL_DATA(partial) \
  BENCH_KERNEL_ALL_DATA(likelihood)



#define BENCH_SEARCH(VEC, REPEATS, MODEL, DATA)\
  BENCHMARK_CAPTURE(BM_search,\
      SEARCH_ ## DATA  ## _ ## MODEL ## _ ## VEC ## _Repeat ## REPEATS ,\
      GET_TREE(DATA), \
      GET_MSA(DATA),\
      STRINGIFY(MODEL),\
      REPEATS,\
      CORAX_ATTRIB_ARCH_ ## VEC \
      );




BENCH_KERNEL_ALL();

BENCH_TREE(128)
BENCH_TREE(286)
BENCH_TREE(10575)

BENCH_SEARCH(AVX, true, GTR+G, 128);
BENCH_SEARCH(AVX, false, GTR+G, 128);
BENCH_SEARCH(SSE, true, GTR+G, 128);
BENCH_SEARCH(SSE, false, GTR+G, 128);

BENCHMARK_MAIN();

Tutorial: How to run an SPR round on a tree.
===============================================================================

1. Read a MSA in
2. Set up the `corax_partition_t`.
    1. (Optional) compress the site patterns
    2. Create partition with appropriate variables.
    3. Set the tip states.
    4. Set initial model parameters
        1. Frequencies
        2. Substitution parameters
        3. `p_inv`
        4. Rate categories, and associated rates
3. Make an initial tree
4. Make a `corax_treeinfo_t`
5. Initialize the partitions for the treeinfo struct
6. Run an SPR round

## 1: Read MSA

There is a data structure to hold an MSA. It has the following definition

```
struct corax_msa_s{
    int count;
    int length;

    char** sequence;
    char** label;
}corax_msa_t;
```

Where the parameters are as follows:

- `count`: The number of taxa in the MSA.
- `length`: the length of each sequence in the alignment. More literally, the
    length of the character arrays indexed by `sequence`.
- `sequence`: The sequences. Not null terminated.
- `label`: sequence label. Is null terminated.

To produce this structure, `corax_phylip_parse_interleaved`, `corax_phylip_parse_sequential` or `corax_fasta_load` should be
used. These functions all take a filename and handle the parsing for you.

## 2: Create a Partition

A partition is a set of sites that all evolve using the same model. This means that the partition structure in `coraxlib`
contains all the parameters _except_ for the tree and branch lengths. For many analyses, a single partition will contain
all the sites in the alignment, but in others multiple will be required. Therefore, it is always important to write
whatever tool to support multiple partitions.

To accomplish steps 2.1 to 2.4, I would recommend reading the [overview](corax_partition_t.md) of the partition
struct,and the [function](@ref corax_partition_t) documentation.

### 3: Make an initial tree

There are two easy methods to create an initial tree. The first is to make a random tree using

```
corax_utree_t* coraxmod_utree_create_random(unsigned int taxa_count,
                           const char* const* names,
                           unsigned int random_seed);
```

This will create a bifurcated tree at random. It does this by creating a minimal 3 tip tree, selecting a branch at
random, and then inserting a new tip off of that. The algorithm repeats this until all the tips have been added to the
tree.

- `taxa_count`: The number of tips that the resulting tree will have.
- `names`: The labels that will be assigned to the tips. Represented as null
    terminated `char*`.
- `random_seed`: The seed to pass to the random number generator.

The other method is to use a parsimony tree. This can be created using

```
corax_utree_t* coraxmod_utree_create_parsimony(unsigned int taxon_count,
                                           unsigned int seq_length,
                                           char* const* names,
                                           char* const* sequences,
                                           const unsigned int* site_weights,
                                           const corax_state_t* map,
                                           unsigned int states,
                                           unsigned int attributes,
                                           unsigned int random_seed,
                                           unsigned int* score);
```

This will create a tree under a parsimony method. The arguments are

- `taxon_count`: The number of tips on the tree.
- `seq_length`: The length of the sequences in the alignment.
- `names`: An array with length `taxon_count` of `char*` null terminated
    strings.
- `sequences`: An array with length `taxon_count` of `char*` that all have
    length `seq_length`.
- `site_weights`: An array that is `seq_length` long of site weights. Can be
    `nullptr` which means that all sites have equal weight.
- `map`: A predefined map from `char` to `int`. The choices are:
    - `corax_map_bin`: For an alphabet of ${0,1}$.
    - `corax_map_nt`: For nucleotide data.
    - `corax_map_aa`: For amino acid data.
- `states`: Number of states in the data.
- `attributes`: Exactly one of the following.
    - `PLL_ATTRIB_ARCH_CPU`: No special extensions,
    - `PLL_ATTRIB_ARCH_SSE`: Use the SSE3 extensions,
    - `PLL_ATTRIB_ARCH_AVX`: Use the AVX extensions,
    - `PLL_ATTRIB_ARCH_AVX2`: Use the AVX2 extensions,
    - `PLL_ATTRIB_ARCH_AVX512`: Use the AVX512 extensions,
    In addition, set the `PLL_ATTRIB_SITE_REPEATS` to enable site repeats.
- `random_seed`: Seed to pass to the random number generator
- `score`: output parameter that contains the score of the tree.

### 4: Make a Treeinfo Struct

Use the function `corax_treeinfo_create()` to create a `corax_treeinfo_t`.

### 5: Initialize the partitions 

Use the function `corax_treeinfo_init_partition()` to initialize a partition. For a basic, single partition inference
with no rate categories, call this function like so:

```{.c}
corax_treeinfo_init_partition(my_treeinfo,            /* pointer to treeinfo */
                              0,                      /* partition index */
                              my_partition,           /* pointer to an initialized partition */
                              CORAX_OPT_PARAM_ALL,    /* definition that controls which parameters to optimize */
                              CORAX_GAMMA_RATES_MEAN, /* definition that controls which gamma mode to use */
                              1.0,                    /* Initial alpha */
                              0,                      /* We only have one set of parameter indicies */
                              my_subst_syms,          /* A list of symmetries for the substitution matrix */
                              )
```


### 6: Perform an SPR round

The relevant function is `corax_algo_spr_round()`, and the documentation there is much more extensive. However, if all
that is required is a function call to start with, here is a quick example:

```{.c}
corax_algo_spr_round(my_treeinfo, 
                     1,    /* Consider reinsertion points at least 1 branch away from the prune point */
                     10,   /* Consider reinsertion points at most 10 branches away from the prune point */
                     5,    /* ntopol_keep, the number of topologies to keep while performing the SPR round */
                     true, /* Perform a "slow" mode round */
                     CORAX_OPT_MIN_BRANCH_LEN, /* Default minimum branch length */
                     CORAX_OPT_MAX_BRANCH_LEN, /* Default maxiumum branch length */
                     1,                  /* Number of BL optimization rounds */
                     CORAX_MISC_EPSILON, /* Default optimization threshold */
                     my_cutoff_info,     /* Pointer to some cutoff info */
                     1.0  /* Default subtree cutoff threshold */
                     )
```

The cutoff information should be initalized as:

```{.c}
cutoff_info_t* my_cutoff_info = (cutoff_info_t*)calloc(sizeof(cutoff_info_t), 1);
my_cutoff_info->lh_cutoff = 1.0;
```

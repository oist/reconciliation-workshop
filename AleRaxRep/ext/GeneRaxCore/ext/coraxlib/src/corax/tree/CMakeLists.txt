include(targetSourcesLocal)

target_sources_local(
    ${CORAX_LIB_TARGET}
    PRIVATE hashtable.c
            hashtable.h
            treeinfo.c
            treeinfo.h
            utree_compare.c
            utree_compare.h
            utree_moves.c
            utree_moves.h
            utree_ops.c
            utree_ops.h
            utree_parsimony.c
            utree_parsimony.h
            utree_random.c
            utree_random.h
            utree_split.c
            utree_split.h
            utree_tbe.c
            utree_tbe.h
            utree_traverse.c
            utree_traverse.h
            utree.c
            utree.h
)

#ifndef CORAX_CORE_PARSIMONY_H_
#define CORAX_CORE_PARSIMONY_H_

#include "corax/core/common.h"

/* structures for parsimony */

typedef struct corax_parsimony_s
{
  /* common information */
  unsigned int tips;
  unsigned int inner_nodes;
  unsigned int sites;
  unsigned int states;
  unsigned int attributes;
  size_t       alignment;

  /* fast unweighted parsimony */
  unsigned int **packedvector;
  unsigned int * node_cost;
  unsigned int   packedvector_count;
  unsigned int   const_cost;
  int *          informative;
  unsigned int   informative_count;

  /* weighted parsimony */
  unsigned int   score_buffers;
  unsigned int   ancestral_buffers;
  double *       score_matrix;
  double **      sbuffer;
  unsigned int **anc_states;
} corax_parsimony_t;

typedef struct corax_pars_buildop_s
{
  unsigned int parent_score_index;
  unsigned int child1_score_index;
  unsigned int child2_score_index;
} corax_pars_buildop_t;

typedef struct corax_pars_recop_s
{
  unsigned int node_score_index;
  unsigned int node_ancestral_index;
  unsigned int parent_score_index;
  unsigned int parent_ancestral_index;
} corax_pars_recop_t;


#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in parsimony.c */

  CORAX_EXPORT int corax_set_parsimony_sequence(corax_parsimony_t *  pars,
                                                unsigned int         tip_index,
                                                const corax_state_t *map,
                                                const char *         sequence);

  CORAX_EXPORT corax_parsimony_t *
               corax_parsimony_create(unsigned int  tips,
                                      unsigned int  states,
                                      unsigned int  sites,
                                      const double *score_matrix,
                                      unsigned int  score_buffers,
                                      unsigned int  ancestral_buffers);

  CORAX_EXPORT void corax_parsimony_destroy(corax_parsimony_t *pars);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_CORE_PARSIMONY_H_ */

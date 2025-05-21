#include "corax/corax.h"

CORAX_EXPORT int corax_update_prob_matrices(corax_partition_t * partition,
                                            const unsigned int *params_indices,
                                            const unsigned int *matrix_indices,
                                            const double *      branch_lengths,
                                            unsigned int        count)
{
  unsigned int n;

#ifdef CORAX_NONREV
  if (partition->attributes & CORAX_ATTRIB_NONREV)
  {
    return corax_core_update_pmatrix_nonrev(partition->pmatrix,
                                            partition->states,
                                            partition->rate_cats,
                                            partition->rates,
                                            branch_lengths,
                                            matrix_indices,
                                            params_indices,
                                            partition->prop_invar,
                                            partition->subst_params,
                                            count,
                                            partition->attributes);
  }
#endif

  /* check whether we have cached an eigen decomposition. If not, compute it */
  for (n = 0; n < partition->rate_cats; ++n)
  {
    if (!partition->eigen_decomp_valid[params_indices[n]])
    {
      if (!corax_update_eigen(partition, params_indices[n]))
        return CORAX_FAILURE;
    }
  }

  return corax_core_update_pmatrix(partition->pmatrix,
                                   partition->states,
                                   partition->rate_cats,
                                   partition->rates,
                                   branch_lengths,
                                   matrix_indices,
                                   params_indices,
                                   partition->prop_invar,
                                   partition->eigenvals,
                                   partition->eigenvecs,
                                   partition->inv_eigenvecs,
                                   count,
                                   partition->attributes);
}

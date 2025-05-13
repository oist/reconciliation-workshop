#include "opt_generic.h"

/******************************************************************************/
/* EXPECTATION-MAXIMIZATION (EM)     */
/* Wang, Li, Susko, and Roger (2008) */
/******************************************************************************/
CORAX_EXPORT void
corax_opt_minimize_em(double *             w,
                      unsigned int         w_count,
                      double *             sitecat_lh,
                      const unsigned int * site_w,
                      unsigned int         l,
                      void *               params,
                      double (*update_sitecatlk_funk)(void *, double *))
{
  unsigned int i, c;
  unsigned int max_steps   = 10;
  int          converged   = 0;
  int          ratio_scale = 0;

  double *new_prop   = (double *)malloc(sizeof(double) * w_count);
  double *ratio_prop = (double *)malloc(sizeof(double) * w_count);

  while (!converged && max_steps--)
  {
    /* update site-cat LK */
    update_sitecatlk_funk(params, sitecat_lh);

    // Expectation
    double *this_lk_cat = sitecat_lh;
    if (ratio_scale)
    {
      for (i = 0; i < l; ++i)
      {
        for (c = 0; c < w_count; c++) { this_lk_cat[c] *= ratio_prop[c]; }
        this_lk_cat += w_count;
      }
    }
    else
      ratio_scale = 1;

    memset(new_prop, 0, w_count * sizeof(double));

    this_lk_cat = sitecat_lh;
    for (i = 0; i < l; ++i)
    {
      // TODO: Check for p_invar
      double lk_ptn = 0;
      for (c = 0; c < w_count; c++) { lk_ptn += this_lk_cat[c]; }
      lk_ptn = site_w[i] / lk_ptn;
      for (c = 0; c < w_count; c++) { new_prop[c] += this_lk_cat[c] * lk_ptn; }
      this_lk_cat += w_count;
    }

    // Maximization
    converged = 1;
    for (c = 0; c < w_count; c++)
    {
      new_prop[c] /= l;

      // check for convergence
      converged     = converged && (fabs(w[c] - new_prop[c]) < 1e-4);
      ratio_prop[c] = new_prop[c] / w[c];
      w[c]          = new_prop[c];
    }
  }

  free(ratio_prop);
  free(new_prop);
}

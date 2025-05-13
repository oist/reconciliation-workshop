#include "opt_generic.h"

#include "lbfgsb/lbfgsb.h"

static inline int is_nan(double v) { return v != v; }

static inline int d_equals(double a, double b) { return (fabs(a - b) < 1e-10); }

/******************************************************************************/
/* L-BFGS-B OPTIMIZATION */
/******************************************************************************/

/**
 * Minimize a multi-parameter function using L-BFGS-B.
 *
 * The target function must compute the score at a certain state,
 * and it requires 2 parameters: (1) custom data (if needed),
 * and (2) the values at which score is computed.
 *
 * @param  x[in,out]   first guess and result of the minimization process
 * @param  xmin        lower bound for each of the variables
 * @param  xmax        upper bound for each of the variables
 * @param  bound       bound type (CORAX_LBFGSB_BOUND_[NONE|LOWER|UPPER|BOTH]
 * @param  n           number of variables
 * @param  factr       convergence tolerance for L-BFGS-B relative to machine
 * epsilon
 * @param  pgtol       absolute gradient tolerance for L-BFGS-B
 * @param  params      custom parameters required by the target function
 * @param  target_funk target function
 *
 * `factr` is a double precision variable. The iteration will stop when
 * (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
 * where epsmch is the machine epsilon
 *
 * `pgtol` is a double precision variable. The iteration will stop when
 * max{|proj g_i | i = 1, ..., n} <= pgtol
 * where pg_i is the ith component of the projected gradient.
 *
 * @return             the minimal score found
 */
CORAX_EXPORT double corax_opt_minimize_lbfgsb(double *     x,
                                              double *     xmin,
                                              double *     xmax,
                                              int *        bound,
                                              unsigned int n,
                                              double       factr,
                                              double       pgtol,
                                              void *       params,
                                              double (*target_funk)(void *,
                                                                    double *))
{
  unsigned int i;

  /* L-BFGS-B parameters */
  //  double initial_score;
  int     max_corrections;
  double  score = 0;
  double *g, *wa;
  int *   iwa;

  int  taskValue;
  int *task = &taskValue;

  int     csaveValue;
  int *   csave = &csaveValue;
  double  dsave[29];
  int     isave[44];
  logical lsave[4];

  int iprint = -1;

  max_corrections = 5;

  /* reset errno */
  corax_errno = 0;

  g = (double *)calloc((size_t)n, sizeof(double));

  /*     We start the iteration by initializing task. */
  *task = (int)START;

  iwa = (int *)calloc(3 * (size_t)n, sizeof(int));
  wa  = (double *)calloc((2 * (size_t)max_corrections + 5) * (size_t)n
                            + 12 * (size_t)max_corrections
                                  * ((size_t)max_corrections + 1),
                        sizeof(double));

  if (!(wa && iwa && g))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for l-bfgs-b variables");
    if (g) free(g);
    if (iwa) free(iwa);
    if (wa) free(wa);
    return (double)-INFINITY;
  }

  double lbfgsb_error = CORAX_LBFGSB_ERROR;
  for (unsigned int i = 0; i < n; i++)
    if (xmin[i] < lbfgsb_error) lbfgsb_error = xmin[i];

  //  double initial_score = target_funk (params, x);
  int continue_opt = 1;
  while (continue_opt)
  {
    /*     This is the call to the L-BFGS-B code. */
    setulb((int *)&n,
           &max_corrections,
           x,
           xmin,
           xmax,
           bound,
           &score,
           g,
           &factr,
           &pgtol,
           wa,
           iwa,
           task,
           &iprint,
           csave,
           lsave,
           isave,
           dsave);
    if (IS_FG(*task))
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */

      score = target_funk(params, x);

      if (is_nan(score) || d_equals(score, (double)-INFINITY)) break;

      double h, temp;
      for (i = 0; i < n; i++)
      {
        temp = x[i];
        h    = lbfgsb_error * fabs(temp);
        if (h < 1e-12) h = lbfgsb_error;

        x[i]           = temp + h;
        h              = x[i] - temp;
        double lnderiv = target_funk(params, x);

        g[i] = (lnderiv - score) / h;

        /* reset variable */
        x[i] = temp;
      }
    }
    else if (*task != NEW_X)
      continue_opt = 0;
  }

  /* fix optimal parameters */
  score = target_funk(params, x);

  free(iwa);
  free(wa);
  free(g);

  if (is_nan(score))
  {
    score = (double)-INFINITY;
    /* set errno only if it was not set by some inner function */
    if (!corax_errno)
    {
      corax_set_error(CORAX_OPT_ERROR_LBFGSB_UNKNOWN, "Unknown LBFGSB error");
    }
  }

  return score;
} /* corax_opt_minimize_lbfgsb */

struct bfgs_multi_opt
{
  unsigned int n;
  double *     x;
  double *     xmin;
  double *     xmax;
  int *        bound;
  double       factr;
  double       pgtol;

  int iprint;
  int max_corrections;
  int task;

  double *g, *wa;
  int *   iwa;
  int     csave;
  double  dsave[29];
  int     isave[44];
  logical lsave[4];

  double score;
  double h, temp;
};

static int init_bfgs_opt(struct bfgs_multi_opt *opt,
                         unsigned int           n,
                         double *               x,
                         double *               xmin,
                         double *               xmax,
                         int *                  bound,
                         double                 factr,
                         double                 pgtol)
{
  /* We start the iteration by initializing task. */
  opt->task  = (int)START;
  opt->score = 0;

  // some magic numbers
  opt->iprint          = -1;
  opt->max_corrections = 5;

  opt->n     = n;
  opt->x     = x;
  opt->xmin  = xmin;
  opt->xmax  = xmax;
  opt->bound = bound;
  opt->factr = factr;
  opt->pgtol = pgtol;

  /* allocate memory */
  opt->g = (double *)calloc((size_t)n, sizeof(double));

  opt->iwa = (int *)calloc(3 * (size_t)n, sizeof(int));
  opt->wa  = (double *)calloc((2 * (size_t)opt->max_corrections + 5) * (size_t)n
                                 + 12 * (size_t)opt->max_corrections
                                       * ((size_t)opt->max_corrections + 1),
                             sizeof(double));

  if (!opt->g || !opt->iwa || !opt->wa)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for l-bfgs-b variables");
    return CORAX_FAILURE;
  }

  return CORAX_SUCCESS;
}

static void destroy_bfgs_opt(struct bfgs_multi_opt *opt)
{
  if (opt)
  {
    if (opt->g) free(opt->g);
    if (opt->iwa) free(opt->iwa);
    if (opt->wa) free(opt->wa);
    free(opt);
  }
}

static int setulb_multi(struct bfgs_multi_opt *opt)
{
  return setulb((int *)&opt->n,
                &opt->max_corrections,
                opt->x,
                opt->xmin,
                opt->xmax,
                opt->bound,
                &opt->score,
                opt->g,
                &opt->factr,
                &opt->pgtol,
                opt->wa,
                opt->iwa,
                &opt->task,
                &opt->iprint,
                &opt->csave,
                opt->lsave,
                opt->isave,
                opt->dsave);
}

CORAX_EXPORT double corax_opt_minimize_lbfgsb_multi(
    unsigned int  xnum,
    double **     x,
    double **     xmin,
    double **     xmax,
    int **        bound,
    unsigned int *n,
    unsigned int  nmax,
    double        factr,
    double        pgtol,
    void *        params,
    double (*target_funk)(void *, double **, double *, int *))
{
  unsigned int i, p;

  double score = (double)-INFINITY;

  double *lh_old    = (double *)calloc((size_t)xnum, sizeof(double));
  double *lh_new    = (double *)calloc((size_t)xnum, sizeof(double));
  int *   converged = (int *)calloc((size_t)xnum + 1, sizeof(int));
  int *   skip      = (int *)calloc((size_t)xnum + 1, sizeof(int));

  struct bfgs_multi_opt **opts = (struct bfgs_multi_opt **)calloc(
      (size_t)xnum, sizeof(struct bfgs_multi_opt *));

  if (!lh_old || !lh_new || !converged || !skip || !opts)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for l-bfgs-b variables");
    goto cleanup;
  }

  for (p = 0; p < xnum; p++)
  {
    if (!x[p])
    {
      /* this is a remote partition - forget about it */
      converged[p] = skip[p] = 1;
      continue;
    }

    opts[p] = (struct bfgs_multi_opt *)calloc(1, sizeof(struct bfgs_multi_opt));
    if (!opts[p])
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for l-bfgs-b variables");
      goto cleanup;
    }

    if (!init_bfgs_opt(
            opts[p], n[p], x[p], xmin[p], xmax[p], bound[p], factr, pgtol))
      goto cleanup;
  }

  /* reset errno */
  corax_errno = 0;

  int continue_opt = 1;
  while (continue_opt)
  {
    continue_opt = 0;
    int all_skip = 1;
    for (p = 0; p < xnum; p++)
    {
      if (converged[p])
      {
        skip[p] = 1;
        continue;
      }

      assert(opts[p]);

      /*     This is the call to the L-BFGS-B code. */
      setulb_multi(opts[p]);

      skip[p]      = !IS_FG(opts[p]->task);
      converged[p] = skip[p] && (opts[p]->task != NEW_X);
      all_skip &= skip[p];
    }

    /* check if ALL partitions have converged, including remote ones */
    target_funk(params, NULL, NULL, converged);
    continue_opt = !converged[xnum];
    target_funk(params, NULL, NULL, skip);
    all_skip = converged[xnum];

    if (!all_skip && continue_opt)
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */
      score = target_funk(params, x, lh_old, skip);

      if (is_nan(score) || d_equals(score, (double)-INFINITY)) break;

      for (p = 0; p < xnum; p++)
      {
        if (!skip[p]) opts[p]->score = lh_old[p];
      }

      for (i = 0; i < nmax; i++)
      {
        for (p = 0; p < xnum; p++)
        {
          if (!skip[p] && i >= n[p])
          {
            /* this partition has less parameters than max -> skip it */
            skip[p] = 1;
          }

          if (skip[p]) continue;

          opts[p]->temp = x[p][i];
          opts[p]->h    = CORAX_LBFGSB_ERROR * fabs(opts[p]->temp);
          if (opts[p]->h < 1e-12) opts[p]->h = CORAX_LBFGSB_ERROR;

          x[p][i]    = opts[p]->temp + opts[p]->h;
          opts[p]->h = x[p][i] - opts[p]->temp;
        }

        score = target_funk(params, x, lh_new, skip);

        for (p = 0; p < xnum; p++)
        {
          if (skip[p]) continue;

          assert(opts[p]);

          /* compute partial derivative */
          opts[p]->g[i] = (lh_new[p] - lh_old[p]) / opts[p]->h;

          /* reset variable */
          x[p][i] = opts[p]->temp;

          DBG("P%u  lh_old: %f, lh_new: %f\n", p, lh_old[p], lh_new[p]);
        }
      }
    }
  }

  /* fix optimal parameters */
  score = target_funk(params, x, NULL, NULL);

cleanup:
  if (lh_old) free(lh_old);
  if (lh_new) free(lh_new);
  if (converged) free(converged);
  if (skip) free(skip);
  if (opts)
  {
    for (p = 0; p < xnum; p++) destroy_bfgs_opt(opts[p]);
    free(opts);
  }

  if (is_nan(score))
  {
    score = (double)-INFINITY;
    /* set errno only if it was not set by some inner function */
    if (!corax_errno)
    {
      corax_set_error(CORAX_OPT_ERROR_LBFGSB_UNKNOWN, "Unknown LBFGSB error");
    }
  }

  return score;
} /* corax_opt_minimize_lbfgsb */

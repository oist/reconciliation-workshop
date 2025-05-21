#include "opt_generic.h"

/******************************************************************************/
/* NEWTON-RAPHSON OPTIMIZATION */
/******************************************************************************/
typedef struct
{
  void (*derive_func)(void *, double, double *, double *);
  void *params;
} newton_wrapper_params;

static void
newton_wrapper_func(void *params, double *proposal, double *df, double *ddf)
{
  newton_wrapper_params *wrapper_params = (newton_wrapper_params *)params;

  wrapper_params->derive_func(wrapper_params->params, proposal[0], df, ddf);
}

/**
 * Minimize a function using Newton-Raphson.
 *
 * This is ultimately a wrapper around `corax_opt_minimize_newton_multi`. Please
 * see that function for a more detailed description of the `deriv_func`.
 *
 * @param  x1         lower bound
 * @param  xguess     first guess for the free variable
 * @param  x2         upper bound
 * @param  tolerance  tolerance of the minimization method
 * @param  max_iters  maximum number of iterations (bounds the effect of slow
 * convergence)
 * @param  params     custom parameters required by the target function
 * @param  deriv_func derivative function
 *
 * @return            the parameter value that minimizes the function in [x1,x2]
 */
CORAX_EXPORT double corax_opt_minimize_newton(
    double       xmin,
    double       xguess,
    double       xmax,
    double       tolerance,
    unsigned int max_iters,
    void *       params,
    void (*deriv_func)(void *, double, double *, double *))
{
  newton_wrapper_params wrapper_params;
  wrapper_params.params      = params;
  wrapper_params.derive_func = deriv_func;

  double xres = xguess;

  int retval = corax_opt_minimize_newton_multi(1,
                                               xmin,
                                               &xres,
                                               xmax,
                                               tolerance,
                                               max_iters,
                                               NULL,
                                               (void *)&wrapper_params,
                                               newton_wrapper_func);

  if (retval)
    return xres;
  else
    return CORAX_FAILURE;
}

/**
 * Minimize multiple functions in parallel using Newton-Raphson.
 * (e.g. unlinked branch lengths for a partitioned alignment)
 *
 * The derivative function must compute _both_ of the first and second
 * derivatives at a given point, and it requires 4 parameters: (1) custom data
 * (if needed), (2) the value at which derivatives are to be computed. The
 * remaining 2 parameters (3, 4) are out parameters for the first and second
 * derivatives.
 *
 * @param  xnum       number of functions/variables to optimize
 * @param  xmin       lower bound
 * @param  xguess     in=first guess for the free variables, out=optimized
 * values
 * @param  xmax       upper bound
 * @param  tolerance  tolerance of the minimization method
 * @param  max_iters  maximum number of iterations (bounds the effect of slow
 * convergence)
 * @param  converged  (optional) in: 1/0=skip/optimize, out: 1/0=converged
 * yes/no
 * @param  params     custom parameters required by the derivative function
 * @param  deriv_func derivative function
 *
 * @return            CORAX_FAILURE on error, CORAX_SUCCESS otherwise
 */
CORAX_EXPORT int corax_opt_minimize_newton_multi(
    unsigned int xnum,
    double       xmin,
    double *     xguess,
    double       xmax,
    double       tolerance,
    unsigned int max_iters,
    int *        converged,
    void *       params,
    void(deriv_func)(void *, double *, double *, double *))
{
  unsigned int i;
  unsigned int iter          = 0;
  int          all_converged = 0;
  int          error_flag    = 0;

  double dxmax = xmax / max_iters;

  double *xl            = (double *)calloc(xnum, sizeof(double));
  double *xh            = (double *)calloc(xnum, sizeof(double));
  double *f             = (double *)calloc(xnum, sizeof(double));
  double *df            = (double *)calloc(xnum, sizeof(double));
  int *   int_converged = NULL;
  double *x             = xguess;

  /* reset errno */
  corax_errno = 0;

  if (!converged)
  {
    int_converged = (int *)calloc(xnum, sizeof(int));
    converged     = int_converged;
  }

  for (i = 0; i < xnum; i++)
  {
    x[i] = CORAX_MAX(CORAX_MIN(x[i], xmax), xmin);

    xl[i] = xmin;
    xh[i] = xmax;
  }

  while (!all_converged && !error_flag)
  {
    if (iter++ > max_iters)
    {
      corax_set_error(CORAX_OPT_ERROR_NEWTON_LIMIT,
                      "Exceeded maximum number of iterations");
      error_flag = 1;
      break;
    }

    /* compute derivative for *all* functions in parallel */
    deriv_func((void *)params, x, f, df);

    /* for every function: check convergence and make the next guess */
    all_converged = 1;
    for (i = 0; i < xnum; i++)
    {
      double dx;

      if (converged[i]) continue;

      if (!isfinite(f[i]) || !isfinite(df[i]))
      {
        DBG("[it=%u][NR deriv][p=%u] BL=%.9f   f=%.12f  df=%.12f\n",
            iter,
            i,
            x[i],
            f[i],
            df[i]);
        corax_set_error(CORAX_OPT_ERROR_NEWTON_DERIV,
                        "Wrong likelihood derivatives");
        error_flag = 1;
        break;
      }

      if (df[i] > 0.0)
      {
        if (fabs(f[i]) < tolerance)
        {
          converged[i] = 1;
          continue;
        }

        if (f[i] < 0.0)
          xl[i] = x[i];
        else
          xh[i] = x[i];

        dx = -1 * f[i] / df[i];
      }
      else
      {
        dx = -1 * f[i] / fabs(df[i]);
      }

      dx = CORAX_MAX(CORAX_MIN(dx, dxmax), -dxmax);

      if (x[i] + dx < xl[i]) dx = xl[i] - x[i];
      if (x[i] + dx > xh[i]) dx = xh[i] - x[i];

      if (fabs(dx) < tolerance)
      {
        converged[i] = 1;
        continue;
      }

      DBG("[it=%u][NR deriv][p=%u] BL=%.9f   f=%.12f  df=%.12f  nextBL=%.9f\n",
          iter,
          i,
          x[i],
          f[i],
          df[i],
          x[i] + dx);

      x[i] += dx;

      x[i] = CORAX_MAX(CORAX_MIN(x[i], xmax), xmin);

      all_converged &= converged[i];
    }
  }

  free(xl);
  free(xh);
  free(f);
  free(df);
  if (int_converged) free(int_converged);

  if (all_converged && !error_flag)
    return CORAX_SUCCESS;
  else
    return CORAX_FAILURE;
}

/******************************************************************************/
/* BRENT'S OPTIMIZATION */
/******************************************************************************/
/* most of the code for Brent optimization taken from IQ-Tree
 * http://www.cibiv.at/software/iqtree
 * --------------------------------------------------------------------------
 * Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013)
 * Ultrafast approximation for phylogenetic bootstrap.
 * Mol. Biol. Evol., 30:1188-1195. (free reprint, DOI: 10.1093/molbev/mst024) */

#include "opt_generic.h"

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-7
#define SHFT(a, b, c, d)                                                       \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

typedef struct
{
  double startx;
  double fstartx;

  double  tol;
  double *foptx;
  double *f2optx;

  int    iter;
  double a;
  double b;
  double d;
  double etemp;
  double fu;
  double fv;
  double fw;
  double fx;
  double p;
  double q;
  double r;
  double tol1;
  double tol2;
  double u;
  double v;
  double w;
  double x;
  double xm;
  double xw;
  double wv;
  double vx;
  double e;
} opt_params;

typedef struct
{
  double (*target_funk)(void *, double);
  void *params;
} brent_wrapper_params;

static int brent_opt_pre_loop(opt_params *bp)
{
  double etemp;

  bp->xm   = 0.5 * (bp->a + bp->b);
  bp->tol2 = 2.0 * (bp->tol1 = bp->tol * fabs(bp->x) + ZEPS);
  if (fabs(bp->x - bp->xm) <= (bp->tol2 - 0.5 * (bp->b - bp->a)))
  {
    if (bp->foptx) *bp->foptx = bp->fx;
    bp->xw = bp->x - bp->w;
    bp->wv = bp->w - bp->v;
    bp->vx = bp->v - bp->x;
    if (bp->f2optx)
    {
      *bp->f2optx = 2.0 * (bp->fv * bp->xw + bp->fx * bp->wv + bp->fw * bp->vx)
                    / (bp->v * bp->v * bp->xw + bp->x * bp->x * bp->wv
                       + bp->w * bp->w * bp->vx);
    }
    return CORAX_FAILURE;
  }

  if (fabs(bp->e) > bp->tol1)
  {
    bp->r = (bp->x - bp->w) * (bp->fx - bp->fv);
    bp->q = (bp->x - bp->v) * (bp->fx - bp->fw);
    bp->p = (bp->x - bp->v) * bp->q - (bp->x - bp->w) * bp->r;
    bp->q = 2.0 * (bp->q - bp->r);
    if (bp->q > 0.0) bp->p = -bp->p;
    bp->q = fabs(bp->q);
    etemp = bp->e;
    bp->e = bp->d;
    if (fabs(bp->p) >= fabs(0.5 * bp->q * etemp)
        || bp->p <= bp->q * (bp->a - bp->x) || bp->p >= bp->q * (bp->b - bp->x))
      bp->d =
          CGOLD * (bp->e = (bp->x >= bp->xm ? bp->a - bp->x : bp->b - bp->x));
    else
    {
      bp->d = bp->p / bp->q;
      bp->u = bp->x + bp->d;
      if (bp->u - bp->a < bp->tol2 || bp->b - bp->u < bp->tol2)
        bp->d = SIGN(bp->tol1, bp->xm - bp->x);
    }
  }
  else
  {
    bp->d = CGOLD * (bp->e = (bp->x >= bp->xm ? bp->a - bp->x : bp->b - bp->x));
  }

  bp->u =
      (fabs(bp->d) >= bp->tol1 ? bp->x + bp->d : bp->x + SIGN(bp->tol1, bp->d));

  return CORAX_SUCCESS;
}

static int brent_opt_init(double      ax,
                          double      bx,
                          double      cx,
                          double      tol,
                          double *    foptx,
                          double *    f2optx,
                          double      fax,
                          double      fbx,
                          double      fcx,
                          opt_params *bp)
{
  memset(bp, 0, sizeof(opt_params));

  bp->tol    = tol;
  bp->foptx  = foptx;
  bp->f2optx = f2optx;

  bp->a      = (ax < cx ? ax : cx);
  bp->b      = (ax > cx ? ax : cx);
  bp->startx = bp->x = bx;
  bp->fstartx = bp->fx = fbx;
  if (fax < fcx)
  {
    bp->w  = ax;
    bp->fw = fax;
    bp->v  = cx;
    bp->fv = fcx;
  }
  else
  {
    bp->w  = cx;
    bp->fw = fcx;
    bp->v  = ax;
    bp->fv = fax;
  }

  /* pre-loop iteration 0 */
  bp->iter = 1;

  return brent_opt_pre_loop(bp);
}

static int brent_opt_post_loop(opt_params *bp)
{
  /* post-loop iteration i */

  if (bp->fu <= bp->fx)
  {
    if (bp->u >= bp->x)
      bp->a = bp->x;
    else
      bp->b = bp->x;

    SHFT(bp->v, bp->w, bp->x, bp->u)
    SHFT(bp->fv, bp->fw, bp->fx, bp->fu)
  }
  else
  {
    if (bp->u < bp->x)
      bp->a = bp->u;
    else
      bp->b = bp->u;
    if (bp->fu <= bp->fw || bp->w == bp->x)
    {
      bp->v  = bp->w;
      bp->w  = bp->u;
      bp->fv = bp->fw;
      bp->fw = bp->fu;
    }
    else if (bp->fu <= bp->fv || bp->v == bp->x || bp->v == bp->w)
    {
      bp->v  = bp->u;
      bp->fv = bp->fu;
    }
  }

  /* pre-loop iteration i+1 */

  ++bp->iter;

  return brent_opt_pre_loop(bp);
}

double
target_funk_wrapper(void *params, double *xopt, double *fxopt, int *converged)
{
  CORAX_UNUSED(converged);
  brent_wrapper_params *wrap_params = (brent_wrapper_params *)params;

  *fxopt = wrap_params->target_funk(wrap_params->params, *xopt);

  return *fxopt;
}

static int
brent_opt_alt(unsigned int xnum,
              int *        opt_mask,
              double *     xmin,
              double *     xguess,
              double *     xmax,
              double       xtol,
              double *     xopt,
              double *     fx,
              double *     f2x,
              void *       params,
              double (*target_funk)(void *, double *, double *, int *),
              int global_range)
{
  opt_params *brent_params = (opt_params *)calloc(xnum, sizeof(opt_params));

  double *ax        = (double *)calloc(xnum, sizeof(double));
  double *cx        = (double *)calloc(xnum, sizeof(double));
  double *fa        = (double *)calloc(xnum, sizeof(double));
  double *fb        = (double *)calloc(xnum, sizeof(double));
  double *fc        = (double *)calloc(xnum, sizeof(double));
  double *fxmin     = (double *)calloc(xnum, sizeof(double));
  double *fxmax     = (double *)calloc(xnum, sizeof(double));
  int *   converged = (int *)calloc(xnum + 1, sizeof(int));

  double *l_xmin = NULL;
  double *l_xmax = NULL;

  unsigned int i;
  int          iterate = 1;

  if (global_range)
  {
    l_xmin = (double *)calloc(xnum, sizeof(double));
    l_xmax = (double *)calloc(xnum, sizeof(double));
    for (i = 0; i < xnum; ++i)
    {
      l_xmin[i] = *xmin;
      l_xmax[i] = *xmax;
    }
  }
  else
  {
    l_xmin = xmin;
    l_xmax = xmax;
  }

  /* this function is a refactored version of brent_opt */
  /* if we consider the following structure:
   *
   *    (a) initialization block
   *    (b) main loop:
   *      (b.1) loop pre-score
   *      (b.2) loop score (call to target)
   *      (b.3) loop post-score
   *
   * the loop has been removed and it has been split into 2 functions:
   *
   * 1. brent_opt_init, that covers (a) and (b.1)
   * 2. brent_opt_post_loop, that covers (b.3) and (b.1)
   *
   * This way, brent_opt can be refactored as follows:
   *
   * (A) brent_opt_init
   * (B) main loop
   *    (B.1) loop score (call to target)
   *    (B.2) brent_opt_post_loop
   *
   * For parallel execution, the synchronization point occurs right at the
   * beginning of the loop.
   *
   * Moreover, the optimization state is encapsulated in a `struct opt_params *`
   * that can be defined as an array holding the optimization state for each
   * local partition.
   *
   * I observed no impact in results nor in runtime
   */

  for (i = 0; i < xnum; ++i)
  {
    double eps;
    int    outbounds_ax, outbounds_cx;

    /* skip params that do not need to be optimized */
    if (opt_mask && !opt_mask[i]) continue;

    /* first attempt to bracketize minimum */
    if (xguess[i] < l_xmin[i]) xguess[i] = l_xmin[i];
    if (xguess[i] > l_xmax[i]) xguess[i] = l_xmax[i];

    // TODO: this is a quick workaround to make it work for xguess==0
    // But we should double-check this bracketing heuristic! (alexey)
    eps = xguess[i] > 0 ? xguess[i] * xtol * 50.0 : 2. * xtol;

    ax[i]        = xguess[i] - eps;
    outbounds_ax = ax[i] < l_xmin[i];
    if (outbounds_ax) ax[i] = l_xmin[i];
    cx[i]        = xguess[i] + eps;
    outbounds_cx = cx[i] > l_xmax[i];
    if (outbounds_cx) cx[i] = l_xmax[i];
  }

  target_funk(params, ax, fa, NULL);
  target_funk(params, xguess, fb, NULL);
  target_funk(params, cx, fc, NULL);
  target_funk(params, l_xmin, fxmin, NULL);
  target_funk(params, l_xmax, fxmax, NULL);

  /* check if this works */
  for (i = 0; i < xnum; ++i)
  {
    /* skip params that do not need to be optimized */
    if (opt_mask && !opt_mask[i]) continue;

    /* if it works use these borders else be conservative */
    if ((fa[i] < fb[i]) || (fc[i] < fb[i]))
    {
      fa[i] = fxmin[i];
      fc[i] = fxmax[i];
      ax[i] = l_xmin[i];
      cx[i] = l_xmax[i];
    }

    DBG("param[%d] a x c / fa fx fc: %lf, %lf, %lf / %lf, %lf, %lf / %lf\n",
        i,
        ax[i],
        xguess[i],
        cx[i],
        fa[i],
        fb[i],
        fc[i],
        l_xmax[i]);

    if (!brent_opt_init(ax[i],
                        xguess[i],
                        cx[i],
                        xtol,
                        fx,
                        f2x,
                        fa[i],
                        fb[i],
                        fc[i],
                        &brent_params[i]))
    {
      DBG("BRENT: converged for partition: %d \n", i);
      converged[i] = 1;
    }
  }

  free(ax);
  free(cx);
  free(fa);
  free(fb);
  free(fc);
  free(fxmin);
  free(fxmax);

  if (global_range)
  {
    free(l_xmin);
    free(l_xmax);
  }

  double *u  = (double *)calloc(xnum, sizeof(double));
  double *fu = (double *)calloc(xnum, sizeof(double));

  int iter_num = 0;
  while (iterate)
  {
    for (i = 0; i < xnum; ++i) u[i] = brent_params[i].u;

    target_funk(params, u, fu, converged);

    DBG("iter: %d, u: %lf, fu: %lf\n", iter_num, u[2], fu[2]);

    /* last element in converged[] array is "all converged" flag */
    iterate = !converged[xnum];
    for (i = 0; i < xnum; ++i)
    {
      /* skip params that do not need to be optimized */
      if (opt_mask && !opt_mask[i]) continue;

      if (!converged[i])
      {
        brent_params[i].fu = fu[i];
        converged[i]       = !brent_opt_post_loop(&brent_params[i]);
      }
    }

    iter_num++;
    iterate &= (iter_num <= ITMAX);
  }

  /* if new score is worse, return initial value */
  for (i = 0; i < xnum; ++i)
  {
    xopt[i] = (brent_params[i].fx > brent_params[i].fstartx)
                  ? brent_params[i].startx
                  : brent_params[i].x;

    DBG("xopt_PRE: %lf, LH_PRE: %lf\n", xopt[i], brent_params[i].fx);
  }

  target_funk(params, xopt, fx, NULL);

  //  for (i = 0; i < xnum; ++i)
  //    DBG("xopt: %lf, LH: %lf\n", xopt[i], fx ? fx[i] : NAN);

  free(u);
  free(fu);
  free(converged);
  free(brent_params);
  return CORAX_SUCCESS;
}

/**
 * Minimize a single-variable function using Brent.
 *
 * The target function must evaluate the function at a certain point,
 * and it requires 2 parameters: (1) custom data (if needed),
 * and (2) the value where the function is evaluated.
 *
 * @param  xmin       lower bound
 * @param  xguess     first guess for the free variable
 * @param  xmax       upper bound
 * @param  xtol       tolerance of the minimization method
 * @param  params     custom parameters required by the target function
 * @param  deriv_func target function
 *
 * @return            the parameter value that minimizes the function in
 * [xmin,xmax]
 */
CORAX_EXPORT double corax_opt_minimize_brent(double  xmin,
                                             double  xguess,
                                             double  xmax,
                                             double  xtol,
                                             double *fx,
                                             double *f2x,
                                             void *  params,
                                             double (*target_funk)(void *,
                                                                   double))
{
  double optx = xguess;

  brent_wrapper_params wrap_params;
  wrap_params.target_funk = target_funk;
  wrap_params.params      = params;

  brent_opt_alt(1,
                NULL,
                &xmin,
                &xguess,
                &xmax,
                xtol,
                &optx,
                fx,
                f2x,
                &wrap_params,
                target_funk_wrapper,
                1);

  //  optx = brent_opt(xmin, xguess, xmax, xtol, fx, f2x, params, target_funk);

  return optx; /* return optimal x */
}

/**
 * Run brent optimization for multiple variables in parallel
 * (e.g., unlinked alphas for multiple partitions)
 *
 * @param xnum number of variables/partitions
 * @param opt_mask opt_mask[i]=1/0: optimize/skip parameter i
 * @param xmin minimum value(s) (see @param global_range)
 * @param xguess array of staring values
 * @param xmax maximum value(s) (see @param global_range)
 * @param xtol tolerance
 * @param xopt optimal variable values [out]
 * @param fx auxiliary variable to keep state between runs [out]
 * @param f2x auxiliary variable to keep state between runs [out]
 * @param params parameters to be passed to target_funk
 * @param target_funk target function, parameters: 1) params, 2) array of x
 * values, 3) array of scores [out], 4) convergence flags
 * @param global_range 0=xmin/xmax point to arrays of size xnum with individual
 * per-variable ranges; 1=xmin/xmax is a global range for all variables
 */
CORAX_EXPORT int corax_opt_minimize_brent_multi(
    unsigned int xnum,
    int *        opt_mask,
    double *     xmin,
    double *     xguess,
    double *     xmax,
    double       xtol,
    double *     xopt,
    double *     fx,
    double *     f2x,
    void *       params,
    double (*target_funk)(void *, double *, double *, int *),
    int global_range)
{
  return brent_opt_alt(xnum,
                       opt_mask,
                       xmin,
                       xguess,
                       xmax,
                       xtol,
                       xopt,
                       fx,
                       f2x,
                       params,
                       target_funk,
                       global_range);
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN

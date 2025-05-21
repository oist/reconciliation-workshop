/*
    Copyright (C) 2015-2021 Diego Darriba, Alexey Kozlov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#ifndef CORAX_OPTIMIZE_GENERIC_H_
#define CORAX_OPTIMIZE_GENERIC_H_

#include "corax/corax_core.h"
#include "corax/corax_tree.h"

#ifdef DEBUG
#define DBG(fmt, ...)                                                          \
  do {                                                                         \
    printf(fmt, ##__VA_ARGS__);                                                \
  } while (0)
#else
#define DBG(fmt, ...)
#endif

/** @ingroup corax_defines
 * @{
 */
// it's actually defined in lbfgsb.h, but not exported from the optimize module
#define CORAX_ALGO_LBFGSB_ERROR 1.0e-4

/* Parameters mask */
#define CORAX_OPT_PARAM_ALL (~0)
#define CORAX_OPT_PARAM_SUBST_RATES (1 << 0)
#define CORAX_OPT_PARAM_ALPHA (1 << 1)
#define CORAX_OPT_PARAM_PINV (1 << 2)
#define CORAX_OPT_PARAM_FREQUENCIES (1 << 3)
#define CORAX_OPT_PARAM_BRANCHES_SINGLE (1 << 4)
#define CORAX_OPT_PARAM_BRANCHES_ALL (1 << 5)
#define CORAX_OPT_PARAM_BRANCHES_ITERATIVE (1 << 6)
#define CORAX_OPT_PARAM_TOPOLOGY (1 << 7)
#define CORAX_OPT_PARAM_FREE_RATES (1 << 8)
#define CORAX_OPT_PARAM_RATE_WEIGHTS (1 << 9)
#define CORAX_OPT_PARAM_BRANCH_LEN_SCALER (1 << 10)
/* !!! NOTE: all params in user code must be defined as
 *           CORAX_OPT_PARAM_USER<<0, CORAX_OPT_PARAM_USER<<1 etc. !!! */
#define CORAX_OPT_PARAM_USER (1 << 16)

/* L-BFGS-B bound type */
#define CORAX_OPT_LBFGSB_BOUND_NONE 0
#define CORAX_OPT_LBFGSB_BOUND_LOWER 1
#define CORAX_OPT_LBFGSB_BOUND_BOTH 2
#define CORAX_OPT_LBFGSB_BOUND_UPPER 3

/* Parameter defaults */
#define CORAX_OPT_DEFAULT_RATE_RATIO 1
#define CORAX_OPT_DEFAULT_FREQ_RATIO 1
#define CORAX_OPT_DEFAULT_PINV 0.01
#define CORAX_OPT_DEFAULT_ALPHA 0.5
#define CORAX_OPT_DEFAULT_BRANCH_LEN 0.1
#define CORAX_OPT_DEFAULT_SMOOTHINGS 32
#define CORAX_OPT_DEFAULT_EPSILON 1.0e-1

/* Default parameter limits */
#define CORAX_OPT_MIN_BRANCH_LEN 1.0e-4
#define CORAX_OPT_MAX_BRANCH_LEN 100.
#define CORAX_OPT_TOL_BRANCH_LEN 1.0e-4
#define CORAX_OPT_MIN_SUBST_RATE 1.0e-3
#define CORAX_OPT_MAX_SUBST_RATE 1000.
#define CORAX_OPT_MIN_FREQ 1.0e-3
#define CORAX_OPT_MAX_FREQ 100.
#define CORAX_OPT_MIN_ALPHA 0.0201 //+ CORAX_LBFGSB_ERROR
#define CORAX_OPT_MAX_ALPHA 100.
#define CORAX_OPT_MIN_PINV 0
#define CORAX_OPT_MAX_PINV 0.99
#define CORAX_OPT_LNL_UNLIKELY -1e+80

/* mixture models limits */
#define CORAX_OPT_MIN_RATE 0.02
#define CORAX_OPT_MAX_RATE 100.
#define CORAX_OPT_MIN_RATE_WEIGHT 1.0e-3
#define CORAX_OPT_MAX_RATE_WEIGHT 100.

/* branch length optimization methods */
#define CORAX_OPT_BLO_NEWTON_FAST 0 /* standard Newton-Raphson (NR) */
#define CORAX_OPT_BLO_NEWTON_SAFE 1 /* NR with per-branch LH check */
#define CORAX_OPT_BLO_NEWTON_FALLBACK                                          \
  2 /* NR-FAST with fallback to NR-SAFE                                        \
     */
#define CORAX_OPT_BLO_NEWTON_GLOBAL                                            \
  3 /* NR variant which looks for local optima */

#define CORAX_OPT_BLO_NEWTON_OLDFAST 4 /* standard Newton-Raphson (NR) - old*/
#define CORAX_OPT_BLO_NEWTON_OLDSAFE 5 /* NR with per-branch LH check - old*/


/** @} */

/** @ingroup corax_errors
 * @{
 */
/* error codes (for this module, 2000-3000) */
#define CORAX_OPT_ERROR_PARAMETER 2000
#define CORAX_OPT_ERROR_TAXA_MISMATCH 2010
#define CORAX_OPT_ERROR_SEQLEN_MISMATCH 2020
#define CORAX_OPT_ERROR_ALIGN_UNREADABLE 2030
#define CORAX_OPT_ERROR_LBFGSB_UNKNOWN 2100
#define CORAX_OPT_ERROR_NEWTON_DERIV 2210
#define CORAX_OPT_ERROR_NEWTON_LIMIT 2220
#define CORAX_OPT_ERROR_NEWTON_UNKNOWN 2230
#define CORAX_OPT_ERROR_NEWTON_WORSE_LK 2240
#define CORAX_OPT_ERROR_NEWTON_BAD_RADIUS 2250
#define CORAX_OPT_ERROR_BRENT_INIT 2310
/** @} */

/* special options */
#define CORAX_OPT_BRLEN_OPTIMIZE_ALL -1

/* Structure with information necessary for evaluating the likelihood */

/* Custom parameters structures provided by CORAX for the
 * high level optimization functions (L-BFGS-B + Brent). */
typedef struct
{
  corax_partition_t  *partition;
  corax_operation_t  *operations;
  double             *branch_lengths;
  unsigned int       *matrix_indices;
  int                 rooted;
  const unsigned int *params_indices;
  union
  {
    struct
    {
      unsigned int root_clv_index;
      int          scaler_index;
    } rooted_t;
    struct
    {
      unsigned int parent_clv_index;
      int          parent_scaler_index;
      unsigned int child_clv_index;
      int          child_scaler_index;
      unsigned int edge_pmatrix_index;
    } unrooted_t;
  } where;

  char   __padding__[4];
  double alpha_value;
} corax_likelihood_info_t;

typedef struct
{
  corax_likelihood_info_t lk_params;
  unsigned int            highest_freq_state;
  unsigned int            highest_weight_state;
  // const unsigned int * params_indices;     /* indices according to rate cats
  // */
  unsigned int params_index; /* individual index to optimize */
  unsigned int which_parameters;
  int         *subst_params_symmetries;
  double       factr;
  double       pgtol;

  double *sumtable;
} corax_optimize_options_t;

/******************************************************************************/

#ifdef __cplusplus
extern "C"
{
#endif


  /* functions in newtom.c */

  /* core Newton-Raphson optimization function (multiple variables) */
  CORAX_EXPORT int corax_opt_minimize_newton_multi(
      unsigned int xnum,
      double       xmin,
      double      *xguess,
      double       xmax,
      double       tolerance,
      unsigned int max_iters,
      int         *converged,
      void        *params,
      void(deriv_func)(void *, double *, double *, double *));

  /* core Newton-Raphson optimization function */
  CORAX_EXPORT double
  corax_opt_minimize_newton(double       xmin,
                            double       xguess,
                            double       xmax,
                            double       tolerance,
                            unsigned int max_iters,
                            void        *params,
                            void(deriv_func)(void *, double, double *, double *));
  /* functions in bfgs.c */

  /* core L-BFGS-B optimization function */
  CORAX_EXPORT double corax_opt_minimize_lbfgsb(double      *x,
                                                double      *xmin,
                                                double      *xmax,
                                                int         *bound,
                                                unsigned int n,
                                                double       factr,
                                                double       pgtol,
                                                void        *params,
                                                double (*target_funk)(void *,
                                                                      double *));

  CORAX_EXPORT double corax_opt_minimize_lbfgsb_multi(
      unsigned int  xnum,
      double      **x,
      double      **xmin,
      double      **xmax,
      int         **bound,
      unsigned int *n,
      unsigned int  nmax,
      double        factr,
      double        pgtol,
      void         *params,
      double (*target_funk)(void *, double **, double *, int *));

  /* functions in brent.c */

  /* core Brent optimization function */
  CORAX_EXPORT double corax_opt_minimize_brent(double  xmin,
                                               double  xguess,
                                               double  xmax,
                                               double  xtol,
                                               double *fx,
                                               double *f2x,
                                               void   *params,
                                               double (*target_funk)(void *,
                                                                     double));

  CORAX_EXPORT int corax_opt_minimize_brent_multi(
      unsigned int xnum,
      int         *opt_mask,
      double      *xmin,
      double      *xguess,
      double      *xmax,
      double       xtol,
      double      *xopt,
      double      *fx,
      double      *f2x,
      void        *params,
      double (*target_funk)(void *, double *, double *, int *),
      int global_range);

  /* functions in em.c */

  /* core Expectation-Maximization (EM) function */
  CORAX_EXPORT void
  corax_opt_minimize_em(double             *w,
                        unsigned int        w_count,
                        double             *sitecat_lh,
                        const unsigned int *site_w,
                        unsigned int        l,
                        void               *params,
                        double (*update_sitecatlk_funk)(void *, double *));

  /* functions in opt_generic.c */

  CORAX_EXPORT double
  corax_opt_optimize_onedim(corax_optimize_options_t *p, double min, double max);

  CORAX_EXPORT double corax_opt_optimize_multidim(corax_optimize_options_t *p,
                                                  const double             *umin,
                                                  const double             *umax);

  CORAX_EXPORT double corax_opt_compute_lk(corax_partition_t  *partition,
                                           corax_unode_t      *tree,
                                           const unsigned int *params_indices,
                                           int                 update_pmatrices,
                                           int                 update_partials);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_OPTIMIZE_GENERIC_H_ */

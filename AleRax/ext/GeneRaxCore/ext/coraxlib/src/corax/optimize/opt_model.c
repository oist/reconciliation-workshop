/*
 Copyright (C) 2016 Diego Darriba, Alexey Kozlov

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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include "opt_model.h"
#include "callback.h"
#include "corax/corax.h"

static void fill_rates(double *     rates,
                       double *     x,
                       int *        bt,
                       double *     lb,
                       double *     ub,
                       double       min_rate,
                       double       max_rate,
                       unsigned int n_rates);

static void fill_weights(double *      weights,
                         unsigned int *fixed_weight_index,
                         double *      x,
                         int *         bt,
                         double *      lb,
                         double *      ub,
                         unsigned int  n_weights);

CORAX_EXPORT double corax_algo_opt_frequencies(corax_partition_t *  partition,
                                               corax_unode_t *      tree,
                                               unsigned int         params_index,
                                               const unsigned int * params_indices,
                                               double               bfgs_factor,
                                               double               tolerance)
{
  double              cur_logl;
  double *            x, *lb, *ub;
  int *               bt;
  unsigned int        i;
  struct freqs_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  opt_params.params_index   = params_index;

  double *     frequencies = partition->frequencies[params_index];
  unsigned int states      = partition->states;
  unsigned int cur_index;

  const double factor = bfgs_factor > 0. ? bfgs_factor : CORAX_ALGO_BFGS_FACTR;

  x  = (double *)malloc(sizeof(double) * (states - 1));
  lb = (double *)malloc(sizeof(double) * (states - 1));
  ub = (double *)malloc(sizeof(double) * (states - 1));
  bt = (int *)malloc(sizeof(int) * (states - 1));

  /* find first state with frequency > min_freq, and use it as fixed freq */
  opt_params.fixed_freq_state = states;
  for (i = 0; i < states; ++i)
  {
    if (frequencies[i] > CORAX_OPT_MIN_FREQ)
    {
      opt_params.fixed_freq_state = i;
      break;
    }
  }

  assert(opt_params.fixed_freq_state < states);
  assert(frequencies[opt_params.fixed_freq_state] > 0.);

  cur_index = 0;
  for (i = 0; i < states; ++i)
  {
    if (i != opt_params.fixed_freq_state)
    {
      x[cur_index]  = frequencies[i] / frequencies[opt_params.fixed_freq_state];
      lb[cur_index] = CORAX_OPT_MIN_FREQ;
      ub[cur_index] = CORAX_OPT_MAX_FREQ;
      bt[cur_index] = CORAX_OPT_LBFGSB_BOUND_BOTH;
      cur_index++;
    }
  }

  cur_logl = corax_opt_minimize_lbfgsb(x,
                                       lb,
                                       ub,
                                       bt,
                                       states - 1,
                                       factor,
                                       tolerance,
                                       (void *)&opt_params,
                                       &target_freqs_func);

  /* update frequencies */
  target_freqs_func((void *)&opt_params, x);

  free(x);
  free(lb);
  free(ub);
  free(bt);

  return cur_logl;
}

CORAX_EXPORT double corax_algo_opt_subst_rates(corax_partition_t *  partition,
                                               corax_unode_t *      tree,
                                               unsigned int         params_index,
                                               const unsigned int * params_indices,
                                               const int *          symmetries,
                                               double               min_rate,
                                               double               max_rate,
                                               double               bfgs_factor,
                                               double               tolerance)
{
  double       cur_logl;
  double *     x, *lb, *ub;
  int *        bt;
  unsigned int i, j, k;

  double *     subst_rates  = partition->subst_params[params_index];
  unsigned int states       = partition->states;
  unsigned int subst_params = (states * (states - 1)) / 2;
  unsigned int subst_free_params;

  const double factor = bfgs_factor > 0. ? bfgs_factor : CORAX_ALGO_BFGS_FACTR;

  if (!symmetries) { subst_free_params = subst_params - 1; }
  else
  {
    subst_free_params = 0;
    for (i = 0; i < subst_params; ++i)
    {
      if ((unsigned int)symmetries[i] > subst_free_params)
      {
        /* check that symmetries vector is correctly formatted */
        assert((unsigned int)symmetries[i] == (subst_free_params + 1));
        ++subst_free_params;
      }
    }
  }

  struct algo_subst_params opt_params;
  opt_params.partition         = partition;
  opt_params.tree              = tree;
  opt_params.params_index      = params_index;
  opt_params.params_indices    = params_indices;
  opt_params.symmetries        = symmetries;
  opt_params.subst_free_params = subst_free_params;

  x  = (double *)malloc(sizeof(double) * (subst_free_params));
  lb = (double *)malloc(sizeof(double) * (subst_free_params));
  ub = (double *)malloc(sizeof(double) * (subst_free_params));
  bt = (int *)malloc(sizeof(int) * (subst_free_params));

  k = 0;
  for (i = 0; i < subst_free_params; ++i)
  {
    bt[i] = CORAX_OPT_LBFGSB_BOUND_BOTH;
    lb[i] = min_rate;
    ub[i] = max_rate;

    if (symmetries)
    {
      if ((unsigned int)symmetries[subst_params - 1] == k) ++k;

      for (j = 0; j < subst_params; ++j)
      {
        if ((unsigned int)symmetries[j] == k)
        {
          x[i] = subst_rates[j];
          break;
        }
      }
      ++k;
    }
    else
    {
      x[i] = subst_rates[i];
    }

    if (!x[i])
    {
      /* initialize to interval center */
      x[i] = (min_rate + max_rate) / 2.0;
    }
    if (x[i] < min_rate)
    {
      /* set to lower bound */
      x[i] = min_rate;
    }
    else if (x[i] > max_rate)
    {
      /* set to upper bound */
      x[i] = max_rate;
    }
  }

  cur_logl = corax_opt_minimize_lbfgsb(x,
                                       lb,
                                       ub,
                                       bt,
                                       subst_free_params,
                                       factor,
                                       tolerance,
                                       (void *)&opt_params,
                                       target_subst_params_func);

  free(x);
  free(lb);
  free(ub);
  free(bt);

  return cur_logl;
}

CORAX_EXPORT double corax_algo_opt_alpha(corax_partition_t  * partition,
                                         corax_unode_t      * tree,
                                         const unsigned int * params_indices,
                                         double               min_alpha,
                                         double               max_alpha,
                                         double *             alpha,
                                         double               tolerance)
{
  double cur_logl;
  double f2x;
  double xres;

  struct default_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  opt_params.gamma_mode     = CORAX_GAMMA_RATES_MEAN; // for now

  xres = corax_opt_minimize_brent(min_alpha,
                                  *alpha,
                                  max_alpha,
                                  tolerance,
                                  &cur_logl,
                                  &f2x,
                                  (void *)&opt_params,
                                  &target_alpha_func);

  cur_logl = target_alpha_func(&opt_params, xres);
  *alpha   = xres;

  return cur_logl;
}

CORAX_EXPORT double corax_algo_opt_pinv(corax_partition_t *  partition,
                                        corax_unode_t *      tree,
                                        const unsigned int * params_indices,
                                        double               min_pinv,
                                        double               max_pinv,
                                        double               tolerance)
{
  double                cur_logl;
  double                f2x;
  double                xres;
  double                start_pinv;
  struct default_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  start_pinv                = partition->prop_invar[params_indices[0]];

  xres = corax_opt_minimize_brent(min_pinv,
                                  start_pinv,
                                  max_pinv,
                                  tolerance,
                                  &cur_logl,
                                  &f2x,
                                  (void *)&opt_params,
                                  &target_pinv_func);

  cur_logl = target_pinv_func(&opt_params, xres);

  return cur_logl;
}

CORAX_EXPORT double corax_algo_opt_alpha_pinv(corax_partition_t *  partition,
                                              corax_unode_t *      tree,
                                              const unsigned int * params_indices,
                                              double               min_alpha,
                                              double               max_alpha,
                                              double *             alpha,
                                              double               min_pinv,
                                              double               max_pinv,
                                              double               bfgs_factor,
                                              double               tolerance)
{
  double cur_logl;
  double x[2], lb[2], ub[2];
  int    bt[2];

  const double factor = bfgs_factor > 0. ? bfgs_factor : CORAX_ALGO_BFGS_FACTR;

  struct default_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  opt_params.gamma_mode     = CORAX_GAMMA_RATES_MEAN; // for now

  /* init alpha */
  x[0]  = *alpha;
  lb[0] = min_alpha > 0. ? min_alpha : CORAX_OPT_MIN_ALPHA;
  ub[0] = max_alpha > 0. ? max_alpha : CORAX_OPT_MAX_ALPHA;
  bt[0] = CORAX_OPT_LBFGSB_BOUND_BOTH;

  /* init p-inv */
  x[1]  = partition->prop_invar[params_indices[0]];
  lb[1] = min_pinv > CORAX_ALGO_LBFGSB_ERROR
              ? min_pinv
              : CORAX_OPT_MIN_PINV + CORAX_ALGO_LBFGSB_ERROR;
  ub[1] = max_pinv > 0. ? max_pinv : CORAX_OPT_MAX_PINV;
  bt[1] = CORAX_OPT_LBFGSB_BOUND_BOTH;

  cur_logl = corax_opt_minimize_lbfgsb(x,
                                       lb,
                                       ub,
                                       bt,
                                       2,
                                       factor,
                                       tolerance,
                                       (void *)&opt_params,
                                       &target_alpha_pinv_func);

  /* save optimal alpha (p-inv is stored in the partition) */
  *alpha = x[0];

  return cur_logl;
}

CORAX_EXPORT double corax_algo_opt_rates_weights(corax_partition_t *  partition,
                                                 corax_unode_t *      tree,
                                                 const unsigned int * params_indices,
                                                 double               min_rate,
                                                 double               max_rate,
                                                 double               bfgs_factor,
                                                 double               tolerance,
                                                 double *             brlen_scaler,
                                                 int                  scale_branches)
{
  double       cur_logl, prev_logl;
  double       sum_weightrates, rate_scaler;
  double *     x, *lb, *ub;
  int *        bt;
  unsigned int i;

  double *     rates     = partition->rates;
  double *     weights   = partition->rate_weights;
  unsigned int rate_cats = partition->rate_cats;

  struct rate_weights_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;

  const double factor = bfgs_factor > 0. ? bfgs_factor : CORAX_ALGO_BFGS_FACTR;

  x  = (double *)malloc(sizeof(double) * (rate_cats));
  lb = (double *)malloc(sizeof(double) * (rate_cats));
  ub = (double *)malloc(sizeof(double) * (rate_cats));
  bt = (int *)malloc(sizeof(int) * (rate_cats));

  /* 2 step BFGS */

  cur_logl = 0;
  do {
    prev_logl = cur_logl;

    /* optimize mixture weights */

    fill_weights(
        weights, &(opt_params.fixed_weight_state), x, bt, lb, ub, rate_cats);

    cur_logl = 1
               * corax_opt_minimize_lbfgsb(x,
                                           lb,
                                           ub,
                                           bt,
                                           rate_cats - 1,
                                           factor,
                                           tolerance,
                                           (void *)&opt_params,
                                           target_weights_func);

    /* optimize mixture rates */

    fill_rates(rates, x, bt, lb, ub, min_rate, max_rate, rate_cats);

    cur_logl = corax_opt_minimize_lbfgsb(x,
                                         lb,
                                         ub,
                                         bt,
                                         rate_cats,
                                         factor,
                                         tolerance,
                                         (void *)&opt_params,
                                         target_rates_func);

  } while (!prev_logl || prev_logl - cur_logl > tolerance);

  /* force constraint sum(weights x rates) = 1.0 */
  sum_weightrates = 0.0;
  for (i = 0; i < rate_cats; ++i) sum_weightrates += rates[i] * weights[i];
  rate_scaler = 1.0 / sum_weightrates;

  for (i = 0; i < rate_cats; ++i) rates[i] *= rate_scaler;

  *brlen_scaler = sum_weightrates;

  if (scale_branches)
  {
    /* scale branch lengths such that likelihood is conserved */
    corax_utree_scale_branches_all(tree, sum_weightrates);

    /* update pmatrices and partials according to the new branches */
    cur_logl = -1
               * corax_opt_compute_lk(partition,
                                      tree,
                                      params_indices,
                                      1,  /* update pmatrices */
                                      1); /* update partials */
  }

  free(x);
  free(lb);
  free(ub);
  free(bt);

  return cur_logl;
}

CORAX_EXPORT double corax_algo_opt_brlen_scaler(corax_partition_t *  partition,
                                                corax_unode_t *      root,
                                                const unsigned int * params_indices,
                                                double *             scaler,
                                                double               min_scaler,
                                                double               max_scaler,
                                                double               tolerance)
{
  double                     cur_logl;
  double                     f2x;
  double                     xres;
  struct brlen_scaler_params opt_params;

  /* create a temporary tree with the scaled branches */
  corax_unode_t *scaled_tree = corax_utree_graph_clone(root);
  corax_utree_scale_branches_all(scaled_tree, *scaler);

  opt_params.partition      = partition;
  opt_params.tree           = scaled_tree;
  opt_params.params_indices = params_indices;
  opt_params.old_scaler     = *scaler;

  xres = corax_opt_minimize_brent(min_scaler,
                                  *scaler,
                                  max_scaler,
                                  tolerance,
                                  &cur_logl,
                                  &f2x,
                                  (void *)&opt_params,
                                  &target_brlen_scaler_func);

  cur_logl = target_brlen_scaler_func(&opt_params, xres);

  corax_utree_graph_destroy(scaled_tree, NULL);

  *scaler = xres;

  return cur_logl;
}

/* STATIC FUNCTIONS */

static void fill_rates(double *     rates,
                       double *     x,
                       int *        bt,
                       double *     lb,
                       double *     ub,
                       double       min_rate,
                       double       max_rate,
                       unsigned int n_rates)
{
  unsigned int i;

  assert(min_rate > 1e-4 && max_rate > min_rate);

  for (i = 0; i < n_rates; ++i)
  {
    bt[i] = CORAX_OPT_LBFGSB_BOUND_BOTH;
    lb[i] = min_rate;
    ub[i] = max_rate;

    if (rates[i] < min_rate)
      x[i] = min_rate;
    else if (rates[i] > max_rate)
      x[i] = max_rate;
    else
      x[i] = rates[i];
  }
}

static void fill_weights(double *      weights,
                         unsigned int *fixed_weight_index,
                         double *      x,
                         int *         bt,
                         double *      lb,
                         double *      ub,
                         unsigned int  n_weights)
{
  unsigned int i, cur_index = 0;

  *fixed_weight_index = n_weights;
  for (i = 0; i < n_weights; ++i)
  {
    if (weights[i] > CORAX_OPT_MIN_FREQ)
    {
      *fixed_weight_index = i;
      break;
    }
  }

  assert(*fixed_weight_index < n_weights);
  assert(weights[*fixed_weight_index] > 0.);

  for (i = 0; i < n_weights; ++i)
  {
    if (i != *fixed_weight_index)
    {
      bt[cur_index] = CORAX_OPT_LBFGSB_BOUND_BOTH;

      double r      = weights[i] / weights[*fixed_weight_index];
      lb[cur_index] = CORAX_ALGO_MIN_WEIGHT_RATIO;
      ub[cur_index] = CORAX_ALGO_MAX_WEIGHT_RATIO;
      if (r < lb[cur_index])
        x[cur_index] = lb[cur_index];
      else if (r > ub[cur_index])
        x[cur_index] = ub[cur_index];
      else
        x[cur_index] = r;

      cur_index++;
    }
  }
}

/*
    Copyright (C) 2016 Diego Darriba

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
#ifndef CORAX_OPTIMIZE_MODEL_H_
#define CORAX_OPTIMIZE_MODEL_H_

#include "opt_generic.h"

#define CORAX_ALGO_MIN_WEIGHT_RATIO 0.001
#define CORAX_ALGO_MAX_WEIGHT_RATIO 10
#define CORAX_ALGO_BFGS_FACTR 1e9

/*
 * Optimize stationary frequencies for parameters `params_index`.
 */
CORAX_EXPORT double corax_algo_opt_frequencies(corax_partition_t *  partition,
                                               corax_unode_t *      tree,
                                               unsigned int         params_index,
                                               const unsigned int * params_indices,
                                               double               bfgs_factor,
                                               double               tolerance);

/*
 * Optimize substitution rates for parameters `params_index`.
 * symmetries is an array with as many positions as substitution parameters,
 *            or 'NULL' for GTR (e.g., for DNA, 012345 and NULL are equivalent)
 *            Must be sorted and start with '0'.
 *            e.g., 000000 = JC/F81, 010010 = K80/HKY, 012314 = TrN
 */
CORAX_EXPORT double corax_algo_opt_subst_rates(corax_partition_t *  partition,
                                               corax_unode_t *      tree,
                                               unsigned int         params_index,
                                               const unsigned int * params_indices,
                                               const int *          symmetries,
                                               double               min_rate,
                                               double               max_rate,
                                               double               bfgs_factor,
                                               double               tolerance);

CORAX_EXPORT double corax_algo_opt_alpha(corax_partition_t  * partition,
                                         corax_unode_t      * tree,
                                         const unsigned int * params_indices,
                                         double               min_alpha,
                                         double               max_alpha,
                                         double *             alpha,
                                         double               tolerance);

CORAX_EXPORT double corax_algo_opt_pinv(corax_partition_t *  partition,
                                        corax_unode_t *      tree,
                                        const unsigned int * params_indices,
                                        double               min_pinv,
                                        double               max_pinv,
                                        double               tolerance);

CORAX_EXPORT double corax_algo_opt_alpha_pinv(corax_partition_t *  partition,
                                              corax_unode_t *      tree,
                                              const unsigned int * params_indices,
                                              double               min_alpha,
                                              double               max_alpha,
                                              double *             alpha,
                                              double               min_pinv,
                                              double               max_pinv,
                                              double               bfgs_factor,
                                              double               tolerance);

/*
 * Optimize free rates and rate weights together, linked to
 * `partition->rate_cats`. Uses 2 step L-BFGS-B algorithm.
 */
CORAX_EXPORT double corax_algo_opt_rates_weights(corax_partition_t *  partition,
                                                 corax_unode_t *      tree,
                                                 const unsigned int * params_indices,
                                                 double               min_rate,
                                                 double               max_rate,
                                                 double               bfgs_factor,
                                                 double               tolerance,
                                                 double *             brlen_scaler,
                                                 int                  scale_branches);

CORAX_EXPORT double corax_algo_opt_brlen_scaler(corax_partition_t *  partition,
                                                corax_unode_t *      tree,
                                                const unsigned int * params_indices,
                                                double *             scaler,
                                                double               min_scaler,
                                                double               max_scaler,
                                                double               tolerance);

#endif /* CORAX_OPTIMIZE_MODEL_H_ */

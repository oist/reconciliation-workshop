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

#ifndef CORAX_TOPOLOGICAL_TESTS_H_
#define CORAX_TOPOLOGICAL_TESTS_H_

#ifdef __cplusplus
extern "C"
{
#endif
 
 #include "corax/corax_tree.h"
 #include "corax/corax_core.h"

 #define CORAX_KH_ORDER_ERROR        6006
 
 /**
   * RELL Bootstrap version of the Kishino-Hasegawa test (KH test)
   * Statistical test between two (different) trees, each one having its own parameters, but assuming the same MSA.
   *  
   * @param  persite_lnl_1          The per-site likelihoods calculated for the tree with the highest likelihood score. 
   *                                The array can be calculated using corax_treeinfo_compute_loglh_persite() function.
   * @param  persite_lnl_2          The per-site likelihoods calculated for the tree with the lowest likelihood score. 
   *                                The array can be calculated using corax_treeinfo_compute_loglh_persite() function.
   * @param  partition_lens         An array where the lengths (number of sites) for each partition are stored.
   * @param  partition_count        Number of partitions of the MSA
   * @param  nBootstrap             Number of bootstrap replicates (e.g. 1000)
   * @param  p_value_threshold      A relatively small threshold (e.g. 0.05)
   * @return                        CORAX_SUCCESS or CORAX_FAILURE, depending on whether the second tree passed or failed the test, respectively.
   *                                Success indicates that the two trees are not significantly different. 
    */
 CORAX_EXPORT int corax_KH_RELL_test(double **persite_lnl_1, 
                                    double **persite_lnl_2,
                                    int *partition_lens, 
                                    int partition_count, 
                                    int nBootstrap,
                                    double p_value_threshold);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_TOPOLOGICAL_TESTS_H_ */
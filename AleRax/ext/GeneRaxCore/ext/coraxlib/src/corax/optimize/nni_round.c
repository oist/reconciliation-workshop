/*
 Copyright (C) 2015-21 Diego Darriba, Alexey Kozlov

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

#include "corax/corax.h"
#include "opt_treeinfo.h"
#include "opt_branches.h"
#include <time.h>

#define DEBUG_MODE 0

// move structure
typedef struct corax_nni_move_s
{
  int           type;
  double      **initial_branch_lengths;
  unsigned int  initial_pmatrix_indices[5];
  unsigned int  partition_num;
  int           collect_branches;

} corax_nni_move_t;

static corax_nni_move_t *create_move(corax_treeinfo_t *treeinfo)
{
  corax_nni_move_t *move = (corax_nni_move_t *) malloc(sizeof(corax_nni_move_t));

  move->initial_branch_lengths =
      (double **) malloc((treeinfo->partition_count) * sizeof(double *));

  move->collect_branches = 0;
  if (treeinfo->partition_count == 1
      || treeinfo->brlen_linkage == CORAX_BRLEN_LINKED
      || treeinfo->brlen_linkage == CORAX_BRLEN_SCALED)
  {
    move->collect_branches = 1;
  }

  move->partition_num = treeinfo->partition_count;

  for (unsigned int part = 0; part < treeinfo->partition_count; part++)
  {
    move->initial_branch_lengths[part] = (double *) malloc(5 * sizeof(double));
  }

  return move;
}

static int update_move(corax_treeinfo_t *treeinfo,
                       corax_nni_move_t *move,
                       corax_unode_t    *node,
                       int               type)
{
  move->initial_pmatrix_indices[0] = node->next->pmatrix_index;
  move->initial_pmatrix_indices[1] = node->next->next->pmatrix_index;
  move->initial_pmatrix_indices[2] = node->pmatrix_index;
  move->initial_pmatrix_indices[3] = node->back->next->pmatrix_index;
  move->initial_pmatrix_indices[4] = node->back->next->next->pmatrix_index;

  for (unsigned int part = 0; part < treeinfo->partition_count; part++)
  {
    if (move->collect_branches)
    {
      move->initial_branch_lengths[part][0] = node->next->length;
      move->initial_branch_lengths[part][1] = node->next->next->length;
      ;
      move->initial_branch_lengths[part][2] = node->length;
      move->initial_branch_lengths[part][3] = node->back->next->length;
      move->initial_branch_lengths[part][4] = node->back->next->next->length;
    }
    else
    {
      move->initial_branch_lengths[part][0] =
          treeinfo->branch_lengths[part][move->initial_pmatrix_indices[0]];
      move->initial_branch_lengths[part][1] =
          treeinfo->branch_lengths[part][move->initial_pmatrix_indices[1]];
      move->initial_branch_lengths[part][2] =
          treeinfo->branch_lengths[part][move->initial_pmatrix_indices[2]];
      move->initial_branch_lengths[part][3] =
          treeinfo->branch_lengths[part][move->initial_pmatrix_indices[3]];
      move->initial_branch_lengths[part][4] =
          treeinfo->branch_lengths[part][move->initial_pmatrix_indices[4]];
    }
  }

  move->type = type;

  return CORAX_SUCCESS;
}

static void delete_move(corax_nni_move_t *move)
{
  for (unsigned int part = 0; part < move->partition_num; part++)
  {
    free(move->initial_branch_lengths[part]);
  }

  free(move->initial_branch_lengths);
//  free(move->initial_pmatrix_indices);
  free(move);
}

// This function is stolen from spr round
int algo_update_pmatrix_nni(corax_treeinfo_t *treeinfo, corax_unode_t *edge)
{
  unsigned int p;
  unsigned int updated       = 0;
  unsigned int pmatrix_index = edge->pmatrix_index;

  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* only selected partitioned will be affected */
    if (treeinfo->partitions[p])
    {
      if (treeinfo->pmatrix_valid[p][pmatrix_index]) continue;

      double p_brlen = treeinfo->branch_lengths[p][pmatrix_index];
      if (treeinfo->brlen_linkage == CORAX_BRLEN_SCALED)
        p_brlen *= treeinfo->brlen_scalers[p];

      int ret = corax_update_prob_matrices(treeinfo->partitions[p],
                                           treeinfo->param_indices[p],
                                           &pmatrix_index,
                                           &p_brlen,
                                           1);

      if (!ret) return CORAX_FAILURE;

      treeinfo->pmatrix_valid[p][pmatrix_index] = 1;
      updated++;
    }
  }

  return CORAX_SUCCESS;
}

static double optimize_branch_lengths_quartet(corax_treeinfo_t *treeinfo,
                                              double            lh_epsilon,
                                              int    brlen_opt_method,
                                              double bl_min,
                                              double bl_max,
                                              int    smoothings)
{
  double logl;

  logl = corax_treeinfo_compute_loglh_flex(treeinfo, 1, 0);

  // new optimization function - we need to hide that from the user I guess (or
  // maybe create a wrapper function)
  logl = -1
         * corax_opt_optimize_branch_lengths_local_multi_quartet(
             treeinfo->partitions,
             treeinfo->partition_count,
             treeinfo->root,
             treeinfo->param_indices,
             treeinfo->deriv_precomp,
             treeinfo->branch_lengths,
             treeinfo->brlen_scalers,
             bl_min,
             bl_max,
             lh_epsilon,
             smoothings,
             1, // keep_update
             brlen_opt_method,
             treeinfo->brlen_linkage,
             treeinfo->parallel_context,
             treeinfo->parallel_reduce_cb);

  if (logl == CORAX_FAILURE)
  {
    printf("Branch Length Optimization Error\n");
    printf("Corax Error num = %d\n", corax_errno);
    printf("Corax Error msg = %s\n", corax_errmsg);
    return CORAX_FAILURE;
  }

  return logl;
}

static double apply_move(corax_treeinfo_t *treeinfo,
                         corax_unode_t    *node,
                         corax_nni_move_t *move,
                         int               brlen_opt_method,
                         double            bl_min,
                         double            bl_max,
                         int               smoothings,
                         double            lh_epsilon)
{
  if (move->type == CORAX_UTREE_MOVE_NNI_LEFT ||
      move->type == CORAX_UTREE_MOVE_NNI_RIGHT)
  {
    if (CORAX_UTREE_IS_TIP(node) || CORAX_UTREE_IS_TIP(node->back))
    {
      /* invalid move */
      corax_set_error(CORAX_TREE_ERROR_NNI_LEAF,
                      "Attempting to apply NNI on a leaf branch\n");
      return CORAX_FAILURE;
    }

    corax_utree_nni(node, move->type, NULL);
    corax_treeinfo_invalidate_clv(treeinfo, node);
    corax_treeinfo_invalidate_clv(treeinfo, node->back);

    double loglh = optimize_branch_lengths_quartet(
        treeinfo, lh_epsilon, brlen_opt_method, bl_min, bl_max, smoothings);

    if (loglh == CORAX_FAILURE)
    {
      printf("Something went wrong with the branch-length optimization. Exit..\n");
      return CORAX_FAILURE;
    }
    return loglh;
  }
  else
  {
    /* invalid move */
    corax_set_error(CORAX_TREE_ERROR_NNI_INVALID_MOVE,
                    "Invalid NNI move type\n");
    return CORAX_FAILURE;
  }
}

static int undo_move(corax_treeinfo_t *treeinfo,
                     corax_unode_t    *node,
                     corax_nni_move_t *move)
{
  int retval = CORAX_SUCCESS;

  retval = corax_utree_nni(node, move->type, NULL);

  node->next->pmatrix_index = node->next->back->pmatrix_index =
      move->initial_pmatrix_indices[0];
  node->next->next->pmatrix_index = node->next->next->back->pmatrix_index =
      move->initial_pmatrix_indices[1];
  node->pmatrix_index = node->back->pmatrix_index =
      move->initial_pmatrix_indices[2];
  node->back->next->pmatrix_index = node->back->next->back->pmatrix_index =
      move->initial_pmatrix_indices[3];
  node->back->next->next->pmatrix_index =
      node->back->next->next->back->pmatrix_index =
          move->initial_pmatrix_indices[4];

  for (unsigned int part = 0; part < move->partition_num; part++)
  {
    treeinfo->branch_lengths[part][node->next->pmatrix_index] =
        move->initial_branch_lengths[part][0];
    treeinfo->branch_lengths[part][node->next->next->pmatrix_index] =
        move->initial_branch_lengths[part][1];
    treeinfo->branch_lengths[part][node->pmatrix_index] =
        move->initial_branch_lengths[part][2];
    treeinfo->branch_lengths[part][node->back->next->pmatrix_index] =
        move->initial_branch_lengths[part][3];
    treeinfo->branch_lengths[part][node->back->next->next->pmatrix_index] =
        move->initial_branch_lengths[part][4];
  }

  if (move->collect_branches)
  {

    node->next->length = node->next->back->length =
        move->initial_branch_lengths[0][0];
    node->next->next->length = node->next->next->back->length =
        move->initial_branch_lengths[0][1];
    node->length = node->back->length = move->initial_branch_lengths[0][2];
    node->back->next->length          = node->back->next->back->length =
        move->initial_branch_lengths[0][3];
    node->back->next->next->length = node->back->next->next->back->length =
        move->initial_branch_lengths[0][4];
  }

  corax_treeinfo_invalidate_pmatrix(treeinfo, node);
  corax_treeinfo_invalidate_pmatrix(treeinfo, node->next);
  corax_treeinfo_invalidate_pmatrix(treeinfo, node->next->next);
  corax_treeinfo_invalidate_pmatrix(treeinfo, node->back->next);
  corax_treeinfo_invalidate_pmatrix(treeinfo, node->back->next->next);

  retval &= algo_update_pmatrix_nni(treeinfo, node);
  retval &= algo_update_pmatrix_nni(treeinfo, node->next);
  retval &= algo_update_pmatrix_nni(treeinfo, node->next->next);
  retval &= algo_update_pmatrix_nni(treeinfo, node->back->next);
  retval &= algo_update_pmatrix_nni(treeinfo, node->back->next->next);

  corax_treeinfo_invalidate_clv(treeinfo, node);
  corax_treeinfo_invalidate_clv(treeinfo, node->back);

  return retval;
}

static int shSupport(corax_treeinfo_t *treeinfo,
                     double           *shSupportValues,
                     double          **persite_lnl_0,
                     int               nBootstrap,
                     double            shEpsilon,
                     double            tolerance,
                     int               brlen_opt_method,
                     double            bl_min,
                     double            bl_max,
                     int               smoothings,
                     double            lh_epsilon,
                     bool             *warning_printed)
{
  double LNL0, LNL1, LNL2, _LNL0, _LNL1, _LNL2, second_best_logl;
  int    nSupport          = 0;
  bool   shittySplit       = false;
  bool   non_optimal_split = false;
  double aLRT              = 0;

  double **persite_lnl_1, **persite_lnl_2;

  persite_lnl_1 =
      (double **)malloc(sizeof(double *) * (treeinfo->partition_count));
  persite_lnl_2 =
      (double **)malloc(sizeof(double *) * (treeinfo->partition_count));

  for (unsigned int i = 0; i < treeinfo->partition_count; i++)
  {
    persite_lnl_1[i] =
        (double *)malloc(sizeof(double) * treeinfo->partitions[i]->sites);
    persite_lnl_2[i] =
        (double *)malloc(sizeof(double) * treeinfo->partitions[i]->sites);
  }

  corax_unode_t    *q    = treeinfo->root;
  corax_nni_move_t *move = create_move(treeinfo);

  // update clvs and matrices
  LNL0 = corax_treeinfo_compute_loglh_flex(treeinfo, 1, 0);
  if (DEBUG_MODE) printf("Loglh = %f, ", LNL0);

  // First NNI topology
  update_move(treeinfo, move, q, CORAX_UTREE_MOVE_NNI_LEFT);
  LNL1 = apply_move(treeinfo,
                    q,
                    move,
                    brlen_opt_method,
                    bl_min,
                    bl_max,
                    smoothings,
                    lh_epsilon);

  if (LNL1 == CORAX_FAILURE)
  {
    printf("Failed to apply NNI move. \n");
    printf("Corax error number %d \n", corax_errno);
    printf("Corax error msg: %s\n", corax_errmsg);
    delete_move(move);
    return CORAX_FAILURE;
  }

  // calculate persite
  assert(
      fabs(corax_treeinfo_compute_loglh_persite(treeinfo, 1, 0, persite_lnl_1)
           - LNL1)
      < 1e-5);

  second_best_logl = LNL1;

  // undo + testing
  if (!undo_move(treeinfo, q, move))
  {
    corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                    "Error in undoing move that generates a tree of worst "
                    "likelihood ...\n");
    return CORAX_FAILURE;
  }

  // second NNI topology
  update_move(treeinfo, move, q, CORAX_UTREE_MOVE_NNI_RIGHT);
  LNL2 = apply_move(treeinfo,
                    q,
                    move,
                    brlen_opt_method,
                    bl_min,
                    bl_max,
                    smoothings,
                    lh_epsilon);

  if (LNL2 > second_best_logl) second_best_logl = LNL2;

  // calculate persite
  assert(
      fabs(corax_treeinfo_compute_loglh_persite(treeinfo, 1, 0, persite_lnl_2)
           - LNL2)
      < 1e-5);

  // undo + testing
  if (!undo_move(treeinfo, q, move))
  {
    corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                    "Error in undoing move that generates a tree of worst "
                    "likelihood ...\n");
    return CORAX_FAILURE;
  }

  // assert(fabs(compute_local_likelihood_treeinfo(treeinfo, q, 1) - LNL0) <
  // 1e-5);
  delete_move(move);

  if (second_best_logl > LNL0)
  {

    if (second_best_logl - LNL0 > tolerance)
    {

      if (!(*warning_printed))
      {
        printf("\nWARNING!\nYou are trying to calculate the SH-like aLRT "
               "support values for a non NNI-optimal tree topology."
               " For the internal branches that aren't NNI optimal, the "
               "SH-like aLRT metric will be equal to -inf.\n"
               "RECOMMENDED: Since the SH-like aLRT statistics makes more "
               "sense for NNI-optimal topologies, we recommend the user to "
               "first run the "
               "corax_algo_nni_round() function, with a relatively small "
               "tolerance value (e.g. 0.1 or 0.01) to optimize the tree "
               "topology.\n"
               "In case you have already called the NNI-optimization function, "
               "try again with a smaller tolerance value.\n\n");
        *warning_printed = true;
      }
      non_optimal_split = true;
    }
    else { shittySplit = true; }
  }
  else { aLRT = 2 * (LNL0 - second_best_logl); }

  if (shittySplit)
  {
    shSupportValues[q->pmatrix_index] = 0;
    if (DEBUG_MODE)
      printf("Branch %d, shitty split, Support = 0\n", q->pmatrix_index);
    return CORAX_SUCCESS;
  }

  if (non_optimal_split)
  {
    shSupportValues[q->pmatrix_index] = -INFINITY;
    if (DEBUG_MODE) printf("Branch %d, Support = -inf\n", q->pmatrix_index);
    return CORAX_SUCCESS;
  }

  unsigned int nSites = 0;
  for (unsigned int i = 0; i < treeinfo->partition_count; i++)
  {
    nSites += treeinfo->partitions[i]->sites;
  }

  unsigned int        _seed  = (unsigned int)rand();
  corax_random_state *rstate = corax_random_create(_seed);

  for (int i = 0; i < nBootstrap; i++)
  {

    double cs[3];
    double tmp, diff;
    int    pIndex, sIndex;
    _LNL0 = 0;
    _LNL1 = 0;
    _LNL2 = 0;

    for (unsigned int j = 0; j < nSites; j++)
    {
      if (treeinfo->partition_count > 1)
      {
        pIndex = corax_random_getint(rstate, (int)treeinfo->partition_count);
      }
      else { pIndex = 0; }

      sIndex =
          corax_random_getint(rstate, (int)treeinfo->partitions[pIndex]->sites);

      _LNL0 += persite_lnl_0[pIndex][sIndex];
      _LNL1 += persite_lnl_1[pIndex][sIndex];
      _LNL2 += persite_lnl_2[pIndex][sIndex];
    }

    cs[0] = _LNL0 - LNL0;
    cs[1] = _LNL1 - LNL1;
    cs[2] = _LNL2 - LNL2;

    if (cs[0] < cs[1])
    {
      tmp   = cs[0];
      cs[0] = cs[1];
      cs[1] = tmp;
    }

    if (cs[1] < cs[2])
    {
      tmp   = cs[1];
      cs[1] = cs[2];
      cs[2] = tmp;
    }

    diff = fabs(cs[0] - cs[1]);
    if (aLRT > 2 * diff + shEpsilon) { nSupport++; }
  }

  // sh support value
  double support                    = (nSupport / (double)nBootstrap);
  shSupportValues[q->pmatrix_index] = support * 100.0;
  if (DEBUG_MODE)
    printf("Branch %d, n_support = %d, Support = %f\n",
           q->pmatrix_index,
           nSupport,
           shSupportValues[q->pmatrix_index]);

  // free heap
  corax_random_destroy(rstate);

  for (unsigned i = 0; i < treeinfo->partition_count; i++)
  {
    free(persite_lnl_1[i]);
    free(persite_lnl_2[i]);
  }

  free(persite_lnl_1);
  free(persite_lnl_2);

  return CORAX_SUCCESS;
}

/* do a single NNI move around the root branch */
CORAX_EXPORT double corax_algo_nni_single(corax_treeinfo_t *treeinfo,
                                          corax_unode_t    *edge,
                                          corax_nni_move_t *move,
                                          int               nni_type,
                                          double            old_logl,
                                          int               brlen_opt_method,
                                          double            bl_min,
                                          double            bl_max,
                                          int               smoothings,
                                          double            lh_epsilon)
{
  double new_logl;

  // setup move
  update_move(treeinfo, move, edge, nni_type);

  // apply nni move
  new_logl = apply_move(treeinfo,
                        edge,
                        move,
                        brlen_opt_method,
                        bl_min,
                        bl_max,
                        smoothings,
                        lh_epsilon);

  if (new_logl == CORAX_FAILURE)
  {
    printf("Failed to apply NNI move. \n");
    printf("Corax error number %d \n", corax_errno);
    printf("Corax error msg: %s\n", corax_errmsg);
    return CORAX_FAILURE;
  }
  else if (new_logl <= old_logl)
  {
    if (!undo_move(treeinfo, edge, move))
    {
      corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                      "Error in undoing move that generates a tree of worse "
                      "likelihood ...\n");
      return CORAX_FAILURE;
    }
  }

  return new_logl;
}

/* do a pair of NNI moves around the root branch */
CORAX_EXPORT double corax_algo_nni_local(corax_treeinfo_t *treeinfo,
                                         int               brlen_opt_method,
                                         double            bl_min,
                                         double            bl_max,
                                         int               smoothings,
                                         double            lh_epsilon)
{
  double tree_logl, new_logl;
  int nni_type;
  bool check_cons = true;

  corax_unode_t    *q    = treeinfo->root;
  corax_nni_move_t *move = create_move(treeinfo);

  // setup move
  tree_logl = optimize_branch_lengths_quartet(treeinfo, lh_epsilon, brlen_opt_method,
                                              bl_min, bl_max, smoothings);
  if (tree_logl == CORAX_FAILURE)
  {
    printf(
        "Something went wrong with the branch-length optimization. Exit..\n");
    return CORAX_FAILURE;
  }

  /* LEFT nni move */
  nni_type = CORAX_UTREE_MOVE_NNI_LEFT;
  /* check if resulting tree would contradict the constraint */
  if (!check_cons || corax_treeinfo_constraint_check_nni(treeinfo, q, nni_type))
  {
    // test nni move
    new_logl = corax_algo_nni_single(treeinfo,
                                     q,
                                     move,
                                     nni_type,
                                     tree_logl,
                                     brlen_opt_method,
                                     bl_min,
                                     bl_max,
                                     smoothings,
                                     lh_epsilon);

    if (new_logl > tree_logl)
    {
      corax_treeinfo_constraint_update_splits(treeinfo);
      tree_logl = new_logl;
    }
    else if (new_logl == CORAX_FAILURE)
    {
      delete_move(move);
      return CORAX_FAILURE;
    }
  }
  else
  {
    DBG("SKIP incompatible LEFT NNI: %u\n", q->clv_index);
  }

  /* RIGHT nni move */
  nni_type = CORAX_UTREE_MOVE_NNI_RIGHT;
  /* check if resulting tree would contradict the constraint */
  if (!check_cons || corax_treeinfo_constraint_check_nni(treeinfo, q, nni_type))
  {
    // test nni move
    new_logl = corax_algo_nni_single(treeinfo,
                                     q,
                                     move,
                                     nni_type,
                                     tree_logl,
                                     brlen_opt_method,
                                     bl_min,
                                     bl_max,
                                     smoothings,
                                     lh_epsilon);

    if (new_logl > tree_logl)
    {
      corax_treeinfo_constraint_update_splits(treeinfo);
      tree_logl = new_logl;
    }
    else if (new_logl == CORAX_FAILURE)
    {
      delete_move(move);
      return CORAX_FAILURE;
    }
  }
  else
  {
    DBG("SKIP incompatible RIGHT NNI: %u\n", q->clv_index);
  }

  delete_move(move);

  return tree_logl;
}

static int nni_recursive(corax_treeinfo_t *treeinfo,
                         corax_unode_t    *node,
                         unsigned int     *interchagnes,
                         bool              compute_aLRT,
                         double           *shSupportValues,
                         double          **persite_lnl,
                         int               nBootstrap,
                         double            shEpsilon,
                         double            tolerance,
                         int               brlen_opt_method,
                         double            bl_min,
                         double            bl_max,
                         int               smoothings,
                         double            lh_epsilon,
                         bool             *warning_printed)
{

  int            retval = CORAX_SUCCESS;
  corax_unode_t *q      = node->back;
  corax_unode_t *pb1    = node->next->back, *pb2 = node->next->next->back;

  double tree_logl;
  double new_logl;

  if (!CORAX_UTREE_IS_TIP(q))
  {
    corax_treeinfo_set_root(treeinfo, q);
    corax_treeinfo_invalidate_clv(treeinfo, q);
    tree_logl = corax_treeinfo_compute_loglh_flex(treeinfo, 1, 0);

    if (!compute_aLRT)
    {
      corax_unode_t *test_node = q->next->back;
      new_logl                 = corax_algo_nni_local(treeinfo, brlen_opt_method, bl_min, bl_max,
                                                      smoothings, lh_epsilon);

      if (new_logl == CORAX_FAILURE) return CORAX_FAILURE;

      if (DEBUG_MODE)
        printf("Old logl = %f and new logl = %f\n", tree_logl, new_logl);

      assert(new_logl - tree_logl > -1e-5); // to avoid numerical errors

      if (q->next->back != test_node) (*interchagnes)++;
    }
    else
    {
      retval = shSupport(treeinfo,
                         shSupportValues,
                         persite_lnl,
                         nBootstrap,
                         shEpsilon,
                         tolerance,
                         brlen_opt_method,
                         bl_min,
                         bl_max,
                         smoothings,
                         lh_epsilon,
                         warning_printed);
    }
  }

  // pb2->back is gonna be the next root
  // so the second time we call the function, we have to update CLVs from the
  // previous root up to the new root
  if (!CORAX_UTREE_IS_TIP(pb1) && retval)
  {
    retval = nni_recursive(treeinfo,
                           pb1,
                           interchagnes,
                           compute_aLRT,
                           shSupportValues,
                           persite_lnl,
                           nBootstrap,
                           shEpsilon,
                           tolerance,
                           brlen_opt_method,
                           bl_min,
                           bl_max,
                           smoothings,
                           lh_epsilon,
                           warning_printed);
  }

  if (!CORAX_UTREE_IS_TIP(pb2) && retval)
  {
    retval = nni_recursive(treeinfo,
                           pb2,
                           interchagnes,
                           compute_aLRT,
                           shSupportValues,
                           persite_lnl,
                           nBootstrap,
                           shEpsilon,
                           tolerance,
                           brlen_opt_method,
                           bl_min,
                           bl_max,
                           smoothings,
                           lh_epsilon,
                           warning_printed);
  }

  return retval;
}

static int sh_support(corax_treeinfo_t       *treeinfo,
                      double                 *sh_support_values,
                      double                **persite_lnl[3],
                      unsigned int            num_bsrep,
                      const unsigned int    **bsrep_site_weights,
                      double                  sh_epsilon,
                      double                  tolerance,
                      int                     brlen_opt_method,
                      double                  bl_min,
                      double                  bl_max,
                      int                     smoothings,
                      double                  lh_epsilon,
                      bool                   *warning_printed)
{
  double LNL0, LNL1, LNL2, second_best_logl;
  int    nSupport          = 0;
  bool   shittySplit       = false;
  bool   non_optimal_split = false;
  double aLRT              = 0;

  double **persite_lnl_0 = persite_lnl[0];
  double **persite_lnl_1 = persite_lnl[1];
  double **persite_lnl_2 = persite_lnl[2];

  corax_unode_t    *q    = treeinfo->root;
  corax_nni_move_t *move = create_move(treeinfo);

  // update clvs and matrices
  LNL0 = corax_treeinfo_compute_loglh_flex(treeinfo, 1, 0);
  if (DEBUG_MODE) printf("Loglh = %f\n", LNL0);

  // First NNI topology
  update_move(treeinfo, move, q, CORAX_UTREE_MOVE_NNI_LEFT);
  LNL1 = apply_move(treeinfo,
                    q,
                    move,
                    brlen_opt_method,
                    bl_min,
                    bl_max,
                    smoothings,
                    lh_epsilon);

  if (LNL1 == CORAX_FAILURE)
  {
    printf("Failed to apply NNI move. \n");
    printf("Corax error number %d \n", corax_errno);
    printf("Corax error msg: %s\n", corax_errmsg);
    delete_move(move);
    return CORAX_FAILURE;
  }

  // calculate persite
  assert(
      fabs(corax_treeinfo_compute_loglh_persite(treeinfo, 1, 0, persite_lnl_1)
           - LNL1)
      < 1e-5);

  second_best_logl = LNL1;

  // undo + testing
  if (!undo_move(treeinfo, q, move))
  {
    corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                    "Error in undoing move that generates a tree of worst "
                    "likelihood ...\n");
    return CORAX_FAILURE;
  }

  // second NNI topology
  update_move(treeinfo, move, q, CORAX_UTREE_MOVE_NNI_RIGHT);
  LNL2 = apply_move(treeinfo,
                    q,
                    move,
                    brlen_opt_method,
                    bl_min,
                    bl_max,
                    smoothings,
                    lh_epsilon);

  if (LNL2 > second_best_logl) second_best_logl = LNL2;

  // calculate persite
  assert(
      fabs(corax_treeinfo_compute_loglh_persite(treeinfo, 1, 0, persite_lnl_2)
           - LNL2)
      < 1e-5);

  // undo + testing
  if (!undo_move(treeinfo, q, move))
  {
    corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                    "Error in undoing move that generates a tree of worst "
                    "likelihood ...\n");
    return CORAX_FAILURE;
  }

  // assert(fabs(compute_local_likelihood_treeinfo(treeinfo, q, 1) - LNL0) <
  // 1e-5);
  delete_move(move);

  if (second_best_logl > LNL0)
  {

    if (second_best_logl - LNL0 > tolerance)
    {

      if (!(*warning_printed))
      {
        printf("\nWARNING!\nYou are trying to calculate the SH-like aLRT "
               "support values for a non NNI-optimal tree topology."
               " For the internal branches that aren't NNI optimal, the "
               "SH-like aLRT metric will be equal to -inf.\n"
               "RECOMMENDED: Since the SH-like aLRT statistics makes more "
               "sense for NNI-optimal topologies, we recommend the user to "
               "first run the "
               "corax_algo_nni_round() function, with a relatively small "
               "tolerance value (e.g. 0.1 or 0.01) to optimize the tree "
               "topology.\n"
               "In case you have already called the NNI-optimization function, "
               "try again with a smaller tolerance value.\n\n");
        *warning_printed = true;
      }
      non_optimal_split = true;
    }
    else { shittySplit = true; }
  }
  else { aLRT = 2 * (LNL0 - second_best_logl); }

  if (shittySplit)
  {
    if (sh_support_values) sh_support_values[q->pmatrix_index] = 0;
    if (DEBUG_MODE)
      printf("Branch %d, shitty split, Support = 0\n", q->pmatrix_index);
    return CORAX_SUCCESS;
  }

  if (non_optimal_split)
  {
    if (sh_support_values) sh_support_values[q->pmatrix_index] = -INFINITY;
    if (DEBUG_MODE) printf("Branch %d, Support = -inf\n", q->pmatrix_index);
    return CORAX_SUCCESS;
  }

  for (unsigned int i = 0; i < num_bsrep; i++)
  {
    double cs[3];
    double _LNL[3] = {0};
    double tmp, diff;

    for (unsigned int p = 0; p < treeinfo->partition_count; p++)
    {
      if (treeinfo->partitions[p])
      {
        const unsigned int *old_wgt = treeinfo->partitions[p]->pattern_weights;
        const unsigned int *new_wgt = bsrep_site_weights[i * treeinfo->partition_count + p];
        for (unsigned int s = 0; s < treeinfo->partitions[p]->sites; s++)
        {
          double ratio = new_wgt[s] / old_wgt[s];
          _LNL[0] += persite_lnl_0[p][s] * ratio;
          _LNL[1] += persite_lnl_1[p][s] * ratio;
          _LNL[2] += persite_lnl_2[p][s] * ratio;
        }
      }
    }

    /* sum up likelihood from all threads */
    corax_treeinfo_parallel_reduce(treeinfo, _LNL, 3, CORAX_REDUCE_SUM);

    cs[0] = _LNL[0] - LNL0;
    cs[1] = _LNL[1] - LNL1;
    cs[2] = _LNL[2] - LNL2;

    if (cs[0] < cs[1])
    {
      tmp   = cs[0];
      cs[0] = cs[1];
      cs[1] = tmp;
    }

    if (cs[1] < cs[2])
    {
      tmp   = cs[1];
      cs[1] = cs[2];
      cs[2] = tmp;
    }

    diff = fabs(cs[0] - cs[1]);
    if (aLRT > 2 * diff + sh_epsilon) { nSupport++; }
  }

  // sh support value
  if (sh_support_values)
  {
    double support                      = (nSupport / (double)num_bsrep);
    sh_support_values[q->pmatrix_index] = support * 100.0;
    if (DEBUG_MODE)
      printf("Branch %d, n_support = %d, Support = %f\n",
             q->pmatrix_index,
             nSupport,
             sh_support_values[q->pmatrix_index]);
  }

  return CORAX_SUCCESS;
}

static int nni_recursive_sh(corax_treeinfo_t       *treeinfo,
                            corax_unode_t          *node,
                            double                 *sh_support_values,
                            double                **persite_lnl[3],
                            unsigned int            num_bsrep,
                            const unsigned int    **bsrep_site_weights,
                            double                  sh_epsilon,
                            double                  tolerance,
                            int                     brlen_opt_method,
                            double                  bl_min,
                            double                  bl_max,
                            int                     smoothings,
                            double                  lh_epsilon,
                            bool                    *warning_printed)
{

  int            retval = CORAX_SUCCESS;
  corax_unode_t *q      = node->back;
  corax_unode_t *pb1    = node->next->back, *pb2 = node->next->next->back;

  if (!CORAX_UTREE_IS_TIP(q))
  {
    corax_treeinfo_set_root(treeinfo, q);
    corax_treeinfo_invalidate_clv(treeinfo, q);
    double tree_logl = corax_treeinfo_compute_loglh_flex(treeinfo, 1, 0);

    retval = sh_support(treeinfo,
                       sh_support_values,
                       persite_lnl,
                       num_bsrep,
                       bsrep_site_weights,
                       sh_epsilon,
                       tolerance,
                       brlen_opt_method,
                       bl_min,
                       bl_max,
                       smoothings,
                       lh_epsilon,
                       warning_printed);
  }

  // pb2->back is gonna be the next root
  // so the second time we call the function, we have to update CLVs from the
  // previous root up to the new root
  if (!CORAX_UTREE_IS_TIP(pb1) && retval)
  {
    retval = nni_recursive_sh(treeinfo,
                           pb1,
                           sh_support_values,
                           persite_lnl,
                           num_bsrep,
                           bsrep_site_weights,
                           sh_epsilon,
                           tolerance,
                           brlen_opt_method,
                           bl_min,
                           bl_max,
                           smoothings,
                           lh_epsilon,
                           warning_printed);
  }

  if (!CORAX_UTREE_IS_TIP(pb2) && retval)
  {
    retval = nni_recursive_sh(treeinfo,
                           pb2,
                           sh_support_values,
                           persite_lnl,
                           num_bsrep,
                           bsrep_site_weights,
                           sh_epsilon,
                           tolerance,
                           brlen_opt_method,
                           bl_min,
                           bl_max,
                           smoothings,
                           lh_epsilon,
                           warning_printed);
  }

  return retval;
}


static double algo_nni_round(corax_treeinfo_t *treeinfo,
                             double            tolerance,
                             int               brlen_opt_method,
                             double            bl_min,
                             double            bl_max,
                             int               smoothings,
                             double            lh_epsilon,
                             bool              print_in_console)
{

  int          retval             = CORAX_SUCCESS;
  unsigned int nniRounds          = 0;
  unsigned int total_interchanges = 0;
  unsigned int interchanges;
  double       tree_logl, diff, new_logl;
  unsigned int tip_index;

  corax_utree_t *tree         = treeinfo->tree;
  corax_unode_t *initial_root = treeinfo->root;
  corax_unode_t *start_node;

  if (print_in_console || DEBUG_MODE) printf("\n\nNNI Round:\n");

  clock_t begin = clock(); // calculate nni time

  // define start node (the back node from leaf 0)
  tip_index  = 0;
  start_node = tree->nodes[tip_index]->back;
  if (!CORAX_UTREE_IS_TIP(start_node->back))
  {
    corax_set_error(
        CORAX_NNI_ROUND_LEAF_ERROR,
        "The NNI round is defined to start from a leaf node, and in particular "
        "this node is the first one in the `**nodes` list of your tree. "
        "The first elements of this list are expected to be leaf nodes, but in "
        "this case, this does not seem to hold.\n");
    printf("Eroor: %s\n", corax_errmsg);
    return CORAX_FAILURE;
  }

  // set root, calculate likelihood, update CLVs in general
  treeinfo->root = start_node;
  tree_logl      = corax_treeinfo_compute_loglh(treeinfo, 0);
  if (DEBUG_MODE) printf("Start logl = %f\n", tree_logl);

  // update splits for constraint check
  corax_treeinfo_constraint_update_splits(treeinfo);

  do {
    interchanges = 0;

    // perform nni moves
    retval = nni_recursive(treeinfo,
                           start_node,
                           &interchanges,
                           false,
                           NULL,
                           NULL,
                           0,
                           0.,
                           tolerance,
                           brlen_opt_method,
                           bl_min,
                           bl_max,
                           smoothings,
                           lh_epsilon,
                           NULL);

    if (retval == CORAX_FAILURE)
    {
      printf("\nSomething went wrong, the return value of NNI round is "
             "CORAX_FAILURE. Exit...\n");
      return CORAX_FAILURE;
    }
    total_interchanges += interchanges;

    // during the nni round, the root of the tree has changed
    // so we first traverse tree to its root, update clvs, and then set again
    // root equal to start node.
    start_node = tree->nodes[tip_index]->back;
    corax_treeinfo_set_root(treeinfo, start_node);
    new_logl = corax_treeinfo_compute_loglh_flex(treeinfo, 1, 0);

    diff = new_logl - tree_logl;

    // we check if (diff < -1e-7) and not (diff < 0) to avoid numerical errors
    if (diff < -1e-7)
    {
      printf("ERROR! It seems that tree likelihood after the NNI round is "
             "smaller than the initial likelihood\n");
      printf("Init logl = %.5f and new logl = %.5f \n", tree_logl, new_logl);
      printf("Round Interchanges = %d\n", interchanges);
      corax_set_error(CORAX_NNI_DIFF_NEGATIVE_ERROR,
                      "The logl after an NNI round was smaller than the "
                      "initial likelihood.\n");
      return CORAX_FAILURE;
    }

    tree_logl = new_logl;
    nniRounds++;

    if (print_in_console || DEBUG_MODE)
    {
      printf("Round %d, interchanges = %d, logl = %.5f\n",
             nniRounds,
             interchanges,
             tree_logl);
    }
  } while ((diff > tolerance || nniRounds < 10) && interchanges != 0);

  clock_t end        = clock();
  double  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  if (print_in_console || DEBUG_MODE)
  {
    printf("\nCompleted.\nTotal interchanges = %d, Time = %.3f secs\n",
           total_interchanges,
           time_spent);
    printf("Final logl = %.5f\n\n", tree_logl);
  }

  if (start_node != initial_root)
  {

    treeinfo->root = initial_root;
    new_logl       = corax_treeinfo_compute_loglh(treeinfo, 0);
    // new_logl = compute_local_likelihood_treeinfo(treeinfo, initial_root, 0);
    assert(fabs(new_logl - tree_logl) < 1e-5);
  }

  return tree_logl;
}

CORAX_EXPORT double corax_algo_nni_round(corax_treeinfo_t *treeinfo,
                                         double            tolerance,
                                         int               brlen_opt_method,
                                         double            bl_min,
                                         double            bl_max,
                                         int               smoothings,
                                         double            lh_epsilon,
                                         bool              print_in_console)
{
  double loglh;

  // Integrity check
  if (!corax_utree_check_integrity(treeinfo->tree))
  {
    corax_set_error(CORAX_NNI_ROUND_INTEGRITY_ERROR,
                    "Tree is not consistent...\n");
    return CORAX_FAILURE;
  }

  // If the tree is a triplet, NNI moves cannot be done
  if (treeinfo->tip_count == 3)
  {
    corax_set_error(CORAX_NNI_ROUND_TRIPLET_ERROR,
                    "Tree is a triplet, NNI moves are not allowed ...\n");
    return CORAX_FAILURE;
  }

  /* make sure initial topology is compatible with constraint */
  if (!corax_treeinfo_constraint_check_current(treeinfo))
  {
    corax_set_error(CORAX_ERROR_INVALID_TREE,
                    "Constraint check failed before SPR round!");
    return CORAX_FAILURE;
  }

  loglh = algo_nni_round(treeinfo,
                         tolerance,
                         brlen_opt_method,
                         bl_min,
                         bl_max,
                         smoothings,
                         lh_epsilon,
                         print_in_console);

  /* make sure final topology is compatible with constraint */
  if (!corax_treeinfo_constraint_check_current(treeinfo))
  {
#if DEBUG_MODE
    corax_utree_show_ascii(treeinfo->root, CORAX_UTREE_SHOW_LABEL | CORAX_UTREE_SHOW_BRANCH_LENGTH |
                                           CORAX_UTREE_SHOW_CLV_INDEX);
#endif
    corax_set_error(CORAX_ERROR_INVALID_TREE,
                     "Constraint check failed after SPR round!");
    return CORAX_FAILURE;
  }

  return loglh;
}

CORAX_EXPORT int corax_shSupport_values(corax_treeinfo_t *treeinfo,
                                        double            tolerance,
                                        double           *shSupportValues,
                                        int               nBootstrap,
                                        double            shEpsilon,
                                        int               brlen_opt_method,
                                        double            bl_min,
                                        double            bl_max,
                                        int               smoothings,
                                        double            lh_epsilon,
                                        bool              print_in_console)
{

  // Integrity check
  if (!corax_utree_check_integrity(treeinfo->tree))
  {
    corax_set_error(CORAX_NNI_ROUND_INTEGRITY_ERROR,
                    "Tree is not consistent...\n");
    return CORAX_FAILURE;
  }

  // If the tree is a triplet, NNI moves cannot be done
  if (treeinfo->tip_count == 3)
  {
    corax_set_error(CORAX_NNI_ROUND_TRIPLET_ERROR,
                    "Tree is a triplet, NNI moves are not allowed ...\n");
    return CORAX_FAILURE;
  }

  int          retval = CORAX_SUCCESS;
  double       tree_logl, test_logl;
  unsigned int tip_index = 0;

  corax_utree_t *tree         = treeinfo->tree;
  corax_unode_t *initial_root = treeinfo->root;
  corax_unode_t *start_node;

  if (print_in_console || DEBUG_MODE)
    printf("Calculating SH-like aLRT statistics ... \n");

  int shalrt_size = 2 * ((int)treeinfo->tip_count) - 3;
  for (int i = 0; i < shalrt_size; i++) { shSupportValues[i] = -INFINITY; }

  double **persite_lnl;
  persite_lnl =
      (double **)malloc(sizeof(double *) * (treeinfo->partition_count));

  for (unsigned int i = 0; i < treeinfo->partition_count; i++)
    persite_lnl[i] =
        (double *)malloc(sizeof(double) * treeinfo->partitions[i]->sites);

  start_node = tree->nodes[tip_index]->back;
  corax_treeinfo_set_root(treeinfo, start_node);
  tree_logl = corax_treeinfo_compute_loglh_persite(treeinfo, 0, 0, persite_lnl);

  // perform nni moves
  bool warning_printed = false;
  retval               = nni_recursive(treeinfo,
                         start_node,
                         NULL,
                         true,
                         shSupportValues,
                         persite_lnl,
                         nBootstrap,
                         shEpsilon,
                         tolerance,
                         brlen_opt_method,
                         bl_min,
                         bl_max,
                         smoothings,
                         lh_epsilon,
                         &warning_printed);

  for (unsigned i = 0; i < treeinfo->partition_count; i++) free(persite_lnl[i]);

  free(persite_lnl);

  if (retval == CORAX_FAILURE)
  {
    printf("\nSomething went wrong in the calculations of SH-like aLRT "
           "statistics. Exit...\n");
    return CORAX_FAILURE;
  }

  if (print_in_console || DEBUG_MODE) printf("Done!\n\n");

  if (start_node != initial_root)
  {

    corax_treeinfo_set_root(treeinfo, initial_root);

    test_logl = corax_treeinfo_compute_loglh(treeinfo, 0);

    if (fabs(tree_logl - test_logl) > 1e-5)
    { // this if is to get rid of "unused variable" warning
      printf("\nLikelihoods before and after the calculation of SH-like aLRT "
             "statistics are different. Exit ... \n");
      return CORAX_FAILURE;
    }
  }

  return retval;
}


CORAX_EXPORT int corax_algo_sh_support(corax_treeinfo_t        *treeinfo,
                                        double                  tolerance,
                                        double                 *sh_support_values,
                                        unsigned int            num_bsrep,
                                        const unsigned int    **bsrep_site_weights,
                                        double                  sh_epsilon,
                                        int                     brlen_opt_method,
                                        double                  bl_min,
                                        double                  bl_max,
                                        int                     smoothings,
                                        double                  lh_epsilon)
{
  // Integrity check
  if (!corax_utree_check_integrity(treeinfo->tree))
  {
    corax_set_error(CORAX_NNI_ROUND_INTEGRITY_ERROR,
                    "Tree is not consistent...\n");
    return CORAX_FAILURE;
  }

  // If the tree is a triplet, NNI moves cannot be done
  if (treeinfo->tip_count == 3)
  {
    corax_set_error(CORAX_NNI_ROUND_TRIPLET_ERROR,
                    "Tree is a triplet, NNI moves are not allowed ...\n");
    return CORAX_FAILURE;
  }

  int          retval = CORAX_SUCCESS;
  double       tree_logl, test_logl;
  unsigned int tip_index = 0;

  corax_utree_t *tree         = treeinfo->tree;
  corax_unode_t *initial_root = treeinfo->root;
  corax_unode_t *start_node;

  if (DEBUG_MODE)
    printf("Calculating SH-like aLRT statistics ... \n");

  if (sh_support_values)
  {
    int shalrt_size = 2 * ((int)treeinfo->tip_count) - 3;
    for (int i = 0; i < shalrt_size; i++) { sh_support_values[i] = -INFINITY; }
  }

  double **persite_lnl[3];

  for (unsigned int k = 0; k < 3; ++k)
  {
    persite_lnl[k] =
        (double **)malloc(sizeof(double *) * (treeinfo->partition_count));

    for (unsigned int i = 0; i < treeinfo->partition_count; i++)
    {
      persite_lnl[k][i] = treeinfo->partitions[i] ?
          (double *)malloc(sizeof(double) * treeinfo->partitions[i]->sites) : NULL;
    }
  }

  start_node = tree->nodes[tip_index]->back;
  corax_treeinfo_set_root(treeinfo, start_node);
  tree_logl = corax_treeinfo_compute_loglh_persite(treeinfo, 0, 0, persite_lnl[0]);

  // perform nni moves
  bool warning_printed = false;
  retval = nni_recursive_sh(treeinfo,
                         start_node,
                         sh_support_values,
                         persite_lnl,
                         num_bsrep,
                         bsrep_site_weights,
                         sh_epsilon,
                         tolerance,
                         brlen_opt_method,
                         bl_min,
                         bl_max,
                         smoothings,
                         lh_epsilon,
                         &warning_printed);

  // cleanup
  for (unsigned int k = 0; k < 3; ++k)
  {
    for (unsigned i = 0; i < treeinfo->partition_count; i++)
      free(persite_lnl[k][i]);

    free(persite_lnl[k]);
  }

  if (retval == CORAX_FAILURE)
  {
    printf("\nSomething went wrong in the calculations of SH-like aLRT "
           "statistics. Exit...\n");
    return CORAX_FAILURE;
  }

  if (DEBUG_MODE) printf("Done!\n\n");

  if (start_node != initial_root)
  {

    corax_treeinfo_set_root(treeinfo, initial_root);

    test_logl = corax_treeinfo_compute_loglh(treeinfo, 0);

    if (fabs(tree_logl - test_logl) > 1e-5)
    { // this if is to get rid of "unused variable" warning
      printf("\nLikelihoods before and after the calculation of SH-like aLRT "
             "statistics are different. Exit ... \n");
      return CORAX_FAILURE;
    }
  }

  return retval;
}

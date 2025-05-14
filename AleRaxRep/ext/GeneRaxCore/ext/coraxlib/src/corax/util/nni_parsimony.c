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

#include "nni_parsimony.h"

#define DEBUG_MODE 0

typedef struct corax_nni_move_s
{

  int           type;
  unsigned int  partition_count;

} corax_nni_move_t;

typedef struct
{
  int clv_valid;
} node_info_t;

/* a callback function for performing a full traversal */
static int cb_full_traversal(corax_unode_t *node)
{
  CORAX_UNUSED(node);
  return CORAX_SUCCESS;
}

static int undo_move(corax_unode_t    *node,
                     corax_nni_move_t *move)
{

    int retval = CORAX_SUCCESS;

    retval = corax_utree_nni(node, move->type, NULL);

    /* check this 
    corax_treeinfo_invalidate_pmatrix(treeinfo, node);
    corax_treeinfo_invalidate_pmatrix(treeinfo, node->next);
    corax_treeinfo_invalidate_pmatrix(treeinfo, node->next->next);
    corax_treeinfo_invalidate_pmatrix(treeinfo, node->back->next);
    corax_treeinfo_invalidate_pmatrix(treeinfo, node->back->next->next);

    
    corax_treeinfo_invalidate_clv(treeinfo, node);
    corax_treeinfo_invalidate_clv(treeinfo, node->back);
    */
    return retval;
}

static unsigned int compute_parsimony_score(corax_parsimony_t **list,
                                            corax_unode_t* node,
                                            unsigned int tips_count,
                                            unsigned int partition_count,
                                            int (*cbtrav)(corax_unode_t *))
{
    corax_unode_t **      travbuffer;
    corax_pars_buildop_t *parsops;
    unsigned int ops_count = 0, 
                traversal_size = 0;

    unsigned int j;

    travbuffer =
      (corax_unode_t **)malloc((2 * tips_count - 2) * sizeof(corax_unode_t *));
    
    parsops = (corax_pars_buildop_t *)malloc((tips_count - 2)
                                    * sizeof(corax_pars_buildop_t));
    
    if (!corax_utree_traverse(node,
                              CORAX_TREE_TRAVERSE_POSTORDER,
                              cbtrav,
                              travbuffer,
                              &traversal_size))
      assert(0);

    /* create parsimony operations */
    corax_utree_create_pars_buildops(
        travbuffer, traversal_size, parsops, &ops_count);

    /* compute the costs for each parsimony partition */
    unsigned int score = 0;
    for (j = 0; j < partition_count; ++j)
    {
        /* update parsimony vectors */
        corax_fastparsimony_update_vectors(list[j], parsops, ops_count);

        /* get parsimony score */
        score += corax_fastparsimony_edge_score(
            list[j], node->node_index, node->back->node_index);
    }

    free(travbuffer);
    free(parsops);

    return score;
}

static unsigned int apply_move(corax_unode_t    *node,
                            corax_parsimony_t **list,
                            corax_nni_move_t *move,
                            unsigned int tips_count,
                            unsigned int partition_count)
{

    if (move->type == CORAX_UTREE_MOVE_NNI_LEFT
        || move->type == CORAX_UTREE_MOVE_NNI_RIGHT)
    {

        if (CORAX_UTREE_IS_TIP(node) || CORAX_UTREE_IS_TIP(node->back))
        {
            /* invalid move */
            corax_set_error(CORAX_TREE_ERROR_NNI_LEAF,
                            "Attempting to apply NNI on a leaf branch\n");
            return CORAX_FAILURE;
        }

        corax_utree_nni(node, move->type, NULL);
        //corax_treeinfo_invalidate_clv(treeinfo, node);
        //corax_treeinfo_invalidate_clv(treeinfo, node->back);

        unsigned int new_score = compute_parsimony_score(list, 
                                                        node, 
                                                        tips_count, 
                                                        partition_count, 
                                                        cb_full_traversal);

        if (new_score == CORAX_FAILURE)
        {
            printf(
                "Something went wrong with the branch-length optimization. Exit..\n");
            return CORAX_FAILURE;
        }

        return new_score;
    
    } else {
    
        /* invalid move */
        corax_set_error(CORAX_TREE_ERROR_NNI_INVALID_MOVE,
                        "Invalid NNI move type\n");
        return CORAX_FAILURE;
    }
}

CORAX_EXPORT unsigned int
    corax_algo_nni_parsimony_local(corax_unode_t* q,
                                    corax_parsimony_t **list,
                                    unsigned int pars_score,
                                    unsigned int partition_count,
                                    unsigned int tips_count)
{
    unsigned int new_score;
    corax_nni_move_t *move = (corax_nni_move_t *)malloc(sizeof(corax_nni_move_t));

    // set up move
    move->partition_count = partition_count;
    move->type = CORAX_UTREE_MOVE_NNI_LEFT;
    
    new_score = apply_move(q, 
                        list, 
                        move,
                        tips_count,
                        partition_count);
    
    if (new_score == CORAX_FAILURE)
    {
        printf("Failed to apply NNI move. \n");
        printf("Corax error number %d \n", corax_errno);
        printf("Corax error msg: %s\n", corax_errmsg);
        free(move);
        return CORAX_FAILURE;
    }

    if (new_score < pars_score)
    {

        pars_score = new_score;

        // update move structure
        move->type = CORAX_UTREE_MOVE_NNI_RIGHT;
        
        new_score = apply_move(q, 
                            list, 
                            move,
                            tips_count,
                            partition_count);

        if (new_score == CORAX_FAILURE)
        {
            printf("Failed to apply NNI move. \n");
            printf("Error number %d \n", corax_errno);
            printf("Error msg: %s\n", corax_errmsg);
            free(move);
            return CORAX_FAILURE;
        }

        if (new_score < pars_score)
        {
            pars_score = new_score;
            free(move);
        }
        else
        {

            if (!undo_move(q, move))
            {
                corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                                "Error in undoing move that generates a tree of worst "
                                "likelihood ...\n");
                return CORAX_FAILURE;
            }

            // assert(fabs(compute_local_likelihood_treeinfo(treeinfo, q, 1) -
            // tree_logl) < 1e-5);
            free(move);
        }
    }
    else
    {

        if (!undo_move(q, move))
        {
            corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                            "Error in undoing move that generates a tree of worst "
                            "likelihood ...\n");
            return CORAX_FAILURE;
        }

        // if(DEBUG_MODE) assert(fabs(compute_local_likelihood_treeinfo(treeinfo, q,
        // 1) - tree_logl) < 1e-5);

        // setup_move_info(move, q, new_logl, CORAX_UTREE_MOVE_NNI_RIGHT);
        move->type = CORAX_UTREE_MOVE_NNI_RIGHT;
        new_score   = apply_move(q, 
                                list, 
                                move, 
                                tips_count, 
                                partition_count);

        if (new_score == CORAX_FAILURE)
        {
            printf("Failed to apply NNI move. \n");
            printf("Error number %d \n", corax_errno);
            printf("Error msg: %s\n", corax_errmsg);
            free(move);
            return CORAX_FAILURE;
        }

        if (new_score < pars_score)
        {
            pars_score = new_score;
            free(move);
        }
        else
        {
            if (!undo_move(q, move))
            {
                corax_set_error(CORAX_NNI_ROUND_UNDO_MOVE_ERROR,
                                "Error in undoing move that generates a tree of worst "
                                "likelihood ...\n");
                return CORAX_FAILURE;
            }

            // assert(fabs(compute_local_likelihood_treeinfo(treeinfo, q, 1) -
            // tree_logl) < 1e-5);
            free(move);
        }
    }

    return pars_score;
}

static int nni_parsimony_recursive(corax_parsimony_t ** list,
                                    corax_unode_t    *node,
                                    unsigned int     *interchagnes,
                                    unsigned int      partition_count,
                                    unsigned int      tip_count)
{

    int            retval = CORAX_SUCCESS;
    corax_unode_t *q      = node->back;
    corax_unode_t *pb1 = node->next->back, *pb2 = node->next->next->back;

    int pars_score;
    int new_score;

    if (!CORAX_UTREE_IS_TIP(q))
    {

        corax_unode_t *test_node = q->next->back;
        pars_score = compute_parsimony_score(list, q, tip_count, partition_count, cb_full_traversal);

        new_score = corax_algo_nni_parsimony_local(q, 
                                                list, 
                                                pars_score, 
                                                partition_count, 
                                                tip_count);

        if (new_score == CORAX_FAILURE) return CORAX_FAILURE;
        if (DEBUG_MODE) printf("Old score = %d and new score = %d\n", pars_score, new_score);

        assert(pars_score - new_score >= 0); // to avoid numerical errors
        if (q->next->back != test_node) (*interchagnes)++;

    }
    // pb2->back is gonna be the next root
    // so the second time we call the function, we have to update CLVs from the
    // previous root up to the new root
    if (!CORAX_UTREE_IS_TIP(pb1) && retval)
        retval = nni_parsimony_recursive(list,
                                pb1,
                                interchagnes,
                                partition_count,
                                tip_count);

    if (!CORAX_UTREE_IS_TIP(pb2) && retval)
        retval = nni_parsimony_recursive(list,
                                pb2,
                                interchagnes,
                                partition_count,
                                tip_count);

    return retval;
}

CORAX_EXPORT unsigned int
    corax_algo_nni_round_parsimony(corax_utree_t * tree,
                                corax_parsimony_t **list,
                                unsigned int partition_count,
                                unsigned int * score)
{
    
    // If the tree is a triplet, NNI moves cannot be done
    if (tree->tip_count == 3)
    {
        corax_set_error(CORAX_NNI_ROUND_TRIPLET_ERROR,
                        "Parsimony tree is a triplet, NNI moves are not allowed ...\n");
        return CORAX_FAILURE;
    }

    // corax_unode_t *initial_root = tree->vroot;
    corax_unode_t *start_node;
    // define start node (the back node from leaf 0)
    unsigned int tip_index  = 0;
    start_node = tree->nodes[tip_index]->back;
    assert(!CORAX_UTREE_IS_TIP(start_node));

    // first thing -> complete the function
    // tree->vroot = start_node;
    unsigned int pars_score = compute_parsimony_score(list, start_node, tree->tip_count, partition_count, cb_full_traversal);
    if (DEBUG_MODE) printf("Start score = %d and init score = %d \n", pars_score, *score);
    assert(pars_score == (*score));

    int retval;
    unsigned int interchanges, 
                new_score,  
                total_interchanges = 0,
                nniRounds = 0;
    
    int diff = 0;
    do {

        interchanges = 0;

        // perform nni moves
        retval = nni_parsimony_recursive(list, 
                                        start_node,
                                        &interchanges,
                                        partition_count,
                                        tree->tip_count);

        if (retval == CORAX_FAILURE)
        {
            printf("\nSomething went wrong, the return value of parsimony-NNI round is "
                    "CORAX_FAILURE. Exit...\n");
            return CORAX_FAILURE;
        }
       
        total_interchanges += interchanges;

        // during the nni round, the root of the tree has changed
        // so we first traverse tree to its root, update clvs, and then set again
        // root equal to start node.
        start_node = tree->nodes[tip_index]->back;
        new_score = compute_parsimony_score(list, start_node, tree->tip_count, partition_count, cb_full_traversal);
        diff = pars_score - new_score;

        // we check if (diff < -1e-7) and not (diff < 0) to avoid numerical errors
        if (diff < 0)
        {
            printf("ERROR! \n");
            // write something in the error here
            return CORAX_FAILURE;
        }

        pars_score = new_score;
        nniRounds++;

        if (DEBUG_MODE)
            printf("Round %d, interchanges = %d, score = %d\n",
                    nniRounds,
                    interchanges,
                    pars_score);
        
        //getchar();

    } while (diff > 0);

    return 0;
}
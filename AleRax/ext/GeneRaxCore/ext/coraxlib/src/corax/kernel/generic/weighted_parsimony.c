/*
    Copyright (C) 2015 Tomas Flouri

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

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "corax/corax.h"

CORAX_EXPORT double
corax_parsimony_build(corax_parsimony_t *         pars,
                      const corax_pars_buildop_t *operations,
                      unsigned int                count)
{
  unsigned int                i, j, k, n;
  const corax_pars_buildop_t *op;

  unsigned int sites  = pars->sites;
  unsigned int states = pars->states;
  double       minimum;

  double *score_buffer;
  double *child1_score_buffer;
  double *child2_score_buffer;

  double *score_matrix = pars->score_matrix;

  /* Implementation of the 'minimum mutation trees' algorithm by David Sankoff.
     For more information see:

  Sankoff D: Minimal Mutation Trees of Sequences. SIAM Appl Math 28:35-32, 1975
  doi: 10.1137/0128004

  Sankoff D, Rousseau P: Locating the Vertices of a Steiner Tree in Arbitrary
  Metric Space. Math Program 1975, 9:240-246.
  doi: 10.1007/BF01681346

  */

  /* iterate through inner nodes (post-order traversal) */
  for (i = 0; i < count; ++i)
  {
    op = &(operations[i]);

    /* get parent score buffer */
    score_buffer = pars->sbuffer[op->parent_score_index];

    /* get child1 score buffer if it's not a tip */
    child1_score_buffer = pars->sbuffer[op->child1_score_index];

    /* get child2 score buffer if it's not a tip */
    child2_score_buffer = pars->sbuffer[op->child2_score_index];

    /* iterate through sites */
    for (j = 0; j < sites; ++j)
    {

      for (n = 0; n < states; ++n)
      {
        /* process child 1 */

        minimum = child1_score_buffer[0] + score_matrix[n];
        for (k = 1; k < states; ++k)
          minimum = fmin(child1_score_buffer[k] + score_matrix[k * states + n],
                         minimum);

        score_buffer[n] = minimum;

        /* process child 2 */

        minimum = child2_score_buffer[0] + score_matrix[n];
        for (k = 1; k < states; ++k)
          minimum = fmin(child2_score_buffer[k] + score_matrix[k * states + n],
                         minimum);

        score_buffer[n] += minimum;
      }

      score_buffer += states;
      child1_score_buffer += states;
      child2_score_buffer += states;
    }
  }

  /* get the sum of minimum scores for each site of the last node in the
     postorder traversal */
  op = &(operations[count - 1]);

  return corax_parsimony_score(pars, op->parent_score_index);
}

CORAX_EXPORT double corax_parsimony_score(corax_parsimony_t *pars,
                                          unsigned int       score_buffer_index)
{
  unsigned int i, j, k;
  unsigned int states = pars->states;
  unsigned int sites  = pars->sites;
  double       sum    = 0;
  double       minimum;

  double *score_buffer = pars->sbuffer[score_buffer_index];

  for (k = 0, i = 0; i < sites; ++i)
  {
    minimum = score_buffer[k++];
    for (j = 1; j < states; ++j) minimum = fmin(score_buffer[k++], minimum);

    sum += minimum;
  }

  return sum;
}

CORAX_EXPORT void
corax_parsimony_reconstruct(corax_parsimony_t *       pars,
                            const corax_state_t *     map,
                            const corax_pars_recop_t *operations,
                            unsigned int              count)
{
  unsigned int i, j, n;
  unsigned int revmap[256];

  double *      score_buffer;
  unsigned int *ancestral_buffer;
  double *      parent_score_buffer;
  unsigned int *parent_ancestral_buffer;
  unsigned int  minindex;

  unsigned int states = pars->states;

  const corax_pars_recop_t *op;

  for (i = 0; i < 256; ++i) revmap[i] = 0;
  for (i = 0; i < 256; ++i)
  {
    if (CORAX_STATE_POPCNT(map[i]) == 1)
    {
      revmap[CORAX_STATE_CTZ(map[i])] = i;
    }
  }

  /* start from root of given subtree */
  op               = &(operations[0]);
  score_buffer     = pars->sbuffer[op->node_score_index];
  ancestral_buffer = pars->anc_states[op->node_ancestral_index];
  for (n = 0; n < pars->sites; ++n)
  {
    minindex = 0;
    for (i = 1; i < pars->states; ++i)
    {
      if (score_buffer[n * states + i] < score_buffer[n * states + minindex])
        minindex = i;
    }
    ancestral_buffer[n] = revmap[minindex];
  }

  /* continue with preorder traversal */
  for (i = 1; i < count; ++i)
  {
    op = &(operations[i]);

    /* get node score and ancestral buffer */
    parent_score_buffer     = pars->sbuffer[op->parent_score_index];
    parent_ancestral_buffer = pars->anc_states[op->parent_ancestral_index];

    /* get node score and ancestral buffer */
    score_buffer     = pars->sbuffer[op->node_score_index];
    ancestral_buffer = pars->anc_states[op->node_ancestral_index];

    for (n = 0; n < pars->sites; ++n)
    {
      /* compute index with minimum value */
      minindex = 0;
      for (j = 1; j < pars->states; ++j)
      {
        if (score_buffer[n * states + j] < score_buffer[n * states + minindex])
          minindex = j;
      }

      double parent_val = parent_score_buffer
          [n * states + CORAX_STATE_CTZ(map[parent_ancestral_buffer[n]])];

      if (score_buffer[n * states + minindex] + 1 > parent_val)
        ancestral_buffer[n] = parent_ancestral_buffer[n];
      else
        ancestral_buffer[n] = revmap[minindex];
    }
  }
}

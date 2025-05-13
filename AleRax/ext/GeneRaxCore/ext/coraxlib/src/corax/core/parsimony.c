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

CORAX_EXPORT int corax_set_parsimony_sequence(corax_parsimony_t *  pars,
                                              unsigned int         tip_index,
                                              const corax_state_t *map,
                                              const char *         sequence)
{
  corax_state_t c;
  unsigned int  i, j;

  unsigned int states   = pars->states;
  double *     tipstate = pars->sbuffer[tip_index];

  double *score_matrix = pars->score_matrix;

  /* set infinity as the highest score in the matrix plus one */
  double inf = score_matrix[0];
  for (i = 1; i < states * states; ++i)
    if (score_matrix[i] > inf) inf = score_matrix[i];
  inf++;

  for (i = 0; i < pars->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      corax_set_error(CORAX_ERROR_TIPDATA_ILLEGALSTATE,
                      "Illegal state code in tip \"%c\"",
                      sequence[i]);
      printf("%s\n", corax_errmsg);
      return CORAX_FAILURE;
    }

    for (j = 0; j < states; ++j)
    {
      if (c & 1)
        tipstate[j] = 0;
      else
        tipstate[j] = inf;
      c >>= 1;
    }

    tipstate += states;
  }

  return CORAX_SUCCESS;
}

CORAX_EXPORT corax_parsimony_t *
             corax_parsimony_create(unsigned int  tips,
                                    unsigned int  states,
                                    unsigned int  sites,
                                    const double *score_matrix,
                                    unsigned int  score_buffers,
                                    unsigned int  ancestral_buffers)
{
  unsigned int i;

  /* create parsimony instance */
  corax_parsimony_t *pars =
      (corax_parsimony_t *)calloc(1, sizeof(corax_parsimony_t));
  if (!pars)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  /* store passed parameters */
  pars->tips              = tips;
  pars->states            = states;
  pars->sites             = sites;
  pars->score_buffers     = score_buffers;
  pars->ancestral_buffers = ancestral_buffers;

  /* create and copy a states*states scoring matrix */
  pars->score_matrix = (double *)calloc(states * states, sizeof(double));
  if (!pars->score_matrix)
  {
    corax_parsimony_destroy(pars);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for scoring matrix.");
    return NULL;
  }
  memcpy(pars->score_matrix, score_matrix, states * states * sizeof(double));

  /* create requested score buffers */
  pars->sbuffer = (double **)calloc(score_buffers + tips, sizeof(double *));
  if (!pars->sbuffer)
  {
    corax_parsimony_destroy(pars);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for score buffers.");
    return NULL;
  }
  for (i = 0; i < score_buffers + tips; ++i)
  {
    pars->sbuffer[i] = (double *)calloc(sites * states, sizeof(double *));
    if (!pars->sbuffer[i])
    {
      corax_parsimony_destroy(pars);
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory for score buffers.");
      return NULL;
    }
  }

  /* create ancestral buffers */
  pars->anc_states =
      (unsigned int **)calloc(tips + ancestral_buffers, sizeof(unsigned int *));
  if (!pars->anc_states)
  {
    corax_parsimony_destroy(pars);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for score buffers.");
    return NULL;
  }
  for (i = tips; i < ancestral_buffers + tips; ++i)
  {
    pars->anc_states[i] = (unsigned int *)calloc(sites, sizeof(unsigned int));
    if (!pars->anc_states[i])
    {
      corax_parsimony_destroy(pars);
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory for score buffers.");
      return NULL;
    }
  }

  return pars;
}

CORAX_EXPORT void corax_parsimony_destroy(corax_parsimony_t *parsimony)
{
  unsigned int i;
  unsigned int nodes_count = 0;

  if (!parsimony) return;

  nodes_count = parsimony->tips + 3 * parsimony->inner_nodes;

  /* deallocate fast parsimony structures */
  if (parsimony->packedvector)
  {
    for (i = 0; i < nodes_count; ++i)
      corax_aligned_free(parsimony->packedvector[i]);
    free(parsimony->packedvector);
  }

  if (parsimony->node_cost) free(parsimony->node_cost);

  if (parsimony->informative) free(parsimony->informative);

  /* if available, deallocate structures for weighted parsimony */

  /* score buffers */
  if (parsimony->sbuffer)
  {
    for (i = 0; i < parsimony->score_buffers + parsimony->tips; ++i)
      free(parsimony->sbuffer[i]);
    free(parsimony->sbuffer);
  }

  /* ancestral state buffers */
  if (parsimony->anc_states)
  {
    for (i = parsimony->tips;
         i < parsimony->ancestral_buffers + parsimony->tips;
         ++i)
      free(parsimony->anc_states[i]);
    free(parsimony->anc_states);
  }

  /* scoring matrix */
  if (parsimony->score_matrix) free(parsimony->score_matrix);

  free(parsimony);
}


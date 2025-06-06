/*
 Copyright (C) 2016 Alexey Kozlov

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

#include <string.h>

#include "corax/model/modutil.h"

/*                                       AA CC GG TT AC AG AT CG CT GT */
static const double gt_rates_equal_sm[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, /* AA */
                                           0, 0, 1, 0, 0, 1, 1, 0,    /* CC */
                                           0, 0, 1, 0, 1, 0, 1,       /* GG */
                                           0, 0, 1, 0, 1, 1,          /* TT */
                                           1, 1, 1, 1, 0,             /* AC */
                                           1, 1, 0, 1,                /* AG */
                                           0, 1, 1,                   /* AT */
                                           1, 1,                      /* CG */
                                           1};                        /* CT */

/*                                      AA CC GG TT AC AG AT CG CT GT */
static const double gt_rates_equal[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, /* AA */
                                        1, 1, 1, 1, 1, 1, 1, 1,    /* CC */
                                        1, 1, 1, 1, 1, 1, 1,       /* GG */
                                        1, 1, 1, 1, 1, 1,          /* TT */
                                        1, 1, 1, 1, 1,             /* AC */
                                        1, 1, 1, 1,                /* AG */
                                        1, 1, 1,                   /* AT */
                                        1, 1,                      /* CG */
                                        1};                        /* CT */

/*                                      AA   CC   GG   TT   AC   AG   AT   CG
 * CT   GT */
static const double gt_freqs_equal[] = {
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

/*                                 A  C  G  T              */
// static int gt_sym_freq_equal[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
// static int gt_sym_freq_free[]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

/*                                  AA CC GG TT AC AG AT CG CT GT          */
static int gt_sym_rate_free_sm[] = {0,  0,  0,  1,  2,  3,  0, 0, 0, /* AA */
                                    0,  0,  4,  0,  0,  5,  6, 0,    /* CC */
                                    0,  0,  7,  0,  8,  0,  9,       /* GG */
                                    0,  0,  10, 0,  11, 12,          /* TT */
                                    13, 14, 15, 16, 0,               /* AC */
                                    17, 18, 0,  19,                  /* AG */
                                    0,  20, 21,                      /* AT */
                                    22, 23,                          /* CG */
                                    24};                             /* CT */

/* A-C: 1, A-G: 2, A-T: 3, C-G: 4, C-T: 5, G-T: 6, others: 0 */
/*                               AA CC GG TT AC AG AT CG CT GT          */
static int gt_sym_rate_dna4[] = {0, 0, 0, 1, 2, 3, 0, 0, 0, /* AA */
                                 0, 0, 1, 0, 0, 4, 5, 0,    /* CC */
                                 0, 0, 2, 0, 4, 0, 6,       /* GG */
                                 0, 0, 3, 0, 5, 6,          /* TT */
                                 4, 5, 2, 3, 0,             /* AC */
                                 6, 1, 0, 3,                /* AG */
                                 0, 1, 2,                   /* AT */
                                 6, 5,                      /* CG */
                                 4};                        /* CT */

/* A-C = A-T = C-G = G-T: 1, A-G = C-T: 2, others: 0 */
/*                               AA CC GG TT AC AG AT CG CT GT          */
static int gt_sym_rate_hky4[] = {0, 0, 0, 1, 2, 1, 0, 0, 0, /* AA */
                                 0, 0, 1, 0, 0, 1, 2, 0,    /* CC */
                                 0, 0, 2, 0, 1, 0, 1,       /* GG */
                                 0, 0, 1, 0, 2, 1,          /* TT */
                                 1, 2, 2, 1, 0,             /* AC */
                                 1, 1, 0, 1,                /* AG */
                                 0, 1, 2,                   /* AT */
                                 1, 2,                      /* CG */
                                 1};                        /* CT */

/* A-C: 1, A-G: 2, A-T: 3, C-G: 4, C-T: 5, G-T: 6, others: 0 */
/*                                    AA CC GG TT AC AG AT CG CT GT CA GA TA GC
 * TC TG          */
static int gt16_sym_rate_dna4[] = {
    0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 3, 0, 0, 0, /* AA */
    0, 0, 1, 0, 0, 4, 5, 0, 1, 0, 0, 4, 5, 0,    /* CC */
    0, 0, 2, 0, 4, 0, 6, 0, 2, 0, 4, 0, 6,       /* GG */
    0, 0, 3, 0, 5, 6, 0, 0, 3, 0, 5, 6,          /* TT */
    4, 5, 2, 3, 0, 0, 0, 0, 2, 3, 0,             /* AC */
    6, 1, 0, 3, 0, 0, 0, 0, 0, 3,                /* AG */
    0, 1, 2, 0, 0, 0, 0, 0, 0,                   /* AT */
    6, 5, 2, 0, 0, 0, 0, 5,                      /* CG */
    4, 3, 0, 0, 0, 0, 0,                         /* CT */
    0, 3, 0, 5, 0, 0,                            /* GT */
    4, 5, 0, 0, 0,                               /* CA */
    6, 1, 0, 0,                                  /* GA */
    0, 1, 2,                                     /* TA */
    6, 0,                                        /* GC */
    4                                            /* TC */
};

static const corax_subst_model_t gt_model_list[] = {
    /*  name    states  model rates         model freqs   rate symmetries
       freq. sym.           */
    {"GTJC-SM", 10, gt_rates_equal_sm, gt_freqs_equal, NULL, NULL, 0},
    {"GTJC", 10, gt_rates_equal, gt_freqs_equal, NULL, NULL, 0},
    {"GTGTR-SM", 10, NULL, NULL, gt_sym_rate_free_sm, NULL, 0},
    {"GTGTR4", 10, NULL, NULL, gt_sym_rate_dna4, NULL, 0},
    {"GTHKY4", 10, NULL, NULL, gt_sym_rate_hky4, NULL, 0},
    {"GTGTR", 10, NULL, NULL, NULL, NULL, 0},
    {"GPGTR4", 16, NULL, NULL, gt16_sym_rate_dna4, NULL, 0}};

const int GT_MODELS_COUNT = sizeof(gt_model_list) / sizeof(corax_subst_model_t);

static int get_model_index(const char *model_name)
{
  int i;
  for (i = 0; i < GT_MODELS_COUNT; ++i)
    if (strcasecmp(model_name, gt_model_list[i].name) == 0) return i;

  /* model not found*/
  return -1;
}

/**
 * @brief Returns number of available built-in protein evolution models
 */
CORAX_EXPORT unsigned int corax_util_model_count_genotype()
{
  return GT_MODELS_COUNT;
}

/**
 * @brief Returns list of available built-in protein evolution models (names)
 */
CORAX_EXPORT char **corax_util_model_names_genotype()
{
  char **names = calloc(GT_MODELS_COUNT, sizeof(char *));

  int i;
  for (i = 0; i < GT_MODELS_COUNT; ++i)
  {
    const char *model_name = gt_model_list[i].name;
    names[i]               = malloc(strlen(model_name) + 1);
    strcpy(names[i], model_name);
  }

  return names;
}

/**
 * @brief Returns 1 if built-in genotype models with a given name exists and 0
 * otherwise
 */
CORAX_EXPORT int corax_util_model_exists_genotype(const char *model_name)
{
  return get_model_index(model_name) >= 0 ? 1 : 0;
}

CORAX_EXPORT int corax_util_model_exists_genotype10(const char *model_name)
{
  int idx = get_model_index(model_name);
  if (idx >= 0)
    return gt_model_list[idx].states == 10 ? 1 : 0;
  else
    return 0;
}

CORAX_EXPORT int corax_util_model_exists_genotype16(const char *model_name)
{
  int idx = get_model_index(model_name);
  if (idx >= 0)
    return gt_model_list[idx].states == 16 ? 1 : 0;
  else
    return 0;
}

/**
 * @brief Returns properties of the specified AA evolution model
 *
 * See corax_subst_model_t definition for details
 *
 * @param model_name name of the AA model
 *
 * @return model info structure, or NULL if model doesn't exist
 */
CORAX_EXPORT corax_subst_model_t *
             corax_util_model_info_genotype(const char *model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
  {
    return corax_util_model_clone(&gt_model_list[model_index]);
  }
  else
  {
    corax_set_error(CORAX_UTIL_ERROR_MODEL_UNKNOWN,
                    "Genotype model not found: %s",
                    model_name);
    return NULL;
  }
}

/**
 * @brief Set given protein model to a corax_partition_t instance
 *
 * @param partition partition instance
 * @param model_name name of the protein model
 * @param model_freqs 0: set model rate matrices only, 1: set model AA
 * frequencies as well
 *
 * @return CORAX_SUCCESS on success, CORAX_FAILURE on error (check corax_errmsg
 * for details)
 */
CORAX_EXPORT int corax_util_model_set_genotype(corax_partition_t *partition,
                                               const char *       model_name,
                                               int                model_freqs)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
  {
    corax_set_subst_params(partition, 0, gt_model_list[model_index].rates);
    if (model_freqs)
    {
      corax_set_frequencies(partition, 0, gt_model_list[model_index].freqs);
    }
    return CORAX_SUCCESS;
  }
  else
  {
    corax_set_error(CORAX_UTIL_ERROR_MODEL_UNKNOWN,
                    "Genotype model not found: %s",
                    model_name);
    return CORAX_FAILURE;
  }
}

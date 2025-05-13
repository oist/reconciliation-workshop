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

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include "corax/corax_core.h"

#ifndef CORAX_MODEL_MODUTIL_H_
#define CORAX_MODEL_MODUTIL_H_


/* error codes for UTIL coraxlib module (5001-6000)*/
#define CORAX_UTIL_ERROR_MODEL_UNKNOWN 5001
#define CORAX_UTIL_ERROR_MODEL_INVALID_DEF 5002
#define CORAX_UTIL_ERROR_MODEL_INVALID_MAPSTRING 5003
#define CORAX_UTIL_ERROR_MODEL_INVALID_MAPFILE 5004
#define CORAX_UTIL_ERROR_MIXTURE_INVALID_SIZE 5011
#define CORAX_UTIL_ERROR_MIXTURE_INVALID_COMPONENT 5012

/* rate heterogeneity mode */
#define CORAX_UTIL_MIXTYPE_FIXED (0)
#define CORAX_UTIL_MIXTYPE_GAMMA (1 << 0)
#define CORAX_UTIL_MIXTYPE_FREE (1 << 1)

/* Substitution model definition */
typedef struct subst_model
{
  const char *  name;   /* name of the model */
  unsigned int  states; /* number of states in this model */
  const double *rates;  /* model substitution rates; NULL = optimize */
  const double *freqs;  /* model base frequencies; NULL = optimize */
  const int *rate_sym;  /* substitution matrix symmetries: AC AG AT CG CT GT */
  const int *freq_sym;  /* base frequencies symmetries: A C G T */
  char       dynamic_malloc;
} corax_subst_model_t;

/* Substitution model definition */
typedef struct mixture_model
{
  char *                name;      /* name of the model                      */
  unsigned int          ncomp;     /* number of mixture components           */
  corax_subst_model_t **models;    /* list of mixture components             */
  double *              mix_rates; /* fixed mixture rates; NULL = optimize   */
  double *mix_weights;             /* fixed mixture weights; NULL = optimize */
  int     mix_type;                /* component rates: fixed, gamma or free  */
} corax_mixture_model_t;

/* Model alias name definition */
typedef struct model_alias
{
  char *alias;        /* model alias name */
  char *primary_name; /* primary name used in the model definition */
} corax_subst_model_alias_t;

#ifdef __cplusplus
extern "C"
{
#endif

/* general model management functions */
CORAX_EXPORT double *corax_util_get_equal_freqs(unsigned int states);
CORAX_EXPORT double *corax_util_get_equal_rates(unsigned int states);

CORAX_EXPORT corax_state_t *corax_util_charmap_create(unsigned int states,
                                                      const char * statechars,
                                                      const char * gapchars,
                                                      int case_sensitive);

CORAX_EXPORT corax_state_t *corax_util_charmap_parse(unsigned int states,
                                                     const char * fname,
                                                     int    case_sensitive,
                                                     char **state_names);

CORAX_EXPORT      corax_subst_model_t *
                  corax_util_model_create_custom(const char *  name,
                                                 unsigned int  states,
                                                 const double *rates,
                                                 const double *freqs,
                                                 const char *  rate_sym_str,
                                                 const char *  freq_sym_str);
CORAX_EXPORT void corax_util_model_destroy(corax_subst_model_t *model);
CORAX_EXPORT      corax_subst_model_t *
                  corax_util_model_clone(const corax_subst_model_t *src);
CORAX_EXPORT int *corax_util_model_string_to_sym(const char *s);

CORAX_EXPORT corax_mixture_model_t *
             corax_util_model_mixture_create(const char *                name,
                                             unsigned int                ncomp,
                                             corax_subst_model_t **const models,
                                             const double *              mix_rates,
                                             const double *              mix_weights,
                                             int                         mix_type);
CORAX_EXPORT void
             corax_util_model_mixture_destroy(corax_mixture_model_t *mixture);
CORAX_EXPORT corax_mixture_model_t *
             corax_util_model_mixture_clone(const corax_mixture_model_t *src);

/* functions for working with built-in DNA models */
CORAX_EXPORT unsigned int corax_util_model_count_dna();
CORAX_EXPORT char **      corax_util_model_names_dna();
CORAX_EXPORT int          corax_util_model_exists_dna(const char *model_name);
CORAX_EXPORT              corax_subst_model_t *
                          corax_util_model_info_dna(const char *model_name);

/* functions for working with built-in protein models */
CORAX_EXPORT unsigned int corax_util_model_count_protein();
CORAX_EXPORT char **      corax_util_model_names_protein();
CORAX_EXPORT int corax_util_model_exists_protein(const char *model_name);
CORAX_EXPORT     corax_subst_model_t *
                 corax_util_model_info_protein(const char *model_name);
CORAX_EXPORT int corax_util_model_set_protein(corax_partition_t *partition,
                                              const char *       model_name,
                                              int                model_freqs);

CORAX_EXPORT int corax_util_model_exists_protmix(const char *model_name);
CORAX_EXPORT     corax_mixture_model_t *
                 corax_util_model_info_protmix(const char *model_name);
CORAX_EXPORT int corax_util_model_set_protmix(corax_partition_t *partition,
                                              const char *       model_name,
                                              int                model_freqs);

/* functions for working with multistates models */
CORAX_EXPORT int corax_util_model_exists_mult(const char *model_name);
CORAX_EXPORT unsigned int
             corax_util_model_numstates_mult(const char *model_name);
CORAX_EXPORT corax_state_t *corax_util_model_charmap_mult(unsigned int states);
CORAX_EXPORT                corax_subst_model_t *
                            corax_util_model_info_mult(const char *model_name);

/* functions for working with built-in genotype models */
CORAX_EXPORT unsigned int corax_util_model_count_genotype();
CORAX_EXPORT char **      corax_util_model_names_genotype();
CORAX_EXPORT int corax_util_model_exists_genotype(const char *model_name);
CORAX_EXPORT     corax_subst_model_t *
                 corax_util_model_info_genotype(const char *model_name);

CORAX_EXPORT int corax_util_model_exists_genotype10(const char *model_name);
CORAX_EXPORT int corax_util_model_exists_genotype16(const char *model_name);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif

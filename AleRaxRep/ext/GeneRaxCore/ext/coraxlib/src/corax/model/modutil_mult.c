/*
 Copyright (C) 2017 Alexey Kozlov, Diego Darriba

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

const char mult_statechars[] =
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ!\"#$%&'()*+,/:;<=>@[\\]^_{|}~";
const char mult_gapchars[] = "-?.";

const char *skip_datatype(const char *full_model_name)
{
  const char *sep = strchr(full_model_name, '_');
  return sep ? sep + 1 : full_model_name;
}

/**
 * @brief Returns 1 if built-in MULTISTATE model with a given name exists and 0
 * otherwise
 */
CORAX_EXPORT int corax_util_model_exists_mult(const char *model_name)
{
  return strncasecmp("MULTI", model_name, 4) == 0 ? 1 : 0;
}

/**
 * @brief Parses model string (MULTIxx) and returns the number of state (xx)
 */
CORAX_EXPORT unsigned int
corax_util_model_numstates_mult(const char *model_name)
{
  unsigned int states;
  if (sscanf(model_name, "MULTI%u", &states) == 1)
    return states;
  else
    return 0;
}

/**
 * @brief Returns a character map for a MULTISTATE model with a given number of
 * states
 *
 * @param states number of states
 *
 * @return array of 256 bit-encoded state identifiers (corax_state_t) indexed by
 * ASCII code
 */
CORAX_EXPORT corax_state_t *corax_util_model_charmap_mult(unsigned int states)
{
  return corax_util_charmap_create(
      states, mult_statechars, mult_gapchars, 0 /* case_sensitive */
  );
}

/**
 * @brief Returns properties of the specified MULTISTATE evolution model
 *
 * See corax_model_t definition for details
 *
 * @param model_name name of the MULTISTATE model
 *
 * @return model info structure, or NULL if model doesn't exist
 */
CORAX_EXPORT corax_subst_model_t *
             corax_util_model_info_mult(const char *model_name)
{
  unsigned int states = corax_util_model_numstates_mult(model_name);
  if (!states)
  {
    corax_set_error(CORAX_UTIL_ERROR_MODEL_UNKNOWN,
                    "Unknown number of states in a MULTISTATE model: %s",
                    model_name);
    return NULL;
  }

  static const unsigned int maxstates = sizeof(corax_state_t) * 8;
  if (states > maxstates)
  {
    corax_set_error(
        CORAX_UTIL_ERROR_MODEL_INVALID_DEF,
        "The specified number of states (%u) exceeds the allowed maximum (%u)",
        states,
        maxstates);
    return NULL;
  }

  const char *subst_model_name = skip_datatype(model_name);
  if (strcasecmp("GTR", subst_model_name) == 0)
  {
    return corax_util_model_create_custom(
        model_name, states, NULL, NULL, NULL, NULL);
  }
  else if (strcasecmp("MK", subst_model_name) == 0
           || strcasecmp("JC", subst_model_name) == 0)
  {
    return corax_util_model_create_custom(model_name,
                                          states,
                                          corax_util_get_equal_rates(states),
                                          corax_util_get_equal_freqs(states),
                                          NULL,
                                          NULL);
  }
  else if (strncasecmp("USER", subst_model_name, 4) == 0)
  {
    return corax_util_model_create_custom(
        model_name, states, NULL, NULL, subst_model_name + 4, NULL);
  }
  else
  {
    corax_set_error(CORAX_UTIL_ERROR_MODEL_UNKNOWN,
                    "MULTISTATE model not found: %s",
                    subst_model_name);
    return NULL;
  }
}

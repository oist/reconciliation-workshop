#ifndef CORAX_UTIL_RANDOM_H_
#define CORAX_UTIL_RANDOM_H_

#include "corax/core/common.h"

/* Reentrant versions of the `random' family of functions.
   These functions all use the following data structure to contain
   state, rather than global state variables. Taken and modified from
   glibc 2.23 */

struct corax_random_data
{
  int32_t *fptr;      /* Front pointer.  */
  int32_t *rptr;      /* Rear pointer.  */
  int32_t *state;     /* Array of state values.  */
  int      rand_type; /* Type of random number generator.  */
  int      rand_deg;  /* Degree of random number generator.  */
  int      rand_sep;  /* Distance between front and rear.  */
  int32_t *end_ptr;   /* Pointer behind state table.  */
};

typedef struct corax_random_state_s
{
  struct corax_random_data rdata;
  char * state_buf; /* Buffer to store state */
} corax_random_state;


#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in random.c */

  CORAX_EXPORT extern int corax_random_r(struct corax_random_data *__buf,
                                         int32_t *                 __result);

  CORAX_EXPORT extern int corax_srandom_r(unsigned int              __seed,
                                          struct corax_random_data *__buf);

  CORAX_EXPORT extern int corax_initstate_r(unsigned int __seed,
                                            char *       __statebuf,
                                            size_t       __statelen,
                                            struct corax_random_data *__buf);

  CORAX_EXPORT extern int corax_setstate_r(char *                    __statebuf,
                                           struct corax_random_data *__buf);

  CORAX_EXPORT corax_random_state *corax_random_create(unsigned int seed);

  CORAX_EXPORT int corax_random_getint(corax_random_state *rstate, uint32_t maxval);

  CORAX_EXPORT void corax_random_destroy(corax_random_state *rstate);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_UTIL_RANDOM_H_ */

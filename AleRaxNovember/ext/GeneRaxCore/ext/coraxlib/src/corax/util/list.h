#ifndef CORAX_UTIL_LIST_H_
#define CORAX_UTIL_LIST_H_

#include "corax/core/common.h"

/* Doubly-linked list */

typedef struct corax_dlist
{
  struct corax_dlist *next;
  struct corax_dlist *prev;
  void *              data;
} corax_dlist_t;


#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in list.c */

  CORAX_EXPORT int corax_dlist_append(corax_dlist_t **dlist, void *data);
  CORAX_EXPORT int corax_dlist_remove(corax_dlist_t **dlist, void *data);
  CORAX_EXPORT int corax_dlist_prepend(corax_dlist_t **dlist, void *data);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_UTIL_LIST_H_ */

#include <stdarg.h>

#include "corax/corax.h"

__thread int  corax_errno;
__thread char corax_errmsg[200] = {0};

/**
 * @brief Set corax error (corax_errno and corax_errmsg)
 *
 * @param[in] errno the error code
 * @param[in] errmsg_fmt formatted error message
 */
__attribute__((format(printf, 2, 3))) void
corax_set_error(int _errno, const char *errmsg_fmt, ...)
{
  corax_errno = _errno;

  va_list args;
  va_start(args, errmsg_fmt);
  vsnprintf(corax_errmsg, CORAX_ERRMSG_LEN, errmsg_fmt, args);
  va_end(args);
}

/**
 * Reset corax error and messages.
 *
 * Call this function within operations whose error status depends on
 * `corax_errno` such that no error leaks in from previous operations.
 */
void corax_reset_error()
{
  corax_errno = 0;
  strcpy(corax_errmsg, "");
}

CORAX_EXPORT void *corax_aligned_alloc(size_t size, size_t alignment)
{
  void *mem;

#if (defined(__WIN32__) || defined(__WIN64__))
  mem = _aligned_malloc(size, alignment);
#else
  if (posix_memalign(&mem, alignment, size)) mem = NULL;
#endif

  return mem;
}

CORAX_EXPORT void corax_aligned_free(void *ptr)
{
#if (defined(__WIN32__) || defined(__WIN64__))
  _aligned_free(ptr);
#else
  free(ptr);
#endif
}


#pragma once
#ifndef _CORAX_ABSDIFF_H

/*
 * @brief Returns the difference of two unsigned int values without risking
 * underflow.
 * @param a Value A
 * @param b Value B
 */
#ifdef __cplusplus
extern "C"
{
#endif
  unsigned int absdiff(unsigned int a, unsigned int b);
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif

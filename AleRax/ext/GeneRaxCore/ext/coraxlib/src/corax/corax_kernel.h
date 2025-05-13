#ifndef CORAX_KERNEL_H_
#define CORAX_KERNEL_H_

#include "corax/kernel/clv.h"
#include "corax/kernel/derivatives.h"
#include "corax/kernel/likelihood.h"
#include "corax/kernel/pmatrix.h"

#include "corax/kernel/kernel_generic.h"

#ifdef HAVE_SSE3
#include "corax/kernel/kernel_sse.h"
#endif

#ifdef HAVE_AVX
#include "corax/kernel/kernel_avx.h"
#endif

#ifdef HAVE_AVX2
#include "corax/kernel/kernel_avx2.h"
#endif


#endif /* CORAX_KERNEL_H_ */

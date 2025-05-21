/*
    Copyright (C) 2017 Tomas Flouri

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

/*
    Apple machines should always default to assembly code due to
    inconsistent versioning in LLVM/clang, see issue #138

    https://github.com/xflouris/libpll/issues/138

*/

#if (defined(__APPLE__) && !defined(__aarch64__)) || \
    (!defined(__clang__) && defined(__GNUC__) && (__GNUC__ < 4 || \
      (__GNUC__ == 4 && __GNUC_MINOR__ < 8))) || \
    (defined(__clang__) && (__clang_major__ < 3 || \
      (__clang_major__ == 3 && __clang_minor__ < 9)))
  
  #if defined(__i386__) && defined(__PIC__)
    #if (defined(__GNUC__) && __GNUC__ < 3)
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("xchgl\t%%ebx, %k1\n\t"                      \
           "cpuid\n\t"                                  \
           "xchgl\t%%ebx, %k1\n\t"                      \
           : "=a" (a), "=&r" (b), "=c" (c), "=d" (d)    \
           : "0" (level), "2" (count))
    #else
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("xchg{l}\t{%%}ebx, %k1\n\t"                  \
           "cpuid\n\t"                                  \
           "xchg{l}\t{%%}ebx, %k1\n\t"                  \
           : "=a" (a), "=&r" (b), "=c" (c), "=d" (d)    \
           : "0" (level), "2" (count))
    #endif
  #elif defined(__x86_64__) && (defined(__code_model_medium__) || \
        defined(__code_model_large__)) && defined(__PIC__)
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("xchg{q}\t{%%}rbx, %q1\n\t"                  \
           "cpuid\n\t"                                  \
           "xchg{q}\t{%%}rbx, %q1\n\t"                  \
           : "=a" (a), "=&r" (b), "=c" (c), "=d" (d)    \
           : "0" (level), "2" (count))
  #else
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("cpuid\n\t"                                  \
           : "=a" (a), "=b" (b), "=c" (c), "=d" (d)     \
           : "0" (level), "2" (count))
  #endif


static void cpu_features_detect()
{
  unsigned int a, b, c, d;

  memset(&corax_hardware, 0, sizeof(corax_hardware_t));

  corax_hardware.is_initialized = 1;

#if defined(__PPC__)
  corax_hardware.altivec_present = 1;
#else

  cpuid(0, 0, a, b, c, d);
  unsigned int maxlevel = a & 0xff;

  if (maxlevel >= 1)
  {
    cpuid(1, 0, a, b, c, d);
    corax_hardware.mmx_present    = (d >> 23) & 1;
    corax_hardware.sse_present    = (d >> 25) & 1;
    corax_hardware.sse2_present   = (d >> 26) & 1;
    corax_hardware.sse3_present   = (c >> 0) & 1;
    corax_hardware.ssse3_present  = (c >> 9) & 1;
    corax_hardware.sse41_present  = (c >> 19) & 1;
    corax_hardware.sse42_present  = (c >> 20) & 1;
    corax_hardware.popcnt_present = (c >> 23) & 1;
    corax_hardware.avx_present    = (c >> 28) & 1;

    if (maxlevel >= 7)
    {
      cpuid(7, 0, a, b, c, d);
      corax_hardware.avx2_present = (b >> 5) & 1;
    }
  }
#endif
}

#else

static void cpu_features_detect()
{
  memset(&corax_hardware, 0, sizeof(corax_hardware_t));

  corax_hardware.is_initialized  = 1;
#if defined(__PPC__)
  corax_hardware.altivec_present = __builtin_cpu_supports("altivec");
#elif defined(__x86_64__) || defined(__i386__)
  corax_hardware.mmx_present    = __builtin_cpu_supports("mmx");
  corax_hardware.sse_present    = __builtin_cpu_supports("sse");
  corax_hardware.sse2_present   = __builtin_cpu_supports("sse2");
  corax_hardware.sse3_present   = __builtin_cpu_supports("sse3");
  corax_hardware.ssse3_present  = __builtin_cpu_supports("ssse3");
  corax_hardware.sse41_present  = __builtin_cpu_supports("sse4.1");
  corax_hardware.sse42_present  = __builtin_cpu_supports("sse4.2");
  corax_hardware.popcnt_present = __builtin_cpu_supports("popcnt");
  corax_hardware.avx_present    = __builtin_cpu_supports("avx");
  corax_hardware.avx2_present   = __builtin_cpu_supports("avx2");
#endif
}

#endif

static void cpu_features_show()
{
  fprintf(stderr, "CPU features:");
  if (corax_hardware.altivec_present) fprintf(stderr, " altivec");
  if (corax_hardware.mmx_present) fprintf(stderr, " mmx");
  if (corax_hardware.sse_present) fprintf(stderr, " sse");
  if (corax_hardware.sse2_present) fprintf(stderr, " sse2");
  if (corax_hardware.sse3_present) fprintf(stderr, " sse3");
  if (corax_hardware.ssse3_present) fprintf(stderr, " ssse3");
  if (corax_hardware.sse41_present) fprintf(stderr, " sse4.1");
  if (corax_hardware.sse42_present) fprintf(stderr, " sse4.2");
  if (corax_hardware.popcnt_present) fprintf(stderr, " popcnt");
  if (corax_hardware.avx_present) fprintf(stderr, " avx");
  if (corax_hardware.avx2_present) fprintf(stderr, " avx2");
  fprintf(stderr, "\n");
}

CORAX_EXPORT int corax_hardware_probe()
{
  /* probe cpu features */
  cpu_features_detect();

  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_hardware_dump()
{
  if (!corax_hardware.is_initialized) { corax_hardware_probe(); }

  cpu_features_show();
}

CORAX_EXPORT void corax_hardware_ignore()
{
  corax_hardware.is_initialized  = 1;
  corax_hardware.altivec_present = 1;
  corax_hardware.mmx_present     = 1;
  corax_hardware.sse_present     = 1;
  corax_hardware.sse2_present    = 1;
  corax_hardware.sse3_present    = 1;
  corax_hardware.ssse3_present   = 1;
  corax_hardware.sse41_present   = 1;
  corax_hardware.sse42_present   = 1;
  corax_hardware.popcnt_present  = 1;
  corax_hardware.avx_present     = 1;
  corax_hardware.avx2_present    = 1;
}

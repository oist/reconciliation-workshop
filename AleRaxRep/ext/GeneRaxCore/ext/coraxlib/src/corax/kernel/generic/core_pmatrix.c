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

#ifdef CORAX_NONREV
#include <cblas.h>
#include <lapacke.h>

int corax_core_update_pmatrix_nonrev_ss(
    double *A, size_t n, size_t lda, double t, double *P)
{
  const size_t matrix_size = n * lda;

  double *X = (double *)calloc(matrix_size, sizeof(double));
  double *N = (double *)calloc(matrix_size, sizeof(double));
  double *D = (double *)calloc(matrix_size, sizeof(double));

  /* Compute the INF norm, which we use to compute scaling factor */
  double inf_norm  = LAPACKE_dlange(CblasRowMajor, 'I', n, n, A, lda);
  int    At_norm   = (int)(inf_norm * t);
  int    scale_exp = CORAX_MIN(30, CORAX_MAX(0, 1 + At_norm));

  double Ascal = t / pow(2.0, scale_exp);

  cblas_dscal(matrix_size, Ascal, A, 1);

  /* q is a magic parameter that controls the number of iterations of the loop
   * higher is more accurate, with each increase of q decreasing error by 4
   * orders of magnitude. Anything above 12 is probably snake oil. Experiments
   * have show that 3 seems to be sufficent.
   */
  const int q    = 3;
  double    c    = 0.5;
  double    sign = -1.0;

  cblas_dcopy(matrix_size, A, 1, X, 1);

  for (size_t i = 0; i < n; ++i) { N[i * lda + i] = 1.0; }
  for (size_t i = 0; i < n; ++i) { D[i * lda + i] = 1.0; }

  /* Using fortran indexing, and we started an iteration ahead to skip some
   * setup. Furhthermore, we are going to unroll the loop to allow us to skip
   * some assignments.
   */

  cblas_daxpy(matrix_size, c, X, 1, N, 1);
  cblas_daxpy(matrix_size, sign * c, X, 1, D, 1);

  for (int i = 2; i <= q; ++i)
  {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                n,
                n,
                n,
                1.0,
                A,
                lda,
                X,
                lda,
                0.0,
                X,
                lda);
    cblas_daxpy(matrix_size, c, X, 1, N, 1);
    cblas_daxpy(matrix_size, sign * c, X, 1, D, 1);
  }

  /* Solve the equation X = N/D or DX = N */
  {
    int *ipiv = (int *)malloc(sizeof(int) * n);

    LAPACKE_dgesv(CblasRowMajor, n, n, D, lda, ipiv, N, lda);

    free(ipiv);
  }

  /*Square until we "fix" the earlier scaling */

  double *r1 = N;
  double *r2 = D;

  for (int i = 0; i < scale_exp; ++i)
  {
    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                n,
                n,
                n,
                1.0,
                r1,
                lda,
                r1,
                lda,
                0.0,
                r2,
                lda);
    double *tmp = r1;
    r1          = r2;
    r2          = tmp;
  }

  for (size_t i = 0; i < matrix_size; ++i) { P[i] = r2[i]; }

  free(X);
  free(N);
  free(D);

  return CORAX_SUCCESS;
}

int setup_ratematrix_nonrev(double *params, size_t n, size_t lda, double *rm)
{
  size_t k = 0;
  for (size_t i = 0; i < n; ++i)
  {
    double row_sum = 0.0;
    for (size_t j = 0; j < n; ++j)
    {
      if (i == j) { continue; }
      double tmp = params[k++];
      row_sum += tmp;
      rm[i * lda + j] = tmp;
    }
    rm[i * lda + i] = -row_sum;
  }
  return CORAX_SUCCESS;
}

CORAX_EXPORT int
corax_core_update_pmatrix_nonrev(double            **pmatrix,
                                 size_t              states,
                                 size_t              rate_cats,
                                 const double       *rates,
                                 const double       *branch_lengths,
                                 const unsigned int *matrix_indices,
                                 const unsigned int *params_indices,
                                 const double       *prop_invar,
                                 double *const      *params,
                                 unsigned int        count,
                                 unsigned int        attrib)
{
  double *tmp_rm = (double *)malloc(states * states * sizeof(double));
  for (size_t i = 0; i < count; ++i)
  {
    for (size_t j = 0; j < rate_cats; ++j)
    {
      double *cur_pmat   = pmatrix[matrix_indices[i]] + j * states * states;
      double *cur_params = params[params_indices[j]];
      double  cur_pinv   = prop_invar[params_indices[j]];
      double  cur_brlen  = branch_lengths[i];
      double  cur_rate   = rates[j];
      double  t          = cur_rate * cur_brlen / (1.0 - cur_pinv);
      assert(cur_pinv < 1.0);
      setup_ratematrix_nonrev(cur_params, states, states, tmp_rm);

      corax_core_update_pmatrix_nonrev_ss(tmp_rm, states, states, t, cur_pmat);
    }
  }
  free(tmp_rm);
  return CORAX_SUCCESS;
}
#endif // CORAX_NONREV

CORAX_EXPORT int corax_core_update_pmatrix(double            **pmatrix,
                                           unsigned int        states,
                                           unsigned int        rate_cats,
                                           const double       *rates,
                                           const double       *branch_lengths,
                                           const unsigned int *matrix_indices,
                                           const unsigned int *params_indices,
                                           const double       *prop_invar,
                                           double *const      *eigenvals,
                                           double *const      *eigenvecs,
                                           double *const      *inv_eigenvecs,
                                           unsigned int        count,
                                           unsigned int        attrib)
{
  unsigned int i, n, j, k, m;
  unsigned int states_padded = states;
  double      *expd;
  double      *temp;

  double  pinvar;
  double *evecs;
  double *inv_evecs;
  double *evals;
  double *pmat;

#ifdef HAVE_SSE3
  if (attrib & CORAX_ATTRIB_ARCH_SSE && CORAX_HAS_CPU_FEATURE(sse3_present))
  {
    if (states == 4)
    {
      return corax_core_update_pmatrix_4x4_sse(pmatrix,
                                               rate_cats,
                                               rates,
                                               branch_lengths,
                                               matrix_indices,
                                               params_indices,
                                               prop_invar,
                                               eigenvals,
                                               eigenvecs,
                                               inv_eigenvecs,
                                               count);
    }
    else if (states == 20)
    {
      return corax_core_update_pmatrix_20x20_sse(pmatrix,
                                                 rate_cats,
                                                 rates,
                                                 branch_lengths,
                                                 matrix_indices,
                                                 params_indices,
                                                 prop_invar,
                                                 eigenvals,
                                                 eigenvecs,
                                                 inv_eigenvecs,
                                                 count);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states + 1) & 0xFFFFFFFE;
  }
#endif
#ifdef HAVE_AVX
  if (attrib & CORAX_ATTRIB_ARCH_AVX && CORAX_HAS_CPU_FEATURE(avx_present))
  {
    if (states == 4)
    {
      return corax_core_update_pmatrix_4x4_avx(pmatrix,
                                               rate_cats,
                                               rates,
                                               branch_lengths,
                                               matrix_indices,
                                               params_indices,
                                               prop_invar,
                                               eigenvals,
                                               eigenvecs,
                                               inv_eigenvecs,
                                               count);
    }
    if (states == 20)
    {
      return corax_core_update_pmatrix_20x20_avx(pmatrix,
                                                 rate_cats,
                                                 rates,
                                                 branch_lengths,
                                                 matrix_indices,
                                                 params_indices,
                                                 prop_invar,
                                                 eigenvals,
                                                 eigenvecs,
                                                 inv_eigenvecs,
                                                 count);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states + 3) & 0xFFFFFFFC;
  }
#endif
#ifdef HAVE_AVX2
  if (attrib & CORAX_ATTRIB_ARCH_AVX2 && CORAX_HAS_CPU_FEATURE(avx2_present))
  {
    if (states == 4)
    {
      /* use AVX version here since FMA doesn't make much sense */
      return corax_core_update_pmatrix_4x4_avx(pmatrix,
                                               rate_cats,
                                               rates,
                                               branch_lengths,
                                               matrix_indices,
                                               params_indices,
                                               prop_invar,
                                               eigenvals,
                                               eigenvecs,
                                               inv_eigenvecs,
                                               count);
    }
    else if (states == 20)
    {
      return corax_core_update_pmatrix_20x20_avx2(pmatrix,
                                                  rate_cats,
                                                  rates,
                                                  branch_lengths,
                                                  matrix_indices,
                                                  params_indices,
                                                  prop_invar,
                                                  eigenvals,
                                                  eigenvecs,
                                                  inv_eigenvecs,
                                                  count);
    }
    else
    {
      return corax_core_update_pmatrix_avx2(pmatrix,
                                            states,
                                            rate_cats,
                                            rates,
                                            branch_lengths,
                                            matrix_indices,
                                            params_indices,
                                            prop_invar,
                                            eigenvals,
                                            eigenvecs,
                                            inv_eigenvecs,
                                            count);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states + 3) & 0xFFFFFFFC;
  }
#endif

  expd = (double *)malloc(states * sizeof(double));
  temp = (double *)malloc(states * states * sizeof(double));

  if (!expd || !temp)
  {
    if (expd) free(expd);
    if (temp) free(temp);

    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pmat = pmatrix[matrix_indices[i]] + n * states * states_padded;

      pinvar    = prop_invar[params_indices[n]];
      evecs     = eigenvecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals     = eigenvals[params_indices[n]];

      if (branch_lengths[i] > 0.)
      {
        /* NOTE: in order to deal with numerical issues in cases when Qt -> 0,
         * we use a trick suggested by Ben Redelings and explained here:
         * https://github.com/xflouris/libpll/issues/129#issuecomment-304004005
         * In short, we use expm1() to compute (exp(Qt) - I), and then correct
         * for this by adding an identity matrix I in the very end */

        /* exponentiate eigenvalues */
        if (pinvar > CORAX_MISC_EPSILON)
        {
          for (j = 0; j < states; ++j)
            expd[j] =
                expm1(evals[j] * rates[n] * branch_lengths[i] / (1.0 - pinvar));
        }
        else
        {
          for (j = 0; j < states; ++j)
            expd[j] = expm1(evals[j] * rates[n] * branch_lengths[i]);
        }

        for (j = 0; j < states; ++j)
          for (k = 0; k < states; ++k)
            temp[j * states + k] = inv_evecs[j * states_padded + k] * expd[k];

        for (j = 0; j < states; ++j)
        {
          for (k = 0; k < states; ++k)
          {
            pmat[j * states_padded + k] = (j == k) ? 1.0 : 0;
            for (m = 0; m < states; ++m)
            {
              pmat[j * states_padded + k] +=
                  temp[j * states + m] * evecs[m * states_padded + k];
            }
          }
        }
      }
      else
      {
        /* if branch length is zero then set the p-matrix to identity matrix
         */
        for (j = 0; j < states; ++j)
          for (k = 0; k < states; ++k)
            pmat[j * states_padded + k] = (j == k) ? 1 : 0;
      }

#ifdef DEBUG
      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k) assert(pmat[j * states_padded + k] >= 0);
#endif
    }
  }

  free(expd);
  free(temp);
  return CORAX_SUCCESS;
}

#include "corax/corax.h"

static int mytqli(double *d, double *e, const unsigned int n, double **z)
{
  unsigned int m, l, iter, i, k;
  double       s, r, p, g, f, dd, c, b;

  for (i = 2; i <= n; i++) e[i - 2] = e[i - 1];

  e[n - 1] = 0.0;

  for (l = 1; l <= n; l++)
  {
    iter = 0;
    do {
      for (m = l; m <= n - 1; m++)
      {
        dd = fabs(d[m - 1]) + fabs(d[m]);
        if (fabs(e[m - 1]) + dd == dd) break;
      }
      if (m != l)
      {
        assert(iter < 30);

        g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);
        r = sqrt((g * g) + 1.0);
        g = d[m - 1] - d[l - 1]
            + e[l - 1]
                  / (g + ((g < 0) ? -fabs(r) : fabs(r))); /*MYSIGN(r, g));*/
        s = c = 1.0;
        p     = 0.0;

        for (i = m - 1; i >= l; i--)
        {
          f = s * e[i - 1];
          b = c * e[i - 1];
          if (fabs(f) >= fabs(g))
          {
            c    = g / f;
            r    = sqrt((c * c) + 1.0);
            e[i] = f * r;
            c *= (s = 1.0 / r);
          }
          else
          {
            s    = f / g;
            r    = sqrt((s * s) + 1.0);
            e[i] = g * r;
            s *= (c = 1.0 / r);
          }
          g    = d[i] - p;
          r    = (d[i - 1] - g) * s + 2.0 * c * b;
          p    = s * r;
          d[i] = g + p;
          g    = c * r - b;
          for (k = 1; k <= n; k++)
          {
            f               = z[i][k - 1];
            z[i][k - 1]     = s * z[i - 1][k - 1] + c * f;
            z[i - 1][k - 1] = c * z[i - 1][k - 1] - s * f;
          }
        }

        d[l - 1] = d[l - 1] - p;
        e[l - 1] = g;
        e[m - 1] = 0.0;
      }
    } while (m != l);
  }

  return (1);
}

static void mytred2(double **a, const unsigned int n, double *d, double *e)
{
  unsigned int l, k, j, i;
  double       scale, hh, h, g, f;

  for (i = n; i > 1; i--)
  {
    l     = i - 1;
    h     = 0.0;
    scale = 0.0;

    if (l > 1)
    {
      for (k = 1; k <= l; k++) scale += fabs(a[k - 1][i - 1]);
      if (scale == 0.0)
        e[i - 1] = a[l - 1][i - 1];
      else
      {
        for (k = 1; k <= l; k++)
        {
          a[k - 1][i - 1] /= scale;
          h += a[k - 1][i - 1] * a[k - 1][i - 1];
        }
        f        = a[l - 1][i - 1];
        g        = ((f > 0) ? -sqrt(h) : sqrt(h)); /* diff */
        e[i - 1] = scale * g;
        h -= f * g;
        a[l - 1][i - 1] = f - g;
        f               = 0.0;
        for (j = 1; j <= l; j++)
        {
          a[i - 1][j - 1] = a[j - 1][i - 1] / h;
          g               = 0.0;
          for (k = 1; k <= j; k++) g += a[k - 1][j - 1] * a[k - 1][i - 1];
          for (k = j + 1; k <= l; k++) g += a[j - 1][k - 1] * a[k - 1][i - 1];
          e[j - 1] = g / h;
          f += e[j - 1] * a[j - 1][i - 1];
        }
        hh = f / (h + h);
        for (j = 1; j <= l; j++)
        {
          f        = a[j - 1][i - 1];
          g        = e[j - 1] - hh * f;
          e[j - 1] = g;
          for (k = 1; k <= j; k++)
            a[k - 1][j - 1] -= (f * e[k - 1] + g * a[k - 1][i - 1]);
        }
      }
    }
    else
      e[i - 1] = a[l - 1][i - 1];
    d[i - 1] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;

  for (i = 1; i <= n; i++)
  {
    l = i - 1;
    if (d[i - 1] != 0.0)
    {
      for (j = 1; j <= l; j++)
      {
        g = 0.0;
        for (k = 1; k <= l; k++) g += a[k - 1][i - 1] * a[j - 1][k - 1];
        for (k = 1; k <= l; k++) a[j - 1][k - 1] -= g * a[i - 1][k - 1];
      }
    }
    d[i - 1]        = a[i - 1][i - 1];
    a[i - 1][i - 1] = 1.0;
    for (j = 1; j <= l; j++) a[i - 1][j - 1] = a[j - 1][i - 1] = 0.0;
  }
}

/* TODO: Add code for SSE/AVX. Perhaps allocate qmatrix in one chunk to avoid
the complex checking when to dealloc */
static double **create_ratematrix(const double *params,
                                  const double *freqs,
                                  unsigned int  states)
{
  unsigned int i, j, k, success;

  double **qmatrix;

  /* normalize substitution parameters */
  unsigned int params_count = CORAX_SUBST_RATE_COUNT(states);
  double *params_normalized = (double *)malloc(sizeof(double) * params_count);
  if (!params_normalized) return NULL;

  memcpy(params_normalized, params, params_count * sizeof(double));

  if (params_normalized[params_count - 1] > 0.0)
  {
    for (i = 0; i < params_count; ++i)
      params_normalized[i] /= params_normalized[params_count - 1];
  }

  /* allocate qmatrix */
  qmatrix = (double **)malloc(states * sizeof(double *));
  if (!qmatrix)
  {
    free(params_normalized);
    return NULL;
  }

  success = 1;
  for (i = 0; i < states; ++i)
    if (!(qmatrix[i] = (double *)malloc(states * sizeof(double)))) success = 0;

  if (!success)
  {
    for (i = 0; i < states; ++i) free(qmatrix[i]);
    free(qmatrix);
    free(params_normalized);
    return NULL;
  }

  /* construct a matrix equal to sqrt(pi) * Q sqrt(pi)^-1 in order to ensure
     it is symmetric */

  for (i = 0; i < states; ++i) qmatrix[i][i] = 0;

  k = 0;
  for (i = 0; i < states; ++i)
  {
    for (j = i + 1; j < states; ++j)
    {
      double factor =
          (freqs[i] <= CORAX_EIGEN_MINFREQ || freqs[j] <= CORAX_EIGEN_MINFREQ)
              ? 0
              : params_normalized[k];
      k++;
      qmatrix[i][j] = qmatrix[j][i] = factor * sqrt(freqs[i] * freqs[j]);
      qmatrix[i][i] -= factor * freqs[j];
      qmatrix[j][j] -= factor * freqs[i];
    }
  }

  double mean = 0;
  for (i = 0; i < states; ++i) mean += freqs[i] * (-qmatrix[i][i]);
  for (i = 0; i < states; ++i)
  {
    for (j = 0; j < states; ++j) qmatrix[i][j] /= mean;
  }

  free(params_normalized);

  return qmatrix;
}

static unsigned int eliminate_zero_states(double **     mat,
                                          const double *forg,
                                          unsigned int  states,
                                          double *      new_forg)
{
  unsigned int i, j, inew, jnew;
  unsigned int new_states = 0;
  for (i = 0; i < states; i++)
  {
    if (forg[i] > CORAX_EIGEN_MINFREQ) new_forg[new_states++] = forg[i];
  }

  assert(new_states <= states);

  if (new_states < states)
  {
    for (i = 0, inew = 0; i < states; i++)
    {
      if (forg[i] > CORAX_EIGEN_MINFREQ)
      {
        for (j = 0, jnew = 0; j < states; j++)
        {
          if (forg[j] > CORAX_EIGEN_MINFREQ)
          {
            mat[inew][jnew] = mat[i][j];
            jnew++;
          }
        }
        inew++;
      }
    }
  }

  return new_states;
}

CORAX_EXPORT int corax_update_eigen(corax_partition_t *partition,
                                    unsigned int       params_index)
{
  unsigned int i, j;
  double *     e, *d;
  double **    a;

  double *eigenvecs     = partition->eigenvecs[params_index];
  double *inv_eigenvecs = partition->inv_eigenvecs[params_index];
  double *eigenvals     = partition->eigenvals[params_index];
  double *freqs         = partition->frequencies[params_index];
  double *subst_params  = partition->subst_params[params_index];

  unsigned int states        = partition->states;
  unsigned int states_padded = partition->states_padded;

  unsigned int inew, jnew;
  unsigned int new_states;
  double *     new_freqs = NULL;

  a = create_ratematrix(subst_params, freqs, states);
  if (!a)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  d         = (double *)malloc(states * sizeof(double));
  e         = (double *)malloc(states * sizeof(double));
  new_freqs = (double *)malloc(states * sizeof(double));
  if (!d || !e || !new_freqs)
  {
    if (d) free(d);
    if (e) free(e);
    if (new_freqs) free(new_freqs);
    for (i = 0; i < states; ++i) free(a[i]);
    free(a);
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  /* Here we use a technical trick to reduce rate matrix if some states
   * have (near) zero frequencies. Code adapted from IQTree, see:
   * https://github.com/Cibiv/IQ-TREE/commit/f222d317af46cf6abf8bcdb70d4db22475e9a7d2
   */
  new_states = eliminate_zero_states(a, freqs, states, new_freqs);

  mytred2(a, new_states, d, e);
  mytqli(d, e, new_states, a);

  for (i = 0, inew = 0; i < states; i++)
    eigenvals[i] = (freqs[i] > CORAX_EIGEN_MINFREQ) ? d[inew++] : 0;

  assert(inew == new_states);

  /* pre-compute square roots of frequencies */
  for (i = 0; i < new_states; i++) new_freqs[i] = sqrt(new_freqs[i]);

  if (new_states < states)
  {
    /* initialize eigenvecs and inv_eigenvecs with diagonal matrix */
    memset(eigenvecs, 0, states_padded * states * sizeof(double));
    memset(inv_eigenvecs, 0, states_padded * states * sizeof(double));

    for (i = 0; i < states; i++)
    {
      eigenvecs[i * states_padded + i]     = 1.;
      inv_eigenvecs[i * states_padded + i] = 1.;
    }

    for (i = 0, inew = 0; i < states; i++)
    {
      if (freqs[i] > CORAX_EIGEN_MINFREQ)
      {
        for (j = 0, jnew = 0; j < states; j++)
        {
          if (freqs[j] > CORAX_EIGEN_MINFREQ)
          {
            /* multiply the eigen vectors from the right with sqrt(pi) */
            eigenvecs[i * states_padded + j] = a[inew][jnew] * new_freqs[jnew];
            /* multiply the inverse eigen vectors from the left with sqrt(pi)^-1
             */
            inv_eigenvecs[i * states_padded + j] =
                a[jnew][inew] / new_freqs[inew];
            jnew++;
          }
        }
        inew++;
      }
    }
  }
  else
  {
    for (i = 0; i < states; i++)
    {
      for (j = 0; j < states; j++)
      {
        /* multiply the eigen vectors from the right with sqrt(pi) */
        eigenvecs[i * states_padded + j] = a[i][j] * new_freqs[j];
        /* multiply the inverse eigen vectors from the left with sqrt(pi)^-1 */
        inv_eigenvecs[i * states_padded + j] = a[j][i] / new_freqs[i];
      }
    }
  }

  partition->eigen_decomp_valid[params_index] = 1;

  free(d);
  free(e);
  free(new_freqs);
  for (i = 0; i < states; ++i) free(a[i]);
  free(a);

  return CORAX_SUCCESS;
}


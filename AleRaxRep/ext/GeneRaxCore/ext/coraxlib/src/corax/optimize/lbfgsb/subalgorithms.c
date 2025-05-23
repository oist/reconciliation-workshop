/*
 * L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
 * or "3-clause license")
 * Please read attached file License.txt
 *
 * ===========   L-BFGS-B (version 3.0.  April 25, 2011  ===================
 *
 *    This is a modified version of L-BFGS-B. Minor changes in the updated
 *    code appear preceded by a line comment as follows
 *
 *    c-jlm-jn
 *
 *    Major changes are described in the accompanying paper:
 *
 *        Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
 *        L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
 *        Optimization"  (2011). To appear in  ACM Transactions on
 *        Mathematical Software,
 */
#include "lbfgsb.h"

static int c__1 = 1;
static int c__11 = 11;

int bmv (int *m, double *sy, double *wt, int *col, double *v, double *p,
         int *info);
int hpsolb (int *n, double *t, int *iorder, int *iheap);

/*     This subroutine initializes iwhere and projects the initial x to */
/*       the feasible set if necessary. */
/*     iwhere is an int array of dimension n. */
/*       On entry iwhere is unspecified. */
/*       On exit iwhere(i)=-1  if x(i) has no bounds */
/*                         3   if l(i)=u(i) */
/*                         0   otherwise. */
/*       In cauchy, iwhere is given finer gradations. */
int active (int *n, double *l, double *u, int *nbd, double *x, int *iwhere,
            int *iprint, logical *prjctd, logical *cnstnd, logical *boxed)
{
  CORAX_UNUSED(iprint);
  int i, nbdd;

  nbdd = 0;
  *prjctd = FALSE_;
  *cnstnd = FALSE_;
  *boxed = TRUE_;

  /*     Project the initial x to the easible set if necessary. */
  for (i = 0; i < *n; ++i)
  {
    if (nbd[i] > 0)
    {
      if (nbd[i] <= 2 && x[i] <= l[i])
      {
        if (x[i] < l[i])
        {
          *prjctd = TRUE_;
          x[i] = l[i];
        }
        ++nbdd;
      }
      else if (nbd[i] >= 2 && x[i] >= u[i])
      {
        if (x[i] > u[i])
        {
          *prjctd = TRUE_;
          x[i] = u[i];
        }
        ++nbdd;
      }
    }
  }
  /*     Initialize iwhere and assign values to cnstnd and boxed. */
  for (i = 0; i < *n; ++i)
  {
    if (nbd[i] != 2)
    {
      *boxed = FALSE_;
    }
    if (nbd[i] == 0)
    {
      /*                                this variable is always free */
      iwhere[i] = -1;
      /*           otherwise set x(i)=mid(x(i), u(i), l(i)). */
    }
    else
    {
      *cnstnd = TRUE_;
      if (nbd[i] == 2 && u[i] - l[i] <= 0.)
      {
        /*                   this variable is always fixed */
        iwhere[i] = 3;
      }
      else
      {
        iwhere[i] = 0;
      }
    }
  }

#ifdef DEBUG
  if (*prjctd)
  {
    DBG("[L-BFGS-B] The initial X is infeasible. Restart with its projection\n");
  }
  if (!(*cnstnd))
  {
    DBG("[L-BFGS-B] This problem is unconstrained\n");
  } DBG("[L-BFGS-B] At X0, %d variables are exactly at the bounds\n", nbdd);
#endif

  return 0;
} /* active */

/*     This subroutine computes the product of the 2m x 2m middle matrix */
/*       in the compact L-BFGS formula of B and a 2m vector v; */
/*       it returns the product in p. */

/*     m is an int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     sy is a double precision array of dimension m x m. */
/*       On entry sy specifies the matrix S'Y. */
/*       On exit sy is unchanged. */

/*     wt is a double precision array of dimension m x m. */
/*       On entry wt specifies the upper triangular matrix J' which is */
/*         the Cholesky factor of (thetaS'S+LD^(-1)L'). */
/*       On exit wt is unchanged. */

/*     col is an int variable. */
/*       On entry col specifies the number of s-vectors (or y-vectors) */
/*         stored in the compact L-BFGS formula. */
/*       On exit col is unchanged. */

/*     v is a double precision array of dimension 2col. */
/*       On entry v specifies vector v. */
/*       On exit v is unchanged. */

/*     p is a double precision array of dimension 2col. */
/*       On entry p is unspecified. */
/*       On exit p is the product Mv. */

/*     info is an int variable. */
/*       On entry info is unspecified. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return when the system */
/*                                to be solved by dtrsl is singular. */
int bmv (int *m, double *sy, double *wt, int *col, double *v, double *p,
         int *info)
{
  /* System generated locals */
  int sy_dim1, sy_offset, wt_dim1, wt_offset;

  /* Local variables */
  int i, k, i2;
  double sum;

  /* Parameter adjustments */
  wt_dim1 = *m;
  wt_offset = 1 + wt_dim1;
  wt -= wt_offset;
  sy_dim1 = *m;
  sy_offset = 1 + sy_dim1;
  sy -= sy_offset;
  --p;
  --v;

  /* Function Body */
  if (*col == 0)
  {
    return 0;
  }
  /*     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ] */
  /*                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ]. */
  /*       solve Jp2=v2+LD^(-1)v1. */
  p[*col + 1] = v[*col + 1];
  for (i = 2; i <= *col; ++i)
  {
    i2 = *col + i;
    sum = 0.;
    for (k = 1; k <= i - 1; ++k)
    {
      sum += sy[i + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
    }
    p[i2] = v[i2] + sum;
  }
  /*     Solve the triangular system */
  dtrsl (&wt[wt_offset], m, col, &p[*col + 1], &c__11, info);
  if (*info != 0)
  {
    return 0;
  }
  /*       solve D^(1/2)p1=v1. */
  for (i = 1; i <= *col; ++i)
  {
    p[i] = v[i] / sqrt (sy[i + i * sy_dim1]);
  }
  /*     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ] */
  /*                    [  0         J'           ] [ p2 ]   [ p2 ]. */
  /*       solve J^Tp2=p2. */
  dtrsl (&wt[wt_offset], m, col, &p[*col + 1], &c__1, info);
  if (*info != 0)
  {
    return 0;
  }
  for (i = 1; i <= *col; ++i)
  {
    p[i] = -p[i] / sqrt (sy[i + i * sy_dim1]);
  }
  for (i = 1; i <= *col; ++i)
  {
    sum = 0.;
    for (k = i + 1; k <= *col; ++k)
    {
      sum += sy[k + i * sy_dim1] * p[*col + k] / sy[i + i * sy_dim1];
    }
    p[i] += sum;
  }
  return 0;
} /* bmv */


/*     For given x, l, u, g (with sbgnrm > 0), and a limited memory */
/*       BFGS matrix B defined in terms of matrices WY, WS, WT, and */
/*       scalars head, col, and theta, this subroutine computes the */
/*       generalized Cauchy point (GCP), defined as the first local */
/*       minimizer of the quadratic */

/*                  Q(x + s) = g's + 1/2 s'Bs */

/*       along the projected gradient direction P(x-tg,l,u). */
/*       The routine returns the GCP in xcp. */
int cauchy(int *n, double *x, double *l,
	double *u, int *nbd, double *g, int *iorder, int *
	iwhere, double *t, double *d__, double *xcp, int *m,
	double *wy, double *ws, double *sy, double *wt,
	double *theta, int *col, int *head, double *p,
	double *c__, double *wbp, double *v, int *nseg,
	int *iprint, double *sbgnrm, int *info, double *
	epsmch)
{
    CORAX_UNUSED(iprint);

    /* System generated locals */
    int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset,
	    wt_dim1, wt_offset, i__1, i__2;
    double d__1;


    /* Local variables */
    int i, j;
    double f1, f2, dt, tj, tl, tu, tj0;
    int ibp;
    double dtm;
    double wmc, wmp, wmw;
    int col2;
    double dibp;
    int iter;
    double zibp, tsum, dibp2;
    logical bnded;
    double neggi;
    int nfree;
    double bkmin;
    int nleft;
    double f2_org__;
    int nbreak, ibkmin;
    int pointr;
    logical xlower, xupper;

    tu = tl = 0.0;

    /* Parameter adjustments */
    --xcp;
    --d__;
    --t;
    --iwhere;
    --iorder;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --v;
    --wbp;
    --c__;
    --p;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;

    /* Function Body */
    if (*sbgnrm <= 0.) {
        DBG("[L-BFGS-B] Subnorm = 0. GCP = X.\n");
        dcopy(n, &x[1], &c__1, &xcp[1], &c__1);
        return 0;
    }
    bnded = TRUE_;
    nfree = *n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = *col << 1;
    f1 = 0.;
    DBG("[L-BFGS-B] CAUCHY entered\n");

    /*     We set p to zero and build it up as we determine d. */
    i__1 = col2;
    for (i = 1; i <= i__1; ++i) {
        p[i] = 0.;
        /* L20: */
    }
    /*     In the following loop we determine for each variable its bound */
    /*        status and its breakpoint, and update p accordingly. */
    /*        Smallest breakpoint is identified. */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
        neggi = -g[i];
        if (iwhere[i] != 3 && iwhere[i] != -1) {
            /*             if x(i) is not a constant and has bounds, */
            /*             compute the difference between x(i) and its bounds. */
            if (nbd[i] <= 2) {
                tl = x[i] - l[i];
            }
            if (nbd[i] >= 2) {
                tu = u[i] - x[i];
            }
            /*           If a variable is close enough to a bound */
            /*             we treat it as at bound. */
            xlower = nbd[i] <= 2 && tl <= 0.;
            xupper = nbd[i] >= 2 && tu <= 0.;
            /*              reset iwhere(i). */
            iwhere[i] = 0;
            if (xlower) {
                if (neggi <= 0.) {
                    iwhere[i] = 1;
                }
            } else if (xupper) {
                if (neggi >= 0.) {
                    iwhere[i] = 2;
                }
            } else {
                if (abs(neggi) <= 0.) {
                    iwhere[i] = -3;
                }
            }
        }
        pointr = *head;
        if (iwhere[i] != 0 && iwhere[i] != -1) {
            d__[i] = 0.;
        } else {
            d__[i] = neggi;
            f1 -= neggi * neggi;
            /*             calculate p := p - W'e_i* (g_i). */
            i__2 = *col;
            for (j = 1; j <= i__2; ++j) {
                p[j] += wy[i + pointr * wy_dim1] * neggi;
                p[*col + j] += ws[i + pointr * ws_dim1] * neggi;
                pointr = pointr % *m + 1;
                /* L40: */
            }
            if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.) {
                /*                                 x(i) + d(i) is bounded; compute t(i). */
                ++nbreak;
                iorder[nbreak] = i;
                t[nbreak] = tl / (-neggi);
                if (nbreak == 1 || t[nbreak] < bkmin) {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
            } else if (nbd[i] >= 2 && neggi > 0.) {
                /*                                 x(i) + d(i) is bounded; compute t(i). */
                ++nbreak;
                iorder[nbreak] = i;
                t[nbreak] = tu / neggi;
                if (nbreak == 1 || t[nbreak] < bkmin) {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
            } else {
                /*                x(i) + d(i) is not bounded. */
                --nfree;
                iorder[nfree] = i;
                if (abs(neggi) > 0.) {
                    bnded = FALSE_;
                }
            }
        }
        /* L50: */
    }
    /*     The indices of the nonzero components of d are now stored */
    /*       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
    /*       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
    if (*theta != 1.) {
        /*                   complete the initialization of p for theta not= one. */
        dscal(col, theta, &p[*col + 1], &c__1);
    }
    /*     Initialize GCP xcp = x. */
    dcopy(n, &x[1], &c__1, &xcp[1], &c__1);
    if (nbreak == 0 && nfree == *n + 1) {
        /*                  is a zero vector, return with the initial xcp as GCP. */
#ifdef DEBUG
            DBG("[L-BFGS-B] Cauchy X = ");
            i__1 = *n;
            for (i = 1; i <= i__1; ++i) {
                DBG("%5.2e ", xcp[i] );
            }
            DBG("\n");
#endif
        return 0;
    }
    /*     Initialize c = W'(xcp - x) = 0. */
    i__1 = col2;
    for (j = 1; j <= i__1; ++j) {
        c__[j] = 0.;
        /* L60: */
    }
    /*     Initialize derivative f2. */
    f2 = -(*theta) * f1;
    f2_org__ = f2;
    if (*col > 0) {
        bmv(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
        if (*info != 0) {
            return 0;
        }
        f2 -= ddot(&col2, &v[1], &c__1, &p[1], &c__1);
    }
    dtm = -f1 / f2;
    tsum = 0.;
    *nseg = 1;
    DBG("[L-BFGS-B] There are %d breakpoints\n", nbreak );

    /*     If there are no breakpoints, locate the GCP and return. */
    if (nbreak == 0) {
        goto L888;
    }
    nleft = nbreak;
    iter = 1;
    tj = 0.;
    /* ------------------- the beginning of the loop ------------------------- */
L777:
    /*     Find the next smallest breakpoint; */
    /*       compute dt = t(nleft) - t(nleft + 1). */
    tj0 = tj;
    if (iter == 1) {
        /*         Since we already have the smallest breakpoint we need not do */
        /*         heapsort yet. Often only one breakpoint is used and the */
        /*         cost of heapsort is avoided. */
        tj = bkmin;
        ibp = iorder[ibkmin];
    } else {
        if (iter == 2) {
            /*             Replace the already used smallest breakpoint with the */
            /*             breakpoint numbered nbreak > nlast, before heapsort call. */
            if (ibkmin != nbreak) {
                t[ibkmin] = t[nbreak];
                iorder[ibkmin] = iorder[nbreak];
            }
            /*        Update heap structure of breakpoints */
            /*           (if iter=2, initialize heap). */
        }
        i__1 = iter - 2;
        hpsolb(&nleft, &t[1], &iorder[1], &i__1);
        tj = t[nleft];
        ibp = iorder[nleft];
    }
    dt = tj - tj0;
#ifdef DEBUG
    if (dt != 0.) {
        DBG("[L-BFGS-B] Piece %d --f1, f2 at start point %.2e %.2e\n", *nseg, f1, f2 );
        DBG("[L-BFGS-B] Distance to the next break point = %.2e\n", dt );
        DBG("[L-BFGS-B] Distance to the stationary point = %.2e\n", dtm );
    }
#endif
    /*     If a minimizer is within this interval, locate the GCP and return. */
    if (dtm < dt) {
        goto L888;
    }
    /*     Otherwise fix one variable and */
    /*       reset the corresponding component of d to zero. */
    tsum += dt;
    --nleft;
    ++iter;
    dibp = d__[ibp];
    d__[ibp] = 0.;
    if (dibp > 0.) {
        zibp = u[ibp] - x[ibp];
        xcp[ibp] = u[ibp];
        iwhere[ibp] = 2;
    } else {
        zibp = l[ibp] - x[ibp];
        xcp[ibp] = l[ibp];
        iwhere[ibp] = 1;
    }
    DBG("[L-BFGS-B] Variable %d is fixed\n", ibp );

    if (nleft == 0 && nbreak == *n) {
        /*                                             all n variables are fixed, */
        /*                                                return with xcp as GCP. */
        dtm = dt;
        goto L999;
    }
    /*     Update the derivative information. */
    ++(*nseg);
    /* Computing 2nd power */
    d__1 = dibp;
    dibp2 = d__1 * d__1;
    /*     Update f1 and f2. */
    /*        temporarily set f1 and f2 for col=0. */
    f1 = f1 + dt * f2 + dibp2 - *theta * dibp * zibp;
    f2 -= *theta * dibp2;
    if (*col > 0) {
        /*                          update c = c + dt*p. */
        daxpy(&col2, &dt, &p[1], &c__1, &c__[1], &c__1);
        /*           choose wbp, */
        /*           the row of W corresponding to the breakpoint encountered. */
        pointr = *head;
        i__1 = *col;
        for (j = 1; j <= i__1; ++j) {
            wbp[j] = wy[ibp + pointr * wy_dim1];
            wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
            pointr = pointr % *m + 1;
            /* L70: */
        }
        /*           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
        bmv(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
        if (*info != 0) {
            return 0;
        }
        wmc = ddot(&col2, &c__[1], &c__1, &v[1], &c__1);
        wmp = ddot(&col2, &p[1], &c__1, &v[1], &c__1);
        wmw = ddot(&col2, &wbp[1], &c__1, &v[1], &c__1);
        /*           update p = p - dibp*wbp. */
        d__1 = -dibp;
        daxpy(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
        /*           complete updating f1 and f2 while col > 0. */
        f1 += dibp * wmc;
        f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }
    /* Computing MAX */
    d__1 = *epsmch * f2_org__;
    f2 = max(d__1,f2);
    if (nleft > 0) {
        dtm = -f1 / f2;
        goto L777;
        /*                 to repeat the loop for unsearched intervals. */
    } else if (bnded) {
        f1 = 0.;
        f2 = 0.;
        dtm = 0.;
    } else {
        dtm = -f1 / f2;
    }
    /* ------------------- the end of the loop ------------------------------- */
L888:
    DBG("[L-BFGS-B] \nGCP found in this segment. Piece %d --f1, f2 at start point %.2e %.2e\n", *nseg, f1, f2 );
    DBG("[L-BFGS-B] Distance to the stationary point = %.2e\n", dtm );
    if (dtm <= 0.) {
        dtm = 0.;
    }
    tsum += dtm;
    /*     Move free variables (i.e., the ones w/o breakpoints) and */
    /*       the variables whose breakpoints haven't been reached. */
    daxpy(n, &tsum, &d__[1], &c__1, &xcp[1], &c__1);
L999:
    /*     Update c = c + dtm*p = W'(x^c - x) */
    /*       which will be used in computing r = Z'(B(x^c - x) + g). */
    if (*col > 0) {
        daxpy(&col2, &dtm, &p[1], &c__1, &c__[1], &c__1);
    }
#ifdef DEBUG
  DBG("[L-BFGS-B] Cauchy X = ");
  i__1 = *n;
  for (i = 1; i <= i__1; ++i)
  {
    DBG("%5.2e ", xcp[i] );
  }
  DBG("\n[L-BFGS-B] -------------- exit CAUCHY -----------\n");
#endif
    return 0;
} /* cauchy */


/*       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using */
/*         wa(2m+1)=W'(xcp-x) from subroutine cauchy. */
int cmprlb(int *n, int *m, double *x,
	double *g, double *ws, double *wy, double *sy,
	double *wt, double *z__, double *r__, double *wa,
	int *index, double *theta, int *col, int *head,
	int *nfree, logical *cnstnd, int *info)
{
    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset,
	    wt_dim1, wt_offset, i__1, i__2;

    int i, j, k;
    double a1, a2;
    int pointr;

/*     ************ */
    /* Parameter adjustments */
    --index;
    --r__;
    --z__;
    --g;
    --x;
    --wa;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;

    /* Function Body */
    if (! (*cnstnd) && *col > 0) {
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    r__[i] = -g[i];
/* L26: */
	}
    } else {
	i__1 = *nfree;
	for (i = 1; i <= i__1; ++i) {
	    k = index[i];
	    r__[i] = -(*theta) * (z__[k] - x[k]) - g[k];
/* L30: */
	}
	bmv(m, &sy[sy_offset], &wt[wt_offset], col, &wa[(*m << 1) + 1], &wa[
		1], info);
	if (*info != 0) {
	    *info = -8;
	    return 0;
	}
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    a1 = wa[j];
	    a2 = *theta * wa[*col + j];
	    i__2 = *nfree;
	    for (i = 1; i <= i__2; ++i) {
		k = index[i];
		r__[i] = r__[i] + wy[k + pointr * wy_dim1] * a1 + ws[k +
			pointr * ws_dim1] * a2;
/* L32: */
	    }
	    pointr = pointr % *m + 1;
/* L34: */
	}
    }
    return 0;
} /* cmprlb */


/* ======================= The end of cmprlb ============================= */
int formk(int *n, int *nsub, int *ind, int *
	nenter, int *ileave, int *indx2, int *iupdat, logical *
	updatd, double *wn, double *wn1, int *m, double *ws,
	double *wy, double *sy, double *theta, int *col,
	int *head, int *info)
{
    /* System generated locals */
    int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset,
	    wy_dim1, wy_offset, sy_dim1, sy_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, k, k1, m2, is, js, iy, jy, is1, js1, col2, dend, pend;
    int upcl;
    double temp1, temp2, temp3, temp4;
    int ipntr, jpntr, dbegin, pbegin;
/*     ************ */

/*     Subroutine formk */

/*     This subroutine forms  the LEL^T factorization of the indefinite */

/*       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */
/*                                                    where E = [-I  0] */
/*                                                              [ 0  I] */
/*     The matrix K can be shown to be equal to the matrix M^[-1]N */
/*       occurring in section 5.1 of [1], as well as to the matrix */
/*       Mbar^[-1] Nbar in section 5.3. */

/*     n is an int variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     nsub is an int variable */
/*       On entry nsub is the number of subspace variables in free set. */
/*       On exit nsub is not changed. */

/*     ind is an int array of dimension nsub. */
/*       On entry ind specifies the indices of subspace variables. */
/*       On exit ind is unchanged. */

/*     nenter is an int variable. */
/*       On entry nenter is the number of variables entering the */
/*         free set. */
/*       On exit nenter is unchanged. */

/*     ileave is an int variable. */
/*       On entry indx2(ileave),...,indx2(n) are the variables leaving */
/*         the free set. */
/*       On exit ileave is unchanged. */

/*     indx2 is an int array of dimension n. */
/*       On entry indx2(1),...,indx2(nenter) are the variables entering */
/*         the free set, while indx2(ileave),...,indx2(n) are the */
/*         variables leaving the free set. */
/*       On exit indx2 is unchanged. */

/*     iupdat is an int variable. */
/*       On entry iupdat is the total number of BFGS updates made so far. */
/*       On exit iupdat is unchanged. */

/*     updatd is a logical variable. */
/*       On entry 'updatd' is true if the L-BFGS matrix is updatd. */
/*       On exit 'updatd' is unchanged. */

/*     wn is a double precision array of dimension 2m x 2m. */
/*       On entry wn is unspecified. */
/*       On exit the upper triangle of wn stores the LEL^T factorization */
/*         of the 2*col x 2*col indefinite matrix */
/*                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*     wn1 is a double precision array of dimension 2m x 2m. */
/*       On entry wn1 stores the lower triangular part of */
/*                     [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a+R_z   S'AA'S   ] */
/*         in the previous iteration. */
/*       On exit wn1 stores the corresponding updated matrices. */
/*       The purpose of wn1 is just to store these inner products */
/*       so they can be easily updated and inserted into wn. */

/*     m is an int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     ws, wy, sy, and wtyy are double precision arrays; */
/*     theta is a double precision variable; */
/*     col is an int variable; */
/*     head is an int variable. */
/*       On entry they store the information defining the */
/*                                          limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         sy(m,m) stores S'Y; */
/*         wtyy(m,m) stores the Cholesky factorization */
/*                                   of (theta*S'S+LD^(-1)L') */
/*         theta is the scaling factor specifying B_0 = theta I; */
/*         col is the number of variable metric corrections stored; */
/*         head is the location of the 1st s- (or y-) vector in S (or Y). */
/*       On exit they are unchanged. */

/*     info is an int variable. */
/*       On entry info is unspecified. */
/*       On exit info =  0 for normal return; */
/*                    = -1 when the 1st Cholesky factorization failed; */
/*                    = -2 when the 2st Cholesky factorization failed. */

/*     Subprograms called: */

/*       Linpack ... dcopy, dpofa, dtrsl. */


/*     References: */
/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a */
/*       limited memory FORTRAN code for solving bound constrained */
/*       optimization problems'', Tech. Report, NAM-11, EECS Department, */
/*       Northwestern University, 1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Form the lower triangular part of */
/*               WN1 = [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a+R_z   S'AA'S   ] */
/*        where L_a is the strictly lower triangular part of S'AA'Y */
/*              R_z is the upper triangular part of S'ZZ'Y. */
    /* Parameter adjustments */
    --indx2;
    --ind;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    wn1_dim1 = 2 * *m;
    wn1_offset = 1 + wn1_dim1;
    wn1 -= wn1_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;

    /* Function Body */
    if (*updatd) {
	if (*iupdat > *m) {
/*                                 shift old part of WN1. */
	    i__1 = *m - 1;
	    for (jy = 1; jy <= i__1; ++jy) {
		js = *m + jy;
		i__2 = *m - jy;
		dcopy(&i__2, &wn1[jy + 1 + (jy + 1) * wn1_dim1], &c__1, &wn1[
			jy + jy * wn1_dim1], &c__1);
		i__2 = *m - jy;
		dcopy(&i__2, &wn1[js + 1 + (js + 1) * wn1_dim1], &c__1, &wn1[
			js + js * wn1_dim1], &c__1);
		i__2 = *m - 1;
		dcopy(&i__2, &wn1[*m + 2 + (jy + 1) * wn1_dim1], &c__1, &wn1[
			*m + 1 + jy * wn1_dim1], &c__1);
/* L10: */
	    }
	}
/*          put new rows in blocks (1,1), (2,1) and (2,2). */
	pbegin = 1;
	pend = *nsub;
	dbegin = *nsub + 1;
	dend = *n;
	iy = *col;
	is = *m + *col;
	ipntr = *head + *col - 1;
	if (ipntr > *m) {
	    ipntr -= *m;
	}
	jpntr = *head;
	i__1 = *col;
	for (jy = 1; jy <= i__1; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
/*             compute element jy of row 'col' of Y'ZZ'Y */
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
/* L15: */
	    }
/*             compute elements jy of row 'col' of L_a and S'AA'S */
	    i__2 = dend;
	    for (k = dbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L16: */
	    }
	    wn1[iy + jy * wn1_dim1] = temp1;
	    wn1[is + js * wn1_dim1] = temp2;
	    wn1[is + jy * wn1_dim1] = temp3;
	    jpntr = jpntr % *m + 1;
/* L20: */
	}
/*          put new column in block (2,1). */
	jy = *col;
	jpntr = *head + *col - 1;
	if (jpntr > *m) {
	    jpntr -= *m;
	}
	ipntr = *head;
	i__1 = *col;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is = *m + i__;
	    temp3 = 0.;
/*             compute element i of column 'col' of R_z */
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L25: */
	    }
	    ipntr = ipntr % *m + 1;
	    wn1[is + jy * wn1_dim1] = temp3;
/* L30: */
	}
	upcl = *col - 1;
    } else {
	upcl = *col;
    }
/*       modify the old parts in blocks (1,1) and (2,2) due to changes */
/*       in the set of free variables. */
    ipntr = *head;
    i__1 = upcl;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *m + iy;
	jpntr = *head;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    temp4 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
/* L35: */
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
/* L36: */
	    }
	    wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
	    wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
	    jpntr = jpntr % *m + 1;
/* L40: */
	}
	ipntr = ipntr % *m + 1;
/* L45: */
    }
/*       modify the old parts in block (2,1). */
    ipntr = *head;
    i__1 = *m + upcl;
    for (is = *m + 1; is <= i__1; ++is) {
	jpntr = *head;
	i__2 = upcl;
	for (jy = 1; jy <= i__2; ++jy) {
	    temp1 = 0.;
	    temp3 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L50: */
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L51: */
	    }
	    if (is <= jy + *m) {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] + temp1 -
			temp3;
	    } else {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] - temp1 +
			temp3;
	    }
	    jpntr = jpntr % *m + 1;
/* L55: */
	}
	ipntr = ipntr % *m + 1;
/* L60: */
    }
/*     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] */
/*                                     [-L_a +R_z        S'AA'S*theta] */
    m2 = *m << 1;
    i__1 = *col;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *col + iy;
	is1 = *m + iy;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *col + jy;
	    js1 = *m + jy;
	    wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
	    wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
/* L65: */
	}
	i__2 = iy - 1;
	for (jy = 1; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
/* L66: */
	}
	i__2 = *col;
	for (jy = iy; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
/* L67: */
	}
	wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
/* L70: */
    }
/*     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] */
/*                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
/*        first Cholesky factor (1,1) block of wn to get LL' */
/*                          with L' stored in the upper triangle of wn. */
    dpofa(&wn[wn_offset], &m2, col, info);
    if (*info != 0) {
	*info = -1;
	return 0;
    }
/*        then form L^-1(-L_a'+R_z') in the (1,2) block. */
    col2 = *col << 1;
    i__1 = col2;
    for (js = *col + 1; js <= i__1; ++js) {
	dtrsl(&wn[wn_offset], &m2, col, &wn[js * wn_dim1 + 1], &c__11, info);
/* L71: */
    }
/*     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
/*        upper triangle of (2,2) block of wn. */
    i__1 = col2;
    for (is = *col + 1; is <= i__1; ++is) {
	i__2 = col2;
	for (js = is; js <= i__2; ++js) {
	    wn[is + js * wn_dim1] += ddot(col, &wn[is * wn_dim1 + 1], &c__1,
		    &wn[js * wn_dim1 + 1], &c__1);
/* L74: */
	}
/* L72: */
    }
/*     Cholesky factorization of (2,2) block of wn. */
    dpofa(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
    if (*info != 0) {
	*info = -2;
	return 0;
    }
    return 0;
} /* formk */

/* ======================= The end of formk ============================== */



int formt(int *m, double *wt, double *sy,
	double *ss, int *col, double *theta, int *info)
{
    /* System generated locals */
    int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1,
	    i__2, i__3;

    /* Local variables */
    int i__, j, k, k1;
    double ddum;

/*     ************ */

/*     Subroutine formt */

/*       This subroutine forms the upper half of the pos. def. and symm. */
/*         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle */
/*         of the array wt, and performs the Cholesky factorization of T */
/*         to produce J*J', with J' stored in the upper triangle of wt. */

/*     Subprograms called: */

/*       Linpack ... dpofa. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Form the upper half of  T = theta*SS + L*D^(-1)*L', */
/*        store T in the upper triangle of the array wt. */
    /* Parameter adjustments */
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;

    /* Function Body */
    i__1 = *col;
    for (j = 1; j <= i__1; ++j) {
	wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
/* L52: */
    }
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *col;
	for (j = i__; j <= i__2; ++j) {
	    k1 = min(i__,j) - 1;
	    ddum = 0.;
	    i__3 = k1;
	    for (k = 1; k <= i__3; ++k) {
		ddum += sy[i__ + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k +
			k * sy_dim1];
/* L53: */
	    }
	    wt[i__ + j * wt_dim1] = ddum + *theta * ss[i__ + j * ss_dim1];
/* L54: */
	}
/* L55: */
    }
/*     Cholesky factorize T to J*J' with */
/*        J' stored in the upper triangle of wt. */
    dpofa(&wt[wt_offset], m, col, info);
    if (*info != 0) {
	*info = -3;
    }
    return 0;
} /* formt */

/* ======================= The end of formt ============================== */






int freev(int *n, int *nfree, int *index,
	int *nenter, int *ileave, int *indx2, int *iwhere,
	logical *wrk, logical *updatd, logical *cnstnd, int *iprint,
	int *iter)
{
    CORAX_UNUSED(iprint);
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i__, k, iact;


/*     ************ */

/*     Subroutine freev */

/*     This subroutine counts the entering and leaving variables when */
/*       iter > 0, and finds the index set of free and active variables */
/*       at the GCP. */

/*     cnstnd is a logical variable indicating whether bounds are present */

/*     index is an int array of dimension n */
/*       for i=1,...,nfree, index(i) are the indices of free variables */
/*       for i=nfree+1,...,n, index(i) are the indices of bound variables */
/*       On entry after the first iteration, index gives */
/*         the free variables at the previous iteration. */
/*       On exit it gives the free variables based on the determination */
/*         in cauchy using the array iwhere. */

/*     indx2 is an int array of dimension n */
/*       On entry indx2 is unspecified. */
/*       On exit with iter>0, indx2 indicates which variables */
/*          have changed status since the previous iteration. */
/*       For i= 1,...,nenter, indx2(i) have changed from bound to free. */
/*       For i= ileave+1,...,n, indx2(i) have changed from free to bound. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --iwhere;
    --indx2;
    --index;

    /* Function Body */
    *nenter = 0;
    *ileave = *n + 1;
    if (*iter > 0 && *cnstnd) {
        /*                           count the entering and leaving variables. */
        i__1 = *nfree;
        for (i__ = 1; i__ <= i__1; ++i__) {
            k = index[i__];
            /*            write(6,*) ' k  = index(i) ', k */
            /*            write(6,*) ' index = ', i */
            if (iwhere[k] > 0) {
                --(*ileave);
                indx2[*ileave] = k;
                DBG("[L-BFGS-B] Variable %d leaves the set of free variables\n", k );
            }
            /* L20: */
        }
        i__1 = *n;
        for (i__ = *nfree + 1; i__ <= i__1; ++i__) {
            k = index[i__];
            if (iwhere[k] <= 0) {
                ++(*nenter);
                indx2[*nenter] = k;
                DBG("[L-BFGS-B] Variable %d leaves the set of free variables\n", k );
            }
            /* L22: */
        }
        DBG("[L-BFGS-B] %d variables leave; %d variables enter\n", *n + 1 - *ileave, *nenter );
    }
    *wrk = *ileave < *n + 1 || *nenter > 0 || *updatd;
    /*     Find the index set of free and active variables at the GCP. */
    *nfree = 0;
    iact = *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (iwhere[i__] <= 0) {
            ++(*nfree);
            index[*nfree] = i__;
        } else {
            --iact;
            index[iact] = i__;
        }
        /* L24: */
    }
    DBG("[L-BFGS-B] %d variables are free at GCP iter %d\n", *nfree, *iter + 1);
    return 0;
} /* freev */

/* ======================= The end of freev ============================== */







int hpsolb (int *n, double *t, int *iorder, int *iheap)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i__, j, k;
    double out, ddum;
    int indxin, indxou;

/*     ************ */

/*     Subroutine hpsolb */

/*     This subroutine sorts out the least element of t, and puts the */
/*       remaining elements of t in a heap. */

/*     n is an int variable. */
/*       On entry n is the dimension of the arrays t and iorder. */
/*       On exit n is unchanged. */

/*     t is a double precision array of dimension n. */
/*       On entry t stores the elements to be sorted, */
/*       On exit t(n) stores the least elements of t, and t(1) to t(n-1) */
/*         stores the remaining elements in the form of a heap. */

/*     iorder is an int array of dimension n. */
/*       On entry iorder(i) is the index of t(i). */
/*       On exit iorder(i) is still the index of t(i), but iorder may be */
/*         permuted in accordance with t. */

/*     iheap is an int variable specifying the task. */
/*       On entry iheap should be set as follows: */
/*         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap, */
/*         iheap .ne. 0 if otherwise. */
/*       On exit iheap is unchanged. */


/*     References: */
/*       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT. */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */

/*     ************ */
    /* Parameter adjustments */
    --iorder;
    --t;

    /* Function Body */
    if (*iheap == 0) {
/*        Rearrange the elements t(1) to t(n) to form a heap. */
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    ddum = t[k];
	    indxin = iorder[k];
/*           Add ddum to the heap. */
	    i__ = k;
L10:
	    if (i__ > 1) {
		j = i__ / 2;
		if (ddum < t[j]) {
		    t[i__] = t[j];
		    iorder[i__] = iorder[j];
		    i__ = j;
		    goto L10;
		}
	    }
	    t[i__] = ddum;
	    iorder[i__] = indxin;
/* L20: */
	}
    }
/*     Assign to 'out' the value of t(1), the least member of the heap, */
/*        and rearrange the remaining members to form a heap as */
/*        elements 1 to n-1 of t. */
    if (*n > 1) {
	i__ = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[*n];
	indxin = iorder[*n];
/*        Restore the heap */
L30:
	j = i__ + i__;
	if (j <= *n - 1) {
	    if (t[j + 1] < t[j]) {
		++j;
	    }
	    if (t[j] < ddum) {
		t[i__] = t[j];
		iorder[i__] = iorder[j];
		i__ = j;
		goto L30;
	    }
	}
	t[i__] = ddum;
	iorder[i__] = indxin;
/*     Put the least member in t(n). */
	t[*n] = out;
	iorder[*n] = indxou;
    }
    return 0;
} /* hpsolb */

/* ====================== The end of hpsolb ============================== */





int matupd(int *n, int *m, double *ws,
	double *wy, double *sy, double *ss, double *d__,
	double *r__, int *itail, int *iupdat, int *col,
	int *head, double *theta, double *rr, double *dr,
	double *stp, double *dtd)
{
    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset,
	    ss_dim1, ss_offset, i__1, i__2;

    /* Local variables */
    int j;
    int pointr;

/*     ************ */

/*     Subroutine matupd */

/*       This subroutine updates matrices WS and WY, and forms the */
/*         middle matrix in B. */

/*     Subprograms called: */

/*       Linpack ... dcopy, ddot. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Set pointers for matrices WS and WY. */
    /* Parameter adjustments */
    --r__;
    --d__;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;

    /* Function Body */
    if (*iupdat <= *m) {
        *col = *iupdat;
        *itail = (*head + *iupdat - 2) % *m + 1;
    } else {
        *itail = *itail % *m + 1;
        *head = *head % *m + 1;
    }
    /*     Update matrices WS and WY. */
    dcopy(n, &d__[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
    dcopy(n, &r__[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
    /*     Set theta=yy/ys. */
    *theta = *rr / *dr;
    /*     Form the middle matrix in B. */
    /*        update the upper triangle of SS, */
    /*                                         and the lower triangle of SY: */
    if (*iupdat > *m) {
        /*                              move old information */
        i__1 = *col - 1;
        for (j = 1; j <= i__1; ++j) {
            dcopy(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1, &ss[j * ss_dim1 + 1]
                    , &c__1);
            i__2 = *col - j;
            dcopy(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1, &sy[j + j *
                    sy_dim1], &c__1);
            /* L50: */
        }
    }
    /*        add new information: the last row of SY */
    /*                                             and the last column of SS: */
    pointr = *head;
    i__1 = *col - 1;
    for (j = 1; j <= i__1; ++j) {
        sy[*col + j * sy_dim1] = ddot(n, &d__[1], &c__1, &wy[pointr *
                wy_dim1 + 1], &c__1);
        ss[j + *col * ss_dim1] = ddot(n, &ws[pointr * ws_dim1 + 1], &c__1, &
                d__[1], &c__1);
        pointr = pointr % *m + 1;
        /* L51: */
    }
    if (*stp == 1.) {
        ss[*col + *col * ss_dim1] = *dtd;
    } else {
        ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
    }
    sy[*col + *col * sy_dim1] = *dr;
    return 0;
} /* matupd */

/* ======================= The end of matupd ============================= */






int projgr(int *n, double *l, double *u,
        int *nbd, double *x, double *g, double *sbgnrm)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Local variables */
    int i__;
    double gi;

/*     ************ */

/*     Subroutine projgr */

/*     This subroutine computes the infinity norm of the projected */
/*       gradient. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --g;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    *sbgnrm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        gi = g[i__];
        if (nbd[i__] != 0) {
            if (gi < 0.) {
                if (nbd[i__] >= 2) {
                    /* Computing MAX */
                    d__1 = x[i__] - u[i__];
                    gi = max(d__1,gi);
                }
            } else {
                if (nbd[i__] <= 2) {
                    /* Computing MIN */
                    d__1 = x[i__] - l[i__];
                    gi = min(d__1,gi);
                }
            }
        }
        /* Computing MAX */
        d__1 = *sbgnrm, d__2 = abs(gi);
        *sbgnrm = max(d__1,d__2);
        /* L15: */
    }
    return 0;
} /* projgr */

/* ======================= The end of projgr ============================= */









int subsm(int *n, int *m, int *nsub, int *
	ind, double *l, double *u, int *nbd, double *x,
	double *d__, double *xp, double *ws, double *wy,
	double *theta, double *xx, double *gg, int *col,
	int *head, int *iword, double *wv, double *wn, 
	int *iprint, int *info)
{
    CORAX_UNUSED(iprint);

    /* System generated locals */
    int ws_dim1, ws_offset, wy_dim1, wy_offset, wn_dim1, wn_offset, i__1,
	    i__2;
    double d__1, d__2;

    /* Local variables */
    int i__, j, k, m2;
    double dk;
    int js, jy;
    double xk;
    int ibd, col2;
    double dd_p__, temp1, temp2, alpha;
    int pointr;

/*     ********************************************************************** */

/*     This routine contains the major changes in the updated version. */
/*     The changes are described in the accompanying paper */

/*      Jose Luis Morales, Jorge Nocedal */
/*      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale */
/*       Bound Constrained Optimization". Decemmber 27, 2010. */

/*             J.L. Morales  Departamento de Matematicas, */
/*                           Instituto Tecnologico Autonomo de Mexico */
/*                           Mexico D.F. */

/*             J, Nocedal    Department of Electrical Engineering and */
/*                           Computer Science. */
/*                           Northwestern University. Evanston, IL. USA */

/*                           January 17, 2011 */

/*      ********************************************************************** */


/*     Subroutine subsm */

/*     Given xcp, l, u, r, an index set that specifies */
/*       the active set at xcp, and an l-BFGS matrix B */
/*       (in terms of WY, WS, SY, WT, head, col, and theta), */
/*       this subroutine computes an approximate solution */
/*       of the subspace problem */

/*       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp) */

/*             subject to l<=x<=u */
/*                       x_i=xcp_i for all i in A(xcp) */

/*       along the subspace unconstrained Newton direction */

/*          d = -(Z'BZ)^(-1) r. */

/*       The formula for the Newton direction, given the L-BFGS matrix */
/*       and the Sherman-Morrison formula, is */

/*          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r. */

/*       where */
/*                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*     Note that this procedure for computing d differs */
/*     from that described in [1]. One can show that the matrix K is */
/*     equal to the matrix M^[-1]N in that paper. */

/*     n is an int variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     m is an int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     nsub is an int variable. */
/*       On entry nsub is the number of free variables. */
/*       On exit nsub is unchanged. */

/*     ind is an int array of dimension nsub. */
/*       On entry ind specifies the coordinate indices of free variables. */
/*       On exit ind is unchanged. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is a int array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x specifies the Cauchy point xcp. */
/*       On exit x(i) is the minimizer of Q over the subspace of */
/*                                                        free variables. */

/*     d is a double precision array of dimension n. */
/*       On entry d is the reduced gradient of Q at xcp. */
/*       On exit d is the Newton direction of Q. */

/*    xp is a double precision array of dimension n. */
/*       used to safeguard the projected Newton direction */

/*    xx is a double precision array of dimension n */
/*       On entry it holds the current iterate */
/*       On output it is unchanged */
/*    gg is a double precision array of dimension n */
/*       On entry it holds the gradient at the current iterate */
/*       On output it is unchanged */

/*     ws and wy are double precision arrays; */
/*     theta is a double precision variable; */
/*     col is an int variable; */
/*     head is an int variable. */
/*       On entry they store the information defining the */
/*                                          limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         theta is the scaling factor specifying B_0 = theta I; */
/*         col is the number of variable metric corrections stored; */
/*         head is the location of the 1st s- (or y-) vector in S (or Y). */
/*       On exit they are unchanged. */

/*     iword is an int variable. */
/*       On entry iword is unspecified. */
/*       On exit iword specifies the status of the subspace solution. */
/*         iword = 0 if the solution is in the box, */
/*                 1 if some bound is encountered. */

/*     wv is a double precision working array of dimension 2m. */

/*     wn is a double precision array of dimension 2m x 2m. */
/*       On entry the upper triangle of wn stores the LEL^T factorization */
/*         of the indefinite matrix */

/*              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                  [L_a -R_z           theta*S'AA'S ] */
/*                                                    where E = [-I  0] */
/*                                                              [ 0  I] */
/*       On exit wn is unchanged. */

/*     info is an int variable. */
/*       On entry info is unspecified. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return */
/*                                  when the matrix K is ill-conditioned. */

/*     Subprograms called: */

/*       Linpack dtrsl. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */



/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */

    /* Parameter adjustments */
    --gg;
    --xx;
    --xp;
    --d__;
    --x;
    --nbd;
    --u;
    --l;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    --wv;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    --ind;

    /* Function Body */
    if (*nsub <= 0) {
        return 0;
    }
    DBG("[L-BFGS-B] ---------------SUBSM entered---------\n");

    /*     Compute wv = W'Zd. */
    pointr = *head;
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
        temp1 = 0.;
        temp2 = 0.;
        i__2 = *nsub;
        for (j = 1; j <= i__2; ++j) {
            k = ind[j];
            temp1 += wy[k + pointr * wy_dim1] * d__[j];
            temp2 += ws[k + pointr * ws_dim1] * d__[j];
            /* L10: */
        }
        wv[i__] = temp1;
        wv[*col + i__] = *theta * temp2;
        pointr = pointr % *m + 1;
        /* L20: */
    }
    /*     Compute wv:=K^(-1)wv. */
    m2 = *m << 1;
    col2 = *col << 1;
    dtrsl(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
    if (*info != 0) {
        return 0;
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
        wv[i__] = -wv[i__];
        /* L25: */
    }
    dtrsl(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
    if (*info != 0) {
        return 0;
    }
    /*     Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
    pointr = *head;
    i__1 = *col;
    for (jy = 1; jy <= i__1; ++jy) {
        js = *col + jy;
        i__2 = *nsub;
        for (i__ = 1; i__ <= i__2; ++i__) {
            k = ind[i__];
            d__[i__] = d__[i__] + wy[k + pointr * wy_dim1] * wv[jy] / *theta
                + ws[k + pointr * ws_dim1] * wv[js];
            /* L30: */
        }
        pointr = pointr % *m + 1;
        /* L40: */
    }
    d__1 = 1. / *theta;
    dscal(nsub, &d__1, &d__[1], &c__1);

    /* ----------------------------------------------------------------- */
    /*     Let us try the projection, d is the Newton direction */
    *iword = 0;
    dcopy(n, &x[1], &c__1, &xp[1], &c__1);

    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
        k = ind[i__];
        dk = d__[i__];
        xk = x[k];
        if (nbd[k] != 0) {

            if (nbd[k] == 1) {
                /* lower bounds only */
                /* Computing MAX */
                d__1 = l[k], d__2 = xk + dk;
                x[k] = max(d__1,d__2);
                if (x[k] == l[k]) {
                    *iword = 1;
                }
            } else {

                if (nbd[k] == 2) {
                    /* upper and lower bounds */
                    /* Computing MAX */
                    d__1 = l[k], d__2 = xk + dk;
                    xk = max(d__1,d__2);
                    /* Computing MIN */
                    d__1 = u[k];
                    x[k] = min(d__1,xk);
                    if (x[k] == l[k] || x[k] == u[k]) {
                        *iword = 1;
                    }
                } else {

                    if (nbd[k] == 3) {
                        /* upper bounds only */
                        /* Computing MIN */
                        d__1 = u[k], d__2 = xk + dk;
                        x[k] = min(d__1,d__2);
                        if (x[k] == u[k]) {
                            *iword = 1;
                        }
                    }
                }
            }

        } else {
            /* free variables */
            x[k] = xk + dk;
        }
        /* L50: */
    }

    if (*iword == 0) {
        goto L911;
    }

    /*     check sign of the directional derivative */

    dd_p__ = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dd_p__ += (x[i__] - xx[i__]) * gg[i__];
        /* L55: */
    }
    if (dd_p__ > 0.) {
        dcopy(n, &xp[1], &c__1, &x[1], &c__1);
        DBG("[L-BFGS-B] Positive dir derivative in projection \n");
        DBG("[L-BFGS-B] Using the backtracking step\n");
    } else {
        goto L911;
    }

    /* ----------------------------------------------------------------- */

    alpha = 1.;
    temp1 = alpha;
    ibd = 0;
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
        k = ind[i__];
        dk = d__[i__];
        if (nbd[k] != 0) {
            if (dk < 0. && nbd[k] <= 2) {
                temp2 = l[k] - x[k];
                if (temp2 >= 0.) {
                    temp1 = 0.;
                } else if (dk * alpha < temp2) {
                    temp1 = temp2 / dk;
                }
            } else if (dk > 0. && nbd[k] >= 2) {
                temp2 = u[k] - x[k];
                if (temp2 <= 0.) {
                    temp1 = 0.;
                } else if (dk * alpha > temp2) {
                    temp1 = temp2 / dk;
                }
            }
            if (temp1 < alpha) {
                alpha = temp1;
                ibd = i__;
            }
        }
        /* L60: */
    }
    if (alpha < 1.) {
        dk = d__[ibd];
        k = ind[ibd];
        if (dk > 0.) {
            x[k] = u[k];
            d__[ibd] = 0.;
        } else if (dk < 0.) {
            x[k] = l[k];
            d__[ibd] = 0.;
        }
    }
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
        k = ind[i__];
        x[k] += alpha * d__[i__];
        /* L70: */
    }
    /* ccccc */
L911:
    DBG("[L-BFGS-B] ----------------- exit SUBSM --------------\n");
    return 0;
} /* subsm */

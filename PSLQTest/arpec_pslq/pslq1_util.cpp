#include <iostream>
#include <cfloat>
//#include <cmath>

//#include "float_precision.h"
//#include "mp_int.h"

#include "../precision/fprecision.h"

#include "pslq1.h"

//#ifdef HAVE_FP_H
//#include <fp.h>
//#endif

//#ifndef HAVE_COPYSIGN
//#define copysign(x, y) ( ((y) != 0.0) ? \
                         ( ((y) > 0.0) ? (x) : -(x) ) : \
                         ( ((1.0 / y) > 0.0) ? (x) : -(x) ) \
                       )
//#endif

using std::cerr;
using std::endl;

/* Global variables */
int debug_level = 2;
int pslq_iter   = 0;
double timers[NR_TIMERS];

void init_pslq(const matrix<float_precision> &x, matrix<float_precision> &y, 
    matrix<float_precision> &h, matrix<float_precision> &b) {

  int i, j;
  int n = x.size();
  matrix<float_precision> s(n);
  float_precision t = 0.0;
  float_precision u;

  /* Compute partial sums. */
  for (i = n-1; i >= 0; i--) {
      t += x(i) * x(i);//sqr(x(i)); ######
      s(i) = sqrt(t);
  }
  t = 1.0 / s(0);

  /* Normalize vector x and put it into y.  Normalize s as well. */
  for (i = 0; i < n; i++) {
    y(i) = x(i) * t;
    s(i) *= t;
  }

  /* Set matrix B to the identity. */
  b.identity();
    
  /* Set H matrix to lower trapezoidal basis of perp(x). */
  for (j = 0; j < n-1; j++) {
    for (i = 0; i < j; i++) {
      h(i, j) = 0.0;
    }

    t = y(j) / (s(j) * s(j+1));
    h(j, j) = s(j+1) / s(j);

    for (i = j+1; i < n; i++) {
      h(i, j) = -y(i) * t;
    }
  }
}

inline double anint(double a) { return a > 0 ? ceil(a - 0.5) : floor(a + 0.5); }

int reduce_pslq(matrix<float_precision> &h, matrix<float_precision> &y, matrix<float_precision> &b, 
                const float_precision &eps) {
  int i, j, k;
  int n = y.size();
  float_precision t;

  for (i = 1; i < n; i++) {
    for (j = i-1; j >= 0; j--) {
      t = anint(h(i, j) / h(j, j));
      if (t == float_precision(0.0))
        continue;

      y(j) += t * y(i);

      for (k = i; k < n; k++) {
        b(k, j) += t * b(k, i);
      }

      for (int k = 0; k <= j; k++) {
        h(i, k) -= t * h(j, k);
      }
    }
  }

  matrix_min(y, t);
  return (t < eps) ? RESULT_RELATION_FOUND : RESULT_CONTINUE;
}

int iterate_pslq(double gamma, matrix<float_precision> &y, 
                 matrix<float_precision> &h, matrix<float_precision> &b, const float_precision &eps, 
                 const float_precision &teps) {
  int i, j, k;
  int n = y.size();
  float_precision t1, t2, t3, t4;
  double d;
  int im, im1;

  /* Find the diagonal element h(j, j) such that 
     |gamma^j h(j, j)| is maximized.                   */
  t1 = 0.0;
  im = -1;
  d = gamma;
  for (i = 0; i < n-1; i++, d *= gamma) {
    t2 = d * abs(h(i, i));
    if (t2 > t1) {
      im = i;
      t1 = t2;
    }
  }

  if (im == -1) {
    cerr << "ERROR: Invalid index." << endl;
    exit(-1);
  }

  /* Exchange the im and im+1 entries of y, rows of h, and columns of b. */
  im1 = im + 1;

  t1 = y(im);
  y(im) = y(im1);
  y(im1) = t1;

  for (i = 0; i < n-1; i++) {
    t1 = h(im, i);
    h(im, i) = h(im1, i);
    h(im1, i) = t1;
  }

  for (i = 0; i < n; i++) {
    t1 = b(i, im);
    b(i, im) = b(i, im1);
    b(i, im1) = t1;
  }

  /* Update H with permutation produced above. */
  if (im <= n-3) {
    t1 = h(im, im);
    t2 = h(im, im1);
    t3 = 1.0 / sqrt( t1 * t1 + t2 * t2 );
    t1 *= t3;
    t2 *= t3;

    for (i = im; i < n; i++) {
      t3 = h(i, im);
      t4 = h(i, im1);
      h(i, im) = t1 * t3 + t2 * t4;
      h(i, im1) = t1 * t4 - t2 * t3;
    }
  }

  /* Reduce H, updating y, B, and H. */
  for (i = im1; i < n; i++) {
    int j1 = (i == im1) ? i-1 : im1;
    
    for (j = j1; j >= 0; j--) {
      t1 = anint(h(i, j) / h(j, j));
      if (t1 == float_precision(0.0))
        continue;

      y(j) += t1 * y(i);

      for (k = 0; k < n; k++) {
        b(k, j) += t1 * b(k, i);
      }

      for (k = 0; k <= j; k++) {
        h(i, k) -= t1 * h(j, k);
      }
    }
  }

  /* Find the min of |y|. */
  matrix_min(y, t1);
  
  int result = RESULT_CONTINUE;
  if (t1 < teps) {
    if (t1 < eps) {
      result = RESULT_RELATION_FOUND;
    } else {
      result = RESULT_PRECISION_EXHAUSTED;
    }
  }

  return result;

}


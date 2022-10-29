#ifndef PSLQ1_H
#define PSLQ1_H

//#include "tictoc.h"
#include "matrix.h"
#include "../precision/fprecision.h"

#undef inline

/* Mode codes */
#define NR_MODES                   1
#define MODE_ALGEBRAIC_TEST        0

/* Result codes */
#define RESULT_CONTINUE            0
#define RESULT_RELATION_FOUND      1
#define RESULT_PRECISION_EXHAUSTED 2

#define RESULT_RANGE_TOO_LARGE     3
#define RESULT_LARGE_VALUE         4
#define RESULT_VERY_LARGE_VALUE    5

#define RESULT_DUPLICATE           6

/* Constants */
#define LOG_2_BASE_10              3.01029995663981e-01
#define DEFAULT_GAMMA              1.154700538379252

/* Timer index constants */
#define NR_TIMERS                  7
#define TIMER_MP_UPDATE            0
#define TIMER_MPM_UPDATE           1
#define TIMER_MPM_LQ               2
#define TIMER_MP_INIT              3
#define TIMER_MPM_INIT             4
#define TIMER_PSLQ_TOTAL           5
#define TIMER_MPM_ITERATE          6

/* Timer macros */
#define TIMER_BEGIN(n)  { tictoc_t tv;  tic(&tv); 
#define TIMER_END(n)      timers(n) += toc(&tv); }

/* Precision control macros */
#define SET_PREC(n) mp::mpsetprecwords(n)
#define PREC_START  int old_nw = 0; \
                    if (nr_words) { \
                      old_nw = mp::mpgetprecwords(); \
                      mp::mpsetprecwords(nr_words); \
                    }
#define PREC_END    if (nr_words) { \
                      mp::mpsetprecwords(old_nw); \
                    }

/* Global variables */
extern int debug_level;
extern int pslq_iter;
extern double timers[];

/* Level 1 routines.  Found in pslq1_util.cpp and pslq1_templates.cpp. */
//void clear_timers();
//void report_timers();
void init_data(int mode, int n, int r, int s, 
               matrix<float_precision> &x, matrix<int_precision> &ans);

int reduce_pslq(matrix<float_precision> &h, matrix<float_precision> &y, 
                matrix<float_precision> &b, const float_precision &eps);
void init_pslq(const matrix<float_precision> &x, matrix<float_precision> &y, 
               matrix<float_precision> &h, matrix<float_precision> &b);
int iterate_pslq(double gamma, matrix<float_precision> &y, matrix<float_precision> &h, 
                 matrix<float_precision> &b, const float_precision &eps, const float_precision &teps);
int pslq1(const matrix<float_precision> &x, matrix<int_precision> &rel, const float_precision &eps, 
          double gamma = DEFAULT_GAMMA); 

/* Swaps the two elements x and y. */
/*
template <class T>
inline void swap(T &x, T &h) {
  T t = x;
  x = y;
  y = t;
}
*/

/* Swaps the two elementx x and y. 
   Specialization for float_precision type. */
/*
inline void swap(float_precision &x, float_precision &y) {
  float_precision::swap(x, y);
}
*/

template <class T>
int  matrix_minmax(const matrix<T> &v, T &v_min, T &v_max, bool isMin=true);
template <class T>
void lq_decomp(int n, int m, matrix<T> &a);
//template <class T>
//void matmul_left(const matrix<T> &a, matrix<float_precision> &b);
template <class T> 
void bound_pslq(const matrix<T> &h, T &r);

#include "pslq1_templates.cpp"

#endif


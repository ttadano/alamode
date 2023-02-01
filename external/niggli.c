/* Copyright (C) 2015 Atsushi Togo */
/* All rights reserved. */

/* This file is part of niggli. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the niggli project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "niggli.h"

#define NIGGLI_MAX_NUM_LOOP 10000
#define DELAUNAY_MAX_NUM_LOOP 1000000
#define DELAUNAY_CHECK_FREQUENCY 1000

typedef struct {
  double A;
  double B;
  double C;
  double eta;
  double xi;
  double zeta;
  double eps;
  int l;
  int m;
  int n;
  long *tmat;
  long *total_tmat;
  double *lattice;
  double *lattice_orig;
} NiggliParams;

static int run_niggli_reduce(double *lattice_, const double eps_);
static NiggliParams * initialize(const double *lattice_, const double eps_);
static void finalize(double *lattice_, NiggliParams *p);
static int reset(NiggliParams *p);
static int step1(NiggliParams *p);
static int step2(NiggliParams *p);
static int step3(NiggliParams *p);
static int step4(NiggliParams *p);
static int step5(NiggliParams *p);
static int step6(NiggliParams *p);
static int step7(NiggliParams *p);
static int step8(NiggliParams *p);
static int set_parameters(NiggliParams *p);
static void set_angle_types(NiggliParams *p);
static int update_lattice(NiggliParams *p);
static int update_total_tmat(NiggliParams *p);

static int run_delaunay_reduction(double *lattice_, const double eps_);
static int purify_basis(double basis[4][3],
                        const double *lattice_,
                        const double eps_);
static int delaunay_reduce_basis(double basis[4][3],
                                 const double symprec);
static void get_delaunay_exteneded_basis(double basis[4][3],
                                         const double *lattice);
static int get_delaunay_shortest_vectors(double basis[4][3],
                                         const double eps_);
static double * get_tmat(const double * orig_lat,
                         const double * cur_lat,
                         const double eps_);
static double * transpose(const double *M);
static double * get_metric(const double *M);
static double norm_squared(const double a[3]);
static void swap_vectors(double a[3], double b[3]);
static double * inverse(const double *m, const double eps_);
static double determinant(const double *m);
static double * multiply(const double *L, const double *R);
static long nint(const double a);

#ifdef NIGGLI_DEBUG
#define debug_print(...) printf(__VA_ARGS__)
static void debug_show(const int j, const NiggliParams *p);
static void debug_show(const int j, const NiggliParams *p)
{
  /* int i; */

  if (j < 0) {
    printf("Finish: ");
  } else {
    printf("Step %d: ", j);
  }
  printf("%f %f %f %f %f %f\n", p->A, p->B, p->C, p->xi, p->eta, p->zeta);

  /* printf("%d %d %d\n", p->l, p->m, p->n); */
  /* for (i = 0; i < 3; i++) { */
  /*   printf("%f %f %f\n", */
  /*       p->lattice[i * 3], p->lattice[i * 3 + 1], p->lattice[i * 3 + 2]); */
  /* } */
}
#else
#define debug_print(...)
#define debug_show(...)
#endif

#ifdef NIGGLI_WARNING
#define warning_print(...) fprintf(stderr,__VA_ARGS__)
#else
#define warning_print(...)
#endif

/*--------------------------------------------*/
/* Version: niggli-[major].[minor].[micro] */
/*--------------------------------------------*/
int niggli_get_major_version(void)
{
  return NIGGLI_MAJOR_VERSION;
}

int niggli_get_minor_version(void)
{
  return NIGGLI_MINOR_VERSION;
}

int niggli_get_micro_version(void)
{
  return NIGGLI_MICRO_VERSION;
}

/* return 0 if failed */
int niggli_reduce(double *lattice_, const double eps_)
{
  int succeeded;

  succeeded = 0;
  succeeded = run_niggli_reduce(lattice_, eps_);

  if (! succeeded) {
    /* Fallback to Delaunay reduction. */
    /* This may induce round-off error when the number of iterations */
    /* is large. */
    run_delaunay_reduction(lattice_, eps_);
    succeeded = run_niggli_reduce(lattice_, eps_);
  }

  return succeeded;
}

static int run_niggli_reduce(double *lattice_, const double eps_)
{
  int i, j, succeeded;
  NiggliParams *p;
  int (*steps[8])(NiggliParams *p) = {step1, step2, step3, step4,
                                      step5, step6, step7, step8};

  p = NULL;
  succeeded = 0;

  if ((p = initialize(lattice_, eps_)) == NULL) {
    return 0;
  }

  /* Step 0 */
  if (! set_parameters(p)) {
    goto err;
  }

  for (i = 0; i < NIGGLI_MAX_NUM_LOOP; i++) {
    for (j = 0; j < 8; j++) {
      debug_show(j + 1, p);
      if ((*steps[j])(p)) {
        if (! reset(p)) {goto err;}
        if (j == 1 || j == 4 || j == 5 || j == 6 || j == 7) {break;}
      }
    }
    if (j == 8) {
      succeeded = 1;
      break;
    }
  }

  debug_show(-1, p);

  if (succeeded) {
    memcpy(lattice_, p->lattice, sizeof(double) * 9);
  }

 err:
  finalize(lattice_, p);
  return succeeded;
}

static NiggliParams * initialize(const double *lattice_, const double eps_)
{
  NiggliParams * p;

  p = NULL;

  if ((p = (NiggliParams*)malloc(sizeof(NiggliParams))) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    return NULL;
  }

  p->A = 0;
  p->B = 0;
  p->C = 0;
  p->eta = 0;
  p->xi = 0;
  p->zeta = 0;
  p->eps = 0;
  p->l = 0;
  p->m = 0;
  p->n = 0;
  p->tmat = NULL;
  p->total_tmat = NULL;
  p->lattice = NULL;
  p->lattice_orig = NULL;

  if ((p->tmat = (long*)malloc(sizeof(long) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    free(p);
    p = NULL;
    return NULL;
  }

  if ((p->total_tmat = (long*)malloc(sizeof(long) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    free(p->tmat);
    p->tmat = NULL;
    free(p);
    p = NULL;
    return NULL;
  }

  if ((p->lattice = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    free(p->total_tmat);
    p->total_tmat = NULL;
    free(p->tmat);
    p->tmat = NULL;
    free(p);
    p = NULL;
    return NULL;
  }

  if ((p->lattice_orig = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    free(p->lattice);
    p->lattice = NULL;
    free(p->total_tmat);
    p->total_tmat = NULL;
    free(p->tmat);
    p->tmat = NULL;
    free(p);
    p = NULL;
    return NULL;
  }

  p->eps = eps_;
  memcpy(p->lattice, lattice_, sizeof(double) * 9);
  memcpy(p->lattice_orig, lattice_, sizeof(double) * 9);
  p->total_tmat[0] = 1;
  p->total_tmat[1] = 0;
  p->total_tmat[2] = 0;
  p->total_tmat[3] = 0;
  p->total_tmat[4] = 1;
  p->total_tmat[5] = 0;
  p->total_tmat[6] = 0;
  p->total_tmat[7] = 0;
  p->total_tmat[8] = 1;

  return p;
}

static void finalize(double *lattice_, NiggliParams *p)
{
  free(p->total_tmat);
  p->tmat = NULL;
  free(p->tmat);
  p->tmat = NULL;
  free(p->lattice_orig);
  p->lattice_orig = NULL;
  free(p->lattice);
  p->lattice = NULL;
  free(p);
  p = NULL;
}

static int reset(NiggliParams *p)
{
  if (! update_total_tmat(p)) {return 0;}
  if (! update_lattice(p)) {return 0;}
  return set_parameters(p);
}

static int set_parameters(NiggliParams *p)
{
  double *G;

  G = NULL;

  if ((G = get_metric(p->lattice)) == NULL) {return 0;}

  p->A = G[0];
  p->B = G[4];
  p->C = G[8];
  p->xi = G[5] * 2;
  p->eta = G[2] * 2;
  p->zeta = G[1] * 2;

  free(G);
  G = NULL;

  set_angle_types(p);

  return 1;
}

static void set_angle_types(NiggliParams *p)
{
  p->l = 0;
  p->m = 0;
  p->n = 0;
  if (p->xi < -p->eps) {p->l = -1;}
  if (p->xi > p->eps) {p->l = 1;}
  if (p->eta < -p->eps) {p->m = -1;}
  if (p->eta > p->eps) {p->m = 1;}
  if (p->zeta < -p->eps) {p->n = -1;}
  if (p->zeta > p->eps) {p->n = 1;}
}

static int step1(NiggliParams *p)
{
  if (p->A > p->B + p->eps ||
      (! (fabs(p->A - p->B) > p->eps) &&
       fabs(p->xi) > fabs(p->eta) + p->eps)) {
    p->tmat[0] = 0,  p->tmat[1] = -1, p->tmat[2] = 0;
    p->tmat[3] = -1, p->tmat[4] = 0,  p->tmat[5] = 0;
    p->tmat[6] = 0,  p->tmat[7] = 0,  p->tmat[8] = -1;
    return 1;
  }
  else {return 0;}
}

static int step2(NiggliParams *p)
{
  if (p->B > p->C + p->eps ||
      (! (fabs(p->B - p->C) > p->eps)
       && fabs(p->eta) > fabs(p->zeta) + p->eps)) {
    p->tmat[0] = -1, p->tmat[1] = 0,  p->tmat[2] = 0;
    p->tmat[3] = 0,  p->tmat[4] = 0,  p->tmat[5] = -1;
    p->tmat[6] = 0,  p->tmat[7] = -1, p->tmat[8] = 0;
    return 1;
  }
  else {return 0;}
}

static int step3(NiggliParams *p)
{
  int i, j, k;
  if (p->l * p->m * p->n == 1) {
    if (p->l == -1) {i = -1;} else {i = 1;}
    if (p->m == -1) {j = -1;} else {j = 1;}
    if (p->n == -1) {k = -1;} else {k = 1;}
    p->tmat[0] = i, p->tmat[1] = 0, p->tmat[2] = 0;
    p->tmat[3] = 0, p->tmat[4] = j, p->tmat[5] = 0;
    p->tmat[6] = 0, p->tmat[7] = 0, p->tmat[8] = k;
    return 1;
  }
  else {return 0;}
}

static int step4(NiggliParams *p)
{
  int i, j, k, r;

  if (p->l == -1 && p->m == -1 && p->n == -1) {
    return 0;
  }

  if (p->l * p->m * p->n == 0 || p->l * p->m * p->n == -1) {
    i = 1;
    j = 1;
    k = 1;
    r = -1; /* 0: i, 1: j, 2: k */
    if (p->l == 1) {i = -1;}
    if (p->l == 0) {r = 0;}
    if (p->m == 1) {j = -1;}
    if (p->m == 0) {r = 1;}
    if (p->n == 1) {k = -1;}
    if (p->n == 0) {r = 2;}

    if (i * j * k == -1) {
      if (r == 0) {i = -1;}
      if (r == 1) {j = -1;}
      if (r == 2) {k = -1;}
    }

    p->tmat[0] = i, p->tmat[1] = 0, p->tmat[2] = 0;
    p->tmat[3] = 0, p->tmat[4] = j, p->tmat[5] = 0;
    p->tmat[6] = 0, p->tmat[7] = 0, p->tmat[8] = k;
    return 1;
  }
  else {return 0;}
}

static int step5(NiggliParams *p)
{
  if (fabs(p->xi) > p->B + p->eps ||
      (! (fabs(p->B - p->xi) > p->eps) && 2 * p->eta < p->zeta - p->eps) ||
      (! (fabs(p->B + p->xi) > p->eps) && p->zeta < -p->eps)) {
    p->tmat[0] = 1, p->tmat[1] = 0, p->tmat[2] = 0;
    p->tmat[3] = 0, p->tmat[4] = 1, p->tmat[5] = 0;
    p->tmat[6] = 0, p->tmat[7] = 0, p->tmat[8] = 1;
    if (p->xi > 0) {p->tmat[5] = -1;}
    if (p->xi < 0) {p->tmat[5] = 1;}
    return 1;
  }
  else {return 0;}
}

static int step6(NiggliParams *p)
{
  if (fabs(p->eta) > p->A + p->eps ||
      (! (fabs(p->A - p->eta) > p->eps) && 2 * p->xi < p->zeta - p->eps) ||
      (! (fabs(p->A + p->eta) > p->eps) && p->zeta < -p->eps)) {
    p->tmat[0] = 1, p->tmat[1] = 0, p->tmat[2] = 0;
    p->tmat[3] = 0, p->tmat[4] = 1, p->tmat[5] = 0;
    p->tmat[6] = 0, p->tmat[7] = 0, p->tmat[8] = 1;
    if (p->eta > 0) {p->tmat[2] = -1;}
    if (p->eta < 0) {p->tmat[2] = 1;}
    return 1;
  }
  else {return 0;}
}

static int step7(NiggliParams *p)
{
  if (fabs(p->zeta) > p->A + p->eps ||
      (! (fabs(p->A - p->zeta) > p->eps) && 2 * p->xi < p->eta - p->eps) ||
      (! (fabs(p->A + p->zeta) > p->eps) && p->eta < -p->eps)) {
    p->tmat[0] = 1, p->tmat[1] = 0, p->tmat[2] = 0;
    p->tmat[3] = 0, p->tmat[4] = 1, p->tmat[5] = 0;
    p->tmat[6] = 0, p->tmat[7] = 0, p->tmat[8] = 1;
    if (p->zeta > 0) {p->tmat[1] = -1;}
    if (p->zeta < 0) {p->tmat[1] = 1;}
    return 1;
  }
  else {return 0;}
}

static int step8(NiggliParams *p)
{
  if (p->xi + p->eta + p->zeta + p->A + p->B < -p->eps ||
      (! (fabs(p->xi + p->eta + p->zeta + p->A + p->B) > p->eps) &&
       2 * (p->A + p->eta) + p->zeta > p->eps)) {
    p->tmat[0] = 1, p->tmat[1] = 0, p->tmat[2] = 1;
    p->tmat[3] = 0, p->tmat[4] = 1, p->tmat[5] = 1;
    p->tmat[6] = 0, p->tmat[7] = 0, p->tmat[8] = 1;
    return 1;
  }
  else {return 0;}
}

static int update_lattice(NiggliParams *p)
{
  int i, j, k;
  double M[9];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M[i * 3 + j] = 0;
      for (k = 0; k < 3; k++) {
        M[i * 3 + j] += p->lattice_orig[i * 3 + k] * p->total_tmat[k * 3 + j];
      }
    }
  }

  memcpy(p->lattice, M, sizeof(double) * 9);

  return 1;
}

static int update_total_tmat(NiggliParams *p)
{
  int i, j, k;
  long M[9];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M[i * 3 + j] = 0;
      for (k = 0; k < 3; k++) {
        M[i * 3 + j] += p->total_tmat[i * 3 + k] * p->tmat[k * 3 + j];
      }
    }
  }

  memcpy(p->total_tmat, M, sizeof(long) * 9);

  return 1;
}


static int run_delaunay_reduction(double *lattice_, const double eps_)
{
  int i, j, attempt, succeeded;
  double basis[4][3];

  succeeded = 0;
  get_delaunay_exteneded_basis(basis, lattice_);

  for (attempt = 1; attempt < DELAUNAY_MAX_NUM_LOOP; attempt++) {
    succeeded = delaunay_reduce_basis(basis, eps_);
    if (succeeded) {
      break;
    }
    if (attempt % DELAUNAY_CHECK_FREQUENCY == 0) {
      if (!purify_basis(basis, lattice_, eps_)) {
        goto err;
      }
    }
  }

  if (!succeeded) {
    goto err;
  }

  if (!purify_basis(basis, lattice_, eps_)) {
    succeeded = 0;
    goto err;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lattice_[i * 3 + j] = basis[j][i];
    }
  }

err:
  return succeeded;
}

static double * get_tmat(const double * orig_lat,
                         const double * cur_lat,
                         const double eps_)
{
  double *inv, *tmat;

  inv = NULL;
  tmat = NULL;

  if ((inv = inverse(orig_lat, eps_)) == NULL) {
    goto err;
  }

  if ((tmat = multiply(inv, cur_lat)) == NULL) {
    goto err;
  }

err:
  if (inv != NULL) {
    free(inv);
    inv = NULL;
  }
  return tmat;
}

static int purify_basis(double basis[4][3],
                        const double *lattice_,
                        const double eps_)
{
  int i, j, succeeded;
  double *inv, *tmat, *p_cur_lat;
  double cur_lat[9];

  succeeded = 0;
  inv = NULL;
  tmat = NULL;
  p_cur_lat = NULL;

  if (!get_delaunay_shortest_vectors(basis, eps_)) {
    goto err;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      cur_lat[i * 3 + j] = basis[j][i];
    }
  }

  if ((tmat = get_tmat(lattice_, cur_lat, eps_)) == NULL) {
    goto err;
  }

  if (fabs(determinant(tmat) - 1) > eps_) {
    goto err;
  }

  for (i = 0; i < 9; i++) {
    tmat[i] = nint(tmat[i]);
  }

  if (fabs(determinant(tmat) - 1) > eps_) {
    goto err;
  }

  if ((p_cur_lat = multiply(lattice_, tmat)) == NULL) {
    goto err;
  }

  get_delaunay_exteneded_basis(basis, p_cur_lat);
  succeeded = 1;

err:

  if (p_cur_lat != NULL) {
    free(p_cur_lat);
    p_cur_lat = NULL;
  }
  if (tmat != NULL) {
    free(tmat);
    tmat = NULL;
  }
  if (inv != NULL) {
    free(inv);
    inv = NULL;
  }

  return succeeded;
}

static int delaunay_reduce_basis(double basis[4][3],
                                 const double symprec)
{
  int i, j, k, l;
  double dot_product;

  for (i = 0; i < 4; i++) {
    for (j = i+1; j < 4; j++) {
      dot_product = 0.0;
      for (k = 0; k < 3; k++) {
        dot_product += basis[i][k] * basis[j][k];
      }
      if (dot_product > symprec) {
        for (k = 0; k < 4; k++) {
          if (! (k == i || k == j)) {
            for (l = 0; l < 3; l++) {
              basis[k][l] += basis[i][l];
            }
          }
        }
        for (k = 0; k < 3; k++) {
          basis[i][k] = -basis[i][k];
        }
        return 0;
      }
    }
  }

  return 1;
}

static void get_delaunay_exteneded_basis(double basis[4][3],
                                         const double *lattice)
{
  int i, j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      basis[i][j] = lattice[j * 3 + i];
    }
  }

  for (i = 0; i < 3; i++) {
    basis[3][i] = -lattice[i * 3] -lattice[i * 3 + 1] -lattice[i * 3 + 2];
  }
}

static int get_delaunay_shortest_vectors(double basis[4][3],
                                         const double eps_)
{
  int i, j, k, succeeded;
  double det;
  double b[7][3];

  succeeded = 0;

  /* Search in the set {b1, b2, b3, b4, b1+b2, b2+b3, b3+b1} */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      b[i][j] = basis[i][j];
    }
  }

  for (i = 0; i < 3; i++) {
    b[4][i] = basis[0][i] + basis[1][i];
  }
  for (i = 0; i < 3; i++) {
    b[5][i] = basis[1][i] + basis[2][i];
  }
  for (i = 0; i < 3; i++) {
    b[6][i] = basis[2][i] + basis[0][i];
  }

  /* Bubble sort */
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      if (norm_squared(b[j]) > norm_squared(b[j + 1]) + eps_) {
        swap_vectors(b[j], b[j + 1]);
      }
    }
  }

  for (i = 1; i < 6; i++) {
    for (j = i + 1; j < 7; j++) {
      det = (b[0][0] * (b[i][1] * b[j][2] - b[i][2] * b[j][1]) +
             b[0][1] * (b[i][2] * b[j][0] - b[i][0] * b[j][2]) +
             b[0][2] * (b[i][0] * b[j][1] - b[i][1] * b[j][0]));
      if (fabs(det) > eps_) {
        for (k = 0; k < 3; k++) {
          if (det > 0) {
            basis[0][k] = b[0][k];
          } else {
            basis[0][k] = -b[0][k];
          }
          basis[1][k] = b[i][k];
          basis[2][k] = b[j][k];
        }
        succeeded = 1;
        goto ret;
      }
    }
  }

ret:
  return succeeded;
}

static double * get_metric(const double *M)
{
  double *G, *M_T;

  G = NULL;
  M_T = NULL;

  if ((M_T = transpose(M)) == NULL) {
    return NULL;
  }

  if ((G = multiply(M_T, M)) == NULL) {
    free(M_T);
    M_T = NULL;
    return NULL;
  }

  free(M_T);
  M_T = NULL;

  return G;
}

static double * multiply(const double *L, const double *R)
{
  int i, j, k;
  double *M;

  M = NULL;

  if ((M = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    return NULL;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M[i * 3 + j] = 0;
      for (k = 0; k < 3; k++) {
        M[i * 3 + j] += L[i * 3 + k] * R[k * 3 + j];
      }
    }
  }

  return M;
}

static double * transpose(const double *M)
{
  int i, j;
  double *M_T;

  M_T = NULL;

  if ((M_T = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    return NULL;
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M_T[i * 3 + j] = M[j * 3 + i];
    }
  }

  return M_T;
}

static double norm_squared(const double a[3])
{
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

static void swap_vectors(double a[3], double b[3])
{
  int i;
  double tmp;

  for (i = 0; i < 3; i++) {
    tmp = a[i];
    a[i] = b[i];
    b[i] = tmp;
  }
}

static double * inverse(const double *m, const double eps_)
{
  double det;
  double *inv;

  inv = NULL;

  det = determinant(m);
  if (fabs(det) < eps_) {
    goto err;
  }

  if ((inv = (double*)malloc(sizeof(double) * 9)) == NULL) {
    warning_print("niggli: Memory could not be allocated.");
    goto err;
  }

  inv[0] = (m[4] * m[8] - m[5] * m[7]) / det;
  inv[1] = (m[7] * m[2] - m[8] * m[1]) / det;
  inv[2] = (m[1] * m[5] - m[2] * m[4]) / det;
  inv[3] = (m[5] * m[6] - m[3] * m[8]) / det;
  inv[4] = (m[8] * m[0] - m[6] * m[2]) / det;
  inv[6] = (m[3] * m[7] - m[4] * m[6]) / det;
  inv[5] = (m[2] * m[3] - m[0] * m[5]) / det;
  inv[7] = (m[6] * m[1] - m[7] * m[0]) / det;
  inv[8] = (m[0] * m[4] - m[1] * m[3]) / det;

err:
  return inv;
}

static double determinant(const double *m)
{
  return (m[0] * (m[4] * m[8] - m[5] * m[7]) +
          m[1] * (m[5] * m[6] - m[3] * m[8]) +
          m[2] * (m[3] * m[7] - m[4] * m[6]));
}

static long nint(const double a)
{
  if (a < 0.0)
    return (long) (a - 0.5);
  else
    return (long) (a + 0.5);
}

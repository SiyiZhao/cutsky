/*******************************************************************************
* survey_geom.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "survey_geom.h"
#include "read_data.h"
#include "cspline.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/*============================================================================*\
                     Functions for applying survey geometry
\*============================================================================*/

/******************************************************************************
Function `load_nz`:
  Load n(z) from file and prepare for the interpolation.
Arguments:
  * `geom`:     interface for survey geometry;
  * `fname`:    name of the file storing n(z);
  * `zmin`:     minimum redshift of interest;
  * `zmax`:     maximum redshift of interest.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int load_nz(GEOM *geom, const char *fname, const double zmin,
    const double zmax) {
  size_t num = 0;
  if (read_ascii_twocol(fname, &geom->z, &geom->nz, &num)) {
    P_ERR("failed to read the n(z) file: `%s'\n", fname);
    return CUTSKY_ERR_FILE;
  }

  if (num > INT_MAX) {
    P_ERR("too many entries in file: `%s'\n", fname);
    return CUTSKY_ERR_FILE;
  }

  /* Validate the redshifts in file. */
  for (size_t i = 1; i < num; i++) {
    if (geom->z[i] <= geom->z[i - 1]) {
      P_ERR("redshifts must be in ascending order in file: `%s'\n", fname);
      return CUTSKY_ERR_FILE;
    }
  }
  if (geom->z[0] > zmin || geom->z[num - 1] < zmax) {
    P_ERR("redshifts in file must cover the range (" OFMT_DBL ", " OFMT_DBL
        "): `%s'\n", zmin, zmax, fname);
    return CUTSKY_ERR_FILE;
  }

  /* Compute the second derivative of n(z) for interpolation. */
  if (!(geom->nzpp = malloc(num * 2 * sizeof(double)))) {
    P_ERR("failed to allocate memory for n(z) interpolation\n");
    return CUTSKY_ERR_MEMORY;
  }

  cspline_ypp(geom->z, geom->nz, num, geom->nzpp);
  geom->nsp = num;

  return 0;
}

/******************************************************************************
Function `bin_search`:
  Binary search the x coordinate for interpolation.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `n`:        number of the sample points;
  * `xv`:       x coordinate of the value to be evaluated;
Return:
  Index of the value in the sample to be evaluated.
******************************************************************************/
static inline int bin_search(const double *x, const int n, const double xv) {
  int l = 0;
  int u = n - 1;
  while (l <= u) {
    int i = ((unsigned int) l + (unsigned int) u) >> 1;
    if (i >= n - 1) {
      if (x[n - 1] == xv) return n - 1;
      else return INT_MAX;
    }
    if (x[i + 1] <= xv) l = i + 1;
    else if (x[i] > xv) u = i - 1;
    else return i;
  }
  return INT_MAX;
}

/*============================================================================*\
                    Interfaces for applying survey geometry
\*============================================================================*/

/******************************************************************************
Function `geom_init`:
  Initialise the interface for applying survey geometry.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Interface for applying survey geometry.
******************************************************************************/
GEOM *geom_init(const CONF *conf) {
  printf("Setting up survey geometry ...");
  if (!conf) {
    P_ERR("configurations are not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

  GEOM *geom = malloc(sizeof *geom);
  if (!geom) {
    P_ERR("failed to allocate memory for survey geometry\n");
    return NULL;
  }
  geom->foot[0] = geom->foot[1] = NULL;
  geom->rng = NULL;
  geom->z = geom->nz = geom->nzpp = NULL;
  geom->nsp = 0;
  geom->infoot = CUTSKY_BITCODE_INFOOT;
  geom->rad_sel = CUTSKY_BITCODE_RAD_SEL;

  /* Process Mangle polygon-format footprints. */
  int err = 0;
  geom->foot[0] = mangle_init(conf->foot_all, CUTSKY_WMIN_FOOT_ALL, &err);
  if (!(geom->foot[0]) || err) {
    P_ERR("failed to process the entire footprint: %s\n", mangle_errmsg(err));
    geom_destroy(geom);
    return NULL;
  }
  if (conf->verbose)
    printf("  The entire DESI footprint is loaded from `%s'\n", conf->foot_all);

  if (conf->foot) {
    geom->foot[1] = mangle_init(conf->foot, CUTSKY_WMIN_FOOT, &err);
    if (!(geom->foot[1]) || err) {
      P_ERR("failed to process the footprint of interest: %s\n",
          mangle_errmsg(err));
      geom_destroy(geom);
      return NULL;
    }
    if (conf->verbose)
      printf("  The footprint of interest is loaded from `%s'\n", conf->foot);
  }

  /* Load the n(z) file and prepare for the interpolation. */
  if (conf->fnz) {
    if (load_nz(geom, conf->fnz, conf->zmin, conf->zmax)) {
      geom_destroy(geom);
      return NULL;
    }

    /* Initialise random number generator. */
    int err = 0;
    for (int i = 0; i < conf->ncap; i++) geom->seed[i] = conf->seed[i];
#ifdef OMP
    geom->rng = prand_init(conf->rng, geom->seed[0], conf->nthread, 0, &err);
#else
    geom->rng = prand_init(conf->rng, geom->seed[0], 1, 0, &err);
#endif
    if (PRAND_IS_ERROR(err) || PRAND_IS_WARN(err)) {
      P_ERR("failed to initialise the random number generator\n");
      geom_destroy(geom);
      return NULL;
    }

    if (conf->verbose)
      printf("  %d n(z) samples read from file `%s'\n", geom->nsp, conf->fnz);
  }

  printf(FMT_DONE);
  return geom;
}

/******************************************************************************
Function `geom_destroy`:
  Deconstruct the interface for applying survey geometry.
Arguments:
  * `geom`:     interface for survey geometry.
******************************************************************************/
void geom_destroy(GEOM *geom) {
  if (!geom) return;
  if (geom->foot[0]) mangle_destroy(geom->foot[0]);
  if (geom->foot[1]) mangle_destroy(geom->foot[1]);
  if (geom->rng) prand_destroy(geom->rng);
  if (geom->z) free(geom->z);
  if (geom->nz) free(geom->nz);
  if (geom->nzpp) free(geom->nzpp);
  free(geom);
}

/******************************************************************************
Function `geom_get_nz`:
  Compute the expected comoving density given redshift by interpolating n(z).
Arguments:
  * `geom`:     interface for survey geometry;
  * `z`:        the given redshift.
Return:
  The comoving density on success; HUGE_VAL on error.
******************************************************************************/
double geom_get_nz(const GEOM *geom, const double z) {
  int idx = bin_search(geom->z, geom->nsp, z);
  if (idx == INT_MAX) return HUGE_VAL;

  double nz = (idx == geom->nsp - 1) ? geom->nz[idx] :
      cspline_eval(geom->z, geom->nz, geom->nzpp, z, idx);
  return nz;
}

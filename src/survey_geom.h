/*******************************************************************************
* survey_geom.h: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __SURVEY_GEOM_H__
#define __SURVEY_GEOM_H__

#include "load_conf.h"
#include "mangle.h"
#include "prand.h"
#include <stdio.h>

/*============================================================================*\
                       Data structure for survey geometry
\*============================================================================*/

typedef struct {
  MANGLE *foot[2];      /* DESI entire and current footprints      */
  prand_t *rng;         /* interface of random number generator    */
  size_t seed[2];       /* random seeds                            */
  double *z;            /* array for redshift values               */
  double *nz;           /* array for comoving number densities     */
  double *nzpp;         /* second derivative of comoving densities */
  int nsp;              /* number of n(z) samples                  */
  uint8_t infoot;       /* bitcode for the current footprint       */
  uint8_t rad_sel;      /* bitcode for radial selection            */
} GEOM;

/*============================================================================*\
                    Interfaces for applying survey geometry
\*============================================================================*/

/******************************************************************************
Macro `geom_infoot`:
  Check if a coordinate is inside a given footprint.
Arguments:
  * `foot`:     the footprint;
  * `ra`:       right ascension of the object;
  * `dec`:      declination of the object.
Return:
  NULL if the coordinate is NOT inside the footprint.
******************************************************************************/
#define geom_infoot(foot, ra, dec)      mangle_query(foot, ra, dec)

/******************************************************************************
Function `geom_init`:
  Initialise the interface for applying survey geometry.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Interface for applying survey geometry.
******************************************************************************/
GEOM *geom_init(const CONF *conf);

/******************************************************************************
Function `geom_destroy`:
  Deconstruct the interface for applying survey geometry.
Arguments:
  * `geom`:     interface for survey geometry.
******************************************************************************/
void geom_destroy(GEOM *geom);

/******************************************************************************
Function `geom_get_nz`:
  Compute the expected comoving density given redshift by interpolating n(z).
Arguments:
  * `geom`:     interface for survey geometry;
  * `z`:        the given redshift.
Return:
  The comoving density on success; HUGE_VAL on error.
******************************************************************************/
double geom_get_nz(const GEOM *geom, const double z);

#endif

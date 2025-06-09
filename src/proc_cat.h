/*******************************************************************************
* proc_cat.h: this file is part of the cutsky program.

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

#ifndef __PROC_CAT_H__
#define __PROC_CAT_H__

#include "load_conf.h"
#include "convert_z.h"
#include "survey_geom.h"
#include <stdio.h>
#include <stdint.h>

/*============================================================================*\
                    Data structures for the cutsky catalogue
\*============================================================================*/

typedef struct {
  size_t n;             /* number of objects in the catalog           */
  size_t max;           /* capacity of the catalog                    */
  float *x[4];          /* coordinates: RA, DEC, Z_RSD, Z_REAL        */
  float *nz;            /* comoving number density                    */
  float *ran;           /* random number for radial selection         */
  uint8_t *status;      /* bitcode for footprint and radial selection */
} DATA;

/*============================================================================*\
                     Interface for the main cutsky process
\*============================================================================*/

/******************************************************************************
Function `process`:
  Read the input catalogue, apply cutsky geometry, and save the results.
Arguments:
  * `conf`:     structure for storing configurations;
  * `zcvt`:     interface for redshift conversion;
  * `geom`:     interface for survey geometry.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int process(const CONF *conf, const ZCVT *zcvt, const GEOM *geom);

#endif

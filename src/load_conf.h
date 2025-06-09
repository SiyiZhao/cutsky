/*******************************************************************************
* load_conf.h: this file is part of the cutsky program.

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

#ifndef __LOAD_CONF_H__
#define __LOAD_CONF_H__

#include <stdbool.h>

/*============================================================================*\
                   Data structure for storing configurations
\*============================================================================*/

typedef struct {
  char *fconf;          /* name of the configuration file */
  char *input;          /* INPUT           */
  int ifmt;             /* INPUT_FORMAT    */
  char **inputs;        /* Input catalogues. */
  int ninput;           /* number of input catalogues */
  char comment;         /* COMMENT         */
  double Lbox;          /* BOX_SIZE        */
  long ndata;           /* NUMBER          */
  double omega_m;       /* OMEGA_M         */
  double omega_l;       /* OMEGA_LAMBDA    */
  double omega_k;       /* 1 - OMEGA_M - OMEGA_LAMBDA */
  double eos_w;         /* DE_EOS_W        */
  char *fzcnvt;         /* Z_CMVDST_CNVT   */
  char *foot_all;       /* DESI_ALL_TILES  */
  char *foot;           /* DESI_TILES      */
  char *gcap;           /* GALACTIC_CAP    */
  int ncap;             /* number of galactic caps */
  char *fnz;            /* NZ_FILE         */
  double zmin;          /* ZMIN            */
  double zmax;          /* ZMAX            */
  int rng;              /* RAND_GENERATOR  */
  long *seed;           /* RAND_SEED       */
  char **output;        /* OUTPUT          */
  int ofmt;             /* OUTPUT_FORMAT   */
  int ovwrite;          /* OVERWRITE       */
  bool verbose;         /* VERBOSE         */
#ifdef OMP
  int nthread;
#endif
} CONF;


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv);

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf);

#endif

/*******************************************************************************
* define.h: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>
 
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

#ifndef __DEFINE_H__
#define __DEFINE_H__

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/
#define SPEED_OF_LIGHT  299792.458
#ifndef M_PI
#define M_PI            0x1.921fb54442d18p+1    /* PI */
#endif
#ifndef M_E
#define M_E             0x1.5bf0a8b145769p+1    /* e */
#endif

#define DEGREE_2_RAD    0x1.1df46a2529d39p-6    /* M_PI / 180 */
#define RAD_2_DEGREE    0x1.ca5dc1a63c1f8p+5    /* 180 / M_PI */
#define DOUBLE_EPSILON  1e-16   /* ~ machine epsilon for double numbers */
#define DOUBLE_TOL      1e-8    /* tolerance for double number comparison */

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
/* Default value for unset parameters. */
#define DEFAULT_CONF_FILE               "cutsky.conf"
#define DEFAULT_INPUT_FORMAT            CUTSKY_FFMT_ASCII
#define DEFAULT_OUTPUT_FORMAT           CUTSKY_FFMT_ASCII
#define DEFAULT_ASCII_COMMENT           '\0'
#define DEFAULT_NDATA                   (-1)
#define DEFAULT_DE_EOS_W                (-1)
#define DEFAULT_RNG                     PRAND_RNG_MT19937
#define DEFAULT_OVERWRITE               0
#define DEFAULT_VERBOSE                 true

/* Priority of parameters from different sources. */
#define CUTSKY_PRIOR_CMD                5
#define CUTSKY_PRIOR_FILE               1

/*============================================================================*\
                            Definitions for file IO
\*============================================================================*/
#define CUTSKY_PATH_SEP         '/'     /* separator for file paths         */
#define CUTSKY_FILE_CHUNK       4194304 /* chunk size for ASCII file IO     */
#define CUTSKY_MAX_CHUNK        INT_MAX /* maximum allowed chunk size       */
#define CUTSKY_MAX_LINES        INT_MAX /* maximum line number read at once */
#define CUTSKY_READ_COMMENT     '#'     /* comment symbol for reading       */
#define CUTSKY_SAVE_COMMENT     '#'     /* comment symbol for writing       */
#define CUTSKY_DATA_INIT_NUM    128     /* initial number of input data     */
#define CUTSKY_SPACE_ESCAPE     '\\'    /* escape character for spaces      */
#define CUTSKY_DATA_CHUNK       4096    /* number of data processed at once */

/* Enumeration of formats for input files. */
typedef enum {
  CUTSKY_FFMT_ASCII     = 0,
  CUTSKY_FFMT_FITS      = 1,
  CUTSKY_FFMT_FITS_LIST = 2
} CUTSKY_FFMT;

/*============================================================================*\
                            Other runtime constants
\*============================================================================*/
#define CUTSKY_CODE_NAME        "CUTSKY"        /* name of the program */

/* Parameters for coordinate conversion. */
#define CUTSKY_ZCNVT_DZ         1e-3    /* interval of redshift (z) samples */
#define CUTSKY_ZCNVT_MIN_NSP    100     /* minimum number of z samples      */
#define CUTSKY_ZCNVT_MAX_NSP    100000  /* maximum number of z samples      */
#define CUTSKY_ZCNVT_ORDER      10      /* order for Gauss integration of z */
#define CUTSKY_ZCNVT_EXT        10      /* number of bins extended on edges */
#define CUTSKY_ZCNVT_MAX_V      3000    /* maximum peculiar velocity        */

/* Parameters for survey geometry */
#define CUTSKY_BITCODE_INFOOT   1       /* code for inside the current foot */
#define CUTSKY_BITCODE_RAD_SEL  2       /* code for passing n(z) selection  */
#define CUTSKY_WMIN_FOOT_ALL    0       /* minimum weight for entire foot   */
#define CUTSKY_WMIN_FOOT        0       /* minimum weight for current foot  */
/* Right ascension range that distinguishes NGC and SGC. */
#define DESI_NGC_RA_MIN         90
#define DESI_NGC_RA_MAX         300
/* Rotation parameters for DESI volume. */
#define DESI_NGC_RA_SHIFT       60
#define DESI_SGC_RA_SHIFT       60

/* Threshold for throwing the warning of data number mismatch. */
#define CUTSKY_NDATA_MISMATCH   0.1

/*============================================================================*\
                     Definitions for the format of outputs
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m"               /* Red "Exit"        */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"  /* Red "FAIL"        */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL "%.10lg"       /* Output format for double parameters */
#define OFMT_FLT "%.8g"         /* Output format for float parameters  */

/*============================================================================*\
                          Definitions for error codes
\*============================================================================*/
#define CUTSKY_ERR_MEMORY       (-1)
#define CUTSKY_ERR_ARG          (-2)
#define CUTSKY_ERR_FILE         (-3)
#define CUTSKY_ERR_CFG          (-4)
#define CUTSKY_ERR_INIT         (-5)
#define CUTSKY_ERR_RAND         (-6)
#define CUTSKY_ERR_ZCVT         (-7)
#define CUTSKY_ERR_GEOM         (-8)
#define CUTSKY_ERR_CUTSKY       (-10)
#define CUTSKY_ERR_ASCII        (-11)
#define CUTSKY_ERR_SAVE         (-12)
#define CUTSKY_ERR_UNKNOWN      (-99)

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_EXIT " " __VA_ARGS__)

#endif


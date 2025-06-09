/*******************************************************************************
* load_conf.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "libcfg.h"
#include "prand.h"
#include "read_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
/* Check existence of configuration parameters. */
#define CHECK_EXIST_PARAM(name, cfg, var)                       \
  if (!cfg_is_set((cfg), (var))) {                              \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return CUTSKY_ERR_CFG;                                      \
  }
#define CHECK_EXIST_ARRAY(name, cfg, var, num)                  \
  if (!(num = cfg_get_size((cfg), (var)))) {                    \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return CUTSKY_ERR_CFG;                                      \
  }

/* Print the warning and error messages. */
#define P_CFG_WRN(cfg)  cfg_pwarn(cfg, stderr, FMT_WARN);
#define P_CFG_ERR(cfg)  {                                       \
  cfg_perror(cfg, stderr, FMT_ERR);                             \
  cfg_destroy(cfg);                                             \
  return NULL;                                                  \
}


/*============================================================================*\
                    Functions called via command line flags
\*============================================================================*/

/******************************************************************************
Function `usage`:
  Print the usage of command line options.
******************************************************************************/
static void usage(void *args) {
  char *pname = (char *) args;
  if (!pname || !(*pname)) pname = CUTSKY_CODE_NAME;
  else {
    /* Get the basename of the executable. */
    char *end = strrchr(pname, CUTSKY_PATH_SEP);
    if (end) pname = end + 1;
  }

  printf("Usage: %s [OPTION]\n\
Generate the cutsky catalogue from a cubic mock catalogue.\n\
  -h, --help\n\
        Display this message and exit\n\
  -t, --template\n\
        Print a template configuration file to the standard output and exit\n\
  -c, --conf            " FMT_KEY(CONFIG_FILE) "     String\n\
        Specify the configuration file (default: `%s')\n\
  -i, --input           " FMT_KEY(INPUT) "           String\n\
        Specify the filename of the input catalog\n\
  -f, --input-format    " FMT_KEY(INPUT_FORMAT) "    Integer\n\
        Specify the format of the input catalog\n\
      --comment         " FMT_KEY(COMMENT) "         Character\n\
        Specify the comment symbol for ASCII-format input catalog\n\
  -b, --box             " FMT_KEY(BOX_SIZE) "        Double\n\
        Set the side length of the cubic simulation box\n\
  -n, --number          " FMT_KEY(NUMBER) "          Long integer\n\
        Set the number of objects in the input catalog\n\
  -m, --omega-m         " FMT_KEY(OMEGA_M) "         Double\n\
        Set the density parameter of matter at z = 0\n\
      --omega-l         " FMT_KEY(OMEGA_LAMBDA) "    Double\n\
        Set the density parameter of Lambda at z = 0\n\
      --de-w            " FMT_KEY(DE_EOS_W) "        Double\n\
        Set the dark energy equation of state\n\
      --cmvdst-file     " FMT_KEY(Z_CMVDST_CNVT) "   String\n\
        Specify the file with redshift to comoving distance conversion table\n\
  -a, --foot-trim       " FMT_KEY(FOOTPRINT_TRIM) "  String\n\
        Set the Mangle polygon file for the footprint to be trimmed\n\
  -A, --foot-mark       " FMT_KEY(FOOTPRINT_MARK) "      String\n\
        Set the Mangle polygon file for the footprint to be marked\n\
  -C, --cap             " FMT_KEY(GALACTIC_CAP) "    Character array\n\
        Specify the galactic caps ('N' or 'S') to be produced\n\
  -N, --nz-file         " FMT_KEY(NZ_FILE) "         String\n\
        Specify the file for radial number density distribution\n\
  -z, --z-min           " FMT_KEY(ZMIN) "            Double\n\
        Set the minimum redshift of the output catalog\n\
  -Z, --z-max           " FMT_KEY(ZMAX) "            Double\n\
        Set the maximum redshift of the output catalog\n\
  -r, --rng             " FMT_KEY(RAND_GENERATOR) "  Integer\n\
        Specify the random number generation algorithm\n\
  -s, --seed            " FMT_KEY(RAND_SEED) "       Long integer array\n\
        Set seeds for random number generation in different galactic caps\n\
  -o, --output          " FMT_KEY(OUTPUT) "          String array\n\
        Specify the output catalogs for different galactic caps\n\
  -F, --output-format   " FMT_KEY(OUTPUT_FORMAT) "   Integer\n\
        Specify the format of output catalogs\n\
  -w, --overwrite       " FMT_KEY(OVERWRITE) "       Integer\n\
        Indicate whether to overwrite existing output files\n\
  -v, --verbose         " FMT_KEY(VERBOSE) "         Boolean\n\
        Indicate whether to display detailed standard outputs\n\
Consult the -t option for more information on the parameters.\n\
Github repository: https://github.com/cheng-zhao/cutsky.\n\
Licence: GPLv3.\n",
    pname, DEFAULT_CONF_FILE);
  exit(0);
}

/******************************************************************************
Function `conf_template`:
  Print a template configuration file.
******************************************************************************/
static void conf_template(void *args) {
  (void) args;
  printf("# Configuration file for cutsky (default: `%s').\n\
# Format: keyword = value # comment\n\
#     or: keyword = [element1, element2]\n\
#    see: https://github.com/cheng-zhao/libcfg for details.\n\
# For supported random number generation algorithms, see\n\
#         https://github.com/cheng-zhao/prand\n\
# NOTE that command line options have priority over this file.\n\
# Unnecessary entries can be left unset.\n\
\n\
######################################################\n\
#  Specifications of the input cubic simulation box  #\n\
######################################################\n\
\n\
INPUT           = \n\
    # String, filename of the input cubic simulation catalog.\n\
INPUT_FORMAT    = \n\
    # Integer, format of the input catalog (unset: %d). Allowed values are:\n\
    # * %d: ASCII file, with the leading 6 columns being (x,y,z,vx,vy,vz);\n\
    # * %d: FITS table, with case-insensitive column names of\n\
    #      ('x','y','z','vx','vy','vz');\n\
    # * %d: an ASCII file with lines being FITS filenames for subvolumes.\n\
COMMENT         = \n\
    # Character, indicate comments of ASCII-format `INPUT` (unset: '%c%s.\n\
    # Empty character ('') means disabling comments.\n\
BOX_SIZE        = \n\
    # Double-precision number, side length of the periodic box.\n\
NUMBER          = \n\
    # Long integer, number of objects in the input catalog.\n\
    # If unset, the number will be read from file.\n\
    # It is used only if radial selection is enabled (when `NZ_FILE` is set).\n\
    # Set this to the number of data objects when generating a random catalog.\n\
\n\n\
##################################################\n\
#  Fiducial cosmology for coordinate conversion  #\n\
##################################################\n\
\n\
OMEGA_M         = \n\
    # Double-precision number, density parameter of matter at z = 0.\n\
OMEGA_LAMBDA    = \n\
    # Double-precision number, density parameter of Lambda at z = 0.\n\
    # (unset: 1 - OMEGA_M).\n\
DE_EOS_W        = \n\
    # Double-precision number, dark energy equation of state: w (unset: "
  OFMT_DBL ").\n\
Z_CMVDST_CNVT   = \n\
    # Filename of a table for redshift to comoving distance conversion.\n\
    # It must be a text file with two columns: (redshift, comoving distance).\n\
    # If this file is set, the cosmological parameters above are omitted.\n\
    # Lines starting with '%c' are omitted.\n\
\n\n\
############################################\n\
#  Configurations for the survey geometry  #\n\
############################################\n\
\n\
FOOTPRINT_TRIM  = \n\
    # String, filename of the Mangle polygon file for the entire footprint.\n\
    # The output catalogs will be trimmed based on this footprint.\n\
FOOTPRINT_MARK  = \n\
    # String, filename of the polygon file for the footprint to be marked.\n\
    # If set, objects inside this footprint will be indicated by bitcode %d\n\
    # in the \"STATUS\" column.\n\
GALACTIC_CAP    = \n\
    # Character array, 'N' for northern galactic cap and 'S' for southern cap.\n\
NZ_FILE         = \n\
    # String, filename of the ASCII file for the radial selection function.\n\
    # The first two columns must be (redshift, comoving number density).\n\
    # If it is unset, then no radial selection is applied.\n\
    # Selected objects are indicated by bitcode %d in the \"STATUS\" column.\n\
    # Lines starting with '%c' are omitted.\n\
ZMIN            = \n\
ZMAX            = \n\
    # Double-precision numbers, minimum and maximum redshifts of the outputs.\n\
RAND_GENERATOR  = \n\
    # Random number generation algorithm (unset: %d).\n\
    # Integer, allowed values are:\n\
    # * 0: MRG32k3a\n\
    # * 1: Mersenne Twister 19937\n\
    # See https://github.com/cheng-zhao/prand for details.\n\
RAND_SEED       = \n\
    # Long integer array, random seeds for different galactic caps.\n\
    # Same dimension as `GALACTIC_CAP`.\n\
\n\n\
##############################\n\
#  Settings for the outputs  #\n\
##############################\n\
\n\
OUTPUT          = \n\
    # String array, name of the output catalogues.\n\
    # Same dimension as `GALACTIC_CAP`. One catalogue per galactic cap.\n\
OUTPUT_FORMAT   = \n\
    # Integer, format of the output catalog (unset: %d). Allowed values are:\n\
    # * %d: ASCII file;\n\
    # * %d: FITS table.\n\
OVERWRITE       = \n\
    # Integer, indicate whether to overwrite existing files (unset: %d).\n\
    # Allowed values are:\n\
    # * 0: quit the program when an output file exist;\n\
    # * positive: force overwriting output files whenever possible;\n\
    # * negative: notify at most this number of times for existing files.\n\
VERBOSE         = \n\
    # Boolean option, indicate whether to show detailed outputs (unset: %c).\n",
  DEFAULT_CONF_FILE, DEFAULT_INPUT_FORMAT, CUTSKY_FFMT_ASCII, CUTSKY_FFMT_FITS,
  CUTSKY_FFMT_FITS_LIST, DEFAULT_ASCII_COMMENT ? DEFAULT_ASCII_COMMENT : '\'',
  DEFAULT_ASCII_COMMENT ? "')" : ")", (double) DEFAULT_DE_EOS_W,
  CUTSKY_READ_COMMENT, CUTSKY_BITCODE_INFOOT, CUTSKY_BITCODE_RAD_SEL,
  CUTSKY_READ_COMMENT, DEFAULT_RNG, DEFAULT_OUTPUT_FORMAT,
  CUTSKY_FFMT_ASCII, CUTSKY_FFMT_FITS, DEFAULT_OVERWRITE,
  DEFAULT_VERBOSE ? 'T' : 'F');
  exit(0);
}


/*============================================================================*\
                      Function for reading configurations
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing configurations.
Return:
  Address of the structure.
******************************************************************************/
static CONF *conf_init(void) {
  CONF *conf = calloc(1, sizeof *conf);
  if (!conf) return NULL;
  conf->fconf = conf->input = conf->fzcnvt = conf->fnz = NULL;
  conf->foot_all = conf->foot = NULL;
  conf->gcap = NULL;
  conf->seed = NULL;
  conf->inputs = conf->output = NULL;
  return conf;
}

/******************************************************************************
Function `conf_read`:
  Read configurations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  Interface of libcfg.
******************************************************************************/
static cfg_t *conf_read(CONF *conf, const int argc, char *const *argv) {
  if (!conf) {
    P_ERR("the structure for configurations is not initialised\n");
    return NULL;
  }
  cfg_t *cfg = cfg_init();
  if (!cfg) P_CFG_ERR(cfg);

  /* Functions to be called via command line flags. */
  const cfg_func_t funcs[] = {
    {   'h',        "help",             usage,          NULL},
    {   't',    "template",     conf_template,          NULL}
  };

  /* Configuration parameters. */
  const cfg_param_t params[] = {
    {'c', "conf"         , "CONFIG_FILE"    , CFG_DTYPE_STR , &conf->fconf   },
    {'i', "input"        , "INPUT"          , CFG_DTYPE_STR , &conf->input   },
    {'f', "input-format" , "INPUT_FORMAT"   , CFG_DTYPE_INT , &conf->ifmt    },
    { 0 , "comment"      , "COMMENT"        , CFG_DTYPE_CHAR, &conf->comment },
    {'b', "box"          , "BOX_SIZE"       , CFG_DTYPE_DBL , &conf->Lbox    },
    {'n', "number"       , "NUMBER"         , CFG_DTYPE_LONG, &conf->ndata   },
    {'m', "omega-m"      , "OMEGA_M"        , CFG_DTYPE_DBL , &conf->omega_m },
    { 0 , "omega-l"      , "OMEGA_LAMBDA"   , CFG_DTYPE_DBL , &conf->omega_l },
    { 0 , "de-w"         , "DE_EOS_W"       , CFG_DTYPE_DBL , &conf->eos_w   },
    { 0 , "cmvdst-file"  , "Z_CMVDST_CNVT"  , CFG_DTYPE_STR , &conf->fzcnvt  },
    {'a', "foot-trim"    , "FOOTPRINT_TRIM" , CFG_DTYPE_STR , &conf->foot_all},
    {'A', "foot-mark"    , "FOOTPRINT_MARK" , CFG_DTYPE_STR , &conf->foot    },
    {'C', "cap"          , "GALACTIC_CAP"   , CFG_ARRAY_CHAR, &conf->gcap    },
    {'N', "nz-file"      , "NZ_FILE"        , CFG_DTYPE_STR , &conf->fnz     },
    {'z', "z-min"        , "ZMIN"           , CFG_DTYPE_DBL , &conf->zmin    },
    {'Z', "z-max"        , "ZMAX"           , CFG_DTYPE_DBL , &conf->zmax    },
    {'r', "rng"          , "RAND_GENERATOR" , CFG_DTYPE_INT , &conf->rng     },
    {'s', "seed"         , "RAND_SEED"      , CFG_ARRAY_LONG, &conf->seed    },
    {'o', "output"       , "OUTPUT"         , CFG_ARRAY_STR , &conf->output  },
    {'F', "output-format", "OUTPUT_FORMAT"  , CFG_DTYPE_INT , &conf->ofmt    },
    {'w', "overwrite"    , "OVERWRITE"      , CFG_DTYPE_INT , &conf->ovwrite },
    {'v', "verbose"      , "VERBOSE"        , CFG_DTYPE_BOOL, &conf->verbose }
  };

  /* Register functions and parameters. */
  if (cfg_set_funcs(cfg, funcs, sizeof(funcs) / sizeof(funcs[0])))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);
  if (cfg_set_params(cfg, params, sizeof(params) / sizeof(params[0])))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read configurations from command line options. */
  int optidx;
  if (cfg_read_opts(cfg, argc, argv, CUTSKY_PRIOR_CMD, &optidx)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read parameters from configuration file. */
  if (cfg_is_set(cfg, &conf->fconf)) {
    if (cfg_read_file(cfg, conf->fconf, CUTSKY_PRIOR_FILE)) P_CFG_ERR(cfg);
  }
  else {
    conf->fconf = NULL;
    if (cfg_read_file(cfg, DEFAULT_CONF_FILE, CUTSKY_PRIOR_FILE))
      P_CFG_ERR(cfg);
  }
  P_CFG_WRN(cfg);

  return cfg;
}


/*============================================================================*\
                      Functions for parameter verification
\*============================================================================*/

/******************************************************************************
Function `check_input`:
  Check whether an input file can be read.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_input(const char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR("the input " FMT_KEY(%s) " is not set\n", key);
    return CUTSKY_ERR_CFG;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot access " FMT_KEY(%s) ": `%s'\n", key, fname);
    return CUTSKY_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `check_output`:
  Check whether an output file can be written.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file;
  * `ovwrite`:  option for overwriting exisiting files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_output(char *fname, const char *key, const int ovwrite) {
  if (!fname || *fname == '\0') {
    P_ERR("the output " FMT_KEY(%s) " is not set\n", key);
    return CUTSKY_ERR_CFG;
  }

  /* Check if the file exists. */
  if (!access(fname, F_OK)) {
    /* not overwriting */
    if (ovwrite == 0) {
      P_ERR("the output " FMT_KEY(%s) " exists: `%s'\n", key, fname);
      return CUTSKY_ERR_FILE;
    }
    /* force overwriting */
    else if (ovwrite > 0) {
      P_WRN("the output " FMT_KEY(%s) " will be overwritten: `%s'\n",
          key, fname);
    }
    /* ask for decision */
    else {
      P_WRN("the output " FMT_KEY(%s) " exists: `%s'\n", key, fname);
      char confirm = 0;
      for (int i = 0; i != ovwrite; i--) {
        fprintf(stderr, "Are you going to overwrite it? (y/n): ");
        if (scanf("%c", &confirm) != 1) continue;
        int c;
        while((c = getchar()) != '\n' && c != EOF) continue;
        if (confirm == 'n') {
          P_ERR("cannot write to the file\n");
          return CUTSKY_ERR_FILE;
        }
        else if (confirm == 'y') break;
      }
      if (confirm != 'y') {
        P_ERR("too many failed inputs\n");
        return CUTSKY_ERR_FILE;
      }
    }

    /* Check file permission for overwriting. */
    if (access(fname, W_OK)) {
      P_ERR("cannot write to file `%s'\n", fname);
      return CUTSKY_ERR_FILE;
    }
  }
  /* Check the path permission. */
  else {
    char *end;
    if ((end = strrchr(fname, CUTSKY_PATH_SEP)) != NULL) {
      *end = '\0';
      if (access(fname, X_OK)) {
        P_ERR("cannot access the directory `%s'\n", fname);
        return CUTSKY_ERR_FILE;
      }
      *end = CUTSKY_PATH_SEP;
    }
  }
  return 0;
}

/******************************************************************************
Function `check_cosmo`:
  Verify cosmological parameters for coordinate conversion.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_cosmo(const cfg_t *cfg, CONF *conf) {
  /* Check OMEGA_M. */
  CHECK_EXIST_PARAM(OMEGA_M, cfg, &conf->omega_m);
  if (conf->omega_m <= 0 || conf->omega_m > 1) {
    P_ERR(FMT_KEY(OMEGA_M) " must be > 0 and <= 1\n");
    return CUTSKY_ERR_CFG;
  }

  /* Check OMEGA_LAMBDA */
  if (!cfg_is_set(cfg, &conf->omega_l)) {
    conf->omega_l = 1 - conf->omega_m;
    conf->omega_k = 0;
  }
  else if (conf->omega_l < 0) {
    P_ERR(FMT_KEY(OMEGA_LAMBDA) " must be >= 0\n");
    return CUTSKY_ERR_CFG;
  }
  else conf->omega_k = 1 - conf->omega_m - conf->omega_l;

  /* Check DE_EOS_W. */
  if (!cfg_is_set(cfg, &conf->eos_w)) conf->eos_w = DEFAULT_DE_EOS_W;
  else if (conf->eos_w > -1 / (double) 3) {
    P_ERR(FMT_KEY(DE_EOS_W) " must be <= -1/3\n");
    return CUTSKY_ERR_CFG;
  }

  /* Finally, make sure that H^2 (z) > 0. */
  double w3 = conf->eos_w * 3;
  double widx = w3 + 1;
  if (conf->omega_k * pow(conf->omega_l * (-widx), widx / w3) <=
      conf->omega_l * w3 * pow(conf->omega_m, widx / w3)) {
    P_ERR("negative H^2 given the cosmological parameters\n");
    return CUTSKY_ERR_CFG;
  }
  return 0;
}

/******************************************************************************
Function `conf_verify`:
  Verify configuration parameters.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int conf_verify(const cfg_t *cfg, CONF *conf) {
  int e;

  /* Check INPUT. */
  CHECK_EXIST_PARAM(INPUT, cfg, &conf->input);
  if ((e = check_input(conf->input, "INPUT"))) return e;

  /* Check INPUT_FORMAT. */
  if (!cfg_is_set(cfg, &conf->ifmt)) conf->ifmt = DEFAULT_INPUT_FORMAT;
  switch (conf->ifmt) {
    case CUTSKY_FFMT_ASCII:
      /* Check COMMENT. */
      if (!cfg_is_set(cfg, &conf->comment))
        conf->comment = DEFAULT_ASCII_COMMENT;
      break;
    case CUTSKY_FFMT_FITS:
#ifdef WITH_CFITSIO
      /* Allocate memory for the list anyway. */
      if (!(conf->inputs = malloc(sizeof(char *)))) {
        P_ERR("failed to allocate memory for the input catalog\n");
        return CUTSKY_ERR_MEMORY;
      }
      size_t len = strlen(conf->input) + 1;
      if (!(conf->inputs[0] = malloc(sizeof(char) * len))) {
        P_ERR("failed to allocate memory for the input catalog\n");
        return CUTSKY_ERR_MEMORY;
      }
      strncpy(conf->inputs[0], conf->input, len);
      conf->ninput = 1;
      break;
#else
      P_ERR("FITS format is not enabled\n"
          "Please re-compile the code with option -DWITH_CFITSIO");
      return CUTSKY_ERR_CFG;
#endif
    case CUTSKY_FFMT_FITS_LIST:
#ifdef WITH_CFITSIO
      /* Read filenames from list. */
      if (conf->ifmt == CUTSKY_FFMT_FITS_LIST) {
        if (read_filelist(conf->input, &conf->inputs, &conf->ninput))
          return CUTSKY_ERR_FILE;

        for (int i = 0; i < conf->ninput; i++) {
          if ((e = check_input(conf->inputs[i], "INPUT"))) return e;
        }
      }
      break;
#else
      P_ERR("FITS format is not enabled\n"
          "Please re-compile the code with option -DWITH_CFITSIO");
      return CUTSKY_ERR_CFG;
#endif
    default:
      P_ERR("invalid " FMT_KEY(INPUT_FORMAT) ": %d\n", conf->ifmt);
      return CUTSKY_ERR_CFG;
  }

  /* Check BOX_SIZE. */
  CHECK_EXIST_PARAM(BOX_SIZE, cfg, &conf->Lbox);
  if (conf->Lbox <= 0) {
    P_ERR(FMT_KEY(BOX_SIZE) " must be > 0\n");
    return CUTSKY_ERR_CFG;
  }

  /* Check the fidual cosmology */
  if (cfg_is_set(cfg, &conf->fzcnvt)) {
    if ((e = check_input(conf->fzcnvt, "Z_CMVDST_CNVT"))) return e;
  }
  else if ((e = check_cosmo(cfg, conf))) return e;

  /* Check FOOTPRINT_TRIM. */
  CHECK_EXIST_PARAM(FOOTPRINT_TRIM, cfg, &conf->foot_all);
  if ((e = check_input(conf->foot_all, "FOOTPRINT_TRIM"))) return e;

  /* Check FOOTPRINT_MASK. */
  if (cfg_is_set(cfg, &conf->foot)) {
    if ((e = check_input(conf->foot, "FOOTPRINT_MASK"))) return e;
  }

  /* Check GALACTIC_CAP. */
  CHECK_EXIST_ARRAY(GALACTIC_CAP, cfg, &conf->gcap, conf->ncap);
  if (conf->ncap > 2) {
    P_ERR("at most 2 elements for " FMT_KEY(GALACTIC_CAP) "\n");
    return CUTSKY_ERR_CFG;
  }
  for (int i = 0; i < conf->ncap; i++) {
    if (conf->gcap[i] != 'N' && conf->gcap[i] != 'S') {
      P_ERR(FMT_KEY(GALACTIC_CAP) " must be 'N' or 'S'\n");
      return CUTSKY_ERR_CFG;
    }
  }
  if (conf->ncap == 2 && conf->gcap[0] == conf->gcap[1]) {
    P_ERR("duplicate " FMT_KEY(GALACTIC_CAP) " element: '%c'\n", conf->gcap[0]);
    return CUTSKY_ERR_CFG;
  }

  /* Check NZ_FILE. */
  if (cfg_is_set(cfg, &conf->fnz)) {
    if ((e = check_input(conf->fnz, "NZ_FILE"))) return e;

    /* Check NUMBER. */
    if (!cfg_is_set(cfg, &conf->ndata)) conf->ndata = DEFAULT_NDATA;

    /* Check RAND_GENERATOR. */
    if (!cfg_is_set(cfg, &conf->rng)) conf->rng = DEFAULT_RNG;
    switch (conf->rng) {
      case PRAND_RNG_MT19937:
      case PRAND_RNG_MRG32K3A:
        break;
      default:
        P_ERR("invalid " FMT_KEY(RAND_GENERATOR) ": %d\n", conf->rng);
        return CUTSKY_ERR_CFG;
    }

    /* Check RAND_SEED. */
    int num = cfg_get_size(cfg, &conf->seed);
    if (num < conf->ncap) {
      P_ERR("too few elements of " FMT_KEY(RAND_SEED) "\n");
      return CUTSKY_ERR_CFG;
    }
    if (num > conf->ncap) {
      P_WRN("omitting the following " FMT_KEY(RAND_SEED) ":");
      for (int i = conf->ncap; i < num; i++)
        fprintf(stderr, " %ld", conf->seed[i]);
      fprintf(stderr, "\n");
    }
    for (int i = 0; i < conf->ncap; i++) {
      if (conf->seed <= 0) {
        P_ERR(FMT_KEY(RAND_SEED) " must be > 0\n");
        return CUTSKY_ERR_CFG;
      }
    }
  }
  else conf->fnz = NULL;

  /* Check ZMIN and ZMAX. */
  CHECK_EXIST_PARAM(ZMIN, cfg, &conf->zmin);
  CHECK_EXIST_PARAM(ZMAX, cfg, &conf->zmax);
  if (conf->zmin < 0 || conf->zmin >= conf->zmax) {
    P_ERR(FMT_KEY(ZMIN) " must be >= 0 and < " FMT_KEY(ZMAX) "\n");
    return CUTSKY_ERR_CFG;
  }

  /* Check OVERWRITE. */
  if (!cfg_is_set(cfg, &conf->ovwrite)) conf->ovwrite = DEFAULT_OVERWRITE;

  /* Check OUTPUT. */
  if (cfg_get_size(cfg, &conf->output) != conf->ncap) {
    P_ERR("Lengths of " FMT_KEY(GALACTIC_CAP) " and " FMT_KEY(OUTPUT)
        " must be equal\n");
    return CUTSKY_ERR_CFG;
  }
  for (int i = 0; i < conf->ncap; i++) {
    if ((e = check_output(conf->output[i], "OUTPUT", conf->ovwrite))) return e;
  }

  /* Check OUTPUT_FORMAT. */
  if (!cfg_is_set(cfg, &conf->ofmt)) conf->ofmt = DEFAULT_OUTPUT_FORMAT;
  switch (conf->ofmt) {
    case CUTSKY_FFMT_ASCII:
      break;
    case CUTSKY_FFMT_FITS:
#ifdef WITH_CFITSIO
      break;
#else
      P_ERR("FITS format is not enabled\n"
          "Please re-compile the code with option -DWITH_CFITSIO\n");
      return CUTSKY_ERR_CFG;
#endif
    default:
      P_ERR("invalid " FMT_KEY(OUTPUT_FORMAT) ": %d\n", conf->ofmt);
      return CUTSKY_ERR_CFG;
  }

  /* Check VERBOSE. */
  if (!cfg_is_set(cfg, &conf->verbose)) conf->verbose = DEFAULT_VERBOSE;

#ifdef OMP
  conf->nthread = omp_get_max_threads();
  if (conf->nthread <= 0) {
    P_ERR("invalid number of OpenMP threads\n");
    return CUTSKY_ERR_CFG;
  }
#endif

  return 0;
}


/*============================================================================*\
                      Function for printing configurations
\*============================================================================*/

/******************************************************************************
Function `conf_print`:
  Print configuration parameters.
Arguments:
  * `conf`:     structure for storing configurations.
******************************************************************************/
static void conf_print(const CONF *conf) {
  /* Configuration file */
  if (conf->fconf)
    printf("\n  CONFIG_FILE     = %s", conf->fconf);
  else
    printf("\n  CONFIG_FILE     = %s", DEFAULT_CONF_FILE);

  /* Input settings. */
  printf("\n  INPUT           = %s", conf->input);
  const char *fmt_name[3] = {"ASCII", "FITS", "FITS_LIST"};
  printf("\n  INPUT_FORMAT    = %d (%s)", conf->ifmt, fmt_name[conf->ifmt]);
  if (conf->ifmt == CUTSKY_FFMT_ASCII) {
    if (conf->comment == '\0') printf("\n  COMMENT         = ''");
    else printf("\n  COMMENT         = '%c'", conf->comment);
  }
  printf("\n  BOX_SIZE        = " OFMT_DBL, conf->Lbox);
  if (conf->fnz && conf->ndata != DEFAULT_NDATA)
    printf("\n  NUMBER          = %ld", conf->ndata);

  /* Fiducial cosmology. */
  if (conf->fzcnvt) printf("\n  Z_CMVDST_CNVT   = %s", conf->fzcnvt);
  else {
    printf("\n  OMEGA_M         = " OFMT_DBL, conf->omega_m);
    printf("\n  OMEGA_LAMBDA    = " OFMT_DBL, conf->omega_l);
    if (conf->eos_w != -1)
      printf("\n  DE_EOS_W        = " OFMT_DBL, conf->eos_w);
  }

  /* Survey geometry. */
  printf("\n  FOOTPRINT_TRIM  = %s", conf->foot_all);
  if (conf->foot) printf("\n  FOOTPRINT_MASK  = %s", conf->foot);
  if (conf->ncap == 1)
    printf("\n  GALACTIC_CAP    = %c", conf->gcap[0]);
  else
    printf("\n  GALACTIC_CAP    = [%c,%c]", conf->gcap[0], conf->gcap[1]);
  if (conf->fnz) printf("\n  NZ_FILE         = %s", conf->fnz);
  printf("\n  ZMIN            = " OFMT_DBL, conf->zmin);
  printf("\n  ZMAX            = " OFMT_DBL, conf->zmax);
  if (conf->fnz) {
    const char *rng_name[2] = {"MRG32K3A", "MT19937"};
    printf("\n  RAND_GENERATOR  = %d (%s)", conf->rng, rng_name[conf->rng]);
    if (conf->ncap == 1)
      printf("\n  RAND_SEED       = %ld", conf->seed[0]);
    else
      printf("\n  RAND_SEED       = [%ld,%ld]", conf->seed[0], conf->seed[1]);
  }

  /* Output. */
  printf("\n  OUTPUT          = %s", conf->output[0]);
  if (conf->ncap == 2) printf("\n                    %s", conf->output[1]);
  printf("\n  OUTPUT_FORMAT   = %d (%s)", conf->ofmt, fmt_name[conf->ofmt]);
  printf("\n  OVERWRITE       = %d\n", conf->ovwrite);
#ifdef OMP
  printf("  OMP_NUM_THREADS = %d\n", conf->nthread);
#endif
}


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
CONF *load_conf(const int argc, char *const *argv) {
  CONF *conf = conf_init();
  if (!conf) return NULL;

  cfg_t *cfg = conf_read(conf, argc, argv);
  if (!cfg) {
    conf_destroy(conf);
    return NULL;
  }

  printf("Loading configurations ...");
  fflush(stdout);

  if (conf_verify(cfg, conf)) {
    if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
    conf_destroy(conf);
    cfg_destroy(cfg);
    return NULL;
  }

  if (conf->verbose) conf_print(conf);

  if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
  cfg_destroy(cfg);

  printf(FMT_DONE);
  return conf;
}

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf) {
  if (!conf) return;
  if (conf->input) free(conf->input);
  if (conf->inputs) {
    if (*(conf->inputs)) free(*(conf->inputs));
    free(conf->inputs);
  }
  if (conf->fzcnvt) free(conf->fzcnvt);
  if (conf->fnz) free(conf->fnz);
  if (conf->foot_all) free(conf->foot_all);
  if (conf->foot) free(conf->foot);
  if (conf->gcap) free(conf->gcap);
  if (conf->seed) free(conf->seed);
  if (conf->output) {
    if (*(conf->output)) free(*(conf->output));
    free(conf->output);
  }
  free(conf);
}


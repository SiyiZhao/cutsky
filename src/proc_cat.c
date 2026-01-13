/*******************************************************************************
* proc_cat.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "proc_cat.h"
#include "read_file.h"
#include "write_file.h"
#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

#ifdef OMP

#include <omp.h>
#include <limits.h>

/* Data structure for maintaining the order of the data for reproductivity. */
typedef struct {
  size_t n;             /* number of data chunks                           */
  size_t max;           /* capacity of this structure                      */
  size_t *start;        /* starting index of each chunk in private catalog */
  size_t *length;       /* length of each chunk in private catalog         */
  size_t *iglobal;      /* global starting index of the current chunk      */
} DATA_CHUNK;

/* Shortcut for garbage collection. */
#define DATA_CLEAN_OMP                                                  \
  for (int ii = 0; ii < conf->ncap; ii++) {                             \
    for (int jj = 0; jj < conf->nthread; jj++) {                        \
      cutsky_destroy(pdata[ii][jj]); chunk_destroy(pchunk[ii][jj]);     \
    }                                                                   \
    free(pdata[ii]); free(pchunk[ii]);                                  \
  }

/*============================================================================*\
              Functions for processing the data chunk information
\*============================================================================*/

/******************************************************************************
Function `chunk_destroy`:
  Deconstruct the instance for storing data chunk information.
Arguments:
  * `chunk`:    instance of the data chunk structure.
******************************************************************************/
static void chunk_destroy(DATA_CHUNK *chunk) {
  if (!chunk) return;
  if (chunk->start) free(chunk->start);
  if (chunk->length) free(chunk->length);
  if (chunk->iglobal) free(chunk->iglobal);
  free(chunk);
}

/******************************************************************************
Function `chunk_init`:
  Initialise the structures for storing data chunk information.
Argument:
  * `nthread`:  number of threads to be run in parallel.
Return:
  Instances of the structures for each thread on success; NULL on error.
******************************************************************************/
static DATA_CHUNK **chunk_init(const int nthread) {
  DATA_CHUNK **chunk = malloc(nthread * sizeof(DATA_CHUNK *));
  if (!chunk) {
    P_ERR("failed to allocate memory for data chunks\n");
    return NULL;
  }
  for (int i = 0; i < nthread; i++) {
    if (!(chunk[i] = malloc(sizeof(DATA_CHUNK)))) {
      P_ERR("failed to allocate memory for data chunks\n");
      for (int j = 0; j < i; j++) chunk_destroy(chunk[j]);
      free(chunk);
      return NULL;
    }
    chunk[i]->start = chunk[i]->length = chunk[i]->iglobal = NULL;
    chunk[i]->n = 0;
    chunk[i]->max = CUTSKY_DATA_INIT_NUM;

    if (!(chunk[i]->start = malloc(chunk[i]->max * sizeof(size_t)))) {
      P_ERR("failed to allocate memory for data chunks\n");
      for (int j = 0; j < i; j++) chunk_destroy(chunk[j]);
      free(chunk);
      return NULL;
    }
  }

  return chunk;
}

/******************************************************************************
Function `chunk_append`:
  Append information of a chunk to the structure storing chunk information.
Arguments:
  * `chunk`:    interface of the data chunks;
  * `start`:    starting index of the current chunk.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int chunk_append(DATA_CHUNK *chunk, const size_t start) {
  /* Enlarge the memory if necessary. */
  if (chunk->n == chunk->max) {
    if (INT_MAX / 2 < chunk->max) {
      P_ERR("too many data chunks\n");
      return CUTSKY_ERR_CUTSKY;
    }
    chunk->max <<= 1;
    size_t *tmp;
    if (!(tmp = realloc(chunk->start, chunk->max * sizeof(size_t)))) {
      P_ERR("failed to allocate memory for data chunks\n");
      return CUTSKY_ERR_MEMORY;
    }
    chunk->start = tmp;
  }

  chunk->start[chunk->n++] = start;
  return 0;
}

#endif          /* OMP */

/*============================================================================*\
                 Functions for processing the cut-sky catalogue
\*============================================================================*/

/******************************************************************************
Function `cutsky_init`:
  Initialise the cut-sky catalogue.
Return:
  Instance of the cut-sky catalogue on success; NULL on error.
******************************************************************************/
static DATA *cutsky_init(void) {
  DATA *data = malloc(sizeof *data);
  if (!data) {
    P_ERR("failed to allocate memory for the cut-sky catalog\n");
    return NULL;
  }
  data->x[0] = data->x[1] = data->x[2] = data->x[3] = NULL;
  data->nz = data->ran = NULL;
  data->status = NULL;
  data->n = 0;
  data->max = CUTSKY_DATA_CHUNK;

  if (data->max) {
    for (int i = 0; i < 4; i++) {
      if (!(data->x[i] = malloc(data->max * sizeof(float)))) {
        P_ERR("failed to allocate memory for the cut-sky catalog\n");
        for (int j = 0; j < i; j++) free(data->x[j]);
        free(data);
        return NULL;
      }
    }
  }

  return data;
}

/******************************************************************************
Function `cutsky_destroy`:
  Deconstruct the cut-sky catalogue.
Arguments:
  * `data`:     instance of the cut-sky catalogue.
******************************************************************************/
static void cutsky_destroy(DATA *data) {
  if (!data) return;
  for (int i = 0; i < 4; i++) {
    if (data->x[i]) free(data->x[i]);
  }
  if (data->nz) free(data->nz);
  if (data->ran) free(data->ran);
  if (data->status) free(data->status);
  free(data);
}

/******************************************************************************
Function `cutsky_append`:
  Append a set of coordinates to the cut-sky catalogue.
Arguments:
  * `data`:     the cut-sky catalogue;
  * `ra`:       the right-acension of the tracer;
  * `dec`:      the declination of the tracer;
  * `z`:        the redshift-space redshift of the tracer;
  * `z_cosmo`:  the real-space redshift of the tracer.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cutsky_append(DATA *data, const float ra, const float dec,
    const float z, const float z_cosmo) {
  /* Enlarge the catalogue if necessary. */
  if (data->n == data->max) {
    if (SIZE_MAX / 2 < data->max) {
      P_ERR("too many objects in the cut-sky catalog\n");
      return CUTSKY_ERR_CUTSKY;
    }
    data->max <<= 1;
    for (int i = 0; i < 4; i++) {
      float *tmp = realloc(data->x[i], data->max * sizeof(float));
      if (!tmp) {
        P_ERR("failed to allocate memory for the cut-sky catalog\n");
        return CUTSKY_ERR_MEMORY;
      }
      data->x[i] = tmp;
    }
  }

  data->x[0][data->n] = ra;
  data->x[1][data->n] = dec;
  data->x[2][data->n] = z;
  data->x[3][data->n] = z_cosmo;
  data->n += 1;
  return 0;
}

/******************************************************************************
Function `cutsky_infoot`:
  Push objects passing the survey geometry test to the cut-sky catalogs.
Arguments:
  * `zcvt`:     interface for distance to redshift conversion;
  * `geom`:     interface for survey geometry;
  * `x`, `y`, `z`:      comoving coordinates;
  * `vx`, `vy`, `vz`:   peculiar velocities;
  * `ncap`:     number of galactic caps to be considered;
  * `ra_shift`: shift of right ascension for box rotation;
  * `is_ngc`:   indicate if the test is for ngc;
  * `data`:     cut-sky catalogs.
Return:
  True if the object is inside the footprint; false otherwise.
******************************************************************************/
static inline int cutsky_infoot(const ZCVT *zcvt, const GEOM *geom,
    const double x, const double y, const double z, const double vx,
    const double vy, const double vz, const int ncap, const double ra_shift[2],
    const bool is_ngc[2], DATA *data[2]) {
  /* Loops for box duplicates. */
  for (int i = -zcvt->ndup; i < zcvt->ndup; i++) {
    double xx = x + i * zcvt->Lbox;
    for (int j = -zcvt->ndup; j < zcvt->ndup; j++) {
      double yy = y + j * zcvt->Lbox;
      for (int k = -zcvt->ndup; k < zcvt->ndup; k++) {
        double zz = z + k * zcvt->Lbox;

        /* Compute the trim and radial distance with tolerance for RSD. */
        double d2 = xx * xx + yy * yy + zz * zz;
        if (d2 > zcvt->d2max || d2 < zcvt->d2min) continue;

        /* Compute the line-of-sight velocity. */
        double d_inv = 1 / sqrt(d2);
        double vel = (vx * xx + vy * yy + vz * zz) * d_inv;

        /* Convert squared distance to redshift. */
        double z_real = convert_z(zcvt, d2);
        double z_red = z_real + vel * (1 + z_real) / SPEED_OF_LIGHT;
        if (z_red < zcvt->zmin || z_red > zcvt->zmax) continue;

        /* Compute sky coordinates. */
        for (int n = 0; n < ncap; n++) {
          double ra, dec;
          if (d_inv > 1 / DOUBLE_TOL) ra = dec = 0;
          else {
            dec = asin(zz * d_inv) * RAD_2_DEGREE;
            ra = atan2(yy, xx) * RAD_2_DEGREE + ra_shift[n];
            if (ra < 0) ra += 360;
          }

          /* Pre-select NGC/SGC. */
          if ((ra > DESI_NGC_RA_MIN && ra < DESI_NGC_RA_MAX) != is_ngc[n])
            continue;

          /* Trim survey footprint. */
          if (geom_infoot(geom->foot[0], ra, dec)) {
            if (cutsky_append(data[n], ra, dec, z_red, z_real))
              return CUTSKY_ERR_CUTSKY;
          }
        }
      }
    }
  }

  return 0;
}

/******************************************************************************
Function `cutsky_save`:
  Save the cut-sky catalogue to an output file.
Arguments:
  * `fname`:    name of the output file;
  * `fmt`:      format of the output file;
  * `data`:     array of cut-sky catalogues to be saved;
  * `ncat`:     number of cut-sky catalogues.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cutsky_save(const char *fname, const CUTSKY_FFMT fmt,
    DATA **data, const int ncat) {
#ifdef WITH_CFITSIO
  if (fmt == CUTSKY_FFMT_ASCII) {        /* ASCII file */
#endif

    /* Open the output file for writing. */
    OFILE *ofile = output_init();
    if (!ofile || output_newfile(ofile, fname)) {
      output_destroy(ofile);
      return CUTSKY_ERR_FILE;
    }

    if (!data[0]->status && !data[0]->nz) {     /* no bitcode and nz */
      /* Header. */
      if (output_writeline(ofile, "%c RA(1) DEC(2) Z(3) Z_COSMO(4)\n",
          CUTSKY_SAVE_COMMENT)) return CUTSKY_ERR_FILE;

      for (int i = 0; i < ncat; i++) {
        for (size_t j = 0; j < data[i]->n; j++) {
          if (output_writeline(ofile, OFMT_FLT " " OFMT_FLT " " OFMT_FLT " "
              OFMT_FLT "\n", data[i]->x[0][j], data[i]->x[1][j],
              data[i]->x[2][j], data[i]->x[3][j]))
            return CUTSKY_ERR_FILE;
        }
      }
    }
    else if (!data[0]->nz) {                    /* bitcode only */
      /* Header. */
      if (output_writeline(ofile, "%c RA(1) DEC(2) Z(3) Z_COSMO(4) "
          "STATUS(5)\n", CUTSKY_SAVE_COMMENT)) return CUTSKY_ERR_FILE;

      for (int i = 0; i < ncat; i++) {
        for (size_t j = 0; j < data[i]->n; j++) {
          if (output_writeline(ofile, OFMT_FLT " " OFMT_FLT " " OFMT_FLT " "
              OFMT_FLT " %" PRId8 "\n", data[i]->x[0][j], data[i]->x[1][j],
              data[i]->x[2][j], data[i]->x[3][j], data[i]->status[j]))
            return CUTSKY_ERR_FILE;
        }
      }
    }
    else {                                      /* both bitcode and nz */
      /* Header. */
      if (output_writeline(ofile, "%c RA(1) DEC(2) Z(3) Z_COSMO(4) "
          "NZ(5) STATUS(6) RAN_NUM_0_1(7)\n", CUTSKY_SAVE_COMMENT))
        return CUTSKY_ERR_FILE;

      for (int i = 0; i < ncat; i++) {
        for (size_t j = 0; j < data[i]->n; j++) {
          if (output_writeline(ofile, OFMT_FLT " " OFMT_FLT " " OFMT_FLT " "
              OFMT_FLT " " OFMT_FLT " %" PRId8 " " OFMT_FLT "\n",
              data[i]->x[0][j], data[i]->x[1][j], data[i]->x[2][j],
              data[i]->x[3][j], data[i]->nz[j], data[i]->status[j],
              data[i]->ran[j]))
            return CUTSKY_ERR_FILE;
        }
      }
    }

    /* Close file. */
    output_destroy(ofile);

#ifdef WITH_CFITSIO
  }
  else {                                        /* FITS file */

    /* Allocate memory for the filename. */
    const size_t len = strlen(fname) + 1;
    char *fitsname = malloc(len + 1);
    if (!fitsname) {
      P_ERR("failed to allocate memory for catalog writing\n");
      return CUTSKY_ERR_MEMORY;
    }
    fitsname[0] = '!';
    strncpy(fitsname + 1, fname, len);

    /* Open the output file for writing. */
    OFFILE *ofile = ofits_init();
    if (!ofile) {
      free(fitsname);
      return CUTSKY_ERR_FILE;
    }

    /* Setup columns. */
    if (!data[0]->status && !data[0]->nz) {     /* no bitcode and nz */
      const int ncol = 4;
      char *names[] = {"RA", "DEC", "Z", "Z_COSMO"};
      char *units[] = {"deg", "deg", NULL, NULL};
      int dtypes[] = {TFLOAT, TFLOAT, TFLOAT, TFLOAT};

      if (ofits_newfile(ofile, fitsname, ncol, names, units, dtypes)) {
        free(fitsname);
        return CUTSKY_ERR_FILE;
      }

      for (int i = 0; i < ncat; i++) {
        for (size_t j = 0; j < data[i]->n; j++) {
          if (ofits_writeline(ofile, data[i]->x[0][j], data[i]->x[1][j],
              data[i]->x[2][j], data[i]->x[3][j])) {
            free(fitsname);
            return CUTSKY_ERR_FILE;
          }
        }
      }
    }
    else if (!data[0]->nz) {                    /* bitcode only */
      const int ncol = 5;
      char *names[] = {"RA", "DEC", "Z", "Z_COSMO", "STATUS"};
      char *units[] = {"deg", "deg", NULL, NULL, NULL};
      int dtypes[] = {TFLOAT, TFLOAT, TFLOAT, TFLOAT, TBYTE};

      if (ofits_newfile(ofile, fitsname, ncol, names, units, dtypes)) {
        free(fitsname);
        return CUTSKY_ERR_FILE;
      }

      for (int i = 0; i < ncat; i++) {
        for (size_t j = 0; j < data[i]->n; j++) {
          if (ofits_writeline(ofile, data[i]->x[0][j], data[i]->x[1][j],
              data[i]->x[2][j], data[i]->x[3][j], data[i]->status[j])) {
            free(fitsname);
            return CUTSKY_ERR_FILE;
          }
        }
      }
    }
    else {                                      /* both bitcode and nz */
      const int ncol = 7;
      char *names[] = {
          "RA", "DEC", "Z", "Z_COSMO", "NZ", "STATUS", "RAN_NUM_0_1"};
      char *units[] = {"deg", "deg", NULL, NULL, NULL, NULL, NULL};
      int dtypes[] = {TFLOAT, TFLOAT, TFLOAT, TFLOAT, TFLOAT, TBYTE, TFLOAT};

      if (ofits_newfile(ofile, fitsname, ncol, names, units, dtypes)) {
        free(fitsname);
        return CUTSKY_ERR_FILE;
      }

      for (int i = 0; i < ncat; i++) {
        for (size_t j = 0; j < data[i]->n; j++) {
          if (ofits_writeline(ofile, data[i]->x[0][j], data[i]->x[1][j],
              data[i]->x[2][j], data[i]->x[3][j], data[i]->nz[j],
              data[i]->status[j], data[i]->ran[j])) {
            free(fitsname);
            return CUTSKY_ERR_FILE;
          }
        }
      }
    }

    free(fitsname);
    /* Close file. */
    ofits_destroy(ofile);

  }
#endif

  return 0;
}

/******************************************************************************
Function `process_serial`:
  Process the input data catalogue serially.
Arguments:
  * `conf`:     structure for storing configurations;
  * `zcvt`:     interface for redshift conversion;
  * `geom`:     interface for survey geometry.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int process_serial(const CONF *conf, const ZCVT *zcvt, const GEOM *geom)
{
  /* Process NGC and SGC individually. */
  DATA *data[2] = {NULL, NULL};
  const double ra_shift[2] = {DESI_NGC_RA_SHIFT, DESI_SGC_RA_SHIFT};
  bool is_ngc[2] = {false, false};

  for (int i = 0; i < conf->ncap; i++) {
    is_ngc[i] = (conf->gcap[i] == 'N');
    if (!(data[i] = cutsky_init())) return CUTSKY_ERR_CUTSKY;
  }

  size_t nline = CUTSKY_DATA_CHUNK;
  size_t nbox = 0;      /* number of objects in the simulation box */

#ifdef WITH_CFITSIO
  if (conf->ifmt == CUTSKY_FFMT_ASCII) {        /* ASCII file */
#endif

    /* Open the file for reading. */
    IFILE *ifile = input_init();
    if (!ifile || input_newfile(ifile, conf->input)) {
      cutsky_destroy(data[0]); cutsky_destroy(data[1]); input_destroy(ifile);
      return CUTSKY_ERR_FILE;
    }

    /* Read the input file by chunk. */
    for (;;) {
      if (input_readlines(ifile, nline)) {
        cutsky_destroy(data[0]); cutsky_destroy(data[1]); input_destroy(ifile);
        return CUTSKY_ERR_FILE;
      }
      if (!ifile->nline) break;

      for (size_t i = 0; i < ifile->nline; i++) {
        char *line = ifile->chunk + ifile->lines[i];
        if (!line) {
          P_ERR("failed to read line from the input catalog\n");
          cutsky_destroy(data[0]); cutsky_destroy(data[1]);
          input_destroy(ifile);
          return CUTSKY_ERR_FILE;
        }

        while (isspace(*line)) ++line;  /* omit leading whitespaces */
        if (*line == conf->comment || *line == '\0') continue;

        /* Parse the line. */
        double x, y, z, vx, vy, vz;
        int ncol = sscanf(line, "%lf %lf %lf %lf %lf %lf",
            &x, &y, &z, &vx, &vy, &vz);
        if (ncol != 6) {
          if (ncol == 3) {
            vx = vy = vz = 0;
          }
          else {
            P_ERR("failed to read data from line: %s\n", line);
            cutsky_destroy(data[0]); cutsky_destroy(data[1]);
            input_destroy(ifile);
            return CUTSKY_ERR_FILE;
          }
        }
        nbox += 1;

        /* Apply coordinate conversion and survey geometry. */
        if (cutsky_infoot(zcvt, geom, x, y, z, vx, vy, vz, conf->ncap,
            ra_shift, is_ngc, data)) {
          cutsky_destroy(data[0]); cutsky_destroy(data[1]);
          input_destroy(ifile);
          return CUTSKY_ERR_CUTSKY;
        }
      }
    }

    /* Close the input file. */
    input_destroy(ifile);

#ifdef WITH_CFITSIO
  }
  else {                                        /* FITS file(s) */

    /* Open the file for reading. */
    IFFILE *ifile = ifits_init();
    if (!ifile ||
        ifits_newfiles(ifile, (const char **) conf->inputs, conf->ninput)) {
      cutsky_destroy(data[0]); cutsky_destroy(data[1]); ifits_destroy(ifile);
      return CUTSKY_ERR_FILE;
    }

    /* Read the input file by chunk. */
    for (;;) {
      if (ifits_readlines(ifile, nline)) {
        cutsky_destroy(data[0]); cutsky_destroy(data[1]); ifits_destroy(ifile);
        return CUTSKY_ERR_FILE;
      }
      if (!ifile->ndata) break;
      nbox += ifile->ndata;

      /* Apply coordinate conversion and survey geometry. */
      for (size_t i = 0; i < ifile->ndata; i++) {
        if (cutsky_infoot(zcvt, geom, ifile->data[0][i], ifile->data[1][i],
            ifile->data[2][i], ifile->data[3][i], ifile->data[4][i],
            ifile->data[5][i], conf->ncap, ra_shift, is_ngc, data)) {
          cutsky_destroy(data[0]); cutsky_destroy(data[1]);
          ifits_destroy(ifile);
          return CUTSKY_ERR_CUTSKY;
        }
      }
    }

    /* Close the input file. */
    ifits_destroy(ifile);

  }
#endif

  if (!nbox) {
    P_ERR("no data in the input catalog\n");
    cutsky_destroy(data[0]); cutsky_destroy(data[1]);
    return CUTSKY_ERR_FILE;
  }

  if (conf->fnz && conf->ndata != DEFAULT_NDATA) {
    if (fabs((double) conf->ndata - nbox) >
        (double) nbox * CUTSKY_NDATA_MISMATCH) {
      P_WRN("mismatched number of objects: %ld from the configuration, "
          "while %zu in the catalog\n"
          "Using the number in the configuration anyway\n",
          conf->ndata, nbox);
    }
    nbox = conf->ndata;
  }
  if (conf->verbose)
    printf("  %zu objects read from the input catalog\n", nbox);

  /* Reduce memory cost if applicable. */
  for (int i = 0; i < conf->ncap; i++) {
    if (!data[i]->n) {
      P_WRN("no data after footprint trimming for %cGC\n", conf->gcap[i]);
      continue;
    }
    if (data[i]->n < data[i]->max) {
      for (int j = 0; j < 4; j++) {
        float *tmp = realloc(data[i]->x[j], data[i]->n * sizeof(float));
        if (tmp) data[i]->x[j] = tmp;
      }
    }
  }

  if (conf->fnz) {      /* add columns for radial selection */
    for (int i = 0; i < conf->ncap; i++) {
      if (!data[i]->n) continue;

      /* Allocate memory. */
      if (!(data[i]->nz = malloc(data[i]->n * sizeof(float))) ||
          !(data[i]->ran = malloc(data[i]->n * sizeof(float))) ||
          !(data[i]->status = calloc(data[i]->n, sizeof(uint8_t)))) {
        P_ERR("failed to allocate memory for the cut-sky catalog\n");
        cutsky_destroy(data[0]); cutsky_destroy(data[1]);
        return CUTSKY_ERR_MEMORY;
      }

      /* Compute the comoving number density. */
      double dens_sim = nbox / pow(conf->Lbox, 3);

      /* Reset random number generator. */
      int err = 0;
      prand_t *rng = geom->rng;
      rng->reset(rng->state, geom->seed[i], 0, &err);
      if (PRAND_IS_ERROR(err) || PRAND_IS_WARN(err)) {
        P_ERR("failed to set the random number generator\n");
        cutsky_destroy(data[0]); cutsky_destroy(data[1]);
        return CUTSKY_ERR_RAND;
      }

      /* Apply radial selection. */
      for (size_t j = 0; j < data[i]->n; j++) {
        data[i]->nz[j] = geom_get_nz(geom, data[i]->x[2][j]);
        data[i]->ran[j] = rng->get_double(rng->state);
        double prop = data[i]->nz[j] / dens_sim;
        if (data[i]->ran[j] < prop) data[i]->status[j] = geom->rad_sel;

        /* Apply the current footprint of interest. */
        if (conf->foot &&
            geom_infoot(geom->foot[1], data[i]->x[0][j], data[i]->x[1][j]))
          data[i]->status[j] += geom->infoot;
      }
    }
  }     /* if (conf->fnz) */
  else if (conf->foot) {        /* add bitcodes for footprint */
    for (int i = 0; i < conf->ncap; i++) {
      if (!data[i]->n) continue;

      /* Allocate memory. */
      if (!(data[i]->status = calloc(data[i]->n, sizeof(uint8_t)))) {
        P_ERR("failed to allocate memory for the cut-sky catalog\n");
        cutsky_destroy(data[0]); cutsky_destroy(data[1]);
        return CUTSKY_ERR_MEMORY;
      }

      /* Check footprint. */
      for (size_t j = 0; j < data[i]->n; j++) {
        if (geom_infoot(geom->foot[1], data[i]->x[0][j], data[i]->x[1][j]))
          data[i]->status[j] = geom->infoot;
      }
    }
  }

  /* Save the catalogues. */
  for (int i = 0; i < conf->ncap; i++) {
    if (!data[i]->n) continue;
    if (cutsky_save(conf->output[i], conf->ofmt, &(data[i]), 1)) {
      for (int j = i; j < conf->ncap; j++) cutsky_destroy(data[j]);
      return CUTSKY_ERR_FILE;
    }

    if (conf->verbose) printf("  %zu objects saved to the output for %cGC\n",
        data[i]->n, conf->gcap[i]);
    cutsky_destroy(data[i]);
  }
  return 0;
}

#ifdef OMP

/******************************************************************************
Function `process_omp`:
  Process the input data catalogue with OpenMP parallelisation.
Arguments:
  * `conf`:     structure for storing configurations;
  * `zcvt`:     interface for redshift conversion;
  * `geom`:     interface for survey geometry.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int process_omp(const CONF *conf, const ZCVT *zcvt, const GEOM *geom) {
  /* Allocate memory for NGC and SGC. */
  DATA **pdata[2] = {NULL, NULL};
  DATA_CHUNK **pchunk[2] = {NULL, NULL};
  const double ra_shift[2] = {DESI_NGC_RA_SHIFT, DESI_SGC_RA_SHIFT};
  bool is_ngc[2] = {false, false};

  for (int i = 0; i < conf->ncap; i++) {
    is_ngc[i] = (conf->gcap[i] == 'N');
    if (!(pdata[i] = malloc(conf->nthread * sizeof(DATA *)))) {
      P_ERR("failed to allocate memory for the cut-sky catalog\n");
      for (int ii = 0; ii < i; ii++) {
        for (int j = 0; j < conf->nthread; j++) cutsky_destroy(pdata[ii][j]);
        free(pdata[ii]);
      }
      return CUTSKY_ERR_MEMORY;
    }

    for (int j = 0; j < conf->nthread; j++) {
      if (!(pdata[i][j]= cutsky_init())) {
        P_ERR("failed to callocate memory for the cut-sky catalog\n");
        for (int ii = 0; ii < i; ii++) {
          for (int jj = 0; jj < conf->nthread; jj++)
            cutsky_destroy(pdata[ii][jj]);
          free(pdata[ii]);
        }
        for (int jj = 0; jj < j; jj++) cutsky_destroy(pdata[i][jj]);
        free(pdata[i]);
        return CUTSKY_ERR_MEMORY;
      }
    }
  }

  for (int i = 0; i < conf->ncap; i++) {
    if (!(pchunk[i] = chunk_init(conf->nthread))) {
      for (int ii = 0; ii < conf->ncap; ii++) {
        for (int j = 0; j < conf->nthread; j++) cutsky_destroy(pdata[ii][j]);
        free(pdata[ii]);
      }
      for (int ii = 0; ii < i; ii++) {
        for (int j = 0; j < conf->nthread; j++) chunk_destroy(pchunk[ii][j]);
        free(pchunk[ii]);
      }
      return CUTSKY_ERR_MEMORY;
    }
  }

  /* Number of lines to be read at once. */
  size_t nline = (size_t) conf->nthread * CUTSKY_DATA_CHUNK;
  size_t nbox = 0;      /* number of objects in the simulation box */

#ifdef WITH_CFITSIO
  if (conf->ifmt == CUTSKY_FFMT_ASCII) {        /* ASCII file */
#endif

    /* Open the file for reading. */
    IFILE *ifile = input_init();
    if (!ifile || input_newfile(ifile, conf->input)) {
      DATA_CLEAN_OMP;
      return CUTSKY_ERR_FILE;
    }

    /* Read the input file by chunk. */
    for (;;) {
      if (input_readlines(ifile, nline)) {
        DATA_CLEAN_OMP; input_destroy(ifile);
        return CUTSKY_ERR_FILE;
      }
      if (!ifile->nline) break;         /* reading completed */

      /* Distribute lines to OpenMP threads. */
      const size_t pnum = ifile->nline / conf->nthread;
      const int rem = ifile->nline % conf->nthread;

#pragma omp parallel num_threads(conf->nthread)
      {
        const int tid = omp_get_thread_num();
        const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
        const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
        const size_t iend = istart + pcnt;
        DATA *data[2] = {pdata[0][tid], NULL};
        if (conf->ncap == 2) data[1] = pdata[1][tid];

        /* Save the starting index of the chunk in the cut-sky catalog. */
        for (int i = 0; i < conf->ncap; i++) {
          if (chunk_append(pchunk[i][tid], pdata[i][tid]->n)) {
            DATA_CLEAN_OMP; input_destroy(ifile);
            exit(CUTSKY_ERR_FILE);
          }
        }

        size_t pnbox = 0;
        for (size_t i = istart; i < iend; i++) {
          char *line = ifile->chunk + ifile->lines[i];
          if (!line) {
            P_ERR("failed to read line from the input catalog\n");
            DATA_CLEAN_OMP; input_destroy(ifile);
            exit(CUTSKY_ERR_FILE);
          }

          while (isspace(*line)) ++line;        /* omit leading whitespaces */
          if (*line == conf->comment || *line == '\0') continue;

          /* Parse the line. */
          double x, y, z, vx, vy, vz;
          int ncol = sscanf(line, "%lf %lf %lf %lf %lf %lf",
              &x, &y, &z, &vx, &vy, &vz);
          if (ncol != 6) {
            if (ncol == 3) {
              vx = vy = vz = 0;
            }
            else {
              P_ERR("failed to read data from line: %s\n", line);
              DATA_CLEAN_OMP; input_destroy(ifile);
              exit(CUTSKY_ERR_FILE);
            }
          }
          pnbox += 1;

          /* Apply coordinate conversion and survey geometry. */
          if (cutsky_infoot(zcvt, geom, x, y, z, vx, vy, vz, conf->ncap,
                            ra_shift, is_ngc, data)) {
            DATA_CLEAN_OMP; input_destroy(ifile);
            exit(CUTSKY_ERR_CUTSKY);
          }
          pdata[0][tid] = data[0];
          if (conf->ncap == 2) pdata[1][tid] = data[1];
        }

#pragma omp critical
        nbox += pnbox;
      } /* omp parallel */
    }

    /* Close the input file. */
    input_destroy(ifile);

#ifdef WITH_CFITSIO
  }
  else {                                        /* FITS file */

    /* Open the file for reading. */
    IFFILE *ifile = ifits_init();
    if (!ifile ||
        ifits_newfiles(ifile, (const char **) conf->inputs, conf->ninput)) {
      DATA_CLEAN_OMP;
      return CUTSKY_ERR_FILE;
    }

    /* Read the input file by chunk. */
    for (;;) {
      if (ifits_readlines(ifile, nline)) {
        DATA_CLEAN_OMP; ifits_destroy(ifile);
        return CUTSKY_ERR_FILE;
      }
      if (!ifile->ndata) break;         /* reading completed */
      nbox += ifile->ndata;

      /* Distribute lines to OpenMP threads. */
      const size_t pnum = ifile->ndata / conf->nthread;
      const int rem = ifile->ndata % conf->nthread;

#pragma omp parallel num_threads(conf->nthread)
      {
        const int tid = omp_get_thread_num();
        const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
        const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
        const size_t iend = istart + pcnt;
        DATA *data[2] = {pdata[0][tid], NULL};
        if (conf->ncap == 2) data[1] = pdata[1][tid];

        /* Save the starting index of the chunk in the cut-sky catalog. */
        for (int i = 0; i < conf->ncap; i++) {
          if (chunk_append(pchunk[i][tid], pdata[i][tid]->n)) {
            DATA_CLEAN_OMP; ifits_destroy(ifile);
            exit(CUTSKY_ERR_FILE);
          }
        }

        /* Apply coordinate conversion and survey geometry. */
        for (size_t i = istart; i < iend; i++) {
          if (cutsky_infoot(zcvt, geom, ifile->data[0][i], ifile->data[1][i],
              ifile->data[2][i], ifile->data[3][i], ifile->data[4][i],
              ifile->data[5][i], conf->ncap, ra_shift, is_ngc, data)) {
            DATA_CLEAN_OMP; ifits_destroy(ifile);
            exit(CUTSKY_ERR_CUTSKY);
          }
        }
        pdata[0][tid] = data[0];
        if (conf->ncap == 2) pdata[1][tid] = data[1];
      } /* omp parallel */
    }

    /* Close the input file. */
    ifits_destroy(ifile);

  }
#endif

  if (!nbox) {
    P_ERR("no data in the input catalog\n");
    DATA_CLEAN_OMP;
    return CUTSKY_ERR_FILE;
  }
  if (conf->verbose)
    printf("  %zu objects read from the input catalog\n", nbox);

  if (conf->fnz && conf->ndata != DEFAULT_NDATA) {
    if (fabs((double) conf->ndata - nbox) >
        (double) nbox * CUTSKY_NDATA_MISMATCH) {
      P_WRN("mismatched number of objects: %ld from the configuration, "
          "while %zu in the catalog\n"
          "Using the number in the configuration anyway\n",
          conf->ndata, nbox);
    }
    nbox = conf->ndata;
  }

#pragma omp parallel num_threads(conf->nthread)
  {
    const int tid = omp_get_thread_num();

    for (int i = 0; i < conf->ncap; i++) {
      DATA *data = pdata[i][tid];
      if (!data->n) continue;

      /* Reduce memory cost if applicable. */
      if (data->n < data->max) {
        for (int k = 0; k < 4; k++) {
          float *tmp = realloc(data->x[k], data->n * sizeof(float));
          if (tmp) data->x[k] = tmp;
        }
      }

      if (conf->fnz) {          /* apply radial selection */
        DATA_CHUNK *chunk = pchunk[i][tid];
        if (!chunk->n) {
          P_ERR("unexpected empty data chunk\n");
          DATA_CLEAN_OMP;
          exit(CUTSKY_ERR_UNKNOWN);
        }

        if (chunk->n < chunk->max) {
          size_t *tmp = realloc(chunk->start, chunk->n * sizeof(size_t));
          if (tmp) chunk->start = tmp;
        }

        /* Allocate memory for the cut-sky catalog and data chunks. */
        if (!(data->nz = malloc(data->n * sizeof(float))) ||
            !(data->ran = malloc(data->n * sizeof(float))) ||
            !(data->status = calloc(data->n, sizeof(uint8_t))) ||
            !(chunk->length = malloc(chunk->n * sizeof(size_t))) ||
            !(chunk->iglobal = malloc(chunk->n * sizeof(size_t)))) {
          P_ERR("failed to allocate memory for the cut-sky catalog\n");
          DATA_CLEAN_OMP;
          exit(CUTSKY_ERR_MEMORY);
        }

        /* Compute the lengths of different chunks. */
        for (int k = 1; k < chunk->n; k++)
          chunk->length[k - 1] = chunk->start[k] - chunk->start[k - 1];
        chunk->length[chunk->n - 1] = data->n - chunk->start[chunk->n - 1];

#pragma omp barrier
        /* Compute the global starting indices of different chunks. */
#pragma omp single
        {
          DATA_CHUNK **c = pchunk[i];
          c[0]->iglobal[0] = 0;
          for (int j = 1; j < conf->nthread; j++)
            c[j]->iglobal[0] = c[j - 1]->iglobal[0] + c[j - 1]->length[0];

          const int last = conf->nthread - 1;
          const int nlast = c[last]->n;
          for (int k = 1; k < nlast; k++) {
            c[0]->iglobal[k] = c[last]->iglobal[k - 1] + c[last]->length[k - 1];
            for (int j = 1; j < conf->nthread; j++)
              c[j]->iglobal[k] = c[j - 1]->iglobal[k] + c[j - 1]->length[k];
          }

          if (c[0]->n > nlast) {
            c[0]->iglobal[c[0]->n - 1] =
                c[last]->iglobal[nlast] + c[last]->length[nlast];
            for (int j = 1; j < conf->nthread; j++) {
              if (c[j]->n == nlast) break;
              c[j]->iglobal[nlast] =
                  c[j - 1]->iglobal[nlast] + c[j - 1]->length[nlast];
            }
          }
        }
#pragma omp barrier

        /* Compute the comoving number density. */
        double dens_sim = nbox / pow(conf->Lbox, 3);

        for (int k = 0; k < chunk->n; k++) {
          /* Reset random number generator. */
          int err = 0;
          prand_t *rng = geom->rng;
          rng->reset(rng->state_stream[tid], geom->seed[i], chunk->iglobal[k],
              &err);
          if (PRAND_IS_ERROR(err) || PRAND_IS_WARN(err)) {
            P_ERR("failed to set the random number generator\n");
            DATA_CLEAN_OMP;
            exit(CUTSKY_ERR_RAND);
          }

          for (size_t n = chunk->start[k];
              n < chunk->start[k] + chunk->length[k]; n++) {
            /* Apply radial selection. */
            data->nz[n] = geom_get_nz(geom, data->x[2][n]);
            data->ran[n] = rng->get_double(rng->state_stream[tid]);
            double prop = data->nz[n] / dens_sim;
            if (data->ran[n] < prop) data->status[n] = geom->rad_sel;

            /* Apply the current footprint of interest. */
            if (conf->foot &&
                geom_infoot(geom->foot[1], data->x[0][n], data->x[1][n]))
              data->status[n] += geom->infoot;
          }
        }
      }         /* if (conf->fnz) */
      else if (conf->foot) {            /* add bitcodes for footprint */
        /* Allocate memory. */
        if (!(data->status = calloc(data->n, sizeof(uint8_t)))) {
          P_ERR("failed to allocate memory for the cut-sky catalog\n");
          DATA_CLEAN_OMP;
          exit(CUTSKY_ERR_MEMORY);
        }

        /* Check footprint. */
        for (size_t n = 0; n < data->n; n++) {
          if (geom_infoot(geom->foot[1], data->x[0][n], data->x[1][n]))
            data->status[n] = geom->infoot;
        }
      }
    }
  }     /* omp parallel */

  for (int i = 0; i < conf->ncap; i++) {
    for (int j = 0; j < conf->nthread; j++) chunk_destroy(pchunk[i][j]);
    free(pchunk[i]);
  }

  /* Save the catalogues. */
  for (int i = 0; i < conf->ncap; i++) {
    if (cutsky_save(conf->output[i], conf->ofmt, pdata[i], conf->nthread)) {
      for (int ii = i; ii < conf->ncap; ii++) {
        for (int jj = 0; jj < conf->nthread; jj++)
          cutsky_destroy(pdata[ii][jj]);
        free(pdata[ii]);
      }
      return CUTSKY_ERR_FILE;
    }

    if (conf->verbose) {
      size_t ndata = 0;
      for (int j = 0; j < conf->nthread; j++) ndata += pdata[i][j]->n;
      printf("  %zu objects saved to the output for %cGC\n",
          ndata, conf->gcap[i]);
    }

    for (int j = 0; j < conf->nthread; j++) cutsky_destroy(pdata[i][j]);
    free(pdata[i]);
  }

  return 0;
}

#endif          /* OMP */

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
int process(const CONF *conf, const ZCVT *zcvt, const GEOM *geom) {
  printf("Running the cutsky process ...");
  if (!conf) {
    P_ERR("configurations are not loaded\n");
    return CUTSKY_ERR_ARG;
  }
  if (!zcvt) {
    P_ERR("redshift to distance interpolation is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (!geom) {
    P_ERR("survey geometry is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

#ifdef OMP
  if (conf->nthread > 1) {
    if (process_omp(conf, zcvt, geom)) return CUTSKY_ERR_CUTSKY;
  }
  else {
#endif
    if (process_serial(conf, zcvt, geom)) return CUTSKY_ERR_CUTSKY;
#ifdef OMP
  }
#endif

  printf(FMT_DONE);
  return 0;
}

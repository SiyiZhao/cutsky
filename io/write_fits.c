/*******************************************************************************
* write_fits.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#ifdef WITH_CFITSIO

#include "define.h"
#include "write_file.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>

#define FITS_ERROR                                                      \
  P_ERR("cfitsio error: ");                                             \
  fits_report_error(stderr, status);

/*============================================================================*\
                        Interfaces for FITS file writing
\*============================================================================*/

/******************************************************************************
Function `ofits_init`:
  Initialise the interface for FITS file writing.
Return:
  Address of the interface.
******************************************************************************/
OFFILE *ofits_init(void) {
  OFFILE *ofile = malloc(sizeof *ofile);
  if (!ofile) {
    P_ERR("failed to initialize the interface for file writing\n");
    return NULL;
  }

  ofile->fp = NULL;
  ofile->data = ofile->nulval = NULL;
  ofile->cnum = NULL;
  ofile->dtypes = NULL;
  ofile->ncol = 0;
  ofile->nrow = 0;

  return ofile;
}

/******************************************************************************
Function `ofits_destroy`:
  Deconstruct the interface for writing FITS files.
Arguments:
  * `ofile`:    interface for FITS file writing.
******************************************************************************/
void ofits_destroy(OFFILE *ofile) {
  if (!ofile) return;
  if (ofits_flush(ofile))
    P_WRN("closing the file with unsaved buffer\n");
  if (ofile->fp) {
    int status = 0;
    if (fits_close_file(ofile->fp, &status)) {
      P_WRN("failed to close FITS file: ");
      fits_report_error(stderr, status);
      fits_clear_errmsg();
    }
  }
  if (ofile->data) {
    for (int i = 0; i < ofile->ncol; i++) {
      if (ofile->data[i]) free(ofile->data[i]);
    }
    free(ofile->data);
  }
  if (ofile->nulval) free(ofile->nulval);
  if (ofile->cnum) free(ofile->cnum);
  free(ofile);
}

/******************************************************************************
Function `ofits_newfile`:
  Flush the buffer to the existing FITS file and open a new file.
Arguments:
  * `ofile`:    interface for FITS file writing;
  * `fname`:    name of the file to be written to;
  * `ncol`:     number of columns to be written;
  * `names`:    names of the columns;
  * `units`:    units of the columns;
  * `dtypes`:   data types of the columns.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ofits_newfile(OFFILE *ofile, const char *fname, const int ncol,
    char **names, char **units, int *dtypes) {
  /* Validate arguments. */
  if (!fname || !(*fname)) {
    P_ERR("invalid output file name\n");
    return CUTSKY_ERR_ARG;
  }
  if (ncol <= 0) {
    P_ERR("non-positive number of FITS columns: %d\n", ncol);
    return CUTSKY_ERR_ARG;
  }
  if (!names || !units || !dtypes) {
    P_ERR("invalid column settings\n");
    return CUTSKY_ERR_ARG;
  }
  for (int i = 0; i < ncol; i++) {
    if (!names[i] || !(*names[i])) {
      P_ERR("invalid column names\n");
    return CUTSKY_ERR_ARG;
    }
    switch (dtypes[i]) {
      case TFLOAT:
      case TDOUBLE:
      case TBYTE:
      case TUSHORT:
      case TINT32BIT:
      case TULONGLONG:
        break;
      default:
        P_ERR("unsupported data type for FITS table: %d\n", dtypes[i]);
        return CUTSKY_ERR_ARG;
    }
  }

  /* Flush the buffer and close the previously opened FITS file. */
  int status = 0;
  if (ofile->fp) {
    if (ofits_flush(ofile)) {
      P_ERR("failed to flush unsaved buffer to the opened file\n");
      return CUTSKY_ERR_FILE;
    }

    if (fits_close_file(ofile->fp, &status)) {
      P_WRN("failed to close file: ");
      fits_report_error(stderr, status);
      status = 0;
      fits_clear_errmsg();
    }
    ofile->fp = NULL;

    /* Release memory for columns. */
    if (ofile->data) {
      for (int i = 0; i < ofile->ncol; i++) {
        if (ofile->data[i]) free(ofile->data[i]);
      }
      free(ofile->data);
      ofile->data = NULL;
    }
    if (ofile->nulval) {
      free(ofile->nulval);
      ofile->nulval = NULL;
    }
    if (ofile->cnum) {
      free(ofile->cnum);
      ofile->cnum = NULL;
    }
  }

  /* Create the new FITS file. */
  if (fits_create_file(&ofile->fp, fname, &status)) {
    FITS_ERROR;
    return CUTSKY_ERR_FILE;
  }

  /* Create an empty table. */
  char **cfmt = malloc(ncol * sizeof(char *));
  if (!cfmt) {
    P_ERR("failed to allocate memory for FITS columns\n");
    return CUTSKY_ERR_MEMORY;
  }

  for (int i = 0; i < ncol; i++) {
    switch (dtypes[i]) {
      case TFLOAT:
        cfmt[i] = "1E";
        break;
      case TDOUBLE:
        cfmt[i] = "1D";
        break;
      case TBYTE:
        cfmt[i] = "1B";
        break;
      case TUSHORT:
        cfmt[i] = "1I";
        break;
      case TINT32BIT:
        cfmt[i] = "1J";
        break;
      case TULONGLONG:
        cfmt[i] = "1K";
        break;
    }
  }

  if (fits_create_tbl(ofile->fp, BINARY_TBL, 0, ncol, names, cfmt, units, NULL,
      &status)) {
    FITS_ERROR;
    free(cfmt);
    return CUTSKY_ERR_FILE;
  }
  free(cfmt);

  /* Get the optimal chunk size for FITS file writing. */
  if (fits_get_rowsize(ofile->fp, &ofile->nrow, &status)) {
    FITS_ERROR;
    return CUTSKY_ERR_FILE;
  }
  if (ofile->nrow <= 0) {
    P_ERR("failed to determine the chunk size for FITS file writing\n");
    return CUTSKY_ERR_FILE;
  }

  /* Allocate memory for the chunks. */
  void **data = malloc(ncol * sizeof(void *));
  if (!data) {
    P_ERR("failed to allocate memory for FITS file writing\n");
    return CUTSKY_ERR_MEMORY;
  }
  for (int i = 0; i < ncol; i++) data[i] = NULL;
  for (int i = 0; i < ncol; i++) {
    switch (dtypes[i]) {
      case TFLOAT:
        data[i] = malloc(ofile->nrow * sizeof(float));
        break;
      case TDOUBLE:
        data[i] = malloc(ofile->nrow * sizeof(double));
        break;
      case TBYTE:
        data[i] = malloc(ofile->nrow * sizeof(uint8_t));
        break;
      case TUSHORT:
        data[i] = malloc(ofile->nrow * sizeof(uint16_t));
        break;
      case TINT32BIT:
        data[i] = malloc(ofile->nrow * sizeof(uint32_t));
        break;
      case TULONGLONG:
        data[i] = malloc(ofile->nrow * sizeof(uint64_t));
        break;
    }

    if (!data[i]) {
      P_ERR("failed to allocate memory for FITS file writing\n");
      for (int j = 0; j < i; j++) {
        if (data[j]) free(data[j]);
      }
      free(data);
      return CUTSKY_ERR_MEMORY;
    }
  }

  ofile->data = data;
  ofile->ncol = ncol;
  ofile->dtypes = dtypes;
  ofile->ndata = 0;

  if (!(ofile->nulval = malloc(ncol * sizeof(void *))) ||
      !(ofile->cnum = malloc(ncol * sizeof(int)))) {
    P_ERR("failed to allocate memory for FITS file writing\n");
    return CUTSKY_ERR_MEMORY;
  }
  for (int i = 0; i < ncol; i++) {
    ofile->nulval[i] = NULL;
    ofile->cnum[i] = i + 1;
  }

  return 0;
}

/******************************************************************************
Function `ofits_flush`:
  Write the buffer string to file.
Arguments:
  * `ofile`:    interface for FITS file writing.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ofits_flush(OFFILE *ofile) {
  if (!ofile) {
    P_ERR("the interface for FITS file writing is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (!ofile->fp || !ofile->data || !ofile->nulval || !ofile->dtypes ||
      ofile->ncol <= 0 || ofile->nrow <= 0) {
    P_ERR("the interface for FITS file writing is not initialized correctly\n");
    return CUTSKY_ERR_ARG;
  }

  if (!ofile->ndata) return 0;

  int status = 0;
  if (fits_write_cols(ofile->fp, ofile->ncol, ofile->dtypes, ofile->cnum,
      1, ofile->ndata, ofile->data, ofile->nulval, &status)) {
    FITS_ERROR;
    return CUTSKY_ERR_FILE;
  }

  ofile->ndata = 0;
  return 0;
}

/******************************************************************************
Function `ofits_writeline`:
  Write a line to the buffer and save it to the FITS file if necessary.
Arguments:
  * `ofile`:    interface for FITS file writing;
  * `...`:      arguments specifying data to be saved.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ofits_writeline(OFFILE *ofile, ...) {
  if (!ofile) {
    P_ERR("the interface for FITS file writing is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (!ofile->fp || !ofile->data || !ofile->nulval || !ofile->dtypes ||
      ofile->ncol <= 0 || ofile->nrow <= 0) {
    P_ERR("the interface for FITS file writing is not initialized correctly\n");
    return CUTSKY_ERR_ARG;
  }

  if (ofile->ndata == ofile->nrow) {
    if (ofits_flush(ofile)) return CUTSKY_ERR_FILE;
  }

  va_list args;
  va_start(args, ofile);

  for (int i = 0; i < ofile->ncol; i++) {
    switch (ofile->dtypes[i]) {
      case TFLOAT:
        ((float *) (ofile->data[i]))[ofile->ndata] = va_arg(args, double);
        break;
      case TDOUBLE:
        ((double *) (ofile->data[i]))[ofile->ndata] = va_arg(args, double);
        break;
      case TBYTE:
        ((uint8_t *) (ofile->data[i]))[ofile->ndata] = va_arg(args, uint32_t);
        break;
      case TUSHORT:
        ((uint16_t *) (ofile->data[i]))[ofile->ndata] = va_arg(args, uint32_t);
        break;
      case TINT32BIT:
        ((uint32_t *) (ofile->data[i]))[ofile->ndata] = va_arg(args, uint32_t);
        break;
      case TULONGLONG:
        ((uint64_t *) (ofile->data[i]))[ofile->ndata] = va_arg(args, uint64_t);
        break;
      default:
        P_ERR("unsupported data type for FITS table: %d\n", ofile->dtypes[i]);
        va_end(args);
        return CUTSKY_ERR_ARG;
    }
  }

  ofile->ndata += 1;
  va_end(args);
  return 0;
}

#endif

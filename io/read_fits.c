/*******************************************************************************
* read_fits.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#ifdef WITH_CFITSIO

#include "define.h"
#include "read_file.h"
#include <stdlib.h>
#include <string.h>

#define FITS_ABORT {                                                    \
  P_ERR("cfitsio error: ");                                             \
  fits_report_error(stderr, status);                                    \
  return CUTSKY_ERR_FILE;                                               \
}

/*============================================================================*\
                        Functions for FITS file reading
\*============================================================================*/

/******************************************************************************
Function `ifits_newfile`:
  Open a new file for reading.
Arguments:
  * `ifile`:    interface for file reading;
  * `fname`:    name of the file to be read from.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int ifits_newfile(IFFILE *ifile, const char *fname) {
  int status = 0;

  /* Close the previous file if needed. */
  if (ifile->fp && fits_close_file(ifile->fp, &status)) {
    P_WRN("failed to close file: ");
    fits_report_error(stderr, status);
    status = 0;
    fits_clear_errmsg();
  }
  ifile->fp = NULL;
  ifile->ntotal = ifile->nread = ifile->nstep = 0;

  /* Open the new file. */
  if (fits_open_data(&ifile->fp, fname, READONLY, &status)) FITS_ABORT;

  /* Check if the coordinate and velocity columns are in the file. */
  char *cname[6] = {"x", "y", "z", "vx", "vy", "vz"};
  for (int i = 0; i < 6; i++) {
    if (fits_get_colnum(ifile->fp, CASEINSEN, cname[i], ifile->col + i,
        &status)) FITS_ABORT;
  }

  /* Get the number of objects. */
  if (fits_get_num_rows(ifile->fp, &ifile->ntotal, &status)) FITS_ABORT;
  if (fits_get_rowsize(ifile->fp, &ifile->nstep, &status)) FITS_ABORT;

  if (!ifile->ntotal) {
    P_ERR("no data in the file: `%s'\n", fname);
    return CUTSKY_ERR_FILE;
  }
  if (!ifile->nstep) {
    P_ERR("failed to read data from file: `%s'\n", fname);
    return CUTSKY_ERR_FILE;
  }

  return 0;
}

/******************************************************************************
Function `readchunk_fits`:
  Read a chunk of data from a FITS file.
Arguments:
  * `ifile`:    interface for file reading;
  * `nline`:    number of lines to be reported.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int readchunk_fits(IFFILE *ifile, const long nline) {
  long nread = 0;
  const long ngoal = nline - ifile->nsave;

  while (nread < ngoal) {
    long nrest = ifile->ntotal - ifile->nread;  /* remaining objects in file */
    long nrow = (ifile->nstep < nrest) ? ifile->nstep : nrest;
    long nsave = ifile->nsave + nrow;

    /* Enlarge the chunk size if needed. */
    while (ifile->maxdata < nsave) {
      if (CUTSKY_MAX_CHUNK / 2 < ifile->maxdata) {
        P_ERR("too many objects to be read at once: %ld\n", nsave);
        return CUTSKY_ERR_FILE;
      }
      const long size = ifile->maxdata * 2;
      for (int i = 0; i < 6; i++) {
        double *tmp = realloc(ifile->data[i], size * sizeof(double));
        if (!tmp) {
          P_ERR("failed to allocate memory for read file by chunk\n");
          return CUTSKY_ERR_MEMORY;
        }
        ifile->data[i] = tmp;
      }
      ifile->maxdata = size;
    }

    for (int i = 0; i < 6; i++) {
      int status = 0;
      int anynul = 0;
      if (nrow && fits_read_col_dbl(ifile->fp, ifile->col[i], ifile->nread + 1,
          1, nrow, 0, ifile->data[i] + ifile->nsave, &anynul, &status))
        FITS_ABORT;
    }

    nread += nrow;
    ifile->nsave += nrow;
    ifile->nread += nrow;

    /* Open a new file if necessary. */
    if (ifile->nread >= ifile->ntotal) {
      if (ifile->current >= ifile->ninput - 1) {        /* no more file */
        ifile->ndata = nread;
        return 0;
      }

      ifile->current += 1;
      if (ifits_newfile(ifile, ifile->fnames[ifile->current]))
        return CUTSKY_ERR_FILE;
    }
  }

  ifile->ndata = nline;
  return 0;
}

/*============================================================================*\
                        Interfaces for file reading
\*============================================================================*/

/******************************************************************************
Function `ifits_init`:
  Initialise the interface for file reading.
Return:
  Address of the interface.
******************************************************************************/
IFFILE *ifits_init(void) {
  IFFILE *ifile = calloc(1, sizeof *ifile);
  if (!ifile) {
    P_ERR("failed to initialize the interface for file reading\n");
    return NULL;
  }

  ifile->fnames = NULL;
  ifile->fp = NULL;
  for (int i = 0; i < 6; i++) ifile->data[i] = NULL;

  return ifile;
}

/******************************************************************************
Function `ifits_destroy`:
  Deconstruct the interface for file reading.
Arguments:
  * `ifile`:    structure for file reading.
******************************************************************************/
void ifits_destroy(IFFILE *ifile) {
  if (!ifile) return;
  if (ifile->fp) {
    int status = 0;
    if (fits_close_file(ifile->fp, &status)) {
      P_WRN("failed to close file: ");
      fits_report_error(stderr, status);
      fits_clear_errmsg();
    }
  }
  free(ifile);
}

/******************************************************************************
Function `ifits_newfiles`:
  Open a set of new files for reading.
Arguments:
  * `ifile`:    interface for file reading;
  * `fnames`:   names of the file to be read from;
  * `ninput`:   number of input files to be read from.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ifits_newfiles(IFFILE *ifile, const char **fnames, const int ninput) {
  /* Validate arguments. */
  if (!ifile) {
    P_ERR("the interface for file reading is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (ninput <= 0) {
    P_ERR("invalid number of input files: %d\n", ninput);
    return CUTSKY_ERR_ARG;
  }
  if (!fnames) {
    P_ERR("invalid input file names\n");
    return CUTSKY_ERR_ARG;
  }
  for (int i = 0; i < ninput; i++) {
    if (!fnames[i] || !(*fnames[i])) {
      P_ERR("invalid input file names\n");
      return CUTSKY_ERR_ARG;
    }
  }

  /* Open the first file. */
  if (ifits_newfile(ifile, fnames[0])) return CUTSKY_ERR_FILE;

  ifile->fnames = fnames;
  ifile->ninput = ninput;
  ifile->current = 0;
  return 0;
}

/******************************************************************************
Function `ifits_readlines`:
  Read multiple records (lines) from the input file.
Arguments:
  * `ifile`:    interface for file reading;
  * `nline`:    number of lines to be read at once.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ifits_readlines(IFFILE *ifile, const size_t nline) {
  if (!ifile) {
    P_ERR("the interface for file reading is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  long nl = (long) nline;
  if (nl <= 0) {
    P_ERR("number of lines read from file should be positive\n");
    return CUTSKY_ERR_ARG;
  }

  /* Allocate memory if needed. */
  if (ifile->maxdata < nl) {
    ifile->maxdata = nl;
    for (int i = 0; i < 6; i++) {
      double *tmp = realloc(ifile->data[i], ifile->maxdata * sizeof(double));
      if (!tmp) {
        P_ERR("failed to allocate memory for the data\n");
        return CUTSKY_ERR_MEMORY;
      }
      ifile->data[i] = tmp;
    }
  }

  /* Remove previously reported data. */
  if (ifile->ndata) {
    if (ifile->nsave < ifile->ndata) {
      P_ERR("unexpected number of objects: %ld cached, but %ld reported\n",
          ifile->nsave, ifile->ndata);
      return CUTSKY_ERR_FILE;
    }
    ifile->nsave -= ifile->ndata;
    for (int i = 0; i < 6; i++)
      memmove(ifile->data[i], ifile->data[i] + ifile->ndata,
          ifile->nsave * sizeof(double));
  }

  /* Early return if there are already enough objects in the buffer. */
  if (nl <= ifile->nsave) {
    ifile->ndata = nl;
    return 0;
  }

  /* Check if all files are read completely. */
  if (ifile->current >= ifile->ninput - 1 && ifile->nread >= ifile->ntotal) {
    ifile->ndata = ifile->nsave;
    return 0;
  }

  /* Read file. */
  return readchunk_fits(ifile, nl);
}

#endif

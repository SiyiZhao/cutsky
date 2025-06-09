/*******************************************************************************
* read_data.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_data.h"
#include "read_file.h"
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>

/*============================================================================*\
                        Functions for reading data files
\*============================================================================*/

/******************************************************************************
Function `read_ascii_twocol`:
  Read the first two columns of an ASCII file as double arrays.
Arguments:
  * `fname`:    filename of the input catalog;
  * `x`:        array for the first column;
  * `y`:        array for the second column;
  * `num`:      number of lines read successfully.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_twocol(const char *fname, double **x, double **y, size_t *num) {
  if (!fname || !(*fname)) {
    P_ERR("invalid name of input file\n");
    return CUTSKY_ERR_ARG;
  }
  if (!x || !y || !num) {
    P_ERR("variables for storing the input data are not initialized\n");
    return CUTSKY_ERR_ARG;
  }

  /* Open the input file for reading. */
  IFILE *ifile = input_init();
  if (!ifile) return CUTSKY_ERR_FILE;
  if (input_newfile(ifile, fname)) {
    input_destroy(ifile);
    return CUTSKY_ERR_FILE;
  }

  /* Allocate memory for the data. */
  size_t nmax = CUTSKY_DATA_INIT_NUM;
  double *nx = malloc(nmax * sizeof(double));
  double *ny = malloc(nmax * sizeof(double));
  if (!nx || !ny) {
    P_ERR("failed to allocate memory for reading file: `%s'\n", fname);
    input_destroy(ifile);
    if (nx) free(nx);
    if (ny) free(ny);
    return CUTSKY_ERR_MEMORY;
  }

  size_t n = 0;
  for (;;) {
    if (input_readline(ifile)) {
      input_destroy(ifile); free(nx); free(ny);
      return CUTSKY_ERR_FILE;
    }

    const char *p = ifile->line;
    if (!p) break;

    while (isspace(*p)) ++p;    /* omit leading whitespaces */
    if (*p == CUTSKY_READ_COMMENT || *p == '\0') continue;

    /* Parse the line. */
    if (sscanf(p, "%lf %lf", nx + n, ny + n) != 2) {
      P_ERR("failed to parse line: %s\n", p);
      input_destroy(ifile); free(nx); free(ny);
      return CUTSKY_ERR_FILE;
    }

    /* Enlarge memory for the data if necessary. */
    if (++n >= nmax) {
      if (SIZE_MAX / 2 < nmax) {
        P_ERR("too many records in file: `%s'\n", fname);
        input_destroy(ifile); free(nx); free(ny);
        return CUTSKY_ERR_FILE;
      }

      nmax <<= 1;
      double *tmp = realloc(nx, nmax * sizeof(double));
      if (!tmp) {
        P_ERR("failed to allocate memory for the input data\n");
        input_destroy(ifile); free(nx); free(ny);
        return CUTSKY_ERR_MEMORY;
      }
      nx = tmp;
      if (!(tmp = realloc(ny, nmax * sizeof(double)))) {
        P_ERR("failed to allocate memory for the input data\n");
        input_destroy(ifile); free(nx); free(ny);
        return CUTSKY_ERR_MEMORY;
      }
      ny = tmp;
    }
  }

  input_destroy(ifile);

  if (!n) {
    P_ERR("no data is read from file: `%s'\n", fname);
    free(nx); free(ny);
    return CUTSKY_ERR_FILE;
  }

  /* Reduce the memory usage if possible. */
  if (n != nmax) {
    double *tmp = realloc(nx, n * sizeof(double));
    if (tmp) nx = tmp;
    tmp = realloc(ny, n * sizeof(double));
    if (tmp) ny = tmp;
  }

  *x = nx;
  *y = ny;
  *num = n;
  return 0;
}

/******************************************************************************
Function `read_filelist`:
  Read filenames from a list.
Arguments:
  * `fname`:    filename of the input list;
  * `list`:     list of filenames read from file;
  * `num`:      number of lines read successfully.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_filelist(const char *fname, char ***list, int *num) {
  if (!fname || !(*fname)) {
    P_ERR("invalid name of input file\n");
    return CUTSKY_ERR_ARG;
  }
  if (!list || !num) {
    P_ERR("variables for storing the input data are not initialized\n");
    return CUTSKY_ERR_ARG;
  }

  IFILE *ifile = input_init();
  if (!ifile) return CUTSKY_ERR_FILE;
  if (input_newfile(ifile, fname)) {
    input_destroy(ifile);
    return CUTSKY_ERR_FILE;
  }

  /* Allocate memory for the filenames. */
  size_t nmax = CUTSKY_DATA_INIT_NUM;
  size_t imax = CUTSKY_DATA_INIT_NUM;
  char *names = NULL;           /* a single array for all filenames */
  size_t *idx = NULL;           /* indices of filenames in the array */

  if (!(names = malloc(nmax * sizeof(char))) ||
      !(idx = malloc(imax * sizeof(size_t)))) {
    P_ERR("failed to allocate memory for filenames\n");
    input_destroy(ifile);
    if (names) free(names);
    return CUTSKY_ERR_MEMORY;
  }

  int n = 0;
  size_t size = 0;
  for (;;) {
    if (input_readline(ifile)) {
      input_destroy(ifile); free(names); free(idx);
      return CUTSKY_ERR_FILE;
    }

    const char *p = ifile->line;
    if (!p) break;

    while (isspace(*p)) ++p;    /* omit leading whitespaces */
    if (*p == CUTSKY_READ_COMMENT || *p == '\0') continue;

    /* Parse the line. */
    idx[n] = size;
    do {
      if (isspace (*p)) {
        /* Overwrite the escape character if applicable. */
        if (*(p - 1) == CUTSKY_SPACE_ESCAPE) size--;
        else break;
      }
      names[size] = *p;
      /* Enlarge memory for the filenames if necessary. */
      if (++size >= nmax - 1) {         /* reserve one byte for '\0' */
        if (SIZE_MAX / 2 < nmax) {
          P_ERR("too many characters in file: `%s'\n", fname);
          input_destroy(ifile); free(names); free(idx);
          return CUTSKY_ERR_MEMORY;
        }
        nmax <<= 1;
        char *tmp = realloc(names, nmax * sizeof(char));
        if (!tmp) {
          P_ERR("failed to allocate memory for filenames\n");
          input_destroy(ifile); free(names); free(idx);
          return CUTSKY_ERR_MEMORY;
        }
        names = tmp;
      }
    }
    while (*(++p) != '\0');
    names[size++] = '\0';               /* null termination */

    /* Enlarge memory for the indices if necessary. */
    if (++n >= imax) {
      if (INT_MAX / 2 < imax) {
        P_ERR("too many entries in file: `%s'\n", fname);
        input_destroy(ifile); free(names); free(idx);
        return CUTSKY_ERR_MEMORY;
      }
      imax <<= 1;
      size_t *tmp = realloc(idx, imax * sizeof(size_t));
      if (!tmp) {
        P_ERR("failed to allocate memory for filename indices\n");
        input_destroy(ifile); free(names); free(idx);
        return CUTSKY_ERR_MEMORY;
      }
      idx = tmp;
    }
  }

  input_destroy(ifile);

  if (!n || !size) {
    P_ERR("no valid filename found in file: `%s'\n", fname);
    free(names); free(idx);
    return CUTSKY_ERR_FILE;
  }

  /* Reduce the memory usage if applicable. */
  if (size != nmax) {
    char *tmp = realloc(names, size * sizeof(char));
    if (tmp) names = tmp;
  }

  /* Allocate memory for the filename array. */
  if (!(*list = malloc(n * sizeof(char *)))) {
    P_ERR("failed to allocate memory for filenames\n");
    free(names); free(idx);
    return CUTSKY_ERR_MEMORY;
  }
  **list = names;
  for (int i = 1; i < n; i++) (*list)[i] = names + idx[i];
  *num = n;

  free(idx);
  return 0;
}

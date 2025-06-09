/*******************************************************************************
* read_ascii.c: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>     /* IWYU pragma: keep */

/*============================================================================*\
                          Function for reading chunks
\*============================================================================*/

/******************************************************************************
Function `readchunk_ascii`:
  Read a chunk from an ASCII file.
Arguments:
  * `ifile`:    interface for file reading.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int readchunk_ascii(IFILE *ifile) {
  /* Copy unused bytes to the beginning of the buffer if needed. */
  if (ifile->used)
    memmove(ifile->chunk, ifile->chunk + ifile->used, ifile->rest);

  /* Check the number of bytes to be read and enlarge the chunk if needed. */
  size_t nread = ifile->size - ifile->rest;
  if (!nread) {
    if (CUTSKY_MAX_CHUNK / 2 < ifile->size) {
      P_ERR("line of ASCII file exceeding %d bytes\n", CUTSKY_MAX_CHUNK);
      return CUTSKY_ERR_FILE;
    }
    const size_t size = ifile->size * 2;

    char *tmp = realloc(ifile->chunk, size * sizeof(char));
    if (!tmp) {
      P_ERR("failed to allocate memory for reading file by chunk\n");
      return CUTSKY_ERR_MEMORY;
    }

    nread = ifile->size;
    ifile->chunk = tmp;
    ifile->size = size;
  }

  nread = fread(ifile->chunk + ifile->rest, sizeof(char), nread, ifile->fp);
  ifile->used = 0;
  ifile->rest += nread;

  if (!nread) return 0;         /* end-of-file */

  /* Append '\n' to the last line in case the file is not ending with '\n'. */
  if (ifile->rest < ifile->size) ifile->chunk[ifile->rest] = '\n';
  return 0;
}


/*============================================================================*\
                          Interfaces for file reading
\*============================================================================*/

/******************************************************************************
Function `input_init`:
  Initialise the interface for file reading.
Return:
  Address of the interface.
******************************************************************************/
IFILE *input_init(void) {
  IFILE *ifile = malloc(sizeof *ifile);
  if (!ifile) {
    P_ERR("failed to initialize the interface for file reading\n");
    return NULL;
  }

  ifile->fname = ifile->line = NULL;
  ifile->fp = NULL;
  ifile->lines = NULL;
  ifile->used = ifile->rest = 0;
  ifile->size = CUTSKY_FILE_CHUNK;
  ifile->nline = ifile->maxline = 0;
  ifile->chunk = calloc(ifile->size, sizeof(char));
  if (!ifile->chunk) {
    P_ERR("failed to allocate memory for reading file by chunk\n");
    free(ifile);
    return NULL;
  }

  return ifile;
}

/******************************************************************************
Function `input_destroy`:
  Deconstruct the interface for file reading.
Arguments:
  * `ifile`:    structure for file reading.
******************************************************************************/
void input_destroy(IFILE *ifile) {
  if (!ifile) return;
  free(ifile->chunk);
  if (ifile->lines) free(ifile->lines);
  if (ifile->fp) {
    if (fclose(ifile->fp))
      P_WRN("failed to close file: `%s'\n", ifile->fname);
  }
  free(ifile);
}

/******************************************************************************
Function `input_newfile`:
  Open a new file for reading.
Arguments:
  * `ifile`:    interface for file reading;
  * `fname`:    name of the file to be read from.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int input_newfile(IFILE *ifile, const char *fname) {
  if (!ifile) {
    P_ERR("the interface for file reading is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (!fname || !(*fname)) {
    P_ERR("invalid input file name\n");
    return CUTSKY_ERR_ARG;
  }

  /* Close the previous file and open the current one. */
  if (ifile->fp && fclose(ifile->fp))
    P_WRN("failed to close file: `%s'\n", ifile->fname);

  ifile->fp = NULL;
  if (!(ifile->fp = fopen(fname, "r"))) {
    P_ERR("failed to open file for reading: `%s'\n", fname);
    return CUTSKY_ERR_FILE;
  }

  /* Initialise the first chunk. */
  ifile->used = 0;
  ifile->rest = fread(ifile->chunk, sizeof(char), ifile->size, ifile->fp);
  if (!ifile->rest) {
    P_ERR("empty file: `%s'\n", fname);
    return CUTSKY_ERR_FILE;
  }
  /* Append '\n' to the last line. */
  if (ifile->rest < ifile->size) ifile->chunk[ifile->rest] = '\n';

  ifile->fname = fname;
  return 0;
}

/******************************************************************************
Function `input_readline`:
  Read a record (line) from the input file.
Arguments:
  * `ifile`:    interface for file reading.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int input_readline(IFILE *ifile) {
  if (!ifile) {
    P_ERR("the interface for file reading is not initialized\n");
    return CUTSKY_ERR_ARG;
  }

  char *begin = ifile->chunk + ifile->used;
  char *end;

  /* Read the file if no line is found. */
  while (!(end = memchr(begin, '\n', ifile->rest))) {
    if (ifile->used + ifile->rest < ifile->size) {      /* end-of-file */
      if (!feof(ifile->fp)) {
        P_ERR("unexpected end of file: `%s'\n", ifile->fname);
        return CUTSKY_ERR_FILE;
      }
      ifile->line = NULL;
      return 0;
    }

    int ecode;
    if ((ecode = readchunk_ascii(ifile))) {
      ifile->line = NULL;
      return ecode;
    }

    begin = ifile->chunk + ifile->used;
  }

  *end = '\0';          /* replace '\n' by string terminator '\0' */
  ifile->line = begin;
  const size_t len = end + 1 - begin;
  ifile->used += len;
  ifile->rest -= len;
  return 0;
}

/******************************************************************************
Function `input_readlines`:
  Read multiple records (lines) from the input file.
Arguments:
  * `ifile`:    interface for file reading;
  * `nline`:    number of lines to be read at once.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int input_readlines(IFILE *ifile, const size_t nline) {
  if (!ifile) {
    P_ERR("the interface for file reading is not initialized\n");
    return CUTSKY_ERR_ARG;
  }
  if (!nline) {
    P_ERR("number of lines read from file should be positive\n");
    return CUTSKY_ERR_ARG;
  }

  /* Allocate memory if needed. */
  if (ifile->maxline < nline) {
    ifile->maxline = nline;
    size_t *tmp = realloc(ifile->lines, ifile->maxline * sizeof(size_t));
    if (!tmp) {
      P_ERR("failed to allocate memory for file lines\n");
      return CUTSKY_ERR_MEMORY;
    }
    ifile->lines = tmp;
  }

  char *begin = ifile->chunk + ifile->used;
  char *end;
  size_t processed = 0;

  /* Read `nline` lines. */
  for (ifile->nline = 0; ifile->nline < nline; ifile->nline++) {
    /* Read the file if no line is found. */
    while (!(end = memchr(begin, '\n', ifile->rest - processed))) {
      if (ifile->used + ifile->rest < ifile->size) {    /* end-of-file */
        if (!feof(ifile->fp)) {
          P_ERR("unexpected end of file: `%s'\n", ifile->fname);
          return CUTSKY_ERR_FILE;
        }
        /* Process the remaining lines. */
        for (size_t i = 0; i < ifile->nline; i++)
          ifile->lines[i] += ifile->used;
        ifile->used += processed;
        ifile->rest -= processed;
        return 0;
      }

      int ecode;
      if ((ecode = readchunk_ascii(ifile))) return ecode;
      begin = ifile->chunk + ifile->used + processed;
    }

    *end = '\0';        /* replace '\n' by string terminator '\0' */

    ifile->lines[ifile->nline] = processed;
    const size_t len = end + 1 - begin;
    begin += len;
    processed += len;
  }

  for (size_t i = 0; i < ifile->nline; i++) ifile->lines[i] += ifile->used;
  ifile->used += processed;
  ifile->rest -= processed;
  return 0;
}
